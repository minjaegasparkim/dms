"""Running many DMS models: input files in a queue folder, several in parallel.

Typical use (from Python or a notebook)::

    import dms
    from dms import batch

    base = dms.Model.from_file("examples/debris_ring.nml")
    jobs = batch.grid(base, {"disk.inclination": [0, 45, 70],
                             "dust.q": [3.0, 3.5, 4.0]})
    batch.write_queue(jobs, "queue")                 # one input file per model
    batch.run("queue", "runs", workers=8)            # run them all (returns when the queue is empty)
    table = batch.collect("runs")                    # one row per model

From a terminal (the same functions)::

    python -m dms.batch run queue runs -j 8          # run everything that is in queue/
    python -m dms.batch run queue runs -j 8 --watch  # keep running and pick up new files
    python -m dms.batch status queue

Queue folder layout: new input files are put into ``queue/`` (any number, at any time; with
--watch they are picked up as they arrive). A file being computed is moved to
``queue/running/``, then to ``queue/done/`` or ``queue/failed/``. Several runner processes (even
on different machines sharing the folder) can work on the same queue: a file is claimed by an
atomic rename. A file ``queue/STOP`` ends a --watch run after the running models.

Each model writes into ``runs/<input name>/``: the output files of DMS, the input file as run
(``input.nml``), the screen output (``dms.log``) and ``result.json`` (wavelengths, SED, degree
of polarisation against scattering angle, disk-integrated Q_phi / I if images were made).

Optical-property caches: with ``optics_cache_name = 'auto'`` (the default set by write_queue) the
cache name is derived from the dust, wavelength and angle parameters, so models that differ only
in geometry share one Mie calculation. The runner never starts two models that would compute
the same cache at the same time; the second one waits until the first has saved it. A cache file
can be 100-200 MB (all optical-data wavelengths, all sizes, all angles), so the cache folder is
kept below a size limit (``cache_max_gb``, default 10 GB): after each model the least recently
used files are deleted (never one that a running model needs). ``cache_max_gb=0`` switches the
caches off.

After an interruption (power, full disk): ``python -m dms.batch requeue queue --running --runs runs``
puts back the models that were running and those marked done whose results are incomplete.
"""

import argparse
import concurrent.futures as cf
import copy
import hashlib
import itertools
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

from .namelist import GROUPS, read_namelist, write_namelist

_PKG = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_PKG, "..", ".."))


def default_executable():
    """$DMS_EXE, else bin/dms of the repository."""
    return os.environ.get("DMS_EXE") or os.path.join(_REPO, "bin", "dms")


# ------------------------------------------------------------------------------------------
# making input files
# ------------------------------------------------------------------------------------------

def _set(params, key, value):
    group, name = key.split(".", 1)
    if group not in GROUPS:
        raise KeyError(f"unknown group in {key!r}")
    params.setdefault(group, {})[name] = value


def _fmt_value(v):
    if isinstance(v, (list, tuple)):
        return "-".join(_fmt_value(x) for x in v)
    if isinstance(v, float):
        return f"{v:g}"
    s = str(v)
    return "".join(ch if ch.isalnum() or ch in ".-" else "_" for ch in os.path.basename(s))


def grid(base, axes, name=None, linked=None):
    """All combinations of the parameter values in ``axes``.

    base:   dms.Model (or dict of groups) with the common parameters
    axes:   {"group.name": [values, ...], ...}
    linked: optional {"group.name": {"other.name": [...]}}: parameters that change together with
            an axis (same length as that axis), e.g. a_min that belongs to each star
    name:   function(dict of axis values) -> job name; default "<name>=<value>_..."
    Returns a list of (job name, params dict).
    """
    params0 = base.params if hasattr(base, "params") else base
    keys = list(axes)
    jobs = []
    for combo in itertools.product(*(range(len(axes[k])) for k in keys)):
        p = copy.deepcopy(params0)
        vals = {}
        for k, i in zip(keys, combo):
            v = axes[k][i]
            _set(p, k, v)
            vals[k] = v
            for lk, lv in (linked or {}).get(k, {}).items():
                _set(p, lk, lv[i])
        if name:
            jn = name(vals)
        else:
            jn = "_".join(f"{k.split('.', 1)[1]}={_fmt_value(v)}" for k, v in vals.items()) or "model"
        jobs.append((jn, p))
    names = [j[0] for j in jobs]
    if len(set(names)) != len(names):
        raise ValueError("job names are not unique; give name=...")
    return jobs


def write_queue(jobs, queue_dir, overwrite=False, auto_cache=True):
    """Write one input file per job into queue_dir.

    jobs: [(name, params or Model), ...] or [(name, params, meta), ...]; meta is a dict of labels
    (e.g. {"star": "A", "f_ice": 0.5}) that is kept with the results (meta.json, collect()).
    """
    os.makedirs(queue_dir, exist_ok=True)
    written = []
    for job in jobs:
        jn, p = job[0], job[1]
        meta = job[2] if len(job) > 2 else None
        p = copy.deepcopy(p.params if hasattr(p, "params") else p)
        if auto_cache:
            p.setdefault("run", {}).setdefault("optics_cache_name", "auto")
        path = os.path.join(queue_dir, f"{jn}.nml")
        if not overwrite and any(os.path.exists(os.path.join(queue_dir, sub, f"{jn}.nml"))
                                 for sub in ("", "running", "done")):
            raise FileExistsError(f"{jn}.nml is already in the queue (use overwrite=True)")
        tmp = path + ".part"
        write_namelist(p, tmp)
        if meta:
            with open(tmp) as f:
                body = f.read()
            with open(tmp, "w") as f:
                f.write("! meta: " + json.dumps(meta, default=str) + "\n" + body)
        os.replace(tmp, path)        # a watching runner never sees a half-written file
        written.append(path)
    return written


# ------------------------------------------------------------------------------------------
# running
# ------------------------------------------------------------------------------------------

_MIE_KEYS = ("n_components", "optical_data", "size_distribution", "a_min", "a_max", "n_sizes",
             "size_file", "emt", "emt_matrix", "emt_inclusion", "emt_volume_fraction", "emt_rule",
             "emt_output", "tdep_component", "tdep_temperatures", "tdep_files",
             "grain_shape", "dhs_fmax", "dhs_nf")


def _cache_key(p):
    """Name of the optical-property cache, from the parameters the Mie calculation depends on
    (grain radii, optical data, wavelengths, angles). Models that differ only in geometry, star,
    size-distribution index or abundances share it. DMS checks radii, wavelengths and n, k when
    it loads a cache, so a cache is never used for other dust."""
    d = p.get("dust", {})
    sub = {k: d[k] for k in _MIE_KEYS if k in d}
    sub["wavelengths"] = p.get("wavelengths", {})
    sub["n_angles"] = p.get("scattering", {}).get("n_angles", 181)
    sub["data_dir"] = p.get("run", {}).get("data_dir", "")
    h = hashlib.sha1(json.dumps(sub, sort_keys=True, default=str).encode()).hexdigest()[:16]
    return "auto_" + h


def _prepare(nml_path, run_dir, cache_dir, data_dir):
    p = read_namelist(nml_path)
    r = p.setdefault("run", {})
    r["output_dir"] = run_dir.rstrip("/") + "/"
    r["cache_dir"] = cache_dir.rstrip("/") + "/"
    r.setdefault("data_dir", data_dir)
    r["write_files"] = True
    key = None
    if str(r.get("optics_cache_name", "")).lower() == "auto":
        key = _cache_key(p)
        r["optics_cache_name"] = key
        r["load_optics"] = True
        r["save_optics"] = True
    elif r.get("save_optics") and r.get("optics_cache_name"):
        key = str(r["optics_cache_name"])
    return p, key


def read_meta(nml_path):
    """The '! meta: {...}' line of an input file written by write_queue (or {})."""
    with open(nml_path) as f:
        for line in f:
            if line.startswith("! meta:"):
                return json.loads(line[len("! meta:"):])
            if line.strip() and not line.startswith("!"):
                break
    return {}


def _run_one(nml_running, run_dir, exe, omp_threads, cache_dir, data_dir, timeout, use_cache=True):
    os.makedirs(run_dir, exist_ok=True)
    meta = read_meta(nml_running)
    if meta:
        with open(os.path.join(run_dir, "meta.json"), "w") as f:
            json.dump(meta, f)
    p, _ = _prepare(nml_running, run_dir, cache_dir, data_dir)
    if not use_cache:
        p["run"]["load_optics"] = False
        p["run"]["save_optics"] = False
    inp = os.path.join(run_dir, "input.nml")
    write_namelist(p, inp)
    env = dict(os.environ, OMP_NUM_THREADS=str(omp_threads))
    t0 = time.time()
    with open(os.path.join(run_dir, "dms.log"), "w") as log:
        try:
            rc = subprocess.run([exe, inp], stdout=log, stderr=subprocess.STDOUT, env=env,
                                timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            log.write(f"\nbatch: stopped after the time limit of {timeout} s\n")
            rc = -1
    dt = time.time() - t0
    if rc == 0:
        try:
            summarise(run_dir, p)
        except Exception as exc:                      # e.g. full disk: the model counts as failed
            with open(os.path.join(run_dir, "dms.log"), "a") as log:
                log.write(f"\nbatch: result.json not written: {exc}\n")
            rc = -3
    return rc, dt


def _cache_bytes(cache_dir):
    total = 0
    for f in os.listdir(cache_dir):
        try:
            total += os.path.getsize(os.path.join(cache_dir, f))
        except OSError:
            pass
    return total


def prune_cache(cache_dir, max_bytes, keep=()):
    """Delete the least recently used cache files until the folder is below max_bytes; files whose
    key is in keep (used by running models) are not deleted. Returns the number deleted."""
    files = []
    for f in os.listdir(cache_dir):
        path = os.path.join(cache_dir, f)
        try:
            st = os.stat(path)
        except OSError:
            continue
        files.append((st.st_mtime, st.st_size, f, path))
    total = sum(x[1] for x in files)
    n = 0
    for mtime, size, f, path in sorted(files):
        if total <= max_bytes:
            break
        if f.split(".")[0] in keep or ".tmp" in f:
            continue
        try:
            os.remove(path)
            total -= size
            n += 1
        except OSError:
            pass
    return n


def _touch_cache(cache_dir, key):
    path = os.path.join(cache_dir, key + ".miex")
    if os.path.exists(path):
        try:
            os.utime(path)                  # marks it as recently used (for prune_cache)
        except OSError:
            pass


def run(queue_dir, runs_dir, workers=None, omp_threads=1, watch=False, exe=None,
        cache_dir=None, data_dir=None, timeout=None, poll=2.0, quiet=False, cache_max_gb=10.0):
    """Run all input files in queue_dir (see the module docstring). Returns {name: status}.

    cache_max_gb: size limit of the optical-property cache folder (least recently used files are
    deleted); 0 switches the caches off."""
    exe = os.path.abspath(exe or default_executable())
    if not os.path.isfile(exe):
        raise FileNotFoundError(f"DMS executable not found: {exe} (run 'make' first)")
    workers = workers or max(1, (os.cpu_count() or 1) // max(1, omp_threads))
    cache_dir = os.path.abspath(cache_dir or os.path.join(runs_dir, "_cache"))
    if data_dir is None:
        from .model import default_data_dir
        data_dir = default_data_dir()
    data_dir = os.path.abspath(data_dir).rstrip("/") + "/"
    for sub in ("running", "done", "failed"):
        os.makedirs(os.path.join(queue_dir, sub), exist_ok=True)
    os.makedirs(runs_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)

    status, active, busy_keys = {}, {}, set()

    def say(msg):
        if not quiet:
            print(time.strftime("%H:%M:%S"), msg, flush=True)

    def pending():
        return sorted(f for f in os.listdir(queue_dir)
                      if f.endswith(".nml") and os.path.isfile(os.path.join(queue_dir, f)))

    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        while True:
            # start new jobs while workers are free
            for f in pending():
                if len(active) >= workers:
                    break
                src = os.path.join(queue_dir, f)
                try:
                    _, key = _prepare(src, "x", cache_dir, data_dir)
                    if cache_max_gb is not None and cache_max_gb <= 0:
                        key = None
                except Exception as exc:
                    say(f"failed (cannot read input): {f}: {exc}")
                    os.replace(src, os.path.join(queue_dir, "failed", f))
                    status[f[:-4]] = "failed"
                    continue
                if key and key in busy_keys and not os.path.exists(os.path.join(cache_dir, key + ".miex")):
                    continue                     # same Mie calculation is running: wait for its cache
                dst = os.path.join(queue_dir, "running", f)
                try:
                    os.replace(src, dst)         # claim (atomic; other runners skip it)
                except FileNotFoundError:
                    continue
                name = f[:-4]
                if key:
                    _touch_cache(cache_dir, key)
                fut = pool.submit(_run_one, dst, os.path.join(runs_dir, name), exe, omp_threads,
                                  cache_dir, data_dir, timeout, cache_max_gb is None or cache_max_gb > 0)
                active[fut] = (name, dst, key)
                if key:
                    busy_keys.add(key)
                say(f"start {name}")
            if not active:
                if watch and not os.path.exists(os.path.join(queue_dir, "STOP")):
                    time.sleep(poll)
                    continue
                break
            done, _ = cf.wait(list(active), timeout=poll, return_when=cf.FIRST_COMPLETED)
            for fut in done:
                name, dst, key = active.pop(fut)
                if key:
                    busy_keys.discard(key)
                try:
                    rc, dt = fut.result()
                except Exception as exc:
                    rc, dt = -2, 0.0
                    with open(os.path.join(runs_dir, name, "dms.log"), "a") as log:
                        log.write(f"\nbatch: {exc}\n")
                ok = (rc == 0)
                os.replace(dst, os.path.join(queue_dir, "done" if ok else "failed", os.path.basename(dst)))
                status[name] = "done" if ok else "failed"
                say(f"{'done  ' if ok else 'FAILED'} {name} ({dt:.1f} s)"
                    + ("" if ok else f" - see {os.path.join(runs_dir, name, 'dms.log')}"))
                if cache_max_gb and cache_max_gb > 0:
                    prune_cache(cache_dir, cache_max_gb * 1e9, keep={a[2] for a in active.values() if a[2]})
            if watch and os.path.exists(os.path.join(queue_dir, "STOP")) and not active:
                break
    n_ok = sum(s == "done" for s in status.values())
    say(f"finished: {n_ok} done, {len(status) - n_ok} failed")
    return status


def queue_status(queue_dir):
    """Number of input files waiting, running, done and failed."""
    out = {}
    for sub, label in (("", "waiting"), ("running", "running"), ("done", "done"), ("failed", "failed")):
        d = os.path.join(queue_dir, sub)
        out[label] = len([f for f in os.listdir(d) if f.endswith(".nml")]) if os.path.isdir(d) else 0
    return out


def requeue_failed(queue_dir, running=False, runs_dir=None):
    """Move the failed input files back into the queue. running=True: also the files left in
    queue/running/ (only when no runner is active). runs_dir: also the done models without a
    readable result.json. Returns the number of files moved."""
    n = 0
    subs = ["failed"] + (["running"] if running else [])
    for sub in subs:
        d = os.path.join(queue_dir, sub)
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if f.endswith(".nml"):
                os.replace(os.path.join(d, f), os.path.join(queue_dir, f))
                n += 1
    if runs_dir:
        d = os.path.join(queue_dir, "done")
        for f in (os.listdir(d) if os.path.isdir(d) else []):
            if not f.endswith(".nml"):
                continue
            res = os.path.join(runs_dir, f[:-4], "result.json")
            try:
                with open(res) as fh:
                    json.load(fh)
                continue
            except Exception:
                os.replace(os.path.join(d, f), os.path.join(queue_dir, f))
                n += 1
    return n


# ------------------------------------------------------------------------------------------
# results
# ------------------------------------------------------------------------------------------

def _ext(p):
    name = p.get("run", {}).get("output_name", "")
    return f".{name}" if name else ""


def summarise(run_dir, params=None):
    """Write run_dir/result.json from the output files of one model and return it as a dict.

    Contents: wavelength [micron]; sed (thermal, scattered, Q, U, polarised, star; Jy);
    dop_theta: degree of polarisation of the scattered light against the scattering angle
    (F_pol / F_sca per 1-degree bin; positive = perpendicular to the scattering plane; NaN where no
    dust is seen at that angle), shape [n_wavelength, 180]; and, if images were made,
    qphi_over_i = sum(Q_phi) / sum(I_scattered) over the image and pi_over_i = sum(PI) / sum(I_sca).
    """
    if params is None:
        params = read_namelist(os.path.join(run_dir, "input.nml"))
    ext = _ext(params)
    sed_file = os.path.join(run_dir, f"SED{ext}.dat")
    if os.path.exists(sed_file):              # modes 'sed' and 'both'
        sed = np.loadtxt(sed_file, comments="#", ndmin=2)
        res = {"wavelength": sed[:, 0].tolist(),
               "sed": {k: sed[:, i].tolist() for k, i in
                       (("thermal", 1), ("scattered", 2), ("Q", 3), ("U", 4), ("polarised", 5), ("star", 7))}}
    else:                                     # mode 'image': wavelengths from the input file
        res = {"wavelength": _wavelengths(params)}
    pf = os.path.join(run_dir, f"phase_function{ext}.dat")
    if os.path.exists(pf):
        ph = np.loadtxt(pf, comments="#", ndmin=2)
        fs, fp = ph[:, 1::2].T, ph[:, 2::2].T
        with np.errstate(invalid="ignore", divide="ignore"):
            dop = np.where(fs > 0, fp / fs, np.nan)
        res["theta"] = ph[:, 0].tolist()
        res["dop_theta"] = [[None if not np.isfinite(v) else float(v) for v in row] for row in dop]
    # disk-integrated polarisation from the images (FITS or ASCII)
    try:
        q, i_ = _image_sum(run_dir, "Qphi", ext), _image_sum(run_dir, "sca", ext)
        pi = _image_sum(run_dir, "pol", ext)
        if q is not None and i_ is not None:
            with np.errstate(invalid="ignore", divide="ignore"):
                res["qphi_over_i"] = (q / i_).tolist()
                if pi is not None:
                    res["pi_over_i"] = (pi / i_).tolist()
    except Exception:
        pass
    od = os.path.join(run_dir, f"optical_depth{ext}.dat")
    if os.path.exists(od):
        keys = (("fractional luminosity", "f_d"), ("mean optical depth", "tau_mean"),
                ("largest radial optical depth", "tau_radial_max"))
        with open(od) as f:
            for line in f:
                for text, key in keys:
                    if text in line and "=" in line:
                        try:
                            res[key] = float(line.rsplit("=", 1)[1].split()[0])
                        except (ValueError, IndexError):
                            pass
    with open(os.path.join(run_dir, "result.json"), "w") as f:
        json.dump(res, f)
    return res


def _wavelengths(params):
    w = params.get("wavelengths", {})
    if "lambda" in w:
        lam = w["lambda"]
        return [float(x) for x in (lam if isinstance(lam, list) else [lam])]
    lmin, lmax, n = float(w["lambda_min"]), float(w["lambda_max"]), int(w["n_lambda"])
    return np.geomspace(lmin, lmax, n).tolist() if n > 1 else [lmin]


def _image_sum(run_dir, q, ext):
    fits_name = os.path.join(run_dir, f"image_{q}{ext}.fits")
    if os.path.exists(fits_name):
        from astropy.io import fits
        with fits.open(fits_name) as h:
            d = np.array(h[0].data, dtype=float)
        if d.ndim == 2:
            d = d[None]
        return np.nansum(d, axis=(1, 2))
    files = sorted((f for f in os.listdir(run_dir) if f.startswith(f"image_{q}{ext}_#") and f.endswith(".dat")),
                   key=lambda s: int(s.rsplit("#", 1)[1].split(".")[0]))
    if files:
        return np.array([np.nansum(np.loadtxt(os.path.join(run_dir, f), comments="#")) for f in files])
    return None


def collect(runs_dir, with_params=True):
    """Results of all finished models: list of dicts (name, params, result.json contents)."""
    rows = []
    for name in sorted(os.listdir(runs_dir)):
        d = os.path.join(runs_dir, name)
        rj = os.path.join(d, "result.json")
        if not os.path.isfile(rj):
            continue
        with open(rj) as f:
            row = {"name": name, **json.load(f)}
        if with_params and os.path.isfile(os.path.join(d, "input.nml")):
            row["params"] = read_namelist(os.path.join(d, "input.nml"))
        mj = os.path.join(d, "meta.json")
        row["meta"] = json.load(open(mj)) if os.path.isfile(mj) else {}
        rows.append(row)
    return rows


def dop_at(row, angle, il=None):
    """Degree of polarisation at a scattering angle [deg] (1-degree bin) for one result row."""
    j = min(int(angle), 179)
    vals = [r[j] for r in row["dop_theta"]]
    return vals if il is None else vals[il]


# ------------------------------------------------------------------------------------------
# command line
# ------------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m dms.batch", description=__doc__.split("\n\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run", help="run the input files of a queue folder")
    r.add_argument("queue")
    r.add_argument("runs")
    r.add_argument("-j", "--workers", type=int, default=None, help="models in parallel (default: CPUs / threads)")
    r.add_argument("-t", "--threads", type=int, default=1, help="OpenMP threads per model (default 1)")
    r.add_argument("--watch", action="store_true", help="keep running and pick up new input files (stop: create queue/STOP)")
    r.add_argument("--exe", default=None, help="DMS executable (default bin/dms)")
    r.add_argument("--cache", default=None, help="optical-property cache folder (default runs/_cache)")
    r.add_argument("--data", default=None, help="optical-data folder (default data/optical_constants)")
    r.add_argument("--timeout", type=float, default=None, help="time limit per model [s]")
    r.add_argument("--cache-max-gb", type=float, default=10.0,
                   help="size limit of the cache folder in GB (least recently used files are deleted; 0: no caches)")
    s = sp.add_parser("status", help="number of waiting / running / done / failed models")
    s.add_argument("queue")
    q = sp.add_parser("requeue", help="put the failed input files back into the queue")
    q.add_argument("queue")
    q.add_argument("--running", action="store_true",
                   help="also the files left in queue/running (only when no runner is active)")
    q.add_argument("--runs", default=None, help="also the done models of this runs folder without a readable result.json")
    a = ap.parse_args(argv)
    if a.cmd == "run":
        st = run(a.queue, a.runs, workers=a.workers, omp_threads=a.threads, watch=a.watch, exe=a.exe,
                 cache_dir=a.cache, data_dir=a.data, timeout=a.timeout, cache_max_gb=a.cache_max_gb)
        return 0 if all(v == "done" for v in st.values()) else 1
    if a.cmd == "status":
        print(queue_status(a.queue))
        return 0
    if a.cmd == "requeue":
        print(f"{requeue_failed(a.queue, running=a.running, runs_dir=a.runs)} input files moved back into the queue")
        return 0


if __name__ == "__main__":
    sys.exit(main())
