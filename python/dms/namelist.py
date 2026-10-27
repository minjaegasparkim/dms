"""Reading and writing DMS input files (Fortran namelists) as Python dictionaries.

The input file has the groups &run, &star, &disk, &image, &wavelengths, &dust and &scattering
(see docs/input_reference.md). Here a namelist is a dict of dicts::

    {"star": {"teff": 5800.0, "radius": 1.0}, "disk": {"r_in": [50.0], ...}, ...}

Only the subset of the namelist syntax used by DMS input files is supported by the reader:
``name = value[, value ...]`` entries, '!' comments, strings in quotes, logicals .true./.false.,
and several entries per line separated by commas.
"""

import re

GROUPS = ("run", "star", "disk", "image", "wavelengths", "dust", "scattering")

_token = re.compile(r"""'[^']*'|"[^"]*"|[^,\s]+""")


def _fmt(v):
    import numbers
    if isinstance(v, bool) or getattr(getattr(v, "dtype", None), "kind", "") == "b":
        return ".true." if bool(v) else ".false."
    if isinstance(v, str):
        return "'" + v.replace("'", "''") + "'"
    if isinstance(v, numbers.Integral):
        return str(int(v))
    if isinstance(v, numbers.Real):            # also numpy floats (repr would give 'np.float64(...)')
        return repr(float(v))
    raise TypeError(f"cannot write value {v!r} to a namelist")


def format_value(value):
    """Namelist representation of a scalar or a list."""
    if isinstance(value, (list, tuple)) or (hasattr(value, "__len__") and not isinstance(value, str)):
        return ", ".join(_fmt(v) for v in value)
    return _fmt(value)


def write_namelist(groups, path):
    """Write a dict of groups to a namelist file."""
    with open(path, "w") as f:
        for g in GROUPS:
            entries = groups.get(g)
            if not entries:
                continue
            f.write(f"&{g}\n")
            for k, v in entries.items():
                f.write(f"  {k} = {format_value(v)}\n")
            f.write("/\n\n")
        unknown = set(groups) - set(GROUPS)
        if unknown:
            raise KeyError(f"unknown namelist group(s): {', '.join(sorted(unknown))}")


def _parse_value(tok):
    t = tok.strip()
    if t[:1] in "'\"":
        return t[1:-1]
    tl = t.lower()
    if tl in (".true.", "t", ".t."):
        return True
    if tl in (".false.", "f", ".f."):
        return False
    try:
        if re.fullmatch(r"[+-]?\d+", t):
            return int(t)
        return float(tl.replace("d", "e"))
    except ValueError:
        return t


def _strip_comment(line):
    out, quote = [], None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def read_namelist(path):
    """Read a DMS input file into a dict of groups (lists for multi-valued entries)."""
    groups, current, text = {}, None, []
    with open(path) as f:
        for raw in f:
            line = _strip_comment(raw).strip()
            if not line:
                continue
            if line.startswith("&"):
                current = line[1:].split()[0].lower()
                groups[current] = {}
                text = []
                line = line[1 + len(current):].strip()
            if current is None:
                continue
            end = line.endswith("/")
            if end:
                line = line[:-1]
            text.append(line)
            if end:
                body = " ".join(text)
                # split at "name =" boundaries
                parts = re.split(r"(?:^|,|\s)\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", " " + body)
                for name, val in zip(parts[1::2], parts[2::2]):
                    vals = []
                    for t in _token.findall(val):
                        rep = re.fullmatch(r"(\d+)\*(.+)", t)          # Fortran repeat count n*value
                        if rep and t[:1] not in "'\"":
                            vals += [_parse_value(rep.group(2))] * int(rep.group(1))
                        else:
                            vals.append(_parse_value(t))
                    groups[current][name.lower()] = vals[0] if len(vals) == 1 else vals
                current = None
    return groups
