"""DMS - Debris disks around Main-sequence Stars.

Thermal emission, scattered light and polarisation of optically thin debris disks
(Kim et al. 2018, A&A 618, A38). Fortran core, Python interface via numpy.f2py.

    import dms
    res = dms.Model.from_file("examples/debris_ring.nml").run()
"""

from .model import Model, Result, DMSError, SizeGridWarning, default_data_dir
from .namelist import read_namelist, write_namelist
from . import io, plot, stellar, batch, emc

__version__ = "2.0.0"

__all__ = ["Model", "Result", "DMSError", "SizeGridWarning", "default_data_dir", "read_namelist", "write_namelist",
           "io", "plot", "stellar", "batch", "emc", "__version__"]
