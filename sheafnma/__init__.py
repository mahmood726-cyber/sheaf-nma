"""SheafNMA — cellular-sheaf inconsistency localization for NMA."""
__version__ = "0.2.0"

from .core import (
    build_coboundary,
    sheaf_laplacian,
    edge_residuals,
    gii,
)
from .io import load_csv, save_csv
from .validate import (
    validate_network,
    is_connected,
    connected_components,
    NetworkValidationError,
)
from .inference import gii_permutation_test
from .comparators import (
    design_by_treatment_chi2,
    bucher_loop_closure,
    enumerate_independent_loops,
)

__all__ = [
    "__version__",
    "build_coboundary",
    "sheaf_laplacian",
    "edge_residuals",
    "gii",
    "load_csv",
    "save_csv",
    "validate_network",
    "is_connected",
    "connected_components",
    "NetworkValidationError",
    "gii_permutation_test",
    "design_by_treatment_chi2",
    "bucher_loop_closure",
    "enumerate_independent_loops",
]
