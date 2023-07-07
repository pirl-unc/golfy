import importlib.metadata
from .initialization import init
from .optimization import optimize
from .solution import Solution
from .validity import is_valid, count_violations, violations_per_replicate

__version__ = importlib.metadata.version("golfy")

__all__ = [
    "init",
    "optimize",
    "count_violations",
    "is_valid",
    "violations_per_replicate",
]
