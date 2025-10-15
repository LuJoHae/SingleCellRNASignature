from importlib import import_module
from typing import Any

__all__ = ["adata", "raw"]

def __getattr__(name: str) -> Any:
    if name in __all__:
        mod = import_module(__name__ + "." + name)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
