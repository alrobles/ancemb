"""Small, versioned Python adapter used by ancemb through reticulate."""

from __future__ import annotations

import numpy as np


__version__ = "0.1.0"


def load_array(path: str, key: str | None = None):
    """Load one NumPy array from NPY or NPZ without changing its values."""
    loaded = np.load(path, allow_pickle=False)
    if isinstance(loaded, np.lib.npyio.NpzFile):
        names = list(loaded.files)
        selected = key or (names[0] if len(names) == 1 else None)
        if selected is None:
            raise ValueError("an NPZ with multiple arrays requires key")
        if selected not in loaded.files:
            raise KeyError(selected)
        return loaded[selected]
    return loaded
