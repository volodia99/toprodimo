__all__ = ["F", "FArray2D"]

from typing import TypeAlias, TypeVar

import numpy as np

f32 = np.float32
f64 = np.float64
F = TypeVar("F", f32, f64)

FArray2D: TypeAlias = np.ndarray[tuple[int, int], np.dtype[F]]
