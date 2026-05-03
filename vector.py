import math
import numpy as np

def magnitude(x: float, y: float, z: float) -> float:
    return math.sqrt(x*x + y*y + z*z)

def a_magnitude(v) -> float:
    if isinstance(v, np.ndarray):
        return magnitude(v[0], v[1], v[2])
    return magnitude(v[0], v[1], v[2])

def distance_between(
    x1: float, y1: float, z1: float,
    x2: float, y2: float, z2: float,
) -> float:
    return magnitude(x1-x2, y1-y2, z1-z2)

def a_distance_between(u, v) -> float:
    return distance_between(u[0], u[1], u[2], v[0], v[1], v[2])

def cross_product(
    x1: float, y1: float, z1: float,
    x2: float, y2: float, z2: float,
):
    return np.array([
        y1 * z2 - z1 * y2,
        z1 * x2 - x1 * z2,
        x1 * y2 - y1 * x2,
    ])

def a_cross_product(u, v) -> np.ndarray:
    return cross_product(u[0], u[1], u[2], v[0], v[1], v[2])
