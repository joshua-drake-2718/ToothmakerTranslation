import numpy as np
from typing import Self

# override of range
class range:
    def __init__(self, start, stop=None, step=1, /):
        if step == 0:
            raise ValueError("range() arg 3 must not be zero")
        elif stop is None:
            stop = start
            start = 1
        self.start = start
        self.stop = stop
        self.step = step

    def __len__(self):
        return max(0, (self.stop - self.start + abs(self.step)) // self.step)

    def __getitem__(self, index):
        if not 0 <= index < len(self):
            raise IndexError("range object index out of range")
        return index*self.step + self.start

    def __iter__(self):
        value = self.start
        while value <= self.stop if self.step > 0 else value >= self.stop:
            yield value
            value += self.step

# transform index or tuple of indices to be one-based
def transform_indices(indices):
    if isinstance(indices, tuple):
        # iterate over tuple
        return tuple(transform_index(index) for index in indices)
    else:
        # single index
        return transform_index(indices)

# transform single index to be one-based
def transform_index(index):
    if isinstance(index, (int, np.integer)):
        # index is integer (Python int or numpy scalar)
        return int(index) - 1
    elif isinstance(index, slice):
        # index is tuple
        start = index.start
        if isinstance(start, int):
            start -= 1
        return slice(start, index.stop, index.step)
    else:
        return transform_generator(index)

# transform generator function to be one-based
def transform_generator(generator):
    for index in generator:
        yield transform_index(index)


# Base class providing 1-based indexing semantics over a numpy ndarray.
# Subclasses set _dtype.
class _array_base:
    _dtype = None

    def __init__(self, shape=None, buffer=None):
        if shape is not None:
            self.inner = np.ndarray(shape, buffer=buffer, dtype=self._dtype)
        elif buffer is not None:
            self.inner = np.array(buffer, dtype=self._dtype)

    def __getitem__(self, indices):
        value = self.inner[transform_indices(indices)]
        if isinstance(value, np.ndarray):
            wrapped = type(self).__new__(type(self))
            wrapped.inner = value
            return wrapped
        return value

    def __setitem__(self, indices, value):
        if isinstance(value, _array_base):
            value = value.inner
        self.inner[transform_indices(indices)] = value

    def copy(self) -> Self:
        wrapped = type(self).__new__(type(self))
        wrapped.inner = self.inner.copy()
        return wrapped

    def __array__(self, dtype=None, copy=None):
        if copy:
            return self.inner.astype(dtype) if dtype else self.inner.copy()
        return self.inner.astype(dtype) if dtype else self.inner

    def __getattr__(self, name):
        # Forward unknown attributes (e.g. .sum, .max, .shape) to the inner ndarray.
        # __getattr__ is only called when normal attribute lookup fails, so it
        # never shadows real attributes like `.inner` on initialised objects.
        return getattr(self.__dict__['inner'], name)

    # Numeric protocol — delegate to numpy and return a plain ndarray.
    # The wrapper is a 1-based indexing aid; once numpy math is involved we
    # drop back to plain ndarrays. Assignments back through __setitem__ accept
    # ndarrays via the else branch.
    def _binop(self, other, op):
        if isinstance(other, _array_base):
            other = other.inner
        return op(self.inner, other)

    def __add__(self, other):     return self._binop(other, np.add)
    def __radd__(self, other):    return self._binop(other, lambda a, b: np.add(b, a))
    def __sub__(self, other):     return self._binop(other, np.subtract)
    def __rsub__(self, other):    return self._binop(other, lambda a, b: np.subtract(b, a))
    def __mul__(self, other):     return self._binop(other, np.multiply)
    def __rmul__(self, other):    return self._binop(other, lambda a, b: np.multiply(b, a))
    def __truediv__(self, other): return self._binop(other, np.true_divide)
    def __rtruediv__(self, other):return self._binop(other, lambda a, b: np.true_divide(b, a))
    def __neg__(self):            return -self.inner
    def __abs__(self):            return np.abs(self.inner)
    def __pos__(self):            return +self.inner

    def __iadd__(self, other):
        if isinstance(other, _array_base):
            other = other.inner
        self.inner += other
        return self
    def __isub__(self, other):
        if isinstance(other, _array_base):
            other = other.inner
        self.inner -= other
        return self
    def __imul__(self, other):
        if isinstance(other, _array_base):
            other = other.inner
        self.inner *= other
        return self
    def __itruediv__(self, other):
        if isinstance(other, _array_base):
            other = other.inner
        self.inner /= other
        return self


class int_array(_array_base):
    _dtype = np.int32

class float_array(_array_base):
    _dtype = np.float64

class bool_array(_array_base):
    _dtype = np.bool_

class str_array(_array_base):
    _dtype = np.str_
