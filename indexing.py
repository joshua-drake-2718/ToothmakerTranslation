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
         return (self.stop - self.start + abs(self.step)) // self.step

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
        return tuple(transform_index(index for index in indices))
    else:
        # single index
        return transform_index(indices)

# transform single index to be one-based
def transform_index(index):
    if isinstance(index, int):
        # index is integer
        return index - 1
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

# wrapper of numpy.NDArray of int32s
class int_array():
    def __init__(self, shape=None, buffer=None):
        if not shape is None:
            self.inner = np.ndarray(shape, buffer=buffer, dtype=np.int32)
        elif not buffer is None:
            self.inner = np.array(buffer, dtype=np.int32)
    
    def __getitem__(self, indices) -> int | Self:
        value = self.inner[transform_indices(indices)]
        if isinstance(value, np.ndarray):
            array = int_array()
            array.inner = value
            return array
        return value
    
    def __setitem__(self, indices, value: int | Self):
        if isinstance(value, int_array):
            self.inner[transform_indices(indices)] = value.inner
        else:
            self.inner[transform_indices(indices)] = value

    def copy(self) -> Self:
        copy = int_array()
        copy.inner = self.inner.copy()
        return copy

# wrapper of numpy.NDArray of float64s
class float_array():
    def __init__(self, shape=None, buffer=None):
        if not shape is None:
            self.inner = np.ndarray(shape, buffer=buffer, dtype=np.float64)
        elif not buffer is None:
            self.inner = np.array(buffer, dtype=np.float64)
    
    def __getitem__(self, indices) -> float | Self:
        value = self.inner[transform_indices(indices)]
        if isinstance(value, np.ndarray):
            array = float_array()
            array.inner = value
            return array
        return value
    
    def __setitem__(self, indices, value: float | Self):
        if isinstance(value, float_array):
            self.inner[transform_indices(indices)] = value.inner
        else:
            self.inner[transform_indices(indices)] = value

    def copy(self) -> Self:
        x = float_array()
        x.inner = self.inner.copy()
        return x

# wrapper of numpy.NDArray of bools
class bool_array():
    def __init__(self, shape=None, buffer=None):
        if not shape is None:
            self.inner = np.ndarray(shape, buffer=buffer, dtype=np.bool)
        elif not buffer is None:
            self.inner = np.array(buffer, dtype=np.bool)
    
    def __getitem__(self, indices) -> bool | Self:
        value = self.inner[transform_indices(indices)]
        if isinstance(value, np.ndarray):
            array = bool_array()
            array.inner = value
            return array
        return value
    
    def __setitem__(self, indices, value: bool | Self):
        if isinstance(value, bool_array):
            self.inner[transform_indices(indices)] = value.inner
        else:
            self.inner[transform_indices(indices)] = value

    def copy(self) -> Self:
        x = bool_array()
        x.inner = self.inner.copy()
        return x

# wrapper of numpy.NDArray of strings
class str_array():
    def __init__(self, shape=None, buffer=None):
        if not shape is None:
            self.inner = np.ndarray(shape, buffer=buffer, dtype=np.str_)
        elif not buffer is None:
            self.inner = np.array(buffer, dtype=np.str_)
    
    def __getitem__(self, indices) -> str | Self:
        value = self.inner[transform_indices(indices)]
        if isinstance(value, np.ndarray):
            array = str_array()
            array.inner = value
            return array
        return value
    
    def __setitem__(self, indices, value: str | Self):
        if isinstance(value, str_array):
            self.inner[transform_indices(indices)] = value.inner
        else:
            self.inner[transform_indices(indices)] = value

    def copy(self) -> Self:
        x = str_array()
        x.inner = self.inner.copy()
        return x