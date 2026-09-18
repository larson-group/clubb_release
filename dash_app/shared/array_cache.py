"""A small LRU with a hard byte budget for plot-ready arrays."""

from collections import OrderedDict
import sys

import numpy as np


def value_bytes(value):
    """Conservatively count arrays and their surrounding Python containers."""
    if isinstance(value, np.ndarray):
        return value.nbytes + 128
    if isinstance(value, dict):
        return sys.getsizeof(value) + sum(value_bytes(k) + value_bytes(v) for k, v in value.items())
    if isinstance(value, (list, tuple)):
        return sys.getsizeof(value) + sum(value_bytes(v) for v in value)
    return sys.getsizeof(value)


class ArrayCache(OrderedDict):
    def __init__(self, max_bytes, max_entries=256):
        super().__init__()
        self.max_bytes = max_bytes
        self.max_entries = max_entries
        self.total_bytes = 0
        self._sizes = {}

    def __setitem__(self, key, value):
        if key in self:
            del self[key]
        size = value_bytes(value)
        if size > self.max_bytes:
            return  # A one-off large image must not evict every useful profile.
        super().__setitem__(key, value)
        self._sizes[key] = size
        self.total_bytes += size
        while self.total_bytes > self.max_bytes or len(self) > self.max_entries:
            self.popitem(last=False)

    def __delitem__(self, key):
        super().__delitem__(key)
        self.total_bytes -= self._sizes.pop(key)

    def popitem(self, last=True):
        if not self:
            raise KeyError("cache is empty")
        key = next(reversed(self)) if last else next(iter(self))
        value = self[key]
        del self[key]
        return key, value

    def clear(self):
        super().clear()
        self._sizes.clear()
        self.total_bytes = 0
