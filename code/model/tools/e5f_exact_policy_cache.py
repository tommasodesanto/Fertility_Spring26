"""Opt-in, in-memory exact cache for repeated dated household policy solves."""
from __future__ import annotations

from collections import OrderedDict
from contextlib import contextmanager
import hashlib
import inspect
import pickle
import struct
import threading
from types import SimpleNamespace

import numpy as np

DEFAULT_MAX_BYTES = 6 * 1024**3


def _field(digest, payload):
    payload = bytes(payload)
    digest.update(struct.pack("!Q", len(payload)))
    digest.update(payload)


def _update(digest, value, active):
    if value is None:
        digest.update(b"N"); return
    if isinstance(value, np.generic):
        digest.update(b"G"); _field(digest, pickle.dumps(value.dtype, protocol=5))
        _field(digest, np.asarray(value).tobytes(order="C")); return
    if isinstance(value, bool):
        digest.update(b"B1" if value else b"B0"); return
    if isinstance(value, int) and not isinstance(value, np.integer):
        digest.update(b"I"); _field(digest, str(value).encode()); return
    if isinstance(value, float):
        digest.update(b"F"); digest.update(struct.pack("!d", value)); return
    if isinstance(value, complex):
        digest.update(b"C"); digest.update(struct.pack("!dd", value.real, value.imag)); return
    if isinstance(value, str):
        digest.update(b"S"); _field(digest, value.encode("utf-8")); return
    if isinstance(value, (bytes, bytearray, memoryview)):
        digest.update(b"Y"); _field(digest, bytes(value)); return
    if isinstance(value, np.ndarray):
        if type(value) is not np.ndarray:
            raise TypeError("array subclasses may carry additional solver state")
        digest.update(b"A"); _field(digest, pickle.dumps(value.dtype, protocol=5))
        _field(digest, pickle.dumps(tuple(int(x) for x in value.shape), protocol=5))
        if value.dtype.hasobject:
            _field(digest, pickle.dumps(value, protocol=5))
        else:
            _field(digest, np.ascontiguousarray(value).tobytes(order="C"))
        return
    identity = id(value)
    if identity in active:
        raise TypeError("cyclic policy argument state is not cacheable")
    active.add(identity)
    try:
        if isinstance(value, (tuple, list)):
            digest.update(b"T" if isinstance(value, tuple) else b"L")
            digest.update(struct.pack("!Q", len(value)))
            for item in value: _update(digest, item, active)
            return
        if isinstance(value, dict):
            digest.update(b"D")
            entries = []
            for key, item in value.items():
                key_digest = hashlib.sha256(); _update(key_digest, key, set())
                entries.append((key_digest.digest(), key, item))
            digest.update(struct.pack("!Q", len(entries)))
            for _, key, item in sorted(entries, key=lambda row: row[0]):
                _update(digest, key, active); _update(digest, item, active)
            return
        if isinstance(value, (set, frozenset)):
            digest.update(b"E" if isinstance(value, set) else b"R")
            entries = []
            for item in value:
                item_digest = hashlib.sha256(); _update(item_digest, item, set())
                entries.append((item_digest.digest(), item))
            digest.update(struct.pack("!Q", len(entries)))
            for _, item in sorted(entries, key=lambda row: row[0]): _update(digest, item, active)
            return
        if type(value) is SimpleNamespace:
            digest.update(b"O")
            _field(digest, f"{type(value).__module__}.{type(value).__qualname__}".encode())
            _update(digest, vars(value), active)
            return
        raise TypeError("unsupported policy argument state is not cacheable")
    finally:
        active.remove(identity)


def exact_call_key(function, args, kwargs):
    """Hash the complete bound call, independent of keyword insertion order."""
    bound = inspect.signature(function).bind(*args, **kwargs)
    bound.apply_defaults()
    digest = hashlib.sha256()
    _update(digest, dict(bound.arguments), set())
    return digest.digest()


class PolicyCacheStats:
    def __init__(self, maximum_bytes):
        self.maximum_bytes = int(maximum_bytes)
        self.actual_solves = 0
        self.hits = 0
        self.misses = 0
        self.bytes = 0
        self.evictions = 0
        self.oversize_bypasses = 0
        self.serialization_bypasses = 0
        self.exceptions = 0

    def snapshot(self):
        return dict(maximum_bytes=self.maximum_bytes, actual_solves=self.actual_solves,
                    hits=self.hits, misses=self.misses, bytes=self.bytes,
                    evictions=self.evictions, oversize_bypasses=self.oversize_bypasses,
                    serialization_bypasses=self.serialization_bypasses,
                    exceptions=self.exceptions)


@contextmanager
def policy_cache(module, max_bytes=DEFAULT_MAX_BYTES):
    """Temporarily cache ``module.solve_date_policy`` by its complete call state."""
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes < 0:
        raise ValueError("max_bytes must be a nonnegative integer")
    original = module.solve_date_policy
    entries = OrderedDict()
    stats = PolicyCacheStats(max_bytes)
    lock = threading.RLock()

    def wrapped(*args, **kwargs):
        try:
            key = exact_call_key(original, args, kwargs)
        except Exception:
            with lock:
                stats.serialization_bypasses += 1
                stats.misses += 1
                stats.actual_solves += 1
            try:
                return original(*args, **kwargs)
            except Exception:
                with lock: stats.exceptions += 1
                raise
        with lock:
            payload = entries.pop(key, None)
            if payload is not None:
                entries[key] = payload
                stats.hits += 1
                return pickle.loads(payload)
            stats.misses += 1
            stats.actual_solves += 1
        try:
            result = original(*args, **kwargs)
        except Exception:
            with lock: stats.exceptions += 1
            raise
        try:
            payload = pickle.dumps(result, protocol=5)
        except Exception:
            with lock: stats.serialization_bypasses += 1
            return result
        size = len(payload)
        with lock:
            if size > max_bytes:
                stats.oversize_bypasses += 1
                return result
            while entries and stats.bytes + size > max_bytes:
                _, evicted = entries.popitem(last=False)
                stats.bytes -= len(evicted); stats.evictions += 1
            entries[key] = payload
            stats.bytes += size
        return result

    module.solve_date_policy = wrapped
    try:
        yield stats
    finally:
        module.solve_date_policy = original
