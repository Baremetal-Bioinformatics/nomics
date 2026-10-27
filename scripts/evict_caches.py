#!/usr/bin/env python3
"""Drop files under the given directories from the page cache (posix_fadvise
DONTNEED), so the next run's huge-page allocations find free memory.

usage: evict_caches.py DIR...   (no root needed; only affects clean pages)"""
import os
import sys

for top in sys.argv[1:]:
    for root, _, files in os.walk(top):
        for name in files:
            path = os.path.join(root, name)
            if os.path.islink(path):
                continue
            try:
                fd = os.open(path, os.O_RDONLY)
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
                os.close(fd)
            except OSError:
                pass
