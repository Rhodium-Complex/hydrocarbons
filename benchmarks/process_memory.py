"""Sample simultaneous RSS across the benchmark process and its children."""
from threading import Event, Thread

import psutil


class ProcessTreeMemory:
    """Sample every 20 ms; RSS includes shared pages in each process."""

    def __init__(self):
        self.peak_bytes = 0
        self._process = psutil.Process()
        self._stop = Event()
        self._thread = Thread(target=self._monitor, daemon=True)

    def sample(self):
        total = 0
        try:
            processes = [self._process, *self._process.children(recursive=True)]
        except psutil.NoSuchProcess:
            return
        for process in processes:
            try:
                total += process.memory_info().rss
            except psutil.NoSuchProcess:
                pass
        self.peak_bytes = max(self.peak_bytes, total)

    def _monitor(self):
        while not self._stop.wait(0.02):
            self.sample()

    def __enter__(self):
        self.sample()
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._thread.join()
        self.sample()
