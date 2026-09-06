"""A transient SVG progress line that stays out of redirected logs."""
from collections.abc import Callable
import shutil
import sys
import time


class SvgProgressBar:
    """Render SVG page progress while allowing callers to supply step logging."""

    def __init__(self, stream=None, print_step: Callable[[object], None] | None = None):
        self._stream = sys.stdout if stream is None else stream
        self._print_step = print_step
        self._enabled = self._stream.isatty()
        self._width = 0
        self._last_update = 0.0

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.clear()

    def clear(self):
        if self._width:
            self._stream.write("\r" + " " * self._width + "\r")
            self._stream.flush()
            self._width = 0

    def update(self, completed, submitted):
        if not self._enabled:
            return
        if completed == submitted:
            self.clear()
            return
        now = time.monotonic()
        if self._width and now - self._last_update < 0.1:
            return
        self._last_update = now
        filled = 20 * completed // max(1, submitted)
        line = f"SVG [{'#' * filled}{'-' * (20 - filled)}] {completed}/{submitted} submitted pages"
        line = line[:max(1, shutil.get_terminal_size().columns - 1)]
        self.clear()
        self._stream.write("\r" + line)
        self._stream.flush()
        self._width = len(line)

    def log_step(self, result):
        self.clear()
        if self._print_step is not None:
            self._print_step(result)
