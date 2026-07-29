"""Small Qt worker primitives for cancellable-by-generation analytics."""

from __future__ import annotations

from PyQt5 import QtCore


class TaskSignals(QtCore.QObject):
    finished = QtCore.pyqtSignal(int, object)
    failed = QtCore.pyqtSignal(int, str)


class AnalysisTask(QtCore.QRunnable):
    """Run a pure analytical callable outside the GUI thread."""

    def __init__(self, token, function, *args, **kwargs):
        super().__init__()
        self.token = int(token)
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.signals = TaskSignals()

    @QtCore.pyqtSlot()
    def run(self):
        try:
            value = self.function(*self.args, **self.kwargs)
        except Exception as exc:  # noqa: BLE001 - worker signal contains safe text
            self.signals.failed.emit(self.token, type(exc).__name__)
            return
        self.signals.finished.emit(self.token, value)


class AnalysisTaskRunner(QtCore.QObject):
    """Own worker lifetimes and ignore results from invalidated UI contexts."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._pool = QtCore.QThreadPool(self)
        self._generation = 0
        self._serial = 0
        self._active = {}

    def invalidate(self):
        self._generation += 1
        records = tuple(self._active.values())
        self._active.clear()
        for _task, _generation, _success, _error, on_settled in records:
            if on_settled is not None:
                on_settled()

    def start(
        self,
        function,
        *args,
        on_success,
        on_error=None,
        on_settled=None,
        **kwargs,
    ):
        self._serial += 1
        token = self._serial
        task = AnalysisTask(token, function, *args, **kwargs)
        self._active[token] = (
            task,
            self._generation,
            on_success,
            on_error,
            on_settled,
        )
        task.signals.finished.connect(self._finished)
        task.signals.failed.connect(self._failed)
        self._pool.start(task)
        return token

    def _finished(self, token, value):
        record = self._active.pop(token, None)
        if record is None:
            return
        _task, generation, on_success, _on_error, on_settled = record
        if generation == self._generation:
            on_success(value)
            if on_settled is not None:
                on_settled()

    def _failed(self, token, error_name):
        record = self._active.pop(token, None)
        if record is None:
            return
        _task, generation, _on_success, on_error, on_settled = record
        if generation == self._generation:
            if on_error is not None:
                on_error(error_name)
            if on_settled is not None:
                on_settled()

    def wait_for_done(self, milliseconds=5000):
        self.invalidate()
        self._pool.waitForDone(milliseconds)


__all__ = ["AnalysisTask", "AnalysisTaskRunner", "TaskSignals"]
