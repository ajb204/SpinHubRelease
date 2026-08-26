"""Authoritative application boundary for the Full Peak List.

The Full Peak List is the canonical peak collection for the complete spectrum.
Legacy ``conn_data`` is intentionally absent from this service: connectivity is
an annotation/relationship concern and must never become a second peak store.
"""


class FullPeakListService:
    def __init__(self, legacy_workspace):
        self._legacy = legacy_workspace

    def payload(self):
        getter = getattr(self._legacy, 'get_full_peak_payload', None)
        return getter() if callable(getter) else None

    def open(self, event=None):
        opener = getattr(self._legacy, 'OnButtonFullPeakList', None)
        return opener(event) if callable(opener) else None

    def refresh_viewers(self):
        refresh = getattr(self._legacy, 'refresh_full_peak_list_viewers', None)
        return refresh() if callable(refresh) else None

    def save(self, *, peaks, records, rows, dimension, source_path=None):
        store = getattr(self._legacy, 'store', None) or getattr(self._legacy, 'data_store', None)
        if store is None:
            raise RuntimeError('FullPeakListService requires the application data store')
        store.save_peak_list('full', peaks=peaks, records=records, rows=rows,
                             dimension=dimension, source_path=source_path)
        rebuild = getattr(self._legacy, '_rebuild_projected_peak_lists', None)
        if callable(rebuild):
            rebuild()
        notify = getattr(self._legacy, '_notify_analysis_changed', None)
        if callable(notify):
            notify()
        self.refresh_viewers()
        status = getattr(self._legacy, 'Status', None)
        if callable(status):
            status()
        return peaks
