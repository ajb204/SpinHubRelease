from __future__ import annotations

from .models import DatasetCapabilities, DatasetStatus, NMRDataset, ResourceState


def capabilities_for(dataset: NMRDataset) -> DatasetCapabilities:
    project = dataset.project
    acquisition = dataset.acquisition
    if project is None:
        return DatasetCapabilities(acquisition is not None, False, False, False, False, False, False)
    if not project.valid:
        return DatasetCapabilities(False, False, False, False, False, False, False)

    r = project.resources
    return DatasetCapabilities(
        can_create_project=False, can_open_project=True,
        can_process_raw=r.raw_available, can_reprocess=r.raw_available,
        can_view_spectrum=r.spectrum_available, can_deconvolve=r.spectrum_available,
        can_view_peaks=r.any_peaks_available,
    )


def status_for(dataset: NMRDataset) -> tuple[DatasetStatus, str]:
    if dataset.project is None:
        return DatasetStatus.ACQUISITION_ONLY, 'Acquisition only - not configured for Decon'

    project = dataset.project
    if not project.valid:
        return DatasetStatus.INVALID_PROJECT, 'Invalid Decon project - ' + (project.error or 'could not read deconParFile')

    r = project.resources
    if r.any_peaks_available:
        if r.raw_available:
            return DatasetStatus.PEAKS_AVAILABLE, 'Peaks available'
        return DatasetStatus.PEAKS_AVAILABLE, 'Peaks available - raw data missing'
    if r.spectrum_available:
        if r.raw_available:
            return DatasetStatus.PROCESSED, 'Spectrum available'
        return DatasetStatus.SPECTRUM_ONLY, 'Spectrum available - raw data missing'
    if r.raw_available:
        if r.spectrum_state is ResourceState.NOT_CONFIGURED:
            return DatasetStatus.READY_TO_PROCESS, 'Raw data available - spectrum not configured'
        return DatasetStatus.READY_TO_PROCESS, 'Raw data available - configured spectrum missing'
    if r.raw_state is ResourceState.NOT_CONFIGURED and r.spectrum_state is ResourceState.NOT_CONFIGURED:
        return DatasetStatus.RESOURCES_UNAVAILABLE, 'Project found - raw data and spectrum not configured'
    return DatasetStatus.RESOURCES_UNAVAILABLE, 'Project found - configured resources missing'


def annotate_dataset(dataset: NMRDataset) -> NMRDataset:
    dataset.capabilities = capabilities_for(dataset)
    dataset.status, dataset.status_text = status_for(dataset)
    return dataset
