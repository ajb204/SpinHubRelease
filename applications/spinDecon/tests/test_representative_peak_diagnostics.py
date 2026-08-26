import numpy as np
from spinDecon.analysis.peak_picker import PeakPicker, PeakPickerSettings


def test_representative_picker_reports_population_widths_and_selected_widths():
    x = np.arange(80.0)
    data = np.zeros((80, 80))
    for cx, cy, sx, sy, amp in [(15,15,3,4,8),(35,18,3,4,9),(18,50,3,4,10),(50,50,3,4,9),(65,30,1,1,30)]:
        data += amp*np.exp(-((x[:,None]-cx)**2/(2*sx*sx) + (x[None,:]-cy)**2/(2*sy*sy)))
    settings = PeakPickerSettings(threshold_fraction=.08, max_peaks=3, isolation_radius=8,
                                  representative_low_percentile=0, representative_high_percentile=90)
    result = PeakPicker(data, settings).run()
    assert result.status == 'complete'
    assert result.representative_count >= 3
    assert result.representative_widths.shape[1] == 2
    assert result.selected_widths.shape == (3, 2)
    # The anomalously sharp strongest peak should not define the representative fit.
    assert np.median(result.selected_widths[:, 0]) > 3
