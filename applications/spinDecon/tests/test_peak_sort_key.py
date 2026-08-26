from spinDecon.domain.peaks import peak_sort_key


def test_peak_sort_uses_longest_digit_run():
    labels = ['A8-N15', 'A12-N9', 'Peak105-H2', 'Peak7-H1234']
    assert sorted(labels, key=peak_sort_key) == ['A12-N9', 'A8-N15', 'Peak105-H2', 'Peak7-H1234']


def test_peak_sort_is_stable_for_unnumbered_and_ties():
    labels = ['beta', 'A12-X34', 'A34-X12', 'alpha']
    assert sorted(labels, key=peak_sort_key) == ['A12-X34', 'A34-X12', 'alpha', 'beta']
