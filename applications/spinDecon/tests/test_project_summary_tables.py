import io

from spinDecon.project.summary import _write_matrix_table


def test_write_matrix_table_accepts_flat_pseudo3d_rows():
    out = io.StringIO()
    _write_matrix_table(
        out,
        'Fitting table',
        ['Peak', 'grp', '%err'],
        [['P1', '1', '2.5'], ['P2', '-', '3.0']],
    )
    tex = out.getvalue()
    assert r'\textbf{Peak} & \textbf{grp} & \textbf{\%err}' in tex
    assert 'P1 & 1 & 2.5' in tex


def test_write_matrix_table_keeps_legacy_parameter_matrix_shape():
    out = io.StringIO()
    _write_matrix_table(
        out,
        '',
        ['F2', 'F1'],
        [('Phase p0', ['1', '2'])],
    )
    tex = out.getvalue()
    assert r'\textbf{Parameter} & \textbf{F2} & \textbf{F1}' in tex
    assert 'Phase p0 & 1 & 2' in tex

from spinDecon.project.summary import _write_pseudo_report


def test_single_peak_pseudo_group_keeps_card_and_plot_side_by_side():
    out = io.StringIO()
    pseudo = {
        'columns': ['Peak', 'grp', '%err', 'f01(ppm)', 'w1(Hz)', 'g1', 'f02(ppm)', 'w2(Hz)', 'g2'],
        'rows': [['P1', '-', '2.0', '8.1', '12', '0.5', '120.0', '15', '0.5']],
        'figures': [('pseudo3d_fit_001.pdf', {'group': None, 'peaks': ['P1']})],
        'analysis': {'columns': ['Peak'], 'rows': [], 'peak_figures': {}},
    }
    _write_pseudo_report(out, pseudo)
    tex = out.getvalue()
    assert r'\begin{minipage}[t]{0.30\linewidth}' in tex
    assert r'\begin{minipage}[t]{0.67\linewidth}\vspace{0pt}\centering' in tex
    assert r'pseudo3d\_fit\_001.pdf' in tex

from spinDecon.project.summary import _append_fitting_to_full_peaks, _write_long_table_with_split


def test_physical_2d_full_peak_list_appends_fitting_values_with_separator():
    pseudo = {
        'kind': '2d',
        'columns': ['Peak', 'grp', '%err', 'f01(ppm)', 'f02(ppm)', 'w1(Hz)', 'w2(Hz)'],
        'rows': [['P1', '7', '1.2', '8.101', '120.2', '11', '14']],
    }
    columns, rows, split = _append_fitting_to_full_peaks(
        ['Peak', 'F2(ppm)', 'F1(ppm)', 'Intensity'],
        [['P1', '8.10', '120.1', '42'], ['P2', '7.20', '115.0', '12']],
        pseudo,
    )
    assert split == 4
    assert columns == ['Peak', 'F2(ppm)', 'F1(ppm)', 'Intensity', 'grp', '%err', 'f01(ppm)', 'f02(ppm)', 'w1(Hz)', 'w2(Hz)']
    assert rows[0][-6:] == ['7', '1.2', '8.101', '120.2', '11', '14']
    assert rows[1][-6:] == ['', '', '', '', '', '']

    out = io.StringIO()
    _write_long_table_with_split(out, columns, rows, split)
    assert r'@{}llll|llllll@{}' in out.getvalue()


def test_physical_2d_fitting_report_uses_group_panel_without_pseudo_slice_gallery():
    out = io.StringIO()
    pseudo = {
        'kind': '2d',
        'columns': ['Peak', 'grp', '%err', 'f01(ppm)', 'w1(Hz)', 'g1', 'f02(ppm)', 'w2(Hz)', 'g2'],
        'rows': [
            ['P1', '2', '1.0', '8.1', '12', '0.5', '120.0', '15', '0.5'],
            ['P2', '2', '1.5', '8.2', '13', '0.6', '121.0', '16', '0.6'],
        ],
        'figures': [('pseudo3d_fit_001.pdf', {
            'group': '2', 'peaks': ['P1', 'P2'],
            'slice_figures': ['pseudo3d_fit_001_slice_001.pdf'],
        })],
        'analysis': None,
    }
    _write_pseudo_report(out, pseudo)
    tex = out.getvalue()
    assert r'\section*{2D fitting results}' in tex
    assert r'\subsection*{Group 2: P1, P2}' in tex
    assert r'pseudo3d\_fit\_001.pdf' in tex
    assert r'pseudo3d\_fit\_001\_slice\_001.pdf' not in tex
    assert 'Fit across pseudo-axis slices' not in tex

from spinDecon.project.summary import _write_joined_results_table, _peak_tables


def test_fitted_parameters_master_table_is_page_breakable_and_repeats_header():
    out = io.StringIO()
    _write_joined_results_table(out, ['Peak', '%err', 'R'], [['P1', '1.2', '3.4']], 2)
    tex = out.getvalue()
    assert r'\begin{longtable}' in tex
    assert r'\endfirsthead' in tex
    assert r'\endhead' in tex
    assert r'\tiny' in tex
    assert r'll|l' in tex


def test_2d_full_peak_table_reports_snr_and_keeps_header_row_alignment():
    class Frame:
        dim = 2
        def get_full_peak_payload(self):
            return {'rows': [['P1', '8.1239', '120.9876', '50.0']]}
        def get_full_peak_headers(self, row_width=None):
            return ['Name', 'F2 (ppm)', 'F1 (ppm)', 'Intensity']
        def get_reference_peak_headers(self):
            return []
        def get_noise_sigma(self):
            return 5.0

    _refs, _ref_rows, headers, rows = _peak_tables(Frame())
    assert headers == ['Name', 'F2 (ppm)', 'F1 (ppm)', 'SNR']
    assert len(headers) == len(rows[0])
    assert rows[0] == ['P1', '8.123', '120.987', '10']


def test_fitting_append_normalises_legacy_row_width_to_header_width():
    pseudo = {
        'kind': '2d',
        'columns': ['Peak', 'grp', '%err'],
        'rows': [['P1', '7', '1.2', 'unexpected legacy value']],
    }
    columns, rows, split = _append_fitting_to_full_peaks(
        ['Name', 'SNR'], [['P1', '10']], pseudo)
    assert split == 2
    assert columns == ['Name', 'SNR', 'grp', '%err']
    assert rows == [['P1', '10', '7', '1.2']]
    assert len(columns) == len(rows[0])

from spinDecon.project.summary import _write_analysis_summary


def test_cpmg_analysis_summary_reports_screen_and_shared_global_parameters():
    out = io.StringIO()
    _write_analysis_summary(out, {'analysis': {
        'name': 'CPMG',
        'screen': {'Rex_threshold': 2.0, 'n_total': 12, 'n_significant': 3},
        'global': {'success': True, 'kex': 850.0, 'pb': 0.035, 'chi2': 1.2, 'n_peaks': 3},
        'summary_figures': [],
    }})
    tex = out.getvalue()
    assert 'CPMG Rex screen' in tex
    assert '3 of 12 peaks significant' in tex
    assert '$k_{ex}$=850.0' in tex
    assert '$p_b$=0.035' in tex
