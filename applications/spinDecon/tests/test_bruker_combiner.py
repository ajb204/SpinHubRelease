from pathlib import Path

from spinDecon.processing.bruker_combiner import (
    combine_bruker_experiments,
    discover_numbered_experiments,
    inspect_combination,
)


def _write_exp(root: Path, number: int, d20: float, p1: float = 10.0, cnst3: float = 2.0):
    exp = root / str(number)
    exp.mkdir()
    (exp / 'acqus').write_text(
        '##TITLE= test\n'
        '##$TD= 256\n##$DTYPA= 0\n##$BYTORDA= 0\n##$AQ_mod= 3\n'
        '##$SW_h= 10000\n##$SFO1= 600.1\n##$BF1= 600.0\n##$O1= 2820\n'
        '##$NUC1= <1H>\n##$DECIM= 32\n##$DSPFVS= 20\n##$GRPDLY= 67.98\n'
        '##$D= (0..20)\n' + ' '.join(['0'] * 20 + [str(d20)]) + '\n'
        '##$P= (0..1)\n0 ' + str(p1) + '\n'
        '##$CNST= (0..3)\n0 0 0 ' + str(cnst3) + '\n'
    )
    # TD=256 int32 words is exactly one 1024-byte Bruker record.
    (exp / 'fid').write_bytes(bytes([number % 256]) * 1024)


def test_combine_numbered_fids(tmp_path):
    _write_exp(tmp_path, 101, 0.1)
    _write_exp(tmp_path, 102, 0.2, p1=11.0, cnst3=3.0)
    _write_exp(tmp_path, 104, 0.5, p1=12.0, cnst3=4.0)
    exps = discover_numbered_experiments(tmp_path, 101, 104)
    assert [x.number for x in exps] == [101, 102, 104]
    info = inspect_combination(exps)
    assert not info.errors
    assert info.numeric_varying_parameters['D20'] == ['0.1', '0.2', '0.5']
    assert info.numeric_varying_parameters['P1'] == ['10.0', '11.0', '12.0']
    assert info.numeric_varying_parameters['CNST3'] == ['2.0', '3.0', '4.0']
    combine_bruker_experiments(tmp_path, exps)
    assert (tmp_path / 'ser').stat().st_size == 3 * 1024
    assert (tmp_path / 'acqus').is_file()
    assert (tmp_path / 'acqu2s').is_file()
    assert (tmp_path / 'pulseprogram').is_file()
    assert (tmp_path / 'decon_expno.list').read_text().split() == ['101', '102', '104']
    assert (tmp_path / 'decon_d20.list').read_text().split() == ['0.1', '0.2', '0.5']
    assert (tmp_path / 'decon_p1.list').read_text().split() == ['10.0', '11.0', '12.0']
    assert (tmp_path / 'decon_cnst3.list').read_text().split() == ['2.0', '3.0', '4.0']

    # The normal conversion detector must read the synthetic lists exactly as
    # it would for an acquired Bruker pseudo-axis and preserve row ordering.
    import sys, types
    sys.modules.setdefault('nmrglue', types.ModuleType('nmrglue'))
    from spinDecon.processing.vpar_decon import detect_bruker_pseudo_axis
    axis = detect_bruker_pseudo_axis(tmp_path)
    assert axis.size == 3
    assert axis.columns == ['expno', 'cnst3', 'd20', 'p1']
    assert axis.rows == [
        ['101', '2.0', '0.1', '10.0'],
        ['102', '3.0', '0.2', '11.0'],
        ['104', '4.0', '0.5', '12.0'],
    ]
    assert (tmp_path / 'combine_manifest.json').is_file()
