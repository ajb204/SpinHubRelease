from pathlib import Path
from types import SimpleNamespace

import pytest
np = pytest.importorskip("numpy")
ng = pytest.importorskip("nmrglue")

from spinDecon.processing.vpar_decon import vpar


def _write_pipe_2d(path, data):
    udic = ng.fileiobase.create_blank_udic(2)
    udic[0]['size'] = data.shape[0]
    udic[1]['size'] = data.shape[1]
    dic = ng.pipe.create_dic(udic)
    ng.pipe.write(str(path), dic, np.asarray(data, dtype=np.float32), overwrite=True)


def test_pseudo_axis_writes_matching_tsv_and_csv(tmp_path):
    inst = vpar.__new__(vpar)
    inst.outdir = str(tmp_path)
    inst.pseudo_axis_info = SimpleNamespace(
        columns=['expno', 'd20'], rows=[['10', '0.1'], ['20', '0.3']]
    )
    inst._write_pseudo_axis_table()
    assert (tmp_path / 'pseudo_axis.tsv').read_text().splitlines() == [
        'spectrum\texpno\td20', '1\t10\t0.1', '2\t20\t0.3'
    ]
    assert (tmp_path / 'pseudo_axis.csv').read_text().splitlines() == [
        'spectrum,expno,d20', '1,10,0.1', '2,20,0.3'
    ]


def test_pseudo2d_first_and_summed_slice(tmp_path):
    data = np.array([[1, 2, 3, 4], [10, 20, 30, 40], [100, 200, 300, 400]], dtype=np.float32)
    _write_pipe_2d(tmp_path / 'test.fid', data)
    inst = vpar.__new__(vpar)
    inst.outdir = str(tmp_path)
    inst.parent = SimpleNamespace(phaseSliceModeValue='First', state=None)
    inst._prepare_pseudo2d_slice_input()
    _, first = ng.pipe.read(str(tmp_path / 'slice.fid'))
    np.testing.assert_allclose(first, data[0])

    inst.parent.phaseSliceModeValue = 'Summed'
    inst._prepare_pseudo2d_slice_input()
    dic, summed = ng.pipe.read(str(tmp_path / 'slice.fid'))
    np.testing.assert_allclose(summed, data.sum(axis=0))
    assert int(float(dic['FDDIMCOUNT'])) == 1
