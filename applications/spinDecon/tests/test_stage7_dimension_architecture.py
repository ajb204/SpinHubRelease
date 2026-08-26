from pathlib import Path


def test_removed_legacy_gui_topology_detectors_do_not_return():
    root = Path(__file__).resolve().parents[1]
    text = "\n".join(
        p.read_text(errors="ignore") for p in root.rglob("*.py")
        if "tests" not in p.parts
    )
    assert "_is_physical_3p" not in text
    assert "canonical_physical_pseudo" not in text


def test_process_frame_does_not_encode_pseudo_in_dim():
    root = Path(__file__).resolve().parents[1]
    text = (root / "gui" / "dialogs" / "processing" / "process.py").read_text()
    assert "self.dim='2p'" not in text
    assert "self.dim='3p'" not in text
    assert "self.dim = '2p'" not in text
    assert "self.dim = '3p'" not in text


def test_migrated_consumers_do_not_branch_on_compatibility_dim_alias():
    """Migrated consumers must branch on explicit spectral/topology fields."""
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    for rel in (
        'gui/workspaces/peak_review.py',
        'gui/dialogs/processing/process.py',
        'gui/dialogs/processing/settings.py',
        'gui/workspaces/slice1d.py',
    ):
        text = (root / rel).read_text()
        assert 'if self.dim' not in text
        assert 'if(self.dim' not in text
        assert 'self.tabOne.dim' not in text


def test_pseudo3d_load_opens_projection_and_uses_topology_for_cache():
    """2 spectral + pseudo must not be excluded by legacy dim>=3 logic."""
    from pathlib import Path
    text = (Path(__file__).resolve().parents[1] / 'gui' / 'workspaces' / 'nmr.py').read_text()
    assert 'needs_projection_cache' in text
    assert 'topology.physical_dim_count == 3' in text
    marker = '# 2 spectral + pseudo has three physical axes'
    assert marker in text
    pseudo_branch = text[text.index(marker):text.index(marker) + 500]
    assert 'self.parent.AddTabTwo(True, self)' in pseudo_branch


def test_decon_dimension_selection_uses_project_state_topology():
    """The NMR tab consumes the canonical Workflow/ProjectState topology."""
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def _selected_topology(self):')
    end = text.index('\n    def ', start + 8)
    block = text[start:end]
    assert 'return self.state.topology()' in block
    assert 'DatasetTopology.from_counts(' not in block


def test_fixed_fit_controls_use_spectral_not_physical_dimensionality():
    """Fixed FIT radii are exposed according to canonical spectral topology."""
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    block = text[text.index('    def SetDim(self):'):text.index('    def DoDim(self, dim):')]
    assert 'show_f1_fit = (topology.spectral_dim_count == 2) or (topology.spectral_dim_count == 1 and topology.has_pseudo_axis)' in block
    assert 'show_f2_fit = (topology.spectral_dim_count == 2)' in block
    assert 'pseudoBox' not in block


def test_decon_load_does_not_infer_scientific_dimensionality_from_array_shape():
    """Loaded ndarray dimensionality is validation input, never the spectral source of truth."""
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def makeinp(self,indir,infile):')
    end = text.index('\n    def ', start + 8)
    makeinp = text[start:end]
    assert 'self.dim=len(self.data.shape)' not in makeinp
    assert 'self.dim = len(self.data.shape)' not in makeinp
    assert 'physical_dim = self.data.ndim' in makeinp
    assert 'state.canonicalize_loaded_dimensions(' in makeinp
    assert 'self.dim = state.spectral_dimensions' in makeinp


def test_decon_load_canonicalizes_nonpseudo_and_pseudo_through_same_state_boundary():
    """Pseudo selection must not choose whether ProjectState participates in load validation."""
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def makeinp(self,indir,infile):')
    end = text.index('\n    def ', start + 8)
    makeinp = text[start:end]
    assert "if(self.pseudoBox.GetValue()==False)" not in makeinp
    assert 'state.canonicalize_loaded_dimensions(' in makeinp
    assert '.pseudoBox' not in makeinp


def test_decon_run_uses_explicit_topology_not_mutable_dim_for_scientific_branches():
    """Launch configuration must distinguish spectral count and pseudo status explicitly."""
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def OnButtonDecon(')
    end = text.index('\n    #Analyse results from decon', start)
    block = text[start:end]
    assert 'topology = self._selected_topology()' in block
    assert 'spectral_dim_count = topology.spectral_dim_count' in block
    assert 'self.dim=(self.dimBox.GetSelection()+1)' not in block
    assert 'if(self.dim' not in block
    assert 'if self.dim' not in block
    assert 'not topology.has_pseudo_axis' in block


def test_decon_analysis_uses_committed_topology_not_gui_dimension_alias():
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def OnButtonAnalyse(')
    end = text.index('\n    def ', start + 8)
    block = text[start:end]
    assert 'topology = self._active_topology()' in block
    assert 'spectral_dim_count = topology.spectral_dim_count' in block
    assert 'self.dim=(self.dimBox.GetSelection()+1)' not in block
    assert 'if(self.dim' not in block
    assert 'if self.dim' not in block


def test_decon_finish_does_not_replace_project_dimension_with_job_dimension():
    """A 2D projection job from a 3D project must not mutate project topology."""
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def _finish_decon_run(')
    end = text.index('\n    def ', start + 8)
    block = text[start:end]
    assert 'topology = self._active_topology()' in block
    assert 'self.dim = run_dim' not in block
    assert 'topology.spectral_dim_count == 3' in block


def test_full_peak_list_uses_canonical_spectral_dimension():
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def load_full_peak_list(')
    end = text.index('\n    def ', start + 8)
    block = text[start:end]
    assert 'spectral_dim_count = self._active_topology().spectral_dim_count' in block
    assert 'dim=int(self.dim)' not in block
    assert 'dimension=int(self.dim)' not in block


def test_decon_output_loading_uses_topology_for_scientific_branches():
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def _load_decon_outputs(')
    end = text.index('\n    def ', start + 8)
    block = text[start:end]
    assert 'topology = self._active_topology()' in block
    assert 'spectral_dim_count = topology.spectral_dim_count' in block
    assert 'not topology.has_pseudo_axis' in block
    assert 'if self.dim' not in block


def test_true_2d_fitting_view_distinguishes_topology_from_array_shape():
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def get_pseudo3d_view(')
    end = text.index('\n    def ', start + 8)
    block = text[start:end]
    assert 'topology.spectral_dim_count == 2' in block
    assert 'not topology.has_pseudo_axis' in block
    assert 'data_ndim == topology.physical_dim_count == 2' in block
    assert "int(getattr(self, 'dim'" not in block


def test_projection_and_four_dimensional_reorganisation_use_spectral_count():
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def OnButtonXA(')
    end = text.index('\n    def CheckFiles(', start)
    block = text[start:end]
    assert 'self._active_topology().spectral_dim_count != 4' in block
    project = block[block.index('    def OnButtonProject('):]
    assert 'spectral_dim_count = self._active_topology().spectral_dim_count' in project
    assert 'if(self.dim' not in project


def test_full_spectrum_setup_branches_on_canonical_topology():
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def makeinp(')
    end = text.index('\n    def SetDmax(', start)
    block = text[start:end]
    setup = block[block.index('        topology = self._active_topology()'):]
    assert 'spectral_dim_count = topology.spectral_dim_count' in setup
    assert 'if(self.dim==' not in setup
    assert 'elif(self.dim==' not in setup
    assert "if lab1 == 'ID'" not in setup
    assert "if 'time_T2' in lab" not in setup
    assert 'self.pseudo = topology.has_pseudo_axis' in setup
    assert 'dim=spectral_dim_count' in setup


def test_reference_peak_loading_uses_topology_not_frame_dim_or_gui_pseudo():
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def ReadPeakListFile(')
    end = text.index('\n    ###############################################\n    #Automatic setting', start)
    block = text[start:end]
    assert 'topology = self._active_topology()' in block
    assert 'spectral_dim_count = topology.spectral_dim_count' in block
    assert 'self.dim==' not in block
    assert 'self.dim ==' not in block
    assert 'self.dim <' not in block
    assert 'self.pseudoBox.GetValue()' not in block
    assert 'not topology.has_pseudo_axis' in block


def test_pseudo3d_recon_uses_physical_dim_for_spinunidec_protocol_dispatch():
    """2 spectral + pseudo must dispatch spinUnidec Protocol3P as dim=3."""
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def OnButtonDecon(')
    end = text.index('\n    #Analyse results from decon', start)
    block = text[start:end]
    assert 'topology.physical_dim_count == 3' in block
    assert "dimVal=str(topology.physical_dim_count)" in block
    assert "decset['pseudo3D']='1'" in block
    assert "decset['FIT']='1'" in block


def test_pseudo3d_protocol_dimension_does_not_redefine_gui_spectral_dimension():
    """The external dim=3 selector must not leak back into canonical spectral state."""
    root = Path(__file__).resolve().parents[1]
    text = (root / 'gui' / 'workspaces' / 'nmr.py').read_text()
    start = text.index('    def OnButtonDecon(')
    end = text.index('\n    #Analyse results from decon', start)
    block = text[start:end]
    assert 'spectral_dim_count = topology.spectral_dim_count' in block
    assert 'self.state.spectral_dimensions = spectral_dim_count' in block
    assert 'self.state.spectral_dimensions = topology.physical_dim_count' not in block
