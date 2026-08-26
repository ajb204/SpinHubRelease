"""wx orchestration for executing NMRPipe scripts.

Scientific script rendering remains in :mod:`decon.processing.nmrpipe_scripts`;
this module owns progress windows, wx thread handoff, and GUI refresh callbacks.
"""
import logging
import threading
import wx

from spinDecon.gui.dialogs.shell_output import ShellOutputFrame, run_command_with_output
from spinDecon.processing.nmrpipe_scripts import pipefile_for

def execute_process_script(frame, script_path: str, lp: str = 'n', on_finish=None, title: str = 'Processing Output') -> str:
    """Run processing, projection generation, and preview refresh in one output window."""
    pipefile = pipefile_for(frame)
    ensure_dir = getattr(frame, '_spec_output_dir', None)
    if callable(ensure_dir):
        try:
            ensure_dir()
        except Exception:
            logging.exception('Could not ensure spec output directory before running processing script')

    frame.ResetReads()
    output = ShellOutputFrame(frame, title=title)
    processing_steps = ['Process spectrum']
    if lp == 'y':
        processing_steps.append('SMILE reconstruction')
    elif lp == 'm':
        processing_steps.append('MDDNMR reconstruction')
    processing_steps.extend(['Generate projections', 'Refresh preview', 'Refresh display'])
    output.set_workflow(processing_steps, 0)
    output.append_text('Processing spectrum\n')
    output.append_text('Input script: %s\n' % script_path)
    if lp == 'y':
        output.append_text('SMILE reconstruction is enabled; progress will be called out when SMILE starts.\n')
    elif lp == 'm':
        output.append_text('MDDNMR reconstruction is enabled; the platform runtime will be selected at runtime (Apple Silicon uses binMAC_ARM_compat via Rosetta).\n')
    output.Show()

    def finish_ui():
        output.start_step('Refresh display')
        try:
            frame.SetLab(refresh=False)
            frame._update_nmrpipe_file_box(pipefile)
            helper = getattr(frame, '_maybe_invert_test_ft_after_processing', None)
            if callable(helper):
                helper(script_path=script_path, pipefile=pipefile)
            frame.UpdateLampLights()
        except Exception:
            logging.exception('Could not refresh processing state after processing script')
        output.append_text('\nProcessing workflow complete.\n')
        output.finish_workflow(True)
        output.set_status('Complete')
        if on_finish is not None:
            try:
                on_finish()
            except Exception:
                logging.exception('Processing completion callback failed')

    def refresh_slice():
        output.start_step('Refresh preview')
        try:
            frame.RefreshDirectSlice(output_frame=output, on_finish=finish_ui)
        except Exception as exc:
            logging.exception('Could not refresh direct-dimension slice after processing script')
            output.append_text('\nDirect-dimension preview could not be refreshed: %s\n' % exc)
            finish_ui()

    def projection_worker():
        try:
            wx.CallAfter(output.start_step, 'Generate projections')
            wx.CallAfter(output.set_status, 'Generating projections...')
            wx.CallAfter(output.append_text, '\n=== Generating projections ===\n')
            frame.DoProjections(pipefile, output_frame=output)
            wx.CallAfter(refresh_slice)
        except Exception as exc:
            logging.exception('Could not generate projections after processing script')
            wx.CallAfter(output.append_text, '\nProjection generation failed: %s\n' % exc)
            wx.CallAfter(refresh_slice)

    def after_main(rc=0):
        if rc not in (0, None):
            output.finish_workflow(False)
            output.set_status('Failed')
            if on_finish is not None:
                on_finish()
            return
        threading.Thread(target=projection_worker, daemon=True).start()

    run_command_with_output(
        ['csh', script_path], parent=frame, title=title, output_frame=output,
        on_finish=after_main, final=False, label='Run NMRPipe processing script')
    return script_path
