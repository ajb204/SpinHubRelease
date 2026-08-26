"""GUI adapter for the headless NMRPipe script helper."""
from spinDecon.processing.nmrpipe_scripts import nmrPipe as _ProcessingNmrPipe
from spinDecon.gui.dialogs.processing.execution import execute_process_script


class NmrPipeGuiAdapter(_ProcessingNmrPipe):
    """Bind script-generation services to wx execution orchestration."""

    def execute_process_script(self, frame=None, script_path='', lp='n', on_finish=None, title='Processing Output'):
        return execute_process_script(
            self._frame(frame), script_path, lp=lp, on_finish=on_finish, title=title
        )

    ExecuteProcessScript = execute_process_script
    run_process_script = execute_process_script
    RunProcessScript = execute_process_script


# Historical GUI spelling.
nmrPipe = NmrPipeGuiAdapter
