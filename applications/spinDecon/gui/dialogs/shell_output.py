"""Reusable shell-output window and background command runner."""

from __future__ import annotations

import logging
import os
import subprocess
import threading
from typing import Optional, Sequence, Union, List

import wx


class NMRProgressBar(wx.Control):
    """Compact stage progress with a spin-half attracted to a bar magnet."""

    def __init__(self, parent):
        super().__init__(parent, size=(-1, 38), style=wx.BORDER_NONE)
        self._range = 1
        self._value = 0
        self._phase = 0
        self._running = True
        self.SetMinSize((300, 38))
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda evt: None)

    def SetRange(self, value):
        self._range = max(1, int(value))
        self.Refresh(False)

    def SetValue(self, value):
        self._value = max(0, min(int(value), self._range))
        self.Refresh(False)

    def SetAnimationPhase(self, phase):
        self._phase = int(phase)
        self.Refresh(False)

    def SetRunning(self, running):
        self._running = bool(running)
        self.Refresh(False)

    def _on_paint(self, event):
        dc = wx.AutoBufferedPaintDC(self)
        rect = self.GetClientRect()
        dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
        dc.Clear()
        w, h = rect.width, rect.height
        if w < 120 or h < 20:
            return

        # Leave a fixed magnet bay on the right; the spin travels through the
        # remaining track according to *real stage progress* only.
        mag_w = 62
        margin = 9
        track_left = margin + 10
        track_right = max(track_left + 30, w - mag_w - 18)
        cy = h // 2
        track_h = 10
        dc.SetPen(wx.Pen(wx.Colour(166, 177, 190), 1))
        dc.SetBrush(wx.Brush(wx.Colour(235, 239, 243)))
        dc.DrawRoundedRectangle(track_left, cy-track_h//2, track_right-track_left, track_h, 5)
        frac = float(self._value) / float(max(1, self._range))
        fill_w = int((track_right-track_left) * frac)
        if fill_w:
            dc.SetPen(wx.TRANSPARENT_PEN)
            dc.SetBrush(wx.Brush(wx.Colour(82, 145, 203)))
            dc.DrawRoundedRectangle(track_left, cy-track_h//2, fill_w, track_h, 5)

        # The spin's slow bob/tilt is cosmetic.  Its horizontal position never
        # advances except when the workflow stage advances.
        import math
        travel = max(0, track_right-track_left-18)
        sx = track_left + 9 + int(travel * frac)
        bob = int(round(1.5 * math.sin(self._phase * math.pi / 4.0))) if self._running else 0
        sy = cy + bob
        r = 7
        dc.SetPen(wx.Pen(wx.Colour(60, 65, 72), 1))
        dc.SetBrush(wx.Brush(wx.Colour(220, 64, 70)))
        dc.DrawCircle(sx, sy, r)
        # Spin axis/arrow, with a small alternating lean.
        lean = (-2, -1, 0, 1, 2, 1, 0, -1)[self._phase % 8] if self._running else 0
        ax, ay = sx + lean, sy - r - 6
        dc.SetPen(wx.Pen(wx.Colour(45, 48, 54), 2))
        dc.DrawLine(sx-r//2, sy+r+4, ax, ay)
        dc.SetBrush(wx.Brush(wx.Colour(45, 48, 54)))
        dc.DrawPolygon([wx.Point(ax, ay-3), wx.Point(ax-3, ay+3), wx.Point(ax+3, ay+3)])
        # A subtle precession arc provides motion without visual noise.
        dc.SetPen(wx.Pen(wx.Colour(92, 145, 190), 1))
        dc.SetBrush(wx.TRANSPARENT_BRUSH)
        dc.DrawEllipse(sx-r-3, sy-3, 2*r+6, 6)

        # Cartoon bar magnet, intentionally simple enough to remain legible at
        # native wx sizes and on high-DPI displays.
        mx = w - mag_w - margin
        my = cy - 10
        half = mag_w // 2
        dc.SetPen(wx.Pen(wx.Colour(80, 86, 94), 1))
        dc.SetBrush(wx.Brush(wx.Colour(215, 60, 65)))
        dc.DrawRoundedRectangle(mx, my, half+3, 20, 4)
        dc.SetBrush(wx.Brush(wx.Colour(55, 112, 190)))
        dc.DrawRoundedRectangle(mx+half-2, my, half+2, 20, 4)
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        font.SetPointSize(max(7, font.GetPointSize()-1))
        dc.SetFont(font)
        dc.SetTextForeground(wx.WHITE)
        dc.DrawText('N', mx+9, my+2)
        dc.DrawText('S', mx+half+8, my+2)

        # Two faint field lines nearest the spin/magnet reinforce the visual
        # metaphor without adding another animation or expensive drawing.
        dc.SetPen(wx.Pen(wx.Colour(174, 204, 229), 1))
        for off in (-5, 5):
            dc.DrawLine(min(sx+r+4, mx-2), sy+off, mx-2, cy+off//2)


class ShellOutputFrame(wx.Frame):
    """Terminal output window with an optional lightweight convergence plot."""

    def __init__(self, parent, title: str = 'Shell Output', convergence_file: Optional[str] = None, decon_profile=None):
        if convergence_file:
            # Match the decon progress window height to its owning/main window.
            # The historical main decon frame is 585 px high; use that only as
            # a fallback when the parent size is not yet available.
            main_height = 585
            try:
                if parent is not None:
                    ph = int(parent.GetSize().GetHeight())
                    if ph > 0:
                        main_height = ph
            except Exception:
                pass
            size = (720, main_height)
        else:
            # Keep utility windows within the owning window vertically.
            main_height = 585
            try:
                if parent is not None:
                    ph = int(parent.GetSize().GetHeight())
                    if ph > 0:
                        main_height = ph
            except Exception:
                pass
            size = (680, min(520, main_height))
        super().__init__(parent, title=title, size=size)

        # Create frame-managed chrome before the sole client-area panel.
        # This mirrors peakFrame: wx.Frame then reserves the status-bar area
        # and automatically keeps the panel fitted to the remaining client area.
        self.statusbar = self.CreateStatusBar(1)
        self.statusbar.SetStatusText('Queued...')

        panel = wx.Panel(self)
        panel.SetBackgroundColour(wx.Colour(255, 255, 255))

        self.convergence_file = convergence_file
        self._decon_profile = decon_profile or {}
        self._decon_run_index = 0
        self._decon_run_start_iter = 0
        self._decon_run_active = False
        self._conv_offset = 0
        self._conv_partial = ''
        self._conv_x = []
        self._conv_y = []
        self._conv_timer = None
        self._conv_line = None
        self._conv_canvas = None
        self._conv_axes = None
        self._conv_events = {}
        self._conv_event_artists = []
        self._conv_conv = []
        self._conv_axes2 = None
        self._conv_line2 = None
        self.detailRadio = None
        self.fullRadio = None
        self.stageRadio = None
        self._peak_axes = None
        self._peak_line = None
        self._toolbar = None
        self._smoothing_note = None

        root = wx.BoxSizer(wx.VERTICAL)

        # Compact, user-facing workflow progress.  This deliberately reports
        # coarse stages rather than pretending NMRPipe exposes a meaningful
        # percentage within each numerical operation.
        self._workflow_steps = []
        self._workflow_index = -1
        self._workflow_done = set()
        self._smile_announced = False
        self._progress_anim_phase = 0
        self._progress_anim_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self._on_progress_animation, self._progress_anim_timer)
        self.progressPanel = wx.Panel(panel)
        ps = wx.BoxSizer(wx.VERTICAL)
        self.progressTitle = wx.StaticText(self.progressPanel, label='Preparing...')
        f = self.progressTitle.GetFont(); f.SetWeight(wx.FONTWEIGHT_BOLD); self.progressTitle.SetFont(f)
        # Custom NMR progress strip: a spin-half is gently pulled toward a
        # bar magnet as workflow stages complete.  Animation is timer-driven
        # on the GUI thread only; it never polls or touches the calculation.
        self.progressGauge = NMRProgressBar(self.progressPanel)
        self.progressHint = wx.StaticText(self.progressPanel, label='Progress is shown by calculation stage.')
        ps.Add(self.progressTitle, 0, wx.BOTTOM, 4)
        ps.Add(self.progressGauge, 0, wx.EXPAND | wx.BOTTOM, 3)
        ps.Add(self.progressHint, 0)
        self.progressPanel.SetSizer(ps)
        root.Add(self.progressPanel, 0, wx.ALL | wx.EXPAND, 8)

        # Deconvolution runs use one main display area.  The graph is shown by
        # default; radio buttons swap that same area between the live graph and
        # the continuously collected shell output.
        self._view_panel = None
        self._view_sizer = None
        self.graphRadio = None
        self.shellRadio = None
        if convergence_file:
            # Decon starts compact.  Progress figure and shell transcript are
            # independent disclosures so routine users see only the stage bar.
            self.graphRadio = None
            self.shellRadio = None
            self._view_panel = wx.Panel(panel)
            self._view_sizer = wx.BoxSizer(wx.VERTICAL)
            self._view_panel.SetSizer(self._view_sizer)
            self._create_convergence_plot(self._view_panel, self._view_sizer)
            self.output = wx.TextCtrl(
                self._view_panel,
                style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH | wx.TE_DONTWRAP,
            )
            self._view_sizer.Add(self.output, 1, wx.EXPAND)
            if self._conv_canvas is not None:
                self._conv_canvas.Hide()
            if self._toolbar is not None:
                self._toolbar.Hide()
            self.output.Hide()
            self._view_panel.Hide()
            root.Add(self._view_panel, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 8)
        else:
            self.output = wx.TextCtrl(
                panel,
                style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH | wx.TE_DONTWRAP,
            )
            root.Add(self.output, 1, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 8)
            # Routine users see only the stage progress; the complete transcript
            # remains one click away for diagnostics and power users.
            self.output.Hide()

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        self.progressViewBtn = None
        self.outputViewBtn = None
        if convergence_file:
            self.progressViewBtn = wx.ToggleButton(panel, label='Show progress')
            self.outputViewBtn = wx.ToggleButton(panel, label='Show output')
            btn_row.Add(self.progressViewBtn, 0, wx.RIGHT, 8)
            btn_row.Add(self.outputViewBtn, 0)
            btn_row.AddStretchSpacer(1)
            self.progressViewBtn.Bind(wx.EVT_TOGGLEBUTTON, self._on_decon_disclosure)
            self.outputViewBtn.Bind(wx.EVT_TOGGLEBUTTON, self._on_decon_disclosure)
        self.detailsBtn = None
        if not convergence_file:
            self.detailsBtn = wx.Button(panel, label='Show details')
            btn_row.Add(self.detailsBtn, 0, wx.RIGHT, 8)
            self.detailsBtn.Bind(wx.EVT_BUTTON, self._on_details)
        self.clearBtn = wx.Button(panel, label='Clear')
        self.closeBtn = wx.Button(panel, label='Close')
        btn_row.Add(self.clearBtn, 0, wx.RIGHT, 8)
        btn_row.Add(self.closeBtn, 0)

        self.clearBtn.Bind(wx.EVT_BUTTON, self._on_clear)
        self.closeBtn.Bind(wx.EVT_BUTTON, self._on_close)
        self.Bind(wx.EVT_CLOSE, self._on_close)

        # The row must expand across the panel: the stretch spacer then keeps
        # view controls on the left and Clear/Close on the right.
        root.Add(btn_row, 0, wx.ALL | wx.EXPAND, 8)

        panel.SetSizer(root)
        self._details_height = size[1]
        if convergence_file:
            self.SetMinSize((560, 180))
            self._details_height = size[1]
            self.SetSize((size[0], min(230, size[1])))
            self.configure_decon_progress(self._decon_profile)
        else:
            self.SetMinSize((560, 180))
            self.SetSize((size[0], min(220, size[1])))

    def configure_decon_progress(self, profile):
        """Configure weighted decon progress; expensive DoRun phases own 82%."""
        self._decon_profile = dict(profile or {})
        dim = int(self._decon_profile.get('dim', 1) or 1)
        enhance = str(self._decon_profile.get('enhance', '0')) == '1'
        recon = str(self._decon_profile.get('recon', '0')) == '1'
        bore = str(self._decon_profile.get('dec3d', '0')) not in ('0', 'False', 'false', '')
        maxiter = int(float(self._decon_profile.get('maxIter', 10000) or 10000))
        maxiter3d = int(float(self._decon_profile.get('maxIter3D', maxiter) or maxiter))
        if enhance:
            runs = [maxiter] if dim < 3 or not bore else [10000, maxiter3d]
            names = ['Enhancing spectrum'] if len(runs)==1 else ['Optimising 1D traces', 'Optimising 3D spectrum']
        elif recon:
            runs = [maxiter3d if dim >= 3 else maxiter]
            names = ['Restricted reconstruction']
        elif dim == 1:
            runs = [10000, 10000, 10000]; names = ['Initial optimisation', 'Refining peak shape', 'Final optimisation']
        elif dim == 2:
            runs = [maxiter, maxiter]; names = ['Initial 2D optimisation', 'Refining 2D solution']
        elif dim >= 3 and bore:
            runs = [10000, 10000, 10000]
            names = ['Initial 1D optimisation', 'Refining 1D traces', 'Final 1D optimisation']
            if maxiter3d > 0:
                runs += [maxiter3d, maxiter3d]; names += ['Building 3D solution', 'Refining 3D spectrum']
        else:
            runs = [maxiter3d, maxiter3d]; names = ['Initial 3D optimisation', 'Refining 3D spectrum']
        self._decon_profile['_run_max'] = [max(1,int(x)) for x in runs]
        self._decon_profile['_run_names'] = names
        self._decon_profile['_run_weights'] = [0.82/len(runs)]*len(runs) if runs else []
        self._decon_run_index = 0; self._decon_run_start_iter = 0; self._decon_run_active = False
        self.progressGauge.SetRange(1000); self.progressGauge.SetValue(25)
        self.progressTitle.SetLabel('Preparing deconvolution')
        self.progressHint.SetLabel('Slow optimisation stages occupy most of the progress bar.')
        if not self._progress_anim_timer.IsRunning(): self._progress_anim_timer.Start(700)

    def _update_decon_progress(self, system_iter=None, event_label=None):
        if not self._decon_profile: return
        runs=self._decon_profile.get('_run_max', []); names=self._decon_profile.get('_run_names', [])
        weights=self._decon_profile.get('_run_weights', [])
        if event_label:
            low=event_label.lower()
            if low.startswith('dorun start'):
                self._decon_run_active=True
                self._decon_run_start_iter=int(system_iter or 0)
                idx=min(self._decon_run_index, max(0,len(names)-1))
                if names: self.progressTitle.SetLabel(names[idx])
            elif low.startswith('dorun end'):
                self._decon_run_active=False
                self._decon_run_index=min(self._decon_run_index+1, len(runs))
                if self._decon_run_index < len(names):
                    self.progressTitle.SetLabel('Preparing — ' + names[self._decon_run_index])
                else: self.progressTitle.SetLabel('Generating outputs')
            elif low.startswith('squash'):
                self.progressTitle.SetLabel('Refining candidate peaks')
            elif low.startswith('cull'):
                self.progressTitle.SetLabel('Removing weak candidates')
            elif low.startswith('setpeaks') or 'peak-shape' in low:
                self.progressTitle.SetLabel('Updating peak shape')
        base=0.08
        done=sum(weights[:min(self._decon_run_index,len(weights))])
        within=0.0
        if self._decon_run_active and self._decon_run_index < len(runs) and system_iter is not None:
            local=max(0,int(system_iter)-self._decon_run_start_iter)
            within=weights[self._decon_run_index]*min(1.0, local/float(runs[self._decon_run_index]))
            self.progressHint.SetLabel('Iteration %d of up to %d in this optimisation.' % (local, runs[self._decon_run_index]))
        frac=min(0.90, base+done+within)
        self.progressGauge.SetValue(int(round(frac*1000)))
        self.progressPanel.Layout()

    def complete_decon_progress(self, success=True, label=None):
        if not self._decon_profile: return
        if success:
            self.progressGauge.SetValue(1000)
            self.progressTitle.SetLabel(label or 'Complete — outputs generated')
            self.progressHint.SetLabel('Calculation and output preparation complete.')
            self.set_status('Complete')
        else:
            self.progressTitle.SetLabel('Calculation stopped — inspect output')
            self.progressHint.SetLabel('The shell output may contain details of the problem.')
            self.set_status('Failed')
        self.progressGauge.SetRunning(False)
        if self._progress_anim_timer.IsRunning(): self._progress_anim_timer.Stop()
        self.progressPanel.Layout()

    def _on_decon_disclosure(self, event=None):
        show_graph=bool(self.progressViewBtn and self.progressViewBtn.GetValue())
        show_output=bool(self.outputViewBtn and self.outputViewBtn.GetValue())
        if self._conv_canvas is not None: self._conv_canvas.Show(show_graph)
        if self._toolbar is not None: self._toolbar.Show(show_graph)
        self.output.Show(show_output)
        self._view_panel.Show(show_graph or show_output)
        self.progressViewBtn.SetLabel('Hide progress' if show_graph else 'Show progress')
        self.outputViewBtn.SetLabel('Hide output' if show_output else 'Show output')
        try:
            w=self.GetSize().GetWidth()
            self.SetSize((w, self._details_height if (show_graph or show_output) else min(230,self._details_height)))
        except Exception: pass
        self.Layout(); self._view_panel.Layout()
        if show_graph and self._conv_canvas is not None: self._conv_canvas.draw_idle()

    def _create_convergence_plot(self, panel, root):
        """Create linked system-state and relative-change progress plots."""
        try:
            from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg
            from matplotlib.figure import Figure

            figure = Figure(figsize=(6.8, 5.2))
            # Equal-height panels with a deliberately small gap.  sharex links
            # toolbar zoom/pan and programmatic x limits between both plots.
            gs = figure.add_gridspec(2, 1, height_ratios=(1, 1), hspace=0.08,
                                     left=0.10, right=0.88, top=0.96, bottom=0.10)
            self._conv_axes = figure.add_subplot(gs[0])
            self._conv_axes2 = figure.add_subplot(gs[1], sharex=self._conv_axes)
            self._peak_axes = self._conv_axes2.twinx()
            # Keep axis/title typography at the same compact scale as tick labels.
            plot_fontsize = 8
            self._plot_fontsize = plot_fontsize
            self._conv_axes.set_ylabel('System state, tack (x1e12)', fontsize=plot_fontsize)
            self._conv_axes2.set_ylabel('Relative change', fontsize=plot_fontsize)
            self._peak_axes.set_ylabel('Peak count', fontsize=plot_fontsize)
            self._conv_axes2.set_xlabel('System iteration', fontsize=plot_fontsize)
            for ax in (self._conv_axes, self._conv_axes2, self._peak_axes):
                ax.tick_params(axis='both', labelsize=plot_fontsize)
            # The plots share X, so only the lower panel needs iteration labels.
            self._conv_axes.tick_params(axis='x', labelbottom=False)
            self._conv_axes2.set_yscale('log')
            self._conv_axes.grid(True, alpha=0.18)
            self._conv_axes2.grid(True, which='major', alpha=0.18)
            (self._conv_line,) = self._conv_axes.plot([], [], linewidth=1.25, label='System state')
            (self._conv_line2,) = self._conv_axes2.plot([], [], linewidth=0.95, label='Relative change')
            # Explicitly distinct from the relative-change trace.
            (self._peak_line,) = self._peak_axes.plot([], [], linewidth=1.15,
                                                       color='tab:orange', label='Peak count')
            self._peak_axes.legend(loc='lower left', fontsize=8, framealpha=0.85)
            self._conv_canvas = FigureCanvasWxAgg(panel, -1, figure)
            root.Add(self._conv_canvas, 1, wx.EXPAND)
            self._toolbar = NavigationToolbar2WxAgg(self._conv_canvas)
            self._toolbar.Realize()
            root.Add(self._toolbar, 0, wx.EXPAND)
            self._conv_timer = wx.Timer(self)
            self.Bind(wx.EVT_TIMER, self._on_convergence_timer, self._conv_timer)
            self._conv_timer.Start(1000)
        except Exception:
            logging.exception('Could not create convergence plot')
            self._conv_timer = None

    def _on_view_change(self, event=None):
        if self._conv_canvas is None or self.output is None:
            return
        show_graph = bool(self.graphRadio and self.graphRadio.GetValue())
        self._conv_canvas.Show(show_graph)
        if self._toolbar is not None:
            self._toolbar.Show(show_graph)
        self.output.Show(not show_graph)
        if self._view_panel is not None:
            self._view_panel.Layout()
        if show_graph:
            self._conv_canvas.draw_idle()

    def _on_scale_change(self, event=None):
        if self._conv_x:
            self._redraw_convergence()

    @staticmethod
    def _median(values):
        vals = sorted(values)
        n = len(vals)
        if not n: return 0.0
        m = n // 2
        return vals[m] if n % 2 else 0.5 * (vals[m-1] + vals[m])

    @classmethod
    def _mad_limits(cls, values, sigma=8.0):
        """Robust limits: median +/- sigma*MAD, with a small visual pad."""
        vals = [v for v in values if v == v and abs(v) != float('inf')]
        if not vals: return None
        med = cls._median(vals)
        mad = cls._median([abs(v-med) for v in vals])
        if mad <= 0:
            spread = max(max(vals)-min(vals), abs(med)*1e-8, 1.0)
            return med-spread*.55, med+spread*.55
        # 1.4826 converts MAD to a normal-distribution sigma estimate.
        half = sigma * 1.4826 * mad
        lo, hi = med-half, med+half
        pad = .06 * (hi-lo)
        return lo-pad, hi+pad

    def _current_stage_start(self):
        starts = []
        for iteration, labels in self._conv_events.items():
            if any(label.lower().startswith('dorun start') for label in labels):
                starts.append(iteration)
        return max(starts) if starts else (self._conv_x[0] if self._conv_x else 0)

    @classmethod
    def _adaptive_smooth(cls, values):
        """Return display values, window and instability score.

        Instability is estimated robustly as iteration-to-iteration motion relative
        to the stage's robust amplitude.  Stable stages are left untouched.  For
        unstable stages the window grows with both noise and sample count, but is
        capped at 15% of the stage so real protocol-scale trends remain visible.
        The moving average is edge-corrected and never crosses a DoRun boundary.
        """
        import math
        vals=[float(v) for v in values]
        n=len(vals)
        if n < 9:
            return vals, 1, 0.0
        diffs=[abs(vals[i]-vals[i-1]) for i in range(1,n)]
        med=cls._median(vals)
        mad=cls._median([abs(v-med) for v in vals])
        robust_amp=max(1.4826*mad, abs(med)*1e-12, 1e-30)
        step=cls._median(diffs)
        score=step/robust_amp
        # A score below ~0.12 means successive motion is small compared with
        # the stage envelope; smoothing would hide rather than clarify detail.
        if score < 0.12:
            return vals, 1, score
        # Choose roughly 3-15% of the stage. More erratic traces receive a
        # wider window. Keep it odd for a centred average.
        frac=min(0.15, max(0.03, 0.03 + 0.12*min(1.0, (score-0.12)/0.9)))
        window=max(5, int(round(n*frac)))
        if window % 2 == 0: window += 1
        maxw=max(5, int(n*0.15))
        if maxw % 2 == 0: maxw -= 1
        window=min(window, maxw, 151)
        if window < 5:
            return vals, 1, score
        half=window//2
        # Prefix sums make this O(n), important during live redraws.
        prefix=[0.0]
        for v in vals: prefix.append(prefix[-1]+v)
        out=[]
        for i in range(n):
            a=max(0,i-half); b=min(n,i+half+1)
            out.append((prefix[b]-prefix[a])/(b-a))
        return out, window, score

    def _stage_ranges(self, first, last):
        """Split an index interval at DoRun starts for stage-local smoothing."""
        starts=[]
        for iteration, labels in self._conv_events.items():
            if any(label.lower().startswith('dorun start') for label in labels):
                starts.append(iteration)
        cuts=[first]
        for it in sorted(starts):
            j=first
            while j < last and self._conv_x[j] < it: j += 1
            if first < j < last and j not in cuts: cuts.append(j)
        cuts.append(last)
        return [(cuts[i],cuts[i+1]) for i in range(len(cuts)-1) if cuts[i+1]>cuts[i]]

    def _on_convergence_timer(self, event=None):
        if not self.convergence_file or self._conv_line is None: return
        try:
            if not os.path.exists(self.convergence_file): return
            size = os.path.getsize(self.convergence_file)
            if size < self._conv_offset:
                self._conv_offset = 0; self._conv_partial = ''
                self._conv_x.clear(); self._conv_y.clear(); self._conv_conv.clear(); self._conv_events.clear()
                for artist in self._conv_event_artists:
                    try: artist.remove()
                    except Exception: pass
                self._conv_event_artists = []
            with open(self.convergence_file, 'r') as handle:
                handle.seek(self._conv_offset); chunk = handle.read(); self._conv_offset = handle.tell()
            if not chunk: return
            text = self._conv_partial + chunk; lines = text.split('\n'); self._conv_partial = lines.pop(); added=False
            for line in lines:
                line=line.strip()
                if not line: continue
                if line.startswith('# EVENT'):
                    fields=line.split('\t')
                    if len(fields)>=3:
                        try: event_iter=int(fields[1])
                        except ValueError: continue
                        label=fields[2]; metadata=' '.join(fields[3:])
                        if metadata: label='%s (%s)' % (label, metadata)
                        labels=self._conv_events.setdefault(event_iter, [])
                        if label not in labels: labels.append(label)
                        self._update_decon_progress(event_iter, label)
                        added=True
                    continue
                if line.startswith('#'): continue
                fields=line.split()
                if len(fields)<2: continue
                try: iteration=int(fields[0]); tack=float(fields[1])
                except ValueError: continue
                self._conv_x.append(iteration); self._conv_y.append(tack)
                self._update_decon_progress(iteration)
                if len(self._conv_y)>1 and tack != 0.0:
                    self._conv_conv.append(abs(tack-self._conv_y[-2])/abs(tack))
                else: self._conv_conv.append(float('nan'))
                added=True
            if added: self._redraw_convergence()
        except (OSError, IOError): return
        except Exception: logging.exception('Could not update progress plot')

    def _redraw_convergence(self):
        if not self._conv_x: return
        import math
        import re

        # Always show the complete protocol.  Adaptive smoothing is stage-local,
        # so violent exploratory motion does not obscure later refinement while
        # the raw .conv data and relative-change trace remain untouched.
        xx=self._conv_x[:]; yy=self._conv_y[:]; cc=self._conv_conv[:]
        display_y=list(yy); smoothing=[]
        for a,b in self._stage_ranges(0, len(self._conv_x)):
            smoothed, window, score=self._adaptive_smooth(yy[a:b])
            display_y[a:b]=smoothed
            if window > 1:
                smoothing.append((self._conv_x[a], self._conv_x[b-1], window, score))

        step=max(1,(len(xx)+4999)//5000)
        x=xx[::step]; y=[v/1e12 for v in display_y[::step]]; c=cc[::step]
        if xx and x[-1] != xx[-1]:
            x.append(xx[-1]); y.append(display_y[-1]/1e12); c.append(cc[-1])
        self._conv_line.set_data(x,y); self._conv_line2.set_data(x,c)
        self._conv_axes.set_xlim(x[0], max(x[-1], x[0]+1))

        # Full range means the complete smoothed system-state range, not a
        # percentile crop.  A small pad keeps extrema away from the frame.
        lo=min(y); hi=max(y); pad=max((hi-lo)*.05, abs(hi)*1e-6, 1e-6)
        self._conv_axes.set_ylim(lo-pad,hi+pad)

        finite_c=[v for v in c if v==v and v>0 and abs(v)!=float('inf')]
        if finite_c:
            logs=[math.log10(v) for v in finite_c]
            llim=self._mad_limits(logs, sigma=6.0)
            if llim: self._conv_axes2.set_ylim(10**llim[0],10**llim[1])

        # Extract sparse peak counts from event metadata, then hold the most
        # recent count between events.  This produces an interpretable step-like
        # model-complexity trace without adding any C++ I/O.
        peak_events=[]
        for iteration in sorted(self._conv_events):
            for label in self._conv_events[iteration]:
                m=re.search(r'\bpeaks=(\d+)', label)
                if m:
                    peak_events.append((iteration, int(m.group(1))))
        if peak_events and self._peak_line is not None:
            px=[]; py=[]
            current=None; ei=0
            for iteration in x:
                while ei < len(peak_events) and peak_events[ei][0] <= iteration:
                    current=peak_events[ei][1]; ei += 1
                if current is not None:
                    px.append(iteration); py.append(current)
            self._peak_line.set_data(px,py)
            if py:
                plo=min(py); phi=max(py); ppad=max(1.0,(phi-plo)*.08)
                self._peak_axes.set_ylim(max(0,plo-ppad),phi+ppad)
        elif self._peak_line is not None:
            self._peak_line.set_data([],[])

        if smoothing:
            wins=sorted(set(w for _,_,w,_ in smoothing))
            smooth_text=(str(wins[0]) if len(wins)==1 else ', '.join(map(str,wins)))
            self._conv_axes.set_title('System state - all stages (adaptive moving average: %s)' % smooth_text,
                                      fontsize=getattr(self, '_plot_fontsize', 8))
        else:
            self._conv_axes.set_title('System state - all stages (raw; stable)', fontsize=getattr(self, '_plot_fontsize', 8))

        for artist in self._conv_event_artists:
            try: artist.remove()
            except Exception: pass
        self._conv_event_artists=[]

        # Structural/model changes only.  DoRun start/end events are used to
        # discover mode transitions, avoiding a forest of redundant labels.
        annotations=[]
        last_mode=None
        for iteration in sorted(self._conv_events):
            labels=self._conv_events[iteration]
            mode_here=None
            for label in labels:
                mm=re.search(r'\bmode=(\d+)', label)
                if mm: mode_here=mm.group(1)
            for label in labels:
                low=label.lower()
                name=None
                if 'squash' in low: name='Squash'
                elif 'cull' in low: name='Cull'
                elif 'setpeaks' in low or 'peak-shape' in low: name='Peak contraction'
                if name:
                    annotations.append((iteration, name, mode_here))
            if mode_here is not None and mode_here != last_mode:
                # The first mode establishes context; later values are actual
                # protocol transitions worth marking.
                if last_mode is not None:
                    annotations.append((iteration, 'Mode change', mode_here))
                last_mode=mode_here

        # De-duplicate coincident semantic events while preserving order.
        clean=[]; seen=set()
        for item in annotations:
            key=(item[0],item[1],item[2])
            if key not in seen:
                clean.append(item); seen.add(key)

        # Matplotlib's categorical cycle provides visually distinct colours.
        # Lines use the same colour on both plots; text exists only on top.
        cycle=__import__('matplotlib').rcParams['axes.prop_cycle'].by_key().get('color', ['C0'])
        # Keep the protocol labels in a compact key down the left-hand side of
        # the System State panel.  Their matching coloured vertical lines still
        # mark the actual iteration on both plots, so labels no longer obscure
        # the traces or depend on where the transition occurs horizontally.
        level=0
        for idx,(iteration,name,mode) in enumerate(clean):
            if iteration < x[0] or iteration > x[-1]: continue
            color=cycle[idx % len(cycle)]
            for ax in (self._conv_axes,self._conv_axes2):
                self._conv_event_artists.append(ax.axvline(iteration, linestyle='--', linewidth=1.0,
                                                           alpha=.78, color=color))
            text=name + ((' - mode=%s' % mode) if mode is not None else '')
            txt=self._conv_axes.text(0.012, 0.985-level*0.105, text,
                                     rotation=0, va='top', ha='left', fontsize=7.5,
                                     color=color, transform=self._conv_axes.transAxes,
                                     clip_on=True,
                                     bbox=dict(facecolor='white', edgecolor='none', alpha=.62, pad=1.0))
            self._conv_event_artists.append(txt)
            level = (level + 1) % 8

        self._conv_canvas.draw_idle()

    def finish_convergence(self):
        """Read the final flushed data and stop polling once the process ends."""
        try:
            self._on_convergence_timer()
        finally:
            if self._conv_timer is not None:
                self._conv_timer.Stop()


    def set_workflow(self, steps, current=0):
        """Configure a small coarse-grained progress display."""
        self._workflow_steps = [str(x) for x in (steps or []) if str(x)]
        self._workflow_done = set()
        self.progressGauge.SetRange(max(1, len(self._workflow_steps)))
        self._progress_anim_phase = 0
        if not self._progress_anim_timer.IsRunning():
            self._progress_anim_timer.Start(700)
        self.start_step(current)

    def start_step(self, step):
        if not self._workflow_steps:
            return
        if isinstance(step, str):
            try:
                idx = self._workflow_steps.index(step)
            except ValueError:
                return
        else:
            idx = max(0, min(int(step), len(self._workflow_steps)-1))
        if self._workflow_index >= 0 and idx > self._workflow_index:
            self._workflow_done.update(range(self._workflow_index, idx))
        self._workflow_index = idx
        self.progressGauge.SetValue(min(idx, len(self._workflow_steps)))
        self.progressTitle.SetLabel('Step %d of %d — %s' % (idx+1, len(self._workflow_steps), self._workflow_steps[idx]))
        self.progressPanel.Layout()
        self.set_status(self._workflow_steps[idx] + '...')

    def finish_workflow(self, success=True):
        if not self._workflow_steps:
            return
        if success:
            self.progressGauge.SetValue(len(self._workflow_steps))
            self.progressTitle.SetLabel('Complete — %s' % self._workflow_steps[-1])
            self.progressHint.SetLabel('Calculation complete. Full text output is available below.')
        else:
            self.progressTitle.SetLabel('Calculation stopped — inspect the output below')
            self.progressHint.SetLabel('The text output may contain details of the problem.')
        if self._progress_anim_timer.IsRunning():
            self._progress_anim_timer.Stop()
        self.progressGauge.SetRunning(False)
        self.progressPanel.Layout()

    def _on_progress_animation(self, event=None):
        # Deliberately slow and tiny: visual reassurance, not fake numerical
        # progress.  Refreshing one small control every 700 ms is negligible
        # compared with the external NMR calculations.
        self._progress_anim_phase = (self._progress_anim_phase + 1) % 8
        self.progressGauge.SetAnimationPhase(self._progress_anim_phase)

    def observe_output(self, text):
        """Use stable NMRPipe messages for useful detail without parsing percentages."""
        low = str(text).lower()
        if self._decon_profile:
            import re
            if 'restricted reconstruction' in low:
                self.progressTitle.SetLabel('Preparing restricted reconstruction')
            elif 'first run of unidec' in low:
                self.progressTitle.SetLabel('Preparing initial 1D optimisation')
            elif 'second run of unidec' in low:
                self.progressTitle.SetLabel('Preparing refined 1D optimisation')
            elif 'contracting peak shape' in low:
                self.progressTitle.SetLabel('Updating peak shape')
            elif 'squash function' in low:
                self.progressTitle.SetLabel('Refining candidate peaks')
            elif 'fit fixed radii' in low or 'fit radii' in low:
                self.progressGauge.SetValue(max(self.progressGauge._value, 910))
                self.progressTitle.SetLabel('Fitting peak positions and widths')
            elif 'writing to' in low or 'written correlate' in low or 'final cross peaks' in low:
                self.progressGauge.SetValue(max(self.progressGauge._value, 930))
                self.progressTitle.SetLabel('Generating outputs')
            elif 'exiting cleanly' in low:
                self.progressGauge.SetValue(max(self.progressGauge._value, 970))
                self.progressTitle.SetLabel('Finalising results')
            # DoRun prints local iteration values at iterShow intervals.  Use
            # these immediately between the one-second .conv tail updates.
            m=re.match(r'^\s*(\d+)\s+[-+0-9.eE]+\s+[-+0-9.eE]+(?:\s+[-+0-9.eE]+)?\s*$', str(text))
            if m and self._decon_run_active:
                idx=self._decon_run_index; runs=self._decon_profile.get('_run_max', [])
                weights=self._decon_profile.get('_run_weights', [])
                if idx < len(runs):
                    local=int(m.group(1)); base=.08+sum(weights[:idx])
                    frac=min(.90, base+weights[idx]*min(1.0,local/float(runs[idx])))
                    self.progressGauge.SetValue(int(round(frac*1000)))
                    self.progressHint.SetLabel('Iteration %d of up to %d in this optimisation.' % (local,runs[idx]))
        if 'smile' in low and not self._smile_announced:
            self._smile_announced = True
            self.append_text('\n--- SMILE reconstruction has begun ---\n')
            if 'SMILE reconstruction' in self._workflow_steps:
                self.start_step('SMILE reconstruction')
        # Dimension messages are useful to power users and make the status bar
        # feel live without adding extra visual clutter.
        if 'processing z dimension' in low:
            self.set_status('Processing Z dimension...')
        elif 'processing a dimension' in low:
            self.set_status('Processing A dimension...')
        elif 'processing y dimension' in low:
            self.set_status('Processing Y dimension...')
        elif 'processing xy dimensions' in low:
            self.set_status('Processing XY dimensions...')

    def append_text(self, text: str):
        if not text:
            return
        try:
            self.output.AppendText(text)
            self.output.ShowPosition(self.output.GetLastPosition())
        except Exception:
            pass

    def set_status(self, text: str):
        try:
            self.statusbar.SetStatusText(text)
        except Exception:
            pass

    def _on_details(self, event=None):
        if self.detailsBtn is None:
            return
        showing = not self.output.IsShown()
        self.output.Show(showing)
        self.detailsBtn.SetLabel('Hide details' if showing else 'Show details')
        try:
            w = self.GetSize().GetWidth()
            self.SetSize((w, self._details_height if showing else min(220, self._details_height)))
        except Exception:
            pass
        self.Layout()

    def _on_clear(self, event):
        try:
            self.output.Clear()
        except Exception:
            pass

    def _on_close(self, event):
        try:
            if self._conv_timer is not None:
                self._conv_timer.Stop()
            if self._progress_anim_timer.IsRunning():
                self._progress_anim_timer.Stop()
            self.Destroy()
        except Exception:
            pass


def _as_command(command: Union[str, Sequence[str]]) -> Sequence[str]:
    if isinstance(command, (list, tuple)):
        return list(command)
    return [str(command)]


def run_command_with_output(
    command: Union[str, Sequence[str]],
    parent=None,
    title: str = 'Shell Output',
    cwd: Optional[str] = None,
    on_finish=None,
    output_frame: Optional[ShellOutputFrame] = None,
    final: bool = True,
    label: Optional[str] = None,
    convergence_file: Optional[str] = None,
    workflow_steps=None,
    workflow_step=None,
    completion_text: Optional[str] = None,
    decon_profile=None,
) -> ShellOutputFrame:
    """Run a command in the background and stream combined stdout/stderr."""

    frame = output_frame or ShellOutputFrame(parent, title=title, convergence_file=convergence_file, decon_profile=decon_profile)
    if workflow_steps is not None:
        frame.set_workflow(workflow_steps, 0 if workflow_step is None else workflow_step)
    elif workflow_step is not None:
        frame.start_step(workflow_step)
    try:
        if isinstance(command, (list, tuple)):
            display_cmd = ' '.join(str(x) for x in command)
        else:
            display_cmd = str(command)
        frame.set_status('Running: ' + (label or display_cmd))
        if label:
            frame.append_text('\n=== ' + label + ' ===\n')
    except Exception:
        pass
    if output_frame is None:
        frame.Show()

    def _invoke_finish(callback, rc):
        try:
            callback(rc)
        except TypeError:
            callback()

    def finalize(rc: int):
        try:
            frame.finish_convergence()
            if rc != 0:
                frame.append_text(f'\n[command failed with exit code {rc}]\n')
            elif final and completion_text:
                frame.append_text('\n' + completion_text.rstrip() + '\n')
            if final:
                if convergence_file and frame._decon_profile:
                    frame.progressGauge.SetValue(1000 if rc == 0 else frame.progressGauge._value)
                    frame.progressTitle.SetLabel('Complete — outputs generated' if rc == 0 else 'Calculation stopped — inspect output')
                    frame.progressHint.SetLabel('Deconvolution complete.' if rc == 0 else 'The shell output may contain details of the problem.')
                    frame.progressGauge.SetRunning(False)
                    if frame._progress_anim_timer.IsRunning(): frame._progress_anim_timer.Stop()
                else:
                    frame.finish_workflow(success=(rc == 0))
            frame.set_status((f'Complete (exit code {rc})' if final else 'Continuing...'))
        except Exception:
            pass
        if on_finish is not None:
            try:
                wx.CallAfter(_invoke_finish, on_finish, rc)
            except Exception:
                logging.exception('Shell output completion callback failed')

    def worker():
        try:
            if isinstance(command, (list, tuple)):
                popen_args = dict(args=_as_command(command), shell=False)
            else:
                popen_args = dict(args=str(command), shell=True)
            proc = subprocess.Popen(
                **popen_args,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                cwd=cwd,
                text=True,
                bufsize=1,
            )
            assert proc.stdout is not None
            for chunk in iter(proc.stdout.readline, ''):
                if not chunk:
                    break
                wx.CallAfter(frame.append_text, chunk)
                wx.CallAfter(frame.observe_output, chunk)
            rc = proc.wait()
            wx.CallAfter(finalize, rc)
        except Exception as exc:
            logging.exception('Could not run shell command')
            wx.CallAfter(frame.append_text, f'\n[error] {exc}\n')
            wx.CallAfter(frame.set_status, 'Failed')

    threading.Thread(target=worker, daemon=True).start()
    return frame
