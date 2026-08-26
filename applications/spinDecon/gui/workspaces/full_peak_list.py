import os
import re
import wx
import math

from spinDecon.processing.peak_list_operations import alias_peak_coordinate, transpose_2d_peaks


class ReferenceAliasFrame(wx.Frame):
    """Modeless alias controls for the currently selected reference peak."""

    def __init__(self, viewer):
        super().__init__(viewer, title=('Alias Full Peak' if viewer.mode == 'full' else 'Alias Reference Peak'), style=wx.DEFAULT_FRAME_STYLE & ~wx.RESIZE_BORDER)
        self.viewer = viewer
        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)
        self.axis_rows = []
        for meta in (self.viewer.controller.get_full_peak_axis_metadata() if self.viewer.mode == 'full' else self.viewer.controller.get_reference_peak_axis_metadata()):
            row = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(panel, label=str(meta['label']), size=(80, -1))
            minus = wx.Button(panel, label='-', size=(42, -1))
            plus = wx.Button(panel, label='+', size=(42, -1))
            minus.Bind(wx.EVT_BUTTON, lambda evt, axis=meta['axis']: self.viewer.alias_selected(axis, -1))
            plus.Bind(wx.EVT_BUTTON, lambda evt, axis=meta['axis']: self.viewer.alias_selected(axis, +1))
            row.Add(label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
            row.Add(minus, 0, wx.RIGHT, 5)
            row.Add(plus, 0)
            outer.Add(row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 10)
            self.axis_rows.append((label, minus, plus))
        close = wx.Button(panel, label='Close')
        close.Bind(wx.EVT_BUTTON, lambda evt: self.Close())
        outer.Add(close, 0, wx.ALIGN_RIGHT | wx.ALL, 10)
        panel.SetSizer(outer)
        outer.Fit(self)
        self.CentreOnParent()


class SaveReferencePeakListDialog(wx.Dialog):
    """Choose a reference peak-list destination below the project's SpecPath."""

    def __init__(self, viewer):
        super().__init__(viewer, title='Save Reference Peak List')
        self.viewer = viewer
        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)
        outer.Add(wx.StaticText(panel, label='File name (relative to SpecPath):'),
                  0, wx.LEFT | wx.RIGHT | wx.TOP, 10)

        row = wx.BoxSizer(wx.HORIZONTAL)
        current = str(viewer.controller.referencePeakBox.GetValue() or '').strip()
        self.path_box = wx.TextCtrl(panel, value=current, size=(390, -1))
        browse = wx.Button(panel, label='Browse...')
        browse.Bind(wx.EVT_BUTTON, self.on_browse)
        row.Add(self.path_box, 1, wx.EXPAND | wx.RIGHT, 6)
        row.Add(browse, 0)
        outer.Add(row, 0, wx.EXPAND | wx.ALL, 10)

        # Create the standard dialog buttons with ``panel`` as their parent.
        # ``Dialog.CreateStdDialogButtonSizer`` parents its buttons to the
        # dialog itself; putting that sizer inside a wx.Panel sizer triggers
        # wxWidgets' parent/sizer assertion (notably with wxPython 4.3).
        buttons = wx.StdDialogButtonSizer()
        ok_button = wx.Button(panel, wx.ID_OK)
        cancel_button = wx.Button(panel, wx.ID_CANCEL)
        buttons.AddButton(ok_button)
        buttons.AddButton(cancel_button)
        buttons.Realize()
        outer.Add(buttons, 0, wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        panel.SetSizer(outer)
        outer.Fit(self)
        self.CentreOnParent()

    def on_browse(self, event):
        controller = self.viewer.controller
        controller._sync_directory_state_only()
        spec_dir = os.path.abspath(controller.state.spec_dir())
        current = str(self.path_box.GetValue() or '').strip()
        default_dir = spec_dir
        default_file = ''
        if current:
            candidate = current if os.path.isabs(current) else os.path.join(spec_dir, current)
            default_dir = os.path.dirname(candidate) or spec_dir
            default_file = os.path.basename(candidate)
            if not os.path.isdir(default_dir):
                default_dir = spec_dir
        with wx.FileDialog(
            self, 'Select reference peak-list file', defaultDir=default_dir,
            defaultFile=default_file, wildcard='Peak list files (*.list)|*.list|All files (*.*)|*.*',
            style=wx.FD_SAVE
        ) as picker:
            if picker.ShowModal() != wx.ID_OK:
                return
            try:
                relative, _ = controller.reference_peak_save_destination(picker.GetPath())
            except ValueError as exc:
                wx.MessageBox(str(exc), 'Invalid save location', wx.OK | wx.ICON_ERROR, self)
                return
            self.path_box.SetValue(relative)

    def value(self):
        return str(self.path_box.GetValue() or '').strip()


class SaveFullPeakListDialog(wx.Dialog):
    """Choose a destination for the authoritative Full nD peak list."""

    def __init__(self, viewer):
        super().__init__(viewer, title='Save Full Peak List')
        self.viewer = viewer
        panel = wx.Panel(self)
        outer = wx.BoxSizer(wx.VERTICAL)
        outer.Add(wx.StaticText(panel, label='File name (relative to SpecPath):'),
                  0, wx.LEFT | wx.RIGHT | wx.TOP, 10)
        row = wx.BoxSizer(wx.HORIZONTAL)
        current = str(viewer.controller.fullPeakBox.GetValue() or '').strip()
        self.path_box = wx.TextCtrl(panel, value=current, size=(390, -1))
        browse = wx.Button(panel, label='Browse...')
        browse.Bind(wx.EVT_BUTTON, self.on_browse)
        row.Add(self.path_box, 1, wx.EXPAND | wx.RIGHT, 6)
        row.Add(browse, 0)
        outer.Add(row, 0, wx.EXPAND | wx.ALL, 10)
        buttons = wx.StdDialogButtonSizer()
        ok_button = wx.Button(panel, wx.ID_OK)
        cancel_button = wx.Button(panel, wx.ID_CANCEL)
        buttons.AddButton(ok_button)
        buttons.AddButton(cancel_button)
        buttons.Realize()
        outer.Add(buttons, 0, wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        panel.SetSizer(outer)
        outer.Fit(self)
        self.CentreOnParent()

    def on_browse(self, event):
        controller = self.viewer.controller
        controller._sync_directory_state_only()
        spec_dir = os.path.abspath(controller.state.spec_dir())
        current = str(self.path_box.GetValue() or '').strip()
        default_dir, default_file = spec_dir, ''
        if current:
            candidate = current if os.path.isabs(current) else os.path.join(spec_dir, current)
            default_dir = os.path.dirname(candidate) or spec_dir
            default_file = os.path.basename(candidate)
            if not os.path.isdir(default_dir):
                default_dir = spec_dir
        with wx.FileDialog(
            self, 'Select full peak-list file', defaultDir=default_dir,
            defaultFile=default_file,
            wildcard='Peak list files (*.list)|*.list|All files (*.*)|*.*',
            style=wx.FD_SAVE
        ) as picker:
            if picker.ShowModal() != wx.ID_OK:
                return
            try:
                relative, _ = controller.full_peak_save_destination(picker.GetPath())
            except ValueError as exc:
                wx.MessageBox(str(exc), 'Invalid save location', wx.OK | wx.ICON_ERROR, self)
                return
            self.path_box.SetValue(relative)

    def value(self):
        return str(self.path_box.GetValue() or '').strip()


class PeakListFrame(wx.Frame):
    """Datastore-backed peak-list viewer."""

    def __init__(self, controller, mode='reference'):
        title = 'Reference 2D Peak List' if mode == 'reference' else 'Full Peak List'
        super().__init__(controller, title=title, size=(760, 360))
        self.controller = controller
        self.mode = mode
        self.alias_frame = None
        # Canonical full-peak names for the rows currently visible in the list.
        # This keeps selection independent of presentation/filtering.
        self._displayed_full_names = []
        self._sort_column = None
        self._sort_ascending = True
        panel = wx.Panel(self)
        self.list = wx.ListCtrl(panel, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        buttons = wx.BoxSizer(wx.VERTICAL)
        self.refresh_button = wx.Button(panel, label='Refresh')
        buttons.Add(self.refresh_button, 0, wx.EXPAND | wx.BOTTOM, 5)
        self.show_button = wx.Button(panel, label='Show')
        buttons.Add(self.show_button, 0, wx.EXPAND | wx.BOTTOM, 5)
        self.show_button.Bind(wx.EVT_BUTTON, self.on_show)
        if mode == 'reference':
            self.classify_button = wx.Button(panel, label='Classify')
            buttons.Add(self.classify_button, 0, wx.EXPAND | wx.BOTTOM, 5)
            self.classify_button.Bind(wx.EVT_BUTTON, self.on_classify)
            self.remove_button = wx.Button(panel, label='Remove')
            buttons.Add(self.remove_button, 0, wx.EXPAND | wx.BOTTOM, 5)
            self.remove_button.Bind(wx.EVT_BUTTON, self.on_remove)
            self.transpose_button = wx.Button(panel, label='Transpose')
            buttons.Add(self.transpose_button, 0, wx.EXPAND | wx.BOTTOM, 5)
            self.transpose_button.Bind(wx.EVT_BUTTON, self.on_transpose)
            self.alias_button = wx.Button(panel, label='Alias')
            buttons.Add(self.alias_button, 0, wx.EXPAND | wx.BOTTOM, 5)
            self.alias_button.Bind(wx.EVT_BUTTON, self.on_alias)
            self.save_button = wx.Button(panel, label='Save')
            buttons.Add(self.save_button, 0, wx.EXPAND | wx.BOTTOM, 5)
            self.save_button.Bind(wx.EVT_BUTTON, self.on_save)
        else:
            self.snr_button = wx.ToggleButton(panel, label='SNR')
            self.snr_button.SetToolTip('Toggle the Intensity column between authoritative intensity and signal-to-noise display')
            self.snr_button.Bind(wx.EVT_TOGGLEBUTTON, self.on_snr_toggle)
            buttons.Add(self.snr_button, 0, wx.EXPAND | wx.BOTTOM, 5)
            self.alias_button = wx.Button(panel, label='Alias')
            buttons.Add(self.alias_button, 0, wx.EXPAND | wx.BOTTOM, 5)
            self.alias_button.Bind(wx.EVT_BUTTON, self.on_alias)
            self.save_button = wx.Button(panel, label='Save')
            buttons.Add(self.save_button, 0, wx.EXPAND | wx.BOTTOM, 5)
            self.save_button.Bind(wx.EVT_BUTTON, self.on_save)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_show)
        if mode == 'full':
            self.list.Bind(wx.EVT_LIST_COL_CLICK, self.on_column_click)
            viewers = getattr(controller, '_full_peak_list_viewers', None)
            if viewers is None:
                viewers = []
                controller._full_peak_list_viewers = viewers
            viewers.append(self)
        self.close_button = wx.Button(panel, label='Close')
        buttons.Add(self.close_button, 0, wx.EXPAND)
        self.refresh_button.Bind(wx.EVT_BUTTON, self.on_refresh)
        self.close_button.Bind(wx.EVT_BUTTON, lambda evt: self.Close())
        self.Bind(wx.EVT_CLOSE, self.on_close)

        list_area = wx.BoxSizer(wx.VERTICAL)
        if mode == 'full':
            search_row = wx.BoxSizer(wx.HORIZONTAL)
            self.search_box = wx.SearchCtrl(panel, style=wx.TE_PROCESS_ENTER)
            self.search_box.ShowSearchButton(True)
            self.search_box.ShowCancelButton(True)
            self.search_box.SetDescriptiveText('Search Name (* and ? supported)')
            search_row.Add(self.search_box, 1, wx.EXPAND)
            list_area.Add(search_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 8)
            self.search_box.Bind(wx.EVT_TEXT_ENTER, self.on_search)
            self.search_box.Bind(wx.EVT_SEARCHCTRL_SEARCH_BTN, self.on_search)
            self.search_box.Bind(wx.EVT_SEARCHCTRL_CANCEL_BTN, self.on_cancel_search)
        list_area.Add(self.list, 1, wx.EXPAND | wx.ALL, 8)

        layout = wx.BoxSizer(wx.HORIZONTAL)
        layout.Add(list_area, 1, wx.EXPAND)
        layout.Add(buttons, 0, wx.TOP | wx.RIGHT, 8)
        panel.SetSizer(layout)
        self.on_refresh(None)
        self.CentreOnParent()
        self.Show(True)

    def _selected_name(self):
        index = self.list.GetFirstSelected()
        if index < 0:
            return None
        if self.mode == 'reference':
            return self.list.GetItemText(index, 1)
        if index < len(self._displayed_full_names):
            return self._displayed_full_names[index]
        return None

    def _selected_reference_peak_index(self):
        if self.mode != 'reference':
            return None
        index = self.list.GetFirstSelected()
        if index < 0 or index >= len(self.controller.get_reference_peaks()):
            return None
        return index

    def _selected_reference_peak(self):
        index = self._selected_reference_peak_index()
        if index is None:
            return None
        return self.controller.get_reference_peaks()[index]

    def _restore_selection(self, name):
        if not name:
            return
        for row in range(self.list.GetItemCount()):
            if self.mode == 'reference':
                row_name = self.list.GetItemText(row, 1)
            else:
                row_name = self._displayed_full_names[row] if row < len(self._displayed_full_names) else None
            if row_name == name:
                self.list.Select(row)
                self.list.EnsureVisible(row)
                return

    @staticmethod
    def _longest_integer(text):
        """Return the digit run containing the most digits (first wins ties)."""
        matches = re.findall(r'\d+', str(text))
        if not matches:
            return None
        return max(enumerate(matches), key=lambda item: (len(item[1]), -item[0]))[1]

    @classmethod
    def _split_2d_generated_name(cls, name):
        """Extract ``(nResID, Number)`` from a 2D-generated 3D peak name.

        The final underscore is the structural separator.  ``nResID`` keeps
        the complete string before that separator (including any characters
        or embedded underscores), while ``Number`` is the integer with the
        most digits in the final component.
        """
        text = str(name)
        left, separator, right = text.rpartition('_')
        if not separator or not left:
            return None
        number = cls._longest_integer(right)
        if number is None:
            return None
        return left, number

    def _uses_split_full_names(self, rows):
        return (int(getattr(self.controller, 'dim', 0) or 0) == 3 and bool(rows) and
                all(fields and self._split_2d_generated_name(fields[0]) is not None for fields in rows))

    @staticmethod
    def _wildcard_match(pattern, value):
        """Case-insensitive whole-field glob match supporting only ``?`` and ``*``.

        ``?`` matches exactly one character and ``*`` matches zero or more
        characters.  Every other character is interpreted literally.
        """
        pattern = str(pattern)
        value = str(value)
        regex_parts = []
        for char in pattern:
            if char == '*':
                regex_parts.append('.*')
            elif char == '?':
                regex_parts.append('.')
            else:
                regex_parts.append(re.escape(char))
        regex = '^' + ''.join(regex_parts) + '$'
        return re.match(regex, value, flags=re.IGNORECASE | re.DOTALL) is not None


    @staticmethod
    def _scientific_notation(value):
        """Format a numeric value as A x 10^B using Unicode superscripts."""
        try:
            number = float(value)
        except (TypeError, ValueError):
            return str(value)
        if number == 0:
            return '0'
        import math
        exponent = int(math.floor(math.log10(abs(number))))
        mantissa = number / (10.0 ** exponent)
        supers = str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻')
        return ('%.4g × 10%s' % (mantissa, str(exponent).translate(supers)))

    @classmethod
    def _identifier_sort_key(cls, value):
        text = str(value)
        digit = cls._longest_integer(text)
        if digit is None:
            return (1, 0, text.casefold())
        return (0, int(digit), text.casefold())

    @staticmethod
    def _numeric_sort_key(value):
        try:
            return (0, float(value))
        except (TypeError, ValueError):
            return (1, str(value).casefold())

    def on_column_click(self, event):
        column = int(event.GetColumn())
        if self._sort_column == column:
            self._sort_ascending = not self._sort_ascending
        else:
            self._sort_column = column
            self._sort_ascending = True
        self.on_refresh(None)

    def on_snr_toggle(self, event):
        # Presentation-only switch.  The Full peak payload/store remains in
        # authoritative intensity units.
        self.on_refresh(None)

    def on_search(self, event):
        self.on_refresh(None)

    def on_cancel_search(self, event):
        self.search_box.SetValue('')
        self.on_refresh(None)
        self.search_box.SetFocus()

    def on_refresh(self, event):
        selected = self._selected_name() if self.list.GetColumnCount() else None
        self.list.ClearAll()
        if self.mode == 'reference':
            headings = self.controller.get_reference_peak_headers()
            for col, heading in enumerate(headings): self.list.InsertColumn(col, heading)
            for i, pk in enumerate(self.controller.get_reference_peaks()):
                match = re.findall(r'[0-9]+', pk.name)
                values = [match[0] if match else str(i + 1), pk.name, str(pk.x), str(pk.y)]
                row = self.list.InsertItem(self.list.GetItemCount(), values[0])
                for col, value in enumerate(values[1:], 1): self.list.SetItem(row, col, value)
        else:
            payload = self.controller.get_full_peak_payload()
            rows = payload.get('rows') or payload.get('peaks') or []
            split_names = self._uses_split_full_names(rows)
            self.search_box.SetDescriptiveText(
                'Search nResID (* and ? supported)' if split_names
                else 'Search Name (* and ? supported)')
            width = max([len(row) for row in rows], default=0)
            headings = self.controller.get_full_peak_headers(row_width=width or (int(getattr(self.controller, 'dim', 0) or 0) + 2))
            if split_names and headings:
                # Replace the displayed Name with the two derived identifiers;
                # the original canonical name remains untouched in memory.
                headings = ['nResID', 'Number'] + headings[1:]
            snr_display = bool(getattr(self, 'snr_button', None) and self.snr_button.GetValue())
            intensity_col = next((i for i, h in enumerate(headings)
                                  if 'intensity' in str(h).lower()), None)
            if snr_display and intensity_col is not None:
                headings[intensity_col] = 'SNR'
            for col, heading in enumerate(headings): self.list.InsertColumn(col, heading)

            query = str(self.search_box.GetValue() or '').strip()
            display_rows = []
            for fields in rows:
                if not fields:
                    continue
                canonical_name = str(fields[0])
                if split_names:
                    nresid, number = self._split_2d_generated_name(canonical_name)
                    raw_values = [nresid, number] + list(fields[1:])
                else:
                    raw_values = list(fields)
                    nresid = None
                search_value = nresid if split_names else canonical_name
                if query and not self._wildcard_match(query, search_value):
                    continue
                display_rows.append((canonical_name, raw_values))

            # Sorting is purely a view operation.  Numeric scientific data use
            # their authoritative numeric values; identifiers use the longest
            # digit run from Name/nResID as requested.
            if self._sort_column is not None and self._sort_column < len(headings):
                col = self._sort_column
                heading = str(headings[col]).strip().lower()
                if heading in ('name', 'nresid'):
                    key = lambda item: self._identifier_sort_key(item[1][col])
                elif heading == 'number' or 'intensity' in heading or 'ppm' in heading:
                    key = lambda item: self._numeric_sort_key(item[1][col])
                else:
                    key = lambda item: self._numeric_sort_key(item[1][col])
                display_rows.sort(key=key, reverse=not self._sort_ascending)

            self._displayed_full_names = []
            ppm_cols = [i for i, h in enumerate(headings)
                        if 'ppm' in str(h).lower()]
            for canonical_name, raw_values in display_rows:
                values = [str(v) for v in raw_values]
                # Display chemical shifts to three decimal places without
                # changing the authoritative values held in raw_values/store.
                # Use truncation (towards zero), as requested, rather than
                # rounding the underlying value for presentation.
                for ppm_col in ppm_cols:
                    if ppm_col < len(raw_values):
                        try:
                            number = float(raw_values[ppm_col])
                            truncated = math.trunc(number * 1000.0) / 1000.0
                            values[ppm_col] = '%.3f' % truncated
                        except (TypeError, ValueError, OverflowError):
                            pass
                if intensity_col is not None and intensity_col < len(values):
                    if snr_display:
                        sigma = self.controller.get_noise_sigma()
                        if sigma is not None:
                            try:
                                values[intensity_col] = '%.4g' % (float(raw_values[intensity_col]) / sigma)
                            except (TypeError, ValueError, OverflowError):
                                values[intensity_col] = str(raw_values[intensity_col])
                        else:
                            values[intensity_col] = 'N/A'
                    else:
                        values[intensity_col] = self._scientific_notation(raw_values[intensity_col])
                row = self.list.InsertItem(self.list.GetItemCount(), values[0])
                for col, value in enumerate(values[1:], 1):
                    self.list.SetItem(row, col, value)
                self._displayed_full_names.append(canonical_name)
        for col in range(self.list.GetColumnCount()): self.list.SetColumnWidth(col, wx.LIST_AUTOSIZE_USEHEADER)
        self._restore_selection(selected)

    def focus_peak_name(self, name):
        """Select a canonical Full-list name after an external plot selection."""
        if self.mode != 'full':
            return
        self._restore_selection(str(name))

    def on_show(self, event):
        name = self._selected_name()
        if not name: return
        if self.mode == 'reference': self.controller.select_reference_peak(name)
        else: self.controller.select_full_peak(name)

    def on_classify(self, event):
        index = self._selected_reference_peak_index()
        if index is None:
            wx.MessageBox('Select a peak in the Reference 2D Peak List first.',
                          'No peak selected', wx.OK | wx.ICON_ERROR, self)
            return

        peak = self.controller.get_reference_peaks()[index]
        dialog = wx.TextEntryDialog(
            self, 'Enter Residue ID (1 letter):', 'Set residue')
        dialog.CentreOnParent()
        try:
            if dialog.ShowModal() != wx.ID_OK:
                print('Residue classification cancelled.')
                return
            residue_type = dialog.GetValue().strip()
        finally:
            dialog.Destroy()

        if len(residue_type) != 1:
            wx.MessageBox('Enter exactly one residue identifier character.',
                          'Invalid residue ID', wx.OK | wx.ICON_ERROR, self)
            return

        old_name, new_name = self.controller.classify_reference_peak(index, residue_type)
        print('Renamed peak: %s -> %s' % (old_name, new_name))
        self.on_refresh(None)
        self._restore_selection(new_name)

    def on_remove(self, event):
        name = self._selected_name()
        if not name: return
        peaks = [pk for pk in self.controller.get_reference_peaks() if pk.name != name]
        self.controller.set_reference_peaks(peaks)
        self.controller.refresh_reference_peak_views()
        frame = getattr(self.controller, 'peak_frame', None)
        if frame is not None:
            try: frame.draw_figure()
            except Exception: pass
        self.on_refresh(None)

    def on_transpose(self, event):
        peaks = self.controller.get_reference_peaks()
        transpose_2d_peaks(peaks)
        for peak in peaks:
            self.controller.refresh_reference_peak_indices(peak)
        self.controller.set_reference_peaks(peaks)
        self.controller.refresh_reference_peak_views()
        frame = getattr(self.controller, 'peak_frame', None)
        if frame is not None:
            try: frame.draw_figure()
            except Exception: pass
        self.on_refresh(None)


    def on_save(self, event):
        if self.mode == 'full':
            return self._on_save_full(event)
        dialog = SaveReferencePeakListDialog(self)
        try:
            while dialog.ShowModal() == wx.ID_OK:
                value = dialog.value()
                if not value:
                    wx.MessageBox('Enter a file name before saving.', 'No file name',
                                  wx.OK | wx.ICON_ERROR, dialog)
                    continue
                try:
                    relative, destination = self.controller.reference_peak_save_destination(value)
                except ValueError as exc:
                    wx.MessageBox(str(exc), 'Invalid save location', wx.OK | wx.ICON_ERROR, dialog)
                    continue

                if os.path.exists(destination):
                    answer = wx.MessageBox(
                        'The file already exists:\n%s\n\nOverwrite it?' % destination,
                        'Confirm overwrite', wx.YES_NO | wx.NO_DEFAULT | wx.ICON_WARNING, dialog)
                    if answer != wx.YES:
                        continue
                try:
                    self.controller.save_reference_peak_list(relative)
                except (OSError, ValueError) as exc:
                    wx.MessageBox('Could not save the reference peak list:\n%s' % exc,
                                  'Save failed', wx.OK | wx.ICON_ERROR, dialog)
                    continue
                break
        finally:
            dialog.Destroy()

    def _on_save_full(self, event):
        dialog = SaveFullPeakListDialog(self)
        try:
            while dialog.ShowModal() == wx.ID_OK:
                value = dialog.value()
                if not value:
                    wx.MessageBox('Enter a file name before saving.', 'No file name',
                                  wx.OK | wx.ICON_ERROR, dialog)
                    continue
                try:
                    relative, destination = self.controller.full_peak_save_destination(value)
                except ValueError as exc:
                    wx.MessageBox(str(exc), 'Invalid save location', wx.OK | wx.ICON_ERROR, dialog)
                    continue
                if os.path.exists(destination):
                    answer = wx.MessageBox(
                        'The file already exists:\n%s\n\nOverwrite it?' % destination,
                        'Confirm overwrite', wx.YES_NO | wx.NO_DEFAULT | wx.ICON_WARNING, dialog)
                    if answer != wx.YES:
                        continue
                try:
                    self.controller.save_full_peak_list(relative)
                except (OSError, ValueError) as exc:
                    wx.MessageBox('Could not save the full peak list:\n%s' % exc,
                                  'Save failed', wx.OK | wx.ICON_ERROR, dialog)
                    continue
                break
        finally:
            dialog.Destroy()

    def on_alias(self, event):
        if self.mode == 'full':
            if not self._selected_name():
                wx.MessageBox('Select a peak before opening the Alias controls.', 'No peak selected', wx.OK | wx.ICON_ERROR, self)
                return
        elif self._selected_reference_peak() is None:
            wx.MessageBox('Select a peak before opening the Alias controls.', 'No peak selected', wx.OK | wx.ICON_ERROR, self)
            return
        if self.alias_frame is not None:
            try:
                if self.alias_frame.IsShown():
                    self.alias_frame.Raise()
                    return
            except Exception:
                pass
        self.alias_frame = ReferenceAliasFrame(self)
        self.alias_frame.Bind(wx.EVT_CLOSE, self._on_alias_close)
        self.alias_frame.Show(True)

    def _on_alias_close(self, event):
        frame = self.alias_frame
        self.alias_frame = None
        if frame is not None:
            frame.Destroy()

    def alias_selected(self, axis, direction):
        if self.mode == 'full':
            name = self._selected_name()
            if not name:
                wx.MessageBox('Select a peak in the Full Peak List first.', 'No peak selected', wx.OK | wx.ICON_ERROR, self)
                return
            try:
                self.controller.alias_full_peak(name, axis, direction)
            except ValueError as exc:
                wx.MessageBox(str(exc), 'Alias unavailable', wx.OK | wx.ICON_ERROR, self)
                return
            self.on_refresh(None)
            self._restore_selection(name)
            return
        peak = self._selected_reference_peak()
        if peak is None:
            wx.MessageBox('Select a peak in the Reference 2D Peak List first.', 'No peak selected', wx.OK | wx.ICON_ERROR, self)
            return
        meta = next((m for m in self.controller.get_reference_peak_axis_metadata() if m['axis'] == axis), None)
        if meta is None:
            wx.MessageBox('Spectrum metadata for this axis is unavailable.', 'Alias unavailable', wx.OK | wx.ICON_ERROR, self)
            return
        alias_peak_coordinate(peak, axis, direction, meta['width_ppm'])
        self.controller.refresh_reference_peak_indices(peak)
        name = peak.name
        self.controller.set_reference_peaks(self.controller.get_reference_peaks())
        self.controller.refresh_reference_peak_views(selected_name=name)
        frame = getattr(self.controller, 'peak_frame', None)
        if frame is not None:
            try: frame.draw_figure()
            except Exception: pass
        self.on_refresh(None)
        self._restore_selection(name)

    def on_close(self, event):
        if self.mode == 'full':
            viewers = getattr(self.controller, '_full_peak_list_viewers', [])
            try:
                viewers.remove(self)
            except ValueError:
                pass
        if self.alias_frame is not None:
            try: self.alias_frame.Destroy()
            except Exception: pass
            self.alias_frame = None
        self.Destroy()


# FUTURE FULL-PEAK / NOE MIGRATION MARKER
# conn_data is legacy and is not an authoritative peak store.  If NOE/connectivity
# features are restored, port them onto the canonical Full Peak List plus a
# dedicated connection model/service rather than restoring Slice2D peak ownership.
_FUTURE_FULL_PEAK_NOE_LEGACY = r"""
def NOETest(refy, maxy): ...
def AddNOE(refy, maxy, i1, i2): ...
def OnRemove(event): ...
def OnSave(): ...
entry.f4
entry.distScore
entry.distppm
reciprocated
"""
