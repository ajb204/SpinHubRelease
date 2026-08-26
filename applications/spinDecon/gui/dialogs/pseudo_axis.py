"""GUI presentation for pseudo-axis tables."""

def show_pseudo_axis_table(parent, table, title="Pseudo-axis table"):
    """Show a pseudo-axis table in a small native wxPython viewer."""
    import wx

    frame = wx.Frame(parent, title=title, style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
    panel = wx.Panel(frame)
    outer = wx.BoxSizer(wx.VERTICAL)

    listing = wx.ListCtrl(panel, style=wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.BORDER_SUNKEN)
    for col_index, header in enumerate(table.headers):
        listing.InsertColumn(col_index, str(header))

    for row_index, row in enumerate(table.rows):
        first = str(row.get(table.headers[0], '') or '')
        item = listing.InsertItem(row_index, first)
        for col_index, header in enumerate(table.headers[1:], start=1):
            listing.SetItem(item, col_index, str(row.get(header, '') or ''))

    # Size columns to their contents, with modest bounds so large tables remain usable.
    for col_index in range(len(table.headers)):
        listing.SetColumnWidth(col_index, wx.LIST_AUTOSIZE_USEHEADER)
        header_width = listing.GetColumnWidth(col_index)
        listing.SetColumnWidth(col_index, wx.LIST_AUTOSIZE)
        content_width = listing.GetColumnWidth(col_index)
        listing.SetColumnWidth(col_index, max(70, min(220, max(header_width, content_width))))

    close_button = wx.Button(panel, wx.ID_CLOSE, label="Close")
    close_button.Bind(wx.EVT_BUTTON, lambda event: frame.Close())

    outer.Add(listing, 1, wx.EXPAND | wx.ALL, 8)
    outer.Add(close_button, 0, wx.ALIGN_RIGHT | wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
    panel.SetSizer(outer)

    width = max(440, min(850, sum(listing.GetColumnWidth(i) for i in range(len(table.headers))) + 40))
    height = max(300, min(600, 120 + 24 * min(len(table.rows), 18)))
    frame.SetClientSize((width, height))
    frame.CentreOnParent()
    frame.Show()
    return frame
