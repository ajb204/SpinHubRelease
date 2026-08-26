"""Embedded, high-resolution PDF viewer used by decon reports.

PDF rendering is handled directly by PyMuPDF.  This deliberately avoids
``wx.lib.pdfviewer`` and its optional/legacy PDF backends (including the old
``fitz`` import path).
"""
from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path

import pymupdf
import wx
import wx.lib.sized_controls as sc


class _PDFCanvas(wx.ScrolledWindow):
    """A small PyMuPDF-backed PDF canvas with lazy, zoom-aware rendering."""

    GAP = 16

    def __init__(self, parent):
        super().__init__(parent, style=wx.HSCROLL | wx.VSCROLL | wx.SUNKEN_BORDER)
        # AutoBufferedPaintDC requires BG_STYLE_PAINT on wxWidgets (notably macOS).
        # Set it before the first EVT_PAINT can be delivered.
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(wx.Colour(96, 96, 96))
        self.SetScrollRate(20, 20)
        self._document = None
        self._path = ''
        self._zoom = 1.0
        self._page_sizes = []
        self._placements = []
        self._cache = {}
        self._content_scale = self._get_content_scale()
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_SIZE, self._on_size)

    def close(self):
        self._cache.clear()
        if self._document is not None:
            self._document.close()
            self._document = None

    def load_file(self, path):
        self.close()
        self._path = str(Path(path).resolve())
        self._document = pymupdf.open(self._path)
        if self._document.page_count == 0:
            raise ValueError('The PDF contains no pages.')
        self._page_sizes = []
        for page in self._document:
            rect = page.rect
            self._page_sizes.append((float(rect.width), float(rect.height)))
        self._zoom = 1.0
        self._layout()
        self.Scroll(0, 0)
        self.Refresh()

    def set_zoom(self, zoom):
        zoom = max(0.25, min(float(zoom), 5.0))
        if abs(zoom - self._zoom) < 0.001:
            return
        self._zoom = zoom
        self._cache.clear()
        self._layout()
        self.Refresh()

    def zoom_in(self):
        self.set_zoom(self._zoom * 1.25)

    def zoom_out(self):
        self.set_zoom(self._zoom / 1.25)

    def actual_size(self):
        self.set_zoom(1.0)

    def fit_width(self):
        if not self._page_sizes:
            return
        client_width = max(1, self.GetClientSize().width - self.GAP * 2)
        widest = max(width for width, _ in self._page_sizes)
        self.set_zoom(client_width / widest)

    def fit_page(self):
        if not self._page_sizes:
            return
        width, height = self._page_sizes[0]
        client = self.GetClientSize()
        zx = max(1, client.width - self.GAP * 2) / width
        zy = max(1, client.height - self.GAP * 2) / height
        self.set_zoom(min(zx, zy))

    def _layout(self):
        self._placements = []
        if not self._page_sizes:
            self.SetVirtualSize((1, 1))
            return
        max_width = max(int(round(w * self._zoom)) for w, _ in self._page_sizes)
        y = self.GAP
        for width, height in self._page_sizes:
            w = max(1, int(round(width * self._zoom)))
            h = max(1, int(round(height * self._zoom)))
            x = self.GAP + (max_width - w) // 2
            self._placements.append(wx.Rect(x, y, w, h))
            y += h + self.GAP
        self.SetVirtualSize((max_width + self.GAP * 2, y))

    def _get_content_scale(self):
        # wx uses logical pixels while Retina/HiDPI displays have more physical
        # pixels. Render at the physical backing resolution and let wx display
        # the bitmap at its logical size.
        try:
            scale = float(self.GetContentScaleFactor())
        except (AttributeError, TypeError, ValueError):
            scale = 1.0
        return max(1.0, scale)

    def _bitmap(self, page_number):
        content_scale = self._get_content_scale()
        key = (page_number, round(self._zoom, 4), round(content_scale, 3))
        bitmap = self._cache.get(key)
        if bitmap is not None:
            return bitmap
        page = self._document.load_page(page_number)
        # Render at the display's physical backing resolution. On a Retina
        # display content_scale is normally 2, so a page occupying 800 logical
        # pixels is rasterised at about 1600 physical pixels.
        render_scale = self._zoom * content_scale
        pix = page.get_pixmap(matrix=pymupdf.Matrix(render_scale, render_scale),
                              colorspace=pymupdf.csRGB, alpha=False)
        image = wx.Image(pix.width, pix.height, bytes(pix.samples))
        # wxPython versions differ in support for scale-aware Bitmap
        # constructors.  Use the universally supported Image -> Bitmap path;
        # the high-resolution raster is scaled to the logical page rectangle
        # explicitly in the paint handler.
        bitmap = wx.Bitmap(image)
        self._cache[key] = bitmap
        return bitmap

    def _on_paint(self, event):
        dc = wx.AutoBufferedPaintDC(self)
        self.PrepareDC(dc)
        dc.SetBackground(wx.Brush(self.GetBackgroundColour()))
        dc.Clear()
        if self._document is None:
            return

        # Work in virtual coordinates and rasterise only pages intersecting the
        # visible viewport. This keeps large summary PDFs responsive.
        vx, vy = self.GetViewStart()
        sx, sy = self.GetScrollPixelsPerUnit()
        client = self.GetClientSize()
        visible = wx.Rect(vx * sx, vy * sy, client.width, client.height)
        # GraphicsContext.DrawBitmap accepts an explicit logical width and
        # height.  This lets us keep the PyMuPDF raster at Retina/HiDPI
        # resolution while displaying it at the intended wx logical size,
        # without relying on version-specific wx.Bitmap(scale=...) APIs.
        gc = wx.GraphicsContext.Create(dc)
        for page_number, rect in enumerate(self._placements):
            if not rect.Intersects(visible):
                continue
            bitmap = self._bitmap(page_number)
            if gc is not None:
                gc.DrawBitmap(bitmap, rect.x, rect.y, rect.width, rect.height)
            else:
                # Conservative fallback for platforms where a graphics
                # context cannot be created.
                source = wx.MemoryDC()
                source.SelectObject(bitmap)
                dc.StretchBlit(rect.x, rect.y, rect.width, rect.height,
                               source, 0, 0, bitmap.GetWidth(),
                               bitmap.GetHeight())
                source.SelectObject(wx.NullBitmap)

    def _on_size(self, event):
        event.Skip()
        content_scale = self._get_content_scale()
        if abs(content_scale - self._content_scale) > 0.001:
            self._content_scale = content_scale
            self._cache.clear()
        self.Refresh()


class PDFViewer(sc.SizedFrame):
    def __init__(self, parent, **kwds):
        super().__init__(parent, **kwds)
        self.pdf_path = ''
        pane = self.GetContentsPane()

        controls = wx.Panel(pane)
        row = wx.BoxSizer(wx.HORIZONTAL)
        for label, handler in (
            ('Fit Width', self.OnFitWidth),
            ('Fit Page', self.OnFitPage),
            ('100%', self.OnActualSize),
            ('-', self.OnZoomOut),
            ('+', self.OnZoomIn),
        ):
            button = wx.Button(controls, label=label)
            button.Bind(wx.EVT_BUTTON, handler)
            row.Add(button, 0, wx.ALL, 3)
        row.AddStretchSpacer(1)
        controls.SetSizer(row)
        controls.SetSizerProps(expand=True)

        self.viewer = _PDFCanvas(pane)
        self.viewer.SetSizerProps(expand=True, proportion=1)

        actions = wx.Panel(pane)
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.AddStretchSpacer(1)
        external = wx.Button(actions, label='External')
        close = wx.Button(actions, label='Close')
        external.Bind(wx.EVT_BUTTON, self.OnExternal)
        close.Bind(wx.EVT_BUTTON, lambda evt: self.Close())
        row.Add(external, 0, wx.ALL, 5)
        row.Add(close, 0, wx.ALL, 5)
        actions.SetSizer(row)
        actions.SetSizerProps(expand=True)
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def LoadFile(self, path):
        self.pdf_path = os.path.abspath(path)
        try:
            self.viewer.load_file(self.pdf_path)
            wx.CallAfter(self.viewer.fit_width)
        except Exception as exc:
            wx.MessageBox('Could not load PDF:\n%s' % exc,
                          'PDF viewer', wx.OK | wx.ICON_ERROR, parent=self)

    def OnFitWidth(self, event=None):
        self.viewer.fit_width()

    def OnFitPage(self, event=None):
        self.viewer.fit_page()

    def OnActualSize(self, event=None):
        self.viewer.actual_size()

    def OnZoomIn(self, event=None):
        self.viewer.zoom_in()

    def OnZoomOut(self, event=None):
        self.viewer.zoom_out()

    def OnClose(self, event=None):
        self.viewer.close()
        self.Destroy()

    def OnExternal(self, event=None):
        if not self.pdf_path or not os.path.isfile(self.pdf_path):
            wx.MessageBox('The PDF file is not available.', 'PDF viewer',
                          wx.OK | wx.ICON_ERROR, parent=self)
            return
        try:
            system = platform.system()
            if system == 'Windows':
                os.startfile(self.pdf_path)
            elif system == 'Darwin':
                subprocess.Popen(['open', self.pdf_path])
            else:
                subprocess.Popen(['xdg-open', self.pdf_path])
        except Exception as exc:
            wx.MessageBox('Could not open external PDF viewer:\n%s' % exc,
                          'PDF viewer', wx.OK | wx.ICON_ERROR, parent=self)


if __name__ == '__main__':
    app = wx.App(False)
    pdfV = PDFViewer(None, size=(800, 600))
    pdfV.Show()
    app.MainLoop()
