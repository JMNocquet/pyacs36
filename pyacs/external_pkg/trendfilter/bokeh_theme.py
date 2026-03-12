"""
Set the Bokeh plotting defaults
"""
import os
from bokeh.io import curdoc
from bokeh.themes import Theme

_here = os.path.dirname(os.path.abspath(__file__))
_theme_path = os.path.join(_here, "bokeh_theme.yaml")
if os.path.isfile(_theme_path):
    curdoc().theme = Theme(filename=_theme_path)
