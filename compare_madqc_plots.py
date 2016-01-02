from __future__ import print_function, division
import pandas
import numpy
from IPython.display import display
import bokeh
from bokeh import mpl
from bokeh.plotting import figure, show, ColumnDataSource, output_file
from bokeh.io import output_notebook
import bokeh.io
import bokeh.resources
from bokeh.models import HoverTool
from bokeh.palettes import *

def create_figure(xname, yname, extra_title='', **kwargs):
    hover = HoverTool(
        tooltips = [
            (xname, '@'+xname),
            (yname, '@'+yname),
            ('accession', '@experiment'),
            # description is too long for my tooltip box
            #('description', '@description'),
            ('organism', '@organism'),
            ('biosample', '@biosample'),
            ('source', '@biosample_lab'),
            ('starting', '@starting'),
            ('age', '@age'),
            ('lab', '@lab'),
            ('rfa', '@rfa'),
        ]
    )    

    p = figure(
        title = "{} vs {} {}".format(xname, yname, extra_title),
        tools=['box_zoom', 'wheel_zoom', 'pan', hover, 'save', 'reset']
    )
    p.xaxis.axis_label = xname
    p.yaxis.axis_label = yname
    return p

def setdefault_style(**kwargs):
    extra = kwargs.copy()
    extra.setdefault('fill_alpha', 0.4)
    extra.setdefault('size', 7)
    extra.setdefault('line_color', 'black')
    extra.setdefault('line_alpha', 0.4)
    return extra
   

def experiment_scatter(matrix, xname, yname, extra_title='', **kwargs):
    """Make a scatter plot comparing two numeric columns from our experiment matrix
    
    Provides a hover tool tip tool.
    """
    p = create_figure(xname, yname, extra_title)

    extra = setdefault_style(**kwargs)
    
    p.circle(xname, 
             yname, 
             source=ColumnDataSource(matrix),
             **extra)
    #show(p)
    return p

