import numpy
from bokeh.plotting import figure, ColumnDataSource

def rafa_plot(matrix, rep1, rep2, labels=None, spike_prefix='gSpike'):
    if labels is not None:
        xname = labels[rep1]
        yname = labels[rep2]
    else:
        xname = rep1
        yname = rep2

    eitherzero = (matrix[rep1] == 0) | (matrix[rep2] == 0)
    replz1 = numpy.log2(matrix[rep1][eitherzero != True])
    replz2 = numpy.log2(matrix[rep2][eitherzero != True])

    M = replz1 - replz2
    A = (replz1 + replz2) / 2.0
    Acutoff = 0

    spikes = numpy.asarray([ x.startswith(spike_prefix) for x in replz1.index ])
    filtered = A < Acutoff
    normal = (spikes == False) & (filtered == False)

    f = figure(title="{} v {}".format(xname, yname))
    f.cross(replz1[normal], replz2[normal], color='blue', fill_alpha=0.2, line_alpha=0.4)
    f.circle(replz1[spikes], replz2[spikes], color='black', fill_alpha=0.8, line_alpha=0.4)
    f.cross(replz1[filtered], replz2[filtered], color='red', fill_alpha=0.2, line_alpha=0.4)
    return f
