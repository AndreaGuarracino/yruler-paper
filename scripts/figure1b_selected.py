#!/usr/bin/env python3
"""Fig. 1b panels (draw), used by figures_1kgp.py: individual CV and lineage errors, depth-error distributions."""
import matplotlib.pyplot as plt
import validation_plots as data
from palette import TITLE, NAMES, FS, letter


def draw(fig, compact=False):
    from validation_ratio import view
    data.CV,data.LO,data.DEP,files=data.load()
    with view(data):
        return _draw(fig,compact,files)


def _draw(fig,compact,files):
    if compact:
        a=fig.add_axes([.06,.55,.40,.275]); b=fig.add_axes([.60,.55,.39,.275])
        letter(fig,'b',.005,.84)
    else:
        a=fig.add_axes([.085,.19,.37,.70]); b=fig.add_axes([.615,.19,.37,.70])
    data.individual_boxes(a,'cv',vertical=True)
    a.set_title('')
    data.depth_boxes(b,compact)
    # Retain every individual error, without reserving space for obsolete >100% tails.
    import math
    maximum=max(r['value'] for k in data.KEYS for f in data.DEPTHS
                for r in data.depth_individual_errors(k,f))
    maximum=max(maximum,max(r['value'] for k in data.KEYS for r in data.individual_errors(k,'cv')))
    limit=max(30,10*math.ceil(maximum*1.06/10))
    b.sharey(a)
    a.set_ylim(0,limit)
    ticks=[t for t in (0,2,5,10,20,30,40,50,100) if t<=limit]
    a.set_yticks(ticks,[str(t) for t in ticks])
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([],[],color=c,marker='o',markersize=2.5,linewidth=0.75,label=data.LABELS[k])
                        for k,c in zip(data.KEYS,data.COLORS)],
               loc='center',bbox_to_anchor=(.53,.68 if compact else .54),ncol=1,
               labelspacing=0.8,borderaxespad=0,title=NAMES['class'],title_fontproperties=TITLE(FS['legend']))
    data.annotate_vertical_outliers(a,compact)
    return a,b,[b],files
