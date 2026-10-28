#!/usr/bin/env python3
"""Fig. 1b: individual-level CV/lineage summaries and depth-error distributions."""
import json
import hashlib
import matplotlib.pyplot as plt
import validation_plots as data


def draw(fig, compact=False):
    from validation_ratio import view
    data.CV,data.LO,data.DEP,files=data.load()
    with view(data):
        return _draw(fig,compact,files)


def _draw(fig,compact,files):
    if compact:
        a=fig.add_axes([.085,.53,.37,.285]); b=fig.add_axes([.615,.53,.37,.285])
        fig.text(.055,.82,'b',fontsize=15,fontweight='bold')
    else:
        a=fig.add_axes([.085,.19,.37,.70]); b=fig.add_axes([.615,.19,.37,.70])
    data.individual_boxes(a,'cv',vertical=True)
    a.set_title('')
    a.tick_params(labelsize=7 if compact else 8)
    a.tick_params(axis='x',labelsize=5.7 if compact else 6.1)
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
    fig.legend(handles=[Line2D([],[],color=c,marker='o',markersize=3,linewidth=1.2,label=data.LABELS[k])
                        for k,c in zip(data.KEYS,data.COLORS)],
               loc='center',bbox_to_anchor=(.505,.672 if compact else .54),ncol=1,
               frameon=False,fontsize=6.6 if compact else 7,handlelength=1,
               handletextpad=.4,labelspacing=1.0,borderaxespad=0)
    data.annotate_vertical_outliers(a,compact)
    return a,b,[b],files


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    out=data.ROOT/'manuscript/figures'
    fig=plt.figure(figsize=(10.4,3.3));a,b,axes,files=draw(fig)
    fig.text(.085,.025,'CV and depth: ten repeats. Boxes: quartiles and median; whiskers: 10th–90th. Both error axes logarithmic above 5%.',fontsize=7.5,color=data.INK)
    data.save(fig,out,'figure1b_selected')
    from validation_ratio import view
    with view(data):
        data.fold_supplement(out)
        data.lineage_supplement(out)
        top=data.write_top_samples(out);data.write_depth_individuals(out)
        manifest={'design':'One point per individual and one box per class/condition; CV means over ten absolute errors, lineage single errors, depth means over ten repeats at each of seven fractions. Main: vertical CV boxes by class and grouped vertical coverage boxes with matching left error axes versus median autosomal coverage (×), with lines connecting the boxplot medians; central shared legend and unframed top-three CV labels placed near their points; lineage comparison with top-three labels in Supplementary Fig. S7.',
                  'colors':dict(zip(data.KEYS,data.COLORS)),
                  'cv_predictions':sum(len(rows) for x in data.CV.values() for rows in x.values()),
                  'lineage_out_predictions':sum(len(rows) for x in data.LO.values() for rows in x.values()),
                  'depth_predictions':sum(len(rows) for x in data.DEP.values() for y in x.values() for rows in y.values()),
                  'base_cv_predictions':7000,'derived_ratio_cv_predictions':820,'base_depth_predictions':49000,'derived_ratio_depth_predictions':5740,'base_lineage_out_predictions':700,'derived_ratio_lineage_out_predictions':82,'cv_repeats':10,'cv_boxplots':9,'lineage_out_boxplots_main':0,'lineage_out_boxplots_supplement':9,'depth_boxplots':63,'fold_supplement_boxplots':450,
                  'box_percentiles':[25,50,75],'whisker_percentiles':[10,90],'tails':'all individuals plotted',
                  'annotated_samples':top,'lineages':data.GROUPS,'depth_repeats':data.depth_reps(),
                  'source_files':[{'path':str(p.relative_to(data.ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]}
    (out/'figure1b_selected_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':main()
