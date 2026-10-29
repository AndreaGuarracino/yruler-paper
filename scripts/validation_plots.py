"""Validation data and plotting helpers for the selected manuscript Fig. 1b."""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from frozen_snapshot import ROOT, FROZEN, TARGETS, DEPTHS

BLUE="#287da7"; GOLD="#c08332"; INK="#263f4b"
COLORS=["#287da7","#338565","#638590","#4c687d","#b09b4a","#b67442","#8f7258","#9962a9"]
GROUPS=["A","B","C","D","E","F","G","H","I","J","L","N","O","Q","R1a","R1b","R","T"]
KEYS=[row[0] for row in TARGETS]
LABELS={row[0]:row[2] for row in TARGETS}


def load():
    from depth_benchmark import load as load_depth
    depth, depth_files, _ = load_depth()
    cv = depth["1"]
    actual={(key,rep,r["sample"]) for key,groups in cv.items() for rep,rows in groups.items() for r in rows}
    assert len(actual)==7120
    lo=defaultdict(lambda:defaultdict(list))
    source=ROOT/"results/model/lineage_out.tsv"
    files=[*depth_files,source]
    seen=set()
    with source.open() as stream:
        for row in csv.DictReader(stream,delimiter="\t"):
            key,sample,group=row["class"],row["sample"],row["lineage"]
            assert (key,sample) not in seen, (key,sample)
            seen.add((key,sample))
            lo[key][group].append({"sample":sample,"signed":float(row["error_pct"]),"reference_bp":float(row["reference_bp"]),"estimate_bp":float(row["estimate_bp"])})
    assert len(seen)==712
    assert {g for groups in lo.values() for g in groups}==set(GROUPS)
    for key in KEYS:
        assert sorted(cv[key])==list(range(10))
        samples={r["sample"] for r in cv[key][0]}
        for rep,rows in cv[key].items():
            assert {r["sample"] for r in rows}==samples
            assert {r["fold"] for r in rows}==set(range(5))
        assert {r["sample"] for group in lo[key].values() for r in group}==samples
        for fraction in DEPTHS:
            assert set(depth[fraction][key]) == set(range(10))
            for rep,rows in depth[fraction][key].items():
                assert {r["sample"] for r in rows}==samples
                assert {r["fold"] for r in rows}==set(range(5))
                expected = {r["sample"]: r for r in cv[key][rep]}
                for row in rows:
                    assert row["fold"] == expected[row["sample"]]["fold"]
                    if fraction == "1":
                        assert np.isclose(row["signed"], expected[row["sample"]]["signed"], atol=.0001)
    return cv,lo,depth,files

def errors(rows):return np.abs([r["signed"] for r in rows])

def all_lo(lo,key):return [r for group in lo[key].values() for r in group]

def reps(key):return sorted(CV[key])

def depth_reps():return sorted(DEP["1"][KEYS[0]])

def counts(key):return len(CV[key][0])

def rowlabels():return [f"{LABELS[k]}  (n={counts(k)})" for k in KEYS]

def style(ax):
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.set_axisbelow(True)

def save(fig,out,stem):
    for ext in ("pdf","svg","png"):
        fig.savefig(out/f"{stem}.{ext}",dpi=190,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def depth_curves(ax,legend=True):
    fracs=sorted(DEPTHS,key=lambda f:DEPTHS[f])
    for key,color in zip(KEYS,COLORS):
        for r in depth_reps():
            vals=[np.median(errors(DEP[f][key][r])) for f in fracs]
            ax.plot([DEPTHS[f] for f in fracs],vals,color=color,alpha=.42,linewidth=.85)
        vals=np.array([[np.median(errors(DEP[f][key][r])) for r in depth_reps()] for f in fracs]).mean(axis=1)
        ax.plot([DEPTHS[f] for f in fracs],vals,"o",color=color,markersize=3.3,label=LABELS[key])
    ax.set_xscale("log");ax.set_xlim(.78,33);ax.set_ylim(0,8)
    ax.set_xticks([.9,2,5,10,30],["0.9","2","5","10","30"])
    ax.set_ylabel("Absolute error (%)",fontsize=8.5)
    ax.set_xlabel("Median autosomal depth (×)",fontsize=8.5)
    ax.set_title("Depth: every repeat, all targets",loc="left",fontsize=10.5,pad=10)
    ax.grid(color="#e8edf0",linewidth=.6)
    if legend:ax.legend(loc="lower left",bbox_to_anchor=(0,1.02),ncol=4,frameon=False,
                        fontsize=6.8,columnspacing=.7,handlelength=.7,handletextpad=.3)
    style(ax)

def error_cloud(ax,kind,target_colors=False):
    rng=np.random.default_rng(710)
    maximum=0
    for i,key in enumerate(KEYS):
        point_color=COLORS[i] if target_colors else (BLUE if kind=="cv" else GOLD)
        if kind=="cv":
            data=[(r,CV[key][r]) for r in reps(key)]
            for rep,rows in data:
                vals=errors(rows);maximum=max(maximum,float(max(vals)))
                y=i+(rep-4.5)*.037+rng.uniform(-.009,.009,len(vals))
                ax.scatter(vals,y,s=3.0,color=point_color,alpha=1,linewidths=0,rasterized=True)
            medians=[np.median(errors(CV[key][r])) for r in reps(key)]
            ax.plot(medians,np.full(10,i),"|",color=INK,markersize=6,markeredgewidth=.7)
        else:
            for j,g in enumerate(GROUPS):
                if g not in LO[key]:continue
                vals=errors(LO[key][g]);maximum=max(maximum,float(max(vals)))
                y=i+(j-8.5)*.02+rng.uniform(-.012,.012,len(vals))
                ax.scatter(vals,y,s=9,color=point_color,alpha=1,linewidths=.15,edgecolors=INK,rasterized=True)
            ax.plot(np.median(errors(all_lo(LO,key))),i,"D",color=INK,markersize=3.3)
    ax.set_xscale("symlog",linthresh=5,linscale=1)
    ax.set_xlim(0,max(120,maximum*1.05))
    ticks=[0,2,5,10,20,50,100]
    ax.set_xticks(ticks,[str(v) for v in ticks])
    ax.set_ylim(len(KEYS)-.5,-.5);ax.set_yticks(range(len(KEYS)))
    ax.set_xlabel("Individual absolute error (%)",fontsize=8.5)
    ax.grid(axis="x",color="#e8edf0",linewidth=.6);style(ax)


def error_box_stats(rows):
    """Exact empirical percentiles; tail values remain visible beyond the whiskers."""
    values=errors(rows)
    q10,q25,median,q75,q90=np.percentile(values,[10,25,50,75,90])
    return dict(q1=q25,med=median,q3=q75,whislo=q10,whishi=q90,
                fliers=values[(values<q10)|(values>q90)])


def draw_error_boxes(ax,stats,positions,color,width,vertical=False,showfliers=True):
    return ax.bxp(stats,positions=positions,widths=width,vert=vertical,
                  patch_artist=True,manage_ticks=False,showcaps=False,showfliers=showfliers,
                  boxprops=dict(facecolor=color,edgecolor=color,alpha=.48,linewidth=.45),
                  whiskerprops=dict(color=color,alpha=.65,linewidth=.55),
                  medianprops=dict(color=INK,linewidth=.7),
                  flierprops=dict(marker='o',markersize=2.0,markeredgewidth=.15,
                                  markeredgecolor=INK,markerfacecolor=color,alpha=1))


def individual_errors(key,kind):
    """One observation per individual: mean absolute CV error or lineage-out error."""
    by_sample=defaultdict(list)
    source=CV[key] if kind=='cv' else LO[key]
    for rows in source.values():
        for row in rows:by_sample[row['sample']].append(row['signed'])
    rng=np.random.default_rng(710+KEYS.index(key))
    records=[]
    for sample,values in sorted(by_sample.items()):
        absolute=np.abs(values)
        records.append(dict(sample=sample,value=float(np.mean(absolute)),
                            mean_signed_error_pct=float(np.mean(values)),n_predictions=len(values),
                            min_absolute_error_pct=float(min(absolute)),max_absolute_error_pct=float(max(absolute)),
                            y=KEYS.index(key)+rng.uniform(-.11,.11)))
    return records


def individual_boxes(ax,kind,vertical=False):
    maximum=0
    for i,(key,color) in enumerate(zip(KEYS,COLORS)):
        records=individual_errors(key,kind)
        values=[r['value'] for r in records]
        stats=error_box_stats([{'signed':v} for v in values])
        artists=draw_error_boxes(ax,[stats],[i],color,.56,vertical=vertical,showfliers=False)
        from matplotlib.colors import to_rgba
        for box in artists['boxes']:
            box.set_alpha(None);box.set_facecolor(to_rgba(color,.18));box.set_edgecolor(color);box.set_linewidth(1.05)
        for median in artists['medians']:
            median.set_linewidth(1.3);median.set_zorder(4)
        for whisker in artists['whiskers']:
            whisker.set_alpha(1);whisker.set_linewidth(.85)
        jitter=[r['y'] for r in records]
        ax.scatter(jitter if vertical else values,values if vertical else jitter,
                   s=5,color=color,edgecolors=INK,linewidths=.1,alpha=.55,zorder=3)
        maximum=max(maximum,max(values))
    if vertical:
        ax.set_yscale('symlog',linthresh=5,linscale=1)
        ax.set_ylim(0,max(130,maximum*1.05));ax.set_xlim(-.6,len(KEYS)-.4)
        ax.set_yticks([0,2,5,10,20,50,100],['0','2','5','10','20','50','100'])
        ax.set_xticks(range(len(KEYS)),[f'{LABELS[k]}\n(n={counts(k)})' for k in KEYS])
        ax.set_ylabel('Mean absolute error (%)' if kind=='cv' else 'Absolute error (%)',fontsize=8)
        ax.grid(axis='y',color='#e8edf0',linewidth=.6);style(ax)
        return
    ax.set_xscale('symlog',linthresh=5,linscale=1)
    ax.set_xlim(0,max(30,maximum*1.05));ax.set_ylim(len(KEYS)-.5,-.5)
    ax.set_xticks([0,2,5,10,20,30],['0','2','5','10','20','30'])
    ax.set_yticks(range(len(KEYS)))
    ax.set_xlabel('Mean absolute CV error (%)' if kind=='cv' else 'Absolute error (%)',fontsize=8.5)
    ax.grid(axis='x',color='#e8edf0',linewidth=.6);style(ax)


def annotate_vertical_outliers(ax,compact=False):
    """Place unframed top-three labels close to points, avoiding label collisions."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.transforms import Bbox
    ax.figure.canvas.draw()
    renderer=ax.figure.canvas.get_renderer()
    fontsize=5.3 if compact else 5.7
    font=FontProperties(size=fontsize)
    bounds=ax.get_window_extent(renderer)
    column_width=bounds.width/(len(KEYS)+.2)
    candidates=[(i,key,color,r) for i,(key,color) in enumerate(zip(KEYS,COLORS))
                for r in top_samples(key,'cv')]
    candidates.sort(key=lambda item:(-item[3]['value'],item[0],item[3]['sample']))
    anchors=[ax.transData.transform((r['y'],r['value'])) for _,_,_,r in candidates]
    occupied=[]
    for i,key,color,record in candidates:
        x,y=ax.transData.transform((record['y'],record['value']))
        width,height,descent=renderer.get_text_width_height_descent(record['sample'],font,False)
        width+=5; height=max(height+descent+4,fontsize*ax.figure.dpi/72+3)
        options=[]
        for level in range(9):
            for direction in (1,-1):
                for shift in (0,-.25,.25,-.45,.45):
                    cx=np.clip(x+shift*column_width,bounds.x0+width/2+2,bounds.x1-width/2-2)
                    cy=np.clip(y+direction*(height*.75+level*height),bounds.y0+height/2+2,bounds.y1-height/2-2)
                    box=Bbox.from_bounds(cx-width/2,cy-height/2,width,height)
                    collisions=sum(box.overlaps(other) for other in occupied)
                    covered=sum(box.contains(px,py) for px,py in anchors)
                    score=collisions*1e6+covered*1e4+(cx-x)**2+(cy-y)**2+(3 if direction<0 else 0)
                    options.append((score,cx,cy,box))
        _,cx,cy,box=min(options,key=lambda item:item[0])
        occupied.append(box)
        position=ax.transAxes.inverted().transform((cx,cy))
        ax.annotate(record['sample'],xy=(record['y'],record['value']),
                    xytext=position,textcoords=ax.transAxes,
                    ha='center',va='center',fontsize=fontsize,color=INK,zorder=8,
                    bbox=dict(boxstyle='square,pad=0.1',facecolor='white',edgecolor='none'),
                    arrowprops=dict(arrowstyle='-',color=color,linewidth=.55,shrinkA=2,shrinkB=2))
        ax.scatter([record['y']],[record['value']],s=12,facecolors='none',
                   edgecolors=INK,linewidths=.65,zorder=7)


def fold_supplement(out):
    """Two measurements per page; retain all 50 test-fold distributions."""
    import csv
    from matplotlib.backends.backend_pdf import PdfPages
    summaries=[]
    with PdfPages(out/'figureS6_cv_folds.pdf') as multipage:
        for page in range((len(KEYS)+1)//2):
            fig,axes=plt.subplots(1,2,figsize=(10.8,12.6),sharex=True,sharey=True)
            if len(KEYS[2*page:2*page+2])==1:axes[1].set_visible(False)
            for ax,key,color,letter in zip(axes,KEYS[2*page:2*page+2],COLORS[2*page:2*page+2],'abcdefghijklmnopqrstuvwxyz'[2*page:2*page+2]):
                stats=[];positions=[]
                for rep in reps(key):
                    rows=CV[key][rep]
                    for fold in [None,*range(5)]:
                        selected=rows if fold is None else [r for r in rows if r['fold']==fold]
                        values=errors(selected);q=np.percentile(values,[0,10,25,50,75,90,100])
                        summaries.append(dict(target=key,repeat=rep+1,fold='all' if fold is None else fold+1,n=len(values),
                                              **dict(zip(['min','p10','p25','median','p75','p90','max'],q))))
                        if fold is not None:
                            stats.append(error_box_stats(selected));positions.append(rep*6+fold)
                assert len(stats)==50
                draw_error_boxes(ax,stats,positions,color,.62)
                ax.set_xscale('symlog',linthresh=5,linscale=1);ax.set_xlim(0,30);ax.set_ylim(59,-1)
                ax.set_xticks([0,2,5,10,20,30],['0','2','5','10','20','30'])
                ax.set_yticks([r*6+2 for r in range(10)],[f'Repeat {r+1}' for r in range(10)])
                ax.tick_params(labelleft=True)
                for rep in range(0,10,2):ax.axhspan(rep*6-.5,rep*6+4.5,color='#f2f5f7',zorder=-2)
                ax.grid(axis='x',color='#e8edf0',linewidth=.6);style(ax)
                ax.set_title(f'{letter}   {LABELS[key]}',loc='left',fontsize=12,fontweight='bold',color=color)
                ax.set_xlabel('Individual absolute error (%)',fontsize=9)
            fig.subplots_adjust(left=.10,right=.985,top=.94,bottom=.085,wspace=.30)
            fig.suptitle('Prediction-error distributions in every test fold',fontsize=14,y=.98,color=INK)
            fig.text(.10,.035,'Within each repeat, folds 1–5 run top to bottom. Boxes: 25th–75th percentiles; median; whiskers: 10th–90th.',fontsize=8,color=INK)
            fig.text(.10,.015,'Points retain errors beyond the whiskers. Axes are logarithmic above 5%. Repeated individuals are not independent.',fontsize=8,color=INK)
            multipage.savefig(fig,bbox_inches='tight',facecolor='white')
            for ext in ['pdf','png','svg']:
                fig.savefig(out/f'figureS6_cv_folds_page{page+1}.{ext}',dpi=190,bbox_inches='tight',facecolor='white')
            if page==0:fig.savefig(out/'figureS6_cv_folds.png',dpi=190,bbox_inches='tight',facecolor='white')
            plt.close(fig)
    old_svg=out/'figureS6_cv_folds.svg'
    if old_svg.exists():old_svg.unlink()
    with (out/'figure1b_cv_distribution_summary.tsv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summaries[0]),delimiter='\t')
        writer.writeheader();writer.writerows(summaries)


def top_samples(key,kind):
    return sorted(individual_errors(key,kind),key=lambda r:(-r['value'],r['sample']))[:3]


def annotate_top_samples(ax,kind,compact=False):
    for i,key in enumerate(KEYS):
        for rank,r in enumerate(top_samples(key,kind),start=1):
            label_y=i+(rank-2)*.25
            ax.annotate(f"{rank} {r['sample']}",xy=(r['value'],r['y']),xycoords='data',
                        xytext=(1.045,label_y),textcoords=ax.get_yaxis_transform(),
                        ha='left',va='center',fontsize=5.8 if compact else 7.5,color=INK,
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle='-',color='#52646e',linewidth=.45,
                                        shrinkA=1.5,shrinkB=2,
                                        connectionstyle='angle3,angleA=0,angleB=90'))
            ax.scatter([r['value']],[r['y']],s=11 if compact else 16,
                       facecolors='none',edgecolors=INK,linewidths=.65,zorder=5)


def write_top_samples(out):
    import csv
    rows=[];individual_rows=[]
    for kind in ['cv','lineage-out']:
        for key in KEYS:
            records=individual_errors(key,'cv' if kind=='cv' else 'lo')
            for r in records:
                individual_rows.append(dict(evaluation=kind,target=LABELS[key],sample=r['sample'],
                                            plotted_absolute_error_pct=r['value'],
                                            **{k:r[k] for k in ['n_predictions','mean_signed_error_pct','min_absolute_error_pct','max_absolute_error_pct']}))
            for rank,r in enumerate(top_samples(key,'cv' if kind=='cv' else 'lo'),start=1):
                rows.append(dict(evaluation=kind,target=LABELS[key],rank=rank,sample=r['sample'],
                                 plotted_absolute_error_pct=r['value'],mean_signed_error_pct=r['mean_signed_error_pct'],
                                 n_predictions=r['n_predictions']))
    for filename,records in [('figure1b_top3_samples.tsv',rows),('figure1b_individual_errors.tsv',individual_rows)]:
        with (out/filename).open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(records[0]),delimiter='\t')
            writer.writeheader();writer.writerows(records)
    return rows


def depth_individual_errors(key,fraction):
    groups=defaultdict(list)
    for rows in DEP[fraction][key].values():
        for r in rows:groups[r['sample']].append(abs(r['signed']))
    return [{'sample':sample,'value':float(np.mean(values)),'n_predictions':len(values)}
            for sample,values in sorted(groups.items())]


def depth_boxes(ax,compact=False):
    """One axis, seven median-autosomal-coverage groups, eight coloured target boxes per group."""
    from matplotlib.colors import to_rgba
    fracs=sorted(DEPTHS,key=float)
    for i,(key,color) in enumerate(zip(KEYS,COLORS)):
        positions=[];medians=[]
        for j,fraction in enumerate(fracs):
            records=depth_individual_errors(key,fraction)
            values=[r['value'] for r in records]
            spacing=.84/len(KEYS)
            pos=j+(i-(len(KEYS)-1)/2)*spacing
            positions.append(pos);medians.append(float(np.median(values)))
            artists=draw_error_boxes(ax,[error_box_stats([{'signed':v} for v in values])],[pos],color,spacing*.88,vertical=True,showfliers=False)
            for box in artists['boxes']:
                box.set_alpha(None);box.set_facecolor(to_rgba(color,.10));box.set_edgecolor(to_rgba(color,.45));box.set_linewidth(.6)
            for whisker in artists['whiskers']:whisker.set_alpha(.35)
            for median in artists['medians']:median.set_zorder(4);median.set_linewidth(.85)
            rng=np.random.default_rng(915+i*10+j)
            ax.scatter(pos+rng.uniform(-.035,.035,len(values)),values,s=1.6 if compact else 2,
                       color=color,alpha=.25,linewidths=0,zorder=3,rasterized=True)
        ax.plot(positions,medians,color='white',linewidth=2.5,zorder=5)
        ax.plot(positions,medians,'o-',color=color,linewidth=1.35,markersize=2.8,
                markeredgecolor='white',markeredgewidth=.35,zorder=6)
    ax.set_xlim(-.55,6.55);ax.set_ylim(0,130);ax.set_yscale('symlog',linthresh=5,linscale=1)
    ax.set_yticks([0,2,5,10,20,50,100],['0','2','5','10','20','50','100'])
    from decimal import Decimal, ROUND_HALF_UP
    labels=[str(Decimal(str(DEPTHS[f])).quantize(Decimal('0.1'),rounding=ROUND_HALF_UP)) for f in fracs]
    ax.set_xticks(range(7),labels)
    ax.set_xlabel('Median autosomal coverage (×)',fontsize=8)
    ax.set_ylabel('Mean absolute error (%)',fontsize=8)
    ax.set_title('')
    ax.yaxis.tick_left();ax.yaxis.set_label_position('left')
    ax.grid(axis='y',color='#e8edf0',linewidth=.5);style(ax)
    ax.spines['left'].set_visible(True);ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=6.8 if compact else 8)


def lineage_supplement(out):
    fig,axes=plt.subplots(1,2,figsize=(10.4,5.7))
    fig.subplots_adjust(left=.15,right=.85,bottom=.13,top=.83,wspace=.62)
    for ax,kind,title in zip(axes,['cv','lo'],['Random five-fold CV, ten repeats','Leave one lineage out, 18 lineages']):
        individual_boxes(ax,kind);annotate_top_samples(ax,kind)
        ax.set_title(title,loc='left',fontsize=10,pad=12)
        ax.set_yticklabels(rowlabels() if kind=='cv' else [],fontsize=9)
    fig.text(.15,.94,'One point per individual: mean of ten CV errors (left); one excluded-lineage error (right).',fontsize=9,color=INK)
    fig.text(.15,.035,'Boxes: 25th–75th percentiles; median; whiskers: 10th–90th. Labels identify the three highest errors per target.',fontsize=8,color=INK)
    save(fig,out,'figureS7_lineage_validation')


def write_depth_individuals(out):
    import csv
    rows=[]
    for key in KEYS:
        for fraction in sorted(DEPTHS,key=lambda f:DEPTHS[f]):
            for r in depth_individual_errors(key,fraction):
                rows.append(dict(target=LABELS[key],fraction=fraction,median_autosomal_depth=DEPTHS[fraction],
                                 sample=r['sample'],mean_absolute_error_pct=r['value'],cv_repeats=r['n_predictions']))
    with (out/'figure1b_depth_individual_errors.tsv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
