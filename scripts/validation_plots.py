"""Validation data and plotting helpers for the selected manuscript Fig. 1b."""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from figures_1kgp import ROOT, FROZEN, TARGETS, DEPTHS

from palette import INK, MUTED
BLUE="#287da7"; GOLD="#c08332"   # error_cloud only
GROUPS=["A","B","C","D","E","F","G","H","I","J","L","N","O","Q","R1a","R1b","R","T"]
KEYS=[row[0] for row in TARGETS]
from palette import KEY as _K, TITLE, NAMES, FS, MM
COLORS=[_K[k] for k in KEYS]   # one palette (palette.py)
LABELS={row[0]:row[2] for row in TARGETS}


def load():
    from depth_benchmark import load as load_depth
    depth, depth_files, _ = load_depth()
    cv = depth["1"]
    actual={(key,rep,r["sample"]) for key,groups in cv.items() for rep,rows in groups.items() for r in rows}
    assert len(actual)==sum(len(rows) for groups in cv.values() for rows in groups.values())
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
    assert len(seen)==sum(len({r["sample"] for r in cv[k][0]}) for k in cv)   # one lineage-out prediction per CV man and class
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
    ax.set_axisbelow(True)

def save(fig,out,stem):
    from palette import save as put
    put(fig,out/f"{stem}.pdf");plt.close(fig)

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
    ax.grid(True)
    if legend:ax.legend(loc="lower left",bbox_to_anchor=(0,1.02),ncol=4,frameon=False,
                        fontsize=6.8,title=NAMES['class'],title_fontproperties=TITLE(6.8),alignment='left')
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
    ax.grid(axis="x");style(ax)


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
                  flierprops=dict(marker='o',markersize=2.0,markeredgewidth=.25,
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
            box.set_alpha(None);box.set_facecolor(to_rgba(color,.18));box.set_edgecolor(color);box.set_linewidth(.6)
        for median in artists['medians']:
            median.set_linewidth(1);median.set_zorder(4)
        for whisker in artists['whiskers']:
            whisker.set_alpha(1);whisker.set_linewidth(.5)
        jitter=[r['y'] for r in records]
        ax.scatter(jitter if vertical else values,values if vertical else jitter,
                   s=2.5,color=color,linewidths=0,alpha=.55,zorder=3)
        maximum=max(maximum,max(values))
    if vertical:
        ax.set_yscale('symlog',linthresh=5,linscale=1)
        ax.set_ylim(0,max(130,maximum*1.05));ax.set_xlim(-.6,len(KEYS)-.4)
        ax.set_yticks([0,2,5,10,20,50,100],['0','2','5','10','20','50','100'])
        ax.set_xticks(range(len(KEYS)),[f"{LABELS[k].replace(':',':\n')}\n(n={counts(k)})" for k in KEYS])
        ax.set_ylabel('Mean absolute error (%)' if kind=='cv' else 'Absolute error (%)')
        ax.grid(axis='y');style(ax)
        return
    ax.set_xscale('symlog',linthresh=5,linscale=1)
    print(f'{kind} axis maximum {maximum:.2f}')
    ax.set_xlim(0,max(25,maximum*1.02));ax.set_ylim(len(KEYS)-.5,-.5)   # 25 covers every point (max 23.93)
    ax.set_xticks([0,2,5,10,25],['0','2','5','10','25'])
    ax.set_yticks(range(len(KEYS)))
    ax.set_xlabel('Mean absolute CV error (%)' if kind=='cv' else 'Absolute error (%)')
    ax.grid(axis='x');style(ax)


def annotate_vertical_outliers(ax,compact=False):
    """Place unframed top-three labels close to points, avoiding label collisions."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.transforms import Bbox
    ax.figure.canvas.draw()
    renderer=ax.figure.canvas.get_renderer()
    fontsize=FS['point']
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
                    arrowprops=dict(arrowstyle='-',color=color,linewidth=.4,shrinkA=1,shrinkB=1.5))
        ax.scatter([record['y']],[record['value']],s=6,facecolors='none',
                   edgecolors=INK,linewidths=.4,zorder=7)


def top_samples(key,kind):
    return sorted(individual_errors(key,kind),key=lambda r:(-r['value'],r['sample']))[:3]


def annotate_top_samples(ax,kind,compact=False):
    for i,key in enumerate(KEYS):
        for rank,r in enumerate(top_samples(key,kind),start=1):
            label_y=i+(rank-2)*.25
            ax.annotate(f"{rank} {r['sample']}",xy=(r['value'],r['y']),xycoords='data',
                        xytext=(1.045,label_y),textcoords=ax.get_yaxis_transform(),
                        ha='left',va='center',fontsize=FS['point'],color=INK,
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle='-',color=MUTED,linewidth=.45,
                                        shrinkA=1.5,shrinkB=2,
                                        connectionstyle='angle3,angleA=0,angleB=90'))
            ax.scatter([r['value']],[r['y']],s=6,
                       facecolors='none',edgecolors=INK,linewidths=.4,zorder=5)


def depth_individual_errors(key,fraction):
    groups=defaultdict(list)
    for rows in DEP[fraction][key].values():
        for r in rows:groups[r['sample']].append(abs(r['signed']))
    return [{'sample':sample,'value':float(np.mean(values)),'n_predictions':len(values)}
            for sample,values in sorted(groups.items())]


def depth_boxes(ax,compact=False):
    """One axis, one group per depth level (median autosomal coverage), one coloured box per target."""
    from matplotlib.colors import to_rgba
    fracs=sorted(DEPTHS,key=lambda f:DEPTHS[f])
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
            ax.scatter(pos+rng.uniform(-.035,.035,len(values)),values,s=.8,
                       color=color,alpha=.25,linewidths=0,zorder=3,rasterized=True)
        ax.plot(positions,medians,color='white',linewidth=1,zorder=5)
        ax.plot(positions,medians,'o-',color=color,linewidth=.75,markersize=2,
                markeredgecolor='white',markeredgewidth=.35,zorder=6)
    ax.set_xlim(-.55,len(fracs)-.45);ax.set_ylim(0,130);ax.set_yscale('symlog',linthresh=5,linscale=1)
    ax.set_yticks([0,2,5,10,20,50,100],['0','2','5','10','20','50','100'])
    from decimal import Decimal, ROUND_HALF_UP
    labels=[str(Decimal(str(DEPTHS[f])).quantize(Decimal('0.1'),rounding=ROUND_HALF_UP)) for f in fracs]
    ax.set_xticks(range(len(fracs)),labels)
    ax.set_xlabel('Median autosomal coverage (×)')
    ax.set_ylabel('Mean absolute error (%)')
    ax.set_title('')
    ax.yaxis.tick_left();ax.yaxis.set_label_position('left')
    ax.grid(axis='y');style(ax)
    ax.spines['left'].set_visible(True);ax.spines['right'].set_visible(False)


def lineage_supplement(out):
    fig,ax=plt.subplots(figsize=(180*MM,80*MM))
    fig.subplots_adjust(left=.1,right=.82,bottom=.12,top=.98)
    individual_boxes(ax,'lo');annotate_top_samples(ax,'lo')
    for i,key in enumerate(KEYS):   # CV reference: median of the per-male mean CV errors
        ax.plot(np.median([r['value'] for r in individual_errors(key,'cv')]),i,'s',color=INK,markersize=3,zorder=6)   # CV median; diamonds mean reference males
    ax.set_yticklabels(rowlabels())
    save(fig,out,'ed/ed6_lineage_validation')


