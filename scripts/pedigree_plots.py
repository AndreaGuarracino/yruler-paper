"""All-target father–son comparisons, retaining batch and reference status."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import figures_1kgp as f
from validation_plots import COLORS,INK
from palette import FS, MM, save, letter, REF, REF_HANDLE, BATCH, MUTED, head, top_legend
from validation_ratio import KEY as RATIO, LABEL as RATIO_LABEL, COLOR as RATIO_COLOR

MEASUREMENTS=[(key,col,label,scale,unit,color) for (key,col,label,scale,unit),color in zip(f.TARGETS,COLORS)]
MEASUREMENTS.append((RATIO,'ratio_bp',RATIO_LABEL,1,'',RATIO_COLOR))
LOOKUP={row[0]:row for row in MEASUREMENTS}
GROUPS=['Cross batch','Same batch','Reverse cross batch']
MARKERS=[BATCH[g] for g in GROUPS]


def collect(valid):
    index={r['sample']:r for r in valid};rows=[];seen=set()
    for father,son in csv.reader((f.ROOT/'data/catalog_d376/fs_pairs.tsv').open(),delimiter='\t'):
        if father not in index or son not in index:continue
        assert (father,son) not in seen;seen.add((father,son))
        a,b=index[father],index[son]
        batches=(a['batch'].removesuffix('_related'),b['batch'].removesuffix('_related'))
        group=('Same batch' if batches[0]==batches[1] else
               'Cross batch' if batches==('2504','698') else 'Reverse cross batch')
        for key,column,label,scale,unit,color in MEASUREMENTS:
            x,y=a[column],b[column]
            assert x>0 and y>0 and np.isfinite(x+y)
            rows.append(dict(father=father,son=son,family_father=a['family'],family_son=b['family'],
                batch_group=group,target=key,label=label,unit=unit or 'ratio',father_value=x/scale,son_value=y/scale,
                father_reference=a['reference_male'],son_reference=b['reference_male'],
                father_y_status=a['y_status'],son_y_status=b['y_status'],
                father_population_included=a['population_included'],son_population_included=b['population_included'],
                signed_difference_pct=100*(y/x-1)))
    return rows


def stats(rows):
    output=[]
    for key,*_ in MEASUREMENTS:
        for group in ['All pairs',*GROUPS]:
            s=[r for r in rows if r['target']==key and (group=='All pairs' or r['batch_group']==group)]
            delta=np.array([r['signed_difference_pct'] for r in s]);q=np.percentile(delta,[25,50,75])
            output.append(dict(target=key,batch_group=group,n_pairs=len(s),n_fathers=len({r['father'] for r in s}),
                n_reference_pairs=sum(r['father_reference']=='yes' or r['son_reference']=='yes' for r in s),
                median_signed_pct=float(q[1]),median_absolute_pct=float(np.median(abs(delta))),
                signed_q25_pct=float(q[0]),signed_q75_pct=float(q[2]),
                maximum_absolute_pct=float(max(abs(delta)))))
    return output


def style(ax):
    ax.spines[['right','top']].set_visible(False)
    ax.grid(True);ax.set_axisbelow(True)


def scatter(ax,rows,key,compact=False):
    _,_,label,scale,unit,color=LOOKUP[key];s=[r for r in rows if r['target']==key]
    upper=1.08*max(max(r['father_value'],r['son_value']) for r in s)
    ax.plot([0,upper],[0,upper],'--',color=MUTED,linewidth=.65,zorder=1)
    for group,marker in zip(GROUPS,MARKERS):
        for ref in (False,True):
            subset=[r for r in s if r['batch_group']==group and
                    (r['father_reference']=='yes' or r['son_reference']=='yes')==ref]
            ax.scatter([r['father_value'] for r in subset],[r['son_value'] for r in subset],zorder=3,
                **(dict(edgecolors=color,**REF) if ref else dict(marker=marker,s=(4 if group=='Cross batch' else 9) if compact else (6 if group=='Cross batch' else 12),
                color=color,linewidths=0,alpha=.5)))
    worst=sorted(s,key=lambda r:-abs(r['signed_difference_pct']))[:1 if compact else 2]
    if compact:
        for r in s:
            if r['son']=='HG02966' and abs(r['signed_difference_pct'])>=20 and r not in worst:worst.append(r)
    n={True:0,False:0}   # labels per empty corner
    for r in worst:
        if abs(r['signed_difference_pct'])<10:continue
        x,y=r['father_value'],r['son_value'];up=y>x;k=n[up];n[up]+=1   # above the diagonal: upper left, else lower right
        ax.annotate(r['son'],(x,y),xytext=(.04,.8-k*.1) if up else (.96,.05+k*.1),textcoords='axes fraction',
                    fontsize=FS['point'],ha='left' if up else 'right',va='top' if up else 'bottom',color=INK,bbox=dict(facecolor='none',edgecolor='none',pad=.3),   # line starts at the text edge
                    arrowprops=dict(arrowstyle='-',color=color,linewidth=.4,shrinkA=1,shrinkB=2))
    med=np.median([abs(r['signed_difference_pct']) for r in s])
    if compact: ax.text(0,1.03,f'Median |Δ| {med:.2f}%',transform=ax.transAxes,va='bottom',ha='left',fontsize=FS['point'],color=INK)   # above the axes, off the data
    else: ax.text(.035,.965,f'Median |Δ| {med:.2f}%',transform=ax.transAxes,va='top',fontsize=FS['point'],color=INK)
    ax.set_xlim(0,upper);ax.set_ylim(0,upper);ax.set_aspect('equal',adjustable='box')
    ax.xaxis.set_major_locator(MaxNLocator(3 if compact else 4));ax.yaxis.set_major_locator(MaxNLocator(3 if compact else 4))
    suffix=f' ({unit})' if unit else ''
    ax.set_xlabel('Parent'+suffix)
    ax.set_ylabel('Offspring'+suffix,labelpad=1)
    style(ax)
    return med


def differences(ax,rows,key,compact=False,reverse=True):
    color=LOOKUP[key][-1];groups=GROUPS if reverse else GROUPS[:2]
    selected=[r for r in rows if r['target']==key and r['batch_group'] in groups]
    maximum=max(10,1.10*max(abs(r['signed_difference_pct']) for r in selected))
    rng=np.random.default_rng(20261028+list(LOOKUP).index(key))
    for i,group in enumerate(groups):
        s=[r for r in selected if r['batch_group']==group]
        values=np.array([r['signed_difference_pct'] for r in s]);ys=i+rng.uniform(-.17,.17,len(s))
        for ref in (False,True):
            mask=np.array([(r['father_reference']=='yes' or r['son_reference']=='yes')==ref for r in s])
            ax.scatter(values[mask],ys[mask],zorder=3,**(dict(edgecolors=color,**REF) if ref else dict(s=4 if compact else 6,marker=MARKERS[i],color=color,linewidths=0,alpha=.35)))
        lo,median,hi=np.percentile(values,[25,50,75])
        ax.barh(i,hi-lo,left=lo,height=.14,color=color,linewidth=0,zorder=4)
        ax.plot(median,i,'|',color=INK,markersize=6,markeredgewidth=1,zorder=5)
        ax.annotate(f'{median:+.2f}%',(1,i),xycoords=('axes fraction','data'),xytext=(3,0),textcoords='offset points',fontsize=FS['point'],color=color,fontweight='bold',va='center',ha='left')   # median right of the panel, off the data
    ax.axvline(0,color=MUTED,linestyle='--',linewidth=.65)
    ax.set_xscale('symlog',linthresh=5);ax.set_xlim(-maximum,maximum)
    tick_pool=([-100,0,100] if maximum>=100 else [-20,0,20] if maximum>=20 else [-10,0,10]) if compact else [-100,-20,-5,0,5,20,100]
    ticks=[x for x in tick_pool if abs(x)<=maximum]
    ax.set_xticks(ticks,[str(x) for x in ticks]);ax.set_ylim(len(groups)-.55,-.55)
    short=['Cross\nbatch','Same\nbatch','Reverse\ncross batch'] if compact else GROUPS
    ax.set_yticks(range(len(groups)),short[:len(groups)])
    ax.tick_params(axis='y',length=0);ax.set_xlabel('Offspring/parent difference (%)')
    style(ax);ax.spines['left'].set_visible(False)


def legend(fig,rows):
    first=[r for r in rows if r['target']=='dyz3'];counts=Counter(r['batch_group'] for r in first)
    h=[Line2D([],[],marker=m,linestyle='none',color=MUTED,markersize=3,
       label=f'{g} (n={counts[g]})') for g,m in zip(GROUPS,MARKERS)]
    h.append(Line2D([],[],markeredgecolor=MUTED,label='Pair with a truth-set male',**REF_HANDLE))
    top_legend(fig,handles=h,ncol=4)   # same gap as ED7 and ED8


def supplement(rows,out):
    fig=plt.figure(figsize=(180*MM,170*MM))   # one page: scatter and differences per measurement
    g=fig.add_gridspec(4,4,left=.065,right=.93,bottom=.04,top=.92,width_ratios=[1,1.2,1,1.2],wspace=.62,hspace=.42)
    pairs,diffs=[],[]
    for i,(key,_,label,_,_,color) in enumerate(MEASUREMENTS):
        r,c=divmod(i,2);a,b=fig.add_subplot(g[r,2*c]),fig.add_subplot(g[r,2*c+1])
        scatter(a,rows,key,compact=True);differences(b,rows,key,compact=True);pairs.append((a,c));diffs.append((a,b,label,color))
    fig.canvas.draw()   # equal-aspect axes settle
    th=[]   # class title centred over each scatter + difference pair
    for a,b,label,color in diffs:   # difference panel at the height of its scatter: one x axis level per row
        pa,pb=a.get_position(),b.get_position();b.set_position([pb.x0,pa.y0,pb.width,pa.height])
        y=pa.y1+4.2*MM/fig.get_figheight();th.append(y);fig.text((pa.x0+pb.x1)/2,y,label,ha='center',va='baseline',fontsize=FS['title'],fontweight='bold',color=color)
    for i,(a,c) in enumerate(pairs):fig.text([.002,.5][c],th[i],chr(97+i),fontsize=FS['letter'],fontweight='bold',va='baseline')   # letter on the title baseline
    legend(fig,rows)
    save(fig,out/'ed/ed6_pedigrees.pdf');plt.close(fig)


def main():
    _,valid,_=f.load_snapshot();rows=collect(valid);out=f.ROOT/'manuscript/figures'
    assert len(rows)%len(MEASUREMENTS)==0; npairs=len(rows)//len(MEASUREMENTS); print('father-son pairs',npairs,flush=True)
    summary=stats(rows);f.write_rows(f.ROOT/'results/figures/pedigree_pairs.tsv',rows);f.write_rows(f.ROOT/'results/figures/pedigree_summary.tsv',summary)
    supplement(rows,out)
    inputs={'cohort':str(f.DEFAULT_SNAPSHOT.relative_to(f.ROOT)),'cohort_sha256':hashlib.sha256(f.DEFAULT_SNAPSHOT.read_bytes()).hexdigest(),
           'pairs_per_measurement':npairs,'measurement_records':len(rows),'measurements':[m[2] for m in MEASUREMENTS],
           'batch_counts':dict(Counter(r['batch_group'] for r in rows if r['target']=='dyz3')),
           'reference_involved_pairs':sum(r['father_reference']=='yes' or r['son_reference']=='yes' for r in rows if r['target']=='dyz3'),
          'profile_policy':'All complete profiles retained, including flagged arm imbalance; no dosage cutoff; no batch correction.'}
    (f.ROOT/'results/figures/pedigree_inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    print(json.dumps([r for r in summary if r['batch_group']=='All pairs'],indent=2))


if __name__=='__main__':main()
