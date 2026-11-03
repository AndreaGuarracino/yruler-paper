"""All-target father–son comparisons, retaining batch and reference status."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
import numpy as np
import figures_1kgp as f
from validation_plots import COLORS,INK
from validation_ratio import KEY as RATIO, LABEL as RATIO_LABEL, COLOR as RATIO_COLOR

MEASUREMENTS=[(key,col,label,scale,unit,color) for (key,col,label,scale,unit),color in zip(f.TARGETS,COLORS)]
MEASUREMENTS.append((RATIO,'ratio_bp',RATIO_LABEL,1,'',RATIO_COLOR))
LOOKUP={row[0]:row for row in MEASUREMENTS}
GROUPS=['Cross batch','Same batch','Reverse cross batch']
MARKERS=['o','^','D']


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
    ax.grid(color='#edf1f3',linewidth=.5);ax.set_axisbelow(True)
    ax.tick_params(labelsize=6.5,length=2.5)


def scatter(ax,rows,key,compact=False):
    _,_,label,scale,unit,color=LOOKUP[key];s=[r for r in rows if r['target']==key]
    upper=1.08*max(max(r['father_value'],r['son_value']) for r in s)
    ax.plot([0,upper],[0,upper],'--',color='#6f8490',linewidth=.65,zorder=1)
    for group,marker in zip(GROUPS,MARKERS):
        for ref in (False,True):
            subset=[r for r in s if r['batch_group']==group and
                    (r['father_reference']=='yes' or r['son_reference']=='yes')==ref]
            ax.scatter([r['father_value'] for r in subset],[r['son_value'] for r in subset],
                marker=marker,s=(7 if group=='Cross batch' else 16) if compact else (12 if group=='Cross batch' else 25),
                facecolors='none' if ref else color,edgecolors=color if ref else 'white',
                linewidths=.5 if ref else .2,alpha=.75 if ref else .50,zorder=3)
    worst=sorted(s,key=lambda r:-abs(r['signed_difference_pct']))[:1 if compact else 2]
    if compact:
        for r in s:
            if r['son']=='HG02966' and abs(r['signed_difference_pct'])>=20 and r not in worst:worst.append(r)
    for rank,r in enumerate(worst):
        if abs(r['signed_difference_pct'])<10:continue
        x,y=r['father_value'],r['son_value'];right=x>upper*.55
        offset=(-10-rank*8) if y>upper*.83 else (4+rank*9)
        ax.annotate(r['son'],(x,y),xytext=(-3 if right else 3,offset),textcoords='offset points',
                    fontsize=5.5 if compact else 6,ha='right' if right else 'left',color=INK,
                    bbox=dict(facecolor='white',edgecolor='none',pad=.1),
                    arrowprops=dict(arrowstyle='-',color=color,linewidth=.4))
    med=np.median([abs(r['signed_difference_pct']) for r in s])
    ax.text(.035,.965,f'Median |Δ| {med:.2f}%',transform=ax.transAxes,va='top',fontsize=6 if compact else 7.5,color=INK)
    ax.set_xlim(0,upper);ax.set_ylim(0,upper);ax.set_aspect('equal',adjustable='box')
    ax.xaxis.set_major_locator(MaxNLocator(3 if compact else 4));ax.yaxis.set_major_locator(MaxNLocator(3 if compact else 4))
    suffix=f' ({unit})' if unit else ''
    ax.set_xlabel('Parent'+suffix,fontsize=7 if compact else 8)
    ax.set_ylabel('Offspring'+suffix,fontsize=7 if compact else 8,labelpad=1)
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
            ax.scatter(values[mask],ys[mask],s=7 if compact else 13,marker=MARKERS[i],
                facecolors='none' if ref else color,edgecolors=color if ref else 'none',
                linewidths=.45,alpha=.65 if ref else .35,zorder=3)
        lo,median,hi=np.percentile(values,[25,50,75])
        ax.plot([lo,hi],[i,i],color=color,linewidth=3,zorder=4)
        ax.plot(median,i,'|',color=INK,markersize=8,zorder=5)
        ax.annotate(f'{median:+.2f}%',(median,i),xytext=(3,6),textcoords='offset points',fontsize=6 if compact else 7,color=INK)
    ax.axvline(0,color='#738793',linestyle='--',linewidth=.65)
    ax.set_xscale('symlog',linthresh=5);ax.set_xlim(-maximum,maximum)
    tick_pool=([-100,0,100] if maximum>=100 else [-20,0,20] if maximum>=20 else [-10,0,10]) if compact else [-100,-20,-5,0,5,20,100]
    ticks=[x for x in tick_pool if abs(x)<=maximum]
    ax.set_xticks(ticks,[str(x) for x in ticks]);ax.set_ylim(len(groups)-.55,-.55)
    short=['Cross\nbatch','Same\nbatch','Reverse\ncross batch'] if compact else GROUPS
    ax.set_yticks(range(len(groups)),short[:len(groups)],fontsize=6 if compact else 7)
    ax.tick_params(axis='y',length=0);ax.set_xlabel('Offspring/parent difference (%)',fontsize=7 if compact else 8)
    style(ax);ax.spines['left'].set_visible(False)


def legend(fig,rows,y=1.0,compact=False):
    first=[r for r in rows if r['target']=='dyz3'];counts=Counter(r['batch_group'] for r in first)
    h=[Line2D([],[],marker=m,linestyle='none',color='#58717e',markersize=4,
       label=f'{g} (n={counts[g]})') for g,m in zip(GROUPS,MARKERS)]
    h.append(Line2D([],[],marker='o',linestyle='none',markerfacecolor='none',markeredgecolor='#58717e',markersize=4,
        label='Open: reference involved'))
    fig.legend(handles=h,loc='upper center',bbox_to_anchor=(.525,y),ncol=4,frameon=False,
               fontsize=6.5 if compact else 8,handletextpad=.4,columnspacing=1.1)


def supplement(rows,out):
    with PdfPages(out/'figureS8_all_pedigrees.pdf') as pdf:
        for kind in ('scatter','differences'):
            fig,axes=plt.subplots(3,3,figsize=(11,10 if kind=='scatter' else 8.5))
            for i,(ax,meta) in enumerate(zip(axes.flat,MEASUREMENTS)):
                key,_,label,_,_,color=meta
                if kind=='scatter':scatter(ax,rows,key)
                else:differences(ax,rows,key)
                ax.set_title(f'{chr(97+i)}   {label}',loc='left',color=color,fontsize=11,pad=7)
            for ax in axes.flat[len(MEASUREMENTS):]:ax.set_visible(False)   # unused grid cells
            legend(fig,rows,y=.985)
            fig.subplots_adjust(left=.075,right=.98,bottom=.12,top=.91,wspace=.42,hspace=.49)
            fig.text(.5,.018,f'{len(rows)//len(MEASUREMENTS)} parent–offspring pairs; all profiles retained. Ratio is DYZ1/DYZ2. No batch correction.',ha='center',fontsize=8,color=INK)
            if kind=='differences':fig.text(.5,.04,'Median and interquartile range; difference axes linear within ±5% and logarithmic outside.',ha='center',fontsize=8,color=INK)
            fig.savefig(out/f'figureS8_pedigree_{kind}.pdf',dpi=170,bbox_inches='tight')
            pdf.savefig(fig,bbox_inches='tight');plt.close(fig)


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    _,valid,_=f.load_snapshot();rows=collect(valid);out=f.ROOT/'manuscript/figures'
    assert len(rows)%len(MEASUREMENTS)==0; npairs=len(rows)//len(MEASUREMENTS); print('father-son pairs',npairs,flush=True)
    # Verify the two existing main-panel comparisons before extending to other measurements.
    previous_path=out/'figure1_pedigree_pairs.tsv'
    previous=list(csv.DictReader(previous_path.open(),delimiter='\t')) if previous_path.exists() else []
    lookup={(r['father'],r['son'],r['target']):r for r in rows}
    for r in previous:
        new=lookup[(r['father'],r['son'],r['target'].removesuffix('_bp'))]
        assert abs(float(r['signed_difference_pct'])-new['signed_difference_pct'])<1e-10
    summary=stats(rows);f.write_rows(out/'figureS8_pedigree_pairs.tsv',rows);f.write_rows(out/'figureS8_pedigree_summary.tsv',summary)
    supplement(rows,out)
    inputs={'cohort':str(f.DEFAULT_SNAPSHOT.relative_to(f.ROOT)),'cohort_sha256':hashlib.sha256(f.DEFAULT_SNAPSHOT.read_bytes()).hexdigest(),
           'pairs_per_measurement':npairs,'measurement_records':len(rows),'measurements':[m[2] for m in MEASUREMENTS],
           'batch_counts':dict(Counter(r['batch_group'] for r in rows if r['target']=='dyz3')),
           'reference_involved_pairs':sum(r['father_reference']=='yes' or r['son_reference']=='yes' for r in rows if r['target']=='dyz3'),
           'existing_main_differences_matched':len(previous),'profile_policy':'All complete profiles retained, including flagged arm imbalance; no dosage cutoff; no batch correction.'}
    (out/'figureS8_pedigree_inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    print(json.dumps([r for r in summary if r['batch_group']=='All pairs'],indent=2))


if __name__=='__main__':main()
