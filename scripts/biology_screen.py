#!/usr/bin/env python3
"""Exploratory, rerunnable biology screens; no publication-level state/selection calls.

Example: --cohort v3a=results/prelim_1kgp_2026-10-28.tsv
         --cohort candidate_B=results/modello_pulito_regolaB_4su4_2026-10-28/cohort_1kgp.tsv
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

import frozen_snapshot as f
import pedigree_plots as ped
from validation_ratio import KEY as RATIO
from validation_plots import COLORS

ROOT=f.ROOT


def read(path):
    with path.open() as stream:return list(csv.DictReader(stream,delimiter='\t'))


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def fa_length(path):
    return sum(len(line.strip()) for line in path.read_text().splitlines() if not line.startswith('>'))


def write(path,rows):
    if rows:f.write_rows(path,rows)


def scopes(pop,kin):
    family={};component={}
    for r in sorted(pop,key=lambda r:r['sample']):
        family.setdefault(r['family'],r)
        component.setdefault(kin.get(r['sample'],r['sample']),r)
    return {'population':pop,'outside_reference':[r for r in pop if not r['reference']],
            'one_per_family':list(family.values()),'one_per_kinship_component':list(component.values())}


def values(rows,key):
    if key==RATIO:return np.array([r['ratio_bp'] for r in rows])
    column=next(t[1] for t in f.TARGETS if t[0]==key)
    return np.array([r[column] for r in rows])


def summary_stats(v):
    q=np.percentile(v,[5,25,50,75,95])
    return dict(n=len(v),p5=float(q[0]),p25=float(q[1]),median=float(q[2]),p75=float(q[3]),p95=float(q[4]),
                minimum=float(min(v)),maximum=float(max(v)))


def cv_ratio_screen(package):
    groups=defaultdict(dict)
    for r in read(package/'depth/depth_f1.tsv'):
        if r['class'] not in {'dyz1','dyz2'}:continue
        key=(r['sample'],int(r['repeat']))
        assert r['class'] not in groups[key], (key,r['class'])
        groups[key][r['class']]=r
    result=[];unpaired=[]
    for (sample,rep),pair in sorted(groups.items()):
        if set(pair)!={'dyz1','dyz2'}:
            unpaired.append(dict(sample=sample,repeat=rep,available_classes=sorted(pair)))
            continue
        a,b=pair['dyz1'],pair['dyz2'];assert a['fold']==b['fold'],(sample,rep)
        e1,e2=float(a['error_pct'])/100,float(b['error_pct'])/100
        result.append(dict(sample=sample,repeat=rep,fold=int(a['fold']),
             dyz1_signed_error_pct=100*e1,dyz2_signed_error_pct=100*e2,
             log2_ratio_error=float(np.log2((1+e1)/(1+e2)))))
    repeats={r['repeat'] for r in result};samples={r['sample'] for r in result}
    assert repeats==set(range(10)), 'Only complete ten-repeat CV packages are accepted.'
    assert len(result)==len(samples)*len(repeats)
    assert all({r['sample'] for r in result if r['repeat']==rep}==samples for rep in repeats)
    error=np.array([r['log2_ratio_error'] for r in result])
    a=np.array([r['dyz1_signed_error_pct'] for r in result]);b=np.array([r['dyz2_signed_error_pct'] for r in result])
    stats=dict(individuals=len(samples),predictions=len(result),repeats=len(repeats),
        incomplete_component_pairs_omitted=unpaired,
        empirical_abs_log2_ratio_error_p95=float(np.quantile(abs(error),.95)),
        empirical_abs_log2_ratio_error_p99=float(np.quantile(abs(error),.99)),
        signed_log2_ratio_error_median=float(np.median(error)),
        component_error_spearman=float(spearmanr(a,b).statistic),
        interpretation='Exploratory empirical error envelope conditional on the matched model references and libraries; repeated predictions are dependent. Not a calibrated prediction interval or biological exception test.')
    return result,stats


def cluster_bootstrap_median(rows,rng,repeats):
    blocks=defaultdict(list)
    for r in rows:blocks[r['component']].append(r['delta_absolute_log2_balance'])
    blocks=[np.array(v) for v in blocks.values()]
    if len(blocks)<3:return None,None
    draws=[np.median(np.concatenate([blocks[i] for i in rng.integers(0,len(blocks),len(blocks))])) for _ in range(repeats)]
    lo,hi=np.percentile(draws,[2.5,97.5]);return float(lo),float(hi)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort',action='append',required=True,help='label=TSV; outputs are kept separate by model')
    parser.add_argument('--out',type=Path,default=ROOT/'results/biology_screen_2026-10-28')
    parser.add_argument('--cv-package',type=Path,default=ROOT/'results/final_model_short_v3_2026-10-28')
    parser.add_argument('--bootstrap',type=int,default=1000)
    args=parser.parse_args();args.out=args.out.absolute();args.cv_package=args.cv_package.absolute();args.out.mkdir(parents=True,exist_ok=True)
    unit_paths=[ROOT/'data/method/hallast_motifs/DYZ1_Yq.fasta',ROOT/'data/method/hallast_motifs/DYZ2_Con.fasta']
    units=list(map(fa_length,unit_paths));unit_ratio=units[0]/units[1]
    kinpath=ROOT/'results/prelim_1kgp_2026-10-28.relatedness_samples.tsv'
    kin={r['sample']:r['kinship_component_id'] for r in read(kinpath)}
    cv_rows,cv_stats=cv_ratio_screen(args.cv_package)
    write(args.out/'joint_cv_error.tsv',cv_rows)
    dyz19_predictions=[r for r in read(args.cv_package/'depth/depth_f1.tsv') if r['class']=='dyz19']
    dyz19_accuracy=[]
    for state in ['reference_below_100kb','reference_at_least_100kb']:
        selected=[r for r in dyz19_predictions if (float(r['reference_bp'])<100000)==(state=='reference_below_100kb')]
        if not selected:continue
        errors=np.abs([float(r['error_pct']) for r in selected])
        medians=[np.median([abs(float(r['error_pct'])) for r in selected if int(r['repeat'])==rep]) for rep in range(10)]
        dyz19_accuracy.append(dict(validation_package=str(args.cv_package.relative_to(ROOT)),reference_state=state,
            n_individuals=len({r['sample'] for r in selected}),predictions=len(selected),
            mean_absolute_pct=float(np.mean(errors)),mean_repeat_median_absolute_pct=float(np.mean(medians)),
            p90_absolute_pct=float(np.percentile(errors,90)),maximum_absolute_pct=float(max(errors))))
    write(args.out/'dyz19_cv_short_vs_other.tsv',dyz19_accuracy)
    source_cohort_hash=digest(args.cv_package/'cohort_1kgp.tsv')
    lineage_stats=[];corrs=[];composition=[];short=[];short_summary=[];transmissions=[];transmission_summary=[];model_compare=[]
    model_rows={};inputs=[];counts={}
    metrics=[t[0] for t in f.TARGETS]+[RATIO]
    rng=np.random.default_rng(20261028)
    for spec in args.cohort:
        label,source=spec.split('=',1);path=Path(source).absolute()
        assert label not in model_rows,label
        raw,valid,_=f.load_snapshot(path);pop=f.population_rows(valid);model_rows[label]=pop
        inputs.append(dict(model=label,path=str(path.relative_to(ROOT)),sha256=digest(path)))
        matched_error_model=digest(path)==source_cohort_hash
        groups=scopes(pop,kin);counts[label]={k:len(v) for k,v in groups.items()}
        assert all(r['sample'] in kin for r in valid),'New individuals need updated independent kinship metadata.'
        for scope,rows in groups.items():
            clades=defaultdict(list)
            for r in rows:clades[r['lineage_display']].append(r)
            clades['ALL']=rows
            for clade,selected in clades.items():
                for key in metrics:
                    lineage_stats.append(dict(model=label,scope=scope,lineage=clade,measurement=key,**summary_stats(values(selected,key))))
            for i,a in enumerate(metrics):
                for b in metrics[i+1:]:
                    x,y=np.log2(values(rows,a)),np.log2(values(rows,b));rho=spearmanr(x,y).statistic
                    # Descriptive within-lineage association: centre both traits in each broad lineage.
                    xr=x.copy();yr=y.copy()
                    for g in {r['lineage_display'] for r in rows}:
                        ix=np.array([r['lineage_display']==g for r in rows]);xr[ix]-=np.median(x[ix]);yr[ix]-=np.median(y[ix])
                    corrs.append(dict(model=label,scope=scope,measurement_a=a,measurement_b=b,n=len(rows),
                        spearman=float(rho),spearman_after_within_lineage_log_centering=float(spearmanr(xr,yr).statistic)))
        for r in pop:
            ratio=r['ratio_bp'];balance=float(np.log2(ratio/unit_ratio));n1,n2=r['DYZ1_bp']/units[0],r['DYZ2_bp']/units[1]
            envelope=cv_stats['empirical_abs_log2_ratio_error_p95']
            flagged=(abs(balance)>envelope) if matched_error_model else ''
            composition.append(dict(model=label,sample=r['sample'],lineage=r['lineage_display'],haplogroup=r['haplogroup'],
                family=r['family'],kinship_component=kin[r['sample']],reference=r['reference'],batch=r['batch'],
                Yq12_bp=r['Yq12_bp'],DYZ1_bp=r['DYZ1_bp'],DYZ2_bp=r['DYZ2_bp'],ratio_bp=ratio,
                dyz1_consensus_unit_equivalents=n1,dyz2_consensus_unit_equivalents=n2,
                log2_unit_equivalent_ratio=balance,unit_equivalent_signed_excess=n1-n2,
                unit_equivalent_relative_excess=(n1-n2)/(n1+n2),
                joint_cv_error_model_matched=matched_error_model,
                screening_deviation_from_consensus_unit_ratio_exceeds_cv95=flagged,
                screening_stable_under_1pct_change_in_unit_ratio='',
                interpretation='Descriptive unit equivalents under the declared consensus convention; neither direct unit counts nor confirmed biological exceptions.'))
            composition[-1]['screening_stable_under_1pct_change_in_unit_ratio']=(min(abs(np.log2(ratio/(unit_ratio*t))) for t in [.99,1,1.01])>envelope) if matched_error_model else ''
        for cutoff in [100000,125000,150000,175000]:
            selected=[r for r in pop if r['DYZ19_bp']<cutoff]
            for group in sorted({r['lineage_display'] for r in selected}):
                members=[r for r in selected if r['lineage_display']==group]
                short_summary.append(dict(model=label,descriptive_cutoff_bp=cutoff,lineage=group,n=len(members),
                    outside_reference=sum(not r['reference'] for r in members),
                    distinct_kinship_components=len({kin[r['sample']] for r in members}),
                    minimum_bp=min(r['DYZ19_bp'] for r in members),maximum_bp=max(r['DYZ19_bp'] for r in members)))
            if cutoff==100000:
                for r in selected:short.append(dict(model=label,sample=r['sample'],haplogroup=r['haplogroup'],lineage=r['lineage_display'],
                    DYZ19_bp=r['DYZ19_bp'],Yq12_bp=r['Yq12_bp'],y_dose=r['y_dose'],reference=r['reference'],
                    family=r['family'],kinship_component=kin[r['sample']],y_status=r['y_status'],
                    interpretation='Below a declared descriptive 100 kb screen; not a deletion/breakpoint call.'))
        index={r['sample']:r for r in valid};pair_rows=[r for r in ped.collect(valid) if r['target']==RATIO]
        local=[]
        for p in pair_rows:
            a,b=index[p['father']],index[p['son']]
            ab,bb=np.log2(a['ratio_bp']/unit_ratio),np.log2(b['ratio_bp']/unit_ratio)
            assert kin[a['sample']]==kin[b['sample']]
            local.append(dict(model=label,father=p['father'],son=p['son'],batch_group=p['batch_group'],component=kin[a['sample']],
                father_log2_unit_balance=float(ab),son_log2_unit_balance=float(bb),
                delta_absolute_log2_balance=float(abs(bb)-abs(ab)),delta_log2_ratio=float(bb-ab),
                both_population_included=a['population_included'] and b['population_included'],
                reference_involved=a['reference'] or b['reference']))
        transmissions.extend(local)
        for scope in ['all_profiles','population_profiles','population_outside_reference']:
            selected=[r for r in local if scope=='all_profiles' or
                      (r['both_population_included'] and (scope!='population_outside_reference' or not r['reference_involved']))]
            for group in ['All pairs',*ped.GROUPS]:
                s=[r for r in selected if group=='All pairs' or r['batch_group']==group]
                if not s:continue
                delta=np.array([r['delta_absolute_log2_balance'] for r in s]);lo,hi=cluster_bootstrap_median(s,rng,args.bootstrap)
                transmission_summary.append(dict(model=label,scope=scope,batch_group=group,n_pairs=len(s),
                    independent_kinship_components=len({r['component'] for r in s}),
                    median_delta_absolute_log2_balance=float(np.median(delta)),fraction_closer_to_consensus_unit_ratio=float(np.mean(delta<0)),
                    family_bootstrap_median_ci_low=lo,family_bootstrap_median_ci_high=hi,
                    median_signed_log2_ratio_change=float(np.median([r['delta_log2_ratio'] for r in s])),
                    interpretation='Exploratory paired contrast; biological convergence is not established without a measurement-error/regression-to-mean null and independent batch calibration.'))
    labels=list(model_rows)
    if len(labels)>1:
        base={r['sample']:r for r in model_rows[labels[0]]}
        for name in labels[1:]:
            other={r['sample']:r for r in model_rows[name]};ids=sorted(base.keys()&other.keys())
            for key in metrics:
                x=values([base[i] for i in ids],key);y=values([other[i] for i in ids],key)
                delta=100*(y/x-1)
                model_compare.append(dict(baseline=labels[0],candidate=name,measurement=key,n_matched=len(ids),
                    median_signed_change_pct=float(np.median(delta)),p95_absolute_change_pct=float(np.quantile(abs(delta),.95)),
                    maximum_absolute_change_pct=float(max(abs(delta)))))
    for name,rows in [('lineage_distributions',lineage_stats),('trait_correlations',corrs),('composition_screen',composition),
                      ('dyz19_below_100kb',short),('dyz19_cutoff_sensitivity',short_summary),('father_son_balance_pairs',transmissions),
                      ('father_son_balance_summary',transmission_summary),('model_sensitivity',model_compare)]:write(args.out/(name+'.tsv'),rows)
    # A compact first look: repeat analyses from input tables, without changing the manuscript.
    fig,axes=plt.subplots(1,3,figsize=(12,3.6))
    first=labels[0];rows=model_rows[first]
    axes[0].scatter([r['Yq12_bp']/1e6 for r in rows],[np.log2(r['ratio_bp']/unit_ratio) for r in rows],s=7,alpha=.4,color='#557c90')
    axes[0].axhline(0,color='#697782',ls='--',lw=.7);axes[0].set_xlabel('Yq12 (Mb)');axes[0].set_ylabel('log₂ unit-equivalent DYZ1:DYZ2')
    axes[0].set_title('Composition: declared unit convention',fontsize=10)
    common_edges=np.histogram_bin_edges(np.concatenate([np.log2([r['ratio_bp']/unit_ratio for r in model_rows[label]]) for label in labels]),bins=45)
    for i,label in enumerate(labels):
        v=np.log2([r['ratio_bp']/unit_ratio for r in model_rows[label]])
        axes[1].hist(v,bins=common_edges,density=True,histtype='step',lw=1.3,label=label)
    axes[1].axvline(0,color='#697782',ls='--',lw=.7);axes[1].set_xlabel('log₂ unit-equivalent DYZ1:DYZ2');axes[1].set_ylabel('Density');axes[1].legend(frameon=False,fontsize=8)
    axes[1].set_title('Sensitivity to model',fontsize=10)
    for i,group in enumerate(ped.GROUPS):
        selected=[r for r in transmissions if r['model']==first and r['batch_group']==group and r['both_population_included']]
        vals=[r['delta_absolute_log2_balance'] for r in selected]
        axes[2].scatter(i+rng.uniform(-.15,.15,len(vals)),vals,s=7,alpha=.4,color=COLORS[i])
        axes[2].plot([i-.18,i+.18],[np.median(vals)]*2,color='#253f4b',lw=2)
    axes[2].axhline(0,color='#697782',ls='--',lw=.7);axes[2].set_xticks(range(3),['Cross batch','Same batch','Reverse'],rotation=15)
    axes[2].set_ylabel('Δ |log₂ unit-equivalent ratio|');axes[2].set_title('Father–son: descriptive contrast',fontsize=10)
    for ax in axes:ax.spines[['right','top']].set_visible(False);ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=.065,right=.99,bottom=.21,top=.87,wspace=.36)
    fig.text(.5,.025,'Exploratory screens; unit equivalents use consensus lengths, not directly annotated unit counts. No selection or convergence claim.',ha='center',fontsize=8)
    for ext in ['pdf']:fig.savefig(args.out/('biology_first_look.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    manifest=dict(inputs=inputs,counts=counts,unit_convention=dict(dyz1_bp=units[0],dyz2_bp=units[1],
        sources=[dict(path=str(p.relative_to(ROOT)),sha256=digest(p)) for p in unit_paths],
        caveat='Consensus-length unit equivalents. Direct complete-unit counts/boundary sensitivity are still needed to validate an equilibrium-at-one interpretation.'),
        joint_error_screen=cv_stats,auxiliary_sources=[dict(path=str(p.relative_to(ROOT)),sha256=digest(p)) for p in [kinpath,args.cv_package/'depth/depth_f1.tsv',ROOT/'data/catalog_d376/fs_pairs.tsv',ROOT/'data/catalog_d376/ped.txt',ROOT/'scripts/biology_screen.py']],bootstrap=dict(repeats=args.bootstrap,resampling_unit='public kinship connected component',seed=20261028),
        readiness=dict(descriptive_population='computed; rerunnable',composition_screen='computed; exploratory convention and empirical model-matched CV error envelope only',
        dyz19_states='computed descriptive screens with cutoff sensitivity; independent sequence confirmation pending',
        father_son_balance='computed exploratory paired contrasts and kinship bootstrap; noise-adjusted biological inference pending',
        phylogenetic_models='Provisional v2 tree supports exploratory relative comparisons now. Robust inference needs final topology/support; time calibration is needed specifically for dates and rates per generation.',
        cross_cohort_replication='partial cohort tables and matched-sample summaries available; full cohort, per-library provenance and model-matched validation pending',
        structural_mechanism='cannot establish breakpoints or donor-versus-cell-line origin from abundance tables'),
        notes=['No dosage cutoff; existing arm/Yp flags excluded from population summaries and retained in pedigree all-profile sensitivity.',
               'Candidate model point estimates are sensitivity analyses; no error envelope is transferred across mismatched model/cohort packages.',
               'No evolutionary, structural or selection conclusions have been inserted into the manuscript.'])
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(output=str(args.out.relative_to(ROOT)),counts=counts,unit_convention=units,joint_error_screen=cv_stats,readiness=manifest['readiness']),indent=2))


if __name__=='__main__':main()
