#!/usr/bin/env python3
"""Exploratory covariance profiles and minimum DYZ19 state changes on the Y tree.

No dates, mutation rates, branch-specific deletion mechanisms, or corrected
measurement-error inference are estimated from this provisional SNP tree.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from threadpoolctl import threadpool_limits
import frozen_snapshot as f
from biology_screen import read,scopes,values,write
from validation_ratio import KEY as RATIO

METRICS=[r[0] for r in f.TARGETS]+[RATIO]
LABELS={r[0]:r[2] for r in f.TARGETS}|{RATIO:'DYZ1:DYZ2'}
CAP=.9999


def covariance(tree,names):
    index={s:i for i,s in enumerate(names)};c=np.zeros((len(names),len(names)))
    branches=[]
    def visit(node,height=0,root=False):
        length=0 if root else node['length'];assert length>=0
        height+=length;branches.append(length)
        if not node['children']:
            if node['name'] not in index:return []
            i=index[node['name']];c[i,i]=height;return [i]
        children=[visit(ch,height) for ch in node['children']]
        for j,a in enumerate(children):
            for b in children[:j]:
                c[np.ix_(a,b)]=height;c[np.ix_(b,a)]=height
        return [i for ids in children for i in ids]
    seen=visit(tree,root=True);assert set(seen)==set(range(len(names)))
    diag=np.diag(c).copy();assert np.all(diag>0)
    z=c/np.sqrt(np.outer(diag,diag));eig,q=eigh(z,check_finite=False)
    assert eig.min()>-1e-8
    return diag,np.maximum(eig,0),q,dict(tips=len(names),minimum_eigenvalue_before_roundoff_clipping=float(eig.min()),zero_length_source_edges=sum(b==0 for b in branches),
        min_root_tip_substitutions_per_site=float(diag.min()),max_root_tip_substitutions_per_site=float(diag.max()))


def profile(diag,eig,q,y):
    """ML of V(lambda)=lambda*C+(1-lambda)*diag(C); spectral whitening."""
    n=len(y);u=q.T@(1/np.sqrt(diag));z=q.T@(y/np.sqrt(diag));ld0=np.log(diag).sum()
    def fit(lam,return_values=False):
        d=1+lam*(eig-1);assert np.all(d>0)
        mu=np.sum(u*z/d)/np.sum(u*u/d);res=z-mu*u
        var=np.sum(res*res/d)/n
        ll=-.5*(n*(np.log(2*np.pi)+1+np.log(var))+np.log(d).sum()+ld0)
        return (ll,mu,var) if return_values else -ll
    grid=np.linspace(0,CAP,41);scores=[fit(x) for x in grid];i=int(np.argmin(scores))
    lo,hi=grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]
    opt=minimize_scalar(fit,bounds=(lo,hi),method='bounded',options={'xatol':1e-8})
    candidates=[(fit(0),0),(fit(CAP),CAP),(opt.fun,opt.x)]
    _,lam=min(candidates);ll,mu,var=fit(lam,True)
    null=fit(0,True)[0];near_bm=fit(CAP,True)[0]
    return dict(lambda_ml=float(lam),mean_log2=float(mu),variance_per_shared_branch_unit=float(var),
        log_likelihood=float(ll),log_likelihood_lambda0=float(null),log_likelihood_lambda_cap=float(near_bm),
        delta_aic_independent_minus_lambda=float(2*(ll-null)-2),
        at_upper_bound=bool(lam>CAP-1e-6),measurement_error_adjusted=False),[(float(x),float(-fit(x))) for x in grid]


def bivariate_profile(diag,eig,q,y):
    n,p=y.shape;u=q.T@(1/np.sqrt(diag));z=q.T@(y/np.sqrt(diag)[:,None]);ld0=np.log(diag).sum()
    def fit(lam,details=False):
        d=1+lam*(eig-1);mu=np.sum(u[:,None]*z/d[:,None],axis=0)/np.sum(u*u/d)
        res=z-u[:,None]*mu;s=(res/d[:,None]).T@res/n
        sign,ld=np.linalg.slogdet(s)
        if sign<=0:return np.inf
        ll=-.5*(n*p*(np.log(2*np.pi)+1)+n*ld+p*(np.log(d).sum()+ld0))
        rho=s[0,1]/np.sqrt(s[0,0]*s[1,1])
        return (ll,rho,s) if details else -ll
    grid=np.linspace(0,CAP,31);i=int(np.argmin([fit(x) for x in grid]))
    opt=minimize_scalar(fit,bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),method='bounded')
    _,lam=min([(fit(0),0),(fit(CAP),CAP),(opt.fun,opt.x)])
    ll,rho,cross=fit(lam,True)
    return dict(shared_lambda_ml=float(lam),residual_correlation=float(rho),
        residual_variance_a=float(cross[0,0]),residual_variance_b=float(cross[1,1]),residual_covariance=float(cross[0,1]),slope_b_on_a=float(cross[0,1]/cross[0,0]),
        correlation_lambda0=float(fit(0,True)[1]),correlation_lambda_cap=float(fit(CAP,True)[1]),
        log_likelihood=float(ll),interpretation='Shared-branch covariance sensitivity; no measurement-error correction and no proof of selection or evolutionary coupling.')


def sankoff(tree,states):
    """Unit-cost, two-state parsimony, valid for multifurcations and unknown roots."""
    nodes=[]
    def visit(node,path='root'):
        if not node['children']:
            if node['name'] not in states:return None
            state=states[node['name']];cost=np.array([0 if s==state else np.inf for s in [0,1]])
            item=dict(path=path,name=node['name'],length=node['length'],children=[],cost=cost,tips=[node['name']])
        else:
            children=[ch for i,c in enumerate(node['children']) if (ch:=visit(c,path+'.'+str(i))) is not None]
            if not children:return None
            cost=np.zeros(2)
            for state in [0,1]:cost[state]=sum(min(ch['cost'][0]+(state!=0),ch['cost'][1]+(state!=1)) for ch in children)
            item=dict(path=path,name=node['name'],length=node['length'],children=children,cost=cost,tips=[t for c in children for t in c['tips']])
        nodes.append(item);return item
    root=visit(tree);minimum=min(root['cost']);root_states=set(np.flatnonzero(root['cost']==minimum));edges=[]
    def forward(node,allowed):
        for child in node['children']:
            pairs=set()
            for state in allowed:
                costs=child['cost']+np.array([state!=0,state!=1])
                for target in np.flatnonzero(costs==min(costs)):pairs.add((int(state),int(target)))
            change=[a!=b for a,b in pairs]
            if any(change):edges.append(dict(child_path=child['path'],branch_length_substitutions_per_site=child['length'],
                child_source_label=child['name'],descendant_tips=len(child['tips']),short_descendants=sum(states[t] for t in child['tips']),
                forced_in_all_optimal_mappings=all(change),allowed_parent_child_states=';'.join(f'{a}>{b}' for a,b in sorted(pairs)),
                example_tips=';'.join(child['tips'][:8])))
            forward(child,{b for a,b in pairs})
    forward(root,root_states)
    return dict(n_tips=len(states),n_short=sum(states.values()),minimum_changes=int(minimum),
                cost_if_root_ordinary=int(root['cost'][0]),cost_if_root_short=int(root['cost'][1]),
                optimal_root_states=';'.join(map(str,sorted(root_states))),
                possibly_changing_edges=len(edges),forced_changing_edges=sum(r['forced_in_all_optimal_mappings'] for r in edges)),edges


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cohort',action='append',required=True)
    p.add_argument('--out',type=Path,default=f.ROOT/'results/biology_followup_2026-10-28/phylogeny')
    p.add_argument('--tree',type=Path,default=f.TREE_SOURCE);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    kin={r['sample']:r['kinship_component_id'] for r in read(f.ROOT/'results/estimates_1kgp.relatedness_samples.tsv')}
    tree=f.parse_newick(args.tree);cache={};summary=[];profiles=[];pairs=[];states_summary=[];states_edges=[];audits=[];inputs=[]
    with threadpool_limits(limits=2):
        for spec in args.cohort:
            label,path=spec.split('=',1);path=Path(path);_,valid,_=f.load_snapshot(path);pop=[r for r in f.population_rows(valid) if r['y_status']!='par2_gain']
            inputs.append(dict(model=label,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            group=scopes(pop,kin)
            unrelated_out=scopes([r for r in pop if not r['reference']],kin)['one_per_kinship_component']
            group['one_per_kinship_outside_reference']=unrelated_out
            for scope in ['one_per_kinship_component','one_per_kinship_outside_reference']:
                rows=sorted(group[scope],key=lambda r:r['sample']);names=tuple(r['sample'] for r in rows)
                if names not in cache:cache[names]=covariance(tree,names)
                diag,eig,q,audit=cache[names];audits.append(dict(model=label,scope=scope,**audit))
                y=np.column_stack([np.log2(values(rows,key)) for key in METRICS])
                for j,key in enumerate(METRICS):
                    s,grid=profile(diag,eig,q,y[:,j]);summary.append(dict(model=label,scope=scope,measurement=key,n=len(rows),**s))
                    profiles.extend(dict(model=label,scope=scope,measurement=key,lambda_value=x,log_likelihood=ll) for x,ll in grid)
                for i,a in enumerate(METRICS):
                    for j in range(i):pairs.append(dict(model=label,scope=scope,measurement_a=METRICS[j],measurement_b=a,n=len(rows),**bivariate_profile(diag,eig,q,y[:,[j,i]])))
            for scope in ['population','one_per_kinship_component']:
                for cutoff in [100000,125000,150000]:
                    mapping={r['sample']:int(r['DYZ19_bp']<cutoff) for r in group[scope]}
                    result,edges=sankoff(tree,mapping)
                    states_summary.append(dict(model=label,scope=scope,descriptive_cutoff_bp=cutoff,**result))
                    states_edges.extend(dict(model=label,scope=scope,descriptive_cutoff_bp=cutoff,**r) for r in edges)
            print(label,'complete',flush=True)
    for stem,rows in [('phylogenetic_signal',summary),('lambda_profiles',profiles),('bivariate_covariance',pairs),('dyz19_parsimony',states_summary),('dyz19_state_edges',states_edges)]:write(args.out/(stem+'.tsv'),rows)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    labels=list(dict.fromkeys(r['model'] for r in summary))
    for i,label in enumerate(labels):
        s=[r for r in summary if r['model']==label and r['scope']=='one_per_kinship_component']
        axes[0].plot([r['lambda_ml'] for r in s],np.arange(len(s))+(i-.5)*.12,'o',label=label,ms=4)
    axes[0].set_yticks(range(len(METRICS)),[LABELS[k] for k in METRICS]);axes[0].invert_yaxis();axes[0].set_xlim(-.02,1.02)
    axes[0].set_xlabel('Profiled shared-branch covariance λ');axes[0].legend(frameon=False,fontsize=8)
    s=[r for r in states_summary if r['scope']=='population' and r['descriptive_cutoff_bp']==100000]
    axes[1].bar(range(len(s)),[r['minimum_changes'] for r in s],color=['#557d92','#c48449'][:len(s)])
    axes[1].set_xticks(range(len(s)),[r['model'] for r in s]);axes[1].set_ylabel('Minimum state changes');axes[1].set_title('DYZ19 below 100 kb; provisional topology',fontsize=10)
    axes[1].set_ylim(0,max(r['minimum_changes'] for r in s)+1)
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=.12,right=.98,bottom=.2,top=.93,wspace=.48)
    fig.text(.5,.035,'Exploratory: provisional SNP tree; uncorrected measurement error; substitutions/site, not time. State changes are not deletion events.',ha='center',fontsize=8)
    for ext in ['pdf']:fig.savefig(args.out/('phylogeny_first_look.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    manifest=dict(inputs=inputs,tree=dict(path=str(args.tree),sha256=hashlib.sha256(args.tree.read_bytes()).hexdigest(),status='A0-rooted v4, branch lengths retained; no clock calibration'),
        covariance_audits=audits,lambda_range=[0,CAP],method='Gaussian ML with Cij=shared root-to-MRCA branch length, diagonal=root-to-tip length; V(lambda)=lambda*C+(1-lambda)*diag(C). Spectral whitening; all means/variances reprofiled.',
        numerical_note='Upper bound below one avoids singular covariance for identical SNP haplotypes. Lambda=1 is not asserted to have been fitted.',
        limitations=['No separation of biological residual variance from measurement error.','No dates, rates per generation, selection claims, or formal tests of satellite turnover.','Parsimony treats the declared length screen as a binary character with equal costs. The minimum refers only to this tree and is not a count of independently proven deletions.'],
        primary_method_reference='https://www.phytools.org/phytools-manual.pdf',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':main()
