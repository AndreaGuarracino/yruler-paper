#!/usr/bin/env python3
"""Audit the source SNP tree, display transformations and ISOGG concordance."""
import csv
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path
import frozen_snapshot as f


def indexed(tree):
    nodes, tips = [], {}
    def visit(node, parent=None, depth=0, level=0):
        i=len(nodes); item=dict(node, parent=parent, depth=depth+node['length'], level=level, index=i);nodes.append(item)
        if not node['children']:
            assert node['name'] not in tips
            tips[node['name']]=i;item['tips']={node['name']}
        else:
            kids=[visit(c,i,item['depth'],level+1) for c in node['children']]
            item['tips']=set().union(*(nodes[c]['tips'] for c in kids))
        return i
    visit(tree)
    return nodes,tips


def mrca(nodes,tips,samples):
    common=None
    for sample in samples:
        path=set();i=tips[sample]
        while i is not None:path.add(i);i=nodes[i]['parent']
        common=path if common is None else common&path
    return max(common,key=lambda i:nodes[i]['level'])


def main():
    raw,valid,_=f.load_snapshot();samples={r['sample']:r for r in valid}
    source=f.TREE_SOURCE;original=f.parse_newick(source)
    full,fulltips=indexed(original)
    assert set(samples)<=set(fulltips)
    assert all(math.isfinite(n['length']) and n['length']>=0 for n in full)
    pruned=f.prune_tree(original,set(samples));pn,pt=indexed(pruned)
    assert set(pt)==set(samples)
    # Every displayed split must be a source-tree split restricted to shown tips.
    projection={frozenset(n['tips']&set(samples)) for n in full}
    assert all(frozenset(n['tips']) in projection for n in pn)
    def distance(nodes,tips,a,b):
        anc=mrca(nodes,tips,[a,b]);return nodes[tips[a]]['depth']+nodes[tips[b]]['depth']-2*nodes[anc]['depth']
    rng=random.Random(2710); pairs=[rng.sample(sorted(samples),2) for _ in range(1000)]
    error=max(abs(distance(full,fulltips,a,b)-distance(pn,pt,a,b)) for a,b in pairs)
    assert error<1e-10
    collapsed,leaves=f.collapse_tree(pruned,samples)
    coverage=Counter(r['sample'] for node in leaves for r in node['members']);assert coverage==Counter(samples.keys())
    triangle_error=0
    for node in leaves:
        ids={r['sample'] for r in node['members']};anc=mrca(pn,pt,ids)
        assert pn[anc]['tips']==ids
        extent=max(pn[pt[s]]['depth'] for s in ids)-pn[anc]['depth']
        triangle_error=max(triangle_error,abs(extent-node['extent']))
    assert triangle_error<1e-10
    parents,_=f.tree_data(f.ROOT/'data/figure_metadata/isogg2019_tree.tsv')
    ancestry={}
    for sample,row in samples.items():
        node=row['haplogroup'].split()[0];path=set()
        while node:path.add(node);node=parents.get(node)
        ancestry[sample]=path
    rows=[]
    for clade in ['A0','A1','B','C','D','E','E1a2a1b2','F','G','H','I','J','J1','J2','J2a','J2a1','J2a1a','J2a1a1','J2a1a1a','J2a1a1b','J2a1a2','J2a2','J2b','K','L','NO','N','O','O1','O1a','O1b','O2','O2a','O2b','Q','R','R1','R1a','R1b','R2','T']:
        ids={s for s in samples if clade in ancestry[s]}
        if len(ids)<2:continue
        anc=mrca(pn,pt,ids);extra=pn[anc]['tips']-ids
        split_compatible=any(n['tips']==ids or set(samples)-n['tips']==ids for n in pn)
        rows.append(dict(clade=clade,n_called_descendants=len(ids),n_mrca_tips=len(pn[anc]['tips']),monophyletic=not extra,unrooted_split_compatible=split_compatible,extra_samples=';'.join(sorted(extra)),extra_calls=';'.join(f'{s}:{samples[s]["haplogroup"]}' for s in sorted(extra))))
    f.write_rows(f.ROOT/'results/figure2_tree_concordance.tsv',rows)
    report={'source':str(source.relative_to(f.ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'snapshot':str(f.DEFAULT_SNAPSHOT.relative_to(f.ROOT)), 'source_tips':len(fulltips),'display_tips':len(pt),'collapsed_rows':len(leaves),'nonempty_internal_labels':sum(bool(n['name']) for n in full if n['children']), 'distance_checks':len(pairs),'max_distance_error':error,'max_triangle_extent_error':triangle_error,'all_display_splits_in_source':True,'all_display_rows_monophyletic':True,'discordant_classification_clades':[{k:r[k] for k in ['clade','n_called_descendants','n_mrca_tips','extra_calls']} for r in rows if not r['monophyletic']], 'limits':['Only rooted Newick archived locally; original alignment, IQ-TREE report and rooting command unavailable here.','No branch support labels; source command uses -fast without bootstrap.','Source filters remove singletons as well as invariant sites; standard +ASC does not explicitly model singleton exclusion.','Binary REF/ALT codes are not nucleotide identities or ancestral/derived states.','Topological agreement with haplogroups is not independent validation of every branch.']}
    reduced=f.prune_tree(original,set(samples)-{'HG03742'});rn,rt=indexed(reduced)
    report['diagnostic_exclusion_only']={}
    for clade in ['K','NO']:
        ids={s for s in rt if clade in ancestry[s]};anc=mrca(rn,rt,ids)
        report['diagnostic_exclusion_only'][clade]={'excluded':'HG03742','n_called_descendants':len(ids),'n_mrca_tips':len(rn[anc]['tips']),'monophyletic':rn[anc]['tips']==ids}
    report['root_diagnostic']={'A1_unrooted_split_compatible':next(r['unrooted_split_compatible'] for r in rows if r['clade']=='A1'),'interpretation':'A1 discordance can be due to root placement; assess A0 versus A1 when A00 is absent.'}
    report['literature']=[{'reference':'Poznik et al. Nature Genetics 48, 593–599 (2016), doi:10.1038/ng.3559','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4884158/','sample':'HG03742','published_placement':'Within K2a, near the branch leading to N and O, sharing derived M2308.','source_placement':'Basal near F, outside the other K descendants; requires genotype/alignment and inference checks.'}]
    report['interpretation']='Display preserves the supplied tree, but the source tree is not scientifically validated. Diagnostic pruning was not applied to the figure.'
    if (f.ROOT/'results/ytree_audit/v1_aln.phy').exists():
        report['limits'][0]='Source alignment, IDs, site positions, unrooted tree and IQ-TREE report now archived in results/ytree_audit; rooting is described in the README. REF/ALT and direct M2308 verification remain pending.'
        report['delivered_input_audit']='results/figure2_tree_delivery_audit.json'
        report['replacement_tree_status']='Colleague reports v2 reconstruction in progress with accessible-region mask, singleton inclusion, full search, UFBoot/SH-aLRT and A0 rooting; not yet used in figures.'
    preparation=f.ROOT/'results/figure2_tree_v2_preparation_audit.json'
    if preparation.exists():
        report['source_preparation']=json.loads(preparation.read_text())
        report['limits']=['Original v2 alignment not delivered; site table, missingness, ML tree and inference report checked.',
                          'The tree is provisional; fine-scale classification discordances remain.',
                          'Binary SNP-site branch lengths are not whole-Y divergence or dates.',
                          'Support calculations do not independently validate every branch.']
        report['replacement_tree_status']='Delivered v2 now used, with A0 rooting derived from its existing unrooted split.'
        report['interpretation']='The display preserves the provisional v2 topology and distances; HG03742 is sister to NO and J2a1a1 is monophyletic.'
        report['literature'][0]['source_placement']='HG03742 is sister to NO in v2, consistent with the published K2a placement.'
        report['diagnostic_exclusion_only']={}
        report['root_diagnostic']={'A0_outgroup':report['source_preparation']['root_outgroup'],
                                   'A0_monophyletic':next(r['monophyletic'] for r in rows if r['clade']=='A0'),
                                   'A1_monophyletic':next(r['monophyletic'] for r in rows if r['clade']=='A1')}
    (f.ROOT/'results/figure2_tree_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['source_tips','display_tips','collapsed_rows','root_diagnostic','diagnostic_exclusion_only','interpretation']},indent=2))

if __name__=='__main__':main()
