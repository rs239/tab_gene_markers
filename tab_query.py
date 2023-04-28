#!/usr/bin/env python
"""
Author: Rohit Singh (rohitsingh@gmail.com)
Copyright (c) 2022 
MIT License


# The script takes as input a Fly Cell Atlas cell-type cluster, along with the fly organ the cluster belongs to 
#  and returns two gene markers for the cluster. It returns upto 3 such marker sets. 
# The code is implemented currently only for Drosophila melanogaster (species: "dm")
# 
# Separately, there are also private functions to process raw bulk RNAseq data and prepare it for downstream analysis. 
# It also prepares data obtained from a separate FCA dataset.
"""

import pandas as pd
import numpy as np
import scipy, os, sys, string, fileinput, glob, re, math, itertools, functools, csv, json, traceback, gc, argparse
import copy, multiprocessing, traceback, logging, pickle, time
import scipy.stats
from scipy.stats import describe
from scipy import sparse
import os.path
import scipy.sparse
from scipy.sparse import csr_matrix, csc_matrix
from collections import defaultdict
import anndata
import scanpy as sc
import urllib.request

import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")

# write a function dprint that prints only when debug is set to True
debug = False #True
def dprint(*args, **kwargs):
    if debug:
        print(*args, **kwargs)  


# Define the available choices for 'species', 'organ', and 'cluster'
VALID_SPECIES = ['dm']

VALID_ORGANS = {}
VALID_ORGANS['dm'] = "Other,AccessoryGland,Brain,Crop,Eye,FatBody,Gut,Head,Heart,MalpighianTubule,Ovary,RectalPad,SalivaryGland,Spermathecum,Testis,ThoracicoAbdominalGanglion".split(",")

VALID_CLUSTERS = {}
VALID_CLUSTERS['dm'] = "unannotated|16-cell germline cyst in germarium region 2a and 2b|adult abdominal pericardial cell|adult alary muscle|adult antenna glial cell|adult brain cell body glial cell|adult brain perineurial glial cell|adult differentiating enterocyte|adult enterocyte|adult esophagus|adult fat body|adult glial cell|adult heart|adult heart ventral longitudinal muscle|adult hindgut|adult hindgut*|adult lamina epithelial/marginal glial cell|adult Malpighian tubule|adult Malpighian tubule bar-shaped cell of initial segment|adult Malpighian tubule principal cell|adult Malpighian tubule principal cell of initial segment|adult Malpighian tubule principal cell of lower segment|adult Malpighian tubule principal cell of lower ureter|adult Malpighian tubule stellate cell of main segment|adult midgut*|adult midgut enterocyte|adult midgut-hindgut hybrid zone|adult neuron|adult oenocyte|adult olfactory receptor neuron acid-sensing, Ir64a|adult olfactory receptor neuron acid-sensing, Ir75a/b/c, Ir64a|adult olfactory receptor neuron Gr21a/63a|adult olfactory receptor neuron Ir56a+, Orco-|adult olfactory receptor neuron Ir75d|adult olfactory receptor neuron Ir84a, Ir31a, Ir76a, Ir76b, Ir8a, Or35a|adult olfactory receptor neuron Or13a|adult olfactory receptor neuron Or22a, Or42b, Or59b|adult olfactory receptor neuron Or47a, Or56a and likely other ORN types|adult olfactory receptor neuron Or47b|adult olfactory receptor neuron Or65|adult olfactory receptor neuron Or67a and likely other unknown ORN types|adult olfactory receptor neuron Or67d|adult olfactory receptor neuron Or85a, Or43b|adult olfactory receptor neuron Or88a|adult olfactory receptor neuron Or92a|adult olfactory receptor neuron unknown type, Orco-|adult optic chiasma glial cell|adult ostium|adult peripheral nervous system|adult peripheral neuron of the heart|adult pylorus|adult renal stem cell|adult reticular neuropil associated glial cell|adult salivary gland|adult tracheal cell|adult ventral nervous system|alpha'/beta' Kenyon cell|alpha/beta Kenyon cell|antennal lobe projection neuron|antennal trichoid sensillum at4|anterior ejaculatory duct|antimicrobial peptide-producing cell|arista and sacculus thermosensory and hygrosensory neuron Ir21a, Ir40a, Gr28b|artefact|auditory sensory neuron|bitter-sensitive labellar taste bristle|btl-GAL4 positive female cell, cluster 1, likely to be ovary cell|btl-GAL4 positive female cell, cluster 2, likely to be ovary cell|btl-GAL4 positive female cell, cluster 3, likely to be ovary cell|btl-GAL4 positive female cell, likely to be ovary cell, sim+|btl-GAL4 positive female cell, likely to be ovary cell, sim+, H15+|cardia (1)|cardia (2)|cardiomyocyte, working adult heart muscle (non-ostia)|cell body glial cell|central main body follicle cell ca. St. 6-8|centrifugal neuron C2|centrifugal neuron C3|choriogenic main body follicle cell and corpus luteum|choriogenic main body follicle cell St. 12|choriogenic main body follicle cell St. 14|CNS surface associated glial cell|columnar neuron T1|cone cell|copper cell|crop|crystal cell|cyst cell branch a|cyst cell branch b|cyst cell intermediate|cyst stem cell|distal medullary amacrine neuron Dm10|distal medullary amacrine neuron Dm11|distal medullary amacrine neuron Dm12|distal medullary amacrine neuron Dm3|distal medullary amacrine neuron Dm8|distal medullary amacrine neuron Dm9|dopaminergic neuron|dopaminergic PAM neuron|dorsal appendage forming follicle cell|dorsal rim area|early cyst cell 1|early cyst cell 2|early elongation stage spermatid|early-mid elongation-stage spermatid|ejaculatory bulb|ejaculatory bulb epithelium|ensheathing glial cell|enteroblast|enterocyte-like|enterocyte of anterior adult midgut epithelium|enterocyte of posterior adult midgut epithelium|enteroendocrine cell|eo support cell|epidermal cell of the abdominal posterior compartment|epidermal cell that specialized in antimicrobial response|epithelial cell|escort cell|eye photoreceptor cell|female reproductive system|follicle cell|follicle cell St. 9+|follicle stem cell and prefollicle cell|gamma Kenyon cell|germ cell stage 4 and later|germline cell|germline cell, unknown stage|gustatory receptor neuron|gustatory receptor neuron of the labellum|head cyst cell|hemocyte|indirect flight muscle|intestinal stem cell|Johnston organ neuron|Kenyon cell|labral sense organ mechanosensory neuron|lamina intrinsic amacrine neuron Lai|lamina monopolar neuron L1|lamina monopolar neuron L2|lamina monopolar neuron L3|lamina monopolar neuron L4|lamina monopolar neuron L5|lamina wide-field 1 neuron|lamina wide-field 2 neuron|late cyst cell branch a|late cyst cell branch b|late primary spermatocyte|leg muscle motor neuron|leg taste bristle chemosensory neuron|lobula columnar neuron LC10|lobula columnar neuron LC12|lobula columnar neuron LC17|main body follicle cell ca. until St. 5|male accessory gland|male accessory gland main cell|male accessory gland secondary cell|male germline differentiating cell|male gonad associated epithelium|male reproductive tract muscle|maxillary palp olfactory receptor neuron|maxillary palpus|mechanosensory neuron|mechanosensory neuron of haltere|mechanosensory neuron of leg chordotonal organ|medullary intrinsic neuron Mi1|medullary intrinsic neuron Mi15|medullary intrinsic neuron Mi4|medullary intrinsic neuron Mi9|midgut|midgut large flat cell|mid-late elongation-stage spermatid|mid-late proliferating spermatogonia|multidendritic neuron|muscle cell|neuron of haltere|nociceptive neuron|ocellus retinula cell|octopaminergic/tyraminergic neuron|olfactory receptor neuron|olfactory receptor neuron, coeloconics|optic-lobe-associated cortex glial cell|outer photoreceptor cell|ovarian sheath muscle|ovary cell|oviduct|pericerebral adult fat mass|perineurial glial sheath|peripheral glial cell|pheromone-sensing neuron|photoreceptor|photoreceptor cell R7|photoreceptor cell R8|photoreceptor-like|pigment cell|polar follicle cell|posterior midgut*|posterior terminal follicle cell ca. St. 5-8|post-mitotic endocycling nurse cell|post-mitotic germ cell early 16-cell cyst|Poxn neuron|prefollicle cell/stalk follicle cell|principal cell*|proximal medullary amacrine neuron Pm2|proximal medullary amacrine neuron Pm4|sacculus/arista neuron|scolopidial neuron|secretory cell of the male reproductive tract|seminal vesicle|sensory neuron|skeletal muscle of head|spermatid|spermatocyte|spermatocyte 0|spermatocyte 1|spermatocyte 2|spermatocyte 3|spermatocyte 4|spermatocyte 5|spermatocyte 6|spermatocyte 7a|spermatocyte cyst cell branch a|spermatocyte cyst cell branch b|spermatogonium|spermatogonium-spermatocyte transition|stalk follicle cell|stretch follicle cell|subperineurial glial cell|tendon cell|testis|T neuron T2|T neuron T2a|T neuron T3|T neuron T4/T5|T neuron T4/T5a-b|T neuron T4/T5c-d|tormogen cell|transmedullary neuron Tm1|transmedullary neuron Tm2|transmedullary neuron Tm20|transmedullary neuron Tm29|transmedullary neuron Tm3a|transmedullary neuron Tm4|transmedullary neuron Tm5c|transmedullary neuron Tm9|transmedullary Y neuron TmY14|transmedullary Y neuron TmY4|transmedullary Y neuron TmY5a|transmedullary Y neuron TmY8|visceral muscle|visceral muscle of the crop|visceral muscle of the midgut|young germ cell".split("|")

DEFAULT_PARAMS = {
    'prep_data_dir': '/afs/csail.mit.edu/u/r/rsingh/work/perrimonlab-marker-gene/data/',
    'datadir':None,
    'bulk_rnaseq_tau_thresh': 0.5,
    'scrnaseq_inv_dispersion_quantile': 0.25,
    'scrnaseq_largeset_diffexp_min_final_size': 20,
    'scrnaseq_largeset_diffexp_max_final_size': 200,
    'hits_sort_style': 'comb_diffexp_wilcoxon_mean',
    'num_hits': 3,
}



def __process_bulk_data(args):
    prep_data_dir = params_dict['prep_data_dir']
    df_bulk_meta = pd.read_table(f"{prep_data_dir}/raw/Metadata_bulk_RNAseq_FB202201_for_Rohit.txt", delimiter='\t')
    df_bulk_seq = pd.read_table(f"{prep_data_dir}/raw/bulk_RNASeq_matrix_FB202201_for_Rohit.txt", delimiter='\t')

    # ignore 1d adults
    idx = (df_bulk_meta.Stage == 'adult') & (df_bulk_meta['Stage_more'].apply(lambda s: int(s[:-1])) > 1)
    df1 = df_bulk_meta[idx].reset_index(drop=True)
    
    t2dset = defaultdict(list)
    for ds,t in zip(df1["DatasetID"], df1["Tissue"]):
        t2dset[t].append(ds)


    dfx = df_bulk_seq.iloc[:,:2]
    
    df_bulk_seq.iloc[:,2:] *= 1e6/df_bulk_seq.iloc[:,2:].sum()
    
    dd = {}
    for t,dset in t2dset.items():
        df2b = df_bulk_seq.loc[:,[c for c in df_bulk_seq.columns if c in dset]]
        dd[t] = np.ravel(df2b.mean(axis=1).values)
    df3 = pd.DataFrame(dd, index = df_bulk_seq["Gene"])
    df4 = np.log1p(df3)

    v_maxval = np.ravel(df4.max(axis=1).values)
    v_maxval_tissue = np.ravel(df4.idxmax(axis=1).values)
    v_maxval_tstat = np.ravel( ((df4.max(axis=1) -  df4.mean(axis=1))*np.sqrt(df4.shape[1])/(df4.std(axis=1))).values)

    # tau is a measure of specificity from https://doi.org/10.1093/bib/bbw008
    v_tau = np.ravel( ((1 - (df4/df4.max(axis=1).values[:,None])).sum(axis=1)/(df4.shape[1]-1)).values )

    dfx["tau"] = v_tau
    dfx["max_exp_logcpm"] = v_maxval
    dfx["max_exp_tissue"] = v_maxval_tissue
    dfx["max_exp_tstat"] = v_maxval_tstat

    dfx.to_csv(f"{prep_data_dir}/processed/bulk_data_gene_summary.csv", index=False)
    print(dfx)



def __prep_fca_data(params_dict):
    cells, genes, v = [],[],[]
    prep_data_dir = params_dict['prep_data_dir']
    for i,row in enumerate(csv.reader(open(f'{prep_data_dir}/raw/FCA_full_matrix.txt','rt'), delimiter=' ')):
        if i==0:
            cells=row
        else:
            genes.append(row[0])
            v.append(np.array([float(a) for a in row[1:]]))
        if i%100==0: print(i)
        
    m = csr_matrix(np.array(v))
    outpfx = f"{prep_data_dir}/processed/FCA_"
    dprint("Flag 298.30 ")
    scipy.sparse.save_npz(f"{outpfx}transcript_counts.npz", m)
    dprint("Flag 298.40 ")
    fh1 = open(f"{outpfx}genes.txt","wt")
    fh1.write("\n".join(genes))
    fh2 = open(f"{outpfx}cells.txt","wt")
    fh2.write("\n".join(hdr))

    dprint("Flag 298.50 ", m.shape)
    m2 = scipy.sparse.csr_matrix(m.T)
    dprint("Flag 298.52 ", m2.shape)

    df_meta = pd.read_table(f"{prep_data_dir}/raw/FCA_metadata_sample.tsv")
    df_meta.columns = ['cells','cell_type']
    df_meta.set_index('cells')

    dprint("Flag 298.55 ", df_meta.shape)
    adata = anndata.AnnData(X=m2)
    adata.obs_names = cells
    adata.var_names = genes
    adata.obs = df_meta
    
    dprint("Flag 298.60 ", adata.shape)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    dprint("Flag 298.65 ", adata.shape)
    sc.tl.rank_genes_groups(adata, 'cell_type', n_genes=500)

    dprint("Flag 298.70 ", adata.shape)
    adata.write(f'{outpfx}adata.h5')
    
    
def __prep_processed_data():
    """
    Internal function meant to process raw data (bulk as well as single-cell) into adata objects and 
    related entities for use. Not intended for broad usage
    """

    params_dict = copy.deepcopy(DEFAULT_PARAMS)
    __prep_processed_data(params_dict)
    __prep_fca_data(params_dict)


def download_data(params_dict):
    """
    Downloads data. If the file was already downloaded before, use the cached copy.
    """

    datadir = 'data' if params_dict['datadir'] is None else params_dict['datadir']

    for f in ['bulk_data_gene_summary.csv', 'FCA_adata.h5.gz']:
        url = 'http://cb.csail.mit.edu/cb/tab-data/' + f

        f_unzipped = f[:-3] if f.endswith('.gz') else f
        local_path = os.path.join(datadir, f_unzipped)
        if os.path.exists(local_path):
            print(f'Using cached copy of {local_path}')
        else:   
            os.makedirs(datadir, exist_ok=True)
            print(f'Downloading data from {url}...')
            dl_path = os.path.join(datadir, f)
            urllib.request.urlretrieve(url, dl_path)
            if f.endswith('.gz'):
                os.system('gunzip ' + dl_path)
            print(f'File details:\n\t')
            os.system('ls -l ' + local_path)


def __get_pairexp_stats(xdf, query_celltype):
    l = []
    xdf1 = xdf[xdf["cell_type"] == query_celltype].reset_index(drop=True)
    xdf2 = xdf[xdf["cell_type"] != query_celltype].reset_index(drop=True)
    for c in ["mean","min"]:
        v_base = xdf1[f"g_{c}"].values
        b = xdf2.groupby('cell_type')[f"g_{c}"].apply(lambda v: scipy.stats.ranksums(v_base, v).statistic).to_frame().reset_index()
        b = b.sort_values(f"g_{c}").reset_index(drop=True)
        l.append( (c, b["cell_type"].iloc[0], b[f"g_{c}"].iloc[0]))
        l.append( (c, b["cell_type"].iloc[1], b[f"g_{c}"].iloc[1]))

    return l



def query_one_cluster(species, organ, cluster, params_dict):
    """
    Query for the two genes for a cluster from the specified species and organ.
    returns a dataframe with the following columns: cluster, g1, g2, ... 
    The first row is the strongest hit, the second row is the second strongest hit, etc.
    The additional columns contain scores for the potentially conflicting clusters 
    """

    datadir = 'data' if params_dict['datadir'] is None else params_dict['datadir']
    adata = sc.read_h5ad(f"{datadir}/FCA_adata.h5")
    dprint("Flag 295.02 ")
    adata_X_t = scipy.sparse.csr_matrix(adata.X.T)
    dprint("Flag 295.03 ")
    
    genes = adata.var_names.tolist()
    cells = adata.obs_names.tolist()
    
    dprint("Flag 295.25 ", adata.shape)
    
    df_bulk = pd.read_csv(f"{datadir}/bulk_data_gene_summary.csv")

    bulk_rnaseq_tau_thresh = float(params_dict.get("bulk_rnaseq_tau_thresh", 0.5))
    bulk_rnaseq_tissue = organ
    query_celltype = cluster

    if bulk_rnaseq_tissue != "Other":
        gL1 = []
        for btissue in bulk_rnaseq_tissue.split(","):
            tL = (df_bulk['Gene'][ (df_bulk.tau > bulk_rnaseq_tau_thresh) &
                                   (df_bulk.max_exp_tissue == btissue) ]).values
            gL1 += list(tL)
        gL1 = list(set(gL1))
        
    else:
        bt_all_genes = set(df_bulk["Gene"].tolist())
        bt_regular_tissue_genes = set()
        for bt in df_bulk["max_exp_tissue"].unique():
            gL2 = (df_bulk['Gene'][ (df_bulk.tau > bulk_rnaseq_tau_thresh) &
                                    (df_bulk.max_exp_tissue == bt) ]).values
            bt_regular_tissue_genes.update(set(gL2))

        dprint("Flag 295.27 ", len(bt_all_genes), len(df_bulk["max_exp_tissue"].unique()), len(bt_regular_tissue_genes))
        gL1 = list( bt_all_genes - bt_regular_tissue_genes)
            

    adata2b = adata[ adata.obs['cell_type'] == query_celltype ].copy()
    adata2b = adata2b[:, adata2b.var_names.isin(gL1)]

    dprint(f"Flag 295.29 for celltype {query_celltype} in bulk tissue {bulk_rnaseq_tissue} matrix {adata2b.X.shape}")
    
    adata2 = anndata.AnnData(X = adata2b.X.todense(), obs=adata2b.obs, var=adata2b.var)
    
    adata3b = adata[:, adata.var_names.isin(gL1)].copy()
    adata3 = anndata.AnnData(X = adata3b.X.todense(), obs=adata3b.obs, var=adata3b.var)
    
    dprint("Flag 295.30 ", adata.shape, adata2.shape, adata3.shape, len(gL1), gL1[:5], adata2.var_names)


    adata2.var['inv_dispersion'] = np.ravel(np.mean(adata2.X,axis=0)/np.var(adata2.X,axis=0))
    inv_disp_thresh = adata2.var['inv_dispersion'].quantile( float(params_dict.get("scrnaseq_inv_dispersion_quantile",0.25)))

    gL2  = list(adata2.var_names[(adata2.var['inv_dispersion'] > inv_disp_thresh)])
    dprint("Flag 295.35 ", len(gL1), len(gL2))
    
    adata3b = adata3[:,adata3.var_names.isin(gL2)]

    xctlist = adata.uns['rank_genes_groups']["names"].dtype.names
    dprint("Flag 295.60 ", adata.shape, len(xctlist))

    K1 = int(params_dict.get("scrnaseq_largeset_diffexp_min_final_size",20))
    K2 = int(params_dict.get("scrnaseq_largeset_diffexp_max_initial_size",200))
    
    xdf = pd.DataFrame({ "cell_type": adata.obs['cell_type'].values})
    xdf["g_mean"] = 0.0
    xdf["g_min"] = 0.0

    
    xd = defaultdict(list) #map from celltypes to list of diffexp genes
    for j,l in enumerate(adata.uns['rank_genes_groups']["names"]):
        for k,b in enumerate(l):
            xd[xctlist[k]].append(b)

    #dprint("Flag 295.62 ", xd[query_celltype].index('esg'), len(xd[query_celltype])) #, xd[query_celltype])
    dprint("Flag 295.63 ", len(gL2), xd.keys(), len(xd), [len(a) for a in xd.values()] )#, [a for a in xd[query_celltype] if a in gL2])
    dprint("Flag 295.635 ", len(adata2.var_names)) #, adata2.var_names.tolist())

    #making the search for diffexp names a bit more generous
    for k in range(10, min(2*K2, len(xd[query_celltype]))):
        #we skip lncRNA and asRNA genes
        xw = [a for a in xd[query_celltype][:k] if a in gL2 and "RNA:" not in a] #top diffexp genes for the querycell type
        if len(xw)>=K1: break
        if len(xw)>=0.5*K1 and k>K2: break #if we have blown past K2, exit when  more than half of K1 is done

        
    dprint("Flag 295.64 ", xw, k, K1)
    
    xd2 = defaultdict(set) #map of diffexp gene from xw to other celltypes it is marked as diffexp in
    for c in xw:
        for u,uL in xd.items():
            if u==query_celltype: continue
            if c in uL: xd2[c].add(u)
    xd3 = []
    for i in range(len(xw)):
        ii = adata.var_names.tolist().index(xw[i])
        vi = np.ravel(adata_X_t[ii,:].todense())
        for j in range(i+1, len(xw)):
            jj = adata.var_names.tolist().index(xw[j])
            vj = np.ravel(adata_X_t[jj,:].todense())
            
            s = (xd2[xw[i]] & xd2[xw[j]]) #count how many other celltypes have both these genes as markers
            xdf["g_mean"] = (vi+vj)/2
            xdf["g_min"] = np.minimum(vi, vj)
            dprint("Flag 295.80 ", i, j, xw[i], xw[j], ii, jj, xdf.shape)
            
            s2 = __get_pairexp_stats(xdf, query_celltype) 
            xd3.append((xw[i], xw[j], len(s), s, *s2))

    dprint("Flag 295.90 ", len(xd3), xd3)
    sort_style = params_dict.get("hits_sort_style", "comb_diffexp_wilcoxon_mean")
    assert sort_style in ["diffexp_other", "wilcoxon_pair_mean", "wilcoxon_pair_min", "avg_wilcoxon_pair_mean_min","comb_diffexp_wilcoxon_mean"]

    f_sortkey = lambda s: s
    if sort_style == "diffexp_other":
        f_sortkey = lambda s: s[2]
    elif sort_style == "wilcoxon_pair_mean":
        f_sortkey = lambda s: -s[4][2]
    elif sort_style == "wilcoxon_pair_min":
        f_sortkey = lambda s: -s[6][2]
    elif sort_style == "avg_wilcoxon_pair_mean_min":
        f_sortkey = lambda s: -0.5*(s[4][2] + s[6][2])
    elif sort_style == "comb_diffexp_wilcoxon_mean":
        f_sortkey = lambda s: -s[4][2]  + 3*s[2]
        
    xd3b = sorted(xd3, key=f_sortkey)[:int(params_dict.get("num_hits", 3))]
    L = []
    for i in range(len(xd3b)):
        g1,g2,overlap1,s1,mean1,_,min1,_ = xd3b[i]
        L.append([query_celltype, g1,g2,overlap1, ";".join(s1), mean1[2], mean1[1], min1[2], min1[1], json.dumps(params_dict)])
    ydf = pd.DataFrame(L, columns="cluster,g1,g2,score_other_clusters_diffexp,list_other_clusters_diffexp,score_other_clusters_wilcoxon_mean,list_other_clusters_wilcoxon_mean,score_other_clusters_wilcoxon_min,list_other_clusters_wilcoxon_min,settings".split(","))
    return ydf


def main():
    # Parse the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--species', choices=VALID_SPECIES, help='the species to query', default='dm')
    parser.add_argument('--describe', action='store_true', help='describe the valid organ and cluster choices')
    parser.add_argument('--organ', help='the organ to query')
    parser.add_argument('--cluster', help='the cluster to query')
    parser.add_argument('--params', nargs='+', help='a list of key:value pairs of parameters')
    args = parser.parse_args()

    if args.organ and args.organ not in VALID_ORGANS[args.species]:
        raise ValueError(f'Invalid organ choice: {args.organ}. Use --describe to see available choices.')
    if args.cluster and args.cluster not in VALID_CLUSTERS[args.species]:     
        raise ValueError(f'Invalid cluster choice: {args.cluster}. Use --describe to see available choices.')
    

    params_dict = copy.deepcopy(DEFAULT_PARAMS)
    if args.params:
        # Create a dictionary from the key-value pairs in the '--params' argument
        params_dict.update( dict([tuple(w.split(':')[:2]) for w in args.params]) )
        print(f'Params: {params_dict}')

    if args.describe:
        # Return the valid organ and cluster choices for the specified species
        s_organs = "\n\t".join([""]+list(VALID_ORGANS[args.species]))
        s_clusters = "\n\t".join([""]+list(VALID_CLUSTERS[args.species]))
        print(f'Valid organ choices for species {args.species}: {s_organs}')
        print()
        print(f'Valid cluster choices for species {args.species}: {s_clusters}')
    else:
        if not args.organ or not args.cluster:
            raise ValueError('Must specify either --describe or both --organ and --cluster')
            
        # Download the data for the specified species
        download_data(params_dict)

        ydf = query_one_cluster(args.species, args.organ, args.cluster, params_dict)
        if ydf is None or ydf.shape[0]==0:
            print("""No results found, unfortunately. Try changing the parameters. 
The default param values are   
    bulk_rnaseq_tau_thresh:0.5  
    scrnaseq_inv_dispersion_quantile:0.25
Try decreasing both, possibly all the way to zero.""")
    
        else:
            ydf.to_csv(sys.stdout, index=False, float_format='%g')



if __name__ == '__main__':
    main()
