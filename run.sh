#!/bin/bash
set -e

########## GCN_old ##########
DATASET="OGBG-MOLTOX21" # OGBG-MOLTOX21, OGBG-MOLHIV, OGBN-ARXIV, OGBN-PROTEINS, OGBG-MOLTOX21-01
OUTPATH="./result/GCN_old/${DATASET}"
mkdir -p ${OUTPATH}
pypy3 gcn_old.py --ds_fn ${DATASET} | tee ${OUTPATH}/result.log

########## GCN ##########
# DATASET="OGBG-MOLTOX21" # OGBG-MOLTOX21, OGBG-MOLHIV, OGBN-ARXIV, OGBN-PROTEINS, OGBG-MOLTOX21-01
# OUTPATH="./result/GCN/${DATASET}"
# mkdir -p ${OUTPATH}
# pypy3 gcn.py --ds_fn ${DATASET} | tee ${OUTPATH}/result.log

########## GraphSage ##########
# DATASET="OGBG-MOLTOX21" # OGBG-MOLTOX21, OGBG-MOLHIV, OGBN-ARXIV, OGBN-PROTEINS, OGBG-MOLTOX21-01
# OUTPATH="./result/GraphSage/${DATASET}"
# mkdir -p ${OUTPATH}
# pypy3 graphsage.py --ds_fn ${DATASET} | tee ${OUTPATH}/result.log

######### GIN ##########
# DATASET="OGBG-MOLTOX21" # OGBG-MOLTOX21, OGBG-MOLHIV, OGBN-ARXIV, OGBN-PROTEINS, OGBG-MOLTOX21-01
# OUTPATH="./result/GIN/${DATASET}"
# mkdir -p ${OUTPATH}
# pypy3 gin.py --ds_fn ${DATASET} | tee ${OUTPATH}/result.log

######### GAT ##########
# DATASET="OGBG-MOLTOX21" # OGBG-MOLTOX21, OGBG-MOLHIV, OGBN-ARXIV, OGBN-PROTEINS, OGBG-MOLTOX21-01
# OUTPATH="./result/GAT/${DATASET}"
# mkdir -p ${OUTPATH}
# pypy3 gat.py --ds_fn ${DATASET} | tee ${OUTPATH}/result.log

########## MoNet ##########
# DATASET="OGBG-MOLTOX21" # OGBG-MOLTOX21, OGBG-MOLHIV, OGBN-ARXIV, OGBN-PROTEINS, OGBG-MOLTOX21-01
# OUTPATH="./result/MoNet/${DATASET}"
# mkdir -p ${OUTPATH}
# pypy3 monet.py --ds_fn ${DATASET} | tee ${OUTPATH}/result.log

########## GatedGCN ##########
# DATASET="OGBG-MOLTOX21" # OGBG-MOLTOX21, OGBG-MOLHIV, OGBN-ARXIV, OGBN-PROTEINS, OGBG-MOLTOX21-01
# OUTPATH="./result/GatedGCN/${DATASET}"
# mkdir -p ${OUTPATH}
# pypy3 gatedgcn.py --ds_fn ${DATASET} | tee ${OUTPATH}/result.log
# # pypy3 gatedgcn_16.py --ds_fn ${DATASET} | tee ${OUTPATH}/result.log

