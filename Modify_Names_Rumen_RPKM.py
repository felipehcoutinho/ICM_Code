#! /usr/bin/env python3

from collections import defaultdict
import pandas as pd
import re

def infer_oid(nid):
    oid = nid
    oid = re.sub("Rumen_Virus_Seq_(\d)+(\s)","",oid)
    #print("Infer_OID","NID",nid,"OID:",oid)
    oid2nid[oid] = nid
    return oid
    
def infer_nid(oid):
    #print("Infer_NID","OID:",oid)
    #sub_info_df = info_df[info_df["Original_ID"] == oid]
    #nid = sub_info_df.index.values[0]
    #print(Infer_NID","NID:",nid)
    nid = oid2nid[oid]
    short_nid = nid.split(" ")[0]
    #print("Oid",oid,"NID",nid,"Short",short_nid)
    return short_nid
    #return "Dummy"

oid2nid = defaultdict()

info_df = pd.read_csv("/mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Seq_Info.tsv", sep="\t",index_col="Sequence",header=0)
info_df["Original_ID"] = info_df["Description"].apply(infer_oid)

abund_df = pd.read_csv("/mnt/lustre/scratch/fcoutinho/Rumen/Abundance/Viral_rpkm.tsv", sep="\t",index_col="Sequence",header=0)
abund_df["New_ID"] = abund_df.index.map(infer_nid)

abund_df.set_index("New_ID",drop=True,inplace=True,verify_integrity=True)
#abund_df.drop("New_ID",axis=1,inplace=True)

abund_df.to_csv("/mnt/lustre/scratch/fcoutinho/Rumen/Abundance/Renamed_Viral_rpkm.tsv",sep="\t",na_rep=0)