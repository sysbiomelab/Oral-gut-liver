#/usr/bin/python3
#The target ouputs are
#A) Reaction richness by Sample (A.1 limited to total number of reactions in global model
# or A.2 the total number of KO->Reactions found in genecount table, irrespective of MSPs)
#B) Reaction relative abundnace table, transform MSP (into genes into  KO) into Reactions
# values in output reflect the sum of relative abundances of genes falling into the Reaction annotation
#C) Reactobiome, presence absece, fraction of species in sample with a funcion in relation to the total
# number of species in the sample.

import argparse
import numpy as np
import pandas as pd
import argparse
import os


def main():
    #####################
    #####Reaction Relative Abundance and Reactobiome Table
    # Output Reactions use KEGG_reactions IDs
    script_dir = os.path.abspath( os.path.dirname( __file__ ) )
    parser=argparse.ArgumentParser()
    parser.add_argument("--msp_table",  metavar='FILE', help="CSV file msp Relative abundance table ouptut form Meteor", default="")
    parser.add_argument("--gct_table",  metavar='FILE', help="TSV file with IGC2 gene ids, and GeneID output by meteor", default="")

    parser.add_argument("--igcmspko", metavar='FILE', help="TSV file with IGC2 gene ids, associated MSP and associated KO",
                        default=os.path.join(script_dir, "data/IGC2_msp_ko.tsv"))
    parser.add_argument("--igc_id",  metavar='FILE', help="CSV file with IGC2 gene ids, and GeneID from gct meteor output",
                        default=os.path.join(script_dir, "data/igc2.meteorids.csv"))
    parser.add_argument("--ko2rn",  metavar='FILE', help="TSV  file with KO and associated RN",
                        default=os.path.join(script_dir, "data/ko2rn.tsv"))

    args=parser.parse_args()


    #Read data tables for conversion from IGC2 catalogue to Kegg Reactions
    igcmspko = pd.read_csv("/Users/jinyi/Documents/github_reactobiome/reactobiome/data/IGC2_msp_ko.tsv", sep="\t", index_col=0, header=None)
    igcmspko.columns = ["msp", "ko"]
    igcmeteorid = pd.read_csv("/Users/jinyi/Documents/github_reactobiome/reactobiome/data/igc2.meteorids.csv", header=None, index_col=0)
    igcmeteorid.columns = ["geneid"]
    ko2rn = pd.read_csv("/Users/jinyi/Documents/github_reactobiome/reactobiome/data/ko2rn.tsv", sep="\t", header=None)
    ko2rn.columns=["ko", "rn"]

    if len(args.msp_table):
        msp_table=pd.read_csv("rytypePJ_merge_gut_lc_and_h.csv", index_col=0)
        msp_table=msp_table*10e6
        msp_presence= msp_table>0
        msp_presence= msp_presence.astype(int)
        #convert MSP table into KO annotation
        ko_table=pd.merge(igcmspko, msp_table,left_on="msp", right_index=True , how="left", sort=False)
        ko_table=ko_table.dropna().groupby("ko").sum()
        ko_table.drop(columns="msp")
        ko_table.to_csv("ko_relab_frommsp.csv")
        
        rn_table=pd.merge(ko2rn, ko_table,left_on="ko", right_index=True , how="left", sort=False)
        rn_table=rn_table.dropna().groupby("rn").sum()
        rn_table.to_csv("reaction_relab_frommsp.csv")
        


        #Reactobiome from MSP table
        msprn = pd.merge( ko2rn,igcmspko, left_on="ko", right_on="ko", how="inner", sort=False ).dropna()
        msprn = msprn.drop(columns="ko").drop_duplicates()
        mspreactobiome = pd.merge(msprn, msp_presence, left_on="msp", right_index=True, how="left", sort=False).dropna()
        reactobiome = mspreactobiome.dropna().groupby("rn").sum()
        reactobiome = reactobiome.drop(columns="msp")
        reactobiome = (reactobiome.divide(msp_presence.sum(), axis=1)*500).astype(int)
        reactobiome.to_csv("reactobiome_frommsp.csv")
        
        mspkobiome = pd.merge(igcmspko, msp_presence, left_on="msp", right_index=True, how="left", sort=False).dropna()
        kobiome =  mspkobiome.dropna().groupby("ko").sum()
        kobiome = kobiome.drop(columns="msp")
        kobiome = (kobiome.divide(msp_presence.sum(), axis=1)*500).astype(int)
        kobiome.to_csv("kobiome_frommsp.csv")






    #Convert genecount table into KO count table (KO against samples, sum of relative ab values)
    #using IGC2 gene catalogue KO annotation (and MSP content for )
    if len(args.gct_table) >0:
        gct_table=pd.read_csv(args.gct_table, index_col=0, sep="\t")
        gct_table=gct_table.loc[(gct_table.sum(axis=1) != 0), :]
        gct_presence=(gct_table>0).astype(int)
        gctmspko=pd.merge(igcmeteorid, igcmspko, left_on="geneid", right_index=True, how="inner", sort=False)
        ko_table=gctmspko.join(gct_table, how="inner")
        ko_table=ko_table.dropna().groupby("ko").sum()
        rn_table=pd.merge(ko2rn, ko_table,left_on="ko", right_index=True , how="left", sort=False)
        rn_table=rn_table.dropna().groupby("rn").sum()
        rn_table.to_csv("reaction_relab.csv")
        ###### Reactobiome table ########
        # Count of reaction present per 500 hundred MSP,
        gctmsprn=pd.merge(igcmspko.reset_index(), ko2rn, left_on="ko", right_on="ko", how="inner", sort=False)
        gctmsprn=gctmsprn.set_index(0, drop=True)
        gctmsprn = pd.merge(igcmeteorid, gctmsprn, left_on="geneid", right_index=True, how="inner", sort=False)
        #gctmsprn = gctmsprn.drop(columns="ko").drop_duplicates(subset=["msp","rn"])

        gctreactobiome = gctmsprn.join(gct_presence, how="inner")
        gctreactobiome = (gctreactobiome.groupby(["msp","rn"]).sum() > 0).astype(int).reset_index()
        reactobiome = gctreactobiome.groupby("rn").sum()
        mspinsample = (gctreactobiome.groupby("msp").sum() >0 ).sum()
        reactobiome = (reactobiome.divide(mspinsample, axis=1)*500).astype(int)
        reactobiome.to_csv("reactobiome_table.csv")

    # From the genecount table sum number genes present with reaction,
    # divide by the total number of MSP present(get this number) from MSPtable trasnformed into MSP Presence table

if __name__== "__main__":
    main()
