# Written by Laura Markey, Dec 11 2023
# Edited by Veda Khadka
# Lieberman Lab, MIT IMES
# Part of Khadka and Markey et al, 2023. "Commensal skin bacteria exacerbate inflammation and delay healing"


# python input: 
# -sample metadata file (provided)
# -counts table of normalized TPM per gene per sample from kallisto and sleuth (part1)
# -list of genes annotated as immune by Jax (external, provided)
# -immune enrichment and not-immune enrichment output from clusterProfiler (part1)

# python output:
# -volcano plots showing immune genes as pdfs
# -filtered by gene count and fold change GO process enrichment csvs
# -heatmaps for GO processes and genes of interest


#%% import packages
import os
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import zscore
from scipy import cluster

#%%set up directory and read in files from R
os.chdir('/Users/liebermanlab/Dropbox (MIT)/Lieberman Lab/Projects/Commensal_memory_mouse/data_and_code')
counts=pd.read_csv("data/transcriptomics/allsamples_genecounts.csv")
meta=pd.read_csv("data/transcriptomics/longmeta.csv")
#differential expression analysis performed in R
d0lrt=pd.read_csv("data/transcriptomics/d0lrtresults.csv").dropna(how="any")
d1lrt=pd.read_csv("data/transcriptomics/d1lrtresults.csv").dropna(how="any")
d4lrt=pd.read_csv("data/transcriptomics/d4lrtresults.csv").dropna(how="any")
#go enrichment analysis performed in R
immune_enrichment=pd.read_csv("data/transcriptomics/full_GO_enrichment_upregulated_immunegenes_bytimepoint.csv")
notimmune_enrichment=pd.read_csv("data/transcriptomics/full_GO_enrichment_nonimmune_upregulated_bytimepoint.csv")
#all immune genes needed as a background list for coloring and analysis
jaximmune=pd.read_csv("data/transcriptomics/jax_immune_list.csv")
jax_immune_genes=jaximmune["Symbol"]
#%% Calculate fold change and add to df containing gene,adjusted pvalue and immune category for volcano plot
#also add immune category
#manually calculated fold change from counts
d0countsm=counts[counts["endpoint"]==0].reset_index(drop=True)
d1countsm=counts[counts["endpoint"]==1].reset_index(drop=True)
latecountsm=counts[counts["endpoint"]==4].reset_index(drop=True)
d0pbs=d0countsm[d0countsm["group"]=="pbs"].iloc[:,23:-1].mean(axis=0)
d0sepi=d0countsm[d0countsm["group"]=="sepi44"].iloc[:,23:-1].mean(axis=0)
d0fc=np.log2(d0sepi.div(d0pbs))
d1pbs=d1countsm[d1countsm["group"]=="pbs"].iloc[:,23:-1].mean(axis=0)
d1sepi=d1countsm[d1countsm["group"]=="sepi44"].iloc[:,23:-1].mean(axis=0)
d1fc=np.log2(d1sepi.div(d1pbs))
d4pbs=latecountsm[latecountsm["group"]=="pbs"].iloc[:,23:-1].mean(axis=0)
d4sepi=latecountsm[latecountsm["group"]=="sepi44"].iloc[:,23:-1].mean(axis=0)
d4fc=np.log2(d4sepi.div(d4pbs))
d0pvalonly=d0lrt[["target_id", "qval", "ext_gene"]]
d0pvalonly.index=d0pvalonly["target_id"]
d0summ=pd.concat([d0fc,d0pvalonly],axis=1)
d0summ["log2fc"]=d0summ[0]
d0summ["immuneornot"]=np.where(d0summ["ext_gene"].isin(jax_immune_genes),"immune", "other")
d1pvalonly=d1lrt[["target_id", "qval", "ext_gene"]]
d1pvalonly.index=d1pvalonly["target_id"]
d1summ=pd.concat([d1fc,d1pvalonly],axis=1)
d1summ["log2fc"]=d1summ[0]
d1summ["immuneornot"]=np.where(d1summ["ext_gene"].isin(jax_immune_genes),"immune", "other")
d4pvalonly=d4lrt[["target_id", "qval", "ext_gene"]]
d4pvalonly.index=d4pvalonly["target_id"]
d4summ=pd.concat([d4fc,d4pvalonly],axis=1)
d4summ["log2fc"]=d4summ[0]
d4summ["immuneornot"]=np.where(d4summ["ext_gene"].isin(jax_immune_genes),"immune", "other")
#%%MAKE VOLCANO PLOTS HIGHLIGHTING IMMUNE GENES FOR SUPPLEMENT
#add color palette for coloring not significant genes based on immune or not
immunepal=dict(immune="red", other="gray")
plt.rcParams['figure.dpi'] = 300
#%%day0
#first plot not sig stuff
sns.scatterplot(data=d0summ[d0summ["qval"]>0.05], x="log2fc", y=-np.log(d0summ[d0summ["qval"]>0.05]["qval"]), hue="immuneornot", palette=immunepal, alpha=0.2,s=30)
#then plot sig genes not immune
sns.scatterplot(data=d0summ[(d0summ["qval"]<0.05)&(d0summ["immuneornot"]!="immune")], x="log2fc", y=-np.log(d0summ[(d0summ["qval"]<0.05)&(d0summ["immuneornot"]!="immune")]["qval"]), color= "gray", alpha=0.2, s=30)
#then plot sig immune genes
sns.scatterplot(data=d0summ[(d0summ["qval"]<0.05)&(d0summ["immuneornot"]=="immune")], x="log2fc", y=-np.log(d0summ[(d0summ["qval"]<0.05)&(d0summ["immuneornot"]=="immune")]["qval"]), color= "red", s=30, linewidth=0.5, edgecolor="black")
plt.legend("")
plt.hlines(-np.log(0.05), -5,5, color="black", linestyles="dashed", alpha=0.5,linewidth=0.5)
#add text to indicate points of interest
#include all d0 genes on d0 immune volcano
d0immunegenes=d0summ[(d0summ["qval"]<0.05)&(d0summ["immuneornot"]=="immune")]["ext_gene"]
for g in d0immunegenes:
    x=d0summ[d0summ["ext_gene"]==g]["log2fc"]
    y=-np.log(d0summ[d0summ["ext_gene"]==g]["qval"])
    plt.text(x,y,str(g),color="black", fontsize=12)
plt.title("Day 0 Immune Genes")
plt.ylabel("-log(qval)")
plt.xlim(-5,5)
plt.ylim(0,20)
plt.savefig("Figs/supplement/d0volcanoalltext.pdf", format="pdf",dpi=300)
#%%DAY 1
#first plot not sig stuff
sns.scatterplot(data=d1summ[d1summ["qval"]>0.05], x="log2fc", y=-np.log(d1summ[d1summ["qval"]>0.05]["qval"]), hue="immuneornot", palette=immunepal, alpha=0.2,s=30)
#then plot sig genes not immune
sns.scatterplot(data=d1summ[(d1summ["qval"]<0.05)&(d1summ["immuneornot"]!="immune")], x="log2fc", y=-np.log(d1summ[(d1summ["qval"]<0.05)&(d1summ["immuneornot"]!="immune")]["qval"]), color= "gray", alpha=0.2, s=30)
#then plot sig immune genes
sns.scatterplot(data=d1summ[(d1summ["qval"]<0.05)&(d1summ["immuneornot"]=="immune")], x="log2fc", y=-np.log(d1summ[(d1summ["qval"]<0.05)&(d1summ["immuneornot"]=="immune")]["qval"]), color= "red", s=30, linewidth=0.5, edgecolor="black")
plt.legend("")
plt.hlines(-np.log(0.05), -5,5, color="black", linestyles="dashed",alpha=0.5,linewidth=0.5)
#plt.vlines(-2, 0, 25, color="black", linestyles="dashed")
#plt.vlines(2, 0, 25, color="black", linestyles="dashed")
plt.ylim(0,25)
plt.xlim(-5,5)
plt.title("Day 1 Immune Genes")
plt.ylabel("-log(qval)")
plt.legend("")
plt.text(2.03,3.4, "Mmp3", color="black", fontweight="bold")
plt.text(2.9,4.5, "Il6", color="black", fontweight="bold")
plt.text(2.9,4.5, "Lcn2", color="black", fontweight="bold")
plt.text(1.4,4.5, "Il1r1", color="black", fontweight="bold")
plt.text(1.1,4.4, "Cxcl12", color="black", fontweight="bold")
plt.text(-3,3.2,"Dapl1", color="black", fontweight="bold")
#plt.show()
plt.savefig("Figs/supplement/d1volcanoalltext.pdf", format="pdf",dpi=300)
plt.show()
#%%day 4
#first plot not sig stuff
sns.scatterplot(data=d4summ[d4summ["qval"]>0.05], x="log2fc", y=-np.log(d4summ[d4summ["qval"]>0.05]["qval"]), hue="immuneornot", palette=immunepal, alpha=0.2,s=30)
#then plot sig genes not immune
sns.scatterplot(data=d4summ[(d4summ["qval"]<0.05)&(d4summ["immuneornot"]!="immune")], x="log2fc", y=-np.log(d4summ[(d4summ["qval"]<0.05)&(d4summ["immuneornot"]!="immune")]["qval"]), color= "gray", alpha=0.2, s=30)
#then plot sig immune genes
sns.scatterplot(data=d4summ[(d4summ["qval"]<0.05)&(d4summ["immuneornot"]=="immune")], x="log2fc", y=-np.log(d4summ[(d4summ["qval"]<0.05)&(d4summ["immuneornot"]=="immune")]["qval"]), color= "red", s=30, linewidth=0.5, edgecolor="black")
plt.legend("")
plt.hlines(-np.log(0.05), -7,7, color="black", linestyles="dashed",alpha=0.5,linewidth=0.5)
#plt.vlines(-2, 0, 20, color="black", linestyles="dashed")
#plt.vlines(2, 0, 20, color="black", linestyles="dashed")
plt.ylim(0,20)
plt.xlim(-7,7)
plt.title("Day 4 Immune Genes")
plt.ylabel("-log(qval)")
plt.text(6.3,3.9, "Lcn2", color="black", fontweight="bold")
plt.text(2.3,4.3, "Ifnar2", color="black", fontweight="bold")
plt.text(2.2,3.3, "Ifitm2", color="black", fontweight="bold")
plt.text(4.2,3.2, "Ccl7", color="black", fontweight="bold")
plt.text(-1,10.2,"Pecam1", color="black", fontweight="bold")
plt.text(1.85,10.2,"Apod", color="black", fontweight="bold")
plt.legend("")
plt.savefig("Figs/supplement/d4volcanoalltext.pdf", format="pdf",dpi=300)
plt.show()
#%%filtering immune gene clusterProfiler outputs and adding "pathway fold change" information for plotting 
#only stuff with more than 5 genes
go_bp_immune_filter=immune_enrichment[immune_enrichment["Count"]>5]
#want to filter by average gene fold change in pathway so need to calculate that
#for every gene in a GO process, calculate the average TPM in PBS and sepi and thus ratio; then average that ratio for all genes in GO process
#generate a df of just the fc for each timepoint for each pathway
go_bp_immune_filter["expandlist"]=go_bp_immune_filter["geneID"].str.split("/")
pathwaylist=[]
d1fc_g=[]
d0fc_g=[]
d4fc_g=[]
for p in go_bp_immune_filter["ID"].unique():
    for g in go_bp_immune_filter[go_bp_immune_filter["ID"]==p]["expandlist"]:
        pbsd0=d0countsm[d0countsm["group"]=="pbs"][g].mean()
        sepid0=d0countsm[d0countsm["group"]=="sepi44"][g].mean()
        d0fc=sepid0.div(pbsd0)
        pbsd1=d1countsm[d1countsm["group"]=="pbs"][g].mean()
        sepid1=d1countsm[d1countsm["group"]=="sepi44"][g].mean()
        d1fc=sepid1.div(pbsd1)
        pbsd4=latecountsm[latecountsm["group"]=="pbs"][g].mean()
        sepid4=latecountsm[latecountsm["group"]=="sepi44"][g].mean()
        d4fc=sepid4.div(pbsd4)
    pathwaylist.append(p)
    d1fc_g.append(d1fc.mean())
    d0fc_g.append(d0fc.mean())
    d4fc_g.append(d4fc.mean())
combined_fc=pd.DataFrame({"ID":pathwaylist,"d0fc":d0fc_g, "d1fc":d1fc_g, "d4fc":d4fc_g })
#save this df to manually add fold change to csv in order to plot later
combined_fc.to_csv("immune_enrich_GO_process_ave_foldchange.csv")
#%% continued immune gene enrichement filtering: using the fold change reference make a list of GO processes to show
#only include stuff with >2 fold change at at least one timepoint
filter2_gobp_immune=combined_fc[(combined_fc["d0fc"]>5)|(combined_fc["d1fc"]>5)|(combined_fc["d4fc"]>5)]["ID"]
#there are 9 pathways with this filter on d0 so include all of those in table
d0_go_bp_include=go_bp_immune_filter[(go_bp_immune_filter["Cluster"]=="d0")&(go_bp_immune_filter["ID"].isin(filter2_gobp_immune))]["ID"]
#there are 70 pathways on d1 so include the top 10
d1bpimmune=go_bp_immune_filter[(go_bp_immune_filter["Cluster"]=="d1")&(go_bp_immune_filter["ID"].isin(filter2_gobp_immune))].sort_values(by="p.adjust").iloc[:10,:]["ID"]
#there are 1367 pathwys on d4 so include the top 10
d4bpimmune=go_bp_immune_filter[(go_bp_immune_filter["Cluster"]=="d4")&(go_bp_immune_filter["ID"].isin(filter2_gobp_immune))].sort_values(by="p.adjust").iloc[:10,:]["ID"]
#combine these series into a list to filter the table
top_bp_immune=pd.concat([d0_go_bp_include, d1bpimmune, d4bpimmune])
#filter combined df by this and then use for plotting because pvalues already calculated
alltp_filteredgenes_bp=immune_enrichment[immune_enrichment["ID"].isin(top_bp_immune)]
alltp_filteredgenes_bp["expandlist"]=alltp_filteredgenes_bp["geneID"].str.split("/")
#this filtered list of enriched immune GO processes was then manually curated for processes of interest and fold-change added for visualization
alltp_filteredgenes_bp.to_csv("data/transcriptomics/filtered_withfoldchange_immune_enriched_GO_processes.csv")
#%%#only stuff with more than 5 genes
notimmune_filter=notimmune_enrichment[notimmune_enrichment["Count"]>5]
#add average gene fold change in pathway
notimmune_filter["expandlist"]=notimmune_filter["geneID"].str.split("/")
pathwaylist_notimmune=[]
d1fc_g_n=[]
d4fc_g_n=[]
for p in notimmune_filter["ID"].unique():
    for g in notimmune_filter[notimmune_filter["ID"]==p]["expandlist"]:
        pbsd1_n=d1countsm[d1countsm["group"]=="pbs"][g].mean()
        sepid1_n=d1countsm[d1countsm["group"]=="sepi44"][g].mean()
        d1fc_n=sepid1_n.div(pbsd1_n)
        pbsd4_n=latecountsm[latecountsm["group"]=="pbs"][g].mean()
        sepid4_n=latecountsm[latecountsm["group"]=="sepi44"][g].mean()
        d4fc_n=sepid4_n.div(pbsd4_n)
    pathwaylist_notimmune.append(p)
    d1fc_g_n.append(d1fc_n.mean())
    d4fc_g_n.append(d4fc_n.mean())
combined_fc_notimmune=pd.DataFrame({"ID":pathwaylist_notimmune, "d1fc_n":d1fc_g_n, "d4fc_n":d4fc_g_n })
notimmune_filter.to_csv("data/transcriptomics/genecount_only_filtered_notimmune.csv")
combined_fc_notimmune.to_csv("data/transcriptomics/notimmune_foldchange_forfigure.csv")
#%% Pulling out genes from processes of interest and averaging TPM per group per timepoint and zscoring for heatmaps
#first need to build a dataframe that has the average gene expression per group at each timepoint
#immune figure has IL-1 and neutrophil migration pathways
d0=d0countsm[(d0countsm['group']=="pbs")|(d0countsm['group']=="sepi44")]
d1=d1countsm[(d1countsm['group']=="pbs")|(d1countsm['group']=="sepi44")]
d4=latecountsm[(latecountsm["group"]=="sepi44")|(latecountsm["group"]=="pbs")]
il1=alltp_filteredgenes_bp[alltp_filteredgenes_bp["ID"]=="GO:0070555"]["expandlist"].explode()
neut=alltp_filteredgenes_bp[alltp_filteredgenes_bp["ID"]=="GO:1990266"]["expandlist"].explode()
immuneheatgenes=pd.concat([il1,neut]).unique()
#d0gene0
d0genes=d0[immuneheatgenes]
d0genes["group"]=d0["group"]
d0ave_0=d0genes.groupby("group").mean()
#d1 gene0
d1genes=d1[immuneheatgenes]
d1genes["group"]=d1["group"]
d1ave_0=d1genes.groupby("group").mean()
#d4 gene0
d4genes=d4[immuneheatgenes]
d4genes["group"]=d4["group"]
d4ave_0=d4genes.groupby("group").mean()
immuneheat=pd.concat([d0ave_0,d1ave_0,d4ave_0])
#making new dataframe to organize genes for heatmap
immune_heat_genes=d1lrt[d1lrt["target_id"].isin(immuneheatgenes)][["target_id","ext_gene"]]
immune_heat_genes["category"]=np.where(immune_heat_genes["target_id"].isin(il1), "il1", "neut")
sorted_immune=immune_heat_genes.sort_values(["category", "ext_gene"])
#reorder ave tpm data by the categories found in reference df
colsort_data=immuneheat[sorted_immune["target_id"]]
colsort_data["group"]=immuneheat.index
#save to csv to visualize in R
data_zscale=colsort_data.iloc[:,:-1].astype("float64").apply(zscore)
#make a list of the genenames that correspond to the ensembl ids used to make dataframe
genename=[]
for c in data_zscale.columns:
    newname=immune_heat_genes[immune_heat_genes["target_id"]==c]["ext_gene"].unique()
    mod=str(newname)[2:-2]
    genename.append(mod)
data_zscale.columns=genename
data_zscale["group"]=colsort_data["group"]
data_zscale["endpoint"]=[0,0,1,1,4,4]
lut=dict(zip(data_zscale["group"].unique(), ["#808080", "#0073C2FF", "red", "orange", "green", "white"]))
row_colors=data_zscale["group"].map(lut)
sns.clustermap(data_zscale.iloc[:,:-2],row_colors=row_colors,xticklabels=genename,col_cluster=False, row_cluster=False,yticklabels=False,linewidths=0.5, figsize=(40,10), linecolor="black", cmap='Purples')
plt.xticks(fontsize=14)
plt.savefig("Figs/Fig3g_heatmap_immune.svg", format="svg")
data_zscale.to_csv("data/transcriptomics/zscale_immune_heatmap.csv")
#%%extracellular matrix/proliferation pathway heatmap
#made a list of endopeptidase genes to avoid plotting very rare or otherwise uninteresting genes
endopeptid=pd.Series(['ENSMUSG00000030022','ENSMUSG00000057530','ENSMUSG00000050063','ENSMUSG00000044430','ENSMUSG00000044737','ENSMUSG00000024713','ENSMUSG00000045027','ENSMUSG00000050762','ENSMUSG00000000957','ENSMUSG00000025355'],index=range(0,10))
ecm=notimmune_filter[notimmune_filter["ID"]=="GO:0005201"]["expandlist"].explode()
prolifheatgenes=pd.concat([endopeptid,ecm]).unique()
#making new dataframe to organize genes for heatmap
prolif_heat_genes=d1lrt[d1lrt["target_id"].isin(prolifheatgenes)][["target_id","ext_gene"]]
prolif_heat_genes["category"]=np.where(prolif_heat_genes["target_id"].isin(endopeptid), "endopeptid", "extra")
prolifsort=prolif_heat_genes.sort_values(["category", "ext_gene"])
#making averaged dataframe
d0genes2=d0[prolifheatgenes]
d0genes2["group"]=d0["group"]
d0ave_2=d0genes2.groupby("group").mean()
d0ave_2["tp"]=0
d1genes2=d1[prolifheatgenes]
d1genes2["group"]=d1["group"]
d1ave_2=d1genes2.groupby("group").mean()
d1ave_2["tp"]=1
d4genes2=d4[prolifheatgenes]
d4genes2["group"]=d4["group"]
d4ave_2=d4genes2.groupby("group").mean()
d4ave_2["tp"]=4
#d4ave_2["path"]="loc2"
prolif_pathway_sum=pd.concat([d0ave_2,d1ave_2,d4ave_2])
#plotting proliferation heatmap
#only plotting day 1 and day 4 
data=prolif_pathway_sum.iloc[2:,:]
colsort_data=data[prolifsort["target_id"]]
colsort_data["group"]=data.index
data_zscale=colsort_data.iloc[:,:-1].astype("float64").apply(zscore)
genename2=[]
for c in data_zscale.columns:
    newname=d0lrt[d0lrt["target_id"]==c]["ext_gene"].unique()
    mod=str(newname)[2:-2]
    genename2.append(mod)
data_zscale.columns=genename2
data_zscale["group"]=colsort_data["group"]
data_zscale["endpoint"]=[1,1,4,4]
lut=dict(zip(data_zscale["group"].unique(), ["#808080", "#0073C2FF", "red", "orange", "green", "white"]))
row_colors=data_zscale["group"].map(lut)
sns.clustermap(data_zscale.iloc[:,:-2],row_colors=row_colors,xticklabels=genename2,col_cluster=False, row_cluster=False,yticklabels=False,linewidths=0.5, figsize=(40,10), linecolor="black", cmap='Greens')
plt.xticks(fontsize=14)
plt.savefig("Figs/Fig4b_heatmap_proliferation.svg", format="svg")
data_zscale.to_csv("data/transcriptomics/zscale_proliferation_heatmap.csv")