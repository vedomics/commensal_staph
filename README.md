# Commensal skin bacteria exacerbate inflammation and delay skin healing
#### Veda D. Khadka*, Laura Markey*, Magalie Boucher, Tami D. Lieberman
'*' Indicates co-first authorship

The skin microbiome can both trigger beneficial immune stimulation and pose a potential infection threat. Previous studies have shown that colonization of mouse skin with the model human skin commensal Staphylococcus epidermidis is protective against subsequent excisional wound or pathogen challenge. However, less is known about concurrent skin damage and exposure to commensal microbes, despite growing interest in interventional probiotic therapy. Here, we address this open question by applying commensal skin bacteria at a high dose to abraded skin. While depletion of the skin microbiome via antibiotics delayed repair from damage, application of commensals-- including the mouse commensal Staphylococcus xylosus, three distinct isolates of S. epidermidis, and all other tested human skin commensals-- also significantly delayed barrier repair. Increased inflammation was observed within four hours of S. epidermidis exposure and persisted through day four, at which point the skin displayed a chronic-wound-like inflammatory state with increased neutrophil infiltration, increased fibroblast activity, and decreased monocyte differentiation. Transcriptomic analysis suggested that the prolonged upregulation of early canonical proliferative pathways inhibited the progression of barrier repair. These results highlight the nuanced role of members of the skin microbiome in modulating barrier integrity and indicate the need for caution in their development as probiotics.

## Included: 

Code and data archive for the publication, "Commensal skin bacteria exacerbate inflammation and delay skin healing" [(Khadka and Markey, Biorxiv, 2023)](https://www.biorxiv.org/content/10.1101/2023.12.04.569980v1). 

This repository contains all the data (with the exception of raw sequencing data, which can be found at BioProject PRJNA1047182) and code necessary to reproduce figures and tables shown in the publication. 

1. The <strong> data </strong> folder contains all raw data accessed by code distributed throughout the repository, including: <p>
   <strong>16s_data </strong>: DADA2 output and metadata from running <a href = "https://qiime2.org/"> QIIME2 </a>  on 16s amplicon sequencing data of mouse skin and gut <br> 
    <strong> supplement </strong>: Data used for generating supplemental figures 1-6 and supplemental tables 1 and 2 <br>
    <strong> transcriptomics </strong>: Processed RNASeq data used for gene and pathway enrichment from damaged mouse skin at three timepoints (4h, 24h, and 4d) <br>
    As well as TEWL, severity score, and skin CFU data from mouse experiments shown in figures 1-4 <br>

2. The <strong> Data_Visualisation </strong> folder contains an R markdown script that creates all figures in the paper with the exception of Figure 3b and Figure 4b <p>
3. The <strong> Figs </strong> folder contains raw figures produces by the R markdown script above <p>
4. The <strong> QIIME </strong> folder contains the QIIME2 workflow used for analysing amplicon sequencing data, as well as the classifier used. Note: the classifier used was created by Alex Poret for the publication, ["Cutaneous surgical wounds have distinct microbiomes from intact skin"](https://pubmed.ncbi.nlm.nih.gov/36541798/) and can also be found at the associated [github repository](https://github.com/ajporet/cutaneous_wound_microbiome) <p>
5. The <strong> Transcriptomics </strong> folder contains code used to analyse RNASeq data, perform differential gene expression analyses, and pathway enrichment analyses.
6. the <strong> Epidermal Thickness </strong> folder contains code used to calculate epidermal thickness from histology images (distance between points)
