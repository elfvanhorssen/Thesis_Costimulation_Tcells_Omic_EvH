Master thesis:
Eralin van Horssen; september 2024 - juni 2025

**The Power of Two:** Combinatorial Co-Stimulation of Naïve Human CD8⁺ T Cells Revealed by Transcriptomics, Proteomics, and Epigenomics. 
Bioinformatics and Systems Biology at VU Amsterdam. 

The research was conducted at the Immunology department in the T cell regulation group of Leids Universitair Medisch Centrum.

# Abstract
Activation of naive human CD8+ T cells requires both T cell receptor (TCR) stimulation and co-stimulatory signals. While individual co-stimulatory receptors such as CD28, CD27, and CD137 (4-1BB) have known roles in proliferation, survival, and memory formation, the molecular signature of combined co-stimulatory signals remains unknown. This study addresses that gap in knowledge, by examining how na¨ ıve human CD8+ T cells respond to combinatorial costimulation using a multi-omics approach. Na¨ ıve human CD8+ T cells were stimulated for 72 hours under six co-stimulatory conditions-αCD3, αCD3+αCD27, αCD3+αCD137, αCD3+αCD28, αCD3+αCD27+αCD28, αCD3+αCD28+αCD137- Here, an anti-CD3 antibody was used to mimic TCRengagement, in combination with monoclonal antibodies against CD27, CD28, and/or CD137. Transcriptomic changes were profiled by RNA-sequencing, and proteomic changes were assessed via Tandem Mass Tag-based mass spectrometry. Data from combined stimulation—particularly the αCD3+αCD28+αCD137 condition— showed the strongest molecular response, >4,000 differentially expressed genes and > 800 differentially expressed proteins compared to the αCD3-alone condition. These changes in expression exceeded the changes of expression in the single co-stimulated conditions, showing a synergistic effect. Principal Component Analysis showed that data from the αCD3+αCD28+αCD137 condition diverged furthest from αCD3-alone condition at both transcriptome and proteome levels. Gene set enrichment analysis revealed that αCD3+αCD28+αCD137 co-stimulation enriched the G2/M checkpoint, cholesterol biosynthesis, mTORC1 signaling, and fatty acid metabolism pathways. Ingenuity Pathway Analysis of the αCD3+αCD28+αCD137condition confirmed activation of these pathways. Key features—SLC2A1 and LDHA—were more than doubled by adding αCD137 to αCD3+αCD28 at the transcript level, proteomics showed an additive effect, indicating synergy. Comparison of transcriptomic and proteomic data showed concordant gene–protein changes, with some discordance suggesting post-transcriptional regulation. Together, these results demonstrate that combinatorial co-stimulation results in a strengthened activation state in na¨ ıve human CD8+ Tcells. Our findings on combined co-stimulatory signaling may help identify promising therapeutic targets



# Directories and files 

The structure consists of multiple directories including the R files used.

Output files: the output of the analysis, see last section of the readme. 
Scripts: containing scripts used in this study, seperated in transcriptomics, proteomics, atac, integration and IPA. 
Transcriptomics: containing the codes for quality control, Diferentially expression analysis using DEseq2 and downstream analysis such as the volcano plots. In addition the calculation of z-scores are provided.
Proteomics: containing the codes for quality control, Diferentially expression analysis using Limma and downstream analysis such as the volcano plots. In addition the calculation of z-scores are provided. Most of them are written in a R notebook. 
Atac: containing the code for enrichment analysis and downstream analysis (PCA). 
Integration: containing the scrips for making the heatmap with the expression z-score, venn diagrams and density plots. 
Remaining: how I made the datafiles which were used as input for IPA. 


# Packages
![image](https://github.com/user-attachments/assets/771643b8-1d7d-4178-878f-9785388f9082)
 

# Data availability 
The datasets used in this study are available upon request from the T Cell Regulation Group, Department of Immunology, LUMC.

# Workflow
![image](https://github.com/user-attachments/assets/0f19444d-6289-44c9-9384-a947d24439c5)


# Output
![image](https://github.com/user-attachments/assets/3e107f9c-f903-4df3-9246-7ac44b6f5235)

