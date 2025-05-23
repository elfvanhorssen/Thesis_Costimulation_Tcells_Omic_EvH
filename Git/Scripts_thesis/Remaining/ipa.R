library(DESeq2)
library(openxlsx)  # of gebruik library(writexl) als je dat prefereert

# Definieer de groepen
groupscsp <- c("Unstim", "aCD3", "aCD3_aCD27", "aCD3_a4_1BB", "aCD3_aCD28", "aCD3_aCD27_aCD28", "aCD3_aCD28_a4_1BB")

# Maak een lege lijst voor de resultaten (optioneel)
results_listcsp <- list()

# Loop over alle vergelijkingen
for (i in 1:(length(groupscsp) - 1)) {
  for (j in (i + 1):length(groupscsp)) {
    
    # Definieer het contrast voor DESeq2
    contrast <- c("condition", groupscsp[i], groupscsp[j])
    res_csp <- results(dds_csp, contrast = contrast)
    
    # Geef de vergelijking een naam, bv. "Unstim_vs_aCD3"
    comparison_name <- paste0(groupscsp[i], "_vs_", groupscsp[j])
    
    # Zet het resultaat om naar een data frame
    res_csp_df <- as.data.frame(res_csp)
    
    # Voeg de genen (rownames) toe als een kolom
    res_csp_df$gene <- rownames(res_csp_df)
    
    # (Optioneel) Zet de kolom met genen als eerste kolom
    res_csp_df <- res_csp_df[, c("gene", setdiff(colnames(res_csp_df), "gene"))]
    
    # Bewaar het resultaat in de lijst (optioneel)
    results_listcsp[[comparison_name]] <- res_csp_df
    
    # Stel de bestandsnaam samen; pas pad en naam aan naar wens
    file_name <- paste0("I:/Research/TCR/2.Lab Members/Eralin van Horssen/co_stim_project/1.unbiased_rnaseq/results/2.DEG/result_dds_", comparison_name, ".xlsx")
    
    # Schrijf het resultaat weg als een Excel-bestand
    write.xlsx(res_csp_df, file = file_name)
  }
}
