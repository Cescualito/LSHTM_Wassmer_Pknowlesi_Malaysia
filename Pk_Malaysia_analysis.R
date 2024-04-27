##############################################################
#  Assessing plasma biomarkers of brain injury in Malaysian  #
#  patients with Plasmodium knowlesi infection               #
#  Analysis of data by Cesc Bertran-Cobo, RA LSHTM           #
##############################################################

# Load packages

  library(data.table)
  library(plyr)
  library(tidyverse)
  library(readr)
  library(readxl)
  library(haven)
  library(ggplot2)
  library(lavaan)
  library(cowplot)
  library(rmarkdown)
  library(gtsummary)
  library(pheatmap)
  library(ggpubr)
  library(knitr)
  library(kableExtra)
  library(plotly)
  library(remotes)
  library(leaflet)
  library(sf)
  library(webshot2)
  library(htmlwidgets)
  library(Hmisc)

# Read Data

# Merge based on sample_ID

# Codebook

  Pk_Malaysia_final_codebook <- describe(Pk_Malaysia_final[, -which(names(Pk_Malaysia_final) == "sample_ID")])
  Pk_Malaysia_final_codebook

# Reshape dataset for figures

  Pk_Malaysia_plots <- Pk_Malaysia_final %>% gather(Biomarker, Levels, 11:63)
  
# Median and average values
  
  Pk_infected <- Pk_Malaysia_final %>% filter(status == "Pk infection")
  summary(Pk_infected)
  Pk_control <- Pk_Malaysia_final %>% filter(status == "Control")
  summary(Pk_control)
  
# Map: Sample collection sites
  
  Pk_Malaysia_provinces <- st_read("/home/cescualito/Dropbox/LABORAL/LSHTM/Wassmer Lab - Cesc Bertran-Cobo/Luminex/2023-05-09_Malaysia_samples/Analysis/Databases/gadm41_MYS_shp/gadm41_MYS_1.shp")
  
  Pk_Malaysia_provinces_CTRL <- Pk_Malaysia_provinces %>%
    filter(NAME_1 %in% c("Johor", "Selangor", "Negeri Sembilan", "Kedah"))
  
  Pk_Malaysia_map <- leaflet(data) %>% # Create a leaflet map 
    addTiles() %>%
    setView(lng = 101.9758, lat = 4.2105, zoom = 6)  # Centered on Malaysia
  
  Pk_Malaysia_map <- Pk_Malaysia_map %>% # Highlight provinces where controls were recruited
    addPolygons(
      data = Pk_Malaysia_provinces_CTRL,
      fillColor = "yellow",  
      fillOpacity = 0.5, 
      stroke = FALSE,
      popup = ~NAME_1
    )
  
  Pk_Malaysia_map_Pk_infected <- Pk_Malaysia_final %>% # Now filter out healthy control rows 
    filter(!is.na(latitude) & !is.na(longitude) & !is.na(NAME_1) & !is.na(hospital))
  
  samples_per_hospital <- Pk_Malaysia_final %>% # This is to add number of samples collected per hospital
    count(hospital) %>%
    rename(samples_collected = n)
  Pk_Malaysia_map_Pk_infected <- left_join(Pk_Malaysia_map_Pk_infected, samples_per_hospital, by = "hospital")
  
  Pk_Malaysia_map <- Pk_Malaysia_map %>% # Add markers for hospitals with infected patients
    addCircleMarkers(
      data = Pk_Malaysia_map_Pk_infected,
      lng = ~longitude,
      lat = ~latitude,
      color = "red",
      radius = 5,
      popup = ~paste("Hospital:", hospital, "<br>Province:", NAME_1, "<br>Samples collected:", samples_collected)
    )
  
  Pk_Malaysia_map
  
# Data normality
  
  Pk_Malaysia_histograms <- ggplot(Pk_Malaysia_plots, aes(x=Levels, fill = status)) +
    geom_histogram(alpha = 0.5) +
    theme_classic() + 
    ggtitle("Histograms: Biomarker levels") + 
    facet_wrap(~Biomarker, scales = "free") + ylab("Count") +
    scale_fill_manual(values=c("black", "red"), name = "Status") 
  
  Pk_Malaysia_histograms
  
  Pk_Malaysia_density <- ggplot(Pk_Malaysia_plots, aes(x=Levels, fill = status)) + 
    geom_density(aes(y=1.1*after_stat(count)), na.rm = TRUE, alpha = 0.75) +
    theme_classic() + 
    ggtitle("Density plots: Biomarker levels") + 
    facet_wrap(~Biomarker, scales = "free") + ylab("Count") +
    scale_fill_manual(values=c("black", "red"), name = "Status")   
  
  Pk_Malaysia_density 
  
# Simple analysis
  
  Pk_Malaysia_summary <- Pk_Malaysia_plots %>%
    group_by(Biomarker, status) %>%
    summarise(
      Participants = sum(!is.na(Levels)),  
      Median = format(median(Levels, na.rm = TRUE), scientific = FALSE),  
      IQR = format(IQR(Levels, na.rm = TRUE), scientific = FALSE)  
    )
  
  Pk_Malaysia_summary_wide <- Pk_Malaysia_summary %>%
    pivot_wider(names_from = status, values_from = c(Participants, Median, IQR)) 
  
  names(Pk_Malaysia_summary_wide) <- c("Biomarker", 
                                       "n CTRL", "n Pk", 
                                       "Median_Control", "Median_Pk_infection", 
                                       "IQR_Control", "IQR_Pk_infection")
  
  Pk_Malaysia_summary_wide$Median_Control <- paste(Pk_Malaysia_summary_wide$Median_Control, 
                                                   "(", Pk_Malaysia_summary_wide$IQR_Control, ")", sep = " ")
  Pk_Malaysia_summary_wide$Median_Pk_infection <- paste(Pk_Malaysia_summary_wide$Median_Pk_infection, 
                                                        "(", Pk_Malaysia_summary_wide$IQR_Pk_infection, ")", sep = " ")
  
  Pk_Malaysia_summary_wide <- Pk_Malaysia_summary_wide[, -c(6, 7)]
  
  names(Pk_Malaysia_summary_wide)[names(Pk_Malaysia_summary_wide) == "Median_Control"] <- "Median (IQR) CTRL"
  names(Pk_Malaysia_summary_wide)[names(Pk_Malaysia_summary_wide) == "Median_Pk_infection"] <- "Median (IQR) Pk"
  
  Pk_Malaysia_final_renamed <- Pk_Malaysia_final
  names(Pk_Malaysia_final_renamed) <- make.names(names(Pk_Malaysia_final_renamed)) # Ooopsies had to correct variable names with symbols such as "-" or else cannot use the function below
  
  biomarkers <- names(Pk_Malaysia_final_renamed)[11:63]
  
  Pk_Malaysia_Bonferroni <- purrr::map_dfr(biomarkers, function(x) {
    f <- as.formula(paste0("`", x, "`", '~ status'))  # Use backticks to handle column names with spaces or symbols
    result <- wilcox.test(formula = f, data = Pk_Malaysia_final_renamed, exact = TRUE)
    return(data.frame(name = x, p.value = result$p.value))
  })
  
  Pk_Malaysia_Bonferroni$p.adj <- p.adjust(Pk_Malaysia_Bonferroni$p,method="bonferroni")
  
  colnames(Pk_Malaysia_Bonferroni)[2:3] <- c("P value", "Corrected P value")
  
  Pk_Malaysia_summary_wide$Biomarker <- make.names(Pk_Malaysia_summary_wide$Biomarker)
  Pk_Malaysia_final_kable <- left_join(Pk_Malaysia_summary_wide, Pk_Malaysia_Bonferroni, by = c("Biomarker" = "name"))
  
  original_biomarker_names <- c("APP", # In the same order as they are in Pk_Malaysia_final_kable
                                "Ang-1", 
                                "Ang-2", 
                                "Ang-2/Ang-1 ratio", 
                                "Aβ(1-42)", 
                                "BDNF", 
                                "BMP-9", 
                                "CCL18", 
                                "CCL2", 
                                "CCL4",
                                "CCL5",
                                "CNTN1",
                                "CRP",
                                "CaBD",
                                "Fetuin A",
                                "GFAP",
                                "GM-CSF",
                                "ICAM-1",
                                "IFNγ",
                                "IL-10",
                                "IL-17",
                                "IL-1RA",
                                "IL-1α",
                                "IL-1β",
                                "IL-2",
                                "IL-4",
                                "IL-6",
                                "IL-8",
                                "KLK6",
                                "Lipocalin 2",
                                "MIF",
                                "MPO",
                                "NCAM-1",
                                "NGF-β",
                                "NRGN",
                                "NSE",
                                "NfL",
                                "OPN",
                                "PDGF AA",
                                "PDGF BB",
                                "Park7",
                                "RAGE",
                                "S100B",
                                "Serpin E1",
                                "TDP-43",
                                "TNF-α",
                                "Tau total",
                                "UCHL-1",
                                "VCAM-1",
                                "VEGF",
                                "YKL40",
                                "vWF",
                                "αSyn")
  Pk_Malaysia_final_kable$Biomarker <- original_biomarker_names
  
  Pk_Malaysia_final_kable <- Pk_Malaysia_final_kable %>%
    kable(caption = "Biomarker median (pg/mL) and IQR, per group") %>%
    kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
    column_spec(column = 6, bold = Pk_Malaysia_final_kable$`P value` < 0.05) %>%
    column_spec(column = 7, bold = Pk_Malaysia_final_kable$`Corrected P value` < 0.05) 
  
  Pk_Malaysia_final_kable
  
# Hierarchical clustering heatmap  
    
  # Brain biomarker subset
    
    Pk_Malaysia_final_brain_subset_list <- c(
      "αSyn", "APP", "Aβ(1-42)", "BDNF", "CaBD", "CNTN1", "NSE",
      "Fetuin A", "GFAP", "KLK6", "NCAM-1", "Lipocalin 2", "NGF-β", "NfL", "NRGN",
      "Park7", "S100B", "TDP-43", "Tau total", "UCHL-1", "YKL40"
    )
    
    Pk_Malaysia_final_brain_subset <- Pk_Malaysia_final %>% select(Participants, Pk_Malaysia_final_brain_subset_list)
    
    row.names(Pk_Malaysia_final_brain_subset) <- Pk_Malaysia_final_brain_subset$Participants
    Pk_Malaysia_final_brain_subset$Participants <- NULL
    
    Pk_Malaysia_final_brain_subset[is.na(Pk_Malaysia_final_brain_subset)] <- 0 # Replace NAs with zeros
    Pk_Malaysia_final_brain_subset <- scale(Pk_Malaysia_final_brain_subset) # Scale data 
    
    Pk_Malaysia_final_brain_subset_dist <- dist(Pk_Malaysia_final_brain_subset) # Base R code
    Pk_Malaysia_final_brain_subset_hclust <- hclust(Pk_Malaysia_final_brain_subset_dist)
    plot(Pk_Malaysia_final_brain_subset_hclust)
    
    pheatmap(Pk_Malaysia_final_brain_subset, cutree_rows = 2, cutree_cols = 2, scale = "none")
    
  # Immune biomarker subset
    
    Pk_Malaysia_final_immune_subset_list <- c(
      "CRP", "GM-CSF", "IFNγ", "IL-1α", "IL-1β", "IL-2", "IL-6", "IL-8",
      "IL-17", "MIF", "MPO", "TNF-α", "IL-1RA", "IL-4", "IL-10",
      "CCL2", "CCL4", "CCL5", "CCL18", "OPN", "RAGE"
    )
    
    Pk_Malaysia_final_immune_subset <- Pk_Malaysia_final %>% select(Participants, Pk_Malaysia_final_immune_subset_list)
    
    row.names(Pk_Malaysia_final_immune_subset) <- Pk_Malaysia_final_immune_subset$Participants
    Pk_Malaysia_final_immune_subset$Participants <- NULL
    
    Pk_Malaysia_final_immune_subset[is.na(Pk_Malaysia_final_immune_subset)] <- 0 # Replace NAs with zeros
    Pk_Malaysia_final_immune_subset <- scale(Pk_Malaysia_final_immune_subset) # Scale data 
    
    Pk_Malaysia_final_immune_subset_dist <- dist(Pk_Malaysia_final_immune_subset) # Base R code
    Pk_Malaysia_final_immune_subset_hclust <- hclust(Pk_Malaysia_final_immune_subset_dist)
    plot(Pk_Malaysia_final_immune_subset_hclust)
    
    pheatmap(Pk_Malaysia_final_immune_subset, scale = "none")
    
  # Vascular biomarker subset
    
    Pk_Malaysia_final_vascular_subset_list <- c(
      "Ang-1", "Ang-2", "Ang-2/Ang-1 ratio", "BMP-9", "ICAM-1",
      "PDGF AA", "PDGF BB", "Serpin E1", "VCAM-1", "VEGF", "vWF"
    )
    
    Pk_Malaysia_final_vascular_subset <- Pk_Malaysia_final %>% select(Participants, Pk_Malaysia_final_vascular_subset_list)
    
    row.names(Pk_Malaysia_final_vascular_subset) <- Pk_Malaysia_final_vascular_subset$Participants
    Pk_Malaysia_final_vascular_subset$Participants <- NULL
    
    Pk_Malaysia_final_vascular_subset[is.na(Pk_Malaysia_final_vascular_subset)] <- 0 # Replace NAs with zeros
    Pk_Malaysia_final_vascular_subset <- scale(Pk_Malaysia_final_vascular_subset) # Scale data 
    
    Pk_Malaysia_final_vascular_subset_dist <- dist(Pk_Malaysia_final_vascular_subset) # Base R code
    Pk_Malaysia_final_vascular_subset_hclust <- hclust(Pk_Malaysia_final_vascular_subset_dist)
    plot(Pk_Malaysia_final_vascular_subset_hclust)
    
    pheatmap(Pk_Malaysia_final_vascular_subset, cutree_rows = 2, cutree_cols = 2,  scale = "none")
  
# Biomarker correlation with age
  
  Pk_Malaysia_age <- Pk_Malaysia_plots %>%    
    ggplot(aes(y = Levels, x = age, color=status, group = sample_ID)) +
    geom_point(size = 1.5) + 
    geom_smooth(method = lm, formula = y ~ x, se = FALSE) + 
    xlab("Age") + ylab("pg/mL") +
    theme_classic() + 
    ggtitle("") + 
    facet_wrap(~Biomarker, scales = "free") +
    scale_color_manual(values=c("black", "red"), name = "Status") + 
    stat_cor(method = "spearman") + 
    geom_smooth(method = lm, formula = y ~ x, se = FALSE)
  
  Pk_Malaysia_age
  
# Biomarker correlation with parasitaemia
  
  Pk_Malaysia_par <- Pk_Malaysia_plots %>%    
    ggplot(aes(y = Levels, x = parasitaemia, color=status, group = sample_ID)) +
    geom_point(size = 1.5) + 
    geom_smooth(method = lm, formula = y ~ x, se = FALSE) + 
    xlab("Parasites/μL") + ylab("pg/mL") +
    theme_classic() + 
    ggtitle("") + 
    facet_wrap(~Biomarker, scales = "free") +
    scale_color_manual(values=c("black", "red"), name = "Status") + 
    stat_cor(method = "spearman") + 
    geom_smooth(method = lm, formula = y ~ x, se = FALSE) +
    scale_x_continuous(labels = scales::scientific_format(digits = 1, decimal.mark = ",", big.mark = "x10^")) +
    scale_y_continuous(labels = scales::scientific_format(digits = 1, decimal.mark = ",", big.mark = "x10^"))
  
  Pk_Malaysia_par
  