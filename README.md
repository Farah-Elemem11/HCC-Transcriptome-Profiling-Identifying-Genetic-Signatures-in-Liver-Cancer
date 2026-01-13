Transcriptomic Profiling of Hepatocellular Carcinoma (GSE112790)
This project focuses on identifying the molecular signatures and differentially expressed genes (DEGs) in Hepatocellular Carcinoma (HCC).
By analyzing the GSE112790 dataset from the Gene Expression Omnibus (GEO), I aimed to distinguish between tumor and healthy liver tissues to pinpoint potential diagnostic biomarkers.

Workflow & Methodology
The analysis was performed in R using a robust bioinformatics pipeline:

  Data Acquisition: Retrieved raw expression data and metadata using GEOquery.

  Experimental Design: Constructed a design matrix and contrast matrix to compare Tumor vs. Normal samples.

  Statistical Analysis: Utilized the limma package and applied eBayes for empirical Bayes moderation.

  Gene Annotation: Mapped Probe IDs to Gene Symbols using the hgu133plus2.db annotation database.

  Data Visualization: Generated a Volcano Plot for global expression trends and a Heatmap for the Top 50 significant genes using Z-score scaling.


Technical Challenges Faced (Problem Solving) 🛠️
This project was a significant learning journey that involved overcoming several technical hurdles:

   Data Integrity: Navigated through the raw dataset to handle missing values (NAs) and ensured the final results were statistically significant and biologically meaningful.

   Duplicate Gene Symbols: Handled the challenge of multiple probes mapping to the same gene symbol, ensuring the data was properly aggregated for accurate biological interpretation.

   Data Scaling: Fine-tuned the Z-score scaling process for the Heatmap to ensure that the visualization accurately represented the relative expression levels across samples.


Key Scientific Findings 📊

   CYP2C19 (Down-regulated, LogFC ≈ -3.9): Identified as a major biomarker showing a massive decrease in tumor samples, reflecting impaired liver detoxification.

  CAP2 & NUSAP1 (Up-regulated): Pinpointed as potential oncogenic drivers in HCC.

  Clear Clustering: The analysis achieved a 100% distinction between "Normal" and "Tumor" groups in the hierarchical clustering.

Tools & Packages

  Language: R

  BioConductor: limma, GEOquery, hgu133plus2.db

   Visualization: pheatmap, ggplot2, ggrepel

🔗 Connect with Me
LinkedIn: https://www.linkedin.com/in/farah-elemam-107969323



I'm passionate about the intersection of Data Science and Healthcare. Feel free to reach out!

    LinkedIn: [Add your LinkedIn profile link here]
