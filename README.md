# Imp/Syp iCLIP Data Explorer

This repository contains the source code for an interactive Shiny web application designed to explore the data from the publication:

> **Imp/IGF2BP and Syp/SYNCRIP temporal RNA interactomes uncover combinatorial networks of regulators of Drosophila brain development** > Lee JY, et al. (2025) *Science Advances* > DOI: [10.1126/sciadv.adr6682](https://www.science.org/doi/10.1126/sciadv.adr6682)

The application provides a user-friendly interface to browse the complete iCLIP dataset, search for specific genes, and visualize the RNA-binding protein (RBP) binding profiles on the UCSC Genome Browser.

---

## Live Application

**The interactive web application can be accessed here:**

[**>> Launch Imp/Syp iCLIP Data Explorer <<**](https://jefflee1103.github.io/Imp-Syp-iCLIP-shiny-app/)

---

## How to Use the Web App

The application is organized into three main tabs, each providing a different way to interact with the data.

### 1. Introduction

This tab provides the scientific background and rationale for the study. It details the experimental question, the approach using individual-nucleotide resolution UV cross-linking and immunoprecipitation (iCLIP), and presents key figures that summarize the experimental design and core concepts.

### 2. Explore Imp/Syp target table

This tab presents the complete, unabridged dataset of Imp and Syp iCLIP targets identified in the study. The interactive table allows you to:
* **Search and Filter:** Use the search boxes at the top of each column to dynamically filter the data for specific genes, features, or value ranges.
* **Sort Data:** Click on any column header to sort the table accordingly.
* **Download Data:** A "Download Full Table (.csv)" button is provided to save the entire dataset for offline analysis.

### 3. UCSC Genome Browser

This tab provides a direct link to a pre-configured session on the UCSC Genome Browser. Clicking the "Open UCSC Genome Browser Session" button will open a new browser tab displaying all the iCLIP tracks from this study. This allows for in-depth, interactive exploration of the RBP binding sites in their genomic context.

---

## Citation

If you use the data or this application in your research, please cite the original publication:


Jeffrey Y. Lee et al. ,Imp/IGF2BP and Syp/SYNCRIP temporal RNA interactomes uncover combinatorial networks of regulators of Drosophila brain development.Sci. Adv.11,eadr6682(2025).DOI:10.1126/sciadv.adr6682. 

## Source data

Source data (including track bigWig/Bed files) used for the publication is available from [https://github.com/jefflee1103/Lee2024_Imp-Syp-iCLIP](https://github.com/jefflee1103/Lee2024_Imp-Syp-iCLIP).  

Raw and processed FASTQ files are available from Gene Expression Omnibus (GEO) [GSE262499](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262499).




## Contact

For any questions or inquiries, please contact:
* **Jeffrey Y Lee**: [jeff.lee@glasgow.ac.uk](mailto:jeff.lee@glasgow.ac.uk)
* **Ilan Davis**: [ilan.davis@glasgow.ac.uk](mailto:ilan.davis@glasgow.ac.uk)
