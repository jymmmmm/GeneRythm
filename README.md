
# Getting Started
## 1. Introduction
scGeneRhythm is an innovative tool designed for single-cell RNA sequencing (scRNA-seq) analysis. It utilizes deep learning techniques and Fourier transformation to examine gene expression data, aiming to reveal gene relationships and provide valuable biological insights.
<figure>
  <img src="images/F1.jpg" >
  <figcaption>Overview of scGeneRythm. a. The flowchart of scGeneRythm. Colored solid circles represent cells mapped to different key segments of a trajectory tree. Average of gene expressions are performed for each key segment. The subsequent step involved obtaining the green-line frame representing the time vector through the calculation of differences, followed by utilizing the Fast Fourier Transform (FFT) to obtain the blue-line frame representing the frequency vector. b. The structure of the machine learning model comprises Graph Convolutional Networks (GCN) along with separate encoders and decoders designed to handle temporal and frequency information. c. For gene embedding clustering and subsequent analysis, including Gene Ontology (GO) terms analysis. This will identify different time and frequency pattern of each gene clusters.</figcaption>
</figure>

### Key Features and Functions

- **Comprehensive Analysis:** GeneRhythm conducts simultaneous analyses in both the time and frequency domains, offering a more holistic perspective.

- **Uncovering Hidden Relationships:** This tool helps identify genes that may appear distinct in the time domain but exhibit similarities in the frequency domain.

- **Outstanding Performance:** scGeneRhythm has demonstrated exceptional performance in gene clustering and pathway enrichment analysis across various single-cell datasets.

- **In-Depth Gene Rhythmicity:** scGeneRhythm assists in uncovering gene rhythmicity and hidden associations, contributing to advancements in cellular and molecular biology research.

## 2. How to use



## 3. Results

<figure>
  <img src="images/F2.jpg">
  <figcaption>Results of scGeneRythm on Human fetal immune cells dataset. a. UMAP visualization of Human fetal immune cells dataset. b. To process an inferred gene trajectory within a UMAP plot using scGeneRythm and subsequently cluster the resulting gene embeddings. c. Subsequent analysis of gene cluster 7, which includes inspecting its temporal and frequency domain patterns as well as conducting Gene Ontology (GO) term analysis. d. Subsequent analysis of gene cluster 8. e. A specific case in cluster 8 and gene expression heatmap.</figcaption>
</figure>

<figure>
  <img src="images/F3.jpg" alt="图片的替代文本">
  <figcaption>scGeneRhythm effectively captures gene clusters in Mouse embryonic blood cells dataset. a. UMAP visualization of Mouse embryonic blood cells dataset. b. To process an inferred gene trajectory within a UMAP plot using scGeneRythm and subsequently cluster the resulting gene embeddings. c. Subsequent analysis of gene cluster 1, which includes inspecting its temporal and frequency domain patterns as well as conducting Gene Ontology (GO) term analysis. d. Subsequent analysis of gene cluster 4. e. A specific case in cluster 4 and gene expression heatmap.</figcaption>
</figure>

<figure>
  <img src="images/F4.jpg">
  <figcaption>Results of scGeneRythm on Mouse embryonic neural crest cells dataset. a. UMAP visualization of Mouse embryonic neural crest cells. b. To process an inferred gene trajectory within a UMAP plot using scGeneRythm and subsequently cluster the resulting gene embeddings. c. Subsequent analysis of gene cluster 5, which includes inspecting its temporal and frequency domain patterns as well as conducting Gene Ontology (GO) term analysis. d. Subsequent analysis of gene cluster 8. e. A specific case in cluster 5 and gene expression heatmap.</figcaption>
</figure>

<figure>
  <img src="images/F5.jpg">
  <figcaption>Comparison between scGeneRythm and scSTEM. a. Perform GO term analysis on the top 4 clusters identified by the scSTEM, and present the top 10 terms with the lowest p-values separately. b. Perform GO term analysis on the top 4 clusters identified by the scGeneRythm and present the top 10 terms. Annotations labeled as 't' represent the temporal patterns of clusters, while those labeled as 'f' represent the frequency domain patterns of clusters. c. a comparative analysis of p-values in GO term analysis for clusters from both scSTEM and scGeneRythm using bar charts.</figcaption>
</figure>

<figure>
  <img src="images/F6.jpg">
  <figcaption>Ablation experiment for scGeneRythm. a. Present the results of scGeneRythm's ablation experiments using a bar chart, including experiments with time information only, frequency information only, and our proposed method. b. Demonstrate that scSTEM's focus solely on time data while neglecting frequency data is impractical by showcasing the patterns of time and frequency of MTSS1 and SLC9A9. c. Illustrate the patterns of GRB10 and PID1 in both the time domain and the frequency domain, thereby demonstrating the significance of frequency information.</figcaption>
</figure>

