# Merged Reference Papers

Generated on: Mon Dec  1 14:30:54 CST 2025

---

# Paper 1: s41467-023-37168-7

# A comprehensive benchmarking with practical guidelines for cellular deconvolution of spatial transcriptomics 

Received: 30 September 2022

## Accepted: 3 March 2023

Published online: 21 March 2023

Check for updates

Haoyang Li ${ }^{(1,2,6}$, Juexiao Zhou ${ }^{1,2,6}$, Zhongxiao Li ${ }^{1,2}$, Siyuan Chen ${ }^{(1)}{ }^{1,2}$, Xingyu Liao $\boldsymbol{(})^{\mathbf{1} \boldsymbol{,} \mathbf{2}}$, Bin Zhang ${ }^{\mathbf{1} \boldsymbol{,} \mathbf{2}}$, Ruochi Zhang ${ }^{\mathbf{3}}$, Yu Wang ${ }^{\mathbf{3}}$, Shiwei Sun ${ }^{\mathbf{4 , 5}}$ \& Xin Gao ${ }^{(1)}{ }^{1,2}$


#### Abstract

Spatial transcriptomics technologies are used to profile transcriptomes while preserving spatial information, which enables high-resolution characterization of transcriptional patterns and reconstruction of tissue architecture. Due to the existence of low-resolution spots in recent spatial transcriptomics technologies, uncovering cellular heterogeneity is crucial for disentangling the spatial patterns of cell types, and many related methods have been proposed. Here, we benchmark 18 existing methods resolving a cellular deconvolution task with 50 real-world and simulated datasets by evaluating the accuracy, robustness, and usability of the methods. We compare these methods comprehensively using different metrics, resolutions, spatial transcriptomics technologies, spot numbers, and gene numbers. In terms of performance, CARD, Cell2location, and Tangram are the best methods for conducting the cellular deconvolution task. To refine our comparative results, we provide decision-tree-style guidelines and recommendations for method selection and their additional features, which will help users easily choose the best method for fulfilling their concerns.


Spatial transcriptomics technologies, named "Method of the Year 2020" ${ }^{1}$, have undergone rapid development in recent years. They are used to profile spatial locations of all detected mRNAs, providing a new perspective for biologists seeking to understand cells per se as well as their microenvironments. Broadly, spatial transcriptomics technologies can identify undiscovered transcriptional patterns and reconstruct transcriptional panoramas of whole tissues. On a finegrained level, these technologies can be used to explore the interactions among neighboring cells and intracellular and extracellular states, which helps redefine the function of cells and improves our knowledge of diseases ${ }^{2}$. The current spatial transcriptomics technologies can be mainly classified into two categories. The first
category is image-based technologies, including in situ sequencingand in situ hybridization-based methods ${ }^{3}$, which can profile mRNA with high spatial resolution, especially at the subcellular level. However, limitations such as the low number of profiled genes, low sensitivity of mRNA detection, and time-consuming processes impede the broad application of image-based technologies. The second category is sequencing-based spatial transcriptomics technologies, which capture position-barcoded mRNA with nongene-specific probes. These technologies can profile the whole transcriptome of tissue sections of any size, and are more user-friendly and less timeconsuming than image-based technologies ${ }^{4}$. Moreover, spatial transcriptomics technologies are highly applicable and have been used to

[^0]improve our understanding of various species, organs, and tissues, including the brain ${ }^{5}$, liver ${ }^{6}$, and tumors ${ }^{7}$.

One critical issue related to sequencing-based spatial transcriptomics technologies is low-resolution spots containing multiple cells with several blended cell types, which can conceal the genuine transcriptional pattern and lead to biological misunderstanding of the tissue resulting in the distorted cellular-level reconstruction of the tissue. An important task, therefore, is to quantify the proportion of all cell types among captured spots, so-called cellular deconvolution. Following deconvolution, all captured spots can be used to better understand intercellular functions and recover the fine-grained panorama of a heterogeneous tissue.

A recent benchmarking study ${ }^{8}$ was focused on single-cell RNA sequencing (scRNA-seq) and spatial transcriptomics integration methods. There is another recent benchmarking study ${ }^{9}$, in which the number of related methods is limited and scRNA-seq reference-free methods are not considered. Despite these efforts, clear guidelines and solid recommendations for the users are still lacking for the comprehensive coverage of available methods.

In the present study, we conducted a comprehensive benchmarking and provided guidelines for the cellular deconvolution of spatial transcriptomics data. Specifically, we evaluated 18 existing computational methods with 50 simulated and real-world datasets by comprehensively testing the accuracy, robustness, and usability of the methods. These methods could be broadly classified as those with and without scRNA-seq references. Based on their computational techniques, we grouped the methods as follows: probabilistic-based, nonnegative matrix factorization-based (NMF-based), graph-based, optimal-transport (OT)-based and deep learning-based methods. During benchmarking, we used multiple metrics and various data resources with different spatial transcriptomics techniques, spot resolutions, gene numbers, spot numbers, and cell types to ensure our assessment was comprehensive and to deepen our understanding of the cellular deconvolution methods. In addition to the quantification and visualization processes, decision-tree-style guidelines were produced, which included the refinement of the benchmarking results and the collection of respective additional features of the methods detailed in related publications. These guidelines recommend scenario-specific methods for users considering computational efficiency and the characteristics of data resources. The general limitations and future perspectives associated with cellular deconvolution are also discussed to give users a clear picture of the cellular deconvolution field and thus facilitate the improvement of tools for the community.

## Results

## Benchmarking pipeline

To evaluate cellular deconvolution methods comprehensively, we identified 18 existing methods from published and preprint papers as follows: CARD ${ }^{10}$, Cell2location ${ }^{11}$, RCTD ${ }^{12}$, DestVI ${ }^{13}$, stereoscope ${ }^{14}$, SpatialDecon ${ }^{15}$, STRIDE ${ }^{16}$, NMFreg ${ }^{17}$, SpatialDWLS ${ }^{18}$, SPOTlight ${ }^{19}$, DSTG ${ }^{20}$, SD ${ }^{221}$, Tangram ${ }^{22}$, Berglund ${ }^{23}$, SpiceMix ${ }^{24}$, STdeconvolve ${ }^{25}$, SpaOTsc ${ }^{26}$ and novoSpaRc ${ }^{27}$. According to the data resources used, Berglund, SpiceMix, and STdeconvolve were scRNA-seq reference-free methods that identified cell-type-specific spatial patterns using only the information from the spatial locations of spots and their gene expression profiles without any reliance on external scRNA-seq data. The remaining 15 methods required scRNA-seq data from the same tissue as the spatial transcriptomics data. Cell-type annotations and cell-type-specific gene expression profiles from scRNA-seq data can help optimize the proportion of all cell types in the spatial transcriptomics data. The 18 methods were classified based their computational techniques as follows. Probabilistic-based methods: Berglund, Cell2location, DestVI, RCTD, SpatialDecon, stereoscope, STRIDE, and STdeconvolve; NMF-based methods: CARD, NMFreg, SpatialDWLS, SPOTlight, and SpiceMix; graph-based methods: DSTG and SD²; deep
learning-based method: Tangram; and OT-based: SpaOTsc and novoSpaRc. These five computational techniques introduced different methods to formulate the spatial transcriptomics data (sometimes with the scRNA-seq data) and solve the cellular deconvolution problem (Fig. 1A).

We also collected seven image-based and sequencing-based spatial transcriptomics datasets: seqFISH+ ${ }^{28}$, MERFISH ${ }^{29,30}$, Spatial Transcriptomics (ST) ${ }^{31}$, 10 X Visium (Visium) ${ }^{32}$, Slide-seqV2 ${ }^{33}$, and stereoseq ${ }^{34,35}$. Their corresponding scRNA-seq datasets were collected as complementary resources (Supplementary Table 1). Among these data resources, the image-based spatial transcriptomics data (seqFISH+ and MERFISH) contained gene expression profiles, spatial locations, and cell-type annotations of individual cells, which could be used to simulate low-resolution spots by binning the cells with a unified square size, and the ground truth could be calculated according to the number of cells with different cell types in each spot. Simulated data could be used to generate different resolutions of spots by defining the different sizes of the binning squares. For the sequencing-based spatial transcriptomics data (ST, Visium, Slide-seqV2, and stereo-seq), realworld scenarios were included that could be solved by the cellular deconvolution methods. The resolution statistics and the number of genes and spots from the abovementioned data resources are plotted in Fig. 1B. Notably, the spots from stereo-seq were of a subcellular resolution ( 500 nm ), and we binned the stereo-seq spots to reduce the resolution for cellular deconvolution (see "Methods").

After obtaining deconvolution results from all 18 methods on all spatial transcriptomics datasets, we assessed the performances of the methods comprehensively according to the following quantities (Fig. 1C): (1) the accuracy of the deconvolution results, which was evaluated using multiple metrics on all methods and datasets; (2) the robustness of all methods tested on different conditions (spatial transcriptomics techniques, number of genes, number of spots, and number of cell types); (3) the usability of all tools, including computational efficiency, quality of documents, publications, and code. To display the benchmarking results of all methods intuitively, a table containing an evaluation of all metrics according to the three listed qualities is provided (Fig. 2), in which darker dots represent better performance (see "Methods"). Moreover, detailed decision-tree-style guidelines are provided, which include scenariospecific recommendations of the methods and a summary of their respective additional features detailed in their related publications.

## Accuracy

To assess accuracy of methods, we used multiple metrics to quantify the performance of deconvolution, including the Jensen-Shannon divergence (JSD) score, root-mean-square error (RMSE), and Pearson correlation coefficient (PCC). For simulated data (MERFISH and seqFISH+), we used the JSD score and RMSE to measure the distance between the predicted cell-type proportion and ground truth (see "Methods"). Because seqFISH+ and MERFISH have single-cell resolutions, we binned them by size ( $51.5 \mu \mathrm{~m}$ and $100 \mu \mathrm{~m}$, respectively) to simulate low-resolution spots. According to the number of simulated spots and genes detected for seqFISH+ and MERFISH, seqFISH+ ( 71 spots and 10,000 genes) had a low number of spots but a high number of genes, whereas MERFISH ( 3067 spots and 135 genes) exhibited the opposite trend. These two data resources were complementary in terms of the number of spots and genes, and hence provided helpful results to observe the performance of all methods at these two extremes. For each sequencing-based dataset (ST, Visium, Slide-seqV2, and stereo-seq), we chose several cell types and their known marker genes to compare the PCC between the spatial distribution of cell-type proportion and marker-gene expression, which was used to quantify the performance of deconvolution without any ground truth (Supplementary Table 2; see the "Methods").
![](https://cdn.mathpix.com/cropped/2025_11_30_20f0f0568aacc5815037g-03.jpg?height=1252&width=1686&top_left_y=158&top_left_x=199)

Fig. 1 | The summarization of benchmarking pipeline. A Eighteen cellular deconvolution methods, classified based on their data requirements and computational techniques, were evaluated with 50 simulated and real-world spatial transcriptomics datasets. B Datasets from six spatial transcriptomics technologies were
used for benchmarking, and the scatter plot shows the resolution of each spot and number of spots and genes in each technology. C The benchmarking results were measured according to the accuracy, robustness, and usability of the methods.

The performance of each method varied with the data resources used, but some methods performed steadily great accuracy with both simulated and real-world datasets, e.g., Cell2location and DestVI (Fig. 2). Detailed quantifications of JSD, RMSE, and PCC are also provided (Supplementary Figs. 1, 5, 9-11). Using simulated data, most of the methods performed well with MERFISH, but only CARD, DestVI, and SpatialDWLS were high-performing methods with seqFISH+, indicating that they worked well with a low number of spots. When the number of spots was higher (i.e., with MERFISH and Slide-seqV2), Cell2location, SpatialDecon, and Tangram were most capable of performing deconvolution with large views of tissues. In addition, SpatialDWLS performed well with simulated datasets but poorly with all real-world datasets. The whole cell-type proportions and their ground truth for all six cell types in seqFISH+ and 12 samples in MERFISH were visualized among all methods (Fig. 3A, B and Supplementary Figs. 2-4, 11-13). To explore the performances associated with these six cell types, spider plots were produced to show the RMSE of the 18 methods for each cell type (Fig. 3C). Excitatory neurons had the highest cell-type abundance in seqFISH+, and this dataset had the highest RMSE among all methods. The same trend was observed for inhibitory neurons in MERFISH. With the sequencing-based datasets, the PCC of ST did not show a strong relationship and distinct spatial patterns among all methods. In the ST dataset, the lowest resolution of the spots of all existing spatial transcriptomics technologies and a highly heterogenous tissue sample of pancreatic cancer led to the disappearance of
the spatial patterns of individual cell types and blurred the relationships between cell types and selected marker genes. Nevertheless, some of the methods still achieved relatively high PCCs, e.g., CARD, Cell2location, SpatialDecon, and stereoscope. With the Visium dataset, most methods performed well with three paired cell types and marker genes. Spatial patterns were unclear, which was related to the choice of cell types rather than a heterogeneity issue in the datasets. With the Slide-seqV2 and stereo-seq datasets, the spatial patterns were distinct and most of the methods achieved relatively high PCCs, especially Cell2location, STdeconvolve, and RCTD, which were the top-three performing methods. Among the top-three least-performing methods, DestVI tended to output the average cell-type proportions, whereas SpiceMix and SpatialDWLS generated the mapping with much noise which means that they could not distinguish the cell type patterns well in the real-world datasets.

## Robustness

To evaluate the robustness of the 18 tested methods, we designed several experiments with different conditions as follows: (1) 10,000, 6000, or 3000 genes were randomly chosen in the seqFISH+ dataset and $26,365,18,000$, and 9000 in zebrafish embryo dataset by stereoseq; (2) three resolutions were simulated with 12 MERFISH datasets and zebrafish embryo dataset using binning sizes of 20,50 , and $100 \mu \mathrm{~m}$ and 5, 10, and $15 \mu \mathrm{~m}$, respectively, which were also used to test performance with different numbers of spots; (3) 17 original cells types
![](https://cdn.mathpix.com/cropped/2025_11_30_20f0f0568aacc5815037g-04.jpg?height=861&width=1677&top_left_y=165&top_left_x=206)

Fig. 2 | The summary table of the performance of all methods. We visualized their performance in terms of accuracy (red), robustness (blue), and usability (green), and we listed the requirement of spatial location, programming language, and the overall performance (gray) for each method. For all colored dots, a darker
shade represents better performance. The black dots shown in normalization and hyperparameter of robustness meant that the methods required the raw count as input spatial transcriptomics data only or do not have hyperparameters to regulate. Source data are provided as a Source Data file.
and 11 integrated cell types were tested in Slide-seqV2 datasets; (4) two kinds of normalization methods in the Visium dataset to test the effect of normalization of input spatial transcriptomics data; (5) three different values for each chosen hyperparameters in Visium and Slide-seq V2 datasets among five pairs of cell types and their marker genes; and (6) the stability of the performance was assessed by repeating the experiments three times with the seqFISH+ dataset using 10,000 genes per spot.

We found that the general JSD and RMSE of almost all methods were barely changed with the seqFISH+ datasets; the exception was SpiceMix, which performed surprisingly weakly with 3000 genes (Supplementary Fig. 5 and Supplementary Dataset 1). In the tests of different spot numbers and resolutions with MERFISH datasets, the performance of all methods worsened with an increasing number of spots, and this tendency was most evident with SpiceMix and STRIDE (Supplementary Fig. 9 and Supplementary Dataset 1). We visualized the JSD and RMSE using all 12 MERFISH samples with three resolutions and all methods, and the patterns mentioned above were more distinct, although Tangram, stereoscope, and DestVI performed steadily with all 36 datasets (Supplementary Fig. 10). However, on these datasets, DSTG and SpiceMix did not perform well. In particular, DSTG could not identify the clear patterns of different cell types from visualization results and SpiceMix did not perform well possibly due to the misalignment of predicted topics and biological cell types. Using the SlideseqV2 datasets, we chose two paired sub-cell types from the 17 original cell types and combined each pair as new cell types with biological senses (see "Methods"). For example, two sub-cell types, CA1 and CA3, are the main subfields of the hippocampus proper with their own spatial pattern at the start and end of the neural circuit ${ }^{31,32}$. We combined these two sub-cell types under the name "Cornu Ammonis" (CA), which was the former name of the hippocampus. The spatial patterns from the sub-cell types (CA1 and CA3) and integrated cell types were located from the results of 17 and 11 cell types in Slide-seqV2 (Supplementary Figs. 14 and 15). Based on the visualization and PCCs of the Slide-seqV2 datasets, we found that SpatialDecon, Tangram, and RCTD
had the capability to handle datasets with fine-grained sub-cell types (Supplementary Fig. 11). For the zebrafish embryo dataset by stereoseq, SpatialDecon, Tangram and CARD showed great performance in terms of both different spot numbers and different gene numbers among three different kinds of cell types (Supplementary Figs. 17 and 18). The effects of normalization methods to the performance were also evaluated on the Visium dataset by comparing the results from input spatial transcriptomics data as the raw count, and two normalization functions: lognorm and sctransform ${ }^{36}$. These two normalization functions were commonly used methods in the analysis of scRNA-seq data. The results showed that many methods (e.g., DestVI, SpatialDecon and SPOTlight) performed better using raw data than using normalized data by lognorm (Supplementary Fig. 19), mainly because there were default normalization procedures in these methods. Thus, lognorm would repeatedly normalize the data which resulted in worse performance than directly inputting the raw data. On the other hand, some of methods (e.g., SpaOTsc and Tangram) did not have any default normalization procedure in their pipelines which caused better performance when being normalized by lognorm. But all the methods performed worse with the sctransform normalization (Supplementary Fig. 19). In this evaluation, Cell2location and STdeconvolve were not included because they required to use the raw count as input data. With respect to the effect of hyperparameters, three different values of each hyperparameter for individual methods were chosen to calculate the variance of PCCs (Supplementary Table 5). This experiment was conducted on Visium and Slide-seq V2 datasets over five pairs of cell types and their marker genes. The results reflected that most of the methods were stable enough among different hyperparameters except DSTG whose variances were over 0.01 (Supplementary Fig. 20). To test the stability of all methods, we repeated the experiments three times with seqFISH+ and 10,000 genes per spot and found that 13 methods showed a steady performance with three identical results. Of the remaining five methods, $\mathrm{SD}^{2}$ and DSTG exhibited high variance in the JSD and RMSE, which was related to their strategies of using synthesized pseudospots (Supplementary Fig. 21).
![](https://cdn.mathpix.com/cropped/2025_11_30_20f0f0568aacc5815037g-05.jpg?height=1856&width=1686&top_left_y=158&top_left_x=199)

Fig. 3 | The performance of all methods for simulated datasets. A Visualization of the ground truth and predicted the proportions of excitatory neurons for 18 methods with the seqFISH+ datasets and 10,000 genes per spot. B Visualization of the ground truth and predicted results of deconvolution for all methods with MERFISH datasets ( $100 \mu \mathrm{~m}$ resolution per spot). The six cell types are represented
by the six different colors shown in (C). C Spider plots showing the RMSE of the deconvolution results for the 18 methods among 6 cell types from the MERFISH $(100,50$, and $20 \mu \mathrm{~m}$ resolution per spot) and seqFISH $+(10,000,6000$, and 3000 genes per spot) datasets. Source data are provided as a Source Data file.

Generally, CARD, Cell2location, Tangram, and $\mathrm{SD}^{2}$ were the most robust methods according to their performance with different resolutions, number of genes, number of spots, and number of cell types (Fig. 2).

## Usability

Besides testing the performance of the methods with different situations and datasets, we also assessed their computational efficiency and
user-friendliness, which are important factors to users. To fulfill the main concerns of users, we recorded the running time with three different spot numbers in the MERFISH datasets and stereo-seq dataset, and scored several aspects of the tutorials of all methods, including document quality, code quality, installation procedure, compatibility for operating system and example analysis.

According to running time (Supplementary Tables 3 and 4), NMFreg, STRIDE, and Tangram were the most efficient methods. As the

![](https://cdn.mathpix.com/cropped/2025_11_30_20f0f0568aacc5815037g-06.jpg?height=858&width=1684&top_left_y=160&top_left_x=201)
Fig. 4 | Scenario-specific decision-tree-style guidelines for users. Four common scenarios are included, and three methods are recommended for each branch.

default hyperparameter is set when using all methods, hyperparameter selection has a substantial effect on method efficiency. In terms of the quality of tutorials and code, most methods satisfied the basic requirements of users. In particular, CARD, Cell2location, RCTD, and DestVI were highly user-friendly with helpful tutorials and readable code that were easy for users to implement.

## Guidelines

Considering the performance of each method and the features described in their related publications, we provided scenario-specific recommendations and guidelines for the methods according to four crucial scenarios (Fig. 4). Because users usually pay less attention to the computational techniques of methods, all scenarios were related to the characteristics of data and computational efficiency. For example, the first scenario was the absence of scRNA-seq reference data from the same tissue, which could be considered a general scenario for users with several possible concerns: (1) the lack of scRNA-seq data from some markedly heterogeneous tissue sample (e.g., cancerous tissue); (2) missing or inconsistent cell types annotated between scRNA-seq and spatial transcriptomics data; and (3) the platform effects of scRNA-seq and spatial transcriptomics data. For each branch of the decision tree, there were three recommended methods. Tangram and Cell2location succeeded in the most situations with the best performance.

Several additional features were described for some methods in their related publications, and these were assessed for our guidelines for method selection. SpatialDecon and RCTD are claimed to correct the variance in gene expression profiles, which resolves the platform effects between scRNA-seq and spatial transcriptomics data ${ }^{11,14}$. Cell2location and SpiceMix can identify potential fine-grained sub-cell types, and CARD can impute cell-type compositions to construct a refined spatial tissue map ${ }^{9,10,22}$. Most methods utilize discrete cell types, although DestVI has the advantage of identifying the continuous variation within the same cell type, which is useful for studying the same tissue section under different conditions ${ }^{13}$.

## Discussion

Here, we presented a comprehensive benchmarking study and guidelines for 18 existing cellular deconvolution methods used in
spatial transcriptomics. We evaluated these methods using 50 datasets with multiple metrics in terms of accuracy, robustness, and usability. The datasets we used included simulated datasets binned by single-cell resolution datasets (e.g., seqFISH+ and MERFISH) which could be used for quantitative evaluation as the ground-truth is known for such datasets, and real-world datasets generated by sequencing-based technologies (e.g., Slide-seq V2 and stereo-seq) which could be used to mimic the real-world scenarios for cellular deconvolution tasks. Considering the performance and additional features of the methods, decision-tree-based scenario-specific recommendations and guidelines for method selection were proposed for users. We found that the performance of the 18 methods varied among multiple spatial transcriptomics technologies with different experimental conditions. Nevertheless, each method category contained at least one highperforming method. In general, CARD, Cell2location, Tangram, and RCTD were the best performing methods. Compared with the existing benchmarking studies ${ }^{8,9}$, our study included most number of existing methods. More importantly, we provided a full-scale summarization of the performance of all the methods and characterized it as a clear guideline including the solid recommendation of the methods and demonstration of their additional features which would give readers an overall understanding of deconvolution in spatial transcriptomics data analyses.

Following our assessment, we raise two general but crucial limitations that await solutions. First, the platform effect causes two problems as follows: (1) systematic variation in gene expression profiles between scRNA-seq and spatial transcriptomics data; owing to differences in technology-dependent library preparation and sequencing platforms, discrepancies in the detected mRNAs from the same tissue section are inevitable, especially in heterogeneous cancer tissue sections, and (2) variation between two modalities affects the mismatch of cell types between scRNA-seq and spatial transcriptomics data. The prior assumption of integrating these two modalities is to share the same cell types between scRNA-seq and spatial transcriptomics data. Even though RCTD solves the platform effect issue via a normalization strategy among all cell types, cell-type-specific platform effects warrant further exploration ${ }^{12}$. Second, the high dropout rate of spatial transcriptomics is a traditional issue in scRNAseq. Small libraries lead to deficient mRNA detection; thus, marker
genes for rare cell types become undetectable ${ }^{16}$. The biological pipeline of spatial transcriptomics technologies should be improved further, even though this situation has been considered in $\mathrm{SD}^{2}$ and an imputation method has been proposed ${ }^{21,37}$.

We also present possible future directions of the field to shed light on the development in the field. (1) Multimodal learning will likely become a hotspot in the development of cellular deconvolution methods used in spatial transcriptomics and its applications. For instance, bioinformaticians could use histological tissue images with image intensity levels that could improve our understanding of spatial transcriptomics. (2) Three-dimensional deconvolution and mapping of tissues will provide more novel biological insights than are currently provided by two-dimensional deconvolution. More spatial transcriptomics datasets with consecutive tissue slices are emerging, and the spatial context of interslices will provide more patterns that assist the deconvolution process ${ }^{11}$. (3) Through the development of spatial transcriptomics technologies, the resolution of spots becomes higher, up to the subcellular resolution by recent technologies ${ }^{34}$. Although the rapid progress of spatial transcriptomics technologies is exciting, the marginal benefits of higher resolution are outweighed by the booming higher dropout issue which has more risks of losing some valuable information. In the future, under current sub-resolution technologies, the spatial transcriptomics technology for the single biological cells should be more essential to develop. (4) The cellular deconvolution of spatial transcriptomics will not only help biologists study the structure of tissues but also become associated with artificial intelligence-assisted computational pathology and the healthcare system.

## Methods

## Benchmark metrics

We first assume that there are $J$ genes per spot and $I$ captured spots in the whole spatial transcriptomics data. $X_{i j}$ represents the expression value of gene $j$ in the $i$ th spot. $T_{i k}$ and $P_{i k}$ represent the true and predicted proportion of cell type $k$, respectively, in the $i$ th spot through the number of total cell types $K$. To evaluate the performance of the tested methods comprehensively, the main benchmark metrics used were RMSE, JSD, and PCC. The definitions of these metrics, as used in our benchmarking pipeline, are provided below.

1. RMSE was calculated between $T_{i k}$ and $P_{i k}$ of per cell type, normalize them by the sum of calculated proportions among all spots $S_{k}$, and then average them as the final RMSE as the following equation:

$$
\mathrm{RMSE}=\sqrt{\frac{1}{K} \sum_{k=1}^{K} \frac{1}{S_{k}} \sum_{i=1}^{l}\left(P_{i k}-T_{i k}\right)^{2}}
$$

2. JSD was calculated as a score between $T_{k}$ and $P_{k}$ per cell type in all spots. The conception of Kullback-Leibler divergence (KL) is used for calculating JSD. $Q\left(P_{k}\right)$ and $Q\left(T_{k}\right)$ represent algorithmpredicted and true distribution of cell type $k$. We average them as the final JSD as the following equation:

$$
\mathrm{JSD}=\frac{1}{2} \mathrm{KL}\left(Q\left(T_{k}\right) \| \frac{Q\left(P_{k}\right)+Q\left(T_{k}\right)}{2}\right)+\frac{1}{2} \mathrm{KL}\left(Q\left(P_{k}\right) \| \frac{Q\left(P_{k}\right)+Q\left(T_{k}\right)}{2}\right)
$$

$$
\mathrm{KL}\left(Q\left(P_{k}\right) \| Q\left(T_{k}\right)\right)=\sum Q\left(P_{k}\right) \ln \frac{Q\left(P_{k}\right)}{Q\left(T_{k}\right)}
$$

3. Because ground truth does not exist in sequencing-based spatial transcriptomics data, we calculated the PCC between the
predicted proportion of specific cell type $P_{k}$ and the expression profile of its marker gene $E_{g}$ using the following equation:

$$
\operatorname{PCC}\left(P_{k}, E_{g}\right)=\frac{\mathbb{E}\left[P_{k} E_{g}\right]-\mathbb{E}\left[P_{k}\right] \mathbb{E}\left[E_{g}\right]}{\sqrt{\mathbb{E}\left[P_{k}{ }^{2}\right]-\left(\mathbb{E}\left[P_{k}\right]\right)^{2}} \sqrt{\mathbb{E}\left[E_{g}{ }^{2}\right]-\left(\mathbb{E}\left[E_{g}\right]\right)^{2}}}
$$

## Evaluation of methods without a scRNA-seq reference

The methods used without a scRNA-seq reference, also called unsupervised methods, deconvolved low-resolution spots based on the gene expression profile and location of spots from spatial transcriptomics data only. STdeconvolve was inspired by the notion of discovering latent topics in collections of documents, which is a common task in natural language processing, and uses Latent Dirichlet Allocation to infer the proportions of cell types based on gene expression profiles in spatial transcriptomics. Berglund uses Poisson factor analysis and Monte-Carlo Markov Chain sampling to deconvolve spots. SpiceMix uses the locations of spots as an extra input in addition to gene expression profiles and incorporates graph representations of spatial relationships into matrix factorization to deconvolve spots.

In general, unsupervised methods needed user to specify the number of topics which represented the clusters waiting for assigning the names of known cell types. Ideally, the number of topics should be similar to the number of cell types. In our evaluation, we set the number of topics as the true number of cell types manually among all datasets. After the deconvolution, a topic-by-spot matrix was generated, and we multiplied this matrix and inputted a spot-by-gene matrix to obtain the topic-by-gene matrix. We aimed to map the topics to real cell types for further evaluation. First, we summed up the same cell types in an annotated cell-by-gene matrix from the scRNA-seq data as the cell-type-by-gene matrix. For each real cell type from scRNA-seq data, we calculated the PCC between this cell type and all topics, chose the best-paired topic with the highest PCC, and assigned the name of current cell type to chosen topic. After assignment, this chosen topic would be ignored in the future steps. Then, we repeated the aforementioned steps on the next cell type until all cell types were iterated. For now, each topic should be paired with the best suitable cell type without duplication and topic-by-spot matrix could be easily transferred to cell-type-by-spot matrix for evaluating the performance of unsupervised methods further.

## Preprocessing of datasets

Owing to the extremely high resolution ( 500 nm per spot) and dropout rate of stereo-seq data, it was necessary to integrate subcellularresolution spots into low-resolution spots by binning them using a $100 \times 100$ spot square (bin100) with slides of $50 \mu \mathrm{~m}(100 \times 500 \mathrm{~nm})$ and summing their gene expression profiles. The bin100 stereo-seq data had a similar resolution as that of Visium, and it performed reasonably in deconvolution tasks. For the zebrafish embryo dataset by stereo-seq, it was binned by 5,10 and $15 \mu \mathrm{~m}$ to test the robustness.

To evaluate the robustness of the methods, we tested their performance in terms of different cell types. We integrated the 17 original cell types into 11 cell types, thereby combining some sub-cell types. We integrated CA1 and CA3 into "Cornu Ammonis"; Neuron.Slc17a6, Neurogenesis, and Cajal_Retzius into "Neuron"; Endothelial_Stalk and Endothelial_Tip into "Endothelial"; and Oligodendrocyte and Polydendrocyte into "Oligo_Poly" ${ }^{12}$.

For marker gene selection in all datasets, most biological marker genes were chosen from publications. Specifically, we chose the topfive highly variable genes (calculating the fold-change of each gene) for each specific cell type from Slide-seqV2 datasets as the marker genes.

## Construction of a summary table

We constructed a summary table to show the performances of the methods (Fig. 2). Because JSD, RMSE, and running time showed better performance with lower values, we normalized the value $x$ in each column according to minmax $\left(\max \left(x_{\text {col }}\right)-x\right)$, where $x_{\text {col }}$ represents the vector of the column. For the other metrics, we normalized according to minmax directly. Thus, we unified a pattern in which darker dots represent better performance, which is a pattern that users will find easy to identify.

## Implementation of methods

CARD ${ }^{10}$ : we used the code of CARD v1.0.0 from https://github.com/ YingMa0107/CARD. We set minCountGene to 5 and minCountSpot to 5, which are the default parameter settings.

SPOTlight ${ }^{19}$ : we used the code of SPOTlight v0.99.0 from https:// github.com/MarcElosua/SPOTlight. We set cl_n to 10 and hvg to 2000.

DSTG ${ }^{20}$ : we used the code of DSTG from https://github.com/Su-informatics-lab/DSTG. We set learning_rate to 0.01 and epoch to 300 .

SpatialDWLS ${ }^{18}$ : we used the code of SpatialDWLS from https:// github.com/RubD/Giotto/, which integrates the SpatialDWLS method. We set min_genes in findMarkers_one_vs_all to 20 and gene_det_in_min_cells and min_det_genes_per_cell in filterGiotto to 5 and 5, respectively.
$\mathrm{SD}^{221}$ : we used the code of $\mathrm{SD}^{2}$ from https://github.com/ leihouyeung/SD2, with the following settings: spot_num $=1000$, lower_cellnum $=2$, and upper_cellnum $=20$.

NMFreg ${ }^{38}$ : we used the code of NMFreg from https://github.com/ tudaga/NMFreg_tutorial. The NMF function was used with the following parameters: number of components $=30$, random_state $=17$, and init = random.

Stereoscope ${ }^{14}$ : we used the code of stereoscope v. 03 from https:// github.com/almaan/stereoscope. Analysis with stereoscope was conducted on a GPU with the following parameters: number of genes = 5000 , st epochs $=75,000$, st batch size $=1000$, sc epochs $=75,000$, sc batch size $=1000$, and learning rate $=0.01$.

Tangram ${ }^{22}$ : we used the code of Tangram v1.0.3 from https:// github.com/broadinstitute/Tangram. The mapping of cells to space was conducted with the function tg.map_cell_to_space with mode = clusters.

Cell2location ${ }^{11}$ : we used the code of Cell2location v0.1 from https://github.com/BayraktarLab/cell2location. The settings maxepochs $=4000$, batch_size $=$ None, and train_size $=1$ were used.

STdeconvolve ${ }^{25}$ : we used the code of STdeconvolve 1.0.0 from https://github.com/JEFworks-Lab/STdeconvolve. We used the default settings, except that the number of factors was set correctly according to each dataset.

Berglund ${ }^{23}$ : we used the code of Berglund 0.2.0 from https:// github.com/SpatialTranscriptomicsResearch/std-poisson. We set the following parameters: -iter $(=2000)$, -feature_alpha $\arg =1$, --mix_alpha arg (= 0.5), --phi_r_1 arg (=1), --phi_r_2 arg (= 0.001), --phi_p_1 arg $(=2)$, --phi_p_2 arg $(=2)$, --theta_r_1 arg $(=1)$, --theta_r_2 arg $(=1)$, --theta_p_1 arg $(=0.050000),-$-theta_p_2 arg $(=0.950000),-$-spot_1 arg $(=10)$, --spot_2 arg (=10), --sigma arg (=1), --residual arg (=100), --bline1 $\arg (=50)$, and --bline2 $\arg (=50)$ The MCMC inference options were set as follows: --MHiter arg (=100) and --MHtemp arg (=1).

SpiceMix ${ }^{24}$ : we used the code of SpiceMix from https://github. com/ma-compbio/SpiceMIx. We selected the number of factors based on the used dataset and set use spatial to True.

RCTD ${ }^{12}$ : we used the code of RCTD from https://github.com/ dmcable/spacexr, which is integrated into a tool called spacexr (2.0.0). Spacexr (RCTD) was run with following the configuration: (1) create.RCTD was used with the parameter CELL_MIN_INSTANCE = 1; (2) run.RCTD was used in the doublet mode.

SpatialDecon ${ }^{15}$ : we used the code of SpatialDecon from https:// github.com/Nanostring-Biostats/SpatialDecon.git. SpatialDecon was run with the expected background count bg set to 0.01 .

STRIDE ${ }^{16}$ : we used the code of STRIDE from https://github.com/ DongqingSun96/STRIDE. The cell-type-associated topic profiles were obtained using the "STRIDE deconvolve" function. If not specified, STRIDE set the $75 \%$ quantile of nCount as the default scaling factor.

DestVI ${ }^{13}$ : we used the code of DestVI (scvi-tools 0.16.0) from https://github.com/scverse/scvi-tools. DestVI first required genes with $<10$ counts to be filtered using the function "sc.pp.filter_genes." To perform deconvolution, the single-cell model was then trained to learn the basis of gene expression with the scRNA-seq data for 300 epochs, whereas the spatial model was trained for 2500 epochs, with a learning rate of $10^{-3}$.

SpaOTsc ${ }^{26}$ : we conducted the code of SpaOTsc from https:// github.com/zcang/SpaOTsc.
novoSpaRc ${ }^{27}$ : we conducted the code of novoSpaRc 0.4.4 from https://github.com/rajewsky-lab/novosparc. Alpha was set as 0.5.

## Computational resource

The workstation used to test all methods had a $2 \operatorname{Intel}(\mathrm{R}) \operatorname{Xeon}(\mathrm{R})$ CPU E5-2680 v3 @ 2.50 GHz ( $30,720 \mathrm{~KB}$ cache size; 24 cores in total) and 528 GB of memory. The GPUs were two Nvidia Quadro M6000 24 GB ( 48 GB in total). The operating system used was Ubuntu 18.04.

## Reporting summary

Further information on research design is available in the Nature Portfolio Reporting Summary linked to this article.

## Data availability

A summary of the data is shown in Supplementary Table 1. seqFISH+: scRNA-seq and spatial transcriptomics data were both obtained from https://github.com/CaiGroup/seqFISH-PLUS. MERFISH: scRNA-seq data were obtained from https://github.com/rdong08/spatialDWLS_ dataset/tree/main/datasets, and spatial transcriptomics data were obtained from https://datadryad.org/stash/dataset/doi:10.5061/dryad. 8t8s248/. ST: scRNA-seq and spatial transcriptomics data were both obtained from GSE111672. Visium: scRNA-seq data and spatial transcriptomics data were obtained from https://github.com/ BayraktarLab/cell2location. Slide-seqV2: scRNA-seq and spatial transcriptomics data were both obtained from https://github.com/ dmcable/spacexr. Stereo-seq (olfactory bulb): scRNA-seq data were obtained from GSE71585, and spatial transcriptomics data were obtained from https://db.cngb.org/stomics/mosta/. Stereo-seq (zebrafish embryo): the scRNA-seq data and spatial transcriptomics data are from https://db.cngb.org/stomics/datasets/STDS0000057. All the mentioned datasets were also integrated and uploaded to public repository ${ }^{39}$. All other relevant data supporting the key findings of this study are available within the article and its Supplementary Information files or from the corresponding author upon reasonable request. No data were excluded from the analyses; the experiments were not randomized; the Investigators were not blinded to allocation during experiments and outcome assessment. Source data are provided with this paper.

## Code availability

The code is available at https://github.com/leihouyeung/STdeconv_ benchmark ${ }^{39}$.

## References

1. Method of the Year 2020: spatially resolved transcriptomics. Nat. Methods 18, 1 (2021).
2. Williams, C. G., Lee, H. J., Asatsuma, T., Vento-Tormo, R. \& Haque, A. An introduction to spatial transcriptomics for biomedical research. Genome Med. 14, 68 (2022).
3. Petukhov, V. et al. Cell segmentation in imaging-based spatial transcriptomics. Nat. Biotechnol. https://doi.org/10.1038/s41587-021-01044-w (2021).
4. Moses, L. \& Pachter, L. Museum of spatial transcriptomics. Nat. Methods 19, 534-546 (2022).
5. Lein, E., Borm, L. E. \& Linnarsson, S. The promise of spatial transcriptomics for neuroscience in the era of molecular cell typing. Science 358, 64-69 (2017).
6. Hildebrandt, F. et al. Spatial transcriptomics to define transcriptional patterns of zonation and structural components in the mouse liver. Nat. Commun. 12, 7046 (2021).
7. Ji, A. L. et al. Multimodal analysis of composition and spatial architecture in human squamous cell carcinoma. Cell 182, 497-514 (2020).
8. Li, B. et al. Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution. Nat. Methods https://doi.org/10.1038/s41592-022-01480-9 (2022).
9. Chen, J. et al. A comprehensive comparison on cell-type composition inference for spatial transcriptomics data. Brief. Bioinform. 23, bbac245 (2022).
10. Ma, Y. \& Zhou, X. Spatially informed cell-type deconvolution for spatial transcriptomics. Nat. Biotechnol. https://doi.org/10.1038/ s41587-022-01273-7 (2022).
11. Kleshchevnikov, V. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat. Biotechnol. https://doi.org/10.1038/ s41587-021-01139-4 (2022).
12. Cable, D. M. et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nat. Biotechnol. https://doi.org/10.1038/ s41587-021-00830-w (2021).
13. Lopez, R. et al. DestVI identifies continuums of cell types in spatial transcriptomics data. Nat. Biotechnol. https://doi.org/10.1038/ s41587-022-01272-8 (2022).
14. Andersson, A. et al. Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography. Commun. Biol. 3, 565 (2020).
15. Danaher, P. et al. Advances in mixed cell deconvolution enable quantification of cell types in spatial transcriptomic data. Nat. Commun. 13, 385 (2022).
16. Sun, D., Liu, Z., Li, T., Wu, Q. \& Wang, C. STRIDE: accurately decomposing and integrating spatial transcriptomics using singlecell RNA sequencing. Nucleic Acids Res. 50, e42 (2022).
17. Rodriques, S. G. et al. Slide-seq: a scalable technology for measuring genome-wide expression at high spatial resolution. Science 363, 1463-1467 (2019).
18. Dong, R. \& Yuan, G.-C. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biol. 22, 145 (2021).
19. Elosua-Bayes, M., Nieto, P., Mereu, E., Gut, I. \& Heyn, H. SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Res. 49, e50 (2021).
20. Song, Q. \& Su, J. DSTG: deconvoluting spatial transcriptomics data through graph-based artificial intelligence. Brief. Bioinform. 22, bbaa414 (2021).
21. Li, H., Li, H., Zhou, J. \& Gao, X. SD2: Spatially resolved transcriptomics deconvolution through integration of dropout and spatial information. Bioinformatics btac605 https://doi.org/10. 1093/bioinformatics/btac605 (2022).
22. Biancalani, T. et al. Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nat. Methods 18, 1352-1362 (2021).
23. Berglund, E. et al. Spatial maps of prostate cancer transcriptomes reveal an unexplored landscape of heterogeneity. Nat. Commun. 9, 2419 (2018).
24. Chidester, B., Zhou, T., Alam, S. \& Ma, J. SpiceMix enables integrative single-cell spatial modeling of cell identity. Nat. Genet. 55, 78-88 (2023).
25. Miller, B. F., Huang, F., Atta, L., Sahoo, A. \& Fan, J. Reference-free cell type deconvolution of multi-cellular pixel-resolution spatially resolved transcriptomics data. Nat. Commun. 13, 2339 (2022).
26. Cang, Z. \& Nie, Q. Inferring spatial and signaling relationships between cells from single cell transcriptomic data. Nat. Commun. 11, 2084 (2020).
27. Moriel, N. et al. NovoSpaRc: flexible spatial reconstruction of singlecell gene expression with optimal transport. Nat. Protoc. 16, 4177-4200 (2021).
28. Eng, C.-H. L. et al. Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH. Nature 568, 235-239 (2019).
29. Xia, C., Fan, J., Emanuel, G., Hao, J. \& Zhuang, X. Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression. Proc. Natl Acad. Sci. USA 116, 19490-19499 (2019).
30. Moffitt, J. R. et al. Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. Science 362, eaau5324 (2018).
31. Ståhl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353, 78-82 (2016).
32. Hajdarovic, K. H. et al. Single-cell analysis of the aging female mouse hypothalamus. Nat. Aging 2, 662-678 (2022).
33. Stickels, R. R. et al. Highly sensitive spatial transcriptomics at nearcellular resolution with Slide-seqV2. Nat. Biotechnol. 39, 313-319 (2021).
34. Chen, A. et al. Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays. Cell 185, 1777-1792 (2022).
35. Liu, C. et al. Spatiotemporal mapping of gene expression landscapes and developmental trajectories during zebrafish embryogenesis. Dev. Cell 57, 1284-1298 (2022).
36. Wolf, F. A., Angerer, P. \& Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 19, 15 (2018).
37. Li, Z., Song, T., Yong, J. \& Kuang, R. Imputation of spatially-resolved transcriptomes by graph-regularized tensor completion. PLoS Comput. Biol. 17, e1008218 (2021).
38. Rodriques, S. G. et al. Slide-seq: a scalable technology for measuring genome-wide expression at high spatial resolution. Science 363, 1463 LP-1461467 (2019).
39. Li, H. A comprehensive benchmarking with practical guidelines for cellular deconvolution of spatial transcriptomics. https://doi.org/ 10.5281/zenodo. 7674290 (2023).

## Acknowledgements

This publication is based upon work supported by the King Abdullah University of Science and Technology (KAUST) Office of Research Administration (ORA) under Award Nos. FCC/1/1976-44-01, FCC/1/1976-45-01, URF/1/4663-01-01, REI/1/5202-01-01, REI/1/5234-01-01, REI/1/ 4940-01-01, and RGC/3/4816-01-01. We thank Hanmin Li for helping running some experiments.

## Author contributions

X.G. and H.L. conceived and initiated this study. X.L., B.Z., R.Z. and Y.W. prepared all spatial transcriptomics datasets and scRNAseq datasets. H.L., J.Z., Z.L. and S.C. conducted all the experiments of all methods under accuracy, robustness and usability. H.L. and J.Z. outputted the figure and tables. H.L. wrote the manuscript under supervision of X.G. S.S. polished the writing of manuscript. All authors are involved in discussion and finalization of the manuscript.

## Competing interests

The authors declare no competing interests.

## Additional information

## Supplementary information The online version contains

supplementary material available at
https://doi.org/10.1038/s41467-023-37168-7.

Correspondence and requests for materials should be addressed to Xin Gao.

Peer review information Nature Communications thanks the anonymous reviewer(s) for their contribution to the peer review of this work. Peer reviewer reports are available.

Reprints and permissions information is available at
http://www.nature.com/reprints

Publisher's note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.

Open Access This article is licensed under a Creative Commons Attribution 4.0 International License, which permits use, sharing, adaptation, distribution and reproduction in any medium or format, as long as you give appropriate credit to the original author(s) and the source, provide a link to the Creative Commons license, and indicate if changes were made. The images or other third party material in this article are included in the article's Creative Commons license, unless indicated otherwise in a credit line to the material. If material is not included in the article's Creative Commons license and your intended use is not permitted by statutory regulation or exceeds the permitted use, you will need to obtain permission directly from the copyright holder. To view a copy of this license, visit http://creativecommons.org/ licenses/by/4.0/.
(c) The Author(s) 2023


[^0]:    ¹Computational Bioscience Research Center, King Abdullah University of Science and Technology (KAUST), Thuwal, Saudi Arabia. ${ }^{2}$ Computer, Electrical and Mathematical Sciences and Engineering Division, King Abdullah University of Science and Technology (KAUST), Thuwal, Saudi Arabia. ${ }^{3}$ Syneron Technology, Guangzhou 510000, China. ${ }^{4}$ Key Lab of Intelligent Information Processing, Institute of Computing Technology, Chinese Academy of Sciences, 100190 Beijing, China. ${ }^{5}$ University of Chinese Academy of Sciences, 100049 Beijing, China. ${ }^{6}$ These authors contributed equally: Haoyang Li, Juexiao Zhou. -e-mail: xin.gao@kaust.edu.sa


---

# Paper 2: s41467-023-43600-9

# Spatial transcriptomics deconvolution at single-cell resolution using Redeconve 

Received: 26 October 2022
Accepted: 14 November 2023
Published online: 01 December 2023

Zixiang Zhou ${ }^{(1)}{ }^{1,2,3}$, Yunshan Zhong ${ }^{(1,3}$, Zemin Zhang ${ }^{(1,2)}$ \& Xianwen Ren ${ }^{(1)}$ ®

Computational deconvolution with single-cell RNA sequencing data as reference is pivotal to interpreting spatial transcriptomics data, but the current methods are limited to cell-type resolution. Here we present Redeconve, an algorithm to deconvolute spatial transcriptomics data at single-cell resolution, enabling interpretation of spatial transcriptomics data with thousands of nuanced cell states. We benchmark Redeconve with the state-of-the-art algorithms on diverse spatial transcriptomics platforms and datasets and demonstrate the superiority of Redeconve in terms of accuracy, resolution, robustness, and speed. Application to a human pancreatic cancer dataset reveals cancer-clone-specific T cell infiltration, and application to lymph node samples identifies differential cytotoxic T cells between $\operatorname{Ig} \mathrm{A}+$ and $\operatorname{IgG}+$ spots, providing novel insights into tumor immunology and the regulatory mechanisms underlying antibody class switch.

Spatial transcriptomics (ST) technologies provide new tools to identify the cellular organization and interactions of biological samples, which is pivotal to biomedical studies. Multiple ST technologies have been developed and applied to mouse and human brains, lymph node, heart, etc., providing novel insights into cellular communication networks underlying different conditions. However, sequencing-based ST technologies, e.g., the 10x Genomics Visium platform and Slide-seq ${ }^{1}$, are essentially of a spot-by-gene matrix structure, needing additional data to provide the cellular information. While the commercial emergence of imaging-based ST technologies, e.g., seqFISH $+^{2}$, MERFISH ${ }^{3}$, 10x Genomics Xenium ${ }^{4}$, and NanoString CosMx ${ }^{5}$, provides subcellular resolution, these technologies are limited by low gene throughput, with hundreds of customized genes detected, making their discovery potential unparallel to whole transcriptome-wide spatial technologies. Therefore, integrative analysis of whole transcriptome-wide ST data together with matched single-cell RNA sequencing (scRNA-seq) data is of high significance for biological discoveries.

Multiple effective and efficient algorithms have been proposed for integrative analysis of whole-transcriptome ST and scRNA-seq data. The current algorithms can be categorized to two groups: (1) mappingbased methods, e.g., NovospaRc ${ }^{6}$, Tangram ${ }^{7}$, Celltrek ${ }^{8}$, and CytoSPACE ${ }^{9}$, which map single cells to the positions of ST data according to gene expression similarity or related measures; and (2)
deconvolution-based methods, e.g., CARD ${ }^{10}$, RCTD ${ }^{11}$, cell2location ${ }^{12}$, DestVI ${ }^{13}$, SpatialDWLS ${ }^{14}$, SPOTlight ${ }^{15}$, STRIDE ${ }^{16}$, CellDART ${ }^{17}$, Celloscope ${ }^{18}$, DSTG ${ }^{19}$, and Stereoscope ${ }^{20}$, which try to reconstruct the ST observations by modeling the experimental process as sampling from different combinations of single cells. Mapping-based methods are superior to the current deconvolution-based methods regarding their single-cell resolution as the resolution of current deconvolution methods is limited to tens of cell types. However, mapping-based methods may introduce artificial biases during the mapping process due to the absence of strong constraint on the reconstruction accuracy of the ST observations. It is urgently needed to develop a deconvolution-based algorithm with single-cell resolution to fully release the biological information hidden in ST data.

In this study, we develop an algorithm, named as Redeconve ${ }^{21}$, to estimate the cellular composition of ST spots. Different from previous deconvolution-based algorithms, Redeconve introduces a regularizing term to solve the collinearity problem of high-resolution deconvolution, with the assumption that similar single cell states have similar abundance in ST spots. This algorithmic innovation not only improves the deconvolution resolution from tens of cell types to thousands of single cell states, but also greatly improve the reconstruction accuracy of ST data, enabling illustration of the nuanced biological mechanisms hidden in the ST data. Stringent comparison with the state-of-the-art

[^0]algorithms including cell2location, CARD, DestVI, CellTrek, NovoSpaRc, and Tangram demonstrates the superiority of Redeconve in terms of reconstruction accuracy, cell abundance estimation per spot, sparseness of the reconstructed cellular composition, cell state resolution, and computational speed. Application to human pancreatic cancer data reveals novel insights into tumor-infiltrating CD8 + T cells, and application to human lymph node data reveals new clues for the regulatory factors of $\operatorname{Ig} \mathrm{A}+$ and $\operatorname{IgG}+\mathrm{B}$ cells.

## Results

## Redeconve: a quadratic programming model for single-cell deconvolution of ST data

Redeconve uses scRNA-seq or single-nucleus RNA-seq (snRNA-seq) as reference to estimate the abundance of different cell states in each spot of ST data (Fig.1a). Different from previous deconvolution methods, Redeconve does not need to group single cells into clusters and then do deconvolution. Instead, Redeconve treats each cell of the sc/snRNA-seq data as a specific cell state serving as reference to estimate the cellular composition of ST data. The direct usage of sc/ snRNA-seq data as reference is conceptually direct and computationally efficient, with the potential to handle the heterogeneity of ST data. However, direct usage of sc/snRNA-seq data as reference will introduce a new challenge, i.e., collinearity. That is, multiple single cells have similar profiles of gene expression, prohibiting the accurate estimation of the abundance of individual cell states. We introduce a biologically reasonable heuristic by assuming that similar cells have similar abundance within ST spots, and thus mathematically introduce a regularization term in the deconvolution model based on nonnegative least regression. Solving this regularized deconvolution model by quadratic programming will produce robust estimation of the cellular composition at single-cell resolution for each spot of ST data.

## High accuracy, resolution, robustness, efficacy, and scalability of Redeconve

We applied Redeconve to multiple ST datasets from various platforms ( 10 x Visium, Slide-seq v2, ST, etc.) and compared the performance with other methods. We first compared the consistency of results among different methods at the cell-type resolution based on a human breast cancer dataset. The results suggested that deconvolution-based methods including Redeconve had higher consistency with each other than mapping-based methods (Fig. 1b), indicating the relative superiority and robustness of deconvolution-based methods. This observation is confirmed on additional ST datasets (Supplementary Fig. 1). Different from previous deconvolution-based methods which only reported cell-type-level results, Redeconve can further dictate fine-grained cell states at single-cell resolution (Fig. 1c and Supplementary Fig. 2). On a ST dataset from a human breast cancer sample, Redeconve resolved 249 different cell states from 9 major cell types (Fig. 1c). On a ST dataset from mouse cerebellum, Redeconve resolved 1000 different cell states from 14 major cell types (Fig. 1c). In contrast, the resolution of previous deconvolution methods is limited by the clustering results of sc/snRNA-seq data.

In addition to the robustness and resolution superiority, Redeconve also improves the reconstruction accuracy of gene expression per spot, and the improvement is independent on similarity measures such as cosine similarity, Pearson's correlation, and Root Mean Square Error (RMSE) between the true ST gene expression profile and the reconstructed gene expression vector (Fig. 1d, and Supplementary Figs. 3-4). Redeconve also reached high accuracy of estimated cell abundance (based on a ground truth by nucleus counting, Fig. 1e and Supplementary Fig. 5), and superior computational speed (Fig. 1f and Supplementary Fig. 6). When suitable reference is provided, e.g., matched scRNA-seq data, Redeconve can reach $>0.8$ cosine accuracy for most ST spots (Fig. 1d). With no suitable reference available (for
example, only snRNA-seq data are accessible for brain samples), Redeconve still outperforms other methods (Fig. 1d). Pairwise comparison between Redeconve and other methods further shows the superiority of Redeconve on almost all spots regarding the reconstruction accuracy (Supplementary Figs. 7-12). Because Redeconve conducts deconvolution analysis spot by spot, parallel computation is enabled and thus Redeconve demonstrates superior computation speed compared with current deconvolution algorithms (Fig. 1f and Supplementary Fig. 6).

To evaluate the performance of Redeconve in estimating the absolute abundance of cells within ST spots, we applied Redeconve to three datasets: Mouse Brain, PDAC and Human Breast Cancer Xenium, in which the cell counts were obtained by nucleus counting based on image segmentation ${ }^{12,22,23}$. Without any priori information, the results of Redeconve showed high conformity with the "ground-truth" cell counts (Fig. 1e), similar to those methods with cell counts (or cell density) as priori knowledge e.g., cell2location and Tangram (Supplementary Fig. 5). We used Shannon entropy to estimate the potential number of different cell states within each spatial spot (see Methods for details about using perplexity as a metric). Redeconve revealed high spot heterogeneity by showing that some spots had complex cellular composition while others had a relatively simple one. In contrast, the entropy of other methods is uniformly high, showing that each spot had been composed of almost all the cell types in reference, which is unrealistic (Supplementary Fig. 13).

## Single-cell resolution is unique to Redeconve compared with previous deconvolution algorithms

Then we examined whether the current deconvolution-based algorithms could be upgraded to single-cell resolution by switching the required cell types to thousands of single cells as Redeconve does. Among all the methods we evaluated, only cell2location and DestVI completed the task but took a rather long time compared with the celltype inputs (Supplementary Fig. 14) while other algorithms reported errors. Although single-cell inputs improved the reconstruction accuracy of cell2location on the ST data of a human lymph node sample based on the 10x Genomics Visium platform, cell2location did not reach improvement on the human pancreatic tumor and mouse brain datasets, and DestVI failed on all three evaluations (Supplement Fig. 15). In contrast, Redeconve outperformed cell2location and DestVI on almost all spots of the evaluated datasets (Fig. 2a). When switching the inputs from cell types to single cells, DestVI achieved well sparsity regarding the different cell states within each spot (measured by perplexity according to Shannon entropy), similar to the performance of Redeconve. But cell2location reported extremely high perplexity for most spots, indicating overpredicted presence of almost all cell types and thus high false positive rate (Fig. 2b). Therefore, changing inputs from cell types to single cells cannot upgrade the performance of current algorithms to levels parallel to that of Redeconve, and the superiority of Redeconve analysis is mainly derived from algorithmic innovation.

## Evaluating the impact of cell-type resolution on deconvolution by simulation

To evaluate how the cell-type resolution of reference data impacts the deconvolution analysis, we devised a series of simulation experiments to showcase the performance differences of Redeconve and the state-of-the-art algorithms. We constructed three pseudo-bulk RNA-seq datasets by averaging the gene expression data of individual cells based on scRNA-seq data from the PDAC ${ }^{24}$, human lymph node ${ }^{12,25}$ and human testis ${ }^{26}$ datasets separately (Fig. 2c and Methods). Then we applied Redeconve and cell2location, the only alternative method capable of this task. With direct comparison with ground-truth, the results indicate that Redeconve performs substantially better than cell2location, as evidenced by its significantly higher accuracy
![](https://cdn.mathpix.com/cropped/2025_11_30_57cdc3a614a421a1b5dfg-03.jpg?height=2104&width=1514&top_left_y=158&top_left_x=285)
(Supplementary Figs. 16-18). When examining the relationship between accuracy and number of clusters in single-cell reference, Redeconve showed an increase in accuracy when the number of clusters grows, while cell2location experienced a sharp drop (Fig. 2d). This suggests that Redeconve is capable of handling large-scale scRNAseq data more effectively and can use finer-grained clusters to increase accuracy instead of becoming confused. Furthermore, simulation
experiments also corroborate the validity of using perplexity as a metric of sparsity (Supplementary Table 1 and Methods).

## Evaluating the estimating accuracy of cell-type proportion by 10x Genomics Xenium data as ground truth

Single-cell ST platforms, such as MERFISH ${ }^{3}$, Xenium ${ }^{4}$ and CoxMx ${ }^{5}$, are commercially emerging as a powerful tool for the high-resolution

Fig. 1 | Overview of the Redeconve algorithm and benchmark analysis.
a overview of Redeconve workflow for deconvoluting spatial transcriptomics data. Redeconve requires sc/snRNA-seq data together with spatial transcriptomics data as input and performs deconvolution by solving a regularized non-negative least regression model with the aims to estimate cellular composition across spots at single-cell resolution. b heatmap illustrated median spot-level Spearman's correlation of cell type proportions among different algorithms on a human breast cancer dataset. c Sankey diagram demonstrated the cell-type and single-cell resolutions of Redeconve results on human breast cancer and mouse cerebellum datasets, respectively. The bar height of cell types or single cells refer to their estimated abundance after deconvolution. d line chart of cosine similarities
between observed and reconstructed expression profiles per spot based on six ST datasets. $N=4039,2426,428,36550,2987$ and 39431 spots for human lymph nodes, human breast cancer, PDAC (pancreatic ductal adenocarcinoma), human testis, mouse brain and mouse cerebellum respaectively. Spots were sorted by an ascending order of the cosine similarities. e Pearson correlation of cell abundances between Redeconve and the cell counts per spot based on a mouse brain dataset. The ground truth cell counts per spot was obtained by nucleus counting of cell segmentation image ${ }^{12}$. $\mathbf{f}$ computational efficiency of different deconvolution-based and mapping-based algorithms on a human lymph nodes dataset. Source data of 1ce are provided as a Source Data file.
mapping of the precise location of single cells, but are limited by the number of genes profiled during experiments because customized probes specific to target genes need to be designed and synthesized before experiments. The high resolution of these platforms provides natural ground truth to evaluate the performance of Redeconve. Here, we used a human breast cancer Xenium dataset generated by 10x Genomics ${ }^{4}$ to evaluate the performance of Redeconve regarding reconstruction of ST spot expression profiles, cell type proportion predictions and abundance of individual cell states. This dataset encompasses not only Xenium data containing coordinates and expression profiles of segmented single cells, but also matched scRNAseq (including $5^{\prime}, 3^{\prime}$ and scFFPE-seq) and Visium data, enabling us to generate ground truths for Visium spots regarding cell abundances and cell type proportions (See Methods for details). 3906 Visium spots overlapped with the Xenium data were extracted for comparative analysis (Fig. 3a). Compared with the state-of-the-art algorithms including cell2location, DestVI, CARD, NovoSpaRc, CellTrek, and Tangram, Redeconve demonstrated superior cosine similarities between the predicted cell type proportions and the ground truths for most of the Visium spots (Fig. 3b). Specially, Redeconve exhibited superior performance on more than $60 \%$ and $70 \%$ of spots compared to alternative deconvolution-based or mapping-based methods, respectively (Supplementary Fig. 19). Redeconve, cell2location and Tangram demonstrated comparable performance in estimating the absolute cell abundance within Visium spots, as evidenced by high Pearson's correlation with the ground-truth cell counts indicated by the overlapped cell counts according to the Xenium data, but the performance of Redeconve was more robust to the selection of scRNAseq references (Fig. 3c and Supplementary Fig. 20). Similarly, the performance of Redeconve in reconstructing the expression profiles of different Visium spots was also more robust to the selection of different scRNA-seq references compared with the state-of-the-art algorithms (Fig. 3d and Supplementary Fig. 21).

## Single-cell resolution by Redeconve enables identification of pancreatic cancer-clone-specific T cell infiltration

To demonstrate the power of deconvolution at single-cell resolution on solving practical biological problems, we further investigated the Redeconve results of the human pancreatic ST dataset ${ }^{24}$. The ST is from the original ST platform, and scRNA-seq data from the same individual were obtained through InDrop. Redeconve with single cells as reference outperformed other methods regarding the reconstruction accuracy for almost all the spots (Fig. 4b, c and Supplementary Fig. 7). Using cell types as reference and varying the cell-type resolution from 20 to 318 clusters, Redeconve still resulted in stable superior performance compared with other methods (with the same inputs) (Supplementary Fig. 22), suggesting the advantage of Redeconve by excluding the interference of single-cell reference vs cell-type reference, although Redeconve is the only algorithm designed to take single cells as reference as we demonstrated in the previous sections. Benchmark regarding individual cell types again showed the superiority of Redeconve. We identified marker genes for each cell type (Supplementary Table 2), and calculated the expression consistency
between ST observation and reconstructed profiles by different algorithms across all spots (See Methods for details). Redeconve outperformed other algorithms on most cell types ( $13 / 20$ in top one), especially for cancer, ductal, endocrine cells, and demonstrated comparable performance to the best performers on the remaining of cell types ( $20 / 20$ in top three, Supplementary Fig. 23a, b). In addition, the performance of Redeconve, cell2location, and Tangram was robust to cell type abundance variations in scRNA-seq data, while the performances of DestVI, CARD, and NovoSpaRc were positively correlated with cell type abundances ( $p$-value < 0.05) (Supplementary Fig. 23c).

Histological analysis based on H\&E staining identified four tissue regions: pancreatic, cancer, duct epithelium, and stroma ${ }^{24}$ (Fig. 4a). Redeconve, CARD, and DestVI successfully distinguished the four types of tissue regions, consistent with histological analysis (Supplementary Fig. 24,). Meanwhile, cell2location, NovoSpaRc and Tangram failed in several conditions (Fig. 4d and Supplementary Fig. 24). Further inspection into a specific spot in the upper cancer region (Fig. 4d, the upper zoomed-in piechart) shows that deconvolution-based methods (Redeconve, cell2location, DestVI and CARD) are able to detect fibroblast, which is known to be abundant in pancreatic cancer ${ }^{24,27,28}$, while mapping methods (Tangram and NovoSpaRc) fail in this task.

Then we examined the detailed characteristics of tumorinfiltrating T cells based on these results, which is important to understand the tumor immune microenvironment of pancreatic cancers. The results of cell2location, NovoSpaRc, Tangram and DestVI reported T cells in almost all spots (Fig. 5a), inconsistent with the nature of PDAC as cold tumors; Meanwhile, Redeconve and CARD clearly suggested the sparsity of tumor-infiltrating T cells in pancreatic cancer, consistent with the spatial distribution of T cell-related genes (CD3, IL32 and TMSB4X, Fig. 5a, Supplementary Figs. 25-27). As CARD is limited by the cell-type resolution, it is difficult to provide more detailed insights, but Redeconve analysis enables deeper investigation. We identified three T cells in the reference scRNA-seq data that appeared in multiple ST spots, indexed as "T.cell.8", "T.cell.11" and "T.cell.35" separately (Fig. 5b). By examining their expression profiles in the reference scRNA-seq, we identified T cell 11 as regulatory T cell ( CD4 $^{+}$FOXP3 $^{+}$) and 8 and 35 as CD8 ${ }^{+}$cytotoxic T cells. For fair comparison, we further divided T cells in the scRNA-seq reference data into three groups, i.e., cytotoxic, helper and regulatory T cells and used these three T cell types together with other cell types as reference to re-run other deconvolution algorithms (Supplementary Fig. 28). Consistent with the spatial distribution of CD8 and FOXP3, the result of Redeconve is the most reasonable (Supplementary Figs. 25 and 27). According to the Redeconve deconvolution results, almost all the T cells within cancer region were similar to regulatory T cell 11, and T cell states similar to 8 and 35 only appeared outside or at the edge of the cancer region (Fig. 5b, c), consistent with the immune suppressive status of the cancer region of pancreatic tumors ${ }^{24,29}$.

We further conducted co-localization analysis of these three T cell states with the resting cell states by calculating the Pearson correlation coefficient of abundance across all spots based on the Redeconve results (Fig. 5d). The results suggested that the regulatory T cell state similar to T cell 11 mainly co-localized with macrophages similar to
![](https://cdn.mathpix.com/cropped/2025_11_30_57cdc3a614a421a1b5dfg-05.jpg?height=2111&width=1596&top_left_y=158&top_left_x=246)

Fig. 2 | Performance benchmarking with single-cell inputs and simulated datasets. Redeconve, cell2location and DestVI are currently the only three deconvolution-based tools with the ability to handle thousands of cell states. a cosine similarity between true and reconstructed spatial expression profiles based on Redeconve, cell2location and DestVI with 1000 single cells as input. Each dot represents a spot of the ST data. $\mathbf{b}$ the number of different cell states within each spot estimated by the perplexity of cell state composition per spot for results
of Redeconve, cell2location and DestVI with 1000 single cells as input (See Methods for details). c workflow of generating simulation data. ScRNA-seq data were aggregated to a pseudo-bulk, which was then used for deconvolution analysis and the results were used for downstream analyses in (d). d cosine similarity between true and reconstructed spatial expression profiles vs. number of clusters on simulated pseudo-bulk. PDAC, pancreatic ductal adenocarcinoma. Source data of $2 \mathrm{a}, \mathrm{b}$ and d are provided as a Source Data file.
![](https://cdn.mathpix.com/cropped/2025_11_30_57cdc3a614a421a1b5dfg-06.jpg?height=1935&width=1686&top_left_y=158&top_left_x=199)

Fig. 3 | Benchmarking Redeconve performance on a human breast cancer Xenium dataset. a Left: Overlapped Xenium cells and Visium spots were illustrated on H\&E image. Right: the overlapped region was employed for benchmarking Redeconve performance by introducing different single-cell references to predict expression profiles, cell type proportions, and cell abundances. b line chart of cosine similarities of cell type proportions between ground truths and algorithmbased predictions per spot. $N=3906$ spots for the dataset and spots were sorted by an ascending order of the cosine similarities. c Heatmap illustrating the pairwise Pearson's correlation of cell abundances among the ground truth, Redeconve, cell2location and Tangram based on various single cell references. d violin and box
plot of cosine similarities between observed and reconstructed expression profiles for Redeconve and alternative approaches with different single cell references ( $3^{\prime}, 5^{\prime}$ and scFFPE-seq). The number of independent single cells in the references are 5527, 13,808 and 28,180 respectively. The center line and the bounds of box refer to median, Q1 and Q3 of scores and the whisker equal to $1.5^{*}(\mathrm{Q} 3-\mathrm{Q} 1)$. The minimum and maximum scores refer to Q1-whisker and Q3+whisker. GT, ground truth. scFFPE-seq, single-cell Formalin Fixed Paraffin Embedded sequencing. Source data of $3 b-d$ are provided as a Source Data file. Display items in this figure were manually generated in Inkscape by the authors.
![](https://cdn.mathpix.com/cropped/2025_11_30_57cdc3a614a421a1b5dfg-07.jpg?height=1607&width=1686&top_left_y=158&top_left_x=199)

Fig. 4 | Single-cell deconvolution of a human PDAC (pancreatic ductal adenocarcinoma) ST dataset. a four regions were annotated by histological analysis of the original paper: pancreatic, ductal, cancer and stroma regions ${ }^{24}$. b spatial distribution of the cosine similarity between true and reconstructed expression
profiles per spot by different computational methods. $\mathbf{c}$ pie charts displaying the spatial distribution of the estimated cell type proportion per spot by different computational methods. RBC red blood cell. mDC myeloid dendritic cell. pDC plasmacytoid dendritic cell. Source data are provided as a Source Data file.
macrophages B. 6, 8, and 16 together with duct cells of two different states. Interestingly, T cell 8 and 35 were mainly co-localized with cancer cells, indicating dispersed cancer cells outside the cancer region. Although provided scRNA-seq reference with higher T cell resolution (cytotoxic/helper/regulatory T cells), such co-localization was not observed by other methods (Supplementary Fig. 29).

Furthermore, these two T cell states were separately co-localized with different cancer clones, with T cell state 8 co-localized with cancer clone B and 35 with cancer clone A. Differential gene expression analysis based on the reference scRNA-seq data further indicated the differences between these two pairs of T cells and cancer cells (Fig. 5e, f). It is revealed previously that TM4SF1+ cancer cells denoted latestage while S1004A+ cancer cells (clone B) denoted early-stage ${ }^{30-32}$. Our analysis identified the co-existence of TM4SF1+ cancer cells (clone A) and S1004A+ cancer cells (clone B) with different $C D 8^{+} \mathrm{T}$ cells, which is important to understand the interactions between cancer and

T cells. We found that interferon-induced genes (IFIT1 and IFI44L, for example) and HLA-related genes ( $H L A-A, H L A-B$ and $H L A-C$ ) were all upregulated in cancer clone B (Fig. 5f), and correspondingly T cell state 8, which is colocalized with cancer clone B, had high expression of HMGB2, HLA-B and HLA-C (Fig. 5f), indicating well-stimulated T cell response ${ }^{33,34}$. In contrast, T cell state 35 was HMGB2-negative, HLA-low and TMBS10-positive and co-localized with more A-type macrophages, indicating a less efficacy state ${ }^{33,34}$. Therefore, with accurate deconvolution at the single-cell resolution, Redeconve can reveal detailed cellcell interaction at single-cell level and enables discoveries revealing the underlying mechanisms of tumor immunity.

## Redeconve sheds novel insights into the regulatory mechanisms underlying antibody class switch

Redeconve were further applied to analyze an ST data of human secondary lymphoid organs ${ }^{12}$. We again compared Redeconve with other
![](https://cdn.mathpix.com/cropped/2025_11_30_57cdc3a614a421a1b5dfg-08.jpg?height=1098&width=1688&top_left_y=158&top_left_x=199)

Fig. 5 | Cancer-clone-specific CD8 + T cell infiltration revealed by Redeconve in human pancreatic cancer. a abundance of T cells per spot estimated by different methods. $\mathbf{b}$ single-cell identity of infiltrated T cells revealed by Redeconve. The three T cells are indexed as "T.cells.8", "T.cells.11", "T.cells.35" separately. c singlecell identity of different cancer clone cells revealed by Redeconve, together with their abundance difference. d co-localization of the three T cell states with other cellular states. Nodes represent single cells and edges represent co-localization (Pearson correlation of cell abundance $>0.4$ ). Cancer clone-specific CD8 + T cell infiltration was revealed. e dot plot displaying characteristics genes among the three T cell states with different spatial preference with cancer clones A and B .
f volcano plot displaying differentially expressed genes between the two cancer clones. The blue and red points refer to up-regulated genes in clones A and B-enriched spots, respectively. Vertical dashed line shows the cutoff of log fold change $( \pm 0.3)$. Horizontal dashed line shows the threshold of $-\lg p(1.301)$. T cell response-related genes including interferon-stimulating genes and human leukocyte antigens were up-regulated in clone B-enriched cells. The two-side exact test was applied in edgeR for the statistical test and the $p$-values were calculated without adjustments. Treg, T regulatory. Source data are provided as a Source Data file.
methods on this dataset. In terms of cosine similarity-based reconstruction accuracy, Redeconve achieved mean similarities of 0.868 and significantly outperformed other methods (Fig. 1d). Redeconve achieved high reconstruction accuracy for almost all spots, while, as for other methods, low similarities regions were obvious (Supplementary Fig. 30). We further checked the sparsity of the results by calculating LO-norm. LO-norm of Redeconve has a reasonable distribution between 4 and 32, indicating that only dozens of cell states appear in one spot. In contrast, other methods except CellTrek demonstrated results that almost all cell types appeared in every spot. CellTrek, a mapping-based algorithm, reached low level of $L O$-norm by generating many "zero-cell" spots, of which Redeconve successfully reconstructed the cellular composition (Supplementary Fig. 31).

We further characterized the spatial heterogeneity at single cell resolution to explore the potential regulators of antibody class switch based on this human lymph node data. During the antibody maturation, an activated B cell can change its antibody production from IgM to either $\operatorname{IgA}, \operatorname{IgG}$, or $\operatorname{IgE}$ depending on the functional requirements, which is termed as class switching ${ }^{35}$. However, the detailed regulators underlying antibody class switching is unclear. Consistent with previous examples, Redeconve outperformed other methods in reconstructing the ST gene expression profiles for almost all spots (Fig. 1d). Spatial pie chart showed that Redeconve produced obvious regional division, while other methods showed blurred or even no boundaries
(Fig. 6a). CellTrek failed to analyze some of the spots. Furthermore, compared with cell-type deconvolution, Redeconve identified 159 different cell states from 17 cell types (Supplementary Fig. 2). 12 different B plasma cell states were identified in the ST data, which can be further divided into 3 groups ( $\operatorname{Ig} \mathrm{A}+$, $\operatorname{IgG}+$ and negative) based on the expression of IGHA and IGHG genes. Interestingly, we found that IgA+ and IgG+ B plasma cells are spatially mapped to spots in different regions with little overlap, which means that we could define $\operatorname{IgA}+$ and IgG+ spots based on the abundance of those B plasma cells (Fig. 6b). Next, we took one spot in each of the two regions for detailed inspection at the single-cell resolution. The cell proportion of the two spots shows that $C D 8^{+}$T cells account for a large proportion in the IgA+ spot, suggesting latent interactions between $C D 8^{+} \mathrm{T}$ cells and $\operatorname{IgA}+\mathrm{B}$ plasma cells (Fig. 6c). To confirm the universality of such phenomenon, we conducted differential gene expression analysis between $\operatorname{Ig} \mathrm{A}$ + and IgG+ spots to identify up-regulated and down-regulated genes (Fig. 6d). As we expected, IGHA and IGHG were the most differentiallyexpressed genes; Genes associated with T cells (TRAC, TRBC2, CD3D, CD8A for example) were more up-regulated in $\operatorname{IgA}+$ spots, confirming the existence of such interaction. Since lymph node is one of the organs that generate $\operatorname{Ig} \mathrm{A}+$ plasma cells, the $\operatorname{Ig} \mathrm{A}+$ spots might be the potential induction sits for $\operatorname{Ig} \mathrm{A}+$ plasma cells, and $C D 8^{+} \mathrm{T}$ cells may play an important role in such process (Fig. 6d). Further co-localization analysis provides more insights (Fig. 6e). We found co-localization of
![](https://cdn.mathpix.com/cropped/2025_11_30_57cdc3a614a421a1b5dfg-09.jpg?height=2053&width=1694&top_left_y=158&top_left_x=197)

Fig. 6 | Single-cell deconvolution of a human secondary lymphoid organ ST dataset by Redeconve revealed differences between IgA+ and IgG+ spots regarding cellular composition. a pie chart displaying the spatial distribution of the estimated cell type proportion by different methods. $\mathbf{b}$ spatial distribution of $\operatorname{IgA}+$ and IgG+ B plasma cells revealed by Redeconve. c comparison of the cell proportion of two selected spots (the $\lg \mathrm{A}+$ and $\operatorname{lgG}+$ spots in Fig. 6a with green squares). d volcano plots showing the differential gene expression between $\operatorname{Ig} \mathrm{A}+$ and $\operatorname{IgG}+$ spots. The red and blue point refer to up-regulated genes in $\operatorname{IgG}+$ and $\operatorname{IgA}+$ spots respectively. Vertical dashed line shows the cutoff of log fold change ( ± 0.3). Horizontal dashed line shows the threshold of $-\lg (p)$, namely 1.301. The two-side
exact test was applied in edgeR for the statistical test and the $p$-values were calculated without adjustments. e co-localization network of $\operatorname{IgA}+$ and $\operatorname{IgG}+$ B plasma cells within the ST data. Nodes represent single cells and edges represent colocated single cells (Pearson correlation of cell abundance $>0.2$ ). Abbreviations: GC germinal center, DZ dark zone, LZ light zone, prePB preplasmablast, mem memory, cDC classical dendritic cell, Endo endothelial, FDC follicular dendritic cell, ILC innate lymphoid cell, NK natural killer, NKT natural killer T, TfH T follicular helper, Treg T regulatory, VSMC vascular smooth muscle cell. Source data of 6a-d are provided as a Source Data file.

IgA+ plasma cells with $C D 8^{+}$cytotoxic T cells, consistent with previous observation that $C D 8^{+}$cytotoxic T cells can help the formation of $\operatorname{Ig} \mathrm{A}+$ plasma cells ${ }^{36,37}$. Furthermore, co-location of IgG+ plasma cells and macrophages was identified (Fig. 6e), indicating the roles of macrophages during the genesis of IgG+ plasma cells ${ }^{38,39}$. Hence, deconvolution at single cell resolution by Redeconve gains additional insights that may be helpful for uncovering previously opaque biological question.

## Discussion

Integrative analysis of disassociated single-cell and in situ ST data is pivotal to construct a comprehensive map of the cellular composition and interactomes of tissues. However, because of technological limitations, current computational methods for integrative analysis of single-cell and ST data are limited to the cell type resolution. To deep mine the biomedical information hidden in the single-cell and ST data, here we present Redeconve, a single-cell resolution deconvolution algorithm for integrative analysis of ST data with sc/snRNA-seq data as reference based on a quadratic programming model with regularization of cell-cell similarity, which enables building of comprehensive spatial maps at single-cell resolution for diverse tissues.

We performed stringent evaluation on multiple datasets from a diverse set of ST platforms. The results suggested superiority of Redeconve compared with the state-of-the-art deconvolution-based and mapping-based algorithms in terms of resolution, accuracy, sparsity, robustness, and computational speed. Such improvement from cell-type to single-cell resolution unlocks novel biological discoveries as exemplified by applications in human pancreatic cancer and lymph node samples.

While Redeconve enables deconvolution at single-cell resolution and thus will be a powerful tool for biomedical discoveries, matching between scRNA-seq and ST data appears to be an important factor determining the quality of deconvolution analysis as shown by our evaluation on different tissues (Fig. 1d). Therefore, construction and selection of reference scRNA-seq data according to the specific ST data configuration will be critical in future applications.

Although Redeconve demonstrates superior computational efficacy compared with the state-of-the-art deconvolution algorithms, the single-cell resolution may require extensive computational cost for resolving thousands of cellular states, especially when the cellular throughput of scRNA-seq technologies increases exponentially. Because of the computational complexity of quadratic programming, Redeconve can currently resolve thousands of cellular states based on a standard machine. An enhanced version based on algorithmic innovation or hardware acceleration is needed to handle scRNA-seq datasets of tens of thousands of cellular states.

Deconvolution at single-cell resolution unlocked by Redeconve may also benefit the imputation of ST data with the aid of the rich information in scRNA-seq data. Redeconve has implemented a function to reconstruct the gene expression profiles of individual spots based on the single-cell deconvolution results based on a parsimony principle. The imputed ST data may be more informative to dissect the cellular states of specific tissues.

In summary, we present an algorithm named as Redeconve for conducting deconvolution-based analysis of scRNA-seq and ST data at single-cell resolution. The usage of Redeconve is expected to help mapping the cellular architecture at fine granularity across diverse biomedical situations including tumor, immune, development, neurology, and other health and disease conditions. Applications to human pancreatic cancer and lymph nodes showed the potential of Redeconve to bring completely novel insights due to the single-cell resolution unlocked and the superior technical metrics of Redeconve compared to the current state-of-the-art algorithms. We expect Redeconve will be a useful tool to advance the application of scRNAseq and ST technologies in diverse research disciplines.

## Methods

## Algorithm

Model overview. In general, we apply an improved linear regression model to deconvolute ST data at single-cell resolution. Given a singlecell (or single-nucleus) expression matrix $X$ with dimensions $n_{\text {genes }} \times n_{\text {cells }}$ and a ST expression matrix $Y$ with dimensions $n_{\text {genes }} \times n_{\text {spots }}$ as input, Redeconve returns a matrix $\beta$ with dimensions $n_{\text {cells }} \times n_{\text {spots }}$ indicating the estimated number of each cell in each spot. The goal of our model is to optimize the following loss function for each spot separately:

$$
\begin{gathered}
L(\beta):=\sum_{j=1}^{J}\left(y_{j}-\sum_{i=1}^{I} x_{i j} \beta_{i}\right)^{2}+c \cdot \sum_{i_{1} \neq i_{2}} R_{i_{1}, i_{2}}\left(\beta_{i_{1}}-\beta_{i_{2}}\right)^{2} \\
\text { s.t. } \beta_{i} \geq 0 \text { for } i=1,2, \ldots, I
\end{gathered}
$$

Here $i=1,2, \ldots, I$ denotes cells and $j=1,2, \ldots J$ denotes genes. The first term is the traditional Least Square (LS) term and the second term is a regularization term, $c$ is a hyperparameter tuning the weight between the two terms. We will later explain the regularization term in details.

Note that this is a typical quadratic programming problem, so we can rewrite our goal as:

$$
\begin{gathered}
\min _{\beta} \frac{1}{2} \beta^{T} G \beta-d^{T} \beta \\
\text { s.t. } a^{T} \beta \geq b
\end{gathered}
$$

Where $G$ is the Hessian matrix, $d^{T}=\left(2 \sum_{j} y_{j} x_{1 j}, \ldots, 2 \sum_{j} y_{j} x_{l j}\right)$, and $a^{T}, b$ are separately

$$
a^{T}=\left(\begin{array}{cccc}
1 & 0 & \cdots & 0 \\
0 & 1 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 1
\end{array}\right), b=\left(\begin{array}{c}
0 \\
0 \\
\vdots \\
0
\end{array}\right)
$$

So we can efficiently solve this problem with the solve. $Q P$ function in R package "quadprog".

The regularization term. In sc/snRNA-seq data, the collinearity among cells is serious: cells of the same cell type have very similar expression profiles. This problem would lead to instability of coefficients and reduction of efficiency when directly doing linear regression. To solve this collinearity problem, we further include a regularization term into the loss function. By add this term, we aim at stabilizing the coefficients while having minor effect on the residuals.

In the regularization term $c \sum_{i_{1} \neq i_{2}} R_{i_{1}, i_{2}}\left(\beta_{i_{1}}-\beta_{i_{2}}\right)^{2}, R_{i_{1}, i_{2}}$ is a measure of similarity between cell $\dot{i}_{1}$ and $\dot{i}_{2}$, which is

$$
R_{i_{1}, i_{2}}=\left\{\begin{array}{cc}
r_{i_{1}, i_{2}}, & r_{i_{1}, i_{2}}>0 \\
0, & r_{i_{1}, i_{2}} \leq 0
\end{array}\right.
$$

Where $r_{i_{1}, i_{2}}$ is the Pearson correlation coefficient between cell $i_{1}$ and $i_{2}$. Namely, when the Pearson correlation coefficient is greater than zero, $R_{i_{1}, i_{2}}$ is equal to the Pearson correlation coefficient; otherwise $R_{i_{1}, i_{2}}$ is zero. So, we manually bring the coefficients of cells whose expression profile is similar closer. By doing this, we can guarantee the robustness and precision of our result.

Determination of the hyperparameter. A key point of this model is how to select the hyperparameter: an extremely small hyperparameter will make the regularization term ineffective, while an extremely large one will greatly affect the fitting residuals. An ideal hyperparameter should be as large as possible while affecting the fitting residual as little as possible. Here we offer 3 ways to set the hyperparameter:

1. "default": use the default hyperparameter we set according to the number of cells and genes;
2. "customized": set the hyperparameter arbitrarily by the user;
3. "autoselection": automatically calculate and select the optimal hyperparameter.

In mode "default", we use the following formula to set the hyperparameter:

$$
c=c_{0} \cdot n_{\text {genes }} / n_{\text {cells }}^{2}
$$

Where $c_{0}$ is a predetermined constant and is set to $10^{5}$. The idea of this formula is: (1) the LS term is approximately proportional to $n_{\text {genes }}$, so as $n_{\text {genes }}$ increases $c$ should synchronously increase; (2) the regularization term is approximately proportional to the square of $n_{\text {cells }}$, so as $n_{\text {cells }}$ increases $c$ should decrease by $n_{\text {cells }}^{2}$.

In mode "autoselection", we apply the following method to determine the optimal hyperparameter:

1. We first calculate a hyperparameter $c_{d}$ according to the formula in mode "default", and set up a series of hyperparameter $c_{1}, c_{2}, c_{3}, c_{4}, c_{5}$ as $0.01 c_{d}, 0.1 c_{d}, c_{d}, 10 c_{d}, 100 c_{d}$;
2. Then we run deconvolution with these hyperparameters separately, and calculate the residual $\varepsilon_{i}$ for each $c_{i}$;
3. We further calculate:

$$
d_{i}=\frac{\Delta \varepsilon}{\Delta c}=\frac{\varepsilon_{i+1}-\varepsilon_{i}}{c_{i+1}-c_{i}}
$$

4. We check these $d_{i}$, then choose $c_{i}$ that maximizes $d_{i}$ as the optimal hyperparameter (This indicates: if the parameter continues to increase, the residual will increase significantly). Namely, we choose $c_{i}$ that satisfies:

$$
\max _{i \in 1,2, \cdots, I} d_{i}=\frac{\varepsilon_{i+1}-\varepsilon_{i}}{c_{i+1}-c_{i}}
$$

By this procedure, we can get the hyperparameter that maximizes the power of regularization term while having minor effect on the LS term.

We use examples to illustrate the effect of hyperparameters on the results. We applied Redeconve to the human lymph node dataset with a series of different hyperparameters from 0 to 1 e 08 , then calculated the deconvolution residuals (RMSE_normal) to evaluate the effect of hyperparameter (Supplementary Fig. 32). The results showed that an optimal hyperparameter can enhance the deconvolution precision in addition to avoiding co-linearity caused by closely similar cell states. Also, the hyperparameter would also affect the number of cell states selected in the result. A bigger hyperparameter would lead to more cell states selected (Supplementary Fig. 33). We set the hyperparameter as 0 and 1 e 04 separately on the PDAC dataset. With a hyperparameter of 1 e 04 , more T cells were detected than a hyperparameter of zero in the PDAC dataset (Supplementary Fig. 34). Considering the distribution of CD3+ cells (Shown in Supplementary Figs. 25-27), this example clearly illustrates how the hyperparameter enables biological discovery.

## Data preprocessing

To run the deconvolution, the following data preprocessing steps are necessary. Note that some steps are alternative according to users' needs.

1. Get the expression profiles of cell type/Sampling of single cells. If a cell-type deconvolution is to be run, we will estimate the expression profile $\bar{x}_{i j}$ of cell type $i$ and gene $j$ as the average expression of gene $j$ across all cells within cell type $i$. If a single-cell deconvolution is to be run and the number of single cells is overwhelming, we will take stratified samples of cells by cell type to get a rational number of cells.
2. Gene filtering. Deconvoluting with tens of thousands of genes is time-consuming or even misleading, so we select highly variable genes before deconvolution for computational efficacy. Filtering criteria include the following three standards: (1) These genes appear in both sc/snRNA-seq data and ST; (2) The variance of these genes in sc/snRNA-seq data must be larger than a threshold (default is 0.025 ); (3) The average counts per spot must be bigger than a threshold (default is 0.003 ). This finally results in $\sim 8000$ genes for deconvolution. Redeconve allows deconvolution without gene filtering with higher computational cost.
3. Normalization of reference. We add a pseudo-count of 0.5 to the "zeros" in sc/snRNA-seq data, and normalize sc/snRNA-seq data to TPM (transcripts per million). Preprocessing operations are not needed for ST data.

## Real datasets for benchmarking

PDAC. ST data of a human pancreatic ductal adenocarcinomas (PDAC-A) with 438 spots and sample-matched scRNA-seq data (InDrop) with 1926 single cells across 20 cell types were integrated by Moncada et al., and an intersection of 19,736 genes was used in our study. The annotation of four main structural regions based on histological analysis by Moncada et al. was used during our analysis to depict the spatial characteristics of the ST data.

Human lymph node. Human lymph node Visium data were downloaded from the 10x Genomics website (https://www.10xgenomics. com/resources/datasets/human-lymph-node-1-standard-1-1-0), which includes a total number of 4035 spots. ScRNA-seq data were collected from Kleshchevnikov et al, of which 73,260 cells across 34 cell types were collected. Since this scRNA-seq dataset captured a wide spectrum of immune cell states spanning lymph nodes, tonsils and spleen, we used it as reference to reveal the phenotypic diversity of immune cells when deconvoluting at single cell resolution.

Mouse cerebellum. The DropViz scRNA-seq dataset were generated by Saunders A. et al. and were collected by Cable D. M. et al. along with the annotations of the cells. The Slide-seq mouse cerebellum data were collected by Cable D. M. et al. using the Slide-seq v2 protocol ${ }^{11}$. Both of these datasets were downloaded from https://singlecell. broadinstitute.org/single_cell/study/SCP948/robust-decomposition-of-cell-type-mixtures-in-spatial-transcriptomics\#study-download.

Human breast cancer. Human Breast Cancer Visium data related to the Wu et al. study ${ }^{40}$ was available at https://zenodo.org/record/ 4739739\#.Ys0v6jdBy3D. Sample 'CID4290' that includes 2426 in tissue spots was used for deconvolution. ScRNA-seq data that includes 100,064 single cells with annotations (Access number: GSE176078, the NCBI GEO database) served as reference to do deconvolution analysis.

Human testis. The processed Human Testis Slide-seq dataset was download from https://www.dropbox.com/s/q5djhy006dq1yhw/ Human.7z?dl=0 and sample 'Puck5' with 36,591 spots was used for evaluation in this study ${ }^{41}$. The reference scRNA-seq data that includes 6490 single cells was obtained from the NCBI GEO database with access number GSE112013, and the corresponding annotations were available in the supplementary information Table S1 by Guo et al.

Mouse brain. 10x Visium and snRNA-seq data (includes annotation) were available in the ArrayExpress database with accession numbers E-MTAB-11114 and E-MTAB-11115, respectively ${ }^{12}$. Sample 'ST8059048' containing 2987 spots was used for evaluation in this study, and all 40,532 single cells across 59 cell types served as reference. In addition, the corresponding data of nuclei counts estimated by histological image segmentation based on deep learning s was downloaded from
https://github.com/vitkl/cell2location_paper/blob/master/notebooks/ selected_results/mouse_visium_snrna/segmentation/144600.csv.

Human breast cancer xenium. The Human Breast Cancer Xenium dataset is available at https://www.10xgenomics.com/products/ xenium-in-situ/preview-dataset-human-breast. A single FFPE tissue block was analyzed by scFFPE-seq, Visium and Xenium. In addition, $3^{\prime}$ and 5' gene expression data from dissociated tumor cells is also available ${ }^{4}$.

## Comparing Redeconve with alternative methods

We compared Redeconve with recently developed deconvolutionbased methods (cell2location, DestVl ${ }^{13}$ and CARD ${ }^{10}$ ) as well as mapping-based methods (NovoSpaRc, CellTrek ${ }^{8}$ and Tangram ${ }^{7}$ ).

Criteria of selecting alternative methods. In considering which methods to include for the comparison, we required methods that (1) are specifically designed for end-to-end estimating the abundance/ proportion of cells or cell types using scRNA-seq and ST data as input; (2) demonstrate superior performance in the corresponding publications and third-party evaluation papers; and (3) are peer reviewed with a publicly available software implementation before Dec 2022.

Parameter setting. Prediction results for the 6 datasets were obtained by running the corresponding programs of the algorithms aforementioned based on the default settings except some special considerations: (1) 1000 cells were randomly selected in NovoSpacRc to avoid large number of total cells; (2) 1000 stratified samples of cells were used for Redeconve in almost all the datasets except PDAC where we used total 1926 cells; ( 3 ) minCountGene and minCountSpot of the createCARDObject function were set to 0 to prevent unexpected gene or spot filtering in CARD. The output of each method was either a cell-by-spot matrix represented absolute abundance (Redeconve, Tangram) or proportion (NovoSpaRc) of single cells existing at each spot or estimated cell-type abundance (cell2location) or proportion (DestVI, CARD) matrix except CellTrek, of which the outcome was predicted spatial coordinates for individual cells. Hence, for CellTrek, we obtained cell-by-spot abundance matrix by assigning single cells to specific spots according to whether the spot area designed by ST platforms covered the predicted coordinates. We only evaluated CellTrek on the two 10x Genomics Visium-based datasets (human lymph node and human breast cancer) because of running errors on other ST datasets in our computational environment.

Calculating performance metrics. To demonstrate superior performance of Redeconve, we firstly estimated predicted expression profiles for spatial spots. For all datasets, spot-wise cosine similarities, Pearson's correlations and RMSEs between observed and predicted spot-by-gene expression matrix were calculated. In order to compute these metrics based on the output of each algorithm, we calculated the predicted expression matrix through two ways: (1) for Redeconve, NovoSpaRc, CellTrek and Tangram, we multiplied spot-by-cell abundance or proportion matrix by the cell-by-gene sc/snRNA expression matrix; (2) for cell2location, DestVI and CARD, we multiplied the celltype abundance or proportion matrix by the reference cell-type expression matrix, where the reference was generated through averaging sc/snRNA expression data according to cell types. When calculating RMSEs, the total number of UMIs for each spot in both observed and predicted expression profile was normalized to $n_{\text {genes }}$. We then estimated sparsity of the results through calculating cell-type proportion matrices of all programs and comparing the results according to cell-type information entropy and $\mathrm{L}_{\mathrm{O}}$ norm. The $\mathrm{L}_{\mathrm{O}}$-norm represents number of cell types present at each spot (nonzero values). We also evaluated the performance of cell abundance estimation by Pearson's correlation between results of individual methods (Redeconve,
cell2location, CellTrek and Tangram) and the cell numbers estimated by histological image segmentation based on deep learning for the mouse brain dataset. Finally, computational efficiencies were estimated through comparing total time spent by each algorithm on a computer with Intel(R) Xeon(R) Platinum 8253 CPU, where we set the maximum number of cores to 96 . In addition, we tested the run time of these programs on a single NVIDIA A40 card if GPU acceleration supported (cell2location, DestVI, NovoSpaRc, and Tangram).

Assessment at single cell resolution. Cell-by-spot abundance matrix is required for comparison among deconvolution-based methods at single-cell resolution. We, therefore, applied Redeconve with 1000 single cells sampled from the reference scRNA-seq data for the two ST datasets (PDAC and human lymph node) and assigned every single cell a unique cell type since cell2location, DestVI and CARD only support cell-type deconvolution. The result matrices of Redeconve, cell2location and DestVI (no result was available for CARD because of running errors) was obtained according to the corresponding default settings. Cosine similarity, information entropy, perplexity and runtime efficiencies were evaluated as mentioned above.

Information entropy and perplexity. We calculate Information entropy $H$ and perplexity $P$ for each spot separately by the following formula:

$$
H=-\sum_{i} \beta_{i} \log _{2}\left(\beta_{i}\right)
$$

$$
P=2^{H}
$$

where $i=1,2, \ldots, I$ denotes different cell states. $\beta_{i}$ were normalized in advance so that their sum equaled to 1 (i.e., they denote proportion rather than absolute abundance). When $\beta$ is uniformly distributed (namely $\beta_{i}$ is a constant, $\frac{1}{7}$, for all $i$ ), we can know by simple calculation that the perplexity equals to the number of states $I$. This means that perplexity can reveal the number of states when the distribution is uniform. For other distributions, perplexity can also approximately represent the number of states. "Number of states" in the setting of single-cell deconvolution refers to "number of cell states (or types)". Namely, the perplexity of each spot can approximately represent the number of cell states/types occurred in this spot. By calculating perplexity on simulated and real datasets, we have verified that perplexity showed good consistency with number of non-zero cell types/states in the result, but poor consistency with absolute cell abundance (Supplementary Fig. 35 and Supplementary Table 1).

Cell-type level benchmark based on the PDAC dataset. Marker genes were first identified for each cell types (Supplementary Table 2). Then, for each cell type, similarities of marker genes expression between ST observation and reconstructed profiles by different algorithms across all spots were calculated. Ranks of cosine similarities of individual cell types were used as metrics to summarize the overall performance. In addition, linear regression and statistical test were used to show relationship between cell type abundances and performance metrics.

Cell abundance of ST spots on PDAC dataset. To generate ground truth of cell abundance for each ST spot, we first registered H\&E and fluorescent images using Adobe Photoshop CC. Such registration enabled the determination of spatial coordinates for ST spots. After that, Cellpose ${ }^{13}$ was applied through squidpy ${ }^{14}$ to detect cell nuclei from the H\&E image. Finally, we counted the absolute number of nuclei within each spot and referred to these values as cell abundance.

## Generating and analyzing simulation datasets

We used 3 scRNA-seq data to generate simulation data separately: PDAC, human lymph node and human testis. Prior to analysis, all scRNA-seq data were down-sampled to around 1000 cells, with the exception of PDAC which contained a total of 1926 cells. To generate a pseudo-bulk for subsequent deconvolution, all single-cells were aggregated together and assigned an abundance value of 1 . To perform deconvolution, we clustered the scRNA-seq reference with 5 different resolutions using FindCluster() function in Seurat package. Together with directly using all single-cells as input, this results in 6 groups of references. Then the differently annotated references were used for deconvolution by Redeconve and cell2location and the results were used to compare with ground-truth, calculate cosine similarity and perplexity (Fig. 2c, d and Supplementary Figs. 16-18).

## Benchmarking on human breast cancer Xenium dataset

The Human Breast Cancer Xenium dataset contains scRNA-seq, Visium and Xenium data for a single FFPE tissue block. By mapping Xenium cells to Visium spots, it becomes possible to generate ground truth data regarding cell abundances and cell type proportions. To achieve this, we chose Replicate 1 of Xenium data to align spatial locations of Xenium cell centers to corresponding H\&E images through translation and rotation. After that, a key-point registration approach was employed to align H\&E images in Xenium and Visium data based on 155 manually identified landmark features on commonly shared microstructures. Then, FindHomography() function in cv2 package with RANSAC method was applied to transform Xenium to Visium coordinates. Hence, the ground truths of cell abundance were generated through counting the transformed cell centers located within each Visium spot. To further generate ground truths of cell type proportion for Visium spots, we labeled each cluster in scFFPE-seq and Xenium data with a corresponding cell type designation (Supplementary Table. 3-4). The proportions of various types of Xenium cells in Visium spots were considered as ground truth cell type proportions.

Based on the generated ground truths, we computed spot-wise cosine similarities between predicted and ground truth cell type proportions for Redeconve and alternative methods. In this approach, we chose scFFPE-seq data as reference for the deconvolution. In addition, Pearson's correlation was applied to measure the performance of cell abundance estimation for Redeconve, cell2location and Tangram. Finally, a selection of distinct single-cell references (including $5^{\prime}, 3^{\prime}$, and scFFPE-seq) were applied for the purpose of assessing robustness of the computational algorithms.

## Downstream analyses after Redeconve deconvolution

Human lymph node. We firstly ran Redeconve on default setting to obtain deconvolution result at single-cell resolution. Then, we investigated the spatial distribution of plasma cells after grouping these plasma cells into $\operatorname{Ig} \mathrm{A}+$, $\operatorname{Ig} \mathrm{G}+$ and others based on the expression of IGHA1, IGHG1, IGHG3 and IGHG4. Ig $\mathrm{A}+$ and IgG + spots were determined by the following three steps: (1) identifying the top $50 \%$ spots with the highest abundance of $\operatorname{IgA}+$ and $\operatorname{IgG}+$ plasma cell enriched, which were named as spot sets A and G ; (2) identifying the difference sets between A and G, and naming as AD and GD; (3) selecting spots from AD and GD with the top $1 \% \operatorname{IgA}+$ and IgG+ plasma abundance, which were assumed to be $\operatorname{IgA}+$ and IgG+ spots respectively. EdgeR ${ }^{42}$ was applied to perform differential gene expression analysis and identified significantly differential genes between $\operatorname{IgA}+$ and IgG+ spots. Then, we calculated Pearson's correlation coefficient among single cell states in the reference across $\operatorname{Ig} \mathrm{A}+$ and $\operatorname{IgG}+$ spots and took single cells as nodes and correlated cells (Pearson > 0.2) as edges to generate the cell-cell colocation network.

PDAC. We ran Redeconve with all the 1926 single cells as reference, and all the parameters were kept default. For downstream analyses, we
first compared Redeconve with existing tools as described in the aforementioned sections. Then, to study the distribution of T cells, we distinguished from NK cells T cells by the expression of CD3D, CD3E or CD3G in the scRNA-seq data. We further picked out those T cells that frequently appeared in the ST spots (T cells 8, 11, and 35). To study the spatial colocalization of these T cells with other cells, we calculated the Pearson's correlation of cell abundance across spatial spots, and generated a colocalization network of single cell resolution using those cell pairs whose Pearson correlation were greater than 0.4 with the R package igraph ${ }^{43}$.

## Statistics and reproducibility

For all datasets except for PDAC, we down sampled the sc/snRNA-seq reference to around 1000 cells. Stratified sampling was performed when cell types are available, otherwise simple random sampling was performed. The exact number of chosen cells for each dataset are as follows: human breast cancer: 1001, human lymph nodes: 1000, human testis: 999, Mouse Brain: 1003, Mouse cerebellum: 1003, human breast cancer Xenium (scFFPE): 1001, human breast cancer Xenium (3'): 998, human breast cancer Xenium (5'): 1002. The seed was set to 2233. All other parts of this study do not involve randomization. The Investigators were not blinded to allocation during experiments and outcome assessment.

## Reporting summary

Further information on research design is available in the Nature Portfolio Reporting Summary linked to this article.

## Data availability

All relevant data supporting the key findings of this study are available within the article and its Supplementary Information files. The PDAC data used in this study are available in the Gene Expression Omnibus database under accession code GSE111672. The processed human lymph nodes Visium data are available at 10x Genomics website [https://www.10xgenomics.com/resources/datasets/human-lymph-node-1-standard-1-1-0]. The processed human lymph nodes scRNA-seq data are available from Kleshchevnikov et al. [https://cell2location.cog. sanger.ac.uk/browser.html]. The mouse cerebellum data used in this study are available in the Single Cell Portal database under accession code SCP948 [https://singlecell.broadinstitute.org/single_cell/study/ SCP948/robust-decomposition-of-cell-type-mixtures-in-spatial-transcriptomics\#study-download]. The processed human breast cancer Visium data are available at zenodo [https://zenodo.org/record/ 4739739\#.Ys0v6jdBy3D]. The processed human breast cancer scRNAseq data used in this study are available in the Gene Expression Omnibus database under accession code GSE176078. The processed human testis Slide-seq data are available at dropbox [https://www. dropbox.com/s/q5djhy006dq1yhw/Human.7z?dl=0]. The processed human testis scRNA-seq data used in this study are available in the Gene Expression Omnibus database under accession code GSE112013. The processed mouse brain Visium data used in this study are available in the ArrayExpress database under accession code E-MTAB-11114. The processed mouse brain snRNA-seq data used in this study are available in the ArrayExpress database under accession code E-MTAB-11115. The processed Visium, 3’ scRNA-seq, 5’ scRNA-seq and scFFPE-seq for human breast cancer Xenium dataset are available at 10x Genomics website [https://www.10xgenomics.com/products/xenium-in-situ/ preview-dataset-human-breast]. Source data are provided with this paper.

## Code availability

The codes used to generate the figures in this paper is available at https://codeocean.com/capsule/1351962/tree/v1. The package is available on GitHub with detailed documentation at https://github.com/ ZxZhou4150/Redeconve, https://doi.org/10.5281/zenodo.8384152 ${ }^{21}$.

## References

1. Stickels, R. R. et al. Highly sensitive spatial transcriptomics at nearcellular resolution with Slide-seqV2. Nat. Biotechnol. 39, 313-319 (2021).
2. Eng, C. L. et al. Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH. Nature 568, 235-239 (2019).
3. Chen, K. H., Boettiger, A. N., Moffitt, J. R., Wang, S. \& Zhuang, X. RNA imaging. Spatially resolved, highly multiplexed RNA profiling in single cells. Science 348, aaa6090 (2015).
4. Janesick, A. et al. High resolution mapping of the breast cancer tumor microenvironment using integrated single cell, spatial and in situ analysis of FFPE tissue. bioRxiv, 2022.2010.2006.510405, https://doi.org/10.1101/2022.10.06.510405 (2022).
5. He, S. et al. High-plex imaging of RNA and proteins at subcellular resolution in fixed tissue by spatial molecular imaging. Nat Biotechnol 40, 1794-1806 (2022).
6. Moriel, N. et al. NovoSpaRc: flexible spatial reconstruction of singlecell gene expression with optimal transport. Nat. Protoc. 16, 4177-4200 (2021).
7. Biancalani, T. et al. Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nat. Methods 18, 1352-1362 (2021).
8. Wei, R. et al. Spatial charting of single-cell transcriptomes in tissues. Nat. Biotechnol. 40, 1190-1199 (2022).
9. Vahid, M. R. et al. High-resolution alignment of single-cell and spatial transcriptomes with CytoSPACE. Nat Biotechnol 41, 1543-1548 (2023).
10. Ma, Y. \& Zhou, X. Spatially informed cell-type deconvolution for spatial transcriptomics. Nat. Biotechnol. 40, 1349-1359 (2022).
11. Cable, D. M. et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nat. Biotechnol. 40, 517-526 (2022).
12. Kleshchevnikov, V. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat. Biotechnol. 40, 661-671 (2022).
13. Lopez, R. et al. DestVI identifies continuums of cell types in spatial transcriptomics data. Nat. Biotechnol. 40, 1360-1369 (2022).
14. Dong, R. \& Yuan, G. C. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biol. 22, 145 (2021).
15. Elosua-Bayes, M., Nieto, P., Mereu, E., Gut, I. \& Heyn, H. SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Res. 49, e50 (2021).
16. Sun, D., Liu, Z., Li, T., Wu, Q. \& Wang, C. STRIDE: accurately decomposing and integrating spatial transcriptomics using singlecell RNA sequencing. Nucleic Acids Res. 50, e42 (2022).
17. Bae, S. et al. CellDART: cell type inference by domain adaptation of single-cell and spatial transcriptomic data. Nucleic Acids Res. 50, e57 (2022).
18. Geras, A. et al. Celloscope: a probabilistic model for marker-genedriven cell type deconvolution in spatial transcriptomics data. Genome Biol 24, 120 (2023).
19. Song, Q. \& Su, J. DSTG: deconvoluting spatial transcriptomics data through graph-based artificial intelligence. Brief Bioinform. 22, https://doi.org/10.1093/bib/bbaa414 (2021).
20. Andersson, A. et al. Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography. Commun. Biol. 3, 565 (2020).
21. Zhou, Z., Zhong, Y., Zhang, Z. \& Ren, X. Spatial transcriptomics deconvolution at single-cell resolution using Redeconve. Zenodo, https://doi.org/10.5281/zenodo. 8384152 (2023).
22. Stringer, C., Wang, T., Michaelos, M. \& Pachitariu, M. Cellpose: a generalist algorithm for cellular segmentation. Nat. Methods 18, 100-106 (2021).
23. Palla, G. et al. Squidpy: a scalable framework for spatial omics analysis. Nat. Methods 19, 171-178 (2022).
24. Moncada, R. et al. Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas. Nat. Biotechnol. 38, 333-342 (2020).
25. 10x Genomics Support. V1_Human_Lymph_Node - Datasets - Spatial Gene Expression, [https://support.10xgenomics.com/spatial-geneexpression/datasets/1.1.0/V1_Human_Lymph_Node](https://support.10xgenomics.com/spatial-geneexpression/datasets/1.1.0/V1_Human_Lymph_Node) (2020).
26. Guo, J. et al. The adult human testis transcriptional cell atlas. Cell Res 28, 1141-1157 (2018).
27. Hutton, C. et al. Single-cell analysis defines a pancreatic fibroblast lineage that supports anti-tumor immunity. Cancer Cell 39, 1227-1244.e1220 (2021).
28. Elyada, E. et al. Cross-species single-cell analysis of pancreatic ductal adenocarcinoma reveals antigen-presenting cancer-associated fibroblasts. Cancer Discov. 9, 1102-1123 (2019).
29. Bear, A. S., Vonderheide, R. H. \& O'Hara, M. H. Challenges and opportunities for pancreatic cancer immunotherapy. Cancer Cell 38, 788-802 (2020).
30. Zheng, B. et al. TM4SF1 as a prognostic marker of pancreatic ductal adenocarcinoma is involved in migration and invasion of cancer cells. Int J. Oncol. 47, 490-498 (2015).
31. Fu, F. et al. Role of transmembrane 4 L Six Family 1 in the development and progression of cancer. Front. Mol. Biosci. 7, 202 (2020).
32. Xu, D. et al. Lost miR-141 and upregulated TM4SF1 expressions associate with poor prognosis of pancreatic cancer: regulation of EMT and angiogenesis by miR-141 and TM4SF1 via AKT. Cancer Biol. Ther. 21, 354-363 (2020).
33. Guo, X. et al. Global characterization of T cells in non-small-cell lung cancer by single-cell sequencing. Nat. Med. 24, 978-985 (2018).
34. Zhang, Q. et al. Landscape and dynamics of single immune cells in hepatocellular carcinoma. Cell 179, 829-845.e820 (2019).
35. Cooper, M. D., Lawton, A. R. \& Kincade, P. W. A two-stage model for development of antibody-producing cells. Clin. Exp. Immunol. 11, 143-149 (1972).
36. Simonelli, C. et al. Both CD8+ and CD16+ human T cell clones can provide B cell help for immunoglobulin production. Int. J. Clin. Lab. Res. 22, 36-39 (1992).
37. Kawanishi, H., Saltzman, L. \& Strober, W. Mechanisms regulating IgA class-specific immunoglobulin production in murine gutassociated lymphoid tissues. II. Terminal differentiation of postswitch slgA-bearing Peyer's patch B cells. J. Exp. Med. 158, 649-669 (1983).
38. Snapper, C. M. \& Mond, J. J. Towards a comprehensive view of immunoglobulin class switching. Immunol. Today 14, 15-17 (1993).
39. De Becker, G. et al. Immunoglobulin isotype regulation by antigenpresenting cells in vivo. Eur. J. Immunol. 24, 1523-1528 (1994).
40. Wu, S. Z. et al. A single-cell and spatially resolved atlas of human breast cancers. Nat. Genet. 53, 1334-1347 (2021).
41. Chen, H. et al. Dissecting mammalian spermatogenesis using spatial transcriptomics. Cell Rep. 37, 109915 (2021).
42. Robinson, M. D., McCarthy, D. J. \& Smyth, G. K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140 (2010).
43. Nepusz, G. C. A. T. The igraph software package for complex network research. InterJ. Complex Sys. 1695, 1-9 http://igraph.org (2006).

## Acknowledgements

This work was supported by Changping Laboratory, the National Natural Science Foundation of China ( 32022016 X.R., 92159305 X.R., and 31991171X.R.), National Key R\&D Program of China (2020YFE0202200 X.R. and 2022YFC3400904 X.R.).

## Author contributions

X.R. conceived this study, designed the algorithm, supervised the analysis, and wrote the manuscript. Z.X.Z developed the software, conducted the data analysis, and wrote the manuscript. Y.Z. conducted the data analysis and wrote the manuscript. Z.M.Z provided valuable discussion on the data analysis and wrote the manuscript.

## Competing interests

The authors declare no competing interests.

## Additional information

Supplementary information The online version contains supplementary material available at https://doi.org/10.1038/s41467-023-43600-9.

Correspondence and requests for materials should be addressed to Xianwen Ren.

Peer review information Nature Communications thanks Ken Chen, Jing Su and the other, anonymous, reviewer(s) for their contribution to the peer review of this work. A peer review file is available.

Reprints and permissions information is available at http://www.nature.com/reprints

Publisher's note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.

Open Access This article is licensed under a Creative Commons Attribution 4.0 International License, which permits use, sharing, adaptation, distribution and reproduction in any medium or format, as long as you give appropriate credit to the original author(s) and the source, provide a link to the Creative Commons licence, and indicate if changes were made. The images or other third party material in this article are included in the article's Creative Commons licence, unless indicated otherwise in a credit line to the material. If material is not included in the article's Creative Commons licence and your intended use is not permitted by statutory regulation or exceeds the permitted use, you will need to obtain permission directly from the copyright holder. To view a copy of this licence, visit http://creativecommons.org/ licenses/by/4.0/.
(c) The Author(s) 2023


[^0]:    ¹Changping Laboratory, Yard 28, Science Park Road, Changping District, Beijing, China. ²Biomedical Pioneering Innovation Center (BIOPIC), Peking University, 100871 Beijing, China. ${ }^{3}$ These authors contributed equally: Zixiang Zhou, Yunshan Zhong. e-mail: renxwise@cpl.ac.cn


---

# Paper 3: s41592-022-01480-9

# Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution 

Bin $\mathrm{Li}^{1,7}$, Wen Zhang ${ }^{1,2,7}$, Chuang Guo ${ }^{1,7}$, Hao Xu ${ }^{1,2}$, Longfei Li ${ }^{3}$, Minghao Fang ${ }^{3}$, Yinlei Hu ${ }^{4}$, Xinye Zhang ${ }^{3}$, Xinfeng Yao ${ }^{1}$, Meifang Tang ${ }^{1}$, Ke Liu ${ }^{1}$, Xuetong Zhao ${ }^{5}$, Jun Lin ${ }^{1,2}$, Linzhao Cheng ${ }^{3}$, Falai Chen ${ }^{4}$, Tian Xue ${ }^{3}$ and Kun $\mathbf{Q u}^{1,2,6} \boxtimes$


#### Abstract

Spatial transcriptomics approaches have substantially advanced our capacity to detect the spatial distribution of RNA transcripts in tissues, yet it remains challenging to characterize whole-transcriptome-level data for single cells in space. Addressing this need, researchers have developed integration methods to combine spatial transcriptomic data with single-cell RNA-seq data to predict the spatial distribution of undetected transcripts and/or perform cell type deconvolution of spots in histological sections. However, to date, no independent studies have comparatively analyzed these integration methods to benchmark their performance. Here we present benchmarking of 16 integration methods using 45 paired datasets (comprising both spatial transcriptomics and scRNA-seq data) and 32 simulated datasets. We found that Tangram, gimVI, and SpaGE outperformed other integration methods for predicting the spatial distribution of RNA transcripts, whereas Cell2location, SpatialDWLS, and RCTD are the top-performing methods for the cell type deconvolution of spots. We provide a benchmark pipeline to help researchers select optimal integration methods to process their datasets.


Spatial transcriptomics approaches allow us to detect RNA transcripts in space, and these approaches have been used to investigate the spatial distribution of gene expression in various tissues and organs, including the brain ${ }^{1}$, heart ${ }^{2}$, pancreas ${ }^{3}$, and skin ${ }^{4}$. On the one hand, the spatial transcriptomics approaches based on in situ hybridization and fluorescence microscopy (image-based)including seqFISH ${ }^{5}$, osmFISH ${ }^{6}$, and MERFISH ${ }^{7}$-detect the spatial distribution of transcripts with high resolution and accuracy, but they are limited in the total number of RNA transcripts that they can detect. On the other hand, spatial transcriptomics approaches based on next-generation sequencing (seq-based), such as ST ${ }^{8}$, 10X Visium ${ }^{9}$, and Slide-seq ${ }^{10,11}$, can capture expressed RNAs at the whole-transcriptome scale from spots in space, but each spot (radius $10-100 \mu \mathrm{~m}$ ) may contain multiple cells, which limits the spatial resolution of these approaches. The limitations of these spatial transcriptomics approaches hinder their capacity to capture whole-transcriptome-scale data at single-cell resolution in space.

To break through the limitations of spatial transcriptomics approaches, bioinformaticians have proposed and developed various integration methods to combine spatial transcriptomics and single-cell RNA-seq (scRNA-seq) data. For example, gimVI ${ }^{12}$ employs a deep generative model to infer the likely spatial distribution of undetected transcripts; SpaGE ${ }^{13}$ uses the domain adaptation algorithm PRECISE ${ }^{14}$ and $k$-nearest-neighbor regression to predict the spatial distribution of undetected transcripts; Tangram ${ }^{15}$ uses non-convex optimization and a deep learning framework to learn
a spatial alignment for scRNA-seq data; Seurat ${ }^{16}$ applies canonical correlation analysis ${ }^{17}$ to embed spatial and scRNA-seq data into a common latent space, and projects cells from scRNA-seq data to the spots of the spatial transcriptomics data; LIGER ${ }^{18}$ uses both integrative non-negative matrix factorization ${ }^{19}$ and shared factor neighborhood graphs to predict gene expression levels in space; novoSpaRc ${ }^{20}$ and SpaOTsc ${ }^{21}$ each use optimal transport methods ${ }^{22}$ to construct spatial metrics of cells on the basis of scRNA-seq data; stPlus ${ }^{23}$ combines the auto-encoder and weighted $k$-nearest-neighbor methods to predict spatial gene expression. These integration methods enable researchers to predict the spatial distribution of undetected transcripts.

In addition, Seurat, Tangram, novoSpaRc, and SpaOTsc have the capacity to assign cells from scRNA-seq data to spatial locations in histological sections; this is useful for improving the resolution of the spatial transcriptomics data generated using spatial transcriptomics approaches, like ST or 10 X Visium. Moreover, Cell2location ${ }^{24}$ uses the gene expression signature of the cell subpopulations in scRNA-seq data to estimate the abundance of each cell type at each spot; RCTD ${ }^{25}$ applies cell type profiles learned from scRNA-seq data and supervised learning to decompose cell type mixtures; SpatialDWLS ${ }^{26}$ adopts the weighted-least-squares approach to infer cell type composition; Stereoscope ${ }^{27}$ leverages the model-based probabilistic method and scRNA-seq data to deconvolve the cell mixtures in spatial data; SPOTlight ${ }^{28}$ applies the seeded non-negative matrix factorization for the deconvolution

[^0]of spots; DSTG ${ }^{29}$ deconvolutes spatial transcriptomics data using graph-based convolutional networks; STRIDE ${ }^{30}$ uses the topic profiles trained from scRNA-seq data to decompose cell types from spatial mixtures; DestVI ${ }^{31}$ adopts the variational inference and latent variable models to delineate cell type proportions. These integration methods allow researchers to predict the cell type composition of spots in histological sections.

The emergence of these integration methods has undoubtedly deepened our understanding of spatial transcriptomics data and related biological and pathological processes. However, to the best of our knowledge, no independent study has comprehensively compared the performance of these integration methods for the prediction of the spatial distribution of transcripts or for the cell type deconvolution of spots in histological sections. Here, we used multiple metrics to systematically benchmark the performance of 16 integration methods that can predict the spatial distribution of undetected transcripts, or the cell type composition of spots in histological sections (Fig. 1a), on the basis of processing of 45 paired datasets containing both spatial transcriptomics data and scRNA-seq data and 32 simulated datasets (Fig. 1b). We assessed the accuracy of each integration method in predicting the spatial distribution of transcripts, including for sparse spatial transcriptomics data that were down-sampled from the original datasets. We also evaluated the accuracy of the integration methods for the cell type deconvolution of spots in histological sections on the basis of the simulation of datasets wherein each spot could contain multiple cells of various types. Finally, we evaluated the computational resources consumed by each integration method. Our findings can help researchers choose appropriate integration methods for their datasets, and they raise interesting questions about how various processing and dataset-specific attributes influence the integration performance of these tools for spatial transcriptomics research.

## Results

Benchmarking framework and datasets examined. To evaluate the performance of the 16 integration methods, we collected 45 paired spatial transcriptomics and scRNA-seq datasets from published studies ${ }^{4-7,10,15,31-61}$ (Fig. 1 and Supplementary Table 1). The spatial transcriptomic datasets were produced by 13 spatial transcriptomics approaches, including FISH, osmFISH, seqFISH, MERFISH, STARmap, ISS, EXseq, BaristaSeq, ST, 10X Visium, Slide-seq, Seq-scope, and HDST, and the scRNA-seq datasets were obtained by Drop-seq ${ }^{62}$, Smart-seq ${ }^{63}$, and the 10X Chromium platform ${ }^{64}$. We designed a pipeline to evaluate the performance of the integration methods for combining spatial and single-cell transcriptomics datasets (Fig. 1a). During preprocessing of the scRNA-seq datasets, we removed cells with fewer than 200 RNAs. For the spatial transcriptomic datasets, we generated a 'ground truth' using 2 criteria: for samples with $<1,000$ detected RNAs, we used all of the RNAs; for samples with $>1,000$ detected RNAs, a set of 1,000 highly variable RNAs (assessed on the basis of the coefficient of variation of each RNA; Methods) was used.

In addition, we adopted the algorithms proposed by RCTD and Stereoscope and generated 32 simulated 10X Visium datasets from 16 paired scRNA-seq datasets (Supplementary Tables 2 and 3). A simulated spot contains $5-15$ cells randomly sampled from the scRNA-seq datasets (Methods), and the gene expression values of each spot represent the sum of all the cells in that spot.

After collecting the datasets, we first assessed the performance of eight integration methods, including Tangram, gimVI, SpaGE, Seurat, SpaOTsc, novoSpaRc, LIGER, and stPlus, in predicting the spatial distribution of RNA transcripts that remain undetected in spatial transcriptomics datasets. We used the 45 collected paired datasets to evaluate the accuracy of these integration methods for predicting the RNA spatial distribution. Then we down-sampled
the spatial transcriptomics data to test the performance of the integration methods for datasets with sparse expression matrices.

Beyond the prediction of the spatial distribution of RNA transcripts, Tangram, Seurat, SpaOTsc, and novoSpaRc can assign cells from scRNA-seq data to spatial locations in histological sections. Also, Cell2location, SpatialDWLS, RCTD, Stereoscope, DestVI, STRIDE, SPOTlight, and DSTG can be used to predict the cell type composition of spots in histological sections by combining spatial transcriptomics data and scRNA-seq data. All 12 of these integration methods are capable of deconvoluting cell types of the spots in the spatial transcriptomics datasets that were generated using the 10 X Visium or ST platforms. To compare the performance of these integration methods in cell type deconvolution, we used datasets 4 and 10 as the basis to simulate 'grids' representing low-spatial-resolution datasets, and we simulated 32 datasets from the scRNA-seq data as the ground truth (Methods). Briefly, in the simulated low-resolution datasets, each gridded 'spot' contains 1-18 cells, similar to the spatial transcriptome datasets generated by the 10X Visium or ST approaches. Finally, we assessed the computational resources consumed by each integration method.

Methods predicting spatial distribution of RNA transcripts. We used tenfold crossvalidation (Methods) on the 45 paired datasets to evaluate the accuracy of each integration method in predicting the spatial distribution of RNA transcripts. We quantified the prediction performance of each integration method by calculating the Pearson correlation coefficient (PCC) between the expression vector of a gene in the ground truth of the spatial transcriptomics dataset and the expression vector for the same gene in the result predicted by each integration method (Methods). We first examined the prediction results of the spatial distribution for known marker genes. For example, Lein et al. reported that Igsf21 and Rprm are highly expressed in the L5/L6 layers of the cortex ${ }^{34}$. Compared with the ground truth for dataset 4 (seqFISH+; Smart-seq; mouse cortex), Tangram performed the best in predicting the spatial distribution of Igsf21 (PCC = 0.79), and gimVI, SpaGE, and Seurat followed closely behind (PCC $=0.77,0.71$, and 0.70 ) (Fig. 2a). For the spatial distribution of Rprm, the results generated by SpaGE and Seurat had the highest PCC values (PCC $=0.79$ ), followed by SpaOTsc, gimVI, Tangram, and LIGER ( $\mathrm{PCC}=0.78,0.71,0.66,0.65$ ) (Fig. 2b).

We also examined the predicted results of the spatial distribution of COL17A1 in dataset 42 (ST; 10X Chromium; human squamous carcinoma). COL17A1 is a known marker gene for basal cells of squamous carcinoma ${ }^{4}$. Tangram, gimVI, novoSpaRc, and SpaGE successfully predicted that COL17A1 was highly expressed in the basal cells of squamous carcinoma; notably, the PCC values of these four integration methods were, respectively, 0.86 (Tangram), 0.84 (gimVI), 0.76 (novoSpaRc), and 0.70 (SpaGE), higher than the best result of the other integration methods (Seurat, 0.48 ; SpaOTsc, 0.40 ; LIGER, 0.31; stPlus, 0.27) (Extended Data Fig. 1a).

To further quantify the prediction accuracy of each integration method, we adopted three metrics besides PCC: (1) structural similarity index (SSIM), which combines mean value, variance, and covariance to measure the similarity between the predicted result and the ground truth; (2) root mean square error (RMSE), the absolute error between the predicted distribution and the ground truth; and (3) Jensen-Shannon divergence (JS), which uses relative information entropy to gauge the difference between two distributions. For one gene, a higher PCC/SSIM or lower RMSE/JS value indicates better prediction accuracy. We also defined an accuracy score (AS) by aggregating the four metrics (Methods) to simplify the evaluation of the accuracy of each integration method (a higher AS value indicates better performance).

Taking dataset 4 (seqFISH+; Smart-seq; mouse cortex) as an example of the image-based spatial transcriptomics approaches, Tangram, gimVI, and SpaGE clearly outperformed the other

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-03.jpg?height=1935&width=1750&top_left_y=197&top_left_x=158)
Fig. 1 | Benchmarking workflow and summary characteristics of the examined paired datasets. a, Schematic overview of the benchmarking workflow used to compare the performance of the integration methods for paired spatial transcriptomics and scRNA-seq datasets. We used the 16 integration methods to combine the spatial and single-cell transcriptomics data, and then compared their performance for (1) predicting the spatial distribution of RNA transcripts and (2) deconvoluting cell types of a histological spot. We also assessed the computational resources consumed by the integration methods. $\mathbf{b}$, Information for the 45 paired datasets and 32 simulated datasets used in this study: each dataset contains both spatial transcriptomic data and scRNA-seq data for the same tissue. Detailed information for the data source is presented in the Methods and in Supplementary Table 1.

integration methods. Specifically, we found that the average PCC/ SSIM of Tangram, gimVI, and SpaGE were $0.54 / 0.45,0.52 / 0.43$, and 0.49/0.39, higher than the PCC/SSIM values for the other 5 methods,
and the average RMSE/JS of these three methods were 0.94/0.18, 0.97/0.19, and 0.99/0.21, lower than the average RMSE/JS for the others (Fig. 2c). Moreover, the average AS for the Tangram, gimVI,

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-04.jpg?height=2055&width=1739&top_left_y=197&top_left_x=154)
Fig. 2 | Comparing the accuracy of eight integration methods capable of predicting the spatial distribution of RNA transcripts. $\mathbf{a}, \mathbf{b}$, The spatial distribution of Igsf21(a) and Rprm (b) in dataset 4 (seqFISH+; Smart-seq; mouse cortex), including the ground truth and the predicted result from each of the integration methods. PCC, Pearson correlation coefficient between the expression vector of a transcript in the ground truth and that of the predicted result. $\mathbf{c}$, The bar plots of PCC, SSIM, RMSE, and JS of each integration method in predicting the spatial distribution of transcripts in dataset 4 . Data are presented as mean values $\pm 95 \%$ confidence intervals; $n=1,000$ predicted genes. $\mathbf{d}$, The violin plot of AS (which is aggregated from the PCC, SSIM, RMSE, and JS values; see Methods) of the 8 integration methods for transcripts in dataset 4 . Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range; $n=4$ benchmark metrics. $\mathbf{e}$, Boxplots of AS of the 8 integration methods for all 45 paired datasets. Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range; $n=45$ independent datasets.

and SpaGE predictions were $1.0,0.875$, and 0.75 , higher than that of the other methods (Fig. 2d). We also calculated the PCC, SSIM, RMSE, JS, and AS values of the prediction results for all transcripts in dataset 42 (as an example of the seq-based spatial transcriptomics approaches), and found that Tangram and gimVI outperformed the other integration methods on the basis of these metrics (Extended Data Fig. 1b,c).

To systematically assess the accuracy of the eight integration methods' predictions of the spatial distribution of undetected transcripts, we determined the PCC, SSIM, RMSE, JS, and AS values of their prediction results for all 45 paired datasets (Fig. 2e, Extended Data Fig. 2). The average ASs for the Tangram, gimVI, and SpaGE predictions were $0.96,0.84$, and 0.69 , respectively, all of which exceed the AS values for Seurat (0.50), SpaOTsc (0.55), LIGER (0.25), novoSpaRc (0.47), and stPlus (0.31). Note that Tangram was still the best-performing integration method when we separately assessed the image-based datasets, the seq-based datasets, and the 32 simulated datasets, followed by gimVI and SpaGE (Extended Data Fig. 3a-c). Because 10X Visium, seqFISH, MERFISH, and Slide-seq have released more than 3 datasets, we further compared the ASs of the eight integration methods when processing data obtained using these four spatial transcriptomics technologies (Extended Data Fig. 3d-g). We found that Tangram, gimVI, and SpaGE outperformed other integration methods for data generated from 10X Visium, seqFISH, and MERFISH platforms, and Tangram and gimVI are top-ranked methods in processing Slide-seq datasets.

Several integration methods (for example, Seurat, LIGER, SpaGE, and stPlus) normalized the spatial transcriptomics data by default prior to integration. Here, we tested four schemes of input expression matrices: (1) raw expression matrix of spatial data and raw expression matrix of scRNA-seq data (R-R); (2) normalized expression matrix of spatial data and raw expression matrix of scRNA-seq data (N-R); (3) raw expression matrix of spatial data and normalized expression matrix of scRNA-seq data (R-N); and (4) normalized expression matrix of spatial data and normalized expression matrix of scRNA-seq data (N-N).

Interestingly, for 28 paired seq-based datasets, the transcript spatial distributions generated by Tangram, gimVI, SpaGE, Seurat, SpaOTsc, and LIGER have significantly higher PCC values when using an R-R and R-N input scheme than when using an N-R or N-N input scheme, and this trend was observed for 16 of the 28 paired datasets ( $P$ values $<0.01$, paired $t$-test) (Extended Data Figs. 4 and 5a); for SpaGE, Seurat, SpaOTsc, and novoSpaRc, the PCC values of the results with the R-R input scheme were higher than those with the other input schemes in 19 of the 28 paired datasets ( $P$ values $<0.01$ ), and stPlus generated results with higher PCC values when using the R-R input scheme than the N-N input scheme in 18 of the 28 paired datasets ( $P$ value $<0.05$ ). For 15 paired image-based datasets (Extended Data Figs. 4 and 5b), the transcript spatial distributions generated by Tangram, gimVI, SpaGE, and Seurat have higher PCC values when using the R-R or R-N input scheme than when using the N-R or N-N input scheme (11 out of 15 datasets, $P$ value $<0.05$ ); SpaGE, Seurat, and LIGER have higher PCC values when using the R-R input scheme than when using the other input schemes ( 11 out of 15 datasets, $P$ value $<0.05$ ); SpaOTsc has a higher PCC value when using the R-R input schemes than when using the NN input scheme (12 out of 15 datasets, $P$ value < 0.05 ). Nevertheless, it should be emphasized that regardless of what input scheme was used, Tangram invariably outperformed the other integration methods (Extended Data Fig. 5c-f).

Impact of matrix sparsity. Notably, for datasets $12,13,40$, and 44 , all eight integration methods had low accuracy in predicting the spatial distribution of transcripts (that is, average PCC/SSIM < 0.3, Extended Data Fig. 2). We investigated this apparently poor performance of the integration methods for these datasets by calculating
correlation coefficients between the four metrics (PCC, SSIM, RMSE, and JS) and considered several features of the spatial transcriptomics datasets, including the sparsity of the expression matrix (the sparsity of the spatial transcriptomics and scRNA-seq data is defined as the percentage of zero elements in the expression matrix), the number of detected genes, the number of detected spots, and the number of genes per spot. Ultimately, we found that the JS values of the 8 methods all linearly increased along with the rising of the sparsity of expression matrices ( $P$ values $<1 \times 10^{-6}$, coefficient of determination $\left(R^{2}\right) \geq 0.50$ ) (Extended Data Fig. 6).

To further characterize the impact of matrix sparsity, we next evaluated the performance of each integration method when inputting a very sparse spatial expression matrix (down-sampled from high quality datasets where sparsity was lower than 0.7 ). Specifically, we examined spatial transcriptomics datasets that captured $>1,000$ genes from $>100$ spots as high quality. To simulate expression matrices with 'high sparsity,' we adopted Splatter ${ }^{65}$ and Scuttle ${ }^{66}$ to down-sample the non-zero elements from the original expression matrices to varying extents (Methods). We then used the original and down-sampled expression matrices as the inputs for the eight integration methods.

First, we evaluated the impact of the expression matrix sparsity in predicting the spatial distribution of known marker genes (Methods). Drew et al. reported that Cplx1 is highly expressed in layer L5 of the cortex ${ }^{67}$. Examining Cplx1 in both the original and down-sampled data (down-sampling rate $=0.8$ ) of dataset 4 , we observed that the spatial distributions of Cplx1 predicted by Tangram, gimVI, and SpaGE each had PCC values $>0.7$ for both the original and down-sampled data (Extended Data Fig. 7a).

We then assessed the performance of each integration method by counting the proportion of transcripts in a dataset exceeding a PCC threshold of 0.5 for both the original and down-sampled data, which we deemed the 'robustness score' (RS). For dataset 4, Tangram had the highest RS value (0.60), followed by gimVI (0.55) and then SpaGE (0.51) (Fig. 3a). Moreover, we noted that (1) the RS values decreased as the down-sampling rate increased and (2) the RS values for Tangram, gimVI, and SpaGE were consistently higher than those for Seurat, SpaOTsc, LIGER, novoSpaRc, and stPlus (Fig. 3b). A combined analysis which included the down-sampled data from 19 datasets again highlighted the strong performance of Tangram, gimVI, and SpaGE: even when the down-sampling rate reached 0.8 , the average RS values of these three methods remained $>0.50$ (Fig. 3c, down-sampled by Splatter; Extended Data Fig. 7b,c, down-sampled by Scuttle). In summary, Tangram, gimVI, and SpaGE outperformed other integration methods in predicting the spatial distribution of transcripts for highly sparse datasets.

Performance of methods in cell type deconvolution. A common issue encountered when using spatial transcriptomics approaches like 10X Visium and ST is that each spot from a histological section may contain multiple cells, so it can be impossible to correctly assign the cell type composition of each spot. As noted above, Seurat, SpaOTsc, Tangram, and novoSpaRc are capable of assigning each cell from a scRNA-seq analysis to a spot from a spatial transcriptomics analysis, implying that they can be used to deconvolute the cell types of each spot. Moreover, Cell2location, SpatialDWLS, RCTD, Stereoscope, DestVI, STRIDE, SPOTlight, and DSTG were also designed for this purpose.

To compare the performance of the 12 integration methods in predicting the cell type composition of spots, we simulated this ‘multi-cell spot problem’ experienced with ST and 10X Visium datasets by 'gridding' a dataset that did not have this problem (dataset 10, acquired using STARmap; Smart-seq; mouse visual cortex). The cell type composition of each spot in dataset 10 has been reported and can be used as the ground truth when simulating a dataset with potentially ambiguous cell type assignations in each spot (Fig. 4a and

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-06.jpg?height=1059&width=1793&top_left_y=197&top_left_x=141)
Fig. 3 | Comparing the accuracy of the eight integration methods for sparse spatial expression matrices down-sampled from the original datasets
using Splatter. a, PCC of the spatial distribution of transcripts predicted from the original data and down-sampled data from dataset 4. The PCC values of the red-colored transcripts are greater than 0.5 for both the original and the down-sampled data. The proportion of the red-colored transcripts in all transcripts was defined as the RS. $\mathbf{b}$, RS values of the 8 integration methods when processing sparse expression matrices down-sampled from dataset 4 at different down-sampling rates. c, RS values of the eight integration methods when processing the sparse expression matrices of the down-sampled datasets. The original datasets (used to generate the down-sampled datasets) capture $>1,000$ genes from $>100$ spots, and the sparsity of the expression matrices is $<0.7$. Data are presented as mean values $\pm 95 \%$ confidence intervals; $n=19$ independent datasets.

Methods). The original dataset 10 captured 1,549 cells, corresponding to 15 cell types. After gridding, the simulated data had 189 spots, with each spot containing 1-18 cells. We plotted the locations of L4 excitatory neurons and found that RCTD and Stereoscope performed better in terms of the PCC values (0.87), followed by Tangram (0.85), Cell2loacation (0.83), STRIDE (0.80), SPOTlight (0.79), Seurat (0.76), SpaOTsc (0.74), and DSTG (0.71) (Fig. 4b). We then employed PCC, SSIM, RMSE, JS, and AS metrics to quantify the accuracy of the 12 integration methods in predicting the cell type composition of spots in gridded dataset 10 (Fig. 4c and Extended Data Fig. 8a). RCTD had the highest AS score (0.94), followed by Stereoscope (0.92).

We also performed the same analysis on dataset 4 (seqFISH+; Smart-seq; mouse cortex), which contains 524 cells of 14 cell types. After 'gridding', the simulated dataset had 72 spots (Extended Data Fig. 8b). Using the ground truth of the locations for the L5/6 excitatory neurons, we found that SpatialDWLS, RCTD, Tangram, Cell2location, and Stereoscope had PCC values of 0.88 , 0.86, 0.85, 0.83, and 0.81 for the assignations of the L5/6 excitatory neurons, higher than other integration methods (Extended Data Fig. 8c). Moreover, in the prediction results for all cell types of dataset 4 , SpatialDWLS, Tangram, and RCTD had the top 1, 2, and 3 ranking AS values ( $1.0,0.92$, and 0.83 respectively), followed by Cell2location (0.67) and Stereoscope (0.65) (Fig. 4d).

We further quantified the performance of these integration methods in cell type deconvolution of spots in the 32 simulated datasets that were synthesized from scRNA-seq datasets (Supplementary Table 2 and 3 and Methods). As the cell type information of each cell
in these scRNA-seq datasets has been reported by the data source papers, the cell type composition of a simulated spot can be inferred from the cells it contains. Note that novoSpaRc and SpaOTsc require spatial location information for each spot and thereby were excluded because spatial location information was not available in the simulated datasets. We used the 32 simulated datasets as the ground truth to assess the performance of the remaining 10 integration methods (including Seurat, Tangram, Cell2location, SpatialDWLS, RCTD, Stereoscope, DestVI, STRIDE, SPOTlight, and DSTG) in deconvoluting cell types in spots. We found that the average PCC and SSIM values of Cell2location, SpatialDWLS, and STRIDE are 0.83/0.75, $0.78 / 0.71$, and $0.83 / 0.69$, higher than those of the other integration methods, and the average RMSE and JS values of the Cell2location, SpatialDWLS, RCTD, and STRIDE are 0.08/0.33, 0.10/0.32, 0.096/0.37, and $0.11 / 0.37$, lower than those of the other integration methods (Extended Data Fig. 8e). We also used the aggregation of the four metrics (that is, the AS score) to rank the performance of these integration methods in predicting cell type composition of spots, and we found that Cell2location, SpatialDWLS, RCTD, and STRIDE outperformed the other integration methods (Fig. 4e).

Computational resources. We used all the 45 paired datasets to compare the computational resources consumed by the 8 integration methods that can predict the spatial distribution of undetected transcripts (Supplementary Table 4). We used an identical CPU platform ( $2.2 \mathrm{GHz}, 45 \mathrm{MB}$ L3 cache, 144 CPU cores) to test each method. We are aware that gimVI and Tangram can support GPU

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-07.jpg?height=1606&width=1795&top_left_y=193&top_left_x=141)
Fig. 4 | Comparing the performance of the 12 integration methods capable of deconvoluting cell types of each histological spot. a, A STARmap slide of dataset 10 (STARmap; Smart-seq; mouse visual cortex), with cells annotated by cell types. Each grid represents a simulated spot containing multiple cells. L1, $2 / 3,4,5$, and 6: layer 1, 2/3, 4, 5, and 6; Excitatory L2/3, 4, 5, and 6: excitatory neurons in layer $2 / 3,4,5$, and 6; Npy, Pvalb, Sst, and Vip: GABAergic interneuron subtypes marked by Npy, Pvalb, Sst, and Vip; Astro: astrocytes; Endo: endothelia cells; Micro: microglia; Oligo: oligodendrocytes; Smc: smooth muscle cells; Other: other unclassified cells. $\mathbf{b}$, The proportion of L4 excitatory neurons in the spots simulated from dataset 10, including the ground truth and the predicted results of 12 integration methods. c, d, Bar plots of AS (aggregated from PCC, SSIM, RMSE, and JS; see Methods) of the cell type composition of the histological spots simulated from dataset $10(\mathbf{c})$ and dataset $4(\mathbf{d})$, predicted by 12 integration methods. Data are presented as mean values $\pm 95 \%$ confidence intervals; $\mathrm{n}=4$ benchmark metrics. $\mathbf{e}$, Boxplots of AS of the 10 integration methods for all the 32 simulated datasets. SpaOTsc and novoSpaRc are excluded, as they require spatial location information for each spot, which is not available in the simulated datasets. Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range; $n=32$ independent datasets.

processing; however, these two integration methods reported memory errors on our GPU platform (NVIDIA Tesla K80 with 12 GB memory) when processing the largest dataset $40(19,522$ spots in the spatial transcriptomics data and 26,252 cells in the scRNA-seq data). Notably, it took Seurat and LIGER less than 10 minutes of CPU time to process each dataset, and Tangram and LIGER consumed less than 32 GB of memory.

We then assessed the impacts of various data attributes (including the number of cells in scRNA-seq data, the number of spots in spatial
data, and the number of genes used for training) on the computational resources consumed by those eight integration methods. By down-sampling the number of cells and the number of spots in dataset 40 and the number of training genes in dataset 6, we found that Seurat was invariably the most computationally efficient method among the 8 integration methods for the prediction of the spatial distribution of undetected transcripts (Extended Data Fig. 9a-c).

To compare the computational resources consumed by the 10 integration methods that can deconvolute the cell types of spots,
we used a large simulated dataset (Methods) that contains 10,000 cells, 20,000 spots, and 56 cell types. For this dataset, Cell2locations reported memory errors on our GPU platform. Seurat and Tangram took less than 30 minutes of CPU time, and Stereoscope, Tangram, and DestVI consumed less than 8 GB of memory (Extended Data Fig. 9d). We then evaluated the impacts of the number of cells in scRNA-seq data, the number of spots in spatial data, and the number of cell types on computing time consumed by those 10 integration methods, and found that Tangram and Seurat are the top two most-efficient methods for processing cell type deconvolution of spots (Extended Data Fig. 9e-g).

## Discussion

In this study, we benchmarked the performance of 16 integration methods capable of combining spatial transcriptomics data and single-cell transcriptomics data. We found that Tangram, gimVI, and SpaGE outperformed other integration methods for predicting the spatial distribution of transcripts, whereas Cell2location, SpatialDWLS, and RCTD were superior to other integration methods for cell type deconvolution of spots in histological sections. Our study helps researchers to choose appropriate tools and to optimize data-analysis workflows to accurately and efficiently integrate spatial transcriptomics data with scRNA-seq data. We have also provided a benchmark pipeline (https://github.com/QuKunLab/ SpatialBenchmarking) and an instructive table (Supplementary Table 5) summarizing the properties and performance of all the benchmarked methods to guide researchers select suitable tools that match their data combinations.

Methods constructed on the basis of probabilistic models combined with negative binomial or Poisson distributions, such as gimVI, Cell2location, and RCTD, generally perform better at predicting the spatial distribution of transcripts or deconvolving cell types of spots. A deep learning algorithm was also applied in several integration methods, among which Tangram is one of the best-performing methods in predicting spatial distribution of the undetected transcripts. Technically, Tangram employs non-convex optimization in the model and selects only the optimal subset of scRNA-seq observations in the loss function. A combination of these measures may help improve the predictive power of the tools.

One observation from our comparative analysis is that the sparsity of the spatial transcriptomics expression matrix seriously affects the performance of the eight integration methods that predict the spatial distribution of RNA transcripts. There are multiple tactics that can be used to combat this sparsity issue for spatial transcriptomics expression matrices: researchers can increase the depth of sequencing, screen spots and genes with strict cut-off values to reduce the sparsity of the filtered expression matrix, or consider applying imputation algorithms (for example, SAVER ${ }^{68}$, MAGIC ${ }^{69}$, and WEDGE ${ }^{70}$ ) to impute the zero elements in the expression matrix.

Another potential application of spatial transcriptomics is to predict ligand-receptor interactions between two cell types that are spatially close to each other. Many analytical tools have been developed for this task, such as SpaOTsc ${ }^{21}$, Giotto ${ }^{71}$, CellChat ${ }^{72}$, NicheNet ${ }^{73}$, ICELLNET ${ }^{74}$, and SingleCellSignalR ${ }^{75}$. However, the vast discrepancies in the results from different methods make informative comparison difficult. For instance, only a small proportion ( $<5 \%$ ) of the predicted ligand-receptor interactions were shared in the results for >3 methods (https://github.com/QuKunLab/SpatialBenchmarking/ tree/main/FigureData). Benchmarking analysis may thereby rely on more experimental validation of cell-cell interactions in the future.

The wide diversity and fast-moving state of scRNA-seq technologies (Drop-seq, Smart-seq, and 10X Chromium) and spatial sequencing technologies (FISH, osmFISH, seqFISH, MERFISH, STARmap, ISS, EXseq, BaristaSeq, ST, 10X Visium, Slide-seq, Seq-scope, and HDST) complicates the task of this benchmarking analysis. We adopted three distinct perspectives to overcome the
challenges caused by the diversity of sequencing technologies: (1) we divided the spatial transcriptome datasets into two categories (seq-based technologies and image-based technologies); (2) we added an independent comparison of the integration methods for the set of four spatial transcriptomics technologies (MERFISH, seqFISH, Slide-seq, and 10X Visium) that have to date released more than three datasets; and (3) we used several intrinsic parameters (for example, the number of captured genes, the number of captured spots, and the sparsity of the expression matrix) to characterize the datasets generated by different spatial transcriptomics technologies. There are still aspects (that is, different numbers of genes and different spatial organization) that may affect the performance of the integration methods and the user's expectations. Nevertheless, on the basis of the current collection of datasets, we found that the performance rankings of these integration methods are barely affected by the spatial transcriptome technologies that generated these data.

Advances in spatial transcriptomics technologies, such as new versions of 10X Visium and BGI Stereo-seq ${ }^{76}$, may enable detection of transcriptomic information for spots with diameters much smaller than cell size. However, in the near future, each spot still might not correspond exactly to a single cell when using these technologies. Moreover, considering that there is a large amount of publicly available spatial transcriptomics datasets that are highly valuable to different research communities, and the intense research efforts using spatial transcriptomics technologies that are now underway, we contend that there will be strong interest in integrating spatial transcriptomics and scRNA-seq data to determine the spatial distribution of cells or undetected transcripts.

In summary, this study presents a much-needed independent comparison of available integration methods of spatial transcriptomics and scRNA-seq data to determine the cell type deconvolution or the spatial distribution of undetected genes. Our results are also a useful resource for biologists who want to analyze their spatial transcriptomics data or methods developers who want to improve state-of-the-art technology.

## Online content

Any methods, additional references, Nature Research reporting summaries, source data, extended data, supplementary information, acknowledgements, peer review information; details of author contributions and competing interests; and statements of data and code availability are available at https://doi.org/10.1038/s41592-022-01480-9.

Received: 10 June 2021; Accepted: 30 March 2022;
Published online: 16 May 2022

## References

1. Maynard, K. R. et al. Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex. Nat. Neurosci. 24, 425-436 (2021).
2. Asp, M. et al. A spatiotemporal organ-wide gene expression and cell atlas of the developing human heart. Cell 179, 1647-1660 e1619 (2019).
3. Moncada, R. et al. Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas. Nat. Biotechnol. 38, 333-342 (2020).
4. Ji, A. L. et al. Multimodal analysis of composition and spatial architecture in human squamous cell carcinoma. Cell 182, 497-514 e422 (2020).
5. Eng, C. L. et al. Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH. Nature 568, 235-239 (2019).
6. Codeluppi, S. et al. Spatial organization of the somatosensory cortex revealed by osmFISH. Nat. Methods 15, 932-935 (2018).
7. Moffitt, J. R. et al. Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. Science 362, eaau5324 (2018).
8. Stahl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353, 78-82 (2016).
9. Visium spatial gene expression (10x Genomics, 2020).
10. Stickels, R. R. et al. Highly sensitive spatial transcriptomics at near-cellular resolution with Slide-seqV2. Nat. Biotechnol. 39, 313-319 (2021).
11. Rodriques, S. G. et al. Slide-seq: a scalable technology for measuring genome-wide expression at high spatial resolution. Science 363, 1463-1467 (2019).
12. Lopez, R. et al. A joint model of unpaired data from scRNA-seq and spatial transcriptomics for imputing missing gene expression measurements. ICML Workshop on Computational Biology (2019).
13. Abdelaal, T., Mourragui, S., Mahfouz, A. \& Reinders, M. J. T. SpaGE: spatial gene enhancement using scRNA-seq. Nucleic Acids Res. 48, e107 (2020).
14. Mourragui, S., Loog, M., van de Wiel, M. A., Reinders, M. J. T. \& Wessels, L. F. A. PRECISE: a domain adaptation approach to transfer predictors of drug response from pre-clinical models to tumors. Bioinformatics 35, i510-i519 (2019).
15. Biancalani, T. et al. Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nat. Methods 18, 1352-1362 (2021).
16. Stuart, T. et al. Comprehensive integration of single-cell data. Cell 177, 1888-1902 e1821 (2019).
17. Butler, A., Hoffman, P., Smibert, P., Papalexi, E. \& Satija, R. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat. Biotechnol. 36, 411-420 (2018).
18. Welch, J. D. et al. Single-Cell Multi-omic integration compares and contrasts features of brain cell identity. Cell 177, 1873-1887 e1817 (2019).
19. Yang, Z. \& Michailidis, G. A non-negative matrix factorization method for detecting modules in heterogeneous omics multi-modal data. Bioinformatics 32, 1-8 (2016).
20. Nitzan, M., Karaiskos, N., Friedman, N. \& Rajewsky, N. Gene expression cartography. Nature 576, 132-137 (2019).
21. Cang, Z. \& Nie, Q. Inferring spatial and signaling relationships between cells from single cell transcriptomic data. Nat. Commun. 11, 2084 (2020).
22. Villani, C. Optimal Transport: Old and New Vol. 338 (Springer, 2009).
23. Chen, S. Q., Zhang, B. H., Chen, X. Y., Zhang, X. G. \& Jiang, R. stPlus: a reference-based method for the accurate enhancement of spatial transcriptomics. Bioinformatics 37, I299-I307 (2021).
24. Kleshchevnikov, V. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat Biotechnol. 1-11, https://doi.org/10.1038/s41587-021-01139-4 (2022).
25. Cable, D. M. et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nat. Biotechnol. 40, 517-526 (2021).
26. Dong, R. \& Yuan, G. C. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biol. 22, 145 (2021).
27. Andersson, A. et al. Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography. Commun. Biol. 3, 565 (2020).
28. Elosua-Bayes, M., Nieto, P., Mereu, E., Gut, I. \& Heyn, H. SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Res. 49, e50 (2021).
29. Song, Q. Q. \& Su, J. DSTG: deconvoluting spatial transcriptomics data through graph-based artificial intelligence. Brief. Bioinform. 22, bbaa414 (2021).
30. Sun, D., Liu, Z., Li, T., Wu, Q. \& Wang, C. STRIDE: accurately decomposing and integrating spatial transcriptomics using single-cell RNA sequencing. Nucleic Acids Res. gkac150 (2022).
31. Lopez, R. et al. Multi-resolution deconvolution of spatial transcriptomics data reveals continuous patterns of inflammation. Nat. Biotechnol. in press (2022).
32. Karaiskos, N. et al. The Drosophila embryo at single-cell transcriptome resolution. Science 358, 194-199 (2017).
33. Berkeley Drosophila Transcription Network Project. http://bdtnp.lbl.gov:8080/ Fly-Net/.
34. Tasic, B. et al. Shared and distinct transcriptomic cell types across neocortical areas. Nature 563, 72-78 (2018).
35. Shah, S., Lubeck, E., Zhou, W. \& Cai, L. In situ transcription profiling of single cells reveals spatial organization of cells in the mouse hippocampus. Neuron 92, 342-357 (2016).
36. Xia, C., Fan, J., Emanuel, G., Hao, J. \& Zhuang, X. Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression. Proc. Natl Acad. Sci. USA 116, 19490-19499 (2019).
37. Wang, X. et al. Three-dimensional intact-tissue sequencing of single-cell transcriptional states. Science 361, eaat5691 (2018).
38. Joglekar, A. et al. A spatially resolved brain region- and cell type-specific isoform atlas of the postnatal mouse brain. Nat. Commun. 12, 463 (2021).
39. Navarro, J. F. et al. Spatial transcriptomics reveals genes associated with dysregulated mitochondrial functions and stress signaling in alzheimer disease. iScience 23, 101556 (2020).
40. Lohoff, T. et al. Integration of spatial and single-cell transcriptomic data elucidates mouse organogenesis. Nat. Biotechnol. 40, 74-85 (2022).
41. Nowotschin, S. et al. The emergent landscape of the mouse gut endoderm at single-cell resolution. Nature 569, 361-367 (2019).
42. Han, X. et al. Mapping the mouse cell atlas by microwell-Seq. Cell 172, 1091-1107 e1017 (2018).
43. Pijuan-Sala, B. et al. A single-cell molecular map of mouse gastrulation and early organogenesis. Nature 566, 490-495 (2019).
44. Brann, D. H. et al. Non-neuronal expression of SARS-CoV-2 entry genes in the olfactory system suggests mechanisms underlying COVID-19-associated anosmia. Sci. Adv. 6, eabc5801 (2020).
45. Cho, C. S. et al. Microscopic examination of spatial transcriptome using Seq-Scope. Cell 184, 3559-3572 e3522 (2021).
46. Tabula Muris, C. et al. Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature 562, 367-372 (2018).
47. Vickovic, S. et al. High-definition spatial transcriptomics for in situ tissue profiling. Nat. Methods 16, 987-990 (2019).
48. Saunders, A. et al. Molecular diversity and specializations among the cells of the adult mouse brain. Cell 174, 1015-1030 e1016 (2018).
49. McCray, T. et al. Erratum: Vitamin D sufficiency enhances differentiation of patient-derived prostate epithelial organoids. iScience 24, 102640 (2021).
50. Janosevic, D. et al. The orchestrated cellular and molecular responses of the kidney to endotoxin define a precise sepsis timeline. eLife 10, e62270 (2021).
51. Melo Ferreira, R. et al. Integration of spatial and single-cell transcriptomics localizes epithelial cell-immune cross-talk in kidney injury. JCI Insight 6, e147703 (2021).
52. Sanchez-Ferras, O. et al. A coordinated progression of progenitor cell states initiates urinary tract development. Nat. Commun. 12, 2627 (2021).
53. Wu, S. Z. et al. A single-cell and spatially resolved atlas of human breast cancers. Nat. Genet. 53, 1334-1347 (2021).
54. Alon, S. et al. Expansion sequencing: spatially precise in situ transcriptomics in intact biological systems. Science 371, 481 (2021).
55. Chen, X., Sun, Y. C., Church, G. M., Lee, J. H. \& Zador, A. M. Efficient in situ barcode sequencing using padlock probe-based BaristaSeq. Nucleic Acids Res. 46, e22 (2018).
56. Booeshaghi, A. S. et al. Isoform cell-type specificity in the mouse primary motor cortex. Nature 598, 195-199 (2021).
57. Hodge, R. D. et al. Conserved cell types with divergent features in human versus mouse cortex. Nature 573, 61-68 (2019).
58. Tepe, B. et al. Single-cell RNA-seq of mouse olfactory bulb reveals cellular heterogeneity and activity-dependent molecular census of adult-born neurons. Cell Rep. 25, 2689-2703 e2683 (2018).
59. Hunter, M. V., Moncada, R., Weiss, J. M., Yanai, I. \& White, R. M. Spatially resolved transcriptomics reveals the architecture of the tumor-microenvironment interface. Nat. Commun. 12, 6278 (2021).
60. McKellar, D. W. et al. Large-scale integration of single-cell transcriptomic data captures transitional progenitor states in mouse skeletal muscle regeneration. Commun. Biol. 4, 1280 (2021).
61. Ratz, M. et al. Clonal relations in the mouse brain revealed by single-cell and spatial transcriptomics. Nat. Neurosci. 25, 285-294 (2022).
62. Macosko, E. Z. et al. Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. Cell 161, 1202-1214 (2015).
63. Ramskold, D. et al. Full-length mRNA-seq from single-cell levels of RNA and individual circulating tumor cells. Nat. Biotechnol. 30, 777-782(2012).
64. Zheng, G. X. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 (2017).
65. Zappia, L., Phipson, B. \& Oshlack, A. Splatter: simulation of single-cell RNA sequencing data. Genome Biol. 18, 174 (2017).
66. McCarthy, D. J., Campbell, K. R., Lun, A. T. L. \& Wills, Q. F. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics 33, 1179-1186 (2017).
67. Drew, C. J. G., Kyd, R. J. \& Morton, A. J. Complexin 1 knockout mice exhibit marked deficits in social behaviours but appear to be cognitively normal. Hum. Mol. Genet. 16, 2288-2305 (2007).
68. Huang, M. et al. SAVER: gene expression recovery for single-cell RNA sequencing. Nat. Methods 15, 539-542 (2018).
69. van Dijk, D. et al. Recovering gene interactions from single-cell data using data diffusion. Cell 174, 716-729 (2018).
70. Hu, Y. et al. WEDGE: imputation of gene expression values from single-cell RNA-seq datasets using biased matrix decomposition. Brief Bioinform 22, bbab085 (2021).
71. Dries, R. et al. Giotto: a toolbox for integrative analysis and visualization of spatial expression data. Genome Biol. 22, 78 (2021).
72. Jin, S. Q. et al. Inference and analysis of cell-cell communication using CellChat. Nat. Commun. 12, 1088 (2021).
73. Browaeys, R., Saelens, W. \& Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat. Methods 17, 159-162 (2020).
74. Noel, F. et al. Dissection of intercellular communication using the transcriptome-based framework ICELLNET. Nat. Commun. 12, 1089 (2021).
75. Cabello-Aguilar, S. et al. SingleCellSignalR: inference of intercellular networks from single-cell transcriptomics. Nucleic Acids Res. 48, e55 (2020).
76. Chen, A. et al. Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball patterned arrays. Preprint at bioRxiv https://doi.org/ 10.1101/2021.01.17.427004 (2021).

Publisher's note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.
© The Author(s), under exclusive licence to Springer Nature America, Inc. 2022

## Methods

Preprocessing of datasets. We preprocessed each dataset with the following steps. (1) Removal of low-quality cells. For scRNA-seq data, we used Seurat with parameters 'min.features $=200$ ' to remove cells for which fewer than 200 RNAs were captured. (2) Normalization of the expression matrix. For spatial transcriptomics datasets, we tested both the non-normalized and normalized expression matrices for input of each integration method. To normalize the expression matrices, we used the following equation:

$$
D_{i j}=\log \left(\bar{N} \times \frac{C_{i j}}{\sum_{j=1}^{M} C_{i j}}+1\right)
$$

where $C_{i j}$ represents the raw read count for gene $i$ in spot $j, D_{i j}$ represents the normalized read count for gene $i$ in spot $j$, and $\bar{N}$ is the median number of detected transcripts per cell. For scRNA-seq datasets, we normalized their expression matrix using the function 'NormalizeData' and default parameters in Seurat 3.2. (3) Selection of highly variable genes. For spatial transcriptomic datasets with more than 1,000 detected genes, we calculated the coefficient of variation of each gene using the following equation:

$$
C V_{i}=\frac{\sigma_{i}}{u_{i}}
$$

where $C V_{i}$ is the coefficient of variation of gene $i ; \sigma_{i}$ is the s.d. of the spatial distribution of gene $i$ in all spots; and $u_{i}$ is the average expression of gene $i$ in all spots. We identified 1,000 genes with the highest $C V_{i}$ values as highly variable genes and used their overlap with detected genes in the corresponding scRNA-seq data to construct the ground truth for each dataset. For spatial transcriptomic datasets with fewer than 1,000 detected genes, we used genes detected in both spatial transcriptomics and scRNA-seq data to build the ground truth of each dataset.

Parameter settings for integration methods. We evaluated the performance of eight integration methods, which can predict the spatial distribution of undetected transcripts, using tenfold crossvalidation. For a set of genes in the spatial data, we divided the genes into 10 portions, and iteratively used nine portions of the genes for integration (that is, the reference gene set used for training); the remaining one portion of genes was used for prediction. The parameters of each integration method were set as described below for each program.
gimVI. We followed the instructions on the gimVI website: https://docs.scvi-tools. org/en/0.8.0/user_guide/notebooks/gimvi_tutorial.html. The spatial distribution of genes was obtained using the model.get_imputed_values function with parameter normalized $=$ False.

SpaGE. We followed the guidelines on the GitHub repository of SpaGE: https:// github.com/tabdelaal/SpaGE/blob/master/SpaGE_Tutorial.ipynb. If the number of genes used for integration was greater than 50 (that is, $N_{\text {gene }}>50$ ), we set the parameter $\mathrm{n} \_\mathrm{pv}=N_{\text {gene }} / 2$.

Tangram. We followed the instructions on the Tangram GitHub repository: https:// github.com/broadinstitute/Tangram. we set the parameters as modes = 'clusters', density = 'rna_count_based'.

Seurat. We followed the instructions on the Seurat 3.2 website: https://satijalab.org/ seurat/archive/v3.2/integration.html. We set the parameter reduction = 'cca', k.filter $=$ NA. If $N_{\text {gene }}>30$, we set dims $=30$, otherwise we set dims $=N_{\text {gene }}$. The predicted spatial distribution of genes was obtained using the Seurat function "TransferData".

SpaOTsc. We followed the instructions on the SpaOTsc GitHub repository: https:// github.com/zcang/SpaOTsc. The spatial distribution of genes was obtained using the function 'issc.transport_plan' with parameters alpha $=0$, rho $=1.0$, epsilon $=$ 0.1 , scaling $=$ False.
novoSpaRc. We followed the guidelines on the GitHub repository of novoSpaRc: https://github.com/rajewsky-lab/novosparc/blob/master/reconstruct_drosophila_ embryo_tutorial.ipynb. We set the parameters as alpha_linea $\mathrm{r}=0.5$, loss_fun $=$ 'square_loss', epsilon $=5 \times 10^{-3}$. We trained novoSpaRc using the expression and spatial information of the training gene set, just as we did for the other methods, to ensure the fairness of the benchmarking study.

LIGER. We followed the instructions on the LIGER GitHub repository: https:// github.com/welch-lab/liger. The predicted spatial distribution of genes was obtained using the function 'imputeKNN' with parameters norm $=$ FALSE, scale $=$ FALSE. If $N_{\text {gene }}>30$, we set $\mathrm{knn} \_\mathrm{k}=30$, otherwise we set $\mathrm{knn} \_\mathrm{k}=N_{\text {gene }}$.
stPlus. We followed the instructions on the LIGER GitHub repository: http:// github.com/xy-chen16/stPlus. we set $\operatorname{tmin}=5$, neighbor $=50$.

We then evaluated the performance of 12 integration methods that can deconvolute the cell types of histological spots. The parameters of each integration method were set as described below for each program.

Cell2location. We followed the guidelines on the Cell2location website: https:// cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html. The single-cell regression model was trained with parameters max_epochs $=250$, $\mathrm{lr}=0.002$. The cell2location model was obtained with parameters max_epochs = 30,000.

RCTD. We followed the guidelines on the RCTD GitHub repository: https://raw. githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html. We set doublet_mode = 'full'.

DestVI. We followed the guidelines on the DestVI website: https://docs.scvi-tools. org/en/stable/tutorials/notebooks/DestVI_tutorial.html. The single-cell model was trained with parameters max_epochs $=250, \mathrm{lr}=0.001$, number of training genes $=$ 7,000 . The spatial model was trained with parameters max_epochs $=2,500$.

Tangram. We followed the instructions on the Tangram GitHub repository: https:// github.com/broadinstitute/Tangram. We set the parameters as modes = 'clusters', density = 'rna_count_based'. To deconvolute the cell types in space, we invoke ‘project_cell_annotation’ to transfer the annotation to space.

Seurat. We followed the instructions on the Seurat 3.2 website: https://satijalab. org/seurat/archive/v3.2/integration.html. We set the parameter $\operatorname{dim}=1: 30$, normalization.method = 'SCT'.

SpatialDWLS. We followed the guidelines on the SpatialDWLS website: https:// rubd.github.io/Giotto_site/articles/tut7_giotto_enrichment.html. We set the parameter as n_cell $=20$.

SPOTlight. We followed the guidelines on the SPOTlight GitHub repository: https://marcelosua.github.io/SPOTlight/. We set the parameter as transf = 'uv', method = 'nsNMF'.

Stereoscope. We followed the guidelines on the website: https://docs.scvi-tools. org/en/stable/user_guide/models/stereoscope.html. The single-cell model was trained with parameters max_epochs $=100$. The spatial model was trained with parameters max_epochs $=10,000$.

STRIDE. We followed the guidelines on the STRIDE website: https://stridespatial. readthedocs.io/en/latest/tutorials/Mouse_embryo.html. We set the parameter as ‘-normalize’.

DSTG. We followed the instructions on the DSTG GitHub repository: https:// github.com/Su-informatics-lab/DSTG.

SpaOTsc. We followed the instructions on the SpaOTsc GitHub repository: https:// github.com/zcang/SpaOTsc. The spatial distribution of genes was obtained using the function 'issc.transport_plan' with parameters alpha $=0$, rho $=1.0$, epsilon $=0.1$, scaling $=$ False.
novoSpaRc. We followed the guidelines on the GitHub repository of novoSpaRc: https://github.com/rajewsky-lab/novosparc/blob/master/reconstruct_drosophila_ embryo_tutorial.ipynb. We set the parameters as alpha_linear $=0.5$, loss_fun $=$ 'square_loss', epsilon $=5 \times 10^{-3}$. We trained novoSpaRc using the expression and spatial information of the training gene set, just as we did for the other methods, to ensure the fairness of the benchmarking study.

Benchmark metrics. We constructed a common pipeline to evaluate the performance of the integration methods for the 45 paired datasets. In the pipeline, we used the following five metrics to assess each integration method.

1. PCC. The PCC value was calculated using the following equation:

$$
P C C=\frac{E\left[\left(\tilde{\boldsymbol{x}}_{i}-\tilde{u}_{i}\right)\left(\boldsymbol{x}_{i}-u_{i}\right)\right]}{\tilde{\sigma}_{i} \sigma_{i}}
$$

where $\boldsymbol{x}_{i}$ and $\tilde{\boldsymbol{x}}_{i}$ are the spatial expression vectors of gene $i$ in the ground truth and the predicted result, respectively; $u_{i}$ and $\tilde{u}_{i}$ are the average expression value of gene $i$ in the ground truth and the predicted result, respectively; and $\sigma_{i}$ and $\tilde{\sigma}_{i}$ are the s.d. of the spatial expression of gene $i$ in the ground truth and the predicted result, respectively. For one gene, a higher PCC value indicates better prediction accuracy.
2. SSIM ${ }^{77}$. We first scaled the expression matrix as follows, so that the expression value of each gene was between 0 and 1 :

$$
x_{i j}^{\prime}=\frac{x_{i j}}{\max \left(\left\{x_{i 1}, \ldots, x_{i M}\right\}\right)}
$$

where $x_{i j}$ denotes the expression of gene $i$ in spot $j$, and $M$ is the total number of spots. Then we used the scaled gene expression and the following equation to calculate the SSIM value of each gene:

$$
\text { SSIM }=\frac{\left(2 \tilde{u}_{i} u_{i}+C_{1}^{2}\right)\left(2 \operatorname{cov}\left(\boldsymbol{x}_{i}^{\prime}, \tilde{\boldsymbol{x}}_{i}^{\prime}\right)+C_{2}^{2}\right)}{\left(\tilde{u}_{i}^{2}+u_{i}^{2}+C_{1}^{2}\right)\left(\tilde{\sigma}_{i}^{2}+\sigma_{i}^{2}+C_{2}^{2}\right)}
$$

where the definitions of $u_{i}, \tilde{u}_{i}, \sigma_{i}$, and $\tilde{\sigma}_{i}$ are similar to those for calculating the PCC value (but for scaled gene expression); $C_{1}$ and $C_{2}$ are 0.01 and 0.03 , respectively; and $\operatorname{cov}\left(\boldsymbol{x}_{i}, \tilde{\boldsymbol{x}}_{i}\right)$ is the covariance between the expression vector of gene $i$ in the ground truth (that is, $\boldsymbol{x}_{i}^{\boldsymbol{\prime}}$ ) and that of the predicted result (that is, $\tilde{\boldsymbol{x}}_{i}^{\boldsymbol{\prime}}$ ). For one gene, a higher SSIM value indicates better prediction accuracy.
3. RMSE. We first calculated the $z$-score of the spatial expression of each gene for all spots, then used the following equation to calculate RMSE:

$$
R M S E=\sqrt{\frac{1}{M} \sum_{j=1}^{M}\left(\tilde{z}_{i j}-z_{i j}\right)^{2}}
$$

where $z_{i j}$ and $\tilde{z}_{i j}$ are the $z$-score of the spatial expression of gene $i$ in spot $j$ in the ground truth and the predicted result, respectively. For one gene, a lower RMSE value indicates better prediction accuracy.
4. JS. JS uses relative information entropy (that is, Kullback-Leibler divergence) to determine the difference between two distributions. We first calculated the spatial distribution probability of each gene as follows:

$$
P_{i j}=\frac{x_{i j}}{\sum_{j=1}^{M} x_{i j}}
$$

where $x_{i j}$ denotes the expression of gene $i$ in spot $j, M$ is the total number of spots, and $P_{i j}$ is the distribution probability of gene $i$ in spot $j$. We then calculated the JS value of each gene using the following equations ${ }^{78}$ :

$$
\begin{gathered}
J S=\frac{1}{2} K L\left(\tilde{\boldsymbol{P}}_{i} \left\lvert\, \frac{\widetilde{\boldsymbol{P}}_{i}+\boldsymbol{P}_{i}}{2}\right.\right)+\frac{1}{2} K L\left(\boldsymbol{P}_{i} \left\lvert\, \frac{\widetilde{\boldsymbol{P}}_{i}+\boldsymbol{P}_{i}}{2}\right.\right) \\
K L\left(\boldsymbol{a}_{i}| | \boldsymbol{b}_{i}\right)=\sum_{j=0}^{M}\left(a_{i j} \times \log \frac{a_{i j}}{b_{i j}}\right)
\end{gathered}
$$

where $\boldsymbol{P}_{\boldsymbol{i}}$ and $\tilde{\boldsymbol{P}}_{i}$ are the spatial distribution probability vectors of gene $i$ in the ground truth and the predicted result, respectively, $K L\left(\boldsymbol{a}_{i} \| \boldsymbol{b}_{i}\right)$ is the Kullback-Leibler divergence between two probability distribution $\boldsymbol{a}_{i}$ and $\boldsymbol{b}_{i}$, and $a_{i j}$ and $b_{i j}$ are the predicted probability and real probability of gene $i$ in spot $j$, respectively. For one gene, a lower JS value indicates better prediction accuracy.
5. AS. We defined AS by aggregating PCC, SSIM, RMSE, and JS to evaluate the relative accuracy of the integration methods for each dataset. For one dataset, we calculated the average PCC, SSIM, RMSE, and JS of all genes predicted by each integration method. Then we sorted the PCC and SSIM values of the integration methods in ascending order to get RANK ${ }_{\text {PCC }}$ and RANK ${ }_{\text {SSIM }}$; the method with the highest PCC/SSIM value will have $\mathbf{R A N K}_{\text {PCC/SSIM }}=\mathbf{N}$, and the method with the lowest PCC/SSIM value will have $\mathbf{R A N K}_{\text {PCC/SSIM }}=1$. We also sorted the RMSE and JS values of the integration methods in descending order to get $\mathbf{R A N K}_{\text {RMSE }}$ and $\mathbf{R A N K}_{\text {JS }}$; the method with the highest RMSE/JS value will have RANK $_{\text {RMSE/JS }}=1$, and the method with the lowest RMSE/JS value will have RANK $_{\text {RMSE/JS }}=\mathbf{N}$. Finally, we calculated the average value of $\mathbf{R A N K}_{\text {PCC }}, \mathbf{R A N K}_{\text {SSIM }}, \mathbf{R A N K}_{\text {RMSE }}$, and $\mathbf{R A N K}_{\text {IS }}$ to obtain the AS value of each integration method, as follows:

$$
\mathbf{A S}=\frac{1}{4}\left(\mathbf{R A N K}_{\mathrm{PCC}}+\mathbf{R A N K}_{\mathrm{SSIM}}+\mathbf{R A N K}_{\mathrm{RMSE}}+\mathbf{R A N K}_{\mathrm{JS}}\right)
$$

For a dataset, the method with the highest AS value had the best performance among the integration methods.

Simulating 'multi-cell-spot problem' datasets. In order to obtain multi-cell-per-spot datasets with known cell compositions at each spot, we gridded single-cell resolution spatial transcriptomics datasets 4 and 10 to simulate datasets with potentially ambiguous cell type assignments per spot. For dataset 4 (seqFISH+; Smart-seq; mouse cortex), we defined a square with $500 \times 500$ pixels ( $\sim 51.5 \mu \mathrm{~m}$ ) as one spot-like region to grid the seqFISH+ slide, referring to the coarse-graining procedure introduced by SpatialDWLS. We summed the expression values of all cells in a grid to simulate a spot that may contain multiple cells, and took the center of the grid as the location of the spot. The simulated data of dataset 4 had 72 spots, and we calculated the percentage of cell types in each spot as the ground truth. For dataset 10 (STARmap; Smart-seq; mouse visual cortex), we used Seurat to cluster cells, and annotated the cell type of each cluster using marker genes ${ }^{79,80}$. We used marker genes Slc17a7 and Gad1 to annotate excitatory neurons and inhibitory neurons, respectively. The L2/3, L4, L5 and L6 excitatory neurons (eL2/3, eL4, eL5, eL6) were annotated, respectively, by marker genes Nov, Rorb, Sulf2, and Pcp4. Moreover, the VIP, SST, and PV inhibitor neurons were annotated, respectively, by marker genes Vip, Sst, and Pvalb. The microglia, astrocytes, oligodendrocytes, smooth-muscle, and endothelial cells were annotated by marker genes Pdgfra, Aqp4, Enpp2, Mgp, and Bsg, respectively. We then used a 750-pixel window to grid the STARmap slide. We summed the expression values of all cells in a grid to simulate a spot that may contain multiple cells, and took the center of the grid as the location of the spot. The simulated data of dataset 10 have 189 spots in total, and we calculated the percentage of cell types in each spot as the ground truth.

We also used PCC, SSIM, RMSE, and JS to assess the accuracy of Seurat, SpaOTsc, Tangram, and novoSpaRc in assigning cells to spatial locations in histological sections. We first counted the proportions of various types of cells in each spot. Then we introduced the cell type proportion of each spot into Eq. 3-Eq. 9 to calculate PCC, SSIM, RMSE, and JS values, which quantified the similarity (PCC/SSIM) or difference (RMSE/JS) between the predicted results and the ground truth. Finally, we used the two-sided Mann-Whitney U test to calculate the statistical significance of the difference in the prediction accuracy between different methods.

Simulating spatial datasets using scRNA-seq datasets. For the 32 simulated datasets, we devised the generation procedure of spatial transcriptomics data simulation by referring to the algorithm introduced by RCTD and Stereoscope. For each simulated spot, we first sampled cell numbers (Nc) in a uniform distribution in the range $5-15$, and sampled the number of cell types ( Nt ) in a uniform distribution in the range $2-6$. Then we assumed that these cell types have equal distribution possibility $\mathrm{P}=1 / \mathrm{Nt}$ in the spot, and randomly assigned cells from each cell types of the scRNA-seq data to the spot. To obtain the gene expression values at each spatial location, we summed the gene expression values of all cells in one spot. Referring to the method for constructing simulated datasets used in RCTD, we used Scuttle (http://bioconductor.org/packages/release/bioc/html/scuttle.html) to down-sample the number of counts per spot to $10 \%$ of the original value. We can obtain the percentage of a cell type at each spot by counting the number of cells corresponding to the cell type. For the large simulated dataset used for assessing the efficiency of each integration method, we used the same algorithm as above, but set the number of spots to 20000 , the number of cells to 10000 , and the number of cell types to 56.

Simulating datasets with high sparsity. We used 19 of the 45 paired datasets (that is, datasets $4,7,18,19,22,23,24,25,27,28,29,30,31,32,33,36,37,38,42$ ) to test the impact of the sparsity of expression matrices for each integration method. To simulate a dataset with high sparsity, we applied the Scuttle package (http://www. bioconductor.org/packages/release/bioc/html/scuttle.html) for the down-sampling of the spatial expression matrices of the datasets, and we also used the Splatter package to down-sample the datasets (down-sample rate $=0.2,0.4,0.6$, and 0.8 ). We down-sampled each dataset 10 times at different rates to avoid errors caused by random selection. To quantify the impact of expression matrix sparsity for each integration method, we counted the percentage of genes whose PCC values of the spatial distribution predicted from the original data and the down-sampled data were both greater than 0.5 , which was defined as the robustness score.

Computer platform. We ran CPU tests of the 16 integration methods on a computer cluster with four Intel Xeon E78860v4 CPUs ( $2.2 \mathrm{GHz}, 45 \mathrm{MB}$ L3 cache, 144 CPU cores in total) and 1 TB memory (DDR4 $2,400 \mathrm{MHz}$ ). The GPU tests for gimVI and Tangram were performed on a computer with Intel Xeon E5-2680v4 CPU ( $2.4 \mathrm{GHz}, 35 \mathrm{MB}$ L3 cache, 14 CPU cores in total), 128 GB memory, and NVIDIA Tesla K80 GPU (12 GB of memory, a total of 2496 CUDA cores).

To assess the impact of various data attributes (including the number of cells in scRNA-seq data, the number of spots in spatial data, and the number of genes used for training) on the computing resources consumed by those 8 integration methods capable of predicting the spatial distribution of undetected transcripts, we down-sampled the number of cells and the number of spots in dataset 40 and down-sampled the number of shared genes in dataset 6 . For the 10 integration methods that can perform cell type composition prediction of spots, we simulated a large dataset ( 10,000 spots, 20,000 cells) to evaluate the computer resources consumed by each method. Then, we down-sampled the number of cells, the number of spots, and the number of cell types in simulated datasets, and then evaluated the impacts of these data attributes on computing resources consumed.

Reporting Summary. Further information on research design is available in the Nature Research Reporting Summary linked to this article.

## Data availability

A summary of the individual accession numbers is given in Supplementary Table 1. The raw data are available from following study:
Dataset 1 (mouse gastrulation): seqFISH, https://content.cruk.cam.ac.uk/jmlab/ SpatialMouseAtlas2020/; 10X Chromium, 'Sample 21' in MouseGastrulationData within the R/Bioconductor data packageMouseGastrulationData.
Dataset 2 (mouse embryonic stem cell): seqFISH, https://zenodo.org/ record/3735329\#.YY69HZMza3J; Microwell-Seq, 'EmbryonicStemCells' in 'MCA_ BatchRemoved_Merge_dge.h5ad' file in https://figshare.com/articles/dataset/ MCA_DGE_Data/5435866.
Dataset 3 (mouse hippocampus): seqFISH, https://ars.els-cdn.com/content/imag e/1-s2.0-S0896627316307024-mmc6.xlsx; 10X Chromium, 'HIPP_sc_Rep1_10X sample' in GSE158450 in the GEO database. Dataset 4 (mouse cortex): seqFISH+, https://github.com/CaiGroup/ seqFISH-PLUS, and the spatial coordinate of each spot was generated using 'stitchFieldCoordinates' function in Giotto; Smart-seq, mouse primary visual
cortex (VISp) in the dataset in https://portal.brain-map.org/atlases-and-data/ rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 5 (mouse olfactory bulb): seqFISH+, https://github.com/CaiGroup/ seqFISH-PLUS; Drop-seq, GSE148360 in the GEO database.
Dataset 6 (mouse hypothalamic preoptic region): MERFISH, the eighteenth female parent mouse (animal ID = 18) in https://datadryad.org/stash/dataset/doi:10.5061/ dryad.8t8s248; 10X Chromium, GSE113576 in the GEO database.
Dataset 7 (human osteosarcoma): MERFISH, the 'B1_cell' used in https://www. pnas.org/doi/suppl/10.1073/pnas.1912459116/suppl_file/pnas.1912459116.sd12.
csv; 10X Chromium, BC22 in GSE152048 in the GEO database.
Dataset 8 (mouse primary motor cortex): MERFISH, 'mousel_slice162' in https://caltech.box.com/shared/static/dzqt6ryytmjbgyai356s1z0phtnsbaol.gz; 10X Chromium, https://data.nemoarchive.org/biccn/lab/zeng/transcriptome/ scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/.
Dataset 9 (mouse VISP): MERFISH, https://github.com/spacetx-spacejam/data/; Smart-seq, mouse primary visual cortex (VISp) in https://portal.brain-map.org/ atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 10 (mouse visual cortex): STARmap, '20180505_BY3_1kgenes' in https://www.starmapresources.com/data; Smart-seq, mouse primary visual cortex (VISp) in https://portal.brain-map.org/atlases-and-data/rnaseq/ mouse-vl-and-alm-smart-seq.
Dataset 11 (mouse prefronatal cortex): STARmap, '20180419_BZ9_control' in https://www.starmapresources.com/data; 10X Chromium, 'PFC_sc_Rep2_10X' in GSE158450 in the GEO database.
Dataset 12 (human middle temporal gyrus): ISS, https://github.com/ spacetx-spacejam/data; Smart-seq, https://portal.brain-map.org/atlases-and-data/ rnaseq/human-mtg-smart-seq.
Dataset 13 (mouse VISP): ISS, https://github.com/spacetx-spacejam/data; Smart-seq, mouse primary visual cortex (VISp) in the dataset in https://portal. brain-map.org/atlases-and-data/rnaseq/mouse-vl-and-alm-smart-seq.
Dataset 14 (Drosophila embryo): FISH, https://github.com/rajewsky-lab/distmap; Drop-seq, GSE95025 in the Gene Expression Omnibus (GEO) database. Dataset 15 (mouse somatosensory cortex): osmFISH, cortical regions in http://linnarssonlab.org/osmFISH/; Smart-seq, mouse somatosensory cortex (SSp) in https://portal.brain-map.org/atlases-and-data/rnaseq/ mouse-whole-cortex-and-hippocampus-smart-seq.
Dataset 16 (mouse VISP): BaristaSeq, https://github.com/spacetx-spacejam/data; Smart-seq, mouse primary visual cortex (VISp) in https://portal.brain-map.org/ atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 17 (mouse VISP): ExSeq, https://github.com/spacetx-spacejam/data; Smart-seq, mouse primary visual cortex (VISp) in https://portal.brain-map.org/ atlases-and-data/rnaseq/mouse-vl-and-alm-smart-seq.
Dataset 18 (mouse hindlimb muscle): 10X Visium, Vis5A in GSE161318 in the GEO database; 10X Chromium, D2_Ev3 in GSE159500 in the GEO database.
Dataset 19 (mouse hindlimb muscle): 10X Visium, Vis9A in GSE161318 in the GEO database; 10X Chromium, D7_Ev3 in GSE159500 in the GEO database.
Dataset 20 (human breast cancer): 10X Visium, 'CID3586' in https://zenodo.org/ record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 21 (human breast cancer): 10X Visium, '1160920F' in https://zenodo.org/ record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 22 (human breast cancer): 10X Visium, 'CID4290' in https://zenodo.org/ record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 23 (human breast cancer): 10X Visium, 'CID4465' in https://zenodo.org/ record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 24 (human breast cancer): 10X Visium, 'CID44971' https://zenodo.org/ record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 25 (human breast cancer): 10X Visium, 'CID4535' in https://zenodo.org/ record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 26 (zebrafish melanoma): 10X Visium, 'Visium-A' in GSE159709 in the GEO database; 10X Chromium, 'SingleCell-E' in GSE159709 in the GEO database.
Dataset 27 (mouse embryo): 10X Visium, 'Visium-A1' in GSE160137 in the GEO database; 10X Chromium, 'Pax2-GFP_SC-2' in GSE143806 in the GEO database. Dataset 28 (human prostate): 10X Visium, 'D25' in GSE159697 in the GEO database; 10X Chromium, 'V8' in GSE142489 in the GEO database.
Dataset 29 (mouse kidney): 10X Visium, Sham Model in GSE171406 in the GEO database; 10X Chromium, wild-type sham mouse in GSE171639 in the GEO database.
Dataset 30 (mouse kidney): 10X Visium, ischemia reperfusion injury model in GSE171406 in the GEO database; 10X Chromium, wild-type ischemic acute kidney injury mouse in GSE171639 in the GEO database.
Dataset 31 (mouse brain): 10X Visium, 'sectionl' in GSE153424 in the GEO database; 10X Chromium, 'brain1_cx' in GSE153424 in the GEO database.

Dataset 32 (mouse prefrontal cortex): 10X Visium, 'Visium_10X' in GSE158450 in the GEO database; 10X Chromium, 'PFC_sc_Rep1_10X' in GSE158450 in the GEO database.
Dataset 33 (mouse hippocampus): 10X Visium, 'Visium_10X' in GSE158450 in the GEO database; 10X Chromium, 'HIPP_sc_Rep1_10X' in GSE158450 in the GEO database.
Dataset 34 (mouse kidney): 10X Visium, GSE154107 in the GEO database; 10X Chromium, sample '(LPS36hr) scRNA-seq' in GSE151658 in the GEO database. Dataset 35 (human prostate): 10X Visium, 'ETOH' in GSE159697 in the GEO database; 10X Chromium, 'V8' in GSE142489 in the GEO database.
Dataset 36 (mouse lymph node): 10X Visium, 'PBS' samples of Tissue 1 in https:// github.com/romain-lopez/DestVI-reproducibility; 10X Chromium, 'PBS' samples in https://github.com/romain-lopez/DestVI-reproducibility.
Dataset 37 (mouse MCA205 tumor): 10X Visium, Tumor A1 of Tissue 1 in https:// github.com/romain-lopez/DestVI-reproducibility; 10X Chromium, https://github. com/romain-lopez/DestVI-reproducibility.
Dataset 38 (mouse primary motor cortex): 10X Visium, https://storage.googleapis. com/tommaso-brain-data/tangram_demo/Allen-Visium_Allen1_cell_count.h5ad; 10X Chromium, 'batch 9' in 'mop_sn_tutorial.h5ad' file from https://console.cloud. google.com/storage/browser/tommaso-brain-data.
Dataset 39 (mouse primary motor cortex): Slide-seq, https://storage.googleapis. com/tommaso-brain-data/tangram_demo/slideseq_MOp_1217.h5ad.gz; 10X Chromium, 'batch 9' in 'mop_sn_tutorial.h5ad' file from https://console.cloud. google.com/storage/browser/tommaso-brain-data.
Dataset 40 (mouse cerebellum): Slide-seqV2, SCP948 in https://singlecell. broadinstitute.org/single_cell/; 10X Chromium, sample M003 of study SCP795 in https://singlecell.broadinstitute.org/single_cell/.
Dataset 41 (mouse hippocampus): Slide-seqV2, 'Puck_200115_08' in https://singlecell.broadinstitute.org/single_cell/study/SCP815/ highly-sensitive-spatial-transcriptomics-at-near-cellular-resolution-with-slide-seq v2\#study-download; Drop-seq, we randomly sampled 10,000 cells from 'GSE116470_F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz' file in GSE116470 in the GEO database.
Dataset 42 (human squamous carcinoma): ST, GSM4284322 in the GEO database; 10X Chromium, 'GSE144236_cSCC_counts.txt.gz' in GSE144236 in the GEO database.
Dataset 43 (mouse hippocampus): ST, wild-type replicate 1 in https://data. mendeley.com/datasets/6s959w2zyr/1; 10X Chromium, GSE116470 in the GEO database.
Dataset 44 (mouse olfactory bulb): HDST, replicate1 in GSE130682 in the GEO database; 10X Chromium, WT1 samples used from GSE121891 in the GEO database.
Dataset 45 (mouse liver): Seq-scope, https://deepblue.lib.umich.edu/data/ downloads/gx41mj14n; Smart-seq2, liver sample in GSE109774 in the GEO database.
We also provide an open source website for users to download all the above datasets: https://drive.google.com/drive/folders/1pHmE9cg_tMcouV1LFJFtbyBJN p7oQo9J?usp=sharing.
Source data for figures and Extended Data Figures are provided with this paper. Source data are provided with this paper.

## Code availability

We uploaded the code and scripts used for the comparative analysis and figure plotting to GitHub: https://github.com/QuKunLab/SpatialBenchmarking. The package can also be used to analyze user's own datasets.

## References

77. Wang, Z., Bovik, A. C., Sheikh, H. R. \& Simoncelli, E. P. Image quality assessment: from error visibility to structural similarity. IEEE Trans. Image Process. 13, 600-612 (2004).
78. Lin, J. Divergence measures based on the Shannon entropy. IEEE Trans. Inf. Theory 37, 145-151 (1991).
79. Tasic, B. et al. Adult mouse cortical cell taxonomy revealed by single cell transcriptomics. Nat. Neurosci. 19, 335-346 (2016).
80. Zeisel, A. et al. Brain structure. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science 347, 1138-1142 (2015).

## Acknowledgements

This work was supported by the National Key R\&D Program of China (2020YFA0112200 to K.Q.), the National Natural Science Foundation of China grants (T2125012, 91940306, 31970858, and 31771428 to K.Q.; 32170668 to B.L.; 81871479 to J.L.), CAS Project for Young Scientists in Basic Research YSBR-005 (to K.Q.), Anhui Province Science and Technology Key Program (202003a07020021 to K.Q.) and the Fundamental Research Funds for the Central Universities (YD2070002019, WK9110000141, and WK2070000158 to K.Q.; WK9100000001 to J.L). We thank the USTC supercomputing center and the School of Life Science Bioinformatics Center for providing computing resources for this project.

## Author contributions

K.Q. and B.L. conceived the project. B.L., W.Z. and C.G. designed the framework and performed data analysis with help from H.X., L.L., M.F, Y.H., X.Y., X.Z., F.C. and T.X., X.Z., M.T., K.L., J.L. and L.C. contributed to revision of the manuscript. B.L., K.Q., C.G., and W.Z. wrote the manuscript with inputs from all authors. K.Q. supervised the entire project. All authors read and approved the final manuscript.

## Competing interests

The authors declare no competing interests.

## Additional information

Extended data is available for this paper at https://doi.org/10.1038/s41592-022-01480-9.
Supplementary information The online version contains supplementary material available at https://doi.org/10.1038/s41592-022-01480-9.
Correspondence and requests for materials should be addressed to Kun Qu.
Peer review information Nature Methods thanks Ahmed Mahfouz and the other, anonymous, reviewer for their contribution to the peer review of this work. Primary Handling editor: Lin Tang, in collaboration with the Nature Methods team. Peer reviewer reports are available.
Reprints and permissions information is available at www.nature.com/reprints.

Data 42 (ST; 10X Chromium; human squamous carcinoma)

a
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-14.jpg?height=334&width=1553&top_left_y=246&top_left_x=274)

b

Data 42 (ST; 10X Chromium; human squamous carcinoma)
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-14.jpg?height=1003&width=1263&top_left_y=661&top_left_x=229)

Data 42 (ST; 10X Chromium; human squamous carcinoma)
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-14.jpg?height=511&width=657&top_left_y=1736&top_left_x=238)

Extended Data Fig. 1 | Comparing the accuracy of eight integration methods in predicting the spatial distribution of RNA transcripts. a, The spatial distribution of COL17A1 in dataset 42 (ST; 10X Chromium; human squamous carcinoma), including the ground truth and prediction results from the integration methods. PCC: Pearson Correlation Coefficient between the expression vector of a transcript in the ground truth and that of the predicted result. b, Bar plots of PCC, SSIM, RMSE, and JS of each integration method in predicting the spatial distribution of transcripts of dataset 42 . SSIM: Structural Similarity Index; RMSE: Root Mean Square Error; JS: Jensen-Shannon divergence. Data are presented as mean values $\pm 95 \%$ confidence intervals; $\mathrm{n}=948$ predicted genes. $\mathbf{c}$, The violin plot of AS (accuracy score, aggregated from PCC, SSIM, RMSE, and JS; see Methods) of the eight integration methods for transcripts in dataset 42 . Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range; $n=4$ benchmark metrics.

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-15.jpg?height=2074&width=1713&top_left_y=188&top_left_x=180)
Extended Data Fig. 2 | The boxplots of PCC, SSIM, RMSE, JS values of each integration method in predicting the spatial distribution of RNA transcripts of $\mathbf{4 5}$ paired spatial transcriptomics and scRNA-seq datasets. The boxplots of PCC, SSIM, RMSE, JS values of each integration method in predicting the spatial distribution of RNA transcripts of 45 paired spatial transcriptomics and scRNA-seq datasets. Center line, median; box limits, upper and lower quartiles; whiskers, $0.5 \times$ interquartile range, the number of genes for each dataset is shown at the top of each panel.

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-16.jpg?height=1505&width=1737&top_left_y=167&top_left_x=165)
Extended Data Fig. 3 | PCC, SSIM, RMSE, JS and AS of spatial distribution of RNA transcripts predicted by each integration method for the 45 paired spatial transcriptomics and scRNA-seq datasets. a-g, Boxplots of AS (accuracy score, aggregated from PCC, SSIM, RMSE, and JS; see Methods) of the integration methods for transcripts in the 17 image-based datasets (a), 28 seq-based datasets (b), 32 simulated datasets (c), 2110 X visium datasets (d), 5 seqFISH datasets (e), 4 MERFISH datasets (f), 3 Slide-seq datasets (g). Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range. Grey dots indicate the prediction result is not available, as the tool made an error when predictions.

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-17.jpg?height=2147&width=1394&top_left_y=184&top_left_x=336)
Extended Data Fig. 4 | The PCC values of each integration method when processing the raw expression matrices and the normalized expression matrices. The PCC values of each integration method when processing the raw expression matrices and the normalized expression matrices. R-R: raw expression matrix of spatial data and raw expression matrix of scRNA-seq data; N-R: normalized expression matrix of spatial data and raw expression matrix of scRNA-seq data; R-N: raw expression matrix of spatial data and normalized expression matrix of scRNA-seq data; N-N: normalized expression matrix of spatial data and normalized expression matrix of scRNA-seq data; $\mathrm{n}=43$ independent datasets. Dataset6 and Dataset8 are excluded, as the normalized expression matrix of spatial data has been normalized.

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-18.jpg?height=1170&width=1761&top_left_y=180&top_left_x=156)
Extended Data Fig. 5 | Impact of normalization on the accuracy of eight integration methods that can predict the spatial distribution of RNA transcripts. $\mathbf{a} \boldsymbol{,} \mathbf{b}$, Boxplots of the PCC values of the eight integration methods for 28 seq-based datasets (a) or 15 image-based datasets (b) when using the four schemes of input expression matrices (that is R-R, R-N, N-R, and N-N, see their definition in the legend of Extended Data Fig. 4). For the genes predicted by each method, we removed outliers using 10\%-90\% confidence interval. Statistical significance was analyzed with two-sided paired t -test, * $\mathrm{P}<0.05$, ${ }^{* *} \mathrm{P}<0.01,{ }^{* * *} \mathrm{P}<0.001$ and ${ }^{* * * *} \mathrm{P}<0.0001$. Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range. c-f, Boxplots of the AS values of the eight integration methods for all the 45 paired datasets when using the four schemes of input expression matrices. For the genes predicted by each method, we removed outliers using $10 \%-90 \%$ confidence interval. Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range; $\mathrm{n}=43$ independent datasets.

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-19.jpg?height=1430&width=793&top_left_y=182&top_left_x=613)

Extended Data Fig. 6 | Correlation between the four metrics (PCC, SSIM, RMSE, and JS) and the sparsity of each examined spatial expression matrix. Correlation between the four metrics (PCC, SSIM, RMSE, and JS) and the sparsity of each examined spatial expression matrix. For all the eight integration methods that can predict the spatial distribution of transcripts, the JS values are linearly positively correlated with the sparsity of expression matrices of the spatial transcriptomics data ( $R 2 \geq 0.50$ ).

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-20.jpg?height=2123&width=1718&top_left_y=182&top_left_x=175)
Extended Data Fig. 7 | See next page for caption.

Extended Data Fig. 7 | Comparing the accuracy of the eight integration methods for sparse expression matrices down-sampled from the original datasets using Scuttle. a, Spatial distribution of Cplx1 expression in dataset 4 (seqFISH+; Smart-seq; mouse cortex), predicted from the original data and down-sampled data (down-sampling rate $=0.8$ ). $\mathbf{b}$, PCC of the spatial distribution of transcripts predicted from the original data and down-sampled data from dataset 4. The PCC values of the red-colored transcripts are greater than 0.5 for both the original and the down-sampled data. The proportion of the red-colored transcripts in all transcripts was defined as the 'robustness score' (RS). $\mathbf{c}$, RS values of the eight integration methods when processing sparse expression matrices down-sampled from dataset 4 at different down-sampling rates. d, RS values of the eight integration methods when processing the sparse expression matrices of the down-sampled datasets. The original datasets (used to generate the down-sampled datasets) capture $>1000$ genes from $>100$ spots, and the sparsity of the expression matrices is $<0.7$. Data are presented as mean values $\pm 95 \%$ confidence intervals; $n=19$ independent datasets.

Dataset 10 (STARmap; Smart-seq; mouse visual cortex)

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-22.jpg?height=2042&width=1761&top_left_y=253&top_left_x=156)
Extended Data Fig. 8 | See next page for caption.

Extended Data Fig. 8 | Comparing the performance of the twelve integration methods in cell type deconvolution. a, PCC, SSIM, RMSE, and JS values for the cell type composition of the spots simulated from dataset 10, generated by twelve integration methods. Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range; $\mathrm{n}=12$ predicted cell types. $\mathbf{b}, \mathrm{A}$ seqFISH + slide of dataset 4 (seqFISH + ; Smart-seq; mouse cortex) with cells annotated by cell type. Each grid represents a simulated spot containing $1 \sim 18$ cells. $\mathbf{c}$, The proportion of L5\&6 excitatory neurons in the spots simulated from dataset 4, including the ground truth and the predicted results of twelve integration methods. d, PCC, SSIM, RMSE, and JS values for the cell type composition of the spots simulated from dataset 4, generated by twelve integration methods. Center line, median; box limits, upper and lower quartiles; whiskers, $1.5 \times$ interquartile range; $\mathrm{n}=8$ predicted cell types. $\mathbf{e}$, PCC, SSIM, RMSE, and JS values for the cell type composition of the spots in all the simulated datasets ( $\mathrm{n}=32$ ), generated by ten integration methods. SpaOTsc and novoSpaRc are excluded, as they require spatial location information for each spot, which is not available in the simulated datasets. Data are presented as mean values $\pm 95 \%$ confidence intervals; $n=32$ independent datasets.

![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-24.jpg?height=1684&width=1722&top_left_y=173&top_left_x=171)
Extended Data Fig. 9 | Computer resources consumed by each integration method. a-c, The impact of the number of cells in scRNA-seq data (a), the number of spots in spatial data (b), and the number of genes used for training (c), on computational resources consumed by the integration methods that can predict the spatial distribution of undetected transcripts. d, The computer time and memory spent by the integration methods that can deconvolute cell types of histological spots, when processing a simulated dataset which contains 20000 spots in its spatial transcriptomics data and 10000 cells in its scRNA-seq data. $\mathbf{e}-\mathbf{g}$, The impacts of the number of cells in scRNA-seq data (e), the number of spots in spatial data (f), and the number of the cell types (g) on computational resources consumed by the integration methods that can deconvolute cell types of histological spots.

## nature research

Corresponding author(s):

Last updated by author(s): Mar 12, 2022

## Reporting Summary

Nature Research wishes to improve the reproducibility of the work that we publish. This form provides structure for consistency and transparency in reporting. For further information on Nature Research policies, see our Editorial Policies and the Editorial Policy Checklist.

## Statistics

For all statistical analyses, confirm that the following items are present in the figure legend, table legend, main text, or Methods section.

## n/a

□
□
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=48&width=61&top_left_y=936&top_left_x=133)
□
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=41&width=57&top_left_y=1014&top_left_x=137)
□
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=49&width=74&top_left_y=1082&top_left_x=120)
□
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=53&width=74&top_left_y=1140&top_left_x=120)
□
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=52&width=61&top_left_y=1220&top_left_x=133)
□
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=55&width=57&top_left_y=1314&top_left_x=137)

For null hypothesis testing, the test statistic (e.g. $F, t, r$ ) with confidence intervals, effect sizes, degrees of freedom and $P$ value noted Give $P$ values as exact values whenever suitable.
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=51&width=57&top_left_y=1400&top_left_x=85)
□ For Bayesian analysis, information on the choice of priors and Markov chain Monte Carlo settings
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=52&width=54&top_left_y=1452&top_left_x=88)
□
Confirmed
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=61&width=57&top_left_y=863&top_left_x=137)

The exact sample size $(n)$ for each experimental group/condition, given as a discrete number and unit of measurement
A statement on whether measurements were taken from distinct samples or whether the same sample was measured repeatedly
The statistical test(s) used AND whether they are one- or two-sided
Only common tests should be described solely by name; describe more complex techniques in the Methods section.
A description of all covariates tested
A description of any assumptions or corrections, such as tests of normality and adjustment for multiple comparisons
A full description of the statistical parameters including central tendency (e.g. means) or other basic estimates (e.g. regression coefficient) AND variation (e.g. standard deviation) or associated estimates of uncertainty (e.g. confidence intervals)
□ For hierarchical and complex designs, identification of the appropriate level for tests and full reporting of outcomes
![](https://cdn.mathpix.com/cropped/2025_11_30_443943e6a7be4b3effb4g-25.jpg?height=59&width=55&top_left_y=1516&top_left_x=141)

Estimates of effect sizes (e.g. Cohen's $d$, Pearson's $r$ ), indicating how they were calculated

Our web collection on statistics for biologists contains articles on many of the points above.

## Software and code

## Policy information about availability of computer code

Data collection
Data analysis

No software was used for data collection.

We compared the performance of 8 integration methods for predicting the spatial distribution of undetected transcripts: gimVI (Version 0.8.0b0), SpaGE(no version), Tangram(Version 1.0.0), Seurat (Version 3.6.3), SpaOTsc (Version 0.2), LIGER (Version 0.5.0), novoSpaRc (Version 0.4 .3 ), stPlus(Version 0.0 .6 ).

We compared the performance of 12 integration methods for predicting the cell type composition of spots: Cell2location(Version 0.6a0), DestVI((Version 0.14.4), SPOTlight((Version 0.1.7), SpatiaIDWLS(Version 1.0.4), Seurat(Version 4.0.5), Tangram((Version 1.0.0), RCTD(Version 1.2.0), Stereoscope(Version 0.14.4), STRIDE(Version 0.0.1b0), DSTG (Version 0.0.1), novoSpaRc (Version 0.4.3), SpaOTsc (Version 0.2). We applied the Scuttle(Version 1.0.4) and Splatte(Version 1.19.3) package for the down-sampling of the spatial expression matrices of the datasets.
The code used in this paper is available at https://github.com/QuKunLab/SpatialBenchmarking.
The file provided at https://github.com/QuKunLab/SpatialBenchmarking/Benchmarkingenvironment.yml lists the software dependencies with version numbers.

For manuscripts utilizing custom algorithms or software that are central to the research but not yet described in published literature, software must be made available to editors and reviewers. We strongly encourage code deposition in a community repository (e.g. GitHub). See the Nature Research guidelines for submitting code \& software for further information.

## Data

## Policy information about availability of data

All manuscripts must include a data availability statement. This statement should provide the following information, where applicable:

- Accession codes, unique identifiers, or web links for publicly available datasets
- A list of figures that have associated raw data
- A description of any restrictions on data availability

Dataset 1 : seqFISH, https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/; 10X Chromium, 'Sample 21' in MouseGastrulationData within the R/ Bioconductor data packageMouseGastrulationData.
Dataset 2 : seqFISH, https://zenodo.org/record/3735329\#.YY69HZMza3J; Microwell-Seq, ‘EmbryonicStemCells’ in ‘MCA_BatchRemoved_Merge_dge.h5ad’ file in https://figshare.com/articles/dataset/MCA_DGE_Data/5435866.
Dataset 3 : seqFISH, https://ars.els-cdn.com/content/image/1-s2.0-S0896627316307024-mmc6.xlsx; 10X Chromium, "HIPP_sc_Rep1_10X sample" in GSE158450 in the GEO database.
Dataset 4 : seqFISH+, https://github.com/CaiGroup/seqFISH-PLUS, and the spatial coordinate of each spot was generated using "stitchFieldCoordinates" function in Giotto; Smart-seq, mouse primary visual cortex (VISp) in the dataset in https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 5 : seqFISH+, https://github.com/CaiGroup/seqFISH-PLUS; Drop-seq, GSE148360 in the GEO database.
Dataset 6 : MERFISH, the 18th female parent mouse (Animal ID = 18) in https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248; 10X Chromium, GSE113576 in the GEO database.
Dataset 7 : MERFISH, the 'B1_cell' used in https://www.pnas.org/doi/suppl/10.1073/pnas.1912459116/suppl_file/pnas.1912459116.sd12.csv; 10X Chromium, BC22 in GSE152048 in the GEO database.
Dataset 8 : MERFISH, 'mouse1_slice162' in https://caltech.box.com/shared/static/dzqt6ryytmjbgyai356s1zOphtnsbaol.gz; 10X Chromium, https:// data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/.
Dataset 9 : MERFISH, https://github.com/spacetx-spacejam/data/; Smart-seq, mouse primary visual cortex (VISp) in https://portal.brain-map.org/atlases-and-data/ rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 10 : STARmap, "20180505_BY3_1kgenes" in https://www.starmapresources.com/data; Smart-seq, mouse primary visual cortex (VISp) in https:// portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 11 : STARmap, "20180419_BZ9_control" in https://www.starmapresources.com/data; 10X Chromium, "PFC_sc_Rep2_10X" in GSE158450 in the GEO database.
Dataset 12 : ISS, https://github.com/spacetx-spacejam/data; Smart-seq, https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq.
Dataset 13 : ISS, https://github.com/spacetx-spacejam/data; Smart-seq, mouse primary visual cortex (VISp) in the dataset in https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 14 : FISH, https://github.com/rajewsky-lab/distmap; Drop-seq, GSE95025 in the Gene Expression Omnibus (GEO) database.
Dataset 15 : osmFISH, cortical regions in http://linnarssonlab.org/osmFISH/; Smart-seq, mouse somatosensory cortex (SSp) in https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq.
Dataset 16 : BARISTASeq, https://github.com/spacetx-spacejam/data; Smart-seq, mouse primary visual cortex (VISp) in https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 17 : ExSeq, https://github.com/spacetx-spacejam/data; Smart-seq, mouse primary visual cortex (VISp) in https://portal.brain-map.org/atlases-and-data/ rnaseq/mouse-v1-and-alm-smart-seq.
Dataset 18 : 10X Visium, Vis5A in GSE161318 in the GEO database; 10X Chromium, D2_Ev3 in GSE159500 in the GEO database.
Dataset 19 : 10X Visium, Vis9A in GSE161318 in the GEO database; 10X Chromium, D7_Ev3 in GSE159500 in the GEO database.
Dataset 20: 10X Visium, "CID3586" in https://zenodo.org/record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 21 : 10X Visium, "1160920F" in https://zenodo.org/record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 22 : 10X Visium, "CID4290" in https://zenodo.org/record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 23 : 10X Visium, "CID4465" in https://zenodo.org/record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 24 : 10X Visium, "CID44971" https://zenodo.org/record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 25 : 10X Visium, "CID4535" in https://zenodo.org/record/4739739\#.YY6N_pMzaWC; 10X Chromium, 'CID3586' in GSE176078 in the GEO database.
Dataset 26 : 10X Visium, "Visium-A" in GSE159709 in the GEO database; 10X Chromium, "SingleCell-E" in GSE159709 in the GEO database.
Dataset 27: 10X Visium, "Visium-A1" in GSE160137 in the GEO database; 10X Chromium, "Pax2-GFP_SC-2" in GSE143806 in the GEO database.
Dataset 28 : 10X Visium, "D25" in GSE159697 in the GEO database; 10X Chromium, "V8" in GSE142489 in the GEO database.
Dataset 29 : 10X Visium, Sham Model in GSE171406 in the GEO database; 10X Chromium, wild-type sham mouse in GSE171639 in the GEO database.
Dataset 30 : 10X Visium, Ischemia Reperfusion Injury Model in GSE171406 in the GEO database; 10X Chromium, wild-type ischemic acute kidney injury mouse in GSE171639 in the GEO database.
Dataset 31 : 10X Visium, "section1" in GSE153424 in the GEO database; 10X Chromium, "brain1_cx" in GSE153424 in the GEO database.
Dataset 32 : 10X Visium, "Visium_10X" in GSE158450 in the GEO database; 10X Chromium, "PFC_sc_Rep1_10X" in GSE158450 in the GEO database.
Dataset 33 : 10X Visium, "Visium_10X" in GSE158450 in the GEO database; 10X Chromium, "HIPP_sc_Rep1_10X" in GSE158450 in the GEO database.
Dataset 34 : 10X Visium, GSE154107 in the GEO database; 10X Chromium, sample "(LPS36hr) scRNAseq" in GSE151658 in the GEO database.
Dataset 35 : 10X Visium, "ETOH" in GSE159697 in the GEO database; 10X Chromium, "V8" in GSE142489 in the GEO database.
Dataset 36 : 10X Visium, 'PBS' samples of Tissue 1 in https://github.com/romain-lopez/DestVI-reproducibility; 10X Chromium, 'PBS' samples in https://github.com/ romain-lopez/DestVI-reproducibility.
Dataset 37 : 10X Visium, Tumor A1 of Tissue 1 in https://github.com/romain-lopez/DestVI-reproducibility; 10X Chromium, https://github.com/romain-lopez/DestVIreproducibility.
Dataset 38 : 10X Visium, https://console.cloud.google.com/storage/browser/tommaso-brain-data; 10X Chromium, 'batch 9' in
'tangram_demo_mop_sn_tutorial.h5ad' file from https://console.cloud.google.com/storage/browser/tommaso-brain-data.
Dataset 39 : Slide-seq, https://console.cloud.google.com/storage/browser/tommaso-brain-data; 10X Chromium, 'batch 9' in
'tangram_demo_mop_sn_tutorial.h5ad' file from https://console.cloud.google.com/storage/browser/tommaso-brain-data.
Dataset 40 : Slide-seqV2, SCP948 in https://singlecell.broadinstitute.org/single_cell/; 10X Chromium, sample M003 of study SCP795 in https:// singlecell.broadinstitute.org/single_cell/.
Dataset 41 : Slide-seqV2, "Puck_200115_08" in https://singlecell.broadinstitute.org/single_cell/study/SCP815/highly-sensitive-spatial-transcriptomics-at-near-cellular-resolution-with-slide-seqv2\#study-download; Drop-seq, we randomly sampled 10,000 cells from
'GSE116470_F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz' file in GSE116470 in the GEO database.
Dataset 42 : ST, GSM4284322 in the GEO database; 10X Chromium, "GSE144236_cSCC" in GSE144236 in the GEO database.
Dataset 43 : ST, wild-type replicate 1 in https://data.mendeley.com/datasets/6s959w2zyr/1; 10X Chromium, GSE116470 in the GEO database.
Dataset 44 : HDST, replicate1 in GSE130682 in the GEO database; 10X Chromium, WT1 samples used from GSE121891 in the GEO database.

The raw and processed data of these datasets are also available at : https://drive.google.com/drive/folders/1pHmE9cg_tMcouV1LFJFtbyBJNp7oQo9J? usp=sharing. A summary of these datasets is given in Supplementary Table 1.

## Field-specific reporting

Please select the one below that is the best fit for your research. If you are not sure, read the appropriate sections before making your selection.
Life sciences □ Behavioural \& social sciences □ Ecological, evolutionary \& environmental sciences

For a reference copy of the document with all sections, see nature.com/documents/nr-reporting-summary-flat.pdf

## Life sciences study design

All studies must disclose on these points even when the disclosure is negative.
Sample size
We used 45 paired spatial transcriptomics and scRNA-seq datasets from published studies. The spatial transcriptomic datasets were produced by 13 different spatial transcriptomics approaches. A summary of the individual accession numbers is given in Supplementary Table 1. The details of these datasets are listed as follows:
Dataset 1 has 8425 spots in spatial transcriptomics data and 4651 cells in scRNA-seq data;
Dataset 2 has 175 spots in spatial transcriptomics data and 9991 cells in scRNA-seq data;
Dataset 3 has 3585 spots in spatial transcriptomics data and 8596 cells in scRNA-seq data;
Dataset 4 has 524 spots in spatial transcriptomics data and 14249 cells in scRNA-seq data;
Dataset 5 has 2050 spots in spatial transcriptomics data and 31217 cells in scRNA-seq data;
Dataset 6 has 4975 spots in spatial transcriptomics data and 31299 cells in scRNA-seq data;
Dataset 7 has 645 spots in spatial transcriptomics data and 9234 cells in scRNA-seq data;
Dataset 8 has 6963 spots in spatial transcriptomics data and 7240 cells in scRNA-seq data;
Dataset 9 has 2399 spots in spatial transcriptomics data and 14249 cells in scRNA-seq data;
Dataset 10 has 1549 spots in spatial transcriptomics data and 14249 cells in scRNA-seq data;
Dataset 11 has 1380 spots in spatial transcriptomics data and 7737 cells in scRNA-seq data;
Dataset 12 has 6000 spots in spatial transcriptomics data and 15928 cells in scRNA-seq data;
Dataset 13 has 6000 spots in spatial transcriptomics data and 14249 cells in scRNA-seq data;
Dataset 14 has 3039 spots in spatial transcriptomics data and 1297 cells in scRNA-seq data;
Dataset 15 has 3405 spots in spatial transcriptomics data and 5613 cells in scRNA-seq data;
Dataset 16 has 11426 spots in spatial transcriptomics data and 14249 cells in scRNA-seq data;
Dataset 17 has 1154 spots in spatial transcriptomics data and 14249 cells in scRNA-seq data;
Dataset 18 has 982 spots in spatial transcriptomics data and 4748 cells in scRNA-seq data;
Dataset 19 has 995 spots in spatial transcriptomics data and 4816 cells in scRNA-seq data;
Dataset 20 has 4784 spots in spatial transcriptomics data and 6178 cells in scRNA-seq data;
Dataset 21 has 4895 spots in spatial transcriptomics data and 6178 cells in scRNA-seq data;
Dataset 22 has 2432 spots in spatial transcriptomics data and 6178 cells in scRNA-seq data;
Dataset 23 has 1211 spots in spatial transcriptomics data and 6178 cells in scRNA-seq data;
Dataset 24 has 1162 spots in spatial transcriptomics data and 6178 cells in scRNA-seq data;
Dataset 25 has 1127 spots in spatial transcriptomics data and 6178 cells in scRNA-seq data;
Dataset 26 has 2425 spots in spatial transcriptomics data and 1911 cells in scRNA-seq data;
Dataset 27 has 198 spots in spatial transcriptomics data and 3415 cells in scRNA-seq data;
Dataset 28 has 277 spots in spatial transcriptomics data and 4740 cells in scRNA-seq data;
Dataset 29 has 1835 spots in spatial transcriptomics data and 10872 cells in scRNA-seq data;
Dataset 30 has 2064 spots in spatial transcriptomics data and 13600 cells in scRNA-seq data;
Dataset 31 has 3805 spots in spatial transcriptomics data and 8798 cells in scRNA-seq data;
Dataset 32 has 3024 spots in spatial transcriptomics data and 3512 cells in scRNA-seq data;
Dataset 33 has 3024 spots in spatial transcriptomics data and 8653 cells in scRNA-seq data;
Dataset 34 has 1888 spots in spatial transcriptomics data and 8346 cells in scRNA-seq data;
Dataset 35 has 744 spots in spatial transcriptomics data and 4561 cells in scRNA-seq data;
Dataset 36 has 369 spots in spatial transcriptomics data and 7268 cells in scRNA-seq data;
Dataset 37 has 2125 spots in spatial transcriptomics data and 7185 cells in scRNA-seq data;
Dataset 38 has 2669 spots in spatial transcriptomics data and 3499 cells in scRNA-seq data;
Dataset 39 has 9852 spots in spatial transcriptomics data and 3499 cells in scRNA-seq data;
Dataset 40 has 41674 spots in spatial transcriptomics data and 26252 cells in scRNA-seq data;
Dataset 41 has 6000 spots in spatial transcriptomics data and 10000 cells in scRNA-seq data;
Dataset 42 has 1145 spots in spatial transcriptomics data and 48164 cells in scRNA-seq data;
Dataset 43 has 604 spots in spatial transcriptomics data and 15095 cells in scRNA-seq data;
Dataset 44 has 6000 spots in spatial transcriptomics data and 10259 cells in scRNA-seq data;
Dataset 45 has 2177 spots in spatial transcriptomics data and 981 cells in scRNA-seq data;
nature research | reporting summary April2020

Data exclusions

Replication
No data was excluded from the study. In the data preprocessing step of scRNA-seq data, we filtered cells that captured >200 RNA transcripts.

Randomization This is not relevant to our study because we reanalyzed publicly available data or generated our own synthetic data.

Blinding
In the comparison, all the methods were blinded to the ground truth of the spatial transcriptomics data. The outputs from the methods were then compared to the ground truth available in the respective datasets.

## Reporting for specific materials, systems and methods

We require information from authors about some types of materials, experimental systems and methods used in many studies. Here, indicate whether each material, system or method listed is relevant to your study. If you are not sure if a list item applies to your research, read the appropriate section before selecting a response.

Materials \& experimental systems
| n/a | Involved in the study | n/a |
| :--- | :--- | :--- |
| X | □ Antibodies | X |
| X | □ Eukaryotic cell lines | X |
| X | □ Palaeontology and archaeology | X |
| X | □ Animals and other organisms |  |
| X | □ Human research participants |  |
| X | □ Clinical data |  |
| X | □ Dual use research of concern |  |


nature research | reporting summary April2020


[^0]:    ${ }^{1}$ Department of Oncology, The First Affiliated Hospital of USTC, School of Basic Medical Sciences, Division of Life Sciences and Medicine, University of Science and Technology of China, Hefei, China. ${ }^{2}$ Institute of Artificial Intelligence, Hefei Comprehensive National Science Center, Hefei, China. ${ }^{3}$ Division of Life Sciences and Medicine, University of Science and Technology of China, Hefei, China. ${ }^{4}$ School of Mathematical Sciences, University of Science and Technology of China, Hefei, China. ${ }^{5}$ CAS Key Laboratory of Microbial Physiological and Metabolic Engineering, Institute of Microbiology, Chinese Academy of Sciences, Beijing, China. ${ }^{6}$ CAS Center for Excellence in Molecular Cell Sciences, the CAS Key Laboratory of Innate Immunity and Chronic Disease, University of Science and Technology of China, Hefei, China. ${ }^{7}$ These authors contributed equally: Bin Li, Wen Zhang, Chuang Guo. ®e-mail: qukun@ustc.edu.cn


---

# Paper 4: SPOTlight

*For correspondence:
yvan.saeys@ugent.be
Competing interest: The authors declare that no competing interests exist.

Funding: See page 17
Preprint posted
24 March 2023
Sent for Review
18 April 2023
Reviewed preprint posted
31 July 2023
Reviewed preprint revised
08 April 2024
Version of Record published 24 May 2024

Reviewing Editor: Luca Pinello, Massachusetts General Hospital, United States
(c) Copyright Sang-aram et al. This article is distributed under the terms of the Creative Commons Attribution License, which permits unrestricted use and redistribution provided that the original author and source are credited.

# Spotless, a reproducible pipeline for benchmarking cell type deconvolution in spatial transcriptomics 

Chananchida Sang-aram ${ }^{1,2}$, Robin Browaeys ${ }^{1,2}$, Ruth Seurinck ${ }^{1,2}$, Yvan Saeys ${ }^{1,2 *}$<br>¹Data Mining and Modelling for Biomedicine, VIB Center for Inflammation Research, Ghent, Belgium; ${ }^{2}$ Department of Applied Mathematics, Computer Science and Statistics, Ghent University, Ghent, Belgium


#### Abstract

Spatial transcriptomics (ST) technologies allow the profiling of the transcriptome of cells while keeping their spatial context. Since most commercial untargeted ST technologies do not yet operate at single-cell resolution, computational methods such as deconvolution are often used to infer the cell type composition of each sequenced spot. We benchmarked 11 deconvolution methods using 63 silver standards, 3 gold standards, and 2 case studies on liver and melanoma tissues. We developed a simulation engine called synthspot to generate silver standards from singlecell RNA-sequencing data, while gold standards are generated by pooling single cells from targeted ST data. We evaluated methods based on their performance, stability across different reference datasets, and scalability. We found that cell 2 location and RCTD are the top-performing methods, but surprisingly, a simple regression model outperforms almost half of the dedicated spatial deconvolution methods. Furthermore, we observe that the performance of all methods significantly decreased in datasets with highly abundant or rare cell types. Our results are reproducible in a Nextflow pipeline, which also allows users to generate synthetic data, run deconvolution methods and optionally benchmark them on their dataset (https://github.com/saeyslab/spotless-benchmark).


## eLife assessment

This study makes a valuable contribution to spatial transcriptomics by rigorously benchmarking cell-type deconvolution methods, assessing their performance across diverse datasets with a focus on biologically relevant, previously unconsidered aspects. The authors demonstrate the strengths of RCTD, cell2location, and SpatialDWLS for their performance, while also revealing the limitations of many methods when compared to simpler baselines. By implementing a full Nextflow pipeline, Docker containers, and a rigorous assessment of the simulator, this work offers robust insights that elevate the standards for future evaluations and provides a resource for those seeking to improve or develop new deconvolution methods. The thorough comparison and analysis of methods, coupled with a strong emphasis on reproducibility, provide solid support for the findings.

## Introduction

Unlike single-cell sequencing, spatial transcriptome profiling technologies can uncover the location of cells, adding another dimension to the data that is essential for studying systems biology, for example cell-cell interactions and tissue architecture (Wagner et al., 2016). These approaches can be categorized into imaging- or sequencing-based methods, each of them offering complementary advantages (Asp et al., 2020). As imaging-based methods use labeled hybridization probes to target specific genes, they offer subcellular resolution and high capture sensitivity, but are limited to a few hundreds
or thousands of genes (Eng et al., 2019; Xia et al., 2019). On the other hand, sequencing-based methods offer an unbiased and transcriptome-wide coverage by using capture oligonucleotides with a polydT sequence (Rodriques et al., 2019; Ståhl et al., 2016). These oligos are printed in clusters, or spots, each with a location-specific barcode that allows identification of the originating location of a transcript. While the size of these spots initially started at $100 \mu \mathrm{~m}$ in diameter, some recent sequencing technologies have managed to reduce their size to the subcellular level, thus closing the resolution gap with imaging technologies (Chen et al., 2022; Cho et al., 2021). However, these spots do not necessarily correspond to individual cells, and therefore, computational methods remain necessary to determine the cell-type composition of each spot.

Deconvolution and mapping are two types of cell type composition inference tools that can be used to disentangle populations from a mixed gene expression profile (Longo et al., 2021). In conjunction with the spatial dataset, a reference scRNA-seq dataset from an atlas or matched sequencing experiment is typically required to build cell-type-specific gene signatures. Deconvolution infers the proportions of cell types in a spot by utilizing a regression or probabilistic framework, and methods specifically designed for ST data often incorporate additional model parameters to account for spot-to-spot variability. On the other hand, mapping approaches score a spot for how strongly its expression profile corresponds to those of specific cell types. As such, deconvolution returns the proportion of cell types per spot, and mapping returns the probability of cell types belonging to a spot. Unless otherwise specified, we use the term 'deconvolution' to refer to both deconvolution and mapping algorithms in this study.

Although there have recently been multiple benchmarking studies (Li et al., 2022; Yan and Sun, 2023; Li et al., 2023b), several questions remain unanswered. First, the added value of algorithms specifically developed for the deconvolution of ST data has not been evaluated by comparing them to a baseline or bulk deconvolution method. Second, it is unclear which algorithms are better equipped to handle challenging scenarios, such as the presence of a highly abundant cell type throughout the tissue or the detection of a rare cell type in a single region of interest. Finally, the stability of these methods to variations in the reference dataset arising from changes in technology or protocols has not been assessed.

In this study, we address these gaps in knowledge and provide a comprehensive evaluation of 11 deconvolution methods in terms of performance, stability, and scalability (Figure 1). The tools include eight spatial deconvolution methods (cell2location Kleshchevnikov and Shmatko, 2020, DestVI Lopez et all., 2022, DSTG Song and Su, 2021, RCTD Cable et al., 2022, SpatialDWLS Dong and Yuan, 2021, SPOTlight Elosua-Bayes et al., 2021, stereoscope Andersson et al., 2020, and STRIDE Sun et al., 2022), one bulk deconvolution method (MuSiC Wang et al., 2019), and two mapping methods (Seurat Stuart et al., 2019 and Tangram Biancalani et al., 2021). For all methods compared, we discussed with the original authors in order to ensure that their method was run in an appropriate setting and with good parameter values. We also compared method performance with two baselines: a 'null distribution' based on random proportions drawn from a Dirichlet distribution, and predictions from the non-negative least squares (NNLS) algorithm. We evaluated method performance on a total of 66 synthetic datasets ( 63 silver standards and 3 gold standards) and two application datasets. Our benchmarking pipeline is completely reproducible and accessible through Nextflow (github.com/ saeyslab/spotless-benchmark). Furthermore, each method is implemented inside a Docker container, which enables users to run the tools without requiring prior installation.

## Results

## Synthspot allows simulation of artificial tissue patterns

Synthetic spatial datasets are commonly generated by developers of deconvolution methods as part of benchmarking their algorithms against others. However, these synthetic datasets typically have spots with random compositions that do not reflect the reality of tissue regions with distinct compositions, such as layers in the brain. On the other hand, general-purpose simulators are more focused on other inference tasks, such as spatial clustering and cell-cell communication, and are usually unsuitable for deconvolution. For instance, generative models and kinetic models like those of scDesign3 (Song and Wang, 2024) and scMultiSim (Li et al., 2023a) are computationally intensive and unable to model entire transcriptomes. SRTSim (Zhu, 2023) focuses on modeling gene expression trends and does not

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-03.jpg?height=1567&width=1627&top_left_y=182&top_left_x=122)
Figure 1. Overview of the benchmark. (a) The datasets used consist of silver standards generated from single-cell RNA-seq data, gold standards from imaging-based data, and two case studies on liver and melanoma. Our simulation engine synthspot enables the creation of artificial tissue patterns. (b) We evaluated deconvolution methods on three overall performance metrics (RMSE, AUPR, and JSD), and further checked specific aspects of performance, that is how well methods detect rare cell types and handle reference datasets from different sequencing technologies. For the case studies, the AUPR and stability are only evaluated on the liver dataset. (c) Our benchmarking pipeline is entirely accessible and reproducible through the use of Docker containers and Nextflow. (d) To evaluate performance on the liver case study, we leveraged prior knowledge of the localization and composition of cell types to calculate the AUPR and JSD. We also investigated method performance on three different sequencing protocols.

The online version of this article includes the following figure supplement(s) for figure 1:
Figure supplement 1. Overview of synthspot abundance patterns used in the study.
Figure supplement 2. Data generation scheme of silver standards.
Figure supplement 3. UMAP and violin plots of three of the seven scRNA-seq datasets used to generate silver standards.
Figure supplement 4. (Continuation of the previous figure supplement) UMAP and violin plots of four of the seven scRNA-seq datasets used to generate silver standards.
explicitly model tissue composition, while spaSim (Feng et al., 2023) only models tissue composition without gene expression. To overcome these limitations, we developed our own simulator called synthspot that can generate synthetic tissue patterns with distinct regions, allowing for more realistic simulations (https://github.com/saeyslab/synthspot; Browaeys and Sang-aram, 2024). We validated that our simulation procedure accurately models real data characteristics and that method performances are comparable between synthetic and real tissue patterns (Appendix 1).

Within a synthetic dataset, synthspot creates artificial regions in which all spots belonging to the same region have the same frequency priors. Frequency priors correspond to the likelihood in which cells from a cell type will be sampled, and therefore, spots within the same region will have similar compositions. These frequency priors are influenced by the chosen artificial tissue pattern, or abundance pattern, which determines the uniformity, distinctness, and rarity of cell types within and across regions (Figure 1a, Figure 1-figure supplement 1). For example, the uniform characteristic will sample the same number of cells for all cell types in a spot, while diverse samples differing number of cells. The distinct characteristic constrains a cell type to only be present in one region, while overlap allows it to be present in multiple regions. Additionally, the dominant cell type characteristic randomly assigns a dominant cell type that is at least $5-15$ times more abundant than others in each spot, while rare cell type does the opposite to create a cell type that is $5-15$ times less abundant. The different characteristics can be combined in up to nine different abundance patterns, each representing a plausible biological scenario.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-04.jpg?height=1056&width=1597&top_left_y=1177&top_left_x=135)
Figure 2. Overall results of the benchmark. (a) Methods ordered according to their overall rankings (d), determined by the aggregated rankings of performance and scalability. (b) Performance of each method across metrics, artificial abundance patterns in the silver standard, and data sources. The ability to detect rare cell types and stability against different reference datasets are also included. (c) Average runtime across silver standards and scalability on increasing dimensions of the spatial dataset.
The online version of this article includes the following source data for figure 2:

Source data 1. Raw data table of Figure 2.

## Cell2location and RCTD are the top performers in synthetic data

Our synthetic spatial datasets consist of 63 silver standards generated from synthspot and 3 gold standards generated from imaging data with single-cell resolution (Supplementary file 1a-b). We generated the silver standards using seven publicly available scRNA-seq datasets and nine abundance patterns. The scRNA-seq datasets consisted of four mouse brain tissues (cortex, hippocampus, and two cerebellum), mouse kidney, mouse skin cancer (melanoma), and human skin cancer (squamous cell carcinoma). Half of the cells from each scRNA-seq dataset were used to generate the synthetic spots and the other half used as the reference for deconvolution. This split was stratified by cell type. We generated 10 replicates for each of the 63 silver standards, with each replicate containing around 750 spots (Figure 1-figure supplement 2). For the gold standard, we used two seqFISH+ sections of mouse brain cortex and olfactory bulb ( 63 spots with 10,000 genes each) and one STARMap section of mouse primary visual cortex ( 108 spots with 996 genes; Eng et al., 2019; Wang et al., 2018). We summed up counts from cells within circles of $55 \mu \mathrm{~m}$ diameter, which are the size of spots in the 10 x Visium commercial platform.

We evaluated method performance with the root-mean-square error (RMSE), area under the precision-recall curve (AUPR), and Jensen-Shannon divergence (JSD; Appendix 2). The RMSE measures how numerically accurate the predicted proportions are, the AUPR measures how well a method can detect whether a cell type is present or absent, and the JSD measure similarities between two probability distributions.

RCTD and cell2location were the top two performers across all metrics in the silver standards, followed by SpatialDWLS, stereoscope, and MuSiC (Figure 2b, Figure 3a). All other methods ranked worse than NNLS in at least one metric. For each silver standard, method rankings were determined using the median value across 10 replicates. We observed that method performances were more consistent between abundance patterns than between datasets (Figure 3-figure supplements 1-3). Most methods had worse performance in the two abundance patterns with a dominant cell type, and there was considerable performance variability between replicates due to different dominant cell types being selected in each replicate. Only RCTD and cell2location consistently outperformed NNLS in all metrics in these patterns (Figure 3-figure supplements 4-5).

For the gold standards, cell2location, MuSiC , and RCTD are the top three performers as well as the only methods to outrank NNLS in all three metrics (Figure 3b). As each seqFISH+ dataset consisted of seven field of views (FOVs), we used the average across FOVs as the representative value to be ranked for that dataset. Several FOVs were dominated by one cell type (Figure 3-figure supplement 6), which was similar to the dominant cell type abundance pattern in our silver standard. Consistent with the silver standard results, half of the methods did not perform well in these FOVs. In particular, SPOTlight, DestVI, stereoscope, and Tangram tended to predict less variation between cell type abundances. DSTG, SpatialDWLS, and Seurat predicted the dominant cell type in some FOVs but did not predict the remaining cell type compositions accurately. Most methods performed worse in the STARMap dataset except for SpatialDWLS, SPOTlight, and Tangram.

## Detecting rare cell types remains challenging even for top-performing methods

Lowly abundant cell types often play an important role in development or disease progression, as in the case of stem cells and progenitor cells or circulating tumor cells (Jindal et al., 2018). As the occurrence of these cell types are often used to create prognostic models of patient outcomes, the accurate detection of rare cell types is a key aspect of deconvolution tools (Ali et al., 2016; Sato et al., 2005). Here, we focus on the two rare cell type patterns in our silver standard (rare cell type diverse and regional rare cell type diverse), in which a rare cell type is $5-15$ times less abundant than other cell types in all or one synthetic region, respectively (Figure 1-figure supplement 1). This resulted in 14 synthetic datasets (seven scRNA-seq datasets × two abundance patterns) with 10 replicates each. We only evaluated methods based on the AUPR of the rare cell type, using the average AUPR across the 10 replicates as the representative value for each dataset. We did not include the RMSE and JSD or consider other cell types, because in prognostic models, the presence of rare cell types is often of more importance than the magnitude of the abundance itself. Therefore, it is more relevant that methods are able to correctly rank spots with and without the rare cell type.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-06.jpg?height=1157&width=1881&top_left_y=180&top_left_x=122)
Figure 3. Method performance on synthetic datasets, evaluated using root-mean-squared error (RMSE), area under the precision-recall curve (AUPR), and Jensen-Shannon divergence (JSD). Non-negative least squares (NNLS) is shaded as a baseline algorithm. Methods are ordered based on the summed ranks across all 63 and three datasets, respectively. (a) The rank distribution of each method across all 63 silver standards, based on the best median value across ten replicates for that standard. (b) Gold standards of two seqFISH+ datasets and one STARMap dataset. We took the average over seven field of views for the seqFISH+ dataset.

The online version of this article includes the following source data and figure supplement(s) for figure 3:
Source data 1. Raw data table of Figure 3.
Figure supplement 1. Boxplots of root-mean-squared error (RMSE) across ten replicates for each silver standard dataset (row) and abundance pattern (column).

Figure supplement 2. Boxplots of the area under the precision-recall curve (AUPR) across ten replicates for each silver standard dataset (row) and abundance pattern (column).
Figure supplement 3. Boxplots of Jensen-Shannon divergence across ten replicates for each silver standard dataset (row) and abundance pattern (column).

Figure supplement 4. Summed rank plots across all silver standard datasets for each abundance pattern (column) and metric (row).
Figure supplement 5. (Continuation of the previous figure supplement) Summed rank plots across all silver standard datasets for each abundance pattern (column) and metric (row).

Figure supplement 6. Summed abundances across all spots for each gold standard dataset.

In line with our previous analysis, RCTD and cell2location were also the best at predicting the presence of lowly abundant cell types (Figure 4a). There is a strong correlation between cell type abundance and AUPR, as clearly demonstrated when plotting the precision-recall curves of cell types with varying proportions (Figure 4b). While most methods can detect cell types with moderate or high abundance, they usually have lower sensitivity for rare cell types, and hence lower AUPRs. Upon visual inspection of precision-recall curves at decreasing abundance levels, we found a similar pattern across all silver standards (Figure 4-figure supplement 1). Nonetheless, we also observe that the AUPR is

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-07.jpg?height=860&width=1627&top_left_y=182&top_left_x=122)
Figure 4. Detection of the rare cell type in the two rare cell type abundance patterns. (a) Area under the precision-recall curve (AUPR) across the seven scRNA-seq datasets, averaged over 10 replicates. Methods generally have better AUPR if the rare cell type is present in all regions compared to just one region. (b) An example on one silver standard replicate demonstrates that most methods can detect moderately and highly abundant cells, but their performance drops for lowly abundant cells.

The online version of this article includes the following source data and figure supplement(s) for figure 4:
Source data 1. Raw data table of Figure 4.
Figure supplement 1. Evaluating the AUPR as a function of cell type abundance.
substantially lower if the rare cell type is only present in one region and not across the entire tissue, indicating that prevalence is also an important factor in detection.

## Technical variability between reference scRNA-seq and spatial datasets can significantly impact predictions

Since deconvolution predictions exclusively rely on signatures learned from the scRNA-seq reference, it should come as no surprise that the choice of the reference dataset has been shown to have the largest impact on bulk deconvolution predictions (Vallania et al., 2018). Hence, deconvolution methods should ideally also account for platform effects, that is the capture biases that may occur from differing protocols and technologies being used to generate scRNA-seq and ST data.

To evaluate the stability of each method against reference datasets from different technological platforms, we devised the inter-dataset scenario, where we provided an alternative reference dataset to be used for deconvolution, in contrast to the intra-dataset analysis done previously, where the same reference dataset was used for both spot generation and deconvolution. We tested this on the brain cortex (SMART-seq) and two cerebellum (Drop-seq and 10x Chromium) silver standards. For the brain cortex, we used an additional 10x Chromium reference from the same atlas (Yao et al., 2021). To measure stability, we computed the JSD between the proportions predicted from the intra- and inter-dataset scenario.

Except for MuSiC , we see that methods with better performance-cell2location, RCTD, SpatialDWLS, and stereoscope-were also more stable against changing reference datasets (Figure 5). Out of these four, only SpatialDWLS did not account for a platform effect in its model definition. Cell2location had the most stable predictions, ranking first in all three datasets, while both NNLS and MuSiC were consistently in the bottom three. For the rest of the methods, stability varied between datasets. DestVI performed well in the cerebellum datasets but not the brain cortex, and SPOTlight had the opposite

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-08.jpg?height=630&width=1883&top_left_y=180&top_left_x=120)
Figure 5. Prediction stability when using different reference datasets. For each synthetic dataset (total $\mathrm{n}=90$, from nine abundance patterns with ten replicates each), we computed the Jensen-Shannon divergence between cell type proportions obtained from two different reference datasets.
The online version of this article includes the following source data and figure supplement(s) for figure 5:

Source data 1. Raw data table of Figure 5.
Figure supplement 1. Changes in performance metrics when using a different reference dataset, from a matched to an unmatched reference (i.e., intradataset vs inter-dataset scenario).
pattern. As expected, deconvolution performance also generally worsens in the inter-dataset scenario (Figure 5-figure supplement 1).

## Evaluation of methods on public Visium datasets validate results on synthetic data

While synthetic datasets provide valuable benchmarks, it is crucial to validate deconvolution methods using real sequencing-based ST data, as they may exhibit distinct characteristics. To achieve this, we leveraged 10x Visium datasets from two mouse tissues, namely liver and skin cancer (melanoma). These datasets were chosen due to the availability of ground truth knowledge regarding cell proportions and localization, as elaborated further below.

## Liver

The liver is a particularly interesting case study due to its tissue pattern, where hepatocytes are highly abundant and constitute more than $60 \%$ of the tissue. This characteristic allows us to compare method performance with those of the dominant cell type abundance pattern from the silver standard, which was challenging for most methods. Here, we use the four Visium slides and single-cell dataset from the liver cell atlas of Guilliams et al., 2022; Supplementary file 1 c. The single-cell dataset was generated from three different experimental protocols-scRNA-seq following ex vivo digestion, scRNA-seq following in vivo liver perfusion, and single-nucleus RNA-seq (snRNA-seq) on frozen liver-additionally enabling us to assess method stability on different reference datasets.

We assessed method performance using AUPR and JSD by leveraging prior knowledge of the localization and composition of cell types in the liver. Although the true composition of each spot is not known, we can assume the presence of certain cell types in certain locations in the tissue due to the zonated nature of the liver. Specifically, we calculated the AUPR of portal vein and central vein endothelial cells (ECs) under the assumption that they are present only in their respective venous locations. Next, we calculated the JSD of the expected and predicted cell type proportions in a liver sample. The expected cell type proportions were based on snRNA-seq data (Figure 6-figure supplement 1), which has been shown to best recapitulate actual cell compositions observed by confocal microscopy (Guilliams et al., 2022). We ran each method on a combined reference containing all three experimental protocols (ex vivo scRNA-seq, in vivo scRNA-seq, and snRNA-seq), as well as on each protocol separately. To ensure consistency, we filtered each reference dataset to only include the nine cell types that were common to all three protocols.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-09.jpg?height=1430&width=1629&top_left_y=180&top_left_x=120)
Figure 6. Method performance on two Visium case studies. (a) In the liver case study, the AUPR was calculated using the presence of portal/central vein endothelial cells in portal and central veins, and the JSD was calculated by comparing predicted cell type proportions with those from snRNA-seq. All reference datasets contain nine cell types. Biological variation refers to the average pairwise JSD between four snRNA-seq samples. Methods are ordered based on the summed rank of all data points. (b) For melanoma, the JSD was calculated between the predicted cell type proportions and those from Molecular Cartography (bold). Biological variation refers to the JSD between the two Molecular Cartography samples. (c) Relationship between the proportions of endothelial cells predicted per spot and their distance to the nearest blood vessel (in arbitrary units, AU), where zero denotes a spot annotated as a vessel. An inverse correlation can be observed more clearly in better-performing methods.

The online version of this article includes the following source data and figure supplement(s) for figure 6:
Source data 1. Raw data table of Figure 6.
Figure supplement 1. Comparison of cell type compositions between three sequencing protocols in the mouse liver atlas from Guilliams et al., 2022.
Figure supplement 2. Predicted cell type abundances averaged across all four Visium slides from the liver atlas.
Figure supplement 3. Concordance of method performance between the synthetic dataset (generated using synthspot's dominant cell type abundance pattern) and the Visium dataset.
Figure supplement 4. The predicted abundance of central vein and portal vein endothelial cells (ECs) for each spot in one liver Visium slide.
Figure supplement 5. Stability of predicted proportions when using three different protocols from the liver atlas as reference for deconvolution.
Figure supplement 6. The ground truth and predicted cell type proportions of the melanoma dataset.

RCTD and cell2location were the top performers for both AUPR and JSD (Figure 6a), achieving a JSD comparable to those of biological variation, that is the average pairwise JSD between four snRNA-seq samples. In contrast, Tangram, DSTG, SPOTlight, and stereoscope had higher JSD values than those of NNLS. Except for SPOTlight, these three methods were not able to accurately infer the overabundance of the dominant cell type, with stereoscope and Tangram predicting a uniform distribution of cell types (Figure 6-figure supplement 2). This corresponds to what we observed in the dominant cell type pattern of the silver standard and certain FOVs in the gold standard. To further substantiate this, we used the ex vivo scRNA-seq dataset to generate 10 synthetic datasets with the dominant cell type abundance pattern. Then, we used the remaining two experimental protocols separately as the reference for deconvolution. Indeed, the method rankings followed the overall pattern, and a comparison of the JSD values between synthetic and real Visium data revealed a strong correlation (Figure 6-figure supplement 3).

Nonetheless, there were some inconsistencies between the AUPR and JSD rankings, specifically with Seurat and SpatialDWLS ranking highly for JSD and lowly for AUPR. This is because the JSD is heavily influenced by the dominant cell type in the dataset, such that even when predicting only the dominant cell type per spot, Seurat still performed well in terms of JSD. However, it was unable to predict the presence of portal and central vein ECs in their respective regions (Figure 6-figure supplement 4). Therefore, complementary metrics like AUPR and JSD must both be considered when evaluating methods.

We observed that using the combination of the three experimental protocols as the reference dataset did not necessarily result in the best performance, and it was often possible to achieve better or similar results by using a reference dataset from a single experimental protocol. The best reference varied between methods, and most methods did not exhibit consistent performance across all references. Interestingly, Cell2location, MuSiC, and NNLS had much higher JSD when using the snRNA-seq data as the reference, while RCTD and Seurat had the lowest JSD on the same reference. To further evaluate the stability of the methods, we calculated the JSD between proportions predicted with different reference datasets. RCTD and Seurat showed the lowest JSD, indicating higher stability (Figure 6-figure supplement 5). Finally, we examined the predicted proportions when using the entire atlas without filtering cell types, which contains all three protocols and 23 cell types instead of only the nine common cell types (Figure 6-figure supplement 2). The additional 14 cell types made up around $20 \%$ of the ground truth proportions. While RCTD, Seurat, SpatialDWLS, and MuSiC retained the relative proportions of the nine common cell types, the rest predicted substantially different cell compositions.

## Melanoma

Melanoma poses a significant challenge in both therapeutic and research efforts due to its high degree of heterogeneity and plasticity. In a recent study, Karras et al., 2022 investigated the cell state diversity of melanoma by generating single-cell and ST data from a mouse melanoma model (Supplementary file 1d). Among others, the spatial data consists of three individual tumor sections profiled by $10 \times$ Visium, as well as 33 regions of interest from two tumor samples profiled by Molecular Cartography (MC), an imaging-based technology that can profile up to 100 genes at subcellular resolution. Using a custom processing pipeline, we obtained cell type annotations-and subsequently the cell type proportions-for each of the MC samples. These cell type proportions were consistent across different sections and samples, and were used as the ground truth for deconvolution (Figure 6figure supplement 6). We aggregated the predicted proportions of the seven malignant cell states in the Visium slides, as we could not reliably annotate these cell states in the MC dataset.

To assess method performance, we calculated the average pairwise JSDs between the two MC samples and three Visium slides. Cell2location, SPOTlight, and RCTD were the best performers with JSDs of around 0.01 (Figure 6b). With the exception of SPOTlight performing well in this dataset, the rankings of the remaining methods followed the usual trend in the silver standards, with NNLS outperforming Seurat, DestVI, STRIDE, and DSTG. Additionally, we sought to corroborate these findings through a more qualitative approach by leveraging the blood vessel annotations provided in the original study. Given that by definition, endothelial cells form the linings of blood vessels, we visualized the relationship between EC abundance and the distance to the nearest vessel (Figure 6c). Although the majority of spots were predicted to have no ECs, a fraction of spots exhibited the expected

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-11.jpg?height=897&width=1879&top_left_y=180&top_left_x=122)
Figure 7. Runtime and scalability of each method. (a) Runtime over the 63 silver standards (three replicates each). Methods are ordered by total runtime. Asterisks indicate when GPU acceleration has been used. Cell2location, stereoscope, DestVI, and STRIDE first build a model for each single-cell reference (red points), which can be reused for all synthetic datasets derived from that reference. (b) Method scalability on increasing dimensions of the spatial dataset. For model-based methods, the model building and fitting time were summed. Methods are ordered based on total runtime.

The online version of this article includes the following source data for figure 7:
Source data 1. Raw data table of Figure 7.
trend where ECs were more abundant the closer the spot was to a blood vessel. This trend was more discernible for higher ranked methods, while the lower-ranked ones either showed no correlation (NNLS, DestVI, DSTG) or had noisy predictions (Seurat, Tangram).

## Runtime and scalability

Most methods are able to deconvolve the silver standard datasets in less than 30 min and only had slight variability in the runtime (Figure 7a). Model-based methods-DestVI, stereoscope, cell2location, and STRIDE-have the advantage that a model built from a single-cell reference can be reused for all synthetic datasets derived from that reference (i.e. the nine abundance patterns $\times 10$ replicates). This typically reduces the runtime needed to fit the model on the synthetic datasets. However, these methods tend to be more computationally intensive and are strongly recommended to be run with GPU acceleration. This was not implemented in STRIDE and could explain its longer runtime during model building.

Next, we investigated method scalability by varying the dimensions of the spatial dataset (Figure 7b). We tested 16 combinations of spots ( $100,1 \mathrm{k}, 5 \mathrm{k}, 10 \mathrm{k}$ ) and genes ( $5 \mathrm{k}, 10 \mathrm{k}, 20 \mathrm{k}, 30 \mathrm{k}$ ), while the single-cell dataset was kept at 5 k cells and 30 k genes. Here, we considered both the model building and fitting step in the runtime (when applicable). Tangram, Seurat, and SPOTlight had a constant runtime across all dimensions, each deconvolving the largest dataset ( $3 \times 10^{8}$ elements in total) in less than ten minutes. Other methods had increased runtime both when spots and genes were increased, except for DestVI which was only affected by the number of spots, and STRIDE by the number of genes. This was because by default, DestVI only uses 2000 highly variable genes for deconvolution. Similarly, the scalability of RCTD, SpatialDWLS, Tangram, and SPOTlight are due to the fact that they only make use of differentially expressed or cell type-specific marker genes. Stereoscope was the least scalable by far, taking over 8 hr for the largest dataset, and 3.5 hr for the smallest one. Note that many methods allow for parameters that can greatly affect runtime, such as the number of training epochs and/or the number of genes to use (Supplementary file 2). For instance, here we
have used all genes to run stereoscope (default parameters), but the authors have suggested that it is possible to use only 5000 genes in the model and maintain the same performance.

## Discussion

In this study, we performed a thorough investigation of various biologically relevant and previously unconsidered aspects related to deconvolution. We evaluated the performance of 11 deconvolution methods on 63 silver standards, 3 gold standards, and 2 Visium datasets using three complementary metrics. We also incorporated two baselines for every analysis, the NNLS algorithm and null distribution proportions from a Dirichlet distribution. In the case studies, we demonstrated two approaches for evaluating deconvolution methods in datasets lacking an absolute ground truth. These approaches include using proportions derived from another sequencing or spatial technology as a proxy, and leveraging spot annotations, for example zonation or blood vessel annotations, that typically have already been generated for a separate analysis.

Our findings indicate that RCTD, cell2location, and SpatialDWLS were the highest ranked methods in terms of performance, consistent with previous studies (Li et al., 2022; Yan and Sun, 2023). However, we also found that over half of the methods did not outperform the baseline (NNLS) and bulk deconvolution method ( MuSiC ). These results were consistent across the silver and gold standards, as well as the liver and melanoma case studies, demonstrating the generalizability and applicability of our simulation and benchmarking framework. We also found that the abundance pattern of the tissue and the reference dataset used had the most significant impact on method performance. Even top-performing methods struggled with challenging abundance patterns, such as when a rare cell type was present in only one region, or when a highly dominant cell type masks the signature of less abundant ones. Furthermore, different reference datasets could result in substantially different predicted proportions. Methods that accounted for technical variability in their models, such as cell2location and RCTD, were more stable to changes in the reference dataset than those that did not, such as SpatialDWLS.

Regarding the reference dataset, the number of genes per cell type (which is generally correlated to the sequencing depth) seems to have a significant impact on the performance of deconvolution methods. We observed that methods were less accurate in datasets with fewer genes per cell type (Figure 1-figure supplements 3-4). For example, all methods performed better in the snRNA-seq cerebellum dataset, which had the same number of cell types as the scRNA-seq cerebellum dataset, but on average 1000 more genes per cell type. The kidney dataset was the most challenging for all methods, with most of its 16 cell types having less than 1000 genes. This was evident from the RMSE and JSD scores that were relatively closer to the null distribution than in other datasets (Figure 3figure supplements 1 and 3). In contrast, the 18 cell types in the brain cortex dataset had an average of 3000-9000 features, leading to much better performance for most methods compared to the kidney dataset despite having more cell types. This trend was also observed in the STARMap gold standard dataset, which consisted of only 996 genes. Most methods performed worse in the STARMap dataset except for SpatialDWLS, SPOTlight, and Tangram. Since these three methods only use marker genes for deconvolution, this may explain why the small number of genes (most of which were already marker genes) did not affect them as much.

In addition to performance, the runtime and scalability of a method is also a factor to consider. Although most runtimes were comparable on our silver standards, we have made use of GPU acceleration for Tangram, DestVI, stereoscope, and cell2location. As this might not be available for some users, using these methods on the CPU might require training on fewer epochs or selecting only a subset of genes. With all these factors in consideration, we recommend RCTD as a good starting point. In addition to being one of the best and fastest methods, it also allows CPU parallelization (i.e. using multiple cores) for users that may not have access to a GPU.

As a general guideline, we recommend comparing the result of multiple deconvolution methods, especially between cell2location and RCTD. If the predictions are highly contradictory, the reference dataset may not be of sufficiently good quality. We recommend obtaining scRNA-seq and spatial data from the same sample to reduce biological variability, as there will always be technical variability across platforms. This can also ensure that the same cell types will be present in both datasets (Longo et al., 2021). If that is not possible, use a reference dataset that has sufficient sequencing depth (at least more than 1000 genes per cell type), preferably from a single platform. In addition to checking
the sequencing depth, our simulator synthspot can also be used to evaluate the quality of the reference dataset. As we have demonstrated in the liver case study, users can generate synthetic spatial datasets with an abundance pattern that best resembles their tissue of interest. With a high-quality reference, both cell2location and RCTD should be able to achieve an AUPR close to one and JSD close to zero. Our Nextflow pipeline seamlessly integrates the complete workflow of synthetic data generation, deconvolution, and performance evaluation.

As spatial omics is still an emerging field, the development of new deconvolution methods can be anticipated in the future. Our benchmark provides a reliable and reproducible evaluation framework for current and upcoming deconvolution tools, making it a valuable resource for researchers in the field.

## Methods

## Synthspot

We generated synthetic spatial datasets using the synthspot R package, whose spot generation process is modeled after that of SPOTlight (Elosua-Bayes et al., 2021). This simulator generates synthetic spot data by considering the gene expression of one spot as a mixture of the expression of $n$ cells, with $n$ being a random number between 2 and 10 that is sampled from a uniform distribution. To generate one spot, the simulator samples $n$ cells from the input scRNA-seq data, sums their counts gene by gene, and downsamples the counts. The target number of counts for downsampling is picked from a normal distribution with mean and standard deviation of $20,000 \pm 5000$ by default, but these values can be changed by the user. To mimic biological tissue, synthspot generates artificial regions, or groups of spots with similar cell type compositions (Figure 1-figure supplement 1). The prior distribution of each cell type in each region is influenced by the selected abundance pattern, called dataset type in the package (Appendix 1).

## Method execution

An overview of the 11 methods can be found in Supplementary file 2. As input, all methods require a reference scRNA-seq dataset with cell type annotations along with a spatial dataset. We first ran the methods based on the provided tutorials, and later reached out to the original authors of each method for additional recommendations. For most methods, we explored different parameters and selected ones that resulted in the best performance. The specifics on how each method was run and the final parameters can be found in Appendix 3. Unless stated otherwise, these parameters were applied on all datasets used in this study.

We implemented a Nextflow pipeline (Di Tommaso et al., 2017) to run the methods concurrently and compute their performance. For reproducibility, each method was installed inside a Docker container and can be run either using Docker or Apptainer (formerly known as Singularity). Our pipeline can be run in three modes: (1) run_standard to reproduce our benchmark, (2) run_dataset to run methods on a given scRNA-seq and spatial dataset, and (3) generate_and_run to also generate synthetic spatial data from the given scRNA-seq dataset. Our pipeline can accept both Seurat (.rds) and AnnData (.h5ad) objects as input, and the proportions are output as a tab-separated file. Parameter tuning has also been implemented.

All methods were deployed on a high-performance computing cluster with an Intel Xeon Gold 6140 processor operating at 2.3 GHz , running a RHEL8 operating system. We made use of GPU acceleration whenever possible using a NVIDIA Volta V100 GPU. We ran all methods with one core and 8 GB memory, dynamically increasing the memory by 8 GB if the process failed.

## Datasets

## Gold

The seqFISH+ dataset consists of two mouse brain tissue slices (cortex and olfactory bulb) with seven field of views (FOVs) per slice (Eng et al., 2019). Ten thousand genes were profiled for each FOV. Each FOV had a dimension of $206 \mu \mathrm{~m} \times 206 \mu \mathrm{~m}$ and a resolution of 103 nm per pixel. We simulated Visium spots of $55 \mu \mathrm{~m}$ diameter and disregarded the spot-to-spot distance, resulting in nine spots per FOV and 126 spots for the entire dataset. The FOVs had varying number of cells and cell types per simulated spot (Supplementary file 1a). We created a reference dataset for each tissue slice using the
combined single-cell expression profiles from all seven FOVs. There were 17 cell types in the cortex tissue section and nine cell types in the olfactory bulb tissue section.

The STARMap dataset of mouse primary visual cortex (VISp) is $1.4 \mathrm{~mm} \times 0.3 \mathrm{~mm}$ and contains 1020 genes and 973 cells (Wang et al., 2018). We generated Visium-like spots as above, assuming that 1 pixel $=810 \mathrm{~nm}$, as this information was not provided in the original paper. We removed any spot with hippocampal cell types, and removed unannotated and Reln cells from spots, resulting in 108 spots comprising 569 cells and 12 cell types in total. We used the VISp scRNA-seq dataset from the Allen Brain Atlas as the reference (Tasic et al., 2016). After filtering for common genes, 996 genes remained.

## Silver

We used six scRNA-seq datasets and one snRNA-seq dataset for the generation of silver standards (Supplementary file 1b). All datasets are publicly available and already contain cell type annotations. We downsized each dataset by filtering out cells with an ambiguous cell type label and cell types with less than 25 cells, and kept only highly variable genes (HVGs) that are expressed in at least $10 \%$ of cells from one cell type. We set this threshold at $25 \%$ for the brain cortex dataset. We further downsampled the kidney dataset to 15,000 cells, and the melanoma dataset to 20,000 cells. For the two cerebellum datasets, we integrated them following the method of Stuart et al., 2019, performed the same gene filtering approach on the integrated dataset, and only retained common cell types.

We generated the synthetic datasets using nine out of the 17 possible abundance patterns in synthspot. For each scRNA-seq dataset, we generated 10 replicates of each abundance pattern. For each replicate, we ran generate_synthetic_visium with five regions (n_regions), with minimally 100 spots and maximally 200 spots per region (n_spots_min, n_spots_max), and target counts of $20,000 \pm 5,000$ per spot (visium_mean, and visium_sd). There were on average 750 spots per replicate.

## Liver

We downloaded single-cell and spatial datasets of healthy mice from the liver cell atlas (Guilliams et al., 2022; Supplementary file 1c). The single-cell data contains 185,894 cells and 31,053 genes from three experimental protocols: scRNA-seq following ex vivo digestion, scRNA-seq following in vivo liver perfusion, and snRNA-seq on frozen liver. We used the finer annotation of CD45 cells, which among others, subclustered endothelial cells into portal vein, central vein, lymphatic, and liver sinusoidal endothelial cells. We only retained cell types where at least 50 cells are present in each protocol, resulting in nine common cell types. The spatial data consisted of four Visium slides, containing on average 1440 spots. Each spot has been annotated by the original authors as belonging to either the central, mid, periportal, or portal zone. This zonation trajectory was calculated based on hepatocyte zonation markers.

We deconvolved the Visium slides using five variations of the single-cell reference: the entire dataset, the dataset filtered to nine cell types, and each experimental protocol separately (also filtered to nine cell types). The proportions obtained from using the entire dataset was not used to compute evaluation metrics but only to visualize method stability (Figure 6-figure supplement 2). Due to the large size of the reference, we ran some methods differently. We ran stereoscope with 5000 HVGs and subsampled each cell type to a maximum of 250 cells (-sub). We gave STRIDE the number of topics to test (--ntopics 2333435363 ) instead of the default range that goes up to triple the number of cell types. We ran DestVI with 2500 training epochs and batch size of 128. For Seurat, we ran FindTransferAnchors with the 'rpca' (reciprocal PCA) option instead of 'pcaproject'. MuSiC and SpatialDWLS use dense matrices in their code (in contrast to sparse matrices in other methods), resulting in a size limit of $2^{31}$ elements in a matrix. Whenever this is exceeded, we downsampled the reference dataset by keeping maximally 10,000 cells per cell type and keeping 3000 HVGs along with any genes that are at least $10 \%$ expressed in a cell type. When a reference with multiple protocols was used, we provided this information to cell2location and Seurat.

To compute the AUPR, we considered only spots annotated as being in the central or portal zone. We considered the following ground truth for portal vein and central vein ECs:

## Continued on next page

|  | Portal vein EC | Central Vein EC |
| :--- | :--- | :--- |
| Portal spot | 1 | 0 |
| Central spot | 0 | 1 |

For the JSD ground truth, we only retained samples from the snRNA-seq protocol where all nine cell types were present, resulting in four samples (ABU11, ABU13, ABU17, and ABU20). For both the ground truth and predicted proportions, we averaged the abundance of each cell type in per slide or sample. Then, we calculated the pairwise JSD between each of the four slides and four snRNA-seq samples and reported the average JSD. The biological variation was obtained by averaging the pairwise JSD between the snRNA-seq samples.

## Melanoma

The scRNA-seq and spatial datasets were downloaded from the original study (Karras et al., 2022; Supplementary file 1d), with the scRNA-seq dataset being the same one used in the silver standard before the split (Supplementary file 1b). We preprocessed the Visium slides by filtering out spots with fewer than 100 features, as their inclusion led to errors when running marker-based deconvolution methods. This filtering removed 26 spots from the third slide. The blood vessel annotation was provided by the pathologist from the original publication. We ran stereoscope and DestVI with the same parameters as the liver case study due to the size of the reference dataset.

The Molecular Cartography (MC) dataset consists of two samples, six sections, and 33 regions of interest (ROI). For each ROI, a DAPI-stained image and a text file containing transcript locations are provided. We obtained cell type proportions using an in-house image processing pipeline. First, the image was cleaned by eliminating tiling effects and background debris using BaSic (Peng et al., 2017) and scikit-image (van der Walt et al., 2014), respectively. Subsequently, DAPI-stained images were segmented, and locations of individual nuclei were obtained using CellPose (Stringer et al., 2021). The transcripts of all measured genes were then assigned to their corresponding cells, resulting in the creation of cell-by-gene count matrices. These matrices were normalized based on the size of segmented nuclei and preprocessed in ScanPy (Wolf et al., 2018). Specifically, counts were logtransformed, scaled, and genes expressed in fewer than five cells and cells with less than 10 transcripts were filtered out. Leiden clustering (Traag et al., 2019) was performed using 17 principal components and 35 neighbors, and cells were annotated using scanpy.tl.score_genes with a curated marker gene list. Finally, the counts of each cell type were aggregated across samples to obtain cell type proportions.

As dendritic cells and melanoma cell states could not be annotated in the MC dataset, we adjusted the predicted proportions from Visium by removing dendritic cells proportions (pDCs and DCs) and recalculating the relative proportions, and aggregating the proportions of the seven malignant cell states.

## Scalability

We generated a synthetic dataset with 10,000 spots and 31,053 genes using only the snRNA-seq protocol from the liver atlas. We then used the remaining two protocols (combined) as the reference dataset. Genes of the synthetic dataset were downsampled randomly based on a uniform distribution, while genes of the reference data were downsampled based on marker genes and HVGs. The spots/ cells of both datasets were downsampled randomly, and this was done in a stratified manner for the reference dataset.

## Evaluation metrics and baselines

The root-mean-squared error between the known and predicted proportions ( $p$ ) of a spot $s$ for cell type $z$, in a reference dataset with $Z$ cell types in total, is calculated as

$$
\operatorname{RMSE}\left(s_{\text {known }}, s_{\text {predicted }}\right)=\sqrt{\frac{1}{Z} \sum_{z}^{Z}\left(p_{s z, \text { known }}-p_{s z, \text { predicted }}\right)^{2}}
$$

We calculated the JSD and AUPR using the R packages philentropy and precrec, respectively (Saito and Rehmsmeier, 2017; Drost, 2018). The JSD is a smoothed and symmetric version of the KullbackLeibler divergence (KL). It is calculated as

$$
J S D\left(s_{\text {known }} \| s_{\text {predicted }}\right)=\frac{1}{2} K L\left(p_{s, \text { known }} \| M\right)+\frac{1}{2} K L\left(p_{s, \text { predicted }} \| M\right),
$$

where

$$
M=\frac{1}{2}\left(p_{s, \text { known }}+p_{s, \text { predicted }}\right)
$$

To calculate the AUPR, we binarized the known proportions by considering a cell type to be present in a spot if its proportion is greater than zero, and absent if it is equal to zero. Then, we compute the micro-averaged AUPR by aggregating all cell types together.

For the silver standards, the RMSE and JSD values across all spots are averaged and used as the representative value for a replicate $k$ for each dataset-abundance pattern combination. In contrast, only one AUPR value is obtained per replicate. Note that the RMSE calculation takes the number of cell types into account and hence should not be compared between datasets.

We created the performance summary plot using the funkyheatmap package in R . To aggregate per data source and per abundance pattern, we calculated the grand mean of each metric across datasets, applied min-max scaling on each metric, and then computed the geometric mean of the scaled metrics (RMSE, AUPR, and JSD for gold and silver standards; AUPR and JSD for liver; and only JSD for melanoma). The aggregation per metric is the arithmetic mean of all datasets evaluated on that metric. Finally, we determined the overall rankings based on the weighted ranks of the following criteria: silver standard, gold standard, liver case study, melanoma case study, rare cell type detection, stability, and runtime. We assigned a weight of 0.5 to each silver standard abundance pattern and a weight of one to the remaining criteria.

To provide reference values for each metric, we used (1) random proportions based on probabilities drawn from a Dirichlet distribution and (2) predictions from the NNLS algorithm. For the former, we used the DirichletReg package (Maier, 2014) to generate reference values for all 63 silver standards using the average value across 100 iterations. The dimension of the $\boldsymbol{\alpha}$ vector was equal to the number of cell types in the corresponding dataset, and all concentration values were set to one. For the NNLS baseline, we used the Lawson-Hanson NNLS implementation from the nnls R package. We solve for $\beta$ in the equation $\boldsymbol{Y}=\boldsymbol{X} \beta$, with $\boldsymbol{Y}$ the spatial expression matrix and $\boldsymbol{X}$ the average gene expression profile per cell type from the scRNA-seq reference. We obtained proportion estimates by dividing each element of $\beta$ with its total sum.

## Acknowledgements

We thank all the authors of the methods who provided valuable feedback on how to optimally run their algorithms. Their input has been crucial in enabling us to compare the methods in a fair manner. We thank Lotte Pollaris for processing the Molecular Cartography dataset. CS is funded by the Ghent University Special Research Fund [grant number BOF21-DOC-105], RB is funded by The Research Foundation - Flanders [grant number 1181318 N], RS is funded by the Flemish Government under the Flanders AI Research Program, and YS is funded by Ghent University Special Research Fund [grant number BOF18-GOA-024], and The Research Foudation - Flanders [Excellence of Science (EOS) program and SBO project, grant number S001121N].

## Additional information

Funding
| Funder | Grant reference number | Author |
| :--- | :--- | :--- |
| Bijzonder Onderzoeksfonds UGent | BOF21-DOC-105 | Chananchida Sang-aram |
| Fonds Wetenschappelijk Onderzoek | 1181318N | Robin Browaeys |
| Vlaamse Overheid | Onderzoeksprogramma Artificiele Intelligentie | Ruth Seurinck |
| Fonds Wetenschappelijk Onderzoek | EOS (Excellence of Science) | Yvan Saeys |
| Bijzonder Onderzoeksfonds UGent | BOF18-GOA-024 | Yvan Saeys |
| Fonds Wetenschappelijk Onderzoek | SBO - S001121N | Yvan Saeys |


The funders had no role in study design, data collection and interpretation, or the decision to submit the work for publication.

## Author contributions

Chananchida Sang-aram, Conceptualization, Data curation, Software, Formal analysis, Validation, Investigation, Visualization, Methodology, Writing - original draft, Writing - review and editing; Robin Browaeys, Conceptualization, Data curation, Software, Methodology, Writing - review and editing; Ruth Seurinck, Conceptualization, Supervision, Methodology, Project administration, Writing - review and editing; Yvan Saeys, Conceptualization, Resources, Supervision, Funding acquisition, Project administration, Writing - review and editing

## Author ORCIDs

Chananchida Sang-aram (BD) https://orcid.org/0000-0002-0922-0822
Robin Browaeys (iD) https://orcid.org/0000-0003-2934-5195
Yvan Saeys ([10) https://orcid.org/0000-0002-0415-1506

## Peer review material

Reviewer \#1 (Public Review): https://doi.org/10.7554/eLife.88431.3.sa1
Reviewer \#2 (Public Review): https://doi.org/10.7554/eLife.88431.3.sa2
Author response https://doi.org/10.7554/eLife.88431.3.sa3

## Additional files

## Supplementary files

- Supplementary file 1. An overview of all datasets used in this study, including the gold standards, silver standards, liver case study, and melanoma case study.
- Supplementary file 2. An overview of the deconvolution tools benchmarked in this study, including the algorithm and usage information.
- MDAR checklist


## Data availability

All datasets used in this article, including the silver standards, gold standards, and case studies, are available on Zenodo at https://doi.org/10.5281/zenodo.5727613. Original download links and accession numbers of individual studies can be found at Supplementary file 1. The synthspot R package can be downloaded from https://github.com/saeyslab/synthspot, copy archived at Browaeys and Sang-aram, 2024. The Nextflow pipeline along with analysis scripts can be found at https://github. com/saeyslab/spotless-benchmark, copy archived at Sang-aram, 2023.

The following dataset was generated:
| Author(s) | Year | Dataset title | Dataset URL | Database and Identifier |
| :--- | :--- | :--- | :--- | :--- |
| Sang-aram C | 2024 | Benchmark datasets for Spotless | https://doi.org/10. 5281/zenodo. 5727613 | Zenodo, 10.5281/ zenodo. 5727613 |


The following previously published datasets were used:

| Author(s) | Year | Dataset title | Dataset URL | Database and Identifier |
| :--- | :--- | :--- | :--- | :--- |
| Wang X, Allen WE, Wright MA, Sylwestrak EL, Samusik N, Vesuna S, Evans K, Liu C, Ramakrishnan C, Liu J, Nolan GP, Bava FA, Deisseroth K | 2018 | Gold standard STARMap | https://www. dropbox.com/sh/ f7ebheru1lbz91s/ AACIAqjvDv--mhnE-PmSXB41a/ visual_1020?dl=0 | Dropbox, visual_1020 |
| Tasic B, Yao Z, Graybuck LT, Smith KA, Nguyen TN, Bertagnolli D, Goldy J, Garren E, Economo MN, Viswanathan S, Penn O | 2018 | Shared and distinct transcriptomic cell types across neocortical areas | https://www.ncbi. nlm.nih.gov/geo/ query/acc.cgi?acc= GSE115746 | NCBI Gene Expression Omnibus, GSE115746 |
| Saunders A, Macosko EZ, Wysoker A, Goldman M, Krienen FM, de Rivera H, Bien E, Baum M, Bortolin L, Wang S, Goeva A | 2018 | A Single-Cell Atlas of Cell Types, States, and Other Transcriptional Patterns from Nine Regions of the Adult Mouse Brain | https://www.ncbi. nlm.nih.gov/geo/ query/acc.cgi?acc= GSE116470 | NCBI Gene Expression Omnibus, GSE116470 |
| Kozareva V, Martin C, Osorno T, Rudolph S, Guo C, Vanderburg C, Nadaf N, Regev A, Regehr WG, Macosko E | 2021 | A transcriptomic atlas of mouse cerebellar cortex reveals novel cell types | https://www.ncbi. nlm.nih.gov/geo/ query/acc.cgi?acc= GSE165371 | NCBI Gene Expression Omnibus, GSE165371 |
| Zeisel A, Hochgerner H, Lönnerberg P, Johnsson A, Memic F, Van Der Zwan J, Häring M, Braun E, Borm LE, La Manno G, Codeluppi S | 2018 | Silver standard mouse hippocampus | http://mousebrain. org/adolescent/ downloads.html | Linarsson Lab Mouse Brain Atlas, adolescent |
| Park J, Shrestha R, Qiu C, Kondo A, Huang S, Werth M, Li M, Barasch J, Suszták K | 2018 | Comprehensive single cell RNAseq analysis of the kidney reveals novel cell types and unexpected cell plasticity | https://www.ncbi. nlm.nih.gov/geo/ query/acc.cgi?acc= GSE107585 | NCBI Gene Expression Omnibus, GSE107585 |
| Continued on next page |  |  |  |  |


Continued
| Author(s) | Year | Dataset title | Dataset URL | Database and Identifier |
| :--- | :--- | :--- | :--- | :--- |
| Karras P, Bordeu I, Pozniak J, Nowosad A, Pazzi C, Van Raemdonck N, Landeloos E, Van Herck Y, Pedri D, Bervoets G, Makhzami S | 2022 | A cellular hierarchy in melanoma uncouples growth and metastasis | https://www.ncbi. nlm.nih.gov/geo/ query/acc.cgi?acc= GSE207592 | NCBI Gene Expression Omnibus, GSE207592 |
| Ji AL, Rubin AJ, Thrane K, Jiang S, Reynolds DL, Meyers RM, Guo MG, George BM, Mollbrink A, Bergenstråhle J, Larsson L | 2022 | Single Cell and Spatial Analysis of Human Squamous Cell Carcinoma [single-cell RNA-seq] | https://www.ncbi. nlm.nih.gov/geo/ query/acc.cgi?acc= GSE144236 | NCBI Gene Expression Omnibus, GSE144236 |
| Yao Z, Van Velthoven CT, Nguyen TN, Goldy J, Sedeno-Cortes AE, Baftizadeh F, Bertagnolli D, Casper T, Chiang M, Crichton K, Ding SL | 2021 | Stability analysis mouse brain cortex | https://portal. brain-map.org/ atlases-and-data/ rnaseq/mouse-whole-cortex-and-hippocampus-10x | Allen Brain Atlas, Mouse-whole-cortex-and-hippocampus-10x |
| Guilliams M, Bonnardel J, Haest B, Vanderborght B, Wagner C, Remmerie A, Bujko A, Martens L, Thoné T, Browaeys R, De Ponti FF | 2022 | Case study - mouse liver | https://www. livercellatlas.org/ download.php | Scott Lab \& Guilliams Lab Liver Cell Atlas, Mouse-StSt |
| Karras P, Bordeu I, Pozniak J, Nowosad A, Pazzi C, Van Raemdonck N, Landeloos E, Van Herck Y, Pedri D, Bervoets G, Makhzami S | 2022 | Case study - mouse melanoma (10x Visium) | https://drive. google.com/drive/ folders/1poq4Lo5 AxVp0WpG1EMg ljleDR4q98zcA | Google Drive, 1poq4Lo5AxVp0WpG1EMgljleDR4q98zcA |
| Karras P, Bordeu I, Pozniak J, Nowosad A, Pazzi C, Van Raemdonck N, Landeloos E, Van Herck Y, Pedri D, Bervoets G, Makhzami S | 2022 | A cellular hierarchy in melanoma uncouples growth and metastasis | https://doi.org/ 10.5281/zenodo. 6856193 | Zenodo, 10.5281/zenodo. 6856193 |


## References

Ali HR, Chlon L, Pharoah PDP, Markowetz F, Caldas C. 2016. Patterns of immune infiltration in breast cancer and their clinical implications: A gene-expression-based retrospective study. PLOS Medicine 13:e1002194. DOI: https://doi.org/10.1371/journal.pmed.1002194, PMID: 27959923

Andersson A, Bergenstråhle J, Asp M, Bergenstråhle L, Jurek A, Fernández Navarro J, Lundeberg J. 2020. Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography. Communications Biology 3:565. DOI: https://doi.org/10.1038/s42003-020-01247-y, PMID: 33037292
Asp M, Bergenstråhle J, Lundeberg J. 2020. Spatially resolved transcriptomes-next generation tools for tissue exploration. BioEssays 42:e1900221. DOI: https://doi.org/10.1002/bies.201900221, PMID: 32363691
Biancalani T, Scalia G, Buffoni L, Avasthi R, Lu Z, Sanger A, Tokcan N, Vanderburg CR, Segerstolpe Å, Zhang M, Avraham-Davidi I, Vickovic S, Nitzan M, Ma S, Subramanian A, Lipinski M, Buenrostro J, Brown NB, Fanelli D, Zhuang X, et al. 2021. Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nature Methods 18:1352-1362. DOI: https://doi.org/10.1038/s41592-021-01264-7, PMID: 34711971
Browaeys R, Sang-aram C. 2024. Synthspot. swh:1:rev:72df76f823f7fd60ad91a1afcf55f5bd7dfebe14. Software Heritage. https://archive.softwareheritage.org/swh:1:dir:21f1a0307d3e3e548198f54b91b78db000d205b0; origin=https://github.com/saeyslab/synthspot;visit=swh:1:snp:eec161d2cb85b6a1048298928ef75b41 dbcd62ad;anchor=swh:1:rev:72df76f823f7fd60ad91a1afcf55f5bd7dfebe14
Cable DM, Murray E, Zou LS, Goeva A, Macosko EZ, Chen F, Irizarry RA. 2022. Robust decomposition of cell type mixtures in spatial transcriptomics. Nature Biotechnology 40:517-526. DOI: https://doi.org/10.1038/s41587-021-00830-w, PMID: 33603203
Chen A, Liao S, Cheng M, Ma K, Wu L, Lai Y, Qui X, Yang J, Xu J, Hao S, Wang X, Lu H, Chen X, Liu X, Huang X, Li Z, Hong Y, Jiang Y, Peng J, Liu S, et al. 2022. Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays. Cell 185:1777-1792. DOI: https://doi.org/10.1016/j.cell.2022.04.003, PMID: 35512705
Cho C-S, Xi J, Si Y, Park S-R, Hsu J-E, Kim M, Jun G, Kang HM, Lee JH. 2021. Microscopic examination of spatial transcriptome using Seq-Scope. Cell 184:3559-3572. DOI: https://doi.org/10.1016/j.cell.2021.05.010, PMID: 34115981
Davis J. 2006. The relationship between Precision-Recall and ROC curves. the 23rd international conference. 233-240. DOI: https://doi.org/10.1145/1143844.1143874
Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. 2017. Nextflow enables reproducible computational workflows. Nature Biotechnology 35:316-319. DOI: https://doi.org/10.1038/nbt.3820, PMID: 28398311
Dong R, Yuan GC. 2021. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biology 22:145. DOI: https://doi.org/10.1186/s13059-021-02362-7, PMID: 33971932
Drost HG. 2018. Philentropy: information theory and distance quantification with R. Journal of Open Source Software 3:765. DOI: https://doi.org/10.21105/joss. 00765
Elosua-Bayes M, Nieto P, Mereu E, Gut I, Heyn H. 2021. SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Research 49:e50. DOI: https://doi. org/10.1093/nar/gkab043, PMID: 33544846
Eng C-HL, Lawson M, Zhu Q, Dries R, Koulena N, Takei Y, Yun J, Cronin C, Karp C, Yuan G-C, Cai L. 2019. Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH. Nature 568:235-239. DOI: https://doi. org/10.1038/s41586-019-1049-y, PMID: 30911168
Feng Y, Yang T, Zhu J, Li M, Doyle M, Ozcoban V, Bass GT, Pizzolla A, Cain L, Weng S, Pasam A, Kocovski N, Huang Y-K, Keam SP, Speed TP, Neeson PJ, Pearson RB, Sandhu S, Goode DL, Trigos AS. 2023. Spatial analysis with SPIAT and spaSim to characterize and simulate tissue microenvironments. Nature Communications 14:2697. DOI: https://doi.org/10.1038/s41467-023-37822-0, PMID: 37188662
Guilliams M, Bonnardel J, Haest B, Vanderborght B, Wagner C, Remmerie A, Bujko A, Martens L, Thoné T, Browaeys R, De Ponti FF, Vanneste B, Zwicker C, Svedberg FR, Vanhalewyn T, Gonçalves A, Lippens S, Devriendt B, Cox E, Ferrero G, et al. 2022. Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches. Cell 185:379-396. DOI: https://doi.org/10.1016/j.cell.2021.12.018, PMID: 35021063
Jindal A, Gupta P, Sengupta D. 2018. Discovery of rare cells from voluminous single cell expression data. Nature Communications 9:4719. DOI: https://doi.org/10.1038/s41467-018-07234-6, PMID: 30413715
Karras P, Bordeu I, Pozniak J, Nowosad A, Pazzi C, Van Raemdonck N, Landeloos E, Van Herck Y, Pedri D, Bervoets G, Makhzami S, Khoo JH, Pavie B, Lamote J, Marin-Bejar O, Dewaele M, Liang H, Zhang X, Hua Y, Wouters J, et al. 2022. A cellular hierarchy in melanoma uncouples growth and metastasis. Nature 610:190-198. DOI: https://doi.org/10.1038/s41586-022-05242-7, PMID: 36131018
Kleshchevnikov V, Shmatko A. 2020. Comprehensive Mapping of Tissue Cell Architecture via Integrated Single Cell and Spatial Transcriptomics. bioRxiv. DOI: https://doi.org/10.1101/2020.11.15.378125
Li B, Zhang W, Guo C, Xu H, Li L, Fang M, Hu Y, Zhang X, Yao X, Tang M, Liu K, Zhao X, Lin J, Cheng L, Chen F, Xue T, Qu K. 2022. Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution. Nature Methods 19:662-670. DOI: https://doi.org/10. 1038/s41592-022-01480-9
Li H, Zhang Z, Squires M, Chen X, Zhang X. 2023a. scMultiSim: Simulation of Multi-Modality Single Cell Data Guided by Cell-Cell Interactions and Gene Regulatory Networks. Research Square. DOI: https://doi.org/10. 21203/rs.3.rs-2675530/v1
Li H, Zhou J, Li Z, Chen S, Liao X, Zhang B, Zhang R, Wang Y, Sun S, Gao X. 2023b. A comprehensive benchmarking with practical guidelines for cellular deconvolution of spatial transcriptomics. Nature Communications 14:1548. DOI: https://doi.org/10.1038/s41467-023-37168-7, PMID: 36941264

Longo SK, Guo MG, Ji AL, Khavari PA. 2021. Integrating single-cell and spatial transcriptomics to elucidate intercellular tissue dynamics. Nature Reviews. Genetics 22:627-644. DOI: https://doi.org/10.1038/s41576-021-00370-8, PMID: 34145435
Lopez R, Li B, Keren-Shaul H, Boyeau P, Kedmi M, Pilzer D, Jelinski A, Yofe I, David E, Wagner A, Ergen C, Addadi Y, Golani O, Ronchese F, Jordan MI, Amit I, Yosef N. 2022. DestVI identifies continuums of cell types in spatial transcriptomics data. Nature Biotechnology 40:1360-1369. DOI: https://doi.org/10.1038/s41587-022-01272-8, PMID: 35449415
Lun ATL, Riesenfeld S, Andrews T, Dao TP, Gomes T, participants in the 1st Human Cell Atlas Jamboree, Marioni JC. 2019. EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome Biology 20:63. DOI: https://doi.org/10.1186/s13059-019-1662-y, PMID: 30902100
Maier M. 2014. DirichletReg: Dirichlet Regression for Compositional Data in R. Wirtschaftsuniversität Wien. DOI: https://doi.org/10.57938/ad3142d3-2fcd-4c37-aec6-8e0bd7d077e1
Peng T, Thorn K, Schroeder T, Wang L, Theis FJ, Marr C, Navab N. 2017. A BaSiC tool for background and shading correction of optical microscopy images. Nature Communications 8:14836. DOI: https://doi.org/10.1038/ ncomms14836, PMID: 28594001
Rodriques SG, Stickels RR, Goeva A, Martin CA, Murray E, Vanderburg CR, Welch J, Chen LM, Chen F, Macosko EZ. 2019. Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution. Science 363:1463-1467. DOI: https://doi.org/10.1126/science.aaw1219, PMID: 30923225
Rossi ND, Chen JG. 2022. Analyzing spatial transcriptomics data using giotto. Current Protocols 2:e405. DOI: https://doi.org/10.1002/cpz1.405, PMID: 35384407
Saito T, Rehmsmeier M. 2017. Precrec: fast and accurate precision-recall and ROC curve calculations in R. Bioinformatics 33:145-147. DOI: https://doi.org/10.1093/bioinformatics/btw570, PMID: 27591081
Sang-aramC. 2023. Spotless-benchmark. swh:1:rev:b6c3a9a47f8f6586fb495ff051f5fba35057ca25. Software Heritage. https://archive.softwareheritage.org/swh:1:dir:e35f7ed0519fca5557d49e8cde7a1b50029c146b; origin=https://github.com/saeyslab/spotless-benchmark;visit=swh:1:snp:51e34cc80c889ef7556b962231cbca54 db105410;anchor=swh:1:rev:b6c3a9a47f8f6586fb495ff051f5fba35057ca25
Sato E, Olson SH, Ahn J, Bundy B, Nishikawa H, Qian F, Jungbluth AA, Frosina D, Gnjatic S, Ambrosone C, Kepner J, Odunsi T, Ritter G, Lele S, Chen YT, Ohtani H, Old LJ, Odunsi K. 2005. Intraepithelial CD8+ tumor-infiltrating lymphocytes and a high CD8+/regulatory T cell ratio are associated with favorable prognosis in ovarian cancer. PNAS 102:18538-18543. DOI: https://doi.org/10.1073/pnas.0509182102, PMID: 16344461
SonesonC2018. Towards unified quality verification of synthetic count data with countsimQC. Bioinformatics 34:691-692. DOI: https://doi.org/10.1093/bioinformatics/btx631, PMID: 29028961
Song Q, Su J. 2021. DSTG: deconvoluting spatial transcriptomics data through graph-based artificial intelligence. Briefings in Bioinformatics 22:bbaa414. DOI: https://doi.org/10.1093/bib/bbaa414, PMID: 33480403
Song D, Wang G. 2024. scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics. Nature Biotechnology 42:247-252. DOI: https://doi.org/10.1038/s41587-023-01772-1
Ståhl PL, Salmén F, Vickovic S, Lundmark A, Navarro JF, Magnusson J, Giacomello S, Asp M, Westholm JO, Huss M, Mollbrink A, Linnarsson S, Codeluppi S, Borg Å, Pontén F, Costea PI, Sahlén P, Mulder J, Bergmann O, Lundeberg J, et al. 2016. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353:78-82. DOI: https://doi.org/10.1126/science.aaf2403
Stringer C, Wang T, Michaelos M, Pachitariu M. 2021. Cellpose: A generalist algorithm for cellular segmentation. Bioinformatics 01:931238. DOI: https://doi.org/10.1101/2020.02.02.931238
Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, Hao Y, Stoeckius M, Smibert P, Satija R. 2019. Comprehensive integration of single-cell data. Cell 177:1888-1902. DOI: https://doi.org/10.1016/j.cell.2019. 05.031, PMID: 31178118

Sun D, Liu Z, Li T, Wu Q, Wang C. 2022. STRIDE: accurately decomposing and integrating spatial transcriptomics using single-cell RNA sequencing. Nucleic Acids Research 50:e42. DOI: https://doi.org/10.1093/nar/gkac150, PMID: 35253896
Tasic B, Menon V, Nguyen TN, Kim TK, Jarsky T, Yao Z, Levi B, Gray LT, Sorensen SA, Dolbeare T, Bertagnolli D, Goldy J, Shapovalova N, Parry S, Lee C, Smith K, Bernard A, Madisen L, Sunkin SM, Hawrylycz M, et al. 2016. Adult mouse cortical cell taxonomy revealed by single cell transcriptomics. Nature Neuroscience 19:335-346. DOI: https://doi.org/10.1038/nn.4216, PMID: 26727548
Traag VA, Waltman L, van Eck NJ. 2019. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports 9:5233. DOI: https://doi.org/10.1038/s41598-019-41695-z, PMID: 30914743
Vallania F, Tam A, Lofgren S, Schaffert S, Azad TD, Bongen E, Haynes W, Alsup M, Alonso M, Davis M, Engleman E, Khatri P. 2018. Leveraging heterogeneity across multiple datasets increases cell-mixture deconvolution accuracy and reduces biological and technical biases. Nature Communications 9:4735. DOI: https://doi.org/ $10.1038 / \mathrm{s} 41467-018-07242-6$, PMID: 30413720
van der Walt S, Schönberger JL, Nunez-lglesias J, Boulogne F, Warner JD, Yager N, Gouillart E, Yu T, scikit-image contributors. 2014. scikit-image: image processing in Python. PeerJ 2:e453. DOI: https://doi.org/10.7717/ peerj.453, PMID: 25024921
Wagner A, Regev A, Yosef N. 2016. Revealing the vectors of cellular identity with single-cell genomics. Nature Biotechnology 34:1145-1160. DOI: https://doi.org/10.1038/nbt.3711, PMID: 27824854
Wang X, Allen WE, Wright MA, Sylwestrak EL, Samusik N, Vesuna S, Evans K, Liu C, Ramakrishnan C, Liu J, Nolan GP, Bava FA, Deisseroth K. 2018. Three-dimensional intact-tissue sequencing of single-cell transcriptional states. Science 361:eaat5691. DOI: https://doi.org/10.1126/science.aat5691, PMID: 29930089

Wang X, Park J, Susztak K, Zhang NR, Li M. 2019. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications 10:380. DOI: https://doi.org/10.1038/s41467-018-08023-x, PMID: 30670690
Wolf FA, Angerer P, Theis FJ. 2018. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19:15. DOI: https://doi.org/10.1186/s13059-017-1382-0, PMID: 29409532
Xia C, Fan J, Emanuel G, Hao J, Zhuang X. 2019. Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression. PNAS 116:19490-19499. DOI: https:// doi.org/10.1073/pnas.1912459116, PMID: 31501331
Yan L, Sun X. 2023. Benchmarking and integration of methods for deconvoluting spatial transcriptomic data. Bioinformatics 39:btac805. DOI: https://doi.org/10.1093/bioinformatics/btac805, PMID: 36515467
Yao Z, van Velthoven CTJ, Nguyen TN, Goldy J, Sedeno-Cortes AE, Baftizadeh F, Bertagnolli D, Casper T, Chiang M, Crichton K, Ding SL, Fong O, Garren E, Glandon A, Gouwens NW, Gray J, Graybuck LT, Hawrylycz MJ, Hirschstein D, Kroll M, et al. 2021. A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. Cell 184:3222-3241. DOI: https://doi.org/10.1016/j.cell.2021.04.021, PMID: 34004146
ZhuJ2023. SRTsim: spatial pattern preserving simulations for spatially resolved transcriptomics. Genome Biology 24:39. DOI: https://doi.org/10.1186/s13059-023-02879-z, PMID: 36869394

## Appendix 1

## Description and validation of the simulation procedure

To generate synthetic spots, synthspot by default samples $2-10$ cells from input scRNA-seq data. The counts from these cells are then summed up per gene then downsampled using the downsampleMatrix function from DropletUtils (Lun et al., 2019) to be within the given mean and standard deviation $(20,000 \pm 5000$ counts by default). The uniqueness of synthspot lies in the variable cell type frequency priors between abundance patterns, which determine the probability that a cell type will be sampled during spot generation. The nine abundance patterns used in our benchmark are made up of three characteristics, (1) the uniformity of each region, the (2) distinctness of cell types within each region, and (3) whether or not there are missing, dominant, or rare cell types (Figure 1-figure supplement 1). Appendix 1-figure 1 depicts the simulation process in detail. For simplicity purposes, we have excluded the partially dominant, regionally rare, and missing cell types abundance patterns from this flowchart. For the partially dominant pattern, we would also randomly select a region the dominant cell type would be absent (cell type prior is then set to zero), and another region where it is as equally abundant as other cell types (prior is sampled from uniform distribution of Wagner et al., 2016; Li et al., 2022). Similarly, for the regionally rare pattern, we would select a select where the rare cell type will only be present in, and then the prior of the rare cell type would be set to zero for all other regions. For the missing cell types pattern, we would first randomly select four cell types to be removed. Note that it is also possible for synthspot to use the cell type composition of regionally annotated scRNA-seq data as frequency priors. In that case, it will create the number of regions equal to the scRNA-seq data and use cell type frequencies in the corresponding real region.

As an example, let us say we want to create a synthetic dataset containing two artificial regions which follow the dominant cell type pattern. Consider an input scRNA-seq dataset with eight cell types called $\mathrm{A}, \mathrm{B}, \mathrm{C}, \ldots, \mathrm{H}$. In this case, a dominant cell type is randomly selected which will be present in all regions $(=\mathrm{H})$. Assume that region 1 contains cell types $[\mathrm{A}, \mathrm{B}, \mathrm{C}, \mathrm{E}, \mathrm{H}]$. The abundance of $H$ is sampled from a uniform distribution from 75 to 100 , while the rest will be sampled from 1 to 10 . Suppose we obtain the following priors: $\mathrm{H}=80, \mathrm{~A}=1, \mathrm{~B}=2, \mathrm{C}=3$, and $\mathrm{E}=4$. The summed abundance is 90 , and the cell type frequencies are now: $\mathrm{H}=0.9, \mathrm{~A}=0.01, \mathrm{~B}=0.02, \mathrm{C}=0.03$ and $\mathrm{E}=0.04$. For all spots generated in region 1, these are the probabilities in which the cell types will be sampled.

We validated that our synthetic data and its abundance patterns sufficiently matches real Visium data in two ways, comparing the distributions with countsimQC (Soneson, 2018) and using frequency priors based on real data. We compared synthspot with the algorithms to generate synthetic data used by cell2location, stereoscope, and SPOTlight. We used brain and kidney Visium datasets as the reference, and generated synthetic data using scRNA-seq data from the respective organs. These were the same scRNA-seq datasets used in our silver standard.

The counts per gene of each dataset seem to be representative of each algorithm's performance for other metrics, for example expression distribution, dispersion, mean-variance trend, and fraction of zeros (Appendix 1-figures 2 and 3). For the brain dataset, synthspot has the most resemblance with real data, followed by SPOTlight (Appendix 1-figure 2). Since cell2location and stereoscope did not implement a downsampling step in their simulation, their synthetic brain datasets had overly abundant counts, a result of the plate-based scRNA-seq dataset (SMART-seq). SPOTlight downsamples each spot to have a total UMI count of 20,000 , so the count distribution becomes uniform, unlike real data. For the kidney, all algorithms except stereoscope's were able to generate synthetic data that resembled real data (Appendix 1-figure 3). For most measures, cell2location's algorithm had the most resemblance with real data. Nonetheless, the robustness of synthspot towards sequencing technologies of the input dataset along with the flexibility in finetuning properties of the synthetic data makes it the preferred tool to aid in benchmarking. Finally, we verified that the distributions of the nine synthspot abundance patterns did not differ from one another visually or from the real abundance pattern, which uses real annotations as the frequency priors (Appendix 1-figure 4).

As a second verification, we used the regional annotation in the brain cortex dataset with the real abundance pattern to generate synthetic data with five brain regions (L1, L2/3, L4, L5, and L6), with each region having the same composition as the real layer. We then compared method performance between the real and diverse overlap patterns (Appendix 1-figure 5). Although the
spot compositions between the patterns are different, method performances are similar, validating that artificial patterns can be used to evaluate model performance.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-24.jpg?height=1756&width=1419&top_left_y=296&top_left_x=580)
Appendix 1-figure 1. Schematic of the synthspot simulation algorithm.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-25.jpg?height=1646&width=1441&top_left_y=174&top_left_x=573)

Appendix 1-figure 2. Plots comparing the characteristics of real Visium data from mouse brain and synthetic datasets generated from brain scRNA-seq data using different algorithms. (a) Average abundance values (log counts per million) per gene. (b) Association between average abundance and the dispersion. (c) Distribution of effective library sizes, or the total count per sample multiplied by the corresponding TMM normalization factor calculated by edgeR. (d-e) Distribution of pairwise Spearman correlation coefficients for 500 randomly selected spots (d) and genes (e), calculated from log CPM values. Only non-constant genes are considered. (f) Distribution of the fraction of zeros observed per spot. $(\mathbf{g}-\mathbf{h})$ The association between fraction zeros and average gene abundance (g) and total counts per spot (h).
![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-26.jpg?height=1625&width=1436&top_left_y=174&top_left_x=573)

Appendix 1-figure 3. Plots comparing the characteristics of real Visium data from mouse kidney and synthetic datasets generated from kidney scRNA-seq data using different algorithms. (a) Average abundance values (log counts per million) per gene. (b) Association between average abundance and the dispersion. (c) Distribution of effective library sizes, or the total count per sample multiplied by the corresponding TMM normalization factor calculated by edgeR. (d-e) Distribution of pairwise Spearman correlation coefficients for 500 randomly selected spots (d) and genes (e), calculated from log CPM values. Only non-constant genes are considered. (f) Distribution of the fraction of zeros observed per spot. $(\mathbf{g}-\mathbf{h})$ The association between fraction zeros and average gene abundance (g) and total counts per spot (h).

(a)
![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-27.jpg?height=1615&width=1430&top_left_y=216&top_left_x=577)

Appendix 1-figure 4. Plots comparing the characteristics of real Visium data from mouse brain and the eight synthetic abundance patterns from synthspot generated from brain scRNA-seq data. (a) Association between average abundance and the dispersion. (b) Average abundance values (log counts per million) per gene. (c) Distribution of effective library sizes, or the total count per sample multiplied by the corresponding TMM normalization factor calculated by edgeR. (d-e) Distribution of pairwise Spearman correlation coefficients for 500 randomly selected spots (d) and genes (e), calculated from log CPM values. Only non-constant genes are considered. (f) Distribution of the fraction of zeros observed per spot. $(\mathbf{g}-\mathbf{h})$ The association between fraction zeros and average gene abundance (g) and total counts per spot (h).
![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-28.jpg?height=1328&width=1443&top_left_y=176&top_left_x=571)

Appendix 1-figure 5. Method performance on synthetic data generated from completely synthetic or annotated brain regions. When using an artificial abundance pattern (diverse overlap) to create synthetic spatial data, method rankings remain almost identical as when using a real abundance pattern. The real pattern uses regional annotations from the scRNA-seq input to create regions with the same cell type frequencies.

## Appendix 2

## Issues with threshold-based classification metrics

In addition to the area under the precision-recall curve (AUPR), we initially included the (balanced) accuracy, specificity, recall (sensitivity), precision, and F1 score in the evaluation. These metrics evaluate the classification capabilities of each method, that is the correctness of cell type presence and absence prediction. Briefly, accuracy is the percentage of correctly classified cell types, specificity measures how many of the cell types predicted as absent are truly absent, sensitivity measures how well a method can detect a cell type within a spot, precision measures how many cell types predicted as present are truly present, and the F1 score integrates sensitivity and precision.

The first issue was that methods that use probabilistic models (e.g. cell2location, stereoscope, RCTD and DestVI) do not return proportions that are exactly zero but instead negligible values as low as $10^{-9}$. This made an unbiased evaluation difficult since a fixed threshold for cell type presence/ absence must be selected to calculate classification metrics. In particular, different methods, datasets and abundance patterns have different thresholds for which the classification metrics are at a maximum.

The second issue stems from the class imbalance of our datasets, in which more cell types are absent than present in a spot (more negative than positive classes). In general, around 15\% of the proportion matrix are positives classes, which made the specificity and precision particularly uninformative. This can be seen by how Seurat was the best performer on both specificity and precision despite having low sensitivity (Appendix 2-figure 1). By mostly predicting cell types to be absent in a spot, there are more false negatives (FNs), but specificity and precision do not take FNs into account. The balanced accuracy and F1 score were also unable to entirely correct for this class imbalance.

Given these two issues, we decided to use the precision-recall curve (PR) for evaluation instead and not include these five classification metrics (although we discuss its calculation below). The PR curve plots precision against recall at different thresholds, and the threshold to distinguish cell type absence/presence is varied from zero to one instead of fixing it at a certain value. Hence, the proportions are used as-is without rounding or binarization. It is also recommended for use with imbalanced datasets (Davis, 2006).

## Calculation of threshold-based classification metrics

First, we rounded the predicted proportion matrices to two decimal points, so a proportion of 0.005 and under was rounded to zero. We calculated the micro-average of each classification metric. This is a global metric where the contributions of all classes are considered. We essentially treat the proportion matrix as in a binary classification problem and go through each element individually. As an example, the micro-precision is calculated as

$$
\text { Precision }_{\text {micro }}=\frac{\sum_{z}^{Z} T P_{z}}{\sum_{z}^{Z} T P_{z}+\sum_{z}^{Z} F P_{z}}
$$

where TP stands for true positive and FP for false positive. This is in contrast to the macro-average, where the metric is computed independently for each cell type using a one-vs-all approach, and the average is taken across all cell types (Precision macro $=\frac{1}{Z} \sum_{z}^{Z}$ Precision $_{z}$ ).
![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-30.jpg?height=581&width=995&top_left_y=182&top_left_x=580)

Appendix 2-figure 1. The relative frequency in which a method performs best in the silver standard, based on the best median value across ten replicates for that combination. 'Tie' means that two or more methods score the same up to the third decimal point. RMSE: root-mean-square error; AUPR: area under the precision-recall curve.

## Appendix 3

## Method execution and parameter choice

In this appendix, we briefly describe each method, how we ran them, and the parameters that we evaluated. As most methods contain adjustable parameters that can affect its performance, we tested a range of parameter options to ensure optimal performance for each method. For reproducibility, users can find the exact parameters we have used for each analysis under the 'conf/' folder in our Github repository.

## cell2location (v0.06a)

Cell2location models the transcripts with a negative binomial distribution. It uses variational inference to estimate all the parameters. One advantage of cell2location is that users can shape model priors by providing hyperparameters that correspond to their prior knowledge. We used default priors from the tutorial of cell2location version 0.3. We filtered genes from the reference and spatial datasets as suggested by the tutorial. Model fitting was performed with cross-validation stratified by cell type annotation (stratify_cv). Sample information was not given to the model. The number of training iterations was set to 30,000 (n_iter).

In a newer tutorial, the default value of detection_alpha hyperparameter changed from 200 to 20. We compared the old and new values on the brain cortex seqFISH+ dataset (gold standard) and the kidney dataset (silver standard) and did not find a difference in performance (Appendix 3-figure 1a). Hence, we kept detection_alpha $=200$ for our benchmark. We also varied the number of cells per location from 10 to 50 cells in the gold standard but also did not find a noticeable change in performance (Appendix 3-figure 1b).

## DestVI (v0.16.0)

DestVI uses latent variable models (i.e. variational autoencoders) for the single-cell and spatial data. Unlike other methods, it also models a continuous estimate of cell state for every cell type in every spot. With the default of 2500 epochs, training did not converge for many of the silver standard datasets. We followed the author's recommendations and increased to 5000 training epochs and reduce the minibatch size (batch_size) to 64, which improved performance for all silver standard datasets (Appendix 3-figure 3). For the gold standards, we used 2500 training epochs and batch_ size $=4$.

## DSTG (v0.0.1)

DSTG uses a graph convolutional neural network to learn the composition of real spots from simulated spots. It performs a joint dimensionality reduction (canonical correlation analysis, CCA) and identifies the mutual nearest neighbors between the real and simulated spots. As many of the parameters were hardcoded, for example, 200 nearest neighbors and 30 canonical vectors, the algorithm did not run on our gold standards where there were only nine spots per FOV. We adjusted the source code to change $k \_$filter, $k$, num_cc, and dims to be equal to the number of spots in such cases, otherwise the default parameters were used.

## MuSiC (v0.2.0)

MuSiC is a bulk deconvolution method developed to handle multiple scRNA-seq references from multiple subjects. It employs weighted non-negative least squares (NNLS) regression. In case there are scRNA-seq reference datasets from multiple samples, the between-subject variance is used as weights for each gene. Genes with consistent expression among subjects are considered more informative and will be given higher weights during regression. We ran MuSiC without pre-grouping of cell types. We did not provide subject information to the model but instead considered each cell as a separate subject. This is not all our silver standards contained sample information, and we observed that the performance remained the same or worse when the sample information was given (Appendix 3-figure 2a). We also experimented with creating 'pseudosamples' for datasets without sample information, where each cell was randomly assigned to one out of three artificial samples. This returned a worse performance than when single cells were used as samples (Appendix 3figure 2b).

RCTD (v1.2.0)
RCTD models the transcripts as being Poisson distributed and uses maximum likelihood estimation to infer cell type proportions. We ran RCTD with doublet_mode="full", indicating that many cell types per spot were to be expected.

## Seurat integration (v4.1.0)

Using joint dimensionality reduction between the scRNA-seq (reference) and spatial (query) data, we defined compatible reference-query pairs and used them as anchors for the label transfer procedure. We obtain a probability for a cell type being in the spot and use this as proxy for the abundance. We compared two normalization methods (SCTransform and vst) and two projection methods (PCA and CCA) and found that using SCTransform with PCA gave the best results (Appendix 3-figure 4). As with DSTG, we changed the number of neighbors and vectors ( $k$. score, $k$.weight, dims) used in the gold standard as equal to the number of spots, otherwise the default was used.

## SpatialDWLS (v1.1.0)

SpatialIDWLS performs cell-type enrichment analysis for each spot, then uses down-weighted least squares on marker genes. We followed the protocol described in Rossi and Chen, 2022, where the makeSignMatrixDWLS function was used with top 100 marker genes, instead of following the online vignette where the makeSignMatrixPAGE function was used with all marker genes (Appendix 3figure 5).

## SPOTlight (v0.1.7)

SPOTlight is the only method based on non-negative matrix factorization (NMF) and NNLS. This method computes topics from the gene expression profile with NMF instead of using the expression values directly. NNLS is used with the topics to obtain the cell type proportions for each spot. We followed the vignette of version 0.1.5, as the recommended parameters have changed slightly between each version of SPOTlight. In this version, SPOTlight uses Seurat's FindAllMarkers functions to calculate marker genes between cell types, and the logfc.threshold (limit testing to genes which show at least X-fold difference) and min.pct (only test genes that are expressed in at least X\% of cells) parameters have a huge impact on the resulting list of marker genes. Furthermore, the parameters cl_n (number of cell types to use) and min_cont (only keep cells that are at least X\% present in a spot) within the deconvolution function itself affects the predictions. We tested three sets of parameters and in the end, ran FindAllMarkers with only.pos $=$ TRUE, logfc.threshold $=1$, and min.pct $=0.9$, and the deconvolution with cl_n=50 and min_cont $=0.09$ (Appendix 3-figure 6). This was also the parameter set that gave the shortest runtime. We also normalized both the reference and spatial data with SCTransform.

## Stereoscope (v0.2.0)

Stereoscope models the transcripts with a negative binomial distribution. It uses maximum likelihood estimation to infer the rate and overdispersion parameters from the scRNA-seq data and then uses maximum a posterior (MAP) estimation to infer cell type proportions. We ran stereoscope using all genes for the silver standard but only the top 5000 most highly expressed genes for the gold standard, as that gave the best results for both cases. Although the authors have noted in their paper that choosing the 5000 most expressed genes is sufficient, we saw that using all genes for the silver standards still gave slightly better performance (Appendix 3-figure 7a). We implemented an option to use the HVGs instead of top expressed genes, but this did not consistently result in better performance (Appendix 3-figure 7b). Finally, we tested the sub parameter, which subsamples each cell type in the single-cell reference to at most X cells, but did not see any improvement in the kidney (silver standard) or liver dataset (Appendix 3-figure 7c). We verified that both the training and test models have converged.

## STRIDE (v0.0.2)

STRIDE trains a topic model from scRNA-seq data, then applies this model to the spatial data. We ran it with --normalize and we let STRIDE automatically select the optimal topic number for each dataset. We did not find a pattern when comparing results from raw and normalized counts (Appendix 3-figure 8).

## Tangram (v1.0.3)

Tangram learns a spatial alignment of scRNA-seq from the spatial data via nonconvex optimization. Like Seurat integration, Tangram returns the probabilistic counts for each cell type in each cell voxel. We tested the three mapping modes: cells (maps single cells to spots), clusters (maps average of cells per cell type instead of single cells), and constrained (constrain the number of mapped single cell profiles). Although the constrained mode was recommended for deconvolution, we found that using the clusters resulted in the best performance (Appendix 3-figure 9). When running the constrained mode, we also provided the ground truth number of cells per spot and total number of cells in the density_prior and target_counts parameters, respectively. Increasing the training epochs did not have an effect on the performance, and we verified that the models have converged. The final parameters we used were map_cells_to_space with mode="clusters" and density_prior="rna_ count_based". We used the top 100 marker genes for each cell type.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-33.jpg?height=652&width=1432&top_left_y=683&top_left_x=575)
Appendix 3-figure 1. Changing hyperparameters in the cell2location model. There is almost no performance difference when changing the (a) detection alpha and (b) number of cells per spot.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-33.jpg?height=744&width=1424&top_left_y=1497&top_left_x=577)
Appendix 3-figure 2. Comparing MuSiC performance when the model was given sample information. MuSiC seems to perform best when single cells were used as samples ('None'), as compared to (a) when the real sample information was given and (b) when pseudosamples were created.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-34.jpg?height=845&width=1423&top_left_y=182&top_left_x=580)

Appendix 3-figure 3. Comparing parameters of DestVI. Compared to the default parameters ( 2500 epochs), DestVI has better performance with 5000 training epochs and batch size of 64 .
![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-34.jpg?height=845&width=1428&top_left_y=1192&top_left_x=577)

Appendix 3-figure 4. Comparing different data transformation and dimensionality reduction methods in Seurat. Seurat has the best performance when the data was transformed using SCTransform and the dimensionality reduction method is PCA. CCA = canonical correlation analysis; PCA = principal component analysis; VST = variance stabilizing transformation.
![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-35.jpg?height=843&width=1428&top_left_y=182&top_left_x=577)

Appendix 3-figure 5. Comparing SpatialDWLS between two signature matrix creation functions.
SpatialDWLS has better performance when the makeSignMatrixDWLS was used to create the signature matrix (as described in the Current Protocols paper), instead of the makeSignMatrixPAGE function (described in the online vignette).

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-35.jpg?height=1108&width=1426&top_left_y=1286&top_left_x=577)
Appendix 3-figure 6. Comparing SPOTlight performance on three sets of parameters. We used Set1 parameters in our benchmark.

![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-36.jpg?height=463&width=1423&top_left_y=182&top_left_x=580)

Appendix 3-figure 7. Comparing stereoscope performance to using only highly variable genes (HVGs). There is no consistent performance difference between using all genes of stereoscope and using the 5000 HVGs (with or without subsampling the scRNA-seq reference).
![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-36.jpg?height=216&width=1421&top_left_y=848&top_left_x=580)

Appendix 3-figure 8. Comparing STRIDE performance on normalized and raw counts. For STRIDE, there is no consistent performance difference between normalizing or using the raw counts.
![](https://cdn.mathpix.com/cropped/2025_11_28_9d7879124f7d67f07ab9g-36.jpg?height=843&width=1423&top_left_y=1233&top_left_x=580)

Appendix 3-figure 9. Comparing the mapping modes in Tangram. Although the constrained mapping mode was recommended in the Tangram vignette, we found that the clusters mode achieve better performance.


