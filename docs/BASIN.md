# BASIN: Bayesian mAtrix variate normal model with Spatial and sparsIty priors in Non-negative deconvolution 

Jiasen Zhang ${ }^{1}$, Xi Qiao ${ }^{2}$, Liangliang Zhang ${ }^{2, ~}{ }^{*}$, and Weihong Guo ${ }^{1, ~}{ }^{*}$<br>${ }^{1}$ Department of Mathematics, Applied Mathematics and Statistics, Case Western Reserve University, Cleveland, OH<br>${ }^{2}$ Department of Population and Quantitative Health Sciences, Case Western Reserve University, Cleveland, OH<br>*Corresponding author: lxz716@case.edu, weihong.guo@case.edu


#### Abstract

Spatial transcriptomics allows researchers to visualize and analyze gene expression within the precise location of tissues or cells. It provides spatially resolved gene expression data but often lacks cellular resolution, necessitating cell type deconvolution to infer cellular composition at each spatial location. In this paper we propose BASIN for cell type deconvolution, which models deconvolution as a nonnegative matrix factorization (NMF) problem incorporating graph Laplacian prior. Rather than find a deterministic optima like other recent methods, we propose a matrix variate Bayesian NMF method with nonnegativity and sparsity priors, in which the variables are maintained in their matrix form to derive a more efficient matrix normal posterior. BASIN employs a Gibbs sampler to approximate the posterior distribution of cell type proportions and other parameters, offering a distribution of possible solutions, enhancing robustness and providing inherent uncertainty quantification. The performance of BASIN is evaluated on different spatial transcriptomics datasets and outperforms other deconvolution methods in terms of accuracy and efficiency. The results also show the effect of the incorporated priors and reflect a truncated matrix normal distribution as we expect.


Keywords: spatial transcriptomics, deconvolution, Bayesian NMF, matrix normal distribution

## 1 Introduction

Spatial transcriptomics (ST) is a cutting-edge technique in the field of genomics that allows researchers to analyze the gene expression profile of tissues while preserving spatial information [1]. It offers insights into the spatial distribution of gene expression within individual cells or regions of a tissue, and it has many exciting applications in fields such as developmental biology, cancer research, neurobiology, nephrology, and immunology 2. 6. ST technology can be performed at either single-cell or spot-level resolution. Obviously, single-cell resolution ST provides higher granularity for downstream analyses and facilitates more accurate biological interpretations, but it is more expensive. Spot resolution ST technology, such as 10x Visium, is cheaper and more prevalent. Typically, 10x Visium data includes three primary types of information: gene expression data, spatial information encoded by barcodes, and histological images, which enables simultaneous integration of molecular and imaging data generated from Formalin-Fixed Paraffin Embedded (FFPE) or frozen tissue sections. However, a main limitation of spot resolution ST data is its inability to precisely map molecular profiles to individual cells, as each spot captures the transcriptomes of roughly 10-200 cells depending on the tissue context [7. To address this, inferring the cell-type composition at each spatial location typically requires a process known as cell type deconvolution 8 or spatial decomposition 9 . It involves employing computational and statistical models to decipher the intricate patterns of gene expression within tissues and infer specific cell types accurately.

To achieve single-cell level recognition within each spatial spot, researchers typically integrate paired dropletbased single-cell RNA-seq (scRNA-seq) profiles into their studies. It is crucial for unraveling the intricacies inherent in ST as scRNA-seq data provides a more detailed and comprehensive view of individual cell profiles and thus bridges the gap between spatial information and detailed cellular characterization. This enables the identification of cell types and facilitates the learning of the distribution of signals at a granular level. With
different ways to represent and analyze ST data, recent deconvolution methods can be broadly grouped as follows: nonnegative matrix factorization (NMF)-based methods such as SPOTlight 10 and NMFreg 11, regression and Bayesian model-based methods such as CARD [12, SpatialDWLS [13, MuSiC [14, RCTD [15], Cell2location [16, DestVI 17, SpatialDecon [18, Stereoscope [19] and STRIDE [20, optimal transport (OT)based methods such as SpaOTsc 21 and novoSpaRc 22, graph convolutional network (GCN)-based methods such as DSTG 23 and SD ${ }^{2}$ 24, and deep learning-based methods such as TANGRAM 25 and SpaDecon 26 . There are also reference-free methods that do not require scRNA-seq data as a reference, but try to simultaneously find the cell type composition and gene expression profiles of each cell type [12] 27-32]. Although there is no need for external data, reference-free methods are naturally less accurate and can produce ambiguous or biologically uninterpretable components.

In this paper, we propose BASIN for spatial transcriptomics cell type deconvolution using scRNA-seq data as a reference. Our contributions and novelties can be summarized as follows:

1. First, we represent cell type deconvolution in NMF model and incorporate graph Laplacian regularity utilizing both spatial information and histology images that are often collected together with ST data. It has been shown that in spatial transcriptomics, similar cell types have been shown to co-localize, with the co-localization pattern decaying as distance increases 12 . However, most of the methods mentioned above did not consider this important spatial information. In addition, there is also no thorough study on how to better use the histological intensity images that come with ST for better performance. High-quality histology images clearly show different regions of the tissue and provide informative clustering guidance.
2. Second, unlike other NMF-based and probabilistic model-based methods that find deterministic solutions based on methods such as MAP estimation and gradient descent, we adopt a Bayesian framework and derive a Gibbs sampler to approximate the posterior distribution of the parameters and variables involved. As a result, BASIN adds an interpretable randomness to the model and helps enhance robustness to local optima and enable quantification of uncertainty.
3. Third, none of the methods mentioned above incorporates sparsity in their models. As a result, they may return too uniformly distributed proportions and fail to distinguish different regions on the tissue. As a solution, we exploit Laplacian distribution as a prior to induce sparsity in our Bayesian model. When replaced with an exponential prior, it retains the sparsity-inducing property while inherently enforcing non-negativity.
4. In computational aspects, we propose a matrix variate Bayesian NMF method to derive the posteriors of the parameters and variables. In previous Bayesian NMF methods 33, 34, the matrices are vectorized to derive a multivariate normal posterior distribution element-wisely. However, we keep all the variables in their matrix form and derive a matrix normal posterior distribution, which enables the combination of Bayesian model and graph Laplacian prior, and improve sampling efficiency.
We evaluate the performance of BASIN on different simulation studies and several published real spatial transcriptomics datasets, and compared it with some other recently proposed methods with high citations. We also illustrate the correctness of the derived posterior distributions and how it benefits the results.

## 2 Results

### 2.1 Simulation study

The workflow of BASIN is summarized in Fig. 1 The details of BASIN, and data processing are discussed in Methods. To evaluate the performance and capability of BASIN, we perform four numerical experiments using simulated data and compare our results with five deconvolution methods: Stereoscope, SPOTlight, SpatialDWLS, RCTD and CARD. In summary, the simulation is based on spatial transcriptomics and scRNA sequence data from a pancreatic ductal adenocarcinoma tumor sample (PDAC-A) 35, with three dominant cell types: acinar cells, cancer cluster 1 , and terminal ductal cells, each predominantly localized in a distinct region. In the different cases, we assume that the cell type proportions are piecewise constant (Simulation 1), piecewise smooth (Simulation 2 and 3), and Dirichlet distributed (Simulation 4). We evaluated the results in three metrics: root mean squared error (RMSE), structural similarity index measure (SSIM), and Jensen-Shannon distance (JSD). The details of the simulation strategy, the compared methods, and the evaluation metrics are described in Methods.

The results shown in Fig. 2 indicate that the proposed method consistently outperforms other methods in all simulations and metrics. In simulation 1, we set the cell type proportions to be piecewise constant to accurately evaluate the fitting ability. In this simple case, the proportions of all the locations consist of $0.7,0.15$ and 0.15 with clear but complicated boundaries among different regions (Supplementary Fig. 1). BASIN achieves the

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-03.jpg?height=717&width=1569&top_left_y=244&top_left_x=277)
Figure 1: Workflow of BASIN. BASIN starts from a spatial transcriptomics data (the first row), which is the target data to be deconvolved, and a scRNA-seq data (the third row), which is a reference to extract cell type information. To utilize the spatial structure in the ST data, we calculate the Graph Laplacian based on the coordinate information and the histology image (the second row). Besides, the mean gene expression for each cell type is computed from the scRNA-seq data as another input, which largely improve the computational efficiency and reduce memory usage. The outputs (right) are sampled from the posterior distribution, indicating the cell type proportions of each location.

lowest RMSE (0.028), compared with Stereoscope (0.119), SPOTlight (0.176), SpatialDWLS (0.156), RCTD (0.189) and CARD (0.06). BASIN also performs the best in terms of SSIM (0.967) compared with Stereoscope (0.724), SPOTlight (0.627), SpatialDWLS (0.614), RCTD (0.681) and CARD (0.805). Besides, BASIN has much lower JSD (0.034) than Stereoscope (0.118), SPOTlight (0.158), SpatialDWLS (0.257), RCTD (0.249) and CARD (0.071).

For the other three simulations, we add more conditions to simulate different properties that could happen in real data. In simulation 2, we simulate the spatial correlation and the cell type proportions are no longer piecewise constant but change smoothly. In simulation 3, we add noise to the cell type proportions based on simulation 2 to evaluate robustness and the ability of recovering more complicated spatial structure. In simulation 4, we implement a new simulation assuming that the cell type proportions follow Dirichlet distributions to simulate data with less obvious spatial structures. For these simulations, all the methods have similar performance as simulation 1 but show different adaptations. For instance, comparing simulation 2 with simulation 1 , BASIN achieves $0.7 \%$ higher SSIM while RMSE and JSD are almost unchanged. SPOTlight and CARD gain obvious improvement in all of the three metrics, and SpatialDWLS performs worse on the three metrics. Overall, BASIN still performs better than other methods for all the simulations. On the other hand, Stereoscope and SPOTlight have relatively larger errors in acinar and cancer regions (Supplementary Fig. 3). SpatialDWLS has relatively larger errors in ductal region, and RCTD has relatively larger errors in acinar and ductal regions (Supplementary Fig. 32$).$

By comparing the outputs with the ground truth (Supplementary Fig. 2), one can find that SpatialDWLS and RCTD tend to make each spot dominated with single cell type, which is also reflected in the results of the real data later. Although such a property can be helpful to distinguish different regions, it generates unrealistic proportions when there are multiple cell types in one spot. In contrast, Stereoscope and SPOTlight tend to generate more uniformly distributed results, making it difficult to distinguish different regions when the spatial structure is complicated (for example, the results of simulation 4 in Supplementary Fig. 2). CARD and BASIN tend to generate proportions that match the best data, while BASIN has fewer errors overall (Supplementary Fig. 3$).$

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-04.jpg?height=824&width=1571&top_left_y=244&top_left_x=277)
Figure 2: Performance of different methods on the simulated data. Simulation 1: Each region has a piecewise constant cell type composition. Simulation 2: The simulated cell type proportion data is smoothed based on simulation 1. Simulation 3: The cell type proportions follow Gaussian distributions smoothed with a Gaussian kernel. Simulation 4: The three cell types follow different Dirichlet distributions smoothed by a Gaussian kernel. RMSE: Root mean square error, the the lower the better. SSIM: Structural similarity index measure, the higher the better. JSD: Jensen-Shannon distance, the lower the better.

### 2.2 Application in human pancreatic ductal adenocarcinoma data A (PDACA)

We also evaluated BASIN with real-world spatial transcriptomics and compared with other methods. The accessing and processing of these public world data are described in Method. We first evaluate BASIN with human pancreatic ductal adenocarcinoma (PDAC) data 35 and utilize scRNA-seq data with twenty cell types generated from the same resource. We tested both tumor cryosections PDAC-A and PDAC-B. The H\&E staining image of PDAC-A in Fig. 3a can be roughly divided into four regions: cancer cells (red), duct epithelium (yellow), pancreatic tissue (blue) and stroma (the rest area). There is a clear difference in the proportions of cell types among the cander region, duct epithelium and pancreatic tissue while the stroma is similar to the duct epithelium 35 .

In Fig. 3b we show the cell type proportions of each spot and compare with other five deconvolution methods. Among the six methods, SPOTlight fails to distinguish the cancer regions, while the other five methods clearly show the cancer region dominant with two kinds of cancer cells. According to the scRNA-seq data from the same resource, both the two cancer clusters are abundant in the cancer region, which is reflected in the results of Stereoscope, RCTD and BASIN but not in SpatialDWLS and CARD. Besides, only BASIN shows the existence of ductal terminal cells in the cancer region 35 . For non-cancer regions, SpatialSWLS, RCTD, CARD and BASIN correctly distinguish the stroma and pancreatic tissue. Only CARD and BASIN reflect that centroacinar cells and terminal duct cells exist at the junction between peripheral acinar cells and the adjacent ductal epithelium 36 . However, only BASIN predicts a significant amount of antigen-presenting ductal cells expressing MHC Class II near the duct epithelium, matching the scRNA-seq data. To further confirm the accuracy of BASIN, in Fig. 3: we compare the distribution of the six cell types having high enrichment in certain regions with corresponding marker genes 35. Although the distribution of cells can be affected by many critical genes rather than a single one, the marker genes can be a reasonable reference to show the accuracy if they match well with the cells.

In Fig. 3d we show the correlations of cell type proportions across the spatial locations as a hierarchicallyclustered heatmap. It shows the spatial correlations for some pairs of cell types inferred by BASIN. Specifically, the two cancer clusters are distinguished from other cells spatially. Acinar and endocrine cells have higher spatial correlation because they have high enrichment only in the pancreatic tissue. On the other hand, cancer and

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-05.jpg?height=1428&width=1556&top_left_y=457&top_left_x=279)
Figure 3: Application in PDAC-A data. a, The histological image 35 with annotations of four regions: cancer (red), pancreatic (blue), ductal (yellow) regions and interstitium (the rest). b, The scatter pie plot of the outputs showing the proportions of the twenty cell types at each spot. We compare BASIN with other five deconvolution methods. c, The first row shows the proportions of six cell types at each spot. The second row shows the amount of the corresponding marker genes at each spot. d, The clustered correlations between each pair of the cell type proportions across all the spots. e, The standard deviations of the proportions of the six cell types at each spot, calculated with 2000 samples. $\mathbf{f}$, The correlations between the predicted cell type proportions and those of the corresponding scRNA-seq data. $\mathbf{g}$, The histograms of the cell type proportions at one spot selected in the cancer region, calculated with 2000 samples. $\mathbf{h}$, The histograms are plotted together.

centroacinar cells show a relatively low spatial correlation, matching the fact that centroacinar cells only have high enrichment in the stroma and duct epithelium 35. In Fig. 3p we show the derived standard deviations of the proportions of the six cell types, calculated with 2000 posterior draws. It indicates how confident we can be in the estimates and the range of plausible values. For instance, the standard deviations of the proportions are lower in the cancer region and pancreatic tissue, but higher in the central region where different regions are mixed (considering that the proportions are ranged from 0 to 1 ). That means the central region is more difficult to distinguish, as shown in the histology image and the results in Fig. 3 p.

Since the scRNA-seq data are obtained in the same tissue as the spatial transcriptomics data, in Fig. 3f we evaluate the Pearson product-moment correlation between the overall mean cell type proportions of the compared methods and that in the scRNA-seq data, showing that stereoscope, SpatialDWLS, CARD and BASIN obtain more accurate overall mean cell type proportions while BASIN performs the best. In Fig. 3g and h, we select a spot from the cancer region and plot histograms of the estimated cell type proportions based on 2,000 posterior samples. All twenty cell types follow truncated normal distributions with varying means and standard deviations. As expected, the two cancer clones are predominant at this spot, further confirming the accuracy of our method.

### 2.3 Application in human pancreatic ductal adenocarcinoma data B (PDACB)

In this section we evaluate the performance of our method on PDAC-B that is from a different patient than of PDAC-A data. The H\&E staining image of PDAC-B in Fig. 4a shows annotations of three regions: cancer region (red), duct epithelium (yellow) and interstitium (green), corresponding to thirteen cell types in its scRNA-seq data. In this data, each region has one or several dominant cell types 35: cancer cluster A (red), ductal cells (yellow) and endothelial cells (green). Because these cell types are also the most abundant ones in the scRNA-seq data from the same resource, we mainly consider them in this section.

We compare the predicted proportions of BASIN with other five deconvolution methods in Fig. 4b, showing that our proposed method can best match each region with its correct dominant cell types. Among the six methods, Stereoscope and SPOTlight fail to distinguish the cancer and ductal cells. CARD generates excessively high proportion of RBCs (red blood cells) and fails to reflect the region of duct epithelium clearly. RCTD identifies the cancer region and duct epithelium, but the proportion of endothelial cells in the duct epithelium is too high. SpatialDWLS and our method can distinguish the three regions correctly, while BASIN can also distinguish ductal cells expressing major histocompatibility complex (MHC) class II from ductal centroacinar cells within the duct epithelium region. In Fig. 4c we show that the distributions of the four cell types match well with their highly expressed genes as well as the histology image. Specifically, the distribution of ductal cells with MHC class II corresponds to that of CD74, a MHC class II gene 35. In Fig. 4d, we show the correlations of cell type proportions across the spatial locations as a hierarchically-clustered heatmap. The endothelial, ductal centroacinar, cancer cells and ductal MHC class II are divided into different clusters according to the heatmap. In Fig. 4 we show the derived standard deviations of the proportions of these four cell types. It shows that the top and lower-right parts of the tissue have relatively higher uncertainty, where different regions are mixed together. In the histology image Fig. 4a, it corresponds to the left and bottom boundary of the cancer region.

Similar to PDAC-A, in Fig. 4f, we evaluate the correlation between the overall mean cell type proportions of the compared methods and that in the scRNA-seq data. Since all the compared methods fail to distinguish the terminal ductal cells, the cell type correlations are lower than 0.3 overall. But it still show that Stereoscope and BASIN obtain the most accurate overall mean cell type proportions. In Fig. 4 g and Fig. 4 h , we select a spot in the cancer region and plot the histograms of 2000 samples separately and together. All the thirteen cell types follow the truncated normal distributions, and the cancer cell is the most abundant at this spot.

### 2.4 Application in mouse olfactory bulb (MOB) data

We evaluate the replicate 12 of mouse olfactory bulb (MOB) data from 37 and utilize the scRNA-seq data of GSE121891 38 as the reference. The H\&E staining image in Fig. 5k shows an obvious layered structure dividing the tissue into four regions: granule cell layer, mitral cell and outer plexiform layer, glomerular layer, and olfactory nerve layer. And there are five groups of cells in the scRNA-seq data: external plexiform layer interneurons (EPL-IN), granule cells (GC), mitral/tufted cells (M/TC), olfactory sensory neurons (OSs) and periglomerular cells (PGC). Except EPL-IN, the other four groups of cells in the scRNA-seq data correspond to the four layers in the ST data and the deconvolution results should reflect such structures too.

We compare BASIN with other five deconvolution methods in Fig. 5b, which shows that our method predicts the most accurate cell-type proportions. In comparison, Stereoscope and SPOTlight predict the correct structure and but the effect of spatial correlation is too strong and the cell type proportions are inaccurate. For SpatialDWLS and RCTD also show the layers clearly but fail to reflect the complete PGC layer. Besides,

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-07.jpg?height=1544&width=1565&top_left_y=397&top_left_x=279)
Figure 4: Application in PDAC-B data. a, The histological image 35 with annotations of three regions: cancer (red), ductal (yellow) regions and interstitium (green). b, The scatter pie plot of the outputs showing the proportions of the thirteen cell types at each spot. We compare BASIN with other five deconvolution methods. c, The first row shows the proportions of three cell types at each spot. The second row shows the amount of the corresponding marker genes at each spot. d, The clustered correlations between each pair of the cell type proportions across all the spots. $\mathbf{e}$, The standard deviations of the proportions of the three cell types at each spot, calculated with 2000 samples. $\mathbf{f}$, The correlations between the predicted cell type proportions and those of the corresponding scRNA-seq data. $\mathbf{g}$, The histograms of the cell type proportions at one spot selected in the cancer region, calculated with 2000 samples. $\mathbf{h}$, The histograms are plotted together.

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-08.jpg?height=942&width=1578&top_left_y=242&top_left_x=270)
Figure 5: Application in MOB data. a, The histological image from 37 clearly shows a spatial structure with different layers. $\mathbf{b}$, The scatter pie plot of the outputs showing the proportions of the five cell types at each spot. We compare BASIN with other five deconvolution methods. c, The first row shows the proportions of four cell types at each spot. The second row shows the amount of the corresponding marker genes at each spot. d, The standard deviations of the proportions of the five cell types at each spot, calculated with 2000 samples. e, The clustered correlations between each pair of the cell type proportions across all the spots. $\mathbf{f}$, The histograms of the cell type proportions at one spot selected in the granular cell layer, calculated with 2000 samples. g, The histograms are plotted together.

SpatialDWLS shows blurred boundaries of GC, M/TC and OSs, and the overall proportion of MT/C in RCTD is too dominant. The result of CARD and BASIN are overall the best, both of which identify the M/TC and PGC layers accurately, while CARD shows more PGCs in the GC layer. Since EPL-in is trivial in the results, we compare the distribution of the other four cell types with corresponding marker genes in Fig. 5: and show that they match well with each other. Both the cell type proportions and marker genes reflect the correct structures for each layer and show obvious boundaries against each other.

In Fig. 5d we show the standard deviations of each cell type at each location to reflect the uncertainty of our results. Obviously, there are four spots with much higher standard deviations, indicating that they are difficult to deconvolve. For example, the upper right spot is predicted as $100 \%$ EPL-IN by SpatialDWLS, and the middle one is predicted as $100 \%$ EPL-IN by CARD. Besides, the spot at bottom is removed by SPOTlight, SpatialDWLS, RCTD and CARD because none of the cell types match the gene profile at this spot. Therefore, our method can help find such highly ambiguous spots, which may be due to technical artifacts or complex cellular interactions. Fig. 5e reflects that GC, OSs and M/TC layers are identified explicitly. M/TC and PGC have higher correlation because these two layers are adjacent and harder to distinguish. In Fig. 55 and Fig. 5b, we select a spot in the granule cell layer and plot the histograms of 2000 posterior samples separately and together, showing that the GC (granule cell) is predominant at this spot. The EPL-IN, GC, M/TC and PGC are distributed as normal distributions obviously. The OSs is shown as a truncated normal distribution because this cell type is trivial in this region. These results confirm the accuracy of our method and prove that BASIN can output cell type proportions in terms of truncated normal distributions. Some look normal distributed because the tail probability is very close to zero, therefore we didn't draw any sample from there.

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-09.jpg?height=1112&width=1550&top_left_y=242&top_left_x=285)
Figure 6: Application in seqFISH+ data. a, Upper: The spatial locations of the cells are from five field of views (FOVs) on the mouse cortex [39. Lower: The five fields are stitched together with different x and y -offset values for each FOV. The five FOVs are represented with different colors. b, The scatter plot of the outputs showing the dominant cell type of each spot. c, The first row shows the proportions of four excitatory clusters at each spot. The second row shows the amount of the corresponding marker genes at each spot. $\mathbf{d}$, The clustered correlations between each pair of the cell type proportions across all the spots. $\mathbf{e}$, The standard deviations of the proportions of six cell types at each spot, calculated with 2000 samples. $\mathbf{f}$, The probability histograms of the cell type proportions at one spot selected in the granular cell layer, calculated with 2000 samples. Since the standard deviations are too small, they are not plotted together as other data.

### 2.5 Application in seqFISH+ mouse cortex data

We evaluated the performance in seqFISH + cortex data 39 and used GSE102827 scRNA-seq data 40 as a reference. The data is generated from five field of views (FOVs) of the cortex, as shown in Fig. 6a, and the coordinates of the spots within each FOV are adjusted so that the five FOVs are assembled together. The FOVs correspond to different excitatory neuron layers (ExcLs) with different gene expressions. From the left to right, the first three FOVs correspond to ExcL23, ExcL4 and ExcL5, and the last two FOVs are dominant with ExcL6.

BASIN clearly shows the layered structures of the cortex in Fig. 6 p and c which match with the layers in Fig. 6a. In Fig. 6d, one can also observe that the ExcLs predicted by BASIN are located in different clusters with low spatial correlations with each other. By contrast, Sterescope does not show the layers clearly. SPOTlight, SpatialDWLS, RCTD and CARD identify the layers of ExcL23 and ExcL4 well but fail to distinguish the layers of ExcL5 and ExcL6. In our results, the predicted ExcL5 has high proportions in general,l, which is also reflected in Fig. 6e where the ExcL5 proportions have higher standard deviations than the others. All these results show that ExcL5 is relatively harder to distinguish. Finally, Fig. 6f shows the probability histogram of the proportions of cells of a spot selected in the ExcL23 layer, and ExcL23 cells take 63. 5\% on average for the 2000 posterior samples. The standard deviations are close to zero, indicating low gene similarities among cell types. One should note that the proportions of ExcL4, ExcL6, and other excitatory cells appear as small non-zero values due to the smoothing effect of the histogram's bandwidth setting.

## 3 Discussion

In this paper, we have presented BASIN, a cellular deconvolution method for spatial transcriptomics that integrates gene expression information from scRNA-seq data, spatial locations, and histology images to estimate cell-type proportions at each location. Compared with existing approaches, BASIN demonstrates superior performance in multiple simulation settings and real spatial transcriptomics datasets generated from diverse tissues and technologies, highlighting its robustness and ability to preserve spatial structure. In the PDAC data, BASIN shows clear boundaries of subregions on the tissue and accurately distinguishes different cell types within each region. In the MOB data, it clearly shows the layer structures. In the seqFISH+ mouse cortex data, it also predicts different excitatory layers correctly. Despite involving sampling from high-dimensional probability distributions, BASIN is highly efficient (Supplementary Fig. 5) and outperforms most existing methods in computation time.

What makes BASIN unique is that it generates results from a probability distribution, which enables easy parameter tuning and uncertainty quantification. We propose a matrix variate Bayesian method for the NMF problem and combine with graph Laplacian and sparsity priors. Results in different data show that the proportions at each location follow truncated normal distribution as we assume with different uncertainties. We report the standard deviations of the estimated cell type proportions at each spot to quantify the uncertainty arising from gene expression similarity among cell types in the scRNA-seq data, the spatial structure in the transcriptomics data, and the selection of hyperparameters.

Cellular deconvolution methods like BASIN still have room to be improved. First, the histology images that contain spatial information could be better utilized. The histology images are typically not aligned with the ST and have much higher spatial resolutions than ST. In our method, we align the spot coordinates in ST data with the histology image manually and resize the image to a smaller size like $1024 \times 1024$, and we compute the average intensity of a small square area ( $5 \times 5$ ) around each spot location in the ST data. However, better image alignment techniques could be studied to align the two sets of data accurately and further increase the accuracy. Second, sampling from truncated multivariate normal distributions in our method may still fail when the sampling dimension is too high. The issues arise either from the sampling algorithm itself which can be addressed through resampling or from inappropriate hyperparameter choices that render the covariance matrices to be non-invertible. Although the current sampling method 41 is the best for BASIN considering both efficiency and accuracy, better substitute can be used in the future.

BASIN can be extended in several ways in future work. First, while BASIN incorporates gene similarity from both ST and scRNA-seq data, it currently treats all genes equally. However, gene importance varies: some genes are uniformly expressed across the tissue, while others may be uninformative or even misleading for distinguishing cell types. The gene expressions can be weighted by their contribution of differentiating cell types and a weight matrix can be added to the NMF-based models. Second, BASIN is a reference-based method using scRNA-seq data and can also be extended to a reference-free method. We provide a matrix variate method for Bayesian NMF but one factorized matrix is known from scRNA-seq data. If scRNA-seq reference is not available, one needs to derive the posteriors for both of the factorized matrices. In the future, we should think of how to use ST data to help cluster the cell-types across spatical locations 17,32. Third, BASIN can be generalized to tensor methods by involving tensor normal distributions [42, enabling us to deal with spatial dimensions separately and extract meaningful latent structures from high dimensional ST data 43. Combination with tensor decomposition techniques may be needed to reduce computational consumption $44-46$.

## 4 Methods

### 4.1 The BASIN method

The goal of cell type deconvolution is to estimate the proportions of different cell types at each spot of the ST data. Suppose the ST data comprises expressions of $n$ genes across $p$ locations, and the reference scRNA-seq data contains expressions of the same $n$ genes across cells annotated into $c$ distinct cell types. Our objective is to estimate the proportions of these $c$ cell types at each of the $p$ spatial locations. In ST data, we denote $\mathbf{X}$ as the $n \times p$ gene expression matrix, where $n$ informative genes are measured across $p$ spatial locations. From the scRNA-seq data, we derive a matrix $\mathbf{B}$ representing the $n \times c$ cell-type specific expression matrix for the same set of $n$ informative genes. Each element denotes the mean expression level of an informative gene in a specific cell type.

We introduce $\mathbf{V}$ as the $c \times p$ cell type composition matrix, where each column of $\mathbf{V}$ indicates the proportions of the $c$ cell types at each spatial location. Our objective is to estimate $\mathbf{V}$ using both $\mathbf{X}$ (the ST data) and $\mathbf{B}$ (derived from the scRNA-seq data). Naturally, the matrix $\mathbf{V}$ is constrained to be non-negative. Besides, since most cell types at each spot take little proportions, to increase the difference among cell type proportions within
one spot to better distinguish different regions, we add a sparsity prior to $\mathbf{V}$ to keep $\mathbf{V}$ getting too uniformly distributed. Laplacian distribution (also called double exponential distribution) is a Bayesian equivalent of $l_{1}$ norm regularization [47] and is applied to induce sparsity in Bayesian models [48. To induce non-negativity, we consider the positive part of Laplacian distribution and therefore the prior becomes an exponential distribution:

$$
\mathbf{V}_{i j} \mid \eta \sim \operatorname{Exp}(\eta)
$$

where $\mathbf{V}_{i j} \geq 0$ for $1 \leq i \leq c$ and $1 \leq j \leq p$. We represent it in a matrix form to combine with the Equation (6) and make an ablation study to show the effect of the exponential distribution in Supplementary Fig 4. We utilize a non-negative matrix factorization (NMF) model to link these matrices:

$$
\mathbf{X}=\mathbf{B V}+\mathbf{E}
$$

Since $\mathbf{X}$ is integer data, we transform it to continuous-valued data with 'lop1p' and assume the error matrix $\mathbf{E}$ as a Laplacian-structured Gaussian Markov Random Fields [49] whose precision matrix (inverse of covariance matrix) is the graph Laplacian:

$$
\mathbf{E} \sim \mathcal{M} \mathcal{N}\left(0, \sigma^{2} \mathbf{I}_{\mathbf{n}}, \mathbf{L}^{-1}\right)
$$

where $\sigma^{2}$ is a parameter controlling the scale of the error, $\mathbf{I}_{\mathbf{n}}$ is an identity matrix of size $n \times n$ and $\mathbf{L}$ is a $p \times p$ matrix representing the graph Laplacian of the spatial domain. $\mathbf{L}$ is defined as $\mathbf{L}=\mathbf{D}-\mathbf{A}$ where $\mathbf{A}$ is the weight (or adjacent) matrix and $\mathbf{D}$ is a diagonal matrix with $\mathbf{D}_{i i}=\sum_{j} \mathbf{A}_{i j}$. There are different choices for the definition of weight matrix $\mathbf{A}$, e.g., $0-1$ weighting, Gaussian kernel weighting and dot-product weighting [50]. We use Gaussian kernel weighting because it can measure distance between two nodes and is naturally suitable for a 2D spatial domain. We also consider the intensity difference of the aligned coordinates on the histology image (if available). Therefore, the weight matrix is defined as:

$$
\mathbf{A}_{i j}= \begin{cases}\exp \left(-\frac{\left\|s_{i}-s_{j}\right\|^{2}+\mu\left\|I_{i}-I_{j}\right\|^{2}}{\sigma_{A}^{2}}\right) & \text { if } i \neq j, \\ 0 & \text { if } i=j,\end{cases}
$$

where $\left\|s_{i}-s_{j}\right\|^{2}$ is the Euclidean distance between spots $i$ and $j$, and $\left\|I_{i}-I_{j}\right\|^{2}$ is the intensity difference obtained from the histology image. $\mu$ and $\sigma_{A}^{2}$ are constants controlling the weight of intensity difference and adapting spatial coordinates with different magnitudes. By setting $\mathbf{L}^{-1}$ as the column covariance matrix, we assume that the covariance of two spatially closer locations is larger, based on the fact that neighboring spatial locations have more similar gene expression 12. Give above, $\mathbf{X}$ is also following matrix normal distribution:

$$
\mathbf{X} \mid \mathbf{V}, \sigma^{2} \sim \mathcal{M} \mathcal{N}\left(\mathbf{B V}, \sigma^{2} \mathbf{I}_{\mathbf{n}}, \mathbf{L}^{-1}\right)
$$

The priors of the hyperparameters $\sigma^{2}$ and $\eta$ also need to be defined. The choice of prior depends on the conjugacy between the prior and the likelihood of the parameter. Therefore we set their priors as follows:

$$
\begin{aligned}
& \sigma^{2} \sim \operatorname{Inverse-Gamma}\left(a_{1}, b_{1}\right) \\
& \eta \sim \operatorname{Gamma}\left(a_{2}, b_{2}\right)
\end{aligned}
$$

where $a_{1}, a_{2}, b_{1}$ and $b_{2}$ are hyperparameters.

### 4.1.1 Posteriors

By applying Bayes' theorem, the conditional posterior distribution of $\mathbf{V}$ can be calculated by the formula:

$$
P(\mathbf{V} \mid \mathbf{X}) \propto P(\mathbf{X} \mid \mathbf{V}) P(\mathbf{V})
$$

Combining the non-negativity of $\mathbf{V}$ brought by the exponential distributions in the prior, We show that the posterior of $\mathbf{V}$ is a truncated matrix normal (TMN) distribution:

$$
\mathbf{V} \mid \mathbf{X} \sim \mathcal{T} \mathcal{M} \mathcal{N}\left(\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1}\left(\mathbf{B}^{T} \mathbf{X}-\sigma^{2} \eta \mathbf{J}_{\mathbf{V}} \mathbf{L}^{-1}\right),\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1} \sigma^{2}, \mathbf{L}^{-1}\right) \quad(\mathbf{V} \geq 0)
$$

where $\mathbf{J}_{\mathbf{V}}$ is an "all ones matrix" of the same dimension as $\mathbf{V}$. Specifically, the mean matrix $\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1}\left(\mathbf{B}^{T} \mathbf{X}-\right. \sigma^{2} \eta \mathbf{J}_{\mathbf{V}} \mathbf{L}^{-1}$ ) captures the cell type proportions at each location. The column covariance $\mathbf{L}^{-1}$ (spatial dimension) still involves the spatial structure defined by the graph Laplacian, while the row covariance $\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1} \sigma^{2}$ (cell type dimension) is determined by the cosine similarity of the cell types measured by $\mathbf{B}^{T} \mathbf{B}$.

Likewise, the conditional posterior distributions of $\sigma^{2}$ and $\eta$ are calculated as:

$$
\begin{gathered}
\left.\sigma^{2} \mid \mathbf{X}, \mathbf{V} \sim \operatorname{Inverse-Gamma}\left(a_{1}+n p / 2+c p / 2, \quad b_{1}+\frac{\operatorname{tr}\left(\mathbf{L E}^{T} \mathbf{E}\right)}{2}\right)\right) \\
\left.\left.\eta \left\lvert\, \mathbf{V} \sim \operatorname{Gamma}\left(a_{2}+c p, \quad 1 /\left(\operatorname{tr}\left(\mathbf{J}_{\mathbf{V}}{ }^{T} \mathbf{V}\right)\right)+\frac{1}{b_{2}}\right)\right.\right)\right)
\end{gathered}
$$

When choosing the hyperparameters, we want to reduce the influence of them. We set $a_{1,2}, b_{1}$ and $\frac{1}{b_{2}}$ to be really small, e.g., $10^{-4}$ to exploit more information from data.

### 4.1.2 Overview of the algorithm

Starting with an arbitrary initial value of $\mathbf{V}$, we iteratively sample the parameters from their conditional posterior distributions using Gibbs sampling. The proposed algorithm is summarized in Algorithm 1 Unlike general

```
Algorithm 1 BASIN
    Require: \(\mathbf{X}, \mathbf{B}, \sigma_{A}^{2}, \mu, \mathbf{L}\), sample size \(T\), burn-in period \(T_{b}\)
    Initialize: choose \(\mathbf{V}^{0}\)
    for \(t=1\) to \(T\) do
        Sample \(\left(\sigma^{2}\right)^{t}\) given \(\mathbf{V}^{t-1}\) with Equation (10)
        Sample \(\eta^{t}\) given \(\mathbf{V}^{t-1}\) with Equation (11)
        Sample \(\mathbf{V}^{t}\) given \(\left(\sigma^{2}\right)^{t}\) and \(\eta^{t}\) with Equation (9)
        Normalize \(\mathbf{V}\) so that the column sum is one
    end for
    Discard \(\mathbf{V}^{1} \cdots \mathbf{V}^{T_{b}}\)
```

Bayesian NMF models with two unknown matrices, our deconvolution model is reference-based with one unknown matrix. Therefore, in practice, the Gibbs sampler can converge within several steps and the burn-in period can be set within five steps. The value of $\sigma_{A}^{2}$ depends on the magnitude of coordinate. For most ST data like MOB and PDAC where the 2D coordinates lie in between $10^{1}$ and $10^{2}$, we set $\sigma_{A}^{2}=0.1$. For the seqFISH + data where the coordinates have magnitude between $10^{3}$ and $10^{4}$, we set $\sigma_{A}^{2}=10^{4} . \mu$ is a hyperparameter controlling the weight of intensity information in the H\&E image. It is set to be one after normalizing image and can be ignored when image is not available.

### 4.2 Sampling from truncated matrix normal distribution

Sampling from a high dimensional truncated matrix normal distribution can be challenging and time consuming. Sampling from matrix normal distribution typically requires first transforming it to a multivariate normal distribution (MVN):

$$
\mathbf{X} \sim \mathcal{M} \mathcal{N}(\mathbf{M}, \mathbf{A}, \mathbf{B}) \quad \rightarrow \quad \operatorname{vec}(\mathbf{X}) \sim \mathcal{M} \mathcal{V} \mathcal{N}(\operatorname{vec}(\mathbf{M}), \mathbf{B} \otimes \mathbf{A})
$$

and then sample from the multivariate normal distribution. Truncated matrix normal distribution can be sampled in the same way. However, when the dimension gets high, sampling from truncated multivariate normal distribution becomes increasingly difficult because the probabilities of each variate are very small and the acceptance rates are almost zero when using Markov chain Monte Carlo methods. In our case, the sampling dimension is the size of $\mathbf{V}(c \times p)$, which is usually in the magnitude of $10^{3}$ or $10^{4}$. Therefore, sampling from the posterior distribution of $\mathbf{V}$ consumes more than $90 \%$ time in our algorithm.

Although there are many implementation schemes proposed to sample from truncated MVN 51-53, we utilize the sampling method proposed in 41 because it's the only one that can deal with high dimension variate according to our test. In short words, it decomposes the truncated multivariate normal distribution into a set of truncated univariate distributions and samples with Gibbs sampling. Although it's efficient, it's still computationally expensive when the parameter dimension is too high (Supplementary Fig. 5b). To further improve the efficiency, we simplify the sampling method 41 by only considering the truncated range from zero to infinity, while the method originally considers general ranges. Besides, we accelerate the method by removing the unnecessary computations which takes the most time in sampling (details are discussed in Supplementary Information). By simplification, we make the sampling algorithm about $10-1000$ times faster than before in
our method (Supplementary Fig. 5b). We compare the computation times spent on the PDAC-A data in Supplementary Fig 5a, showing that BASIN is still efficient when sampling from more than 8000 dimensions.

However, sampling from a a high-dimensional truncated multivariate normal distribution is essentially challenging due to the fact that the probability of each variate is extremely small as the dimension increases. Therefore, in our algorithm, when sampling dimension (number of locations times number of cell types) is larger than $10^{4}$, we apply rectified matrix normal distribution instead, which is as fast as sampling from regular matrix normal distributions. Specifically, we sample from regular matrix normal distribution first, and then set the negative values to be zeros. We have observed that rectified matrix normal distribution is a good approximation of truncated matrix normal distribution (Supplementary Fig. 6).

### 4.3 Data preprocessing and acquisition

All the ST data and scRNA-seq data we used are publicly available. Each ST dataset includes two data tables: one containing the gene expression abundance for each spot, and another containing the spatial coordinates ( X and Y ) of the spots. The scRNA-seq data includes the abundance of each gene for each single cell and a meta data revealing the type of cells. To make sure the two kind of data corresponding well with each other, we preprocess the data using the similar methods as in 12,14 . First, we filter out the genes not showing in any spots or cells for the two data, and the genes that just show in one of the STdata and scRNA-seq data. Second, the ST and scRNA-seq data are re-ordered by genes, spots, and cells to align perfectly with their corresponding locations and metadata. Third, we construct the cell-type expression matrix using the scRNA-seq data and its meta data, by computing the average gene expression levels for each cell type. Fourth, we calculate the VMR (variance-to-mean ratio) for each gene in each cell type and mean VMR across the cell types. We remove the genes with top $10 \%$ largest mean VMR. Finally, to filter out the genes with low expression heterogeneity and reduce the data size without losing performance, we only keep the genes whose highest mean expression level in one cell type is at least 1.25 -log-fold higher than the mean expression level across all remaining cell types.

The mouse olfactory bulb (MOB) data 37 is downloaded from the website of Spatial Research Lab (https://www.spatialresearch.org/). We choose the replicate 12, containing 16034 genes and 282 spots. The corresponding scRNA-seq data is obtained from the website of GSE121891 [38], containing 18560 genes and 12801 cells including 5 cell types. After preprocessing, there are 4284 genes retained. The human pancreatic ductal adenocarcinoma (PDAC) data is downloaded from the website of GSE111672 35, including both ST data and scRNA-seq data generated from the same samples. In the paper we choose PDAC-A ST1 from patient ID GSM3036911, and PDAC-B ST1 from patient ID GSM3405534 as the ST data, and PDAC-A-indrop and PDAC-B-indrop as the scRNA-seq data. The original PDAC-A ST1 data has 19738 genes and 428 spots, and PDAC-A-indrop containing 19738 genes and 1926 cells from 20 cell types. After preprocessing, there are 10431 genes retained. The original PDAC-B ST1 data 19738 genes and 224 spots, and PDAC-B-indrop contains 19738 genes and 1733 cells from 13 cell types. After preprocessing, there are 9401 genes retained. The seqFISH+ mouse cortex data is downloaded from Cai Lab's Github 39 (https://github.com/CaiGroup/seqFISH-PLUS), which contains 10000 genes and 524 spots. The corresponding scRNA-seq data is obtained from GSE102827 [40, containing 25187 genes and 48266 cells from 8 main types and 33 subtypes. We divide them into 12 types as done in CARD to distinguish different excitatory layers. After preprocessing, there are 9498 genes retained.

### 4.4 Simulation strategy

We use simulated data to evaluate the performance of BASIN and compare with other deconvolution methods. The simulation is based on the real PDAC-A data 35 . We first divide the spatial domain into three regions according to the histological annotation shown in Supplementary Information. We assume there are three cell types (acinar, cancer cluster 1 and terminal ductal cells) whose proportions follow different probability distributions and are dominant respectively in the three regions. Four simulation schemes are designed as below:

In simulation 1, each region has a piecewise constant cell type composition in which the dominant cell type has a proportion of 0.7 and the other two types take 0.15 . Simulation 2 is based on simulation 1 , in which the simulated cell type proportion data is smoothed with a large Gaussian kernel of size $40 \times 40$ with $\sigma=0.5$. In simulation 3, the cell type proportions follow Gaussian distributions whose mean values are equal to those of simulation 1 and the standard deviation is 0.05 . Then the data is smoothed with the same Gaussian kernel used in simulation 2. In simulation 4, the three cell types follow different Dirichlet distributions in the three regions. For each region. the concentration parameter is equal to 3 for the main type and equal to 1 for the other two types.

For each of the simulations, the simulated cell type proportion matrix is limited to be nonnegative, with each column normalized to sum to one. These simulated $V$ 's are treated as ground truth for each of the four simulations and compared with results of various methods. We extract the gene profiles of these three cell types
from the real scRNA-seq data of PDAC-A. Then we obtain the simulated transcriptomic data based on the NMF model $\mathbf{X}=\mathbf{B V}$ plus a zero mean Gaussian noise with $\sigma=0.5$ and limited to be nonnegative again.

### 4.5 Evaluation metrics

Three metrics are used to evaluate the simulations. Root mean squared error (RMSE) estimates the average errors between the values and their ground truth, defined as

$$
\operatorname{RMSE}(\mathbf{a}, \mathbf{b})=\frac{1}{\sqrt{p}}\|\mathbf{a}-\mathbf{b}\|_{F}
$$

where $p$ is the number of entries of $\mathbf{a}$ and $\mathbf{b}$. Structural similarity index measure (SSIM) is used to measure the similarity between two images. Compared with RMSE, SSIM considers not only absolute errors, but also difference of structural information, defined as

$$
\operatorname{SSIM}(\mathbf{a}, \mathbf{b})=\frac{\left(2 \mu_{a} \mu_{b}+c_{1}\right)\left(2 \sigma_{a b}+c_{2}\right)}{\left(\mu_{a}^{2}+\mu_{b}^{2}+c_{1}\right)\left(\sigma_{a}^{2}+\sigma_{b}^{2}+c_{2}\right)}
$$

where $\mu_{a}, \sigma_{a}^{2}$ are the mean and variance of $\mathbf{a}$, and $\sigma_{a b}$ is the covariance measuring the structure similarity. $c_{1}$ and $c_{2}$ are parameters to stabilize the division with weak denominator. In our case, we compute the average SSIM of each cell type. Jensen-Shannon distance (JSD) measures the difference between two probability distributions and is defined as the square root of Jensen-Shannon divergence. The Jensen-Shannon divergence is defined as:

$$
\mathrm{JSD}(\mathbf{a}, \mathbf{b})=\frac{1}{2} \mathrm{KL}\left(\mathbf{a} \| \frac{1}{2}(\mathbf{a}+\mathbf{b})\right)+\frac{1}{2} \mathrm{KL}\left(\mathbf{b} \| \frac{1}{2}(\mathbf{a}+\mathbf{b})\right)
$$

where KL means the Kullback-Leibler divergence. For each location, the cell type can be regarded as following a certain distribution. In this case we compute the average JSD of each location.

### 4.6 Benchmark methods

To evaluate the performance of our proposed method, we compared it with several bulk tissue cell-type deconvolution methods using reference scRNA-seq data: Stereoscope, SPOTlight, SpatialDWLS, RCTD and CARD. Stereoscope 19 models both the scRNAseq and spatial data using negative binomial distribution. The parameters are estimated by MLE using gradient descent for some steps and then the cell type proportions are obtained with MAP estimate in a similar way. In our analysis, we set the learning rate as 0.01 , batch size as 100 and number of steps as 5000 for each data. SPOTlight [10] employs a seeded NMF approach to deconvolute ST spots. The process starts by initializing cell-type marker genes and then applies non-negative least squares (NNLS) to further deconvolute the ST capture locations (spots). In both simulation and real data analysis, all genes are used, with a maximum of 75 cells selected per cell type. All other parameters were set to their default values. SpatialDWLS 13 estimates cell proportions by first identifying the cell types likely present at each location. This is accomplished using marker genes through differential expression analysis with Giotto 54. It then applies dampened weighted least squares (DWLS) 55 to infer the fraction of each selected cell type. An enrichment score cutoff of 2 is used to select the cell types, and all other parameters are set to their default values. RCTD [15] models the observed gene counts in each pixel of ST data using a Poisson distribution. It links individual cell types with a generalized hierarchical linear regression model, and the cell type proportions are inferred using the maximum-likelihood estimation (MLE) method. For both simulation and real data analysis, we used the "doublet" model, which constrains the number of cell types to two per pixel, as recommended in 15 to reduce overfitting and increase power. All other parameters were set to their default values. CARD 12 models the the problem as NMF and assumes the cell type composition as a conditional autoregressive (CAR) model. The parameters are estimated with MAP estimation using a NMF framework. In our analysis all hyperparameters were set to their default values.

## Declarations

- Funding: Dr. Liangliang Zhang is supported by his junior faculty start-up grant BGT630267, funded by School of Medicine at Case Western Reserve University. Additional support is provided by the Data Management and Statistics Core under grant number 5P30AG072959, led by Dr. Jonathan L. Haines at Cleveland Alzheimer's Disease Research Center.
- Competing interests: The authors declare no competing interests.
- Ethics approval and consent to participate: Not applicable.
- Consent for publication: All authors agree
- Data availability: All datasets are publicly available.
- Materials availability: Not applicable.
- Code availability: https://github.com/Jiasen-Zhang/BASIN-Deconvolution
- Author contribution: Guo and Zhang provide the main ideas and guide Zhang who implements the ideas numerically and prepares the draft of the manuscript. Qiao helps with comparison with some of the methods. All coauthors contribute to polish the manuscript.


## References

[1] Luyi Tian, Fei Chen, and Evan Z Macosko. The expanding vistas of spatial transcriptomics. Nature Biotechnology, 41(6):773-782, 2023.
[2] Anjali Rao, Dalia Barkley, Gustavo S França, and Itai Yanai. Exploring tissue architecture using spatial transcriptomics. Nature, 596(7871):211-220, 2021.
[3] Lambda Moses and Lior Pachter. Museum of spatial transcriptomics. Nature methods, 19(5):534-546, 2022.
[4] Yang Jin, Yuanli Zuo, Gang Li, Wenrong Liu, Yitong Pan, Ting Fan, Xin Fu, Xiaojun Yao, and Yong Peng. Advances in spatial transcriptomics and its applications in cancer research. Molecular Cancer, 23(1):129, 2024.
[5] Sanjay Jain and Michael T Eadon. Spatial transcriptomics in health and disease. Nature Reviews Nephrology, 20(10):659-671, 2024.
[6] Ed Lein, Lars E Borm, and Sten Linnarsson. The promise of spatial transcriptomics for neuroscience in the era of molecular cell typing. Science, 358(6359):64-69, 2017.
[7] Manuel Saiselet, Joël Rodrigues-Vitória, Adrien Tourneur, Ligia Craciun, Alex Spinette, Denis Larsimont, Guy Andry, Joakim Lundeberg, Carine Maenhaut, and Vincent Detours. Transcriptional output, cell-type densities, and normalization in spatial transcriptomics. Journal of molecular cell biology, 12(11):906-908, 2020.
[8] Xinyang Zhang, Evan W. Miller, Joshua T. Vogelstein, Mauro Maggioni, and Risi Kondor. Spatial deconvolution of multiplexed transcriptomic data with multimodal learning. Nature Biotechnology, 38(9):1028-1037, September 2020.
[9] Zexian Zeng, Yawei Li, Yiming Li, and Yuan Luo. Statistical and machine learning methods for spatially resolved transcriptomics data analysis. Genome biology, 23(1):83, 2022.
[10] Marc Elosua-Bayes, Paula Nieto, Elisabetta Mereu, Ivo Gut, and Holger Heyn. Spotlight: seeded nmf regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic acids research, 49(9):e50-e50, 2021.
[11] Samuel G Rodriques, Robert R Stickels, Aleksandrina Goeva, Carly A Martin, Evan Murray, Charles R Vanderburg, Joshua Welch, Linlin M Chen, Fei Chen, and Evan Z Macosko. Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution. Science, 363(6434):1463-1467, 2019.
[12] Ying Ma and Xiang Zhou. Spatially informed cell-type deconvolution for spatial transcriptomics. Nature biotechnology, 40(9):1349-1359, 2022.
[13] Rui Dong and Guo-Cheng Yuan. Spatialdwls: accurate deconvolution of spatial transcriptomic data. Genome biology, 22(1):145, 2021.
[14] Xuran Wang, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao Li. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications, 10(1):380, 2019.
[15] Dylan M Cable, Evan Murray, Luli S Zou, Aleksandrina Goeva, Evan Z Macosko, Fei Chen, and Rafael A Irizarry. Robust decomposition of cell type mixtures in spatial transcriptomics. Nature biotechnology, 40(4):517-526, 2022.
[16] Vitalii Kleshchevnikov, Artem Shmatko, Emma Dann, Alexander Aivazidis, Hamish W King, Tong Li, Rasa Elmentaite, Artem Lomakin, Veronika Kedlian, Adam Gayoso, et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nature biotechnology, 40(5):661-671, 2022.
[17] Romain Lopez, Baoguo Li, Hadas Keren-Shaul, Pierre Boyeau, Merav Kedmi, David Pilzer, Adam Jelinski, Ido Yofe, Eyal David, Allon Wagner, et al. Destvi identifies continuums of cell types in spatial transcriptomics data. Nature biotechnology, 40(9):1360-1369, 2022.
[18] Patrick Danaher, Youngmi Kim, Brenn Nelson, Maddy Griswold, Zhi Yang, Erin Piazza, and Joseph M Beechem. Advances in mixed cell deconvolution enable quantification of cell types in spatial transcriptomic data. Nature communications, 13(1):385, 2022.
[19] Alma Andersson, Joseph Bergenstråhle, Michaela Asp, Ludvig Bergenstråhle, Aleksandra Jurek, José Fernández Navarro, and Joakim Lundeberg. Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography. Communications biology, 3(1):565, 2020.
[20] Dongqing Sun, Zhaoyang Liu, Taiwen Li, Qiu Wu, and Chenfei Wang. Stride: accurately decomposing and integrating spatial transcriptomics using single-cell rna sequencing. Nucleic Acids Research, 50(7):e42-e42, 2022.
[21] Zixuan Cang and Qing Nie. Inferring spatial and signaling relationships between cells from single cell transcriptomic data. Nature communications, 11(1):2084, 2020.
[22] Noa Moriel, Enes Senel, Nir Friedman, Nikolaus Rajewsky, Nikos Karaiskos, and Mor Nitzan. Novosparc: flexible spatial reconstruction of single-cell gene expression with optimal transport. Nature protocols, 16(9):4177-4200, 2021.
[23] Qianqian Song and Jing Su. Dstg: deconvoluting spatial transcriptomics data through graph-based artificial intelligence. Briefings in bioinformatics, 22(5):bbaa414, 2021.
[24] Haoyang Li, Hanmin Li, Juexiao Zhou, and Xin Gao. Sd2: spatially resolved transcriptomics deconvolution through integration of dropout and spatial information. Bioinformatics, 38(21):4878-4884, 2022.
[25] Tommaso Biancalani, Gabriele Scalia, Lorenzo Buffoni, Raghav Avasthi, Ziqing Lu, Aman Sanger, Neriman Tokcan, Charles R Vanderburg, Åsa Segerstolpe, Meng Zhang, et al. Deep learning and alignment of spatially resolved single-cell transcriptomes with tangram. Nature methods, 18(11):1352-1362, 2021.
[26] Kyle Coleman, Jian Hu, Amelia Schroeder, Edward B Lee, and Mingyao Li. Spadecon: cell-type deconvolution in spatial transcriptomics with semi-supervised learning. Communications Biology, 6(1):378, 2023.
[27] Brendan F Miller, Feiyang Huang, Lyla Atta, Arpan Sahoo, and Jean Fan. Reference-free cell type deconvolution of multi-cellular pixel-resolution spatially resolved transcriptomics data. Nature communications, 13(1):2339, 2022.
[28] Bárbara Andrade Barbosa, Saskia D van Asten, Ji Won Oh, Arantza Farina-Sarasqueta, Joanne Verheij, Frederike Dijk, Hanneke WM van Laarhoven, Bauke Ylstra, Juan J Garcia Vallejo, Mark A van de Wiel, et al. Bayesian log-normal deconvolution for enhanced in silico microdissection of bulk gene expression data. Nature communications, 12(1):6106, 2021.
[29] Emelie Berglund, Jonas Maaskola, Niklas Schultz, Stefanie Friedrich, Maja Marklund, Joseph Bergenstråhle, Firas Tarish, Anna Tanoglidi, Sanja Vickovic, Ludvig Larsson, et al. Spatial maps of prostate cancer transcriptomes reveal an unexplored landscape of heterogeneity. Nature communications, 9(1):2419, 2018.
[30] Meng Zhang, Joel Parker, Lingling An, Yiwen Liu, and Xiaoxiao Sun. Flexible analysis of spatial transcriptomics data (fast): a deconvolution approach. BMC bioinformatics, 26(1):35, 2025.
[31] Jian Hu, Xiangjie Li, Kyle Coleman, Amelia Schroeder, Nan Ma, David J Irwin, Edward B Lee, Russell T Shinohara, and Mingyao Li. Spagcn: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network. Nature methods, 18(11):13421351, 2021.
[32] Edward Zhao, Matthew R Stone, Xing Ren, Jamie Guenthoer, Kimberly S Smythe, Thomas Pulliam, Stephen R Williams, Cedric R Uytingco, Sarah EB Taylor, Paul Nghiem, et al. Spatial transcriptomics at subspot resolution with bayesspace. Nature biotechnology, 39(11):1375-1384, 2021.
[33] Jun Lu and Xuanyu Ye. Flexible and hierarchical prior for bayesian nonnegative matrix factorization. arXiv preprint arXiv:2205.11025, 2022.
[34] Mikkel N Schmidt, Ole Winther, and Lars Kai Hansen. Bayesian non-negative matrix factorization. In Independent Component Analysis and Signal Separation: 8th International Conference, ICA 2009, Paraty, Brazil, March 15-18, 2009. Proceedings 8, pages 540-547. Springer, 2009.
[35] Reuben Moncada, Dalia Barkley, Florian Wagner, Marta Chiodin, Joseph C Devlin, Maayan Baron, Cristina H Hajdu, Diane M Simeone, and Itai Yanai. Integrating microarray-based spatial transcriptomics and single-cell rna-seq reveals tissue architecture in pancreatic ductal adenocarcinomas. Nature biotechnology, 38(3):333-342, 2020.
[36] Meritxell Rovira, Sherri-Gae Scott, Andrew S Liss, Jan Jensen, Sarah P Thayer, and Steven D Leach. Isolation and characterization of centroacinar/terminal ductal progenitor cells in adult mouse pancreas. Proceedings of the National Academy of Sciences, 107(1):75-80, 2010.
[37] Patrik L Ståhl, Fredrik Salmén, Sanja Vickovic, Anna Lundmark, José Fernández Navarro, Jens Magnusson, Stefania Giacomello, Michaela Asp, Jakub O Westholm, Mikael Huss, et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science, 353(6294):78-82, 2016.
[38] Burak Tepe, Matthew C Hill, Brandon T Pekarek, Patrick J Hunt, Thomas J Martin, James F Martin, and Benjamin R Arenkiel. Single-cell rna-seq of mouse olfactory bulb reveals cellular heterogeneity and activity-dependent molecular census of adult-born neurons. Cell reports, 25(10):2689-2703, 2018.
[39] Chee-Huat Linus Eng, Michael Lawson, Qian Zhu, Ruben Dries, Noushin Koulena, Yodai Takei, Jina Yun, Christopher Cronin, Christoph Karp, Guo-Cheng Yuan, et al. Transcriptome-scale super-resolved imaging in tissues by rna seqfish+. Nature, 568(7751):235-239, 2019.
[40] Sinisa Hrvatin, Daniel R Hochbaum, M Aurel Nagy, Marcelo Cicconet, Keiramarie Robertson, Lucas Cheadle, Rapolas Zilionis, Alex Ratner, Rebeca Borges-Monroy, Allon M Klein, et al. Single-cell analysis of experience-dependent transcriptomic states in the mouse visual cortex. Nature neuroscience, 21(1):120129, 2018.
[41] Yifang Li and Sujit K Ghosh. Efficient sampling methods for truncated multivariate normal and student-t distributions subject to linear inequality constraints. Journal of Statistical Theory and Practice, 9:712-732, 2015.
[42] Ameur M Manceur and Pierre Dutilleul. Maximum likelihood estimation for the tensor normal distribution: Algorithm, minimum sample size, and empirical bias and dispersion. Journal of Computational and Applied Mathematics, 239:37-49, 2013.
[43] Myson Burch, Jiasen Zhang, Gideon Idumah, Hakan Doga, Richard Lartey, Lamis Yehia, Mingrui Yang, Murat Yildirim, Mihriban Karaayvaz, Omar Shehab, et al. Towards quantum tensor decomposition in biomedical applications. arXiv preprint arXiv:2502.13140, 2025.
[44] Tianci Song, Charles Broadbent, and Rui Kuang. Gntd: reconstructing spatial transcriptomes with graphguided neural tensor decomposition informed by spatial and functional relations. Nature communications, 14(1):8276, 2023.
[45] Zhuliu Li, Tianci Song, Jeongsik Yong, and Rui Kuang. Imputation of spatially-resolved transcriptomes by graph-regularized tensor completion. PLoS computational biology, 17(4):e1008218, 2021.
[46] Charles Broadbent, Tianci Song, and Rui Kuang. Deciphering high-order structures in spatial transcriptomes with graph-guided tucker decomposition. Bioinformatics, 40(Supplement_1):i529-i538, 2024.
[47] Pierre Garrigues and Bruno Olshausen. Group sparse coding with a laplacian scale mixture prior. Advances in neural information processing systems, 23, 2010.
[48] Facundo Costa, Hadj Batatia, Lotfi Chaari, and Jean-Yves Tourneret. Sparse eeg source localization using bernoulli laplacian priors. IEEE Transactions on Biomedical Engineering, 62(12):2888-2898, 2015.
[49] Jiaxi Ying, José Vinícius de Miranda Cardoso, and Daniel Palomar. Minimax estimation of laplacian constrained precision matrices. In International Conference on Artificial Intelligence and Statistics, pages 3736-3744. PMLR, 2021.
[50] Deng Cai, Xiaofei He, Jiawei Han, and Thomas S Huang. Graph regularized nonnegative matrix factorization for data representation. IEEE transactions on pattern analysis and machine intelligence, 33(8):1548-1560, 2010.
[51] Stefan Wilhelm and BG Manjunath. tmvtnorm: A package for the truncated multivariate normal distribution. 2010.
[52] Zdravko I Botev. The normal law under linear restrictions: simulation and estimation via minimax tilting. Journal of the Royal Statistical Society Series B: Statistical Methodology, 79(1):125-148, 2017.
[53] Jonathan Goodman and Jonathan Weare. Ensemble samplers with affine invariance. Communications in applied mathematics and computational science, 5(1):65-80, 2010.
[54] Jiaji George Chen, Joselyn Cristina Chávez-Fuentes, Matthew O'Brien, Junxiang Xu, Edward Ruiz, Wen Wang, Iqra Amin, Irzam Sarfraz, Pratishtha Guckhool, Adriana Sistig, Guo-Cheng Yuan, and Ruben Dries. Giotto suite: a multi-scale and technology-agnostic spatial multi-omics analysis ecosystem. bioRxiv, 2023.
[55] Daphne Tsoucas, Rui Dong, Haide Chen, Qian Zhu, Guoji Guo, and Guo-Cheng Yuan. Accurate estimation of cell-type composition from gene expression data. Nature Communications, 10(1):2975, 2019.

## Supplementary Information

## 1 Supplementary figures

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-18.jpg?height=652&width=796&top_left_y=560&top_left_x=674)
Figure 1: Spatial division of our simulated data. The simulation is based on the real PDAC-A data 35. We divide the spatial domain into three regions according to the histological annotation, and assume there are three cell types (acinar, cancer cluster A and terminal ductal cells) whose proportions follow different probability distributions and are dominant respectively in the three regions.

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-19.jpg?height=991&width=1737&top_left_y=247&top_left_x=193)
Figure 2: The results of the four simulation studies compared with the ground truth. Blue: acinar cells. Yellow: terminal ductal cells. Red: cancer cells

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-20.jpg?height=1000&width=1745&top_left_y=242&top_left_x=191)
Figure 3: The spot-wise RMSE of the four simulations between the six compared methods and the ground truth.

## 2 Details of the method

### 2.1 Matrix normal distribution

Matrix normal distribution is a generalization of multivariate normal distribution. Suppose there is a random matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$ whose entries follow normal distributions, and the mean values of $\mathbf{X}$ is a matrix $\mathbf{M}$ of size $n \times p$. Each column of $\mathbf{X}$ follows a multivariate normal distribution with covariance matrix $\mathbf{U}$ of size $n \times n$. Each row of $\mathbf{X}$ follows a multivariate normal distribution with covariance matrix $\mathbf{V}$ of size $p \times p$. Then $\mathbf{X}$ follows a matrix normal distribution $\mathbf{X} \sim \mathcal{M} \mathcal{N}(\mathbf{M}, \mathbf{U}, \mathbf{V})$ where $\mathbf{U}$ and $\mathbf{V}$ are column and row covariances. Its probability density function is

$$
P(\mathbf{X})=\frac{\exp \left(-\frac{1}{2} \operatorname{tr}\left(\mathbf{V}^{-1}(\mathbf{X}-\mathbf{M})^{T} \mathbf{U}^{-1}(\mathbf{X}-\mathbf{M})\right)\right.}{(2 \pi)^{n p / 2} \operatorname{det}(\mathbf{U})^{p / 2} \operatorname{det}(\mathbf{V})^{n / 2}}
$$

On the other hand, if the probability density function of a random matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$ satisfies

$$
P(\mathbf{X}) \propto \exp \left(-\frac{1}{2} \operatorname{tr}\left(\mathbf{A} \mathbf{X}^{T} \mathbf{B X}-\mathbf{2} \mathbf{C X}\right)\right)
$$

where $\mathbf{A}$ and $\mathbf{B}$ are symmetric, then $\mathbf{X}$ can be represented in form of matrix normal distribution:

$$
\begin{aligned}
P(\mathbf{X}) & \propto \exp \left(-\frac{1}{2} \operatorname{tr}\left(\mathbf{A} \mathbf{X}^{T} \mathbf{B} \mathbf{X}-\mathbf{2} \mathbf{C X}\right)\right) \\
& \propto \exp \left(-\frac{1}{2} \operatorname{tr}\left(\mathbf{A} \mathbf{X}^{T} \mathbf{B} \mathbf{X}-\mathbf{2} \mathbf{A} \mathbf{A}^{-1} \mathbf{C B}^{-1} \mathbf{B} \mathbf{X}\right)\right) \\
& \propto \exp \left(-\frac{1}{2} \operatorname{tr}\left(\mathbf{A} \mathbf{X}^{T} \mathbf{B} \mathbf{X}-\mathbf{2} \mathbf{A} \mathbf{A}^{-1} \mathbf{C B}^{-1} \mathbf{B} \mathbf{X}+\mathbf{A}\left(\mathbf{A}^{-1} \mathbf{C B}^{-1}\right) \mathbf{B}\left(\mathbf{B}^{-T} \mathbf{C}^{T} \mathbf{A}^{-T}\right)\right)\right) \\
& \propto \exp \left(-\frac{1}{2} \operatorname{tr}\left(\mathbf { A } \left(\mathbf{X}-\left(\mathbf{B}^{-T} \mathbf{C}^{T} \mathbf{A}^{-T}\right)^{T} \mathbf{B}\left(\mathbf{X}-\left(\mathbf{B}^{-T} \mathbf{C}^{T} \mathbf{A}^{-T}\right)\right)\right.\right.\right.
\end{aligned}
$$

which is equivalent to the pdf of the matrix normal distribution

$$
\mathcal{M} \mathcal{N}\left(\mathbf{B}^{-T} \mathbf{C}^{T} \mathbf{A}^{-T}, \mathbf{B}^{-1}, \mathbf{A}^{-1}\right)
$$

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-21.jpg?height=622&width=1482&top_left_y=268&top_left_x=330)
Figure 4: We make an experiment to show the effect of the exponential distribution in the prior of $V$. It not only constrains the values of $V$ to be nonnegative, but also introduce sparsity to $V$ to reduce noise and better distinguish the cells. We use the MOB data as an example and plot the histograms of the five cell types without (left) and with (right) the exponential distribution in the prior. We can see that by introducing sparsity, there are more spots with trivial proportion of EPL-IN (blue) and OSs (red) on the right. It matches the fact that EPL-IN is trivial in the tissue and OSs only exist in the a small area. On the other hand, there are more spots with high proportions of GC cells (orange) on the right, which is also reasonable since GC cells mainly distributed in the central dominant area. Therefore, the exponential distribution in the prior help us better identify the dominant cell type of each spot, and avoid getting too uniform distributed results.

### 2.2 Posterior of V

In the main text we have defined the following likelihood and priors

$$
\begin{gathered}
\mathbf{X} \mid \mathbf{V}, \sigma^{2} \sim \mathcal{M} \mathcal{N}\left(\mathbf{B V}, \sigma^{2} \mathbf{I}_{\mathbf{n}}, \mathbf{L}^{-1}\right) \\
\mathbf{V}_{i j} \mid \eta \sim \operatorname{Exp}(\eta) \quad 1 \leq i \leq c \quad 1 \leq j \leq p
\end{gathered}
$$

Their probability density function (pdf) of likelihood satisfies

$$
\begin{aligned}
P\left(\mathbf{X} \mid \mathbf{V}, \sigma^{2}\right) & \propto \frac{1}{\left|\sigma^{2}\right|^{n p / 2}} \exp \left(-\frac{1}{2 \sigma^{2}} \operatorname{tr}\left(\mathbf{L}(\mathbf{X}-\mathbf{B} \mathbf{V})^{T}(\mathbf{X}-\mathbf{B} \mathbf{V})\right)\right) \\
& \propto \exp \left(-\frac{1}{2 \sigma^{2}} \operatorname{tr}\left(\mathbf{L}(\mathbf{X}-\mathbf{B} \mathbf{V})^{T}(\mathbf{X}-\mathbf{B} \mathbf{V})\right)\right) \\
& \propto \exp \left(-\frac{1}{2 \sigma^{2}} \operatorname{tr}\left(\mathbf{L V}^{T} \mathbf{B}^{T} \mathbf{B} \mathbf{V}-\mathbf{2} \mathbf{L} \mathbf{X}^{T} \mathbf{B} \mathbf{V}+\mathbf{L} \mathbf{X}^{T} \mathbf{X}\right)\right) \\
& \propto \exp \left(-\frac{1}{2 \sigma^{2}} \operatorname{tr}\left(\mathbf{L} \mathbf{V}^{T} \mathbf{B}^{T} \mathbf{B} \mathbf{V}-\mathbf{2} \mathbf{L} \mathbf{X}^{T} \mathbf{B} \mathbf{V}\right)\right)
\end{aligned}
$$

The pdf of $\mathbf{V}_{i j} \mid \eta$ satisfies

$$
\begin{gathered}
P\left(\mathbf{V}_{i j} \mid \eta\right)=\eta \exp \left(-\eta \mathbf{V}_{i j}\right) \quad \mathbf{V}_{i j} \geq 0 \\
P(\mathbf{V} \mid \eta)=\prod_{i, j} P\left(\mathbf{V}_{i j} \mid \eta\right)=\eta^{c p} \exp \left(-\eta \sum \mathbf{V}_{i j}\right)=\eta^{c p} \exp \left(-\eta \operatorname{tr}\left(\mathbf{J}_{\mathbf{V}}{ }^{T} \mathbf{V}\right)\right)
\end{gathered}
$$

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-22.jpg?height=547&width=1685&top_left_y=240&top_left_x=193)
Figure 5: Computational efficiency. a, We use the PDAC-A data to compare the computational efficiency of the six compared methods. All the methods are implemented with 11th Gen Intel(R) Core(TM) i7-11800H @ 2.30 GHz and Stereoscope also utilizes NVIDIA GeForce RTX 3080 Laptop GPU. We run Stereoscope for 5000 epochs for both the spatial transcriptomics and scRNA-seq data, and the time of BASIN means that for one sample. The other four methods are implemented in the default settings. b, The times of generating one sample from a truncated multivariate normal distribution vs sampling dimensions. Red: the algorithm in 41. Blue: our improved algorithm. In practice we generate one sample by ignoring first five samples and take the sixth sample.

where $\mathbf{J}_{\mathbf{V}}$ is an "all ones matrix" of the same size as $\mathbf{V}$. According to Bayes's formula,

$$
\begin{aligned}
P(\mathbf{V} \mid \mathbf{X}) & \propto P(\mathbf{X} \mid \mathbf{V}) P(\mathbf{V}) \\
& \propto \exp \left(-\frac{1}{2 \sigma^{2}} \operatorname{tr}\left(\mathbf{L V}^{T} \mathbf{B}^{T} \mathbf{B V}-\mathbf{2} \mathbf{L} \mathbf{X}^{T} \mathbf{B V}\right)\right) \exp \left(-\eta \operatorname{tr}\left(\mathbf{J}_{\mathbf{V}}^{T} \mathbf{V}\right)\right) \\
& \propto \exp \left(-\frac{1}{2 \sigma^{2}} \operatorname{tr}\left(\mathbf{L} \mathbf{V}^{T} \mathbf{B}^{T} \mathbf{B V}-\mathbf{2} \mathbf{L} \mathbf{X}^{T} \mathbf{B V}+2 \sigma^{2} \eta \mathbf{J}_{\mathbf{V}}^{T} \mathbf{V}\right)\right) \\
& \propto \exp \left(-\frac{1}{2 \sigma^{2}} \operatorname{tr}\left(\mathbf{L} \mathbf{V}^{T} \mathbf{B}^{T} \mathbf{B V}-\mathbf{2}\left(\mathbf{L} \mathbf{X}^{T} \mathbf{B}-\sigma^{2} \eta \mathbf{J}_{\mathbf{V}}^{T}\right) \mathbf{V}\right)\right)
\end{aligned}
$$

Applying the equation 3 it satisfies a matrix normal distribution with the mean matrix:

$$
\begin{aligned}
\mathbf{M} & =\left(\mathbf{B}^{T} \mathbf{B}\right)^{-T}\left(\mathbf{L} \mathbf{X}^{T} \mathbf{B}-\sigma^{2} \eta \mathbf{J}_{\mathbf{V}}{ }^{T}\right)^{T} \mathbf{L}^{-T} \\
& =\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1}\left(\mathbf{B}^{T} \mathbf{X} \mathbf{L}-\sigma^{2} \eta \mathbf{J}_{\mathbf{V}}\right) \mathbf{L}^{-1} \\
& =\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1}\left(\mathbf{B}^{T} \mathbf{X}-\sigma^{2} \eta \mathbf{J}_{\mathbf{V}} \mathbf{L}^{-1}\right)
\end{aligned}
$$

The two covariance matrices are $\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1}$ and $\mathbf{L}^{-1}$, and the variance $\sigma^{2}$ can be merged with either of them. Considering that the exponential distribution is constraint to be nonnegative, the posterior of $\mathbf{V}$ follows the truncated matrix normal distribution:

$$
\mathbf{V} \mid \mathbf{X} \sim \mathcal{T} \mathcal{M} \mathcal{N}\left(\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1}\left(\mathbf{B}^{T} \mathbf{X}-\sigma^{2} \eta \mathbf{J}_{\mathbf{V}} \mathbf{L}^{-1}\right),\left(\mathbf{B}^{T} \mathbf{B}\right)^{-1} \sigma^{2}, \mathbf{L}^{-1}\right) \quad(\mathbf{V} \geq 0)
$$

### 2.3 Posterior of $\sigma^{2}$

With $\sigma^{2} \sim$ Inverse-Gamma $\left(a_{2}, b_{2}\right)$, the pdf satisfies

$$
P\left(\sigma^{2}\right) \propto\left|\sigma^{2}\right|^{-a_{2}-1} \exp \left(-\frac{b_{2}}{\sigma^{2}}\right)
$$

According to equation 6 we have (let $\mathbf{E}=\mathbf{X}-\mathbf{B V}$ )

$$
P\left(\mathbf{X} \mid \sigma^{2}\right) \propto \frac{1}{\left|\sigma^{2}\right|^{n p / 2}} \exp \left(-\frac{1}{2 \sigma^{2}} \operatorname{tr}\left(\mathbf{L} \mathbf{E}^{T} \mathbf{E}\right)\right)
$$

Then

$$
\begin{aligned}
P\left(\sigma^{2} \mid \mathbf{X}\right) & \propto P\left(\mathbf{X} \mid \sigma^{2}\right) P\left(\sigma^{2}\right) \\
& \propto\left|\sigma^{2}\right|^{-n p / 2-a_{2}-1} \exp \left(\frac{1}{\sigma^{2}}\left(-b_{2}-\frac{\operatorname{tr}\left(\mathbf{L E}^{T} \mathbf{E}\right)}{2}\right)\right)
\end{aligned}
$$

![](https://cdn.mathpix.com/cropped/522e954c-7ef5-4851-8a69-29d8ad2df089-23.jpg?height=798&width=1645&top_left_y=244&top_left_x=210)
Figure 6: We make an experiment to compare rectified matrix normal (RMN) and truncated matrix normal (TMN) distribution. We use the MOB data and plot the probability histograms of the cell type proportions at the same spot as in the main text, still calculated with 2000 samples. We can observe that when the mean values are far from zero (like GC and M/TC), the samples of RMN are really close to those of TMN. When the mean values are close to zero (like EPL-IN, OSs and PGC), the difference between RMN and TMN is getting obvious. However, sampling from RMN is much faster than TMN.

So the posterior of $\sigma^{2}$ follows the inverse Gamma distribution

$$
\left.\sigma^{2} \mid \mathbf{X} \sim \text { Inverse-Gamma }\left(a_{2}+n p / 2, \quad b_{2}+\frac{\operatorname{tr}\left(\mathbf{L E}^{T} \mathbf{E}\right)}{2}\right)\right)
$$

### 2.4 Posterior of $\eta$

With $\eta \sim \operatorname{Gamma}\left(a_{3}, b_{3}\right)$, the pdf satisfies

$$
P(\eta) \propto \eta^{a_{3}-1} \exp \left(-\frac{\eta}{b_{3}}\right)
$$

According to equation 8 we have

$$
\begin{gathered}
P(\mathbf{V} \mid \eta) \propto \eta^{c p} \exp \left(-\eta \operatorname{tr}\left(\mathbf{J}_{\mathbf{V}}{ }^{T} \mathbf{V}\right)\right) \\
P(\eta \mid \mathbf{V}) \propto P(\mathbf{V} \mid \eta) P(\eta) \propto \eta^{c p+a_{3}-1} \exp \left(-\eta\left(\operatorname{tr}\left(\mathbf{J}_{\mathbf{V}}{ }^{T} \mathbf{V}\right)+\frac{1}{b_{3}}\right)\right)
\end{gathered}
$$

So the posterior of $\eta$ follows the Gamma distribution

$$
\eta \left\lvert\, \mathbf{V} \sim \operatorname{Gamma}\left(a_{3}+c p, \quad 1 /\left(\operatorname{tr}\left(\mathbf{J}_{\mathbf{V}}{ }^{T} \mathbf{V}\right)+\frac{1}{b_{3}}\right)\right)\right.
$$

## 3 Sampling from nonnegative multivariate normal distribution

We sample from truncated multivariate normal distribution (TMVN) based on the method proposed in 41 in which a Gibbs sampler is built. To improve the efficiency, we only consider the case that the lower bound is zero and the upper bound is infinity. Here we briefly describe the algorithm and how we improve it. Given an n -dimensional multivariate normal distribution truncated in $[0, \infty)$ :

$$
\mathbf{w} \sim \mathcal{T} \mathcal{M} \mathcal{V} \mathcal{N}(\mu, \boldsymbol{\Sigma}) \quad \mathbf{w} \geq \mathbf{0}
$$

whose covariance $\boldsymbol{\Sigma}$ is positive-definite, one can find the Choleskey decomposition $\boldsymbol{\Sigma}=\mathbf{L L}^{T}$ and introduce the transformation $\mathbf{x}=\mathbf{L}^{-1}(\mathbf{w}-\mu)$. The new random variable $\mathbf{x}$ follows

$$
\mathbf{x} \sim \mathcal{T} \mathcal{M} \mathcal{V} \mathcal{N}(\mathbf{0}, \mathbf{I}) \quad \mathbf{L} \mathbf{x} \geq-\mu
$$

Then it's proved in 41 that the ith variate of $\mathbf{x}$ follows the conditional distribution:

$$
\mathbf{x}_{\mathbf{i}} \mid \mathbf{x}_{-\mathbf{i}} \sim \mathcal{T} \mathcal{N}(0,1) \quad \mathbf{L}_{\mathbf{i}} \mathbf{x}_{\mathbf{i}} \geq-\mu-\mathbf{L}_{-\mathbf{i}} \mathbf{x}_{-\mathbf{i}}
$$

where $\mathbf{x}_{-\mathbf{i}}=\left(\mathbf{x}_{\mathbf{1}}, \cdots, \mathbf{x}_{\mathbf{i}-\mathbf{1}}, \mathbf{x}_{\mathbf{i}+\mathbf{1}}, \cdots, \mathbf{x}_{\mathbf{n}}\right)$ represents the ( $n-1$ ) $\times 1$ vector by removing the ith entry of $\mathbf{x}$, and $\mathbf{L}_{-\mathbf{i}}$ represents the $n \times(n-1)$ matrix by removing the ith column of $\mathbf{L}$. Sampling from the truncated univariate normal distribution in Equation 20 is not discussed here. The original algorithm can be summerized as Algorithm 1

```
Algorithm 1
    Given: \(\mu, \boldsymbol{\Sigma}\), sample size \(T\), burn-in period \(T_{b}\)
    Initialize: choose \(\mathbf{x}^{0}\), compute \(\boldsymbol{\Sigma}=\mathbf{L L}^{T}\)
    for \(t=1\) to \(T\) do
        for \(i=1\) to \(n\) do
            Define \(\mathrm{x}_{-\mathrm{i}}^{\mathrm{t}}=\left(\mathrm{x}_{1}^{\mathrm{t}}, \cdots, \mathrm{x}_{\mathrm{i}-1}^{\mathrm{t}}, \mathrm{x}_{\mathrm{i}+1}^{\mathrm{t}-1}, \cdots, \mathrm{x}_{\mathrm{n}}^{\mathrm{t}-1}\right)\)
            Sample from \(\mathbf{x}_{\mathbf{i}}^{\mathbf{t}} \mid \mathbf{x}_{-\mathbf{i}}^{\mathbf{t}} \sim \mathcal{T} \mathcal{N}(0,1) \quad \mathbf{L}_{\mathbf{i}} \mathbf{x}_{\mathbf{i}}^{\mathbf{t}} \geq-\mu-\mathbf{L}_{-\mathbf{i}} \mathbf{x}_{-\mathbf{i}}^{\mathbf{t}}\)
            Update \(\mathbf{x}^{\mathbf{t}}=\left(\mathbf{x}_{\mathbf{1}}^{\mathbf{t}}, \cdots, \mathbf{x}_{\mathbf{i}}^{\mathbf{t}}, \mathbf{x}_{\mathbf{i}+\mathbf{1}}^{\mathbf{t}-\mathbf{1}}, \cdots, \mathbf{x}_{\mathbf{n}}^{\mathbf{t}-\mathbf{1}}\right)\)
        end for
        \(\mathbf{w}^{\mathbf{t}}=\mathbf{L} \mathbf{x}^{\mathbf{t}}+\mu\)
    end for
    Discard \(\mathbf{w}^{1} \cdots \mathbf{w}^{T_{b}}\)
```

We improve the algorithm by replacing the matrix-vector multiplication in each step with scalar-vector multiplication. Therefore the larger the sampling dimension is, the more computation time we can save (Fig.5b). The improved sampling algorithm is summarized in Algorithm 2.

```
Algorithm 2
    Given: \(\mu, \boldsymbol{\Sigma}\), sample size \(T\), burn-in period \(T_{b}\)
    Initialize: choose \(\mathbf{x}^{0}\), compute \(\boldsymbol{\Sigma}=\mathbf{L L}^{T}, \mathbf{z}=\mathbf{L} \mathbf{x}^{0}\)
    for \(t=1\) to \(T\) do
        for \(i=1\) to \(n\) do
            Define \(\mathrm{x}_{-\mathrm{i}}^{\mathrm{t}}=\left(\mathrm{x}_{1}^{\mathrm{t}}, \cdots, \mathrm{x}_{\mathrm{i}-1}^{\mathrm{t}}, \mathrm{x}_{\mathrm{i}+1}^{\mathrm{t}-1}, \cdots, \mathrm{x}_{\mathrm{n}}^{\mathrm{t}-1}\right)\)
            Sample from \(\mathbf{x}_{\mathbf{i}}^{\mathbf{t}} \mid \mathbf{x}_{-\mathbf{i}}^{\mathbf{t}} \sim \mathcal{T} \mathcal{N}(0,1) \quad \mathbf{L}_{\mathbf{i}} \mathbf{x}_{\mathbf{i}}^{\mathbf{t}} \geq-\mu-\mathbf{z}+\mathbf{L}_{\mathbf{i}} \mathbf{x}_{\mathbf{i}}^{\mathbf{t}-\mathbf{1}}\)
            Update \(\mathbf{x}^{\mathbf{t}}=\left(\mathbf{x}_{1}^{\mathbf{t}}, \cdots, \mathbf{x}_{\mathbf{i}}^{\mathbf{t}}, \mathbf{x}_{\mathbf{i}+\mathbf{1}}^{\mathbf{t}-\mathbf{1}}, \cdots, \mathbf{x}_{\mathbf{n}}^{\mathbf{t}-\mathbf{1}}\right)\)
            Update \(\mathbf{z}=\mathbf{z}+\mathbf{L}_{\mathbf{i}}\left(\mathbf{x}_{\mathbf{i}}^{\mathbf{t}}-\mathbf{x}_{\mathbf{i}}^{\mathbf{t}-\mathbf{1}}\right)\)
        end for
        \(\mathbf{w}^{\mathbf{t}}=\mathbf{L} \mathbf{x}^{\mathbf{t}}+\mu\)
    end for
    Discard \(\mathbf{x}^{0} \cdots \mathbf{x}^{T_{b-1}}\)
```


## 4 Details of datasets

| Dataset | \# genes | \# spots | H\&E image | scRNA-seq data |
| :--- | :--- | :--- | :--- | :--- |
| MOB (Replicate 12) [37] | 16034 | 282 | Yes | GSE121891 [38] |
| Human PDAC-A-1 (GSM3036911) [35] | 19738 | 428 | Yes | GSE111672 (PDAC-A) [35] |
| Human PDAC-B (GSM3405534) [35] | 19738 | 224 | No | GSE111672 (PDAC-B) [35] |
| Mouse brain cortex (seqFISH+) [39] | 10000 | 524 | No | GSE102827 [40] |

Table 1: Spatial transcriptomics data we used in our studies.

| Cell type | \# |
| :--- | :--- |
| All | 12801 |
| EPL-IN (external plexiform layer interneurons) | 161 |
| GC (granule cells) | 8614 |
| MT-C (mitral and tufted cells) | 1133 |
| OSNs (olfactory sensory neurons) | 1200 |
| PGC (periglomerular cells) | 1693 |

Table 2: Detailed cell type information of the scRNA-seq dataset GSE121891.

| Cell type | \# |
| :--- | :--- |
| All | 48266 |
| Astrocytes | 7039 |
| Endothelial cells | 4071 |
| Excitatory cells layer 2\&3 | 2963 |
| Excitatory cells layer 4 | 3198 |
| Excitatory cells layer 5 | 3793 |
| Excitatory cells layer 6 | 3276 |
| Excitatory cells | 1057 |
| Interneurons | 936 |
| Macrophage | 537 |
| Microglia | 10158 |
| Mural | 782 |
| Oligodendrocytes | 10456 |

Table 3: Detailed cell type information of the scRNA-seq dataset GSE102827.

| Cell type | \# |
| :--- | :--- |
| All | 1926 |
| Acinar cells | 13 |
| Cancer clone A | 126 |
| Cancer clone B | 170 |
| Ductal - APOL high/hypoxic | 215 |
| Ductal - CRISP3 high/centroacinar like | 529 |
| Ductal - MHC Class II | 287 |
| Ductal - terminal ductal like | 350 |
| Endocrine cells | 3 |
| Endothelial cells | 11 |
| Fibroblasts | 5 |
| Macrophages A | 21 |
| Macrophages B | 19 |
| Mast cells | 14 |
| mDCs A | 12 |
| mDCs B | 33 |
| Monocytes | 18 |
| pDCs | 13 |
| RBCs | 15 |
| T cells \& NK cells | 40 |
| Tuft cells | 32 |

Table 4: Detailed cell type information of the scRNA-seq dataset GSE111672 (PDAC-A).

| Cell type | \# |
| :--- | :--- |
| All | 1733 |
| Acinar cells | 6 |
| Cancer clone A | 339 |
| Ductal - CRISP3 high/centroacinar like | 152 |
| Ductal - MHC Class II | 211 |
| Ductal - terminal ductal like | 736 |
| Endocrine cells | 13 |
| Endothelial cells | 159 |
| Macrophages | 9 |
| Mast cells | 13 |
| mDCs | 35 |
| Monocytes | 20 |
| RBCs | 3 |
| Tuft cells | 37 |

Table 5: Detailed cell type information of the scRNA-seq dataset GSE111672 (PDAC-B).

