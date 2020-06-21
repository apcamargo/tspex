---
title: 'tspex: a tissue-specificity calculator for gene expression data'
tags:
  - Python
  - bioinformatics
  - gene expression
  - transcriptomics
  - tissue-specificity
authors:
  - name: Antonio Pedro Camargo
    orcid: 0000-0003-3913-2484
    affiliation: 1
  - name: Adrielle Ayumi Vasconcelos
    orcid: 0000-0001-8145-4669
    affiliation: 1
  - name: Mateus Bernabe Fiamenghi
    orcid: 0000-0003-4535-8594
    affiliation: 1
  - name: Gonçalo Amarante Guimarães Pereira
    orcid: 0000-0003-4140-3482
    affiliation: 1
  - name: Marcelo Falsarella Carazzolle
    orcid: 0000-0002-5474-2830
    affiliation: 1
affiliations:
  - name: Department of Genetics, Evolution, Microbiology and Immunology, Institute of Biology, University of Campinas, Campinas, SP, 13083-862, Brazil.
    index: 1
date: 21 Jun 2020
bibliography: paper.bib
---

# Summary

When comparing gene expression data of different tissues it is often interesting to identify tissue-specific genes or transcripts. Even though there are several metrics to measure tissue-specificity, a user-friendly tool that facilitates this analysis is not available yet. We present tspex, a software that allows easy computation of a comprehensive set of different tissue-specificity metrics from gene expression data. tspex can be used through a web interface, command-line or the Python API. Its package version also provides visualization functions that facilitate inspection of results. The documentation and the source code of tspex are available at https://apcamargo.github.io/tspex/ and the web application can be accessed at https://tspex.lge.ibi.unicamp.br/.

# Introduction

High-throughput sequencing technologies have allowed the quantification of gene expression of several tissues from different species, generating huge amounts of data that can be analyzed for discovering novel expression profiles. Regarding expression patterns across tissues, genes can be positioned within a continuous scale that goes from housekeeping genes to tissue-specific genes. Housekeeping genes are responsible for critical cell functions and are ubiquitously expressed [@goldman2001housekeeping], while tissue-specific genes are expressed in a single or a small subset of tissues, suggesting specialized functions.

Detection of tissue-specific genes can be important, for instance, for gene discovery, evolutionary comparisons [@kryuchkova2016tissue; @duret2000determinants], drug target identification [@yang2018systematic], association with diseases [@amazon2018using], and cancer studies [@haigis2019tissue; @kim2017tissgdb]. In this context, there are several initiatives to sequence various tissues of an organism, such as the Genotype-Tissue Expression (GTEx) project [@lonsdale2013genotype], which offers expression data of 53 human tissues sampled from nearly 1,000 individuals, providing a foundation for the identification of new tissue-specific genes.

Although methods to quantify gene tissue-specificity have been extensively used in the literature, there is no available tool that allows easy measurement of tissue-specificity from gene expression data, forcing users to develop their own one-time use solutions. Herein, we present tspex, a tool that allows easy computation of twelve different tissue-specificity metrics from gene expression data. It provides visualization functions that facilitate exploration of results and can be used through different interfaces, including a web version. Additionally, tspex is also useful to measure expression specificity among distinct biological conditions, such as developmental stages, time points, genetic varieties, etc.

# Implementation

tspex is implemented as a Python package and it can be used locally through a Python Application Programming Interface (API), command-line interface or web version. Local installation of tspex is as easy as calling it with pip or conda and requires few dependencies. Refer to the tspex GitHub repository for the most up-to-date source code, dependency details and instructions. An open source web interface (Figure 1A) built with Flask and deployed using Docker containers can be accessed at https://tspex.lge.ibi.unicamp.br/ and its source code is available in https://github.com/apcamargo/tspex-webapp/.

tspex provides twelve distinct tissue-specificity metrics, which differ in their assumptions, scale and properties. Broadly, these metrics can be divided into two groups [@kryuchkova2017benchmark]: (1) general scoring metrics, that summarize in a single value how tissue-specific or ubiquitous is a gene across all tissues and (2) individualized scoring metrics that quantify how specific is the expression of each gene to each tissue.

The general scoring metrics provided by tspex are: Counts [@duret2000determinants], Tau [@yanai2004genome], Gini coefficient [@ceriani2012origins], Simpson index [@simpson1949measurement], Shannon entropy specificity [@schug2005promoter], ROKU specificity [@kadota2006roku], Specificity measure dispersion (SPM DPM) [@pan2012pagefinder], and Jensen-Shannon specificity dispersion (JSS DPM) [@cabili2011integrative]. As for individualized scoring metrics, tspex includes: Tissue-specificity index (TSI) [@julien2012mechanisms], Z-score [@vandenbon2009modeling], Specificity measure (SPM) [@xiao2010tisged], and Jensen-Shannon specificity (JSS) [@cabili2011integrative]. Each metric provides values that range within different scales, thus tspex includes an option to transform tissue-specificity values so that they fall within 0 (ubiquitous expression) and 1 (tissue-specific expression). The equations for all provided metrics as well as their transformations can be found in the documentation (https://apcamargo.github.io/tspex/metrics/).

As input, tspex requires an expression matrix (TSV, CSV or Excel formats) in any appropriate unit, such as TPM, FPKM or CPM. Optionally, tspex allows the expression values to be log-transformed before computation of tissue-specificity, which reduces the dependency between expression variance and expression level, improving the reliability of tissue-specificity measurements [@kryuchkova2017benchmark]. Internally, expression data and the tissue-specificity values are stored in a Python object and can be easily accessed for further investigation through the Python API.

Finally, the tspex package provides built-in functions for data visualization. Specifically, the user can plot histograms of tissue-specificity values (Figure 1B) and heatmaps of the expression of genes whose tissue-specificity is above a chosen value (Figure 1C). These visualizations allow quick inspection of the results and can be helpful for deciding threshold values.

# Results

To address the lack of tools for calculating tissue-specificity metrics from gene expression data, we developed tspex, a tool for easy computation of twelve different tissue-specificity scoring metrics. It is available to be used through three different interfaces: Python API, command-line, and a web version. Detailed tutorials on the usage of tspex interfaces can be found in the tspex documentation.

On the tspex website, calculation is easily done by simply uploading the input file, choosing the metric in the drop-down menu, and submitting it for calculation (Figure 1a). For guidance in choosing a tissue-specific metric for your gene expression data refer to this benchmark [@kryuchkova2017benchmark]. When it is done, the results can be downloaded and/or viewed it on the results page (only for data with up to 5000 genes). This can also be achieved by using tspex locally as a command line tool or in a Python API. Running tspex on a command line doesn’t require prior knowledge in Python and can be useful when adding it into an automated analysis pipelines or to run on multiple files at once. The Python API offers the advantages of tight integration with popular data analysis libraries, such as NumPy [@walt2011numpy], SciPy [@virtanen2020scipy], and pandas [@reback2020pandas], as well as visualization functions to create publication quality figures that aid the inspection of results.

Here, we demonstrate the features of tspex by running the tool in the Python API with real gene expression data from the Genotype-Tissue Expression (GTEx) project [@lonsdale2013genotype], which provides a large catalogue of gene expression across 54 human tissues. To showcase an example analysis, we used gene expression data from only five tissues: Bladder, Liver, Lung, Pancreas and Stomach. After removing genes that are not expressed in any of these tissues, expression values in transcripts per million (TPM) of 31872 genes were used as input for tspex. In order to obtain tissue-specificity values for these genes in the sampled human tissues, we calculated the general scoring metric Tau, which results in a single tissue-specificity score per gene. By running the visualization functions on the TissueSpecificity object, we can have an overview of the tissue-specificity results. The histogram plotting function can be used to verify the distribution of the tissue-specificity values, which is helpful for deciding thresholds for selecting genes (Fig 1b). Whereas, the heatmap plotting function is useful to visualize the amount of genes that are specific for each tissue above a given threshold (Fig 1c). In this heatmap, it is possible to see that lung has the largest number of tissue-specific genes among the five tissues.

The user can further explore and manipulate the tissue-specificity results obtained with tspex depending on what is being investigated. For example, searching for the top most tissue-specific gene for each tissue, generating lists of tissue-specific genes above a threshold for Gene Ontology term enrichment analysis, or filtering it to look at only certain class of genes. It is relevant to note that tissue-specificity values rely and will change according to the number of tissues being analysed and how different they are. Besides its application for finding tissue-specific expression, tspex is also useful to measure expression specificity among distinct biological conditions, such as developmental stages, time points, genetic varieties, etc.

![(a) tspex web interface home page. (b) Histogram and (c) heatmap created with the visualization functions available in the Python API.](figure.pdf)

# Conclusion

tspex is a software for calculating a variety of tissue-specificity metrics from gene expression data, addressing the lack of tools that perform this important task. It is entirely developed in Python to provide integration with an extensive library of data analysis packages. Finally, tspex can be used through three different interfaces, including a web version, providing solutions for different use cases.

# Acknowledgements

This work was supported by The São Paulo Research Foundation (FAPESP) grants #2013/08293-7 and #2016/23218-0, FAPESP fellowships to A.P.C (#2018/04240-0) and A.A.V (#2017/13015-7) and a fellowship from Coordenação de Aperfeiçoamento de Pessoal de Nível Superior – Brasil (CAPES) to M.B.F.

# References