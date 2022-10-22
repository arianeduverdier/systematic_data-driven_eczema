# Systematic literature review on data-driven research on atopic dermatitis and eczema
 
This repository contains the code developed for the article by [Duverdier et al. (2022), "Data-driven research on eczema: Systematic characterization of the field and recommendations for the future"](https://doi.org/10.1002/clt2.12170), published in Clinical and Translational Allergy. 
 
This project performs a systematic literature review on data-driven atopic dermatitis (AD) and eczema research. 
1. The literature search was conducted on the SCOPUS database, retrieving all documents that apply multivariate statistics (MS), machine learning and artificial intelligence (ML&AI), and/or Bayesian statistics (BS) methods to AD and eczema research. 
2. A bibliometric analysis was conducted on the corpus of documents to highlight the publication trends and conceptual knowledge structure of the field of research. 
3. Topic modelling, using the Latent Dirichlet Allocation (LDA) algorithm was applied on the corpus of documents to retrieve the key topics present within the literature.


The code is written in R version 4.0.4.


## Files
The datasets required for this project are found in the folder [data](data). These csv files contain the downloaded SCOPUS entries according to the methodology (MS, ML&AI, or BS) and the term (AD or eczema). The data can be loaded and pre-processed using the functions in [001_data_import.R](001_data_import.R).

The bibliometric analysis is performed in [002_bibliometric_analysis.R](002_bibliometric_analysis.R) using the [bibliometrix R package](https://www.bibliometrix.org) introduced in [Aria & Cuccurullo (2017)](https://doi.org/10.1016/j.joi.2017.08.007).

LDA is applied to the collection of documents in [003_lda.R](003_lda.R) using the [topicmodels R package](https://cran.r-project.org/web/packages/topicmodels/index.html) implementation. 

All functions and examples of how to run them are provided in the files.
	
## Libraries
The project is created with the following libraries and their versions:
* bibliometrix 3.0.4
* tm 0.7-8
* topicmodels 0.2-12
* topicdoc 0.1.0
* RVenn 1.1.0
* ggplot2 3.3.3
* dplyr 1.0.4
* tidyr 1.1.2
* plyr 1.8.6
* stringr 1.4.0
* textmineR 3.0.4
* tidytext 0.3.1
* wordcloud 2.6
* reshape2 1.4.4
* pheatmap 1.0.12

All packages used in this project can be installed using `install.packages("name_of_package")`.
