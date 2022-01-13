
#---------------- About ---------------
# Author: Ariane Duverdier
# Project: AD systematic literature review

# This script imports the data downloaded from the SCOPUS database on the application of multivariate statistics (MS),
# machine learning & artificial intelligence (ML&AI), and Bayesian statistics (BS).
# Data consists of the following important field tags:
# AU = authors
# TI = document title
# SO = publication name (or source)
# DT = document type
# ID = keywords associated by SCOP or ISI database
# DE = authors' keywords
# AB = abstract
# UT = unique article identifier

# Functions
# import_scopus_data() : imports data
# remove_miscellaneous_articles(df) : helper function to clean data

#-----------------------------------------

library("bibliometrix")


remove_miscellaneous_articles <- function(df){
  # Function that removes from the dataframe, articles without an abstract, conference reviews, and erratums.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package
  
  # Output:
  #   df: cleaned dataframe
  
  df <- df[df$AB != "[NO ABSTRACT AVAILABLE]",]
  df <- df[df$DT != "CONFERENCE REVIEW",]
  df <- df[df$DT != "ERRATUM",]
  return(df)
}

import_scopus_data <- function(){
  # Function that imports the SCOPUS data into dataframes.
  
  # Output:
  #   MS_eczema_df, MLAI_eczema_df, BS_eczema_df : dataframes containing SCOPUS document information for MS, ML&AI, and BS methodologies for eczema
  #   MS_AD_df, MLAI_AD_df, BS_AD_df : dataframes containing SCOPUS document information for MS, ML&AI, and BS methodologies for AD
  #   MS_df, MLAI_df, BS_df : dataframes containing SCOPUS document information for MS, ML&AI, and BS methodologies for both AD and eczema
  #   full_df : dataframe containing SCOPUS document information for the full collection
  #   AD_df, eczema_df : dataframes containing SCOPUS document information for documents retrieved using the AD and eczema terms
  
  # A. Multivariate statistics
  MS_eczema_df <- convert2df(file = "data/MS_eczema.csv", dbsource = "scopus", format = "csv")
  MS_AD_df <- convert2df(file = "data/MS_AD.csv", dbsource = "scopus", format = "csv")
  MS_eczema_df <- remove_miscellaneous_articles(MS_eczema_df)
  MS_AD_df <- remove_miscellaneous_articles(MS_AD_df)
  
  # B. Machine learning & artificial intelligence
  MLAI_eczema_df <- convert2df(file = "data/MLAI_eczema.csv", dbsource = "scopus", format = "csv")
  MLAI_AD_df <- convert2df(file = "data/MLAI_AD.csv", dbsource = "scopus", format = "csv")  
  MLAI_eczema_df <- remove_miscellaneous_articles(MLAI_eczema_df)
  MLAI_AD_df <- remove_miscellaneous_articles(MLAI_AD_df)
  
  # C. Bayesian statistics
  BS_eczema_df <- convert2df(file = "data/BS_eczema.csv", dbsource = "scopus", format = "csv")
  BS_AD_df <- convert2df(file = "data/BS_AD.csv", dbsource = "scopus", format = "csv")
  BS_eczema_df <- remove_miscellaneous_articles(BS_eczema_df)
  BS_AD_df <- remove_miscellaneous_articles(BS_AD_df)
  # Remove Arora paper, it compares its method to Bayes methods but is actually deep learning framework
  BS_eczema_df <- BS_eczema_df[BS_eczema_df$UT != "2-S2.0-85077700986",]

  # D. Combine all datasets together to form a full collection
  full_df <- dplyr::bind_rows(MS_eczema_df,MS_AD_df,BS_AD_df,BS_eczema_df,MLAI_AD_df,MLAI_eczema_df)
  full_df <- dplyr::distinct(full_df, UT, .keep_all = TRUE)
  
  # E. Combine AD and eczema for their respective method
  BS_df <- dplyr::bind_rows(BS_AD_df,BS_eczema_df)
  BS_df <- dplyr::distinct(BS_df, UT, .keep_all = TRUE) 
  MS_df <- dplyr::bind_rows(MS_AD_df,MS_eczema_df)
  MS_df <- dplyr::distinct(MS_df, UT, .keep_all = TRUE) 
  MLAI_df <- dplyr::bind_rows(MLAI_AD_df,MLAI_eczema_df)
  MLAI_df <- dplyr::distinct(MLAI_df, UT, .keep_all = TRUE) 
  
  # F. Combine methods to create an AD and eczema df
  AD_df <- dplyr::bind_rows(BS_AD_df,MS_AD_df, MLAI_AD_df)
  AD_df <- dplyr::distinct(AD_df, UT, .keep_all = TRUE) 
  eczema_df <- dplyr::bind_rows(BS_eczema_df,MS_eczema_df, MLAI_eczema_df)
  eczema_df <- dplyr::distinct(eczema_df, UT, .keep_all = TRUE)
  
  return(list(MS_eczema_df,MS_AD_df,MLAI_AD_df,MLAI_eczema_df,BS_eczema_df,BS_AD_df,MS_df,MLAI_df,BS_df,full_df,eczema_df,AD_df))
}



