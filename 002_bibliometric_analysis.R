
#---------------- About ------------------
# Project: Systematic literature review of data-driven research on eczema and atopic dermatitis (AD).

# This script performs a bibliometric analysis on the collection of documents retrieved and imported from the SCOPUS database.
# Uses the bibliometrix R package to do so.
## Aria, M. & Cuccurullo, C. (2017) bibliometrix: An R-tool for comprehensive science mapping analysis, Journal of Informetrics, 11(4), pp 959-975, Elsevier.
## http:\\www.bibliometrix.org

# Functions
# descriptive_analysis(df) : displays key information about the collection of documents
# keyword_coocurences(df, n, edge_min, keyword_type) : produces a keyword co-occurence network
# thematic_map(df, n, freq_min, term_type, n_labels) : produces a thematic map
# venn_methodology(MS_df, MLAI_df, BS_df) and venn_term(AD_df, eczema_df) : produce venn diagrams by separating documents according to methodology (MS,ML&AI,BS) and term (AD, eczema) respectively
# freq_year_methodology(MS_df, BS_df, MLAI_df) and freq_year_term(eczema_df, AD_df, separation) : produce frequency line plots according to methodology (MS,ML&AI,BS) and term (AD, eczema) respectively

#-----------------------------------------


library("bibliometrix")
library(RVenn)
library(ggplot2)
library(dplyr)
library(tidyr)

source("001_data_import.R")



# ---------------- 1. Import SCOPUS data ----------------

# Use function defined in 001_data_import.R to import data.
l_df <- import_scopus_data()
MS_eczema_df <- l_df[[1]]
MS_AD_df <- l_df[[2]]
MLAI_AD_df <- l_df[[3]]
MLAI_eczema_df <- l_df[[4]]
BS_eczema_df <- l_df[[5]]
BS_AD_df <- l_df[[6]]
MS_df <- l_df[[7]]
MLAI_df <- l_df[[8]]
BS_df <- l_df[[9]]
full_df <- l_df[[10]]
eczema_df <- l_df[[11]]
AD_df <- l_df[[12]]




# --------------- 2. Descriptive analysis ---------------

descriptive_analysis <- function(df){
  # Function that displays key information about the dataset.
  # Provides some snapshots about the annual research development, the top “k” productive authors, papers, 
  # countries and most relevant keywords.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package
  
  results <- biblioAnalysis(df)
  
  # Print to the screen main information about the dataset.
  summary(results, k=10, pause=F, width=130)
  
  # Display key figures, hit <Return> to see next plot.
  plot(x=results, k=10, pause=T) 
}

descriptive_analysis(df = full_df)




# ----------------- 3. Co-word analysis ----------------- 

# --- A. Keyword co-occurence network
keyword_coocurences <- function(df, n, edge_min, keyword_type){
  # Function that runs produces a keyword co-occurence network.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package
  #   n : number of vertices to plot in co-occurence network
  #   edge_min: minimum number of co-occurences between vertices for edge to be plotted in the network
  #   keyword_type: "author_keywords" or "keywords" to represent if network uses the authors' keywords or the SCOPUS keywords, respectively
  
  # Output:
  #   cluster_df: dataframe containing the top n keywords and their assigned cluster
  
  NetMatrix <- biblioNetwork(df, analysis = "co-occurrences", network = keyword_type, sep = ";") 
  
  # Display co-occurence network.
  net=networkPlot(NetMatrix, normalize="association", n = n, Title = "Keyword Co-occurrences", cluster = "louvain", type = "auto", size.cex=TRUE, size=20, remove.multiple=T, edgesize = 10, labelsize=3,label.cex=TRUE,label.n=(n-1),edges.min=edge_min)

  # Produce cluster_df (keywords and their assigned cluster number).
  cluster_df <- data.frame(as.character(net$cluster_res$vertex),net$cluster_res$cluster)
  colnames(cluster_df) <- c("keyword", "cluster")
  return(cluster_df)
}

cluster_df <- keyword_coocurences(df = ull_df, n = 50, edge_min = 2, keyword_type = "author_keywords")


# --- B. Thematic map
thematic_map <- function(df, n, freq_min, term_type, n_labels){
  # Function that produces a thematic map using co-word analysis. 
  # Co-word analysis draws clusters of keywords. They are considered as themes, whose density and centrality can be used in classifying themes and mapping in a two-dimensional diagram.
  # Can analyze themes according to the quadrant in which they are placed: (1) upper-right quadrant: motor-themes; (2) lower-right quadrant: basic themes; (3) lower-left quadrant: emerging or disappearing themes; (4) upper-left quadrant: very specialized/niche themes.
  # Cobo, M. J., López-Herrera, A. G., Herrera-Viedma, E., & Herrera, F. (2011). An approach for detecting, quantifying, and visualizing the evolution of a research field: A practical application to the fuzzy sets theory field. Journal of Informetrics, 5(1), 146-166.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package
  #   n : number of terms to include in the analysis
  #   freq_min: minimum number of occurences of a term to be analyzed and plotted
  #   term_type: determines which terms in dataframe to analyze, "ID" for keywords associated by SCOPUS database, "DE" for authors' keywords, TI" for terms extracted from titles, and "AB" for terms extracted from abstracts
  #   n_labels : number of labels to associate to each cluster (top n_labels most important terms for each cluster)
  
  theme_map <- thematicMap(df, field = term_type, n = n, minfreq = freq_min, size = 0.2, repel = TRUE, n.labels = n_labels) 
  
  # Display thematic map.
  plot(theme_map$map)
}

thematic_map(df = full_df, n = 100, freq_min = 5, term_type = "DE", n_labels = 6)




# ---------------- 4. More Visualization ---------------- 

# --- A. Venn diagram by methodology
venn_methodology <- function(MS_df, MLAI_df, BS_df){
  # Function that produces a Venn diagram of the overlap in documents returned for the different
  # methodologies (MS, BS, ML&AI) using the RVenn library.
  
  # Input: dataframes created from upload of SCOPUS download using the bibliometrix package
  #   MS_df : dataframe containing the collection of MS documents
  #   MLAI_df : dataframe containing the collection of ML&AI documents
  #   BS_df : dataframe containing the collection of BS documents
  
  df_lists <- list('MS' = MS_df$UT,
                   'BS' = BS_df$UT,
                   'MLAI' = MLAI_df$UT)
  venn <- Venn(df_lists)
  
  # Plot venn diagram.
  ggvenn(venn)
}

venn_methodology(MS_df, MLAI_df, BS_df)


# --- B. Venn diagram by term
venn_term <- function(AD_df, eczema_df){
  # Function that produces a Venn diagram of the overlap in documents returned for the different
  # terms (eczema and AD) using the RVenn library.
  
  # Input: dataframes created from upload of SCOPUS download using the bibliometrix package
  #   AD_df : dataframe containing the collection of AD documents
  #   eczema_df : dataframe containing the collection of eczema documents
  
  df_lists <- list('AD' = AD_df$UT,
                   'eczema' = eczema_df$UT)
  venn <- Venn(df_lists)
  
  # Plot venn diagram.
  ggvenn(venn) 
}

venn_term(AD_df, eczema_df)
venn_term(MS_AD_df, MS_eczema_df)
venn_term(MLAI_AD_df, MLAI_eczema_df)
venn_term(BS_AD_df, BS_eczema_df)


plot_freq <- function(freq_year, colors){
  # Helper function that produces a line plot of frequency over years.
  
  # Input: 
  #   freq_df : dataframe containing the number of articles published for each year in each group
  #   colors : vector of colors for plotting the groups
  
  freq_year$Year<-as.integer(as.character(freq_year$Year))
  freq_year <- arrange(freq_year, freq_year$Year) 
  min_year <- min(freq_year$Year)
  max_year <- max(freq_year$Year)
  freq_year <- freq_year %>%
    complete(Type, Year = min_year:max_year, 
             fill = list(Articles = 0)) %>%
    as.data.frame()
  
  ggplot(data=freq_year, aes(x=Year, y=Articles, group = Type, color = Type)) + 
    geom_line() + geom_point()+
    scale_color_manual(values=colors)+
    theme_minimal()+
    ggtitle("Annual Publications")+
    ylab("Number of documents")+
    xlab("Years")+
    scale_x_continuous(limits=c(min_year,max_year), breaks=seq(1975,max_year,5))+
    scale_y_continuous(breaks=seq(min(freq_year$Articles), max(freq_year$Articles),5))+
    theme(axis.text=element_text(size=12),
           axis.title=element_text(size=14),
          legend.title = element_text(size=12),
          legend.text = element_text(size=12))
}


# --- C. Frequency plot by methodology
freq_year_methodology <- function(MS_df, BS_df, MLAI_df) {
  # Function that produces a plot of the number of documents published per year, separated by
  # methodology (MS, ML&AI, or BS). 
  
  # Input: dataframes created from upload of SCOPUS download using the bibliometrix package
  #   MS_df : dataframe containing the collection of MS documents
  #   MLAI_df : dataframe containing the collection of ML&AI documents
  #   BS_df : dataframe containing the collection of BS documents
  
  freq_year <- as.data.frame(table(MS_df$PY))
  freq_year <- cbind("Type" = "MS", freq_year)
  colnames(freq_year) <- c("Type","Year", "Articles")
  freq_year_BS <- as.data.frame(table(BS_df$PY))
  freq_year_BS <- cbind("Type" = "BS", freq_year_BS)
  colnames(freq_year_BS) <- c("Type","Year", "Articles")
  freq_year_MLAI <- as.data.frame(table(MLAI_df$PY))
  freq_year_MLAI <- cbind("Type" = "MLAI", freq_year_MLAI)
  colnames(freq_year_MLAI) <- c("Type","Year", "Articles")
  freq_year <- rbind(freq_year,freq_year_BS,freq_year_MLAI)
  colors <- c("lightgoldenrod1","cornflowerblue","hotpink2")
  
  # Remove 2021, since only 3 months in
  freq_year <- freq_year[freq_year$Year != 2021,]
  print(freq_year)

  # Plot frequency.
  plot_freq(freq_year, colors)
}

freq_year_methodology(MS_df, BS_df, MLAI_df)


# --- D. Frequency plot by term
freq_year_term <- function(eczema_df, AD_df, separation="both") {
  # Function that produces a plot of the number of documents published per year, separated by
  # term (eczema, AD, or both) that were used to retrieve them. 
  
  # Input: dataframes created from upload of SCOPUS download using the bibliometrix package
  #   AD_df : dataframe containing the collection of AD documents
  #   eczema_df : dataframe containing the collection of eczema documents
  #   separation : "both" means that publications will be grouped according to whether they were only
  #                 tagged with the term AD, eczema, or if both were used. Any other value will group documents by
  #                 the term (eczema or AD) that were used to retrieve them. Some documents may have been retrieved
  #                 in both searches (tagged with both terms) and these will be counted twice, once in eczema 
  #                 and once in AD.
  
  if(separation == "both"){
    common <- intersect(eczema_df$UT, AD_df$UT)  
    both_df <- eczema_df[eczema_df$UT %in% common,]
    freq_year_both <- as.data.frame(table(both_df$PY))
    freq_year_both <- cbind("Type" = "Both", freq_year_both)
    colnames(freq_year_both) <- c("Type","Year", "Articles")
    
    eczema_only_df <- dplyr::filter(eczema_df,!UT %in% common)
    freq_year <- as.data.frame(table(eczema_only_df$PY))
    freq_year <- cbind("Type" = "Only Eczema", freq_year)
    colnames(freq_year) <- c("Type","Year", "Articles")
    AD_only_df <- dplyr::filter(AD_df,!UT %in% common)
    freq_year_AD <- as.data.frame(table(AD_only_df$PY))
    freq_year_AD <- cbind("Type" = "Only AD", freq_year_AD)
    colnames(freq_year_AD) <- c("Type","Year", "Articles")
    
    freq_year_total <- as.data.frame(table(full_df$PY))
    freq_year_total <- cbind("Type" = "Total", freq_year_total)
    colnames(freq_year_total) <- c("Type","Year", "Articles")
    
    freq_year <- rbind(freq_year,freq_year_AD,freq_year_both,freq_year_total)
    colors <- c('plum',"#00A9FF","#7CAE00",'#A0A0A0')
  }
  else{
    freq_year <- as.data.frame(table(eczema_df$PY))
    freq_year <- cbind("Type" = "eczema", freq_year)
    colnames(freq_year) <- c("Type","Year", "Articles")
    freq_year_AD <- as.data.frame(table(AD_df$PY))
    freq_year_AD <- cbind("Type" = "AD", freq_year_AD)
    colnames(freq_year_AD) <- c("Type","Year", "Articles")
    freq_year <- rbind(freq_year,freq_year_AD)
    colors <- c("#00A9FF","#7CAE00")
  }

  # Remove 2021, since only 3 months in
  freq_year <- freq_year[freq_year$Year != 2021,]
  print(freq_year)
  
  # Plot frequency.
  plot_freq(freq_year, colors)
}

freq_year_term(eczema_df,AD_df,separation="both")



