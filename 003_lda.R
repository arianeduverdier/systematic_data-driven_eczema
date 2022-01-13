
#---------------- About ------------------
# Project: Systematic literature review of data-driven research on eczema and atopic dermatitis (AD).

# This script performs Latent Dirichlet Allocation (LDA) on the collection of documents retrieved and imported from the SCOPUS database.
# Uses the tm and topicmodels R packages to do so.

# Functions
# define_corpus(df) : defines the corpus for the collection of documents
# preprocess_corpus(review_corpus) : pre-processes the corpus
# dtm_from_corpus(review_corpus) : converts the corpus into a document-term matrix
# finding_k(review_dtm, fitting_method, max_k) : calculates statistical measures for a range of k values
# finding_k_plot(mod_results_df, measure) and finding_k_boxplot(df, measure) : plots these statistical measures
# run_lda(review_dtm, k, alpha) : computes the LDA model
# log_likelihood_iterations(topicModel) : plots the log-likelihood of the LDA model at each iteration in its Gibbs sampling inference
# model_diagnostics(topicModel, review_dtm) : produces summary information about the LDA model
# topic_word_cloud(topicModel,path_file), topic_proportions_decade(df, review_dtm, topicModel), topic_numbers_decade(df, review_dtm, topicModel), topic_numbers_years(df, review_dtm, topicModel, type_figure), top_words_topic(topicModel, n), and topic_dendogram(topicModel) : visualize the LDA model
# heat_map() : produces a heatmap of the use of MS, ML&AI, and BS methodologies in the LDA topics
# search_methodologies_MS(bigram_counts), search_methodologies_MLAI(bigram_counts), and search_methodologies_BS(bigram_counts) : return proxies of how many documents use specific methods within the major methodologies using bigram frequencies

#-----------------------------------------

library(tm)
library(stringr)
library(textmineR) 
library(topicmodels)
library(topicdoc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(wordcloud)
library(reshape2)
library(plyr)
library(tidytext)
library(pheatmap)

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

# LDA is performed on a collection of documents. Each document consists of a piece of text.
# For this project, the piece of text consists of the document's title, author keywords, and abstract.
# Combine the keywords, abstract, and title into a DOC column to be used in LDA.
full_df$DOC <- paste(full_df$TI,full_df$DE, full_df$AB)
AD_df$DOC <- paste(AD_df$TI,AD_df$DE, AD_df$AB)
eczema_df$DOC <- paste(eczema_df$TI,eczema_df$DE, eczema_df$AB)
BS_df$DOC <- paste(BS_df$TI,BS_df$DE, BS_df$AB)
MS_df$DOC <- paste(MS_df$TI,MS_df$DE, MS_df$AB)
MLAI_df$DOC <- paste(MLAI_df$TI,MLAI_df$DE, MLAI_df$AB)




# ------- 2. Pre-processing text using tm package -------

# --- First, define corpus.
define_corpus <- function(df){
  # Function that defines the corpus for the collection of documents.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package
  
  # Output:
  #   review_corpus: corpus of the collection of documents in the systematic review
  
  # End of abstract contains a copyright symbol and either authors or name of journal, not of interest so remove.
  df$DOC <- gsub("©.*","",df$DOC) 
  
  # Corpus consists of the text (title, author keywords, and abstract) of each document and its ID.
  corpus_df <- dplyr::select(df, "DOC","UT")
  colnames(corpus_df)<- c("text", "doc_id")
  review_corpus <- Corpus(DataframeSource(corpus_df))
  return(review_corpus)
}

review_corpus <- define_corpus(df = full_df)


# --- Second, pre-process the corpus.
preprocess_corpus <- function(review_corpus){
  # Function that pre-processes the corpus.
  
  # Input:
  #   review_corpus : corpus of the collection of documents in the systematic review
  
  # Output:
  #   review_corpus : pre-processed corpus of the collection of documents in the systematic review
  
  # Convert all characters to lower case.
  review_corpus = tm_map(review_corpus, content_transformer(tolower))
  
  # Remove numbers.
  review_corpus = tm_map(review_corpus, removeNumbers)
  
  # Remove punctuation marks and stopwords. 
  review_corpus = tm_map(review_corpus, removePunctuation)
  review_corpus = tm_map(review_corpus, removeWords, c("the", "and", stopwords("english")))
  
  # Remove extra whitespaces,
  review_corpus =  tm_map(review_corpus, stripWhitespace)
  
  # Remove text quotation marks, -, ~.
  removeSpecialChars <- function(x) gsub("“","",x)
  review_corpus <- tm_map(review_corpus, removeSpecialChars)
  removeSpecialChars <- function(x) gsub("”","",x)
  review_corpus <- tm_map(review_corpus, removeSpecialChars)
  removeSpecialChars <- function(x) gsub(" – "," ",x)
  review_corpus <- tm_map(review_corpus, removeSpecialChars)
  removeSpecialChars <- function(x) gsub(" ∼ "," ",x)
  review_corpus <- tm_map(review_corpus, removeSpecialChars)
  
  # Standardize common words.
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " diseases ", replacement = " disease ")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " factor ", replacement = " factors ")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = "clusters", replacement = "cluster")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = "clustering", replacement = "cluster")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = "identified", replacement = "identify")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " gene ", replacement = " genes ")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = "children", replacement = "childhood")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " cell ", replacement = " cells ")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " level ", replacement = "levels")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = "study", replacement = "studies")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " studied ", replacement = "studies")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = "associated", replacement = "association")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " image ", replacement = "images")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " patient ", replacement = "patients")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " adults ", replacement = "adult")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " expressed ", replacement = "expression")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " assessed ", replacement = "assessment")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " trial ", replacement = "trials")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " effective ", replacement = "effectiveness")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " score ", replacement = "scores")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = "patterns", replacement = "pattern")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " control ", replacement = "controls")
  review_corpus <- tm_map(review_corpus, content_transformer(gsub),
                          pattern = " predictive ", replacement = "prediction")
  return(review_corpus)
}

review_corpus <- preprocess_corpus(review_corpus)


# --- Third, convert corpus into a document-term matrix.
dtm_from_corpus <- function(review_corpus){
  # Function converts the corpus into a document-term matrix.
  
  # Input:
  #   review_corpus : pre-processed corpus of the collection of documents in the systematic review
  
  # Output:
  #   review_dtm : document-term matrix for the collection of documents
  
  # Convert corpus into a DocumentTermMatrix, only keep terms that occur >= 10 times (in other words, in >= 10 documents) 
  review_dtm <- DocumentTermMatrix(review_corpus, control = list(bounds = list(global = c(10, Inf))))
  print(paste0("Number of documents and terms in document-term matrix --before-- removing noisy common words: ",review_dtm$nrow," , ",review_dtm$ncol)) #  620 documents and 1257 terms 
  
  # Remove common words, they create 'noise' in the models.
  terms_remove <- c("can", "used","test","two","three", "significant","significantly","found", "showed","results","result","results","also","background","introduction","secondary", "objective","conclusion","conclusions", "method","methods","based","four","reported","respectively","including","using","included","among","one","performed","use","well","however","first","may","compared","different","aim","revealed","reveal","approach","important","affected","obtained","proposed","general","work","will","higher","low","lower","high","developed","new","items","increased","per","difference","several","provide","show","furthermore","good","useful","version","whether","second","although","copyright","investigated","possible","identify","better","paper","potential") #add terms here like c("can", "used")
  review_dtm <- review_dtm[, !(colnames(review_dtm) %in% terms_remove)]
  
  # Remove the documents that have a frequency of 0 for all terms remaining in matrix (terms that remain are non-'noisy' words).
  # LDA cannot handle empty documents.
  raw.sum <- apply(review_dtm,1,FUN=sum)
  review_dtm <- review_dtm[raw.sum!=0,]
  print(paste0("Number of documents and terms in document-term matrix --after--  removing noisy common words: ",review_dtm$nrow," , ",review_dtm$ncol)) # 619 documents and 1183 terms, only 1 document had a row sum of 0 and had to be deleted
  
  return(review_dtm)
}

review_dtm <- dtm_from_corpus(review_corpus)




# --------- 3. Run LDA using topicmodels package ---------

# --- A. Find the best k using statistical measures.
finding_k <- function(review_dtm, fitting_method, max_k){
  # Function that fits an LDA model for each k in range 2 to max_k. And calculates the 
  # coherence, topic exclusivity, and distance from corpus for each of the models.
  # Choice between VEM and Gibbs sampling for fitting of the models.
  
  # Input:
  #   fitting_method : "Gibbs" or "VEM" method used for fitting
  #   max_k : maximum number of k 
  
  # Output:
  #   mod_avg_coherence_df: dataframe of the average coherence for k = 2 to max_k
  #   mod_coherence_df: dataframe of the coherence of each of the topics in each model from k = 2 to max_k
  #   mod_topic_exclusivity: dataframe of the exclusivity of each of the topics in each model from k = 2 to max_k
  #   mod_dist_corpus_df: dataframe of the distance from corpus of each of the topics in each model from k = 2 to max_k
  
  # Finding the best K: Fit the model for several values of k, Plot the values, Pick the one where improvements are small, Similar to "elbow plot" in k-means clustering
  
  mod_avg_coherence = numeric(max_k-1)
  mod_coherence = vector(mode="list", length=max_k-1)
  mod_diagnostic = vector(mode="list", length=max_k-1)
  topics = seq(from = 2, to = max_k, by= 1)
  for (i in 2:max_k) {
    if(fitting_method == "VEM"){
      mod = LDA(review_dtm, k=i, method="VEM", control=list(seed = list(2003,5,63,100001,765), estimate.alpha = T, nstart = 5))
    }
    else{
      mod = LDA(review_dtm, k=i, method="Gibbs", control=list(iter=2000, seed=123)) 
    }
    mod_avg_coherence[i-1] = mean(topic_coherence(mod, review_dtm,top_n_tokens = 10))
    mod_coherence[[i-1]] = topic_coherence(mod, review_dtm,top_n_tokens = 10)
    mod_diagnostic[[i-1]] = topic_diagnostics(mod, review_dtm)
    tmResult <- topicmodels::posterior(mod)
    theta <- tmResult$topics
    print(paste0("Model fitted for k = ",i))
  }
  
  # Place results into dataframes for easy visualization.
  mod_avg_coherence_df <- data.frame(topics,mod_avg_coherence)
  mod_coherence_df <- data.frame(topics=numeric(0),coherence=numeric(0))
  for (i in 1:(length(mod_coherence))){
    list_coherence <- mod_coherence[[i]]
    for (entry in list_coherence){
      mod_coherence_df <- rbind(mod_coherence_df,c((i+1),entry))
    }
  }
  colnames(mod_coherence_df) <- c("topics","coherence")
  mod_topic_exclusivity <- data.frame(topics=numeric(0),exclusivity=numeric(0))
  for (i in 1:(length(mod_diagnostic))){
    list_exclusivity <- (mod_diagnostic[[i]])$topic_exclusivity
    for (entry in list_exclusivity){
      mod_topic_exclusivity <- rbind(mod_topic_exclusivity,c((i+1),entry))
    }
  }
  colnames(mod_topic_exclusivity) <- c("topics","exclusivity")
  mod_dist_corpus_df <- data.frame(topics=numeric(0),dist_from_corpus=numeric(0))
  for (i in 1:(length(mod_diagnostic))){
    list_dist_corpus <- (mod_diagnostic[[i]])$dist_from_corpus
    for (entry in list_dist_corpus){
      mod_dist_corpus_df <- rbind(mod_dist_corpus_df,c((i+1),entry))
    }
  }
  colnames(mod_dist_corpus_df) <- c("topics","dist_from_corpus")
  
  return(list(mod_avg_coherence_df,mod_coherence_df,mod_topic_exclusivity,mod_dist_corpus_df))
}

measures_list <- finding_k(review_dtm, fitting_method = "Gibbs", max_k = 16) 


finding_k_avg_coherence <- function(mod_avg_coherence_df){
  # Function that plots average coherence of models over a range of k values.
  
  # Input:
  #   mod_avg_coherence_df : output of finding_k, dataframe of the average coherence over the range of k values
  
  ggplot(data=mod_avg_coherence_df, aes(x=topics, y=mod_avg_coherence)) +
    geom_line() + geom_point()+
    scale_color_brewer(palette="Paired")+
    theme_minimal()+
    xlab("Number of Topics, K") + ylab("Average Coherence")+
    scale_x_continuous(breaks=seq(from = 2, to = max(mod_avg_coherence_df$topics), by= 1))
}

finding_k_avg_coherence(mod_avg_coherence_df = measures_list[[1]]) # mod_avg_coherence_df = measures_list[[1]]


finding_k_boxplot<- function(df, measure){
  # Function that produces boxplots of topic coherence, exclusivity, and distance from corpus for all the topics
  # in the models over a range of k values.
  
  # Input:
  #   df : output of finding_k, either mod_coherence_df, mod_topic_exclusivity, or mod_dist_corpus_df
  
  if(measure == "coherence"){
    y <- df$coherence
    ylab <- "Topic Coherence"
  }
  else if(measure == "exclusivity"){
    y <- df$exclusivity
    ylab <- "Topic Exclusivity"
  }
  else{ #assume measure == "distcorpus"
    y <- df$dist_from_corpus
    ylab <- "Topic Distance from Corpus"
  }
  
  df$topics <- as.factor(df$topics)
  ggplot(data = df, aes(x=topics, y=y)) +
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(position=position_jitter(width=.1, height=0)) + ylab(ylab) + xlab("Number of Topics, K")

}
finding_k_boxplot(df = measures_list[[2]], measure = "coherence") # mod_coherence_df = measures_list[[2]]
finding_k_boxplot(df = measures_list[[3]], measure = "exclusivity") # mod_topic_exclusivity = measures_list[[3]]
finding_k_boxplot(df = measures_list[[4]], measure = "distcorpus") # mod_dist_corpus_df = measures_list[[4]]



# --- B. Fit the LDA model.
run_lda <- function(review_dtm, k, alpha){
  # Function that computes the LDA model, inference via 2000 iterations of Gibbs sampling.
  
  # Input:
  #   review_dtm : document-term matrix for the collection of documents
  #   k : number of topics
  #   alpha : document-topic density factor (controls the number of topics expected in the document, low value of alpha
  #           implies fewer number of topics and a higher value implies a higher number topics)
  
  # Output:
  #   topicModel : LDA topic model
  
  topicModel <- LDA(review_dtm, k=k, method="Gibbs", control=list(iter = 2000, seed = 123, verbose = 25, keep = 1, alpha = alpha))
  
  return(topicModel)
}

topicModel <- run_lda(review_dtm, k = 8, alpha = 0.5) # default alpha value for LDA() function is 50 / k


# --- C. Plot topic statistics and diagnostics for the LDA model.
log_likelihood_iterations <- function(topicModel){
  # Function that produces a plot of the log-likelihood of the LDA model over the iterations of the Gibbs sampling.
  
  # Input:
  #   topicModel : LDA topic model
  
  iterations_x <- seq(from = 1, to = 2000, by= 1)
  log_lik_df <- data.frame(topicModel@logLiks,iterations_x)
  ggplot(data=log_lik_df, aes(x=iterations_x, y=topicModel.logLiks)) + 
    geom_line() + geom_point()+
    scale_color_brewer(palette="Paired")+
    theme_minimal()+
    xlab("Number of Iterations") + ylab("Log-Likelihood")
}

log_likelihood_iterations(topicModel)


model_diagnostics <- function(topicModel,review_dtm){
  # Function that produces summary information about the topicModel.
  # It prints the average probability of topics and basic diagnostics of the model.
  # It also produces a plot of these basic diagnostics.
  
  # Input:
  #   review_dtm : document-term matrix for the collection of documents
  #   topicModel : LDA topic model
  
  tmResult <- topicmodels::posterior(topicModel)
  beta <- tmResult$terms
  theta <- tmResult$topics 
  
  # Print average probability of topics over the collection of documents in decreasing order.
  topicProportions <- colSums(theta) / nDocs(review_dtm)  
  print("Topic Proportions in Decreasing Order:")
  print(sort(topicProportions, decreasing = TRUE))
  cat("\n\n")
  
  # Print basic diagnostics of the model.
  diag_df <- topic_diagnostics(topicModel, review_dtm)
  print("Basic Diagnostics of Model per Topic:s")
  print(diag_df)
  
  # Produce a plot of these basic diagnostics.
  colnames(diag_df) <- c("topic_num","Topic Size","Mean Token Length","Distance from Corpus","Distance between Token and Document Frequencies","Number of Documents Topic Appears","Topic Coherence","Topic Exclusivity")
  diag_df %>%
    gather(diagnostic, value, -topic_num) %>%
    ggplot(aes(x = as.factor(topic_num), y = value,
               fill = str_wrap(topic_num, 25))) +
    geom_bar(stat = "identity") +
    facet_wrap(~diagnostic, scales = "free") +
    labs(x = "Topic Number", y = "Diagnostic Value",
         fill = "Topic", title = "All Topic Model Diagnostics")
}

model_diagnostics(topicModel,review_dtm)




# -------- 4.  Visualization of LDA results/model --------

# --- A. Word cloud representation of topics
topic_word_cloud <- function(topicModel,path_file = paste0(getwd(),"/")){
  # Function that produces a word cloud representation of each topic in the LDA model.
  # Plots the top 40 terms in each model, size of the word is proportional to its probability (beta).
  
  # Input:
  #   topicModel : LDA topic model
  #   path_file : file path to directory want to store plots in, by default is current working directory
  
  mycolors <- brewer.pal(8, "Dark2")
  tmResult <- topicmodels::posterior(topicModel)
  for (topicNum in 1:topicModel@k){
    top40terms <- sort(tmResult$terms[topicNum,], decreasing=TRUE)[1:40]
    words <- names(top40terms)
    file_name <- paste0(path_file,"Wordcloud topic = ",topicNum,".png",sep = "")
    png(file_name, width=1000,height=1000) 
    wordcloud(words, top40terms, random.order = FALSE, color = mycolors,scale=c(10,1))
    dev.off()
  }
}

topic_word_cloud(topicModel)


# --- B. Mean topic proportions per decade
topic_proportions_decade <- function(df, review_dtm, topicModel){
  # Function that produces a bar plot of the mean topic proportions per decade of the LDA model.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package, used to build LDA model
  #   review_dtm : document-term matrix for the collection of documents
  #   topicModel : LDA topic model
  
  tmResult <- topicmodels::posterior(topicModel)
  theta <- tmResult$topics 
  
  # Publications are sorted into decades according to the year they were published.
  year_df <- df[df$UT %in% review_dtm$dimnames$Docs,]
  year_df$decade <- paste0(substr(year_df$PY, 0, 3), "0")
  
  # Get mean topic proportions per decade.
  topic_proportion_per_decade <- aggregate(theta, by = list(decade = year_df$decade), mean)
  
  # Reshape data frame.
  vizDataFrame <- melt(topic_proportion_per_decade, id.vars = "decade")
  colnames(vizDataFrame) <- c("decade", "topic","proportion")
  
  ggplot(vizDataFrame, aes(x=decade, y=proportion, fill=topic)) + 
    geom_bar(stat = "identity") +
    labs(x = "Decade", y = "Proportion",
         fill = "Topic")+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

topic_proportions_decade(df = full_df, review_dtm, topicModel)


# --- C. Topic publication numbers per decade
max_three <- function(row){
  # Helper function that returns the topic assignment for a document.
  
  # Input:
  #   row : row of theta dataframe of LDA model
  
  # Output:
  #   document_topic_assignment : dataframe that stores topic assignment for each document in collection
  
  # Top 3 most probable topics are assigned to each document (provided probabilities are above 0.1).
  # Topic assignment is stored in Choice1,Choice2,Choice3 rows of document_topic_assignment dataframe, 
  # will be NA if probability is below 0.1, otherwise will correspond to topic number (1,2,...,k).
  
  # Get the top 3 values in the row
  sorted_row <- sort(row, decreasing = T)[1:3] 
  
  # Set the values that are less than 0.1 to NA.
  sorted_row[sorted_row <0.1] <- NA 
  colnames_na <- names(which(colSums(is.na(sorted_row)) > 0))
  if(length(colnames_na)>0){ 
    for (colname_na in colnames_na){
      colnames(sorted_row)[colnames(sorted_row) == colname_na] <- NA  #change the column names of those with NA values to NA
    }
  }
  doc_id <- rownames(sorted_row)[1]
  
  # Update topic choices in the document_topic_assignment dataframe for that document.
  document_topic_assignment[,colnames(document_topic_assignment) == doc_id] <<- as.integer(colnames(sorted_row))
  return(document_topic_assignment)
}

topic_assignment <- function(topicModel){
  # Helper function that returns the topic assignment for all documents in the collection.
  
  # Input:
  #   topicModel : LDA topic model
  
  # Output:
  #   document_topic_assignment : dataframe that stores topic assignment for each document in collection
  
  # The estimated probabilities (theta) of observing each topic in each document can be used to assign each 
  # document its top three most probable topics. Provided the probabilities are greater than 0.1.
  tmResult <- topicmodels::posterior(topicModel)
  theta <- tmResult$topics 
  
  # Set up empty document_topic_assignment dataframe, to be filled by max_three() helper function.
  document_topic_assignment <<- data.frame(matrix(NA, nrow = nrow(theta), ncol = 3))
  rownames(document_topic_assignment) <<- rownames(theta)
  colnames(document_topic_assignment) <<- c("Choice1","Choice2","Choice3")
  document_topic_assignment <<- t(document_topic_assignment)
  
  # Use max_three() helper function to assign topics to each document.
  by(theta, seq_len(nrow(theta)), function(row) max_three(row))
  document_topic_assignment <<- as.data.frame(t(document_topic_assignment))
  return(document_topic_assignment)
}

topic_numbers_decade <- function(df, review_dtm, topicModel){
  # Function that produces a bar plot of the number of publications per topic of the LDA model grouped by decade.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package, used to build LDA model
  #   review_dtm : document-term matrix for the collection of documents
  #   topicModel : LDA topic model
  
  # Use topic_assignment() helper function to assign top 3 most probable topics for each document.
  # Provided the probabilities are greater than 0.1
  document_topic_assignment <- topic_assignment(topicModel)
  
  # Assign decades to each document according to year published.
  year_df <- full_df[full_df$UT %in% review_dtm$dimnames$Docs,]
  year_df$decade <- paste0(substr(year_df$PY, 0, 3), "0")
  document_topic_assignment$decade <- year_df$decade
  min(document_topic_assignment$decade) #1970
  max(document_topic_assignment$decade) #2020
  
  # Obtain publication numbers for each decade by adding up the topics for each document. 
  number_decades <- length(unique(document_topic_assignment$decade))
  topic_numbers_per_decade <- data.frame(matrix(NA, nrow = 8, ncol = number_decades))
  colnames(topic_numbers_per_decade) <- sort(unique(document_topic_assignment$decade))
  topic_numbers_per_decade <- melt(document_topic_assignment, id.vars="decade")
  colnames(topic_numbers_per_decade) <- c("decade", "variable","topic")
  topic_numbers_per_decade <- ddply(topic_numbers_per_decade,.(decade,topic),nrow)
  topic_numbers_per_decade <- na.omit(topic_numbers_per_decade)
  colnames(topic_numbers_per_decade) <- c("decade","topic","value")
  
  # Fill in the topics with missing decades to have a publication number of 0.
  freq_decade <- topic_numbers_per_decade %>%
    complete(topic, decade = unique(topic_numbers_per_decade$decade), 
             fill = list(value = 0)) %>%
    as.data.frame()
  
  topic_numbers_per_decade$topic <- as.factor(topic_numbers_per_decade$topic)
  ggplot(topic_numbers_per_decade, aes(x=decade, y=value, fill=topic)) + 
    geom_bar(stat = "identity") +
    labs(x = "Decade", y = "Number of Documents",
         fill = "Topic")+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

topic_numbers_decade(df = full_df, review_dtm, topicModel) 


# --- D. Topic publication numbers per year
topic_numbers_years <- function(df, review_dtm, topicModel, type_figure){
  # Function that produces a plot of the number of publications per topic of the LDA model over time.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package, used to build LDA model
  #   review_dtm : document-term matrix for the collection of documents
  #   topicModel : LDA topic model
  #   type_figure : "allin1" or "separate" to determine whether all topics numbers are on one plot or separate subplots 
  
  # Use topic_assignment() helper function to assign top 3 most probable topics for each document.
  # Provided the probabilities are greater than 0.1
  document_topic_assignment <- topic_assignment(topicModel)
  
  # Assign year published to each document.
  year_df <- full_df[full_df$UT %in% review_dtm$dimnames$Docs,]
  document_topic_assignment$year <- year_df$PY
  
  # Obtain publication numbers each year by adding up the topics for each document. 
  number_years <- length(unique(document_topic_assignment$year))
  topic_numbers_per_year <- data.frame(matrix(NA, nrow = topicModel@k, ncol = number_years))
  colnames(topic_numbers_per_year) <- sort(unique(document_topic_assignment$year))
  topic_numbers_per_year <- melt(document_topic_assignment, id.vars="year")
  colnames(topic_numbers_per_year) <- c("year", "variable","topic")
  topic_numbers_per_year <- ddply(topic_numbers_per_year,.(year,topic),nrow)
  topic_numbers_per_year <- na.omit(topic_numbers_per_year)
  colnames(topic_numbers_per_year) <- c("year","topic","value")
  
  min_year <- min(topic_numbers_per_year$year) #1973
  max_year <- max(topic_numbers_per_year$year) #2021
  
  # Fill in the topics with missing decades to have a publication number of 0.
  freq_year <- topic_numbers_per_year %>%
    complete(topic, year = min_year:max_year, 
             fill = list(value = 0)) %>%
    as.data.frame()
  freq_year$topic <- as.factor(freq_year$topic)
  
  if(type_figure=="allin1"){ 
    # Plot topic numbers as line plot all on one figure.
    ggplot(data=freq_year, aes(x=year, y=value, group = topic, color = topic)) + 
      geom_line() + geom_point()+
      labs(x = "Years", y = "Number of Documents",
           color = "Topic")+
      scale_x_continuous(breaks=seq(1970,2021,5))+
      scale_y_continuous(breaks=seq(min(freq_year$value), max(freq_year$value),10))
  }
  else{ #assume type_figure == "separate"
    # Plot topic numbers as separate subplots.
    freq_year %>%
      ggplot(aes(x = year, y = value, color = topic)) +
      geom_line() + geom_point()+
      facet_wrap(~topic, nrow = 2) + 
      labs(x = "Years", y = "Number of Documents", color = "Topic")
  }
}

topic_numbers_years(df = full_df, review_dtm, topicModel, type_figure="allin1")
topic_numbers_years(df = full_df, review_dtm, topicModel, type_figure="separate")


# --- E. Top terms (words) per topic
top_words_topic <- function(topicModel, n = 10){
  # Function that produces a plot of the top 10 terms (words) per topic in the LDA model.
  
  # Input:
  #   topicModel : LDA topic model
  #   n : number of top terms in each topic to plot, default is 10
  
  x_lab <- paste0("Top ",n," Terms",sep="")
  
  # Extract the per-topic-per-word probabilities (beta), from the model.
  text_topics <- tidy(topicModel, matrix = "beta") 
  
  text_top_terms <- text_topics %>%
    group_by(topic) %>%
    top_n(n, beta) %>%
    ungroup() %>%
    arrange(topic, -beta)
  text_top_terms$topic <- as.factor(text_top_terms$topic)
  text_top_terms %>%
    mutate(term = reorder(term, beta)) %>%
    ggplot(aes(term, beta, fill = topic)) +
    geom_col(show.legend = T) +
    facet_wrap(~ topic, scales = "free") +
    labs(x = x_lab, y = "Beta", fill = "Topic")+
    coord_flip()
}

top_words_topic(topicModel, n = 10)


# --- F. Topic dendogram
topic_dendogram <- function(topicModel){
  # Function that plots a topic dendogram for the LDA model using Hellinger distance 
  # (distance between 2 probability vectors) to decide if the topics are closely related.
  
  # Input:
  #   topicModel : LDA topic model
  
  tmResult <- topicmodels::posterior(topicModel)
  beta <- tmResult$terms
  topic_linguistic_dist <- CalcHellingerDist(beta)
  topic_hclust <- hclust(as.dist(topic_linguistic_dist), "ward.D")
  plot(topic_hclust, main = "Topic Dendogram",xlab = "Topics",ylab = "Hellinger Distance", sub = "")
}

topic_dendogram(topicModel)


# --- G. Heatmap of distribution of the application of MS, ML&AI, and BS methodologies in the LDA topics
heat_map <- function(topicModel, topic_assignment_type){
  # Function that produces a heatmap of the distribution of the application of MS, ML&AI, and BS methodologies in the LDA topics.
  
  # Input:
  #   topicModel : LDA topic model
  #   topic_assignment_type : "top3" or "top1" to determine if we consider only the top 1 or 3 topics assigned to each document
  
  # Use topic_assignment() helper function to assign top 3 most probable topics for each document.
  # Provided the probabilities are greater than 0.1
  document_topic_assignment <- topic_assignment(topicModel)
  document_topic_assignment$UT <- rownames(document_topic_assignment)
  
  if(topic_assignment_type == "top3"){
    # Use all top 3 topic assignments.
    document_topic_melted <- melt(document_topic_assignment, id.vars="UT")
  }
  else{ # assume topic_assignment_type == "top1"
    # Use only top 1 topic assignment.
    document_topic_melted <- dplyr::select(document_topic_assignment, "UT","Choice1")
    colnames(document_topic_melted) <- c("UT","value")
  }

  document_topic_assignment_MS <- document_topic_melted[document_topic_melted$UT %in% MS_df$UT,]
  document_topic_assignment_MLAI <- document_topic_melted[document_topic_melted$UT %in% MLAI_df$UT,]
  document_topic_assignment_BS <- document_topic_melted[document_topic_melted$UT %in% BS_df$UT,]
  
  # A. Multivariate statistics
  freq_MS <- as.data.frame(table(document_topic_assignment_MS$value)) #value represents topic
  colnames(freq_MS) <- c("Topic","Freq")
  freq_MS$Topic <- as.numeric(freq_MS$Topic)
  # Fill in all the missing topics with 0 documents with Freq 0.
  freq_MS <- freq_MS %>%
    complete(Topic = 1:8,
             fill = list(Freq = 0)) %>%
    as.data.frame()
  freq_MS$Topic <- as.factor(freq_MS$Topic)
  freq_MS <- cbind(freq_MS, Percentage = numeric(nrow(freq_MS)))
  freq_MS$Percentage <- (freq_MS$Freq/sum(freq_MS$Freq))*100

  # B. Machine learning & artificial intelligence
  freq_MLAI <- as.data.frame(table(document_topic_assignment_MLAI$value)) #value represents topic
  colnames(freq_MLAI) <- c("Topic","Freq")
  freq_MLAI$Topic <- as.numeric(freq_MLAI$Topic)
  # Fill in all the missing topics with 0 documents with Freq 0.
  freq_MLAI <- freq_MLAI %>%
    complete(Topic = 1:8,
             fill = list(Freq = 0)) %>%
    as.data.frame()
  freq_MLAI$Topic <- as.factor(freq_MLAI$Topic)
  freq_MLAI <- cbind(freq_MLAI, Percentage = numeric(nrow(freq_MLAI)))
  freq_MLAI$Percentage <- (freq_MLAI$Freq/sum(freq_MLAI$Freq))*100
  
  # Bayesian statistics.
  freq_BS <- as.data.frame(table(document_topic_assignment_BS$value)) #value represents topic
  colnames(freq_BS) <- c("Topic","Freq")
  freq_BS$Topic <- as.numeric(freq_BS$Topic)
  # Fill in all the missing topics with 0 documents with Freq 0.
  freq_BS <- freq_BS %>%
    complete(Topic = 1:8,
             fill = list(Freq = 0)) %>%
    as.data.frame()
  freq_BS$Topic <- as.factor(freq_BS$Topic)
  freq_BS <- cbind(freq_BS, Percentage = numeric(nrow(freq_BS)))
  freq_BS$Percentage <- (freq_BS$Freq/sum(freq_BS$Freq))*100
  
  
  freq_methodology <- data.frame(Percentage_MS <- freq_MS$Percentage, Percentage_MLAI <- freq_MLAI$Percentage, Percentage_BS <- freq_BS$Percentage)
  colnames(freq_methodology) <- c("MS","ML&AI","BS")
  
  pheatmap(freq_methodology, display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 15,color = colorRampPalette(c('white','lightblue'))(nrow(freq_methodology)*3),legend = T,angle_col = 0)

}
heat_map(topicModel, topic_assignment_type = "top1")
heat_map(topicModel, topic_assignment_type = "top3")




# ------ 5.  Distribution of methods in MS and ML&AI ------
bigram_freq <- function(df){
  # Helper function that produces a dataframe containing the bigram frequencies of the collection of documents.
  
  # Input:
  #   df : dataframe created from upload of SCOPUS download using the bibliometrix package
  
  # Output:
  #   bigram_counts : dataframe containing the bigrams and their frequencies within the collection of documents
  
  # Produce a dataframe that contains the text (title, keywords, abstract) of each document in the collection
  # and its ID.
  df$DOC <- gsub("©.*","",df$DOC) # excluding the (c)
  corpus_df <- dplyr::select(df, "DOC","UT")
  colnames(corpus_df)<- c("text", "doc_id")
  
  # Retrieve all bigrams (every consecutive pair of words) in the texts.
  corpus_bigrams <- corpus_df %>%
    unnest_tokens(bigram, text, token = "ngrams", n = 2)
  
  # Check don't have any duplicate bigrams within the same document (only count a bigram once in each document).
  corpus_bigrams <- corpus_bigrams %>% 
    distinct(doc_id, bigram, .keep_all = TRUE)
  
  bigrams_separated <- corpus_bigrams %>%
    separate(bigram, c("word1", "word2"), sep = " ")
  
  # For each bigram, count the number that appear in the collection of documents.
  bigram_counts <- bigrams_separated %>% 
    dplyr::count(word1, word2, sort = TRUE)
  
  # Standardize analyses.
  word1 <- bigram_counts$word1
  word1 <- word1 %>% str_replace_all("analyses","analysis")
  bigram_counts$word1 <- word1
  word2 <- bigram_counts$word2
  word2 <- word2 %>% str_replace_all("analyses","analysis")
  bigram_counts$word2 <- word2
  
  return(bigram_counts)
}
remove_na_sum <- function(obj){
  # Helper function to return the sum of an object.
  return(sum(obj[!is.na(obj)]))
}

search_methodologies_MS <- function(bigram_counts){
  # Function that prints the count of each bigram of interest, used to see in how many documents the bigram
  # is found. This function searchs for major multivariate statistics methods, this can be used as a proxy
  # of how many documents use the respective methods.
  
  # Input:
  #   bigram_counts : dataframe containing the bigrams and their frequencies within the collection of documents
  
  # factor analysis or factorial analysis
  print(paste0("factor analysis: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="factor"|bigram_counts$word1 == "factorial")&(bigram_counts$word2=="analysis"),]$n)))
  # cluster analysis or clustering
  print(paste0("cluster analysis: n = ",remove_na_sum(bigram_counts[((bigram_counts$word1=="cluster")&(bigram_counts$word2=="analysis"))|(bigram_counts$word1=="clustering"),]$n))) #only check for clustering in one word, not both or else will be double counted, for example if in phrase 'use clustering approach' will get both 'use clustering' and 'clustering approach'
  # markov models
  print(paste0("markov models: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="markov"),]$n))) #only check for clustering in one markov, not both or else will be double counted
  # principal component
  print(paste0("principal component analysis: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="principal")&(bigram_counts$word2=="component"|bigram_counts$word2=="components"),]$n)))
  # discriminant analysis
  print(paste0("discriminant analysis: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="discriminant")&(bigram_counts$word2=="analysis"),]$n)))
  # correspondence analysis
  print(paste0("correspondence analysis: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="correspondence")&(bigram_counts$word2=="analysis"),]$n)))
  # canonical correlation
  print(paste0("canonical correlation: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="canonical")&(bigram_counts$word2=="correlation"),]$n)))
  # latent class/transition models (only 1 is transition model, 0 transition models)
  print(paste0("latent class/transition model: n = ",remove_na_sum(bigram_counts[((bigram_counts$word1=="latent")&(bigram_counts$word2=="class"|bigram_counts$word2=="classes"))|((bigram_counts$word1=="transition")&(bigram_counts$word2=="model")),]$n)))
  # mixture model
  print(paste0("mixture model: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="mixture")&(bigram_counts$word2=="model"|bigram_counts$word2=="models"),]$n)))
  # latent variable
  print(paste0("latent variable model: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="latent")&(bigram_counts$word2=="variable"|bigram_counts$word2=="variables"),]$n)))
  # multidimensional scaling
  print(paste0("multidimensional scaling: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="multidimensional")&(bigram_counts$word2=="scaling"),]$n)))
  # profile regression
  print(paste0("profile regression: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="profile")&(bigram_counts$word2=="regression"),]$n)))
  # latent profile
  print(paste0("latent profile analysis: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="latent")&(bigram_counts$word2=="profile"),]$n)))
  # structural equation modeling (SEM)
  print(paste0("structural equation modeling (SEM): n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="structural")&(bigram_counts$word2=="equation"),]$n)))
  # logistic regression
  print(paste0("logistic regression: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="logistic")&(bigram_counts$word2=="regression"),]$n)))
}

search_methodologies_MLAI <- function(bigram_counts){
  # Function that prints the count of each bigram of interest, used to see in how many documents the bigram
  # is found. This function searchs for major machine learning & artificial intelligence methods, this can be used 
  # as a proxy of how many documents use the respective methods.
  
  # Input:
  #   bigram_counts : dataframe containing the bigrams and their frequencies within the collection of documents
  
  # neural network(s) (this will retrieve neural networks,artificial neural networks, and convolutional neural networks), then also need their acronyms ANN and CNN
  ANN <- remove_na_sum(bigram_counts[(bigram_counts$word1=="neural")&(bigram_counts$word2=="network"|bigram_counts$word2 == "networks"),]$n)
  ANN <- ANN+remove_na_sum(bigram_counts[(bigram_counts$word1=="cnn")|(bigram_counts$word1=="ann"),]$n) #only need to do one word, or else will be counted twice since for a phrase 'developed ann framework' you get the bigrams 'developed ann' and 'ann framework' so only count one
  ANN <- ANN - remove_na_sum(bigram_counts[(bigram_counts$word1=="network")&(bigram_counts$word2=="ann"),]$n) #remove "network ann", I assume these come from the longer phrase 'artificial neural network (ann)' and so are already counted
  ANN <- ANN - remove_na_sum(bigram_counts[(bigram_counts$word1=="network")&(bigram_counts$word2=="cnn"),]$n) #remove "network cnn", I assume these come from the longer phrase 'convolutional neural network (cnn)' and so are already counted
  print(paste0("artificial neural networks (including CNNs): n = ",ANN))
  # machine learning
  print(paste0("machine learning: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="machine")&(bigram_counts$word2=="learning"),]$n)))
  # artificial intelligence
  print(paste0("artificial intelligence: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="artificial")&(bigram_counts$word2=="intelligence"),]$n)))
  # deep learning
  print(paste0("deep learning: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="deep")&(bigram_counts$word2=="learning"),]$n)))
  # unsupervised learning
  print(paste0("unsupervised learning: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="unsupervised")&(bigram_counts$word2=="learning"),]$n)))
  # supervised learning
  print(paste0("supervised learning: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="supervised")&(bigram_counts$word2=="learning"),]$n)))
  # support vector machine
  SVM <- remove_na_sum(bigram_counts[(bigram_counts$word1=="support")&(bigram_counts$word2=="vector"),]$n)
  SVM <- SVM + remove_na_sum(bigram_counts[(bigram_counts$word1=="svm"),]$n) #as for ann and cnn above only need 1 word
  SVM <- SVM - remove_na_sum(bigram_counts[(bigram_counts$word1=="machine")&(bigram_counts$word2=="svm"),]$n) #remove "machine svm", I assume these come from the longer phrase 'support vector machine (svm)' and so are already counted
  print(paste0("support vector machine: n = ",SVM))
  # decision trees, classification trees, regression trees
  print(paste0("decision trees: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="decision"|bigram_counts$word1=="classification"|bigram_counts$word1=="regression")&(bigram_counts$word2=="tree"|bigram_counts$word2=="trees"),]$n))) 
  # random forest
  print(paste0("random forest: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="random")&(bigram_counts$word2=="forest"|bigram_counts$word2=="forests"),]$n)))
  # natural language processing
  NLP <- remove_na_sum(bigram_counts[((bigram_counts$word1=="natural")&(bigram_counts$word2=="language"))|((bigram_counts$word1=="nlp")),]$n) #as for ann, cnn, and svm above only need 1 word don't check both or else will count twice
  NLP <- NLP - remove_na_sum(bigram_counts[(bigram_counts$word1=="processing")&(bigram_counts$word2=="nlp"),]$n) #remove "processing nlp", I assume these come from the longer phrase 'natural language processing (nlp)' and so are already counted
  print(paste0("natural language processing: n = ",NLP))
  # reinforcement learning, reinforced learning
  print(paste0("reinforcement learning: n = ",remove_na_sum(bigram_counts[(bigram_counts$word1=="reinforced"|bigram_counts$word1=="reinforcement"|bigram_counts$word1=="regression")&(bigram_counts$word2=="learning"),]$n))) 
}

bigram_counts <- bigram_freq(df = full_df)
search_methodologies_MS(bigram_counts) 
search_methodologies_MLAI(bigram_counts) 



