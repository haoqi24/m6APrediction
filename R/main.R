
#' DNA sequence encoding
#'
#' Encode DNA sequences into feature matrix for m6A prediction
#'
#' @param dna_strings A vector of DNA sequences
#'
#' @return A data frame of encoded sequence features
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Multiple Sample Prediction for m6A Sites
#'
#' Predict m6A sites for multiple samples from a feature data frame.
#'
#' @param ml_fit A trained random forest model
#' @param feature_df A data frame containing all features for prediction
#' @param positive_threshold A numeric value between 0 and 1
#'
#' @return Input data frame appended with predicted probability and m6A status columns
#'
#' @examples
#' \dontrun{
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
#' }
#' @export
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  dna_enc_df <- dna_encoding(feature_df$DNA_5mer)
  pred_input <- cbind(feature_df[, !colnames(feature_df) %in% "DNA_5mer"], dna_enc_df)
  pred_prob <- predict(ml_fit, newdata = pred_input, type = "prob")[, "Positive"]
  pred_status <- ifelse(pred_prob >= positive_threshold, "Positive", "Negative")
  feature_df$predicted_m6A_prob <- pred_prob
  feature_df$predicted_m6A_status <- pred_status
  return(feature_df)
}

#' Single Sample Prediction for m6A Sites
#'
#' This function takes individual feature values and outputs m6A probability and status for a single sample.
#' It returns the predicted values as a named vector.
#'
#' @param ml_fit A trained random forest model
#' @param gc_content A numeric value for GC content
#' @param RNA_type A character string specifying RNA type
#' @param RNA_region A character string specifying RNA region
#' @param exon_length A numeric value for exon length
#' @param distance_to_junction A numeric value
#' @param evolutionary_conservation A numeric value
#' @param DNA_5mer A character string representing DNA 5-mer
#' @param positive_threshold A numeric value between 0 and 1; if the probability is greater than the threshold, treat it as positive
#'
#' @return A named vector with two elements: predicted_prob and predicted_m6A_status (either "Positive" or "Negative")
#'
#' @examples
#' \dontrun{
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5, evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA", positive_threshold = 0.5)
#' }
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  if(any(length(c(gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer)) == 0)){stop("Empty input")}
  input_df <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
    RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )
  dna_enc <- dna_encoding(DNA_5mer)
  pred_df <- cbind(input_df[, -7], dna_enc)
  pred_prob <- predict(ml_fit, newdata = pred_df, type = "prob")[ , "Positive"]
  pred_status <- ifelse(pred_prob >= positive_threshold, "Positive", "Negative")
  returned_vector <- c(
    "predicted_m6A_prob" = pred_prob,
    "predicted_m6A_status" = pred_status
  )
  return(returned_vector)
}
