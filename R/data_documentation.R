#' URLS to active databases
#'
#' @description A list containing the urls to active databases. Named by organism (eg 'kn99' or 's288cr64')
#'
#' @format A list of lists, organized by organism. Currently contains info for
#'     KN99 and S288C_R64 databases
#' \describe{
#'   \item{host}{host of the database server, eg ec2-18-224-181-136.us-east-2.compute.amazonaws.com}
#'   \item{db_name}{cryptococcus database name. something like kn99_database}
#'   \item{urls}{urls to all tables in database}
#'   ...
#' }
#' @source \url{https://rnaseq-databases-documentation.readthedocs.io/en/latest/}
"database_info"

#'
#' A named list containing a run number without a leading zero, eg 647, with the value being the same runnumber with
#' a leading 0, eg 0647.
#' @description this is best remedied in the database itself by forcing the column to be a string and adding the 0s
#'
"run_numbers_with_leading_zero"

#'
#' the 2016 grant summary represented as a dataframe
#'
#' @description also check the google sheet
#'
"grant_df"

#'
#' genes considered NOT to have 'high dispersion'
#'
#' @description union of old and new protocol
#'
"passing_genes_all"

#'
#' thresholds set for the novo+htseq pipeline
#'
#' @description list of metrics and their thresholds
#'
"kn99_novo_htseq_thresholds"

#'
#' QC status codes for the novo+htseq pipeline
#'
#' @description list of metrics and their statuses (powers of 2)
#'
"kn99_novo_htseq_status"
