#' Query GO Terms for Disease Associations
#'
#' This function checks whether a set of GO terms are associated with a given disease
#' using a user-specified CSV database file.
#'
#' @param goterms A character vector of GO term IDs (e.g., "GO:0008150").
#' @param disease A character string for the disease name or keyword to search.
#' @param database Path to the CSV file containing the GO-disease associations.
#' @param go_column Name of the column in the CSV file that contains GO term IDs. Default is "GOID".
#' @param disease_column Name of the column in the CSV file that contains disease names. Default is "DiseaseName".
#'
#' @return A data frame with two columns: GO term and a logical indicating disease association.
#' @export
#'
#' @examples
#' \dontrun{
#' query_GO_disease(c("GO:0008150", "GO:0003674"), "cancer", "ctd_data.csv")
#' }
query_GO_disease <- function(goterms, disease, database, go_column = "GOID", disease_column = "DiseaseName") {
  # Ensure required package
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }

  # Load data
  databaseC <- read.csv(database, stringsAsFactors = FALSE)

  # Filter for disease of interest and select GO column
  filtered <- databaseC |>
    dplyr::filter(grepl(disease, .data[[disease_column]], ignore.case = TRUE)) |>
    dplyr::pull(.data[[go_column]])

  # Check presence of GO terms
  result <- data.frame(
    GO = goterms,
    disease_association = goterms %in% filtered,
    stringsAsFactors = FALSE
  )

  return(result)
}
