#' Extract TRN(s) from input
#'
#' @param x Character vector.
#' @param registry Character. Registry of TRN. Default to NULL for all registries.
#' @param collapse Character. When a collapse string is given, (by default ";"),
#' then the return vector will be collapsed to a single string with the chosen
#' separator string. If collapse is set to "none", the return is a vector with
#' all extracted TRNs instead.
#' @param clean Logical. By default (TRUE), the individual TRNs will also be
#' cleaned with `trn_clean`. Skip cleaning by setting to FALSE.
#'
#' @return Character. A vector of the same length as the input or
#'  of length 1, depending on the value of collapse.
#'  Returns NA if no trial registration numbers detected.
#'
#' @importFrom rlang .data
#' @export
#' @examples
#' which_trn("NCT00312962")
#' which_trn("NCT00312962 and euctr2020-001808-42")
#' which_trn("NCT00312962 and euctr2020-001808-42", collapse = "none")
#' which_trn("NCT00312962 and euctr2020-001808-42", registry = "EudraCT")
#' which_trn("NCT00312962 and euctr2020-001808-42", registry = "ClinicalTrials.gov")
#' which_trn("NCT 00312962 and euctr2020-001808-42")
#' which_trn("NCT 00312962 and euctr2020-001808-42", clean = FALSE)
#' which_trn("hello") # return NA if nothing was detected
which_trn <- function(x, registry = NULL, collapse = ";", clean = TRUE) {

  if (is.null(registry)) {
    trn_regex <- paste0(ctregistries::registries$trn_regex,
                        collapse = "|")
  } else {

    regexes <- ctregistries::registries |>
      dplyr::select(registry, trn_regex) |>
      tibble::deframe()

    trn_regex <- regexes[names(regexes) == registry]
  }


  trn <-
    stringr::str_extract_all(x,
                             trn_regex #,
                             # simplify = TRUE
    ) |>
    purrr::flatten_chr() |>
    stats::na.omit()

  if (rlang::is_empty(trn)) {
    # Return NA if no trn match
    return(NA_character_)
  } else {

    if (clean == TRUE) {
      trn <- furrr::future_map_chr(trn,
                                   \(x) clean_trn(x, quiet = TRUE))
    }

    if (collapse != "none") {
      trn <- trn |>
        unique() |>
        dplyr::na_if("NA") |>
        dplyr::na_if("") |>
        stats::na.omit() |>
        paste(collapse = collapse)
      if (trn == "") return(NA_character_)
    }
    return(trn)
  }
}


#' Extract TRN(s) from each element of input
#'
#' @param trn_vec A character vector.
#' @param registry Character. Registry of TRN. Default to NULL for all registries.
#' @param collapse A string. If one of the vector elements has multiple TRNS
#' the element of the output vector will be collapsed to a single
#' string with the chosen separator string. To make the length of the return
#' vector the same as the input vector, setting to "none" has been disabled.
#' @param clean Logical. By default (TRUE), the individual TRNs will also be
#' cleaned with `trn_clean`. Skip cleaning by setting to FALSE.
#'
#' @return A character vector, length of input.
#'
#' @importFrom rlang .data
#' @export
#' @examples
#' which_trns(c("NCT00312962", "hello", "euctr2020-001808-42", NA))
#'
#' df <- dplyr::tibble(trn = c("NCT00312962", "hello", "euctr2020-001808-42", NA))
#' dplyr::mutate(df, trn_extract = which_trns(trn))
#'
#' # Cannot set collapse to "none"
#' \dontrun{
#' which_trns(c("NCT00312962 and euctr2020-001808-41", "hello", "euctr2020-001808-42", NA),
#'  collapse = "none")
#' }
#' which_trns(c("NCT 00312962 and euctr2020-001808-41", "hello", "euctr2020-001808-42", NA))
#' which_trns(c("NCT 00312962 and euctr2020-001808-41", "hello", "euctr2020-001808-42", NA),
#'  clean = FALSE)
#' # character(0) and NULL in input vector are omitted
#' # convert NULL values to NAs when importing e.g. with `readr` to preserve
#' messy_input <- c("", NULL, "Nct02529358 and 2020-001934-37-ES", NA,
#' "EudraCT 2004-002714-11", character(0), "NA")
#' which_trns(messy_input)
#' which_trns(messy_input, registry = "EudraCT")
#' # compare to which_trn:
#' which_trn(messy_input)
which_trns <- function(trn_vec, registry = NULL, collapse = ";", clean = TRUE) {
  stopifnot("Collapse argument cannot be set to 'none'. Did you mean to use which_trn?" = collapse != "none")
  furrr::future_map_chr(trn_vec, \(x) which_trn(x, registry = registry, collapse = collapse, clean = clean), .progress = TRUE)
}
