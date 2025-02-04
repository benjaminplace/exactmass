#' Extract individual elements from a formula string
#'
#' @param composition.str character string containing an elemental formula
#'
#' @returns list object containing elements and counts of individual elements
#' @export
#'
#' @examples
#' extract_elements("C2H6O6")
extract_elements <- function(composition.str) {
  single.elem <- gregexpr("[A-Z]", composition.str)[[1]]
  double.elem <- gregexpr("[A-Z][a-z]", composition.str)[[1]]
  single.elem <- single.elem[!single.elem %in% double.elem]
  elements <- c()
  if (length(single.elem) > 0) {
    elements <- substring(composition.str, single.elem, single.elem)
  }
  if (length(double.elem) > 0) {
    if (double.elem[1] != -1) {
      elements <- c(elements, substring(composition.str, double.elem, double.elem + 1))
    }
  }
  ecounts <- rep(1, length(elements))
  nums <- gregexpr("[0-9]+", composition.str)[[1]]
  nums.count <- attr(nums, "match.length")
  for (i in 1:length(nums)) {
    if (substr(composition.str, nums[i] - 1, nums[i] - 1) %in% elements) {
      ecounts[which(elements == substr(composition.str, nums[i] - 1, nums[i] - 1))] <- substr(composition.str, nums[i], nums[i] + nums.count[i] - 1)
    }
    if (substr(composition.str, nums[i] - 2, nums[i] - 1) %in% elements) {
      ecounts[which(elements == substr(composition.str, nums[i] - 2, nums[i] - 1))] <- substr(composition.str, nums[i], nums[i] + nums.count[i] - 1)
    }
  }
  results <- c()
  results$elements <- elements
  results$counts <- as.integer(ecounts)
  results
}

#' Calculate the monoisotopic mass of an elemental formula
#'
#' @param elementcomp character string containing an elemental formula
#' @param ionstate character string containing the ion state of the formula in `[M+H]+` format
#'
#' @returns numeric monoisotopic mass
#' @import tibble
#' @export
#'
#' @examples
#' calculate_monoisotope("C2H6O6", ionstate = "[M+H]+")
calculate_monoisotope <- function(elementcomp, ionstate = "none") {
  exactmasses <- tibble(ELEMENT = exactmasschart$elements, EXACTMASS = exactmasschart$masses[,1])
  elementlist <- extract_elements(elementcomp)
  mass <- 0
  for (i in 1:length(elementlist$elements)) {
    mass <- mass + as.numeric(exactmasses$EXACTMASS[which(exactmasses$ELEMENT == elementlist$elements[i])])*elementlist$counts[i]
  }
  if (ionstate != "none") {
    mass_adj <- get_ionstate_massadj(ionstate)
    mass <- mass + mass_adj
  }
  mass
}


#' Get mass adjustment for compound's ion state
#'
#' @param ionstate character string containing the ion state of the formula in `[M+H]+` format
#'
#' @returns numeric monoisotopic mass
#' @import stringr
#' @export
#'
#' @examples
#' get_ionstate_massadj("[M-H]+")
get_ionstate_massadj <- function(ionstate = "none") {
  exactmasses <- tibble(ELEMENT = exactmasschart$elements, EXACTMASS = exactmasschart$masses[,1])
  massadj <- 0
  if (ionstate != "none") {
    chargestate <- str_sub(ionstate, start = str_locate(ionstate, "]")[1]+1)
      if (chargestate == "+" | chargestate == "-") {chargestate <- paste0(chargestate, "1")}
    ionfunc <- str_replace(str_sub(ionstate, str_locate(ionstate, "M[+-]")), "M", "")
    adduct <- str_sub(ionstate, start = str_locate(ionstate, "M[+-]")[2]+1, end = str_locate(ionstate, "]")[1]-1)
    if (!is.na(adduct)) {
      adductmass <- calculate_monoisotope(adduct, ionstate = "none")
      massadj <- switch(ionfunc, `-` = -adductmass, `+` = adductmass)
    }
    if (is.na(adduct)) {
      massadj <- 0
    }
    electrons <- -0.0005484 * eval(parse(text = chargestate)) #account for electrons, mass is MW of 1 electron
    massadj <- massadj + electrons
  }
  massadj
}
