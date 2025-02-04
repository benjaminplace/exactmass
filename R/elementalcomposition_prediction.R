#' Predict possible elemental compositions for a specific mass
#'
#' @param mass numeric mass to predict elemental composition (in Da)
#' @param error numeric relative error (in ppm)
#' @param min.error numeric minimum absolute error  (in Da)
#' @param elements list of characters containing elements to include in calculation
#' @param counts list of integer vectors of count of elements to include in calculation
#'
#' @returns table containing elemental combinations, calculated mass, ppm error, and double bond equivalents (DBE)
#' @import tibble
#' @export
#'
#' @examples
#' predictcomposition(498.9302, error = 5, min.error = 0.002, elements = c("C", "H", "O", "N", "F", "S"), counts = list(0:39, 0:72, 0:20, 0:20, 0:20, 0:10))
predictcomposition <- function(mass, error, min.error, elements = c("C", "H", "N", "O", "F", "S"), counts = list(0:39, 0:72, 0:20, 0:20, 0:20, 0:10)) {
  exactmasses <- tibble(ELEMENT = exactmasschart$elements, EXACTMASS = exactmasschart$masses[,1])
  exactmasses <- exactmasses[order(exactmasses$ELEMENT),]
  counts <- counts[order(elements)]
  elements <- sort(elements)
  emass <- exactmasses$EXACTMASS[which(exactmasses$ELEMENT %in% elements)]
  for (i in 1:length(counts)) {
    if (max(counts[[i]])*emass[i] > mass + max(mass*error*10^-6,min.error)) {counts[[i]] <- 0:ceiling((mass + max(mass*error*10^-6,min.error))/emass[i])}
  }
  combos <- as.matrix(expand.grid(counts))
  masses <- combos %*% emass
  output <- cbind(combos[which(masses >= mass - max(mass*error*10^-6,min.error) & masses <= mass + max(mass*error*10^-6,min.error)),], masses[which(masses >= mass - max(mass*error*10^-6,min.error) & masses <= mass + max(mass*error*10^-6,min.error))], (10^6)*(masses[which(masses >= mass - max(mass*error*10^-6,min.error) & masses <= mass + max(mass*error*10^-6,min.error))]-mass)/mass)
  colnames(output) <- c(elements, "mass", "error")
  DBE <- cbind(rep(1, nrow(output)))
  try(DBE <- DBE + output[,"C"], silent = TRUE)
  try(DBE <- DBE - output[,"H"]/2, silent = TRUE)
  try(DBE <- DBE + output[,"N"]/2, silent = TRUE)
  try(DBE <- DBE + output[,"P"]/2, silent = TRUE)
  try(DBE <- DBE + output[,"Si"], silent = TRUE)
  if (length(which(colnames(output) %in% c("Cl", "Br", "I", "F"))) > 1) {
    try(DBE <- DBE - rowSums(output[,which(colnames(output) %in% c("Cl", "Br", "I", "F"))])/2, silent = TRUE)
  }
  if (length(which(colnames(output) %in% c("Cl", "Br", "I", "F"))) == 1) {
    try(DBE <- DBE - output[,which(colnames(output) %in% c("Cl", "Br", "I", "F"))]/2, silent = TRUE)
  }
  output <- cbind(output, DBE)
  colnames(output) <- c(elements, "mass", "ppmerror", "DBE")
  output
}
