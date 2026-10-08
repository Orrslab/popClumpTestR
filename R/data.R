#' Bird Location Data
#' 
#' Example bird location data using ATLAS in Harod Valley, Israel. 
#' There are 18 populations of 17 avian species, mostly in the 7-21 range,
#' with one large Spur-Winged Lapwing (Vanellus spinosus) with 63 individuals.
#' Note that the locations are in Israel Transverse Mercator (ITM) - https://en.wikipedia.org/wiki/Israeli_Transverse_Mercator
#'
#' @format ## `bird.location.data`
#' A data frame with 247 rows (of 18 populations of 17 avian species) and 3 columns:
#' \describe{
#'   \item{POPULATION}{Population name, in format "SPECIES AMOUNT"}
#'   \item{X}{X of centroid of the individual (in ITM)}
#'   \item{Y}{Y of centroid of the individual (in ITM)}
#'   ...
#' }
#' @source ATLAS
"bird.location.data"