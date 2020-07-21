#' Madras Longitudinal Schizophrenia Study.
#'
#' Monthly records of presence/abscence of psychiatric symptom 'thought disorder' of 86 patients
#' over the first year after initial hospitalisation for disease.
#'
#' @format A data frame with 1032 observations and 5 variables
#' \describe{
#'   \item{id}{identification number of a patient}
#'   \item{y}{response 'thought disorder': 0 = absent, 1 = present}
#'   \item{month}{month since hospitalisation}
#'   \item{age}{age indicator: 0 = less than 20 years, 1 = 20 or over}
#'   \item{gender}{sex indicator: 0 = male, 1 = female}
#' }
#' @source Diggle PJ, Heagerty P, Liang KY, Zeger SL (2002). The analysis of Longitudinal Data, second ed., pp. 234-43. Oxford University Press, Oxford.\cr
#' <http://faculty.washington.edu/heagerty/Books/AnalysisLongitudinal/datasets.html>
#' @references Jokinen J. Fast estimation algorithm for likelihood-based
#'  analysis of repeated categorical responses.
#'  \emph{Computational Statistics and Data Analysis} 2006; 51:1509-1522.

"madras"
