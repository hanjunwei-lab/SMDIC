#' A large list of 24 immune cells type-specific gene signatures from Bindea et al
#'
#' It's a built-in data. The name of the list represent 24 immune cells, the value of the list are 24 immune cells type-specific gene signatures from Bindea et al
#'
#' @format A list
#' @references
#' Bindea G, Mlecnik B, Tosolini M, Kirilovsky A, Waldner M, Obenauf AC, et al. Spatiotemporal dynamics of intratumoral immune cells reveal the immune landscape in human cancer. Immunity. 2013;39:782–95.
'immunelist'

#' A data.frame of 64 immune cells name from xCell method
#'
#' It's a built-in data. The first column represents the abbreviation of 64 immune cells, the second column represents the full name of 64 immune cells
#'
#' @format A data.frame with 64 rows and 2 column
#' @references
#' Aran D , Hu Z , Butte A J . xCell: digitally portraying the tissue cellular heterogeneity landscape[J]. Genome Biology, 2017, 18(1):220.
'cell64'

#' A data.frame of 24 immune cells name from Bindea et al
#'
#' It's a built-in data. The first column represents the abbreviation of 24 immune cells, the second column represents the full name of 24 immune cells
#'
#' @format A data.frame with 24 rows and 2 column
#' @references
#' Bindea G, Mlecnik B, Tosolini M, Kirilovsky A, Waldner M, Obenauf AC, et al. Spatiotemporal dynamics of intratumoral immune cells reveal the immune landscape in human cancer. Immunity. 2013;39:782–95.
'cell24'

#' @title xCell datasets
#' @description xCell datasets. It's a built-in data.
#' @usage xCell.data
#'
#' @format list:
#' \describe{
#'   \item{spill}{spillover matrix and calibration parameters}
#'   \item{signatures}{the signatures for calculating scores}
#'   \item{genes}{genes to use to calculate xCell}
#' }
"xCell.data"

