#' @title  Access and extract MODIS imagery.
#'
#' @section Introduction:
#'
#' This package complements and replaces some of the functions of the MODIS package (v0.10-31; https://r-forge.r-project.org/R/?group_id=1252)
#' For the functions that are replaced, most of the code is simply copy-pasted from the MODIS package,
#' with a few bugs corrected and some changes to address specific needs.
#'
#' @section Functions that are replaced:
#'
#' \itemize{
#'    \item{\code{\link{getHdf}}}{Modifications are only related to the way MODIS files are downloaded}
#'    \item{\code{\link{runGdal}}}{Use \code{gdalUtils} to simplify the code, with a few bugs corrected}
#'    \item{\code{\link{orgTime}}}{Allows to define user-defined output time stamps with the \code{nDays} argument}
#'    \item{\code{\link{filesUrl}}}{Use \code{rvest} instead of \code{RCurl} package}
#' }
#'
#' @section New functions:
#'
#' \itemize{
#'    \item{\code{\link{convert_modis}}}{Convert MODIS files to a \code{\link[raster]{Raster-class}} object}
#'    \item{\code{\link{convert_dn_modis}}}{Convert DN to physical values}
#'    \item{\code{\link{interpolate_raster}}}{Interpolate Raster* time series, possibly to new time stamps (temporal resampling)}
#'    \item{\code{\link{gapfill_raster}}}{Fill gaps in Raster* time series}
#'    \item{\code{\link{convert_qf_modis}}}{Extract Quality Flags}
#' }
#'
#'
#' @note
#' Function arguments follow the same naming as in MODIS, as much as possible. The package should be attached AFTER the
#' MODIS package.
#'
#' @author Antoine Stevens and Matteo Mattiuzzi (MODIS package)
#' @name MODISExtra
#' @docType package
#' @import MODIS raster foreach
NULL
