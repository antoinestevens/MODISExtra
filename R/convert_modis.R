#' @title Convert MODIS files
#' @description Convert MODIS files into a Raster* object, and possibly convert raw DN to physical values.
#' @usage convert_modis(path = ".", pattern = NULL, type, convertDN = TRUE,
#'                      extractAll = FALSE, filename = rasterTmpFile(), datatype = "FLT4S",
#'                      overwrite = TRUE, ...)
#' @param path Path to the folder where MODIS files are stored. Default is the working directory. See \code{\link[MODIS]{preStack}}.
#' @param pattern optional regular expression, only file names matching the regular expression will be extracted. See \code{\link[MODIS]{preStack}}.
#' @param type character \code{vector} of length 1 giving the type of MODIS data. Should be one of these values:
#' 'Fpar_1km','Lai_1km','FparLai_QC','FparExtra_QC','FparStdDev_1km','LaiStdDev_1km','NDVI','EVI','VI','red','NIR','view_zenith_angle','sun_zenith_angle','relative_azimuth_angle','composite_day_of_the_year','pixel_reliability'
#' @param convertDN a logical value indicating whether MODIS DN values should be converted to physical values.
#' See \code{\link{convert_dn_modis}}. Default is \code{TRUE}.
#' @param extractAll a logical value indicating whether time info, file names and date object should be returned in addition to the processed raster object.
#' Default is \code{FALSE}.
#' @param filename Passed to \code{\link[raster]{writeRaster}}.
#' @param datatype Passed to \code{\link[raster]{writeRaster}}.
#' @param overwrite Passed to \code{\link[raster]{writeRaster}}.
#' @param ... arguments passed to \code{\link{orgTime}} (for instance to restrict the studied period with \code{begin} or \code{end}
#'        or define new time stamps to interpolate to with \code{nDays})
#' @return a \code{list} with the following elements:
#' \itemize{
#'    \item{\code{raster}}{the processed Raster* object}
#'    \item{\code{raster_fill}}{if \code{type} is 'Fpar_1km|Lai_1km|FparStdDev_1km|LaiStdDev_1km', fill values are automatically extracted (https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod15a2)}
#'    \item{\code{timeInfo}}{time info as returned by \code{\link[MODIS]{orgTime}}}
#'    \item{\code{filename}}{file names as returned by \code{\link[MODIS]{preStack}}}
#'    \item{\code{date}}{\code{vector} of dates as returned by \code{\link[MODIS]{extractDate}}}
#' }
#' If there is only one element in the \code{list}, the \code{convert_modis} returns the \code{unlist}-ed Raster object.
#' @author Antoine Stevens
#' @examples
#' \dontrun{
#' # Get and extract NDVI MOD13A2 for JUNE 2009
#' runGdal(product="MOD13A2", begin = "2009.06.01", end = "2009.06.30",
#'        tileH = 19, tileV = 5 ,
#'        SDSstring =  "100000000000",job="H19V5")
#' # process and convert to raster
#' path <- "~/COPERNICUS/MODIS_DATA/PROCESSED/H19V5"
#' convert_modis(path = path,pattern = "16_days_NDVI")
#' }
#'
#' @export
convert_modis <- function(path = ".",pattern = NULL,type,convertDN = TRUE,extractAll = FALSE,filename = rasterTmpFile(),datatype="FLT4S",overwrite = TRUE,...) {

  capture.output(fn <- preStack(path = path, pattern = pattern)) # capture.output to ignore annoying calls to print in preStack

  if(is.null(fn))
    stop("No data found")

  if(missing(type)){
    type <- stringr::str_extract(fn[1],"Fpar_1km|Lai_1km|FparLai_QC|FparExtra_QC|FparStdDev_1km|LaiStdDev_1km|NDVI|EVI|VI_Quality|red|NIR|view_zenith_angle|sun_zenith_angle|relative_azimuth_angle|composite_day_of_the_year|pixel_reliability")
    if(is.na(type))
      stop("type is missing")
  }

  ti <- orgTime(fn,...)  # ... is to possibly restrict the studied period
  fn <- preStack(path = path, pattern = pattern,timeInfo = ti)
  d <- extractDate(fn,asDate=T)

   # create raster
  modis <- stack(fn)
  modis_list <- list() # to store results

  if(stringr::str_detect(type,"Fpar_1km|Lai_1km|FparStdDev_1km|LaiStdDev_1km")){
    # extract fill values for MOD15A2 and create a new raster*
    # see here https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod15a2
    modis_fill <- clamp(modis,lower = 248,useValues=F) - 248 # convert to zero-based indices
    modis_fill_levels <- c("no sd available","unclassified land cover", "urban", "permanent wetlands", "permanent snow/ice", "barren" , "water" , "fill reflectance value")
    # write to disk
    modis_fill <- add_raster_attributes(modis_fill,modis_fill_levels)
    # set time info
    modis_fill <- setZ(modis_fill,d$inputLayerDates)
    # convert to brick and write to disk
    modis_fill <- brick(writeRaster(modis_fill,filename = paste0(tools::file_path_sans_ext(filename),"_fill",extension(filename))))
    # names(modis_fill) <- format(d$inputLayerDates,"%Y_%m_%d")
    # assign to global env
    modis_list[[2]] <- list(modis_fill)
    names(modis_list)[2] <- "raster_fill"
  }

  # Set time info
  modis <- setZ(modis,d$inputLayerDates)
  # names(modis) <- format(d$inputLayerDates,"%Y_%m_%d")

  # convert DN to VI
  if(convertDN)
    modis <- convert_dn_modis(modis,type)

  # Convert to brick and write to disk
  modis <- writeRaster(modis,filename = filename,datatype = datatype,overwrite = overwrite)

  meta_list <- list()
  # assign result to global env
  if(extractAll){
    meta_list[[1]] <- ti
    meta_list[[2]] <- fn
    meta_list[[3]] <- d
    names(meta_list) <- c("timeInfo","filename","date")
  }
  modis_list[[1]] <- modis
  names(modis_list)[1] <- "raster"
  result <- c(modis_list,meta_list)
  if(length(result)==1)
    return(result[[1]])
  else
    return(result)
}

#' @title Add attributes to a Raster* object
#' @description Convert a Raster* object with numeric values to classes
#' @usage add_raster_attributes(r,attribute)
#' @param r A \code{\link[raster]{Raster-class}} object to add attributes (or categories) to
#' @param attribute a character \code{vector} of class names, with the same length as the number
#' of unique numeric values in the input raster object
#' @details
#' The function uses rgdal::writeGDAL with \code{catNames} argument to write a tif file to disk with attributes
#' and read the tif back to R
#' @author Antoine Stevens
add_raster_attributes <- function(r,attribute){
  out_raster <- sub("grd","tif",rasterTmpFile()) # temporary name
  n <- nlayers(r)
  rgdal::writeGDAL(as(r,"SpatialGridDataFrame"), out_raster, catNames=rep(list(attribute),n),mvFlag=255L)
  if(n == 1)
    return(raster(out_raster))
  else
    return(brick(out_raster))
}

#' @title Convert MOD13A2 and MOD15A2 DN values to physical values
#' @description Set values outside valid range to \code{NA} and scale input data using pre-determined gain values
#' @usage convert_dn_modis(r,type)
#' @param r A \code{\link[raster]{Raster-class}} object to convert
#' @param type character \code{vector} of length 1 giving the type of MODIS data. Should be one of these values:
#' 'Fpar_1km','Lai_1km','FparLai_QC','FparExtra_QC','FparStdDev_1km','LaiStdDev_1km','NDVI','EVI','VI','red','NIR','view_zenith_angle','sun_zenith_angle','relative_azimuth_angle','composite_day_of_the_year','pixel_reliability'
#' @return A \code{\link[raster]{Raster-class}} object with converted DN values
#' @note It is used in \code{\link{convert_modis}}. This function is obviously valid for MOD13A2 and MOD15A2 only !
#' @references Solano et al. 2010. MODIS Vegetation Index User's Guide (MOD13 Series), 42 p.
#' http://gis.stackexchange.com/questions/84058/reprojected-modis-ndvi-has-range-from-32768-to-32767-expected-1-to-1
#' @author Antoine Stevens
#' @seealso \code{\link{convert_modis}}
#'
#' @export
convert_dn_modis <- function(r,type = c("Fpar_1km","Lai_1km","FparLai_QC","FparExtra_QC","FparStdDev_1km","LaiStdDev_1km","NDVI","EVI","VI","red","NIR","view_zenith_angle","sun_zenith_angle","relative_azimuth_angle","composite_day_of_the_year","pixel_reliability")){

  if (!inherits(r, "Raster"))
    stop("r should be a Raster* object")

  type <- match.arg(type)
  if(length(type)>1)
    stop("type should be of length 1")

  n <- names(r)
  if(length(getZ(r)))
    d <- getZ(r)
  else
    d <- NULL

  if(missing(type)){
    type <- stringr::str_extract(n[1],"Fpar_1km|Lai_1km|FparLai_QC|FparExtra_QC|FparStdDev_1km|LaiStdDev_1km|NDVI|EVI|VI|red|NIR|view_zenith_angle|sun_zenith_angle|relative_azimuth_angle|composite_day_of_the_year|pixel_reliability")
    if(is.na(type))
      stop("type is missing")
  }
  trans_fac <- switch(type,
                      Fpar_1km = c(0.01,0,100), # gain, min, max
                      Lai_1km = c(0.1,0,100),
                      FparLai_QC = c(NA,0,254),
                      FparExtra_QC = c(NA,0,254),
                      FparStdDev_1km = c(0.01,0,100),
                      LaiStdDev_1km = c(0.1,0,100),
                      NDVI = c(0.0001,-2000,10000),
                      EVI = c(0.0001,-2000,10000),
                      VI = c(NA,0,65534),
                      red = c(0.0001,0,10000),
                      NIR = c(0.0001,0,10000),
                      blue = c(0.0001,0,10000),
                      MIR = c(0.0001,0,10000),
                      view_zenith_angle = c(0.01,-9000,9000), # view zenith angle
                      sun_zenith_angle = c(0.01,-9000,9000), # sun zenith angle
                      relative_azimuth_angle = c(0.1,-3600,3600), # relative azimuth angle
                      composite_day_of_the_year = c(NA,1,366), # composite day of the year
                      pixel_reliability = c(NA,0,4) # pixel reliability
  )
  # make pixel outside valid range to NA
  r <- clamp(r,lower = trans_fac[2], upper = trans_fac[3], useValues = FALSE)

  if(!is.na(trans_fac[1]))
    gain(r) <-  trans_fac[1] # apply scale factor
  names(r) <- n
  if(!is.null(d))
    r <- setZ(r,d)
  r
}
