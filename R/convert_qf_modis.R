#' @title Convert bit values from a QFLAG MODIS product
#' @description Convert a Raster* object representing MODIS QFLAG bit values to a Raster* with two categories representing flagged and non-flagged pixels
#' @usage convert_qf_modis(r,qf,type,cl = NULL,filename = rasterTmpFile(), ...)
#' @param r A \code{\link[raster]{Raster-class}} object
#' @param qf Quality Flag to be extracted. Can be more than one of these:
#' If type = 'MOD13': 'MODLAND_QA','VI_usefulness','aerosol_quantity','adjacent_cloud','atmosphere_brdf','mixed_clouds','land_water_flag','snow','shadow'.
#' If type = 'MOD15':'MODLAND_QC','sensor','dead_detector','cloud','scf_qc'.
#' If type = 'MOD15Extra':'sea', 'snow', 'aerosol', 'cirrus', 'cloud', 'shadow', 'biome'
#' @param type character \code{vector} of length 1 giving the MODIS data type. Should be one of these: 'MOD13', 'MOD15', 'MOD15Extra'
#' @param cl cluster object for parallel processing. Default is \code{NULL}
#' @param filename Passed to \code{\link[raster]{writeRaster}}.
#' @param ... arguments passed to \code{\link[raster]{writeRaster}}.
#' @return A \code{\link[raster]{RasterBrick-class}} object with two values (1,0), representing respectively pixels that are flagged by at least one of the give \code{qf},
#'        and pixels that are not flagged in any of the given \code{qf}
#' @references https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13a2,
#' https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod15a2
#' @seealso \code{\link[MODIS]{extractBits}} and \code{\link{interpolate_raster}} for an example
#' @author Antoine Stevens
#' @export
convert_qf_modis <- function(r,qf,type=c("MOD13","MOD15","MOD15Extra"), cl = NULL, filename = rasterTmpFile(), ...){
  
  if(!is.null(cl))
    if(!"cluster"%in%class(cl))
      stop("cl should be a cluster object. See ?getCluster")
      
  if (!inherits(r, "Raster"))
    stop("r should be a Raster* object")

  type <- match.arg(type)
  if(length(type)>1)
    stop("type should be of length 1")

  if (type=="MOD13") {
    pattern <- c("MODLAND_QA","VI_usefulness","aerosol_quantity","adjacent_cloud","atmosphere_brdf","mixed_clouds","land_water_flag","snow","shadow")
  } else if (type=="MOD15") {
    pattern <- c("MODLAND_QC","sensor","dead_detector","cloud","scf_qc")
  } else {
    pattern <- c("sea", "snow", "aerosol", "cirrus", "cloud", "shadow", "biome")
  }

  q <- stringr::str_detect(paste(qf,collapse="|"),pattern)

  if(!sum(q))
   stop(paste0("qf should match one of these : ", paste(pattern,collapse=", ")))

  b <- list()
  b[[1]] <- brick(r,nl=nlayers(r), values=FALSE)
  names(b[[1]]) <- names(r)
  if(length(getZ(r)))
    b[[1]] <- setZ(b[[1]],getZ(r))
  b[[1]] <- writeStart(b[[1]], filename = filename,...)
  
  if (is.null(cl)) {
    tr <- blockSize(r)
    for ( i in seq_along(tr$row) )
      b[[1]] <- writeValues(b[[1]], .convert_qf_mod(i = i, r = r, row = tr$row, nrows = tr$nrow,pattern = pattern,q = q, type = type), tr$row[i])
  } else {
    
    cores <- length(cl)
    # send expr and data to cluster nodes
    parallel::clusterEvalQ(cl,{library(sfsmisc)})
    # number of blocks
    tr <- blockSize(r, minblocks=cores)
    for (i in 1:cores)
      raster:::.sendCall(cl[[i]],.convert_qf_mod,list(i = i, r = r, row = tr$row, nrows = tr$nrow,pattern = pattern,q = q, type = type),tag=i)

    for (i in 1:tr$n)
    {
      d <- raster:::.recvOneData(cl);
      if (!d$value$success)
        stop("Cluster error in Row: ", tr$row[d$value$tag],"\n")

      b[[1]] <- writeValues(b[[1]], d$value$value, tr$row[d$value$tag])

      ni <- cores + i
      if (ni <= tr$n)
        raster:::.sendCall(cl[[d$node]],.convert_qf_mod,list(i = ni, r = r, row = tr$row, nrows = tr$nrow,pattern = pattern,q = q, type = type),tag=ni)
    }
  }

  for (a in seq_along(b))
    b[[a]] <- writeStop(b[[a]])

  b <- brick(filename)
  b
}

.convert_qf_mod <- function(i,r,row,nrows,pattern,q,type){
  val  <-  raster::getValues(r, row=row[i], nrows=nrows[i])
  n <- nrow(val)
  p <- ncol(val)
  val <- as.factor(val)
  lev <- as.numeric(levels(val))
  # This binary bit-string is parsed from right to left, and the individual bits within a bit-field are read from left to right
  # All	HDF-EOS	products	are	written	in	the	big-endian	referencing	scheme.	The	bits are always	numbered
  # from	right	(least-significant	bit)	to	left	(most-significant	bit).
  # See MODIS_LP_QA_Tutorial-2
  if(type=="MOD13"){
    bits <- sapply(lev,function(x)as.integer(intToBits(x)[1:16])) # intToBits convert to least-significant bit first (big endian)
    MODLAND_QA <- sfsmisc::as.intBase(bits[2:1,])
    VI_usefulness <- sfsmisc::as.intBase(bits[6:3,])
    aerosol_quantity <- sfsmisc::as.intBase(bits[8:7,])
    adjacent_cloud <- sfsmisc::as.intBase(bits[9,,drop=F])
    atmosphere_brdf <- sfsmisc::as.intBase(bits[10,,drop=F])
    mixed_clouds <- sfsmisc::as.intBase(bits[11,,drop=F])
    land_water_flag <- sfsmisc::as.intBase(bits[14:12,])
    snow <- sfsmisc::as.intBase(bits[15,,drop=F])
    shadow <- sfsmisc::as.intBase(bits[16,,drop=F])
  } else if (type=="MOD15") {
    bits <- sapply(lev,function(x)as.integer(intToBits(x)[1:8]))
    MODLAND_QC <- sfsmisc::as.intBase(bits[1,,drop=F])
    sensor <- sfsmisc::as.intBase(bits[2,,drop=F])
    dead_detector <- sfsmisc::as.intBase(bits[3,,drop=F])
    cloud <- sfsmisc::as.intBase(bits[5:4,])
    scf_qc <- sfsmisc::as.intBase(bits[8:6,])
  } else {
    bits <- sapply(lev,function(x)as.integer(intToBits(x)[1:8]))
    sea <- sfsmisc::as.intBase(bits[2:1,])
    snow <- sfsmisc::as.intBase(bits[3,,drop=F])
    aerosol <- sfsmisc::as.intBase(bits[4,,drop=F])
    cirrus <- sfsmisc::as.intBase(bits[5,,drop=F])
    cloud <- sfsmisc::as.intBase(bits[6,,drop=F])
    shadow <- sfsmisc::as.intBase(bits[7,,drop=F])
    biome <- sfsmisc::as.intBase(bits[8,,drop=F])
  }
  if(sum(q)==1)
    pat <- get(pattern[q])
  else
    pat <- as.numeric(as.logical(do.call(pmax,lapply(pattern[q], function(x)get(x)))))
  # replace values
  val <- matrix(pat[as.numeric(val)],nrow=n,ncol=p)
  val
}
