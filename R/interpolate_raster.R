#' @title Interpolate a Raster* time series to new time stamps (temporal resampling)
#' @description Interpolate Raster* times series to new time stamps
#' @usage interpolate_raster(x, w = NULL, t = NULL, timeInfo = orgTime(x),
#'        method = c("spline","linear","constant","nn","whittaker"),
#'        lambda = 5000, nIter= 3, outlierThreshold = NULL,
#'        df = 6,
#'        cores, filename=rasterTmpFile(),...)
#' @param x A \code{\link[raster]{Raster-class}} object to interpolate
#' @param w Optional, a \code{\link[raster]{Raster-class}} object with weight information. See e.g. \code{\link[MODIS]{makeWeights}}.
#' @param t Optional, a \code{\link[raster]{Raster-class}} object with MODIS 'composite day of the year' data.
#' @param timeInfo object returned by \code{orgTime} defining the input and output time stamps to interpolate to.
#' @param method character \code{vector} of length 1 with one of the following value:
#' 'spline' for smoothing splines ,'linear' for linear approximation ,'constant' for constant approximation,
#' 'nn' for nearest neighbour , 'whittaker' for whittaked smoothing.
#' @param lambda Yearly lambda value. Default is 5000. See \code{\link[MODIS]{whittaker.raster}}.
#' @param nIter numeric \code{vector} of length 1 giving the number of iterations for the upper envelope fitting ('whittaker' \code{method}).
#'  Default is 3. See \code{\link[MODIS]{whittaker.raster}}.
#' @param outlierThreshold Numeric. If provided, \code{outlierTreshold} allows to remove outliers before smoothing by fitting a
#' whittaker function and removing values with residuals higher than \code{outlierTreshold}.
#' @param df Yearly degree of freedom for 'spline' \code{method}. Default is 6. See \code{\link[MODIS]{smooth.spline.raster}}
#' @param cores number of cores used when interpolating. Default is 1.
#' @param filename name of the file where the interpolated raster is saved. Default is created through \code{\link[raster]{rasterTmpFile}}
#' @param ... arguments passed to \code{\link[raster]{writeRaster}}
#' @details 'linear' and 'constant' \code{method} use the \code{\link[stats]{approx}} interpolation functions while 'spline' use
#' \code{\link[stats]{splinefun}}.
#' @return A \code{\link[raster]{Raster-class}} object with values interpolated to the output time series (as defined in \code{timeInfo})
#' @author Antoine Stevens & Matteo Mattiuzzi
#' @note Most of the code is borrowed from \code{\link[MODIS]{smooth.spline.raster}} and \code{\link[MODIS]{whittaker.raster}}.
#' Compared to the \code{MODIS} version, code is simplified and made general for reflectance values
#' @seealso \code{gapfill_raster}
#' @references Atzberger, C., and Eilers, P.H.C. (2011). A time series for monitoring vegetation activity and phenology at 10-daily time steps covering large parts of South America. International Journal of Digital Earth 4, 365-386.
#' @examples
#' \dontrun{
#' # Get and extract NDVI MOD13A2 for 2009
#' runGdal(product="MOD13A2", begin = "2009.01.01", end = "2009.12.31",
#'         tileH = 19, tileV = 5 ,
#'         SDSstring =  "100000000000",job="H19V5",
#'         gdalPath = "C:/OSGeo4W64/bin")
#' # process and convert to raster
#' path <- "~/MODIS_DATA/PROCESSED/H19V5"
#' r <- convert_modis(path = path,pattern = "16_days_NDVI",
#'                    begin = "2009.01.01", end = "2009.12.31",extractAll=T)
#' str(r)
#' r_NDVI <- r$raster_NDVI
#' # Smooth by splines
#' r_smooth <- interpolate_raster(x = r_NDVI, method = "spline")
#' # plot
#' x <- r$ti_NDVI$inputLayerDates
#' y1 <- r_NDVI[1000]
#' y2 <- r_smooth[1000]
#' plot(x,y1)
#' lines(x,y2,col = "red")
#' # now, using weights and composite day of the year to smooth
#' # first extract data
#' runGdal(product="MOD13A2", begin = "2009.01.01", end = "2009.12.31",
#'         tileH = 19, tileV = 5 ,
#'         SDSstring =  "001000000010", # extract VI_usefulness and composite_day_of_the_year
#'         job="H19V5",
#'         gdalPath = "C:/OSGeo4W64/bin")
#' cdoy <- convert_modis(path = path ,pattern = "16_days_composite",
#'                       begin = "2009.01.01", end = "2009.12.31")
#' QF <- convert_modis(path = path ,pattern = "16_days_VI_Quality",
#'                    begin = "2009.01.01", end = "2009.12.31", convertDN = FALSE)
#' # Convert QF to 'usefulness'
#' usefulness <- convert_qf_modis(r = QF,qf = "VI_usefulness")
#' # Convert usefulness to weights (see reference section)
#' w <- (15-usefulness)/15
#' r_smooth2 <- interpolate_raster(x = r_NDVI, w = w, t = cdoy, method = "spline")
#' y3 <- r_smooth2[1000]
#' lines(x,y3,col="green")
#' # One could aslo interpolate to new time stamps
#' ti <- orgTime(r_NDVI,nDays = 10,
#'               begin = "2009.01.01", end = "2009.12.31") # resample every 10 days
#' r_smooth_resampled <- interpolate_raster(x = r$raster_NDVI, w = w, t = cdoy,
#'                                          timeInfo = ti, method = "spline")
#' y4 <- r_smooth_resampled[1000]
#' points(ti$outputLayerDates,y4,pch=3)
#' }
#'
#' @export
interpolate_raster <- function(x, w=NULL, t=NULL, timeInfo = orgTime(x),
                               method = c("spline","linear","constant","nn","whittaker"),
                               lambda = 5000, nIter= 3, outlierThreshold = NULL,
                               df = 6,
                               cores, filename=rasterTmpFile(),...)
{
  method <- match.arg(method)

  opts <- MODIS:::combineOptions(...)

  if(is.null(opts$datatype))
    datatype <- dataType(x)
  else
    datatype <- opts$datatype

  if (is.null(opts$minDat))
    minDat <- 3 # 3 is very small!
  else
    minDat <- opts$minDat

  tsLength <- as.numeric(max(timeInfo$inputLayerDates) - (min(timeInfo$inputLayerDates))) - 1
  nx <- nlayers(x)

  if(method == "whittaker"){
    if (is.character(lambda))
    {
      cat("Using fixed 'lambda':",lambda,"\n")
      lambda <- as.numeric(lambda)
    } else
    {
      intlambda <- lambda
      lambda <- lambda*(tsLength/365)
      cat("Yearly 'lambda' is:",intlambda,"\nNow changed with lambda*('length of input data period in days'/365) to:",lambda,"\n")
    }
  }

  if(method == "spline"){
    if (is.character(df))
    {
      cat("Using fixed 'df':",df,"\n")
      df <- as.numeric(df)
    } else
    {
      indf <- df
      df   <- df*(tsLength/365)
      cat("Yearly 'df' is:",indf,"\nNow changed with df*('length of input data period in days'/365) to:",df,"\n")
    }
    if(df>=nx|df<1){
      df <- pmax(pmin(df,nx)-1,1)
      warning(paste0("df should be in 1 < df <= nlayers(x). The parameter is adapted to: ",df))
    }
  }

  .process_raster(x,w,t,timeInfo,
                  fun = .interp,
                  args = list(method = method,lambda = lambda,nIter = nIter,
                              outlierThreshold = outlierThreshold,
                              df = df, minDat = minDat,
                              ws = NULL, gap = FALSE),
                  cores, filename, datatype)
}

#' @title Fill NA values in Raster* time series by interpolation (gap filling)
#' @description Fill \code{NA} values in Raster* time series
#' @usage gapfill_raster(x, w=NULL, t=NULL, timeInfo = orgTime(x), method,
#'                       lambda = 5000, nIter= 3, outlierThreshold = NULL,
#'                       df = 6, ws = 2, cores, filename=rasterTmpFile(),...)
#' @inheritParams interpolate_raster
#' @param ws window size for the `fensholt` \code{method}. Default to 2. See details
#' @author Antoine Stevens & Matteo Mattiuzzi
#' @return A \code{\link[raster]{Raster-class}} object with NA values filled
#' @details The 'fensholt' \code{method} fill gaps in a time serie according to the following rule:
#  Given a pixel at time t (t0) and two previous (t-1, t-2) and subsequent (t+1, t+2) composite period values
#' If t0 is flagged as NA (e.g. cloudy) an average between t-1, t+1 is performed;
#' if only one of t-1 or t+1 is not missing t0 will be replaced by this value.
#' If both t-1 and t+1 are missing then the filter is extended to cover also t-2 and t+2 following the same rules as for
#' t-1 and t+1
#' @note Most of the code is borrowed from \code{\link[MODIS]{smooth.spline.raster}} and \code{\link[MODIS]{whittaker.raster}}
#' @references Fensholt, R., Rassmussen, K., Nielsen, T.T., and Mbow, C. (2009). Evaluation of earth observation based long term vegetation trends - Intercomparing NDVI time series trend analysis consistency of Sahel from AVHRR GIMMS, Terra MODIS and SPOT VGT data. Remote Sensing of Environment 113, 1886-1898.
#' @seealso \code{interpolate_raster}
#' @examples
#' # See interpolate_raster
#' @export
gapfill_raster <- function(x, w=NULL, t=NULL, timeInfo = orgTime(x),
                               method = c("spline","linear","constant","nn","whittaker","fensholt"),
                               lambda = 5000, nIter= 3, outlierThreshold = NULL,
                               df = 6, # parameters for splines
                               ws = 2,
                               cores, filename=rasterTmpFile(),...)
{
  # fensholt

  method <- match.arg(method)

  if(method=="fensholt"&!identical(timeInfo$inputLayerDates,timeInfo$outputLayerDates))
    stop("The Fensholt and Proud method is valid only when timeInfo$inputLayerDates and timeInfo$outputLayerDates are identical")

  opts <- MODIS:::combineOptions(...)

  if(is.null(opts$datatype))
    datatype <- dataType(x)
  else
    datatype <- opts$datatype

  if (is.null(opts$minDat))
    minDat <- 3 # 3 is very small!
  else
    minDat <- opts$minDat

  tsLength <- as.numeric(max(timeInfo$inputLayerDates) - (min(timeInfo$inputLayerDates))) - 1
  nx <- nlayers(x)

  if(method == "whittaker"){
    if (is.character(lambda))
    {
      cat("Using fixed 'lambda':",lambda,"\n")
      lambda <- as.numeric(lambda)
    } else
    {
      intlambda <- lambda
      lambda <- lambda*(tsLength/365)
      cat("Yearly 'lambda' is:",intlambda,"\nNow changed with lambda*('length of input data period in days'/365) to:",lambda,"\n")
    }
  }

  if(method == "spline"){
    if (is.character(df))
    {
      cat("Using fixed 'df':",df,"\n")
      df <- as.numeric(df)
    } else
    {
      indf <- df
      df   <- df*(tsLength/365)
      cat("Yearly 'df' is:",indf,"\nNow changed with df*('length of input data period in days'/365) to:",df,"\n")
    }
    if(df>=nx|df<1){
      df <- pmax(pmin(df,nx)-1,1)
      warning(paste0("df should be in 1 < df <= nlayers(x). The parameter is adapted to: ",df))
    }
  }

  .process_raster(x,w,t,timeInfo,
                  fun = .interp,
                  args = list(method = method,lambda = lambda,nIter = nIter,
                              outlierThreshold = outlierThreshold,
                              df = df, minDat = minDat,
                              ws = ws, gap = TRUE),
                  cores, filename,datatype)
}



.process_raster <- function(x, w, t, timeInfo,
                            fun,args,
                            cores, filename,datatype,...)
{
  # Main workhorse + check raster data consistency

  if(!inherits(x,"Raster"))
    x <- stack(x, quick=TRUE)

  if(!inherits(w,"Raster") & !is.null(w))
    w <- stack(w, quick=TRUE)

  if(!inherits(t,"Raster") & !is.null(t))
    t <- stack(t, quick=TRUE)

  # compare dimensions
  if(!is.null(w)){
    if(!compareRaster(x,w))
      stop("x and w rasters have different geometries")
    if(nlayers(x) != nlayers(w))
      stop("x and w rasters have a different number of bands")
  }

  if(!is.null(t)){
    if(!compareRaster(x,t))
      stop("x and t raster have different geometries")
    if(nlayers(x) != nlayers(t))
      stop("x and t rasters have a different number of bands")
  }

  b <- list()
  b[[1]] <- brick(x,nl=as.integer(length(timeInfo$outSeq)), values=FALSE)
  b[[1]] <- writeStart(b[[1]],filename, datatype=datatype, ...)

  tr <- blockSize(x)

  if (missing(cores))
  {
    for ( i in seq_along(tr$row) )
      b[[1]] <- writeValues(b[[1]], do.call(fun,c(list(i = i, row = tr$row, nrows = tr$nrow, x = x, w = w, t = t, timeInfo = timeInfo), args)), tr$row[i])
  } else
  {

    beginCluster(cores)

    cl <- getCluster()
    on.exit(endCluster())

    # send expr and data to cluster nodes
    parallel::clusterEvalQ(cl,{library(MODIS);library(rgdal);library(raster);library(ptw)})

    # number of blocks
    tr <- blockSize(x, minblocks=cores)
    for (i in 1:cores)
      parallel:::sendCall(cl[[i]],fun,c(list(i = i, row = tr$row, nrows = tr$nrow, x = x, w = w, t = t, timeInfo = timeInfo), args),tag=i)

    for (i in 1:tr$n)
    {
      d <- parallel:::recvOneData(cl);

      if (!d$value$success)
        stop("Cluster error in Row: ", tr$row[d$value$tag],"\n")

      b[[1]] <- writeValues(b[[1]], d$value$value, tr$row[d$value$tag])

      ni <- cores + i
      if (ni <= tr$n)
        parallel:::sendCall(cl[[d$node]],fun,c(list(i = ni, row = tr$row, nrows = tr$nrow, x = x, w = w, t = t, timeInfo = timeInfo), args),tag=ni)
    }
  }

  for (a in seq_along(b))
    b[[a]] <- writeStop(b[[a]])

  b <- brick(filename)
  b <- setZ(b,timeInfo$outputLayerDates)
  names(b) <- format(timeInfo$outputLayerDates,"%Y_%m_%d")
  return(b)
}


.interp <- function(i,row,nrows,x,w,t,timeInfo,method,lambda,nIter,outlierThreshold,df,minDat,ws,gap)
{
  val  <-  raster::getValues(x, row=row[i], nrows=nrows[i])

  yRow <- nrow(val)
  yCol <- ncol(val)
  isna <- is.na(val)

  set0   <- matrix(FALSE,nrow = yRow,ncol = yCol)
  set0[isna] <- TRUE
  set0[rowSums(val,na.rm=TRUE)==0] <- TRUE

  if (!is.null(w))
  {
    wtu <- raster::getValues(w, row=row[i], nrows=nrows[i])

    # is it a weight info?
    if(max(wtu,na.rm=T) > 1)
      stop("max(weights) > 1. Check...")
    set0[wtu==0] <- TRUE

  } else
  {
    # if no weighting info is available then weight = 1
    wtu <- matrix(1,nrow=yRow,ncol=yCol)
  }

  if (inherits(t,"Raster"))
  {
    inTu  <- raster::getValues(t, row=row[i], nrows=nrows[i])
    inTu  <- MODIS::repDoy(inTu,timeInfo,bias=timeInfo$inSeq[1]-1)
    set0[is.na(inTu)] <- TRUE
    set0[inTu <= 0] <- TRUE
    inTu[set0] <- 0
  } else
  {
    inTu <- matrix(timeInfo$inSeq,nrow=yRow,ncol=yCol,byrow=TRUE)
  }

  wtu[set0] <- 0 # the entire info to include or not a pixel is in "wtu"
  if(method %in% c("spline","whittaker"))
    val[set0] <- 0  # missing values are not allowed as input

  # for detecting outliers with whittaker
  if(!is.null(outlierThreshold)){
    kickOutlier <- function(y,w,lambda,threshold)
    {
      fTS <- ptw::whit2(y,w,lambda=lambda)
      w[abs(y-fTS) > threshold] <- 0
      return(w)
    }
  } else
  {
    # if missing(outlierThreshold) generate a fake function to avoid a per pixel "if"
    kickOutlier <- function(y,w,lambda,threshold)
      return(w)
  }

  # generate output matrix
  if (is.null(timeInfo)|gap)
  {
    outTu <- inTu # out time period
    if(gap)
      out <- val
    else
      out <- matrix(nrow=yRow, ncol=yCol)
  } else {
    outTu <- as.matrix(timeInfo$outSeq)
    if (ncol(outTu)==1)
      outTu <- matrix(outTu, nrow=yRow, ncol=length(outTu),byrow=T)
    out <- matrix(nrow=nrow(outTu), ncol=ncol(outTu))
  }

  # minimum "minVal" and "minDat" input values for filtering
  if(method=="spline"){
    Cvec <- rowSums(wtu>0) > df
  }
  else if(method=="whittaker")
    Cvec <- rowSums(wtu>0) >= minDat
  else
    Cvec <- rowSums(!isna) >= minDat

  if(!gap){ # better to place the if outside the loop (not sure this is so much faster, though)
    for (u in (1:yRow)[Cvec])
    {
      if(method=="spline"){
        out[u,] <- stats::predict(smooth.spline(y=val[u,],x=inTu[u,],w=wtu[u,],df=df, tol=1), outTu[u,])$y
      } else if(method %in% c("constant","linear")){
        out[u,] <- stats::approx(y=val[u,],x= inTu[u,] , outTu[u,],rule=2,method = method)$y
      } else if(method == "nn"){
        out[u,] <- nn(val[u,],inTu[u,],outTu[u,]) # nearest neighbour
      } else if (method =="whittaker"){
        wtu[u,] <- kickOutlier(val[u,], wtu[u,], lambda=lambda,threshold=outlierThreshold) # assign w to 0 for outliers
        out[u,] <- whittaker(val[u,], inTu[u,], outTu[u,], wtu[u,],lambda = lambda, nIter = nIter) # whittaker interpolation
      }
    }
  } else {
    for (u in (1:yRow)[Cvec])
    {
      isna <- set0[u,]
      if(any(isna)){
        if(method=="spline"){
          out[u,isna] <- stats::predict(smooth.spline(y=val[u,], x=inTu[u,], wtu[u,], df=df, tol=1), outTu[u,])$y[isna]
        } else if(method %in% c("constant","linear")){
          out[u,isna] <- stats::approx(y=val[u,], x=inTu[u,] , xout = outTu[u,],rule=2,method = method)$y[isna]
        } else if(method == "nn"){
          out[u,isna] <- nn(y=val[u,],x=inTu[u,],xout = outTu[u,])[isna] # nearest neighbour
        } else if (method =="whittaker"){
          wtu[u,] <- kickOutlier(y=val[u,], w=wtu[u,], lambda=lambda,threshold=outlierThreshold) # assign w to 0 for outliers
          out[u,isna] <- whittaker(y=val[u,], x=inTu[u,], xout = outTu[u,], w=wtu[u,],lambda = lambda, nIter = nIter)[isna] # whittaker interpolation
        } else {
          out[u,] <- fnp(val[u,],ws=2) # fensholt and proud
        }
      } else {
        out[u,] <- val[u,] # do nothing if no NA's
      }
    }
  }
  out[rowSums(abs(out))==0,] <- NA
  return(out)
}

# nearest neighbour interpolator
nn <- function(y,x,xout) {
  # find nearest neighbour and assign y values
  if(identical(x,xout)){
    stop("x and xout are identical, use the fensholt method with ws = 1 instead")
  } else {
    m = length(x)
    n = length(xout)
    # squared distance between vectors
    d <- matrix(rep(x^2, n), nrow=m)  +  matrix(rep(xout^2, m), nrow=m, byrow=TRUE) - 2 * tcrossprod(x,xout)
    return(y[apply(d,2,which.min)])
  }
}

# fensholt and proud gap filling algorithm
fnp <- function(y,ws=2) {
  # compute mean y values within a given window size
  # ws = window size
  if(!ws %in% 1:2)
    stop("window size should be 1 or 2")
  l <- length(y)
  nn <- pmin(1, l)
  isna <- is.na(y)
  if(!any(isna)){
    return(y)
  } else {
    # rbind the signal with a lag and lead of 1
    y1 <- rbind(c(y[-seq_len(nn)], NA),c(NA, y[seq_len(l - nn)])) # this is refactored from dplyr::lead and dplyr::lag
    y[isna] <- colMeans(y1[,isna,drop=F],na.rm=T)
    isnan <- is.nan(y)
    if(any(isnan)&ws==2){
      # lag of 2
      y2 <-  rbind(c(y1[1,-seq_len(nn)],NA),c(NA, y1[2,seq_len(l - nn)]))
      y[isnan] <- colMeans(y2[,isnan,drop=F],na.rm=T)
    }
    y[is.nan(y)] <- NA
    return(y)
  }
}

# whittaker interpolator
whittaker <- function(y,x,xout,w,lambda,nIter) {
  # upper enveloppe fitting, see Atzerberg and Eiler, 2003
  # all observed values
  # that lie below the fitted curve (or that were missing/flagged values) are
  # replaced by their fitted value. With the updated values, the smoothing is
  # repeated

  # add a safe length of data (because layer doy + effective composite doy)
  vec0 <- rep(0,max(x,xout) - min(x,xout) - 1  + 30)
  names(vec0) <- 1:length(vec0)
  y_long <- w_long <- vec0
  y_long[x] <- y
  w_long[x] <- w

  for(i in 1:nIter)
  {
    fTS <- ptw::whit2(y_long,w=w_long,lambda=lambda)
    idx <- y_long < fTS
    y_long[idx] <- fTS[idx]
  }
  names(fTS) <- names(vec0)
  fTS[xout]
}
