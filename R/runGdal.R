#' @title Process MODIS hdf with GDAL
#' @description Download MODIS data from archive and process the files with GDAL. This function
#' is basically copy-pasted from MODIS::runGdal, with a few bugs corrected. It uses also now the \code{GdalUtils}
#' package, to simplify the code
#' @usage runGdal(product, collection=NULL, begin=NULL, end=NULL, extent=NULL, tileH=NULL,
#' tileV=NULL, buffer=0, SDSstring=NULL, job=NULL, checkIntegrity=TRUE, wait=0.5,
#' quiet=FALSE,gdalPath = "",...)
#' @inheritParams MODIS::runGdal
#' @param gdalPath Path to the gdal binaries.
#' @return See \code{\link[MODIS]{runGdal}} in the MODIS package
#' @author Matteo Mattiuzzi and Antoine Stevens
#' @examples
#' \dontrun{
#' # Get and extract NDVI MOD13A2 for June 2009
#' runGdal(product="MOD13A2", begin = "2009.06.01", end = "2009.06.30",
#'        tileH = 19, tileV = 5 ,
#'        SDSstring =  "100000000000",job="H19V5",
#'        gdalPath = "C:/OSGeo4W64/bin")
#'}
#'
#' @export
runGdal <- function (product, collection = NULL, begin = NULL, end = NULL,
                     extent = NULL, tileH = NULL, tileV = NULL, buffer = 0, SDSstring = NULL,
                     job = NULL, checkIntegrity = TRUE, wait = 0.5, quiet = FALSE, gdalPath="",
                     ...)
{

  gdalUtils::gdal_setInstallation(search_path = gdalPath)
  opts <- .combineOptions(...)

  product <- getProduct(product, quiet = TRUE)
  product$CCC <- getCollection(product, collection = collection)
  tLimits <- transDate(begin = begin, end = end)
  dataFormat <- toupper(opts$dataFormat)

  if (dataFormat == "RAW BINARY")
    stop("in argument dataFormat='raw binary', format not supported by GDAL (it is MRT specific) type: 'options(\"MODIS_gdalOutDriver\")' (column 'name') to list available inputs")

  if (dataFormat == "HDF-EOS")
    dataFormat <- "HDF4IMAGE"

  if (dataFormat == "GEOTIFF")
    dataFormat <- "GTIFF"

  if (is.null(opts$gdalOutDriver)) {
    opts$gdalOutDriver <- .gdalWriteDriver()
    options(MODIS_gdalOutDriver = opts$gdalOutDriver)
  }

  if (dataFormat %in% toupper(opts$gdalOutDriver$name)) {
    dataFormat <- grep(opts$gdalOutDriver$name, pattern = paste("^", dataFormat, "$", sep = ""), ignore.case = TRUE, value = TRUE)
    of <- dataFormat
    extension <-.getExtension(dataFormat)
  }
  else {
    stop("in argument dataFormat='", opts$dataFormat, "', format not supported by GDAL type: 'gdalWriteDriver()' (column 'name') to list available inputs")
  }

  # Set outProj
  if (product$TYPE[1] == "Tile" | (all(!is.null(extent) | !is.null(tileH) & !is.null(tileV)) & product$TYPE[1] == "CMG")) {
    if(class(extent)%in%c("RasterLayer","RasterBrick","RasterStack"))
      if(is.na(projection(extent)))
        stop("Provide a projection to the raster object")
    requireNamespace("rgeos")
    extent <- getTile(extent = extent, tileH = tileH, tileV = tileV, buffer = buffer)
  }
  else {
    extent <- NULL
  }

  t_srs <- NULL
  cat("########################\n")
  if (!is.null(extent$target$outProj)) {
    outProj <- .checkOutProj(extent$target$outProj, tool = "GDAL")
    cat("outProj          = ", outProj, " (Specified by raster*/spatial* object)\n")
  }
  else {
    outProj <- .checkOutProj(opts$outProj, tool = "GDAL")
    cat("outProj          = ", outProj, "\n")
  }
  if (outProj == "asIn") {
    if (product$SENSOR[1] == "MODIS") {
      if (product$TYPE[1] == "Tile") {
        outProj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
      }
      else {
        outProj <- "+proj=longlat +ellps=clrk66 +no_defs"
      }
    }
    else if (product$SENSOR[1] == "SRTM") {
      outProj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    }
  }
  t_srs <- outProj

  # Set pixel size
  if (!is.null(extent$target$pixelSize)) {
    pixelSize <- extent$target$pixelSize
    cat("pixelSize        = ", pixelSize, " (Specified by raster* object)\n")
  }
  else {
    pixelSize <- opts$pixelSize
    cat("pixelSize        = ", pixelSize, "\n")
  }
  tr <- NULL
  if (pixelSize[1] != "asIn") {
    if (length(pixelSize) == 1) {
      tr <- c(pixelSize, pixelSize)
    }
    else {
      tr <- pixelSize
    }
  }

  # set resampling type
  opts$resamplingType <- .checkResamplingType(opts$resamplingType, tool = "gdal")
  cat("resamplingType   = ", opts$resamplingType, "\n")
  r <- opts$resamplingType

  # set source projection
  if (product$SENSOR[1] == "MODIS") {
    if (product$TYPE[1] == "Tile") {
      s_srs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    }
    else {
      s_srs <- "+proj=longlat +ellps=clrk66 +no_defs"
    }
  }
  else if (product$SENSOR[1] == "SRTM") {
    s_srs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  }

  # Set extent
  te <- NULL
  if (!is.null(extent$target$extent)) {
    if (is.null(extent$target$outProj)) {
      rx <- raster(extent$target$extent, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      rx <- projectExtent(rx, outProj)
      rx <- extent(rx)
    }
    else {
      rx <- extent$target$extent
    }
    te <- c(rx@xmin, rx@ymin, rx@xmax, rx@ymax)
  }
  if (is.null(extent$target)) {
    if (!is.null(extent$extent)) {
      rx <- raster(extent$extent, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      rx <- projectExtent(rx, outProj)
      rx <- extent(rx)
      te <- c(rx@xmin, rx@ymin, rx@xmax, rx@ymax)
    }
  }

  # Set misc options
  if (is.null(opts$blockSize))
    bs <- NULL
  else
    bs <- list(paste0("BLOCKYSIZE=", as.integer(opts$blockSize)))

  if (is.null(opts$compression)|isTRUE(opts$compression))
    cp <- list("compress=lzw","predictor=2")
  else
    cp <- NULL

  bscp <- c(bs,cp)

  if (quiet)
    q <- T
  else
    q <- F
  # Download, extract and project
  for (z in seq_along(product$PRODUCT)) {
    todo <- paste(product$PRODUCT[z], ".", product$CCC[[product$PRODUCT[z]]], sep = "")
    if (z == 1) {
      if (is.null(job)) {
        job <- paste0(todo[1], "_", format(Sys.time(), "%Y%m%d%H%M%S"))
        cat("Output directory = ", paste0(normalizePath(opts$outDirPath, "/", mustWork = FALSE), "/", job), " (no 'job' name specified, generated (date/time based))\n")
      } else {
        cat("Output Directory = ", paste0(normalizePath(opts$outDirPath, "/", mustWork = FALSE), "/", job), "\n")
      }
      cat("########################\n")
      outDir <- file.path(opts$outDirPath, job, fsep = "/")
      dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
    }
    if(!is.null(opts$auxPath)){
      # create dir if necessary
      dir.create(opts$auxPath,showWarnings = F)

    } else {
      stop("Specify an auxPath directory. See .combineOptions()")
    }
    for (u in seq_along(todo)) {
      ftpdirs <- list()
      ftpdirs[[1]] <- as.Date(.getStruc(product = strsplit(todo[u],"\\.")[[1]][1], collection = strsplit(todo[u], "\\.")[[1]][2], begin = tLimits$begin, end = tLimits$end, server = opts$MODISserverOrder[1])$dates)
      prodname <- strsplit(todo[u], "\\.")[[1]][1]
      coll <- strsplit(todo[u], "\\.")[[1]][2]
      avDates <- ftpdirs[[1]]
      avDates <- avDates[avDates != FALSE]
      avDates <- avDates[!is.na(avDates)]
      sel <- as.Date(avDates)
      us <- sel >= tLimits$begin & sel <= tLimits$end
      if (sum(us, na.rm = TRUE) > 0) {
        avDates <- avDates[us]
        for (l in seq_along(avDates)) {
          # download files
          files <- unlist(getHdf(product = prodname,
                                 collection = coll, begin = avDates[l], end = avDates[l],
                                 tileH = extent$tileH, tileV = extent$tileV,
                                 checkIntegrity = checkIntegrity, stubbornness = opts$stubbornness,
                                 MODISserverOrder = opts$MODISserverOrder))
          files <- files[basename(files) != "NA"]
          if (length(files) > 0) {
            w <- getOption("warn")
            options(warn = -1)
            SDS <- list()
            for (z in seq_along(files)) {
              SDS[[z]] <- getSds(HdfName = files[z], SDSstring = SDSstring, method = "GDAL")
            }
            options(warn = w)
            if (!exists("NAS")) {
              NAS <- .getNa(SDS[[1]]$SDS4gdal)
            }
            for (i in seq_along(SDS[[1]]$SDSnames)) {
              outname <- paste0(paste0(strsplit(basename(files[1]), "\\.")[[1]][1:2], collapse = "."), ".", gsub(SDS[[1]]$SDSnames[i], pattern = " ", replacement = "_"), extension)
              gdalSDS <- sapply(SDS, function(x) {x$SDS4gdal[i]})
              # get no data values
              naID <- which(SDS[[1]]$SDSnames == names(NAS)[i])
              if (length(naID) > 0) {
                srcnodata <- as.character(NAS[[naID]]);dstnodata <- as.character(NAS[[naID]])
              } else {
                srcnodata <- NULL ; dstnodata <- NULL
              }

              if (length(grep(todo, pattern = "M.D13C2\\.005")) > 0) {
                if (i == 1) {
                  cat("\n###############\nM.D13C2.005 is likely to have a problem in metadata extent information, it is corrected on the fly\n###############\n")
                }
                ranpat <- .makeRandomString(length = 21)
                randomName <- paste0(outDir, "/deleteMe_", ranpat, ".tif")
                on.exit(unlink(list.files(path = outDir, pattern = ranpat, full.names = TRUE), recursive = TRUE))
                for (ix in seq_along(gdalSDS)) {
                  gdalUtils::gdal_translate(src_dataset = gdalSDS[ix],dst_dataset = randomName[ix],a_nodata=NAS[[naID]],a_ullr = c(-180,90,180,-90))
                }
                gdalSDS <- randomName
              }
              # extract and project
              ifile <-gdalSDS
              ofile <- paste0(outDir, "/", outname)
              args <- list(srcfile = ifile, dstfile = ofile, s_srs = s_srs, t_srs = t_srs,
                           of = of, te = te, tr = tr,
                           co = bscp, r = r, q = q,
                           srcnodata = srcnodata, dstnodata = dstnodata,
                           overwrite = TRUE, output_Raster = FALSE,multi = TRUE) # arguments to gdalwarp
              #remove args not used
              args <- args[!sapply(args,is.null)]
              # perform warping
              do.call(gdalUtils::gdalwarp,args)

              if (length(grep(todo, pattern = "M.D13C2\\.005")) > 0)
                unlink(list.files(path = outDir, pattern = ranpat, full.names = TRUE), recursive = TRUE)

            }
          } else {
            warning(paste0("No file found for date: ", avDates[l]))
          }
        }
      } else {
        warning("No product found within the date range")
      }
    }
  }
}
