#' @title Create or update a local, user-defined subset of the global online MODIS grid data pool
#' @description See \code{\link[MODIS]{getHdf}} in the \code{MODIS} package, which this function is essentially a copy-paste. The only changes are related
#' to the way MODIS files are downloaded which uses now the \code{rvest} package.
#' @usage getHdf(product, begin=NULL, end=NULL, tileH=NULL, tileV=NULL, extent=NULL,
#' collection=NULL, HdfName, quiet=FALSE, wait=0.5, checkIntegrity=FALSE, ...)
#' @inheritParams MODIS::getHdf
#' @author Matteo Mattiuzzi and Antoine Stevens
#' @return See \code{\link[MODIS]{getHdf}} in the MODIS package
#' @examples
#' \dontrun{
#' # Get NDVI MOD13A2 for June 2009
#' getHdf(product="MOD13A2", begin = "2009.06.01", end = "2009.06.30",
#'        tileH = 19, tileV = 5 ,job="H19V5")
#' }
#' @export
getHdf <- function(product, begin=NULL, end=NULL, tileH=NULL, tileV=NULL, extent=NULL,
                   collection=NULL, HdfName, quiet=FALSE, wait=0.5, checkIntegrity=FALSE,...){

  # Author: Matteo Mattiuzzi, Anja Klisch, matteo.mattiuzzi@boku.ac.at
  # Date : July 2011
  # Licence GPL v3

  # product="MOD13Q1"; begin="2010001"; end="2010005"; tileH=NULL; tileV=NULL; extent=NULL; collection=NULL; quiet=FALSE; wait=0.5; checkIntegrity=FALSE; z=1;u=1
  opts <- .combineOptions(...)

  sturheit <- .stubborn(level=opts$stubbornness)
  wait     <- as.numeric(wait)

  # TODO HdfName as regex
  if (!missing(HdfName))
  {
    HdfName <- unlist(HdfName)
    dates <- list()

    for (i in seq_along(HdfName))
    {
      HdfName[i]     <- basename(HdfName[i]) # separate name from path
      path           <- .genString(HdfName[i],...)
      path$localPath <- .setPath(path$localPath)

      if (!file.exists(paste0(path$localPath,"/",HdfName[i])))
      {
        .ModisFileDownloader(HdfName[i],quiet=quiet,...)
      }

      if(checkIntegrity)
      {
        .doCheckIntegrity(HdfName[i], quiet=quiet,...)
      }
      dates[[i]] <- paste0(path$local,"/",HdfName[i])
    }
    return(invisible(unlist(dates)))

  } else
  { # if HdfName isn't provided:

    if (missing(product))
      stop("Please provide the supported-'product'. See in: 'getProduct()'")

    #######
    # check product
    product <- getProduct(x=product,quiet=TRUE)
    # check if missing collection, else bilieve it
    if(is.null(collection))
      product$CCC <- getCollection(product=product,collection=collection,quiet=TRUE)[[1]]
    else
      product$CCC <- sprintf("%03d",as.numeric(unlist(collection)[1]))
    #########

    if (product$SENSOR[1]=="MODIS")
    {
      if (is.null(begin))
        cat("No begin(-date) set, getting data from the beginning\n")
      if (is.null(end))
        cat("No end(-date) set, getting data up to the most actual\n")

      # tranform dates
      tLimits <- transDate(begin=begin,end=end)
    }

    dates  <- list()
    output <- list() # path info for the invisible output
    l=0
    for(z in seq_along(product$PRODUCT))
    { # Platforms MOD/MYD

      if (product$TYPE[z]=="Swath")
      {
        cat("'Swath'-products not supported, jumping to the next.\n")
      } else
      {
        todo <- paste0(product$PRODUCT[z],".",product$CCC)
        for (u in seq_along(todo))
        {
          # tileID
          if (product$TYPE[z]=="CMG")
          {
            tileID="GLOBAL"
            ntiles=1
          } else
          {
            if (!is.null(tileH) & !is.null(tileV))
              extent <- getTile(tileH=tileH,tileV=tileV)
            else
              extent <- getTile(extent=extent)

            tileID <- extent$tile
            ntiles <- length(tileID)
          }

          onlineInfo <- .getStruc(product=product$PRODUCT[z],collection=product$CCC,server=opts$MODISserverOrder[1],begin=tLimits$begin,end=tLimits$end,wait=0)
          if(!is.na(onlineInfo$online))
          {
            if (!onlineInfo$online & length(opts$MODISserverOrder)==2)
            {
              cat(opts$MODISserverOrder[1]," seams not online, trying on '",opts$MODISserverOrder[2],"':\n",sep="")
              onlineInfo <- .getStruc(product=product$PRODUCT[z],collection=product$CCC,begin=tLimits$begin,end=tLimits$end,wait=0,server=opts$MODISserverOrder[2])
            }

            if(is.null(onlineInfo$dates))
              stop("Could not connect to server(s), and no data is available offline!\n")

            if(!is.na(onlineInfo$online))
              if(!onlineInfo$online)
                cat("Could not connect to server(s), data download disabled!\n")
          }

          datedirs <- as.Date(onlineInfo$dates)
          datedirs <- datedirs[!is.na(datedirs)]
          sel <- datedirs
          us  <- sel >= tLimits$begin & sel <= tLimits$end

          if (sum(us,na.rm=TRUE)>0)
          {
            suboutput <- list()
            l=l+1
            dates[[l]] <- datedirs[us]

            dates[[l]] <- cbind(as.character(dates[[l]]),matrix(rep(NA, length(dates[[l]])*ntiles),ncol=ntiles,nrow=length(dates[[l]])))
            colnames(dates[[l]]) <- c("date",tileID)

            for (i in 1:nrow(dates[[l]]))
            { # i=1
              #cat(dates[[l]][i,1],"\n")
              #flush.console()

              year <- format(as.Date(dates[[l]][i,1]), "%Y")
              doy  <- as.integer(format(as.Date(dates[[l]][i,1]), "%j"))
              doy  <- sprintf("%03d",doy)
              mtr  <- rep(1,ntiles) # for file availability flaging
              path <- .genString(x=strsplit(todo[u],"\\.")[[1]][1],collection=strsplit(todo[u],"\\.")[[1]][2],date=dates[[l]][i,1])

              for(j in 1:ntiles)
              {
                dates[[l]][i,j+1] <- paste0(strsplit(todo[u],"\\.")[[1]][1],".",paste0("A",year,doy),".",if (tileID[j]!="GLOBAL") {paste0(tileID[j],".")},strsplit(todo[u],"\\.")[[1]][2],".*.hdf$") # create pattern

                # print(dir(path$localPath,pattern=dates[[l]][i,j+1]))
                if (length(dir(path$localPath,pattern=dates[[l]][i,j+1]))>0)
                { # if available locally
                  HDF <- dir(path$localPath,pattern=dates[[l]][i,j+1]) # extract HDF file

                  if (length(HDF)>1)
                  { # in very recent files sometimes there is more than 1 file/tile/date if so get the most recent processing date
                    select <- list()
                    for (d in 1:length(HDF))
                    {
                      select[[d]]<- strsplit(HDF[d],"\\.")[[1]][5]
                    }
                    HDF <- HDF[which.max(unlist(select))]
                  }
                  dates[[l]][i,j+1] <- HDF
                  mtr[j] <- 0
                }
              }

              if (sum(mtr)!=0 & (onlineInfo$online | is.na(onlineInfo$online)))
              { # if one or more of the tiles in the given date is missing, its necessary to go online

                if(exists("ftpfiles"))
                {
                  rm(ftpfiles)
                }

              for (g in 1:sturheit)
                { # get list of FILES in remote dir
                  # server <- names(path$remotePath)[g%%length(path$remotePath)+1]
                  ftpfiles <- filesUrl(path$remotePath[[which(names(path$remotePath)==onlineInfo$source)]])
                  if(length(ftpfiles))
                    break

                  Sys.sleep(wait)
                }
                if(!length(ftpfiles))
                {
                  stop("Problems with online connections try a little later")
                }
                if (ftpfiles[1] != "total 0")
                {
                  ftpfiles <- unlist(lapply(strsplit(ftpfiles," "),function(x){x[length(x)]})) # found empty dir!

                  for(j in 1:ntiles)
                  { # j=1
                    if(mtr[j]==1)
                    { # if tile is missing get it
                      onFtp <- grep(ftpfiles,pattern=dates[[l]][i,j+1],value=TRUE)
                      HDF   <- grep(onFtp,pattern=".hdf$",value=TRUE)
                      if(length(HDF)>0)
                      {
                        if (length(HDF)>1)
                        { # in very recent files sometimes there is more than 1 file/tile/date if so get the last
                          select <- list()
                          for (d in seq_along(HDF))
                          {
                            select[[d]] <- strsplit(HDF[d],"\\.")[[1]][5]
                          }
                          HDF <- HDF[which.max(unlist(select))]
                        }

                        dates[[l]][i,j+1] <- HDF

                        hdf <- .ModisFileDownloader(HDF, wait=wait, quiet=quiet)
                        mtr[j] <- hdf

                      } else
                      {
                        dates[[l]][i,j+1] <- NA
                      }
                    }
                  }
                } else
                {
                  dates[[l]][i,(j+1):ncol(dates[[l]])] <- NA
                } # on ftp is possible to find empty folders!
              }
              if(checkIntegrity)
              { # after each 'i' do the sizeCheck
                isIn <- .doCheckIntegrity(paste0(path$localPath,dates[[l]][i,-1]), wait=wait, quiet=quiet,...)
              }
              suboutput[[i]] <- paste0(path$localPath,dates[[l]][i,-1])
            } # end i

            output[[l]] <-  as.character(unlist(suboutput))
            names(output)[l] <- todo[u]
          } else
          {
            cat(paste0("No files on ftp in date range for: ",todo[u],"\n\n"))
          }
        }
      }
    }
    return(invisible(output))
  }
} ## END: FTP vs ARC check and download
