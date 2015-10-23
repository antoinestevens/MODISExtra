#' @title Get MODIS file URL's
#' @usage filesUrl(url)
#' @param url an url pointing to MODIS data for a given date
#' @note This produce the same results as \code{MODIS:::filesUrl} but uses the \code{rvest}
#' package instead of \code{RCurl}
#' @return a character \code{vector} of urls pointing to data locations
#' @seealso \code{MODIS:::filesUrl} in the MODIS package
#' @examples
#' \dontrun{
#' filesUrl(url = "http://e4ftl01.cr.usgs.gov/MOLT/MOD13A2.005/2015.05.25")
#  filesUrl(url = "http://e4ftl01.cr.usgs.gov/MOLT/MOD13A2.005/2000.03.18")
#' }
#'
#' @export

filesUrl <- function(url)
{
  # use now rvest instead of old RCurl

  if (substr(url,nchar(url),nchar(url))!="/")
    url <- paste0(url,"/")
  try(co <- rvest::html(url,httr::timeout(5)),silent=TRUE)
  if (class(co)[1]=="try-error"|as(co,"character")=="")
    return(NULL)
  fnames <- character(0)
  if (substring(url,1,4)=="http")
  { 
    co <- rvest::html_nodes(co,"a")
    if(!length(co)){  # try a second time (sometimes, the file list is not returned by the http request the first time it is GET)
      Sys.sleep(0.1)
      try(co <- rvest::html(url,httr::timeout(5)),silent=TRUE)
      if (class(co)[1]=="try-error")
        return(co)
      co <- rvest::html_nodes(co,"a")
    }
    co <- rvest::html_attr(co,"href")
    fnames <- co[stringr::str_detect(co,"hdf$")]
    fnames <- gsub(fnames,pattern="/",replacement="")
  } else
  {
    co <- stringr::str_split(co, ifelse(.Platform$OS.type=="unix","\n","\r\n"))[[1]]
    co <- stringr::str_split(" ",co)
    elim    <- grep(co,pattern="total")
    if(length(elim)==1)
      co <- co[-elim]
    fnames <- basename(sapply(co,function(x){x[length(x)]}))
    fnames <- gsub(fnames,pattern="/",replacement="")
  }

  return(fnames)
}
