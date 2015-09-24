#' @title Handles the in/output date used in the filtering/interpolation
#'
#' @description See \code{\link[MODIS]{orgTime}} in the MODIS package, which this function is a copy-paste, with one modification
#'              of the \code{nDays} argument which accepts now user-defined time series
#' @param files MODIS filenames, eg output of \code{runGdal}
#' @param nDays Integer, \code{Date} vector or character string indicating the time interval/time stamps for the output layers.
#'        Default is "asIn" that includes the exact input dates within the period selected using \code{begin} and \code{end}.
#'        Can also be nDays="1 month" or "1 week"
#' @param begin Default is from earliest input dataset. Here you can supply the begin date of the output
#' @param end Default to the end of the input dataset. Here you can specify the end date of the output
#'        (note, the exact end-date depends on \code{begin} and \code{nDays} argument.
#' @param pillow Number of days added on the beginning and on the end of a time serie.
#' @param pos1 Start position of date in the filename.
#' @param pos2 End position of date in the filename.
#' @param format How is the date formatted in the file, default expects: 'YYYYDDD' ("\%Y\%j").
#' @return See \code{\link[MODIS]{orgTime}} in the MODIS package
#' @author Matteo Mattiuzzi and Antoine Stevens
#' @export
orgTime <- function (files, nDays = "asIn", begin = NULL, end = NULL, pillow = 75,
                     pos1 = 10, pos2 = 16, format = "%Y%j")
{
  if (inherits(files, "Raster"))
    files <- names(files)

  files <- basename(files)

  allDates <- sort(extractDate(files, asDate = TRUE, pos1 = pos1, pos2 = pos2, format = format)$inputLayerDates)
  if(!length(allDates))
    stop("Problems with Date parsing in ORGTime!")
  datLim <- transDate(begin = begin, end = end)
  if (!is.null(begin)) {
    minOUT <- datLim$begin
    minIN <- minOUT - pillow
    minHAVE <- min(allDates[allDates >= minIN])
    if (as.character(nDays[1]) == "asIn")
      minIN <- minHAVE
  }
  else {
    minIN <- minOUT <- minHAVE <- min(allDates)
  }
  if (!is.null(end)) {
    maxOUT <- datLim$end
    maxIN <- maxOUT + pillow
    maxHAVE <- max(allDates[allDates <= maxIN])
    if (as.character(nDays[1]) == "asIn") {
      maxIN <- maxHAVE
    }
  }
  else {
    maxIN <- maxOUT <- maxHAVE <- max(allDates)
  }
  inputLayerDates <- allDates[allDates >= minHAVE & allDates <= maxHAVE]
  inDoys <- as.numeric(format(as.Date(inputLayerDates), "%j"))
  if (FALSE) {
    if (minIN < minHAVE) {
      if (as.numeric(minHAVE - minIN) <= pillow)
        warning("'begin'-date - 'pillow' is earlier by, ",
                as.numeric(minHAVE - minIN), " days, than the available input dates!\nPillow at the start of the time serie is reduced to ",
                pillow - as.numeric(minHAVE - minIN), " days!")

      else if (minOUT == minHAVE)
        warning("Is is not possible to use the pillow at the begin of the time series since there is no data available before 'begin'-date!")

    }
    if (maxIN > maxHAVE)
      warning("'end'-date + 'pillow' is later by, ", as.numeric(maxIN - max(inputLayerDates)), " days, than the available input dates!")
  }
  if (as.character(nDays[1]) == "asIn")
    outputLayerDates <- inputLayerDates[datLim$begin <= inputLayerDates & datLim$end > inputLayerDates]
  else if ((is.numeric(nDays)|is.character(nDays))&length(nDays)==1)
    outputLayerDates <- seq(minOUT, maxOUT, by = nDays)
  else
    outputLayerDates <- nDays

  t0 <- as.numeric(min(outputLayerDates, inputLayerDates)) - 1
  inSeq <- as.numeric(inputLayerDates) - t0
  outSeq <- as.numeric(outputLayerDates) - t0

  return(list(inSeq = inSeq, outSeq = outSeq, inDoys = inDoys,
              inputLayerDates = inputLayerDates, outputLayerDates = outputLayerDates,
              call = list(pos1 = pos1, pos2 = pos2, format = format,
                          asDate = TRUE, nDays = nDays, pillow = pillow)))
}
