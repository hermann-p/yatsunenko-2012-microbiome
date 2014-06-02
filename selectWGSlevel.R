selectWGSlevel <- function(theLevel, theDir="/home/pah11816/Yatsunenko/Data/processed_wgs/") {
  cat("reading ", theLevel, " stats from samples in ", theDir, "...\n")
  sampleList = dir(theDir)
  theData <- data.frame()
  for (sample in sampleList) {
    fileName <- paste(c(theDir, "/", sample, "/processed/999.done.", theLevel, ".stats"), collapse="")
    tmpData <- read.delim(fileName, header=FALSE)
    tmpValues <- tmpData[,2]
    names(tmpValues) <- tmpData[,1]
    i <- nrow(theData) + 1
    for (name in names(tmpValues)) {
      theData[i,name] <- tmpValues[name]
    }
  }
  theData
}
