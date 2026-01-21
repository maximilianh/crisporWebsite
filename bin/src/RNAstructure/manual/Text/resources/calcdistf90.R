calcdistf90 <- function(pairdata) {
  # output: distance matrix 
  # input:  pairing_data matrix
  dyn.load('calcdist.so')
  # get number of nucleotides and structures
  size <- dim(pairdata)
  nucl <- size[2]
  struct <- size[1]
  distsize <- struct*struct
  dist <- .Fortran("calcdist",nuc=as.integer(nucl),struc=as.integer(struct),pairing=as.integer(pairdata),distance=integer(distsize))
  return(matrix(dist$distance,struct,struct))
}