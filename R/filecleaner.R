#' clean o02 file
#'
#' this script is before fetching information out of output file from ensemble
#' currently only deal with RESTART run which cause some block rearrangement in o02 file
#'
#' @param path string. input file. must be provided
#' @return string. the path to the cleaned file
#' @export
#' @import stringr
o02_clean<-function(path=NULL){
  if(is.null(path)){
    stop("please provide input path")
  }
  outputpath=paste0(path,".cleaned")
  lines=readLines(path)
  start_ind=str_which(string=lines,pattern="iout_th")
  RESTART_ind=str_which(string=lines,pattern="RESTART")
  if(length(RESTART_ind)>0){
    bfind=start_ind[start_ind<RESTART_ind]
    afind=start_ind[start_ind>RESTART_ind]
    remrange=c(bfind[length(bfind)],afind[1])
    lines=lines[c(1:(remrange[1]-1),remrange[2]:length(lines))]
  }
  cat(lines,sep="\n",file=outputpath)
  return(outputpath)
}
