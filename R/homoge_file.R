#' modify output script for consistency
#'
#' makesure the output analysis script agreed with the original input file
#'
#' @param localpath string. the local working directory. must be provided
#' @param foldname string. the name of the project folder. must be provided
#' @return print out result no return
#' @export
homoge_file<-function(localpath=NULL,foldname=NULL){
  if(is.null(localpath)){
    stop("please provide directory path")
  }
  if(is.null(foldname)){
    stop("please provide project folder")
  }
  outputfiles=list.files(path=localpath,pattern="\\.+output.R")
  if(length(outputfiles)>1){
    stop("multiple local output files")
  }else{
    lines=readLines(paste0(localpath,outputfiles))
    ind_targ=str_which(string=lines,pattern="^foldname=")
    lines[ind_targ]=paste0("foldname=\"",foldname,"\"")
  }
  cat(lines,sep="\n",file=paste0(localpath,outputfiles))
}
