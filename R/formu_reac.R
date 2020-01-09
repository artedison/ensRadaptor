#' formualte transcription and translation
#'
#' this is funciton is for producing transcription and translation reaction from a gene/enzyme list
#'
#' @param enz array. the enzyme list. must be provided
#' @param path string. the result location. must be provided
#' @return no return
#' @export
formu_reac<-function(enz=NULL,path=NULL){
  if(is.null(enz)){
    stop("please provide enzyme list")
  }
  if(is.null(path)){
    stop("please provide output path")
  }
  lines<-sapply(seq(length(enz)),function(x){
    enzyme=enz[x]
    transcr=paste0("Reaction (",x,") ",enzyme,"_trans: ",
                enzyme,"_1 <-> ",enzyme,"_1 ",enzyme,"_r")
    transl=paste0("Reaction (",x,") ",enzyme,"_transl: ",
                enzyme,"_r <-> ",enzyme,"_r ",enzyme)
    decayr=paste0("Reaction (",x,") ",enzyme,"_decay: ",
                enzyme," <-> ")
    decayl=paste0("Reaction (",x,") ",enzyme,"_r_decay: ",
                enzyme,"_r <-> ")
    c(transcr,transl,decayr,decayl)
  })
  cat(lines,sep="\n",file=path)
}
