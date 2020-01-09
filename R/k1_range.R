#' k1 range production
#'
#' this is for producing a constraint enough k1 range to start the training
#' based on the assumption that each enzyme will have similar region size for kinetic parameter
#'
#' @param list.reac.addon list. the list of reactions. must be provided
#' @param list.kine kinetic list. parameter list. must be provided
#' @param varn numeric. how variance the range can be, the range is defined as mean +- varn*sigma
#'        default 1
#' @return list. formualted k1 range
k1_range<-function(list.reac.addon=NULL,list.kine=NULL,varn=1){
  if(is.null(list.reac.addon)){
    stop("please provide reaction list")
  }
  if(is.null(list.kine)){
    stop("please provide parameter list")
  }
  k1list=sapply(list.reac.addon,simplify=FALSE,function(x){
    name=x[["name"]]
    para=sapply(c("km","kcat"),simplify=FALSE,function(type){
      paratemp=list.kine[[type]][[name]]
      if(is.null(paratemp)){
        temp=list.kine[[type]][["default"]]
      }else{
        temp=paratemp
      }
      if(temp[1]==temp[2]){
        temp[1]=temp[1]/2
        temp[2]=temp[2]*2
      }
      temp[temp==0]=smallvalue
      temp
    })
    k1ragraw=c(min(para[["kcat"]])/max(para[["km"]]),max(para[["kcat"]])/min(para[["km"]]))
  })
  list.res<-densi_region(k1list,TRUE)
  sigma=list.res[["sigma"]]
  k1list_constr=sapply(k1list,simplify=FALSE,function(x){
    logx=log10(x)
    logmean=mean(logx)
    min=max(logmean-sigma*varn,logx[1])
    max=min(logmean+sigma*varn,logx[2])
    newx=c(10^min,10^max)
  })
  return(k1list_constr)
}
