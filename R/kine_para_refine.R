#' refine kinetic parameters
#'
#' this function will take kinetic parameters structure and modify it
#' will calcualte range for all kinetic parameters and give default value for those without enough information
#'
#' @param list.res.refine list. the species non-specific data struct. must be provided
#' @param range.speci list. the species specific data struct. must be provided
#' @param extendrag numeric. extent to extend unknown range in kinetic parameters(when there is only one value in in record). must be provided
#' @return list. refined parameter list
#' @export
kine_para_refine<-function(list.res.refine=NULL,range.speci=NULL,extendrag=NULL){
  if(is.null(list.res.refine)||is.null(range.speci)){
    stop("please provide the two parameter lists")
  }
  if(is.null(extendrag)){
    stop("please provide the extension range")
  }
  ## the enzyme with kinetic parameter
  ecs=unique(c(names(list.res.refine[[1]]),names(list.res.refine[[2]]), ##for all species
               names(range.speci[[1]]),names(range.speci[[2]]))) ## for NC
  ##expected variance and value in general
  general.val=list(km=median(as.numeric(unlist(list.res.refine[["km"]]))),
                   kcat=median(unlist(list.res.refine[["kcat"]])),
                   km.var=extendrag*10^median(sapply(list.res.refine[["km"]],
                        function(x){
                          x[x<=0]=smallvalue
                          abs(log10(x[2]/x[1]))})),
                   kcat.var=extendrag*10^median(sapply(list.res.refine[["kcat"]],
                        function(x){
                          x[x<=0]=smallvalue
                          abs(log10(x[2]/x[1]))})))
  parameter.list=list(km=c(),kcat=c())#,ec=conc.ec)
  for(type in names(parameter.list)[1:2]){
    temp=sapply(ecs,simplify=FALSE,function(x){
            datatemp=range.speci[[type]][[x]]
            datalist=range.speci[[paste0(type,"list")]][[x]]
            ##if not in the nc specific struct then find if rom teh more general one
            if(is.null(datatemp)){
              datatemp=list.res.refine[[type]][[x]]
              datalist=list.res.refine[[paste0(type,"list")]][[x]]
            }
            # if(!is.null(datatemp)){
            #   if(datatemp[1]==datatemp[2]){
            #     datatemp[1]=datatemp[1]/general.val[[paste0(type,".var")]]
            #     datatemp[2]=datatemp[2]*general.val[[paste0(type,".var")]]
            #   }
            #   ## 0 is not allowed for lower limit
            #   lowrag=datatemp[1]
            #   if(lowrag==0){
            #     datatemp[1]=sort(unique(datalist))[2]
            #   }
            #   datatemp*conv.factor[[type]]
            # }
    })
    temp=temp[sapply(temp,length)!=0]
    parameter.list[[type]]=temp
    parameter.list[[type]][["default"]]=c(general.val[[type]]/general.val[[paste0(type,".var")]],general.val[[type]]*general.val[[paste0(type,".var")]])
  }
  return(parameter.list)
}
