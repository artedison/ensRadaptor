#' formualte the dataset (i02)
#'
#' formualte the dataset (i02)
#'
#' @param listdata list. the list contains time and value
#' @param classnum int. the start class number to be considered, a new class number will be return
#' @param output string. ouput the add part of experiment measrument
#' @param dir.data string. the working/data directory
#' @param modified.file string. addon file name
#' @param zspec_wid numeric. the zspec_wid for new input data, the default is 0.4
#' @param sameclass int. whether different scale class for different input measurement,
#'        the default is 0 (use different class number). 1 use the same class number
#' @return int. the new class number to be used
#' @export
template_exp_data<-function(listdata,classnum,output,dir.data,modified.file,zspec_wid=0.4,sameclass=0){
  ##find the used class type
  load(paste0(dir.lib,"default.measurement.RData"))
  classini=classnum
  exp_output=sapply(names(listdata),function(x){
    data=listdata[[x]]
    tempmat=data[,2,drop=FALSE]
    timeseq=data[,1]
    len=dim(tempmat)[1]
    tab=data.frame(ldtin=seq(len),
                   time_xpt=timeseq,
                   zspec_xpt=c(t(tempmat)),
                   zspec_wid=rep(NA,times=len),#it can be changed
                   zscal=rep(1.0,times=len),
                   zproc=rep(2,times=len),
                   lscal=rep(classini,times=len)
                )
    if(length(zspec_wid)==1&&class(zspec_wid)=="numeric"){
      tab[,"zspec_wid"]=rep(zspec_wid,times=len)
    }else{
      vec.rela.uncern=zspec_wid[[x]]
      if(length(vec.rela.uncern)==1){
        tab[,"zspec_wid"]=rep(vec.rela.uncern,times=len)
      }else{
        tab[,"zspec_wid"]=vec.rela.uncern
      }
    }
    if(!sameclass){
      classini<<-classini+1 # uncomment this line if you need all input from the measurement to be different scale classes
    }
    string=c(specline_exp,x,expnumline_exp,len,
             paste0(" ",paste(titledata_exp,collapse=" ")),
             unlist(apply(tab,1,function(y){
               paste(y,collapse=" ")
             })))
     output=paste(string,collapse="\n ")
  })
  # classini=classini+1 #uncomment this line if you need all input from the measrument to be the same class
  cat(exp_output,file=output,sep="\n   \n")
  # this part are deleted as the good thing might be input it into the corresponding part
  return(classini)
}
