#' comb accumulated run
#'
#' this is for submitting multiple accumulation run individually and combine them by one run
#' the function will correct the index number and make the combined file looks like output of one run
#'
#' @param path string. the input o02 file. must be provided
#' @return no return
#' @export
#' @import stringr magrittr
comb_accum<-function(path=NULL){
  if(is.null(path)){
    stop("please provide the path to o02 file")
  }
  lines=readLines(path)
  iout.ind=str_which(string=lines,pattern="iout_th")
  lines[iout.ind+1]=str_replace_all(string=lines[iout.ind+1],pattern="^\\s+\\d+",replacement=paste("     ",seq(from=0,by=1,to=length(iout.ind)-1)))
  listres=change_para("imc\\_rep",NA,infile=path,outfile=NA,type="show")
  numval=as.numeric(listres[["val"]])
  diffval=diff(numval)
  if(length(which(diffval<0))!=0){
    addon=max(numval)
    flagvec=cumsum(c(0,diffval<0))
    numval=numval+flagvec*addon
    imc_rep_ind=listres[["ind"]]
    lines[imc_rep_ind] %>% str_trim(string=.,side="both") %>%
                           str_split(string=.,simplify=TRUE,pattern="\\s+") -> tab
    tab[,2]=as.character(numval)
    lines[imc_rep_ind]=apply(tab,1,function(x){
      paste0("      ",paste(x,collapse="     "))
    })
  }
  cat(lines,sep="\n",file=paste0(path,"_addon"))
}
