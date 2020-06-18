#' Editing or showing ensemble input files (ens.i01) for single term
#'
#' The function will edit or show the value of specific parameters
#'
#' @param term string. the term that need to be changed. must be provided
#' @param value array, dim(array)=1. the set value. default NULL
#' @param infile string. input file path. must be provided
#' @param outfile string. output file path. default NULL
#' @param type string. "show"(show the parameter, default), "edit"(edit the value)
#'         if there is mulitple locations for the corresponding term, only show value is possible
#'         and both location/index and value are returned.
#' @param indrep array. the index of the one that need to be changed if there are multiple (edit). default 1
#' @return If type="edit", no value is returned and just change the file correspondingly
#'         If type="show" OR length(uncomm)!=1&&indrep==1, a list of value and line number will be returned
#' @seealso [change_block()] for modifying one block
#' @export
change_para<-function(term=NULL,value=NULL,infile=NULL,outfile=NULL,type="show",indrep=1){
  if(is.null(term)){
    stop("please provide the term to change")
  }
  if(is.null(infile)){
    stop("please provide full path for input at least")
  }
  lines=readLines(infile)
  comment_ind=str_detect(string=lines,pattern="^#")
  ind=str_which(string=lines,pattern=term)
  uncomm=which(!comment_ind[ind])
  # if(is.na(outfile)||infile==outfile){
  #   print("pay attention to file name!")
  # }
  if(length(uncomm)!=1&&indrep==1){
    print("the term is duplicated!\n")
    print("only \'show\' value is possible\n")
    ind=ind[!comment_ind[ind]]
    titles=str_split(string=str_trim(lines[ind],side="both"),pattern="\\s+")[[1]]
    colind=str_which(string=titles,pattern=term)
    ind=ind+1
    strs=str_split(string=str_trim(lines[ind],side="both"),pattern="\\s+",simplify=TRUE)
    res=strs[,colind]
    list=list(val=res,ind=ind)
    # print(paste0(term," : ",strs[colind]))
    return(list)
  }else{
    ind=ind[!comment_ind[ind]]
    titles=str_split(string=str_trim(lines[ind],side="both"),pattern="\\s+")[[1]]
    colind=str_which(string=titles,pattern=term)
    ind=ind+1
    ind=ind[indrep]
    strs=str_split(string=str_trim(lines[ind],side="both"),pattern="\\s+")[[1]]
    if(type=="show"){
      print(paste0(term," : ",strs[colind]))
      res=strs[colind]
      list=list(val=res,ind=ind)
      return(list)
    }else{
      strs[colind]=value
      lines[ind]=paste(" ",strs,collapse="    ")
      cat(lines,sep="\n",file=outfile)
    }
  }
}
#' Editing or showing specific blocks (ens.i01)
#'
#' The function will edit or show the value of specific blocks
#'
#' @param paternlist list. indicates the location for block changing
#'          ->pattern0: the start of the block to be matched
#'          ->pattern1: the end of the block to be matched.
#'          ->shift: a vector. if length(shift)==2. Index shift will be added on pattern0 and pattern1
#'                   if length(shift)==1. The shift is the same for both pattern0 and pattern1
#'          must be provided.
#' @param contentlist list. the value to be changed
#'          ->file: the file to read from
#'          ->ilinerag: the index range of lines as input
#'          Default NULL
#' @param infile string. input file path. must be provided
#' @param outfile string. output file path. default NULL
#' @param type string. "show"(show the parameter, default), "edit"(edit the value)
#'        if there is mulitple value for the corresponding term, only show value is possible
#' @return If type="edit", no value is returned and just change the file correspondingly
#'         If type="show", return the show part
#' @seealso [change_para()] for modifying one parameter
#' @export
change_block<-function(paternlist=NULL,contentlist=NULL,infile=NULL,outfile=NULL,type="show"){
  if(is.null(paternlist)){
    stop("please provide correct searching pattern")
  }
  if(is.null(infile)){
    stop("please provide full path for input at least")
  }
  # if(is.na(outfile)||infile==outfile){
  #   print("pay attention to file name!")
  # }
  lines=readLines(infile)
  lenlines=length(lines)
  shifts=paternlist[["shift"]]
  ind0=multi_line_map(string=lines,pattern=paternlist[["pattern0"]])+shifts[1]
  ind1=multi_line_map(string=lines,pattern=paternlist[["pattern1"]])+shifts[2]
  if(type=="show"){
    # print(paste0(term," : ",lines[ind0:ind1]))
    return(lines[ind0:ind1])
  }else{
    line_replace=readLines(contentlist[["file"]])
    ilinerag=contentlist[["ilinerag"]]
    linesnew=c()
    if(ind0>1){
      linesnew=c(linesnew,lines[1:(ind0-1)])
    }
    linesnew=c(linesnew,line_replace[ilinerag[1]:ilinerag[2]])
    if(ind1<lenlines){
      linesnew=c(linesnew,lines[(ind1+1):lenlines])
    }
    cat(linesnew,sep="\n",file=outfile)
    # return(linesnew)
  }
}

#' editing or showing ens.def files
#'
#' the function is for editing or showing specific parameters
#'
#' @param term string. the value that need to be changed. must be provided
#' @param value array, dim(array)=1. the value. default NULL
#' @param infile string. input file path. must be provided
#' @param outfile string. output file path. default NULL
#' @param type string. "show"(show the parameter. default), "edit"(edit the value)
#' @return If type="edit", no value is returned and just change the file correspondingly
#'         If type="show", a list of value will be returned
#' @export
change_def<-function(term=NULL,value=NULL,infile=NULL,outfile=NULL,type="show"){
  if(is.null(term)){
    stop("please provide the term to change")
  }
  if(is.null(infile)){
    stop("please provide full path for input at least")
  }
  lines=readLines(infile)
  comment_ind=str_detect(string=lines,pattern="^!")
  patternwrap=paste0("parameter\\(",term,"\\=[^\\)]+\\)")
  term_nonescape=str_replace(string=term,pattern="\\\\",replacement="")
  repalcewrap=paste0("parameter(",term_nonescape,"=",value,")")
  ind=str_which(string=lines,pattern=patternwrap)
  uncomm=which(!comment_ind[ind])
  # if(is.na(outfile)||infile==outfile){
  #   print("pay attention to file name!")
  # }
  ind=ind[!comment_ind[ind]]
  if(type=="show"){
    matchpart=lines[ind]
    print(matchpart)
    matchpart %>% str_extract(string=.,pattern="\\=[^\\)]+") %>%
                  str_replace(string=.,pattern="\\=",replacement="") %>%
                  as.numeric(.) -> paranum
    return(paranum)
  }else{
    lines[ind]=str_replace(string=lines[ind],pattern=patternwrap,replacement=repalcewrap)
    cat(lines,sep="\n",file=outfile)
  }
}

#' editing or showing ens.sh files
#'
#' the function is for editing or showing specific parameters
#'
#' @param term string. the value that need to be changed. must be provided
#' @param value array, dim(array)=1. the value. default NULL
#' @param infile string. input file path.  must be provided
#' @param outfile string. output file path. default NULL
#' @param type string. "show"(show the parameter. default), "edit"(edit the value)
#' @return If type="edit", no value is returned and just change the file correspondingly
#'         If type="show", a list of value will be returned
#' @export
change_sh<-function(term=NULL,value=NULL,infile=NULL,outfile=NULL,type="show"){
  if(is.null(term)){
    stop("please provide the term to change")
  }
  if(is.null(infile)){
    stop("please provide full path for input at least")
  }
  lines=readLines(infile)
  patternwrap=paste0(term,"\\=.+")
  repalcewrap=paste0(term,"=",value)
  ind=str_which(string=lines,pattern=patternwrap)
  # if(is.na(outfile)||infile==outfile){
  #   print("pay attention to file name!")
  # }
  if(type=="show"){
    matchpart=lines[ind]
    print(matchpart)
    matchpart %>% str_extract(string=.,pattern="\\=[^\\)]+") %>%
                  str_replace(string=.,pattern="\\=",replacement="") -> paranum
    return(paranum)
  }else{
    lines[ind]=str_replace(string=lines[ind],pattern=patternwrap,replacement=repalcewrap)
    cat(lines,sep="\n",file=outfile)
  }
}
