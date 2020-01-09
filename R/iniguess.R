#' set initial guess
#'
#' this function will use the vector to initilize the inital guess for parameters
#'
#' @param input string. input file path. must be provided
#' @param output string. output file path. must be provided
#' @param valvec array. named vector of values for initial guess. block of value that wasn't inlcuded in the vector will be initilized as lowval + smallvalue. must be provided
#' @return no return just change the file
iniguess<-function(input=NULL,output=NULL,valvec=NULL){
  if(is.null(input)){
    stop("please provide input path")
  }
  if(is.null(output)){
    stop("please provide output path")
  }
  if(is.null(valvec)){
    stop("please provide the value array")
  }
  lines=readLines(input)
  indspec=str_which(string=lines,pattern="\\#\\s+Species\\s+control\\-\\s+and\\s+\\\\Theta\\-variables")
  lines_abov=lines[1:indspec]
  lines_below=lines[(indspec+1):length(lines)]
  indnames=str_which(string=lines_below,pattern="(namespec)|(namereac)")
  indnamesrag=c(indnames,length(lines))
  names=str_trim(string=lines_below[indnames+1],side="both")
  indtheta=str_which(string=lines_below,pattern="\\s+1\\s+1\\s*$")
  lines_below[indtheta] %>% str_trim(string=.,side="both") %>%
                    str_split(string=.,pattern="\\s+",simplify=TRUE) %>%
                    extract(,3) %>% as.numeric(.) -> lowval
  for(indnamesi in seq(length(indnames))){
    # print(indnamesi)
    name=names[indnamesi]
    ranghere=c(indnamesrag[indnamesi],indnamesrag[indnamesi+1])
    indtheta_here=indtheta[indtheta>ranghere[1]&indtheta<ranghere[2]]
    for(subtheta_indi in seq(length(indtheta_here))){
      # print(subtheta_indi)
      subtheta_ind=indtheta_here[subtheta_indi]
      valname=paste0(name,"_",as.character(subtheta_indi))
      if(valname%in%names(valvec)){
        newval=valvec[valname]
      }else{
        lowvalhere=lowval[indtheta>ranghere[1]&indtheta<ranghere[2]][subtheta_indi]
        newval=lowvalhere+smallvalue
      }
      lines_below[subtheta_ind+3]=str_replace(string=lines_below[subtheta_ind+3],pattern="[\\d\\.]+\\s*$",replacement=as.character(newval))
    }
  }
  cat(c(lines_abov,lines_below),sep="\n",file=output)
}
