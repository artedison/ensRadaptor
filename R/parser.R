#' parser for name information
#'
#' this function parse the part wtih species and reaction information
#' it will return names of theta.
#' only for o02 file 'namespec' and 'namereac' block
#'
#' @param lines array. specific lines to parse. must be provided
#' @param ind array. the index to define start searching point. must be provided
#' @param pattern string. the pattern for titile. 'namespec' or 'namereac'. must be provided
#' @return array. parameter name list
nameparser<-function(lines=NULL,ind=NULL,pattern=NULL){
  if(is.null(lines)){
    stop("please provide lines to parse")
  }
  if(is.null(ind)){
    stop("please provide the searching index array")
  }
  if(is.null(pattern)){
    stop("please provide pattern to match")
  }
  lineshere=lines[(ind[1]+1):(ind[2]-1)]
  lineshere %>% str_detect(string=.,pattern=pattern) %>%
              which(.) -> block_ind
  block_ind=c(block_ind,ind[2]-ind[1]-1)
  unlist(sapply(seq(length(block_ind)-1),function(x){
    loci=block_ind[x]
    locinext=block_ind[x+1]
    name=str_trim(string=lineshere[loci+1],side="both")
    reps=locinext-loci-3
    paste(rep(name,times=reps),seq(reps),sep="_")
  }))
}

#' parse the parameter numbers for theta speceis/parameters
#'
#' o02 file blocks of species parameter, reaction parameter, acceptance rate information and step size information.
#'
#' @param lines_models array. lines for searching. must be provided
#' @param this_ind array. define the range start. must be provided
#' @param next_ind array. define the range end. must be provided
#' @param name array. name of the theta vector. default NULL. no names
#' @return list. the named parameter list
numparse<-function(lines_models=NULL,this_ind=NULL,next_ind=NULL,name=NULL){
  if(is.null(lines_models)){
    stop("please provide lines to parse")
  }
  if(is.null(this_ind)||is.null(next_ind)){
    stop("please provide the searching index range")
  }
  sapply(seq(length(this_ind)),simplify=FALSE,function(x){
                loci=this_ind[x]+1
                lociend=next_ind[x]-1
                linestemp=lines_models[loci:lociend]
                linestemp=linestemp[!str_detect(string=linestemp,pattern="^\\s+$")]
                linepas=paste(linestemp,collapse="")
                linepas=str_trim(linepas,side="both")
                vec=str_split(string=linepas,pattern="\\s+")[[1]]
                vec=as.numeric(str_replace_all(string=vec,pattern="D",replacement="e"))
                names(vec)=name
                vec
  })
}
#' table data parse
#'
#' this funciton parse the table in input and output file and return specific column
#' space will be treated as separater
#'
#' @param lines_models array. the lines to fetch. must be provided
#' @param rowrange array. row range in the table. must be provided
#' @param colrange array. column range in the table. must be provided
#' @return array. the selected part of the table
#' @export
#' @examples
#' lines_models=c("1 2 3 4 5","6 7 8  9 10","11 12 13\t14 15")
#' tabparse(lines_models,c(1,3),c(2,4))
tabparse<-function(lines_models=NULL,rowrange=NULL,colrange=NULL){
  if(is.null(lines_models)){
    stop("please provide lines to parse")
  }
  if(is.null(rowrange)||is.null(colrange)){
    stop("please provide the index range for both row and column")
  }
  lines_models[rowrange]%>% str_trim(string=.,side="both") %>%
            str_replace_all(string=.,pattern="D",replacement="e") %>%
            str_replace_all(string=.,pattern="\\s+",replacement="\t") %>%
            read.table(text=.,sep="\t",header=FALSE) ->tab
  return(tab[,colrange])
}
