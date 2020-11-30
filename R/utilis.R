#' str_which for multiple lines
#'
#' this function performs as str_which (no \n in pattern).
#' If there is \n in pattern, the function will try to match corresponding continuous line and return index of the last line matched.
#'
#' @param string string. input string to match with. array with each element for a line. must be provided
#' @param pattern string. the pattern to find. \\n is in the pattern. must be provided
#' @return array. the index location for the last line of the matching block
#' @export
#' @examples
#' string=c("aabb","cc","eeab","ce")
#' pattern="b\\nc"
#' multi_line_map(string=string,pattern=pattern)
#' pattern="eea"
#' multi_line_map(string=string,pattern=pattern)
#' pattern="eea\\n"
#' multi_line_map(string=string,pattern=pattern)
#' @import stringr
multi_line_map<-function(string=NULL,pattern=NULL){
  if(is.null(string)||is.null(pattern)){
    stop("please provide both the string and searching pattern")
  }
  pattern_parts=str_split(string=pattern,pattern="\\\\n")[[1]]
  list.match=vector(mode="list")
  vecmark=rep(NA,times=length(string))
  lenpart=length(pattern_parts)
  for(parti in seq(lenpart)){
    pattern_part=pattern_parts[parti]
    list.match[[parti]]=str_which(string=string,pattern=pattern_part)
    vecmark[list.match[[parti]]]=parti
  }
  blockind=c()
  if(lenpart==1){
    blockind=list.match[[1]]
  }else{
    vecmarki=1
    parti=1
    while(vecmarki<=length(vecmark)){
      vecmarkele=vecmark[vecmarki]
      if(!is.na(vecmarkele)&&vecmarkele==parti){
        parti=parti+1
      }else{
        parti=1
      }
      if(parti==lenpart+1){
        blockind=c(blockind,vecmarki)
        parti=1
      }
      vecmarki=vecmarki+1
    }
  }
  return(blockind)
}

#' record history in the history file
#'
#' record history in the history file
#'
#' @param histpath string. the history record file. must be provided
#' @param projpath string. the current working project directory. must be provided
#' @return no return just modify the file
#' @export
update_history<-function(histpath=NULL,projpath=NULL){
  if(is.null(histpath)||is.null(projpath)){
    stop("please provide the path for both record file and current working folder")
  }
  time=format(Sys.time(), "%H%M%S_%m%d%Y")
  string=paste0(time,"\t",projpath)
  cat(string,file=histpath,sep="\n",append=TRUE)
}

#' read reactions file
#'
#' load the reaction list file
#'
#' @param path string. the path to the reaction list file. must be provided
#' @return list. return the list containing c(reaction_name, substrates, products)
#' @import stringr magrittr
read_reac<-function(path=NULL){
  if(is.null(path)){
    stop("please provide the path for reaction file")
  }
  lines=readLines(path)
  list.reac.addon=sapply(lines,simplify=FALSE,function(x){
    x %>% str_trim(.) %>%
          str_split(string=.,pattern="\\s+") %>%
          extract2(.,1) -> slots
    name=str_replace_all(string=slots[3],pattern="\\:",replacement="")
    converter_ind=str_which(string=slots,pattern="\\<\\-\\>")
    subs=slots[4:(converter_ind-1)]
    prods=c()
    if(converter_ind+1<=length(slots)){
      prods=slots[(converter_ind+1):length(slots)]
      prods=prods[(!is.na(prods))&(prods!="")]
    }
    list(name=name,subs=subs,prods=prods)
  })
  return(list.reac.addon)
}

#' file based str_which
#'
#' file verion of str_which, used especially for large files
#'
#' @param path string. input file for the string. must be provided
#' @param pattern string. pattern to be matched. must be provided
#' @return array. the index(line number) of matched position
#' @export
#' @import stringr
file_str_which<-function(path=NULL,pattern=NULL){
  if(is.null(path)||is.null(pattern)){
    stop("please provide both the file path and searching pattern")
  }
  outputs=system(paste0("LC_ALL=C grep -n ","\'",pattern,"\' \'",path,"\'"),intern=TRUE)
  ind=as.numeric(str_split(outputs,pattern=":",simplify=TRUE)[,1])
  return(ind)
}

#' use relative locaiton between indx and indy to filter indy
#'
#' use relative locaiton between indx and indy to filter indy
#'
#' @param indx array. the reference index. must be provided
#' @param indy array. the index to be filtered. must be provided
#' @param rela string. the expected relative location between indx and indy. default ">="
#' @return array. the filtered indy
rela_ind_screen<-function(indx=NULL,indy=NULL,rela=">="){
  if(is.null(indx)||is.null(indy)){
    stop("please provide both index arrays")
  }
  if(rela==">="){
    indy_cho=sapply(indx,function(x){
      indy[indy>=x][1]
    })
  }else if(rela==">"){
    indy_cho=sapply(indx,function(x){
      indy[indy>x][1]
    })
  }else{
    indy_cho=sapply(indx,function(x){
      temp=indy[indy<x]
      temp[length(temp)]
    })
  }
  return(indy_cho)
}

#' constructed density frunction from a list of interval
#'
#' constructed density frunction from a list of interval
#'
#' @param list.intev list. list of parameter range. must be provided
#' @param logtrans bool. whether log transform the interval. TRUE log transformed. default FALSE
#' @return list. list of range distribution information
#' @export
densi_region<-function(list.intev=NULL,logtrans=FALSE){
  if(is.null(list.intev)){
    stop("please provide the parameter range list")
  }
  if(logtrans){
    bounds=log10(unlist(list.intev))
  }else{
    bounds=unlist(list.intev)
  }
  boundtype=rep(c("L","U"),times=length(bounds))
  ind=order(bounds)
  bounds.sort=bounds[ind]
  boundtype.sort=boundtype[ind]
  dens.pre=cumsum(ifelse(boundtype.sort=="L",1,-1))
  dens.pre=dens.pre[-length(dens.pre)]
  reglength=bounds.sort[-1]-bounds.sort[-length(bounds.sort)]
  scale=sum(dens.pre*reglength)
  dens=dens.pre/scale
  # plot(log10(bounds.sort[-length(bounds.sort)]),dens.pre,xlab="log10(x)",ylab="density")
  power2=bounds.sort^2
  power3=bounds.sort^3
  power2delta=power2[-1]-power2[-length(power2)]
  power3delta=power3[-1]-power3[-length(power3)]
  mean=sum(power2delta*dens)/2
  meanpower2=sum(power3delta*dens)/3
  sigma=sqrt(meanpower2-mean^2)
  list.res=list(densi=dens,reglength=reglength,bounds.sort=bounds.sort,
                mean=mean,sigma=sigma)
  return(list.res)
}

#' modified version of as.data.frame for rownames
#'
#' modified version of as.data.frame to remove the changed caused by make.name
#'
#' @param list list. list to be converted. must be provided
#' @return dataframe. will return t(as.data.frame(list))
#' @export
#' @examples
#' inputlist=list("1.34"=c(1,2,3),"8.5"=c(3,2,3))
#' as_data_frame_rowname(inputlist)
as_data_frame_rowname<-function(list=NULL){
  if(is.null(list)){
    stop("please provide the list")
  }
  names=names(list)
  dataframe=t(as.data.frame(list))
  rownames(dataframe)=names
  return(dataframe)
}

#' folder creating
#'
#' this function will create the folder if the folder doesn't exist
#'
#' @param dirpath string. path to the folder to create. must be provided
#' @return no return just create the folder
#' @export
foldcreate<-function(dirpath=NULL){
  if(is.null(dirpath)){
    stop("please provide directory path")
  }
  if(!dir.exists(dirpath)){
    system(paste0("mkdir \'",dirpath,"\'"))
  }
}

#' count number of lines in the file
#'
#' this function will count number of lines in the file
#'
#' @param path string. the file location. must be provided
#' @return int. line number
#' @export
#' @import stringr
linecount<-function(path=NULL){
  if(is.null(path)){
    stop("please provide file path")
  }
  linenum=system(paste0("wc -l \"",path,"\""),intern=TRUE)
  linenum=as.numeric(str_trim(string=str_extract(string=linenum,pattern="^\\s*\\d+\\s+"),side="both"))
  return(linenum)
}
