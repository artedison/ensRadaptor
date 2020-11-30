#' summary input files
#'
#' this funciton summary exist content in i01 and i02 file
#' output exist species and reaction
#' print imporant parameters
#' these information can be used for modify ens.def file
#'
#' @param i01 string. the locaiton of i01. must be provided
#' @param i02 string. the locaiton of i02. must be provided
#' @return list. list of information contains in input files
#' @export
#' @import stringr magrittr
summary_input<-function(i01=NULL,i02=NULL){
  if(is.null(i01)||is.null(i02)){
    stop("please provide path to both i01 and i02 files")
  }
  list.exi=vector(mode="list")
  #i01
  nexpt=change_para("nexpt",NA,infile=i01,outfile=NA,type="show")
  ntxpt=change_para("ntxpt",NA,infile=i01,outfile=NA,type="show")
  nspec=change_para("nspec",NA,infile=i01,outfile=NA,type="show")
  nreac=change_para("nreac",NA,infile=i01,outfile=NA,type="show")
  ntime=change_para("ntime",NA,infile=i01,outfile=NA,type="show")
  lines=readLines(i01)
  lenlines=length(lines)
  lines %>% str_which(string=.,pattern="^\\s*namespec") %>%
            add(.,1) %>%
            extract(lines,.) %>%
            str_trim(string=.,side="both") %>%
            str_split(string=.,pattern="\\s+") %>%
            unlist(.) -> species
  list.exi[["species"]]=species
  lines %>% str_which(string=.,pattern="^\\s*namereac") %>%
            add(.,1) %>%
            extract(lines,.) %>%
            str_trim(string=.,side="both") -> reactions
  list.exi[["reactions"]]=reactions
  lines %>% str_which(string=.,pattern="^\\s*namereac") %>%
            add(.,2) %>%
            extract(lines,.) %>%
            str_extract(.,pattern="^\\s+\\d+\\s+\\d+\\s") %>%
            str_trim(string=.,side="both") %>%
            str_split(string=.,pattern="\\s+") %>%
            unlist(.) %>%
            as.numeric(.) %>%
            max(.) -> npart
  list.exi[["npart_x"]]=npart
  print(paste0("npart_x: ",npart))
  ##the table for range in total and initial condition for species
  ind_spec=str_which(string=lines,pattern="\\#\\s+Species\\s+control\\-\\s+and\\s+\\\\Theta\\-variables")
  ind_reac=str_which(string=lines,pattern="#\\s+Reaction\\s+control\\-\\s+and\\s+\\\\Theta\\-variables")
  indbound=c(ind_spec,ind_reac,lenlines)
  ind_theta=str_which(string=lines,pattern="\\s+1\\s+1\\s*$")
  indloc=c(0,1,3)
  blocknamevec=c("namespec","namereac")
  for(blocknamei in seq(length(blocknamevec))){
    blockname=blocknamevec[blocknamei]
    ind_blockname=str_which(string=lines,pattern=blockname)
    ind_blockname=ind_blockname[ind_blockname>indbound[blocknamei]&ind_blockname<indbound[blocknamei+1]]
    ind_blockname=c(ind_blockname,indbound[blocknamei])
    listpara=vector(mode="list")
    matpara_name=c()
    for(indblockele in seq(length(ind_blockname)-1)){
      blockstarti=ind_blockname[indblockele]
      blockendi=ind_blockname[indblockele+1]
      name=str_trim(lines[blockstarti+1],side="both")
      cat(paste0(name,"\n"))
      ind_theta_seles=ind_theta[ind_theta>blockstarti&ind_theta<blockendi]
      for(ind_theta_sele in ind_theta_seles){
        valvec=c()
        for(indloci in seq(length(indloc))){
          indlocele=indloc[indloci]
          lines[ind_theta_sele+indlocele] %>% str_trim(string=.,side="both") %>%
                            str_split(string=.,pattern="\\s+",simplify=TRUE) %>%
                            extract(.,3) %>% as.numeric(.) -> val
          valvec=c(valvec,val)
        }
        cat(paste0(valvec[1],"\t",valvec[2],"\t",valvec[3],"\n"))
        listpara[[length(listpara)+1]]=valvec
        matpara_name=c(matpara_name,name)
      }
    }
    matpara=as.data.frame(Reduce("rbind",listpara))
    matpara=cbind(matpara_name,matpara)
    rownames(matpara)=NULL
    colnames(matpara)=c("name","low","high","ini_guess")
    list.exi[[blockname]]=matpara
  }
  #i02
  ## the structure of experiment data
  lines=readLines(i02)
  exp_ind=str_which(string=lines,pattern="FF\\s+iexpt")
  exp_ind=c(exp_ind,length(lines))
  specname_ind=str_which(string=lines,pattern="FFF\\s+namespec")
  time_ind=str_which(string=lines,pattern="\\s*ndtin_txpt")
  exp_blocks=sapply(seq(length(exp_ind)-1),function(x){
         thisblock=exp_ind[x]
         nextblock=exp_ind[x+1]
         specname_block=specname_ind[specname_ind>thisblock&specname_ind<nextblock]
         time_block=time_ind[time_ind>thisblock&time_ind<nextblock]
         (specname_block+1) %>% extract(lines,.) %>%
                   str_trim(string=.,side="both") %>%
                   str_split(string=.,pattern="\\s+") %>%
                   unlist() -> specs
         (time_block+1) %>% extract(lines,.) %>%
                   str_trim(string=.,side="both") %>%
                   str_split(string=.,pattern="\\s+") %>%
                   unlist(.) %>% as.numeric(.) -> timepointnum
         names(timepointnum)=specs
         print(paste0("BLOCK: ",x))
         print(timepointnum)
         timepointnum
  })
  print(paste0("ALL:",sum(unlist(exp_blocks))))
  return(list.exi)
}

#' reaction summary
#'
#' path the reaction file location
#' warning on input equation file
#'
#' @param path string. path to the reaction list file. must be provided
#' @param enzpattern string. the pattern to search for the enzyme entity. must be provided
#' @return list. list containing enzyme information
#' @export
#' @import stringr
summary_reac<-function(path=NULL,enzpattern=NULL){
  if(is.null(path)){
    stop("please provide the input path")
  }
  if(is.null(enzpattern)){
    stop("please provide the searching pattern for enzyme")
  }
  list.reac.addon=read_reac(path)
  specs=unique(unlist(sapply(list.reac.addon,function(x){
    c(x[[2]],x[[3]])
  })))
  enzind=str_detect(string=specs,pattern=enzpattern)
  enz=unique(specs[enzind])
  print(paste0("enzymes: ",length(enz)))
  print(enz)
  compounds=unique(specs[!enzind])
  print(paste0("compounds: ",length(compounds)))
  print(compounds)
  reacts=unlist(sapply(list.reac.addon,function(x){
    x[[1]]
  }))
  print(paste0("reactions: ",length(reacts)))
  print(reacts)
  reacsdup=reacts[duplicated(reacts)]
  if(length(reacsdup)!=0){
    print("duplicated reactions:")
    print(unique(reacsdup))
  }
  tabenz=table(specs[enzind])
  if(max(tabenz)>2){
    print("duplicated enzymes:")
    print(names(tabenz[tabenz>2]))
  }
  reac.type=vector(mode="list")
  temp=sapply(list.reac.addon,function(x){
    sub=x[["subs"]]
    prod=x[["prods"]]
    name=x[["name"]]
    enzy=unique(c(sub[str_detect(string=sub,pattern=enzpattern)],prod[str_detect(string=prod,pattern=enzpattern)]))
    sub=sub[!str_detect(string=sub,pattern=enzpattern)]
    prod=prod[!str_detect(string=prod,pattern=enzpattern)]
    flag=sapply(reac.type,function(y){
      setequal(sub,y[["sub"]])&&setequal(prod,y[["prod"]])
    })
    ind=NULL
    if(length(flag)!=0){
      ind=which(flag)
    }
    if(length(ind)!=0){
      reac.type[[ind]][["name"]]<<-c(name,reac.type[[ind]][["name"]])
    }else{
      temp.list=list(sub=sub,prod=prod,name=name,enzy=enzy)
      reac.type[[length(reac.type)+1]]<<-temp.list
    }
  })
  list.res=list(enz=enz,compounds=compounds,reacs=reacts,reac.type=reac.type)
  return(list.res)
}

#' summary plot for o02
#'
#' summary plot for o02
#'
#' @param o02.data list. the struct from o02.reader. must be provided
#' @param dir.res string. the location for producing figures. must be provided
#' @param addonname string. the addon name. default ""
#' @param linethick bool. whether thicker lines shoulde be draw for figure. default FALSE
#' @return just plot no return
#' @export
summary_o02<-function(o02.data=NULL,dir.res=NULL,addonname="",linethick=FALSE){
  if(is.null(o02.data)){
    stop("please provide the input path")
  }
  if(is.null(dir.res)){
    stop("please provide the output result path")
  }
  ##length of modeling
  print(paste0("length: ",length(o02.data[["ids"]])))
  ##sweep~t
  tab=cbind(o02.data[["ids"]],o02.data[["chisq"]]*2)
  draw_sweep(tab,ylab="chi^2",
      loci=paste0(dir.res,"chi2-sweep.",addonname,".pdf"),
      linethick=linethick
  )
  ##hisogram of chi^2
  tab=o02.data[["chisq"]]*2
  dim(tab)=c(length(o02.data[["chisq"]]),1)
  draw_hist(tab,loci=paste0(dir.res,"chi2.distr.",addonname,".pdf"),xlab="chi^2")
  print(paste0("mean: ",mean(tab[,1])))
  ##acceptance rates
  vec=unlist(o02.data[["raccp"]])
  tab=vec
  tab=tab[tab>=0]
  tab[tab>1]=1
  dim(tab)=c(length(tab),1)
  draw_hist(tab,loci=paste0(dir.res,"acceptance.rate.distr.",addonname,".pdf"),xlab="acceprate")
  print("acceptance rate:")
  print(summary(vec))
  ##step size
  vec=unlist(o02.data[["fstp"]])
  tab=vec
  dim(tab)=c(length(vec),1)
  draw_hist(tab,loci=paste0(dir.res,"footstepsize.rate.distr.",addonname,".pdf"),xlab="stepsize")
  print("stepsize:")
  print(summary(vec))
}
