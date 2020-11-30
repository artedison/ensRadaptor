#' reader for o02 file
#'
#' file reader for o02 file
#'
#' @param path string. input file. must be provided
#' @param flagsweep bool. TRUE or FALSE.
#'        TRUE: only the end of sweep part are considered for final reasult.
#'        FALSE: output all blocks. Default FALSE
#' @return list. list containing information in o02
#' @export
#' @import stringr
o02_reader<-function(path=NULL,flagsweep=FALSE){
  if(is.null(path)){
    stop("please provide input path")
  }
  lines=readLines(path)
  spec_ind=str_which(string=lines,pattern="^\\s+species\\s+Theta\\s+parameters$")
  reac_ind=str_which(string=lines,pattern="^\\s+reaction\\s+Theta\\s+parameters$")
  output_ind=str_which(string=lines,pattern="^iout_th\\s+jclC\\s+jclO")
  specs=nameparser(lines,ind=c(spec_ind,reac_ind),pattern="namespec")
  reacs=nameparser(lines,ind=c(reac_ind,min(output_ind)),pattern="namereac")
  lines_models=lines[min(output_ind):(length(lines))]
  start_ind=str_which(string=lines_models,pattern="iout_th")
  chisq_ind=str_which(string=lines_models,pattern="chisq_th")
  repacc_ind=str_which(string=lines_models,pattern="imc_rep")
  theta_spec_ind=str_which(string=lines_models,pattern="cspec_th")
  theta_reac_ind=str_which(string=lines_models,pattern="rreac_th")
  raccp_swo_ind=str_which(string=lines_models,pattern="raccp_swo")
  fstp_swo_ind=str_which(string=lines_models,pattern="fstp_swo")
  ch_ind=str_which(string=lines_models,pattern="[:alpha:][:alpha:]+")
  if(flagsweep){
    end_ind=str_which(string=lines_models,pattern="End\\s+Of\\s+Sweep")
    end_ind=end_ind[-1]
    ##use relative location to screen other inds
    indlist=c("start_ind","chisq_ind","theta_spec_ind","theta_reac_ind","raccp_swo_ind","fstp_swo_ind","repacc_ind")
    for(indele in indlist){
      assign(indele,rela_ind_screen(indx=end_ind,indy=get(indele),rela=">="))
    }
  }
  ids=as.numeric(str_split(string=str_trim(lines_models[start_ind+1],side="both"),
                pattern="\\s+",simplify=TRUE)[,1])
  chisq=str_split(string=str_trim(lines_models[chisq_ind+1],side="both"),
                pattern="\\s+",simplify=TRUE)[,1]
  chisq=as.numeric(str_replace_all(string=chisq,pattern="D",replacement="e"))
  numvec=str_split(string=str_trim(lines_models[repacc_ind+1],side="both"),
                pattern="\\s+",simplify=TRUE)
  nrep=as.numeric(str_replace_all(string=numvec[,2],pattern="D",replacement="e"))
  nacc=as.numeric(str_replace_all(string=numvec[,4],pattern="D",replacement="e"))
  theta_spec_ini=numparse(lines_models,this_ind=theta_spec_ind,next_ind=theta_reac_ind,name=specs)
  theta_reac_para=numparse(lines_models,this_ind=theta_reac_ind,next_ind=raccp_swo_ind,name=reacs)
  raccp_swo=numparse(lines_models,this_ind=raccp_swo_ind,next_ind=fstp_swo_ind,name=NULL)
  ch_ind_cho=rela_ind_screen(indx=fstp_swo_ind,indy=ch_ind,rela=">")
  # ch_ind_cho=sapply(fstp_swo_ind,function(x){
  #   ch_ind[ch_ind>x][1]
  # })
  if(is.na(ch_ind_cho[length(ch_ind_cho)])){
    ch_ind_cho[length(ch_ind_cho)]=length(lines_models)
  }
  fstp_swo=numparse(lines_models,this_ind=fstp_swo_ind,next_ind=ch_ind_cho,name=NULL)
  o02.data=list(ids=ids,chisq=chisq,nrep=nrep,nacc=nacc,theta_spec_ini=theta_spec_ini,
                theta_reac_para=theta_reac_para,
                raccp=raccp_swo,fstp=fstp_swo
              )
  return(o02.data)
}
#' reader for o03 file
#'
#' file reader for o03 file
#'
#' @param pathlist list. input file list for o03. must be provided
#' @param ioutselec int. the iout_id to be read. must be provided
#' @return array. array of time series data. time x species x #blocks
#' @export
#' @import stringr magrittr
o03_reader<-function(pathlist=NULL,ioutselec=NULL){
  if(is.null(pathlist)){
    stop("please provide input path")
  }
  if(is.null(ioutselec)){
    stop("please provide the index list")
  }
  ####combine the files in the pathlist
  input=pathlist[[1]]
  ##was used for the case the run was run in multiple section(restart and continued) so there was multiple o03 file
  if(length(pathlist)>1){
    system("rm .tempt.lines")
    paste0("cat \'",paste0(pathlist,sep="\' \'"),"\' > .tempt.lines")
    input=".tempt.lines"
  }
  ####block index and cut files
  blocksize=4000 ##the # of block each time readin
  iout_ind=file_str_which(input,"iout_th")
  ####find the first match of block 1
  prevpath=Sys.getenv("PATH")
  Sys.setenv(PATH=paste(prevpath,"/Users/yuewu/anaconda3/bin/",sep=":"))
  on.exit(Sys.setenv(PATH=prevpath),add=TRUE)
  outputs=system(paste0("LC_ALL=C pcregrep -M -n ","\'","iout_th\\n+\\s+1\\s*$","\' \'",input,"\'"),intern=TRUE) ## pcregrep use perl notation in regex so + instead of \+
  iout_1_ind=as.numeric(str_split(outputs[1],pattern=":",simplify=TRUE)[,1])
  iout_ind=iout_ind[iout_ind>(iout_1_ind-1)]
  nlines=system(paste0("wc -l \'",input,"\'"),intern=TRUE)
  nlines %>% str_trim(string=.,side="both") %>%
             str_split(string=.,pattern="\\s+",simplify=TRUE) %>%
             extract(.,,1) %>% as.numeric(.) -> nlines
  nblock=floor(length(iout_ind)/blocksize)
  blockseq=c(iout_ind[1]-1,iout_ind[seq(from=blocksize,by=blocksize,length.out=nblock-1)]-1,nlines)
  ####calculate the array size and initialization
  namespec_ind=file_str_which(input,"namespec")
  specind=namespec_ind[namespec_ind<iout_ind[2]&namespec_ind>iout_ind[1]]
  nspec=length(specind)
  nsweep=length(ioutselec)
  space_ind=file_str_which(input,"^\\s*$")
  space_ind=space_ind[space_ind>specind[1]][1]
  linesspike=read_lines(file=input,skip=specind[1]+2,n_max=space_ind-specind[1]-2)
  linesspike %>% paste0(.,collapse="") %>%
                 str_trim(string=.,side="both") %>%
                 str_split(string=.,pattern="\\s+") %>%
                 extract2(.,1) %>%
                 length(.) -> ntime
  ####x time, y species, z #blocks
  array_data=array(NA,c(ntime,nspec,nsweep))
  ####array dim name on species name
  lines_block=read_lines(file=input,skip=iout_ind[1],n_max=iout_ind[2]-iout_ind[1]-1)
  spec_ind_block=str_which(lines_block,pattern="namespec")
  specnames=str_trim(string=lines_block[spec_ind_block+1],side="both")
  lines_block[spec_ind_block-1] %>% str_trim(string=.,side="both") %>%
                                    str_split(string=.,pattern="\\s+",simplify=TRUE) %>%
                                    extract(.,,1) -> expid
  specnames=paste(specnames,"_",expid,sep="")
  dimnames(array_data)[[2]]<-specnames
  startind=1
  for(blocki in seq(length(blockseq)-1)){
    range=blockseq[blocki:(blocki+1)]
    endind=range[2]-range[1]+1
    lines=read_lines(file=input,skip=range[1],n_max=endind)
    iout_ind_in=c(str_which(lines,pattern="iout_th"),endind)
    lines[iout_ind_in+1] %>% str_trim(string=.,side="both") %>%
                             as.numeric(.) -> iout_num
    selec_ind=which(iout_num %in% ioutselec)
    iout_ind_in=c(iout_ind_in[selec_ind],iout_ind_in[selec_ind[length(selec_ind)]+1])
    namespec_ind_in=str_which(lines,pattern="namespec")
    # xspec_ind_in=str_which(lines,pattern="xspec_o")
    space_ind_in=sort(str_which(lines,pattern="^\\s*$"))
    temparray=array(NA,c(ntime,nspec,length(iout_ind_in)-1))
    dimnames(temparray)[[2]]<-specnames
    for(iout_i in seq(length(iout_ind_in)-1)){
      blockrange=c(iout_ind_in[iout_i],iout_ind_in[iout_i+1])
      namespec_ind_in_sub=namespec_ind_in[namespec_ind_in>blockrange[1]&namespec_ind_in<blockrange[2]]
      specnames_block=str_trim(string=lines[namespec_ind_in_sub+1],side="both")
      lines[namespec_ind_in_sub-1] %>% str_trim(string=.,side="both") %>%
                                       str_split(string=.,pattern="\\s+",simplify=TRUE) %>%
                                       extract(.,,1) -> expid_in
      specnames_block=paste(specnames_block,"_",expid_in,sep="")
      specnames_block_convind=order(match(specnames,specnames_block))
      namesepc_seq=seq(length(specnames_block_convind))
      specnames_block_matcind=specnames_block_convind[namesepc_seq]
      specnames_from=namespec_ind_in_sub[specnames_block_matcind]+3
      blockrange=range(specnames_from)
      spacblocrange=c(which(space_ind_in>blockrange[1])[1],which(space_ind_in>blockrange[2])[1])
      space_ind_in_block=space_ind_in[spacblocrange[1]:spacblocrange[2]]
      specnames_to=space_ind_in_block[sapply(specnames_from,function(x){
        which(space_ind_in_block>x)[1]
      })]-1
      strinput=sapply(seq(length(specnames_to)),function(x){
        paste0(lines[specnames_from[x]:specnames_to[x]],collapse="")
      })
      strinput %>%  str_trim(string=.,side="both")%>%
                    str_split(string=.,pattern="\\s+",simplify=TRUE) %>%
                    t() -> strinput_tempt
      strinput_tempt %>% str_replace_all(string=.,pattern="D",replacement="E") %>%
          as.numeric(.) -> strinput_tempt2
      dim(strinput_tempt2)=dim(strinput_tempt)
      temparray[,specnames_block_matcind,iout_i]=strinput_tempt2
      # for(namespec_i in seq(length(specnames_block_convind))){
      #   specnames_block_matcind=specnames_block_convind[namespec_i]
      #   specnames_from=namespec_ind_in_sub[namespec_i]+3
      #   specnames_to=space_ind_in[space_ind_in>specnames_from][1]-1
      #   lines[specnames_from:specnames_to] %>% paste0(.,collapse="") %>%
      #                   str_trim(string=.,side="both")%>%
      #                   str_split(string=.,pattern="\\s+") %>%
      #                   extract2(.,1) %>%
      #                   str_replace_all(string=.,pattern="D",replacement="E") %>%
      #                   as.numeric(.) -> temparray[,specnames_block_matcind,iout_i]
      # }
    }
    array_data[,specnames_block_matcind,startind:(startind+length(iout_ind_in)-2)]=temparray
    startind=startind+length(iout_ind_in)-1
  }
  return(array_data)
}

#' read the prior information
#'
#' read prior range of the enyzme initial concentration and kinetic parameters
#' initial value for enzyme is not correct if you are using dependence
#' there is some manual modification in i01 file, and so
#' there can be some differences in space number. so use '\s' to replace ' '
#' this function only works for the simple formulation of i01
#' for later version with much complex dependence formulation it is no longer working.
#'
#' @param path string. the location of i01 for prior information. must be provided
#' @return list. list containing prioir information
#' @export
#' @import stringr magrittr
prior_reader<-function(path=NULL){
  if(is.null(path)){
    stop("please provide input path")
  }
  lines=readLines(path)
  spec_tit_ind=str_which(string=lines,pattern="^\\s*namespec/jfix\\s+npulse\\s+nperiod\\s+jmsspec\\s+xspec_min\\s+xspec_max\\s*$")
  spec_dattit_ind=str_which(string=lines,pattern="^\\s*ipm\\s+ivpm\\s+pmspec\\s+jctspec\\s+jbcspec\\s*$")
  blanck_ind=str_which(string=lines,pattern="\\#|([:alpha:][:alpha:]+)")
  lines[spec_tit_ind+1] %>% str_trim(string=.,side="both") ->names
  conc_init_range=sapply(spec_dattit_ind,function(x){
    start=x+1
    end=blanck_ind[blanck_ind>start][1]-1
    inds=start:end
    range=c()
    if(length(inds)>=2){
      lines[start] %>% str_trim(string=.,side="both") %>%
                  str_split(string=.,pattern="\\s+") -> low
      range[1]=low[[1]][3]
      lines[start+1] %>% str_trim(string=.,side="both") %>%
                  str_split(string=.,pattern="\\s+") -> high
      range[2]=high[[1]][3]
    }
    if(is.null(range)){
      range=c(NA,NA)
    }
    range
  })
  conc_init_range=t(conc_init_range)
  rownames(conc_init_range)=names
  conc_init_range=conc_init_range[!is.na(conc_init_range[,1]),]
  reac_tit_ind=str_which(string=lines,pattern="^\\s*namereac\\s+/nipart\\s+nopart\\s+jkin\\s*$")
  reac_dattit_ind=str_which(string=lines,pattern="^\\s*ipm\\s+ivpm\\s+pmreac\\s+jctreac\\s+jbcreac\\s*$")
  blanck_ind=str_which(string=lines,pattern="\\#|([:alpha:][:alpha:]+)")
  lines[reac_tit_ind+1] %>% str_trim(string=.,side="both") ->names
  conc_kine_range=sapply(reac_dattit_ind,function(x){
    start=x+1
    end=blanck_ind[blanck_ind>start][1]-1
    inds=start:end
    range=c()
    if(length(inds)>=2){
      loci=str_which(string=lines[start:end],pattern="\\s+1\\s+1\\s*$")+start-1
      for(locus in loci){
        lines[locus] %>% str_trim(string=.,side="both") %>%
                    str_split(string=.,pattern="\\s+") -> low
        range=c(range,low[[1]][3])
        lines[locus+1] %>% str_trim(string=.,side="both") %>%
                    str_split(string=.,pattern="\\s+") -> high
        range=c(range,high[[1]][3])
      }
    }
    if(is.null(range)){
      range=c(NA,NA,NA,NA)
    }
    range
  })
  conc_kine_range=t(conc_kine_range)
  # k1=conc_kine_range[,1:2]
  # k2=conc_kine_range[,3:4]
  kcat=conc_kine_range[,1:2]
  km=conc_kine_range[,3:4]
  # rownames(k1)=names
  # rownames(k2)=names
  rownames(kcat)=names
  rownames(km)=names
  km=km[!is.na(km[,1]),]
  kcat=kcat[!is.na(kcat[,1]),]
  # list.res=list(enz=conc_init_range,k1=k1,k2=k2)
  list.res=list(enz=conc_init_range,kcat=kcat,km=km)
  return(list.res)
}
