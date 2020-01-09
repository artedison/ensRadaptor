#' produce i01 file
#'
#' produce the template for species and reactions
#' read the reaction file formated as template/addon.reac.tab
#'
#' @param path string. the reaction file path
#' @param type string. the reaction type that need to be formated
#'        "mr": mass reaction
#'        "mm": Michaelis–Menten kinetics
#' @param dir.data string. the data directory(working folder)
#' @param modified.file string. the addon name
#' @param para.list list. the input parameter list that need to be changed
#'        a list contains parameters related to species and reactions
#'
#'        Attention: if there is no regulation to formualte then regu=NULL
#'        species: obsv(jmsspec), fixed(jctspec), initial(4th, 5th of ivpm) const(jfix)
#'        react: kine(kinetic parameter list), enz(the enzyme in mm reaction, need for ranking), rev(reversible, 1 indicate reversible 0 non-reversible),
#'
#' @param list.exi list. existing species and reactions list
#' @param extend numeric the enlarge amount/by proportion of parameters(k1 k2). default 1
#' @param rand bool. related to generating initial theta. FALSE not randomeness, geometry mean of the boundary. TRUE unif random generating.
#' @param rand.seed numeric. the seed for sampling initial condiiton
#' @return list. parameters and related sizes
#' @export
template_spec_reac<-function(path,type,dir.data,modified.file,para.list,list.exi,extend=1,rand=FALSE,rand.seed=0){
  ## reading source files
  # smallvalue=0.0000000000000000001
  sourcefile=paste0(dir.data,"ens.i01")
  outfile=paste0(dir.data,"ens.modified.i01")
  outfile12=paste0(dir.data,"ens.modified.add.i12")
  modilines=readLines(sourcefile)
  # prioroutput=paste0(dir.data,"prior.",modified.file,".RData")
  list.reac.addon=read_reac(path)
  names(list.reac.addon)=sapply(list.reac.addon,function(x){
    x[["name"]]
  })
  vec.spec.addon=unique(unlist(sapply(list.reac.addon,function(x){
    c(x[["subs"]],x[["prods"]])
  })))
  vec.spec.addon=setdiff(vec.spec.addon,list.exi[["species"]])
  save(list.reac.addon,vec.spec.addon,file=paste(dir.data,modified.file,"addon.reac.speci.RData"))
  templist=vector(mode="list")
  local.para=list(enz=templist,kcat=templist,km=templist,kcatr=templist,kmr=templist)
  ## load default parameters
  if((!para.list[["species"]][["format"]])||(!para.list[["react"]][["format"]])){
    load(paste0(dir.lib,"default.para.mr.RData"))
    # load(paste0(dir.lib,"default.para.mm.RData"))
  }
  if(type=="mr"){
    jkin=1
  }else if(type=="mm"){
    jkin=10011
    # jmsspec="0011"
    jmsspec="0001"
  }
  ###regulation formulation
  reguinfor=vector(mode="list")
  if(!is.null(para.list[["regu"]][["listfile"]])){
    regutab=read.table(para.list[["regu"]][["listfile"]]);
    regutab=regutab[,c(2,3,4)]
    colnames(regutab)=c("compd","regu","enz")
    regutab[,"regu"]=ifelse(regutab[,"regu"]=="->",1,-1)
    reguinfor[["regutab"]]=regutab
    wrongname=setdiff(regutab[,"compd"],vec.spec.addon)
    if(length(wrongname)!=0){
      warning("wrong names for compound in the regulation table")
      warning(paste(wrongname,collapse=","))
    }
    #### load the fomat for regulation to use later
    lines=readLines(para.list[["regu"]][["format"]])
    boundind=str_which(string=lines,pattern="^####")
    reguformat=list(compd_regu=lines[(boundind[2]+1):(boundind[3]-1)],
                  regu_enz=lines[(boundind[3]+1):(boundind[4]-1)])
    reguinfor[["reguformat"]]=reguformat
  }
  ###output the species format
  if(para.list[["species"]][["format"]]){
    # environment(spec_output_format)<-environment()
    spec_output=spec_output_format()
  }else{
    # environment(spec_output)<-environment()
    spec_output=spec_output()
  }
  cat(spec_output,file=paste(dir.data,"species.addon.txt"),sep="\n")
  add_ind=str_which(string=modilines,pattern="\\# add on by yue new species")
  seqfront=seq(add_ind)
  seqafter=seq(from=add_ind+1,by=1,to=length(modilines))
  modilines=c(modilines[seqfront],
              spec_output,
              modilines[seqafter])
  ###output the reaction format
  reac.para.addon.len=0
  reac.count=0
  # if(type=="mm"){
  #   k1rag=k1.range(list.reac.addon,para.list[["react"]][["kine"]])
  # }
  # list.kine.prior=list(km=vector(mode="list"),kcat=vector(mode="list"))
  if(para.list[["react"]][["format"]]){
    # environment(react.output.format)<-environment()
    list.res=react_output_format()
  }else{
    # environment(react.output)<-environment()
    list.res=react_output()
  }
  reac_output=list.res[["out"]]
  reac.para.addon.len=list.res[["len"]]
  reac.count=list.res[["count"]]
  cat(reac_output,file=paste(dir.data,"reactions.addon.txt"),sep="\n")
  add_ind=str_which(string=modilines,pattern="\\# add on by yue new reactions")
  seqfront=seq(add_ind)
  seqafter=seq(from=add_ind+1,by=1,to=length(modilines))
  modilines=c(modilines[seqfront],
              reac_output,
              modilines[seqafter])
  cat(modilines,sep="\n",file=outfile)
  # list.kine.prior.tab[["km"]]=Reduce("rbind",list.kine.prior[["km"]])
  # list.kine.prior.tab[["kcat"]]=Reduce("rbind",list.kine.prior[["kcat"]])
  # save(list.kine.prior,file=prioroutput)
  ## output for the .i12 file
  output=paste("# species addon",
              paste(rep(0.02,times=length(vec.spec.addon)),collapse=" "),
              "# reactions addon",
              paste(rep(0.4,times=reac.para.addon.len),collapse=" "),
              sep="\n")
  cat(output,file=outfile12,sep="\n")
  lenlist=list(spec=length(vec.spec.addon),reac=reac.count,localpara=local.para)
  return(lenlist)
}
