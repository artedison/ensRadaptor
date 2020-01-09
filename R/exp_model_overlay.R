#' overlay experiment measurement with model predicted result in o01
#'
#' measurement is actually "scaled" measurement
#' the range [low high] will also be plot
#'
#' @param o01 string. file path. must be provided
#' @param species string. the species that need overlay. must be provided
#' @param addon string. addon on output file name. default ""
#' @param exp numeric. the experiment number that need to be plot. must be provided
#' @param dir.res string. the output folder. must be provided
#' @param logtrans bool. whether log transformaed is performed on dat. default FALSE
#' @return figure. return the figure itself
#' @export
exp_model_overlay<-function(o01=NULL,species=NULL,addon="",exp=NULL,dir.res=NULL,logtrans=FALSE){
  if(is.null(o01)){
    stop("please provide input path")
  }
  if(is.null(species)){
    stop("please provide species list")
  }
  if(is.null(exp)){
    stop("please provide experiment list")
  }
  if(is.null(dir.res)){
    stop("please provide output path")
  }
  ## reading o01
  lines=readLines(o01)
  start_mod_ind=str_which(string=lines,pattern="^\\s*MC\\s+averaged\\s+model\\s+species\\s+concentrations")
  start_exp_ind=str_which(string=lines,pattern="^\\s*exptl\\.\\s+observed\\s+species\\s+at\\s+exptl\\.\\s+time\\s+points")
  specname_ind=str_which(string=lines,pattern="namespec")
  power_ind=str_which(string=lines,pattern="^\\s*ipow")
  ch_ind=str_which(string=lines,pattern="[:alpha:][:alpha:]+")
  ### experiment result
  tabexp=exp_model_tab(lines,start_exp_ind,start_mod_ind,specname_ind,power_ind,ch_ind,species)
  ### modeled result
  tabmod=mod_model_tab(lines,start_mod_ind,start_mod_ind,specname_ind,power_ind,species)
  ##plot
  exppart=tabexp[tabexp[,"exp"]==exp,c(1,2)]
  modpart=tabmod[tabmod[,"mod"]==exp,c(1,2,3)]
  timerange=c(min(exppart[,"time"]),max(exppart[,"time"]))
  modpart=modpart[modpart[,"time"]>=timerange[1]&modpart[,"time"]<=timerange[2],]
  exppart=cbind(cbind(exppart,exppart[,"val"]),exppart[,"val"])
  colnames(exppart)=c("time","val","low","high")
  modpart2=cbind(cbind(modpart[,c("time","val")],modpart[,"val"]-2*modpart[,"var"]),
                modpart[,"val"]+2*modpart[,"var"])
  colnames(modpart2)=c("time","val","low","high")
  # if(logtrans){
  #   for(col in c("val","low","high")){
  #     exppart[,col]=log10(exppart[,col])
  #     modpart2[,col]=log10(modpart2[,col])
  #   }
  # }
  tab=cbind(rbind(exppart,modpart2),rep(c("exp","mod"),times=c(dim(exppart)[1],dim(modpart2)[1])))
  colnames(tab)=c("time","val","low","high","type")
  tab=as.data.frame(tab)
  tab[,"type"]=as.factor(tab[,"type"])
  p<-ggplot(data=tab,aes(time,val,colour=type,group=type))+
        geom_ribbon(data=tab[tab[,"type"]=="mod",],aes(ymin=low,ymax=high),fill="grey70",colour="grey70")+
        geom_line(data=tab[tab[,"type"]=="mod",],size=1,alpha=0.85)+
        geom_point(data=tab[tab[,"type"]=="exp",],size=2)+
        xlab("time")+
        ylab(species)+
        expand_limits(x=0.0,y=0.0)+
        theme_bw()
  if(logtrans){
    p=p+scale_y_log10(breaks=round(10^seq(log10(min(tab[,"low"])),log10(max(tab[,"high"])),length.out=5),4))
  }
  ggsave(plot=p,filename=paste0(dir.res,species,"_",addon,"_",exp,"concentr.pdf"))
  return(p)
}
