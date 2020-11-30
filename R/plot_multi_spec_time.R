#' plot of time trajectory of multiple species
#'
#' plot of time trajectory of multiple species
#'
#' @param o01 string. the location of o01. must be provided
#' @param species.list array. the species to plot. must be provided
#' @param addon string. addon file name. default ""
#' @param exp int. the experiment number to choose. must be provided
#' @param dir.res string. result locaiton. must be provided
#' @param ncol int. the number of columne in plot_grid plot. default 5
#' @param logtrans bool. whether log transformation is  performed. default FALSE
#' @return figure. the figure itself
#' @export
#' @import stringr magrittr ggplot2 cowplot
plot_multi_spec_time<-function(o01=NULL,species.list=NULL,addon="",exp=NULL,dir.res=NULL,ncol=5,logtrans=FALSE){
  if(is.null(o01)){
    stop("please provide input path")
  }
  if(is.null(species.list)){
    stop("please provide species list")
  }
  if(is.null(exp)){
    stop("please provide the experiment number")
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
  list.p=vector(mode="list")
  for(spec_meta in species.list){
    tabmod=mod_model_tab(lines,start=start_mod_ind,end=start_mod_ind,specname_ind=specname_ind,power_ind=power_ind,species=spec_meta)
    modpart=tabmod[tabmod[,"mod"]==exp,c(1,2,3)]
    tab=modpart
    colnames(tab)=c("time","val","var")
    tab=as.data.frame(tab)
    ##a not perfect solution for start jump
    # if(spec_meta %in% c("pyruvate","Smalate")){
    tab=tab[tab[,"time"]>min(tab[,"time"]),]
    # }
    ##
    p<-ggplot(data=tab,aes(time,val))+
          geom_ribbon(aes(ymin=val-2*var,ymax=val+2*var),fill="grey70",colour="grey70")+
          geom_line(size=1)+
          xlab("time(h)")+
          ylab("")+
          expand_limits(x=0.0,y=0.0)+
          theme_bw()
    if(logtrans){
      p=p+scale_y_log10(breaks=round(10^seq(log10(min(tab[,"val"])),log10(max(tab[,"val"])),length.out=5),4))
    }
    list.p[[spec_meta]]=p
  }
  p<-plot_grid(plotlist=list.p,labels=c(species.list),
              label_size=10,ncol=ncol,label_x=0.3)
  ggsave(plot=p,filename=paste0(dir.res,addon,"_",exp,"multi.concentr.pdf"),
    limitsize=FALSE)
  return(p)
}
