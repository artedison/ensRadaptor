#' boxplot of multiple species and different experiments
#'
#' can be parameters or intial value
#' prior range of species can be also presented
#'
#' @param list.spe list. data.list
#' @param path string. result location
#' @param xlab string. x axis
#' @param ylab string. y axis
#' @param list.prior.part list. the prior range
#' @param logtrans bool. whether use log transform
#' @param rank array. rank of presenting the data
#' @return just plot no return
#' @export
#' @seealso [boxplot_multi_exp_comb()]
boxplot_multi_exp<-function(list.spe,path,xlab,ylab,list.prior.part=NULL,logtrans=FALSE,rank){
    list.tab=sapply(names(list.spe),simplify=FALSE,function(name){
      temptab=list.spe[[name]]
      temp=cbind(unlist(temptab),rep(names(temptab[[1]]),times=length(temptab)))
      cbind(temp,rep(name,times=dim(temp)[1]))
    })
    tab=Reduce(rbind,list.tab)
    temp_pattern=str_replace(string=tab[,2],pattern="\\_[\\dr]+$",replacement="")
    tab_list=sapply(rank,simplify=FALSE,function(x){
      tab[temp_pattern==x,]
    })
    tab=Reduce("rbind",tab_list)
    rownames(tab)=NULL
    tab=as.data.frame(tab)
    colnames(tab)=c("val","type","exp")
    tab[,3]=as.factor(tab[,3])
    tab[,3]=ordered(tab[,3],levels=unique(tab[,3]))
    ecid=str_replace_all(string=tab[,2],pattern="\\_\\d+$",replacement="")
    tab[,2]=factor(ecid,levels=unique(ecid))
    tab[,1]=as.numeric(tab[,1])
    # if(logtrans){
    #   tab[,1]=log10(tab[,1])
    # }
    p<-ggplot(data=tab,aes(x=type,y=val,color=exp))+
          geom_boxplot(outlier.alpha=0.05)+
          xlab(xlab)+
          ylab(ylab)+
          theme_bw()+
          theme(axis.text.x=element_text(angle=90,hjust=1))
    if(logtrans){
      breaks=10^seq(log10(min(tab[,"val"])),log10(max(tab[,"val"])),length.out=8)
      p=p+scale_y_log10(breaks=breaks,labels=scientific)
      # p=p+scale_y_log10(breaks=round(10^seq(log10(min(tab[,"val"])),log10(max(tab[,"val"])),length.out=5),6),labels=scientific)
    }
  if(!is.null(list.prior.part)){
    list.tab.add=sapply(names(list.spe),simplify=FALSE,function(name){
      temptab=list.prior.part[[name]]
      temp=cbind(temptab,rownames(temptab))
      cbind(temp,rep(name,times=dim(temp)[1]))
    })
    tab.add=Reduce(rbind,list.tab.add)
    temp_pattern=str_replace(string=tab[,2],pattern="\\_[\\dr]+$",replacement="")
    tab_list=sapply(rank,simplify=FALSE,function(x){
      tab[temp_pattern==x,]
    })
    tab=Reduce("rbind",tab_list)
    rownames(tab.add)=NULL
    tab.add=as.data.frame(tab.add)
    colnames(tab.add)=c("low","high","type","exp")
    tab.add[,4]=as.factor(tab.add[,4])
    tab.add[,3]=as.factor(str_replace_all(string=tab.add[,3],pattern="\\_\\d+$",replacement=""))
    # tab.add[,1]=log10(as.numeric(tab.add[,1]))
    # tab.add[,2]=log10(as.numeric(tab.add[,2]))
    tab.add[,1]=as.numeric(tab.add[,1])
    tab.add[,2]=as.numeric(tab.add[,2])
    p2<-ggplot(data=tab.add)+
          geom_point(aes(x=type,y=low,color=exp))+
          geom_point(aes(x=type,y=high,color=exp))+
          xlab(xlab)+
          ylab(ylab)+
          theme_bw()+
          theme(axis.text.x=element_text(angle=90,hjust=1))
    if(logtrans){
      # p2=p2+scale_y_log10(breaks=round(10^seq(log10(min(tab.add[,"low"])),log10(max(tab.add[,"high"])),length.out=5),4),labels=scientific)
      breaks=10^seq(log10(min(tab[,"val"])),log10(max(tab[,"val"])),length.out=8)
      p2=p2+scale_y_log10(breaks=breaks,labels=scientific)
    }
    p=plot_grid(p,p2,labels=c("result", "prior"),nrow=2)
  }
  # p
  ggsave(plot=p,file=path,width=14,height=7)
}

#' boxplot of multiple species and different experiments
#'
#' can be parameters or intial value
#' prior range of species can be also presented
#' the differences from boxplot.multi.exp is that this funciton will plot prior and result side by side
#'
#' @param list.spe list. list of data
#' @param path string. result location
#' @param xlab string. x axis
#' @param ylab string. y axis
#' @param list.prior.part list. the prior range
#' @param logtrans bool. whether use log transform
#' @param rank array. rank of presenting the data
#' @return just plot no return
#' @seealso [boxplot_multi_exp()]
boxplot_multi_exp_comb<-function(list.spe,path,xlab,ylab,list.prior.part=NULL,logtrans=FALSE,rank){
    list.tab=sapply(names(list.spe),simplify=FALSE,function(name){
      temptab=list.spe[[name]]
      temp=cbind(unlist(temptab),rep(names(temptab[[1]]),times=length(temptab)))
      cbind(temp,rep(name,times=dim(temp)[1]))
    })
    tab=Reduce(rbind,list.tab)
    temp_pattern=str_replace(string=tab[,2],pattern="\\_[\\dr]+$",replacement="")
    tab_list=sapply(rank,simplify=FALSE,function(x){
      tab[temp_pattern==x,]
    })
    tab=Reduce("rbind",tab_list)
    rownames(tab)=NULL
    tab=as.data.frame(tab)
    colnames(tab)=c("val","type","exp")
    tab[,3]=as.factor(tab[,3])
    tab[,3]=ordered(tab[,3],levels=unique(tab[,3]))
    ecid=str_replace_all(string=tab[,2],pattern="\\_\\d+$",replacement="")
    tab[,2]=factor(ecid,levels=unique(ecid))
    tab[,1]=as.numeric(tab[,1])
    ###prior
    list.tab.add=sapply(names(list.spe),simplify=FALSE,function(name){
      temptab=list.prior.part[[name]]
      namessele=names(list.spe[[1]][[1]])
      namessele=str_replace(string=namessele,pattern="\\_\\d+$",replacement="")
      temptab=temptab[namessele,]
      temp=cbind(temptab,rownames(temptab))
      cbind(temp,rep(name,times=dim(temp)[1]))
    })
    tab.add=Reduce(rbind,list.tab.add)
    rownames(tab.add)=NULL
    tab.add=as.data.frame(tab.add)
    colnames(tab.add)=c("low","high","type","exp")
    tab.add[,4]=as.factor(tab.add[,4])
    tab.add[,3]=as.factor(str_replace_all(string=tab.add[,3],pattern="\\_\\d+$",replacement=""))
    tab.add[,1]=as.numeric(tab.add[,1])
    tab.add[,2]=as.numeric(tab.add[,2])
    p<-ggplot(data=tab,aes(x=type,y=val),color="blue")+
          geom_boxplot(outlier.alpha=0.05,outlier.color="grey70")+
          geom_point(data=tab.add,aes(x=type,y=low),color="red")+
          geom_point(data=tab.add,aes(x=type,y=high),color="red")+
          xlab(xlab)+
          ylab(ylab)+
          theme_bw()+
          theme(axis.text.x=element_text(angle=90,hjust=1))
    if(logtrans){
      breaks=10^seq(log10(min(tab.add[,"low"])),log10(max(tab.add[,"high"])),length.out=8)
      p=p+scale_y_log10(breaks=breaks,labels=scientific)
    }
  # p
  ggsave(plot=p,file=path,width=14,height=7)
}
