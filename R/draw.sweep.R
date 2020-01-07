#' plot var~sweep
#'
#' plot the variable through sweeps
#'
#' @param tab array. the data table
#' @param para string. y axis
#' @param loci string. location for result (add on string for name)
#' @param linethick bool. whether make line thicker
#' @return no return just plot
draw_sweep<-function(tab,para,loci,linethick){
  if(linethick){
    size=1
  }else{
    size=0.5
  }
  tab=as.data.frame(tab)
  colnames(tab)=c("sweep","val")
  tab[,1]=as.numeric(tab[,1])
  tab[,2]=as.numeric(tab[,2])
  p<-ggplot(data=tab,aes(x=sweep,y=val))+
        geom_line(size=size,alpha=0.75)+
        # geom_point(color="red",alpha=0.5,size=0.5)+
        xlab("sweep")+
        ylab(para)+
        theme_bw()
  ggsave(plot=p,file=loci)
}
