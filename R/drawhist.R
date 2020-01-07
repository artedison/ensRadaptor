#' plot distribution histogram for summary
#'
#' plot distribution histogram for summary
#'
#' @param tab array. the data table
#' @param loci string. location for result(add on string for name)
#' @param xlab string. x axis
#' @return no return just plot the histogram
draw_hist<-function(tab,loci,xlab){
  tab=as.data.frame(tab)
  colnames(tab)=c("val")
  tab[,1]=as.numeric(tab[,1])
  p<-ggplot(data=tab,aes(val))+
        geom_histogram(position="stack")+
        xlab(xlab)+
        ylab("frequency")+
        theme_bw()
  ggsave(plot=p,file=loci)
}
