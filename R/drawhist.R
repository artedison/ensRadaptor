#' plot distribution histogram for summary
#'
#' plot distribution histogram for summary
#'
#' @param tab array. the data table. must be provided
#' @param loci string. location for result(add on string for name). default ""
#' @param xlab string. x axis. default "xlab"
#' @return no return just plot the histogram
draw_hist<-function(tab=NULL,loci="",xlab="xlab"){
  if(is.null){
    stop("please provide input data")
  }
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
