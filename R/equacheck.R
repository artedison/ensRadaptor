#' check equalabration
#'
#' this function check the trend of chi^2 and parameters
#' this is prefered before accumulation run
#' it can check the chi^2~sweeps, initial_concentration~time, and kinetic parameters~time
#'
#' @param o02data list. o02 data struct
#' @param sweeps numeric. the last n sweeps that should be included in the drawing
#' @param addon string. addon string for names
#' @param comp string. whether plot parameter with sweep?
#'       "chi2": only plot chi2
#'       "all": plot all
#' @return no return just plot
#' @export
#' @examples
#'
equa_check<-function(o02data,sweeps,addon,comp){
  len=length(o02data[["ids"]])
  rang=(len-sweeps+1):len
  ##chi^2
  tab=cbind(o02data[["ids"]][rang],2*o02data[["chisq"]][rang])
  draw_sweep(tab,"chi^2",
      paste0(dir.res,"chi2-sweep.check",addon,".pdf"),FALSE)
  if(comp=="all"){
    ##ini concentration
    inis=names(o02data[["theta_spec_ini"]][[1]])
    tabchose=Reduce(rbind,o02data[["theta_spec_ini"]][rang])
    for(ini in inis){
      tab=cbind(o02data[["ids"]][rang],tabchose[,ini])
      draw_sweep(tab,ini,
          paste0(dir.res,"iniconec-sweep.check",ini,".",addon,".pdf"),FALSE)
    }
    ##kinetic parameter
    paras=names(o02data[["theta_reac_para"]][[1]])
    tabchose=Reduce(rbind,o02data[["theta_reac_para"]][rang])
    for(para in paras){
      tab=cbind(o02data[["ids"]][rang],tabchose[,para])
      draw_sweep(tab,para,
          paste0(dir.res,"kinepara-sweep.check",para,".",addon,".pdf"),FALSE)
    }
  }
}
