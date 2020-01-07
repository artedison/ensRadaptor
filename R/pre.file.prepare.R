#' prepare input files
#'
#' thie function will copy all template file needed into the temptory directory
#'
#' @param dir.tempt string. the temptory directory
#' @param templatepath string. the template path
#' @return just move the file no return value
#' @export
pre_file_prepare<-function(dir.tempt,templatepath){
  system(paste0("cp \"",templatepath[["input.template"]],"ens.i01\" \"",dir.tempt,"\""))
  system(paste0("cp \"",templatepath[["input.template"]],"ens.i02\" \"",dir.tempt,"\""))
  system(paste0("cp \"",templatepath[["enscode"]],"ens.f90\" \"",dir.tempt,"\""))
  system(paste0("cp \"",templatepath[["enscode"]],"ens.def\" \"",dir.tempt,"\""))
  system(paste0("cp \"",templatepath[["input.template"]],"ens.sh\" \"",dir.tempt,"\""))
  system(paste0("cp \"",templatepath[["input.template"]],"submit.sh\" \"",dir.tempt,"\""))
}
