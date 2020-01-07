#' initial guess generate
#'
#' this function will generate the initial guess for all parameters by lhs (Latin Hypercube sample) in log-uniform distribution
#' the lhs matrix will be calcualted for isample and used afterwards
#'
#' @param nsample int. number of replicate to produce
#' @param isample int. the sample index
#' @param input string. input file path
#' @param output string. output file path
#' @return just change the file no return
#' @export
randloglhs<-function(nsample,isample,input,output){
  lines=readLines(input)
  indspec=str_which(string=lines,pattern="\\#\\s+Species\\s+control\\-\\s+and\\s+\\\\Theta\\-variables")
  lines_abov=lines[1:indspec]
  lines_below=lines[(indspec+1):length(lines)]
  indtheta=str_which(string=lines_below,pattern="\\s+1\\s+1\\s*$")
  ntheta=length(indtheta)
  if(isample==1){
    set.seed(1)## this seed is just set here.
    lhsmatrix<<-randomLHS(nsample,ntheta)
  }
  lines_below[indtheta] %>% str_trim(string=.,side="both") %>%
                    str_split(string=.,pattern="\\s+",simplify=TRUE) %>%
                    extract(,3) %>% as.numeric(.) -> lowval
  lines_below[indtheta+1] %>% str_trim(string=.,side="both") %>%
                    str_split(string=.,pattern="\\s+",simplify=TRUE) %>%
                    extract(,3) %>% as.numeric(.) -> highval
  iniguess=lowval*(highval/lowval)^lhsmatrix[isample,]
  lines_below[indtheta+3]=str_replace(string=lines_below[indtheta+3],pattern="[\\d\\.]+\\s*$",replacement=as.character(iniguess))
  cat(c(lines_abov,lines_below),sep="\n",file=output)
}
