#' experiment specie concentration fetching(o01)
#'
#' experiment specie concentration fetching from the table in file o01
#'
#' @param lines array. lines to parse
#' @param start array. start of the index range
#' @param end array. end of the index range
#' @param specname_ind array. the index of species part
#' @param power_ind array. the index of power line
#' @param ch_ind array. the index of the number line
#' @param species string. the specie to fetch
#' @return array. the informaiton table
#' @seealso [mod_model_tab()] for similar function on model estimation information
exp_model_tab<-function(lines,start,end,specname_ind,power_ind,ch_ind,species){
  specname_ind_exp=specname_ind[specname_ind>start&specname_ind<end]
  titltable_ind_exp=str_which(string=lines,pattern="^\\s+itxpt")
  ch_ind_exp=ch_ind[ch_ind>start&ch_ind<=end]
  exp_str=sapply(seq(length(specname_ind_exp)),simplify=FALSE,function(x){
          specname_block=specname_ind_exp[x]
          (specname_block+1) %>% extract(lines,.) %>%
                    str_trim(string=.,side="both") %>%
                    str_split(string=.,pattern="\\s+") %>%
                    unlist() -> specs
          if(specs==species){
            tabstart=titltable_ind_exp[titltable_ind_exp>specname_block][1]+1
            tabend=ch_ind_exp[ch_ind_exp>tabstart][1]-1
            tab_block=tabstart:tabend
            tab_block=tab_block[!str_detect(string=tab_block,pattern="^\\s+$")]
            temp=lines[tab_block]
            temp=temp[!str_detect(string=temp,pattern="^\\s*$")]
            if(length(temp)>0){
              tab=tabparse(lines,tab_block,c(2,4))
            }else{
              NULL
            }
         }else{
           NULL
         }
  })
  exp_str=exp_str[sapply(exp_str,function(x){!is.null(x)})]
  len=sapply(exp_str,function(x){dim(x)[1]})
  tabexp=cbind(Reduce(rbind,exp_str),rep(seq(length(len)),times=unlist(len)))
  colnames(tabexp)=c("time","val","exp")
  return(tabexp)
}
#' model specie concentration fetching(o01)
#'
#' model specie concentration fetching from the table in file o01
#'
#' @param lines array. lines to parse
#' @param start array. start of the index range
#' @param end array. end of the index range
#' @param specname_ind array. the index of species part
#' @param ch_ind array. the index of the number line
#' @param species string. the specie to fetch
#' @return array. the informaiton table
#' @seealso [exp_model_tab()] for similar function on model estimation information
mod_model_tab<-function(lines,start,end,specname_ind,power_ind,species){
  specname_ind_mod=specname_ind[specname_ind>start]
  power_ind_mod=power_ind[power_ind>start]
  titltable_ind_mod=str_which(string=lines,pattern="^\\s+io_time")
  mod_str=sapply(seq(length(specname_ind_mod)),simplify=FALSE,function(x){
          specname_block=specname_ind_mod[x]
          power_block=power_ind_mod[power_ind_mod>specname_block][1]
          (specname_block+1) %>% extract(lines,.) %>%
                    str_trim(string=.,side="both") %>%
                    str_split(string=.,pattern="\\s+") %>%
                    unlist() -> specs
          (power_block+1) %>% extract(lines,.) %>%
                    str_trim(string=.,side="both") %>%
                    str_split(string=.,pattern="\\s+") %>%
                    unlist() -> power
          if(specs==species && power=="1"){
            tabstart=titltable_ind_mod[titltable_ind_mod>specname_block][1]+1
            tabend=power_ind_mod[power_ind_mod>tabstart][1]-1
            tab_block=tabstart:tabend
            tab=tabparse(lines,tab_block,c(4,7,8))
         }else{
           NULL
         }
  })
  mod_str=mod_str[sapply(mod_str,function(x){!is.null(x)})]
  len=sapply(mod_str,function(x){dim(x)[1]})
  tabmod=cbind(Reduce(rbind,mod_str),rep(seq(length(len)),times=unlist(len)))
  colnames(tabmod)=c("time","val","var","mod")
  return(tabmod)
}