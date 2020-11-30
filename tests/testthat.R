library(testthat)
library(ensRadaptor)

# unit test
test_check("ensRadaptor")

# integration test
currdir=getwd()
cat(currdir)
system("rm ./temp/input_test/history/history.record.tab",ignore.stderr=TRUE)
system("rm -r ./temp/input_test/testmodel/*",ignore.stderr=TRUE)
source(paste0(currdir,"/workflowtest_open/script.input.standalone.R"))
dir_ext_data=system.file("extdata","", package = "ensRadaptor")
if(!str_detect(string=dir_ext_data,pattern="\\/$")){
  dir_ext_data=paste0(dir_ext_data,"/")
}
predatadir=paste0(dir_ext_data,"internal_data/pre_result_1/")
filelist=c("ens.sh","ens.i01","ens.i02","ens.def","ens.f90")#
cat("Test: differences between the produced files and stored templates\n")
for(file in filelist){
  output=system(paste0("diff \'",dir.input,"1/",file,"\' \'",predatadir,file,"\'"),intern=TRUE)
  if(length(output)==0){
    cat(paste0(file,"\tPASSED\n"))
  }else{
    cat(paste0(file,"\n"))
    cat(output)
  }
}
listpre1=mget(load(paste0(predatadir,"local.prior.pre.170008_01102020.RData")))
listnew=mget(load(savedrdata))
comapared_names=names(listpre1)
# the data that will be definitely different (false positive) as they depends on the locations.
comapared_names=comapared_names[!comapared_names%in%c("modified.file","foldname","dir.res")]
flagmatch=sapply(comapared_names,function(namedata){
  identical(listpre1[[namedata]],listnew[[namedata]])
})
cat(paste0("Test: difference in stored data \n"))## there can be differences in folder
cat(length(flagmatch[!flagmatch]),"\n")
cat(paste0(names(flagmatch[!flagmatch]),"\n"))
# system("rm ./temp/input_test/history/history.record.tab")
# system("rm -r ./temp/input_test/testmodel/*")

# # the following block will take a few minutes to run and produce a lot of files. run this part only when changes has been done on script.output.standalone.R
# # this test is dependent on some big files which are not shared thorough github
# rm(list=ls())
# # # del folder
# currdir=getwd()
# system("rm -r ./temp/output_test/*",ignore.stderr=TRUE)
# source(paste0(currdir,"/workflowtest_open/script.output.standalone.R"))
# #compare expected one difference in each comparison as it is from time
# #skip first 100 bytes in pdf file
# predatadir=paste0(dir_ext_data,"internal_data/pre_result_testplotting/")
# pdflist=list.files(predatadir,pattern="*.pdf")
# cat("Test: compare each produced figures")
# ## test of identical in figure files are hard
# # for(file in pdflist){
# #   if(!file.exists(paste0(dir.res,file))){
# #     cat(paste0(file,"\tNOT EXIST\n"))
# #   }
# #   output=system(paste0("cmp -i 100 \'",dir.res,file,"\' \'",predatadir,file,"\'"),intern=TRUE)
# #   if(length(output)==0){
# #     cat(paste0(file,"\tPASSED\n"))
# #   }else{
# #     cat(paste0(file,"\n"))
# #     cat(output)
# #   }
# # }
# # system("rm -r ./temp/output_test/*")
