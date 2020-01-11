## test workflow for ensRadaptor
## this testing script analyze output files from modeling central metabolism in Neurospora crassa based on experimental data.
## the output file should be fed to ens modeling for ODE parameter estimation
## works with sapelo2, template files need to be modified for adaption to other sever systems
## unit time: hour h
## unit concentration: mM
## kinetic form: MM (assumed)
## optimized parameter: kcat, km, E
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(R.matlab)
require(ggplot2)
require(xml2)
require(cowplot)
require(scales)
require(lhs)
require(foreach)
require(doMC)
require(ensRadaptor)
registerDoMC(cores=10)

#------------------------------path setup------------------------------#
compdir="/Users/yuewu/"
dirpack=paste0(compdir,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/package.formulate/ensRadaptor/")
dir=paste0(dirpack,"temp/testworkflow/")
dir_ext_data=paste0(dirpack,"inst/extdata/")## the open external data directory
dir.lib=paste0(dir_ext_data,"template_format/")##this path is used within many functions
dir.data=paste0(dirpack,"internal_data/pre_input_testplotting/")
# dir.template.format=paste0(compdir,"Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/ensemble_infor/doc/template.file/")
foldname="1"
dir.res=paste0(dir,"testmodel/",foldname,"/res/")

load(paste0(dir.data,"local.prior.pre.122657_06242019.RData")) ## this need setup manually
foldname="1"
dir.res=paste0(dir,"testmodel/",foldname,"/res/")
##fetch of information from output file
name=modified.file
# replicateseq=c("1","2","3","4","5","8","10","11","12","13")
# not run as data not provided for file size
# info<-foreach(replicatei=replicateseq)%dopar%{
#   path.equala=paste0(dir.res,"equala/",replicatei,"/ens.o02")
#   o02.data.equ=o02_reader(path.equala)
#   summary_o02(o02.data=o02.data.equ,dir.res=dir.res,addonname=paste0("equala.",name,".",replicatei),linethick=TRUE)
#   equa_check(o02.data.equ,sweeps=1000,dir.res=dir.res,addon=paste0("equala.",name,".",replicatei),comp="chi2")
# }
### accumulation combine
# comb.accum("/Users/mikeaalv/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/nc_model/model/nc.glucose.fermetation/inputens/addon_model_measurement_new_cor_speed3_inivar_rem_aa_anaerobic/accumu.comb/ens.o02")

#------------------------------o02 summary------------------------------#
replicatei=1
path.accumu=paste0(dir.data,"/ens.o02")
o02.data.accumu=o02_reader(path.accumu)
summary_o02(o02.data=o02.data.accumu,dir.res=dir.res,addonname=paste0("accumu.",name,".",replicatei),linethick=TRUE)
equa_check(o02.data.accumu,sweeps=1000,dir.res=dir.res,
          addon=paste0("accumu.",name,".",replicatei),comp="chi2")

##new method return while output ens.i01
list.prior.pre=vector(mode="list")
list.prior.pre[["enz"]][["comb"]]=as_data_frame_rowname(len.list[["localpara"]][["enz"]])
list.prior.pre[["kcat"]][["comb"]]=as_data_frame_rowname(len.list[["localpara"]][["kcat"]])
list.prior.pre[["km"]][["comb"]]=as_data_frame_rowname(len.list[["localpara"]][["km"]])
list.prior.pre[["kcatr"]][["comb"]]=as_data_frame_rowname(len.list[["localpara"]][["kcatr"]])
list.prior.pre[["kmr"]][["comb"]]=as_data_frame_rowname(len.list[["localpara"]][["kmr"]])
rownames(list.prior.pre[["kcat"]][["comb"]])=paste0(rownames(list.prior.pre[["kcat"]][["comb"]]),"_S")
rownames(list.prior.pre[["kcatr"]][["comb"]])=paste0(rownames(list.prior.pre[["kcatr"]][["comb"]]),"_P")
list.prior.pre[["kcat"]][["comb"]]=rbind(list.prior.pre[["kcat"]][["comb"]],list.prior.pre[["kcatr"]][["comb"]])
rownames(list.prior.pre[["km"]][["comb"]])=paste0(rownames(list.prior.pre[["km"]][["comb"]]),"_S")
rownames(list.prior.pre[["kmr"]][["comb"]])=paste0(rownames(list.prior.pre[["kmr"]][["comb"]]),"_P")
list.prior.pre[["km"]][["comb"]]=rbind(list.prior.pre[["km"]][["comb"]],list.prior.pre[["kmr"]][["comb"]])
# list.prior.pre[["kmr"]][["comb"]]["comb_glycolysis_P",]=c(0.0018,7.5000)

#------------------------------o01 summary------------------------------#

##overlap with experiment data
paths.accu.01=paste0(dir.data,"/ens.o01")
for(exp in seq(length(experiments))){
  for(comp in obsv){
    exp_model_overlay(paths.accu.01,species=comp,
                    addon=paste0("_",name,"_repli_",replicatei),
                    exp=exp,dir.res=dir.res)
  }
}

##metabolite concentration
# allspe=list.exi[["species"]]
# meta.list=allspe[!str_detect(string=allspe,pattern="^NCU")]
names(meta.list)=meta.list
path.accumu.o01=paste0(dir.data,"/ens.o01")
meta.list.use=meta.list
for(exp in seq(length(experiments))){
  plot_multi_spec_time(path.accumu.o01,
            species.list=meta.list.use,
            addon=paste0(name,".multisub",".",replicatei),
            exp=exp,dir.res=dir.res)
}

# distribution of estimated enzyme
# compare between different extend range in parameters
reac.rank=unique(str_replace(string=rownames(list.prior.pre[["kcat"]][["comb"]]),pattern="\\_[SP]+$",replacement=""))
enz.all.list=vector(mode="list")
inilist=o02.data.accumu[["theta_spec_ini"]]
namelist=names(inilist[[1]])
for(exp in seq(length(experiments))){
  enz.list=namelist[str_detect(string=namelist,pattern=paste0("((^NCU)|(comp)).*\\_",exp))]
  data.list=sapply(inilist,simplify=FALSE,function(x){
              x[enz.list]
  })
  enz.all.list[[paste(exp)]]=data.list
}
enzrealval.list=enz.all.list
for(y in seq(length(enz.all.list))){
  names=names(enz.all.list[[1]][[y]])
  names=str_replace(string=names,pattern="\\_\\d+$",replacement="")
  for(x in seq(length(enz.all.list[[y]]))){
    tempvec=enz.all.list[[1]][[x]]*enz.all.list[[y]][[x]][[paste0("NCU05627_",y)]]
    names(tempvec)=paste0(names,"_",y)
    enzrealval.list[[y]][[x]]=tempvec
  }
}
rank=vector(length=length(list.res[[4]]))
temp=sapply(list.res[[4]],function(x){
    name=x[["name"]]
    ind=str_which(string=reac.rank,pattern=fixed(name))
    rank[ind]<<-x[["enzy"]]
})
boxplot_multi_exp(enzrealval.list,paste0(dir.res,"boxplot.enz.ini.all.multiexp",".",replicatei,".pdf"),
            rank=rank,
            xlab="enzymes",ylab="log_concentration",
            list.prior.part=NULL,logtrans=TRUE)

temptest=enzrealval.list[1]
names(temptest)="comb"
list.prior.pre[["enz"]][["comb"]]=list.prior.pre[["enz"]][["comb"]][enz,]
list.prior.pre[["enz"]][["comb"]][,1]=0.0000000001
list.prior.pre[["enz"]][["comb"]][,2]=0.1
boxplot_multi_exp_comb(temptest,paste0(dir.res,"boxplot.enz.ini.all.multiexp.single",".",replicatei,".pdf"),
            rank=rank,
            xlab="enzymes",ylab="log_concentration",
            list.prior.part=list.prior.pre[["enz"]],logtrans=TRUE)

# distribution of estimated kinetic parameters
# compare between different extend range in parameters
kine.all.list=list(km=vector(mode="list"),kcat=vector(mode="list"))
inikinelist=o02.data.accumu[["theta_reac_para"]]
# names=unique(str_replace_all(string=names(inikinelist[[1]]),pattern="\\_\\d+$",replacement=""))
names=reac.rank
len=length(names)
kine.km.list=vector(mode="list")
kine.kcat.list=vector(mode="list")
nameexp="comb"
kine.rank=c()
temp=sapply(inikinelist,simplify=FALSE,function(x){
        kcatvec=c()
        kmvec=c()
        indloc=1
        nameskine=c()
        sapply(seq(len),simplify=FALSE,function(y){
          # cat(y)
          if(rev[rank[y]]==1){
            kcatP=x[indloc]
            kcatS=x[indloc+1]
            kmS=x[indloc+2]
            kmP=x[indloc+3]
            kcatvec<<-c(kcatvec,kcatS,kcatP)
            kmvec<<-c(kmvec,kmS,kmP)
            indloc<<-indloc+4
            nameskine<<-c(nameskine,paste0(reac.rank[y],"_S"),paste0(reac.rank[y],"_P"))
          }else{
            kcat=x[indloc]
            km=x[indloc+1]
            kcatvec<<-c(kcatvec,kcat)
            kmvec<<-c(kmvec,km)
            indloc<<-indloc+2
            nameskine<<-c(nameskine,paste0(reac.rank[y],"_S"))
          }
        })
        names(kcatvec)=nameskine
        names(kmvec)=nameskine
        kine.km.list[[length(kine.km.list)+1]]<<-kmvec
        kine.kcat.list[[length(kine.kcat.list)+1]]<<-kcatvec
        kine.rank<<-nameskine
})
kine.all.list[["km"]][[nameexp]]=kine.km.list
kine.all.list[["kcat"]][[nameexp]]=kine.kcat.list
boxplot_multi_exp_comb(kine.all.list[["km"]],paste0(dir.res,"boxplot.km.all",".",replicatei,".pdf"),
            rank=kine.rank,
            xlab="enzymes",ylab="log_km",
            list.prior.part=list.prior.pre[["km"]],logtrans=TRUE)
boxplot_multi_exp_comb(kine.all.list[["kcat"]],paste0(dir.res,"boxplot.kcat.all",".",replicatei,".pdf"),
            rank=kine.rank,
            xlab="enzymes",ylab="log_kcat",
            list.prior.part=list.prior.pre[["kcat"]],logtrans=TRUE)

#vmax distribution
enz.list=enz.all.list
vmax.list=vector(mode="list")
model_len=length(kine.kcat.list)
rep_len=3
enz.rank=c()
temp=sapply(kine.rank,function(x){
  kine.inds=str_replace(string=x,pattern="_[SP]+$",replacement="")
  enz.rank<<-c(enz.rank,rank[which(reac.rank==kine.inds)])
})
temp=sapply(seq(model_len),function(x){
  sapply(seq(rep_len),function(y){
    vmax.list[[length(vmax.list)+1]]<<-kine.kcat.list[[x]][kine.rank]*enz.list[[1]][[x]][paste0(enz.rank,"_1")]*enz.list[[y]][[x]][[paste0("NCU05627_",y)]]
  })
})
##prior range
tablist=sapply(seq(dim(list.prior.pre[["kcat"]][["comb"]])[1]),simplify=FALSE,function(x){
  as.numeric(list.prior.pre[["kcat"]][["comb"]][x,])*c(0.0000000001,0.1)
})
list.prior.pre[["vmax"]][["comb"]]=Reduce("rbind",tablist)
rownames(list.prior.pre[["vmax"]][["comb"]])=rownames(list.prior.pre[["kcat"]][["comb"]])
vmax.list.list=vector(mode="list")
vmax.list.list[["comb"]]=vmax.list
boxplot_multi_exp_comb(vmax.list.list,paste0(dir.res,"boxplot.vmax.all",".",replicatei,".pdf"),
            rank=kine.rank,
            xlab="reaction",ylab="log_vmax",
            list.prior.part=list.prior.pre[["vmax"]],logtrans=TRUE)

#compare expected one difference in each comparison as it is from time
#skip first 100 bytes in pdf file
predatadir=predatadir=paste0(dirpack,"internal_data/pre_result_testplotting/")
pdflist=list.files(dir.res,pattern="*.pdf")
for(file in pdflist){
  cat(paste0(file,"\n"))
  system(paste0("cmp -i 100 \'",dir.res,file,"\' \'",predatadir,file,"\'"))
}
