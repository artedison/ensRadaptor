## test workflow for ensRadaptor
## this testing script construct input files for modeling central metabolism in Neurospora crassa based on experimental data.
## the output file should be fed to ens modeling for ODE parameter estimation
## works with sapelo2, template files need to be modified for adaption to other sever systems
## unit time: hour h
## unit concentration: mM
## kinetic form: MM (assumed)
## optimized parameter: kcat, km, E
## multiple folders will be created and "input" is the converted input files for ens and "res" is the folder for output from ens. In "input", mulitple different random seed is used
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15,scipen=999)
require(stringr)
require(magrittr)
# require(ggplot2)
# require(cowplot)
# require(scales)
# require(lhs)
require(ensRadaptor)

#------------------------------path setup-------------------#
currdir=getwd()
dir=paste0(currdir,"/temp/input_test/")## the working folder. User defined
dir_ext_data=system.file("extdata","", package = "ensRadaptor")#"../inst/extdata/"## the folder for needed data files. User defined
if(!str_detect(string=dir_ext_data,pattern="\\/$")){
  dir_ext_data=paste0(dir_ext_data,"/")
}
dir.lib=paste0(dir_ext_data,"template_format/")##this path is used within many functions. User defined
dir.his=paste0(dir,"history/history.record.tab")##recording for running history. The time of running and project folder. User defined
if(!dir.exists("./temp")){
  foldcreate("./temp")
}
if(!dir.exists(dir)){
  foldcreate(dir)
}
if(!dir.exists(paste0(dir,"history/"))){
  foldcreate(paste0(dir,"history/"))
}
if(!file.exists(dir.his)){
  file.create(dir.his)
}
dir.data=paste0(dir,"testmodel/")
foldcreate(paste0(dir.data))
dir.template.format=paste0(dir_ext_data,"template_format/")
foldname="1" ##the model directory
foldcreate(paste0(dir.data,foldname))
update_history(histpath=dir.his,projpath=paste0(dir.data,foldname,"/"))
dir.res=paste0(dir.data,foldname,"/res/")
foldcreate(dir.res)
dir.tempt=paste0(dir.data,foldname,"/tempt/")
foldcreate(dir.tempt)
dir.input=paste0(dir.data,foldname,"/input/")
foldcreate(dir.input)
dir.equala=paste0(dir.res,"equala/")
foldcreate(dir.equala)
dir.accumu=paste0(dir.res,"accumu/")
foldcreate(dir.accumu)
##create the multiple run folder for independent start point
nreplicate=10##number of model replicate to run in the model (independent MCMC runs). User defined
for(replicatei in seq(nreplicate)){
  foldcreate(paste0(dir.input,replicatei,"/"))
}
## make sure the output anlaysis file agree with the input production file
# localpath=paste0(dir.data,foldname,"/")
# homoge_file(localpath,foldname)

#------------------------------parameter setup--------------#
modified.file="testmodel"#model name. User defined
##model parameter
obsv=c("gluc_all","ethanol","lactate","succinate","citrate","glu_1_phos","fumarate") ##observable compounds (in the measurement). all other compounds defined in the pathway file will be unobservable. User defined
time1=12# the end time for simulation, time0 is assumed to be 0 but you can always shift your measuremnt to make start time to 0. User defined
extend=1#extend of all parameters range. For those kinetic parameters that a known range exist, the range will be extended by dividing (multiplying) this value for lower (upper) range. User defined
extendrag=100#extend of parameters without range (no informaiton from database or just one value). User defined
smallvalue<-0.0000000000000000001#used to replace 0 when math on 0 is not working. used within many functions. should be fine under different models
refenz="NCU05627"#the enzyme used to scale different experimental globally(as they are with different concentration). Choose the first enzyme as the reference but it doesn't matter. User defined

##measurement parameter
classnum=1# start scale class for measurement (concept in ens program). the class number for measurement. if not 0 the scale class, if 0 not using scale factor. The same classs number define a group of measurements that change relative to each other. User defined
conv.factor=c(1,3600)#unit converting factor for km and kcat. User defined
names(conv.factor)=c("km","kcat")
zspec_wid=0.3# the default relative uncerntainty for measurement. User defined
experiments=c(4:6)#4,5,6anaerobic, 1,2,3 aerobic. User defined

lhsmatrix<-c()#global variable as used by a function multiple times. used for randloglhs function
# lowthrespec=0.00000000000000001
# highthrespec=100.00006

#----------------------construt template setup--------------#
###the template files are used for easy construction of different reactions and species
templatepath=list(input.template=paste0(dir_ext_data,"template_input/"),enscode=paste0(dir_ext_data,"enscode/"))#User defined
pre_file_prepare(dir.tempt,templatepath)
specformat=list(metabolites=paste0(dir.template.format,"spec.metabolite.tab"),enzyme=paste0(dir.template.format,"spec.enz.massscal.tab"))#User defined
reacformat=list(rev=paste0(dir.template.format,"reac.rev.mm.tab"),irrev=paste0(dir.template.format,"reac.irrev.mm.tab"))#User defined

#----------------------pathway setup -----------------------#
pathfile=paste0(dir_ext_data,"pathway/testpathway.tab")#pathway decide the species and reactions. User defined
# the wrapper invisible(capture.output()) to stop print out. In practice, reading the output might be really helpful to obtain information on the model
invisible(capture.output(list.res<-summary_reac(pathfile,"(^NCU)|(comp$)")))#User defined
enz=list.res$enz

meta.real=list.res$compounds##metabolites
meta.comb=c("gluc_all")## addon species to represent a combination of intracellular and extracellular part. User defined
meta.list=c(meta.real,meta.comb)

sub_var=c(enz,meta.real)##those species that do not depends on other species

reacs_dependent=list.res[["reacs"]]
names(reacs_dependent)=NULL
const=enz##not change through time

##reversibility list.  User defined
rev=unlist(list(
  "NCU05627"=0,"NCU00575"=0,"NCU07281"=1,"NCU00629"=0,
  "NCU06075"=0,"NCU02193"=0,"NCU02476"=1,"NCU02179"=0,
  "NCU07659"=0,"NCU01692"=0,"NCU00775"=0,"2oxodehycomp"=0,
  "NCU09810"=1,"NCU00959"=0,"NCU10008"=1,"NCU04899"=1,
  "NCU09873"=1,"NCU02505"=1,"NCU02906"=0,"NCU10058"=1,
  "glycogen_syn_comp"=0,"glycogen_deg_comp"=0,"fattyacid_syn_comp"=0,"fattyacid_deg_comp"=0,"NCU04230"=0,"NCU10007"=0,
  "enz_glycolysis_comp"=1,"enz_tca_comp"=1,"NCU04797"=0,
  "oxy_phos_enz_comp"=0,"gluc1p_phosase_enz_comp"=0
))

#----------------------kinetic parameters-------------------#
##stored kinetic parameters for enzyme. User defined
load(paste0(dir_ext_data,"kin_infor/Neurospora crassa.kine.new.RData"))
load(paste0(dir_ext_data,"kin_infor/all.kine.RData"))
parameter.list=kine_para_refine(list.res.refine=list.res.refine,range.speci=range.speci,extendrag=extendrag)

##further refinement of the value range. User defined
##for km: 2.7.1.1 is proper choice as it is with larger range than 2.7.1.2. Both are used in glucose->glu-6-p
parameter.list[["kcat"]][["2.7.1.1"]]=parameter.list[["kcat"]][["2.7.1.2"]]##choose larger one
parameter.list[["kcat"]][["comb_glycolysis"]]=parameter.list[["kcat"]][["4.1.2.13"]]#choose the slowest one kcat
parameter.list[["km"]][["comb_glycolysis"]]=parameter.list[["km"]][["4.1.2.13"]]#the first reaction km for forward km and last reaction km for backward km
parameter.list[["km"]][["comb_tca"]]=parameter.list[["km"]][["4.2.1.3"]]##the same enzyme for two step and use the same km

##for lactate fermentation there are both R and S type with different EC, the larger range is chosen here
parameter.list[["km"]][["1.1.1.28"]]=c(2.0e-04,6.7e+02)
parameter.list[["kcat"]][["1.1.1.28"]]=c(3.6e+02,1.8e+07)

##for malate->pyruvate there might be NADPH or NADH use the larger ones
parameter.list[["km"]][["1.1.1.38"]]=c(0.00139,45.00000)
parameter.list[["kcat"]][["1.1.1.38"]]=c(16.56,6631200000)

##for oxidative phosphorylation reaction
parameter.list[["km"]][["oxy_phos"]]=c()
parameter.list[["kcat"]][["oxy_phos"]]=c()

len.metab=length(sub_var)

##dependence
###the dependence from different experiments
# specdep=data.frame(name=sub_var,y=rep(1,times=len.metab),
#                    x1=rep(2,times=len.metab),x2=rep(3,times=len.metab),
#                    x3=rep(4,times=len.metab),
#                    func=rep("LOOK4_PM_IEXPT",times=len.metab))
# len.reacs=length(reacs_dependent)

#-----------------species parameter range-------------------#
changedrag=c(obsv,"gluc_in","glycogen","fattyacid","Gluc_ex")#User defined
ini.names=c(meta.list,enz)
meta.val=rep(0.00000001,times=length(meta.list))
names(meta.val)=meta.list
meta.val[changedrag]=rep(0.00001,times=length(changedrag))
meta.val[c("Gluc_ex","gluc_all")]=c(200,200)
enz.val=rep(0.0000000001,times=length(enz))
names(enz.val)=enz
low.ini.spec=c(unlist(meta.val),enz.val)#User defined
#
meta.val=rep(1,times=length(meta.list))
names(meta.val)=meta.list
meta.val[changedrag]=rep(2000,times=length(changedrag))
meta.val[c("gluc_in")]=c(20)
enz.val=rep(0.1,times=length(enz))
names(enz.val)=enz
high.ini.spec=c(unlist(meta.val),enz.val)#User defined
#
low.all.spec=c(rep(0.000000000001,times=length(meta.list)),
              rep(0.000000000001,times=length(enz)))#User defined
names(low.all.spec)=ini.names
#
meta.val=rep(2,times=length(meta.list))
names(meta.val)=meta.list
meta.val[changedrag]=rep(2000,times=length(changedrag))
enz.val=rep(0.1,times=length(enz))
names(enz.val)=enz
high.all.spec=c(unlist(meta.val),enz.val)#User defined

#---------------------regulation values---------------------#
reg_rang=c(0.000000000001,2000)#c(low high). User defined
alpha_vec=c(-5,5,1)#c(low,high,value). User defined

#---------------------pack up values-----------------------#
regupara=list(reg_low=reg_rang[1],reg_high=reg_rang[2],alpha_low=alpha_vec[1],alpha_high=alpha_vec[2],alpha_val=alpha_vec[3])
spe.list=list(obsv=obsv,const=const,
              low.ini.spec=low.ini.spec,high.ini.spec=high.ini.spec,
              low.all.spec=low.all.spec,high.all.spec=high.all.spec,
              format=TRUE,
              formatfile=specformat,
              refenz=refenz)#,fixed=fixed,initial=initial)
reac.list=list(kine=parameter.list,enz=enz,rev=rev,
               format=TRUE,
               formatfile=reacformat)

##this is previous setting parameters for dependence, which turn out not flexible enough for new condition, so not used anymore. This block of code and function are still in script (NULL) and functions but will later be discarded.
reactdep=NULL
specdep=NULL
depend=list(spec=specdep,react=reactdep)
allosteric.regu=list(listfile=NULL,###or NULL
                     format=NULL,
                     parameter=regupara)
para.list=list(species=spe.list,react=reac.list,depend=depend,regu=allosteric.regu)

## file path
len.list=c()
temp.change.file1=paste0(dir.tempt,"ens.modified.i01")
temp.change.file2=paste0(dir.tempt,"ens.modified2.i01")
temp.change.file3=paste0(dir.tempt,"ens.modified3.i01")
ens.def=paste0(dir.tempt,"ens.def")
ens.f90=paste0(dir.tempt,"ens.f90")
ens.sh=paste0(dir.tempt,"ens.sh")
submit.sh=paste0(dir.tempt,"submit.sh")

#----------------------i01 construct------------------------#
len.list.temp<-template_spec_reac(pathfile,
                   type="mm",##Michaelis Menten reaction
                   dir.data=dir.tempt,modified.file=modified.file,
                   para.list=para.list,list.exi=c(),extend=extend)
##only len.list produced by one loop is needed as they are the same
if(is.null(len.list)){
  len.list=len.list.temp
}

##addon blocks User defined
##!change backward km (ipm=6) for comb_glycolysis as 4.2.1.11
# 0.0018 7.5000 0.116189500386223
change_block(list(pattern0="namereac\\s+\\/nipart\\s+nopart\\sjkin\\n\\s+comb_glycolysis",pattern1="namereac\\s+\\/nipart\\s+nopart\\s+jkin\\n\\s+2.7.1.40",shift=c(26,-31)),
             list(file=paste0(dir.template.format,"comb_glycolysis.tab"),ilinerag=c(31,35)),
             infile=temp.change.file1,outfile=temp.change.file1,type="edit")
## add on the combined measrument{addon.comb.meas.tab}
### gluc_all=Gluc_ex + gluc_in
change_block(list(pattern0="^\\#\\s+Reaction\\s+control\\-\\s+and\\s+\\\\Theta\\-variables",pattern1="^\\#\\s+Reaction\\s+control\\-\\s+and\\s+\\\\Theta\\-variables",shift=c(-1,-2)),
             list(file=paste0(dir.template.format,"addon.comb.meas.tab"),ilinerag=c(76,99)),
             infile=temp.change.file1,outfile=temp.change.file1,type="edit")
## modify the first enzyme to be scale factor {addon.scalenz.tab}
change_block(list(pattern0="namespec\\/jfix\\s+npulse\\s+nperiod\\s+jmsspec\\s+xspec\\_min xspec\\_max\\n\\s+NCU05627",pattern1="namespec\\/jfix\\s+npulse\\s+nperiod\\s+jmsspec\\s+xspec\\_min\\s+xspec\\_max\\n\\s+gluc_in",shift=c(-1,-3)),
              list(file=paste0(dir.template.format,"addon.scalenz.tab"),ilinerag=c(66,127)),
              infile=temp.change.file1,outfile=temp.change.file1,type="edit")
## add differences on oxidative phosphorylation for anaerobic and aerobic conditions {addon.aero.vs.anaero.tab}
change_block(list(pattern0="namereac\\s+\\/nipart\\s+nopart\\s+jkin\\n\\s+oxy\\_phos\\_enz1",pattern1="namereac\\s+\\/nipart\\s+nopart\\s+jkin\\n\\s+gluc1p\\_phosase\\_enz",shift=c(-1,-2)),## this line need modification
             list(file=paste0(dir.template.format,"addon.aero.vs.anaero.tab"),ilinerag=c(245,392)),
             infile=temp.change.file1,outfile=temp.change.file1,type="edit")

## change dimension parameters. User defined
nspec=length(change_para("namespec",NA,infile=temp.change.file1,outfile=NA,type="show")[[1]])
ndepn=length(change_para("iparn",NA,infile=temp.change.file1,outfile=NA,type="show")[[1]])
list.i01.para=list(nspec=nspec,nreac=len.list[["reac"]],time1=time1,nexpt=length(experiments),"th\\_opt1"=102010001,##supress unuseful information in o02 file
                   "jmc\\_ini"=40)
system(paste0("cp \"",temp.change.file1,"\" \"",temp.change.file2,"\""))
for(para.name in names(list.i01.para)){
  change_para(para.name,list.i01.para[[para.name]],infile=temp.change.file2,outfile=temp.change.file2,type="edit")
}

#----------------------def modificate-----------------------#
# User defined
list.def.para=list("mp\\_dim\\_yspec\\_x"=120000,"iline\\_x"="200+2*nspec_x+2*nreac_x+5500",## the two array size parameter that seldom need modification
                   "npart\\_x"=8,"nms\\_spec\\_x"=18,
                   "nspec\\_x"=nspec,"nreac\\_x"=len.list[["reac"]],"ndpen\\_tot\\_x"=ndepn)
for(paraname in names(list.def.para)){
  change_def(paraname,list.def.para[[paraname]],
             infile=ens.def,outfile=ens.def,type="edit")
}

#-------------------submit.sh modificate--------------------#
# User defined
lines=readLines(submit.sh)
lines=str_replace(string=lines,pattern="N\\_REPLICATE",replacement=as.character(nreplicate))
cat(lines,sep="\n",file=submit.sh)
system(paste0("cp \"",submit.sh,"\" \"",dir.input,"submit.sh\""))

#-------------------replicates construct---------------------#
for(replicatei in seq(nreplicate)){
  print(paste0("replicate: ",replicatei))
  ##Latin Hypercube log uniform initial guess reformulating
  randloglhs(nreplicate,replicatei,input=temp.change.file2,output=temp.change.file3)
  change_para("imc_ran",9518730-replicatei+1,infile=temp.change.file3,outfile=temp.change.file3,type="edit")
  change_para("iseed2",replicatei+1-1,infile=temp.change.file3,outfile=temp.change.file3,type="edit")
  system(paste0("cp \"",ens.def,"\" \"",dir.input,replicatei,"/ens.def\""))
  system(paste0("cp \"",ens.f90,"\" \"",dir.input,replicatei,"/ens.f90\""))
  system(paste0("cp \"",temp.change.file3,"\" \"",dir.input,replicatei,"/ens.i01\""))
  system(paste0("cp \"",ens.sh,"\" \"",dir.input,replicatei,"/ens.sh\""))
}
savedrdata=paste0(dir.res,"local.prior.pre.",format(Sys.time(), "%H%M%S_%m%d%Y"),".RData")
save(len.list,rev,meta.list,list.res,experiments,obsv,modified.file,foldname,dir.res,enz,refenz,file=savedrdata)

#-----------------------i02 construct-----------------------#
load(paste0(dir_ext_data,"measurements/compound.quan.record.absmore.upd.RData"))
## manually add the experiment part
classnumfix=classnum
names=names(list.quan)
names[names=="glucose-1-phosphate"]="glu_1_phos"
names[names=="glucose"]="gluc_all"
# setdiff(obsv,names)
names(list.quan)=names
if(exists("list.rela.uncerntainty")){
  names=names(list.rela.uncerntainty)
  names[names=="glucose-1-phosphate"]="glu_1_phos"
  names[names=="glucose"]="gluc_all"
  names(list.rela.uncerntainty)=names
  zspec_wid=list.rela.uncerntainty
}

i02path=paste0(dir.tempt,"ens.i02")
for(exp in experiments){
  list.nmr=sapply(obsv,simplify=FALSE,function(comp){
    tempmat=list.quan[[comp]]
    tempmat[tempmat[,"exp"]==exp,]
  })
  list.nmr=list.nmr[!sapply(list.nmr,function(y){dim(y)[1]==0})]
  # print(sapply(list.nmr,dim))
  localblockpath=paste0(dir.tempt,"addon.measure.nmr.",exp,".tab")
  classnum<-template_exp_data(list.nmr,classnum=classnumfix,
                              output=localblockpath,
                              dir.data=dir.tempt,modified.file=modified.file,
                              zspec_wid=zspec_wid,sameclass=1)
  patternmatch=paste0("^\\s+",exp-3,"\\s+18\\s+$")
  nlineblock=linecount(localblockpath)
  change_block(list(pattern0=patternmatch,pattern1=patternmatch,shift=c(2,2)),
               list(file=localblockpath,ilinerag=c(1,nlineblock)),
               infile=i02path,outfile=i02path,type="edit")
}

##modify i02 measurement number
lines=readLines(i02path)
indnum=str_which(string=lines,pattern="^\\s+\\d+\\s+\\d+\\s+$")
lines[indnum] %>% str_split(string=.,pattern="\\s+",simplify=TRUE) -> numtab
numtab[,3]=as.character(length(obsv))
lines[indnum]=apply(numtab,1,function(x){paste(x,collapse="     ")})
cat(lines,sep="\n",file=i02path)

for(replicatei in seq(nreplicate)){
  system(paste0("cp \"",i02path,"\" \"",dir.input,replicatei,"/ens.i02\""))
}

#-----------------------summary input-----------------------#
# the wrapper invisible(capture.output()) to stop print out. In practice, reading the output might be really helpful to obtain information on the model
invisible(capture.output(list.exi<-summary_input(i01=paste0(dir.tempt,"ens.modified3.i01"),i02=paste0(dir.tempt,"ens.i02"))))
##transfer the folder to sever, cd into the folder
##!a small local interactive run(if debug clear the whole folder)
##chmod +x submit.sh
##./submit.sh
