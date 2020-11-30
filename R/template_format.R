#' formulate the reaction part of i01 based on format files
#'
#' this function is indicated to depends on exsting format/template for each part
#' lines start with more than two '#' will be removed
#' template is needed
#' no input arguments as arguments are transferred by environment
#'
#' @return list. input on file information and size
#' @seealso [spec_output_format()] for producing species block in i01 for default format
#'          [react_output()] for similar function based on default format
#' @import stringr
react_output_format<-function(){
  oldstr<-options(stringsAsFactors=FALSE)
  on.exit(options(oldstr),add=TRUE)
  ##environment change for digit
  oldscipen<-options(scipen=20)
  on.exit(options(oldscipen),add=TRUE)
  vals=para.list[["react"]][["initial"]]
  lines_exm=sapply(c("rev","irrev"),simplify=FALSE,function(x){
    lines=readLines(para.list[["react"]][["formatfile"]][[x]])
    lines[!str_detect(string=lines,pattern="^##")]
  })
  reac_output=sapply(list.reac.addon,function(x){
      flagmmrev=0
      rev=para.list[["react"]][["rev"]]
      infor=vector(mode="list")
      # infor[["nsubs"]]=length(x[["subs"]])
      # infor[["nprods"]]=length(x[["prods"]])
      # infor[["jkin"]]=jkin
      # infor[["name"]]=x[["name"]]
      name=x[["name"]]
      # print(name)
      ##reactants formulation
      reactants=c(x[["subs"]],x[["prods"]])
      reactants=reactants[!is.null(reactants)&reactants!=""]
      subs=x[["subs"]][x[["subs"]]%in%reactants]
      enzyme=c()
      indsub=subs%in%para.list[["react"]][["enz"]]
      enzyme=c(enzyme,subs[indsub])
      # subs=c(subs[indsub],subs[!indsub])
      # infor[["subs"]]=subs[!indsub]
      prods=x[["prods"]][x[["prods"]]%in%reactants]
      indprod=prods%in%para.list[["react"]][["enz"]]
      enzyme=c(enzyme,prods[indprod])
      # prods=c(prods[indprod],prods[!indprod])
      # infor[["prods"]]=prods[!indprod]
      enzyme=unique(enzyme)
      if(rev[[enzyme]]){
        template=lines_exm[["rev"]]
      }else{
        template=lines_exm[["irrev"]]
      }
      enzymeoutput=enzyme
      if(enzyme%in%reguinfor[["regutab"]][,"enz"]){
        enzymeoutput=paste0("regulated_",enzyme)
      }
      template=str_replace_all(string=template,pattern="EC",replacement=x[["name"]])
      template=str_replace_all(string=template,pattern="NSUBSTRATE",replacement=as.character(length(x[["subs"]])))
      template=str_replace_all(string=template,pattern="NPRODUCT",replacement=as.character(length(x[["prods"]])))
      template=str_replace_all(string=template,pattern="XSUBSTRATE",replacement=paste(subs[!indsub],collapse="\n   "))
      template=str_replace_all(string=template,pattern="XPRODUCT",replacement=paste(prods[!indprod],collapse="\n         "))
      template=str_replace_all(string=template,pattern="ENZYME",replacement=enzymeoutput)
      if(type=="mm"){
        parameter.kine.list=para.list[["react"]][["kine"]]
        kmdefault=0
        para=sapply(c("km","kcat"),simplify=FALSE,function(x){
          nameec=str_replace_all(string=name,pattern="\\_\\d+\\w+$",replacement="")
          tempara=parameter.kine.list[[x]][[nameec]]
          if(is.null(tempara)){
            temp=parameter.kine.list[[x]][["default"]]
            if(x=="km"){
              kmdefault<<-1
            }
          }else{
            temp=parameter.kine.list[[x]][[nameec]]
            # if(x=="km"&&length(unique(subs))>=3){
            #   temp=temp^2
            # }
          }
          temp[temp==0]=smallvalue
          temp
        })
        # kcat=para[["kcat"]]
        # km=para[["km"]]
        ## power the km to the extend of substrate number - 1
        if(kmdefault==0){#only power up km for multiple compound reaction for reactions with known km parameter
          para[["km"]]=para[["km"]]^(length(subs)-1)
        }
        if(str_detect(string=name,pattern="\\d+\\.\\d+\\.\\d+")||str_detect(string=name,pattern="^trans\\_")||str_detect(string=name,pattern="\\_enz$")||str_detect(string=name,pattern="^comb")){
          #kcat
          for(typeval in c("km","kcat")){
            template=str_replace_all(string=template,pattern=paste0(typeval,".*_min"),replacement=as.character(para[[typeval]][1]/extend))
            template=str_replace_all(string=template,pattern=paste0(typeval,".*_max"),replacement=as.character(para[[typeval]][2]*extend))
            local.para[[typeval]][[name]]<<-para[[typeval]]#store the value without extend
            local.para[[paste0(typeval,"r")]][[name]]<<-para[[typeval]]#store the value without extend
            namepara=paste0(name,".",typeval)
            if(namepara%in%names(vals)){
              val1=vals[namepara]
            }else if(rand){
              set.seed(rand.seed)
              val1=runif(1,min=para[[typeval]][1],max=para[[typeval]][2])
            }else{
              val1=sqrt(para[[typeval]][1]*para[[typeval]][2])#kcat[1]^(2/3)*kcat[2]^(1/3)#median(kcat)#sqrt(kcat[1]*kcat[2])
            }
            template=str_replace_all(string=template,pattern=paste0(typeval,".*_val\\d+"),replacement=as.character(val1))
          }
          reac.para.addon.len<<-reac.para.addon.len+5
          if(rev[[enzyme]]){
            reac.para.addon.len<<-reac.para.addon.len+2
          }
        }else{
          reac.para.addon.len<<-reac.para.addon.len+1
        }
      }else{
        ###the format need to be set here for the specific mass reaction pattern
        warning("the format of mass reaction need to be set!")
        if((str_detect(string=name,pattern="\\d+\\.\\d+\\.\\d+\\.\\d+")||str_detect(string=name,pattern="^trans\\_"))&&rev[enzyme]){
          reac.para.addon.len<<-reac.para.addon.len+2
        }else{
          reac.para.addon.len<<-reac.para.addon.len+1
        }
      }
      reac.count<<-reac.count+1
      template=template[!str_detect(string=template,pattern="^\\s+$")]
      output=paste(template,collapse="\n")
    })
    list.res=list(out=reac_output,len=reac.para.addon.len,count=reac.count)
    return(list.res)
}

#' formulate the species part of i01 based on format files
#'
#' this function is indicated to depends on exsting format/template for each part
#' template is needed
#' no input arguments as arguments are transferred by environment
#'
#' @return array. input on file information
#' @seealso [react_output_format()] for producing reaction block in i01 for default format
#'          [spec_output()] for similar function based on default format
#' @import stringr
spec_output_format<-function(){
  # print("****")
  # print(search())
  # vec_spec_addon
  oldstr<-options(stringsAsFactors=FALSE)
  on.exit(options(oldstr),add=TRUE)
  ##environment change for digit
  oldscipen<-options(scipen=20)
  on.exit(options(oldscipen),add=TRUE)
  list.inival=vector(mode="list")
  list.inival[["low.all.spec"]]=para.list[["species"]][["low.all.spec"]]
  list.inival[["high.all.spec"]]=para.list[["species"]][["high.all.spec"]]
  list.inival[["low.ini.spec"]]=para.list[["species"]][["low.ini.spec"]]
  list.inival[["high.ini.spec"]]=para.list[["species"]][["high.ini.spec"]]
  vals=para.list[["species"]][["initial"]]
  lines_exm=sapply(c("metabolites","enzyme"),simplify=FALSE,function(x){
    lines=readLines(para.list[["species"]][["formatfile"]][[x]])
    lines[!str_detect(string=lines,pattern="^##")]
  })
  parareg=para.list[["regu"]][["parameter"]]
  spec_output=sapply(vec_spec_addon,function(x){
      # print(x)
      ###regulation format?
      templatetrans=c()
      if(x%in%reguinfor[["regutab"]][,"enz"]){
        repcompreg=reguinfor[["reguformat"]][["compd_regu"]]
        enzreg=reguinfor[["reguformat"]][["regu_enz"]]
        lenenzreg=length(enzreg)
        repenzreg=enzreg[(lenenzreg-2):(lenenzreg-1)]
        enzreg=enzreg[1:(lenenzreg-3)]
        compds=reguinfor[["regutab"]][which(reguinfor[["regutab"]][,"enz"]==x),"compd"]
        compdi=2
        sumtempl=enzreg
        sumtempl=str_replace_all(string=sumtempl,pattern="regu_enz",replacement=paste0("regulated_",x))
        sumtempl=str_replace_all(string=sumtempl,pattern="reg_low",replacement=paste0(parareg[["reg_low"]]))
        sumtempl=str_replace_all(string=sumtempl,pattern="reg_high",replacement=paste0(parareg[["reg_high"]]))
        comptempl=c()
        for(compd in compds){
          reppartch=repenzreg
          reppartch=str_replace_all(string=reppartch,pattern="compd_regu",replacement=paste0(compd,"_regu_",x))
          reppartch=str_replace_all(string=reppartch,pattern="^\\s+2\\s+",replacement=paste0("    ",compdi," "))
          compdi=compdi+1
          sumtempl=c(sumtempl,reppartch)
          repparttransch=repcompreg
          repparttransch=str_replace_all(string=repparttransch,pattern="compd_regu",replacement=paste0(compd,"_regu_",x))
          repparttransch=str_replace_all(string=repparttransch,pattern="compd",replacement=compd)
          repparttransch=str_replace_all(string=repparttransch,pattern="reg_low",replacement=paste0(parareg[["reg_low"]]))
          repparttransch=str_replace_all(string=repparttransch,pattern="reg_high",replacement=paste0(parareg[["reg_high"]]))
          regutabhere=reguinfor[["regutab"]]
          regustat=regutabhere[regutabhere[,"enz"]==x&regutabhere[,"compd"]==compd,"regu"]
          if(regustat==1){
            alpharange=c(0,parareg[["alpha_high"]])
          }else{
            alpharange=c(parareg[["alpha_low"]],0)
          }
          repparttransch=str_replace_all(string=repparttransch,pattern="alpha_low",replacement=paste0(alpharange[1]))
          repparttransch=str_replace_all(string=repparttransch,pattern="alpha_high",replacement=paste0(alpharange[2]))
          alphaval_ini=parareg[["alpha_val"]]
          if(!(alphaval_ini>alpharange[1]&&alphaval_ini<alpharange[2])){
            if((-alphaval_ini>alpharange[1]&&(-alphaval_ini<alpharange[2]))){
              alphaval_ini=-alphaval_ini
            }else{
              alphaval_ini=mean(alpharange)
            }
          }
          repparttransch=str_replace_all(string=repparttransch,pattern="alpha_val\\d",replacement=paste0(alphaval_ini))
          comptempl=c(comptempl,repparttransch)
        }
        # template=c(template,templatetrans)
        templatetrans=c(comptempl,sumtempl,"#")
        templatetrans=str_replace_all(string=templatetrans,pattern="NPART",replacement=paste0(compdi-1))
        templatetrans=str_replace_all(string=templatetrans,pattern="ENZYME",replacement=x)
      }
      if(x%in%para.list[["react"]][["enz"]]){
        template=lines_exm[["enzyme"]]
        template=str_replace_all(string=template,pattern="SCALENZ",replacement=paste0(para.list[["species"]][["refenz"]]))
        template=str_replace_all(string=template,pattern="ENZYME",replacement=paste0(x))
      }else{
        template=lines_exm[["metabolites"]]
        template=str_replace_all(string=template,pattern="SPEC",replacement=paste0(x))
      }
      template=c(template,templatetrans)
      ##obsv
      obj=0
      if(x%in%para.list[["species"]][["obsv"]]){
        if(!is.null(jmsspec)){
          obj=jmsspec
        }else{
          obj=1
        }
      }
      template=str_replace_all(string=template,pattern="OBSERVE",replacement=paste0(obj))
      ##const
      const=0
      if(x%in%para.list[["species"]][["const"]]){
        const=1
      }
      parts=str_split(template[str_which(template,pattern="^\\s+\\d+\\s+")[1]],pattern=" ")[[1]]
      replapart=parts[parts!=""][1]
      const=as.numeric(replapart)+const
      template=str_replace_all(string=template,pattern=paste0("^\\s+",replapart),replacement=paste0("      ",const))
      ##initial value
      template=str_replace_all(string=template,pattern="ranglow",replacement=paste0(list.inival[["low.all.spec"]][x]))
      template=str_replace_all(string=template,pattern="ranghigh",replacement=paste0(list.inival[["high.all.spec"]][x]))
      template=str_replace_all(string=template,pattern="expenz\\d_high",replacement=paste0(list.inival[["high.ini.spec"]][x]))
      template=str_replace_all(string=template,pattern="expenz\\d_low",replacement=paste0(list.inival[["low.ini.spec"]][x]))
      local.para[["enz"]][[x]]<<-c(list.inival[["low.ini.spec"]][x],list.inival[["high.ini.spec"]][x])
      if(x%in%names(vals)){
        template=str_replace_all(string=template,pattern="expenz\\d_val\\d",replacement=paste0(vals[x]))
        # if(fixed==1){
        #   ###
        # }
      }else if(rand){
        set.seed(rand.seed)
        valrand=runif(1,min=list.inival[["low.ini.spec"]][x],max=list.inival[["high.ini.spec"]][x])
        template=str_replace_all(string=template,pattern="expenz\\d_val\\d",replacement=paste0(valrand))
      }else{
        valrand=sqrt(list.inival[["low.ini.spec"]][x]*list.inival[["high.ini.spec"]][x])
        template=str_replace_all(string=template,pattern="expenz\\d_val\\d",replacement=paste0(valrand))
      }
      template=template[!str_detect(string=template,pattern="^\\s+$")]
      output=paste(template,collapse="\n")
    })
    return(spec_output)
}
