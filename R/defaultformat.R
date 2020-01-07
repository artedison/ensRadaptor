#' formulate the species part of i01(in a default format)
#'
#' formulate the species part for i01 file with default format
#' the argument is not set explicit and they are transferred by environment (the function that call spec_output())
#'
#' @return array. output information
#' @seealso [spec_output_format()] for similar function based on user defined format files
#'          [react_output()] for producing reaction block in i01 for default format
spec_output<-function(){
  ##recover environment
  old<-options(stringsAsFactors=FALSE)
  on.exit(options(old),add=TRUE)
  ##environment change for digit
  options(scipen=20)
  list.inival=vector(mode="list")
  list.inival[["low.all.spec"]]=para.list[["species"]][["low.all.spec"]]
  list.inival[["high.all.spec"]]=para.list[["species"]][["high.all.spec"]]
  list.inival[["low.ini.spec"]]=para.list[["species"]][["low.ini.spec"]]
  list.inival[["high.ini.spec"]]=para.list[["species"]][["high.ini.spec"]]
  vals=para.list[["species"]][["initial"]]
  spec_outputreturn=sapply(vec.spec.addon,function(x){
    lineannonum2=spelineannonum
    ## settings
    if(x%in%para.list[["species"]][["obsv"]]){
      if(!is.null(jmsspec)){
        lineannonum2[4]=jmsspec
      }else{
        lineannonum2[4]=1
      }
    }
    if(x%in%para.list[["species"]][["const"]]){
      lineannonum2[1]=1
    }
    fixed=0
    if(x%in%para.list[["species"]][["fixed"]]){
      spelinenum[[1]][4]=100
      fixed=1
    }
    ###initial value
    if(x %in% names(list.inival[["low.all.spec"]])){
      lineannonum2[5]=list.inival[["low.all.spec"]][x]
    }
    if(x %in% names(list.inival[["high.all.spec"]])){
      lineannonum2[6]=list.inival[["high.all.spec"]][x]
    }
    if(x %in% names(list.inival[["low.ini.spec"]])){
      spelinenum[[1]][3]=list.inival[["low.ini.spec"]][x]
    }
    if(x %in% names(list.inival[["high.ini.spec"]])){
      spelinenum[[2]][3]=list.inival[["high.ini.spec"]][x]
    }
    local.para[["enz"]][[x]]<<-c(spelinenum[[1]][3],spelinenum[[2]][3])#store the real output value
    # init.rag=para.list[["species"]][["initial"]]
    if(x%in%names(vals)){
      spelinenum[[3]][3]=vals[x]
      if(fixed==1){
        spelinenum[[1]][3]=vals[x]
      }
    }else if(rand){
      set.seed(rand.seed)
      spelinenum[[3]][3]=runif(1,min=spelinenum[[1]][3],max=spelinenum[[2]][3])
    }else{
      spelinenum[[3]][3]=sqrt(spelinenum[[1]][3]*spelinenum[[2]][3])
    }
    spelinenum[[4]][3]=spelinenum[[3]][3]
    spelinenum[[5]][3]=spelinenum[[3]][3]
    if(x%in%para.list[["depend"]][["spec"]][,"name"]){
      temptab=para.list[["depend"]][["spec"]]
      spelinenumpre=spelinenum
      spelinenumpre[[1]][4]=10
      # spelinenumpre[[1]][5]=10001
      coln=dim(temptab)[2]
      lineannonum2[1]=as.numeric(lineannonum2[1])+10000*(coln-3)
      temptab=temptab[temptab[,"name"]==x,]
      spelinenumadd=spelinenum
      spelinenum=vector(mode="list")
      spelinenum[[1]]=spelinenumpre[[1]]
      for(i in 1:(coln-3)){
        spelinenumadd=sapply(spelinenumadd,simplify=FALSE,function(y){
          y[1]=y[1]+1
          y
        })
        spelinenum=c(spelinenum,spelinenumadd)
      }
      depenc=depend_for(temptab[,1],unlist(temptab[,3:(coln-1)]),temptab[,coln])
    }
    if(fixed==1){
      string=c(paste(spelineanno1,collapse=" "),paste0(" ",x),
               paste(" ",lineannonum2,collapse=" "),
               paste(" ",spelineanno2,collapse=" "),
               paste(" ",spelinenum[[1]],collapse="   "))
    }else{
      string=c(paste(spelineanno1,collapse=" "),paste0(" ",x),
               paste(" ",lineannonum2,collapse=" "),
               paste(" ",spelineanno2,collapse=" "),
               unlist(sapply(spelinenum,function(y){
                 paste(" ",y,collapse="   ")
               })))
    }
    if(x%in%para.list[["depend"]][["spec"]][,"name"]){
      string=c(string,"#",depenc)
    }
    output=paste(string,collapse="\n")
  })
  return(spec_outputreturn)
}

#' formulate the reaction part of i01 (in a default format)
#'
#' formulate the species part for i01 file with default format
#' the argument is not set explicit and they are transferred by environment (the function that call react_output())
#'
#' @return list. output information
#' @seealso [react_output_format()] for similar function based on user defined format files
#'          [spec_output()] for producing species block in i01 for default format
react_output<-function(){
  ##recover environment
  old<-options(stringsAsFactors=FALSE)
  on.exit(options(old),add=TRUE)
  ##environment change for digit
  options(scipen=20)
  vals=para.list[["react"]][["initial"]]
  reac_outputreturn=sapply(list.reac.addon,function(x){
      flagmmrev=0
      lineannonum2=c(length(x[["subs"]]),length(x[["prods"]]),jkin)
      name=x[["name"]]
      print(name)
      reactants=c(x[["subs"]],x[["prods"]])
      reactants=reactants[!is.null(reactants)&reactants!=""]
      subs=x[["subs"]][x[["subs"]]%in%reactants]
      enzyme=c()
      indsub=subs%in%para.list[["react"]][["enz"]]
      enzyme=c(enzyme,subs[indsub])
      subs=c(subs[indsub],subs[!indsub])
      prods=x[["prods"]][x[["prods"]]%in%reactants]
      indprod=prods%in%para.list[["react"]][["enz"]]
      enzyme=c(enzyme,prods[indprod])
      prods=c(prods[indprod],prods[!indprod])
      enzyme=unique(enzyme)
      reactionnum=c()
      depenc.reac=NULL
      if(type=="mm"){
        parameter.kine.list=para.list[["react"]][["kine"]]
        rev=para.list[["react"]][["rev"]]
        # revflag=para.list[["react"]][["revflag"]]
        para=sapply(c("km","kcat"),simplify=FALSE,function(x){
          nameec=str_replace_all(string=name,pattern="\\_\\d+\\w+$",replacement="")
          tempara=parameter.kine.list[[x]][[nameec]]
          if(is.null(tempara)){
            temp=parameter.kine.list[[x]][["default"]]
          }else{
            temp=parameter.kine.list[[x]][[nameec]]
            # if(x=="km"&&length(unique(subs))>=3){
            #   temp=temp^2
            # }
          }
          temp[temp==0]=smallvalue
          temp
        })
        kcat=para[["kcat"]]
        km=para[["km"]]
        ## power the km to the extend of substrate number - 1
        km=km^(length(subs)-1)
        # list.kine.prior[["kcat"]][[length(list.kine.prior[["kcat"]])+1]]<<-kcat
        # list.kine.prior[["km"]][[length(list.kine.prior[["km"]])+1]]<<-km
        fixedpara=paste(reaclinenum[[6]],collapse=" ")
        temptabreac=para.list[["depend"]][["react"]]
        coln=dim(temptabreac)[2]
        if(str_detect(string=name,pattern="\\d+\\.\\d+\\.\\d+")||str_detect(string=name,pattern="^trans\\_")||str_detect(string=name,pattern="\\_enz$")){
          #the dependent line k1
          reaclinenum[[6]][1]=1
          reaclinenum[[6]][4]=10
          reaclinenum[[6]][3]=0.001
          fixedpara1=paste(reaclinenum[[6]],collapse="     ")
          #k-1
          reaclinenum[[6]][1]=2
          reaclinenum[[6]][4]=100
          reaclinenum[[6]][3]=0
          fixedpara2=paste(reaclinenum[[6]],collapse="     ")
          #k2 kcat
          sapply(1:length(reaclinenum),function(x){
            reaclinenum[[x]][1]<<-3
          })
          reaclinenum[[1]][3]=kcat[1]/extend
          reaclinenum[[2]][3]=kcat[2]*extend
          local.para[["kcat"]][[name]]<<-kcat#store the value without extend
          namepara=paste0(name,".kcat")
          if(namepara%in%names(vals)){
            reaclinenum[[3]][3]=vals[namepara]
          }else if(rand){
            set.seed(rand.seed)
            reaclinenum[[3]][3]=runif(1,min=kcat[1],max=kcat[2])
          }else{
            reaclinenum[[3]][3]=sqrt(kcat[1]*kcat[2])#kcat[1]^(2/3)*kcat[2]^(1/3)#median(kcat)#sqrt(kcat[1]*kcat[2])
          }
          reaclinenum[[4]][3]=reaclinenum[[3]][3]#kcat[1]^(2/3)*kcat[2]^(1/3)#median(kcat)#sqrt(kcat[1]*kcat[2])
          reaclinenum[[5]][3]=reaclinenum[[3]][3]#kcat[1]^(2/3)*kcat[2]^(1/3)#median(kcat)#sqrt(kcat[1]*kcat[2])
          trainedpara2=unlist(sapply(reaclinenum[1:5],function(y){
            paste(y,collapse="     ")
          }))
          #k-2
          reaclinenum[[6]][1]=4
          fixedpara3=paste(reaclinenum[[6]],collapse="     ")
          #km
          sapply(1:length(reaclinenum),function(x){
            reaclinenum[[x]][1]<<-5
          })
          reaclinenum[[1]][3]=km[1]/extend
          reaclinenum[[2]][3]=km[2]*extend
          local.para[["km"]][[name]]<<-km#store the value without extend
          namepara=paste0(name,".km")
          if(namepara%in%names(vals)){
            reaclinenum[[3]][3]=vals[namepara]
          }else if(rand){
            set.seed(rand.seed)
            reaclinenum[[3]][3]=runif(1,min=km[1],max=km[2])
          }else{
            reaclinenum[[3]][3]=sqrt(km[1]*km[2])#km[1]^(2/3)*km[2]^(1/3)#median(km)#sqrt(km[1]*km[2])
          }
          reaclinenum[[4]][3]=reaclinenum[[3]][3]#km[1]^(2/3)*km[2]^(1/3)#sqrt(km[1]*km[2])
          reaclinenum[[5]][3]=reaclinenum[[3]][3]#km[1]^(2/3)*km[2]^(1/3)#sqrt(km[1]*km[2])
          trainedpara1=unlist(sapply(reaclinenum[1:5],function(y){
            paste(y,collapse="     ")
          }))
          tempvecreac=temptabreac[temptabreac[,1]==name,]
          depenc.reac=depend_for(tempvecreac[,1],unlist(tempvecreac[,3:(coln-1)]),tempvecreac[,coln])
          reac.para.addon.len<<-reac.para.addon.len+5
          reactionnum=c(fixedpara1,fixedpara2,trainedpara2,fixedpara3,trainedpara1)
          if(rev[enzyme]){
            flagmmrev=1
            depenc.reac.rev=depend_for(paste0(" ",tempvecreac[,1],"_r"),unlist(tempvecreac[,3:(coln-1)]),tempvecreac[,coln])
          }
        }else{
          trainedpara=unlist(sapply(reaclinenum[1:5],function(y){
            paste(y,collapse="     ")
          }))
          reac.para.addon.len<<-reac.para.addon.len+1
          reactionnum=c(trainedpara,fixedpara)
        }
      }else{
        trainedpara=unlist(sapply(reaclinenum[1:5],function(y){
          paste(y,collapse="     ")
        }))
        fixedpara=paste(reaclinenum[[6]],collapse="     ")
        if((str_detect(string=name,pattern="\\d+\\.\\d+\\.\\d+\\.\\d+")||str_detect(string=name,pattern="^trans\\_"))&&rev[enzyme]){
          reac.para.addon.len<<-reac.para.addon.len+2
          reactionnum=c(trainedpara,trainedpara)
        }else{
          reac.para.addon.len<<-reac.para.addon.len+1
          reactionnum=c(trainedpara,fixedpara)
        }
      }
      string=c(paste(reaclineanno1,collapse=" "),paste0(" ",x[["name"]]),
               paste(" ",lineannonum2,collapse=" "),
               paste(" ",subs),
               paste("           ",prods),
               paste(" ",reaclineanno2,collapse=" "),
               paste(" ",reactionnum,collapse="\n"),"#"
              )
       if(name%in%para.list[["depend"]][["react"]][,"name"]){
         string=c(string,depenc.reac)
       }
       reac.count<<-reac.count+1
       if(flagmmrev){
         lineannonum2=c(length(x[["prods"]]),length(x[["subs"]]),jkin)
         string=c(string,
                  paste(reaclineanno1,collapse=" "),paste0(" ",x[["name"]],"_r"),
                  paste(" ",lineannonum2,collapse=" "),
                  paste(" ",prods),
                  paste("           ",subs),
                  paste(" ",reaclineanno2,collapse=" "),
                  paste(" ",reactionnum,collapse="\n"),"#"
         )
         if(name%in%para.list[["depend"]][["react"]][,"name"]){
           string=c(string,depenc.reac.rev)
         }
         reac.count<<-reac.count+1
       }
       string=string[!str_detect(string=string,pattern="^\\s+$")]
       output=paste(string,collapse="\n")
    })
    list.res=list(out=reac_outputreturn,len=reac.para.addon.len,count=reac.count)
    return(list.res)
}

#' dependences formulation
#'
#' automatic dependence formulation with default format
#'
#' @param y array. the location of the dependent variable
#' @param x array. the independent variable
#' @param func array. the function to be formulated for the dependency
#'
#' @return string. output information
#' @seealso [spec_output()] for producing species block in i01 for default format
#'          [react_output()] for producing reaction block in i01 for default format
depend_for<-function(y,x,func){
  load(paste0(dir.lib,"default.depend.RData"))
  depannonum[2]=length(x)
  dependnum=sapply(seq(length(x)),function(i){
    paste0(c(i,x[i],0.000001,"\n     ",y),collapse=" ")
  })
  string=c(paste(" ",depanno,collapse=" "),
           paste(" ",depannonum,collapse=" "),
           paste(" ",func),
           paste(" ",depannolist,collapse=" "),
           paste("          ",dependnum,collapse="\n"),"#"
  )
  depenc=paste(string,collapse="\n")
  return(depenc)
}
