#Dec 2024 Glen Pridham
#streamlined code for GitHub
#please see https://arxiv.org/abs/2412.07795 for more info and citation info
  #also on GitHub page
library(survival) #needed
library(ggplot2) #needed for plots
#library(metr) #recommended

GenBasisDefault = function(f)
{
  x = matrix(1,nrow=length(f),ncol=2)
  x[,2] = f
  colnames(x)=c("intercept","fi")
  return(x)
}

DiagnoseModel.exp = function(data,
                             fiVar, #fi variable names
                             s, #survival function #should be start-stop format (same size as data)
                             dt,
                             GenBasis=GenBasisDefault, #fi basis functions - generates from fi (function) #shuld return N x Nx matrix or dataframe
                             fix=NULL, #which parameters do you want to fix at default values?
                             fixedValues = NULL,
                             idCol="id",
                             Ntraj = 20, #number of trajectories to plot
                             ageCol = "age",
                             lastMeasurementCol="last_measurement",
                             nboot=100,
                             Ngridf=21, #number of grid points for f-t plane
                             Ngridt=Ngridf,
                             fi_cuts = seq(0,1,by=.1),
                             fi_cuts_surv = c(0,.25,.5,1), #exclusively for survival plots
                             age_cuts = seq(min(data[,ageCol],na.rm=T),max(data[,ageCol]),length=11),
                             Nsim=1000,
                             simStep=median(dt,na.rm=T)/2, #sampling step size
                             simMultiplier=10, #simulation step size is simStep/simMultiplier
                             verbose=FALSE,
                             b0 = rep(0,length(fiVar)),
                             age0 = min(data[,ageCol],na.rm=T), #can be 1 value or vector of length Nsim
                             end_age = max(data[,ageCol],na.rm=T), #must be 1 value
                             par0=NULL,
                             file=NULL,
                             save=F,
                             load=T,
                             skipPlots=F
)
{
  #wrapper function for model selection
  #only fits model of the form:
  #D = exp(basis(f) + basis(f)*t)
  #R = exp(basis(f) + basis(f)*t)
  #h = exp(basis(f) + basis(f)*t)
  #basis variables must provide EVERYTHING, including offsets
  #what it does:
  #1. preprocess data
  #2. estimate initial parameter values
  #3. fit specific model (bootstrap)
  #4. simulate fit parameters
  #5. diagnose model fit
  # -RMSE?
  # -Youden?
  #6. return plots of fit quality
  
  #data:
  #must have columns:
  #fiVar
  #sprintf("%s_next",fiVar)
  #lastMeasurement (before death)
  #age
  #age_next
  
  if(load & !is.null(file))
  {
    if(file.exists(file)) return(readRDS(file))
  }
  
  if(nrow(s)!=nrow(data)) stop("s should be start-stop format and same nrow as data")
  if(!is.null(dim(b0))) b0 = as.matrix(b0) #strip off dataframe if present
  
  Nx = 0
  if(!is.null(GenBasis))
  {
    x = GenBasis(apply(data[,fiVar],1,mean,na.rm=T))
    Nx = ncol(x)
  }
  else stop("Basis functions needed - there are no base parameters in this model")
  
  
  #default values
  if(is.null(fix))  fix=rep(F,6*Nx) #which parameters do you want to fix at default values?
  if(is.null(fixedValues)) fixedValues = rep(0,6*Nx)
  
  
  dr = DamageAndRepair(data[,fiVar],data[,sprintf("%s_next",fiVar)])
  dr[,"ddamagedt"] = dr[,"damage_rate"]/dt
  dr[,"drepairdt"] = dr[,"repair_rate"]/dt
  dr[,"fi"]  = apply(data[,fiVar],1,mean,na.rm=T)
  dr[,"fi_cut"] = as.numeric(as.character(cut(dr[,"fi"],fi_cuts,labels=fi_cuts[-1]/2+fi_cuts[-length(fi_cuts)]/2,include.lowest=T)))
  dr[,"fi_next"]  = apply(data[,sprintf("%s_next",fiVar)],1,mean,na.rm=T)
  
  data[,"fi"] = apply(data[,fiVar],1,mean,na.rm=T)
  data[,"fi_next"] = apply(data[,sprintf("%s_next",fiVar)],1,mean,na.rm=T)
  data[,"age"] = data[,ageCol]
  data[,"age_next"] = data[,sprintf("%s_next",ageCol)]
  data[,"dt"] = data[,"age_next"]-data[,"age"]
  
  #estimate initial parameter values:
  #assume basic exp model to start
  d0 = log(mean(subset(dr,fi < .4)[,"ddamagedt"],na.rm=T))
  r0 = log(mean(subset(dr,fi < .4)[,"drepairdt"],na.rm=T))
  h0 = log(log(2)/quantile(survfit(s~1),probs=.5)[[1]])
  if(is.na(h0)) h0 = log(log(2)/median(data[,ageCol],na.rm=T)) #happends if s never crosses t_1/2
  
  #check for NAs
  d0[is.nan(d0)] = 0
  r0[is.nan(r0)] = 0
  
  if(is.null(par0))
  {
    print(Nx)
    par0 = rep(0,6*Nx)
    par0[1] = d0
    par0[1+2*Nx] = r0
    par0[1+4*Nx] = h0
    
    par0[fix] = fixedValues[fix]
  }
  
  
  if(verbose)
  {
    print("starting guesses:")
    print(par0)
    
    print("starting ll:")
    print(logLik.exp(par0,
                     bn=data[,fiVar],bnp1=data[,sprintf("%s_next",fiVar)],x=x,
                     tn=data[,ageCol],data[,sprintf("%s_next",ageCol)],
                     aliveLogi=rep(T,nrow(data)),lastMeasurementLogi=data[,lastMeasurementCol],fix=fix))
    
    print("starting grad:")
    print(logLik.grad.exp(par0,
                          bn=data[,fiVar],bnp1=data[,sprintf("%s_next",fiVar)],x=x,
                          tn=data[,ageCol],data[,sprintf("%s_next",ageCol)],
                          aliveLogi=rep(T,nrow(data)),lastMeasurementLogi=data[,lastMeasurementCol],fix=fix))
  }
  
  par = list()
  dfdt = list()
  #acc = list() #I'm not sure how to do this without PMF or sim
  success = rep(F,nboot)
  for (i in 1:nboot)
  { 
    cat(".")
    inds = sample(1:nrow(data),replace=T)
    
    xsamp  = x[inds,,drop=F]
    op = tryCatch(optim(par0,fn=logLik.exp,gr=logLik.grad.exp,method="BFGS",hessian=F,control=list(fnscale=-1),
                        bn=data[inds,fiVar],bnp1=data[inds,sprintf("%s_next",fiVar)],x=xsamp,
                        tn=data[inds,ageCol],tnp1=data[inds,sprintf("%s_next",ageCol)],
                        aliveLogi=rep(T,length(inds)),lastMeasurementLogi=data[inds,lastMeasurementCol],fix=fix),
                  error=function(e){return("e")})
    
    
    #if(debug) print(op)
    
    if(all(op=="e"))
    {
      warning("bootstrap failed")
      next
    }
    
    par[[i]] = op$par
    
    dfdt[[i]] = expand.grid(f=seq(0,1,length=Ngridf),age=seq(min(data[,"age"]),max(data[,"age"]),length=Ngridt))
    dfdt[[i]][,"d"] = DamageFun.exp(dfdt[[i]][,"f"],dfdt[[i]][,"age"],par=op$par,x=GenBasis(dfdt[[i]][,"f"]))
    dfdt[[i]][,"r"] = RepairFun.exp(dfdt[[i]][,"f"],dfdt[[i]][,"age"],par=op$par,x=GenBasis(dfdt[[i]][,"f"]))
    dfdt[[i]][,"h"] = SurvivalFun.exp(dfdt[[i]][,"f"],dfdt[[i]][,"age"],par=op$par,x=GenBasis(dfdt[[i]][,"f"]))
    dfdt[[i]][,"S"] = S.exp(dfdt[[i]][,"f"],dfdt[[i]][,"age"],par=op$par,x=GenBasis(dfdt[[i]][,"f"]))
    dfdt[[i]][,"dfdt"] = (1-dfdt[[i]][,"f"])*dfdt[[i]][,"d"]-dfdt[[i]][,"f"]*dfdt[[i]][,"r"]
    
    Ndatum = sum(!is.na(data[inds,fiVar]*data[inds,sprintf("%s_next",fiVar)]))
    #to do: how many datum to include for survival part? #coxph says use number of events
    #eha uses number of data points irrespective of event
    #Volinsky, C. T. & Raftery, A. E. Bayesian information criterion for censored survival models. Biometrics 56, 256–262 (2000)
    #this guy says use number of events
    #I also prefer this for start-stop
    Ndatum = Ndatum + sum(data[inds,lastMeasurementCol],na.rm=T)
    dfdt[[i]][,"ll_train"] = logLik.exp(op$par,
                                        bn=data[inds,fiVar],bnp1=data[inds,sprintf("%s_next",fiVar)],x=xsamp,
                                        tn=data[inds,ageCol],tnp1=data[inds,sprintf("%s_next",ageCol)],
                                        aliveLogi=rep(T,length(inds)),lastMeasurementLogi=data[inds,lastMeasurementCol],fix=fix)
    dfdt[[i]][,"aic_train"] = 2*6*ncol(x) - 2*dfdt[[i]][,"ll_train"]
    dfdt[[i]][,"bic_train"] = 6*ncol(x)*log(Ndatum)-2*dfdt[[i]][,"ll_train"]
    dfdt[[i]][,"ll_ave_train"] = dfdt[[i]][,"ll_train"]/Ndatum
    
    testdata = data[-inds,]
    Ndatumtest = sum(!is.na(testdata[inds,fiVar]*testdata[inds,sprintf("%s_next",fiVar)]))
    Ndatumtest = Ndatumtest + sum(testdata[,lastMeasurementCol],na.rm=T)
    dfdt[[i]][,"ll_test"] = logLik.exp(op$par,
                                       bn=testdata[,fiVar],bnp1=testdata[,sprintf("%s_next",fiVar)],x=x[-inds,,drop=F],
                                       tn=testdata[,ageCol],tnp1=testdata[,sprintf("%s_next",ageCol)],
                                       aliveLogi=rep(T,nrow(testdata)),lastMeasurementLogi=testdata[,lastMeasurementCol],
                                       fix=fix)
    dfdt[[i]][,"aic_test"] = 2*6*ncol(x) - 2*dfdt[[i]][,"ll_test"]
    dfdt[[i]][,"bic_test"] = 6*ncol(x)*log(Ndatumtest)-2*dfdt[[i]][,"ll_test"]
    dfdt[[i]][,"ll_ave_test"] = dfdt[[i]][,"ll_test"]/Ndatumtest
    
    if(any(dfdt[[i]][,"d"] >= Inf))
    {
      warning("Infinite damage")
      dfdt[[i]] = NULL
      next
    }
    if(any(dfdt[[i]][,"r"] >= Inf))
    {
      warning("Infinite repair")
      dfdt[[i]] = NULL
      next
    }
    
    #nullcline estimate
    unt = unique(dfdt[[i]][,"age"])
    for (ii in 1:length(unt))
    {
      logi = dfdt[[i]][,"age"] == unt[ii]
      ind = which.min(abs(dfdt[[i]][logi,"dfdt"]))
      if(length(ind) < 1) dfdt[[i]][logi,"null"] = NA
      else dfdt[[i]][logi,"null"] = dfdt[[i]][ind,"f"]
    }
    
    #error
    #in and out of sample - how????
    #acc[[i]] = data.frame()
    
    success[i] = T
  }
  par = do.call(rbind,par[success])
  
  #newly added, could cause problems
  coef_nms = rep("NA",6*Nx)
  coef_nms[1:Nx] = sprintf("ln_damage_f^%01d",1:Nx-1)
  coef_nms[1:Nx+Nx] = sprintf("ln_damage_f^%01d*t",1:Nx-1)
  coef_nms[1:Nx+2*Nx] = sprintf("ln_repair_f^%01d",1:Nx-1)
  coef_nms[1:Nx+3*Nx] = sprintf("ln_repair_f^%01d*t",1:Nx-1)
  coef_nms[1:Nx+4*Nx] = sprintf("ln_hazard_f^%01d",1:Nx-1)
  coef_nms[1:Nx+5*Nx] = sprintf("ln_hazard_f^%01d*t",1:Nx-1)
  colnames(par) = coef_nms[success]
  
  dfdt = ListMeanSD(dfdt[success],na.rm=T)
  nsuccess = sum(success)
  
  mu = apply(par,2,mean,na.rm=T)
  se = apply(par,2,sd,na.rm=T)
  
  if(verbose)
  {
    print("parameter estimate:")
    print(mu)
  }
  
  
  ################now simulate
  if(verbose) print("Done bootstraps, starting sim...")
  DamageFun = function(f)
  {
    x = GenBasis(f)
    Nx = ncol(x)
    damagePar =  mu[1:Nx]
    lnDamageBasisHazard = x*NA
    
    for (j in 1:Nx)
    {
      lnDamageBasisHazard[,j]    = x[,j]*damagePar [j]
    }
    lnDamageBasisHazard       = apply(lnDamageBasisHazard,1,sum,na.rm=T)
    
    d = exp(lnDamageBasisHazard)
    return(d)
  }
  DamageGompertzFun = function(f)
  {
    x = GenBasis(f)
    Nx = ncol(x)
    damageParT = mu[1:Nx+  Nx]
    damageGompertz = x*NA
    
    for (j in 1:Nx)
    {
      damageGompertz[,j]         = x[,j]*damageParT[j]
    }
    damageGompertz            = apply(damageGompertz,1,sum,na.rm=T)
    
    return(damageGompertz)
  }
  
  RepairFun = function(f)
  {
    x = GenBasis(f)
    Nx = ncol(x)
    repairPar =  mu[1:Nx+2*Nx]
    
    lnRepairBasisHazard = x*NA
    
    for (j in 1:Nx)
    {
      lnRepairBasisHazard[,j]    = x[,j]*repairPar [j]
      
    }
    lnRepairBasisHazard       = apply(lnRepairBasisHazard,1,sum,na.rm=T)
    
    r = exp(lnRepairBasisHazard)
    return(r)
  }
  
  RepairGompertzFun = function(f)
  {
    x = GenBasis(f)
    Nx = ncol(x)
    repairParT = mu[1:Nx+3*Nx]
    repairGompertz = x*NA
    
    for (j in 1:Nx)
    {
      repairGompertz[,j]         = x[,j]*repairParT[j]
      
    }
    repairGompertz            = apply(repairGompertz,1,sum,na.rm=T)
    
    return(repairGompertz)
  }
  SurvivalFun = function(f)
  {
    x = GenBasis(f)
    Nx = ncol(x)
    hazardPar =  mu[1:Nx+4*Nx]
    lnSurvivalBasisHazard = x*NA
    
    for (j in 1:Nx)
    {
      lnSurvivalBasisHazard[,j]  = x[,j]*hazardPar [j]
    }
    lnSurvivalBasisHazard     = apply(lnSurvivalBasisHazard,1,sum,na.rm=T)
    
    h = exp(lnSurvivalBasisHazard)
    return(h)
  }
  
  SurvivalGompertzFun = function(f)
  {
    x = GenBasis(f)
    Nx = ncol(x)
    hazardParT = mu[1:Nx+5*Nx]
    survivalGompertz = x*NA
    
    for (j in 1:Nx)
    {
      survivalGompertz[,j]       = x[,j]*hazardParT[j]
    }
    survivalGompertz          = apply(survivalGompertz,1,sum,na.rm=T)
    
    return(survivalGompertz)
  }
  
  
  if(is.null(dim(b0))) p = length(b0)
  else p = ncol(b0)
  sim = QuickSimExp(N=Nsim,twindow=c(0,end_age-min(age0,na.rm=T)),t0=age0,
                    dt=simStep,innerMultiplier = simMultiplier,
                    DamageFun=DamageFun,RepairFun=RepairFun, SurvivalFun=SurvivalFun,
                    DamageGompertzFun=DamageGompertzFun,RepairGompertzFun=RepairGompertzFun,
                    SurvivalGompertzFun=SurvivalGompertzFun,b0=b0,p=p) 
  #sim0 = QuickSimExp(N=Nsim,twindow=c(min(age0,na.rm=T),end_age),t0=0, #for debug
  #                dt=simStep,innerMultiplier = simMultiplier,
  #                DamageFun=DamageFun,RepairFun=RepairFun, SurvivalFun=SurvivalFun,
  #                DamageGompertzFun=DamageGompertzFun,RepairGompertzFun=RepairGompertzFun,
  #                SurvivalGompertzFun=SurvivalGompertzFun,b0=b0,p=p) 
  sim0 = NULL
  
  if(skipPlots)
  {
    return(list(par=par,sim=sim,mu=mu,sd=sd))
  }
  

  print("Done simulations")
  
  #mean-field FI
  t_mft = seq(mean(age0,na.rm=T),end_age,length=1001)
  fi_mft = EulerFI(t=t_mft,DamageFun=DamageFun,RepairFun=RepairFun,damageGompertz=mu[3],repairGompertz=mu[6],f0=mean(b0,na.rm=T))
  fi_mft = data.frame(t=t_mft,fi=fi_mft)
  
  ###############plots
  plots = DiagPlots(data=data,dr=dr,s=s,sim=sim,dfdt=dfdt,fi_cuts=fi_cuts,fi_cuts_surv=fi_cuts_surv,
                    Ntraj=Ntraj,
                    age_cuts=age_cuts,idCol=idCol,ageCol=ageCol,fix=fix,fixedValues=fixedValues,fiVar=fiVar)
  
  l = list(par=par,se=se,dfdt=dfdt,plots=plots,sim=sim,dr=dr,par0=par0,
           nboot=nboot,Ntraj=Ntraj,Nsim=Nsim,
           DamageFun=DamageFun,RepairFun=RepairFun,SurvivalFun=SurvivalFun,
           DamageGompertzFun=DamageGompertzFun,RepairGompertzFun=RepairGompertzFun,
           SurvivalGompertzFun=SurvivalGompertzFun,GenBasis=GenBasis,
           b0=b0,age0=age0,
           fix=fix,fixedValues=fixedValues,data=data,
           fi_mft=fi_mft,sim0=sim0,nsuccess=nsuccess)
  #agg=agg,fitraj=fitraj,fidata=fidata,
  if(save)
  {
    saveRDS(l,file)
  }
  
  return(l)
}


DiagPlots = function(data,
                     dr,
                     s,
                     sim,
                     dfdt,
                     fi_cuts,
                     fi_cuts_surv,
                     age_cuts,
                     idCol,
                     ageCol,
                     Ntraj,
                     fix,
                     fixedValues,
                     fiVar
)
{
  #to do:
  #add damage and repair rates from sim
  #add PMF
  #auto-regressive stuff to confirm paths loko right
  #keep in mind the time lag may matter, so might be best to fit SF-like model and cmpute those
  #add spaghetti plots with quantile colours
  #can later add transition stuff
  
  library(colorspace)
  
  #survival plot
  plots = list()
  #basic survival plot:
  sf = survfit(s~1)
  sf_sim = survfit(sim$s~1)
  
  sort_age = sort.list(dfdt[[1]][,"age"])
  survplotdata = rbind(data.frame(t=sf$time,f=0,s=sf$surv,slow=sf$lower,shigh=sf$upper,label="GT"),
                       data.frame(t=sf_sim$time,f=0,s=sf_sim$surv,slow=sf_sim$lower,shigh=sf_sim$upper,label="Sim")
                       # data.frame(t=dfdt[[1]][sort_age,"age"],f=dfdt[[2]][sort_age,"f"],s=dfdt[[1]][sort_age,"S"], #problem: need to aggregate by f - or stratify
                       #            slow=dfdt[[1]][sort_age,"S"]-1.96*dfdt[[2]][sort_age,"S"],
                       #           shigh=dfdt[[1]][sort_age,"S"]+1.96*dfdt[[2]][sort_age,"S"],label="Fit")
  )
  
  plots[["s0"]] = ggplot(survplotdata,aes(x=t,y=s,ymin=slow,ymax=shigh,colour=label,fill=label))+
    geom_point(data=subset(survplotdata,label=="GT"))+
    geom_line()+
    geom_ribbon(alpha=.15,colour=NA)+
    theme_minimal()
  if(fix[8] & abs(fixedValues[8]) < 1e-6) #no FI term in survival
  {
    plots[["S"]] = NULL
  } else
  {
    #GT
    fi_cut = cut(dr[,"fi"],fi_cuts_surv,labels=fi_cuts_surv[-1]/2+fi_cuts_surv[-length(fi_cuts_surv)]/2,include.lowest=T)
    sf = survfit(s~fi_cut)
    survplotdata = list()
    indShift = 0
    if(length(sf$strata) > 0)
    {
      for (j in 1:length(sf$strata))
      {
        inds = 1:sf$strata[j]+indShift
        nm = strsplit(names(sf$strata)[j],"fi_cut=")[[1]][2]
        survplotdata[[j]] = data.frame(t=sf$time[inds],f=nm,s=sf$surv[inds],slow=sf$lower[inds],shigh=sf$upper[inds],label="GT")
        indShift = indShift+sf$strata[j]
      }
    }
    
    shift = length(survplotdata)
    
    #Sim
    #simfi = apply(sim$wave,3,mean,na.rm=T)
    fi_cut = cut(sim$stst[,"fi"],fi_cuts_surv,labels=fi_cuts_surv[-1]/2+fi_cuts_surv[-length(fi_cuts_surv)]/2,include.lowest=T)
    sf = survfit(sim$sstst~fi_cut)
    indShift = 0
    if(length(sf$strata) > 0)
    {
      for (j in 1:length(sf$strata))
      {
        inds = 1:sf$strata[j]+indShift
        nm = strsplit(names(sf$strata)[j],"fi_cut=")[[1]][2]
        survplotdata[[j+shift]] = data.frame(t=sf$time[inds],f=nm,s=sf$surv[inds],slow=sf$lower[inds],shigh=sf$upper[inds],label="Sim")
        indShift = indShift+sf$strata[j]
      }
    }
    
    #fc = cut(dfdt[[2]][,"f"],fi_cuts,labels=fi_cuts[-1]/2+fi_cuts[-length(fi_cuts)]/2,include.lowest=T)
    #survplotdata[[length(survplotdata)+1]] = data.frame(t=dfdt[[1]][,"age"],f=fc,s=dfdt[[1]][,"S"],
    #                                                    slow=dfdt[[1]][,"S"]-dfdt[[2]][,"S"],
    #                                                    shigh=dfdt[[1]][,"S"]+dfdt[[2]][,"S"],label="Fit")
    survplotdata=do.call(rbind,survplotdata)
    plots[["S"]] = ggplot(survplotdata,aes(x=t,y=s,ymin=slow,ymax=shigh,colour=f,fill=f,lty=label))+
      geom_point(data=subset(survplotdata,label=="GT"))+
      geom_line()+
      geom_ribbon(alpha=.15,colour=NA)+
      scale_colour_discrete_diverging("Blue-Red")+
      scale_fill_discrete_diverging("Blue-Red")+
      #scale_colour_brewer(palette="RdBu")+
      #scale_fill_brewer(palette="RdBu")+
      #geom_point()+
      theme_minimal()
  }
  
  
  dfdt[[1]][,"dfdt_se"] = dfdt[[2]][,"dfdt"]
  dfdt[[1]][,"null_se"] = dfdt[[2]][,"null"]
  dfdt[[1]][,"null_low"] = dfdt[[1]][,"null"] - dfdt[[2]][,"null"]
  dfdt[[1]][,"null_high"] = dfdt[[1]][,"null"] + dfdt[[2]][,"null"]
  logi = dfdt[[1]][,"null_high"] > 1
  logi[is.na(logi)] = F
  dfdt[[1]][logi,"null_high"] = 1
  
  #vector plot for df/dt
  if (!requireNamespace("metR", quietly = TRUE)) {
    warning(sprintf("metR library not found, skipping nullcline plots"))
    plots[["dfdt"]] = NULL
  }
  else
  {
    library(metR)
    plots[["dfdt"]] = ggplot(dfdt[[1]],aes(x=age,y=f,z=dfdt,colour=dfdt,fill=dfdt))+labs(title="df/dt")+
      geom_vector(aes(dx=1/sqrt(2), dy=dfdt/sqrt(2)))+
      geom_vector(aes(dx=1/sqrt(2), dy=(dfdt-dfdt_se)/sqrt(2)),arrow.length=0)+
      geom_vector(aes(dx=1/sqrt(2), dy=(dfdt+dfdt_se)/sqrt(2)),arrow.length=0)+
      geom_line(aes(x=age,y=null),colour="red",size=1)+
      geom_ribbon(aes(x=age,y=null,ymin=null_low,ymax=null_high),colour=NA,fill="red",alpha=.2)+
      scale_colour_gradient2(mid="gray75",low="red",high="blue")+
      scale_y_continuous(limits=c(0,1))+
      theme_classic()
    
    #include points and/or density
    #plots[["dfdt2"]] = plots[["dfdt"]]+geom_point(data=data,aes(x=age,y=fi),inherit.aes=F) #prob too big for HRDS
    plots[["dfdt2"]] = plots[["dfdt"]]+geom_density2d(data=data,aes(x=age,y=fi),inherit.aes=F)
    plots[["dfdt3"]] = plots[["dfdt"]]+geom_density2d(data=data,aes(x=age,y=fi),inherit.aes=F,breaks=10^seq(-4,0,length=10)) #bins by densities seems to be < 1
    plots[["dfdt4"]] = plots[["dfdt"]]+geom_density_2d_filled(data=data,aes(x=age,y=fi),inherit.aes=F,breaks=10^seq(-4,0,length=10)) #fails?

  }

  #damage and repair contours
  plots[["d"]] = ggplot(dfdt[[1]],aes(x=age,y=f,z=d))+geom_contour_filled()+labs(x="Age",y="FI",title="Damage")+theme_minimal()
  plots[["r"]] = ggplot(dfdt[[1]],aes(x=age,y=f,z=r))+geom_contour_filled()+labs(x="Age",y="FI",title="Repair")+theme_minimal()
  plots[["h"]] = ggplot(dfdt[[1]],aes(x=age,y=f,z=h))+geom_contour_filled()+labs(x="Age",y="FI",title="Survival (hazard)")+theme_minimal()
  
  plots[["logd"]] = ggplot(dfdt[[1]],aes(x=age,y=f,z=log(d,10)))+geom_contour_filled()+labs(x="Age",y="FI",title=bquote("log"[10]*"(Damage)"))+theme_minimal()
  plots[["logr"]] = ggplot(dfdt[[1]],aes(x=age,y=f,z=log(r,10)))+geom_contour_filled()+labs(x="Age",y="FI",title=bquote("log"[10]*"(Repair)"))+theme_minimal()
  plots[["logh"]] = ggplot(dfdt[[1]],aes(x=age,y=f,z=log(h,10)))+geom_contour_filled()+labs(x="Age",y="FI",title=bquote("log"[10]*"(hazard)"))+theme_minimal()
  
  
  #marginal damage and repair plots
  #could stratify by age/fi...
  dam = aggregate(dr[,"ddamagedt",drop=F],by=list(fi=dr[,"fi_cut"]),mean,na.rm=T)
  dam[,"damage_se"] = aggregate(dr[,"ddamagedt",drop=F],by=list(fi=dr[,"fi_cut"]),SEM,na.rm=T)[,"ddamagedt"]
  colnames(dam) = c("fi","rate","se")
  
  funData = data.frame(fi=seq(0,1,length=101))
  funData[,"rate"] = sim$DamageFun(funData[,"fi"])
  funData[,"se"] = NA
  
  plots[["dfi"]] = ggplot(dam,aes(x=fi,y=rate,ymin=rate-se,ymax=rate+se))+
    geom_pointrange()+
    geom_smooth()+
    geom_line(data=funData)+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_classic()
  
  repa = aggregate(dr[,"drepairdt",drop=F],by=list(fi=dr[,"fi_cut"]),mean,na.rm=T)
  repa[,"repair_se"] = aggregate(dr[,"drepairdt",drop=F],by=list(fi=dr[,"fi_cut"]),SEM,na.rm=T)[,"drepairdt"]
  colnames(repa) = c("fi","rate","se")
  
  funData = data.frame(fi=seq(0,1,length=101))
  funData[,"rate"] = sim$RepairFun(funData[,"fi"])
  funData[,"se"] = NA
  
  plots[["rfi"]] = ggplot(repa,aes(x=fi,y=rate,ymin=rate-se,ymax=rate+se))+
    geom_pointrange()+
    geom_smooth()+
    geom_line(data=funData)+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_classic()
  
  
  fidata = list()
  fidata[[1]] = data[,c(idCol,ageCol,"fi")]
  colnames(fidata[[1]])[1:2]=c("id","age")
  fidata[[1]][,"label"] = "GT"
  fi = apply(sim$wave,c(1,3),mean,na.rm=T)
  fidata[[2]] = data.frame(id=c(outer(1:nrow(fi),rep(1,ncol(fi)))),age=c(outer(sim$t0,sim$t,FUN="+")),fi=c(fi), label="Sim")
  #fi_se =apply(fi,2,SEM,na.rm=T),fi_sd =apply(fi,2,sd,na.rm=T),
  fi0 = apply(sim$wave0,c(1,3),mean,na.rm=T)
  fidata[[3]] = data.frame(id=c(outer(1:nrow(fi),rep(1,ncol(fi)))),age=c(outer(sim$t0,sim$t,FUN="+")),fi=c(fi0), 
                           label="Sim - w zombies")
  
  fidata = do.call(rbind,fidata)
  fidata[,"age_cut"] = as.numeric(as.character(cut(fidata[,"age"],age_cuts,labels=age_cuts[-1]/2+age_cuts[-length(age_cuts)]/2,include.lowest=T)))
  
  
  
  plots[["fi"]] = ggplot(fidata,aes(x=age,y=fi,colour=label,fill=label,lty=label))+
    #geom_point(size=.1)+
    geom_smooth(size=.5)+ #se=FALSE,
    stat_summary(mapping=aes(x=age_cut,shape=label),position=position_dodge(1))+
    labs(x="Age",y="FI")+
    theme_minimal()
  
  
  agg = aggregate(fidata[,c("age","fi")],by=list(age_cut=fidata[,"age_cut"],label=fidata[,"label"]),mean,na.rm=T)
  agg[,"fi_sd"] = aggregate(fidata[,c("age","fi")],by=list(age_cut=fidata[,"age_cut"],label=fidata[,"label"]),sd,na.rm=T)[,"fi"]
  plots[["cv"]] =  ggplot(agg,aes(x=age,y=fi_sd/fi,colour=label,fill=label,lty=label))+
    geom_smooth()+
    labs(x="Age",y=expression("coefficient of variation, "~frac(sigma,mu)))+
    theme_minimal()
  
  samp = sample(unique(data[,idCol]),min(c(unlen(data[,idCol]),Ntraj)))
  logi = (fidata[,"id"] %in% samp) & fidata[,"label"] == "GT"
  samp = sample(1:nrow(fi),Ntraj)
  #print(samp)
  logi = logi | ( (fidata[,"id"] %in% samp) & fidata[,"label"] == "Sim" )
  #logi = logi | (fidata[,"id"] %in% sample(1:nrow(fi),Ntraj)) & fidata[,"label"] == "Sim - w zombies"
  logi[is.na(logi)] = F
  fitraj = fidata[logi,]
  plots[["fitraj"]] = ggplot(fitraj,aes(x=age,y=fi,colour=label,fill=label,group=id))+
    #geom_smooth()+
    geom_line(alpha=.5)+
    theme_minimal()
  
  plots[["box"]] = ggplot(subset(fidata,label%in%c("GT","Sim")),aes(x=ordered(round(age_cut,1)),y=fi,colour=label))+
    geom_boxplot()+labs(x="Age",y="FI")+theme_minimal()
  plots[["violin"]] = ggplot(subset(fidata,label%in%c("GT","Sim")),aes(x=ordered(round(age_cut,1)),y=fi,colour=label,fill=label))+
    geom_violin()+labs(x="Age",y="FI")+theme_minimal()
  
  scale = 0.95
  subdata = subset(fidata,label%in%c("GT","Sim")) #make a dataset with same available groups
  subdata[,"age_cut_ord"] = ordered(round(subdata[,"age_cut"],1))
  #exclude any cuts that aren't shared by both GT and Sim
  aggCheck = aggregate(subdata[,"label"],by=list(age=subdata[,"age_cut_ord"]),unlen)
  allowAges = aggCheck[aggCheck[,"x"] == 2,"age"]
  subdata = subdata[subdata[,"age_cut_ord"]%in%allowAges,]
  
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    warning(sprintf("ggridges library not found, skipping nullcline plots"))
  }
  else
  {
    library(ggridges)
  plots[["ridge_gt"]] = ggplot(subset(subdata,label%in%c("GT")),aes(x=fi,y=ordered(round(age_cut,1)),fill=factor(stat(quantile))))+
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,quantiles = 4, quantile_lines = TRUE,scale=scale)+
    scale_fill_viridis_d(name = "Quartiles")+
    geom_density_ridges(inherit.aes=F,aes(x=fi,y=ordered(round(age_cut,1))),stat = "binline", scale = scale, 
                        draw_baseline = FALSE,alpha=.33)+
    scale_x_continuous(limits=c(0,1))+
    labs(y="Age",x="FI",title="GT")+
    theme_minimal()
  plots[["ridge_sim"]] = ggplot(subset(subdata,label%in%c("Sim")),aes(x=fi,y=ordered(round(age_cut,1)),fill=factor(stat(quantile))))+
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,quantiles = 4, quantile_lines = TRUE,scale=scale)+
    scale_fill_viridis_d(name = "Quartiles")+
    geom_density_ridges(inherit.aes=F,aes(x=fi,y=ordered(round(age_cut,1))),stat = "binline", scale = scale, 
                        draw_baseline = FALSE,alpha=.33)+
    scale_x_continuous(limits=c(0,1))+
    labs(y="Age",x="FI",title="Sim")+
    theme_minimal()
  plots[["ridge"]] = ggplot(subdata,aes(x=fi,y=ordered(round(age_cut,1)),fill=label,colour=label))+
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,quantiles = 4, quantile_lines = TRUE,scale=scale,fill=NA)+
    geom_density_ridges(inherit.aes=F,aes(x=fi,y=ordered(round(age_cut,1)),fill=label),stat = "binline", scale = scale, draw_baseline = FALSE,alpha=.15)+
    scale_x_continuous(limits=c(0,.5))+
    labs(y="Age",x="FI")+
    theme_minimal()
  
  plots[["ridge2"]] = ggplot(subset(fidata,label%in%c("GT","Sim")),aes(x=fi,y=ordered(round(age_cut,1)),fill=label,colour=label))+
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,quantiles = 4, quantile_lines = TRUE,scale=scale,fill=NA)+
    geom_density_ridges(inherit.aes=F,aes(x=fi,y=ordered(round(age_cut,1)),fill=label),stat = "binline", scale = scale, draw_baseline = FALSE,alpha=.15)+
    labs(y="Age",x="FI")+
    theme_minimal()
  }
  
  #hazard ratios parameter plot
  #marginal
  
  plotdata = dfdt[[1]]
  plotdata[,"dse"] = dfdt[[2]][,"d"]
  plotdata[,"rse"] = dfdt[[2]][,"r"]
  plotdata[,"hse"] = dfdt[[2]][,"h"]
  sub = plotdata
  sub[,"age"] = round(sub[,"age"],2)
  unage = unique(sub[,"age"])
  unage = unage[!is.na(unage)]
  unage = unage[round(seq(1,length(unage),length=5))]
  sub = subset(sub,age %in% unage)
  plots[["dam_fi"]] = ggplot(sub,aes(x=f,y=d,ymin=d-dse,ymax=d+dse,colour=ordered(round(age,2)),fill=ordered(round(age,2))))+
    geom_pointrange()+
    #geom_ribbon(alpha=.15,colour=NA)+
    scale_colour_discrete_diverging("Blue-Red")+
    scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_minimal()
  plots[["rep_fi"]] = ggplot(sub,aes(x=f,y=r,ymin=r-rse,ymax=r+rse,colour=ordered(round(age,2)),fill=ordered(round(age,2))))+
    geom_pointrange()+
    #geom_ribbon(alpha=.15,colour=NA)+
    scale_colour_discrete_diverging("Blue-Red")+
    scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_minimal()
  plots[["haz_fi"]] = ggplot(sub,aes(x=f,y=h,ymin=h-hse,ymax=h+hse,colour=ordered(round(age,2)),fill=ordered(round(age,2))))+
    geom_pointrange()+
    #geom_ribbon(alpha=.15,colour=NA)+
    scale_colour_discrete_diverging("Blue-Red")+
    scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_minimal()
  
  sub = plotdata
  sub[,"f"] = round(sub[,"f"],2)
  unf = unique(sub[,"f"])
  unf = unf[!is.na(unf)]
  unf = unf[round(seq(1,length(unf),length=5))]
  sub = subset(sub,f %in% unf)
  plots[["dam_age"]] = ggplot(sub,aes(x=age,y=d,ymin=d-dse,ymax=d+dse,colour=ordered(round(f,2)),fill=ordered(round(f,2))))+
    geom_pointrange()+
    #geom_ribbon(alpha=.15,colour=NA)+
    scale_colour_discrete_diverging("Blue-Red")+
    scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_minimal()
  plots[["rep_age"]] = ggplot(sub,aes(x=age,y=r,ymin=r-rse,ymax=r+rse,colour=ordered(round(f,2)),fill=ordered(round(f,2))))+
    geom_pointrange()+
    #geom_ribbon(alpha=.15,colour=NA)+
    scale_colour_discrete_diverging("Blue-Red")+
    scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_minimal()
  plots[["haz_age"]] = ggplot(sub,aes(x=age,y=h,ymin=h-hse,ymax=h+hse,colour=ordered(round(f,2)),fill=ordered(round(f,2))))+
    geom_pointrange()+
    #geom_ribbon(alpha=.15,colour=NA)+
    scale_colour_discrete_diverging("Blue-Red")+
    scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_minimal()
  
  #combined
  hazardData = list()
  hazardData[[1]] = plotdata[,c("f","age","d","dse")]
  colnames(hazardData[[1]])=c("f","age","rate","se")
  hazardData[[1]][,"type"] = "damage"
  hazardData[[2]] = plotdata[,c("f","age","r","rse")]
  colnames(hazardData[[2]])=c("f","age","rate","se")
  hazardData[[2]][,"type"] = "repair"
  hazardData[[3]] = plotdata[,c("f","age","h","hse")]
  colnames(hazardData[[3]])=c("f","age","rate","se")
  hazardData[[3]][,"type"] = "survival"
  hazardData = do.call(rbind,hazardData)
  hazardData[,"age"] = round(hazardData[,"age"],2)
  hazardData[,"f"] = round(hazardData[,"f"],2)
  
  
  logi = hazardData[,"age"]==median(hazardData[,"age"])
  plots[["rates_fi"]] = ggplot(hazardData[logi,],aes(x=f,y=rate,ymin=rate-se,ymax=rate+se,colour=type,fill=type))+
    geom_line(size=1)+
    geom_ribbon(alpha=.15,colour=NA)+
    #scale_colour_discrete_diverging("Blue-Red")+
    #scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_minimal()
  logi = hazardData[,"f"]==median(hazardData[,"f"])
  plots[["rates_age"]] = ggplot(hazardData[logi,],aes(x=age,y=rate,ymin=rate-se,ymax=rate+se,colour=type,fill=type))+
    geom_line(size=1)+
    geom_ribbon(alpha=.15,colour=NA)+
    #scale_colour_discrete_diverging("Blue-Red")+
    #scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    theme_minimal()
  
  #damage and repair transitions
  tr = WavesToTransitions(current_wave=data[,fiVar],next_wave=data[,sprintf("%s_next",fiVar)],
                          age=data[,ageCol],age_next=data[,sprintf("%s_next",ageCol)])
  simtr = WavesToTransitions(wave=sim[["wave"]],age=c(sim[["age"]][,-ncol(sim[["age"]])]),age_next=c(sim[["age"]][,-1]))
  
  #damage 'survival'
  sf = survfit(tr[["Sd"]]~1)
  sf_sim = survfit(simtr[["Sd"]]~1)
  
  Sdplotdata = rbind(data.frame(t=sf$time,f=0,s=sf$surv,slow=sf$lower,shigh=sf$upper,label="GT"),
                     data.frame(t=sf_sim$time,f=0,s=sf_sim$surv,slow=sf_sim$lower,shigh=sf_sim$upper,label="Sim")
  )
  
  plots[["Sd"]] = ggplot(Sdplotdata,aes(x=t,y=s,ymin=slow,ymax=shigh,colour=label,fill=label))+
    geom_point(data=subset(Sdplotdata,label=="GT"))+
    geom_line()+
    geom_ribbon(alpha=.15,colour=NA)+
    labs(x="Age",y="Sd")+
    theme_minimal()
  
  sfmodel = "sfd"
  plotdata = data.frame(t=simtr[[sfmodel]][["time"]],s=simtr[[sfmodel]][["surv"]],slow=simtr[[sfmodel]][["lower"]],shigh=simtr[[sfmodel]][["upper"]])
  plotdata[,"fi"] = NA
  shift = 0
  for (i in 1:length(simtr[[sfmodel]][["strata"]]))
  {
    inds = 1:simtr[[sfmodel]][["strata"]][i]+shift
    plotdata[inds,"fi"] = names(simtr[[sfmodel]][["strata"]])[i]
    shift = max(inds)
  }
  
  simdata = list()
  levs = levels(simtr[["fi"]])
  for (i in 1:length(levs))
  {
    simdata[[i]] = data.frame(t=seq(0,100,length=101),fi=as.numeric(as.character(levs[i])))
    simdata[[i]][,"s"] = exp(-sim$DamageFun(simtr[["fi_cuts"]][i])*simdata[[i]][,"t"])
  }
  simdata = do.call(rbind,simdata)
  simdata[,"fi"] = sprintf("fid=%.02f",simdata[,"fi"])
  simdata[,"slow"] = simdata[,"s"]
  simdata[,"shigh"] = simdata[,"s"]
  
  
  plots[["Sdfi"]] = ggplot(plotdata,aes(x=t,y=s,ymin=slow,ymax=shigh,colour=fi,fill=fi))+
    geom_line()+
    geom_ribbon(alpha=.15,colour=NA)+
    geom_line(data=simdata,lty=3,size=1)+
    #geom_pointrange()+
    scale_colour_discrete_diverging("Blue-Red")+
    scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    labs(x="Age",y="Sd")
  
  #repair 'survival'
  sf = survfit(tr[["Sr"]]~1)
  sf_sim = survfit(simtr[["Sr"]]~1)
  
  Srplotdata = rbind(data.frame(t=sf$time,f=0,s=sf$surv,slow=sf$lower,shigh=sf$upper,label="GT"),
                     data.frame(t=sf_sim$time,f=0,s=sf_sim$surv,slow=sf_sim$lower,shigh=sf_sim$upper,label="Sim")
  )
  
  plots[["Sr"]] = ggplot(Srplotdata,aes(x=t,y=s,ymin=slow,ymax=shigh,colour=label,fill=label))+
    geom_point(data=subset(Srplotdata,label=="GT"))+
    geom_line()+
    geom_ribbon(alpha=.15,colour=NA)+
    labs(x="Age",y="Sr")+
    theme_minimal()
  
  sfmodel = "sfr"
  plotdata = data.frame(t=simtr[[sfmodel]][["time"]],s=simtr[[sfmodel]][["surv"]],slow=simtr[[sfmodel]][["lower"]],shigh=simtr[[sfmodel]][["upper"]])
  plotdata[,"fi"] = NA
  shift = 0
  for (i in 1:length(simtr[[sfmodel]][["strata"]]))
  {
    inds = 1:simtr[[sfmodel]][["strata"]][i]+shift
    plotdata[inds,"fi"] = names(simtr[[sfmodel]][["strata"]])[i]
    shift = max(inds)
  }
  
  simdata = list()
  levs = levels(simtr[["fi"]])
  for (i in 1:length(levs))
  {
    simdata[[i]] = data.frame(t=seq(0,100,length=101),fi=as.numeric(as.character(levs[i])))
    simdata[[i]][,"s"] = exp(-sim$DamageFun(simtr[["fi_cuts"]][i])*simdata[[i]][,"t"])
  }
  simdata = do.call(rbind,simdata)
  simdata[,"fi"] = sprintf("fid=%.02f",simdata[,"fi"])
  simdata[,"slow"] = simdata[,"s"]
  simdata[,"shigh"] = simdata[,"s"]
  
  
  plots[["Srfi"]] = ggplot(plotdata,aes(x=t,y=s,ymin=slow,ymax=shigh,colour=fi,fill=fi))+
    geom_line()+
    geom_ribbon(alpha=.15,colour=NA)+
    geom_line(data=simdata,lty=3,size=1)+
    #geom_pointrange()+
    scale_colour_discrete_diverging("Blue-Red")+
    scale_fill_discrete_diverging("Blue-Red")+
    scale_y_log10()+
    annotation_logticks(sides="l")+
    labs(x="Age",y="Sr")
  
  #auto-regressive stuff
  #SF model
  data[,"intercept"] = 1
  #print(FitSFModel(yn=data[,"fi",drop=F],ynp1=data[,"fi_next",drop=F],xn=data[,c("intercept","age")],dt=data[,"dt"],
  #                 options=list(diagonalW=T,error=T))) #debug
  sfm = tryCatch(FitSFModel(yn=data[,"fi",drop=F],ynp1=data[,"fi_next",drop=F],xn=data[,c("intercept","age")],dt=data[,"dt"],
                            options=list(diagonalW=T,error=T)),error=function(e){return(NULL)})
  if(all(is.null(sfm))) #try fitting again
  {
    warning("sf fit failed for gt, trying again...")
    sfm = tryCatch(FitSFModel(yn=data[,"fi",drop=F],ynp1=data[,"fi_next",drop=F],xn=data[,c("intercept","age")],dt=data[,"dt"],
                              options=list(diagonalW=T,error=F)),error=function(e){return(NULL)})
  }
  
  simdata = list()
  for (k in 2:dim(sim$wave)[3])
  {
    simdata[[k]] = data.frame(fi=apply(sim$wave[,,k-1],1,mean,na.rm=T),fi_next=apply(sim$wave[,,k],1,mean,na.rm=T),
                              id=1:nrow(sim$wave),
                              age=sim$age[,k-1],age_next=sim$age[,k],intercept=1,
                              dt=sim$age[,k]-sim$age[,k-1]
    )
  }
  simdata = do.call(rbind,simdata)
  sfmsim = tryCatch(FitSFModel(yn=simdata[,"fi",drop=F],ynp1=simdata[,"fi_next",drop=F],xn=simdata[,c("intercept","age")],dt=simdata[,"dt"],
                               options=list(diagonalW=T,error=T)),error=function(e){return(NULL)})
  
  
  if(all(is.null(sfm)) | all(is.null(sfmsim))) #fit failed
  {
    warning("sf fit failed")
    sftab = NULL
    sfMean = NULL
    plots[["sf_fi"]] = NULL
    plots[["sf_par"]] = NULL
  }
  else
  {
    #how to visualize?
    #1. pinhead plot
    #problem: each parameter has different scale
    #could to W - W_est
    #2. mean plot - predicted mean for each vs real mean vs age
    #include mean and model mean for both GT and sim so I can simultaneously compare SF model and Sim vs GT
    #3. table of  values and z score comparisons
    gtpar = c(sfm$W[1,1],as.numeric(sfm$lambda[1,]),sqrt(sfm$sigma[1,1]),1/sfm$W[1,1])
    dgtpar = c(sfm$dW[1,1],as.numeric(sfm$dlambda[1,]),NA,sfm$dW[1,1]/sfm$W[1,1]^2)
    simpar = c(sfmsim$W[1,1],as.numeric(sfmsim$lambda[1,]),sqrt(sfmsim$sigma[1,1]),1/sfmsim$W[1,1])
    dsimpar = c(sfmsim$dW[1,1],as.numeric(sfmsim$dlambda[1,]),NA,sfmsim$dW[1,1]/sfmsim$W[1,1]^2)
    sftab = data.frame(GT=gtpar,GTSE=dgtpar, Sim=simpar,SimSE=dsimpar, z=abs(gtpar-simpar)/sqrt(dgtpar^2+dsimpar^2))
    rownames(sftab)=c("W","mu0","mut","sigma","tau")
    sftab[,"p"] = 2*(1-pnorm(sftab[,"z"]))
    
    ttest = seq(min(data[,"age"],na.rm=T),max(data[,"age"],na.rm=T),length=101)
    
    fi0 = mean(sim$wave[,,1],na.rm=T)
    figt = (fi0-sfm$lambda[1,"intercept"]-ttest[1]*sfm$lambda[1,"age"])*exp(sfm$W[1,1]*(ttest-ttest[1]))+
      sfm$lambda[1,"age"]/sfm$W[1,1]*(1-exp(sfm$W[1,1]*(ttest-ttest[1])))+
      sfm$lambda[1,"intercept"]+ttest*sfm$lambda[1,"age"]
    fisim = (fi0-sfmsim$lambda[1,"intercept"]-ttest[1]*sfmsim$lambda[1,"age"])*exp(sfmsim$W[1,1]*(ttest-ttest[1]))+
      sfmsim$lambda[1,"age"]/sfmsim$W[1,1]*(1-exp(sfmsim$W[1,1]*(ttest-ttest[1])))+
      sfmsim$lambda[1,"intercept"]+ttest*sfmsim$lambda[1,"age"]
    
    sfMean = rbind(data.frame(age=ttest,fi=figt,label="GT"),
                   data.frame(age=ttest,fi=fisim,label="Sim")
    )
    
    plots[["sf_fi"]] = ggplot(subset(fidata,label%in%c("GT","Sim")),aes(x=age,y=fi,colour=label,fill=label,lty=label))+
      #geom_point(size=.1)+
      stat_summary(mapping=aes(x=age_cut,shape=label),position=position_dodge(1))+
      geom_line(data=sfMean,mapping=aes(x=age,y=fi))+
      labs(x="Age",y="FI")+
      theme_minimal()
    
    dgtpar[is.na(dgtpar)] = gtpar[is.na(dgtpar)]
    parcomp = rbind(data.frame(par=c("W","mu0","mut","sigma","tau"),value=(gtpar-gtpar)/dgtpar,se=1,label="GT"),
                    data.frame(par=c("W","mu0","mut","sigma","tau"),value=(simpar-gtpar)/dgtpar,se=dsimpar/dgtpar,label="Sim"))
    plots[["sf_par"]] = ggplot(parcomp,aes(x=par,y=value,ymin=value-se,ymax=value+se,colour=label))+
      geom_pointrange(position=position_dodge(.1))+
      labs(y=bquote("(GT - Sim) /"~Delta*"GT)"),x="")+
      theme_minimal()
  }
  
  #spaghetti plots
  plotdata = subset(fidata,label=="GT")
  ids = unique(plotdata[,"id"])
  plotdata[,'firstvalue'] = NA
  for (i in 1:length(ids))
  {
    logi = plotdata[,"id"] == ids[i]
    ind = which.min(plotdata[logi,"age"])
    plotdata[logi,"firstvalue"] = plotdata[logi,"fi"][ind]
  }
  q = quantile(plotdata[,"firstvalue"],probs=seq(0,1,length=11),na.rm=T)
  q = q[!duplicated(q)]
  if(length(q) < 2) q = c(0,.02,.05,.1,.2,1)
  plotdata[,"q"] = cut(plotdata[,"firstvalue"],q,include.lowest=T)
  #plotdata[,"idrank"] = rank(plotdata[,"id"])
  #q = quantile(plotdata[,"idrank"])
  #plotdata[,"q"] = cut(plotdata[,"idrank"],q,include.lowest=T)
  plots[["gt_traj"]] = ggplot(plotdata,aes(x=age,y=fi,colour=q,fill=q,group=id))+
    #geom_smooth()+
    geom_line(alpha=.5)+
    labs(x="Age",y="FI",title="GT")+
    scale_colour_discrete_diverging("Blue-Red")+
    theme_minimal()
  plotdata = subset(fidata,label=="Sim")
  ids = unique(plotdata[,"id"])
  plotdata[,'firstvalue'] = NA
  for (i in 1:length(ids))
  {
    logi = plotdata[,"id"] == ids[i]
    ind = which.min(plotdata[logi,"age"])
    plotdata[logi,"firstvalue"] = plotdata[logi,"fi"][ind]
  }
  q = quantile(plotdata[,"firstvalue"],probs=seq(0,1,length=11),na.rm=T)
  q = q[!duplicated(q)]
  if(length(q) < 2) q = c(0,.02,.05,.1,.2,1)
  plotdata[,"q"] = cut(plotdata[,"firstvalue"],q,include.lowest=T)
  #plotdata[,"idrank"] = rank(plotdata[,"id"])
  #q = quantile(plotdata[,"idrank"])
  #plotdata[,"idrank"] = rank(plotdata[,"id"])
  plots[["sim_traj"]] = ggplot(plotdata,aes(x=age,y=fi,colour=q,fill=q,group=id))+
    #geom_smooth()+
    geom_line(alpha=.5)+
    labs(x="Age",y="FI",title="Sim")+
    scale_colour_discrete_diverging("Blue-Red")+
    theme_minimal()
  
  #`mean' crossing
  #step 1. define a function to track mean
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    warning(sprintf("mgcv library not found, skipping nullcline plots"))
    fidata[,"ddelta_mu"] = NA
    fidata[,"cross_mean"] =NA
    fidata[,"Ncross"] = NA
    fidata[,"observation_time"] = NA
    fidata[,"cross_per_t"] = NA
  }
  else
  {
    library(mgcv)
    gGT = gam(fi~s(age),data=subset(fidata,label=="GT"))
    gS = gam(fi~s(age),data=subset(fidata,label=="Sim"))
    gSz = gam(fi~s(age),data=subset(fidata,label=="Sim - w zombies"))
    #step 2. measure distance from mean at each time point
    fidata[,"delta_mu"] = NA
    logi = fidata[,"label"] == "GT"
    fidata[logi,"delta_mu"] = fidata[logi,"fi"] - predict(gGT,fidata[logi,])
    logi = fidata[,"label"] == "Sim"
    fidata[logi,"delta_mu"] = fidata[logi,"fi"] - predict(gS,fidata[logi,])
    logi = fidata[,"label"] == "Sim - w zombies"
    fidata[logi,"delta_mu"] = fidata[logi,"fi"] - predict(gSz,fidata[logi,])
    #step 3. check how offen each person crosses the mean
    fidata[,"ddelta_mu"] = NA #change in delta_mu
    fidata[,"cross_mean"] = NA #times crossed the mean
    fidata[,"Ncross"] = NA #times crossed the mean
    fidata[,"observation_time"] = NA #times crossed the mean
    fidata[,"cross_per_t"] = NA #times crossed the mean
    ids = unique(fidata[,"id"])
    labs = unique(fidata[,"label"]) #do separately for each label to avoid duplicated ids
    for (i in 1:length(ids))
    {
      for (j in 1:length(labs))
      {
        logi = fidata[,"id"] == ids[i] & fidata[,"label"] == labs[j]
        dmu = fidata[logi,"delta_mu"]
        fidata[logi,"ddelta_mu"] = c(diff(dmu),NA)
        cr = sign(dmu[-1]) != sign(dmu[-length(dmu)])
        fidata[logi,"cross_mean"] = c(cr,NA)
        fidata[logi,"Ncross"] = sum(cr,na.rm=T)
        fidata[logi,"observation_time"] = diff(range(fidata[logi,"age"],na.rm=T))
        fidata[logi,"cross_per_t"] = sum(cr,na.rm=T)/diff(range(fidata[logi,"age"],na.rm=T))
      }
    
    }
  }
  
  plots[["cross"]] = ggplot(fidata,aes(x=label,y=cross_per_t))+geom_boxplot()+theme_minimal()
  plots[["cross_total"]] = ggplot(fidata,aes(x=observation_time,y=Ncross,colour=label))+geom_point()+theme_minimal()
  
  keyplots = c("s0","Sd","Sr","ridge","rates_fi","rates_age","fi","fitraj","dfdt")
  rateplots = c("d","r","h","logd","logr","logh","dam_fi","rep_fi","haz_fi","dam_age","rep_age","haz_age")
  splots = c("s0","S","Sd","Sr","Sdfi","Sdfi")
  ridgeplots = c("ridge","ridge_gt","ridge_sim","ridge2")
  sfplots = c("sf_fi","sf_par")
  trajplots = c("fitraj","gt_traj","sim_traj","cross","cross_total")
  otherplots = setdiff(names(plots),c(keyplots,rateplots,ridgeplots,sfplots,trajplots))
  return(list(key=plots[keyplots],rates=plots[rateplots],splots=plots[splots],ridge=plots[ridgeplots],
              sf=plots[sfplots],traj=plots[trajplots],other=plots[otherplots],fidata=fidata,agg=agg,fitraj=fitraj,
              hazardData,sftab=sftab,sfMean=sfMean,sfm=sfm,sfmsim=sfmsim))
}


logLik.exp = function(par,
                      bn, #binary damage/not damage at tn ( matrix / dataframe)
                      bnp1, #binary damage/not damage at tn+1 ( matrix / dataframe)
                      tn, #vector of times
                      tnp1,
                      x, #covariates at time tn, usually this would be transformed fi matrix / dataframe
                      #f=apply(bn,1,mean,na.rm=T), #fi vector (at tn)
                      aliveLogi, #vector of alive to check for joint events
                      lastMeasurementLogi,
                      eps = 1e-8, #smallest allowed value of gompertz parameter (abs value) before assume 0
                      fix=NA #stub
)
{
  #log - generalized additive model for f with linear model for t
  #joint model
  #damage, repair and survival hazard functions are of the form: exp(basis + basis*t)
  
  #notes:
  #looks right, but there are strong correlations between the f and t-dependent terms, so you'll probably have to fix some of them when fitting
  
  #notes:
  #intercept should be included in x
  
  dt = tnp1-tn
  
  if(length(dim(bn))< 2) stop("bn must be a 2d array or data.frame")
  if(!all(dim(bnp1)==dim(bn))) stop("bnp1 must be same size as bn")
  if(length(dt)!=nrow(bnp1)) stop("dt should be aligned to b")
  
  Nx = ncol(x)
  damagePar =  par[1:Nx]
  damageParT = par[1:Nx+  Nx]
  repairPar =  par[1:Nx+2*Nx]
  repairParT = par[1:Nx+3*Nx]
  hazardPar =  par[1:Nx+4*Nx]
  hazardParT = par[1:Nx+5*Nx]
  
  
  lnDamageBasisHazard = x*NA
  damageGompertz = x*NA
  lnRepairBasisHazard = x*NA
  repairGompertz = x*NA
  lnSurvivalBasisHazard = x*NA
  survivalGompertz = x*NA
  
  for (j in 1:length(damagePar))
  {
    lnDamageBasisHazard[,j]    = x[,j]*damagePar [j]
    damageGompertz[,j]         = x[,j]*damageParT[j]
    lnRepairBasisHazard[,j]    = x[,j]*repairPar [j]
    repairGompertz[,j]         = x[,j]*repairParT[j]
    lnSurvivalBasisHazard[,j]  = x[,j]*hazardPar [j]
    survivalGompertz[,j]       = x[,j]*hazardParT[j]
  }
  lnDamageBasisHazard       = apply(lnDamageBasisHazard,1,sum,na.rm=T)
  damageGompertz            = apply(damageGompertz,1,sum,na.rm=T)
  lnRepairBasisHazard       = apply(lnRepairBasisHazard,1,sum,na.rm=T)
  repairGompertz            = apply(repairGompertz,1,sum,na.rm=T)
  lnSurvivalBasisHazard     = apply(lnSurvivalBasisHazard,1,sum,na.rm=T)
  #lnSurvivalBasisHazardLast = lnSurvivalBasisHazard[lastMeasurementLogi] #not needed
  survivalGompertz          = apply(survivalGompertz,1,sum,na.rm=T)
  #survivalGompertz = survivalGompertz[lastMeasurementLogi] #why is this here?
  
  #the only thing we have to check for is when the time-dependent term is close to 0, then we have to replace it with dt
  SDGompertz =  (exp(damageGompertz*dt)-1)/damageGompertz
  SDGompertz[abs(damageGompertz) < eps] = dt[abs(damageGompertz) < eps]
  Sd = exp(-exp(lnDamageBasisHazard + damageGompertz*tn)*SDGompertz) 
  
  SRGompertz =  (exp(repairGompertz*dt)-1)/repairGompertz
  SRGompertz[abs(repairGompertz) < eps] = dt[abs(repairGompertz) < eps]
  Sr = exp(-exp(lnRepairBasisHazard + repairGompertz*tn)*SRGompertz) 
  
  SGompertz =  (exp(survivalGompertz*dt)-1)/survivalGompertz #survGompertz << 0 -> 1/survGompertz -> 0
  SGompertz[abs(survivalGompertz) < eps] = dt[abs(survivalGompertz) < eps]
  lnS = -exp(lnSurvivalBasisHazard + survivalGompertz*tn)*SGompertz 
  
  
  #hazard when you die:
  deathHazard = exp(lnSurvivalBasisHazard[lastMeasurementLogi] + survivalGompertz[lastMeasurementLogi]*tnp1[lastMeasurementLogi])
  
  #note: if you have no deaths it'll just make survival -> 0 (as it should!)
  #print(sprintf("survival component: %.1e, hazard component: %.1e, par: %.1e",sum(lnS[aliveLogi],na.rm=T),sum(log(deathHazard),na.rm=T),hazardParT[1])) #for debug #if you make it hazardPaT << 0 this goes to 0 which is good
  l = sum(lnS[aliveLogi],na.rm=T) + sum(log(deathHazard),na.rm=T)
  #l = sum(log(deathHazard),na.rm=T) #debug
  for (j in 1:ncol(bn))
  {
    wasRepaired = bn[,j]  == 0
    isRepaired = bnp1[,j] == 0
    wasDamaged = bn[,j]   == 1
    isDamaged = bnp1[,j]  == 1
    #4 cases:
    #1. was damaged, stayed damage
    logi = wasDamaged & isDamaged & aliveLogi
    l = l + sum(log(Sr[logi]),na.rm=T)
    #2. was damaged but repaired
    logi = wasDamaged & isRepaired & aliveLogi
    l = l + sum(log(1-Sr[logi]),na.rm=T)
    #3. was repaired, stayed repaired
    logi = wasRepaired & isRepaired & aliveLogi
    l = l + sum(log(Sd[logi]),na.rm=T)
    #4. was repaired, but damaged
    logi = wasRepaired & isDamaged & aliveLogi
    l = l + sum(log(1-Sd[logi]),na.rm=T) 
  }
  
  return(l)
}

logLik.grad.exp = function(par,
                           bn, #binary damage/not damage at tn ( matrix / dataframe)
                           bnp1, #binary damage/not damage at tn+1 ( matrix / dataframe)
                           tn, #vector of times
                           tnp1,
                           x, #covariates at time tn, usually this would be transformed fi matrix / dataframe
                           #f=apply(bn,1,mean,na.rm=T), #fi vector (at tn)
                           aliveLogi, #vector of alive to check for joint events
                           lastMeasurementLogi,
                           eps = 1e-8, #smallest allowed value of gompertz parameter (abs value) before assume 0
                           fix=rep(F,length(par))
)
{
  
  
  dt = tnp1-tn
  
  if(length(dim(bn))< 2) stop("bn must be a 2d array or data.frame")
  if(!all(dim(bnp1)==dim(bn))) stop("bnp1 must be same size as bn")
  if(length(dt)!=nrow(bnp1)) stop("dt should be aligned to b")
  
  Nx = ncol(x)
  damagePar =  par[1:Nx]
  damageParT = par[1:Nx+  Nx]
  repairPar =  par[1:Nx+2*Nx]
  repairParT = par[1:Nx+3*Nx]
  hazardPar =  par[1:Nx+4*Nx]
  hazardParT = par[1:Nx+5*Nx]
  
  
  lnDamageBasisHazard = x*NA
  damageGompertz = x*NA
  lnRepairBasisHazard = x*NA
  repairGompertz = x*NA
  lnSurvivalBasisHazard = x*NA
  survivalGompertz = x*NA
  
  for (j in 1:length(damagePar))
  {
    lnDamageBasisHazard[,j]    = x[,j]*damagePar [j]
    damageGompertz[,j]         = x[,j]*damageParT[j]
    lnRepairBasisHazard[,j]    = x[,j]*repairPar [j]
    repairGompertz[,j]         = x[,j]*repairParT[j]
    lnSurvivalBasisHazard[,j]  = x[,j]*hazardPar [j]
    survivalGompertz[,j]       = x[,j]*hazardParT[j]
  }
  lnDamageBasisHazard       = apply(lnDamageBasisHazard,1,sum,na.rm=T)
  damageGompertz            = apply(damageGompertz,1,sum,na.rm=T)
  lnRepairBasisHazard       = apply(lnRepairBasisHazard,1,sum,na.rm=T)
  repairGompertz            = apply(repairGompertz,1,sum,na.rm=T)
  lnSurvivalBasisHazard     = apply(lnSurvivalBasisHazard,1,sum,na.rm=T)
  #lnSurvivalBasisHazardLast = lnSurvivalBasisHazard[lastMeasurementLogi] #not needed
  survivalGompertz          = apply(survivalGompertz,1,sum,na.rm=T)
  #survivalGompertz = survivalGompertz[lastMeasurementLogi] #why is this here?
  
  #the only thing we have to check for is when the time-dependent term is close to 0, then we have to replace it with dt
  SDGompertz =  (exp(damageGompertz*dt)-1)/damageGompertz
  SDGompertz[abs(damageGompertz) < eps] = dt[abs(damageGompertz) < eps]
  Sd = exp(-exp(lnDamageBasisHazard + damageGompertz*tn)*SDGompertz) 
  
  damageGompertzGrad = (damageGompertz*tn-1+(1-damageGompertz*tnp1)*exp(damageGompertz*dt))/damageGompertz^2
  damageGompertzGrad[abs(damageGompertz) < eps] = -dt[abs(damageGompertz) < eps]
  
  SRGompertz =  (exp(repairGompertz*dt)-1)/repairGompertz
  SRGompertz[abs(repairGompertz) < eps] = dt[abs(repairGompertz) < eps]
  Sr = exp(-exp(lnRepairBasisHazard + repairGompertz*tn)*SRGompertz) 
  
  repairGompertzGrad = (repairGompertz*tn-1+(1-repairGompertz*tnp1)*exp(repairGompertz*dt))/repairGompertz^2
  repairGompertzGrad[abs(repairGompertz) < eps] = -dt[abs(repairGompertz) < eps]
  
  SGompertz =  (exp(survivalGompertz*dt)-1)/survivalGompertz
  SGompertz[abs(survivalGompertz) < eps] = dt[abs(survivalGompertz) < eps]
  lnS = -exp(lnSurvivalBasisHazard + survivalGompertz*tn)*SGompertz 
  
  hazardGompertzGrad = (survivalGompertz*tn-1+(1-survivalGompertz*tnp1)*exp(survivalGompertz*dt))/survivalGompertz^2
  hazardGompertzGrad[abs(survivalGompertz) < eps] = -dt[abs(survivalGompertz) < eps]
  
  D = exp(lnDamageBasisHazard + damageGompertz*tn)
  R = exp(lnRepairBasisHazard + repairGompertz*tn)
  
  gradDamage = rep(0,length(damagePar))
  gradDamageT = rep(0,length(damageParT))
  gradRepair = rep(0,length(repairPar))
  gradRepairT = rep(0,length(repairParT))
  for (j in 1:ncol(bn))
  {
    wasRepaired = bn[,j]  == 0
    isRepaired  = bnp1[,j] == 0
    wasDamaged  = bn[,j]   == 1
    isDamaged   = bnp1[,j]  == 1
    
    wasRepIsRep = wasRepaired & isRepaired & aliveLogi
    wasRepIsDam = wasRepaired & isDamaged & aliveLogi 
    wasDamIsDam = wasDamaged & isDamaged & aliveLogi
    wasDamIsRep = wasDamaged & isRepaired & aliveLogi
    
    for (jj in 1:ncol(x))
    {
      #damage - Sd
      gradDamage[jj] = gradDamage[jj]+sum(-x[wasRepIsRep,jj]*D[wasRepIsRep]*SDGompertz[wasRepIsRep], na.rm=T)
      gradDamage[jj] = gradDamage[jj]-sum(-x[wasRepIsDam,jj]*D[wasRepIsDam]*SDGompertz[wasRepIsDam]*Sd[wasRepIsDam]/(1-Sd[wasRepIsDam]),na.rm=T) 
      
      gradDamageT[jj] = gradDamageT[jj]+sum(x[wasRepIsRep,jj]*D[wasRepIsRep]*damageGompertzGrad[wasRepIsRep], na.rm=T)
      gradDamageT[jj] = gradDamageT[jj]-sum(x[wasRepIsDam,jj]*D[wasRepIsDam]*damageGompertzGrad[wasRepIsDam]*Sd[wasRepIsDam]/(1-Sd[wasRepIsDam]),na.rm=T) 
      
      #repair - Sr
      gradRepair[jj] = gradRepair[jj]+sum( -x[wasDamIsDam,jj]*R[wasDamIsDam]*SRGompertz[wasDamIsDam], na.rm=T)
      gradRepair[jj] = gradRepair[jj]-sum( -x[wasDamIsRep,jj]*R[wasDamIsRep]*SRGompertz[wasDamIsRep]*Sr[wasDamIsRep]/(1-Sr[wasDamIsRep]),na.rm=T)
      
      
      gradRepairT[jj] = gradRepairT[jj]+sum(x[wasDamIsDam,jj]*R[wasDamIsDam]*repairGompertzGrad[wasDamIsDam], na.rm=T)
      gradRepairT[jj] = gradRepairT[jj]-sum(x[wasDamIsRep,jj]*R[wasDamIsRep]*repairGompertzGrad[wasDamIsRep]*Sr[wasDamIsRep]/(1-Sr[wasDamIsRep]),na.rm=T)
      
    }
  }
  
  
  hazard = exp(lnSurvivalBasisHazard + survivalGompertz*tn)
  
  gradHazard = rep(0,length(hazardPar))
  gradHazardT = rep(0,length(hazardParT))
  
  for (jj in 1:ncol(x))
  {
    gradHazard[jj] = sum(  -x[aliveLogi,jj]*hazard[aliveLogi]*SGompertz[aliveLogi],na.rm=T)+ #from S
      sum( x[lastMeasurementLogi,jj],na.rm=T) #from h
    
    gradHazardT[jj] = sum(  x[aliveLogi,jj]*hazard[aliveLogi]*hazardGompertzGrad[aliveLogi],na.rm=T)+ #from S
      sum( x[lastMeasurementLogi,jj]*tnp1[lastMeasurementLogi],na.rm=T) #from h
  }
  
  
  gr = c(gradDamage,gradDamageT,gradRepair,gradRepairT,gradHazard,gradHazardT)
  gr[fix] = 0
  
  return(gr)
}

QuickSimExp = function(N=1000, #number of individuals
                       p=30, #number of nodes
                       DamageFun=DamageFun0, #depends only on f #hazard for damage
                       DamageGompertzFun=function(f){rep(0,length(f))},
                       RepairFun=RepairFun0, #depends only on f #hazard for repair
                       RepairGompertzFun=function(f){rep(0,length(f))},
                       SurvivalFun=SurvivalFun0, #depends only on f #hazard for death
                       SurvivalGompertzFun=function(f){rep(0.089,length(f))},
                       t0=0, #time or vector of times to shift individuals
                       twindow=c(0,150), #window to use for conversion to uniform sampling data
                       dt=.5, #timesteps to use for conversion to uniform sampling (to look like wave data); doesn't affect base sim #needs to be pretty small
                       innerMultiplier=10, #runs inner loop for dt of step size dt/innerMultiplier -> these data are not saved (except survival)
                       b0=rep(0,p), #initial deficits
                       eps = 1e-8, #for finding gompertz cases
                       saveQALY=T #compute and save QALY
                       #saveFI = F #will keep FI at ALL TIMES
)
{
  #uses discrete sampling
  #assumes no intermediate events between time points dt
  #so dt should be small!
  # only saves a portion of the data
  
  #to do:
  #validate start-stop formatting
  #esp make sure no duplicate deaths
  
  #new approach, for damage/repair just sample directly using 1-survival
  
  
  
  matrixStart = F
  if(is.null(dim(b0)) & length(b0) != p) stop("p and length(b0) must be same (or b0 can be matrix with ncol = p)")
  else if(!is.null(dim(b0)))
  {
    matrixStart=T
    if(nrow(b0) != N | ncol(b0)!= p)  stop("p and ncol(b0) must be same and N and nrow(b0) must be same")
    b0 = as.matrix(b0) #in case was supplied as dataframe
  }
  if(length(t0) != 1 & length(t0) != N) stop("t0 must be length 1 or N")
  if(length(t0)==1) t0 = rep(t0,N)
  
  #I don't think I want or need this
  #if(matrixStart) Ndef = apply(b0,1,sum,na.rm=T)
  #else Ndef = sum(b0)
  #if(any(DamageFun(Ndef)<=1e-30) & !skipDamageCheck) stop("Initial damage rate is 0 for some individuals, this causes a fixed point at 0 and is not acceptable")
  
  
  QALY = NULL
  if(saveQALY) QALY = rep(0,N)
  
  times = seq(twindow[1],twindow[2],by=dt)
  waveAge = array(NA,dim=c(N,length(times))) #age at each sampling wave - or age at death if they are dead at wave
  waveAge[,1] = times[1]+t0
  waveData0 = array(NA,dim=c(N,p,length(times))) 
  #print(dim(waveData0))
  #initialize
  if(matrixStart) waveData0[,,1] = b0
  else for (j in 1:length(b0)) waveData0[,j,1] = b0[j]
  
  alive = array(T,dim=c(N,length(times)))
  tdeath = rep(NA,N)
  
  fi = apply(waveData0[,,1],1,mean,na.rm=T)  
  fi_t = rep(NA,length(times))
  fi_t[1] = mean(fi,na.rm=T)
  #currentWave = waveData0[,,1]
  #prevWave = waveData0[,,1]
  stst = list() #start-stop formatting for survival
  stst[[1]] = data.frame(id=1:N,age = waveAge[,1],fi=fi,status=rep(0,N),alreadyDead=F)
  teventlast = matrix(NA,nrow=N,ncol=p) #keeps track of the time of the last event that occurred
  for (j in 1:p) teventlast[,j] = waveAge[,1]
  teventprev = teventlast #keeps track of the time of the previous event that occurred
  events = data.frame(t0=NA,t=NA,event=NA) #data.frame of events
  tfp = NA*teventlast #first passage times for damage
  for (k in 2:length(times))
  {
    #print(k)
    #print(sum(is.na(waveData0[,,k-1])))
    waveAge[,k] = times[k]+t0
    stst[[k-1]][,"age_next"] = times[k]+t0
    
    #forward carry
    waveData0[,,k] = waveData0[,,k-1]
    alive[,k] = alive[,k-1]
    
    #note: multi-events not allowed, so no damage-repair-damage etc
    if(innerMultiplier >= 1) #changed to >= not 100% sure this is right
    {
      currentWave = waveData0[,,k]
      previousWave = waveData0[,,k]*NA
      dt = (times[k]-times[k-1])/innerMultiplier
      aliveNow = alive[,k-1]
      for (i in 1:innerMultiplier)
      {
        tprev = times[k-1]+t0+dt*(i-1)
        tcur = tprev+dt
        
        
        
        #sample damage
        #directly sample from survival probability (i.e. prob of surviving dt)
        #exponential model with DamageFun for hazard
        damageGompertz = DamageGompertzFun(fi)
        d = DamageFun(fi)*exp(damageGompertz*tprev) 
        gompertzTerm = (exp(damageGompertz*dt)-1)/damageGompertz
        gompertzTerm[abs(damageGompertz) < eps] = dt #special case when gompertz term ->0
        Sd = exp(-d*gompertzTerm)
        #sample all nodes for damage
        for (j in 1:ncol(currentWave))
        {
          logi = currentWave[,j] == 0
          logit = runif(sum(logi),0,1) < 1-Sd[logi] #didn't survive damaging
          currentWave[logi,j][logit] = 1
          #update: keep track of events
          if(sum(logi[logit])>0)
          {
            temp = tprev[logi][logit] + dt/2
            teventlast[logi,j][logit] = temp
            events = rbind(events, data.frame(t0=teventprev[logi,j][logit],t=temp,event="d"))
            virg = is.na(tfp[logi,j][logit])
            tfp[logi,j][logit][virg] = temp[virg]
          }
        }
        
        #sample repair
        #exponential model with DamageFun for hazard
        repairGompertz = RepairGompertzFun(fi)
        r = RepairFun(fi)*exp(repairGompertz*tprev)
        gompertzTerm = (exp(repairGompertz*dt)-1)/repairGompertz
        gompertzTerm[abs(repairGompertz) < eps] = dt #special case when gompertz term ->0
        Sr = exp(-r*gompertzTerm)
        for (j in 1:ncol(currentWave))
        {
          logi = currentWave[,j] == 1
          logit = runif(sum(logi),0,1) < 1-Sr[logi] #didn't survive repairing
          currentWave[logi,j][logit] = 0
          #update: keep track of events
          if(sum(logi[logit])>0)
          {
            temp = tprev[logi][logit] + dt/2
            teventlast[logi,j][logit] = temp
            events = rbind(events, data.frame(t0=teventprev[logi,j][logit],t=temp,event="r"))
          }
        }
        
        
        
        #sample survival
        survivalGompertz = SurvivalGompertzFun(fi[aliveNow])
        s = SurvivalFun(fi[aliveNow])*exp(survivalGompertz*tprev[aliveNow])
        ts =   (- log(runif(sum(aliveNow),0,1))/s) #gompertz = 0 case (default)
        #positiont alpha Gompertz case:
        gompLogi = survivalGompertz > eps
        ts[gompLogi] = log(1-survivalGompertz[gompLogi]/s[gompLogi]*log(runif(sum(gompLogi),0,1)))/survivalGompertz[gompLogi]
        #negative alpha Gompertz case:
        alphaLogi = -survivalGompertz > eps 
        #umin = exp(s[alphaLogi]/survivalGompertz[alphaLogi]) #not right, right shape but wrong scale - use accept/reject instead
        #ts[alphaLogi] = log(1-survivalGompertz[alphaLogi]/s[alphaLogi]*log(runif(sum(alphaLogi),umin,1)))/survivalGompertz[alphaLogi]    
        ts2 = 1-survivalGompertz[alphaLogi]/s[alphaLogi]*log(runif(sum(alphaLogi),0,1))
        ts2[ts2 < 0] = 0 #these will be rejects (ts -> Inf)
        ts[alphaLogi] = log(ts2)/survivalGompertz #negative gompertz
        
        deaths = ts < dt
        tdeath[aliveNow][deaths] = ts[deaths] + tprev[aliveNow][deaths]
        waveAge[aliveNow,k][deaths] = ts[deaths] + tprev[aliveNow][deaths]
        aliveNow[aliveNow][deaths] = F
        
        
        #update fi
        fi0 = fi
        fi = apply(currentWave,1,mean,na.rm=T)
        previousWave = currentWave
        if(saveQALY) QALY = QALY + (1-fi/2-fi0/2)*dt*aliveNow
        
        teventprev = teventlast
      }
    }
    
    waveData0[,,k] = currentWave
    alive[,k] = aliveNow
    fi = apply(waveData0[,,k],1,mean,na.rm=T)  
    fi_t[k] = mean(fi[alive[,k]],na.rm=T)
    stst[[k]] = data.frame(id=1:N,age = waveAge[,k],age_next=NA,fi=fi,status=1-as.integer(alive[,k]),alreadyDead=!alive[,k-1]) 
  }
  stst = do.call(rbind,stst)
  stst = subset(stst,!alreadyDead) #only include people who were alive at start of wave
  
  #survival
  tevent = tdeath
  tmax = twindow[2]+t0
  tevent[is.na(tevent)] = tmax[is.na(tevent)] #censors
  tevent[tevent > tmax] = tmax[tevent > tmax] 
  status = as.integer(tdeath <= tmax) #deaths
  status[is.na(tdeath)] = 0 #censors
  
  s = Surv(twindow[1]+t0,tevent,status)
  sstst = Surv(stst[,"age"],stst[,"age_next"],event=stst[,"status"])
  
  events = events[-1,] #first one is just place-holder
  
  #first passage survival (damage)
  fpevent = !is.na(c(tfp))
  for (j in 1:ncol(tfp)) tfp[is.na(tfp[,j]),j] = tevent[is.na(tfp[,j])]
  Sdfp = Surv(time=c(outer(t0,rep(1,p))),time2=c(tfp),event=fpevent)
  
  #damage 'survival'
  de = subset(events,event=="d")
  Sd = Surv(de[,"t"]-de[,"t0"])
  
  #repair 'survival'
  re = subset(events,event=="r")
  Sr = Surv(re[,"t"]-re[,"t0"])
  
  #make a censored version
  waveData = waveData0
  for (j in 1:ncol(waveData)) waveData[,j,][!alive] = NA
  
  l = list()
  
  l[["wave"]] = waveData #observed events
  l[["wave0"]] = waveData0 #all events even after death
  l[["age"]] = waveAge #times of waves, with time of death if you died instead of being measured
  l[["t"]] = times #time relative to t0
  l[["t0"]] = t0 #initial time
  l[["s"]] = s #observed survival
  l[["fi"]] = fi_t
  l[["deaths"]] = deaths
  l[["alive"]] = alive
  l[["par"]] = c(N=N,p=p,dt=dt,tmin=twindow[1],tmax=twindow[2])
  l[["DamageFun"]]=DamageFun
  l[["RepairFun"]]=RepairFun
  l[["SurvivalFun"]]=SurvivalFun
  l[["DamageGompertzFun"]]=DamageGompertzFun
  l[["RepairGompertzFun"]]=RepairGompertzFun
  l[["SurvivalGompertzFun"]]=SurvivalGompertzFun
  l[["stst"]]=stst #start-stop formatted survival and FI info
  l[["sstst"]]=sstst #survival in start-stop format
  l[["tdeath"]] = tdeath
  l[["tevent"]] = tevent
  l[["q"]] = QALY
  l[["qn"]] = QALY/(s[,"stop"]-s[,"start"])
  l[["events"]] = events[-1,] #first one is just place-holder
  l[["Sd"]] = Sd #unvalidated
  l[["Sr"]] = Sr #unvalidated
  l[["Sdfp"]] = Sdfp #unvalidated #survival of damage/ first passage i.e. first damage
  return(l)
}


DamageAndRepair = function(thiswave,
                           nextwave,
                           rdcols=c("n_damage","damage_rate","n_repair","repair_rate",
                                    "n_damage_se","damage_rate_se","n_repair_se","repair_rate_se",
                                    "can_damage","can_repair")
)
{
  x = data.frame(id=1:nrow(thiswave))
  
  canDamage = thiswave == 0 #at risk
  canDamageSum = apply(canDamage,1,sum,na.rm=T)
  x[,"n_damage"] = apply(nextwave == 1 & canDamage,1,sum,na.rm=T)
  x[,"damage_rate"] = x[,"n_damage"]/canDamageSum
  #standard error is sqrt(n*p*(1-p)/N) #n: number of trials; N: number of samples, p = unknown prob to estimate #N=1
  #mean is n*p
  x[,"n_damage_se"] = sqrt(canDamageSum*x[,"damage_rate"]*(1-x[,"damage_rate"]))
  x[,"damage_rate_se"] = x[,"n_damage_se"]/canDamageSum
  x[,"can_damage"] = canDamageSum
  noDam = apply(canDamage,1,sum,na.rm=T)==0
  noDam[is.na(noDam)] = F
  x[noDam,"damage_rate"] = NA
  x[noDam,"n_damage_se"] = NA
  
  canRepair = thiswave == 1 #at risk
  canRepairSum = apply(canRepair,1,sum,na.rm=T)
  x[,"n_repair"] = apply(nextwave == 0 & canRepair,1,sum,na.rm=T)
  x[,"repair_rate"] = x[,"n_repair"]/canRepairSum
  x[,"n_repair_se"] = sqrt(canRepairSum*x[,"repair_rate"]*(1-x[,"repair_rate"]))
  x[,"repair_rate_se"] = x[,"n_repair_se"]/canRepairSum
  x[,"can_repair"] = canRepairSum
  noRep = apply(canRepair,1,sum,na.rm=T)==0
  noRep[is.na(noRep)] = F
  x[noRep,"repair_rate"] = NA
  
  #check for special case: people who died i.e. will die in next wave
  #add new columns encoded with fi = 1 at/after death
  x[,sprintf("%s_deathfi",rdcols)] = x[,rdcols]
  
  
  cols = setdiff(colnames(x),"id")
  return(x[,cols])
}


ListMeanSD = function(l,sem=F,negToZero=F,na.rm=F,skipsd=F)
{
  #calculates SD of list of objects of same size (e.g. matrices/arrays)
  #preserves dimensions!
  #returns mean and sd or sem
  #negToZero: force negative variance to 0 (happens if mu^2 ~= v)
  
  if(length(l)==1) 
  {
    if(is.null(dim(l[[1]]))) return(list(mean=l[[1]],sd=rep(0,length(l))))
    else return(list(mean=l[[1]],sd=array(0,dim(l[[1]]))))
  }
  else if (length(l)<1) return(list(mean=NA,sd=NA))
  else if(na.rm)
  {
    lnotNA = l
    for (i in 1:length(l)) 
    {
      lnotNA[[i]] = !is.na(l[[i]])
      l[[i]][is.na(l[[i]])] = 0
    }
    numNotNA = Reduce("+",lnotNA)
    
    mu = Reduce("+",l)/numNotNA
    l = lapply(l,function(x){return(x^2)})
    
    if(skipsd) return(list(mean=mu,sd=NULL))
    
    v = Reduce("+",l)/numNotNA
    v = v - mu^2
    v = v*numNotNA/(numNotNA-1) #convert to unbiased estimate
    
    s = sqrt(v)
    
    if(negToZero)
    {
      v[is.nan(v)]=NA
      s[v<0] = 0
    }
    
    if(sem) s = s/sqrt(numNotNA)
    return(list(mean=mu,sd=s))
  }
  else
  {
    mu = Reduce("+",l)/length(l)
    l = lapply(l,function(x){return(x^2)})
    
    if(skipsd) return(list(mean=mu,sd=NULL))
    
    v = Reduce("+",l)/length(l)
    v = v - mu^2
    v = v*length(l)/(length(l)-1) #convert to unbiased estimate
    
    s = sqrt(v)
    
    if(negToZero)
    {
      v[is.nan(v)]=NA
      s[v<0] = 0
    }
    
    if(sem) s = s/sqrt(length(l))
    return(list(mean=mu,sd=s))
  }
  
}



DamageFun.exp = function(f,t,par,x)
{
  Nx = ncol(x)
  damagePar =  par[1:Nx]
  damageParT = par[1:Nx+  Nx]
  lnDamageBasisHazard = x*NA
  damageGompertz = x*NA
  
  for (j in 1:length(damagePar))
  {
    lnDamageBasisHazard[,j]    = x[,j]*damagePar [j]
    damageGompertz[,j]         = x[,j]*damageParT[j]
  }
  lnDamageBasisHazard       = apply(lnDamageBasisHazard,1,sum,na.rm=T)
  damageGompertz            = apply(damageGompertz,1,sum,na.rm=T)
  
  d = exp(lnDamageBasisHazard + damageGompertz*t)
  return(d)
}
Sd.exp  = function(f,t,par,x=NULL,eps=1e-8)
{
  warning("unvalidated")
  Nx = ncol(x)
  damagePar =  par[1:Nx]
  damageParT = par[1:Nx+  Nx]
  lnDamageBasisHazard = x*NA
  damageGompertz = x*NA
  
  for (j in 1:length(damagePar))
  {
    lnDamageBasisHazard[,j]    = x[,j]*damagePar [j]
    damageGompertz[,j]         = x[,j]*damageParT[j]
  }
  lnDamageBasisHazard       = apply(lnDamageBasisHazard,1,sum,na.rm=T)
  damageGompertz            = apply(damageGompertz,1,sum,na.rm=T)
  
  SDGompertz =  (exp(damageGompertz*t)-1)/damageGompertz
  SDGompertz[abs(damageGompertz) < eps] = t[abs(damageGompertz) < eps]
  sd = exp(-exp(lnDamageBasisHazard)*SDGompertz) 
  return(sd)
}
RepairFun.exp  = function(f,t,par,x=NULL)
{
  Nx = ncol(x)
  repairPar =  par[1:Nx+2*Nx]
  repairParT = par[1:Nx+3*Nx]
  
  lnRepairBasisHazard = x*NA
  repairGompertz = x*NA
  
  for (j in 1:length(repairPar))
  {
    lnRepairBasisHazard[,j]    = x[,j]*repairPar [j]
    repairGompertz[,j]         = x[,j]*repairParT[j]
    
  }
  lnRepairBasisHazard       = apply(lnRepairBasisHazard,1,sum,na.rm=T)
  repairGompertz            = apply(repairGompertz,1,sum,na.rm=T)
  
  r = exp(lnRepairBasisHazard + repairGompertz*t)
  return(r)
}
Sr.exp  = function(f,t,par,x=NULL,eps=1e-8)
{
  warning("unvalidated")
  Nx = ncol(x)
  hazardPar =  par[1:Nx+4*Nx]
  hazardParT = par[1:Nx+5*Nx]
  
  
  lnSurvivalBasisHazard = x*NA
  survivalGompertz = x*NA
  
  for (j in 1:length(repairPar))
  {
    lnSurvivalBasisHazard[,j]  = x[,j]*hazardPar [j]
    survivalGompertz[,j]       = x[,j]*hazardParT[j]
  }
  lnSurvivalBasisHazard     = apply(lnSurvivalBasisHazard,1,sum,na.rm=T)
  survivalGompertz          = apply(survivalGompertz,1,sum,na.rm=T)
  
  SRGompertz =  (exp(repairGompertz*t)-1)/repairGompertz
  SRGompertz[abs(repairGompertz) < eps] = t[abs(repairGompertz) < eps]
  sr = exp(-exp(lnRepairBasisHazard )*SRGompertz) 
  
  return(sr)
}

SurvivalFun.exp  = function(f,t,par,x=NULL)
{
  Nx = ncol(x)
  hazardPar =  par[1:Nx+4*Nx]
  hazardParT = par[1:Nx+5*Nx]
  
  
  lnSurvivalBasisHazard = x*NA
  survivalGompertz = x*NA
  
  for (j in 1:length(hazardPar))
  {
    lnSurvivalBasisHazard[,j]  = x[,j]*hazardPar [j]
    survivalGompertz[,j]       = x[,j]*hazardParT[j]
  }
  lnSurvivalBasisHazard     = apply(lnSurvivalBasisHazard,1,sum,na.rm=T)
  survivalGompertz          = apply(survivalGompertz,1,sum,na.rm=T)
  h = exp(lnSurvivalBasisHazard + survivalGompertz*t)
  return(h)
}

S.exp  = function(f,t,par,x=NULL,eps=1e-8)
{
  Nx = ncol(x)
  hazardPar =  par[1:Nx+4*Nx]
  hazardParT = par[1:Nx+5*Nx]
  
  
  lnSurvivalBasisHazard = x*NA
  survivalGompertz = x*NA
  
  for (j in 1:length(hazardPar))
  {
    lnSurvivalBasisHazard[,j]  = x[,j]*hazardPar [j]
    survivalGompertz[,j]       = x[,j]*hazardParT[j]
  }
  lnSurvivalBasisHazard     = apply(lnSurvivalBasisHazard,1,sum,na.rm=T)
  survivalGompertz          = apply(survivalGompertz,1,sum,na.rm=T)
  
  SGompertz =  (exp(survivalGompertz*t)-1)/survivalGompertz
  SGompertz[abs(survivalGompertz) < eps] = t[abs(survivalGompertz) < eps]
  
  s = exp(-exp(lnSurvivalBasisHazard)*SGompertz )
  return(s)
}


EulerFI = function(t, #times to estimate at
                   DamageFun, #function of current FI
                   RepairFun,
                   f0 = 0,
                   damageGompertz=0,
                   repairGompertz=0
)
{
  f = rep(f0,length(t))
  for (i in 2:length(t))
  {
    dt = t[i]-t[i-1]
    f[i] = f[i-1] + (1-f[i-1])*DamageFun(f[i-1])*exp(damageGompertz*t[i-1])*dt-f[i-1]*RepairFun(f[i-1])*exp(repairGompertz*t[i-1])*dt
  }
  return(f)
}

SEM = function(x,na.rm=F)
{
  if(na.rm) x = x[!is.na(x)]
  return(sd(x)/sqrt(length(x)))
}


unlen = function(x,na.rm=T) 
{
  un = unique(x)
  if(na.rm) un = un[!is.na(un)]
  return(length(un))
}


WavesToTransitions = function(current_wave, #matrix
                              next_wave, #matrix
                              age, #vector
                              age_next, #vector
                              baseline_age=NA,
                              fi_cuts = seq(0,1,by=.2),
                              waves = NULL #3d array #optional, provide instead of current_wave/next_wave
)
{
  library(survival)
  #converts waves of binary data into time-to-event transitions
  #uses start-stop formatting
  #note: I've had trouble with this vs tmerge before so this is very unvalidated
  #looks reasonable but hard to fully validate
  #warning("unvalidated")
  
  #not sure if this is exact but its definitely close
  
  if(!is.null(waves))
  {
    
    current_wave = list()
    next_wave = list()
    for (k in 2:dim(waves)[3])
    {
      current_wave[[k-1]] = waves[,,k-1]
      next_wave[[k-1]] = waves[,,k]
    }
    current_wave = do.call(rbind,current_wave)
    next_wave = do.call(rbind,next_wave)
  }
  
  data = list()
  for (j in 1:ncol(current_wave))
  {
    data[[j]] = data.frame(current_wave=current_wave[,j],next_wave=next_wave[,j])
    data[[j]][,"fi"] = apply(current_wave,1,mean,na.rm=T)
    data[[j]][,"age"] = age
    data[[j]][,"age_next"] = age_next
    data[[j]][,"baseline_age"] = baseline_age
  }
  data = do.call(rbind,data)
  data[,"transition"] = NA
  current_one = data[,"current_wave"] == 1
  current_zero = !current_one
  next_one = data[,"next_wave"] == 1
  next_zero = !next_one
  current_one[is.na(current_one)] = F
  current_zero[is.na(current_zero)] = F
  next_one[is.na(next_one)] = F
  next_zero[is.na(next_zero)] = F
  
  #survive repair
  data[current_one & next_one,"transition"] = "sr"
  #survive damage
  data[current_zero & next_zero,"transition"] = "sd"
  #repair
  data[current_one & next_zero,"transition"] = "r"
  #damage
  data[current_zero & next_one,"transition"] = "d"
  
  #time-to-event stats
  #use only at-risk
  Sr = Surv(data[current_one,"age"],data[current_one,"age_next"],event=as.integer(data[current_one,"transition"]=="r"))
  Sd = Surv(data[current_zero,"age"],data[current_zero,"age_next"],event=as.integer(data[current_zero,"transition"]=="d"))
  St = Surv(data[,"age"],data[,"age_next"],event=as.integer( (current_one & data[,"transition"]=="r") | (current_zero & data[,"transition"]=="d"))) #either
  #St unvalidated and potentially wrong
  
  fi = cut(data[,"fi"],fi_cuts,labels=fi_cuts[-1]/2+fi_cuts[-length(fi_cuts)]/2,include.lowest=T)
  fir = fi[current_one]
  fid = fi[current_zero]
  sfr = survfit(Sr~fir)
  sfd = survfit(Sd~fid)
  sft = survfit(St~fi)
  
  #survival terms that start from 0 #doesn't make sense... every step is dt
  #Sr0 = Surv(data[current_one,"age"]-data[current_one,"age"],data[current_one,"age_next"]-data[current_one,"age"],
  #           event=as.integer(data[current_one,"transition"]=="r"))
  #Sd0 = Surv(data[current_zero,"age"]-data[current_zero,"age"],data[current_zero,"age_next"]-data[current_zero,"age"],
  #           event=as.integer(data[current_zero,"transition"]=="d"))
  #sfr0 = survfit(Sr0~fir)
  #sfd0 = survfit(Sd0~fid)
  
  return(list(Sr=Sr,Sd=Sd,St=St,sfr=sfr,sfd=sfd,sft=sft,
              #Sr0=Sr0,Sd0=Sd0,sfr0=sfr0,sfd0=sfd0, #makes no sense
              data=data,current_one=current_one,current_zero=current_zero,
              #datar=data[current_one,], datad=data[current_zero,], #too much data just use prev line
              fi=fi,fir=fir,fid=fid,fi_cuts=fi_cuts))
}