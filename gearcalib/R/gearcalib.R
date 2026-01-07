##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Checks the a gear calib data set has the correct format
##' @param d A 'gear calib list'
##' @return nothing
check.gearcalib.data<-function(d){
    stopifnot(class(d)=="list")
    stopifnot(all(c("N","SweptArea","group","Gear")  %in%  names(d)))
    stopifnot(class(d$N)[1]=="matrix")
    #stopifnot(class(d$SweptArea)=="numeric")
    stopifnot(class(d$SweptArea)[1]=="matrix")
    stopifnot(class(d$group)=="factor")
    stopifnot(class(d$Gear)=="factor")
    no.hauls <- nrow(d$N)
    stopifnot(dim(d$N)==dim(d$SweptArea))
    stopifnot(nlevels(d$group)==no.hauls/2)
    stopifnot(nlevels(d$Gear)==2)
    stopifnot(all(is.finite(d$SweptArea)))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Fit relative catch efficiency by length group between to gears
##' @param d A list with four elements: N (matrix of integers), SweptArea(numeric vector), group(factor vector), and Gear(factor vector). 
##' @param fit0 If TRUE a Chisq-test of no size structure in gear effect is performed.
##' @return A list
gearcalibFit <- function(d,fit0=FALSE,logsd=NA,phi=NA,logsdnug=NA,logsdres=NA,logsdGearRW=NA,logalpha=0)
    {
        check.gearcalib.data(d)
        nsize <- ncol(d$N)
        ngroup <- nlevels(d$group)
        nhaul <- nrow(d$N)
        ngear <- nlevels(d$Gear)

        ## Set default RW order
        rw_order <- c(1,1)

        data <- list(
            N=d$N,
            SweptArea=d$SweptArea,
            group=d$group,
            Gear=d$Gear,
            huge=10,
            tiny=0.01,
            rw_order=rw_order
            )

        parameters=list(
            logspectrum=matrix(0,ngroup,nsize),
            nugget=matrix(0,nhaul,nsize),
            residual=matrix(0,nhaul,nsize),
            loggear=numeric(nsize),
            logsd=-1,
            phi=0.9,
            logsdnug=-1,
            logsdres=-1,
            logsdGearRW=-1,
            logalpha = 0
            )

        random <- c("logspectrum","residual","loggear","nugget")

        map <- list()
        
        ## If any parameters have been specified in the call, do not estimate those
        ## parameters but fix them to the specified value
        setparameter <- function(name) {
            var <- get(name)
            if(!is.na(var))
            {
                map[[name]] <<- factor(NA)
                parameters[[name]] <<- var
            }
        }
        parameternames <- c("logsd","phi","logsdnug","logsdres","logsdGearRW","logalpha")
        sapply(parameternames,setparameter)

        
        obj <- MakeADFun(
            data = data,
            parameters = parameters,
            map = map, 
            random = random,
            DLL="gearcalib"
        )

        
        lower <- 0*obj$par-Inf
        upper <- 0*obj$par+Inf
        if(any("phi" == names(obj$par)))
            {
                lower["phi"] <- 0
                upper["phi"] <- 0.99
            }
                

        obj$env$tracepar <- TRUE

        system.time( opt <- nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper) )
        Pvalue <- NA

        if(fit0)
            {
                parameters0 <- parameters
                parameters0$logsdGearRW <- -10
                
                obj0 <- MakeADFun(
                    data = data,
                    parameters = parameters0,
                    map = c(map,list(logsdGearRW=factor(NA))),
                    random = random,
                    DLL="gearcalib"
                    )

                obj0$env$tracepar <- TRUE
                
                system.time( opt0 <- nlminb(obj0$par,obj0$fn,obj0$gr,lower=lower,upper=upper) )

                Pvalue <- pchisq(2*(opt0$objective-opt$objective),df=1)

                print(paste("Chisq-test of no size structure in gear effect:",Pvalue))

                rm(obj0)
            }

        rep <- sdreport(obj)
        repsum <- summary(rep,"random")
        s <- repsum[grep("gear",rownames(repsum)),]
        est <- 2*s[,1]
        sd <- 2*s[,2]

        ret <- list(d=d,
                    Pvalue=Pvalue,rep=rep,opt=opt,obj=obj,est=est,sd=sd)

        class(ret)<-"gearcalibFit" 

        return(ret)

    }
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Empirical estimate with bootstrap confidence regions
##' @param d gearcalib data
##' @param quantiles report these quantiles
##' @param Nboot number of bootstrap replicates
##' @return a list
boot <- function(d,quantiles = c(0.025,0.10),Nboot=1000){
  GearNames <- levels(d$Gear)

  ## Raw data: Total catches
  I1 <- d$Gear==GearNames[1]
  I2 <- d$Gear==GearNames[2]

  totcatch1 <- apply(d$N[I1,],2,sum)
  totcatch2 <- apply(d$N[I2,],2,sum)

  # Density1 <- totcatch1 / sum(d$SweptArea[I1])
  # Density2 <- totcatch2 / sum(d$SweptArea[I2])
  Density1 <- totcatch1 / apply(d$SweptArea[I1,],2,sum)
  Density2 <- totcatch2 / apply(d$SweptArea[I2,],2,sum)

  tiny <- 1e-6
  
  RawEstimate <- (Density2 + tiny) / (Density1 + tiny)

  BootEstimate <- array(NA,c(Nboot,length(RawEstimate)))

  for(i in 1:Nboot)
      {
          ## Select random groups with replacement
          GroupBoot <- sample(1:length(levels(d$group)),length(levels(d$group)),replace=TRUE)
          
          Haul1 <- as.numeric(sapply(GroupBoot,function(g) which( (as.numeric(d$group) == g) & I1)))
          Haul2 <- as.numeric(sapply(GroupBoot,function(g) which( (as.numeric(d$group) == g) & I2)))

          Nb1 <- d$N[Haul1,]
          Nb2 <- d$N[Haul2,]

          # Db1 <- apply(Nb1,2,sum) / sum(d$SweptArea[Haul1])
          # Db2 <- apply(Nb2,2,sum) / sum(d$SweptArea[Haul2])
          Db1 <- apply(Nb1,2,sum) / apply(d$SweptArea[Haul1,],2,sum)
          Db2 <- apply(Nb2,2,sum) / apply(d$SweptArea[Haul2,],2,sum)

          BootEstimate[i,] <- (Db2 + tiny) / (Db1 + tiny) 
      }

  lst = list(); i=0;
  for (q in quantiles){
    lwr = apply(BootEstimate,2,quantile,probs=q);
    upr = apply(BootEstimate,2,quantile,probs=1-q);
    lst[[i<-i+1]] = tibble::tibble(ci=1-2*q,L=d$L,lwr=lwr,upr=upr);
  }
  BootQuantiles <- dplyr::bind_rows(lst);
  medBoot <- tibble::tibble(L=d$L,val=apply(BootEstimate,2,quantile,probs=0.5));

  return(list(RawEstimate=RawEstimate,
              BootEstimate=BootEstimate,
              quantiles=quantiles,
              medBoot=medBoot,
              BootQuantiles=BootQuantiles,
              Density1=Density1,
              Density2=Density2))

}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plots a gear calibration fit
##' @param fit an onject of class 'gearcalibFit' as fitted by gearcalib()
##' @param select a vector specifying which plots are wanted. Can be "relsel" for "Relative selectivity" or "density".
##' @param boot (optional) list with bootstrap estimates as produced by boot()
##' @param Lvec (optional) vector with labels for each length group
##' @return nothing
plot.gearcalibFit <- function(fit,select=c("relsel","density"),boot=NULL, Lvec=NULL, ymax=NULL){
  
  if(is.null(Lvec)) Lvec <- 1:(ncol(fit$d$N))
  Lmin <- min(Lvec)
  Lmax <- max(Lvec)
  
  if("relsel" %in% select){
       # plot(range(Lvec),c(0,max(exp(est + sd %o% c(0,-2,2)))),type="n",
       #      ylim=c(0,2),
       #      xlim=range(Lvec),
       #      xlab="Length group",
       #      ylab=paste(levels(fit$d$Gear)[2]," vs. ",levels(fit$d$Gear)[1]),
       #      main="Relative selectivity")
       # 
       # lines(range(Lvec),rep(1,2),col="darkgrey",lwd=3,lty="dashed")
       # 
       # polygon(c(Lvec,rev(Lvec)),c(exp(est-2*sd),rev(exp(est+2*sd))),
       #         col="grey",border=NA)
       # lines(Lvec,exp(est),lwd=3)
       # if(!is.null(boot)){
       #     points(Lvec,b$RawEstimate)
       #     apply(b$BootQuantiles,1,function(x)lines(Lvec,x))
       # }
       # grid()
    est = fit$est; sd = fit$sd;
    dfr = tibble::tibble(x=Lvec,
                         y=exp(est),
                         ymin=exp(est-2*sd),
                         ymax=exp(est+2*sd));
    if (!is.null(boot)) dfr$RawEstimate = boot$RawEstimate;
    p1 = ggplot(dfr,aes(x=x,y=y,ymin=ymin,ymax=ymax));
    if (!is.null(boot)){
      cis = sort(unique(boot$BootQuantiles$ci),decreasing=TRUE);
      for (ci_ in cis) {
        p1 = p1 + geom_ribbon(mapping=aes(x=L,ymin=lwr,ymax=upr),
                              data=boot$BootQuantiles |> dplyr::filter(ci==ci_),
                              alpha=0.2,fill="green",colour=NA,inherit.aes=FALSE);
      }
    }
    p1 = p1 + 
           geom_ribbon(colour=NA,alpha=0.2) + 
           geom_line() + 
           geom_hline(yintercept=c(0,1),linetype=c(1,3),alpha=c(1,0.5)) + 
           labs(x="size (mm CW)",y="relative selectivity") + 
           wtsPlots::getStdTheme();
    if (!is.null(boot)){
      if (is.null(ymax)){
        ymax = max(1,dfr$ymax);
        ymx1 = max(1,dfr$ymax,dfr$RawEstimate,boot$BootQuantiles$upr);
        if (ymx1>3*ymax) ymax = 3*ymax;
      }
      p1 = p1 + 
           geom_point(aes(y=RawEstimate),shape=21) + 
           geom_line(aes(x=L,y=val),data=boot$medBoot,linetype=2,colour="green",inherit.aes=FALSE);
    }
    p1 = p1 + scale_y_continuous(limits=c(0,ymax),oob=scales::oob_squish_any);
    print(p1)
  }
       

  if(!is.null(boot) && ("density" %in% select)){
  #     with(boot,{
  #     plot(Lvec,log10(1+Density1),
  #           ylim=log10(1+range(c(Density1,Density2))),
  #           xlim=c(Lmin,Lmax),
  #           type="l",lty="dashed",
  #           xlab="Length [cm]",ylab="Density (log10(N/A+1))")
  #      points(Lvec,log10(1+Density1),pch="o")
  #      lines(Lvec,log10(1+Density2))
  #      points(Lvec,log10(1+Density2),pch="+")
  #      legend("topright",legend=GearNames,lty=c("dashed","solid"),pch=c("o","+"))
  #      
  #      grid()
  #     })
    dfr = tibble::tibble(L=Lvec,
                         Density1=boot$Density1,
                         Density2=boot$Density2,
                         LD1=log10(1+Density1),
                         LD2=log10(1+Density2));
    gcs = c("blue","green");        #--gear colours
    names(gcs) = levels(fit$d$Gear);
    p2 = ggplot(dfr,aes(x=L)) + 
           geom_point(aes(y=LD1,colour=names(gcs)[1]),shape=21) + 
           geom_line( aes(y=LD1,colour=names(gcs)[1])) + 
           geom_point(aes(y=LD2,colour=names(gcs)[2]),shape=22) + 
           geom_line( aes(y=LD2,colour=names(gcs)[2])) + 
           labs(x="size (mm CW)",y="Density (log10(N/A+1))") + 
           scale_color_manual(name="gear",
                              breaks=names(gcs),
                              values=gcs) +
           wtsPlots::getStdTheme();
    print(p2);
  }
  
  return(list(p1=p1,p2=p2));
}
