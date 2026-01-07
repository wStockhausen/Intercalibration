#--fitmodel code for gearcalib0.R
fitmodel <- function(Species,GearNames,d)
    {
        ## Construct table of which gears are used in which groups
        tt <- table(d$group,d$Gear)

        ## Identify those groups where the two gears are used
        UseGroups <- tt[,GearNames[1]] & tt[,GearNames[2]] 

        I <- sapply(d$group,function(g)UseGroups[g])

        ## Use only if SweptArea is finite
        I <- I & is.finite(d$SweptArea)
        
        ## Select only these records
        d$SweptArea <- d$SweptArea[I]
        d$N         <- d$N[I,Lmin:Lmax]
        d$group     <- factor(d$group[I],levels=unique(d$group[I]))
        d$Gear      <- factor(d$Gear[I],levels=GearNames)

        nsize <- ncol(d$N)
        ngroup <- nlevels(d$group)
        nhaul <- nrow(d$N)
        ngear <- nlevels(d$Gear)

        obj <- MakeADFun(
            data=list(
                N=d$N,
                SweptArea=d$SweptArea,
                group=d$group,
                Gear=d$Gear,
                huge=10,
                tiny=0.01
            ),
            parameters=list(
                logspectrum=matrix(0,ngroup,nsize),
                ## nugget=matrix(0,nhaul,nsize),     #<== not estimating nugget
                residual=matrix(0,nhaul,nsize),
                loggear=numeric(nsize),
                logsd=-1,
                phi=.5,
                ## logsdnug=-1,                      #<== not estimating nugget
                logsdres=-1,
                logsdGearRW=-1,
                alpha = 1
            ),
            map = list(alpha=factor(NA)), # phi=factor(NA),logsdnug=factor(NA),logsdGearRW=factor(NA)),
            random=c("logspectrum","residual","loggear"),
            DLL="gearcalib0"
        )

        obj0 <- MakeADFun(
            data=list(
                N=d$N,
                SweptArea=d$SweptArea,
                group=d$group,
                Gear=d$Gear,
                huge=10,
                tiny=0.01
            ),
            parameters=list(
                logspectrum=matrix(0,ngroup,nsize),
                ## nugget=matrix(0,nhaul,nsize),
                residual=matrix(0,nhaul,nsize),
                loggear=numeric(nsize),
                logsd=-1,
                phi=.5,
                ## logsdnug=-1,
                logsdres=-1,
                logsdGearRW=-10,
                alpha = 1
            ),
            ## map = list(phi=factor(NA),logsdnug=factor(NA),logsdGearRW=factor(NA)),    
            map = list(alpha=factor(NA),logsdGearRW=factor(NA)),    
            random=c("logspectrum","residual","loggear"),
            DLL="gearcalib0"
        )

        lower <- 0*obj$par-Inf
        upper <- 0*obj$par+Inf
        lower["phi"] <- 0
        upper["phi"] <- 0.9999

        system.time( opt  <- nlminb(obj$par, obj$fn, obj$gr, lower=lower,upper=upper) )
        system.time( opt0 <- nlminb(obj0$par,obj0$fn,obj0$gr,lower=lower,upper=upper) )

        Pvalue <- pchisq(2*(opt0$objective-opt$objective),df=1)

        print(paste("Chisq-test of no size structure in gear effect:",Pvalue))

        rm(obj0)
        save(list=ls(), file=paste("Image-",Species,"-",GearNames[2],"-vs-",GearNames[1],".Rdata",sep=""))

        rep <- sdreport(obj)

        repsum <- summary(rep)
        s   <- repsum[grep("gear",rownames(repsum)),]
        est <- 2*s[,1]
        sd  <- 2*s[,2]

        # GearNames <- levels(d$Gear)
        # 
        # ## Raw data: Total catches
        # I1 <- d$Gear==GearNames[1]
        # I2 <- d$Gear==GearNames[2]
        # 
        # totcatch1 <- apply(d$N[I1,],2,sum)
        # totcatch2 <- apply(d$N[I2,],2,sum)
        # 
        # Density1 <- totcatch1 / sum(d$SweptArea[I1])
        # Density2 <- totcatch2 / sum(d$SweptArea[I2])

        ## Plots
        # X11()
        # 
        # plot(c(Lmin,Lmax),c(0,max(exp(est + sd %o% c(0,-2,2)))),type="n",
        #      ylim=c(0,2),
        #      xlim=c(Lmin,Lmax),
        #      xlab="Length [cm]",
        #      ylab=paste(GearNames[2]," vs. ",GearNames[1]),
        #      main=paste("Relative selectivity for", Species))
        # 
        # lines(c(Lmin,Lmax),rep(1,2),col="darkgrey",lwd=3,lty="dashed")
        # 
        # polygon(c(Lmin:Lmax,Lmax:Lmin),c(exp(est-2*sd),rev(exp(est+2*sd))),
        #         col="grey",border=NA)
        # lines(Lmin:Lmax,exp(est),lwd=3)
        # 
        # points(Lmin:Lmax,Density2/Density1)
        # grid()
        # 
        # dev.copy2pdf(file=paste("Relsel-",Species,"-",GearNames[2],"-vs-",GearNames[1],".pdf",sep=""))
        # 
        # plot(Lmin:Lmax,log10(1+Density1),
        #      ylim=log10(1+range(c(Density1,Density2))),
        #      xlim=c(Lmin,Lmax),
        #      type="l",lty="dashed",
        #      xlab="Length [cm]",ylab="Density (log10(N/A+1))")
        # points(Lmin:Lmax,log10(1+Density1),pch="o")
        # lines(Lmin:Lmax,log10(1+Density2))
        # points(Lmin:Lmax,log10(1+Density2),pch="+")
        # legend("topright",legend=GearNames,lty=c("dashed","solid"),pch=c("o","+"))
        # 
        # grid()
        # 
        # dev.copy2pdf(file=paste("Density-",Species,"-",GearNames[2],"-vs-",GearNames[1],".pdf",sep=""))
        # 
        ## #############################################################################################
        ## Plot example of a size spectrum
        # pl <- obj$env$parList()
        # spec <- pl$logspectrum
        # 
        # ## Take the station witht the largest catch
        # ii <- which.max(apply(spec,1,sum))
        # 
        # iii <- which(d$group == levels(d$group)[ii])
        # 
        # res1 <- pl$residual[iii,]
        # 
        # plot(Lmin:Lmax,exp(spec[ii,])*d$SweptArea[ii],log="y",type="l",ylim=c(1e-1,max(d$N[iii,])),
        #      xlim=c(Lmin,Lmax),
        #      xlab="Length [cm]",ylab="Size spectrum of population [#/cm]",lwd=3,
        #      main=paste(Species,"at station",d$group[ii]))
        # grid()
        # 
        # lines(Lmin:Lmax,exp(spec[ii,] + res1[1,])*d$SweptArea[ii],lty="solid")
        # lines(Lmin:Lmax,exp(spec[ii,] + res1[2,])*d$SweptArea[ii],lty="dashed")
        # 
        # points(Lmin:Lmax,d$N[iii[1],],pch="+")
        # points(Lmin:Lmax,d$N[iii[2],],pch="o")
        # 
        # legend("bottom",legend=c("Size structure at station",
        #                     paste(GearNames[1],", observed and modelled"),
        #                     paste(GearNames[2],", observed and modelled")),
        #        lwd=c(3,1,1),lty=c("solid","solid","dashed"),pch=c(NA,"+","o"))
        # 
        # grid()
        # 
        # dev.copy2pdf(file=paste("Sample-",Species,"-",GearNames[2],"-vs-",GearNames[1],".pdf",sep=""))
        # 
        res <- data.frame(L=Lmin:Lmax,est=est,sd=sd)

        # write.table(res,file=paste("Results-",Species,"-",GearNames[2],"-vs-",GearNames[1],".csv",sep=""))
        # write.table(res,file=paste("Estimates-",Species,"-",GearNames[2],"-vs-",GearNames[1],".csv",sep=""))

        return(list(Species=Species,GearNames=GearNames,d=d,res=res,Pvalue=Pvalue,rep=rep,opt=opt,obj=obj))
    }

plot_gearcalib0<-function(lst){
  require(ggplot2);
  #--model estimates
  dfrRes = lst$res;
  #--data
  d = lst$d;
  GearNames <- levels(d$Gear);
  ## Raw data: Total catches
  I1 <- d$Gear==GearNames[1];
  I2 <- d$Gear==GearNames[2];
  totcatch1 <- apply(d$N[I1,],2,sum);
  totcatch2 <- apply(d$N[I2,],2,sum);
  Density1 <- totcatch1 / sum(d$SweptArea[I1]);
  Density2 <- totcatch2 / sum(d$SweptArea[I2]);
  dfrData = tibble::tibble(L=min(dfrRes$L):max(dfrRes$L),
                           Density1=Density1,
                           Density2=Density2,
                           obsR=Density2/Density1);

        # plot(c(Lmin,Lmax),c(0,max(exp(est + sd %o% c(0,-2,2)))),type="n",
        #      ylim=c(0,2),
        #      xlim=c(Lmin,Lmax),
        #      xlab="Length [cm]",
        #      ylab=paste(GearNames[2]," vs. ",GearNames[1]),
        #      main=paste("Relative selectivity for", Species))
        # 
        # lines(c(Lmin,Lmax),rep(1,2),col="darkgrey",lwd=3,lty="dashed")
        # polygon(c(Lmin:Lmax,Lmax:Lmin),c(exp(est-2*sd),rev(exp(est+2*sd))),
        #         col="grey",border=NA)
        # lines(Lmin:Lmax,exp(est),lwd=3)
        # points(Lmin:Lmax,Density2/Density1)
  p1 = ggplot(dfrRes,aes(x=L,y=exp(est),ymin=exp(est-2*sd),ymax=exp(est+2*sd))) + 
        geom_ribbon(colour=NA,alpha=0.2) + 
        geom_point(aes(x=L,y=obsR),data=dfrData,inherit.aes=FALSE) + 
        geom_line() + 
        geom_hline(yintercept=c(0,1),linetype=3,alpha=0.8) + 
        labs(x="Length (cm)",y="estimated relative selctivity",
          subtitle=paste0(lst$Species,": ",GearNames[2],"-",GearNames[1])) + 
        wtsPlots::getStdTheme();
  print(p1);
  
        # plot(Lmin:Lmax,log10(1+Density1),
        #      ylim=log10(1+range(c(Density1,Density2))),
        #      xlim=c(Lmin,Lmax),
        #      type="l",lty="dashed",
        #      xlab="Length [cm]",ylab="Density (log10(N/A+1))")
        # points(Lmin:Lmax,log10(1+Density1),pch="o")
        # lines(Lmin:Lmax,log10(1+Density2))
        # points(Lmin:Lmax,log10(1+Density2),pch="+")
        # legend("topright",legend=GearNames,lty=c("dashed","solid"),pch=c("o","+"))
  p2 = ggplot(dfrData,aes(x=L)) + 
         geom_point(aes(y=log10(1+Density1)),shape=21) + 
         geom_point(aes(y=log10(1+Density2)),shape=22) + 
         geom_line(aes(y=log10(1+Density2))) + 
        labs(x="Length (cm)",y="Density (log10(N/A+1))",
          subtitle=paste0(lst$Species,": ",GearNames[2],"-",GearNames[1])) + 
        wtsPlots::getStdTheme();
  print(p2);
  
  ## Plot example of a size spectrum
  pl <- lst$obj$env$parList();
  spec <- pl$logspectrum;
  ### Take the station with the largest catch
  ii <- which.max(apply(spec,1,sum))
  iii <- which(d$group == levels(d$group)[ii])
  res1 <- pl$residual[iii,]

      # plot(Lmin:Lmax,exp(spec[ii,])*d$SweptArea[ii],log="y",type="l",ylim=c(1e-1,max(d$N[iii,])),
      #      xlim=c(Lmin,Lmax),
      #      xlab="Length [cm]",ylab="Size spectrum of population [#/cm]",lwd=3,
      #      main=paste(Species,"at station",d$group[ii]))
      # grid()
      # 
      # lines(Lmin:Lmax,exp(spec[ii,] + res1[1,])*d$SweptArea[ii],lty="solid")
      # lines(Lmin:Lmax,exp(spec[ii,] + res1[2,])*d$SweptArea[ii],lty="dashed")
      # 
      # points(Lmin:Lmax,d$N[iii[1],],pch="+")
      # points(Lmin:Lmax,d$N[iii[2],],pch="o")
      # 
      # legend("bottom",legend=c("Size structure at station",
      #                     paste(GearNames[1],", observed and modelled"),
      #                     paste(GearNames[2],", observed and modelled")),
      #        lwd=c(3,1,1),lty=c("solid","solid","dashed"),pch=c(NA,"+","o"))
  dfrSpec <- tibble::tibble(L=min(dfrRes$L):max(dfrRes$L),
                            y=exp(spec[ii,])*d$SweptArea[ii],
                            y1=exp(spec[ii,] + res1[1,])*d$SweptArea[ii],
                            y2=exp(spec[ii,] + res1[2,])*d$SweptArea[ii])
  p3 = ggplot(dfrSpec,aes(x=L,y=y)) + 
        geom_line() + 
        geom_point(aes(x=L,y=y1),shape=21) + 
        geom_point(aes(x=L,y=y2),shape=22) + 
        geom_line(aes(x=L,y=y1),colour="blue", linetype=3) + 
        geom_line(aes(x=L,y=y2),colour="green",linetype=3) + 
        labs(x="Length (cm)",y="Size spectrum of population [#/cm]",
            subtitle=paste(lst$Species,"at station",d$group[ii])) + 
        wtsPlots::getStdTheme();
  print(p3);    
  

  return(list(p1=p1,p2=p2,p3=p3))
}

