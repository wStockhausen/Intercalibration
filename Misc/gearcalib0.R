library(TMB)
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/04a_ThygesenIntercalibration/Intercalibration/Misc");

#--compile cpp file
setwd(dirThs);
compile("gearcalib0.cpp",CXX="clang++")
if(is.loaded("EvalDoubleFunObject")) dyn.unload("gearcalib0.so") #<-may error out
dyn.load("gearcalib0.so")

## Specify range of length classes (cm) to include in the model
Lmin <- 10
Lmax <- 80

Lvec <- seq(Lmin,Lmax)

Do.Intercal <- TRUE
Do.Verification <- TRUE

source(file.path(dirThs,"fitmodel_for_gearcalib0.R"));

fits <- list()

if(Do.Intercal){
        Dir <- "Repaired-Trawl-Pairs"
        filename <- file.path(dirThs,"..",Dir,"InterCal_Merluccius.RData")
        load(filename,envir=e1);

        for(Species in c("Capensis","Paradoxus"))
            for(gear in c("Afr_New","Afr_Old"))
                {
                    #--testing: Species = "Capensis"; gear = "Afr_New";
                    d <- eval(parse(text=paste0("e1$listIntercalData_",Species)))

                    d$group <- factor(d$group)

                    ## Change order of gear such that Gisund is reference
                    d$Gear <- factor(d$Gear,levels=rev(unique(d$Gear)))

                    d$SweptArea <- d$SweptArea[,1]


                    ## Compare these two gear types
                    GearNames <- c("Gisund",gear)

                    lst  <- fitmodel(Species,GearNames,d);
                    lst$plts <-plot_gearcalib0(lst);
                    fits[[paste(Species,gear)]] = lst;
                    rm(lst);
            }
  if (FALSE){
    #--re-do plots (without refitting models)
    for (case in names(fits)){
      lst = fits[[case]];
      lst$plts = plot_gearcalib0(lst);
      fits[[case]] = lst;
    }
  }
}

if(Do.Verification)
    {
        
        Dir <- "Verification"
        filename <- paste(Dir,"InterCal_Merluccius_1998_1999.RData",sep="/")
        load(filename)


        for(Species in c("Capensis","Paradoxus"))
                {
                    d <- eval(parse(text=paste("listIntercalData_",Species,sep="")))

                    ## Pick only data from the 1990's
                    I <- grep("199",d$group)
                    d$Gear <- d$Gear[I]
                    d$N <- d$N[I,]
                    d$group <- d$group[I]
                    d$SweptArea <- d$SweptArea[I,]
                    
                    d$Gear[d$Gear == "Blue Sea"] <- "Blue_Sea"
                    d$group <- factor(d$group)

                    ## Change order of gear such that Gisund is reference
                    d$Gear <- factor(d$Gear,levels=rev(unique(d$Gear)))

                    d$SweptArea <- d$SweptArea[,1]


                    ## Compare these two gear types
                    GearNames <- c("Nansen","Blue_Sea")

                    fits[[length(fits)+1]] <- fitmodel(Species,GearNames,d)
                }
    }

save(fits,file="fits.Rdata")

casetab <- lapply(fits,function(f)return(
    list(Species=f$Species,
      Gear1=f$GearNames[1],
      Gear2=f$GearNames[2])))

casetab <- do.call(rbind,casetab)
partab <- t(sapply(fits,function(f)f$opt$par))
sdtab <- t(sapply(fits,function(f)sqrt(diag(f$rep$cov.fixed))))

colnames(sdtab) <- paste(colnames(sdtab),".sd",sep="")
partab <- cbind(casetab,partab,sdtab)

