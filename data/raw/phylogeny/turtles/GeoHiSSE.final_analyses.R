#library(devtools)
#devtools::install_github("thej022214/hisse", ref="higeosse_postsub")

library(hisse)
library(diversitree)
library(parallel)

#set up data
setwd("~/Documents/manuscripts/totw/totw.v7/for_Dryad")
trees<-read.tree("GeoHISSEInput.tree") 

#read.tree("pruned_turtles.tree") -> phy
read.nexus.data("GeoHiSSEInput.turtles_coast1_both0_not2.final.nex") -> data
turtle.dat <- data.frame(attributes(data), unlist(data,use.names=FALSE))
as.matrix(turtle.dat)->dat

#set up output
df <- data.frame(iteration=integer(1),
                 IndexBestModel=integer(1),
                 mod1=double(1),
                 mod2=double(1),
                 mod3=double(1),
                 mod4=double(1),
                 mod5=double(1),
                 mod6=double(1),
                 mod7=double(1),
                 mod8=double(1),
                 mod9=double(1),
                 mod10=double(1),
                 stringsAsFactors=FALSE)

for(i in 1:1) {
phy<- trees
outname <- sprintf("Turtle_GeoHiSSE")
dirname <- paste0("./results_", outname)
dir.create(path = dirname )

## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                                    , separate.extirpation=FALSE)
mod1 <- GeoHiSSE(phy, data=dat, f=c(0.81,0.93,0.81), speciation=speciation
                      , extirpation=extirpation, hidden.areas=FALSE
                      , trans.rate=trans.rate, assume.cladogenetic=TRUE)
saveRDS(mod1, file=paste0(dirname, "/", outname, "_mod1.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 1 done.") ) , file=paste0(outname,".log") )

## Model 2 - Dispersal parameters vary only, no range-dependent diversification.
speciation <- c(1,1,1)
extirpation <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=TRUE)
mod2 <- try( GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81), speciation=speciation
                      , extirpation=extirpation, hidden.areas=FALSE
                      , trans.rate=trans.rate, assume.cladogenetic=TRUE) )
saveRDS(mod2, file=paste0(dirname, "/", outname, "_mod2.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 2 done.") ) , file=paste0(outname,".log"), append=TRUE )

## Model 3. GeoHiSSE model with 1 hidden area, no range-dependent diversification.
## Note below how parameters vary among hidden classes but are the same within each ## hidden class.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, make.null=TRUE
                                    , include.jumps=FALSE, separate.extirpation=FALSE)
mod3 <- GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81),
                 speciation=speciation, extirpation=extirpation,
                 hidden.areas=TRUE, trans.rate=trans.rate)
saveRDS(mod3, file=paste0(dirname, "/", outname, "_mod3.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 3 done.") ) , file=paste0(outname,".log"), append=TRUE )

## Model 4. GeoHiSSE model with 1 hidden area, no range-dependent diversification.
## Note below how parameters vary among hidden classes but are the same within each ## hidden class.
speciation <- c(1,1,1,2,2,2)
extirpation <- c(1,1,2,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, make.null=TRUE
                                    , include.jumps=FALSE, separate.extirpation=TRUE)
mod4 <- GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81),
                  speciation=speciation, extirpation=extirpation,
                  hidden.areas=TRUE, trans.rate=trans.rate)
saveRDS(mod4, file=paste0(dirname, "/", outname, "_mod4.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 4 done.") ) , file=paste0(outname,".log"), append=TRUE )

## Model 5. Heterogeneous diversification, not tied to range evolution.
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=4, make.null=TRUE
                                    , include.jumps=FALSE, separate.extirpation=FALSE)
mod5 <- try( GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81), speciation=speciation
                      , extirpation=extirpation, hidden.areas=TRUE
                      , trans.rate=trans.rate, assume.cladogenetic=TRUE) )
saveRDS(mod5, file=paste0(dirname, "/", outname, "_mod5.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 5 done.") ) , file=paste0(outname,".log"), append=TRUE )

## Model 6. Heterogeneous diversification, not tied to range evolution. 
## Assumes 5 distinct diversification rates.
speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=4, make.null=TRUE
                                    , include.jumps=FALSE, separate.extirpation=TRUE)
mod6 <- try( GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81), speciation=speciation
                       , extirpation=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=TRUE) )
saveRDS(mod6, file=paste0(dirname, "/", outname, "_mod6.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 6 done.") ) , file=paste0(outname,".log"), append=TRUE )

## Model 7. Canonical GeoSSE model, range effect on diversification
speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=FALSE)
mod7 <- try( GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81), speciation=speciation
                      , extirpation=extirpation, hidden.areas=FALSE
                      , trans.rate=trans.rate, assume.cladogenetic=TRUE) )
saveRDS(mod7, file=paste0(dirname, "/", outname, "_mod7.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 7 done.") ) , file=paste0(outname,".log"), append=TRUE )


## Model 8. GeoSSE model, with range effect on diversification
speciation <- c(1,2,3)
extirpation <- c(1,2)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, include.jumps=FALSE
                                    , separate.extirpation=TRUE)
mod8 <- try( GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81), speciation=speciation
                      , extirpation=extirpation, hidden.areas=FALSE
                      , trans.rate=trans.rate, assume.cladogenetic=TRUE) )
saveRDS(mod8, file=paste0(dirname, "/", outname, "_mod8.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 8 done.") ) , file=paste0(outname,".log"), append=TRUE )

## Model 9. Heterogeneous diversification, tied to range evolution.
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=FALSE , separate.extirpation=FALSE)
mod9 <- try( GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81), speciation=speciation
                      , extirpation=extirpation, hidden.areas=TRUE
                      , trans.rate=trans.rate, assume.cladogenetic=TRUE) )
saveRDS(mod9, file=paste0(dirname, "/", outname, "_mod9.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 9 done.") ) , file=paste0(outname,".log"), append=TRUE )





## Model 10. Heterogeneous diversification, tied to range evolution.
## Assumes 6 distinct diversification rates.
speciation <- c(1,2,3,4,5,6)
extirpation <- c(1,2,3,4)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, include.jumps=FALSE
                                    , separate.extirpation=TRUE)
mod10 <- try( GeoHiSSE(phy, dat, f=c(0.81,0.93,0.81), speciation=speciation
                       , extirpation=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=TRUE) )
saveRDS(mod10, file=paste0(dirname, "/", outname, "_mod10.rds"))
#capture.output( print( paste0("Analysis - ", outname, " - model 10 done.") ) , file=paste0(outname,".log"), append=TRUE )






file=paste0(dirname, "/", outname, "_mod1.rds")

#use something like the following to reconstrcut states
mod1 <- readRDS(paste0(dirname, "/", outname, "_mod1.rds"))
mod2 <- readRDS(paste0(dirname, "/", outname, "_mod2.rds"))
mod3 <- readRDS(paste0(dirname, "/", outname, "_mod3.rds"))
mod4 <- readRDS(paste0(dirname, "/", outname, "_mod4.rds"))
mod5 <- readRDS(paste0(dirname, "/", outname, "_mod5.rds"))
mod6 <- readRDS(paste0(dirname, "/", outname, "_mod6.rds"))
mod7 <- readRDS(paste0(dirname, "/", outname, "_mod7.rds"))
mod8 <- readRDS(paste0(dirname, "/", outname, "_mod8.rds"))
mod9 <- readRDS(paste0(dirname, "/", outname, "_mod9.rds"))
mod10 <- readRDS(paste0(dirname, "/", outname, "_mod10.rds"))


recon.mod1 <- MarginReconGeoSSE(phy = mod1$phy, data = mod1$data, f = mod1$f,
                                pars = mod1$solution, hidden.areas = mod1$hidden.areas,
                                root.type = mod1$root.type, root.p = mod1$root.p,
                                aic = mod1$AIC, n.cores = 4)
                                
recon.mod2 <- MarginReconGeoSSE(phy = mod2$phy, data = mod2$data, f = mod2$f,
                                pars = mod2$solution, hidden.areas = mod2$hidden.areas,
                                root.type = mod2$root.type, root.p = mod2$root.p,
                                aic = mod2$AIC, n.cores = 4) 

recon.mod3 <- MarginReconGeoSSE(phy = mod3$phy, data = mod3$data, f = mod3$f,
                                 pars = mod3$solution, hidden.areas = mod3$hidden.areas,
                                 root.type = mod3$root.type, root.p = mod3$root.p,
                                 aic = mod3$AIC, n.cores = 4)      

recon.mod4 <- MarginReconGeoSSE(phy = mod4$phy, data = mod4$data, f = mod4$f,
                                 pars = mod4$solution, hidden.areas = mod4$hidden.areas,
                                 root.type = mod4$root.type, root.p = mod4$root.p,
                                 aic = mod4$AIC, n.cores = 4)    


recon.mod5 <- MarginReconGeoSSE(phy = mod5$phy, data = mod5$data, f = mod5$f,
                                pars = mod5$solution, hidden.areas = mod5$hidden.areas,
                                root.type = mod5$root.type, root.p = mod5$root.p,
                                aic = mod5$AIC, n.cores = 4)       

recon.mod6 <- MarginReconGeoSSE(phy = mod6$phy, data = mod6$data, f = mod6$f,
                                 pars = mod6$solution, hidden.areas = mod6$hidden.areas,
                                 root.type = mod6$root.type, root.p = mod6$root.p,
                                 aic = mod6$AIC, n.cores = 4)       

recon.mod7 <- MarginReconGeoSSE(phy = mod7$phy, data = mod7$data, f = mod7$f,
                                pars = mod7$solution, hidden.areas = mod7$hidden.areas,
                                root.type = mod7$root.type, root.p = mod7$root.p,
                                aic = mod7$AIC, n.cores = 4)
     
recon.mod8 <- MarginReconGeoSSE(phy = mod8$phy, data = mod8$data, f = mod8$f,
                                pars = mod8$solution, hidden.areas = mod8$hidden.areas,
                                root.type = mod8$root.type, root.p = mod8$root.p,
                                aic = mod8$AIC, n.cores = 4)       

recon.mod9 <- MarginReconGeoSSE(phy = mod9$phy, data = mod9$data, f = mod9$f,
                                pars = mod9$solution, hidden.areas = mod9$hidden.areas,
                                root.type = mod9$root.type, root.p = mod9$root.p,
                                aic = mod9$AIC, n.cores = 4)       

recon.mod10 <- MarginReconGeoSSE(phy = mod10$phy, data = mod10$data, f = mod10$f,
                                 pars = mod10$solution, hidden.areas = mod10$hidden.areas,
                                 root.type = mod10$root.type, root.p = mod10$root.p,
                                 aic = mod10$AIC, n.cores = 4)       


saveRDS(recon.mod1, file=paste0(dirname, "/", outname, "_recon_mod1.rds"))
saveRDS(recon.mod2, file=paste0(dirname, "/", outname, "_recon_mod2.rds"))
saveRDS(recon.mod3, file=paste0(dirname, "/", outname, "_recon_mod3.rds"))
saveRDS(recon.mod4, file=paste0(dirname, "/", outname, "_recon_mod4.rds"))
saveRDS(recon.mod5, file=paste0(dirname, "/", outname, "_recon_mod5.rds"))
saveRDS(recon.mod6, file=paste0(dirname, "/", outname, "_recon_mod6.rds"))
saveRDS(recon.mod7, file=paste0(dirname, "/", outname, "_recon_mod7.rds"))
saveRDS(recon.mod8, file=paste0(dirname, "/", outname, "_recon_mod8.rds"))
saveRDS(recon.mod9, file=paste0(dirname, "/", outname, "_recon_mod9.rds"))
saveRDS(recon.mod10, file=paste0(dirname, "/", outname, "_recon_mod10.rds"))



recon.mod1 <- readRDS(paste0(dirname, "/", outname, "_recon_mod1.rds"))
recon.mod2 <- readRDS(paste0(dirname, "/", outname, "_recon_mod2.rds"))
recon.mod3 <- readRDS(paste0(dirname, "/", outname, "_recon_mod3.rds"))
recon.mod4 <- readRDS(paste0(dirname, "/", outname, "_recon_mod4.rds"))
recon.mod5 <- readRDS(paste0(dirname, "/", outname, "_recon_mod5.rds"))
recon.mod6 <- readRDS(paste0(dirname, "/", outname, "_recon_mod6.rds"))
recon.mod7 <- readRDS(paste0(dirname, "/", outname, "_recon_mod7.rds"))
recon.mod8 <- readRDS(paste0(dirname, "/", outname, "_recon_mod8.rds"))
recon.mod9 <- readRDS(paste0(dirname, "/", outname, "_recon_mod9.rds"))
recon.mod10 <- readRDS(paste0(dirname, "/", outname, "_recon_mod10.rds"))


#Model set
recon.models<-list(recon.mod1, recon.mod2, recon.mod3, recon.mod4, recon.mod5, recon.mod6, recon.mod7, recon.mod8, recon.mod9,recon.mod10)

model.ave.rates <- GetModelAveRates(recon.models, type = "both", bound.par.matrix=cbind(c(0,0,0,0,0),c(10000000,10000000,10000000,10000000,10000000)) )

model.ave.rates.copy <- model.ave.rates

head(model.ave.rates.copy$tips)

#set up an extra colum using the 1,2,0 coding scheme
x<-1
for (i in model.ave.rates.copy$tips[,2])
{
  if (i == 1){model.ave.rates.copy$tips[x,10] <- "Coastal"} else {model.ave.rates.copy$tips[x,10] <- "NA"}
  x<-x+1
}
x<-1
for (i in model.ave.rates.copy$tips[,3])
{
  if (i == 1){model.ave.rates.copy$tips[x,10] <- "Inland"}
  x<-x+1
}
x<-1
for (i in model.ave.rates.copy$tips[,4])
{
  if (i == 1){model.ave.rates.copy$tips[x,10] <- "Widespread"}
  x<-x+1
}

#save this to file
pdf(paste0(dirname, "/", outname, "_boxplot_diversification.pdf"), width=2.25, height=1)
boxplot(model.ave.rates.copy$tips[,6] ~ model.ave.rates.copy$tips[,10], model.ave.rates.copy$tips, notch=TRUE , las = 1, ylab="Net Diversification", ylim=c(0,0.10), col=c("#00A600FF","#F2F2F2FF","#ECB176FF"))
dev.off()

pdf(paste0(dirname, "/", outname, "_boxplot_speciation.pdf"))
boxplot(model.ave.rates.copy$tips[,7] ~ model.ave.rates.copy$tips[,10], model.ave.rates.copy$tips, notch=FALSE , las = 1, ylab="Speciation", ylim=c(0,0.10), col=c("#00A600FF","#F2F2F2FF","#ECB176FF"))
dev.off()

pdf(paste0(dirname, "/", outname, "_boxplot_extinction.pdf"))
boxplot(model.ave.rates.copy$tips[,9] ~ model.ave.rates.copy$tips[,10], model.ave.rates.copy$tips, notch=FALSE , las = 1, ylab="Extinction", ylim=c(0,0.10), col=c("#00A600FF","#F2F2F2FF","#ECB176FF"))
dev.off()

pdf(paste0(dirname, "/", outname, "_boxplot_turnover.pdf"))
boxplot(model.ave.rates.copy$tips[,5] ~ model.ave.rates.copy$tips[,10], model.ave.rates.copy$tips, notch=FALSE , las = 1, ylab="Turnover", ylim=c(0,0.10), col=c("#00A600FF","#F2F2F2FF","#ECB176FF"))
dev.off()

mods<-list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10)
weights<-GetModelWeight(mods,criterion="AIC")


#write tree number followed by model weights
df[i,]<-c(i,which(weights==max(weights)),weights)
}
write.table(df,file="output.gplates.csv",sep=",",row.names=FALSE)

