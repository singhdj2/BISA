#######################################################################
### objective: create functions for phenotypic data analysis       ####
### Functions: heritability, BLUE/BLUP, mixture normal,            ####
### .........  DEM slicing, genetic correlations                   ####
### author: Daljit Singh (singhdj2@ksu.edu)                        ####
#######################################################################
library(lme4)
library(dplyr)
library(tidyr)
library(mixtools)
library(RMySQL)
library(sommer)
#library(ggplot2)

### heritablity function for single trait based broad-sense heritability calculations
calcH2func <- function(dat) {
  v = data.frame(VarCorr(lmer(value ~ 0 + (1|gid) + rep + (1|rep:block), data=dat)))
  cbind(H2=round(v[1,4]/(v[1,4]+(v[3,4]/2)),2), vG=round(v[1,4],2), vE=round(v[3,4],2))
}
# end heritability function

### Best Linear Unbiased Estimates (BLUE) function (with trial main effect)
calcBlueFunc <- function(data) {
  if(!is.character(data$rep) | !is.character(data$block) | !is.character(data$trial) | !is.character(data$gid)) {
    data$rep = as.character(data$rep)
    data$range = as.character(data$range)
    data$column = as.character(data$trial)
    data$plot_name = as.character(data$gid)
    #data$block = as.character(data$block)
  }
  mod <- try(lmer(value ~ 0 + gid + (1|trial) + (1|trial:rep) + (1|trial:rep:block), data=data))
  blue <- fixef(mod)
  entry <- names(fixef(mod))
  data.frame(entry, blue) #outputs dataframe with blues and entry names
}
# end BLUE function

### Function to calculate univariate normal mixture distributions 
### To extract digital elevation model based lodging ...
# updated function for better error handling and added iterations to 
# ..account for random start values for some troubling plots (non)
mixCalc <- function(x){
  w <- try(suppressWarnings(normalmixEM(x$dem, lambda = c(0.5,0.5), 
                                        mu= c(-60,10), mean.constr = c(NA,0),
                                        maxit = 2000, verb = F, maxrestarts=60)), silent = T)
  if(class(w)=="try-error"){
    cat("Note: Error in model run. \n")
  } else{
    for(i in 1:5){
      z <- try(suppressWarnings(normalmixEM(x$dem, lambda = c(0.5,0.5), 
                                            mu= c(-60,10), mean.constr = c(NA,0),
                                            maxit = 2000, verb = F, maxrestarts=60)), silent = T)
      if (i == 1){ #discard first
        next
      } else {
        if (class(z) =="try-error"){
          next
        } else {
          df <- cbind(DL_nIter=length(w$all.loglik)-1,DL_nRestarts=w$restarts,
                      DL_lamb1=w$lambda[1], DL_lamb2=w$lambda[2], DL_mu1=w$mu[1],
                      DL_mu2=w$mu[2], DL_sig1=w$sigma[1],DL_sig2=w$sigma[2],DL_loglik=w$loglik)
          #cat(df,i,"\n")
        }
      }
    }
    df %>% data.frame() %>% filter(DL_lambda1==max(DL_lambda1))   #summarise_each(funs(max))
  }
}
# end mixnorm function


### Function to slice DEM by percentiles/quantile
## quantile function
qnt <- function(x,a) quantile(x, probs=a, na.rm=TRUE)
ldh.0315minus0302.osav <- ldh.0315minus0302.osav %>% dplyr::rename(dem=height_diff_cm) %>% 
              group_by(plot_id) %>% 
              summarize(num.pix=length(dem),
                        top0.05=mean(dem[dem > qnt(dem,0.95)], na.rm=T),
                        top0.10=mean(dem[dem > qnt(dem,0.90)], na.rm =T),
                        top0.15=mean(dem[dem > qnt(dem,0.85)], na.rm =T),
                        top0.20=mean(dem[dem > qnt(dem,0.80)], na.rm =T),
                        top0.80=mean(dem[dem > qnt(dem,0.20)], na.rm =T),
                        top0.90=mean(dem[dem > qnt(dem,0.10)], na.rm =T),
                        mean=mean(dem, na.rm =T),
                        median=median(dem,na.rm=T),
                        var=var(dem),
                        mix.param=mixCalc(dem)
              ) %>% 
              separate(mix.param, into = c("lamb1","lamb2","mu1","mu2","sig1","sig2","loglik"), sep=",")
# end percentile function 


### example usage of mixnormal function (input is a DEM dataframe with plot attributes)
## run mix normal function on full frame (can parallelize the function for larger sets)
## note: takes approx 30 minutes for >1700 plots with 3-4k data points each)
system.time(DL.mix <- dem.clean %>%
              group_by(plot_id) %>% do(data.frame(mixCalc(.))) %>% data.frame())

### example usage of broad-sense heritability looping through each trait
H2.all.mixLI <- tst %>% dplyr::select(plot_id:DL_mean,DL_mixLI) %>% 
  gather(trait,value,DAYSMT:DL_mixLI) %>% mutate(trial=as.character(trial)) %>% 
  group_by(trait,trial) %>% do(data.frame(calcH2func(.)))

## example usage of BLUE function
pheno.blue <- pheno.mix.new %>% gather(trait,value,DAYSMT:DL_mixLI) %>% 
  group_by(trait) %>% do(calcBlueFunc(.)) %>% data.frame() %>%
  mutate(entry=as.character(entry)) %>% 
  mutate(GID=paste0("GID",substr(entry,4,nchar(entry)))) %>%
  dplyr::select(-entry)

##################################################################
#### genetic correlations across traits using sommer package #####
##################################################################
## Fit Multivariate Linear mixed model on pheno traits
# check genotype order in two dfs
all(names(data.frame(geno01))==blue.rrblup$GID)
# create marker matrix for sommer input
A <- A.mat(t(geno01))
# run Multi-variate model
system.time(ans.A <- mmer2(cbind(LOI,LOS,PH,DTHD,GRYLD,DAYSMT,DL_mean,DL_mixLI)~1, random=~us(trait):g(GID), rcov=~us(trait):units, G=list(GID=A),
               data=blue.rrblup, silent = FALSE))
summary(ans.A)
## genetic variance covariance
gvc <- ans.A$var.comp$`g(GID)`
## extract variances (diagonals) and get standard deviations 
sd.gvc <- as.matrix(sqrt(diag(gvc)))
## get possible products sd(Vgi) * sd(Vgi')
prod.sd <- sd.gvc %*% t(sd.gvc)
## genetic correlations cov(gi,gi')/[sd(Vgi) * sd(Vgi')] 
gen.cor <- gvc/prod.sd
gen.cor

## genetic corr with standard errors
pin(ans.A, gen.corr ~ V2 / sqrt(V1*V9)) # LOI vs LOS
pt(-abs(0.92/0.0235),df=pmin(589,589)-1) # calculate p-values: (estimate/std error), degrees of freedom
pin(ans.A, gen.corr ~ V3 / sqrt(V1*V16)) # LOI vs PH
pin(ans.A, gen.corr ~ V10 / sqrt(V9*V16)) # LOS vs PH
pt(-abs(0.22/0.133),df=pmin(589,589)-1)

## narrow-sense heritability
(h2 <- diag(gvc) / diag(cov(blue.rrblup[,names(diag(gvc))], use = "complete.obs")))
sqrt(h2)

########################################################
#### fetching data from wheatgenetics mySQL database ###
########################################################
cimm= dbConnect(MySQL(),user="user13", dbname='cimmyt', host="apate", password="xpwxxxx") #when on server use 'apate'
#dbListTables(cimm); dbListFields(cimm,'phenotypes') ; dbListFields(cimm,'traits')
# single query combining pheno and plot infromation
pheno <- dbGetQuery(cimm, "SELECT  phenotypes.plot_id,
	                             phenotypes.phenotype_value,
                               phenotypes.phenotype_date,
                               plots.iyear,
                               plots.ilocation,
                               plots.icondition,
                               plots.plot_no,
                               #plots.plot_id,
                               plots.entry,
                               plots.gid,
                               plots.rep,
                               plots.block,
                               plots.col,
                               plots.`row`,
                               plots.trial,
                               plots.location,
                               plots.planting_date
                               FROM plots
                               LEFT JOIN phenotypes ON plots.plot_id = phenotypes.plot_id
                               #AND gbs.plexing LIKE barcodes.`set`
                           WHERE (phenotypes.plot_id LIKE '17-FAS-STN%' OR phenotypes.plot_id LIKE '16-FAS-STN%' OR 
                                phenotypes.plot_id LIKE '15-FAS-STN%' OR phenotypes.plot_id LIKE '14-FAS-STN%') 
                                 AND (trait_id = 'LOI' AND phenotype_date !='2017-03-29')
")

#close DB connection
dbDisconnect(cimm)
##
head(pheno)
pheno <- pheno %>% mutate(loc_date=paste(.$iyear,.$ilocation,.$phenotype_date,sep="_"))

## add genomic prediction functions below....
