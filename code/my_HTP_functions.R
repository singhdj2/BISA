#######################################################################
### objective: create functions for phenotypic data analysis       ####
### Functions: heritability, BLUE/BLUP, mixture Normal,DEM slicing ####
#######################################################################
library(lme4)
library(dplyr)
library(tidyr)
library(mixtools)

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

