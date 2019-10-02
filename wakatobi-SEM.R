# PATH ANALYSIS of Wakatobi SES
rm(list=ls())
library(piecewiseSEM)
library(nlme) # Version 3.1.122
library(ggplot2)
library(corrplot)
library(car) #VIF calculations
library(vegan)


# FIRST: create scatterplots

# Select all response + predictor variables + hierarchical variables
sitecols<-names(alldat.site) %in% c("site_id", "Site.Name", "lat_dd", "long_dd")
scatter.final<-alldat.site[,!sitecols]
loc.col<-grep("location", names(scatter.final))
response.col<-grep(fish.col, names(scatter.final))


scatter.names<-names(scatter.final)[-c(loc.col, response.col)]


########################################################################################
########################################################################################
########################################################################################
# Set graph names here: should match object scatter.names
scatter.titles<-c( # site journal columns
                  "Exposure", "Visibility", "Reef Type",  
                  
                  # Landings dat
                  "Landings per Trip", "Total Landings", 
                  "Total Personal", "Proportion Personal",
                  "Total Pengumpul", "Proportion Pengumpul", "Mean Pengumpul",
                  "Total Papalele", "Proportion Papalele",
                  "Total Market", "Proportion Market",
                  "Total On-Island", "Proportion On-Island", "Mean On-Island",
                  
                  # MSEC columns
                  "Mean NPP", "Min NPP", "Max NPP", "NPP SD", "Interannual NPP SD",
                  "Mean Wave Energy", "SD of Wave Energy", "Interannual SD of Wave Energy", "Wind Fetch", 
                  
                  # getWakatobiPop columns
                  "Population within 5km", "Fishers within 5km", "Rowboats within 5km", "Motorboats within 5km",
                  
                  # benthic columns
                  "Rugosity",
                  "Dead Coral With Algae", "Macroalgae", "Rubble", "Rock", "Soft Coral", "Sand", 
                  "All Hard Coral", "All Abiotic Substrate",
                  
                  # SST columns
                  "SST SD", "SST 50th Percentile", "SST 98th Percentile", "SST 2nd Percentile", "SST Kurtosis", "SST Skewness"
                  )

# Create scatterplots for fish.col
for(i in 1:length(scatter.names))
{
  newfile=paste("plot_scatter_", fish.col, "_vs_", scatter.names[i], ".pdf", sep="")
  p<-ggplot(data=scatter.final, aes(x=get(scatter.names[i]), y=get(fish.col))) + 
    geom_point(aes(x=get(scatter.names[i]), y=get(fish.col), color=location, shape=type_reef), size=2) +
    labs(y=fish.title, x=scatter.titles[i]) +  
    theme_light()+
    theme(text = element_text(size=18), 
          legend.position="right")
  
  pdf(file=newfile)
  print(p)      
  dev.off()
}

##### CALCULATE CORRELATIONS and VIF
##### NOTE: Correlations and VIF are unaffected by scaling (which preserves the shape and spread of data)
##### i.e., scaling can wait until the last step

##### Remove discreet variables
discreetvar<-match(c("wave_wind_fetch", "type_reef", "exposed"), scatter.names)
corr.names<-scatter.names[-discreetvar]
corr.final<-scatter.final[corr.names]


cor.dat<-cor(corr.final)
setwd(outdir)
write.csv(cor.dat, file="_Table_CorrelationsPearson.csv")
cor.test<-abs(cor.dat)>0.5
write.csv(cor.test, file="_Table_CorrelationsTestPearson.csv")

# Spearman's for all ranked, ordered vars (e.g., DACOR data)
cor.spear<-cor(corr.final, method="spearman")
write.csv(cor.spear, file="_Table_CorrelationsSpearman.csv")
cor.spear.test<-abs(cor.spear)>0.5
write.csv(cor.spear.test, file="_Table_CorrelationsTestSpearman.csv")

# Generate p values and confidence intervals for each correlation pair
pvals<-cor.mtest(corr.final, conf.level=0.95)

pdf(file="_Figure_CorrelationVisualPearson.pdf")
corrplot.mixed(cor.dat, upper="circle", lower="number", tl.pos="lt", tl.col="black", tl.cex=0.7, lower.col="black", addCoefasPercent=TRUE, number.cex=0.7, p.mat=pvals$p, sig.level=0.05, insig="blank", diag="n")
dev.off()

pdf(file="_Figure_CorrelationVisualSpearman.pdf")
corrplot.mixed(cor.spear, upper="circle", lower="number", tl.pos="lt", tl.col="black", tl.cex=0.7, lower.col="black", addCoefasPercent=TRUE, number.cex=0.7, p.mat=pvals$p, sig.level=0.05, insig="blank", diag="n")
dev.off()

##### QUESTION: do we have to worry about collinearity / multicollinearity in SEM?
# Supplementary table 6 here only tests correlation of exogenous variables (only those at the very start of a causation path): https://royalsocietypublishing.org/doi/suppl/10.1098/rspb.2016.0561#secSuppl
# Other studies, only address multicollinearity (but not sure which variables were tested): https://www.pnas.org/content/pnas/113/22/6230.full.pdf

##### ANSWER: as per Jonathan Lefcheck's email, both collinearity and multicollinearity are potential issues (but will only affect tests of significance, not parameter estimates)
## i.e., Both Rugosity and Hard Corals should not be used together as predictors of fish biomass

##### List of correlated variables:

## Suites of variables: Choose one of each
## All NPP variables with each other
## All wave variables with each other
## All human population variables with each other
## Some SST variables with each other: 
# SD with 98th and 2nd percentile and kurtosis
# 50th Percentile with 2nd percentile and skewness
# 98th Percentile with kurtosis and skewness
# 2nd percentile with kurtosis and skewness

## Other correlations:
## Visibility vs. 2nd percentile SST
## All wave variables vs. Rugosity, Soft Coral, Hard Coral, and some SST variables (SD, kurtosis, 2nd percentile)
## Rugosity vs. Hard Corals and SST_skew: 
## Soft corals vs Hard Corals and some SST variables (SD, 2nd and 98th percentile, and kurtosis)
## Hard coral vs. All abiotic benthos and SST_kurtosis
## Rubble with all human population variables: Don't use rubble

## Landings: TBD
# Currently, personal use vs. sold to market (papalele + pengumpul); Consider changing to: on-island (papalele + personal) vs off-island (pengumpul)

## STRATEGY: FIRST, Create path diagram, then trim based on correlations

### NOW, test for MULTICOLLINEARITY (calculate VIF)
### Seems like test for multicollinearity should be constructed for EACH linear model found within the SEM

#### Define PSEM equations here
#### Then call them later again for actual PSEM function

### REMINDER: for biomass, response.col indicates columns for raw and logged response variable
analysis.col<-grep(fish.col, names(alldat.site))

form1<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + landings_mean_onisland + Population_2017", sep=""))
form2<-as.formula("landings_mean_onisland ~ Population_2017")
form3<-as.formula("All_HardCoral ~ Population_2017")




#form3<-as.formula("Rugosity ~ All_HardCoral + wave_mean + Population_2017 ")  
#Had to remove rugosity because it was correlated with hard coral

fit1 <- lm(form1, data=alldat.site)
fit2 <- lm(form2, data=alldat.site)
fit3 <- lm(form3, data=alldat.site)


# Note if some path equations only contain two variables, VIF test below is invalid
vif.test1<-vif(fit1)
vif.test2<-vif(fit2)
vif.test3<-vif(fit3)


### SUBSET only model variables from alldat
vars1<-all.vars(form1)
vars2<-all.vars(form2)
vars3<-all.vars(form3)
model.vars<-unique(c(vars1, vars2, vars3))


# Aggregate to fishing ground level, scale and center data:
sem.vars.ground<-aggregate(sem.vars.site[model.vars], FUN=mean, by=list(sem.vars.site$location))
sem.ground.scaled<-as.data.frame(apply(sem.vars.ground[,-1], 2, scale))
sem.ground.scaled<-cbind(sem.vars.ground[,1], sem.ground.scaled)
names(sem.ground.scaled)[1]<- "location"

waka.psem<-psem(lm(form1, data=sem.ground.scaled), 
                lm(form2, data=sem.ground.scaled),
                lm(form3, data=sem.ground.scaled))


basisSet(waka.psem)
# NOTE: So far, there are no independence claims so basis set = 0

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.psem, .progressBar = F))
sink()


# coefficients should already be standardized since data was already scaled
#coefs(waka.psem, standardize="scale")


#### RANDOM EFFECTS:
### Use site-level (unaggregated) data and incorporate random effects by location
# i.e., start with sem.dat.site

sem.site.scaled<-as.data.frame(apply(sem.vars.site[model.vars], 2, scale))
sem.site.scaled<-cbind(sem.vars.site$location, sem.site.scaled)
names(sem.site.scaled)[1]<- "location"

## For examples on hierarchical model specification see:
## http://www.rensenieuwenhuis.nl/r-sessions-21-multilevel-model-specification-nlme/

# Try to get simple models to converge before adding more complexity:
# Only random intercepts by location
wakarandom.psem<-psem(lme(form1, random = ~ 1 | location, data=sem.site.scaled), 
                      lme(form2, random = ~ 1 | location, data=sem.site.scaled),
                      lme(form3, random = ~ 1 | location, data=sem.site.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_", fish.col, "_randomIntercepts.txt", sep="")
sink(txtname)
print(summary(wakarandom.psem, .progressBar = F))
sink()

# For non-convergence problems, https://stats.stackexchange.com/questions/40647/lme-error-iteration-limit-reached
# use lmeControl to help with convergence?
#ctrl<-lmeControl(opt = "optim", maxIter)
#wakarandom.psem<-psem(lme(form1, random = ~ 1 | location, data=sem.site.scaled, method="ML", control=ctrl), 
#                      lme(form2, random = ~ 1 | location, data=sem.site.scaled, method="ML", control=ctrl),
#                      lme(form3, random = ~ 1 | location, data=sem.site.scaled, method="ML", control=ctrl))

## Does removing repeated measures help? NO
## sem.site.scaled2<-sem.site.scaled[!sem.site.scaled$location %in% c("hoga island", "tomia island", "wanci island", "komponaone island"),]
#wakarandom.psem<-psem(lme(form1, random = ~ 1 | location, data=sem.site.scaled2), 
#                      lme(form2, random = ~ 1 | location, data=sem.site.scaled2),
#                      lme(form3, random = ~ 1 | location, data=sem.site.scaled2))


## Random intercepts for site-level + 
## predictor that is allowed to vary between groups (e.g., effect of coral) + 
## group-level predictor (e.g., effect of total landings)
wakarandom.psem<-psem(lme(form1, random = ~ 1 + All_HardCoral | location, data=sem.site.scaled, method="ML"), 
                      #lme(form2, random = ~ 1 + Population_2017 | location, data=sem.site.scaled),
                      lme(form3, random = ~ 1 + All_HardCoral | location, data=sem.site.scaled, method="ML"))
     
wakarandom.psem<-psem(lme(form1, random = ~ 1 + All_HardCoral | location, data=sem.site.scaled, method="ML", control=ctrl), 
                      #lme(form2, random = ~ 1 + Population_2017 | location, data=sem.site.scaled),
                      lme(form3, random = ~ 1 + All_HardCoral | location, data=sem.site.scaled, method="ML", control=ctrl))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_", fish.col, "_randomInterceptsAndSlopes.txt", sep="")
sink(txtname)
print(summary(wakarandom.psem, .progressBar = F))
sink()

