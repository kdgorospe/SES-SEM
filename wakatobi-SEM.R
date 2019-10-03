# PATH ANALYSIS of Wakatobi SES
rm(list=ls())
library(piecewiseSEM)
library(nlme) # Version 3.1.122
library(ggplot2)
library(corrplot)
library(car) #VIF calculations
library(vegan)

# Read-in organized data
setwd("~/Analyses/_RESULTS/SES-SEM")
alldat.site<-read.csv("data_wakatobi_allDataMerged.csv")

# Set output directory
outdir<-"~/Analyses/_RESULTS/SES-SEM/"

################################################################################################################
################################################################################################################
responseDF<-as.data.frame(cbind(fish.response=c("log_biomass_g", "biomass_g", "no_of_species", "shannon", "invsimpson", "SimpsonEvenness"),
                                fish.title=c("log Total Biomass (g)", "Total Biomass (g)", "Richness", "Shannon Diversity (H')", "Inverse Simpson's Diversity (D2)", "Simpson's Evenness (E)" )))


#### FIRST, identify fish response column name, and title
## CHOICES:
## for biomass, set as biomass_g
## for log biomass, set as log_biomass_g
## for species richness, set as no_of_species
## for shannon diversity, set as shannon
## for inverse simpson's, set as invsimpson
## for simpson's evenness, set as SimpsonEvenness
fish.col<-"biomass_g" # Set response here
fish.row<-responseDF$fish.response %in% fish.col
fish.title<-as.character(responseDF[fish.row, "fish.title"])


# Setting "rish.response", "fish.title", and "fish.col" above allows for the remainder of code below to be flexible based on desired response variable


# NEXT: Create df for scatter plots (scatter.final) with all response + predictor variables + hierarchical variables
otherresponse<-responseDF$fish.response[!responseDF$fish.response %in% fish.col]
nonpredictors<-names(alldat.site) %in% c("site_id", "Site.Name", "lat_dd", "long_dd", as.character(otherresponse))
scatter.final<-alldat.site[,!nonpredictors]
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

# Create scatterplots
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


### Test for MULTICOLLINEARITY (calculate VIF)
### Test for multicollinearity should be constructed for EACH linear model found within the SEM
### Equations for MULTICOLLINEARITY can be recycled as PSEM equations
analysis.col<-grep(fish.col, names(alldat.site))

### FIRST, construct simple model that only uses site-level data (i.e., no fish landings data)
form1<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + Rugosity", sep=""))
form2<-as.formula("All_HardCoral ~ Rugosity")

fit1 <- lm(form1, data=alldat.site)
fit2 <- lm(form2, data=alldat.site)

# Note: if some path equations only contain one predictor, VIF test below is invalid
vif.test1<-vif(fit1)
vif.test2<-vif(fit2)

### SUBSET only model variables from alldat
vars1<-all.vars(form1)
vars2<-all.vars(form2)
model.vars_site<-unique(c(vars1, vars2))


sem.vars.site<-alldat.site[c("location", "type_reef", model.vars_site)]
sem.site.scaled<-as.data.frame(apply(sem.vars.site[model.vars_site], 2, scale))
sem.site.scaled<-cbind(sem.vars.site$location, sem.vars.site$type_reef, sem.site.scaled)
names(sem.site.scaled)[1:2]<- c("location", "reef_type")

# Site-level PSEM with no fish landings (i.e., fishing ground-level) data
waka.sitelevel.psem<-psem(lm(form1, data=sem.site.scaled), 
                lm(form2, data=sem.site.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.psem, .progressBar = F))
sink()


## Site-level PSEM with reef_type hierarchy
waka.sitelevel.reeftype.psem<-psem(lme(form1, random = ~ 1 | reef_type, data=sem.site.scaled), 
                      lme(form2, random = ~ 1 | reef_type, data=sem.site.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, "_reefTypeEffects.txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.reeftype.psem, .progressBar = F))
sink()


# Now, include fish landings (i.e., CATCH) data

form_a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + Rugosity + landings_sum_tot", sep=""))
form_b<-as.formula("All_HardCoral ~ Rugosity")

fit_a <- lm(form_a, data=alldat.site)
fit_b <- lm(form_b, data=alldat.site)

# Note: if some path equations only contain one predictor, VIF test below is invalid
vif.test_a<-vif(fit_a)
vif.test_b<-vif(fit_b)

### SUBSET only model variables from alldat
vars_a<-all.vars(form_a)
vars_b<-all.vars(form_b)
model.vars_catch<-unique(c(vars_a, vars_b))


sem.vars.catch<-alldat.site[c("location", "type_reef", model.vars_catch)]
sem.catch.scaled<-as.data.frame(apply(sem.vars.catch[model.vars_catch], 2, scale))
sem.catch.scaled<-cbind(sem.vars.catch$location, sem.vars.catch$type_reef, sem.catch.scaled)
names(sem.catch.scaled)[1:2]<- c("location", "reef_type")


waka.catch.psem<-psem(lm(form_a, data=sem.catch.scaled), 
                lm(form_b, data=sem.catch.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withCatchData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.catch.psem, .progressBar = F))
sink()


## Catch data PSEM with reef_type hierarchy
waka.catch.reeftype.psem<-psem(lme(form_a, random = ~ 1 | reef_type, data=sem.catch.scaled), 
                                   lme(form_b, random = ~ 1 | reef_type, data=sem.catch.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withCatchData_", fish.col, "_reefTypeEffects.txt", sep="")
sink(txtname)
print(summary(waka.catch.reeftype.psem, .progressBar = F))
sink()

### Notice: coefficient signs (positive vs negative) make more sense now


## Catch data PSEM with reef_type AND location hierarchies
#### LEFT OFF HERE: the following takes too long, simplify further...
waka.catch.allgroups.psem<-psem(lme(form_a, random = ~ 1 + All_HardCoral + Rugosity + landings_sum_tot | reef_type/reef_type/reef_type/location, data=sem.catch.scaled), 
                               lme(form_b, random = ~ 1 | reef_type, data=sem.catch.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withCatchData_", fish.col, "_reefTypeAndLocationEffects.txt", sep="")
sink(txtname)
print(summary(waka.catch.allgroups.psem, .progressBar = F))
sink()


# coefficients should already be standardized since data was already scaled
#coefs(waka.psem, standardize="scale")


#### RANDOM EFFECTS:
### Use site-level (unaggregated) data and incorporate random effects by location
# i.e., start with sem.dat.site


## For examples on hierarchical model specification see:
## http://www.rensenieuwenhuis.nl/r-sessions-21-multilevel-model-specification-nlme/

## Random intercepts for site-level + 
## predictor that is allowed to vary by reef type (e.g., effect of coral and landings on fish) + 
## group-level predictor (e.g., effect of total landings by location)
#wakarandom.psem<-psem(lme(form1, random = ~ 1 + All_HardCoral + landings_mean_onisland | reef_type, data=sem.site.scaled), 
#                      lme(form2, random = ~ 1 + Population_2017 | location, data=sem.site.scaled),
#                      lme(form3, random = ~ 1 + Population_2017 | reef_type, data=sem.site.scaled))


#setwd(outdir)
#txtname<-paste("stats_wakatobiSEM_", fish.col, "_randomIntercepts.txt", sep="")
#sink(txtname)
#print(summary(wakarandom.psem, .progressBar = F))
#sink()

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


