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

### NOW, test for MULTICOLLINEARITY (calculate VIF)
### Seems like test for multicollinearity should be constructed for EACH linear model found within the SEM
### Equations for MULTICOLLINEARITY can be recycled as PSEM equations
analysis.col<-grep(fish.col, names(alldat.site))

form1<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + landings_mean_onisland", sep=""))
form2<-as.formula("landings_mean_onisland ~ Population_2017 + landings_onisland_prop")
form3<-as.formula("All_HardCoral ~ Population_2017")

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


sem.vars.site<-alldat.site[c("location", "type_reef", model.vars)]

# Aggregate to fishing ground level, scale and center data:
sem.vars.ground<-aggregate(sem.vars.site[model.vars], FUN=mean, by=list(sem.vars.site$location))
sem.ground.scaled<-as.data.frame(apply(sem.vars.ground[,-1], 2, scale))
sem.ground.scaled<-cbind(sem.vars.ground[,1], sem.ground.scaled)
names(sem.ground.scaled)[1]<-"location"

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
sem.site.scaled<-cbind(sem.vars.site$location, sem.vars.site$type_reef, sem.site.scaled)
names(sem.site.scaled)[1:2]<- c("location", "reef_type")

## For examples on hierarchical model specification see:
## http://www.rensenieuwenhuis.nl/r-sessions-21-multilevel-model-specification-nlme/

## Random intercepts for site-level + 
## predictor that is allowed to vary by reef type (e.g., effect of coral and landings on fish) + 
## group-level predictor (e.g., effect of total landings by location)
wakarandom.psem<-psem(lme(form1, random = ~ 1 + All_HardCoral + landings_mean_onisland | reef_type, data=sem.site.scaled), 
                      lme(form2, random = ~ 1 + Population_2017 | location, data=sem.site.scaled),
                      lme(form3, random = ~ 1 + Population_2017 | reef_type, data=sem.site.scaled))


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


