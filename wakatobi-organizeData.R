# Requires permission (shared Google drive folder) and access (internet) to Wakatobi Data
rm(list=ls())

#update.packages(ask = FALSE, checkBuilt = TRUE)
library(googledrive)
library(tidyverse)
library(vegan)
library(codyn) #Simpson's evenness calculation
library(sf)


# FIRST: setwd for where you want outputs saved: 
outdir<-"~/Analyses/_RESULTS/SES-SEM/"

# Data files are in GoogleDrive
drive_auth() # Will require you to sign into Google account and grant permission to tidyverse for access 

# input / munge fish data
# get file ID from Google Drive's "shareable link" for the file: https://drive.google.com/open?id=11jIP-ZlqgE9q2C046Kc_kg-7T3mjgelD
# masterdat<-read.csv("/Users/KGthatsme/Projects/Google Drive/Wakatobi-SEMAnalysis/_fishData/fish_df.csv")
drive_download(as_id("11jIP-ZlqgE9q2C046Kc_kg-7T3mjgelD"), overwrite=TRUE) # Saves file to working directory 
pauldat<-read.csv("fish_df.csv") 
file.remove("fish_df.csv") # Now that it's loaded into R, can delete file that was just downloaded

### NOTE fish_df is created by mungeing MASTER_fish_data.csv through wakatobi-cleanFishData.R (author: Paul Carvalho)
# Paul's species filters: only include non-pelagic, reef-associated fish (based on Aaron MacNeil's list), minus those species that the Indo team didn't survey well
#pauldat<-masterdat %>%
#  filter(region=="wakatobi") %>%
#  filter(site_name!="Sombano") #%>% # missing benthic data; result: 3892 rows
#  filter(family == "acanthuridae" | family == "aulostomidae" | family == "balistidae" |
#         family == "caesionidae" | family == "cirrhitidae" | family == "diodontidae" |
#         family == "fistulariidae" | family == "haemulidae" | family == "kyphosidae" |
#         family == "labridae" | family == "lutjanidae" | family == "monacanthidae" |
#         family == "mullidae" | family == "nemipteridae" | family == "ostraciidae" |
#         family == "pomacanthidae" | family == "pomacentridae" | family == "priacanthidae" |
#         family == "scaridae" | family == "serranidae" | family == "siganidae" |
#         family == "tetraodontidae" | family == "zanclidae" | family == "chaetodontidae") %>%
#  filter(species != "Acanthurus coeruleus" & species != "Kyphosus sectatrix" &
#           species != "Sectator ocyurus" & species != "Aphareus furca") %>%
#  droplevels()



# BEGIN CREATING RESPONSE VARIABLES FOR ANALYSIS
wakadat<-pauldat %>%
  filter(region=="wakatobi") %>%
  filter(site_name!="Sombano") %>% # missing benthic data; result: 3892 rows
  droplevels()

# 1 - BIOMASS: Calculate total fish biomass at site level (site_id column), averaged across three transects
# Old code using base R functions:
#fish.mass1<-aggregate(biomass_g ~ site_id + transect, data=fishdat, FUN=sum)
#fish.mass<-aggregate(biomass_g ~ site_id, data=fish.mass1, FUN=mean)
# Weight = a * Length^b
fishdat<-wakadat %>%
  mutate(biomass_g=a*size_cm^b) %>%
  arrange(site_name, transect)

plot(fishdat$biomass_g)

fish.mass<-fishdat %>%
  group_by(site_name, transect) %>%
  summarise(transect_mass=sum(biomass_g, na.rm=TRUE)) %>%
  summarise(biomass_g=mean(transect_mass))

plot(fish.mass$biomass_g)


# Calculate different metrics of species diversity and test all of them
# See Morris et al. (Ecology and Evolution) for discussion simultaneously considering analyses
# of multiple indices can provide greater insight

# 2 - richness 
fish.rich<-fishdat %>% 
  group_by(site_name, transect) %>%
  summarise(transect.count=n())  %>%
  summarise(no_of_species=mean(transect.count))
  

### NEXT, DIVERSITY INDICES: for now, maintain focus on "trophic group" diversity" because "functional diversity" is essentially "trait diversity" (see Bellwood et al. 2019)
### FIRST, FILL IN MISSING DATA
fishdat %>%
  filter(is.na(trophic_group)) %>%
  distinct(genus_species)

# Species with missing trophic data:
needs_troph<-c("Hemigymnus melapterus", "Leptojulis urostigma", "Hemigymnus fasciatus", "Choreodon fasciatus", "Halichoeres melanurus",
               "Thalassoma lutescens", "Thalassoma spp.", "Cirrhilabrus filamentosus", "Paracheilinus angulatus", "Pteragogus guttatus")

# See if it's possible to identify trophic group by inspecting the full dataset (i.e., pauldat)
pauldat %>%
  filter(genus_species %in% needs_troph) %>%
  filter(!is.na(trophic_group)) %>%
  select(genus_species, trophic_group) %>%
  distinct(genus_species, trophic_group)
# Use this result to populate NAs in fishdat:
fishdat[(is.na(fishdat$trophic_group) & fishdat$genus_species=="Hemigymnus melapterus"),]$trophic_group<-"Benthic Invertivore" #+2
fishdat[(is.na(fishdat$trophic_group) & fishdat$genus_species=="Hemigymnus fasciatus"),]$trophic_group<-"Benthic Invertivore" #+2
fishdat[(is.na(fishdat$trophic_group) & fishdat$genus_species=="Halichoeres melanurus"),]$trophic_group<-"Benthic Invertivore" #+1
fishdat[(is.na(fishdat$trophic_group) & fishdat$genus_species=="Thalassoma lutescens"),]$trophic_group<-"Benthic Invertivore" #+2

#Still missing trophic groups:
fishdat %>%
  filter(is.na(trophic_group)) 

## FILL IN THE REST BASED ON FISHBASE INFORMATION (DIET TABLE, etc)
fishdat[(is.na(fishdat$trophic_group) & fishdat$genus_species=="Leptojulis urostigma"),]$trophic_group<-"Benthic Invertivore" # feeds on benthic animals, polychaets, crustaceans
fishdat[(is.na(fishdat$trophic_group) & fishdat$genus_species=="Choreodon fasciatus"),]$trophic_group<-"Benthic Invertivore" # feed on mollusks, crustaceans, various worms, and echinoderms
fishdat[(is.na(fishdat$trophic_group) & fishdat$genus_species=="Cirrhilabrus filamentosus"),]$trophic_group<-"Planktivore" # feed above substrate on zooplankton
fishdat[(is.na(fishdat$trophic_group) & fishdat$genus_species=="Paracheilinus angulatus"),]$trophic_group<-"Planktivore" # see "Food Items" info

#These observations should be dropped (no trophic level information):
#Thalassoma spp: can't assume it's trophic level (could be benthic invertivore, planktivore, carnivore...)
#Pteragogus guttatus: no info in fishbase
fishdat<-fishdat %>%
  filter(!is.na(trophic_group))

# NOW THAT MISSING DATA IS FILLED IN, calculate 3 - shannon diversity (aka H') - using "TROPHIC_GROUP" 
trophcount_site<-fishdat %>%
  group_by(site_name, transect, trophic_group) %>%
  summarise(transect.sum=sum(abundance)) %>%
  group_by(site_name, trophic_group) %>%
  summarise(number_of_fish=mean(transect.sum)) 

trophcount_spread<-trophcount_site %>%
  spread(key=trophic_group, value=number_of_fish)

troph_matrix<-data.matrix(trophcount_spread[,-1])
fish.shan<-diversity(troph_matrix, index="shannon")
names(fish.shan)<-trophcount_spread$site_name
fish.shan<-as.data.frame(fish.shan)
fish.shan<-cbind(rownames(fish.shan), fish.shan)
rownames(fish.shan)<-NULL
names(fish.shan)<-c("site_name", "shannon")


# 3 - use INVERSE of simpsons diversity (aka D2), more commonly used than original D1 index (see Morris et al. 2014)
fish.isim<-diversity(troph_matrix, index="invsimpson")
names(fish.isim)<-trophcount_spread$site_name
fish.isim<-as.data.frame(fish.isim)
fish.isim<-cbind(rownames(fish.isim), fish.isim)
rownames(fish.isim)<-NULL
names(fish.isim)<-c("site_name", "invsimpson")


# 4 - calculate Simpson's Evenness (use codyn package) - requires dataframe of counts (not matrix)
fish.even<-community_structure(trophcount_site, abundance.var="number_of_fish", replicate.var="site_name", metric="SimpsonEvenness")
fish.even<-fish.even[,-2] # Remove richness column
names(fish.even)[2]<-"evenness"
### OTHER POTENTIAL RESPONSE VARIABLES: functional trait diversity using "fun.mass" above?


# NEXT: plot histograms of all response variables, consider log transforming some variables

# Total biomass
setwd(outdir)
pdf(file="plot_histogram.totalbiomass.pdf")
p<-ggplot(fish.mass, aes(x=biomass_g))+
  geom_histogram(bins=7)
print(p)
dev.off()

# Try log biomass
pdf(file="plot_histogram.LOGtotalbiomass.pdf")
p<-ggplot(fish.mass, aes(x=log10(biomass_g)))+
  geom_histogram(bins=7)
print(p)
#hist(log10(fish.mass[,"biomass_g"]),xlab="log Biomass", main="Histogram of Site-Level log Biomass")
dev.off()

# Option to use log biomass as response variable:
fish.logmass<-fish.mass %>%
  mutate(log_biomass_g=log10(biomass_g)) %>%
  select(site_name, log_biomass_g)
#fish.logmass$biomass_g<-log10(fish.logmass[,"biomass_g"])

# Richness
setwd(outdir)
pdf(file="plot_histogram.richness.pdf")
p<-ggplot(fish.rich, aes(x=no_of_species))+
  geom_histogram(bins=7)
print(p)
#hist(fish.rich[,"no_of_species"],xlab="Richness", main="Histogram of Site-Level Species Richness")
dev.off()

# Shannon
setwd(outdir)
pdf(file="plot_histogram.shannon.pdf")
p<-ggplot(fish.shan, aes(x=shannon))+
  geom_histogram(bins=7)
print(p)
#hist(fish.shan[,"shannon"],xlab="Shannon Diversity", main="Histogram of Site-Level Diversity")
dev.off()

# inverse Simpson
setwd(outdir)
pdf(file="plot_histogram.invsimpson.pdf")
p<-ggplot(fish.isim, aes(x=invsimpson))+
  geom_histogram(bins=7)
print(p)
#hist(fish.isim[,"invsimpson"],xlab="Inverse Simpson's Diversity", main="Histogram of Site-Level Diversity")
dev.off()

# Simpson's Evenness
setwd(outdir)
pdf(file="plot_histogram.evenness.pdf")
p<-ggplot(fish.even, aes(x=evenness))+
  geom_histogram(bins=7)
print(p)
#hist(fish.even[,"SimpsonEvenness"],xlab="Simpson's Evenness", main="Histogram of Site-Level Diversity")
dev.off()

# merge all fish data together
fish.dat<-full_join(fish.mass, fish.logmass, by="site_name")
fish.dat<-full_join(fish.dat, fish.rich, by="site_name")
fish.dat<-full_join(fish.dat, fish.shan, by="site_name")
fish.dat<-full_join(fish.dat, fish.isim, by="site_name")
fish.dat<-full_join(fish.dat, fish.even, by="site_name")

# input / munge benthic cover data: https://drive.google.com/open?id=1ba04__uY3alCvHNXI1CLmQmInstLwach
drive_download(as_id("1ba04__uY3alCvHNXI1CLmQmInstLwach"), overwrite=TRUE) # Saves file to working directory 
coraldat<-read.csv("raw_coral_data_wakatobi_may_2018.csv") 
file.remove("raw_coral_data_wakatobi_may_2018.csv")

table(coraldat$Dive.Site, coraldat$Transect) # 3 Transects per site (100 measurements per transect)
table(coraldat$Life.Form)

# NEED TO DO SOME CLEANING:
# cm -> CM
# s -> S

agg.clean<-coraldat %>%
  mutate(Life.Form = fct_recode(Life.Form, "S"="s", "CM"="cm")) %>%
  group_by(Dive.Site) %>%
  count(Life.Form)

# OLD BASE-R CODE:
#coraldat$Life.Form[coraldat$Life.Form=="s"]<-"S"
#coraldat$Life.Form[coraldat$Life.Form=="cm"]<-"CM"
#coraldat$Life.Form<-factor(coraldat$Life.Form) # Reset levels
agg.clean.tmp<-aggregate(Life.Form ~ Dive.Site, data=coraldat, FUN = table)

setwd(outdir)
write.csv(agg.clean, file="data_wakatobi_benthicCountsTable-allcategories.csv", row.names=FALSE, quote = FALSE)

# OLD Base-R code
#site.counts<-as.data.frame(table(coraldat$Dive.Site))
#benth.site<-(agg.clean.tmp[,-1] / site.counts$Freq)
#rownames(benth.site)<-as.character(agg.clean.tmp$Dive.Site)

#hardCor<-c("ACB", "ACD", "ACE", "ACT", "CB", "CE", "CF", "CHE", "CM", "CME", "CMR", "CS")
#allAbio<-c("R", "RCK", "S")

benth.site<-agg.clean %>%
  group_by(Dive.Site) %>%
  mutate(site.counts=sum(n)) %>% # total benthic observations per transect
  mutate(benth.cover=n/site.counts) %>% # frequency of each Life.Form type
  select(Dive.Site, Life.Form, benth.cover) %>%
  spread(key=Life.Form, value=benth.cover) %>%
  select(Dive.Site, DCA, MA, R, RCK, SC, S) %>% # FOR NOW, remove: (1) all taxa and morphology-specific coral data; (2) OTHER (very little variation); (3) Sponge - not interested
  replace_na(list(DCA=0, MA=0, R=0, RCK=0, SC=0, S=0))
  
benth.site2<-agg.clean %>%
  group_by(Dive.Site) %>%
  mutate(site.counts=sum(n)) %>%
  mutate(new.benthic=recode_factor(Life.Form, "ACB"="All_HardCoral", 
                                   "ACD"="All_HardCoral",
                                   "ACE"="All_HardCoral",
                                   "ACT"="All_HardCoral",
                                   "CB"="All_HardCoral",
                                   "CE"="All_HardCoral",
                                   "CF"="All_HardCoral",
                                   "CHE"="All_HardCoral",
                                   "CM"="All_HardCoral",
                                   "CME"="All_HardCoral",
                                   "CMR"="All_HardCoral",
                                   "CS"="All_HardCoral",
                                   "R"="All_Abiotic",
                                   "RCK"="All_Abiotic",
                                   "S"="All_Abiotic")) %>% # simplify benthic categories
  group_by(Dive.Site, new.benthic) %>% 
  mutate(new.n=sum(n)) %>% # Sum counts of new categories
  mutate(new.cover=new.n/site.counts) %>% # calculate new frequency of categories
  select(Dive.Site, new.benthic, new.cover) %>%
  distinct() %>%
  filter(new.benthic %in% c("All_HardCoral", "All_Abiotic")) %>% 
  spread(key=new.benthic, value=new.cover)
  
  

# Percent Cov Key:
# ACB -   Acropora branching
# ACD - Acropora digitate
# ACE -    Acropora encrusting
# ACT -   Acropora tabulate
# CB -   Branching coral
# CE -   Coral encrusting
# CF -    Foliose coral
# CHE - Heliopora
# CM - Massive coral
# CME - Millephora
# CMR - Mushroom coral
# CS - Submassive coral
# DCA -   Dead coral with algae
# MA -   Macro Algae
# OT -    Others
# R -    Ruble
# RCK - Rock
# S - Sand
# SC - Soft Coral
# SP - Sponge


benthcov.site<-full_join(benth.site, benth.site2, by="Dive.Site") %>%
  rename("site_name"="Dive.Site")
setwd(outdir)
write.csv(benthcov.site, "data_wakatobi_benthicPercentCover-allcategories.csv", quote = FALSE)


# input / munge rugosity data: https://drive.google.com/open?id=1bB-UTzGzF2CEJhrUsJH9xxZ067khCqgw
drive_download(as_id("1bB-UTzGzF2CEJhrUsJH9xxZ067khCqgw"), overwrite=TRUE)
rugdat<-read.csv("raw_rugosity_data_wakatobi_may_2018.csv")
file.remove("raw_rugosity_data_wakatobi_may_2018.csv")

rug.site<-rugdat %>% 
  group_by(Site.Name) %>%
  summarise(Rugosity = mean(Rugosity),
            n=n()) %>% # Good to include n=n() to see counts or sum(!is.na(x)) to see count of non-missing values for each grouping
  rename("site_name"="Site.Name") %>%
  select(site_name, Rugosity)
setwd(outdir)
write.csv(rug.site, "data_wakatobi_benthicRugosity.csv", quote=FALSE, row.names=FALSE)


# input 5km human population data: https://drive.google.com/open?id=1DcVqeVEx6yGksBqLGzbWhU8WxN6UcYBz

#drive_download(as_id("1DcVqeVEx6yGksBqLGzbWhU8WxN6UcYBz"), overwrite=TRUE) # Saves file to working directory 
#humanDensity.dat<-read.csv("data_wakatobiHumans_areaWeightedDensityMetrics_5_km_buffer.csv") # weights each village's population density by its area to get "total population" within 5km buffer
#file.remove("data_wakatobiHumans_areaWeightedDensityMetrics_5_km_buffer.csv")

# input 2.5km human population data: https://drive.google.com/open?id=1xf91oaXfqp-BDKIPcgK0bW-AfvJ3joHY
#drive_download(as_id("1xf91oaXfqp-BDKIPcgK0bW-AfvJ3joHY"), overwrite=TRUE) # Saves file to working directory 
#humanDensity.dat<-read.csv("data_wakatobiHumans_areaWeightedDensityMetrics_2.5_km_buffer.csv") # weights each village's population density by its area to get "total population" within 5km buffer
#file.remove("data_wakatobiHumans_areaWeightedDensityMetrics_2.5_km_buffer.csv")

# input 10 km human population data: https://drive.google.com/open?id=1c05yh4thZTqhEPMA0SrdRCuM8n86b1G1
drive_download(as_id("1c05yh4thZTqhEPMA0SrdRCuM8n86b1G1"), overwrite=TRUE) # Saves file to working directory 
humanDensity.dat<-read.csv("data_wakatobiHumans_areaWeightedDensityMetrics_10_km_buffer.csv") # weights each village's population density by its area to get "total population" within 5km buffer
file.remove("data_wakatobiHumans_areaWeightedDensityMetrics_10_km_buffer.csv")


# input oceanographic (and other) variables from MSEC: https://drive.google.com/open?id=12CErWykopoj2_gpQI47XYHbUEocOdOSr
drive_download(as_id("12CErWykopoj2_gpQI47XYHbUEocOdOSr"), overwrite=TRUE) # Saves file to working directory 
msec.dat<-read.csv("msec_out_5km.csv")
file.remove("msec_out_5km.csv")

###### FOR NOW, remove all human population data and other unnecessary columns
msec.dat<-subset(msec.dat, select=-c(no, long, lat, 
                                         npp_flag, 
                                         land_area_5km,
                                         wave_ww3_res, 
                                         pop1990_5km, pop2010_5km, pop2015_5km, pop2000_5km, dist_market))


# input SST data: https://drive.google.com/open?id=1ROPUFf6yi6r78vw9eTOa78WyqxFfYuBK
drive_download(as_id("1ROPUFf6yi6r78vw9eTOa78WyqxFfYuBK"), overwrite=TRUE) # Saves file to working directory 
sst.dat<-read.csv("Wakatobi_2018_SSTExtract.csv")
file.remove("Wakatobi_2018_SSTExtract.csv")

sst.dat<-subset(sst.dat, select=c(site_name, SST_stdev, SST_50Perc, SST_98perc, SST_2perc, SST_kurtosis, SST_skewness))


# input cleaned fish flow data: https://drive.google.com/open?id=1PRrdjBQ-aWsjKwO5gHV91VXKZ4cRxE_G
drive_download(as_id("1PRrdjBQ-aWsjKwO5gHV91VXKZ4cRxE_G"), overwrite=TRUE) # Saves file to working directory 
trip.dat<-read.csv("Wakatobi-landings_201909124_TRIP_cleanedFishFlows.csv")
file.remove("Wakatobi-landings_201909124_TRIP_cleanedFishFlows.csv")

# Trim down dataset
trip.dat<-subset(trip.dat, select=c(trip_id, fishing_grnd1,
                                    fishflow_abund,
                                    landings_sold_personally_abund,
                                    landings_sold_Papalele_abund,
                                    landings_sold_Pengumpul_abund,
                                    landings_eaten_abund,
                                    landings_given_abund))




## input aggregation file for landings trips: https://drive.google.com/open?id=1PkaXlA1r1RA6tUWX7Tm3sk3SPm7kJxMf
drive_download(as_id("1PkaXlA1r1RA6tUWX7Tm3sk3SPm7kJxMf"), overwrite=TRUE)
trip.agg<-read.csv("aggregationKey-FishingGround_PC.csv")
file.remove("aggregationKey-FishingGround_PC.csv")


# Aggregate (calculate mean) of groups of fishing grounds based on column: new_fg
# OLD BASE-R CODE
#trip.dat[!(trip.dat$fishing_grnd1 %in% trip.agg$original_fg),]
#trip.dat<-trip.dat[(trip.dat$fishing_grnd1 %in% trip.agg$original_fg),]
trip.dat %>%
  filter(!fishing_grnd1 %in% trip.agg$original_fg) # These are the trips (n=5) with missing fishing ground information (i.e., doesn't match list of fishing grounds in aggregation file)

trip.dat<-trip.dat %>% 
  filter(fishing_grnd1 %in% trip.agg$original_fg) %>% # Only keep trips with fishing ground information
  rename(original_fg="fishing_grnd1") %>%
  left_join(trip.agg, by = "original_fg")

table(trip.dat$new_fg)


# How do fish catch landings affect fish response
landings.dat<-trip.dat %>%
  group_by(new_fg) %>%
  summarise(landings_sum_tot = sum(fishflow_abund),
            landings_sum_personal = sum(landings_sold_personally_abund + landings_eaten_abund + landings_given_abund), # NOTE: PERSONAL USE = eaten, given away, or sold personally
            landings_sum_onisland = sum(landings_sold_personally_abund + landings_eaten_abund + landings_given_abund + landings_sold_Papalele_abund), # NOTE: ONISLAND = "personal use" + sold to papalele (on-island dealer0)
            landings_sum_papalele = sum(landings_sold_Papalele_abund), # NOTE: PAPALELE = on-island dealer
            landings_sum_pengumpul = sum(landings_sold_Pengumpul_abund), # NOTE: PENGUMPUL = off-island dealer
            landings_sum_market = sum(landings_sold_Pengumpul_abund + landings_sold_Papalele_abund), # NOTE: MARKET = sold to any market (PAPALELE + PENGUMPUL)
            
            landings_mean_tot = mean(fishflow_abund),
            landings_mean_personal = mean(landings_sold_personally_abund + landings_eaten_abund + landings_given_abund),
            landings_mean_onisland = mean(landings_sold_personally_abund + landings_eaten_abund + landings_given_abund + landings_sold_Papalele_abund),
            landings_mean_papalele = mean(landings_sold_Papalele_abund),
            landings_mean_pengumpul = mean(landings_sold_Pengumpul_abund),
            landings_mean_market = mean(landings_sold_Pengumpul_abund + landings_sold_Papalele_abund),
            
            landings_prop_personal = landings_sum_personal / landings_sum_tot,
            landings_prop_onisland = landings_sum_onisland / landings_sum_tot,
            landings_prop_papalele = landings_sum_papalele / landings_sum_tot,
            landings_prop_pengumpul = landings_sum_pengumpul / landings_sum_tot, 
            landings_prop_market = landings_sum_market / landings_sum_tot
            ) %>%
            rename("location"="new_fg")


### PLOT Landings Fish Flow data
allfiguretheme<-theme_bw()+
  theme(panel.grid.minor = element_line(colour="transparent"), panel.grid.major = element_line(colour="transparent"),
        panel.border=element_rect(colour="black", size=1.5),
        panel.background = element_rect(fill = "white"),
        axis.text.x=element_text(size=12, colour="black"), # Controls size of x axis tickmark labels
        axis.text.y=element_text(size=12, colour="black"), # Controls size of y axis tickmark labels
        legend.text=element_text(size=8),
        legend.title=element_text(size=12),
        legend.position="bottom",
        plot.title=element_text(size=20))

# read-in village level sf file: get file ID from Google Drive's "shareable link" for the file: https://drive.google.com/open?id=1G4zUP5w2AmdCGRemqIdz_bgSfjjABf49
drive_download(as_id("1G4zUP5w2AmdCGRemqIdz_bgSfjjABf49"), overwrite=TRUE) # Saves file to working directory 
indo_4_sf<-readRDS("gadm36_IDN_4_sf.rds")
file.remove("gadm36_IDN_4_sf.rds") # Now that it's loaded into R, can delete file that was just downloaded


# Subset just Wakatobi (with villages)
wakatobi_sf<-indo_4_sf[indo_4_sf$NAME_2=="Wakatobi",]
wakatobi_sf

names(wakatobi_sf)[grep("^NAME_4$", names(wakatobi_sf))]<-"Village"
names(wakatobi_sf)[grep("^NAME_3$", names(wakatobi_sf))]<-"District"



### Read in site journal metadata
### File path: /Users/KGthatsme/Projects/Google Drive/Wakatobi-SEMAnalysis/site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv
### Google Drive Shareable Link: https://drive.google.com/open?id=1SNHtCmszbl6SYMPng1RLCDQVmap3e27n
drive_download(as_id("1SNHtCmszbl6SYMPng1RLCDQVmap3e27n"), overwrite=TRUE) 
fishsites<-read.csv("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
file.remove("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")


# Convert to an sf object:
fish_sf<-st_as_sf(fishsites, coords=c("long_dd", "lat_dd"), crs="+proj=longlat +datum=WGS84")

# We only need Site.Name and geometry columns
fish_sf<-fish_sf[,c("site_name", "geometry")]


# Read-in fishing grounds shape file (from Melati)
# File path: /Users/KGthatsme/Projects/Google Drive/Wakatobi-SEMAnalysis/_MapData/FishingGrounds/Waka_files_4KG
# Because this is a shapefile, need to download entire shapefile folder
grounds.files<-drive_ls("Waka_files_4KG")
for(i in 1:length(grounds.files$id)){
  drive_download(as_id(grounds.files$id[i]), overwrite=TRUE) 
}

# Read-in shapefile
grounds<-st_read("F_grnd.shp")

# Now delete shapefile and associated files
for(i in 1:length(grounds.files$name)){
  file.remove(grounds.files$name[i]) 
}

# Read-in fishing grounds aggregation file (this one uses capitalized names)
# SHAREABLE LINK: https://drive.google.com/open?id=1llGZzqRbkLssnH3OKEzLUaO4usGsyf-G
drive_download(as_id("1llGZzqRbkLssnH3OKEzLUaO4usGsyf-G"), overwrite=TRUE) 
trip.agg_forMapping<-read.csv("aggregationKey-FishingGround_PC-ForMatchingWithShapefileNames.csv")
file.remove("aggregationKey-FishingGround_PC-ForMatchingWithShapefileNames.csv")

# Which polygons in shapefile do not have a corresponding fishing ground in aggregation file (spelling mistakes?)
shapefile_nomatch<-grounds %>%
  filter(!Name %in% trip.agg_forMapping$original_fg)

# and vice versa, which fishing grounds in aggregation file do not have a corresponding polygon in shapefile?
aggfile_nomatch<-trip.agg_forMapping %>%
  filter(!original_fg %in% grounds$Name)

fishflows<-grounds %>%
  rename("original_fg"="Name") %>%
  left_join(trip.agg_forMapping, by="original_fg") %>%
  na.omit() %>%
  rename("location"="new_fg") %>%
  left_join(landings.dat, by="location")

flowvar_index<-grep("landings", names(fishflows))
for (i in flowvar_index)
{
flowvar<-names(fishflows)[i]
wa <- ggplot() +
  geom_sf(data=fishflows, aes(fill=get(flowvar), color=get(flowvar))) +
  scale_fill_viridis_c(option = "plasma") +
  scale_color_viridis_c(option = "plasma")+
  geom_sf(data = wakatobi_sf, fill="grey45", color="black") +
  geom_sf(data=fish_sf, color="black", size=2) +
  allfiguretheme +
  coord_sf(xlim=c(123.3,124.2), ylim=c(-6.1,-5.2), expand=FALSE) +
  labs(x=NULL, y=NULL) 
newfile<-paste("map-Wakatobi-fishflow-", flowvar, ".pdf", sep="")
pdf(newfile)
print(wa)
dev.off()
}





###### Merge fish, oceanographic (MSEC), human pop data, rugosity, benthic cover, SST AND catch data using "site journal.xlsx" as site key: https://drive.google.com/open?id=1SNHtCmszbl6SYMPng1RLCDQVmap3e27n
drive_download(as_id("1SNHtCmszbl6SYMPng1RLCDQVmap3e27n"), overwrite=TRUE) # Saves file to working directory 
site.key<-read.csv("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
file.remove("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
site.key<-subset(site.key, select=c("site_name", "lat_dd", "long_dd", "exposed", "u_visibility", "type_reef", "location"))


# Merge all data: 
##### Do this in the following order: fish, fishing grounds (catch), oceanographic (MSEC), human pop data, rugosity, benthic cover, SST

# Merge fish data
alldat.site<-site.key %>%
  left_join(fish.dat, by="site_name") %>%
  left_join(landings.dat, by="location") %>% 
  replace(is.na(.), 0) %>% # previous left join now includes fish UVC site Furake (on Hoga Island) which is restricted from fishing; replace NAs in landings data with 0
  left_join(rug.site, by="site_name") %>%
  left_join(benthcov.site, by="site_name") %>%
  left_join(msec.dat, by="site_name") %>%
  left_join(sst.dat, by="site_name") %>%
  left_join(humanDensity.dat, by="site_name") %>%
  mutate_at(vars(Population_2017, No_of_Fishermen, Row_Boats, Total_Motorboats), ~replace_na(., 0)) %>%
  arrange(Population_2017)
  
### DIVIDE human metrics data by reef area
#alldat.site$Population_2017<-alldat.site$Population_2017/alldat.site$reef_area_5km
#alldat.site$No_of_Fishermen<-alldat.site$No_of_Fishermen/alldat.site$reef_area_5km
#alldat.site$Row_Boats<-alldat.site$Row_Boats/alldat.site$reef_area_5km
#alldat.site$Total_Motorboats<-alldat.site$Total_Motorboats/alldat.site$reef_area_5km
#tmp.col<-grep("reef_area_5km", names(alldat.site))
#alldat.site<-alldat.site[,-tmp.col]

setwd(outdir)
write.csv(alldat.site, "data_wakatobi_allDataMerged.csv", quote=FALSE, row.names = FALSE)
