# Get Wakatobi population data

# References: 
# set the ggmap API key: https://www.littlemissdata.com/blog/maps
# ggmap cheatsheet, also how to read-in a .SHP file: https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/ggmap/ggmapCheatsheet.pdf
# intro to sp and sf - two packages at the forefront of handling spatial data in R: https://www.jessesadler.com/post/gis-with-r-intro/
# Tutotials/vignettes for sf: https://cran.r-project.org/web/packages/sf/
# Especially section on "Geometrical operations" for creating buffers around polygons and finding the intersection of polygons (will need this for human popualtion data): https://cran.r-project.org/web/packages/sf/vignettes/sf1.html


# Possible sources of Wakatobi polygons (still not sure if these contain village boundaries):
# https://www.naturalearthdata.com/
# https://www.statsilk.com/maps/download-free-shapefile-maps
# Use: (1) this site to find polygon ID: https://nominatim.openstreetmap.org/ then (2) search by ID# on this site to download polygon as a .POLY, GeoJSON, or WKT file: http://polygons.openstreetmap.fr/index.py
# GADM: https://gadm.org/download_country_v3.html also see: https://gadm.org/formats.html


#if using ggmap package:
#Get the latest Install of ggmap
#if(!requireNamespace("devtools")) install.packages("devtools")
#devtools::install_github("dkahle/ggmap", ref = "tidyup", force=TRUE)
#library(ggmap) #Load the ggmap library
#Set your API Key
#ggmap::register_google(key = "AIzaSyBuvn8zwc5UPfMBQykn0dcOOwjCZ--WHLs") #Notes: If you get still have a failure try restarting R? 
# Get google map:
#ggmap(get_googlemap(center = c(lon = LONGITUDE coords, lat = LATITUDE coords),
#                    zoom = ?, scale = ?,
#                    maptype ='terrain',
#                    color = 'color'))


#if using sp package to read in map data:
#library(sp) # old way is readRDS in sp package
#readRDS("~/Analyses_notGit/fish-otakotak/indo-dat/Wakatobi/_MapData/GADM/R-sf/gadm36_IDN_0_sf.rds") # level 0 thru 4 - different scales of Indonesia?


#here we'll use the sf package
rm(list=ls())
outdir<-"~/Analyses/_RESULTS/SES-SEM/" 


library(sf)
library(ggplot2)
library(googledrive)


# Read in sf objects from GADM
# level 0 is lowest resolution (country-level); level 4 is highest resolution (village-level?)
# File Path: /Users/KGthatsme/Projects/Google Drive/Wakatobi-SEMAnalysis/_MapData/GADM/R-sf/gadm36_IDN_4_sf.rds
drive_auth() # Will require you to sign into Google account and grant permission to tidyverse for access 

# read-in village level sf file: get file ID from Google Drive's "shareable link" for the file: https://drive.google.com/open?id=1G4zUP5w2AmdCGRemqIdz_bgSfjjABf49
drive_download(as_id("1G4zUP5w2AmdCGRemqIdz_bgSfjjABf49"), overwrite=TRUE) # Saves file to working directory 
indo_4_sf<-readRDS("gadm36_IDN_4_sf.rds")
file.remove("gadm36_IDN_4_sf.rds") # Now that it's loaded into R, can delete file that was just downloaded

# read-in regency-level sf file (ie, wakatobi-level but no village info): https://drive.google.com/open?id=1SuZT6iTVVA7mKGaeOd6ubgCES6CInwDB
drive_download(as_id("1SuZT6iTVVA7mKGaeOd6ubgCES6CInwDB"), overwrite=TRUE) # Saves file to working directory 
indo_2_sf<-readRDS("gadm36_IDN_2_sf.rds")
file.remove("gadm36_IDN_2_sf.rds") # Now that it's loaded into R, can delete file that was just downloaded


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
 
# Plot indo_4_sf geometry: super dense village map of entire Indonesia
#setwd(outdir)
#pdf("st_geometry_4.pdf")
#plot(st_geometry(indo_4_sf)) # Plot just the geometry column
#dev.off()

# How to work with sf objects:
# indo_4_sf is an sf object; but behaves like a data.frame - example: table(indo_4_sf$NAME_2)
# str(indo_4_sf)
# print(indo_4_sf[6], n=5) # print column 6, but only the first 5 features
# column 4 "NAME_2": equilvanet to Wakatobi 
# column 6 "NAME_3": equivalent to "Binongko" or "Togo Binongko" (district-level) # binongko.names<-grep("Binongko", names(table(indo_4_sf$NAME_3))); table(indo_4_sf$NAME_3)[binongko.names]
# column 8 "NAME_4": equivalent to "Wanci" (village?) within Wangi-wangi (northern) district 

# Subset just Wakatobi (with villages)
wakatobi_sf<-indo_4_sf[indo_4_sf$NAME_2=="Wakatobi",]
wakatobi_sf

names(wakatobi_sf)[grep("^NAME_4$", names(wakatobi_sf))]<-"Village"
names(wakatobi_sf)[grep("^NAME_3$", names(wakatobi_sf))]<-"District"

# Do the same for Regency-level sf file (ie, Wakatobi but no villages) 
wakatobi_sf2<-indo_2_sf[indo_2_sf$NAME_2=="Wakatobi",]
wakatobi_sf2



# Plot just the geometry column
setwd(outdir)
pdf("map_wakatobi.pdf")
p<-ggplot(data=wakatobi_sf) +
  geom_sf() +
  allfiguretheme
print(p)
dev.off()


# Do the same for regency-level
setwd(outdir)
pdf("map_wakatobi2.pdf")
p<-ggplot(data=wakatobi_sf2) +
  geom_sf() +
  allfiguretheme
print(p)
dev.off()

########################################################################
########################################################################
#### Add dive sites to map (for brainstorming purposes - humans vs fish biomass)

### Read in site journal metadata
### File path: /Users/KGthatsme/Projects/Google Drive/Wakatobi-SEMAnalysis/site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv
### Google Drive Shareable Link: https://drive.google.com/open?id=1SNHtCmszbl6SYMPng1RLCDQVmap3e27n
drive_download(as_id("1SNHtCmszbl6SYMPng1RLCDQVmap3e27n"), overwrite=TRUE) 
fishsites<-read.csv("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
file.remove("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")


# Convert to an sf object:
fish_sf<-st_as_sf(fishsites, coords=c("long_dd", "lat_dd"), crs="+proj=longlat +datum=WGS84")

# We only need Site.Name and geometry columns
fish_sf<-fish_sf[,c("Site.Name", "geometry")]



# Full Wakatobi Map with Fish Site Names
setwd(outdir)
pdf("map_wakatobi_FishSiteNames.pdf")
p<-ggplot() +
  geom_sf(data=wakatobi_sf) +
  geom_sf(data=fish_sf, aes(shape=Site.Name, color=Site.Name), show.legend="point")+
  scale_shape_manual(values = 1:nlevels(fish_sf$Site.Name))+
  allfiguretheme
print(p)
dev.off()


# Full Wakatobi Map with Fish Site (no Names)
setwd(outdir)
pdf("map_wakatobi_FishSiteNoNames.pdf")
p<-ggplot() +
  geom_sf(data=wakatobi_sf) +
  geom_sf(data=fish_sf, show.legend=FALSE)+
  #geom_sf(data=fish_sf, aes(shape=Site.Name, color=location), show.legend="point")+
  #scale_color_viridis_d()+
  allfiguretheme
print(p)
dev.off()


########################################################################
########################################################################
# Plot dive sites and fishing grounds together

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


setwd(outdir)
pdf("map_wakatobi_fishingGrounds.pdf")
p<-ggplot() +
  geom_sf(data=grounds, aes(fill=Name), alpha=0.5, show.legend="polygon")+
  geom_sf(data=wakatobi_sf) +
  geom_sf(data=fish_sf, show.legend=FALSE)+
  coord_sf(xlim = c(123.366, 124.2202))+
  #geom_sf(data=fish_sf, aes(shape=Site.Name, color=location), show.legend="point")+
  #scale_color_viridis_d()+
  allfiguretheme
print(p)
dev.off()

# Repeat but plot just ONE fishing ground at a time
setwd(outdir)
for (i in 1:length(grounds$Name))
{
  pdfname<-paste("map_wakatobi_fishingGrounds-", grounds$Name[i], ".pdf", sep="")
  pdf(pdfname)
  p<-ggplot() +
    ggtitle(grounds$Name[i]) +
    geom_sf(data=wakatobi_sf2) +
    geom_sf(data=grounds[i,], colour="black", fill="black", alpha=0.5, show.legend=FALSE)+
    geom_sf(data=fish_sf, aes(shape=site_name, color=site_name), show.legend="point")+
    scale_shape_manual(values = 1:nlevels(fish_sf$site_name))+
    coord_sf(xlim = c(123.366, 124.2202))+
    #geom_sf(data=fish_sf, aes(shape=Site.Name, color=location), show.legend="point")+
    #scale_color_viridis_d()+
    allfiguretheme
  print(p)
  dev.off()
}



# READ IN CENSUS DATA (CSV FILE)
# File path: /Users/KGthatsme/Projects/Google Drive/Wakatobi-SEMAnalysis/_humanPopData/Wakatobi-HumanPopulationData-CLEANED.csv
# GoogleDrive shareable link: https://drive.google.com/open?id=1ye4n_XBXWO563RqHSh5WXgVOfq2OuV1l
drive_download(as_id("1ye4n_XBXWO563RqHSh5WXgVOfq2OuV1l"), overwrite=TRUE) 
popdata<-read.csv("Wakatobi-HumanPopulationData-CLEANED.csv")
file.remove("Wakatobi-HumanPopulationData-CLEANED.csv")


length(popdata$Village) # 100 villages in csv file
length(names(table(wakatobi_sf$Village))) # 101 unique villages in sf object: NAME_4 == "Waha" occurs twice (village in Tomia and village in Wangi-wangi both have the same name, "Waha")
length(wakatobi_sf$Village) # 102 villages in SF object (two are "extra") 

# Match village names in SF object with village names in CSV file
#name_4<-names(table(wakatobi_sf$NAME_4))
village_sf<-wakatobi_sf$Village
mismatch_csv<-(popdata$Village %in% village_sf)
popdata$Village[!mismatch_csv] # Names in popdata csv that don't match sf object

# Change the following names in CSV file to match SF object
#Laulua -> Lau-lua       
#Lewuto -> Lefuto       
#Kasuari -> Kaswari    
#Kel Sowa -> Sowa     
#Kel Popalia  -> Popalia 
#Desa Oihu -> Oihu    
#Desa Waloindi -> Waloindi 
#Desa Haka -> Haka


# VICE VERSA: change the following names in SF object to match CSV file
mismatch_sf<-(village_sf %in% popdata$Village)
village_sf[!mismatch_sf]
#Mandati Ii -> Mandati II
#Mandati Iii  -> Mandati III     
#Patua Ii -> Patua II 

# The following villages in SF object are "EXTRA" - ie, not found in census data; 
# Pulau        
# Pulau Sampora 
# Where are these?
##### NOTE: Checked google maps satellite view, both places appear uninhabited
wakatobi_sf[grep("^Pulau$", wakatobi_sf$Village),] # Village in Tomia
wakatobi_sf[grep("^Pulau Sampora$", wakatobi_sf$Village),] # Village in Wangi-wangi selatan

# Highlight these on a map
waka_unpopulated<-wakatobi_sf
waka_unpopulated$Status<-"Populated"
waka_unpopulated[grep("^Pulau$", waka_unpopulated$Village),]$Status<-"Unpopulated" # Give Pulau village red color
waka_unpopulated[grep("^Pulau Sampora$", waka_unpopulated$Village),]$Status<-"Unpopulated"


setwd(outdir)
pdf("map_wakatobi_check_noPopulationData.pdf")
#plot(st_geometry(waka_unpopulated), col = waka_unpopulated$population_color)
p<-ggplot(data=waka_unpopulated)+
  geom_sf(aes(fill=Status), color="black")+
  scale_fill_manual(values=c("white", "red"))+
  allfiguretheme
print(p)
dev.off()

####################################################
# Change village name in CSV file to match SF object
popdata$Village<-as.character(popdata$Village) # convert to character so that names can be changed

popdata[grep("^Desa Haka$", popdata$Village),]$Village<-"Haka"
popdata[grep("^Haka$", popdata$Village),]

popdata[grep("^Desa Waloindi$", popdata$Village),]$Village<-"Waloindi"
popdata[grep("^Waloindi$", popdata$Village),]

popdata[grep("^Desa Oihu$", popdata$Village),]$Village<-"Oihu"
popdata[grep("^Oihu$", popdata$Village),]

popdata[grep("^Kel Popalia$", popdata$Village),]$Village<-"Popalia"
popdata[grep("^Popalia$", popdata$Village),]

popdata[grep("^Kel Sowa$", popdata$Village),]$Village<-"Sowa"
popdata[grep("^Sowa$", popdata$Village),]

popdata[grep("^Kasuari$", popdata$Village),]$Village<-"Kaswari"
popdata[grep("^Kaswari$", popdata$Village),]

popdata[grep("^Lewuto$", popdata$Village),]$Village<-"Lefuto"
popdata[grep("^Lefuto$", popdata$Village),]

popdata[grep("^Laulua$", popdata$Village),]$Village<-"Lau-lua"
popdata[grep("^Lau-lua$", popdata$Village),]

popdata$Village<-factor(popdata$Village) # convert back to factor



####################################################
# Change village name in SF object to match CSV file
# for grep: use ^ to signifiy start and $ to signify end

wakatobi_sf[grep("^Mandati Ii$", wakatobi_sf$Village),]$Village<-"Mandati II"
wakatobi_sf[grep("^Mandati II$", wakatobi_sf$Village),] # Second line in each pair checks that replacement worked (search new name: Mandati II)

wakatobi_sf[grep("^Mandati Iii$", wakatobi_sf$Village),]$Village<-"Mandati III"
wakatobi_sf[grep("^Mandati III$", wakatobi_sf$Village),]

wakatobi_sf[grep("^Patua Ii$", wakatobi_sf$Village),]$Village<-"Patua II"
wakatobi_sf[grep("^Patua II$", wakatobi_sf$Village),]


## RETEST names: 
village_sf<-wakatobi_sf$Village
mismatch2_csv<-(popdata$Village %in% village_sf) 
popdata$Village[!mismatch2_csv] # Do any spreadsheet names mismatch with names in SF object

mismatch2_sf<-(village_sf %in% popdata$Village)
wakatobi_sf[!mismatch2_sf,] # Only Pulau in Tomia and Pulau Sampora in Wangi-wangi selatan are not found in csv file (i.e., they're missing population census data)

##### Make separate maps for each district


binongko_sf<-wakatobi_sf[wakatobi_sf$District=="Binongko",]
#st_centroid(binongko_sf) #st_point_on_surface is more appropriate for irregular polygons
binongko_point<-st_point_on_surface(binongko_sf)

kaledupa_sf<-wakatobi_sf[wakatobi_sf$District=="Kaledupa",]
kaledupa_point<-st_point_on_surface(kaledupa_sf)

kaledupa_s_sf<-wakatobi_sf[wakatobi_sf$District=="Kaledupa Selatan",]
kaledupa_s_point<-st_point_on_surface(kaledupa_s_sf)

binongko_t_sf<-wakatobi_sf[wakatobi_sf$District=="Togo Binongko",]
binongko_t_point<-st_point_on_surface(binongko_t_sf)

tomia_sf<-wakatobi_sf[wakatobi_sf$District=="Tomia",]
tomia_point<-st_point_on_surface(tomia_sf)

tomia_t_sf<-wakatobi_sf[wakatobi_sf$District=="Tomia Timur",]
tomia_t_point<-st_point_on_surface(tomia_t_sf)

wangi_sf<-wakatobi_sf[wakatobi_sf$District=="Wangi-Wangi",]
wangi_point<-st_point_on_surface(wangi_sf)

wangi_s_sf<-wakatobi_sf[wakatobi_sf$District=="Wangi-Wangi Selatan",]
wangi_s_point<-st_point_on_surface(wangi_s_sf)



setwd(outdir)
districtnames<-c("binongko_sf", "kaledupa_sf", "kaledupa_s_sf", 
                 "binongko_t_sf", "tomia_sf", "tomia_t_sf",
                 "wangi_sf", "wangi_s_sf")
centroidnames<-c("binongko_point", "kaledupa_point", "kaledupa_s_point",
                 "binongko_t_point", "tomia_point", "tomia_t_point",
                 "wangi_point", "wangi_s_point")
for(i in 1:length(districtnames))
{
  districtpdf<-paste("map_wakatobi_villagenames", districtnames[i], ".pdf", sep="")
  nextdistrict<-get(districtnames[i])
  nextcentroid<-get(centroidnames[i])
  pdf(districtpdf)
  p<-ggplot()+
    geom_sf(data=nextdistrict, (aes(color=Village)), show.legend=FALSE)+
    geom_sf(data=nextcentroid, aes(shape=Village, color=Village), show.legend="point")+
    scale_shape_manual(values=1:length(nextcentroid$Village))+
  allfiguretheme
  print(p)
  dev.off()
}



# CALCULATE DENSITIES for all population metrics: All calculation in km
################################################################################
################################################################################
################################################################################
#################### CHOSE POPULATION METRICS FOR CALCULATIONS and SET GRAPH LEGEND NAMES
popdata$Total_Motorboats<-popdata$Inboard_Motorboats + popdata$Outboard_Motorboats
popcolumns<-c("Population_2017", "No_of_Fishermen", "Row_Boats", "Total_Motorboats")
legendnames<-c(expression(paste("Population (", km^-2, ")", sep="")), 
               expression(paste("No. of Fishers (", km^-2, ")", sep="")),
               expression(paste("No. of Row Boats (", km^-2, ")", sep="")),
               expression(paste("No. of Motorboats (", km^-2, ")", sep="")))
################################################################################
################################################################################

popdensities<-popdata[popcolumns]/popdata$Area_km2
densitynames<-paste(popcolumns, "_per_km2", sep="")
names(popdensities)<-densitynames
popdata<-cbind(popdata, popdensities)



## Merge data and create population maps
wakatobi_pop<-merge(x=wakatobi_sf, y=popdata, by=c("Village", "District"), all.x=TRUE) # Keep all rows from wakatobi_sf even if these are not found in popdata (i.e., keep Pulau and Pulau Sampora in the mix)
wakatobi_pop$Village<-factor(wakatobi_pop$Village)

for (i in 1:length(densitynames))
{
  pdfname<-paste("map_wakatobi_humanData_", densitynames[i], ".pdf", sep="")
  #get(paste("wakatobi_pop$", densitynames[i], sep=""))
  pdf(pdfname)
  #wakatobi_pop$Population_2017_per_km2
  p<-ggplot(data=wakatobi_pop) +
    geom_sf(aes(fill=eval(as.symbol(densitynames[i])))) +
    scale_fill_viridis_c(option = "viridis", na.value="white", direction=-1) +
    labs(fill=legendnames[i])+
    allfiguretheme
  print(p)
  dev.off()
  
  
  # With fish sites:
  pdfname<-paste("map_wakatobi_humanDataWFishSites_", densitynames[i], ".pdf", sep="")
  pdf(pdfname)
  p<-ggplot(data=wakatobi_pop) +
    geom_sf(aes(fill=eval(as.symbol(densitynames[i])))) +
    scale_fill_viridis_c(option = "viridis", na.value="white", direction=-1) +
    geom_sf(data=fish_sf, show.legend=FALSE)+
    labs(fill=legendnames[i])+
    allfiguretheme
  print(p)
  dev.off()
  
}

# Maps of non-density columns
rawnames<-c("Population", "No. of Fishermen", "No. of Row Boats", "No. of Motorboats")
  
 
for (i in 1:length(popcolumns))
{
  pdfname<-paste("map_wakatobi_humanData_", popcolumns[i], ".pdf", sep="")
  #get(paste("wakatobi_pop$", densitynames[i], sep=""))
  pdf(pdfname)
  #wakatobi_pop$Population_2017_per_km2
  p<-ggplot(data=wakatobi_pop) +
    geom_sf(aes(fill=eval(as.symbol(popcolumns[i])))) +
    scale_fill_viridis_c(option = "viridis", na.value="white", direction=-1) +
    labs(fill=rawnames[i])+
    allfiguretheme
  print(p)
  dev.off()
  
  
  # With fish sites:
  pdfname<-paste("map_wakatobi_humanDataWFishSites_", popcolumns[i], ".pdf", sep="")
  pdf(pdfname)
  p<-ggplot(data=wakatobi_pop) +
    geom_sf(aes(fill=eval(as.symbol(popcolumns[i])))) +
    scale_fill_viridis_c(option = "viridis", na.value="white", direction=-1) +
    geom_sf(data=fish_sf, show.legend=FALSE)+
    labs(fill=rawnames[i])+
    allfiguretheme
  print(p)
  dev.off()
  
}

# Two types of coordinate reference systems: (1) georeference system (units of degrees) and (2) projected reference system (linear units, e.g., meters)
# See: https://www.jessesadler.com/post/gis-with-r-intro/
# Transform from geographic coordinate system (units = degrees) to projected coordinate system (units = meters)
# Need to specify projection model (find EPSG number)
fish_sf_merc<-st_transform(fish_sf, "+init=epsg:3857 +units=km") # global mercator projection
fish_sf_itm<-st_transform(fish_sf, "+init=epsg:23840 +units=km") # DGN95 / Indonesia TM-3 zone 51.2 from spatialreference.org
fish_sf_utm<-st_transform(fish_sf, "+init=epsg:23891 +units=km") # ID74 / UTM zone 51S from spatialreference.org
# Test: st_crs(fish_sf_utm) # all transformations result in units = meters

# Decide how big of a buffer to use (in km):
# distance = 2.5 # meant to represent human effects limited to the same side of the island


# Create a buffer around each one of the above projections; see if they map differently
distance=10 # in kilometers
fish_merc_buffer<-st_buffer(fish_sf_merc, distance)
fish_itm_buffer<-st_buffer(fish_sf_itm, distance)
fish_utm_buffer<-st_buffer(fish_sf_utm, distance)


# When plotting make sure the same projection is used for all sf objects 
# Transform wakatobi_sf to corresponding projection
wakatobi_sf_merc<-st_transform(wakatobi_pop, "+init=epsg:3857 +units=km") # global mercator projection
wakatobi_sf_itm<-st_transform(wakatobi_pop, "+init=epsg:23840 +units=km") # DGN95 / Indonesia TM-3 zone 51.2 from spatialreference.org
wakatobi_sf_utm<-st_transform(wakatobi_pop, "+init=epsg:23891 +units=km") # ID74 / UTM zone 51S from spatialreference.org



setwd(outdir)
buffproj<-c("fish_merc_buffer", "fish_itm_buffer", "fish_utm_buffer")
wakaproj<-c("wakatobi_sf_merc", "wakatobi_sf_itm", "wakatobi_sf_utm")
fishproj<-c("fish_sf_merc", "fish_sf_itm", "fish_sf_utm")
for(i in 1:length(buffproj))
{
  pdfname<-paste("map_wakatobi_checkProjections_", fishproj[i], "_", distance, "_km_buffer.pdf", sep="")
  pdf(pdfname)
  nextbuff<-get(buffproj[i])
  nextwaka<-get(wakaproj[i])
  nextfish<-get(fishproj[i])
  p<-ggplot() +
    geom_sf(data=nextwaka) +
    geom_sf(data=nextbuff, fill="transparent")+
    geom_sf(data=nextfish)+
    #geom_sf(data=fish_sf, aes(shape=Site.Name, color=location), show.legend="point")+
    #scale_shape_manual(values = 1:nlevels(fish_sf$Site.Name))+
    #scale_color_viridis_d()+
    allfiguretheme
  print(p)
  dev.off()
}


### Pick ONE projection for the rest of the analysis
# Pick wakatobi_sf_utm and fish_utm_buffer (EPSG=23891) since it's projection is centered on Wakatobi: http://www.spatialreference.org/ref/epsg/23891/

#### NEXT - Calculate intersection between buffers and village polygons, and calculate area-weighted population metric
intersect_dat<-st_intersection(wakatobi_sf_utm, fish_utm_buffer)

# Plot intersection:
pdf("map_wakatobi_intersection_of_villages_and_fishsites.pdf")
p<-ggplot() +
  geom_sf(data=fish_utm_buffer, fill="transparent")+
  geom_sf(data=wakatobi_sf_utm)+
  geom_sf(data=intersect_dat, fill="red") +
  geom_sf(data=fish_sf_utm)+
  #geom_sf(data=fish_sf, aes(shape=Site.Name, color=location), show.legend="point")+
  #scale_shape_manual(values = 1:nlevels(fish_sf$Site.Name))+
  #scale_color_viridis_d()+
  allfiguretheme
print(p)
dev.off()


# Each row in intersect_dat corresponds to a unique polygon intersection between wakatobi_sf_utm and fish_utm_buffer
table(intersect_dat$Site.Name) # Example: the site "Shark Point" intersects with 6 villages

# Calculate the area of each intersection polygon and add this as a data column
intersect_dat$intersection_area_km2<-as.numeric(st_area(intersect_dat))

# Calculate area weighted population metric
intersect_nogeo<-intersect_dat
st_geometry(intersect_nogeo)<-NULL # this converts back to data.frame
popdata_df<-intersect_nogeo$intersection_area_km2 * intersect_nogeo[,densitynames] # Population in each intersecting polygon

# popcolumns is defined above
names(popdata_df)<-popcolumns
popdata_df<-cbind(intersect_nogeo$Site.Name, popdata_df)
names(popdata_df)[1]<-"Site.Name"

intersect_pop<-aggregate(.~Site.Name, data=popdata_df, sum, na.action=na.omit) # Sum all polygons associated with the same fish site

setwd(outdir)
csvname<-paste("data_wakatobiHumans_areaWeightedDensityMetrics_", distance, "_km_buffer.csv", sep="")
write.csv(intersect_pop, file = csvname, quote = FALSE, row.names = FALSE)

# NEXT - Create distance weighted metrics of raw population data
# Use same CRS projections: here, wakatobi_sf_utm and fish_sf_utm
raw_centroid<-st_centroid(wakatobi_sf_utm) # raw_centroid is 102 rows (villages) by 19 columns (fish sites)
raw_dist<-st_distance(raw_centroid, fish_sf_utm)

colnames(raw_dist)<-fish_sf_utm$Site.Name
rownames(raw_dist)<-wakatobi_sf_utm$Village

# Try weighting by: 1 / (distance)^2
sq_dist<-raw_dist^2
units(sq_dist)<-NULL
dist_weights<-1/sq_dist
dist_weights<-as.data.frame(dist_weights)

# 
#dist_weights<-raw_dist / max(raw_dist)
#units(dist_weights)<-NULL


wakatobi_nogeo<-wakatobi_sf_utm
st_geometry(wakatobi_nogeo)<-NULL

global_human_df<-NA
for (i in 1:length(popcolumns))
{
  global_human_column<-NA
  for (j in 1:ncol(dist_weights))
  {
    global_human_effect<-wakatobi_nogeo[popcolumns[i]] * dist_weights[j]
    names(global_human_effect)<-colnames(dist_weights)[j]
    global_human_column<-cbind(global_human_column, global_human_effect)
  }
  global_human_sum<-colSums(global_human_column, na.rm=TRUE)
  global_human_df<-cbind(global_human_df, global_human_sum)
  colnames(global_human_df)[i+1]<-paste("distWeighted_", popcolumns[i], sep="")
}

global_human_df<-as.data.frame(global_human_df)
global_human_df<-global_human_df[-1,] 
global_human_df[,1]<-rownames(global_human_df) # Replace first column with Site.Names  
names(global_human_df)[1]<-"Site.Name"
global_human_df$Site.Name<-as.factor(global_human_df$Site.Name)

setwd(outdir)
csvname<-"data_wakatobiHumans_distanceWeighted.csv"
write.csv(global_human_df, file = csvname, quote = FALSE, row.names = FALSE)

# NEXT - generate dataset of distance from UVCs to nearest fish landing site and nearest export market
# As per Austin, most fish are landed in the "Mola" villages; FOR NOW: choose Mola Samaturu because it is central to all of these
mola_samaturu_centroid<-raw_centroid[grep("Mola Samaturu", raw_centroid$Village),]
dist_landings<-st_distance(mola_samaturu_centroid, fish_sf_utm)
units(dist_landings)<-NULL
dist_landings<-as.vector(dist_landings)
fish_nogeo<-fish_sf_utm
st_geometry(fish_nogeo)<-NULL
dist_landings_df<-as.data.frame(cbind(as.character(fish_nogeo$Site.Name), dist_landings))

names(dist_landings_df)[1]<-"Site.Name"

setwd(outdir)
csvname<-"data_wakatobiHumans_distanceToLandingsSite.csv"
write.csv(dist_landings_df, file = csvname, quote = FALSE, row.names = FALSE)


#### TO DO - create meap of distance weighted human metrics data: color code fish survey sites by human metrics
#### TO DO - add scale bar and north arrow (add this to object "allfiguretheme")
#### TO DO - create maps of distance weighted metrics

