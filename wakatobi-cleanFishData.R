##########################################################################################################################################
# Cleaning data for size spectra analysis
#
# Author: Paul Carvalho
# Date created: September 19, 2019
# Original Filename was: data_prep.R
##########################################################################################################################################

# SETUP ----------------------------------------------------------------------------------------------------------------------------------

# Clean workspace
rm(list=ls())

# Libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(rlang)
library(tidyr)

# Directories
setwd("C:/Users/pgcar/Google Drive/Paul Carvalho/dissertation/chapter 3/analysis/data")

# Load data
fish_df <- read_xlsx("Master_fish_data.xlsx")
benthic_df <- read_xlsx("Master_benthic_data.xlsx")
family_key <- read.csv("family_key.csv")


fg_key <- read.csv("species_fg_key.csv")

# FUNCTIONS ------------------------------------------------------------------------------------------------------------------------------
check.spelling <- function(df_in_1, var_name_1, df_in_2=NULL, var_name_2=NULL, distance_sensitivity){
  # turn bare args into quosures
  quo_var <- enquo(var_name_1)
  if(is.null(df_in_2) == FALSE){
    quo_var_1 <- enquo(var_name_2)  
  }
  
  # quo_text(quo_var) gets the string of the column name from bare parameter input
  char.list <- as.character(df_in_1[[quo_text(quo_var)]])
  if(is.null(df_in_2) == FALSE){
    char.list.1 <- as.character(df_in_2[[quo_text(quo_var_1)]])
  }

  # Remove na's from the list of characters
  char.list <- na.omit(char.list)
  if(is.null(df_in_2) == FALSE){
    char.list.1 <- na.omit(char.list.1)
  }

  # check spelling and obtain a list of words that have similar spelling
  x <- unique(char.list)
  df_list <- data_frame(x=x, name=rep("df_1",length(x)))
  # Only run this code if there is a second dataframe
  if(is.null(df_in_2) == FALSE){
    y <- unique(char.list.1)
    df_x <- data_frame(x=x, name=rep("df_1",length(x)))
    df_y <- data_frame(x=y, name=rep("df_2",length(y)))
    df_list <- rbind(df_x,df_y)
  }
  
  # Need to run this nested for loop to get the length of the data frame
  tmp_count <- 0
  for(i in 1:length(df_list$x)){
    for(j in i:length(df_list$x)){
      dist = adist(df_list$x[i],df_list$x[j])
      if(dist > 0 && dist <= distance_sensitivity){
        tmp_count <- tmp_count+1
      }
    }
  }
  df_dist = data.frame(name_1=character(length=tmp_count),
                       df_name_1=character(length=tmp_count),
                       name_2=character(length=tmp_count),
                       df_name_2=character(length=tmp_count),
                       distance=numeric(length=tmp_count))
  # need to make sure these vectors are characters and not factors
  df_dist$name_1 <- as.character(df_dist$name_1)
  df_dist$df_name_1 <- as.character(df_dist$df_name_1)
  df_dist$name_2 <- as.character(df_dist$name_2)
  df_dist$df_name_2 <- as.character(df_dist$df_name_2)

  count = 0
  for(i in 1:length(df_list$x)){
    for(j in i:length(df_list$x)){
      dist = adist(df_list$x[i],df_list$x[j])
      if(dist > 0 && dist <= distance_sensitivity){
        count <- count+1
        # print(df_dist$name_1[count])
        # print(df_list$x[i])
        df_dist$name_1[count]    <- df_list$x[i]
        df_dist$name_2[count]    <- df_list$x[j]
        df_dist$df_name_1[count] <- df_list$name[i]
        df_dist$df_name_2[count] <- df_list$name[j]
        df_dist$distance[count]  <- dist
        
      }
    }
  }

  # create a new list that omits repeated pairs
  tmp_name_1 <- NULL
  tmp_df_name_1 <- NULL
  tmp_name_2 <- NULL
  tmp_df_name_2 <- NULL
  tmp_distance <- NULL
  # refined_df_dist <- data.frame(name_1=character(),
  #                               df_name_1=character(),
  #                               name_2=character(),
  #                               df_name_2=character(),
  #                               distance=numeric())

  
  count.1 = 0
  for(k in 1:length(df_dist$name_1)){
    if(k == 1){
      count.1 <- count.1+1
      tmp_name_1[count.1]    <- df_dist$name_1[k]
      tmp_name_2[count.1]    <- df_dist$name_2[k]
      tmp_df_name_1[count.1] <- df_dist$df_name_1[k]
      tmp_df_name_2[count.1] <- df_dist$df_name_2[k]
      tmp_distance[count.1]  <- df_dist$distance[k]

    } else {
      flag = 0
      for(m in 1:length(tmp_name_1)){
        if(df_dist$name_1[k] == tmp_name_1[m] && df_dist$name_2[k] == tmp_name_2[m]){
          flag = 1
        } else if (df_dist$name_1[k] == tmp_name_2[m] && df_dist$name_2[k] == tmp_name_1[m]){
          flag = 1
        }
      }
      if(flag == 0){
        count.1 = count.1+1
        tmp_name_1[count.1]    <- df_dist$name_1[k]
        tmp_name_2[count.1]    <- df_dist$name_2[k]
        tmp_df_name_1[count.1] <- df_dist$df_name_1[k]
        tmp_df_name_2[count.1] <- df_dist$df_name_2[k]
        tmp_distance[count.1]  <- df_dist$distance[k]

      }
    }
  }
  refined_df_dist <- data.frame(name_1=tmp_name_1,
                                df_name_1=tmp_df_name_1,
                                name_2=tmp_name_2,
                                df_name_2=tmp_df_name_2,
                                distance=tmp_distance)
  list.order <- order(refined_df_dist$distance)
  refined_df_dist <- refined_df_dist[list.order,]

  # If function does not detect differences in spelling, then do nothing
  if(is.na(refined_df_dist$name_1[1])){
    cat("No differences detected for specified distance sensitivity!")

  # If differences detected, then enter interactive loop
  } else {
    # User interaction to replace spelling
    for(i in 1:(length(refined_df_dist$name_1))){
      cat("\n1.\"", as.character(refined_df_dist$name_1[i]),"\"","[",as.character(refined_df_dist$df_name_1[i]),"]",
          "   2.\"", as.character(refined_df_dist$name_2[i]),"[",as.character(refined_df_dist$df_name_2[i]),"]","\"\n", sep="")
      cat("Distance =", as.character(refined_df_dist$distance[i]))
      cat("\n1. If spelling 1 is correct")
      cat("\n2. If spelling 2 is correct")
      cat("\n3. If both incorrect and you want to enter the correct spelling")
      cat("\n4. If both correct or you do not want to make any substitutions")
      val <- eval(parse(text=readline(prompt="Enter option: ")))
      # Save correct and incorrect spelling
      if(val==1){
        correct.sp    <- refined_df_dist$name_1[i]
        incorrect.sp1 <- refined_df_dist$name_2[i]
        # substitute spelling and adjust factor levels - no changes will be made to df_in_2
        df_in_1[[quo_text(quo_var)]] <- gsub(incorrect.sp1, correct.sp, df_in_1[[quo_text(quo_var)]])
      } else if(val==2) {
        correct.sp    <- refined_df_dist$name_2[i]
        incorrect.sp1 <- refined_df_dist$name_1[i]
        # substitute spelling and adjust factor levels
        df_in_1[[quo_text(quo_var)]] <- gsub(incorrect.sp1, correct.sp, df_in_1[[quo_text(quo_var)]])
      } else if(val==3){
        correct.sp <- readline(prompt="Enter correct spelling: ")
        incorrect.sp1 <- refined_df_dist$name_1[i]
        incorrect.sp2 <- refined_df_dist$name_2[i]
        # substitute spelling and adjust factor levels
        df_in_1[[quo_text(quo_var)]] <- gsub(incorrect.sp1, correct.sp, df_in_1[[quo_text(quo_var)]])
        df_in_1[[quo_text(quo_var)]] <- gsub(incorrect.sp2, correct.sp, df_in_1[[quo_text(quo_var)]])
      } else if(val==4){
        # do nothing
      } else {
        cat("\nInvalid option")
      }
    } # End of interactive for loop
  }

  return(df_in_1)
}

fix_family <- function(df,family_key){
        for(i in 1:length(family_key$species)){
                spp_i <- as.character(family_key$species[i])
                fam_i <- as.character(family_key$family[i])
                df$family[which(df$species == spp_i)] = fam_i
        }
        return(df)
}

bin_sizes_5cm <- function(df){
     df$size_5cm_bin <- cut(df$size_cm, breaks=c(0,seq(5,60,5),Inf),
                            labels=c("0-5","6-10","11-15","16-20","21-25","26-30","31-35","36-40","41-45",
                                     "46-50","51-55","56-60",">60"))
     return(df)
}
 
# std_error <- function(val){
#   error <- (sd(val)) / (sqrt(length(val)))
#   return(error)
# }

# CLEAN MASTER FISH DATAFRAME ------------------------------------------------------------------------------------------------------------
# check spelling of species names
# fish_df1 <- check.spelling(df_in_1 = fish_df, var_name_1 = species, distance_sensitivity = 3)
# fish_df2 <- check.spelling(df_in_1 = fish_df1, var_name_1 = species, distance_sensitivity = 3)
# fish_df3 <- check.spelling(df_in_1 = fish_df2, var_name_1 = species, distance_sensitivity = 3)
# write.csv(fish_df3, "fish_df.csv")

# check the spelling of family names
# fish_df3$family <- tolower(fish_df3$family)
# fish_df4 <- check.spelling(df_in_1 = fish_df3, var_name_1 = family, distance_sensitivity = 3)
# write.csv(fish_df4, "fish_df.csv")

# fix missing family names
# fish_df4$species[(which(fish_df4$species == "Cheili"))] <- "Cheilinus bimaculatus" # fix species name
# fish_df4$family[(which(fish_df4$species == "Cheilinus bimaculatus"))] <- "labridae" # fix family name
# fish_df4$species[(which(fish_df4$species == "Stethojulis spp.."))] = "Stethojulis spp." # fix species name
# fish_df4$family[(which(fish_df4$family == "scarini-labridae"))] = "scaridae" # fix family name
# tmp_df1 <- fish_df4[which(is.na(fish_df4$family)),]
# write.csv(tmp_df1, "family_key.csv")
# fish_df5 <- fish_df4 %>%
#         fix_family(., family_key)
# fish_df6 <- check.spelling(df_in_1 = fish_df5, var_name_1 = family, distance_sensitivity = 3)
# write.csv(fish_df6, "fish_df.csv")
# fish_df6 <- read.csv("fish_df6.csv")

# Get only coral reef fish families for this study
fish_df7 <- fish_df6 %>%
        filter(family == "acanthuridae" | family == "aulostomidae" | family == "balistidae" |
               family == "caesionidae" | family == "cirrhitidae" | family == "diodontidae" |
               family == "fistulariidae" | family == "haemulidae" | family == "kyphosidae" |
               family == "labridae" | family == "lutjanidae" | family == "monacanthidae" |
               family == "mullidae" | family == "nemipteridae" | family == "ostraciidae" |
               family == "pomacanthidae" | family == "pomacentridae" | family == "priacanthidae" |
               family == "scaridae" | family == "serranidae" | family == "siganidae" |
               family == "tetraodontidae" | family == "zanclidae" | family == "chaetodontidae") %>%
        filter(species != "Acanthurus coeruleus" & species != "Kyphosus sectatrix" &
               species != "Sectator ocyurus" & species != "Aphareus furca") %>%
        mutate(genus_species = species) %>%
        separate(species,c("genus","species")) %>%
        filter(genus != "Pseudanthias" & genus != "Pterocaesio")

# Edit the size bins
fish_df8 <- bin_sizes_5cm(fish_df7)
# write.csv(fish_df8,"fish_df.csv")

# Add a and b coefficients for all fishes ====
fish_df <- read.csv("fish_df.csv")

# find row indices for missing information
idx <- which(is.na(fish_df$a))
missing <- fish_df[idx,]

# find unique species that are missing information
missing_uni <- as.character(unique(missing$genus_species))
missing_df <- data.frame(species = missing_uni,
                         family = NA,
                         functional_group = NA,
                         a = NA,
                         b = NA)

# write.csv(missing_df, "missing_coef.csv")
missing_df1 <- read.csv("missing_coef.csv")
missing_df1 <- missing_df1 %>%
    select(species, family, functional_group, a, b)

fish_df$functional_group <- as.character(fish_df$functional_group)

##
for(i in 1:length(idx)){
    j <- idx[i] # index for fish_df
    k <- which(as.character(missing_df1$species) == as.character(fish_df$genus_species[j])) # index for missing_df1
    fish_df$family[j] <- as.character(missing_df1$family[k])
    fish_df$functional_group[j] <- as.character(missing_df1$functional_group[k])
    fish_df$a[j] <- missing_df1$a[k]
    fish_df$b[j] <- missing_df1$b[k]
}


write.csv(fish_df, "fish_df.csv")
