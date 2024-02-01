
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code clean and combine NSFG and CPS data into one tidy data frame.
##
## 
##  Author: Benjamin Schlüter
##  Date: December 2023
##==============================================================================
##
##  Notes
## ------
## 1. 
## 
##==============================================================================



## LOAD PACKAGES ===============================================================

# Install/load packages
packages <- c("tidyverse", "rvest", "here", "survey",
              "SAScii", "readr")
for(p in packages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
}



## ===== FUNCTIONS =============================================================


## ===== LOAD & TIDY CPS DATA ==================================================

# Load CPS data (extract from IPUMS)
# read_csv allows to open compressed .gz 
data.cps <- read_csv(
    here(
        "data_private", 
        "cps_00009.csv.gz"
        )
    ) |> 
    rename_all(
        .funs = tolower
        ) |> 
    # Keep useful variables
    dplyr::select(
        year, month, serial,pernum, # identification
        metarea, metfips, county, # strata for sampling design
        momloc, momloc2, poploc, poploc2, frever,# parenthood
        wtfinl, # weights,
        statefip, # states
        age, sex, race, hispan # demographic
        ) |> 
    filter(
        # Consider year from 2008 (change in data collection)
        year >= 2008,
        # Question on children only asked to women
        sex == 2
    ) |> 
    mutate(
        # Create abridged race (!! NEED TO BE CHECKED !!)
        race_eth = case_when(race == 100 & hispan == 0 ~ "NH-White",
                             race == 200 & hispan == 0 ~ "NH-Black",
                             !(race %in% c(100, 200)) & hispan == 0 ~ "NH-Other",
                             hispan != 0 ~ "Hispanic"),
        # Create HH for survey design
        hhid = paste0(serial, year) |> as.numeric()
        )

# Data before 2012 need to be corrected (see Census Bureau report)
data.cps.2012_22 <- 
    data.cps |>
    filter(
        # Survey conducted every two years 
        year %in% seq(2012, 2022, 2),
        # aged 15-50 not having 999.
        # Only from 2012 CPS that itwees
        # are until 50 yo
        age >= 15,
        age <= 50,
        # Remove when info on frever is not collected
        frever != 999
    ) |> 
    mutate(
        # Prop mother
        mother = ifelse(frever == 0, 0, 1)
    )

# Correct the data for 2008 and 2010
# Use hh info and coresident child info to correct proportion childless
data.cps.2008_10 <- 
    data.cps |> 
    filter(
        # Two years were info of coresidency wasn't use
        year %in% c(2008, 2010),
        # aged 15-45 not having 999.
        # Only in from 2012 that itwees
        # are until 50 yo
        age >= 15,
        age < 45
    ) |> 
    mutate(
        # Set NA because if no child recorded in hh,
        # it can be that the child already left the hh
        # so we can't set default value to 0
        frever_cor = NA
    )

# Correct frever according to relation to child(ren)
correct.frever <- function(d, ...) {
    
    # Who is a mother
    find <- d$momloc != 0 | d$momloc2 != 0 
    # How many time are they listed parent
    n <- table(c(d$momloc[find], d$momloc2[find]))
    # Do not account for 0 (construction of "find" will always contain a 0)
    n <- tail(n, -1)
    id <- as.numeric(names(n))
    # Store
    d$frever_cor[d$pernum %in% id] <- n
    
    return(d)
}

# Apply function and get frever_cor using
# momloc
data.cps.2008_10 <- 
    data.cps.2008_10 |> 
    group_by(
        year, month, serial
    ) |> 
    group_modify(correct.frever) |> 
    mutate(
        frever_merge = case_when(
            # non-response but no child(ren) in hh 
            (frever == 999 & is.na(frever_cor)) ~ 0,
            # non-response but child(ren) in hh
            (frever == 999 & !is.na(frever_cor)) ~ frever_cor,
            # frever == 0 but child(ren) recorded in hh
            (frever == 0 & !is.na(frever_cor)) ~ frever_cor,
            TRUE ~ frever
        )
    ) |> 
    ungroup() |> 
    mutate(
        # Prop mother
        mother = ifelse(frever_merge == 0, 0, 1)
    )

# Combine the two data frames replacing the uncorrected years
# with the corrected ones
data.cps.cor <- 
    bind_rows(
        data.cps.2008_10,
        data.cps.2012_22
    ) |> 
    filter(
        # Consider the most populous
        # subpopulation groups
        race_eth != "NH-Other"
    )

# Compute proportions of mother and SE accounting for weights
# Complex survey design
cps_design <- 
    svydesign( 
        id = ~ hhid , 
        data = data.cps.cor , 
        weights = ~ wtfinl
    )

# Weighted childlessness by age and race
df.cps <- tibble(
    svyby( ~ mother , ~ race_eth + year + age, cps_design , svymean )
) |> 
    mutate(
        l95 = mother - (1.96 * se),
        u95 = mother + (1.96 *se),
        # Harmonize with NSFG data
        period = NA
    )



## LOAD & TIDY NSFG DATA ===================================================================

# Code taken from
# https://asdfree.com/national-survey-of-family-growth-nsfg.html

# Need to add period 2006_2010
dwld.period <- c("2006_2010", "2011_2013", "2013_2015", "2015_2017", "2017_2019")

list.data.nsfg <- lapply(dwld.period,
                       function(p) {
                           
                           # File name changed after 2006_2010
                           if (p != "2006_2010") {
                               file_name <- "_FemRespData.dat"
                           } else { 
                               file_name <- "_FemResp.dat"
                           }
                           # Url to download data
                           dat_url <- paste0(
                               "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NSFG/",
                               p,
                               file_name
                           )
                           
                           sas_url <-
                               file.path(dirname(dat_url) , 
                                         paste0("sas/",
                                                p,
                                                "_FemRespSetup.sas") )
                           
                           sas_positions <-
                               parse.SAScii( sas_url )
                           
                           sas_positions[ , 'varname' ] <-
                               tolower( sas_positions[ , 'varname' ] )
                           
                           sas_positions[ , 'column_types' ] <-
                               ifelse( sas_positions[ , 'char' ] , "c" , "d" )
                           
                           # Load data
                           nsfg_tbl <-
                               read_fwf(
                                   dat_url ,
                                   fwf_widths( 
                                       abs( sas_positions[ , 'width' ] ) , 
                                       col_names = sas_positions[ , 'varname' ] 
                                   ) ,
                                   col_types = paste0( sas_positions[ , 'column_types' ] , collapse = "" ) ,
                                   na = c( "" , "." )
                               )
                           
                           # Weights have different names in 2006_2010
                           if (p != "2006_2010") {
                               wgt_var <- paste0("wgt", p)
                           } else { 
                               wgt_var <- "wgtq1q16"
                           }
                           # Select variables of interest
                           nsfg_df <- data.frame( nsfg_tbl ) |> 
                               dplyr::select(
                                   age = age_r, hisprace2, parity, secu, sest, !!sym(wgt_var)
                               ) |> 
                               rename(
                                   # Harmonize weight var name to row bind later
                                   "wgt" = !!sym(wgt_var)
                               ) |> 
                               mutate(
                                   # Var for period
                                   period = p
                               )
                       }
)
data.nsfg <- do.call("rbind", list.data.nsfg)

# Store data
saveRDS(data.nsfg, 
        here(
            "data_private",
            "2006_2019_NSFG_FemRespData.rds"
        ), 
        compress = "xz" )

# Load data
# data.nsfg <- readRDS(
#     here(
#         "data_private",
#         "2006_2019_NSFG_FemRespData.rds"
#     )
# )

# Tidy data
data.nsfg <- 
    data.nsfg |> 
    mutate(
        ## Create 5 variables for parity (to be used with svyby())
        mother = ifelse(parity == 0, 0, 1),
        # Name abbridged race
        race_eth = case_when(
            hisprace2 == 1 ~ "Hispanic",
            hisprace2 == 2 ~ "NH-White",
            hisprace2 == 3 ~ "NH-Black",
            hisprace2 == 4 ~ "NH-Other"
        )
    ) |> 
    filter(
        # Only consider most populous subpopulations
        race_eth != "NH-Other"
    )

# Complex survey design
nsfg_design <- 
    svydesign( 
        id = ~ secu, 
        strata = ~ sest, 
        data = data.nsfg, 
        weights = ~ wgt, 
        nest = TRUE 
    )

# Weighted prop of mother by age and race
df.nsfg <- tibble(
    svyby( ~ mother , ~ age + race_eth + period, nsfg_design , svymean )
) |> 
    mutate(
        l95 = mother - (1.96 * se),
        u95 = mother + (1.96 *se),
        year = case_when(
            period == "2006_2010" ~ 2008,
            period == "2011_2013" ~ 2012,
            period == "2013_2015" ~ 2014,
            period == "2015_2017" ~ 2016,
            period == "2017_2019" ~ 2018,
        )
    )


## ===== COMBINE THE TWO TIDY DF ===============================================

df <- bind_rows(
    df.cps |> 
        mutate(
            survey = "CPS"
        ),
    df.nsfg |> 
        mutate(
            survey = "NSFG"
        )
)

    

## VISUALIZATION ===============================================================

df.cps |> 
    ggplot(aes(x = age,
               y = mother,
               ymin = l95,
               ymax = u95,
               group = year,
               col = year)) +
    facet_wrap(~ race_eth) +
    geom_line() +
    geom_point() +
    geom_pointrange() +
    theme_bw()


df |> 
    ggplot(aes(x = age,
               y = mother,
               ymin = l95,
               ymax = u95,
               group = survey,
               col = survey)) +
    facet_grid(race_eth ~ year) +
    geom_line() +
    geom_point() +
    geom_pointrange() +
    theme_bw()




    
