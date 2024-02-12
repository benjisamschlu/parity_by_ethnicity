
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



## ===== STATE NAMES & FIP =====================================================

# Get data set with state names and fip codes

# # Wikipedia pages with States names and fip codes
url <- "https://en.wikipedia.org/wiki/Federal_Information_Processing_Standard_state_code"
page = read_html(url)

# Obtain the piece of the web page that corresponds to the "wikitable" node
df.fip <- 
    page |> 
    html_node(".wikitable") |> 
    # Convert the html table element into a data frame
    html_table(fill = TRUE) |> 
    # Tidying
    filter(`Alpha code` != "") |>
    dplyr::select(
        state = Name, statefip = `Numeric code`
        )

# Store data on states
saveRDS(
    df.fip,
    here(
        "data",
        "df_fip.rds"
        )
    )

# Load data on fip code and state names
# df.fip <- readRDS(
#     here(
#         "data",
#         "df_fip.rds"
#     )
# )


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
        ) |> 
    # Add States names
    left_join(
        df.fip,
        by = c("statefip")
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
        p = ifelse(frever == 0, 0, 1),
        # Categorical variable for counts with surey()
        mother = ifelse(frever == 0, "No", "Yes")
    )

# Correct the data for 2008 and 2010
# Utilization of information from the household rosters both to allocate data 
# in cases of non-response, and to validate responses.
data.cps.2008_10 <- 
    data.cps |> 
    filter(
        # Two years were info of coresidency wasn't use
        year %in% c(2008, 2010)
    ) |> 
    mutate(
        # Never frever based on hh roster
        frever_cor = 0
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
            # non-response so use frever_cor
            frever == 999 ~ frever_cor,
            # frever == 0 but possible that child(ren) recorded in hh
            # otherwise frever_cor is also = 0
            frever == 0 ~ frever_cor,
            # frever is > 0 so use it (frever_cor might miss children
            # that left the hh)
            TRUE ~ frever
            )
        ) |> 
    ungroup() |> 
    mutate(
        # Prop mother
        p = ifelse(frever_merge == 0, 0, 1),
        # Categorical variable for counts with surey()
        mother = ifelse(frever_merge == 0, "No", "Yes")
    ) |> 
    filter(
        # Years 2008 and 2010 only interviewed women aged 15-44
        age >= 15, 
        age < 45
    )

# Combine the two data frames replacing the uncorrected years
# with the corrected ones
data.cps.cor <- 
    bind_rows(
        data.cps.2008_10,
        data.cps.2012_22
    ) 

# Compute proportions of mother and SE accounting for weights
# Complex survey design
cps_design <- 
    svydesign( 
        id = ~ hhid , 
        data = data.cps.cor , 
        weights = ~ wtfinl
    )

# Weighted motherhood prop. by age at US level
df.cps.p.us <- tibble(
    svyby( ~ p , ~ year + age, cps_design , svymean )
) |> 
    mutate(
        l95 = p - (1.96 * se),
        u95 = p + (1.96 *se)
    ) |> 
    arrange(
        year, age
    )

# Weighted count of mother by age at US level
df.cps.y.us <- tibble(
    svyby( ~ mother , ~ year + age, cps_design , svytotal)
) |> 
    mutate(
        n = motherYes + `motherNo`
    ) |> 
    rename(
        "y" = motherYes
    ) |> 
    arrange(
        year, age
    ) |> 
    filter(
        # Remove cases where y==0 and y==n:
        # Unobserved & log(p / (1-p)) is undefined
        y != 0,
        y != n
    )

# Combine p, y, and n into one df at US level
df.cps.us <-
    df.cps.y.us |> 
    dplyr::select(
        year, age, y, n 
    ) |> 
    left_join(
        df.cps.p.us,
        by = c("year", "age")
    ) 


# Weighted motherhood prop. by age at RACE level
df.cps.p.race <- tibble(
    svyby( ~ p , ~ race_eth + year + age, cps_design , svymean )
) |> 
    mutate(
        l95 = p - (1.96 * se),
        u95 = p + (1.96 *se)
    ) |> 
    arrange(
        race_eth, year, age
    )

# Weighted count of mother by age at RACE level
df.cps.y.race <- tibble(
    svyby( ~ mother , ~ race_eth + year + age, cps_design , svytotal)
    ) |> 
    mutate(
        n = motherYes + `motherNo`
    ) |> 
    rename(
        "y" = motherYes
    ) |> 
    arrange(
        race_eth, year, age
    ) |> 
    filter(
        # Remove cases where y==0 and y==n:
        # Unobserved & log(p / (1-p)) is undefined
        y != 0,
        y != n
    )

# Combine p, y, and n into one df at RACE level
df.cps.race <-
    df.cps.y.race |> 
    dplyr::select(
        race_eth, year, age, y, n 
        ) |> 
    left_join(
        df.cps.p.race,
        by = c("race_eth", "year", "age")
    ) 

# Weighted motherhood prop. by age at STATE level

# !!! NOT ENOUGH MEMORY !!!
df.cps.p.state <- tibble(
    svyby( ~ p , ~ state + year + age, cps_design , svymean )
) |> 
    mutate(
        l95 = p - (1.96 * se),
        u95 = p + (1.96 *se)
    ) |> 
    arrange(
        state, year, age
    )

# Weighted count of mother by age at STATE level
df.cps.y.state <- tibble(
    svyby( ~ mother , ~ state + year + age, cps_design , svytotal)
) |> 
    mutate(
        n = motherYes + `motherNo`
    ) |> 
    rename(
        "y" = motherYes
    ) |> 
    arrange(
        state, year, age
    )

# Combine p, y, and n into one df at STATE level
df.cps.state <-
    df.cps.y.state |> 
    dplyr::select(
        state, year, age, y, n 
    ) |> 
    left_join(
        df.cps.p.state,
        by = c("state", "year", "age")
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
                           # Different age bounds according to the period
                           if (p %in% c("2006_2010", "2011_2013", "2013_2015")) {
                               age_bnd <- 45
                           } else { 
                               age_bnd <- 50
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
                               filter(
                                   # Remove age that shouldn't be surveyed
                                   age < age_bnd
                                   
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
        # Prop mother
        p = ifelse(parity == 0, 0, 1),
        # Categorical variable for counts with surey()
        mother = ifelse(parity == 0, "No", "Yes"),
        # Name abbridged race
        race_eth = case_when(
            hisprace2 == 1 ~ "Hispanic",
            hisprace2 == 2 ~ "NH-White",
            hisprace2 == 3 ~ "NH-Black",
            hisprace2 == 4 ~ "NH-Other"
        )
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

# Weighted prop of mother by age at US level
df.nsfg.p.us <- tibble(
    svyby( ~ p , ~ period + age, nsfg_design , svymean )
) |> 
    mutate(
        l95 = p - (1.96 * se),
        u95 = p + (1.96 *se)
    ) |> 
    arrange(
        period, age
    )

# Weighted count of mother by age at US level
df.nsfg.y.us <- tibble(
    svyby( ~ mother , ~ period + age, nsfg_design , svytotal)
) |> 
    mutate(
        n = motherYes + `motherNo`
    ) |> 
    rename(
        "y" = motherYes
    ) |> 
    arrange(
        period, age
    ) |> 
    filter(
        # Remove cases where y==0 and y==n:
        # Unobserved & log(p / (1-p)) is undefined
        y != 0,
        y != n
    )

# Combine p, y, and n into one df at US level
df.nsfg.us <-
    df.nsfg.y.us |> 
    dplyr::select(
        period, age, y, n 
    ) |> 
    left_join(
        df.nsfg.p.us,
        by = c("period", "age")
    ) |> 
    mutate(
        year = case_when(
            period == "2006_2010" ~ 2008,
            period == "2011_2013" ~ 2012,
            period == "2013_2015" ~ 2014,
            period == "2015_2017" ~ 2016,
            period == "2017_2019" ~ 2018,
        )
    )

# Weighted prop of mother by age at RACE level
df.nsfg.p.race <- tibble(
    svyby( ~ p , ~ race_eth + period + age, nsfg_design , svymean )
) |> 
    mutate(
        l95 = p - (1.96 * se),
        u95 = p + (1.96 *se)
    ) |> 
    arrange(
        race_eth, period, age
    )

# Weighted count of mother by age at RACE level
df.nsfg.y.race <- tibble(
    svyby( ~ mother , ~ race_eth + period + age, nsfg_design , svytotal)
) |> 
    mutate(
        n = motherYes + `motherNo`
    ) |> 
    rename(
        "y" = motherYes
    ) |> 
    arrange(
        race_eth, period, age
    ) |> 
    filter(
        # Remove cases where y==0 and y==n:
        # Unobserved & log(p / (1-p)) is undefined
        y != 0,
        y != n
    )

# Combine p, y, and n into one df at RACE level
df.nsfg.race <-
    df.nsfg.y.race |> 
    dplyr::select(
        race_eth, period, age, y, n 
    ) |> 
    left_join(
        df.nsfg.p.race,
        by = c("race_eth", "period", "age")
    ) |> 
    mutate(
        year = case_when(
            period == "2006_2010" ~ 2008,
            period == "2011_2013" ~ 2012,
            period == "2013_2015" ~ 2014,
            period == "2015_2017" ~ 2016,
            period == "2017_2019" ~ 2018,
        )
    )

# NSFG data does not allow to get estimates at the state level


## ===== COMBINE THE TWO TIDY DF ===============================================

# US level
df.us <- bind_rows(
    df.cps.us |> 
        mutate(
            survey = "CPS"
        ),
    df.nsfg.us |> 
        mutate(
            survey = "NSFG"
        )
)

# Store the tidy data
saveRDS(
    df.us,
    here(
        "data",
        "df_us.rds"
    ),
    compress = "xz"
)

# Race level
df.race <- bind_rows(
    df.cps.race |> 
        mutate(
            survey = "CPS"
        ),
    df.nsfg.race |> 
        mutate(
            survey = "NSFG"
        )
)

# Store the tidy data
saveRDS(
    df.race,
    here(
        "data",
        "df_race.rds"
    ),
    compress = "xz"
)

# State level !!!!

    

