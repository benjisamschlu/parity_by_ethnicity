
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code loads, cleans up, and tidy the data used to estimate
##  the proportion of childless women.
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
packages <- c("tidyverse", "rvest", "here")
for(p in packages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
}



## FUNCTIONS ===================================================================



## LOAD DATA ===================================================================

# Careful about age of mother: not the same boundaries over the years
# https://cps.ipums.org/cps-action/variables/FREVER#universe_section

# read_csv allows to open compressed .gz 
data <- read_csv(
    here(
        "data_private", 
        "cps_00005.csv.gz"
        )
    ) |> 
    rename_all(.funs = tolower) |> 
    # Keep useful variables
    dplyr::select(
        year, month, serial,pernum, # identification
        momloc, momloc2, poploc, poploc2, frever,# parenthood
        wtfinl, # weights,
        statefip, # states
        age, sex, race, hispan) |> # demographic
    mutate(
        # Create abridged race (!! NEED TO BE CHECKED !!)
        race_eth = case_when(race == 100 & hispan == 0 ~ "NH-White",
                             race == 200 & hispan == 0 ~ "NH-Black",
                             !(race %in% c(100, 200)) & hispan == 0 ~ "NH-Other",
                             hispan != 0 ~ "Hispanic")
    ) |> 
    filter(
        # Survey conducted every two years (2022 has no info on frever)
        year %in% seq(2000, 2020, 2)
    )

# Wikipedia pages with States names and fip codes 
url <- "https://en.wikipedia.org/wiki/Federal_Information_Processing_Standard_state_code"
page = read_html(url)

# Obtain the piece of the web page that corresponds to the "wikitable" node
fip.table = html_node(page, ".wikitable")

# Convert the html table element into a data frame
fip.table = html_table(fip.table, fill = TRUE)

# Tidying 
fip.table <- fip.table |> 
    filter(`Alpha code` != "") |> 
    dplyr::select(name = Name, statefip = `Numeric code`)



## TIDY DATA ===================================================================

# Tidy
df <- data |> 
    filter(
        # Question on children only asked to women
        sex == 2,
        # aged 15-45 not having 999.
        # Only in latest CPS that itwees
        # are until 50 yo
        age >= 15,
        age < 45,
        # Remove when info on frever is not collected
        frever != 999
        )
# Frever has no NA

# Childlessness from variable frever
df.chldness <- df |>  
    mutate(
        # Rescale weights
        w = wtfinl/1000,
        # Create 4+ nber of minor
        childless = ifelse(frever == 0, 1, 0)
    ) |> 
    group_by(
        age, year, statefip, race_eth,
        .drop = F
    ) |> 
    summarise(
        y = sum(w[childless == 1]),
        n = sum(w)
    ) |> 
    ## Add States names
    left_join(
        fip.table,
        by = c("statefip")
    ) |> 
    ungroup()


# Store data for modeling
saveRDS(df.chldness,
        here(
            "data", 
            "df_childness.rds"
            )
        )



## CORRECTION OF 2008, 2010 ====================================================

# Use hh info and coresident child info to correct proportion childless
# for the years 2008 and 2010
df_cor <- data |> 
    filter(
        # Question on children only asked to women
        sex == 2,
        # Two years were info of coresidency wasn't use
        year %in% c(2008, 2010)
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
df.chldness.cor <- df_cor |> 
    group_by(
        year, month, serial
        ) |> 
    group_modify(correct.frever) |> 
    mutate(
        frever_merge = case_when(
            # non-response but not child(ren) in hh 
            (frever == 999 & is.na(frever_cor)) ~ 0,
            # non-response but child(ren) in hh
            (frever == 999 & !is.na(frever_cor)) ~ frever_cor,
            # frever == 0 but child(ren) recorded in hh
            (frever == 0 & !is.na(frever_cor)) ~ frever_cor,
            TRUE ~ frever
        )
    ) |> 
    filter(
        # aged 15-45 not having 999.
        # Only in latest CPS that itwees
        # are until 50 yo
        age >= 15,
        age < 45
    ) |> 
    mutate(
        # Rescale weights
        w = wtfinl/1000,
        # Create 4+ nber of minor
        childless = ifelse(frever_merge == 0, 1, 0)
    ) |> 
    group_by(
        age, year, statefip, race_eth,
        .drop = F
    ) |> 
    summarise(
        y = sum(w[childless == 1]),
        n = sum(w)
    ) |> 
    ## Add States names
    left_join(
        fip.table,
        by = c("statefip")
    ) |> 
    ungroup()

# Store data 
saveRDS(df.chldness.cor,
        here(
            "data", 
            "df_childness_cor_2008_10.rds"
        )
) 

# use nber of children from hh info in other years?
# risky, children might be out of hh
# -> only use when frever==0 but child(ren) in hh

# can frever == 999 and frever_cor>0?
