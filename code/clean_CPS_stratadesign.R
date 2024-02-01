
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code clean CPS data and creates strata and cluster for survey()
##  following recommendations from Davern et al. (2006, 2007) paper.
##  It means scrapping a lot of info on counties and metro area. This
##  code is no longer used because defining strata had no impact on the
##  computation of SE.
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


## ===== DATA ON METRO AREA & COUNTY ===========================================

# Davern et al. (2006, 2007) paper shows that 
# using a sampling design improves the estimation of SE.
# In order to define a sampling design they build strata
# based on state, populous county and biggest metro area.

# Data with county name, code, and population
# Accessed on the 31/1/2024
page <- "https://www2.census.gov/programs-surveys/popest/datasets/2020-2022/counties/totals/co-est2022-alldata.csv"

df.county.pop <- read.csv(
    url(page),
) |> 
    dplyr::select(
        STNAME, STATE, CTYNAME, COUNTY, POPESTIMATE2020
    ) |> 
    rename(
        "state" = STNAME,
        "statefip" = STATE,
        "county" = CTYNAME,
        "countyfip" = COUNTY,
        "pop" = POPESTIMATE2020
    ) |> 
    filter(
        # Consist of States so remove
        countyfip != 0
    ) |> 
    mutate(
        county = sub("\\ County.*", "", county),
        # Harmonize the county code with CPS data
        # All need have 3 digits
        countyfip = as.character(countyfip),
        countyfip = case_when(
            nchar(countyfip) == 1 ~ paste0("00", countyfip),
            nchar(countyfip) == 2 ~ paste0("0", countyfip),
            nchar(countyfip) >= 2 ~ countyfip
        ),
        countyfip_cps = paste0(statefip, countyfip) |> as.numeric()
    ) 

counties.kept <- 
    df.county.pop |> 
    filter(
        # Davern et al. (2006, 2007) consider counties with > 100,000 pop
        # Here, we place the threshold at 135,000 to account for population growth
        # between 1990 (year used in paper) and 2020
        pop > 135000
    ) |> 
    pull(
        countyfip_cps
    )


# CPS data does not provide metro stat area name so we scraped it for
# two periods as fip codes changed. We need the names as
# Davern et al. (2006, 2007) in building their strata only account
# for the biggest 242 metro stat area. Ranking of metro stat area 
# are based on the name, not the fip so we had to get the name
# to associate them with the fip and keep the 242 biggest.

# 2004-2014
page <- read_html("https://cps.ipums.org/cps/codes/metfips_20042014_codes.shtml")

metarea.fip <-
    page %>%
    html_nodes("dt") |> 
    html_text()

metarea.name <-
    page %>%
    html_nodes("dd") |> 
    html_text()

metro.area.fip.2004_14 <-
    tibble(
        metarea = metarea.name,
        metareafip = metarea.fip
    ) |> 
    mutate(
        metarea = sub("\\,.*", "", metarea)
    )


# 2014-
page <- read_html("https://cps.ipums.org/cps/codes/metfips_2014onward_codes.shtml")

metarea.fip <-
    page %>%
    html_nodes("dt") |> 
    html_text()

metarea.name <-
    page %>%
    html_nodes("dd") |> 
    html_text()

metro.area.fip.2014_ <-
    tibble(
        metarea = metarea.name,
        metareafip = metarea.fip
    ) |> 
    mutate(
        metarea = sub("\\,.*", "", metarea)
    )

# Data on the metro stat area ranked by population
url <- "https://en.wikipedia.org/wiki/Metropolitan_statistical_area"
page = read_html(url)

# Obtain the piece of the web page that corresponds to the "wikitable" node
metro.area = html_node(page, ".wikitable")

# Convert the html table element into a data frame
big.metro.area = html_table(metro.area, fill = TRUE) |> 
    dplyr::select(
        "metroarea" = `Metropolitan statistical area`
    ) |> 
    slice(
        # Keep 242 largest (1st row is empty)
        2:243
    ) |> 
    mutate(
        metroarea = sub("\\,.*", "", metroarea)
    ) |> 
    pull(metroarea)

# Get metarea fip of the biggest according to the period
metro.area.fip.2004_14 <-
    metro.area.fip.2004_14 |> 
    filter(
        metarea %in% big.metro.area
    ) |> 
    pull(metareafip)

metro.area.fip.2014_ <-
    metro.area.fip.2014_ |> 
    filter(
        metarea %in% big.metro.area
    ) |> 
    pull(metareafip)




## ===== LOAD & TIDY CPS DATA ==================================================

# Load CPS data (extract from IPUMS)
# read_csv allows to open compressed .gz 
raw.cps <- read_csv(
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
                             hispan != 0 ~ "Hispanic")
    )

# Create Strata following Davern et al. (2006, 2007)
data.cps <-
    raw.cps |> 
    mutate(
        # Consider only county with more than 150,000 population (place bound higher than in paper)
        big_county = ifelse(county %in% counties.kept, county, 0),
        big_metarea = case_when(
            (year < 2014) & (metfips %in% metro.area.fip.2004_14) ~ metfips,
            (year >= 2014) & (metfips %in% metro.area.fip.2014_) ~ metfips,
            TRUE ~ 0
        ),
        strataid = case_when(
            big_county == 0 & big_metarea == 0 ~ statefip,
            big_metarea != 0 ~ (statefip * 10000) + big_metarea,
            big_county != 0 & big_metarea == 0 ~ (statefip * 10000) + big_county
        ),
        # Create HH following Davern et al. (2006, 2007)
        hhid = paste0(serial, year) |> as.numeric(),
        # Some strata have only one cluster, which create issue in survey()
        # later in the code. Merge these into bigger strata
        strataid = case_when(
            year == 2014 & strataid == 143223 ~ 13,
            year == 2014 & strataid == 297171 ~ 27,
            year == 2020 & strataid == 429103 ~ 39,
            year == 2022 & strataid == 407159 ~ 37,
            TRUE ~ strataid
        )
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
        p0 = ifelse(frever == 0, 1, 0)
    )

# Strata with one cluster creates problem in survey()

# data.cps.2012_22 |> 
#     group_by(
#         year, strataid, hhid
#     ) |> 
#     summarise(
#         nber = n()
#     ) |> 
#     ungroup() |> 
#     summarise(
#         .by = c(year, strataid),
#         n_psu = n()
#     ) |> 
#     filter(
#         n_psu == 1
#     )
