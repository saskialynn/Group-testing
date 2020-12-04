#### Preprocessing scripts
library("zoo")
library("reshape2")
source("c3aidatalake.R")
df = read_csv("countries_names.csv", col_types = list(col_character(),
                                                      col_character(),
                                                      col_character(),
                                                      col_character(),
                                                      col_double()))
states_abbreviations = read_csv("state_abbreviations.csv",
                                col_types = list(col_character(),
                                                 col_character()))

test_performed = read_csv("daily-tests-per-thousand-people-smoothed-7-day.csv")
test_performed$Date = as.Date(test_performed$Date, format = "%m/%d/%y")
mymonths <- c("jan","feb","mar",
              "april","may","june",
              "july","aug","sept",
              "oct","nov","dec")



fetch_all_prevalence <- function(country, region,  date){
  date0 <- date
  if (date > as.Date("2020-11-10")){
    date0 = as.Date("2020-11-10")
  }
  database_lookup = df[which((df$country_name == country) & (df$region == region)),]
  pop = database_lookup$pop
  database = database_lookup$database
  name = database_lookup$country
  abbr = ifelse((country == "United States") &(region %in% states_abbreviations$state) , (states_abbreviations %>% filter(state == region))$state_abb, "")
  if ((region == "American Samoa") ||(region== "Puerto Rico")){
      abbr = region
  }
  

  #### Fetch cases on that day
  casecounts <- evalmetrics(
    "outbreaklocation",
    list(
      spec = list(
        ids = list(name),
        expressions = list(paste0(database, "_ConfirmedCases")),
        start = date0 - 31,
        end = date0,
        interval = "DAY"
      )
    )
  )
  casecounts = casecounts %>% mutate(new_cases = data - lag(data, default = 1)) %>% filter(dates >= date0 - 21)
  casecounts = casecounts %>% mutate(smoothed_cases=rollapply(new_cases, 7, mean, align='right',fill=NA))  ### thats the new cases counts smoothed over
  ### Let's now compute the prevalence: active cases over the last 14 days
  prevalence = sum(casecounts %>% dplyr::filter(dates > date0 -14) %>% dplyr::select(smoothed_cases))/ pop
  ### Thats the global prevalence in the population
  ### Now let's fetch the prevalence among the tested
  return(list(prev= prevalence))
}
