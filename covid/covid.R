library(tidyverse)
library(lubridate)
library(cowplot)

# inspired by https://twitter.com/_ansgar/status/1480557441012150272 and  https://bydata.github.io/nyt-corona-spiral-chart/

# get data
owid_url <- "https://github.com/owid/covid-19-data/blob/master/public/data/owid-covid-data.csv?raw=true"

country <- "United States"

covid <- read_csv(owid_url)

covid_cases <- covid %>% 
  filter(location == country) %>% 
  select(date, new_cases, new_cases_smoothed) %>% 
  arrange(date) %>% 
  # Add the dates before the 1st confirmed case
  add_row(date = as_date("2020-01-01"), new_cases = 0, new_cases_smoothed = 0,
          .before = 1) %>% 
  complete(date = seq(min(.$date), max(.$date), by = 1),
          fill = list(new_cases = 0, new_cases_smoothed = 0)) %>% 
  mutate(day_of_year = yday(date),
         year = year(date)
  )

head(covid_cases)
names(covid_cases)


a <- covid_cases %>% 
  ggplot(aes(x = date)) +
  geom_bar(aes(y = new_cases), stat = "identity") +
  geom_line(aes(y = new_cases_smoothed), color = "red") +
  scale_x_date(breaks = "month") +
  scale_y_continuous(labels = scales::label_number_si()) +
  theme_linedraw(base_size = 16) +
  theme(axis.text.x = element_text(angle =45, hjust = 1),) +
  labs(x = " ", y = "New cases",
       subtitle = "New Covid-19 cases in the US")


b <- covid_cases %>% 
  ggplot(aes(x = date)) +
  geom_bar(aes(y = new_tests), stat = "identity") +
  scale_x_date(breaks = "month") +
  scale_y_continuous(labels = scales::label_number_si()) +
  theme_linedraw(base_size = 16) +
  theme(axis.text.x = element_text(angle =45, hjust = 1),) +
  labs(x = " ", y = "New tests",
       subtitle = "Number of Covid-19 tests administered the US")

c <- covid_cases %>%
  mutate(cases_per_test = new_cases/new_tests) %>%
  ggplot(aes(x = date)) +
  geom_bar(aes(y = cases_per_test), stat = "identity") +
  scale_x_date(breaks = "month") +
  #scale_y_continuous(labels = scales::label_number_si()) +
  theme_linedraw(base_size = 16) +
  theme(axis.text.x = element_text(angle =45, hjust = 1),) +
  labs(x = "Date", y = "Cases per tests",
       subtitle = "Rate of Covid-19 positive cases per test in the US",
       caption = "Data from https://github.com/owid/covid-19-data/")

plot_grid(a,b,c, nrow= 3, rel_heights = c(1,1,1.1))
