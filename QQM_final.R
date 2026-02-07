############################################
# W25 SEMESTER PROJECT: ADVANCED METHODS II
# Budgeting for Safety:  Linking Defence and Climate Spending 
# in the EU to Security Perceptions in Austria 
############################################


install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("sf") 
install.packages("eurostat")
install.packages("modelsummary")
install.packages("patchwork")

packages <- c("tidyverse","sf", "leaflet", "ggplot2", "eurostat", "forecast", "modelsummary",
              "patchwork", "rnaturalearth","rnaturalearthdata", "dplyr", "htmltools", "spdep")
sapply(packages, library, character.only = T)


rm(list = ls())
gc()



############################################
# 1) DATA WRANGLING - MAIN DATASET 
############################################

setwd("/Users/filipareis/Desktop/QQM_II")

exp_gdp <- read.csv("/Users/filipareis/Desktop/QQM_II/data/expenditure.csv")

dataset_long <- exp_gdp %>%
  select(cofog99, geo, TIME_PERIOD, OBS_VALUE) %>%
  drop_na() %>%
  mutate(
    TIME_PERIOD = as.integer(TIME_PERIOD), #------------------------------------------> why as integer?
    geo = recode(geo, "European Union - 27 countries (from 2020)" = "EU-27") 
  ) %>%
  filter(TIME_PERIOD >= 1998, TIME_PERIOD <= 2023) %>%
  rename(
    category = cofog99,
    country  = geo,
    year     = TIME_PERIOD,
    value    = OBS_VALUE
  )

dataset <- dataset_long %>%
  pivot_wider(names_from = category, values_from = value) %>%
  mutate(year = as.integer(year)) %>% #------------------------------------------> was ist das?
  arrange(country, year) %>%
  mutate(
    new_env = `Environmental protection` - `Waste management` - `Waste water management`
  )


############################################
# 2) FIT + FORECAST PER COUNTRY 
############################################

fit_forecast_country_nowindow <- function(df, country_name, var,
                                          train_end_target = 2019, 
                                          arima_order = c(1,1,2),
                                          h = 4,
                                          use_auto_arima = FALSE,
                                          min_train_years = 10) {  
  
  d_raw <- df %>%
    filter(country == country_name) %>%
    arrange(year) %>%
    select(year, value = all_of(var))
  
  actual_all_years <- tibble(
    country  = country_name,
    variable = var,
    year     = d_raw$year,
    actual   = d_raw$value
  )
  
  d <- d_raw %>% drop_na(value)
  
  if (nrow(d) < min_train_years) {
    return(list(
      actual = actual_all_years,
      forecast = tibble(),
      compare = tibble(),
      status = tibble(country = country_name, variable = var, status = "SKIP: too few non-missing observations")
    ))
  }
  
  d <- d %>%
    complete(year = full_seq(min(year):max(year), 1)) %>%
    drop_na(value)
  
  if (nrow(d) < min_train_years) {
    return(list(
      actual = actual_all_years,
      forecast = tibble(),
      compare = tibble(),
      status = tibble(country = country_name, variable = var, status = "SKIP: too few after completing years")
    ))
  }
  
  years  <- d$year
  values <- d$value
  
  train_idx <- years <= train_end_target
  if (sum(train_idx) < min_train_years) {
    return(list(
      actual = tibble(country = country_name, variable = var, year = years, actual = values),
      forecast = tibble(),
      compare = tibble(),
      status = tibble(country = country_name, variable = var, status = "SKIP: not enough years <= train_end_target")
    ))
  }
  
  train_years  <- years[train_idx]
  train_values <- values[train_idx]
  
  train_start <- min(train_years)
  train_ts <- ts(train_values, start = train_start, frequency = 1)
  
  fit <- tryCatch(
    {
      if (use_auto_arima) {
        auto.arima(train_ts, seasonal = FALSE, stepwise = TRUE, approximation = FALSE)
      } else {
        arima(train_ts, order = arima_order)
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(
      actual = tibble(country = country_name, variable = var, year = years, actual = values),
      forecast = tibble(),
      compare = tibble(),
      status = tibble(country = country_name, variable = var, status = "FAIL: ARIMA fit error") # status, related to na and missing obsv if years can be skipped
    ))
  }
  
  fc <- forecast(fit, h = h)
  
  last_train_year <- max(train_years)
  fc_years <- (last_train_year + 1):(last_train_year + h)
  
  forecast_tbl <- tibble(
    country  = country_name,
    variable = var,
    year     = fc_years,
    forecast = as.numeric(fc$mean),
    lo80     = as.numeric(fc$lower[, "80%"]),
    hi80     = as.numeric(fc$upper[, "80%"]),
    lo95     = as.numeric(fc$lower[, "95%"]),
    hi95     = as.numeric(fc$upper[, "95%"]),
    train_end_used = last_train_year
  )
  
  actual_h_tbl <- tibble(year = years, actual = values) %>%
    filter(year %in% fc_years)
  
  compare_tbl <- forecast_tbl %>%
    left_join(actual_h_tbl, by = "year") %>%
    mutate(
      error     = actual - forecast,
      abs_error = abs(error)
      )
  
  list(
    actual = tibble(country = country_name, variable = var, year = years, actual = values),
    forecast = forecast_tbl,
    compare = compare_tbl,
    status = tibble(country = country_name, variable = var, status = "OK")
  )
}

run_all_countries_nowindow <- function(df, var, arima_order,
                                       train_end_target = 2019,
                                       h = 4,
                                       use_auto_arima = FALSE,
                                       min_train_years = 10) {
  
  countries <- unique(df$country)
  
  res <- lapply(countries, function(cn) {
    fit_forecast_country_nowindow(
      df = df,
      country_name = cn,
      var = var,
      train_end_target = train_end_target,
      arima_order = arima_order,
      h = h,
      use_auto_arima = use_auto_arima,
      min_train_years = min_train_years)})
  
  list(
    actual   = bind_rows(lapply(res, `[[`, "actual")),
    forecast = bind_rows(lapply(res, `[[`, "forecast")),
    compare  = bind_rows(lapply(res, `[[`, "compare")),
    status   = bind_rows(lapply(res, `[[`, "status")))
}


############################################
# 3) RUN MODELS
############################################

def_out <- run_all_countries_nowindow(
  df = dataset,
  var = "Defence",
  arima_order = c(1,1,2),
  train_end_target = 2019,
  h = 4,
  use_auto_arima = FALSE)

env_out <- run_all_countries_nowindow(
  df = dataset,
  var = "new_env",
  arima_order = c(1,2,2),
  train_end_target = 2019,
  h = 4,
  use_auto_arima = FALSE)

# Diagnostics (skips/failures)
print(def_out$status %>% count(status) %>% arrange(desc(n)))
print(env_out$status %>% count(status) %>% arrange(desc(n)))

# Save forecast tables and model status
write.csv(def_out$compare, "defence_forecast_table.csv", row.names = FALSE)
write.csv(env_out$compare, "environment_forecast_table.csv", row.names = FALSE)
write.csv(def_out$status, "defence_model_status.csv", row.names = FALSE)
write.csv(env_out$status, "environment_model_status.csv", row.names = FALSE)



######################
# 4) FACET PLOTS 
######################

make_facet_plot <- function(actual_df, forecast_df,
                            title, subtitle,
                            interval = "80",
                            facet_cols = 6,
                            y_limits = NULL) {
  
  lo <- if (interval == "95") "lo95" else "lo80"
  hi <- if (interval == "95") "hi95" else "hi80"
  
  # Alphabetical facet order
  actual_df <- actual_df %>%
    mutate(country = factor(country, levels = sort(unique(country))))
  forecast_df <- forecast_df %>%
    mutate(country = factor(country, levels = levels(actual_df$country)))
  
  # Fixed split at 2019 (your requested training cutoff)
  train_end_line <- 2019
  
  # Forecast shading range (if forecasts exist)
  if (nrow(forecast_df) > 0) {
    fc_min <- min(forecast_df$year, na.rm = TRUE)
    fc_max <- max(forecast_df$year, na.rm = TRUE)
  } else {
    fc_min <- train_end_line + 1
    fc_max <- train_end_line + 5
  }
  
  p <- ggplot() +
    # Forecast-period shading
    annotate(
      "rect",
      xmin = fc_min - 0.5, xmax = fc_max + 0.5,
      ymin = -Inf, ymax = Inf,
      fill = "#F4F7FB", alpha = 0.9
    ) +
    # Confidence interval ribbon
    geom_ribbon(
      data = forecast_df,
      aes(x = year, ymin = .data[[lo]], ymax = .data[[hi]], group = 1),
      fill = "#BFD7FF", alpha = 0.35
    ) +
    # Forecast line
    geom_line(
      data = forecast_df,
      aes(x = year, y = forecast, group = 1),
      linewidth = 0.5,
      color = "#1F5AA6"
    ) +
    # Actual line
    geom_line(
      data = actual_df,
      aes(x = year, y = actual, group = 1),
      linewidth = 0.5,
      color = "#1A1A1A"
    ) +
    # Train/forecast split line
    geom_vline(
      xintercept = train_end_line + 0.5,
      linewidth = 0.3,
      linetype = "dashed",
      color = "#6E6E6E"
    ) +
    facet_wrap(~ country, ncol = facet_cols, scales = "fixed") +
    scale_x_continuous(breaks = c(1998, 2010, 2019, 2023)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "% of GDP"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 14, color = "#111111"),
      plot.subtitle = element_text(size = 10, color = "#4A4A4A"),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "#FFFFFF", color = "#DDDDDD"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Fixed y-limits across facets (recommended when you choose "fixed")
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  p
}

# Choose y-limits for comparability (adjust as needed)
y_def <- c(0, 4)
y_env <- c(0, 2)

p_def <- make_facet_plot(
  actual_df = def_out$actual,
  forecast_df = def_out$forecast,
  title = "Defence: Actual vs Forecast (ARIMA per country)",
  subtitle = "Training: years ≤ 2019 where available | Forecast: next 5 years (80% interval)",
  interval = "80",
  facet_cols = 6,
  y_limits = y_def
)

p_env <- make_facet_plot(
  actual_df = env_out$actual,
  forecast_df = env_out$forecast,
  title = "Environment (excl. waste & water): Actual vs Forecast",
  subtitle = "Training: years ≤ 2019 where available | Forecast: next 5 years (80% interval)",
  interval = "80",
  facet_cols = 6,
  y_limits = y_env
)

print(p_def)
print(p_env)

ggsave("defence_facet_forecasts.png", p_def, width = 16, height = 10, dpi = 300)
ggsave("environment_facet_forecasts.png", p_env, width = 16, height = 10, dpi = 300)


############################################
# 5) ACCURACY SUMMARY (MAE/RMSE) OVER AVAILABLE ACTUALS IN FORECAST YEARS
############################################

accuracy_by_country <- function(tbl) {
  tbl %>%
    filter(!is.na(actual)) %>%
    group_by(country, variable) %>%
    summarise(
      n = n(),
      MAE = mean(abs_error, na.rm = TRUE),
      RMSE = sqrt(mean(error^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(variable, RMSE)
}

def_accuracy <- accuracy_by_country(def_out$compare)
env_accuracy <- accuracy_by_country(env_out$compare)

print(def_accuracy)
print(env_accuracy)

write.csv(def_accuracy, "defence_accuracy_summary.csv", row.names = FALSE)
write.csv(env_accuracy, "environment_accuracy_summary.csv", row.names = FALSE)


############################################
# 6) DATA WRANGLING FOR OLS
############################################

def <-  read.csv("/Users/filipareis/Desktop/QQM_II/defence_forecast_table.csv")
env <-  read.csv("/Users/filipareis/Desktop/QQM_II/environment_forecast_table.csv")
dist <- read.csv("/Users/filipareis/Desktop/QQM_II/data/dist_cepii.csv", sep = ";", header = TRUE)
debt <- read.csv("/Users/filipareis/Desktop/QQM_II/debt.csv")
energy <- read.csv("/Users/filipareis/Desktop/QQM_II/energy.csv")

df1 <- def  %>%
  filter (year == 2023) %>%
  mutate(diff_def= actual - forecast) %>%
  select(country, diff_def)

df2<- env  %>%
  filter(year == 2023) %>%
  mutate(diff_env= actual - forecast) %>%
  select(country, diff_env)

df3 <- dist %>%
  drop_na() %>% 
  select(iso_o, iso_d, distcap) %>% 
  filter(iso_o == "RUS") %>% 
  rename(country = iso_d) %>%
  mutate(country = recode(country,
                          "AUT" = "Austria",
                          "BEL" = "Belgium",
                          "HRV" = "Croatia",
                          "CYP" = "Cyprus",
                          "CZE" = "Czechia",
                          "DNK" = "Denmark",
                          "EST" = "Estonia",
                          "FIN" = "Finland",
                          "FRA" = "France",
                          "DEU" = "Germany",
                          "GRC" = "Greece",
                          "HUN" = "Hungary",
                          "ISL" = "Iceland",
                          "IRL" = "Ireland",
                          "ITA" = "Italy",
                          "LVA" = "Latvia",
                          "LTU" = "Lithuania",
                          "LUX" = "Luxembourg",
                          "MLT" = "Malta",
                          "NLD" = "Netherlands",
                          "NOR" = "Norway",
                          "POL" = "Poland",
                          "PRT" = "Portugal",
                          "ROM" = "Romania",
                          "SVK" = "Slovakia",
                          "SVN" = "Slovenia",
                          "ESP" = "Spain",
                          "SWE" = "Sweden",
                          "CHE" = "Switzerland")) %>%
           filter(nchar(country) > 3)%>%
  mutate(distcap = as.numeric(gsub(",", ".", distcap))) %>%
  select(country, distcap)


df4 <- debt %>%
  select(geo, TIME_PERIOD, OBS_VALUE) %>%  
  drop_na() %>%                                     
  filter(TIME_PERIOD == 2023) %>%
  rename(
    country = geo,
    year = TIME_PERIOD,
    debt = OBS_VALUE) %>%
  group_by(country) %>%
  mutate(country = recode(country,
                          "European Union - 27 countries (from 2020)" = "EU-27"))

df5 <- energy %>%
  select(geo, TIME_PERIOD, OBS_VALUE) %>%  
  drop_na() %>%                                     
  filter(TIME_PERIOD == 2023) %>%
  rename(
    country = geo,
    year = TIME_PERIOD,
    energy = OBS_VALUE) %>%
  group_by(country) %>%
  mutate(country = recode(country,
                          "European Union - 27 countries (from 2020)" = "EU-27"))

countries <- dataset_long %>%
  distinct(country) %>%
  pull(country)

nato <- tibble(country = countries,year = 2023) %>%
  mutate(nato = if_else (country %in% c("Austria","Cyprus","Ireland","Malta","Sweden","Finland","Switzerland"),0,1))

df_ols <- df1 %>%
  full_join(df2, by = "country") %>%
  full_join(df3, by = "country") %>%
  full_join(nato, by = "country") %>%
  full_join(df4, by = "country") %>%
  full_join(df5, by = "country")


############################################
# 7) OLS
############################################

ols_env <- lm(diff_env ~ diff_def + log(distcap) + nato + debt + energy, data= df_ols)
summary(ols_env)
ols_def <- lm(diff_def ~ diff_env + log(distcap) + nato + debt + energy, data= df_ols)
summary(ols_def)

modelsummary(list("Model 1: Environmental Spending" = ols_env, 
                  "Model 2: Defense Spending" = ols_def),
             stars = TRUE,
             gof_omit = "Adj|F|Log|AIC|BIC",
             output = "html") 

models <- list(
  "Model 1" = ols_env,
  "Model 2" = ols_def
)

var_labels <- c(
  "(Intercept)" = "Intercept",
  "diff_def"    = "Change in Defense Spending",
  "diff_env"    = "Change in Environmental Spending",
  "log distcap"     = "Distance from Russia",
  "nato"        = "NATO Membership",
  "debt"        = "Debt-to-GDP",
  "energy"      = "Energy Intensity"
)



#ROBBUSTNES

library(lmtest)

#heterosckedascity
bptest(ols_env)
bptest(ols_def)

#multicoliniarity
library(car)

vif(ols_env)
vif(ols_def)


############################################
# 8) MAPS AND VISUALS 
############################################


#### EU: MILITARY VS CLIMATE EXPENDITURES

plots <- lapply(unique(dataset$country), function(c) {
  
    sub_df <- subset(dataset, country == c)
  
    ggplot(sub_df, aes(x = year)) +
    geom_line(aes(y = `Defence`, color = "Defense")) +
    geom_line(aes(y = `new_env`, color = "Environmental Protection without Waste Management")) +
    geom_line(aes(y = `Environmental protection`, color = "Environmental Protection")) +
    labs(
      #without a main title, just the country as title, and no axis names
      title = paste(unique(sub_df$country)),
      x = NULL,
      y = NULL,
      color = NULL) +
    
    #fix the x axis so only shows 3 dates
    scale_x_continuous(breaks = c(1998, 2010, 2023)) +
    
    #fix y scale
    scale_y_continuous(limits = c(0, 4))+
    
    scale_color_manual(values = c(
      "Defense" = "#3A488AFF",
      "Environmental Protection without Waste Management" = "#DABD61FF",
      "Environmental Protection" = "#D95F30FF")) +
    
    theme_minimal(base_size = 5) +
    
    #also no legend
    theme(
      plot.title = element_text(, size = 5, hjust = 0.5),
      legend.text = element_text(size = 8),
      legend.position = "bottom")
})

all_gdp <- wrap_plots(plots, ncol = 7) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Military vs. Climate Spending",
    subtitle = "in Europe 1998-2023",
    theme = theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "gray30", hjust = 0.5),
      legend.text = element_text(size = 8),
      legend.position = "bottom")
  )

all_gdp

ggsave("all_gdp.png", plot = all_gdp, width = 10, height = 10, dpi = 300)


#### ENVIRONMENT BUDGET BREAKDOWN

#EUROPE
sub_env1 <- subset(dataset, country == "EU-27")

eu <- ggplot(sub_env1, aes(x = year)) +
  geom_line(aes(y = `Environmental protection`, color = "Environmental Protection")) +
  geom_line(aes(y = `new_env`, color = "Environmental Protection (without WWM)")) +
  geom_line(aes(y = `Waste management`, color = "Waste Management")) +
  geom_line(aes(y = `Waste water management`, color = "Water Waste Management")) +
  geom_line(aes(y = `Pollution abatement`, color = "Pollution Abatment")) +
  geom_line(aes(y = `Protection of biodiversity and landscape`, color = "Biodiversidy & Lanscape Protection")) +
  geom_line(aes(y = `R&D Environmental protection`, color = "R & D")) +
  geom_line(aes(y = `Environmental protection n.e.c.`, color = "Other")) +
  labs(
    title = "EU: Climate Spending",
    subtitle = "Breakdown by Category 1998–2023",
    x = "Year",
    y = "Spending (% of GDP)",
    color = "Key"
  ) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2020)) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
  scale_color_manual(
    values = c(
      "Environmental Protection" = "#063B41FF",
      "Environmental Protection (without WWM)" = "#009392FF",
      "Waste Management" = "#75D8D5FF",
      "Water Waste Management" = "#9CCB86FF",
      "Pollution Abatment" = "#E9E29CFF",
      "Biodiversidy & Lanscape Protection" = "#EEB479FF",
      "R & D" = "#E88471FF",
      "Other" = "#CF597EFF"
    ),
    breaks = c(
      "Environmental Protection",
      "Environmental Protection (without WWM)",
      "Waste Management",
      "Water Waste Management",
      "Pollution Abatment",
      "Biodiversidy & Lanscape Protection",
      "R & D",
      "Other"
    )
  ) +
  theme_void(base_size = 10) + 
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, face = "bold", size = 10, margin = margin(r = 10)),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

eu

#AUSTRIA
sub_df2 <- subset(dataset, country == "AUSTRIA")

at_env <- ggplot(sub_df2, aes(x = year)) +
  geom_line(aes(y = `Environmental protection`, color = "Environmental Protection")) +
  geom_line(aes(y = `new_env`, color = "Environmental Protection (without WWM)")) +
  geom_line(aes(y = `Waste management`, color = "Waste Management")) +
  geom_line(aes(y = `Waste water management`, color = "Water Waste Management")) +
  geom_line(aes(y = `Pollution abatement`, color = "Pollution Abatment")) +
  geom_line(aes(y = `Protection of biodiversity and landscape`, color = "Biodiversidy & Lanscape Protection")) +
  geom_line(aes(y = `R&D Environmental protection`, color = "R & D")) +
  geom_line(aes(y = `Environmental protection n.e.c.`, color = "Other")) +
  labs(
    title = "Austria: Climate Spending",
    subtitle = "Breakdown by Category 1998–2023",
    x = "Year",
    y = "Spending (% of GDP)",
    color = "Key"
  ) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2020)) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
  scale_color_manual(
    values = c(
      "Environmental Protection" = "#063B41FF",
      "Environmental Protection (without WWM)" = "#009392FF",
      "Waste Management" = "#75D8D5FF",
      "Water Waste Management" = "#9CCB86FF",
      "Pollution Abatment" = "#E9E29CFF",
      "Biodiversidy & Lanscape Protection" = "#EEB479FF",
      "R & D" = "#E88471FF",
      "Other" = "#CF597EFF"
    ),
    breaks = c(
      "Environmental Protection",
      "Environmental Protection (without WWM)",
      "Waste Management",
      "Water Waste Management",
      "Pollution Abatment",
      "Biodiversidy & Lanscape Protection",
      "R & D",
      "Other"
    )
  ) +
  theme_void(base_size = 10) + 
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, face = "bold", size = 10, margin = margin(r = 10)),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

at_env


ggsave("eu_env.png", plot = eu, width = 13, height = 10, dpi = 300)
ggsave("at_env.png", plot = eu, width = 13, height = 10, dpi = 300)



#### MILITARY BUDGET BREAKDOWN


#EUROPE
sub_df3 <- subset(dataset, country == "EU-27")

eu_def <- ggplot(sub_df3, aes(x = year)) +
  geom_line(aes(y = Defence, color = "Defense")) +
  geom_line(aes(y = `Military defence`, color = "Military")) +
  geom_line(aes(y = `Civil defence`, color = "Civil Protection")) +
  geom_line(aes(y = `Foreign military aid`, color = "Foreign Military Aid")) +
  geom_line(aes(y = `R&D Defence`, color = "R & D")) +
  geom_line(aes(y = `Defence n.e.c.`, color = "Other")) +
  guides(size = "none") +
  labs(
    title = "EU-27: Military Spending",
    subtitle = "Breakdown by Category 1998–2023",
    x = "Year",
    y = "Spending (% of GDP)",
    color = "Key"
  ) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2020)) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
  scale_color_manual(
    values = c(
      "Defense" = "#65323EFF",
      "Military" = "#FE7F9DFF", 
      "Civil Protection" = "#FCA315FF",
      "Foreign Military Aid" = "#1BB6AFFF",
      "R & D" = "#088BBEFF",
      "Other" = "#BF9BDDFF"),
    breaks = c("Defense","Military","Civil Protection","Foreign Military Aid","R & D","Other")
  ) +
  theme_void(base_size = 15) + 
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, face = "bold", size = 10, margin = margin(r = 10)),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
eu_def

#AUSTRIA
sub_df4 <- subset(dataset, country == "Austria")

at_def <- ggplot(sub_df4, aes(x = year)) +
  geom_line(aes(y = Defence, color = "Defense")) +
  geom_line(aes(y = `Military defence`, color = "Military")) +
  geom_line(aes(y = `Civil defence`, color = "Civil Protection")) +
  geom_line(aes(y = `Foreign military aid`, color = "Foreign Military Aid")) +
  geom_line(aes(y = `R&D Defence`, color = "R & D")) +
  geom_line(aes(y = `Defence n.e.c.`, color = "Other")) +
  guides(size = "none") +
  labs(
    title = "Austria: Military Spending",
    subtitle = "Breakdown by Category 1998–2023",
    x = "Year",
    y = "Spending (% of GDP)",
    color = "Key"
  ) +
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2020)) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
  scale_color_manual(
    values = c(
      "Defense" = "#65323EFF",
      "Military" = "#FE7F9DFF", 
      "Civil Protection" = "#FCA315FF",
      "Foreign Military Aid" = "#1BB6AFFF",
      "R & D" = "#088BBEFF",
      "Other" = "#BF9BDDFF"),
    breaks = c("Defense","Military","Civil Protection","Foreign Military Aid","R & D","Other")
  ) +
  theme_void(base_size = 15) + 
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, face = "bold", size = 10, margin = margin(r = 10)),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
at_def


ggsave("eu_def.png", plot = at, width = 13, height = 10, dpi = 300)
ggsave("at_def.png", plot = at, width = 13, height = 10, dpi = 300)



#### MAP EU: CHANGE IN EXPENDITURES 


#DATA WRANGLING

diffs <- dataset %>%
  filter(year %in% c(2013, 2023)) %>%
  select(country, year, new_env, Defence) %>%
  pivot_wider(
    names_from = year,
    values_from = c(new_env, Defence)
  ) %>%
  mutate(
    env_dif     = new_env_2023 - new_env_2013,
    defence_dif = Defence_2023 - Defence_2013
  ) %>%
  select(country, env_dif, defence_dif)

europe_sf <- get_eurostat_geospatial(
  resolution = "60",  
  nuts_level = 0)

diffs_clean <- diffs %>%
  mutate(
    CNTR_CODE = case_when(
      country == "Austria"         ~ "AT",
      country == "Belgium"         ~ "BE",
      country == "Croatia"         ~ "HR",
      country == "Cyprus"          ~ "CY",
      country == "Czechia"         ~ "CZ",
      country == "Denmark"         ~ "DK",
      country == "Estonia"         ~ "EE",
      country == "Finland"         ~ "FI",
      country == "France"          ~ "FR",
      country == "Germany"         ~ "DE",
      country == "Greece"          ~ "EL",
      country == "Hungary"         ~ "HU",
      country == "Iceland"         ~ "IS",
      country == "Ireland"         ~ "IE",
      country == "Italy"           ~ "IT",
      country == "Latvia"          ~ "LV",
      country == "Lithuania"       ~ "LT",
      country == "Luxembourg"      ~ "LU",
      country == "Malta"           ~ "MT",
      country == "Netherlands"     ~ "NL",
      country == "Norway"          ~ "NO",
      country == "Poland"          ~ "PL",
      country == "Portugal"        ~ "PT",
      country == "Romania"         ~ "RO",
      country == "Slovakia"        ~ "SK",
      country == "Slovenia"        ~ "SI",
      country == "Spain"           ~ "ES",
      country == "Sweden"          ~ "SE",
      country == "Switzerland"     ~ "CH",
      TRUE ~ NA_character_  
    ))

map_data <- europe_sf %>%
  left_join(diffs_clean, by = c("CNTR_CODE" = "CNTR_CODE"))


#MILITARY MAP 
scale_min <- -0.5
scale_max <- 2.5

my_colors <- paletteer::paletteer_c(
  "grDevices::Geyser",
  n = 100,
  direction = 1
)

# Military Spending Map
def <- ggplot(map_data) +
  geom_sf(aes(fill = defence_dif), color = NA) +  
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(scale_min, scale_max),
    na.value = "lightgrey",
    name = "Percentage Points"
  ) +
  coord_sf(xlim = c(-25, 40), ylim = c(34, 72), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  ) +
  labs(
    title = "Defence Spending Change across Europe",
    subtitle = "2013–2023"
  )
def
ggsave("defence_diff.png", def, width = 16, height = 10, dpi = 300)

# Climate Spending Map
climate <- ggplot(map_data) +
  geom_sf(aes(fill = env_dif), color = NA) +  
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(scale_min, scale_max),  
    na.value = "lightgrey",
    name = "Percentage Points"
  ) +
  coord_sf(xlim = c(-25, 40), ylim = c(34, 72), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  ) +
  labs(
    title = "Environmental Spending Change across Europe",
    subtitle = "2013–2023"
  )

climate
ggsave("defence_diff.png", climate, width = 16, height = 10, dpi = 300)

