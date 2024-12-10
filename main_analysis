
#####################################
Loading necessary libraries
#####################################

# Load necessary libraries
library(data.table)
library(dplyr)
library(lubridate)
library(tidyverse)
library(mice)
library(splines)
library(dlnm)
library(glmmTMB)
library(MuMIn)
library(MASS)
library(ggplot2)
library(sf)
library(weathermetrics)
library(readr)
library(viridis)
library(sf)
library(patchwork)


########################################
# STAGE 1: DATA PREPROCESSING
###################################

# Set working directory
setwd("/Users/fati0049/Desktop/My_data_2024/Research/Zohra1/Rcode/analysis_paper2/Rcode")

# Load data
LBW <- fread("Final_LBW_data.csv", header = TRUE) %>%
  rename(province = Region)

# Ensure the date column is in Date format
LBW$date <- as.Date(LBW$date)

# Extract Year and Month from the date column
LBW$Year <- year(LBW$date)
LBW$Month <- month(LBW$date)

# Load population data and merge with the main dataset
population <- fread("population_data_prov.csv", header = TRUE) %>%
  rename(province = Province)

LBW <- left_join(LBW, population, by = c("province", "Month", "Year"))

# Exclude specific provinces and filter years
LBW <- LBW %>%
  filter(!province %in% c("AJK")) %>%  # Exclude AJK from provinces
  filter(Year >= 2008 & Year <= 2017)   # Exclude years beyond 2008 to 2017

# LBW <- LBW %>%
#   filter(!province %in% c("AJK", "GB")) %>%  # Exclude AJK and GB from provinces
#   filter(Year >= 2008 & Year <= 2017)   # Exclude years beyond 2008 to 2017


# Assign 0 to NA values in LBW_count
LBW$LBW_count[is.na(LBW$LBW_count)] <- 0

# Convert relevant variables to factors
LBW_clean <- LBW %>%
  mutate(
    Education = as.factor(Education),
    `Place of delivery` = as.factor(`Place of delivery`),
    `Wealth Index` = as.factor(`Wealth Index`),
    province = as.factor(province),
    `Mother Age Group` = as.factor(`Mother Age Group`)
  )


# Rename columns to avoid spaces or special characters before imputation
LBW_clean <- LBW_clean %>%
  rename(
    Place_of_delivery = `Place of delivery`,
    Wealth_Index = `Wealth Index`,
    Mother_Age_Group = `Mother Age Group`
  )

# Select only the variables to impute
impute_data <- LBW_clean %>%
  dplyr::select(Education, Place_of_delivery, Wealth_Index, Mother_Age_Group)

# Perform the imputation with appropriate methods for each variable
imputation_methods <- make.method(impute_data)
imputation_methods["Education"] <- "polyreg"             # Polynomial regression for categorical data
imputation_methods["Place_of_delivery"] <- "polyreg"
imputation_methods["Wealth_Index"] <- "polyreg"
imputation_methods["Mother_Age_Group"] <- "polyreg"

# Perform multiple imputations
set.seed(123)  # Set a seed for reproducibility
imputed <- mice(impute_data, method = imputation_methods, m = 5)

# Extract the completed dataset (first imputed dataset)
LBW_clean_imputed <- complete(imputed, 1)

# Replace the original columns with the imputed ones, renaming back if needed
LBW_clean <- LBW_clean %>%
  mutate(
    Education = LBW_clean_imputed$Education,
    `Place of delivery` = LBW_clean_imputed$Place_of_delivery,
    `Wealth Index` = LBW_clean_imputed$Wealth_Index,
    `Mother Age Group` = LBW_clean_imputed$Mother_Age_Group
  )

# Check for any remaining missing values
summary(LBW_clean)

LBW_clean <- LBW_clean %>%
  arrange(province, date)



library(data.table)
library(dplyr)

# Step 1: Load province-level precipitation data
province_precipitation <- fread("/Users/fati0049/Desktop/My_data_2024/Research/Processed_env_data/Rcode/province_precipitation.csv", header = TRUE)

# Ensure the date column is in Date format
province_precipitation$date <- as.Date(province_precipitation$date)

# Check if province column is a factor (convert if necessary)
province_precipitation$ADM1_EN <- as.factor(province_precipitation$ADM1_EN)

# Step 2: Harmonize column names for merging
province_precipitation <- province_precipitation %>%
  rename(province = ADM1_EN)

# Replace province names with desired names
province_precipitation <- province_precipitation %>%
  mutate(
    province = case_when(
      province == "Azad Kashmir" ~ "AJK",           # Replace Azad Kashmir with AJK
      province == "Balochistan" ~ "Baluchistan",    # Ensure consistent spelling
      province == "Gilgit Baltistan" ~ "GB",        # Replace Gilgit Baltistan with GB
      province == "Khyber Pakhtunkhwa" ~ "KPK",     # Replace Khyber Pakhtunkhwa with KPK
      TRUE ~ province                                # Keep other names as is
    )
  )

# Verify the changes
unique_provinces <- unique(province_precipitation$province)
print(unique_provinces)


LBW_clean <- LBW_clean %>%
  left_join(province_precipitation, by = c("province", "date"))


# Create scaled variables directly within LBW_clean
# Create scaled variables directly within LBW_clean and store medians
LBW_clean <- LBW_clean %>%
  mutate(
    scaled_tmean = scale(tmean_C, center = TRUE, scale = TRUE),                  # Scale mean temperature
    scaled_humidity = scale(humidity, center = TRUE, scale = TRUE),             # Scale humidity
    scaled_precipitation = scale(total_precipitation, center = TRUE, scale = TRUE), # Scale total precipitation
    PM25_scaled = scale(PM25, center = TRUE, scale = TRUE)                      # Scale PM2.5
  )

# Calculate medians of scaled variables
scaled_medians <- LBW_clean %>%
  summarize(
    median_scaled_tmean = median(scaled_tmean, na.rm = TRUE),
    median_scaled_humidity = median(scaled_humidity, na.rm = TRUE),
    median_scaled_precipitation = median(scaled_precipitation, na.rm = TRUE),
    median_PM25_scaled = median(PM25_scaled, na.rm = TRUE)
  )

print(scaled_medians)



# Save the LBW_clean dataset as a CSV file
write.csv(LBW_clean, "LBW_clean.csv", row.names = FALSE)

# Print a confirmation message
cat("The LBW_clean dataset has been successfully saved as 'LBW_clean.csv'.\n")




###############################################
#STAGE 2: MODELING PARAMETERS 
############################################


# Lag and polynomial degrees
lag <- 7
degree_exposure <- 2
degree_lag <- 1

# Define function to extract coefficients and variance-covariance matrices
extract_cb_parameters <- function(model, cb_name) {
  coef_full <- fixef(model)$cond
  cb_coef_names <- grep(paste0("^", cb_name), names(coef_full), value = TRUE)
  cb_coefs <- coef_full[cb_coef_names]
  vcov_full <- vcov(model)$cond
  cb_vcov <- vcov_full[cb_coef_names, cb_coef_names]
  list(coefs = cb_coefs, vcov = cb_vcov)
}



# Define function to predict RR values using raw median for centering
predict_rr <- function(cb, coefs, vcov, data, cb_name, scaled_median) {
  # Use median of scaled temperature for centering
  tmean_percentiles <- quantile(data$scaled_tmean, probs = c(0.01, 0.10, 0.90, 0.99), na.rm = TRUE)
  
  # Calculate the median of the scaled temperature
  scaled_median <- median(LBW_clean$scaled_tmean, na.rm = TRUE)
  
  pred <- crosspred(
    cb,
    coef = coefs,
    vcov = vcov,
    at = tmean_percentiles,  # Use percentiles of scaled temperature
    cen = scaled_median      # Center on median of scaled values
  )
  
  rr_results <- data.frame(
    Percentile = names(tmean_percentiles),
    Scaled_Temperature = tmean_percentiles,
    RR = exp(pred$allfit),
    LCI = exp(pred$alllow),
    UCI = exp(pred$allhigh)
  )
  
  cat("\nRelative Risk Predictions for", cb_name, ":\n")
  print(rr_results)
  return(rr_results)
}


# Null Model
null_model <- glmmTMB(
  LBW_count ~ poly(PM25_scaled, 2) + factor(Education) +
    offset(log(WRA_pop)) + (1 | province),
  data = LBW_clean,
  family = nbinom2()
)
null_model_aic <- AIC(null_model)
cat("\nNull Model AIC:", null_model_aic, "\n")



# Linear Model
cb_linear <- crossbasis(
  LBW_clean$scaled_tmean,  # Use scaled temperature
  lag = lag,
  argvar = list(fun = "lin"),
  arglag = list(fun = "poly", degree = degree_lag)
)
linear_model <- glmmTMB(
  LBW_count ~ cb_linear + poly(PM25_scaled, 2) + factor(Education) +
    offset(log(WRA_pop)) + (1 | province),
  data = LBW_clean,
  family = nbinom2()
)
linear_params <- extract_cb_parameters(linear_model, "cb_linear")
linear_rr <- predict_rr(cb_linear, linear_params$coefs, linear_params$vcov, LBW_clean, "Linear Model")

# Polynomial Models
cb_polynomial <- crossbasis(
  LBW_clean$scaled_tmean,  # Use scaled temperature
  lag = lag,
  argvar = list(fun = "poly", degree = degree_exposure),
  arglag = list(fun = "poly", degree = degree_lag)
)

# Model 2: Scaled PM25 & Education
model_2 <- glmmTMB(
  LBW_count ~ cb_polynomial + poly(PM25_scaled, 2) + factor(Education) +
    offset(log(WRA_pop)) + (1 | province),
  data = LBW_clean,
  family = nbinom2()
)
model_2_params <- extract_cb_parameters(model_2, "cb_polynomial")
model_2_rr <- predict_rr(cb_polynomial, model_2_params$coefs, model_2_params$vcov, LBW_clean, "Model 2")

# Model 3: Scaled PM25, Scaled Humidity & Education
model_3 <- glmmTMB(
  LBW_count ~ cb_polynomial + poly(PM25_scaled, 2) + poly(scaled_humidity, 2) + factor(Education) +
    offset(log(WRA_pop)) + (1 | province),
  data = LBW_clean,
  family = nbinom2()
)
model_3_params <- extract_cb_parameters(model_3, "cb_polynomial")
model_3_rr <- predict_rr(cb_polynomial, model_3_params$coefs, model_3_params$vcov, LBW_clean, "Model 3")

# Model 4: Scaled PM25, Scaled Precipitation & Education
model_4 <- glmmTMB(
  LBW_count ~ cb_polynomial + poly(PM25_scaled, 2) + poly(scaled_precipitation, 2) + factor(Education) +
    offset(log(WRA_pop)) + (1 | province),
  data = LBW_clean,
  family = nbinom2()
)
model_4_params <- extract_cb_parameters(model_4, "cb_polynomial")
model_4_rr <- predict_rr(cb_polynomial, model_4_params$coefs, model_4_params$vcov, LBW_clean, "Model 4")

# Model 5: Scaled PM25, Scaled Humidity, Scaled Precipitation & Education
model_5 <- glmmTMB(
  LBW_count ~ cb_polynomial + poly(PM25_scaled, 2) + poly(scaled_humidity, 2) +
    poly(scaled_precipitation, 2) + factor(Education) +
    offset(log(WRA_pop)) + (1 | province),
  data = LBW_clean,
  family = nbinom2()
)
model_5_params <- extract_cb_parameters(model_5, "cb_polynomial")
model_5_rr <- predict_rr(cb_polynomial, model_5_params$coefs, model_5_params$vcov, LBW_clean, "Model 5")

# Summarize AIC for all models
models <- list(
  Null_Model = list(model = null_model, AIC = null_model_aic),
  Linear_Model = list(model = linear_model, AIC = AIC(linear_model)),
  Model_2 = list(model = model_2, AIC = AIC(model_2)),
  Model_3 = list(model = model_3, AIC = AIC(model_3)),
  Model_4 = list(model = model_4, AIC = AIC(model_4)),
  Model_5 = list(model = model_5, AIC = AIC(model_5))
)

# Extract AIC values
aic_values <- sapply(models, function(x) x$AIC)
min_aic <- min(aic_values)
delta_aic <- aic_values - min_aic
weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))

# Model selection summary
model_summary <- data.frame(
  Model = names(models),
  AIC = aic_values,
  Delta_AIC = delta_aic,
  Weight = weights
)

# Display the summary
print(model_summary)
summary(model_2)

cat("\nVariance-Covariance Matrix:\n")
print(vcov(model_2)$cond)

cat("\nRandom Effects for Province 'province_name':\n")
print(coef(model_2)$cond)


############################
#2.1 Model visualisation
########################
generate_plot_data <- function(cb, coefs, vcov, data, model_name, color) {
  # Extract scaling parameters
  mean_tmean <- attr(data$scaled_tmean, "scaled:center")
  sd_tmean <- attr(data$scaled_tmean, "scaled:scale")
  
  # Create a sequence of original temperatures
  temp_seq <- seq(min(data$tmean_C, na.rm = TRUE), max(data$tmean_C, na.rm = TRUE), length.out = 100)
  
  # Transform to the scaled space
  temp_seq_scaled <- (temp_seq - mean_tmean) / sd_tmean
  
  # Calculate the median of scaled temperature for centering
  temp_median_scaled <- median(data$scaled_tmean, na.rm = TRUE)
  
  # Perform prediction using the crossbasis
  pred <- crosspred(
    cb,
    coef = coefs,
    vcov = vcov,
    at = temp_seq_scaled,  # Use scaled temperature sequence
    cen = temp_median_scaled  # Center on scaled median temperature
  )
  
  # Create a data frame for plotting with original temperatures
  data.frame(
    Temperature = temp_seq,  # Original temperature
    RR = exp(pred$allfit),   # Rescaled RR
    LCI = exp(pred$alllow),  # Rescaled lower confidence interval
    UCI = exp(pred$allhigh), # Rescaled upper confidence interval
    Model = model_name,
    Color = color
  )
}

plot_data <- bind_rows(
  generate_plot_data(cb_linear, linear_params$coefs, linear_params$vcov, LBW_clean, "Linear Model", "blue"),
  generate_plot_data(cb_polynomial, model_2_params$coefs, model_2_params$vcov, LBW_clean, "Model 2", "red"),
  generate_plot_data(cb_polynomial, model_3_params$coefs, model_3_params$vcov, LBW_clean, "Model 3", "green"),
  generate_plot_data(cb_polynomial, model_4_params$coefs, model_4_params$vcov, LBW_clean, "Model 4", "purple"),
  generate_plot_data(cb_polynomial, model_5_params$coefs, model_5_params$vcov, LBW_clean, "Model 5", "orange")
)

combined_plot <- ggplot(plot_data, aes(x = Temperature, y = RR, color = Model, fill = Model)) +
  # Main lines for each model
  geom_line(size = 1.2) +
  # Shaded confidence intervals for each model
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2, color = NA) +
  # Reference line at RR = 1
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  # Labels and aesthetics
  labs(
    title = "Combined Exposure-Response Associations for All Models (Original Scale)",
    x = "Temperature (°C)",
    y = "Relative Risk (RR)"
  ) +
  scale_color_manual(values = unique(plot_data$Color)) +
  scale_fill_manual(values = unique(plot_data$Color)) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Save the plot as a PNG file
ggsave("Combined_Exposure_Response_Original_Scale.png", combined_plot, width = 10, height = 8)

# Display the plot
print(combined_plot)

##########################
#STAGE 3: RUNNING ALL MODELING AND ESTIMATING WEIGHTAGES
########################




# Extract coefficients and variance-covariance matrices for all models
model_params <- lapply(models, function(x) {
  extract_cb_parameters(x$model, "cb_polynomial")  # Extract coefficients and vcov
})

# Align the names of coefficients and vcov matrices to ensure compatibility
all_coef_names <- unique(unlist(lapply(model_params, function(x) names(x$coefs))))
aligned_params <- lapply(model_params, function(params) {
  # Align coefficients
  aligned_coefs <- setNames(rep(0, length(all_coef_names)), all_coef_names)
  aligned_coefs[names(params$coefs)] <- params$coefs
  
  # Align vcov matrix
  aligned_vcov <- matrix(0, nrow = length(all_coef_names), ncol = length(all_coef_names),
                         dimnames = list(all_coef_names, all_coef_names))
  aligned_vcov[rownames(params$vcov), colnames(params$vcov)] <- params$vcov
  
  list(coefs = aligned_coefs, vcov = aligned_vcov)
})

# Compute weighted coefficients and variance-covariance matrices
weighted_coefs <- rep(0, length(all_coef_names))
names(weighted_coefs) <- all_coef_names
weighted_vcov <- matrix(0, nrow = length(all_coef_names), ncol = length(all_coef_names),
                        dimnames = list(all_coef_names, all_coef_names))

for (i in seq_along(aligned_params)) {
  params <- aligned_params[[i]]
  weight <- weights[i]
  
  # Add weighted coefficients and vcov
  weighted_coefs <- weighted_coefs + weight * params$coefs
  weighted_vcov <- weighted_vcov + weight * params$vcov
}

# Display results
cat("\nModel Selection Summary:\n")
print(model_summary)

cat("\nWeighted Coefficients:\n")
print(weighted_coefs)

cat("\nWeighted Variance-Covariance Matrix:\n")
print(weighted_vcov)


# Province-Level Predictions with Scaled Variables

# Create results storage
province_full_rr_results <- list()   # For full array
province_percentile_rr_results <- list()  # For specific percentiles
provinces <- unique(LBW_clean$province)   # Get unique provinces

# Extract scaling parameters
mean_tmean <- attr(LBW_clean$scaled_tmean, "scaled:center")
sd_tmean <- attr(LBW_clean$scaled_tmean, "scaled:scale")

for (prov in provinces) {
  cat("\n--- Processing Province:", prov, "---\n")
  
  # Subset data for the current province
  prov_data <- subset(LBW_clean, province == prov)
  
  # Define province-specific median and percentiles
  median_tmean_prov <- median(prov_data$scaled_tmean, na.rm = TRUE)  # Scaled median for centering
  percentiles <- c(0.01, 0.10, 0.90, 0.99)                           # Percentile definitions
  tmean_percentiles_scaled <- quantile(prov_data$scaled_tmean, probs = percentiles, na.rm = TRUE)
  
  # Define a full temperature sequence in scaled space
  temp_seq_scaled <- seq(min(prov_data$scaled_tmean, na.rm = TRUE), 
                         max(prov_data$scaled_tmean, na.rm = TRUE), 
                         length.out = 100)
  
  # Rescale the sequence back to the original temperature scale for output
  temp_seq_original <- temp_seq_scaled * sd_tmean + mean_tmean
  
  # Use weighted coefficients and variance-covariance matrix for full array prediction
  pred_full <- crosspred(
    cb_polynomial,  # Use the crossbasis from the polynomial model
    coef = weighted_coefs,  # Weighted coefficients
    vcov = weighted_vcov,   # Weighted variance-covariance matrix
    at = temp_seq_scaled,   # Predict for the scaled temperature sequence
    cen = median_tmean_prov # Centering at the province-specific scaled median
  )
  
  # Create a data frame with full array results
  rr_full <- data.frame(
    Province = prov,
    Temperature = temp_seq_original,              # Rescaled to original temperature scale
    Percentile = ecdf(prov_data$tmean_C)(temp_seq_original),  # Percentiles for full array
    RR = exp(pred_full$allfit),
    LCI = exp(pred_full$alllow),
    UCI = exp(pred_full$allhigh)
  )
  province_full_rr_results[[prov]] <- rr_full
  
  # Use weighted coefficients for specific percentiles prediction
  pred_percentiles <- crosspred(
    cb_polynomial, 
    coef = weighted_coefs, 
    vcov = weighted_vcov, 
    at = tmean_percentiles_scaled,  # Scaled temperature percentiles
    cen = median_tmean_prov         # Centering at the province-specific scaled median
  )
  
  # Rescale percentiles back to original scale
  tmean_percentiles_original <- tmean_percentiles_scaled * sd_tmean + mean_tmean
  
  # Create a data frame with specific percentile results
  rr_at_percentiles <- data.frame(
    Province = prov,
    Percentile = names(tmean_percentiles_scaled),
    tmean_C = tmean_percentiles_original,  # Rescaled to original temperature
    RR = exp(pred_percentiles$allfit),
    LCI = exp(pred_percentiles$alllow),
    UCI = exp(pred_percentiles$allhigh)
  )
  province_percentile_rr_results[[prov]] <- rr_at_percentiles
  
  # Print previews
  cat("\nFull Array Preview:\n")
  print(head(rr_full))
  
  cat("\nSpecific Percentiles Preview:\n")
  print(rr_at_percentiles)
}

# Combine all results into data frames
all_full_rr_results <- do.call(rbind, province_full_rr_results)
all_percentile_rr_results <- do.call(rbind, province_percentile_rr_results)

# Save the results as CSV files
write.csv(all_full_rr_results, "Province_Specific_Temperature_RR_Full_Array_Scaled.csv", row.names = FALSE)
write.csv(all_percentile_rr_results, "Province_Specific_Temperature_RR_Specific_Percentiles_Scaled.csv", row.names = FALSE)

cat("Results saved as 'Province_Specific_Temperature_RR_Full_Array_Scaled.csv' and 'Province_Specific_Temperature_RR_Specific_Percentiles_Scaled.csv'.\n")

# Extract the 99th percentile RR for each province
rr_99th_percentile <- all_percentile_rr_results %>%
  filter(Percentile == "99%") %>%  # Select only the 99th percentile
  dplyr::select(Province, RR, LCI, UCI)  # Keep relevant columns

# Preview the extracted data
cat("\n99th Percentile Relative Risks (RR):\n")
print(rr_99th_percentile)

# Save the 99th percentile RR as a CSV file
write.csv(rr_99th_percentile, "Province_RR.csv", row.names = FALSE)

cat("\nRR at the 99th percentile saved as 'Province_RR_99th_Percentile_Scaled.csv'.\n")


#####################
#3.1: Data visualisation
#####################


# Ensure proper column naming
all_full_rr_results <- all_full_rr_results %>%
  rename(province = Province)

# Rescale temperatures back to the original scale and align RR ~ 1 with minimal CI range
all_full_rr_results <- all_full_rr_results %>%
  group_by(province) %>%
  mutate(
    AlignTemp = Temperature[which.min(abs(RR - 1) + (UCI - LCI))],  # Closest RR to 1 with narrowest CI
    P1 = quantile(Temperature, 0.01, na.rm = TRUE),  # 1st percentile
    P10 = quantile(Temperature, 0.10, na.rm = TRUE),  # 10th percentile
    P90 = quantile(Temperature, 0.90, na.rm = TRUE),  # 90th percentile
    P99 = quantile(Temperature, 0.99, na.rm = TRUE)   # 99th percentile
  ) %>%
  ungroup()

# Create a dataset for the temperature frequency histogram
hist_data <- LBW_clean %>%
  group_by(province) %>%
  mutate(
    tmean_C_bin = cut(tmean_C, breaks = seq(min(tmean_C, na.rm = TRUE), max(tmean_C, na.rm = TRUE), by = 1))
  ) %>%
  count(province, tmean_C_bin) %>%
  ungroup() %>%
  mutate(
    tmean_C_mid = as.numeric(sub("\\((.+),.*", "\\1", as.character(tmean_C_bin))) + 0.5  # Extract bin midpoints
  )

# Initialize an empty list to store combined plots for each province
combined_plots <- list()

# Create plots for each province
for (prov in unique(all_full_rr_results$province)) {
  # Subset data for the current province
  rr_data <- all_full_rr_results %>% filter(province == prov)
  hist_data_subset <- hist_data %>% filter(province == prov)
  
  # Calculate breaks for histogram y-axis
  max_count <- max(hist_data_subset$n, na.rm = TRUE)
  mid_count <- round(max_count / 2)
  y_breaks <- c(0, mid_count, max_count)  # Min, mid, and max values
  
  # Risk curve plot
  risk_plot <- ggplot(rr_data, aes(x = Temperature, y = RR)) +
    geom_line(data = subset(rr_data, Temperature <= AlignTemp), 
              aes(color = "Cold"), size = 1) +
    geom_line(data = subset(rr_data, Temperature > AlignTemp), 
              aes(color = "Hot"), size = 1) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "grey", alpha = 0.2) +
    geom_vline(aes(xintercept = AlignTemp), linetype = "solid", size = 0.8, color = "black") +
    geom_vline(aes(xintercept = P1), linetype = "dotted", size = 0.5, color = "black") +
    geom_vline(aes(xintercept = P10), linetype = "dotted", size = 0.5, color = "black") +
    geom_vline(aes(xintercept = P90), linetype = "dotted", size = 0.5, color = "black") +
    geom_vline(aes(xintercept = P99), linetype = "dotted", size = 0.5, color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.5, color = "black") +
    labs(
      title = prov,
      x = NULL,
      y = "RR"
    ) +
    scale_color_manual(values = c("Cold" = "blue", "Hot" = "red")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  # Histogram plot
  hist_plot <- ggplot(hist_data_subset, aes(x = tmean_C_mid, y = n)) +
    geom_bar(stat = "identity", fill = "lightblue", alpha = 0.6, color = "lightblue") +
    scale_y_continuous(
      breaks = y_breaks,  # Min, mid, and max
      labels = scales::comma  # Format numbers
    ) +
    labs(
      x = "Temperature (°C)",
      y = "Months"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(color = "black")
    )
  
  # Combine the two plots using `patchwork`
  combined_plot <- risk_plot / hist_plot + plot_layout(heights = c(3, 1))
  combined_plots[[prov]] <- combined_plot
}

# Combine all province-specific plots into a grid
final_combined_plot <- wrap_plots(combined_plots, ncol = 2)

# Save the final plot
ggsave("Figure1_Final_Province_ER_Curves_Scaled.png", final_combined_plot, width = 16, height = 12)

# Display the final plot
print(final_combined_plot)

########################
#STAGE 4: ATTRIBUTABLE FRACTION
#######################


# Set a random seed for reproducibility
set.seed(123)


# Load datasets
#LBW_clean <- read.csv("LBW_clean.csv")  # Baseline dataset
prov_RCP4.5_48 <- read.csv("prov_RCP4.5_48.csv")  # RCP 4.5 (2048-2057)
prov_RCP8.5_48 <- read.csv("prov_RCP8.5_48.csv")  # RCP 8.5 (2048-2057)
prov_RCP4.5_68 <- read.csv("prov_RCP4.5_68.csv")  # RCP 4.5 (2068-2077)
prov_RCP8.5_68 <- read.csv("prov_RCP8.5_68.csv")  # RCP 8.5 (2068-2077)

# Ensure date column is in Date format for all datasets
LBW_clean$date <- as.Date(LBW_clean$date)
prov_RCP4.5_48$date.y <- as.Date(prov_RCP4.5_48$date.y, format = "%d/%m/%Y")
prov_RCP8.5_48$date.y <- as.Date(prov_RCP8.5_48$date.y, format = "%d/%m/%Y")
prov_RCP4.5_68$date.y <- as.Date(prov_RCP4.5_68$date.y, format = "%d/%m/%Y")
prov_RCP8.5_68$date.y <- as.Date(prov_RCP8.5_68$date.y, format = "%d/%m/%Y")


# Function to calculate SD for log-transformed RR
calculate_sd <- function(uci, lci) {
  (log(uci) - log(lci)) / (2 * 1.96)  # For 95% CI
}

# Monte Carlo simulation to estimate AF
simulate_af <- function(rr_values, prop_exposure, n_simulations, sd_rr) {
  sim_af <- numeric(n_simulations)
  for (i in 1:n_simulations) {
    sim_rr <- exp(rnorm(length(rr_values), mean = log(rr_values), sd = sd_rr))
    sim_af[i] <- mean(((sim_rr - 1) / sim_rr) * prop_exposure * 100, na.rm = TRUE)
  }
  list(
    af = mean(sim_af, na.rm = TRUE),
    lci = quantile(sim_af, probs = 0.025, na.rm = TRUE),
    uci = quantile(sim_af, probs = 0.975, na.rm = TRUE)
  )
}

# Function to scale projected temperatures using baseline parameters
scale_projected_temperatures <- function(data, mean_tmean, sd_tmean) {
  data <- data %>%
    mutate(scaled_tmean = (tmean_C - mean_tmean) / sd_tmean)
  return(data)
}


# Calculate AlignTemp using baseline RR results
align_temp_values <- all_full_rr_results %>%
  group_by(province) %>%
  summarize(
    AlignTemp = Temperature[which.min(abs(RR - 1) + (UCI - LCI))]  # Closest RR to 1
  )



# Add AlignTemp to baseline
LBW_clean <- LBW_clean %>%
  left_join(align_temp_values, by = "province")

# Add AlignTemp to projected datasets
datasets <- list(
  prov_RCP4.5_48 = prov_RCP4.5_48,
  prov_RCP8.5_48 = prov_RCP8.5_48,
  prov_RCP4.5_68 = prov_RCP4.5_68,
  prov_RCP8.5_68 = prov_RCP8.5_68
)

for (name in names(datasets)) {
  datasets[[name]] <- datasets[[name]] %>%
    rename(province = ADM1_EN) %>%
    left_join(align_temp_values, by = "province")
}

list2env(datasets, envir = .GlobalEnv)  # Assign back to global environment



# Extract scaling parameters from baseline
mean_tmean <- attr(LBW_clean$scaled_tmean, "scaled:center")
sd_tmean <- attr(LBW_clean$scaled_tmean, "scaled:scale")

# Scale temperatures in projections
prov_RCP4.5_48 <- scale_projected_temperatures(prov_RCP4.5_48, mean_tmean, sd_tmean)
prov_RCP8.5_48 <- scale_projected_temperatures(prov_RCP8.5_48, mean_tmean, sd_tmean)
prov_RCP4.5_68 <- scale_projected_temperatures(prov_RCP4.5_68, mean_tmean, sd_tmean)
prov_RCP8.5_68 <- scale_projected_temperatures(prov_RCP8.5_68, mean_tmean, sd_tmean)



# AF calculation for baseline
heat_af_results <- list()
cold_af_results <- list()

for (prov in unique(LBW_clean$province)) {
  prov_data <- subset(LBW_clean, province == prov)
  prov_rr <- subset(all_full_rr_results, province == prov)
  align_temp <- first(prov_rr$AlignTemp)
  
  # Proportion of exposure
  prop_heat <- mean(prov_data$tmean_C > align_temp, na.rm = TRUE)
  prop_cold <- mean(prov_data$tmean_C < align_temp, na.rm = TRUE)
  
  # Extract RR and SD
  rr_heat <- prov_rr$RR[prov_rr$Temperature > align_temp]
  sd_heat <- calculate_sd(prov_rr$UCI[prov_rr$Temperature > align_temp],
                          prov_rr$LCI[prov_rr$Temperature > align_temp])
  rr_cold <- prov_rr$RR[prov_rr$Temperature < align_temp]
  sd_cold <- calculate_sd(prov_rr$UCI[prov_rr$Temperature < align_temp],
                          prov_rr$LCI[prov_rr$Temperature < align_temp])
  
  # Simulate AF
  heat_results <- simulate_af(rr_heat, prop_heat, 1000, sd_heat)
  cold_results <- simulate_af(rr_cold, prop_cold, 1000, sd_cold)
  
  # Store results
  heat_af_results[[prov]] <- data.frame(
    Province = prov, AF_Heat = heat_results$af,
    LCI_Heat = heat_results$lci, UCI_Heat = heat_results$uci
  )
  cold_af_results[[prov]] <- data.frame(
    Province = prov, AF_Cold = cold_results$af,
    LCI_Cold = cold_results$lci, UCI_Cold = cold_results$uci
  )
}

# Combine AF results
all_heat_af_results <- do.call(rbind, heat_af_results)
all_cold_af_results <- do.call(rbind, cold_af_results)
baseline_af_results <- merge(all_heat_af_results, all_cold_af_results, by = "Province")

# Save baseline AF results
#write.csv(baseline_af_results, "Baseline_AF_Results.csv", row.names = FALSE)



# AF calculation for projections
projected_scenarios <- list(
  RCP4.5_48 = prov_RCP4.5_48,
  RCP8.5_48 = prov_RCP8.5_48,
  RCP4.5_68 = prov_RCP4.5_68,
  RCP8.5_68 = prov_RCP8.5_68
)

projected_af_results <- list()

for (scenario_name in names(projected_scenarios)) {
  scenario_data <- projected_scenarios[[scenario_name]]
  
  scenario_af_results <- list()
  
  for (prov in unique(scenario_data$province)) {
    prov_data <- subset(scenario_data, province == prov)
    prov_rr <- subset(all_full_rr_results, province == prov)
    align_temp <- first(prov_rr$AlignTemp)
    
    # Proportion of exposure
    prop_heat <- mean(prov_data$tmean_C > align_temp, na.rm = TRUE)
    prop_cold <- mean(prov_data$tmean_C < align_temp, na.rm = TRUE)
    
    # Extract RR and SD
    rr_heat <- prov_rr$RR[prov_rr$Temperature > align_temp]
    sd_heat <- calculate_sd(prov_rr$UCI[prov_rr$Temperature > align_temp],
                            prov_rr$LCI[prov_rr$Temperature > align_temp])
    rr_cold <- prov_rr$RR[prov_rr$Temperature < align_temp]
    sd_cold <- calculate_sd(prov_rr$UCI[prov_rr$Temperature < align_temp],
                            prov_rr$LCI[prov_rr$Temperature < align_temp])
    
    # Simulate AF
    heat_results <- simulate_af(rr_heat, prop_heat, 1000, sd_heat)
    cold_results <- simulate_af(rr_cold, prop_cold, 1000, sd_cold)
    
    # Store results
    scenario_af_results[[prov]] <- data.frame(
      Province = prov, Scenario = scenario_name,
      AF_Heat = heat_results$af, LCI_Heat = heat_results$lci, UCI_Heat = heat_results$uci,
      AF_Cold = cold_results$af, LCI_Cold = cold_results$lci, UCI_Cold = cold_results$uci
    )
  }
  
  projected_af_results[[scenario_name]] <- do.call(rbind, scenario_af_results)
}

# Combine all scenario results
all_projected_af_results <- do.call(rbind, projected_af_results)

# Save projected AF results
#write.csv(all_projected_af_results, "Projected_AF_Results.csv", row.names = FALSE)



# Add Period and Scenario details to the baseline dataset
baseline_af_results <- baseline_af_results %>%
  mutate(
    Period = "Baseline",
    Scenario = "NA"
  )

# Add Period and Scenario details to the projected dataset
all_projected_af_results <- all_projected_af_results %>%
  mutate(
    Period = ifelse(grepl("48", Scenario), "2048–2057", "2068–2077"),
    Scenario = gsub("_.*", "", Scenario)  # Clean scenario names (e.g., RCP4.5_48 to RCP4.5)
  )

# Combine baseline and projected results
all_combined_af_results <- bind_rows(
  baseline_af_results %>% dplyr::select(Province, AF_Heat, LCI_Heat, UCI_Heat, AF_Cold, LCI_Cold, UCI_Cold, Period, Scenario),
  all_projected_af_results %>% dplyr::select(Province, AF_Heat, LCI_Heat, UCI_Heat, AF_Cold, LCI_Cold, UCI_Cold, Period, Scenario)
)

# Ensure consistent formatting
all_combined_af_results <- all_combined_af_results %>%
  mutate(
    Period = factor(Period, levels = c("Baseline", "2048–2057", "2068–2077")),
    Scenario = factor(Scenario, levels = c("NA", "RCP4.5", "RCP8.5"))
  )

# Save the combined dataset
write.csv(all_combined_af_results, "Combined_AF_Results.csv", row.names = FALSE)

# Print the combined dataset to verify
cat("Combined AF results (Baseline and Projections):\n")
print(head(all_combined_af_results))


# Remove AJK from the combined dataset and create a heat-only dataset
all_af_combined_heat <- all_combined_af_results %>%
  filter(Province != "AJK") %>%  # Exclude AJK
  dplyr::select(Province, Period, Scenario, AF_Heat, LCI_Heat, UCI_Heat)  # Keep heat-related columns

# Save the heat-only dataset
write.csv(all_af_combined_heat, "Final_AF_Combined_Heat_Only.csv", row.names = FALSE)

# Print the first few rows to verify
cat("Final Heat-Only Dataset (Excluding AJK):\n")
print(head(all_af_combined_heat))


#######################
# 4.1 Data visualisation
####################



# Ensure consistent column formats in `all_af_combined`
all_af_combined_heat <- all_af_combined_heat %>%
  mutate(
    AF_Heat = as.numeric(AF_Heat),       # Ensure numeric
    LCI_Heat = as.numeric(LCI_Heat),    # Ensure numeric
    UCI_Heat = as.numeric(UCI_Heat),    # Ensure numeric
    Period = factor(Period, levels = c("Baseline", "2048–2057", "2068–2077")),
    Scenario = ifelse(Scenario == "NA", "NA", gsub("_.*", "", Scenario)),  # Remove extra characters
    Scenario = factor(Scenario, levels = c("NA", "RCP4.5", "RCP8.5"))
  )

# Rescale temperature back to the original scale
rescale_temperature <- function(scaled_tmean, mean_tmean, sd_tmean) {
  return(scaled_tmean * sd_tmean + mean_tmean)  # Rescale to original temperature
}

# Apply rescaling to the projected datasets for visualization
prov_RCP4.5_48 <- prov_RCP4.5_48 %>%
  mutate(tmean_rescaled = rescale_temperature(scaled_tmean, mean_tmean, sd_tmean))

prov_RCP8.5_48 <- prov_RCP8.5_48 %>%
  mutate(tmean_rescaled = rescale_temperature(scaled_tmean, mean_tmean, sd_tmean))

prov_RCP4.5_68 <- prov_RCP4.5_68 %>%
  mutate(tmean_rescaled = rescale_temperature(scaled_tmean, mean_tmean, sd_tmean))

prov_RCP8.5_68 <- prov_RCP8.5_68 %>%
  mutate(tmean_rescaled = rescale_temperature(scaled_tmean, mean_tmean, sd_tmean))

generate_af_plot <- function(province_name) {
  # Filter data for the given province
  AF_data_filtered <- all_af_combined_heat %>%
    filter(Province == province_name)
  
  # Create a custom factor for Period and Scenario to reorder them
  AF_data_filtered$Period_Scenario <- factor(
    interaction(AF_data_filtered$Period, AF_data_filtered$Scenario),
    levels = c("2068–2077.RCP8.5", "2068–2077.RCP4.5",
               "2048–2057.RCP8.5", "2048–2057.RCP4.5", 
               "Baseline.NA"),
    labels = c("2068–2077 RCP8.5", "2068–2077 RCP4.5",
               "2048–2057 RCP8.5", "2048–2057 RCP4.5", 
               "Baseline")
  )
  
  # Create a text column to combine AF and 95% CI values
  AF_data_filtered$AF_CI_label <- paste0(round(AF_data_filtered$AF_Heat, 2), " (", 
                                         round(AF_data_filtered$LCI_Heat, 2), ", ", 
                                         round(AF_data_filtered$UCI_Heat, 2), ")")
  
  # Create the plot for the given province
  af_plot <- ggplot(AF_data_filtered, 
                    aes(x = AF_Heat, y = Period_Scenario, color = Scenario)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = LCI_Heat, xmax = UCI_Heat), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_text(aes(label = AF_CI_label, x = 19),  # Adjust label position to fit the range
              hjust = 0, size = 3, color = "black") +
    scale_color_manual(values = c("black", "blue", "red"), 
                       labels = c("Baseline", "RCP4.5", "RCP8.5")) +
    labs(x = "Heat-related AF (%)", y = province_name, 
         title = paste("Heat-related Attributable Fraction (AF) with 95% CI for", province_name)) +
    xlim(5, 20) +  # Manually set the x-axis range
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.grid.major.x = element_blank())
  
  return(af_plot)
}

# Generate the plots for provinces
af_plot_baluchistan <- generate_af_plot("Baluchistan")
af_plot_kpk <- generate_af_plot("KPK")
af_plot_punjab <- generate_af_plot("Punjab")
af_plot_sindh <- generate_af_plot("Sindh")

# Arrange the plots in a grid
combined_plot <- grid.arrange(af_plot_baluchistan, af_plot_kpk, af_plot_punjab, af_plot_sindh, 
                              ncol = 1, heights = c(1, 1, 1, 1))

# Save the final combined plot
ggsave(filename = "Figure2_proj_AF_plot_updated.png", plot = combined_plot, width = 10, height = 12, dpi = 300)







################################
#STAGE 5: SNESIVITY ANALYSIS TO COMPARE MEAN ADN MEDIAN TEMPERATURE
###########################

# Load necessary library
library(dplyr)

# Ensure the data is loaded
LBW_clean <- read.csv("LBW_clean.csv")

# Calculate mean and median temperature for each province
province_temp_stats <- LBW_clean %>%
  group_by(province) %>%
  summarize(
    Mean_Temperature = mean(tmean_C, na.rm = TRUE),  # Calculate mean temperature
    Median_Temperature = median(tmean_C, na.rm = TRUE),  # Calculate median temperature
    Count = n()  # Count number of observations per province
  ) %>%
  arrange(province)  # Sort results by province

# Display the results
print(province_temp_stats)



predict_rr_mean <- function(cb, coefs, vcov, data, cb_name) {
  # Calculate the mean of scaled temperature for centering
  scaled_mean <- mean(data$scaled_tmean, na.rm = TRUE)
  
  # Generate predictions
  pred <- crosspred(
    cb,
    coef = coefs,
    vcov = vcov,
    at = seq(min(data$scaled_tmean, na.rm = TRUE), max(data$scaled_tmean, na.rm = TRUE), length.out = 100),
    cen = scaled_mean  # Center on mean of scaled temperature
  )
  
  # Create a data frame with results
  rr_results <- data.frame(
    Temperature = seq(min(data$tmean_C, na.rm = TRUE), max(data$tmean_C, na.rm = TRUE), length.out = 100),
    RR = exp(pred$allfit),
    LCI = exp(pred$alllow),
    UCI = exp(pred$allhigh)
  )
  
  return(rr_results)
}

province_rr_mean <- list()  # Store results for mean-based sensitivity analysis

for (prov in unique(LBW_clean$province)) {
  cat("\n--- Processing Province:", prov, "---\n")
  
  # Subset data for the current province
  prov_data <- subset(LBW_clean, province == prov)
  
  # Generate RR estimates using mean temperature for centering
  rr_mean <- predict_rr_mean(cb_polynomial, weighted_coefs, weighted_vcov, prov_data, prov)
  
  # Add province information
  rr_mean$Province <- prov
  
  # Store the results
  province_rr_mean[[prov]] <- rr_mean
}

# Combine all provinces' results into a single data frame
all_rr_mean_results <- do.call(rbind, province_rr_mean)


# Add a column to differentiate between median and mean
all_full_rr_results <- all_full_rr_results %>%
  mutate(Method = "Median")

all_rr_mean_results <- all_rr_mean_results %>%
  mutate(Method = "Mean")

# Combine the two datasets
combined_rr_results <- bind_rows(all_full_rr_results, all_rr_mean_results)

# Visualize the comparison for each province
library(ggplot2)
comparison_plots <- list()

for (prov in unique(combined_rr_results$Province)) {
  plot_data <- combined_rr_results %>% filter(Province == prov)
  
  p <- ggplot(plot_data, aes(x = Temperature, y = RR, color = Method, fill = Method)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2, color = NA) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    labs(
      title = paste("Exposure-Response Curve for", prov),
      x = "Temperature (°C)",
      y = "Relative Risk (RR)"
    ) +
    theme_classic() +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  comparison_plots[[prov]] <- p
}

# Arrange all province-specific plots in a grid
library(patchwork)
final_comparison_plot <- wrap_plots(comparison_plots, ncol = 2)

# Save and display the plot
ggsave("Province_ER_Curve_Comparison_Median_vs_Mean.png", final_comparison_plot, width = 16, height = 12)
print(final_comparison_plot)


#########################
#STAGE 6: HEAT VULNERABLILITY INDEX
#########################


dist_births <- fread("/Users/fati0049/Desktop/My_data_2024/Research/Zohra1/Rcode/analysis/dist_births_mor_mpi.csv")
# Specify the columns to keep
columns_to_keep <- c("ADM2_EN", "ADM1_EN.x", 
                     "avg_pm", "avg_tmean", "avg_heat_index", 
                     "avg_humidity", "avg_prec", "SUM", 
                     "IHME_LMICS", "mpi")

# Select only the specified columns
dist_births <- dist_births %>%
  dplyr::select(all_of(columns_to_keep))



merged_data <- inner_join(dist_births, rr_99th_percentile, by = c("ADM1_EN.x" = "Province"))
colnames(merged_data)


# Normalize numeric variables for HVI calculation
numeric_vars <- c("avg_pm", "avg_tmean", "avg_heat_index", "avg_humidity", "avg_prec", "IHME_LMICS", "mpi", "SUM", "RR")
merged_data <- merged_data %>%
  mutate(across(all_of(numeric_vars), ~ scale(.x, center = TRUE, scale = TRUE)))

# Step 3: Monte Carlo Simulation for HVI
set.seed(123)  # For reproducibility
num_samples <- 1000  # Number of Monte Carlo simulations

# Generate Monte Carlo samples for RR
rr_samples <- lapply(1:num_samples, function(i) {
  rnorm(nrow(merged_data), mean = merged_data$RR, sd = (merged_data$UCI - merged_data$LCI) / (2 * qnorm(0.975)))
})

# Define weights for HVI factors
weight_temperature <- 0.20
weight_pm <- 0.10
weight_IHME <- 0.20
weight_mpi <- 0.20
weight_rr <- 0.20
weight_other <- 0.10

# Calculate HVI for each simulation
for (i in 1:num_samples) {
  col_name <- paste0("HVI_", i)
  merged_data[[col_name]] <- rowSums(cbind(
    merged_data$avg_heat_index * weight_temperature,
    merged_data$avg_pm * weight_pm,
    merged_data$IHME_LMICS * weight_IHME,
    merged_data$mpi * weight_mpi,
    rr_samples[[i]] * weight_rr
  ))
}


# Step 4: Calculate Mean HVI Across Simulations
# Calculate mean HVI across all simulation columns
hvi_data <- merged_data %>%
  mutate(mean_HVI = rowMeans(across(starts_with("HVI_")), na.rm = TRUE)) %>%  # Compute mean HVI
  dplyr::select(ADM2_EN, mean_HVI)  # Keep only ADM2_EN and mean_HVI



# Specify the columns to remove
districts_shapefile <- st_read("/Users/fati0049/Desktop/My_data_2024/Research/Zohra1/Rcode/analysis/pak_adm_wfp_20220909_shp/pak_admbnda_adm2_wfp_20220909.shp")
#connect back to shapefile
district_data <- merge(districts_shapefile, hvi_data, by.x = "ADM2_EN", by.y = "ADM2_EN", all.x = TRUE)
#connect back to other env variables
final_district_data <- merge(district_data, dist_births, by.x = "ADM2_EN", by.y = "ADM2_EN", all.x = TRUE)
#remove unnecessary columns
columns_to_remove <- c("Shape_Leng", "Shape_Area", "ADM2_PCODE", "ADM2_REF", 
                       "ADM2ALT1EN", "ADM2ALT2EN", "ADM1_EN", "ADM1_PCODE", 
                       "ADM0_EN", "ADM0_PCODE", "date", "validOn", "validTo")

# Remove the specified columns from the dataset
final_district_data <- final_district_data %>%
  dplyr::select(-all_of(columns_to_remove))

colnames(final_district_data)


##########################
#6.1 Data visualisation
#######################




# Ensure the dataset is in sf format
final_district_data <- st_as_sf(final_district_data)

# Load province-level shapefile (if available)
province_shapefile <- st_read("/Users/fati0049/Desktop/My_data_2024/Research/Zohra1/Rcode/analysis/pak_adm_wfp_20220909_shp/pak_admbnda_adm1_wfp_20220909.shp")

# Add Province Boundaries
province_boundaries <- st_as_sf(province_shapefile)

# Calculate centroids for province labels
province_centroids <- province_boundaries %>%
  mutate(centroid = st_centroid(geometry)) %>%
  st_as_sf() %>%
  mutate(
    lon = st_coordinates(centroid)[, 1],
    lat = st_coordinates(centroid)[, 2]
  )

# Update the plot with grid lines on gray background
hvi_plot <- ggplot(data = final_district_data) +
  geom_sf(aes(fill = mean_HVI), color = NA) +  # District-level HVI map
  geom_sf(data = province_boundaries, fill = NA, color = "black", size = 0.8) +  # Province boundaries
  geom_text(
    data = province_centroids %>%
      filter(!ADM1_EN %in% c("Azad Kashmir", "Islamabad")),  # Exclude Azad Kashmir and Islamabad
    aes(
      x = ifelse(ADM1_EN == "Khyber Pakhtunkhwa", lon - 1.0, lon),  # Adjust placement of Khyber Pakhtunkhwa
      y = ifelse(ADM1_EN == "Khyber Pakhtunkhwa", lat - 0.5, lat),  # Adjust placement of Khyber Pakhtunkhwa
      label = ADM1_EN
    ),
    size = 5, fontface = "bold", color = "white"  # All labels in white
  ) +
  scale_fill_viridis_c(
    option = "turbo",
    na.value = "gray",
    name = "heat vulnerability index",  # Set legend title
    breaks = c(-1, 1),
    labels = c("Low", "High")
  ) +
  theme_light() +
  theme(
    panel.background = element_rect(fill = "gray40", color = NA),  # Change background to light gray
    legend.position = c(0.8, 0.2),  # Position legend inside map
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    panel.grid.major = element_line(color = "white", size = 0.2),  # Thin white major grid lines
    panel.grid.minor = element_blank()  # No minor grid lines
  )

# Save the updated map
ggsave(
  filename = "HVI_Map_Final.png",
  plot = hvi_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# Display the updated map
print(hvi_plot)



# Smaller environmental variable plots
plot_pm <- ggplot(data = final_district_data) +
  geom_sf(aes(fill = avg_pm)) +
  scale_fill_viridis_c(option = "plasma", na.value = "gray") +
  labs(title = "Particulate Matter (PM2.5)", fill = "PM2.5") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

plot_precip <- ggplot(data = final_district_data) +
  geom_sf(aes(fill = avg_prec)) +
  scale_fill_viridis_c(option = "viridis", na.value = "gray") +
  labs(title = "Precipitation (mm)", fill = "Precipitation") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

plot_births <- ggplot(data = final_district_data) +
  geom_sf(aes(fill = SUM)) +
  scale_fill_viridis_c(option = "mako", na.value = "gray") +
  labs(title = "Live Births", fill = "Births") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

plot_mpi <- ggplot(data = final_district_data) +
  geom_sf(aes(fill = mpi)) +
  scale_fill_viridis_c(option = "inferno", na.value = "gray") +
  labs(title = "Poverty Index (MPI)", fill = "MPI") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

plot_deaths <- ggplot(data = final_district_data) +
  geom_sf(aes(fill = IHME_LMICS)) +
  scale_fill_viridis_c(option = "cividis", na.value = "gray") +
  labs(title = "Probability of Deaths", fill = "Deaths") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

plot_heat_index <- ggplot(data = final_district_data) +
  geom_sf(aes(fill = avg_heat_index)) +
  scale_fill_viridis_c(option = "turbo", na.value = "gray") +
  labs(title = "Heat Index (°C)", fill = "Heat Index") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

# Combine all smaller plots into a panel
bottom_panel <- (plot_pm | plot_precip | plot_births) / (plot_heat_index | plot_mpi | plot_deaths)

# Save the paneled map
ggsave(
  filename = "Paneled_Environmental_Maps.png",
  plot = bottom_panel,
  width = 18,
  height = 12,
  dpi = 300
)

# Display the paneled map
print(bottom_panel)


#############################
#STAGE 7: SUBGROUP ANALYSIS
#######################

# Calculate global scaling parameters from the full dataset
mean_tmean <- mean(LBW_clean$tmean_C, na.rm = TRUE)
sd_tmean <- sd(LBW_clean$tmean_C, na.rm = TRUE)

cat("Mean temperature (tmean_C):", mean_tmean, "\n")
cat("Standard deviation of temperature (tmean_C):", sd_tmean, "\n")

# Update subgroups in `LBW_clean`
# Recode Education: Combine No Education and Primary as "Low Education", and Secondary and Higher as "Secondary+"
LBW_clean <- LBW_clean %>%
  mutate(Education_Group = case_when(
    Education %in% c("No education", "Primary") ~ "Low Education",
    Education %in% c("Secondary", "Higher") ~ "Secondary+",
    TRUE ~ NA_character_  # Explicitly mark unexpected values as NA
  ))

# Recode Wealth Index: Combine Poorest and Poorer as "Low Wealth", and Middle, Richer, and Richest as "High Wealth"
LBW_clean <- LBW_clean %>%
  mutate(Wealth_Group = case_when(
    Wealth_Index %in% c("Poorest", "Poorer", "Middle") ~ "Low Wealth",
    Wealth_Index %in% c("Richer", "Richest") ~ "High Wealth"
  ))

# Recode PM2.5: <25 as "Fair", >=25 as "Poor/Hazardous"
LBW_clean <- LBW_clean %>%
  mutate(PM25_Group = ifelse(PM25 < 25, "Fair", "Poor/Hazardous"))

# Define updated subgroup analysis function with scaling
run_weighted_subgroup_analysis <- function(subgroup_data, subgroup_name, weighted_coefs, weighted_vcov, mean_tmean, sd_tmean) {
  if (nrow(subgroup_data) > 10) {  # Ensure sufficient data for analysis
    
    # Scale temperatures within the subgroup
    subgroup_data <- subgroup_data %>%
      mutate(scaled_tmean = (tmean_C - mean_tmean) / sd_tmean)
    
    # Define crossbasis for scaled temperature
    cb_subgroup <- crossbasis(
      subgroup_data$scaled_tmean,
      lag = lag,
      argvar = list(fun = "poly", degree = degree_exposure),
      arglag = list(fun = "poly", degree = degree_lag)
    )
    
    # Calculate median temperature for centering
    median_tmean_scaled <- median(subgroup_data$scaled_tmean, na.rm = TRUE)
    
    # Define specific percentiles for predictions (in scaled space)
    percentiles <- c(0.01, 0.10, 0.90, 0.99)
    tmean_C_percentiles_scaled <- quantile(subgroup_data$scaled_tmean, probs = percentiles, na.rm = TRUE)
    
    # Rescale percentiles back to the original temperature
    tmean_C_percentiles <- tmean_C_percentiles_scaled * sd_tmean + mean_tmean
    
    # Use crosspred with weighted coefficients and variance-covariance matrix
    pred_percentiles <- crosspred(
      cb_subgroup, 
      coef = weighted_coefs, 
      vcov = weighted_vcov, 
      at = tmean_C_percentiles_scaled,  # Percentile predictions in scaled space
      cen = median_tmean_scaled         # Centering at median scaled temperature
    )
    
    # Extract the RR estimates for this subgroup
    rr_at_percentiles <- data.frame(
      Subgroup = subgroup_name,
      Percentile = names(tmean_C_percentiles),  # Use original scale for percentiles
      tmean_C = tmean_C_percentiles,           # Original temperature values
      RR = exp(pred_percentiles$allfit),
      LCI = exp(pred_percentiles$alllow),
      UCI = exp(pred_percentiles$allhigh)
    )
    
    return(rr_at_percentiles)
  } else {
    cat("Skipping subgroup:", subgroup_name, "- insufficient data\n")
    return(NULL)
  }
}

# Define subgroups
education_subgroup1 <- LBW_clean %>% filter(Education_Group == "Low Education")
education_subgroup2 <- LBW_clean %>% filter(Education_Group == "Secondary+")

wealth_subgroup1 <- LBW_clean %>% filter(Wealth_Group == "Low Wealth")
wealth_subgroup2 <- LBW_clean %>% filter(Wealth_Group == "High Wealth")

pm25_subgroup1 <- LBW_clean %>% filter(PM25_Group == "Fair")
pm25_subgroup2 <- LBW_clean %>% filter(PM25_Group == "Poor/Hazardous")

urban_subgroup <- LBW_clean %>% filter(Area == "Urban")
rural_subgroup <- LBW_clean %>% filter(Area == "Rural")

# Run subgroup analyses using weighted coefficients and variance-covariance
rr_education_group1 <- run_weighted_subgroup_analysis(
  education_subgroup1, "Low_Education", weighted_coefs, weighted_vcov, mean_tmean, sd_tmean
)
rr_education_group2 <- run_weighted_subgroup_analysis(
  education_subgroup2, "Secondary_Higher", weighted_coefs, weighted_vcov, mean_tmean, sd_tmean
)

rr_wealth_group1 <- run_weighted_subgroup_analysis(
  wealth_subgroup1, "Low_Wealth", weighted_coefs, weighted_vcov, mean_tmean, sd_tmean
)
rr_wealth_group2 <- run_weighted_subgroup_analysis(
  wealth_subgroup2, "High_Wealth", weighted_coefs, weighted_vcov, mean_tmean, sd_tmean
)

rr_pm25_group1 <- run_weighted_subgroup_analysis(
  pm25_subgroup1, "PM25_Fair", weighted_coefs, weighted_vcov, mean_tmean, sd_tmean
)
rr_pm25_group2 <- run_weighted_subgroup_analysis(
  pm25_subgroup2, "PM25_Poor_Hazardous", weighted_coefs, weighted_vcov, mean_tmean, sd_tmean
)

rr_urban <- run_weighted_subgroup_analysis(
  urban_subgroup, "Urban", weighted_coefs, weighted_vcov, mean_tmean, sd_tmean
)
rr_rural <- run_weighted_subgroup_analysis(
  rural_subgroup, "Rural", weighted_coefs, weighted_vcov, mean_tmean, sd_tmean
)

# Combine all results
rr_subgroup_results <- bind_rows(
  rr_education_group1, 
  rr_education_group2, 
  rr_wealth_group1, 
  rr_wealth_group2, 
  rr_pm25_group1, 
  rr_pm25_group2, 
  rr_urban, 
  rr_rural
)

# Format results for reporting
rr_subgroup_formatted <- rr_subgroup_results %>%
  pivot_wider(
    names_from = Percentile,
    values_from = c(RR, LCI, UCI),
    names_sep = "_"
  ) %>%
  mutate(
    `Extreme Cold (1st Percentile)` = paste0(round(`RR_1%`, 2), " (", round(`LCI_1%`, 2), ", ", round(`UCI_1%`, 2), ")"),
    `Moderate Cold (10th Percentile)` = paste0(round(`RR_10%`, 2), " (", round(`LCI_10%`, 2), ", ", round(`UCI_10%`, 2), ")"),
    `Moderate Heat (90th Percentile)` = paste0(round(`RR_90%`, 2), " (", round(`LCI_90%`, 2), ", ", round(`UCI_90%`, 2), ")"),
    `Extreme Heat (99th Percentile)` = paste0(round(`RR_99%`, 2), " (", round(`LCI_99%`, 2), ", ", round(`UCI_99%`, 2), ")")
  ) %>%
  dplyr::select(
    Subgroup, 
    `Extreme Cold (1st Percentile)`, 
    `Moderate Cold (10th Percentile)`, 
    `Moderate Heat (90th Percentile)`, 
    `Extreme Heat (99th Percentile)`
  )

# Print formatted results
print(rr_subgroup_formatted)


# Save the results into a CSV file
write.csv(rr_subgroup_formatted, "Subgroup_Analysis_Results.csv", row.names = FALSE)

# Print a confirmation message
cat("Subgroup analysis results saved as 'Subgroup_Analysis_Results.csv'\n")

# Display the first few rows of the results for verification
head(rr_subgroup_results)
