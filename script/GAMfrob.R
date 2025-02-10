# load libraries

library(mgcv)
library(ggplot2)
library(pROC)
library(scales)

####

setwd("..")  # Set your project root as the working directory

####

# log10_trans() function
log10_trans <- function() {
  trans <- function(x) log10(x)
  inv <- function(x) 10^x
  scales::trans_new(
    name = "log10",
    transform = trans,
    inverse = inv,
    breaks = scales::extended_breaks(),
    minor_breaks = scales::log_breaks(),
    domain = c(1e-100, Inf)
  )
}

# load Data

####
Data=read.csv("data/datafrob.csv")
x=c("bio1","bio12","tconmean","lenferr","lenroadr","shdi","popde","urbr","popdi","numazr","redre","coorx","coory")
var.names=c("bio1","bio12","tconmean","lenferr","lenroadr","shdi","popde","urbr","popdi","numazr","redre","coorx","coory")
n=nrow(Data)
n0=sum(Data$frob==0)
p = length(var.names) # Number of covariates

Data$nonzerofrob=ifelse(Data$frob == 0, 0, 1)
Data$nonzerofrob=as.integer(Data$nonzerofrob)

names(Data)

##########################################
##########################################


Binomial_GAM_W_rob  <- gam(nonzerofrob ~ s(bio1,bs="cr") +  s(bio12,bs="cr")  +  
                       s(tconmean,bs="cr") + s(lenferr,bs="cr") +  
                       s(lenroadr,bs="cr")  +  s(shdi,bs="cr")  + s(popde,bs="cr")+s(urbr,bs="cr")+s(popdi,bs="cr")+s(numazr,bs="cr")+s(redre,bs="cr")+s(coorx,coory, bs = "gp", m = 2),                    
                     data=Data,family= binomial(link="logit"),method="REML",weights=areagis_ha/mean(areagis_ha))
summary(Binomial_GAM_W_rob)

# AUC calculation

# 1. Calculate Predicted Probabilities
predicted_probs <- predict(Binomial_GAM_W_rob, type = "response")  # Replace 'Model' with your GAM model

# 2. Calculate AUC
roc_obj <- roc(Data$nonzerofrob, predicted_probs)  # Replace 'Data' with your dataset

# Display AUC
auc_value <- auc(roc_obj)
cat("AUC:", auc_value, "\n")

# 3. Create ROC Curve (optional)
plot(roc_obj, main = "ROC Curve")

# give predicted value for bio1 

# Create a data frame with the desired values
new_data <- data.frame(
  bio1 = 6,                    # Set bio1 to the desired value
  bio12 = median(Data$bio12),  # Set other variables to their medians
  tconmean = median(Data$tconmean),
  lenferr = median(Data$lenferr),
  lenroadr = median(Data$lenroadr),
  shdi = median(Data$shdi),
  popde = median(Data$popde),
  urbr = median(Data$urbr),
  popdi = median(Data$popdi),
  numazr = median(Data$numazr),
  redre = median(Data$redre),
  coorx = median(Data$coorx),
  coory = median(Data$coory)
)

# Predict and Produce confidence intervals for the specific value of bio1
Preds_nzero_specific <- predict(Binomial_GAM_W_rob, newdata = new_data, type = "link", se = TRUE, unconditional = TRUE)
fit.link_specific <- Preds_nzero_specific$fit
se_specific <- Preds_nzero_specific$se
CI.L_specific <- fit.link_specific - 1.96 * se_specific # Lower limit of the 95% confidence interval
CI.R_specific <- fit.link_specific + 1.96 * se_specific # Upper limit of the 95% confidence interval

# Transform the link value back to the original scale (probability)
predicted_probability_specific <- exp(fit.link_specific) / (1 + exp(fit.link_specific))
CI.L_prob_specific <- exp(CI.L_specific) / (1 + exp(CI.L_specific))
CI.R_prob_specific <- exp(CI.R_specific) / (1 + exp(CI.R_specific))

# Create a data frame with the predicted value and confidence intervals
specific_prediction_df <- data.frame(
  Predicted = predicted_probability_specific,
  CI_Lower = CI.L_prob_specific,
  CI_Upper = CI.R_prob_specific
)

# Print or inspect the specific_prediction_df
print(specific_prediction_df)

num_repetitions <- 30  # You can adjust this as needed

# Initialize a vector to store AUC values
auc_values <- numeric(num_repetitions)

# Perform the procedure multiple times
for (i in 1:num_repetitions) {
  # Step 1: Create a random train-test split (75% train, 25% test)
  set.seed(i)  # Set a random seed for reproducibility
  sample_indices <- sample(nrow(Data), nrow(Data) * 0.75)
  train_data <- Data[sample_indices, ]
  test_data <- Data[-sample_indices, ]
  
  # Step 2: Fit the model M1 to the training data (adjust your model formula as needed)
  M1 <- gam(nonzerofrob ~ s(bio1,bs="cr") +  s(bio12,bs="cr")  +  
              s(tconmean,bs="cr") + s(lenferr,bs="cr") +  
              s(lenroadr,bs="cr")  +  s(shdi,bs="cr")  + s(popde,bs="cr")+s(urbr,bs="cr")+s(popdi,bs="cr")+s(numazr,bs="cr")+s(redre,bs="cr")+s(coorx,coory, bs = "gp", m = 2),                    
            data=train_data,family= binomial(link="logit"),method="REML",weights=areagis_ha/mean(areagis_ha))
  
  # Step 3: Use the fitted model to predict responses in the test data
  predicted_probs <- predict(M1, newdata = test_data, type = "response")
  
  # Step 4: Calculate the AUC and store it in the vector
  roc_obj <- roc(test_data$nonzerofrob, predicted_probs)
  auc_values[i] <- auc(roc_obj)
}

# Calculate the mean AUC value and assess predictive ability
mean_auc <- mean(auc_values)

# prepare for plotting all significant variables (excluding coorx,coory)

model_summary <- summary(Binomial_GAM_W_rob)

# Extract the p-values of the smooth terms
smooth_terms_pvalues <- model_summary$s.table[, "p-value"]

# Define a significance level, e.g., 0.05
significance_level <- 0.05

# Identify significant variables
significant_vars <- names(which(smooth_terms_pvalues < significance_level))

# Remove 's()', 'coorx' and 'coory' from the names
significant_vars <- gsub("s\\((.*?)\\)", "\\1", significant_vars)
significant_vars <- setdiff(significant_vars, c("coorx,coory"))

# Initialize lists for storing plot data
gg_x <- list()
gg_original <- list()

# Number of data points for prediction
nn <- 3 * 10^4

# Create a baseline data frame with median values for all predictors
# Including 'coorx' and 'coory' as they are required for the model prediction
all_vars <- setdiff(names(Data), "nonzerofrob")
baseline_data <- as.data.frame(lapply(Data[, all_vars], median))

# Loop over significant variables
for (var in significant_vars) {
  # Create a copy of the baseline data for modification
  R <- baseline_data[rep(1, nn), ]
  
  # Replace the column of the significant variable with a sequence across its range
  a <- quantile(Data[, var], 0.001)
  b <- quantile(Data[, var], 0.999)
  R[, var] <- seq(a, b, length.out = nn)
  
  # Predict and produce confidence intervals
  Preds_nzero <- predict(Binomial_GAM_W_rob, newdata = R, type = "link", se = T, unconditional = T) 
  fit.link <- Preds_nzero$fit 
  se <- Preds_nzero$se  
  CI.L <- fit.link - 2 * se 
  CI.R <- fit.link + 2 * se 
  CI <- cbind(fit.link, CI.L, CI.R) 
  CI <- exp(CI) / (1 + exp(CI)) # Convert to probabilities
  colnames(CI) <- c("Predictions", "CI_L", "CI_R")
  
  # Store data for plotting
  gg_x[[var]] <- data.frame(cbind(CI, x = R[, var]), var = rep(var, nn))
  
  # Original data for rug plots
  I <- Data[, var] > a & Data[, var] < b
  gg_original[[var]] <- data.frame(x = Data[I, var], y = Data[I, "nonzerofrob"], var = rep(var, sum(I)))
}

# Combine data for facets
ggg_x <- do.call(rbind, gg_x)
ggg_original <- do.call(rbind, gg_original)

# Plot via ggplot2

pdf("Binomial_GAM_rob.pdf", width = 6, height = 6) # Fig. 4

# Adjust font sizes
base_font_size <- 14  # Base font size, adjust as needed

ggplot(ggg_x, aes(x = x, y = Predictions, group = var)) +
  geom_rug(data = filter(ggg_original, y == 0), aes(x = x, y = y), color = "gray75", size = 0.15, alpha = 0.5, sides = "b") +
  geom_rug(data = filter(ggg_original, y == 1), aes(x = x, y = y), color = "gray75", size = 0.15, alpha = 0.5, sides = "t") +
  geom_ribbon(aes(ymin = CI_L, ymax = CI_R), fill = "goldenrod", alpha = 0.5) +
  geom_line(aes(x = x, y = Predictions), col = "black", size = 0.5) +
  ylab(expression(paste("Probability of ", ~f[rob] > 0, sep = ""))) +
  xlab("") +
  theme_bw(base_size = base_font_size) +
  theme(legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title.y = element_text(vjust = 0.1, margin = margin(0, 10, 0, 0)),
        axis.title.x = element_text(vjust = -0.25),
        axis.text = element_text(size = base_font_size),
        axis.title = element_text(size = base_font_size),
        strip.text = element_text(size = base_font_size),
        plot.margin = margin(2, 2, 2, 2, "mm"),  # Reduce plot margins
        panel.spacing = unit(1, "mm")) +  # Reduce spacing between panels
  facet_wrap(~var, scales = "free_x", nrow = 2, labeller = label_parsed, strip.position = "bottom") 
dev.off()


#######################################################################################
## The Second Stage
### Model  frob when frob > 0 throug Beta_GAM with logistic link
###
r <- 35
Beta_GAM_W_rob    <- gam(frob ~ s(bio1,bs="cr",k=r) +  s(bio12,bs="cr",k=r)  +  s(tconmean,bs="cr",k=r) + 
                     s(lenferr,bs="cr",k=r) +
                     s(lenroadr,bs="cr",k=r)  +  s(shdi,bs="cr",k=r)   +  s(popde,bs="cr",k=r) + s(urbr,bs="cr",k=r)+s(popdi,bs="cr",k=r)+
                     s(numazr,bs="cr",k=r)+s(redre,bs="cr",k=r)+s(coorx,coory, bs = "gp", m = 2),weights=areagis_ha/mean(areagis_ha),
                   subset=frob>0,data=Data,family=betar(link="logit"),method="REML")

summary(Beta_GAM_W_rob)


###########################
# prepare for plotting

# prepare for plotting all significant variables (excluding coorx,coory)

model_summary <- summary(Beta_GAM_W_rob)

# Extract the p-values of the smooth terms
smooth_terms_pvalues <- model_summary$s.table[, "p-value"]

# Define a significance level, e.g., 0.05
significance_level <- 0.05

# Identify significant variables
significant_vars <- names(which(smooth_terms_pvalues < significance_level))

# Remove 's()', 'coorx' and 'coory' from the names
significant_vars <- gsub("s\\((.*?)\\)", "\\1", significant_vars)
significant_vars <- setdiff(significant_vars, c("coorx,coory"))

# Initialize lists for storing plot data
gg_x <- list()
gg_original <- list()

# Number of data points for prediction
nn <- 3 * 10^4

# Create a baseline data frame with median values for all predictors
# Including 'coorx' and 'coory' as they are required for the model prediction
all_vars <- setdiff(names(Data), "frob")
baseline_data <- as.data.frame(lapply(Data[, all_vars], median))

# Loop over significant variables
for (var in significant_vars) {
  # Create a copy of the baseline data for modification
  R <- baseline_data[rep(1, nn), ]
  
  # Replace the column of the significant variable with a sequence across its range
  a <- quantile(Data[, var], 0.001)
  b <- quantile(Data[, var], 0.999)
  R[, var] <- seq(a, b, length.out = nn)
  #
  # Predict and Produce confidence intervals:
  #
  Preds_frob  <- predict(Beta_GAM_W_rob,newdata = R,type="link",se=T,unconditional=T) 
  fit.link    <- Preds_frob$fit 
  se          <- Preds_frob$se 
  CI.L        <- fit.link-2*se 
  CI.R        <- fit.link+2*se 
  CI          <- cbind(fit.link,CI.L,CI.R) 
  CI          <- exp(CI)/(1+exp(CI)) # The first column corresponds to the estimated average of frob when frob > 0
  colnames(CI) <- c("Predictions","CI_L","CI_R")
  # Store data for plotting
  gg_x[[var]] <- data.frame(cbind(CI, x = R[, var]), var = rep(var, nn))
  
  # Original data for rug plots
  I <- Data[, var] > a & Data[, var] < b
  gg_original[[var]] <- data.frame(x = Data[I, var], y = Data[I, "frob"], var = rep(var, sum(I)))
}

# Combine data for facets
ggg_x <- do.call(rbind, gg_x)
ggg_original <- do.call(rbind, gg_original)

# Plot via ggplot2

pdf("Beta_GAM_rob.pdf", width = 6, height = 6) # Fig. 6

# Increase font sizes for better readability
base_font_size <- 14  # Base font size, adjust as needed

ggplot(ggg_x, aes(x = x, y = Predictions, group = var)) +
  geom_hex(data = ggg_original, fill = "gray75", size = 0.15, alpha = 0.5, bins = 75,
           aes(x = x, y = y, alpha = ..count..)) +
  scale_alpha_continuous(trans = log10_trans()) +
  geom_ribbon(aes(ymin = CI_L, ymax = CI_R), fill = "goldenrod", alpha = 0.5) +
  geom_line(col = "black", size = 0.5) +
  ylab(expression(paste("Average of ", ~f[rob], " given that ", ~f[rob] > 0, sep = ""))) +
  xlab("") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  theme_bw(base_size = base_font_size) +
  theme(legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title.y = element_text(vjust = 0.1, margin = margin(0, 10, 0, 0)),
        axis.title.x = element_text(vjust = -0.25),
        axis.text = element_text(size = base_font_size),
        axis.title = element_text(size = base_font_size),
        strip.text = element_text(size = base_font_size),
        plot.margin = margin(2, 2, 2, 2, "mm"),
        panel.spacing = unit(1, "mm")) +
  facet_wrap(~var, scales = "free_x", ncol = 2, labeller = label_parsed, strip.position = "bottom")
dev.off()


### Hurdle model:Putting things together
###
### Fitted values:
####
#### Produce hurdle plots, using delta method
#### 

###########################
# prepare for plotting

# Function to extract significant variables from GAM model summary
extract_significant_vars <- function(model_summary) {
  # Extract the p-values of the smooth terms
  smooth_terms_pvalues <- model_summary$s.table[, "p-value"]
  
  # Define a significance level, e.g., 0.05
  significance_level <- 0.05
  
  # Identify significant variables
  significant_vars <- names(which(smooth_terms_pvalues < significance_level))
  
  # Remove 's()', 'coorx' and 'coory' from the names
  significant_vars <- gsub("s\\((.*?)\\)", "\\1", significant_vars)
  significant_vars <- setdiff(significant_vars, c("coorx,coory"))
  
  return(significant_vars)
}

# Extract significant variables from both models
significant_vars_beta <- extract_significant_vars(summary(Beta_GAM_W_rob))
significant_vars_binomial <- extract_significant_vars(summary(Binomial_GAM_W_rob))

# Combine and get unique list of significant variables from both models
all_significant_vars <- unique(c(significant_vars_beta, significant_vars_binomial))


#### 
gg          <-list()
gg_x        <- list()
gg_original <-list()
for(i in 1:p){
  nn     = 3*10^4+1; 
  R      = matrix(apply(Data[,x],2,median),nrow=1); 
  R      = R%x% rep(1,nn);colnames(R) = x; 
  R      = as.data.frame(R)
  a      = quantile(Data[,x[i]],0.001); b= quantile(Data[,x[i]],0.999); I = Data[,x[i]] > a & Data[,x[i]] < b
  R[,i]  = seq(a,b,length=nn)
  #
  fit_Binomial = predict(Binomial_GAM_W_rob,newdata=R,type="response",se=T,unconditional=T) 
  fit_Beta     = predict(Beta_GAM_W_rob,newdata=R,type="response",se=T,unconditional=T)     
  mu_Binomial  = fit_Binomial$fit
  mu_Beta      = fit_Beta$fit
  se_Binomial  = fit_Binomial$se 
  se_Beta      = fit_Beta$se
  ##
  mu_Hurdle    = mu_Binomial*mu_Beta
  se_Hurdle    = sqrt(se_Binomial^2*mu_Beta^2  + mu_Binomial^2*se_Beta^2 + se_Binomial^2*se_Beta^2)
  ## 
  phi <- Beta_GAM_W_rob$family$getTheta(TRUE)
  sd.y       <-  sqrt(mu_Binomial*(mu_Beta*(1-mu_Beta)/(1+phi) + (1-mu_Binomial)*mu_Beta^2))
  #
  CI.L        <- mu_Hurdle-2*se_Hurdle 
  CI.R        <- mu_Hurdle+2*se_Hurdle 
  CI          <- cbind(mu_Hurdle,CI.L,CI.R,sd.y) 
  colnames(CI) <- c("Predictions","CI_L","CI_R","SD")
  ##
  gg_x[[i]] <- data.frame(cbind(CI,x=R[,i]),var = rep(var.names[i],nn))
  gg_original[[i]] <- data.frame(x=Data[I,x[i]],y=Data[I,"frob"],var = rep(var.names[i],sum(I)))
}        

# put altogether for facets
ggg_x<-c()
for(i in 1:p){
  ggg_x <- rbind(ggg_x,gg_x[[i]])
}  
ggg_original <- c()
for(i in 1:p){
  ggg_original <- rbind(ggg_original,gg_original[[i]])
}  


# Plot
pdf("Hurdle_GAM_rob.pdf", width = 6, height = 6)  # Adjust the PDF size # Fig. S1

# Increase font sizes for better readability
base_font_size <- 12  # Base font size, adjust as needed

ggplot(ggg_x[ggg_x$var %in% all_significant_vars, ], aes(x = x, y = Predictions, group = var, fill = var)) +
  geom_hex(data = ggg_original[ggg_original$var %in% all_significant_vars, ], fill = "gray75", size = 0.15, alpha = 0.5, bins = 75, aes(x = x, y = y, alpha = ..count..)) +
  geom_line(aes(x = x, y = SD), size = 0.5, color = "darkgreen") +
  geom_ribbon(aes(ymin = CI_L, ymax = CI_R), fill = "goldenrod", alpha = 0.5) +
  geom_line(color = "black", size = 0.5) +
  ylab(expression(paste(~f[rob], sep = ""))) +
  xlab("") +
  theme(legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title.y = element_text(vjust = 0.1, margin = margin(0, 10, 0, 0)),
        axis.title.x = element_text(vjust = -0.25),
        axis.text = element_text(size = base_font_size),
        axis.title = element_text(size = base_font_size),
        strip.text = element_text(size = base_font_size),
        plot.margin = margin(2, 2, 2, 2, "mm"),
        panel.spacing = unit(1, "mm"),
        # Modify the panel background and border
        panel.background = element_rect(fill = "white", color = "black"),  # White background, black border
        panel.grid.major = element_line(color = "gray75", size = 0.15),  # Set grid color and size
        panel.grid.minor = element_line(color = "gray75", size = 0.15),  # Set grid color and size
        panel.border = element_rect(fill = NA, color = "black"),  # Black border around the panel
        strip.background = element_rect(fill = "lightgray", color = "black")) +  # Lighter grey background with a black border for the variable name panel
  facet_wrap(~var, scales = "free_x", ncol = 2, labeller = label_parsed, strip.position = "bottom") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.25))  # Set limits from 0 to 1

dev.off()

####

# Function to extract chi-square statistics, edf, and calculate p-values for smooth terms
extract_gam_info <- function(model) {
  s_table <- summary(model)$s.table
  chi_sq <- s_table[, "Chi.sq"]
  edf <- s_table[, "edf"]
  # Calculate p-values from chi-square statistics and edf
  p_values <- 1 - pchisq(chi_sq, df = edf)
  return(list(chi_sq = chi_sq, edf = edf, p_values = p_values))
}

# Define the function to filter significant variables

filter_significant <- function(info, threshold = 0.05) {
  # Extract p-values
  p_values <- info$p_values
  
  # Identify significant terms
  significant_vars <- names(p_values)[p_values < threshold]
  
  return(significant_vars)
}

# Extract and filter information from both models
info_M1 <- extract_gam_info(Binomial_GAM_W_rob)
significant_vars_M1 <- filter_significant(info_M1)
info_M2 <- extract_gam_info(Beta_GAM_W_rob)
significant_vars_M2 <- filter_significant(info_M2)

# Combine significant variables from both models
significant_vars <- unique(c(significant_vars_M1, significant_vars_M2))

# Function to create a results dataframe with only significant variables
create_results_df <- function(info, model_name, significant_vars) {
  data.frame(
    Model = model_name,
    Predictor = names(info$chi_sq)[names(info$chi_sq) %in% significant_vars],
    Chi_Square = info$chi_sq[names(info$chi_sq) %in% significant_vars],
    EDF = info$edf[names(info$edf) %in% significant_vars],
    P_Value = info$p_values[names(info$p_values) %in% significant_vars]
  )
}

# Create results dataframes for significant variables of both models
results_M1 <- create_results_df(info_M1, "Binomial_GAM_W_rob", significant_vars)
results_M2 <- create_results_df(info_M2, "Beta_GAM_W_rob", significant_vars)

# Combine chi-square statistics and subtract sum of edf
combined_chi_sq <- info_M1$chi_sq + info_M2$chi_sq
degrees_of_freedom <- info_M1$edf + info_M2$edf
adjusted_stats <- combined_chi_sq - degrees_of_freedom

# Create combined results dataframe
combined_results <- data.frame(
  Predictor = names(adjusted_stats),
  Combined_Chi_Square = adjusted_stats
)

# Print results
print("Results
for Binomial_GAM_W_rob:")
print(results_M1)
print("Results for Beta_GAM_W_rob:")
print(results_M2)
print("Combined Results for Both Models:")
print(combined_results)


# Combine the row names and chi-square values for plotting
plot_data <- rbind(
  data.frame(Model = "Binomial", Variable = rownames(results_M1), Chi_Square = results_M1$Chi_Square),
  data.frame(Model = "Beta", Variable = rownames(results_M2), Chi_Square = results_M2$Chi_Square),
  data.frame(Model = "Hurdle", Variable = rownames(combined_results), Chi_Square = combined_results$Combined_Chi_Square)
)

# Clean variable names by removing "s()" or similar prefixes
clean_variable_name <- function(name) {
  gsub("s\\((.*?)\\)", "\\1", name)
}
plot_data$Variable <- sapply(plot_data$Variable, clean_variable_name)

# Specify the order of facets
plot_data$Model <- factor(plot_data$Model, levels = c("Binomial", "Beta", "Hurdle"))

# Sort the data by Combined_Chi_Square in descending order
plot_data <- plot_data[order(plot_data$Model, -plot_data$Chi_Square), ]


# retain only significant variables from the models
plot_data <- plot_data[plot_data$Variable %in% c(all_significant_vars, "coorx,coory"), ]

# plot

p <- ggplot(plot_data, aes(x = Chi_Square, y = reorder(Variable, Chi_Square))) +
  geom_bar(stat = "identity", fill = "goldenrod", color = "black", alpha = 0.5) +
  facet_wrap(~ Model, scales = "free_x") +
  labs(x = element_text("Test statistic", size = 14)) +  # Modify the size here
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        panel.spacing = unit(0.5, "lines"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(10, 10, 10, 10)) +
  theme(plot.title = element_text(hjust = 0.5))

# Define the file path and name for the plot
file_path <- "Wald_rob.png" # Fig. 8

# Save the plot to the specified file path
ggsave(filename = file_path, plot = p, device = "png", width = 15, height = 6)
 
 # Print a message confirming the save
cat("Plot saved to:", file_path, "\n")





