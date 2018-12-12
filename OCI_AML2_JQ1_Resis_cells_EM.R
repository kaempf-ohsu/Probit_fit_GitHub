## top of script   --------------------------------------------------------------------------------------



## PURPOSE: Use dummy dose-response data to fit probit and compute AUC and then connect to GitHub remote repo



# script start date: 12/11/18




# Working directory -------------------------------------------------------




# Load packages -----------------------------------------------------------

library(tidyverse)




# Global vars -------------------------------------------------------------

inc <- 0.00001  




## Dose-response data ------------------------------------------------------

# Raw MTS absorbance values from OCI-AML2 cell line on Marctest v13 plate (from Kyle Romine) 


#** JQ1   ------------------------------------------------------------------

Raw_response_vec <- c(0.168,0.183,0.188,0.230,0.313,0.344,0.295)

Concentration_uM_vec <- c(10,10/3,10/(3^2),10/(3^3),10/(3^4),10/(3^5),10/(3^6))


# quick check
(Raw_response_vec - ave_FSV) / (PAC - ave_FSV)


# ** CREATE scalars  ****

PAC <- 0.3250217

ave_FSV <- 0.153


# ** CREATE df  ****
DR_df <- tibble(Raw_viab = Raw_response_vec,
                Norm_viab = (Raw_viab - ave_FSV) / (PAC - ave_FSV),
                Norm_trunc_viab = case_when(Norm_viab < 0 ~ 0,
                                            Norm_viab > 1 ~ 1,
                                            Norm_viab >= 0 & Norm_viab <= 1 ~ Norm_viab),
                conc = Concentration_uM_vec,
                log10_conc = log10(conc))


# check
dim(DR_df)  # 7 by 5


#** CPI-0610   ----------------------------------------------------------------------------------

CPI_response_vec <- c(0.189,0.196,0.260,0.260,0.285,0.294,0.289)



## (1) Probit reg   ---------------------------------------------------------------------------------------------


#** JQ1 probit   -------------------------------------------------------------------------------------------

JQ1_Probit_betas <- suppressWarnings(glm(Norm_trunc_viab ~ log10_conc, family=binomial(link="probit"), 
                                     data=DR_df))$coefficients

# check
length(JQ1_Probit_betas)  # 2



#** JQ1 IC50   -----------------------------------------------------------------


# Max dose for JQ1
max_log_dose <- log10(10)

# Min dose for JQ1
min_log_dose <- log10(10/(3^6))

IC50_log10 <- max_log_dose


log_dose_range <- seq(from=min_log_dose, to=max_log_dose, by=inc)


for (i in log_dose_range) {
  curve_height <- pnorm(q=(JQ1_Probit_betas[1] + i*JQ1_Probit_betas[2]), lower.tail=TRUE) * 100 
  
  if (curve_height <= 50 & i < IC50_log10) {
    IC50_log10 <- i
  }
}


# IC50 should be 0.445
(IC50_uM <- 10^IC50_log10)




#** JQ1 AUC -------------------------------------------------------------


curve_height_100 <- pnorm(q=(JQ1_Probit_betas[1] + log_dose_range*JQ1_Probit_betas[2]),lower.tail=TRUE)*100


# 'x' = abscissa values, 'fx' = corresponding values of f(x)

# AUC should be 150
(AUC <- sfsmisc::integrate.xy(x=log_dose_range, fx=curve_height_100))

# 150.21 


## (2) Cubic reg   -----------------------------------------------------------------------------


# checks
DR_df


#** JQ1 cubic   -------------------------------------------------------------------------------------------


# (1) y = proportion (Viability)
JQ1_Cubic_prop_betas <- lm(formula = Norm_trunc_viab ~ poly(log10_conc, 3), data=DR_df)$coefficients



# (2) y = percentage (Viability)
JQ1_Cubic_pct_betas <- lm(formula = (Norm_trunc_viab*100) ~ poly(log10_conc, 3), data=DR_df)$coefficients


# checks
JQ1_Cubic_pct_betas
length(JQ1_Cubic_pct_betas)  # 4



#** JQ1 AUC   --------------------------------------------------------------------------------------
 
# Max dose for JQ1
max_log_dose <- log10(10)

# Min dose for JQ1
min_log_dose <- log10(10/(3^6))


log_dose_range <- seq(from=min_log_dose, to=max_log_dose, by=inc)


JQ1_cubic_y_vals <- JQ1_Cubic_pct_betas[1] + (JQ1_Cubic_pct_betas[2] * log_dose_range) + 
                      (JQ1_Cubic_pct_betas[3] * (log_dose_range^2)) + (JQ1_Cubic_pct_betas[4] * (log_dose_range^3))

# check
identical(length(log_dose_range),length(JQ1_cubic_y_vals)) # TRUE


# compute AUC
(cubic_AUC <- sfsmisc::integrate.xy(x=log_dose_range, 
                                   fx=JQ1_cubic_y_vals))
# 148.45 

# probit AUC was 150.21


## (3) Linear interpolation   ---------------------------------------------------------------------------------------

# li = (fine) linear interpolation

# kulife::auc() --> Compute the AUC using linear interpolation for two vectors where one corresponds
#                   to the x values and the other corresponds to the y values.

(JQ1_li_AUC <- kulife::auc(x=DR_df$log10_conc,y=(DR_df$Norm_trunc_viab)*100))

# AUC = 153.25

# cubic AUC = 148.45
# probit AUC = 150.21
