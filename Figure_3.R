########################################## SURVIVAL ANALYSIS #################################################

#################################################################################
#Survival

library(openxlsx)
library(survival)
library(survminer)
library(dplyr)
###############################################
############Load data##########################

mat <- data
mat_df <- as.data.frame(mat)

# Subset NR and R patents 

NR <- mat_df[mat_df$Response == "NR", ]
R <- mat_df[mat_df$Response == "R", ]

# PARAMETERS
#docDir <- file.path("~/Research/Gross/Inflammasome_070223/doc/")

survDir <- file.path("/Users/...../Survival")
#dir.create(survDir)

kaplanDir <- file.path("/Users/....../Kaplan")
#dir.create(kaplanDir)

coxDir <- file.path("/Users/............/Cox")
#dir.create(coxDir)

# LOAD DATA
#setwd(docDir)
#ann.sample <- read.xlsx("TCGA_LAML_annotation_Q_PAM2.xlsx")

ann.sample <- as.data.frame(mat_df)

ann.sample$TIME <- ann.sample$Overall_survival 

#ann.sample$TIME[is.na(ann.sample$TIME)] <- ann.sample$days_to_last_followup[is.na(ann.sample$TIME)]

ann.sample$STATUS <- ifelse(ann.sample$Alive == "Yes", 1, 2)
#ann.sample$SEX <- ifelse(ann.sample$gender == "FEMALE", 1, 2)
#ann.sample$AGE <- ann.sample$age_at_initial_pathologic_diagnosis
#ann.sample$AGEBIN <- ifelse(ann.sample$age_at_initial_pathologic_diagnosis <= 60, 1, 2)
#zs <- scale(ann.sample[, c(10:47)])
#qt <- ann.sample[, grepl("_Q$", colnames(ann.sample))]

################################################################################
# KAPLAN MEIER

survInput <- ann.sample[, c("TIME", "STATUS", "Therapy", "Imm.Subtype", "Macrophages.M2", "LA_TAM", "NonLA_TAM")]

# SEX, AGEBIN
covariates <- c("Therapy", "Imm.Subtype", "Macrophages.M2", "LA_TAM", "NonLA_TAM")

survInput$TIME <- as.numeric(survInput$TIME)

setwd(kaplanDir)

for(cvt in covariates){
  survInput.sub <-  cbind(CVT = survInput[, cvt],
                          survInput[, ])
  fit <- survfit(Surv(TIME, STATUS) ~ CVT, data = survInput.sub)
  #pv <- surv_pvalue(fit)
  
  p <- ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = survInput.sub,             # data used to fit survival curves.
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = FALSE,         # show confidence intervals for 
    # point estimates of survival curves.
    xlim = c(0,500),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in days",   # customize X axis label.
    break.time.by = 100,     # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )
  pdf(paste0(cvt, "_KaplanNR.pdf"), width = 5, height = 6, onefile = FALSE)
  print(p)
  dev.off()
}


################################################################################

# COXPH

survInput <- ann.sample[, c("TIME", "STATUS", "Therapy", "Imm.Subtype", "Macrophages.M2", "LA_TAM", "NonLA_TAM")] # WO Response

survInput$TIME <- as.numeric(survInput$TIME)
survInput$Response <- as.numeric(factor(survInput$Response))
#survInput$Therapy <- as.numeric(factor(survInput$Therapy))
#survInput$Imm.Subtype <- as.numeric(factor(survInput$Imm.Subtype))

write.xlsx(survInput, "/Users/......../survInput.xlsx")


library(tidyr)
library(dplyr)

Imm.Subtype <- survInput %>% mutate(value = 1)  %>% spread(Imm.Subtype, value,  fill = 0 ) 
Response <- survInput %>% mutate(value = 1) %>% spread(Therapy, value, fill = 0)

Res_col <- Response[, 7:9]

survInput2 <- cbind(Imm.Subtype, Res_col)
survInput <- survInput2

#Change colnames
names(survInput)[names(survInput) == 'CTLA4+PD1'] <- 'CTLA4.PD1'
names(survInput)[names(survInput) == 'PD1+KIR'] <- 'PD1.KIR'

setwd(coxDir)

# UNIVARIATE ANALYSIS

covariates <- c("Macrophages.M2","LA_TAM", "NonLA_TAM", "Subtype1", "Subtype2", "Subtype3", "Subtype4", "Subtype6", "CTLA4", "CTLA4.PD1", "PD1")
#covariates <- colnames(survInput)[-c(1:2)]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(TIME, STATUS)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = survInput)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
res$p.value <- as.numeric(res$p.value)
res <- res[order(res$p.value), ]

write.xlsx(res, "coxph_univariate.xlsx", rowNames = TRUE)

# MULTI-VARIATE ANALYSIS

# Restrict to significant features
mdrcox <- coxph(Surv(TIME, STATUS) ~ CTLA4.PD1 + Subtype4 + PD1 + LA_TAM + NonLA_TAM, 
                data=survInput)
sink("coxph_multivariate_summary_Response.txt")
summary(mdrcox)
sink()

p <- ggforest(mdrcox, data = survInput)
ggsave(plot = p, filename = "coxph_multivariate_forestplot_Response.pdf", width = 7, height = 7)

hzd <- cox.zph(mdrcox)
write.table(hzd$table, "coxph_multivariate_hazard_assumption_Response.txt",
            row.names = TRUE, sep = "\t", quote = FALSE)

# Stratify AGEBIN
mdrcox <- coxph(Surv(TIME, STATUS) ~ strata(AGEBIN),
                data=survInput)
sink("coxph_final_multivariate_summary.txt")
summary(mdrcox)
sink()

hzd <- cox.zph(mdrcox)
write.table(hzd$table, "coxph_final_multivariate_hazard_assumption.txt",
            row.names = TRUE, sep = "\t", quote = FALSE)


################################################################################
# RANDOM FOREST
r_fit <- ranger(Surv(TIME, STATUS) ~ APOE +
                  TRIM2 + IL10,
                data = ann.sample,
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE)

# Average the survival models
death_times <- r_fit$unique.death.times 
surv_prob <- data.frame(r_fit$survival)
avg_prob <- sapply(surv_prob,mean)

# Plot the survival models for each patient
plot(r_fit$unique.death.times,r_fit$survival[1,], 
     type = "l", 
     ylim = c(0,1),
     col = "red",
     xlab = "Days",
     ylab = "survival",
     main = "Patient Survival Curves")

#
cols <- colors()
for (n in sample(c(2:dim(ann.sample)[1]), 20)){
  lines(r_fit$unique.death.times, r_fit$survival[n,], type = "l", col = cols[n])
}
lines(death_times, avg_prob, lwd = 2)
legend(500, 0.7, legend = c('Average = black'))

vi <- data.frame(sort(round(r_fit$variable.importance, 4), decreasing = TRUE))
names(vi) <- "importance"
head(vi)
