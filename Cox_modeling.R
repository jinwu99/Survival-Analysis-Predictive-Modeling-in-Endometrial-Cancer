getwd()
rm(list=ls())
library("survival")
library("survminer")
library(coxphf)

df <- read.csv('stat_final_csv.csv')
df <- df[-c(110:121),]
df$MI.3group. <- df$MI.3group. + 1
df$MI.2group. <- df$MI.2group. + 1
df$Histology_3group <- df$Histology_specific
df[df$Histology_3group>2,'Histology_3group'] <- 3
df$Histology_2group <- df$Histology_specific
df[df$Histology_2group>1,'Histology_2group'] <- 2
df$MI.2group_v2. <- df$MI.3group.
df[df$MI.2group_v2.>1,'MI.2group_v2.'] <- 2

col_use <- c('new_OS','new_PFS','recurrence','new_status',
             'Age','BMI','ultrastaging_lymphnode_metastasis',
             'LVIbinary.3group.','MI.2group.',
             'CA.125','Histology_3group','CSI','Mass.size.binary.2cm.')
data <- df[col_use]
data$Age.2group.51. <- 0; data$Age.2group.61. <- 0
data[data$Age>=51,'Age.2group.51.'] <- 1
data[data$Age>=61,'Age.2group.61.'] <- 1

tmp <- colnames(data)
cat_vars <- tmp[-c(1:6,10)]
for (var in cat_vars) {
  data[[var]] <- as.factor(data[[var]])
}
data$new_OS <- as.numeric(data$new_OS)
data$new_PFS <- as.numeric(data$new_PFS)
data$recurrence <- as.numeric(data$recurrence)
data$new_status <- as.numeric(data$new_status)
data$CA.125 <- as.numeric(data$CA.125)
summary(data)

penalty <- 0.6

# variable selection
tmp <- colnames(data)
covariates <- tmp[-c(1:4,15)]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(new_PFS, recurrence)~', x)))
univ_formulas2 <- sapply(covariates,
                         function(x) as.formula(paste('Surv(new_OS, new_status)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxphf(x, data = data, penalty=penalty)})
univ_models2 <- lapply(univ_formulas2, function(x){coxphf(x, data = data, penalty=penalty)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         name <- names(x$coefficients)
                         HR <- round(exp(x$coefficients),3)
                         p_val <- round(x$prob,3)
                         HR.l <- round(x$ci.lower,3)
                         HR.u <- round(x$ci.upper,3)
                         res <- matrix(c(name,HR,HR.l,HR.u,p_val),nrow=length(HR),byrow=F,
                                       dimnames=list(rownames(x$coef),c('Variable','HR','HR.lower','HR.upper','p_value')))
                         return(res)
                       })
univ_results2 <- lapply(univ_models2,
                        function(x){ 
                          x <- summary(x)
                          name <- names(x$coefficients)
                          HR <- round(exp(x$coefficients),3)
                          p_val <- round(x$prob,3)
                          HR.l <- round(x$ci.lower,3)
                          HR.u <- round(x$ci.upper,3)
                          res <- matrix(c(name,HR,HR.l,HR.u,p_val),nrow=length(HR),byrow=F,
                                        dimnames=list(rownames(x$coef),c('Variable','HR','HR.lower','HR.upper','p_value')))
                          return(res)
                        })
res <- do.call(rbind,univ_results); res2 <- do.call(rbind,univ_results2)
res
res2

valid_cov <- covariates[-c(2,4,6,8:10)]
valid_cov <- c(valid_cov,'new_PFS','recurrence')
cov0 <- valid_cov
data0 <- data[cov0]
cox0 <- coxph(Surv(new_PFS, recurrence)~
                Age +
                ultrastaging_lymphnode_metastasis +
                MI.3group.
              , data=data0)
cox0 <- coxph(Surv(new_PFS, recurrence)~
                MI.3group. +
                ultrastaging_lymphnode_metastasis
              , data=data0)
cox0 <- coxph(Surv(new_PFS, recurrence)~
                MI.3group. +
                Age
              , data=data0)
cox0 <- coxph(Surv(new_PFS, recurrence)~
                MI.3group. +
                Histology_3group
              , data=data0)

names(data0)[2:4] <- c('Ultrastaging','M','H')
cox0 <- coxphf(Surv(new_PFS, recurrence)~., data=data0, penalty=penalty)

# proportional Hazard test
cox0 <- coxph(Surv(new_PFS, recurrence)~., data=data0)
test.ph <- cox.zph(cox0)
test.ph
ggcoxzph(test.ph)
# functional form of each covariate
cox.null <- coxph(Surv(new_PFS, recurrence)~1, data=data0)
rr <- resid(cox.null)
plot(data0$Age, rr); lines(lowess(data0$Age, rr))
# influential : nothing to consider
par(mfrow=c(2,2))
rr <- resid(cox0,type='dfbeta')
for(j in 1:4){
  plot(rr[,j],
       main=colnames(data0)[j],
       xlab='data index',
       ylab='residuals')
}

cox_list <- list(cox0)
cox_summary <- function(cox){
  x <- summary(cox)
  name <- names(x$coefficients)
  HR <- round(exp(x$coefficients),3)
  p_val <- round(x$prob,3)
  HR.l <- round(x$ci.lower,3)
  HR.u <- round(x$ci.upper,3)
  res <- matrix(c(name,HR,HR.l,HR.u,p_val),nrow=length(HR),byrow=F,
                dimnames=list(rownames(x$coef),c('Variable','HR','HR.lower','HR.upper','p_value')))
  return(res)
}

res_list <- list()
for(i in 1:1){
  cat('Cox model',i,' : \n')
  res_list[[i]] <- cox_summary(cox_list[[i]])
  cat('\n')
}
res_list
table(data$recurrence,data$FIGO2023binary)
table(data$recurrence,data$Cytology.2group)


valid_cov <- covariates[-c(2,4,6,8:10)]
valid_cov <- c(valid_cov,'new_OS','new_status')
cov0 <- valid_cov
data0 <- data[cov0]
cox0 <- coxphf(Surv(new_OS, new_status)~., data=data0, penalty=penalty)

# proportional Hazard test
names(data0)[2:4] <- c('Ultrastaging','M','H')
cox0 <- coxph(Surv(new_OS, new_status)~., data=data0)
test.ph <- cox.zph(cox0)
test.ph
ggcoxzph(test.ph)
# functional form of each covariate
cox.null <- coxph(Surv(new_PFS, recurrence)~1, data=data0)
rr <- resid(cox.null)
plot(data0$Age, rr); lines(lowess(data0$Age, rr))
# influential : nothing to consider
par(mfrow=c(2,2))
rr <- resid(cox0,type='dfbeta')
for(j in 1:4){
  plot(rr[,j],
       main=colnames(data0)[j],
       xlab='data index',
       ylab='residuals')
}

cox_list <- list(cox0)

cox_summary <- function(cox){
  x <- summary(cox)
  name <- names(x$coefficients)
  HR <- round(exp(x$coefficients),3)
  p_val <- round(x$prob,3)
  HR.l <- round(x$ci.lower,3)
  HR.u <- round(x$ci.upper,3)
  res <- matrix(c(name,HR,HR.l,HR.u,p_val),nrow=length(HR),byrow=F,
                dimnames=list(rownames(x$coef),c('Variable','HR','HR.lower','HR.upper','p_value')))
  return(res)
}

res_list <- list()
for(i in 1:1){
  cat('Cox model',i,' : \n')
  res_list[[i]] <- cox_summary(cox_list[[i]])
  cat('\n')
}
res_list
