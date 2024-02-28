#####################################################################################################################################################
#                                                                                                                                                   #
#   Replication code for Wilder, Rosalle and Bishop (2023) 'Eco-welfare states and just transitions: a multi-method analysis and research agenda'   #    
#                                                                                                                                                   #
#                                                            Updated 22 January, 2024                                                               #                                                                                                                                                                                                        
#                                                                                                                                                   #
#####################################################################################################################################################

# download and store replication data files to your working directory here:
getwd()

#####################
# load main dataset #
#####################

library(readxl)
data <- read_excel("./dataset.xlsx") 

###########################################################
# create the dependent variable 'welfare state robustness #
###########################################################

# standardize measures to be used in dimensionality reduction 
data$sickgen.scaled = scale(data$sickgen, center= TRUE, scale=TRUE)
data$uegen.scaled = scale(data$uegen, center= TRUE, scale=TRUE)
data$pengen.scaled = scale(data$pengen, center= TRUE, scale=TRUE)
data$almp_pmp.scaled = scale(data$almp_pmp, center= TRUE, scale=TRUE)
data$ud.scaled = scale(data$ud, center= TRUE, scale=TRUE)
data$unioncent.scaled = scale(data$unioncent, center= TRUE, scale=TRUE)
data$wcoord.scaled = scale(data$wcoord, center= TRUE, scale=TRUE)
data$wc_rights.scaled = scale(data$wc_rights, center= TRUE, scale=TRUE)
data$wc.scaled = scale(data$wc, center= TRUE, scale=TRUE)
data$wc_struct.scaled = scale(data$wc_struct, center= TRUE, scale=TRUE)
data$wc_negot.scaled = scale(data$wc_negot, center= TRUE, scale=TRUE)

# subset data to include just those variables used in factor analysis
new_data <- as.data.frame(subset(data, year < 2020, select = c(country, uegen.scaled, sickgen.scaled ,pengen.scaled, 
                                                almp_pmp.scaled, ud.scaled,unioncent.scaled, wcoord.scaled,
                                                wc_rights.scaled, wc.scaled, wc_struct.scaled, wc_negot.scaled)))

# factor analysis

# split the dataset to avoid overfitting 
set.seed(666)
N <- nrow(subset(new_data, select = -c(country)))
indices <- seq(1, N)
indices_efa <- sample(indices, floor((.5*N)))
indices_cfa <- indices[!(indices %in% indices_efa)]
efa_data <- subset(new_data, select = -c(country))[indices_efa, ]
cfa_data <- subset(new_data, select = -c(country))[indices_cfa, ]

# inspect eigen values to calculate dimensionality
library(psych)
efa_data_cor <- cor(efa_data, use = "pairwise.complete.obs") # calculate the correlation matrix
scree(efa_data_cor, factors = FALSE)  # look for components with values > 1                  
               
# run an exploratory factor analysis (efa) with 2 factors (as arguably indicated by scree plot)
efa_model <- fa(efa_data, nfactors = 2) 
efa_model$loadings
fa.diagram(efa_model)

# interpret varimax rotated loadings
efa_rotated <- fa(efa_data, nfactors = 2, rotate = "varimax")
efa_rotated$loadings

# goodness of fit
# ideally, likelihood chi-square would have a non-significant result, meaning observed and expected data are not significantly different, but a large N will return significant values 
# TLI (penalized non-normal fit index) penalizes more complex models for additional parameters: should be > 0.9, measures how well observed data match expected data
# RMSEA quantifies differences between observed and expected data: should be < 0.05 and no greater than 0.10
efa_model

# assess relative fit of 1, 2 and 3 factor models
fa(efa_data, nfactors = 1)
fa(efa_data, nfactors = 2)
fa(efa_data, nfactors = 3)

# although all manifest variables are positively correlated, we find that models including all manifest variables return poor goodness of fit statistics, hence the omission of manifest variables in the cfa below


# confirmatory factor analysis (cfa)
library(lavaan)
cfa_loadings <- cfa(model ='ws_robust =~ 
                            almp_pmp.scaled+
                            uegen.scaled + 
                            sickgen.scaled +
                            unioncent.scaled',       # there is a strong theoretical rationale for including union centrality (Scharpf 1987)
                        #    ud.scaled +             # ud as the lowest loading (and the variable has missing values)
                        #    wc.scaled +             # wc_struct has the strongest theoretical justification for works council variables, but has high mi values vis-a-vis union centrality and the generosity variables and so is omitted to avoid overfitting
                        #    wc_rights.scaled +      
                        #    wc_struct.scaled +
                        #    wcoord.scaled           # ditto for wcoord
                        #    other_factor =~         # it is difficult to discern what the second factor might represent (if anything); loadings are not terribly high
                        #    pengen.scaled +
                        #    wc_negot.scaled', 
                             data = cfa_data)  

summary(cfa_loadings, standardized = TRUE, fit.measures = TRUE) # std.all, lv of 0.3 is an established minimum threshold for loading
# CFI and TLI should be above 0.9
# RMSEA and SRMR should be below 0.1 (ideally below 0.05)

# examine modification indices to improve model fit. High mi between variables indicate correlated errors ( > 10 is a common threshold)
modificationIndices(cfa_loadings, sort. = TRUE)

# plot the results
library(semPlot)
semPaths(object = cfa_loadings,
         whatLabels = "std",
         edge.label.cex = 1,
         rotation = 2,
         edge.color = "black")

# confirm robustness to full data set (near identical loadings)
cfa_loadings_full <- cfa(model ='ws_robust =~ almp_pmp.scaled + uegen.scaled + sickgen.scaled + unioncent.scaled',  data = data)
summary(cfa_loadings_full, standardized = TRUE, fit.measures = TRUE) 
semPaths(object = cfa_loadings_full, whatLabels = "std", edge.label.cex = 1, rotation = 2, edge.color = "blue")

# confirm robustness to static data
loadings_1995 <- cfa(model ='ws_robust =~ almp_pmp.scaled + uegen.scaled + sickgen.scaled + unioncent.scaled',  data = subset(data, year == 1995))
summary(loadings_1995, standardized = TRUE, fit.measures = TRUE) 
semPaths(object = loadings_1995, whatLabels = "std", edge.label.cex = 1, rotation = 2, edge.color = "red3")

loadings_2005 <- cfa(model ='ws_robust =~ almp_pmp.scaled + uegen.scaled + sickgen.scaled + unioncent.scaled',  data = subset(data, year == 2005))
summary(loadings_2005, standardized = TRUE, fit.measures = TRUE) 
semPaths(object = loadings_2005, whatLabels = "std", edge.label.cex = 1, rotation = 2, edge.color = "green4")

loadings_2015 <- cfa(model ='ws_robust =~ almp_pmp.scaled + uegen.scaled + sickgen.scaled + unioncent.scaled',  data = subset(data, year == 2015))
summary(loadings_2015, standardized = TRUE, fit.measures = TRUE) 
semPaths(object = loadings_2015, whatLabels = "std", edge.label.cex = 1, rotation = 2, edge.color = "purple3")


# model fit is preserved using single year data, and is not biased longitudinally
# non-stationarity in time series do not appear to unduly bias factor loadings
# we therefore leave implementation using dynamic factor analysis to future research  


# create the ws_robust variable from factor loadings 
factor_loadings  <- predict(cfa(model ='ws_robust =~ almp_pmp.scaled + uegen.scaled + sickgen.scaled + unioncent.scaled',
                                data = subset(data, select = -country))) 

# append factor loadings object to new_data dataframe 
new_data <- cbind(new_data, factor_loadings)


# convert welfare state robustness to index between 0 and 1
new_data$ws_robust = (new_data$ws_robust-min(new_data$ws_robust))/(max(new_data$ws_robust)-min(new_data$ws_robust))



# join to to main dataframe (data)
library(dplyr)
data <- data %>% 
  left_join(new_data)



###########################################
# create GHG per capita variable, ghg_pop #
###########################################
data$ghg_pop <- (data$ghg*1000)/data$population # change to 'ghg2' for CO2 exluding LULUCF 

# set the color palette
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(8, "Set2"))(22)



###########################################
#    plot the bivariate relationship      #
###########################################

# generate a plot of ghg per capita (ghg_pop) against welfare state robustness (ws_robust) generated from confirmatory factor analysis (cfa)
library(ggplot2)


plot_cfa <- ggplot(aes(x = ws_robust, y = ghg_pop, color = id, label = label3), data = data)+  # do "label = label2" for industrial emissions
                    geom_smooth(method = "lm", se = FALSE, color = "light grey")+
                    geom_smooth(method = "lm", se = FALSE, size = 0.5, linewidth = 0.5,
                                aes(group = id))+
                    geom_point(size = 0.4)+
                    #  ylab(bquote("industrial "*CO[2]*" tons per capita"))+
                    ylab(bquote(CO[2]*" tons per capita"))+
                    xlab("welfare state robustness")+
                    #   scale_y_continuous(limits = c(0, 15),expand = c(0,0))+
                    scale_y_continuous(limits = c(0, 40), expand = c(0,0))+
                    theme_classic() +
                    theme(legend.position="none")+
                    stat_ellipse(geom = "polygon",type = "norm", level = 0.90,  alpha = 0.2,
                                 aes(fill = id,), linetype = 0)+
                    scale_color_manual(values = colors)+
                    scale_fill_manual(values = colors)+
                    geom_text(size = 3, color = "black")

#view the plot
plot_cfa


#assess Pierson correlation coefficient 
cor.test(data$ws_robust, data$ghg_pop)

ggsave('FIG4.png',plot_cfa, height = 4, width = 5, units = "in", dpi = 300)


###############################
#  hierachical mixed effects  #
###############################

# random effects anova (no predictors)
library(lme4)
baseline <- lmer(ghg_pop~1 + (1|id), data = subset(data, year < 2020)) # variables outside the parentheses = fixed effects, inside the parentheses = random effects (fixed slope as well as random slope)
summary(baseline) 
# id variance is the average variance of the country means from the fixed effect mean
# residual variance is the average within country variance from the country mean
# fixed effects estimate is average mean across all countries and years

# ICC (intra-class correlation) test tells us how much of the total variance is due to clustering 
library(performance)
icc(baseline) # 92% of variability is due to clustering, so a hierarchical model is required!

#plot the model to get a sense of what's going on
library(flexplot)  # run devtools::install_github("dustinfife/flexplot", ref="development") if necessary
visualize(baseline)

# random effects (random slopes, random intercepts)
rand.slope.int <- lmer(ghg_pop ~  ws_robust + (ws_robust|id), data = subset(data, year < 2020))
summary(rand.slope.int)
# id variance is variability about the slope AND intercept

# plot the relationship; black line is the fixed effect, colored lines are random effects
visualize(rand.slope.int, plot = "model")

# compare with baseline 
model.comparison(rand.slope.int, baseline)


# bivariate model with regulation instead of ws_robust
# convert environmental policy stringency to a score  between 0 and 1
data$stringency.index = (data$stringency-min(data$stringency, na.rm=T))/(max(data$stringency,na.rm=T)-min(data$stringency, na.rm=T))
rand.slope.int.reg <- lmer(ghg_pop ~  stringency.index + (stringency.index|id), data = subset(data, year < 2020))  # note stringency is a random effect and fixed effect in this model but not in the full model below
summary(rand.slope.int.reg)


# plot the relationship; black line is the fixed effect, colored lines are random effects
visualize(rand.slope.int.reg, plot = "model")

# compare with baseline 
model.comparison(rand.slope.int.reg, baseline)


# bivariate model with post materialist values instead of ws_robust
# convert post materialist values (pmv) score to an index between 0 and 1
data$pmv1w.index <- (data$pmv1w-min(data$pmv1w, na.rm=T))/(max(data$pmv1w, na.rm=T)-min(data$pmv1w, na.rm=T))
rand.slope.int.val <- lmer(ghg_pop ~  pmv1w.index + (pmv1w.index|id), data = subset(data, year < 2020)) 
summary(rand.slope.int.val)


# plot the relationship; black line is the fixed effect, colored lines are random effects
visualize(rand.slope.int.val, plot = "model")

# compare with baseline 
model.comparison(rand.slope.int.val, baseline)


# now add controls (EU and Kyoto and post-materialist values as fixed effects, others as random effects)
# but first get variables on to similar scales...
# convert political constraints (polcon, i.e., veto players) to a score  between 0 and 1
data$polconiii.index = (data$polconiii-min(data$polconiii, na.rm=T))/(max(data$polconiii,na.rm=T)-min(data$polconiii, na.rm=T))
# convert political constraints (polcon, i.e., veto players) to a score  between 0 and 1
data$GVC.index = (data$GVC-min(data$GVC, na.rm=T))/(max(data$GVC,na.rm=T)-min(data$GVC, na.rm=T))
# convert growth to a decimal
data$growth.dec <- data$growth/100
# convert pollution to deaths per 100 million 
data$pollution.hm <- data$pollution/100
# convert government R&D to proportion of GDP 
data$GBARD_gdp <- data$GBARD/(data$GDP/1000)

# for the sake of the simplifaction (which is required for convergence) assess which control variables can be fixed effects only 
# ...starting with post materialist values
fixed_slope_pmv <- lmer(ghg_pop ~ pmv1w.index + (1|id), data = subset(data, year < 2020)) 
random_slope_pmv <- lmer(ghg_pop ~ pmv1w.index + (pmv1w.index|id), data = subset(data, year < 2020))
compare.fits(ghg_pop ~ pmv1w.index|id, data = subset(data, year < 2020), fixed_slope_pmv, random_slope_pmv)
# models differ, so pmv should be modeled as both a random and fixed effect

# now growth
fixed_slope_growth <- lmer(ghg_pop ~ growth.dec + (1|id), data = subset(data, year < 2020))
random_slope_growth <- lmer(ghg_pop ~ growth.dec + (growth.dec|id), data = subset(data, year < 2020))
compare.fits(ghg_pop ~ growth.dec|id, data = subset(data, year < 2020), fixed_slope_growth, random_slope_growth)
# models are highly similar, so growth can just be a fixed effect 

# environmental policy stringency index
fixed_slope_stringent <- lmer(ghg_pop ~ stringency.index + (1|id), data = subset(data, year < 2020))
random_slope_stringent <- lmer(ghg_pop ~ stringency.index + (stringency.index|id), data = subset(data, year < 2020))
compare.fits(ghg_pop ~ stringency.index|id, data = subset(data, year < 2020), fixed_slope_stringent, random_slope_stringent)
# models are highly similar, so stringency can just be a fixed effect

# veto players
fixed_slope_veto <- lmer(ghg_pop ~ polconiii.index + (1|id), data = subset(data, year < 2020))
random_slope_veto <- lmer(ghg_pop ~ polconiii.index + (polconiii.index|id), data = subset(data, year < 2020))
compare.fits(ghg_pop ~ polconiii.index|id, data = subset(data, year < 2020), fixed_slope_veto, random_slope_veto)
# models differ, so veto players should be modeled as both a random and fixed effect

# GVC participation
fixed_slope_GVC <- lmer(ghg_pop ~ GVC.index + (1|id), data = subset(data, year < 2020))
random_slope_GVC <- lmer(ghg_pop ~ GVC.index + (GVC.index|id), data = subset(data, year < 2020))
compare.fits(ghg_pop ~ GVC.index|id, data = subset(data, year < 2020), fixed_slope_GVC, random_slope_GVC)
# models are fairly similar, so GVC index can be modeled as a fixed effect only


# government R&D
fixed_slope_RD <- lmer(ghg_pop ~ GBARD_gdp + (1|id), data = subset(data, year < 2020))
random_slope_RD <- lmer(ghg_pop ~ GBARD_gdp + (GBARD_gdp|id), data = subset(data, year < 2020))
compare.fits(ghg_pop ~ GBARD_gdp|id, data = subset(data, year < 2020), fixed_slope_RD, random_slope_RD)
# models are mostl similatly signed across most, but not all countries; however, obtaining non-singular fit requires this random effect be omitted

# random AND fixed effects for some variables (inside parentheses) but fixed effects only for others (outside parentheses) 
full.rand.slope.int <-  lmer(ghg_pop ~ ws_robust + kyoto +  pmv1w.index + eu + growth.dec + stringency.index + polconiii.index +  GVC.index + GBARD_gdp + (ws_robust + pmv1w.index + polconiii.index|id), data = subset(data, year < 2020)) 
summary(full.rand.slope.int) # in the paper the growth coefficient is reported as a percentage point change, not a 100 percentage point change, so divide the coefficient and standard error by 100

# compare with baseline 
model.comparison(full.rand.slope.int, baseline)

#obtain p-values by refitting models with lmerTest loaded 
library(lmerTest)
baseline <- lmer(ghg_pop~1 + (1|id), data = subset(data, year < 2020)) # variables outside the parentheses = fixed effects, inside the parentheses = random effects (fixed slope as well as random slope)
summary(baseline)

rand.slope.int <- lmer(ghg_pop ~  ws_robust + (ws_robust|id), data = subset(data, year < 2020))
summary(rand.slope.int)

rand.slope.int.reg <- lmer(ghg_pop ~  stringency.index + (stringency.index|id), data = subset(data, year < 2020))  # note stringency is a random effect and fixed effect in this model but not in the full model below
summary(rand.slope.int.reg)

rand.slope.int.val <- lmer(ghg_pop ~  pmv1w.index + (pmv1w.index|id), data = subset(data, year < 2020)) 
summary(rand.slope.int.val)

full.rand.slope.int <-  lmer(ghg_pop ~ ws_robust + kyoto +  pmv1w.index + eu + growth.dec + stringency.index + polconiii.index + GVC.index + GBARD_gdp + (ws_robust + pmv1w.index + polconiii.index|id), data = subset(data, year < 2020)) 
summary(full.rand.slope.int) # in the paper the growth coefficient is reported as a percentage point change, not a 100 percentage point change, so divide the coefficient and standard error by 100

# run detach("package:lme4", unload=TRUE) if necessary

# test for heteroskedasticity 

# plot residuals by country
sub_data <- subset(data, year < 2020)
sub_data$residuals <- residuals(full.rand.slope.int)
sub_data$fitted <- fitted(full.rand.slope.int)

# clouds indicate homoskedasticity, patterns indicate heteroskedasticity
ggplot(sub_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ id) +  
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  labs(x = "fitted values", y = "residuals")

# fit a robust lmer model for model 8
library(robustlmm) 
robust_full.rand.slope.int <- rlmer(ghg_pop ~ ws_robust + kyoto + pmv1w.index + eu + growth.dec +
                               stringency.index + polconiii.index +  GVC.index + GBARD_gdp +
                               (ws_robust + pmv1w.index + polconiii.index | id),
                               data = subset(data, year < 2020))

# compare with original model
summary(robust_full.rand.slope.int) # in the paper the growth coefficient is reported as a percentage point change, not a 100 percentage point change, so divide the coefficient and standard error by 100
summary(full.rand.slope.int)
# coefficients change slightly in the robust model, suggesting that heteroskedacity slightly biases the less conservative model 

# calculate approximate p-values for the robust model
summary_robust <- summary(robust_full.rand.slope.int)
estimates <- summary_robust$coefficients[, "Estimate"]
std_errors <- summary_robust$coefficients[, "Std. Error"]

# Wald statistic (z-value)
z_values <- estimates / std_errors

# p-values from the normal distribution
p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
p_values


# fit robust models for and baseline and (simplified) rand.slope.int models
robust_baseline <- rlmer(ghg_pop~1 + (1|id), data = subset(data, year < 2020))
summary(robust_baseline)

# obtain p-values
summary_robust.bl <- summary(robust_baseline)
estimates <- summary_robust.bl$coefficients[, "Estimate"]
std_errors <- summary_robust.bl$coefficients[, "Std. Error"]
# Wald statistic (z-value)
z_values <- estimates / std_errors
# p-values from the normal distribution
p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
p_values 

# model 5
robust_rsi<- rlmer(ghg_pop ~  ws_robust + (ws_robust|id), data = subset(data, year < 2020))
summary(robust_rsi)

# obtain p-values
summary_robust.rsi <- summary(robust_rsi)
estimates <- summary_robust.rsi$coefficients[, "Estimate"]
std_errors <- summary_robust.rsi$coefficients[, "Std. Error"]
# Wald statistic (z-value)
z_values <- estimates / std_errors
# p-values from the normal distribution
p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
p_values 

# model 6
robust_rsi.r<- rlmer(ghg_pop ~  stringency.index + (stringency.index|id), data = subset(data, year < 2020))
summary(robust_rsi.r)

# obtain p-values
summary_robust.rsi.r <- summary(robust_rsi.r)
estimates <- summary_robust.rsi.r$coefficients[, "Estimate"]
std_errors <- summary_robust.rsi.r$coefficients[, "Std. Error"]
# Wald statistic (z-value)
z_values <- estimates / std_errors
# p-values from the normal distribution
p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
p_values 

# model 7
robust_rsi.v <- rlmer(ghg_pop ~  pmv1w.index + (pmv1w.index|id), data = subset(data, year < 2020)) 
summary(robust_rsi.v)

summary_robust_rsi.v <- summary(robust_rsi.v)
estimates <- summary_robust_rsi.v$coefficients[, "Estimate"]
std_errors <- summary_robust_rsi.v$coefficients[, "Std. Error"]
# Wald statistic (z-value)
z_values <- estimates / std_errors
# p-values from the normal distribution
p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
p_values 


###################################
#          pooled model           #
###################################

mod <- ghg_pop ~ ws_robust + kyoto + eu + pmv1w.index + growth.dec + stringency.index + polconiii.index + GVC.index + GBARD_gdp 


# segment data for a pooled model that avoids autocorrelation by randomly drawing N=50 samples and finding one (or more) without correlated errors
pooldata <- data %>%
        select(country, year, ghg_pop, ws_robust, kyoto, eu, pmv1w.index, growth.dec, stringency.index, polconiii.index, GVC.index, GBARD_gdp) %>% 
        na.omit


set.seed(250) #set.seed(28)

randomsample <- pooldata[sample(nrow(pooldata), size = 50), ]
library(plm)
m.pool <- plm(mod, data = randomsample, model = "pooling", index = c("country","year"))
summary(m.pool) # in the paper the growth coefficient is reported as a percentage point change, not a 100 percentage point change, so divide the coefficient and standard error by 100

# lagrange multiplier test for correlated errors (residuals)
plmtest(m.pool)
# p < 0.05 = autocorrelation (errors are correlated), which biases estimators.


# to err on the conservative side, calculate HAC (robust) standard errors 
library(lmtest)
coeftest(m.pool, vcov=function(x) vcovHC(x, method="arellano",type="HC1"))
coeftest(m.pool, vcov=function(x) vcovBK(m.pool,cluster="time")) # Beck and Katz

# obtain the periods t
length(unique(randomsample$year))

# test for generality of estimates by creating a sampling distribution
# take 10,000 samples of N = 50
f <- function () {
  fit <- plm(mod, data = pooldata, model = "pooling", index = c("country","year"), subset = sample(nrow(pooldata), 50))
  coef(fit) 
}

set.seed(11); pooled_estimates <- t(replicate(10000, f()))
head(pooled_estimates) 

# return the distributions of the coefficients and see where ours lands; presumably, those with higher absolute values are biased whereas those with lower absolute values are not
pooled_estimates <- data.frame(pooled_estimates)
hist(pooled_estimates$X.Intercept)
mean(pooled_estimates$X.Intercept)
var(pooled_estimates$X.Intercept)

hist(pooled_estimates$ws_robust)
mean(pooled_estimates$ws_robust)
var(pooled_estimates$ws_robust)


################################
#  between effects estimation  #
################################

m.be.mod <- plm(mod, data=data, model="between", index = c("country","year")) #index refers to the variables representing the time series cross sectional index, so the time and place 
summary(m.be.mod) # in the paper the growth coefficient is reported as a percentage point change, not a 100 percentage point change, so divide the coefficient and standard error by 100


b.be.mod <- plm(ghg_pop~ws_robust, data=data, model="between", index = c("country","year")) #index refers to the variables representing the time series cross sectional index, so the time and place 
summary(b.be.mod)


r.be.mod <- plm(ghg_pop~stringency.index, data= subset(data, year < 2020), model="between", index = c("country","year")) #index refers to the variables representing the time series cross sectional index, so the time and place 
summary(r.be.mod)

v.be.mod <- plm(ghg_pop~pmv1w.index, data= subset(data, year < 2020), model="between", index = c("country","year")) #index refers to the variables representing the time series cross sectional index, so the time and place 
summary(v.be.mod)


##############################################
#  robustness checks: consumption based CO2  #
##############################################

plot_cons <- ggplot(aes(x = ws_robust, y = ghg_cons, color = id, label = label2), data = data)+  # do "label = label2" for industrial emissions
              geom_smooth(method = "lm", se = FALSE, color = "light grey")+
              geom_smooth(method = "lm", se = FALSE, size = 0.5, linewidth = 0.5,
                          aes(group = id))+
              geom_point(size = 0.4)+
              ylab(bquote("consumption-based "*CO[2]*" tons per capita"))+
              xlab("welfare state robustness")+
              #   scale_y_continuous(limits = c(0, 15),expand = c(0,0))+
              scale_y_continuous(limits = c(0, 40), expand = c(0,0))+
              theme_classic() +
              theme(legend.position="none")+
              stat_ellipse(geom = "polygon",type = "norm", level = 0.90,  alpha = 0.2,
                           aes(fill = id,), linetype = 0)+
              scale_color_manual(values = colors)+
              scale_fill_manual(values = colors)+
              geom_text(size = 3, color = "black")

#view the plot
plot_cons

#assess Pierson correlation coefficient 
cor.test(data$ws_robust, data$ghg_cons)


# random effects 
baseline_cons <- lmer(ghg_cons~1 + (1|id), data = subset(data, year < 2020)) 
summary(baseline_cons) 

# ICC (intra-class correlation) test tells us how much of the total variance is due to clustering 
icc(baseline_cons) # 83% of variability is due to clustering, so a hierarchical model is required!

# random effects (random slopes, random intercepts)
rand.slope.int_cons <- lmer(ghg_cons ~  ws_robust + (ws_robust|id), data = subset(data, year < 2020))
summary(rand.slope.int_cons)

# compare with baseline (you may need to detach packages that are marking functions, e.g., lmerTest for model.comparison to work properly)
model.comparison(rand.slope.int_cons, baseline_cons)

# plot the relationship; black line is the fixed effect, colored lines are random effects
visualize(rand.slope.int_cons, plot = "model")

# stringency index instead 
r.rand.slope.int_cons <- lmer(ghg_cons ~  stringency.index + (stringency.index|id), data = subset(data, year < 2020))
summary(r.rand.slope.int_cons)

visualize(r.rand.slope.int_cons, plot = "model")
model.comparison(r.rand.slope.int_cons, baseline_cons)


# post-materialist values instead 
v.rand.slope.int_cons <- lmer(ghg_cons ~  pmv1w.index+ (pmv1w.index|id), data = subset(data, year < 2020))
summary(v.rand.slope.int_cons)

visualize(v.rand.slope.int_cons, plot = "model")
model.comparison(v.rand.slope.int_cons, baseline_cons)



# add controls 
# random AND fixed effects for some variables (inside parentheses) but fixed effects only for others (outside parentheses) 
full.rand.slope.int_cons <-  lmer(ghg_cons ~ ws_robust + kyoto +  pmv1w.index + eu + growth.dec + stringency.index + polconiii.index + GVC.index + GBARD_gdp + (ws_robust + pmv1w.index + polconiii.index|id), data = subset(data, year < 2020)) 
summary(full.rand.slope.int_cons) # in the paper the growth coefficient is reported as a percentage point change, not a 100 percentage point change, so divide the coefficient and standard error by 100

model.comparison(full.rand.slope.int_cons, baseline_cons)


# pooled model           
mod_cons <- ghg_cons ~ ws_robust + kyoto + eu + pmv1w.index + growth.dec + stringency.index + polconiii.index  + GVC.index + GBARD_gdp 


# segment data for a pooled model that avoids autocorrelation by randomly drawing N=50 samples and finding one (or more) without correlated errors
pooldata_cons <- data %>%
  select(country, year, ghg_cons, ws_robust, kyoto, eu, pmv1w.index, growth.dec, stringency.index, polconiii.index,  GVC.index, GBARD_gdp) %>% 
  na.omit



set.seed(21) # set.seed(37)
 

randomsample <- pooldata_cons[sample(nrow(pooldata_cons), size = 50), ]

m.pool_cons <- plm(mod_cons, data = randomsample, model = "pooling", index = c("country","year"))
summary(m.pool_cons) # pooled model

# lagrange multiplier test for correlated errors (residuals)
plmtest(m.pool_cons)
# p < 0.05 = autocorrelation (errors are correlated), which biases estimators.

# to err on the conservative side, calculate HAC (robust) standard errors 
coeftest(m.pool_cons, vcov=function(x) vcovHC(x, method="arellano",type="HC1"))
coeftest(m.pool_cons, vcov=function(x) vcovBK(m.pool_cons,cluster="time")) # Beck and Katz

# obtain the periods t
length(unique(randomsample$year))

# test for generality of estimates
# take 10,000 samples of N = 50
f <- function () {
  fit <- plm(mod_cons, data = pooldata_cons, model = "pooling", index = c("country","year"), subset = sample(nrow(pooldata_cons), 50))
  coef(fit) 
}

set.seed(11); pooled_estimates <- t(replicate(10000, f()))
head(pooled_estimates) 

# return the distributions of the coefficients and see where ours lands; presumably, those with higher absolute values are biased whereas those with lower absolute values are not
pooled_estimates <- data.frame(pooled_estimates)
hist(pooled_estimates$X.Intercept.)
mean(pooled_estimates$X.Intercept)
var(pooled_estimates$X.Intercept)

hist(pooled_estimates$ws_robust)
mean(pooled_estimates$ws_robust)
var(pooled_estimates$ws_robust)

#  between effects estimation  
m.be.mod_cons <- plm(mod_cons, data=data, model="between", index = c("country","year")) 
summary(m.be.mod_cons)


b.be.mod_cons <- plm(ghg_cons~ws_robust, data=data, model="between", index = c("country","year")) #index refers to the variables representing the time series cross sectional index, so the time and place 
summary(b.be.mod_cons)

r.be.mod_cons <- plm(ghg_cons~stringency.index, data=subset(data, year < 2020), model="between", index = c("country","year")) #index refers to the variables representing the time series cross sectional index, so the time and place 
summary(r.be.mod_cons)

v.be.mod_cons <- plm(ghg_cons~pmv1w.index, data=subset(data, year < 2020), model="between", index = c("country","year")) #index refers to the variables representing the time series cross sectional index, so the time and place 
summary(v.be.mod_cons)


#######################################################################################
#  relationship between welfare state robustness and environmental policy stringency  #
#######################################################################################

plot_stringency <- ggplot(aes(x = ws_robust, y =  stringency, color = id, label = label), data = data)+  # do "label = label2" for industrial emissions
                  geom_smooth(method = "lm", se = FALSE, color = "light grey")+
                  geom_smooth(method = "lm", se = FALSE, size = 0.5, linewidth = 0.5,
                              aes(group = id))+
                  geom_point(size = 0.4)+
                  ylab(bquote("environmental policy stringency"))+
                  xlab("welfare state robustness")+
                  scale_y_continuous(limits = c(-0.2, 5), expand = c(0,0))+
                  theme_classic() +
                  theme(legend.position="none")+
                  stat_ellipse(geom = "polygon",type = "norm", level = 0.90,  alpha = 0.2,
                               aes(fill = id,), linetype = 0)+
                  scale_color_manual(values = colors)+
                  scale_fill_manual(values = colors)+
                  geom_text(size = 3, color = "black")

#view the plot
plot_stringency

#assess Pierson correlation coefficient 
cor.test(data$ws_robust, data$stringency.index)


range(data$stringency, na.rm=T)

# random effects 
baseline_strin <- lmer(stringency~1 + (1|id), data = subset(data, year < 2020)) 
summary(baseline_strin) 

# ICC (intra-class correlation) test tells us how much of the total variance is due to clustering 
icc(baseline_strin) # 25% of variability is due to clustering, so a hierarchical model is required!

# random effects (random slopes, random intercepts)
rand.slope.int_strin <- lmer(stringency ~  ws_robust + (ws_robust|id), data = subset(data, year < 2020))
summary(rand.slope.int_strin)

# compare with baseline 
model.comparison(rand.slope.int_strin, baseline_strin)


# plot the relationship; black line is the fixed effect, colored lines are random effects
visualize(rand.slope.int_strin, plot = "model")

# add controls 
# random AND fixed effects for some variables (inside parentheses) but fixed effects only for others (outside parentheses) 
full.rand.slope.int_strin <-  lmer(stringency ~ ws_robust + kyoto +  pmv1w.index + eu + growth.dec + polconiii.index + GVC.index +GBARD_gdp + (ws_robust + pmv1w.index + polconiii.index|id), data = subset(data, year < 2020)) 
summary(full.rand.slope.int_strin) # in the paper the growth coefficient is reported as a percentage point change, not a 100 percentage point change, so divide the coefficient and standard error by 100

model.comparison(full.rand.slope.int_strin, baseline_strin)


# pooled model           
mod_strin <- stringency ~ ws_robust + kyoto + eu + pmv1w.index + growth.dec +  polconiii.index +  GVC.index + GBARD_gdp 


# segment data for a pooled model that avoids autocorrelation by randomly drawing N=50 samples and finding one (or more) without correlated errors
pooldata_strin <- data %>%
                    select(country, year, stringency, ws_robust, kyoto, eu, pmv1w.index, growth.dec, polconiii.index, GVC.index, GBARD_gdp) %>% 
                    na.omit


set.seed(250)

randomsample <- pooldata_strin[sample(nrow(pooldata_strin), size = 50), ]

m.pool_strin <- plm(mod_strin, data = randomsample, model = "pooling", index = c("country","year"))
summary(m.pool_strin) # pooled model


# lagrange multiplier test for correlated errors (residuals)
plmtest(m.pool_strin)
# p < 0.05 = autocorrelation (errors are correlated), which biases estimators.

# to err on the conservative side, calculate HAC (robust) standard errors 
coeftest(m.pool_strin, vcov=function(x) vcovHC(x, method="arellano",type="HC1"))
coeftest(m.pool_strin, vcov=function(x) vcovBK(m.pool,cluster="time")) # Beck and Katz

# obtain the periods t
length(unique(randomsample$year))

# test for generality of estimates
# take 10,000 samples of N = 50
f <- function () {
  fit <- plm(mod_strin, data = pooldata_strin, model = "pooling", index = c("country","year"), subset = sample(nrow(pooldata_strin), 50))
  coef(fit) 
}

set.seed(11); pooled_estimates <- t(replicate(10000, f()))
head(pooled_estimates) 

# return the distributions of the coefficients and see where ours lands; presumably, those with higher absolute values are biased whereas those with lower absolute values are not
pooled_estimates <- data.frame(pooled_estimates)
hist(pooled_estimates$X.Intercept)
mean(pooled_estimates$X.Intercept)
var(pooled_estimates$X.Intercept)

hist(pooled_estimates$ws_robust)
mean(pooled_estimates$ws_robust)
var(pooled_estimates$ws_robust)


#  between effects estimation  
m.be.mod_strin<- plm(mod_strin, data=data, model="between", index = c("country","year")) 
summary(m.be.mod_strin)


b.be.mod_strin <- plm(stringency~ws_robust, data=data, model="between", index = c("country","year")) #index refers to the variables representing the time series cross sectional index, so the time and place 
summary(b.be.mod_strin)



############################
#  descriptive statistics  #
############################


# labour/employment statistics 
UKdata <- read_excel("./employment statistics.xlsx") 

UK_plotL1 <- ggplot()+ 
             theme_bw() +  
          #  ggtitle("United Kingdom")+
              theme(axis.title.x=element_blank(), 
                  panel.grid = element_blank(),
                  axis.title.y = element_text(size =9),
                  legend.position = "none",
                  plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle =90, vjust = 0.5))+
              geom_line(data = subset (data, country == "United Kingdom" & year > 1995), aes(x=year, y=(emp_industry*1000)/population), size = 0.25) +
              geom_point(data = subset (data, country == "United Kingdom" & year > 1995), aes(x=year, y=(emp_industry*1000)/population), size = 2, shape = 4) +
              geom_line(data = subset (data, country == "United Kingdom" & year > 1995), aes(x=year, y=((emp_services_J+emp_services_K+emp_services_M)*1000)/population), size = 0.25) +
              geom_point(data = subset (data, country == "United Kingdom" & year > 1995), aes(x=year, y=((emp_services_J+emp_services_K+emp_services_M)*1000)/population), size = 2, shape = 21, fill ="white") +
              scale_x_continuous(limits= c(1995,2020),  breaks = seq(1995, 2020, 5), expand = c(0,0))+
              ylab("employment (capita)") +
              scale_y_continuous(limits = c(0, 0.15),breaks = c(0, 0.05, 0.10, 0.15), labels = c(0, 0.05, 0.10, 0.15), expand = c(0,0))+
          #    geom_point(aes(x=2010, y=0.04), shape = 4, size =2, color = "red3")  +  # create legend items
              annotate("text", x=1999, y=0.125, label= "industry", size = 4) +
         #     geom_point(aes(x=2010, y=0.03), shape = 19, size =1.5, color = "lightblue4")  +  
              annotate("text", x=2007, y=0.050, label= "knowledge-based services", size = 4)

UK_plotL1


UK_plotL2 <- ggplot(aes(x=year, y=employment/1000,shape = industry), data = subset(UKdata, industry == "steel" | industry == "chemicals" | industry == "coal" & year > 1985))+
                  geom_point(size = 1.5)+
                  scale_shape_discrete(breaks=c("chemicals","steel","coal"))+
                  geom_line(size = 0.25)+
                  theme_bw() +  
                  theme(axis.title.x=element_blank(), 
                        legend.spacing.x  = unit(-0.1, "cm"),
                        legend.title = element_blank(),
                        legend.text = element_text(size =12),
                        legend.position = c(.74,.8),
                        panel.grid.major.y = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.title.y = element_text(size =9),
                        plot.title = element_text(hjust = 0.5),
                        axis.text.x = element_text(angle =90, vjust = 0.5))+
                  ylab("employment (000)") +
                  scale_x_continuous(limits= c(1990,2020),  breaks = seq(1990, 2020, 5), expand = c(0,0))+
                  scale_y_continuous(limits= c(0, 120), expand = c(0,0))

UK_plotL2


UK_plotL3 <- ggplot(aes(x=year, y=employment/1000,shape = industry), data = subset(UKdata, industry == "refining" | industry == "cement"  & year > 1985))+
                geom_point(size = 1.2, fill = "black")+
                scale_shape_manual(breaks=c("refining","cement"), values=c(23,25))+
                geom_line(size = 0.25)+
                theme_bw() +  
                theme(axis.title.x=element_blank(), 
                      legend.spacing.x  = unit(-0.1, "cm"),
                      legend.title = element_blank(),
                      legend.text = element_text(size =12),
                      legend.position = c(.78,.8),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.title.y = element_text(size =9),
                      plot.title = element_text(hjust = 0.5),
                      axis.text.x = element_text(angle =90, vjust = 0.5))+
                ylab("employment (000)") +
                scale_x_continuous(limits= c(1990,2020),  breaks = seq(1990, 2020, 5), expand = c(0,0))+
                scale_y_continuous(limits= c(0, 25), expand = c(0,0))

UK_plotL3


library('ggpubr')
arrangement <- ggarrange(UK_plotL2, UK_plotL3, UK_plotL1, ncol = 3, nrow = 3)


ggsave('FIG6.png',
       arrangement, height = 9, width = 8.5, units = "in", dpi = 300)


# emissions by source 
denmark_data <- subset(data, country == "Denmark")
uk_data <- subset(data, country == "United Kingdom")

denmark_data$energy_per_capita <- denmark_data$ghg_energy * 1000 / denmark_data$population
denmark_data$industry_per_capita <- (denmark_data$ghg_industry * 1000 + denmark_data$ghg_manufacturing_construction * 1000) / denmark_data$population
uk_data$energy_per_capita <- uk_data$ghg_energy * 1000 / uk_data$population
uk_data$industry_per_capita <- (uk_data$ghg_industry * 1000 + uk_data$ghg_manufacturing_construction * 1000) / uk_data$population


combined_data <- rbind(
  data.frame(country = "Denmark", year = denmark_data$year, 
             source = "energy", value = denmark_data$energy_per_capita),
  data.frame(country = "Denmark", year = denmark_data$year, 
             source = "industry", value = denmark_data$industry_per_capita),
  data.frame(country = "United Kingdom", year = uk_data$year, 
             source = "energy", value = uk_data$energy_per_capita),
  data.frame(country = "United Kingdom", year = uk_data$year, 
             source = "industry", value = uk_data$industry_per_capita)
)


shape_energy <- 17 
shape_industry <- 16 


co2_plot <- ggplot(combined_data, aes(x = year, y = value)) +
              geom_line(aes(group = source), size = 0.25) +
              geom_point(aes(shape = source), size = 1.2) + 
              scale_shape_manual(values=c(energy = shape_energy, industry = shape_industry)) +
              theme_bw() +  
              theme(axis.title.x=element_blank(), 
                    legend.spacing.x  = unit(-0.1, "cm"),
                    legend.title.align = 0.5,
                    legend.title = element_blank(),
                    legend.text = element_text(size = 11),
                    legend.position = c(.88,.8),
                    panel.grid.major.y = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.title.y = element_text(size =9),
                    plot.title = element_text(hjust = 0.5),
                    axis.text.x = element_text(angle =90, vjust = 0.5),
                    strip.background = element_blank(),
                    strip.text = element_text(size = 12, color = "black"))+
              ylab(bquote(""*CO[2]*" tons per capita")) +
              scale_x_continuous(breaks = seq(1990, 2020, by = 5)) +
              scale_y_continuous(limits= c(0, 9), expand = c(0,0)) +
              facet_wrap(~ country, ncol = 2) +
              guides(linetype = guide_legend(override.aes = list(size = 2)))

co2_plot


ggsave('FIG5.png', co2_plot, height = 3, width = 6, units = "in", dpi = 300)


# multi panel plots 
AUS_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Australia"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Australia"), aes(x=year, y=(ws_robust)/.025), size =0.2)+
  geom_point(data = subset (data, country == "Australia"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Australia")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1995, y=37, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=9, label= "WS robustness", size = 3.5) 
AUS_plot1




AUT_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Austria"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Austria"), aes(x=year, y=(ws_robust)/.025), size =0.2)+
  geom_point(data = subset (data, country == "Austria"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Austria")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=12, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=32, label= "WS robustness", size = 3.5) 
AUT_plot1



BEL_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Belgium"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Belgium"), aes(x=year, y=(ws_robust)/.025), size =0.2)+
  geom_point(data = subset (data, country == "Belgium"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Belgium")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=19, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=33, label= "WS robustness", size = 3.5) 
BEL_plot1



CAN_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Canada"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Canada"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Canada"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Canada")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=23, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=6, label= "WS robustness", size = 3.5)
CAN_plot1


DEN_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Denmark"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Denmark"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Denmark"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Denmark")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=20, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=30, label= "WS robustness", size = 3.5)

DEN_plot1


FIN_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Finland"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Finland"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Finland"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Finland")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=14, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=30, label= "WS robustness", size = 3.5)
FIN_plot1



FRA_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "France"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "France"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "France"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("France")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=13, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=31, label= "WS robustness", size = 3.5)
FRA_plot1



GER_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Germany"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Germany"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Germany"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Germany")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=19, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=32, label= "WS robustness", size = 3.5)
GER_plot1



GRE_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Greece"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Greece"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Greece"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Greece")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=5, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=17, label= "WS robustness", size = 3.5)
GRE_plot1



IRE_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Ireland"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Ireland"), aes(x=year, y=(ws_robust)/.025), size =0.2)+
  geom_point(data = subset (data, country == "Ireland"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Ireland")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=22, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=28, label= "WS robustness", size = 3.5)
IRE_plot1


ITA_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Italy"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Italy"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Italy"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Italy")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=12, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=24, label= "WS robustness", size = 3.5)
ITA_plot1



JPN_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Japan"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Japan"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Japan"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Japan")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=7, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=22, label= "WS robustness", size = 3.5)
JPN_plot1


NET_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Netherlands"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Netherlands"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Netherlands"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Netherlands")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=19, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=35, label= "WS robustness", size = 3.5)

NET_plot1



NZ_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "New Zealand"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "New Zealand"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "New Zealand"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("New Zealand")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=2005, y=17, label= "GHG", size = 3.5) +
  annotate("text", x=1997, y=7, label= "WS robustness", size = 3.5)


NZ_plot1


NOR_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Norway"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Norway"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Norway"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Norway")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=13, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=31, label= "WS robustness", size = 3.5)

NOR_plot1


POR_plot1 <- ggplot() + 
  geom_line(data = subset (data, country == "Portugal"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Portugal"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Portugal"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Portugal")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=9, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=29, label= "WS robustness", size = 3.5)


POR_plot1



SPN_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Spain"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Spain"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Spain"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Spain")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=10, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=26, label= "WS robustness", size = 3.5)


SPN_plot1


SWE_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Sweden"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Sweden"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Sweden"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Sweden")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=7, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=30, label= "WS robustness", size = 3.5)


SWE_plot1


SWI_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "Switzerland"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "Switzerland"), aes(x=year, y=ws_robust/.025), size =0.2)+
  geom_point(data = subset (data, country == "Switzerland"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("Switzerland")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=11, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=31, label= "WS robustness", size = 3.5)

SWI_plot1


UK_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "United Kingdom"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "United Kingdom"), aes(x=year, y=(ws_robust)/.025), size =0.2)+
  geom_point(data = subset (data, country == "United Kingdom"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("United Kingdom")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=17, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=21, label= "WS robustness", size = 3.5)

UK_plot1


US_plot1 <- ggplot() +
  geom_line(data = subset (data, country == "United States"), aes(x=year, y=ghg_pop), size = 1) +
  geom_line(data = subset (data, country == "United States"), aes(x=year, y=(ws_robust)/.025), size =0.2)+
  geom_point(data = subset (data, country == "United States"), aes(x=year, y=(ws_robust)/.025), size =1, shape = 4)+
  theme_bw() +  
  scale_color_manual(values = c("0" = "black",
                                "1" ="gray")) + 
  ggtitle("United States")+
  theme(axis.title.x=element_blank(), 
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle =90, vjust = 0.5))+
  scale_x_continuous(breaks= c(1990, 1995, 2000, 2005, 2010, 2015, 2020))+
  ylab("CO2 t/pop")+
  scale_y_continuous(limits = c(0, 40),
                     sec.axis = sec_axis(~.*.025))+
  annotate("text", x=1992, y=26, label= "GHG", size = 3.5) +
  annotate("text", x=2013, y=3, label= "WS robustness", size = 3.5)

US_plot1



FIG3a <- ggarrange(US_plot1, CAN_plot1, AUS_plot1, IRE_plot1, UK_plot1,  NZ_plot1,
                     JPN_plot1, GRE_plot1, ITA_plot1, SPN_plot1, POR_plot1, FRA_plot1,  
                     ncol = 3, nrow = 4)
FIG3a


FIG3b <- ggarrange(GER_plot1, AUT_plot1, SWI_plot1, BEL_plot1,  DEN_plot1, FIN_plot1, 
                     NET_plot1,  SWE_plot1, NOR_plot1,
                     ncol = 3, nrow = 4)
FIG3b


ggsave('FIG3a.png',
       FIG3a, height = 9, width = 8.5, units = "in", dpi = 300)

ggsave('FIG3b.png',
       FIG3b, height = 9, width = 8.5, units = "in", dpi = 300)


