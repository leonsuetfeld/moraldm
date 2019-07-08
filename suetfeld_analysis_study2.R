# ###############################################
# ### load library / top statements #############
# ###############################################

library(sjPlot)
l# ###############################################
# ### load library / top statements #############
# ###############################################

library(brms)
library(tidyverse)

options(mc.cores = parallel::detectCores())
options(max.print = 2000)

# ###############################################
# ### load & pre-process data ###################
# ###############################################
d1 = read.csv("trialmatrix_2.csv")
source('suetfeld_preprocessing.R')
d1.hum.effect = preprocessing(d1,study='2')


# ###############################################
# ### define & fit bayesian models ##############
# ###############################################

# info on priors:
# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-the-regression-coefficients-in-logistic-regression-non-sparse-case
# http://jakewestfall.org/misc/SorensenEtAl.pdf


# full age sex sds vg model (some between others within)
model_agesexsdsvg = choice_left ~ 0+
  (0+Intercept + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction + 
  (0+Intercept + visonset_left + sex_diff + young_diff + elderly_diff)*sex + 
  (0+Intercept + visonset_left + sex_diff + young_diff + elderly_diff)*age_center + 
  (0+Intercept + visonset_left + sex_diff + young_diff + elderly_diff)*sds17_center + 
  (0+Intercept + visonset_left + sex_diff + young_diff + elderly_diff)*vg_center + 
  ((0+Intercept + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction | sn_idx)

# define priors
priors_agesexsdsvg <- c(set_prior("lkj(2)", class = "cor"),
                        set_prior("normal(0,3)", class = "b"),
                        set_prior("cauchy(0,1)", class = "sd", group = "sn_idx"))
mres_agesexsdsvg = brm(model_agesexsdsvg,
                       family=bernoulli,
                       data=d1.hum.effect,
                       control = list(adapt_delta=0.9),sample_prior = 'yes',
                       chains=4,iter=8000,cores=4,warmup=2000,
                       prior = priors_agesexsdsvg)


summary(mres_agesexsdsvg, maxsum=2)
save("mres_agesexsdsvg",file="stan_model_2_withPriors.RData")


# ###############################################
# ### Model Comparison to check effect of Abstraction
# ###############################################

# Model with abstraction fixed at abstraction = 0
model_agesexsdsvg_noabstraction = choice_left ~
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed + 
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*sex + 
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*age_center + 
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*sds17_center + 
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*vg_center + 
  ((1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction | sn_idx)

priors_agesexsdsvg <- c(set_prior("lkj(2)", class = "cor"),
                        set_prior("normal(0,3)", class = "b"),
                        set_prior("cauchy(0,1)", class = "sd", group = "sn_idx"))

mres_agesexsdsvg_noabstraction = brm(model_agesexsdsvg_noabstraction,
                                     family=bernoulli,
                                     data=d1.hum.effect,
                                     control = list(adapt_delta=0.9),
                                     chains=4,iter=8000,cores=4,warmup=2000,
                                     prior = priors_agesexsdsvg)


loo::compare(c(loo::loo(mres_agesexsdsvg),loo::loo(mres_agesexsdsvg_noabstraction)))