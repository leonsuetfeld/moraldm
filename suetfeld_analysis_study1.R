# ###############################################
# ### load library / top statements #############
# ###############################################

library(brms)
library(tidyverse)

options(mc.cores = parallel::detectCores())
options(max.print = 2000)

# ###############################################
# ### load & pre-process data ###################
# ###############################################
d1 = read.csv("trialmatrix.csv")
source('git_moraldm/suetfeld_preprocessing.R')
d1.hum.effect = preprocessing(d1,study='1')


# ###############################################
# ### define & fit bayesian models ##############
# ###############################################

# info on priors:
# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-the-regression-coefficients-in-logistic-regression-non-sparse-case
# http://jakewestfall.org/misc/SorensenEtAl.pdf

# full model
model_full = choice_left ~
  0 + Intercept + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction +
  (0 + Intercept + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction | sn_idx)

priors <- c(set_prior("lkj(2)", class = "cor"),
            set_prior("normal(0,3)", class = "b"),
            set_prior("cauchy(0,1)", class = "sd", group = "sn_idx")) # alternatives: exponential(1) / cauchy(0,1)

mres_full = brm(model_full,
                family=bernoulli,
                data=d1.hum.effect,
                control = list(adapt_delta=0.95),sample_prior = 'yes',
                chains=4,iter=8000,cores=4,warmup=2000,
                prior = priors)
summary(mres_full, maxsum=2)
save("mres_full",file="stan_model_1_withPriors.RData")



model_full_noabstraction = choice_left ~
  1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality +
  (1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction | sn_idx)

priors <- c(set_prior("lkj(2)", class = "cor"),
            set_prior("normal(0,3)", class = "b"),
            set_prior("cauchy(0,1)", class = "sd", group = "sn_idx")) # alternatives: exponential(1) / cauchy(0,1)

mres_full_noabstraction = brm(model_full_noabstraction,
                              family=bernoulli,
                              data=d1.hum.effect,
                              control = list(adapt_delta=0.95),
                              chains=4,iter=8000,cores=4,warmup=2000,
                              prior = priors)

save("mres_full_noabstraction",file="stan_model_1_noabstraction.RData")
