# ###############################################
# ### load library / top statements #############
# ###############################################

library(ggplot2)
library(brms)
library(bayesplot)
library(dplyr)
options(mc.cores = parallel::detectCores())
options(max.print = 2000)
#rstan_options(auto_write = TRUE)

# ###############################################
# ### load & pre-process data ###################
# ###############################################

setwd("/DATEN/PhD/research/VRdriveBA/_data_preStudy/data_csv/")
d1 = read.csv("trialmatrix.csv")

# obstacle.left / obstacle.right coding scheme:
# 1 = girl, 2 = boy, 3 = woman, 4 = man, 5 = oldwoman, 6 = oldman, 7 = deer, 8 = goat, 9 = boar, 10 = nothing
# condition coding scheme:
# 1 = desktop text, 2 = desktop image, 3 = vr text, 4 = vr image

# add dummy variables for age and sex
d1$left_young <- ifelse( d1$obstacle.left %in% c(1,2) , 1 , 0 )
d1$left_adult <- ifelse( d1$obstacle.left %in% c(3,4) , 1 , 0 )
d1$left_elderly <- ifelse( d1$obstacle.left %in% c(5,6) , 1 , 0 )
d1$right_young <- ifelse( d1$obstacle.right %in% c(1,2) , 1 , 0 )
d1$right_adult <- ifelse( d1$obstacle.right %in% c(3,4) , 1 , 0 )
d1$right_elderly <- ifelse( d1$obstacle.right %in% c(5,6) , 1 , 0 )
d1$young_diff <- d1$right_young - d1$left_young # Bradley–Terry model (Codierung mit -1 l 1 r etc.)
d1$adult_diff <- d1$right_adult - d1$left_adult # Bradley–Terry model (Codierung mit -1 l 1 r etc.)
d1$elderly_diff <- d1$right_elderly - d1$left_elderly # Bradley–Terry model (Codierung mit -1 l 1 r etc.)
d1$left_sex <- ifelse( d1$obstacle.left %in% c(1,3,5) , 1 , 0 ) # m = 0, f = 1
d1$right_sex <- ifelse( d1$obstacle.right %in% c(1,3,5) , 1 , 0 ) # m = 0, f = 1
d1$sex_diff <- d1$left_sex - d1$right_sex # 1 = f left m right, 0 = same, -1 = m left f right

# change coding for choice
d1$choice_left <- -1*d1$finish.lane + 2 # left = 1, right = 0

# change coding for start lane and visonset lane (left = 1, right = -1)
d1$start_left <- d1$start.lane
d1$start_left[d1$start.lane==2] <- rep(-1,nrow(d1))[d1$start.lane==2]
d1$visonset_left <- d1$visonset.lane 
d1$visonset_left[d1$visonset.lane==2] <- rep(-1,nrow(d1))[d1$visonset.lane==2]
d1$visonset_interaction <- d1$visonset_left*((d1$condition==4)-0.5) # similar to centering for dummy coding, so the intercept gives the mean of conditions (instead of mean of data)

# calculate percentage progress from trial number
d1$trial_progress <- (d1$trial.number-min(d1$trial.number))/max(d1$trial.number)

# change block to block_id, since block is taken in R
names(d1)[names(d1) == 'block'] <- 'block_id'

# create a subject number index to have a full set of integers
d1$sn_idx <- rethinking::coerce_index(d1$subject.number)

# transform condition into modality / abstract dummy coding
d1$modality <- ifelse( d1$condition %in% c(1,2) , 0 , 1 )
d1$abstraction <- ifelse( d1$condition %in% c(1,3) , 0 , 1 )

# create data frame containing only human-human trials
d1.hum <- d1[d1$obstacle.left<7,]
d1.hum <- d1.hum[d1.hum$obstacle.right<7,]

# create d1.effect data frame
d1.hum.effect <- d1.hum
d1.hum.effect$visonset_left <- d1.hum.effect$visonset_left/2
d1.hum.effect$sex_diff <- d1.hum.effect$sex_diff/2
d1.hum.effect$elderly_diff <- d1.hum.effect$elderly_diff/2
d1.hum.effect$young_diff <- d1.hum.effect$young_diff/2
d1.hum.effect$modality <- d1.hum.effect$modality - 0.5
d1.hum.effect$abstraction <- d1.hum.effect$abstraction - 0.5

d1.hum.effect %>% purrr::map(function(x)unique(x))

# ###############################################
# ### fit models ################################
# ###############################################

# glm fit (no priors for quick check) 
mres_lme = lme4::glmer(choice_left ~ 
                       1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction + (1 | sn_idx),
                       family=binomial, data=d1.hum.effect)

# model for baysian model
modelformula = choice_left ~
               1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction +
              (1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction | sn_idx)

# priors for bayesian model
modelprior = get_prior(modelformula, family=binomial, data=d1.hum.effect)
priors <- c(set_prior("lkj(2)", class = "cor"),
            set_prior("normal(0,3)", class = "b"),
            set_prior("exponential(1)", class = "sd", group = "sn_idx")) # alternative: cauchy(0,1)
# info on priors:
# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-the-regression-coefficients-in-logistic-regression-non-sparse-case
# http://jakewestfall.org/misc/SorensenEtAl.pdf

# fit model
mres = brm(modelformula,
           family=bernoulli,
           data=d1.hum.effect,
           control = list(adapt_delta=0.9),
           prior = priors)

# print results of fit
summary(mres)

# ###############################################
# ### results and interpretation ################
# ###############################################

# results and interpretation:
# Population-Level Effects: 
#                                   Estimate Est.Error l-95% CI u-95% CI
# Intercept                            -0.17      0.17    -0.50     0.17 = no strong general preference for a side
# visonset_left                         1.45      0.56     0.37     2.61 = small omission bias: higher likelihood of staying in the lane
# sex_diff                             -3.02      0.39    -3.83    -2.29 = strong effect for sex: females clearly preferred
# young_diff                            9.53      1.27     7.24    12.12 = extreme effect for age: advantage for children
# elderly_diff                         -7.43      1.03    -9.61    -5.63 = extreme effect for age: disadvantage for elderly

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# the above are main effects based on the mean of the abstraction levels
# and the mean of the modality levels, and they can be interpreted as
# population means.
# ----------------------------------------------------------------------

# abstraction                           0.38      0.30    -0.20     0.95
# visonset_left:abstraction            -1.26      0.83    -2.95     0.33
# sex_diff:abstraction                 -0.30      0.69    -1.68     1.02
# young_diff:abstraction               -0.01      1.34    -2.73     2.57
# elderly_diff:abstraction              1.81      1.09    -0.22     4.07

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# these are the effects of the abstraction level on the parameters.
# in order to make predictions w.r.t. abstraction==1 (image), we would
# add 0.5 of the estimated value to the main effects, while making
# predictions for abstraction==0 (text), we would subtract them.
# ----------------------------------------------------------------------

# modality                              0.16      0.34    -0.50     0.83
# sex_diff:modality                     0.70      0.71    -0.66     2.09
# young_diff:modality                   1.04      1.68    -2.31     4.21
# elderly_diff:modality                 1.99      1.49    -0.86     4.98

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# these are the effects of modality level on parameter. same calculation
# as with the abstraction parameters above.
# ----------------------------------------------------------------------

# modality:abstraction                  0.29      0.58    -0.84     1.43
# sex_diff:modality:abstraction        -0.84      1.30    -3.44     1.66
# young_diff:modality:abstraction      -0.22      2.07    -4.31     3.89
# elderly_diff:modality:abstraction     1.22      1.76    -2.23     4.76

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# these are interaction parameters and would be added (*0.5) if abstr==1 and
# modal==1 at the same time, and subtracted (*0.5) wenn abstr==0 && modal==0.
# ----------------------------------------------------------------------

# questions / to do:
# - we should make predictions with this model and see if they fit the data on a descriptive level
#   to make sure the model is making good predictions.
# - to make a statement about whether or not something is "significant", you suggested bayes factor,
#   but bayes factor is a tool for model comparison, right? so how can we use it to make statements
#   about individual parameters? couldn't we just fit 4 models (one without A/M, one with A, one
#   with M, one with both) and compare WAIC values?

# results with exponential(1) prior on sd for the sn_idx group:
# Population-Level Effects: 
#                                   Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                            -0.16      0.16    -0.47     0.14       4000 1.00
# visonset_left                         1.30      0.51     0.31     2.36       3751 1.00
# sex_diff                             -2.80      0.35    -3.52    -2.15       4000 1.00
# young_diff                            8.34      1.03     6.51    10.59       2133 1.00
# elderly_diff                         -6.72      0.88    -8.56    -5.13       2317 1.00

# abstraction                           0.37      0.27    -0.13     0.90       4000 1.00
# visonset_left:abstraction            -1.20      0.78    -2.77     0.31       4000 1.00
# sex_diff:abstraction                 -0.28      0.64    -1.55     0.98       4000 1.00
# young_diff:abstraction               -0.04      1.22    -2.43     2.29       4000 1.00
# elderly_diff:abstraction              1.61      0.97    -0.27     3.53       4000 1.00

# modality                              0.15      0.30    -0.45     0.75       4000 1.00
# sex_diff:modality                     0.68      0.65    -0.57     1.95       4000 1.00
# young_diff:modality                   1.10      1.50    -1.90     4.07       3513 1.00
# elderly_diff:modality                 1.78      1.33    -0.88     4.42       3487 1.00

# modality:abstraction                  0.23      0.54    -0.84     1.32       4000 1.00
# sex_diff:modality:abstraction        -0.73      1.20    -3.15     1.61       4000 1.00
# young_diff:modality:abstraction      -0.26      1.90    -3.96     3.54       4000 1.00
# elderly_diff:modality:abstraction     1.21      1.66    -2.07     4.54       4000 1.00

# Warning message:
#   There were 17 divergent transitions after warmup. Increasing adapt_delta above 0.9 may help.
# See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 

# ###############################################
# ### posterior predictive check ################
# ###############################################

pp_check(mres, "data", nsamples = NULL, group = "sn_idx")

# steps to manually sample from the distribution
# 1) sample all population-level parameters from the posterior distribution
# 2) sample 85 subjects (i.e., their individual parameter-offsets) from the distribution
# 3) randomly split them up into 42 and 43 and assign them to the modalities
# 4) get the distributions over the number of hum-hum trials, and sample two numbers
#    of hum-hum trials for each participant (one for each abstraction level)
# 5) take the exact trials from the data (obstacle left, obstacle right, starting lane) and encode
#    them in the age_diff / sex_diff coding scheme
# 6) for each trial, add the relevant parameters up and make a prediction
# 7) repeat the process 1000 times and make a distribution of predictions?

# ###############################################
# ### plots #####################################
# ###############################################

# plot posteriors
stanplot(mres, pars=c("b_"), type="areas", exact_match=FALSE)

# adding labels to the d1.hum data frame
d1.hum.plot = d1.hum %>% mutate(sex_diff = factor(sex_diff, levels=c(-1,0,1), labels=c('femaleLeftMaleRight','sameSex','femaleRightMaleLeft')),
                                visonset_left = factor(visonset_left, levels=c(-1,0,1), labels=c('visonsetLeft','none','visonsetRight')),
                                elderly_diff = factor(elderly_diff, levels=c(-1,0,1), labels=c('oldLeftNotRight','noneOrBothOld','oldRightNotLeft')),
                                young_diff = factor(young_diff, levels=c(-1,0,1), labels=c('youngLeftNotRight','noneOrBothYoung','youngRightNotLeft')),
                                abstraction = factor(abstraction, levels=c(0,1), labels=c('text','image')),
                                modality = factor(modality, levels=c(0,1), labels=c('desktop','vr')))

# plots: what are these exactly? data based / descriptive plots? where do the error bars come from?
# error message: "No summary function supplied, defaulting to `mean_se()" -- ?
ggplot(d1.hum.plot,aes(x=modality,y=choice_left,group=abstraction,color=abstraction))+stat_summary()+geom_hline(yintercept=0.5)
# ^ indicating a minimal bias towards the right lane, but probably not significant. in accordance with the model fit.

ggplot(d1.hum.plot,aes(x=sex_diff,y=choice_left,color=young_diff))+stat_summary()+geom_hline(yintercept=0.5)
ggplot(d1.hum.plot,aes(x=sex_diff,y=choice_left,color=elderly_diff))+stat_summary()+geom_hline(yintercept=0.5)
# ^ indicating medium to strong effects for gender, approx. 70/30 chance of saving the female over the male,
# provided they have the same age (and irrespective of other factors)

ggplot(d1.hum.plot,aes(x=visonset_left,y=choice_left,color=modality,shape=factor(abstraction)))+stat_summary()+geom_hline(yintercept=0.5)
# ^ this seems to indicate that only VR-text and not VR-image has a bias here, and it's a "panic reaction bias", not an "omission bias"
# contradicting the model fit, but probably ok since this graph only shows means irrespective of other factors, right?

ggplot(d1.hum.plot,aes(x=elderly_diff,y=choice_left,color=modality,shape=abstraction))+stat_summary()+geom_hline(yintercept=0.5)
ggplot(d1.hum.plot,aes(x=young_diff,y=choice_left,color=modality,shape=abstraction))+stat_summary()+geom_hline(yintercept=0.5)
# ^ indicates that neither modality not abstraction should have a strong effect on elderly_diff, young_diff

ggplot(d1.hum.plot,aes(x=sex_diff,y=choice_left,color=modality,shape=abstraction))+stat_summary()+geom_hline(yintercept=0.5)
# ^ no difference between text and image in the desktop setting, but a bit of weirdness in the VR settings but with large error bars,
# so probably no real effects here