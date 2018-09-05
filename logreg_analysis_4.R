# ###############################################
# ### load library / top statements #############
# ###############################################

library(ggplot2)
library(brms)
library(bayesplot)
library(dplyr)
options(mc.cores = parallel::detectCores())
options(max.print = 2000)

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
# ### define & fit bayesian models ##############
# ###############################################

# single level model
model_single = choice_left ~
               1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction

# full model
model_full = choice_left ~
             1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction +
            (1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction | sn_idx)

# abstraction only model (no modality levels)
model_abst = choice_left ~
             1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*abstraction +
            (1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*abstraction | sn_idx)

# modality only model (no abstraction levels)
model_mod = choice_left ~
            1 + visonset_left + (sex_diff + young_diff + elderly_diff)*modality +
          (1 + visonset_left + (sex_diff + young_diff + elderly_diff)*modality | sn_idx)

# none model (no modality or abstraction levels)
model_none = choice_left ~
             1 + visonset_left + sex_diff + young_diff + elderly_diff +
            (1 + visonset_left + sex_diff + young_diff + elderly_diff | sn_idx)

# priors for bayesian model
modelprior_single = get_prior(model_single, family=binomial, data=d1.hum.effect)
priors_single <- c(set_prior("normal(0,3)", class = "b"))
modelprior = get_prior(model_full, family=binomial, data=d1.hum.effect)
priors <- c(set_prior("lkj(2)", class = "cor"),
            set_prior("normal(0,3)", class = "b"),
            set_prior("cauchy(0,1)", class = "sd", group = "sn_idx")) # alternatives: exponential(1) / cauchy(0,1)
# info on priors:
# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-the-regression-coefficients-in-logistic-regression-non-sparse-case
# http://jakewestfall.org/misc/SorensenEtAl.pdf

# fit single level model
mres_single = brm(model_single,
                  family=bernoulli,
                  data=d1.hum.effect,
                  control = list(adapt_delta=0.9),
                  prior = priors_single)

# fit full model
mres_full = brm(model_full,
                family=bernoulli,
                data=d1.hum.effect,
                control = list(adapt_delta=0.9),
                prior = priors)

# fit abst model
mres_abst = brm(model_abst,
                family=bernoulli,
                data=d1.hum.effect,
                control = list(adapt_delta=0.9),
                prior = priors)

# fit modu model
mres_mod = brm(model_mod,
               family=bernoulli,
               data=d1.hum.effect,
               control = list(adapt_delta=0.9),
               prior = priors)

# fit none model
mres_none = brm(model_none,
                family=bernoulli,
                data=d1.hum.effect,
                control = list(adapt_delta=0.9),
                prior = priors)

# ###############################################
# ### print results #############################
# ###############################################

summary(mres_full, maxsum=2)
# Population-Level Effects: 
#                                   Estimate Est.Error l-95% CI u-95% CI
# Intercept                            -0.17      0.17    -0.49     0.17 = no strong general preference for a side
# visonset_left                         1.47      0.56     0.42     2.65 = small omission bias: higher likelihood of staying in the lane
# sex_diff                             -3.01      0.40    -3.85    -2.30 = strong effect for sex: females clearly preferred
# young_diff                            9.48      1.28     7.19    12.30 = extreme effect for age: advantage for children
# elderly_diff                         -7.42      1.01    -9.54    -5.61 = extreme effect for age: disadvantage for elderly

# modality                              0.16      0.33    -0.51     0.80 => all these encode the total difference in parameter between the two modalities
# sex_diff:modality                     0.68      0.71    -0.68     2.09 => all have their 95% CIs on both sides of 0, i.e., no strong effects of modality
# young_diff:modality                   1.09      1.67    -2.16     4.37
# elderly_diff:modality                 1.99      1.54    -0.92     5.14

# abstraction                           0.37      0.30    -0.22     0.97 => all these encode the total difference in parameter between the two abstraction levels
# visonset_left:abstraction            -1.29      0.85    -3.07     0.27 => all have their 95% CIs on both sides of 0, i.e., no strong effects of abstraction
# sex_diff:abstraction                 -0.28      0.69    -1.66     1.01
# young_diff:abstraction               -0.04      1.38    -2.82     2.64
# elderly_diff:abstraction              1.83      1.07    -0.23     3.98
# modality:abstraction                  0.30      0.59    -0.86     1.48

# sex_diff:modality:abstraction        -0.81      1.25    -3.23     1.64 => no strong interaction effects either
# young_diff:modality:abstraction      -0.21      2.14    -4.41     4.05
# elderly_diff:modality:abstraction     1.25      1.81    -2.28     4.87

# predictions for abstraction==1 (image): main + 0.5*abstraction
# predictions for abstraction==0 (text): main - 0.5*abstraction

# Group-Level Effects: 
#   ~sn_idx (Number of levels: 85) 
#                   Estimate Est.Error l-95% CI u-95% CI
# sd(Intercept)         0.52      0.30     0.03     1.10 = low between-subject variance w.r.t. lane preference
# sd(visonset_left)     2.28      0.76     0.73     3.83 = higher between-subject variance w.r.t. omission bias
# sd(sex_diff)          0.71      0.53     0.02     1.95 = low between-subject variance w.r.t. value of gender / sex
# sd(young_diff)        3.57      2.44     0.10     7.91 = higher between-subject variance w.r.t. value of age
# sd(elderly_diff)      4.09      2.13     0.15     7.35 = higher between-subject variance w.r.t. value of age
# sd(modality)          0.82      0.56     0.04     2.01
# sd(abstraction)       0.43      0.33     0.02     1.23

# summary(mres_abst)
# summary(mres_mod)
# summary(mres_none)

# plot posteriors
stanplot(mres_full, pars=c("b_"), type="areas", exact_match=FALSE)

# ###############################################
# ### model comparison ##########################
# ###############################################

kfold_full <- kfold(mres_full)
kfold_abst <- kfold(mres_abst)
kfold_mod <- kfold(mres_mod)
kfold_none <- kfold(mres_none)

# compare models
compare_ic(kfold_full, kfold_abst, kfold_mod, kfold_none)

# mres_full              649.52 37.71
# mres_abst              644.69 37.76 = shared best models according to WAIC
# mres_mod               644.66 38.94 = shared best models according to WAIC
# mres_none              653.15 37.93
# mres_full - mres_abst    4.84  8.72
# mres_full - mres_mod     4.87 12.25
# mres_full - mres_none   -3.62 12.27
# mres_abst - mres_mod     0.03 10.49
# mres_abst - mres_none   -8.46  8.88
# mres_mod - mres_none    -8.49  8.92

# The ic is -2 * summed expected log pointwise predictive density,
# putting the elpd on the deviance scale so being similar to other
# information criteria metrics such as AIC or DIC.

# ###############################################
# ### posterior predictive check ################
# ###############################################

# I would treat these as "just for us" to make sure our model doesn't predict nonsense.

# all four conditions combined, single-level model vs. multi-level model
# pp_single = posterior_predict(mres_single, nsamples = 1000)
# pp = posterior_predict(mres_full, nsamples = 1000)
# d1.hum.effect$pp_single = apply(pp,2,mean)
# d1.hum.effect$pp = apply(pp,2,mean)
# d1.pp = reshape2::melt(d1.hum.effect,measure.vars=c("pp","pp_single","choice_left"))
# ggplot(d1.pp,aes(x=young_diff,y=value,group=variable,color=factor(variable)))+stat_summary(position=position_dodge(width=0.2))+geom_hline(yintercept=0.5)+facet_wrap(~abstraction)

# by condition, sampling "new random subjects"
dnew = d1.hum.effect %>% group_by(modality,abstraction,sex_diff,visonset_left,elderly_diff,young_diff)%>%summarise(choice_left=mean(choice_left))
pp_newsub = posterior_predict(mres_full,newdata = data.frame(dnew),allow_new_levels = TRUE,sample_new_levels="gaussian")
pp_single_newsub = posterior_predict(mres_single,newdata = data.frame(dnew))
dnew$pp_newsub = apply(pp_newsub,2,mean)
pp_single_newsub = posterior_predict(mres_single,newdata = data.frame(dnew))
dnew$pp_single_newsub = apply(pp_single_newsub,2,mean)
pp_uncert = posterior_predict(mres_full,newdata = data.frame(dnew),allow_new_levels = TRUE,sample_new_levels="uncertainty")
dnew$pp_uncert = apply(pp_uncert,2,mean)

# by condition, using the observed subjects
tmp = d1.hum.effect %>% group_by(modality,abstraction,sex_diff,visonset_left,elderly_diff,young_diff,sn_idx)%>%summarise(choice_left=mean(choice_left))
pp_samesub = posterior_predict(mres_full,newdata = data.frame(tmp),allow_new_levels = FALSE,sample_new_levels="gaussian")
tmp$pp_samesub = apply(pp_samesub,2,mean)
dnew$pp_samesub = tmp %>% group_by(modality,abstraction,sex_diff,visonset_left,elderly_diff,young_diff)%>%summarise(pp_samesub=mean(pp_samesub))%>%.$pp_samesub

# plot
dnew.pp = reshape2::melt(dnew,measure.vars=c("pp_samesub","pp_newsub","pp_single_newsub","choice_left"))
ggplot(dnew.pp,aes(x=young_diff,y=value,group=variable,color=factor(variable)))+stat_summary(position=position_dodge(width=0.2))+geom_hline(yintercept=0.5)+facet_grid(modality~abstraction)+ylim(c(0,1))
ggplot(dnew.pp,aes(x=elderly_diff,y=value,group=variable,color=factor(variable)))+stat_summary(position=position_dodge(width=0.2))+geom_hline(yintercept=0.5)+facet_grid(modality~abstraction)+ylim(c(0,1))
ggplot(dnew.pp,aes(x=sex_diff,y=value,group=variable,color=factor(variable)))+stat_summary(position=position_dodge(width=0.2))+geom_hline(yintercept=0.5)+facet_grid(modality~abstraction)+ylim(c(0,1))

# #############
# ### misc ####
# #############

# also report EDA results