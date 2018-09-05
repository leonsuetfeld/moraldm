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
d1 = read.csv("trialmatrix_2.csv")

str(d1)

# obstacle.left / obstacle.right coding scheme:
# 1 = girl, 2 = boy, 3 = woman, 4 = man, 5 = oldwoman, 6 = oldman, 7 = deer, 8 = goat, 9 = boar, 10 = nothing
# condition coding scheme:
# 1 = VR text fast, 2 = VR text slow, 3 = VR image fast, 4 = VR image slow

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
d1$age_center <- d1$age - mean(d1$age)

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
d1$sn_idx <- rethinking::coerce_index(d1$X..subject.number)

# transform condition into modality / abstract dummy coding
d1$abstraction <- ifelse( d1$condition %in% c(1,2) , 0 , 1 ) # 0 = text, 1 = image
d1$speed <- ifelse( d1$condition %in% c(1,3) , 1 , 0 )       # 0 = slow, 1 = fast

# filter out trials with an empty lane
d1.hum <- d1[d1$obstacle.left<7,]
d1.hum <- d1.hum[d1.hum$obstacle.right<7,]

# create d1.effect data frame
d1.hum.effect <- d1.hum
d1.hum.effect$visonset_left <- d1.hum.effect$visonset_left/2
d1.hum.effect$sex_diff <- d1.hum.effect$sex_diff/2
d1.hum.effect$elderly_diff <- d1.hum.effect$elderly_diff/2
d1.hum.effect$young_diff <- d1.hum.effect$young_diff/2
d1.hum.effect$speed <- d1.hum.effect$speed - 0.5
d1.hum.effect$abstraction <- d1.hum.effect$abstraction - 0.5
d1.hum.effect$sex <- d1.hum.effect$sex - 1.5

d1.hum.effect %>% purrr::map(function(x)unique(x))

# ###############################################
# ### define & fit bayesian models ##############
# ###############################################

# single level model
model_single = choice_left ~
  1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*speed*abstraction
modelprior_single = get_prior(model_single, family=binomial, data=d1.hum.effect)
priors_single <- c(set_prior("normal(0,3)", class = "b"))
mres_single = brm(model_single,
                  family=bernoulli,
                  data=d1.hum.effect,
                  control = list(adapt_delta=0.9),
                  prior = priors_single)

# full model
model_full = choice_left ~
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction +
  ((1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction | sn_idx)
modelprior = get_prior(model_full, family=binomial, data=d1.hum.effect)
priors <- c(set_prior("lkj(2)", class = "cor"),
            set_prior("normal(0,3)", class = "b"),
            set_prior("cauchy(0,1)", class = "sd", group = "sn_idx")) # alternatives: exponential(1) / cauchy(0,1)
mres_full = brm(model_full,
                family=bernoulli,
                data=d1.hum.effect,
                control = list(adapt_delta=0.9),
                prior = priors)

# sex model
model_sex = choice_left ~
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction + 
  (visonset_left + sex_diff + young_diff + elderly_diff)*sex + 
  ((1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction + 
  (visonset_left + sex_diff + young_diff + elderly_diff)*sex | sn_idx)
modelprior_sex = get_prior(model_sex, family=binomial, data=d1.hum.effect)
priors_sex <- c(set_prior("lkj(2)", class = "cor"),
                set_prior("normal(0,3)", class = "b"),
                set_prior("cauchy(0,1)", class = "sd", group = "sn_idx")) # alternatives: exponential(1) / cauchy(0,1)
mres_sex = brm(model_sex,
               family=bernoulli,
               data=d1.hum.effect,
               control = list(adapt_delta=0.9),
               prior = priors_sex)

# agesex model
model_agesex = choice_left ~
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction + 
  (visonset_left + sex_diff + young_diff + elderly_diff)*sex + 
  (visonset_left + sex_diff + young_diff + elderly_diff)*age + 
  ((1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction + 
  (visonset_left + sex_diff + young_diff + elderly_diff)*sex + 
  (visonset_left + sex_diff + young_diff + elderly_diff)*age | sn_idx)
modelprior_agesex = get_prior(model_agesex, family=binomial, data=d1.hum.effect)
priors_agesex <- c(set_prior("lkj(2)", class = "cor"),
                set_prior("normal(0,3)", class = "b"),
                set_prior("cauchy(0,1)", class = "sd", group = "sn_idx")) # alternatives: exponential(1) / cauchy(0,1)
mres_agesex = brm(model_agesex,
                  family=bernoulli,
                  data=d1.hum.effect,
                  control = list(adapt_delta=0.9),
                  prior = priors_agesex)

# ###############################################
# ### print results #############################
# ###############################################

summary(mres_full, maxsum=2)

# Population-Level Effects: 
#                                 Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
# Intercept                           0.05      0.12    -0.18     0.28       4000 1.00
# visonset_left                       0.50      0.29    -0.07     1.06       2652 1.00
# sex_diff                           -2.40      0.32    -3.08    -1.78       2322 1.00
# young_diff                          7.41      0.82     5.91     9.17       1729 1.00
# elderly_diff                       -6.03      0.65    -7.41    -4.84       1009 1.00
# speed                               0.28      0.21    -0.13     0.71       4000 1.00
# abstraction                         0.00      0.24    -0.49     0.48       4000 1.00
# visonset_left:speed                -0.14      0.41    -0.95     0.65       4000 1.00
# sex_diff:speed                      0.91      0.49    -0.05     1.86       4000 1.00
# young_diff:speed                   -2.89      0.80    -4.53    -1.36       4000 1.00
# elderly_diff:speed                  2.33      0.85     0.74     4.05       4000 1.00
# visonset_left:abstraction          -0.50      0.42    -1.34     0.34       4000 1.00
# sex_diff:abstraction                0.27      0.49    -0.70     1.22       4000 1.00
# young_diff:abstraction              1.63      0.81     0.08     3.27       4000 1.00
# elderly_diff:abstraction            4.20      0.78     2.74     5.87       4000 1.00
# speed:abstraction                  -0.04      0.41    -0.83     0.77       4000 1.00
# visonset_left:speed:abstraction    -0.00      0.78    -1.52     1.51       4000 1.00
# sex_diff:speed:abstraction          0.42      0.91    -1.34     2.22       4000 1.00
# young_diff:speed:abstraction        1.33      1.34    -1.27     3.98       4000 1.00
# elderly_diff:speed:abstraction      0.17      1.33    -2.40     2.77       4000 1.00

# predictions for abstraction==1 (image): main + 0.5*abstraction
# predictions for abstraction==0 (text): main - 0.5*abstraction

# plot posteriors
stanplot(mres_full, pars=c("b_"), type="areas", exact_match=FALSE)

# ###############################################
# ### model comparison ##########################
# ###############################################

kfold_full <- kfold(mres_full)

# compare models
# compare_ic(kfold_full, kfold_abst, kfold_mod, kfold_none)

# The ic is -2 * summed expected log pointwise predictive density,
# putting the elpd on the deviance scale so being similar to other
# information criteria metrics such as AIC or DIC.

# ###############################################
# ### posterior predictive check ################
# ###############################################

# I would treat these as "just for us" to make sure our model doesn't predict nonsense.

# THESE ARE STILL THE OLD ONES COPIED OVER FROM THE OTHER CODE, NO NEED TO CHECK HERE.

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