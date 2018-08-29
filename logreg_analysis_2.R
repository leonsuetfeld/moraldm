# ###############################################
# ### load library / top statements #############
# ###############################################

library(magrittr)
library(dplyr)
library(brms)
library(loo)
library(ggplot2)
library(bayesplot)
options(mc.cores = parallel::detectCores())
options(max.print = 2000)

# library(rethinking)
# rstan_options(auto_write = TRUE)

# ###############################################
# ### load & pre-process data ###################
# ###############################################

setwd("/DATEN/PhD/research/VRdriveBA/_data_preStudy/data_csv/")
d1 = read.csv("trialmatrix.csv")
str(d1)

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

# calculate percentage progress from trial number
d1$trial_progress <- (d1$trial.number-min(d1$trial.number))/max(d1$trial.number)

# change block to block_id, since block is taken in R
names(d1)[names(d1) == 'block'] <- 'block_id'

# create a subject number index to have a full set of integers
d1$sn_idx <- coerce_index(d1$subject.number)

# create data frame containing only human-human trials
d1.hum <- d1[d1$obstacle.left<7,]
d1.hum <- d1.hum[d1.hum$obstacle.right<7,]

# convert data frame for bernoulli
d1.hum.sum = d1.hum %>% group_by(sn_idx, condition, visonset_left, sex_diff, left_young, left_adult, left_elderly, right_young, right_adult, right_elderly) %>% summarise(count_left = sum(choice_left), total_left=length(choice_left))
d1.hum.sum = d1.hum.sum %>% mutate(young_diff = right_young - left_young, elderly_diff = right_elderly - left_elderly)

d1.effect = d1.hum.sum
d1.effect$visonset_left = d1.effect$visonset_left - 0.5
d1.effect$sex_diff = d1.effect$sex_diff - 0.5
d1.effect$elderly_diff = d1.effect$elderly_diff - 0.5
d1.effect$young_diff = d1.effect$young_diff - 0.5
d1.effect$visonset_left = d1.effect$visonset_left - 0.5

str(d1.hum.sum)

# ###############################################
# ### 
# ###############################################

# fit models (bene)
mres0 = brm(count_left | trials(total_left) ~ 
            1,
            family=binomial, data=d1.hum.sum, control=list(adapt_delta = 0.9))
mres = brm(count_left | trials(total_left) ~ 
           1 + visonset_left + sex_diff + condition + young + elderly +
          (1 + visonset_left + sex_diff + condition + young + elderly | sn_idx),
           family=binomial, data=d1.hum.sum, control=list(adapt_delta = 0.9))

# fit model without priors (leon)
mres1 = brm(count_left | trials(total_left) ~ 
            1 + visonset_left + sex_diff + young + elderly +
           (1 + visonset_left + sex_diff + young + elderly | sn_idx) +
          (-1 + sex_diff + young + elderly | condition),
            family=binomial, data=d1.hum.sum, control=list(adapt_delta = 0.99))

# fit model with priors (leon)
get_prior(count_left | trials(total_left) ~ 
          1 + visonset_left + sex_diff + young_diff + elderly_diff +
         (1 + visonset_left + sex_diff + young_diff + elderly_diff | sn_idx) +
        (-1 + sex_diff + young_diff + elderly_diff | condition),
          family=binomial, data=d1.hum.sum)
prior_mres2 <- c(set_prior("normal(0,3)", class = "b", coef = "visonset_left"),
                 set_prior("normal(0,3)", class = "b", coef = "sex_diff"),
                 set_prior("normal(0,3)", class = "b", coef = "young"),
                 set_prior("normal(0,3)", class = "b", coef = "elderly"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "Intercept"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "visonset_left"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "sex_diff"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "young"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "elderly"),
                 set_prior("cauchy(0,1)", class = "sd", group = "condition", coef = "sex_diff"),
                 set_prior("cauchy(0,1)", class = "sd", group = "condition", coef = "young"),
                 set_prior("cauchy(0,1)", class = "sd", group = "condition", coef = "elderly"))
mres2 = brm(count_left | trials(total_left) ~ 
            1 + visonset_left + sex_diff + young_diff + elderly_diff +
           (1 + visonset_left + sex_diff + young_diff + elderly_diff | sn_idx) +
          (-1 + sex_diff + young + elderly | condition),
            family=binomial, data=d1.hum.sum, prior=prior_mres2, control=list(adapt_delta = 0.98))

# fit some other model (bene)
prior_mres3 <- c(set_prior("normal(0,3)", class = "b", coef = "visonset_left"),
                 set_prior("normal(0,3)", class = "b", coef = "sex_diff"),
                 set_prior("normal(0,3)", class = "b", coef = "young"),
                 set_prior("normal(0,3)", class = "b", coef = "elderly"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "Intercept"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "visonset_left"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "sex_diff"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "young"),
                 set_prior("cauchy(0,1)", class = "sd", group = "sn_idx", coef = "elderly"),
                 set_prior("cauchy(0,1)", class = "sd", group = "condition", coef = "sex_diff"),
                 set_prior("cauchy(0,1)", class = "sd", group = "condition", coef = "young"),
                 set_prior("cauchy(0,1)", class = "sd", group = "condition", coef = "elderly"))

mres3 = brm(modelformula,
           family=binomial,data=d1.hum.sum,
           control = list(adapt_delta=0.9),
           prior = c(prior(student_t(4,0,3),class=b), prior(lkj(2),class=cor)))

mres3 = brm(count_left | trials(total_left) ~
            1 + visonset_left + visonset_interaction + (sex_diff + young_diff + elderly_diff) * condition +
           (1 + visonset_left + visonset_interaction + (sex_diff + young_diff + elderly_diff) * condition | sn_idx),
            family=binomial, data=d1.hum.sum, prior=prior_mres3, control=list(adapt_delta = 0.8))

stanplot(mres3, pars=c("b_"), type="areas", exact_match=FALSE)
summary(mres3)

# show results
mres0
mres
mres1
mres2

# model comparison
loo(mres0)
loo(mres)
loo(mres1)
loo(mres2)
