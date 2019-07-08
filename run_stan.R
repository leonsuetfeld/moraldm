
#library(sjPlot)
#library(sjlabelled)
#library(sjmisc)
#library(ggplot2)
library(brms)
#library(bayesplot)
library(dplyr)
#library(gridExtra)
options(mc.cores = parallel::detectCores())
options(max.print = 2000)

Sys.setenv(USE_CXX14 = 1)

# ###############################################
# ### load & pre-process data ###################
# ###############################################

setwd("/home/predatt/benehi/projects/moraldm/")
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
d1$sex_diff <- d1$right_sex - d1$left_sex # 1 = m left f right, 0 = same, -1 = f left m right

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
d1$modality <- ifelse( d1$condition %in% c(1,2) , 0 , 1 ) # 0 = desktop, 1 = vr
d1$abstraction <- ifelse( d1$condition %in% c(1,3) , 0 , 1 ) # 0 = slow, 1 = fast

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

# info on priors:
# https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-the-regression-coefficients-in-logistic-regression-non-sparse-case
# http://jakewestfall.org/misc/SorensenEtAl.pdf

# full model
model_full = choice_left ~
             1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction +
            (1 + visonset_left + visonset_left:abstraction + (sex_diff + young_diff + elderly_diff)*modality*abstraction | sn_idx)
modelprior = get_prior(model_full, family=binomial, data=d1.hum.effect)
priors <- c(set_prior("lkj(2)", class = "cor"),
            set_prior("normal(0,3)", class = "b"),
            set_prior("cauchy(0,1)", class = "sd", group = "sn_idx")) # alternatives: exponential(1) / cauchy(0,1)
mres_full = brm(model_full,
                family=bernoulli,
                data=d1.hum.effect,
                control = list(adapt_delta=0.9),chains=4,cores=4,iter=8000,warmup=2000,
                prior = priors)

save('mres_full',file='moral_model.Rdata')
