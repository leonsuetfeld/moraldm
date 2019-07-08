# ###############################################
# ### load library / top statements #############
# ###############################################

library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(brms)
library(bayesplot)
library(dplyr)
library(gridExtra)
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

# subject variables:
# sex: 1 = female, 2 = male

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
d1$age_center <- d1$age - mean(d1$age)
d1$sds17_center <- d1$sds17 - mean(d1$sds17)
d1$vg_center <- d1$hours.video.games - mean(d1$hours.video.games)

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
d1.hum.effect$visonset_left <- d1.hum.effect$visonset_left
d1.hum.effect$sex_diff <- d1.hum.effect$sex_diff
d1.hum.effect$elderly_diff <- d1.hum.effect$elderly_diff
d1.hum.effect$young_diff <- d1.hum.effect$young_diff
d1.hum.effect$abstraction <- d1.hum.effect$abstraction - 0.5

d1.hum.effect$speed <- d1.hum.effect$speed - 0.5
d1.hum.effect$sex <- d1.hum.effect$sex - 1.5

d1.hum.effect %>% purrr::map(function(x)unique(x))

# ###############################################
# ### define & fit bayesian models ##############
# ###############################################

# "full" model
# model_full = choice_left ~
#   (1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction +
#   ((1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction | sn_idx)
# modelprior = get_prior(model_full, family=binomial, data=d1.hum.effect)
# priors <- c(set_prior("lkj(2)", class = "cor"),
#             set_prior("normal(0,3)", class = "b"),
#             set_prior("cauchy(0,1)", class = "sd", group = "sn_idx")) # alternatives: exponential(1) / cauchy(0,1)
# mres_full = brm(model_full,
#                 family=bernoulli,
#                 data=d1.hum.effect,
#                 control = list(adapt_delta=0.9),
#                 chains=4,iter=8000,cores=4,warmup=2000,
#                 
#                 prior = priors)
# 
# summary(mres_full, maxsum=2)

# full age sex sds vg model
model_agesexsdsvg = choice_left ~
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction + 
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*sex + 
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*age_center + 
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*sds17_center + 
  (1 + visonset_left + sex_diff + young_diff + elderly_diff)*vg_center + 
  ((1 + visonset_left + sex_diff + young_diff + elderly_diff)*speed*abstraction | sn_idx)
modelprior_agesexsdsvg = get_prior(model_agesexsdsvg, family=binomial, data=d1.hum.effect)
priors_agesexsdsvg <- c(set_prior("lkj(2)", class = "cor"),
                        set_prior("normal(0,3)", class = "b"),
                        set_prior("cauchy(0,1)", class = "sd", group = "sn_idx"))
mres_agesexsdsvg = brm(model_agesexsdsvg,
                       family=bernoulli,
                       data=d1.hum.effect,
                       control = list(adapt_delta=0.9),
                       chains=4,iter=8000,cores=4,warmup=2000,
                       prior = priors_agesexsdsvg)


summary(mres_agesexsdsvg, maxsum=2)
save("mres_agesexsdsvg",file="stan_model_2.RData")




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

# plot posteriors
stanplot(mres_full, pars=c("b_"), type="areas", exact_match=FALSE)

# ###############################################
# ### plots #####################################
# ###############################################

s2_p1 <- plot_model(mres_agesexsdsvg,
                    order.terms = c(1,2,3,4,5,
                                    7,12,13,14,15,
                                    6,8,9,10,11,
                                    16,17,18,19,20),
                    group.terms = c(1,1,1,1,1,3,2,3,3,3,3,2,2,2,2,4,4,4,4,4),
                    terms = c("Intercept",
                              "visonset_left",
                              "sex_diff",
                              "young_diff",
                              "elderly_diff",
                              "abstraction",
                              "visonset_left.abstraction",
                              "sex_diff.abstraction",
                              "young_diff.abstraction",
                              "elderly_diff.abstraction",
                              "speed",
                              "visonset_left.speed",
                              "sex_diff.speed",
                              "young_diff.speed",
                              "elderly_diff.speed",
                              "speed.abstraction",
                              "visonset_left.speed.abstraction",
                              "sex_diff.speed.abstraction",
                              "young_diff.speed.abstraction",
                              "elderly_diff.speed.abstraction"),
                    bpe.color = "black",
                    dot.size = 2.0,
                    dot.color = "black",
                    bpe = "mean",
                    show.intercept = TRUE,
                    width = 0.3,
                    prob.inner = .5,
                    prob.outer = .95,
                    title = "")

s2_p2 <- plot_model(mres_agesexsdsvg,
                    order.terms = c(21,25,26,27,28,
                                    22,29,30,31,32,
                                    24,37,38,39,40,
                                    23,33,34,35,36),
                    group.terms = c(1,2,4,3,1,1,1,1,2,2,2,2,4,4,4,4,3,3,3,3),
                    rm.terms = c("Intercept",
                                 "visonset_left",
                                 "sex_diff",
                                 "young_diff",
                                 "elderly_diff",
                                 "abstraction",
                                 "visonset_left.abstraction",
                                 "sex_diff.abstraction",
                                 "young_diff.abstraction",
                                 "elderly_diff.abstraction",
                                 "speed",
                                 "visonset_left.speed",
                                 "sex_diff.speed",
                                 "young_diff.speed",
                                 "elderly_diff.speed",
                                 "speed.abstraction",
                                 "visonset_left.speed.abstraction",
                                 "sex_diff.speed.abstraction",
                                 "young_diff.speed.abstraction",
                                 "elderly_diff.speed.abstraction"),
                    bpe.color = "black",
                    dot.size = 2.0,
                    dot.color = "black",
                    bpe = "mean",
                    show.intercept = TRUE,
                    width = 0.3,
                    prob.inner = .5,
                    prob.outer = .95,
                    title = "")

s2_p2a <- plot_model(mres_agesexsdsvg,
                    order.terms = c(21,25,26,27,28),
                    group.terms = c(1,1,1,1,1),
                    terms = c("sex",
                              "visonset_left.sex",
                              "sex_diff.sex",
                              "young_diff.sex",
                              "elderly_diff.sex"),
                    axis.lim = c(0.01,50),
                    bpe.color = "black",
                    dot.size = 2.0,
                    dot.color = "black",
                    bpe = "mean",
                    show.intercept = TRUE,
                    width = 0.3,
                    prob.inner = .5,
                    prob.outer = .95,
                    title = "")

s2_p2b <- plot_model(mres_agesexsdsvg,
                    order.terms = c(22,29,30,31,32,
                                    24,37,38,39,40,
                                    23,33,34,35,36),
                    group.terms = c(1,2,3,1,1,1,1,2,2,2,2,3,3,3,3),
                    rm.terms = c("Intercept",
                                 "visonset_left",
                                 "sex_diff",
                                 "young_diff",
                                 "elderly_diff",
                                 "abstraction",
                                 "visonset_left.abstraction",
                                 "sex_diff.abstraction",
                                 "young_diff.abstraction",
                                 "elderly_diff.abstraction",
                                 "speed",
                                 "visonset_left.speed",
                                 "sex_diff.speed",
                                 "young_diff.speed",
                                 "elderly_diff.speed",
                                 "speed.abstraction",
                                 "visonset_left.speed.abstraction",
                                 "sex_diff.speed.abstraction",
                                 "young_diff.speed.abstraction",
                                 "elderly_diff.speed.abstraction",
                                 "sex",
                                 "visonset_left.sex",
                                 "sex_diff.sex",
                                 "young_diff.sex",
                                 "elderly_diff.sex"),
                    axis.lim = c(0.8,1.5),
                    bpe.color = "black",
                    dot.size = 2.0,
                    dot.color = "black",
                    bpe = "mean",
                    show.intercept = TRUE,
                    width = 0.3,
                    prob.inner = .5,
                    prob.outer = .95,
                    title = "")

# grid.arrange(s2_p1,s2_p2, nrow=1)

# create plots: distribution of random effects (with added fixed effects)
fe_sex_diff <- fixef(mres_agesexsdsvg)[3,1]
sex_coefs <- ranef(mres_agesexsdsvg)$sn_idx[,1,]%>%data.frame()%>%.$sex_diff+fe_sex_diff
# sex_coefs <- exp(sex_coefs)
s2_p3 <- qplot(sex_coefs, xlim = c(0,6), xlab = "value of life (female), relative to male")

fe_young_diff <- fixef(mres_agesexsdsvg)[4,1]
young_coefs <- ranef(mres_agesexsdsvg)$sn_idx[,1,]%>%data.frame()%>%.$young_diff+fe_young_diff
# young_coefs <- exp(young_coefs)
s2_p4 <- qplot(young_coefs, xlim = c(-4,16), xlab = "value of life (young), relative to adults")

fe_elederly_diff <- fixef(mres_agesexsdsvg)[5,1]
elderly_coefs <- ranef(mres_agesexsdsvg)$sn_idx[,1,]%>%data.frame()%>%.$elderly_diff+fe_elederly_diff
# elderly_coefs <- exp(elderly_coefs)
s2_p5 <- qplot(elderly_coefs, xlim = c(-16,4), xlab = "value of life (elderly), relative to adults")

# plot
lay <- rbind(c(1,1,2,2,4),
             c(1,1,3,3,5),
             c(1,1,3,3,6))
grid.arrange(s2_p1,s2_p2a,s2_p2b,s2_p3,s2_p4,s2_p5, nrow=3, ncol=3, layout_matrix = lay)

# ###############################################
# ### hypothesis tests ##########################
# ###############################################

# bayes factor (BF10) calculation for all population-level parameters
# question: since it's a ratio of prior vs. posterior it's very dependent on our prior choice -- how do we make sure the priors are ok?

# Intercept
h1 <- hypothesis(mres_agesexsdsvg, "Intercept = 0")
h1
BF10_1 <- 1/h1$hypothesis$Evid.Ratio
BF10_1

# visonset_left
h2 <- hypothesis(mres_agesexsdsvg, "visonset_left = 0")
h2
BF10_2 <- 1/h2$hypothesis$Evid.Ratio
BF10_2

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