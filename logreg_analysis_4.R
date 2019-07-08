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
d1.hum.effect$visonset_left <- d1.hum.effect$visonset_left
d1.hum.effect$sex_diff <- d1.hum.effect$sex_diff
d1.hum.effect$elderly_diff <- d1.hum.effect$elderly_diff
d1.hum.effect$young_diff <- d1.hum.effect$young_diff
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
                control = list(adapt_delta=0.95),
                chains=4,iter=8000,cores=4,warmup=2000,
                prior = priors)
summary(mres_full, maxsum=2)



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


# ###############################################
# ### print results #############################
# ###############################################

# Population-Level Effects: 
#                                   Estimate Est.Error l-95% CI u-95% CI
# Intercept                            -0.17      0.17    -0.49     0.15
# visonset_left                         1.44      0.56     0.38     2.57 !!!
# sex_diff                             -3.02      0.38    -3.83    -2.32 !!!
# young_diff                            9.51      1.27     7.27    12.19 !!!
# elderly_diff                         -7.41      1.03    -9.55    -5.57 !!!
# modality                              0.17      0.33    -0.48     0.84
# abstraction                           0.38      0.30    -0.20     0.99 !
# visonset_left:abstraction            -1.30      0.84    -3.02     0.31 !
# sex_diff:modality                     0.70      0.70    -0.66     2.04
# young_diff:modality                   1.08      1.71    -2.29     4.44
# elderly_diff:modality                 2.02      1.57    -1.03     5.15 !
# sex_diff:abstraction                 -0.30      0.67    -1.62     1.03
# young_diff:abstraction                0.00      1.32    -2.61     2.64
# elderly_diff:abstraction              1.85      1.10    -0.25     4.13 !
# modality:abstraction                  0.30      0.58    -0.83     1.46
# sex_diff:modality:abstraction        -0.78      1.23    -3.12     1.65
# young_diff:modality:abstraction      -0.18      2.14    -4.35     4.02
# elderly_diff:modality:abstraction     1.25      1.82    -2.26     4.93

# Group-Level Effects: 
#   ~sn_idx (Number of levels: 85) 
#                                     Estimate Est.Error l-95% CI u-95% CI
# sd(Intercept)                           0.53      0.30     0.03     1.13
# sd(visonset_left)                       2.29      0.75     0.91     3.80 !!
# sd(sex_diff)                            0.69      0.51     0.04     1.84
# sd(young_diff)                          4.00      2.49     0.13     8.34 !!!
# sd(elderly_diff)                        3.35      2.28     0.10     7.23 !!!
# sd(modality)                            0.80      0.56     0.03     2.02
# sd(abstraction)                         0.42      0.32     0.02     1.19
# sd(visonset_left:abstraction)           1.03      0.89     0.03     3.30 !
# sd(sex_diff:modality)                   1.02      0.89     0.04     3.35 !
# sd(young_diff:modality)                 6.33      5.25     0.06    16.09 !!!
# sd(elderly_diff:modality)               6.25      4.77     0.10    14.29 !!!
# sd(sex_diff:abstraction)                0.94      0.80     0.03     2.94
# sd(young_diff:abstraction)              2.86      2.55     0.05     8.79 !!
# sd(elderly_diff:abstraction)            1.65      1.55     0.05     5.74 !
# sd(modality:abstraction)                0.68      0.56     0.03     2.05
# sd(sex_diff:modality:abstraction)       1.30      1.34     0.04     4.98 !
# sd(young_diff:modality:abstraction)     3.31      4.25     0.05    15.48 !!!
# sd(elderly_diff:modality:abstraction)   2.06      2.52     0.04     9.14 !!

# predictions for abstraction==1 (image): main + 0.5*abstraction
# predictions for abstraction==0 (text): main - 0.5*abstraction

# plot posteriors
stanplot(mres_full, pars=c("b_"), type="areas", exact_match=FALSE)

# ###############################################
# ### plots #####################################
# ###############################################

# create plot: model posteriors
p1 <- plot_model(mres_full,
                 order.terms = c(1,2,3,4,5,7,8,12,13,14,6,9,10,11,15,16,17,18),
                 group.terms = c(1,1,1,1,1,3,2,2,3,3,3,2,2,2,4,4,4,4),
                 bpe.color = "black",
                 dot.size = 2.0,
                 dot.color = "black",
                 bpe = "mean",
                 show.intercept = TRUE,
                 width = 0.3,
                 prob.inner = .5,
                 prob.outer = .95,
                 axis.labels = c("± diff. elderly [abstr.+mod.]",
                                 "± diff. young [abstr.+mod.]",
                                 "± diff. male [abstr.+mod.]",
                                 "± diff. left lane [abstr.+mod.]",
                                 "± diff. elderly [mod.]",
                                 "± diff. young [mod.]",
                                 "± diff. female [mod.]",
                                 "± diff. left lane [mod.]",
                                 "± diff. elderly [abstr.]",
                                 "± diff. young [abstr.]",
                                 "± diff. female [abstr.]",
                                 "± omission bias [abstr.]",
                                 "± diff. left lane [abstr.]",
                                 "adv. elderly",
                                 "adv. young",
                                 "adv. female",
                                 "omission bias",
                                 "adv. left lane"),
                 title = "")


# create plots: distribution of random effects (with added fixed effects)
fe_sex_diff <- fixef(mres_full)[3,1]
sex_coefs <- ranef(mres_full)$sn_idx[,1,]%>%data.frame()%>%.$sex_diff+fe_sex_diff
# sex_coefs <- exp(sex_coefs)
p2 <- qplot(sex_coefs, xlim = c(0,4), xlab = "value of life (female), relative to male")

fe_young_diff <- fixef(mres_full)[4,1]
young_coefs <- ranef(mres_full)$sn_idx[,1,]%>%data.frame()%>%.$young_diff+fe_young_diff
# young_coefs <- exp(young_coefs)
p3 <- qplot(young_coefs, xlim = c(-2,12), xlab = "value of life (young), relative to adults")

fe_elederly_diff <- fixef(mres_full)[5,1]
elderly_coefs <- ranef(mres_full)$sn_idx[,1,]%>%data.frame()%>%.$elderly_diff+fe_elederly_diff
# elderly_coefs <- exp(elderly_coefs)
p4 <- qplot(elderly_coefs, xlim = c(-12,2), xlab = "value of life (elderly), relative to adults")

# plot
lay <- rbind(c(1,2),
             c(1,3),
             c(1,4))
grid.arrange(p1,p2,p3,p4, nrow=3, layout_matrix = lay, widths = 2:1)

# ###############################################
# ### hypothesis tests ##########################
# ###############################################

# bayes factor (BF10) calculation for all population-level parameters
# question: since it's a ratio of prior vs. posterior it's very dependent on our prior choice -- how do we make sure the priors are ok?

# intercept
h1 <- hypothesis(mres_full, "intercept = 0")
h1
BF10_1 <- 1/h1$hypothesis$Evid.Ratio
BF10_1

# Hypothesis Tests for class b:
#   Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
# 1 (intercept) = 0    -0.17      0.16    -0.49     0.15         NA        NA     
# ---
#   '*': The expected value under the hypothesis lies outside the 95%-CI.
# Posterior probabilities of point hypotheses assume equal prior probabilities.

# visonset_left
h2 <- hypothesis(mres_full, "visonset_left = 0")
h2
BF10_2 <- 1/h2$hypothesis$Evid.Ratio
BF10_2

# ###############################################
# ### accuracy ##################################
# ###############################################

# get 1000 random samples from the posterior (already 0&1's)
p_single = posterior_predict(mres_single,nsamples=1000)
p_full = posterior_predict(mres_full,nsamples=1000)

#calculate mean for each data point over posterior (mean prediction accuracy for each data point)
d1.hum.effect$accuracy_single = apply(apply(p_single,1,function(x)(x==d1.hum.effect$choice_left)),1,mean)
d1.hum.effect$accuracy_full = apply(apply(p_full,1,function(x)(x==d1.hum.effect$choice_left)),1,mean)
mean(d1.hum.effect$accuracy_single)
mean(d1.hum.effect$accuracy_full)

#plot it subjectwise (that is the accuracy each subject has)
ggplot(d1.hum.effect%>%group_by(sn_idx)%>%summarise(acc=mean(accuracy_single)),aes(x=0,y=acc))+ggbeeswarm::geom_quasirandom() + stat_summary(method="mean_cl_boot",color='red',size=1)+xlim(c(-1,1))+ylim(c(0,1))
ggplot(d1.hum.effect%>%group_by(sn_idx)%>%summarise(acc=mean(accuracy_full)),aes(x=0,y=acc))+ggbeeswarm::geom_quasirandom() + stat_summary(method="mean_cl_boot",color='red',size=1)+xlim(c(-1,1))+ylim(c(0,1))

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