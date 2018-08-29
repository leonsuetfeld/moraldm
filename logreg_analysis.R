# ###############################################
# ### load library / top statements #############
# ###############################################

library(rethinking)
options(mc.cores = parallel::detectCores())
options(max.print = 2000)
rstan_options(auto_write = TRUE)

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

# ###############################################
# ### pooled model: all conditions in one #######
# ###############################################

# trim data frame
d1.hum.all <- d1.hum[,c("choice_left","visonset_left","sex_diff",
                        "left_young","left_adult","left_elderly",
                        "right_young","right_adult","right_elderly")]

# fit model
m.pooled <- map2stan(
  alist(
    choice_left ~ dbinom(1, p),
    logit(p) <- lb
              + visonset_left*ob
              - sex_diff*vf
              + (right_young-left_young)*vy 
              + (right_elderly-left_elderly)*ve,
    c(lb, ob, vf, vy, ve) ~ dnorm(0, 3)
  ),
  data=d1.hum.all,
  warmup=500, iter=3500,
  chains=8, cores=3,
  control=list(adapt_delta=0.95))

# results
tracerplot(m.pooled, window=c(50, 3500)) # these look fine
precis(m.pooled)
pairs(m.pooled)
show(m.pooled)

# variable names:
# lb = lane bias (intercept), ob = omission bias, vf = value female
# vy = value young, va = value adult (not used), ve = value elderly 

# ###############################################
# ### partially pooled model: condition #########
# ###############################################

# trim data frame
d1.hum.pp_cond <- d1.hum[,c("choice_left","condition","visonset_left","sex_diff",
                            "left_young","left_adult","left_elderly",
                            "right_young","right_adult","right_elderly")]

# fit full model
m.partial_pooled.cond <- map2stan(
  alist(
    choice_left ~ dbinom(1, p),
    logit(p) <- lb + lb_cond[condition]
    + visonset_left*(ob+ob_cond[condition])
    - sex_diff*(vf+vf_cond[condition])
    + (right_young-left_young)*(vy+vy_cond[condition])
    + (right_elderly-left_elderly)*(ve+ve_cond[condition]),
    c(lb,ob,vf,vy,ve) ~ dnorm(0, 3),
    lb_cond[condition] ~ dnorm(0, sigma_lbc),
    ob_cond[condition] ~ dnorm(0, sigma_obc),
    vf_cond[condition] ~ dnorm(0, sigma_vfc),
    vy_cond[condition] ~ dnorm(0, sigma_vyc),
    ve_cond[condition] ~ dnorm(0, sigma_vec),
    c(sigma_lbc, sigma_obc, sigma_vfc, sigma_vyc, sigma_vec) ~ dcauchy(0, 1)
  ),
  data=d1.hum.pp_cond,
  warmup=500, iter=3500,
  chains=8, cores=3,
  control=list(adapt_delta=0.95))

# fit reduced model
m.partial_pooled.cond.reduced <- map2stan(
  alist(
    choice_left ~ dbinom(1, p),
    logit(p) <- lb
    + visonset_left*ob
    - sex_diff*(vf+vf_cond[condition])
    + (right_young-left_young)*(vy+vy_cond[condition])
    + (right_elderly-left_elderly)*(ve+ve_cond[condition]),
    c(lb,ob,vf,vy,ve) ~ dnorm(0, 3),
    vf_cond[condition] ~ dnorm(0, sigma_vfc),
    vy_cond[condition] ~ dnorm(0, sigma_vyc),
    ve_cond[condition] ~ dnorm(0, sigma_vec),
    c(sigma_vfc, sigma_vyc, sigma_vec) ~ dcauchy(0, 1)
  ),
  data=d1.hum.pp_cond,
  warmup=500, iter=3500,
  chains=8, cores=3,
  control=list(adapt_delta=0.95))

# results full model
tracerplot(m.partial_pooled.cond, window=c(50, 3500)) # are these chains healthy?
precis(m.partial_pooled.cond, depth=2)
pairs(m.partial_pooled.cond)
show(m.partial_pooled.cond)

# results reduced model
tracerplot(m.partial_pooled.cond.reduced, window=c(50, 3500)) # are these chains healthy?
precis(m.partial_pooled.cond.reduced, depth=2)
pairs(m.partial_pooled.cond.reduced)
show(m.partial_pooled.cond.reduced)

# comparison
coeftab(m.partial_pooled.cond, m.partial_pooled.cond.reduced)
compare(m.partial_pooled.cond, m.partial_pooled.cond.reduced)

# ###############################################
# ### partially pooled model: condition & sn ####
# ###############################################

# trim data frame
d1.hum.pp_cond_sn <- d1.hum[,c("sn_idx","choice_left","condition","visonset_left","sex_diff",
                               "left_young","left_adult","left_elderly",
                               "right_young","right_adult","right_elderly")]

# fit full model
m.partial_pooled.cond_sn <- map2stan(
  alist(
    choice_left ~ dbinom(1, p),
    logit(p) <- lb + lb_cond[condition] + lb_sn[sn_idx]
    + visonset_left*(ob+ob_cond[condition]+ob_sn[sn_idx])
    - sex_diff*(vf+vf_cond[condition]+vf_sn[sn_idx])
    + (right_young-left_young)*(vy+vy_cond[condition]+vy_sn[sn_idx])
    + (right_elderly-left_elderly)*(ve+ve_cond[condition]+ve_sn[sn_idx]),
    c(lb,ob,vf,vy,ve) ~ dnorm(0, 3),
    lb_cond[condition] ~ dnorm(0, sigma_lbc),
    ob_cond[condition] ~ dnorm(0, sigma_obc),
    vf_cond[condition] ~ dnorm(0, sigma_vfc),
    vy_cond[condition] ~ dnorm(0, sigma_vyc),
    ve_cond[condition] ~ dnorm(0, sigma_vec),
    c(sigma_lbc, sigma_obc, sigma_vfc, sigma_vyc, sigma_vec) ~ dcauchy(0, 1),
    lb_sn[sn_idx] ~ dnorm(0, sigma_lbsn),
    ob_sn[sn_idx] ~ dnorm(0, sigma_obsn),
    vf_sn[sn_idx] ~ dnorm(0, sigma_vfsn),
    vy_sn[sn_idx] ~ dnorm(0, sigma_vysn),
    ve_sn[sn_idx] ~ dnorm(0, sigma_vesn),
    c(sigma_lbsn, sigma_obsn, sigma_vfsn, sigma_vysn, sigma_vesn) ~ dcauchy(0, 1)
  ),
  data=d1.hum.pp_cond_sn,
  warmup=500, iter=3500,
  chains=8, cores=3,
  control=list(adapt_delta=0.95))

# fit reduced model
m.partial_pooled.cond_sn.reduced <- map2stan(
  alist(
    choice_left ~ dbinom(1, p),
    logit(p) <- lb
    + visonset_left*ob
    - sex_diff*(vf+vf_cond[condition]+vf_sn[sn_idx])
    + (right_young-left_young)*(vy+vy_cond[condition]+vy_sn[sn_idx])
    + (right_elderly-left_elderly)*(ve+ve_cond[condition]+ve_sn[sn_idx]),
    c(lb,ob,vf,vy,ve) ~ dnorm(0, 3),
    vf_cond[condition] ~ dnorm(0, sigma_vfc),
    vy_cond[condition] ~ dnorm(0, sigma_vyc),
    ve_cond[condition] ~ dnorm(0, sigma_vec),
    c(sigma_vfc, sigma_vyc, sigma_vec) ~ dcauchy(0, 1),
    vf_sn[sn_idx] ~ dnorm(0, sigma_vfsn),
    vy_sn[sn_idx] ~ dnorm(0, sigma_vysn),
    ve_sn[sn_idx] ~ dnorm(0, sigma_vesn),
    c(sigma_vfsn, sigma_vysn, sigma_vesn) ~ dcauchy(0, 1)
  ),
  data=d1.hum.pp_cond_sn,
  warmup=500, iter=3500,
  chains=8, cores=3,
  control=list(adapt_delta=0.95))

# results full model
tracerplot(m.partial_pooled.cond_sn, window=c(50, 3500)) # are these chains healthy?
precis(m.partial_pooled.cond_sn, depth=1)
pairs(m.partial_pooled.cond_sn, pars = c("ob", "vf", "vy", "ve", "sigma_obsn", "sigma_vfsn", "sigma_vysn", "sigma_vesn"))
show(m.partial_pooled.cond_sn)

# results reduced model
tracerplot(m.partial_pooled.cond_sn.reduced, window=c(50, 3500)) # are these chains healthy?
precis(m.partial_pooled.cond_sn.reduced, depth=1)
pairs(m.partial_pooled.cond_sn.reduced)
show(m.partial_pooled.cond_sn.reduced)

# comparison
coeftab(m.partial_pooled.cond_sn, m.partial_pooled.cond_sn.reduced)
compare(m.partial_pooled.cond_sn, m.partial_pooled.cond_sn.reduced)

# ###############################################
# ### overall model comparison ##################
# ###############################################


# ###############################################
# ### useful info from the documentation ########
# ###############################################

# If you have many unknowns in your Stan program, then the pairs() plot may be illegible.
# There is a pars argument to the pairs() function that allows you to specify a subset of parameters.
# Also, there is an include argument that can be set to FALSE that will result in the complement
# of the pars argument being included in the pairs() plot. You should definitely include variance
# (or standard deviation, etc.) parameters, since they illustrate the funnel phenomenon most often.
# In addition, you should include at least some lower-level parameters in hierarchical models.

# library(rstan)
# funnel <- stan_demo("funnel", seed = 12345)   # has 5 divergent transitions
# pairs(funnel, pars = c("y", "x[1]", "lp__"), las = 1) # below the diagonal
# funnel_reparam <- stan_demo("funnel_reparam") # has no divergent transitions


# ###############################################
# ### pooled models: all cond in one + separate #
# ###############################################

# the following is my "old way" of doing it (splitting the dataset by condition
# and performing individual fits).

# pooled model, all conditions in one, fixed va (value adult) at 0 (threw it out)
d1.hum.all <- d1.hum[,c("choice_left","visonset_left","sex_diff",
                        "left_young","left_adult","left_elderly",
                        "right_young","right_adult","right_elderly")]
m.pooled <- map2stan(
  alist(
    choice_left ~ dbinom( 1 , p ),
    logit(p) <- lb + visonset_left*ob - sex_diff*vf +
      (right_young-left_young)*vy + (right_elderly-left_elderly)*ve,
    c(lb,ob,vf,vy,ve) ~ dnorm( 0, 3 )
  ),
  data=d1.hum.all, warmup=1000 , iter=6000 , chains=1 , cores=1 )

precis(m.pooled)
pairs(m.pooled)

# pooled model, desktop text
d1.hum.dtt <- d1.hum[d1.hum$condition==1,]
d1.hum.dtt <- d1.hum.dti[,c("choice_left","visonset_left","sex_diff",
                            "left_young","left_adult","left_elderly",
                            "right_young","right_adult","right_elderly")]
m.pooled.dtt <- map2stan(
  alist(
    choice_left ~ dbinom( 1 , p ),
    logit(p) <- lb - sex_diff*vf +
      (right_young-left_young)*vy + (right_elderly-left_elderly)*ve,
    c(lb,vf,vy,ve) ~ dnorm( 0, 3 )
  ),
  data=d1.hum.dtt, warmup=1000 , iter=6000 , chains=1 , cores=1 )

precis(m.pooled.dtt)

# pooled model, desktop image
d1.hum.dti <- d1.hum[d1.hum$condition==2,]
d1.hum.dti <- d1.hum.dti[,c("choice_left","visonset_left","sex_diff",
                            "left_young","left_adult","left_elderly",
                            "right_young","right_adult","right_elderly")]
m.pooled.dti <- map2stan(
  alist(
    choice_left ~ dbinom( 1 , p ),
    logit(p) <- lb - sex_diff*vf +
      (right_young-left_young)*vy + (right_elderly-left_elderly)*ve,
    c(lb,vf,vy,ve) ~ dnorm( 0, 3 )
  ),
  data=d1.hum.dti, warmup=1000 , iter=6000 , chains=1 , cores=1 )

precis(m.pooled.dti)

# pooled model, desktop text
d1.hum.vrt <- d1.hum[d1.hum$condition==3,]
d1.hum.vrt <- d1.hum.dti[,c("choice_left","visonset_left","sex_diff",
                            "left_young","left_adult","left_elderly",
                            "right_young","right_adult","right_elderly")]
m.pooled.vrt <- map2stan(
  alist(
    choice_left ~ dbinom( 1 , p ),
    logit(p) <- lb + visonset_left*ob - sex_diff*vf +
      (right_young-left_young)*vy + (right_elderly-left_elderly)*ve,
    c(lb,ob,vf,vy,ve) ~ dnorm( 0, 3 )
  ),
  data=d1.hum.vrt, warmup=1000 , iter=6000 , chains=1 , cores=1 )

precis(m.pooled.vrt)

# pooled model, desktop image
d1.hum.vri <- d1.hum[d1.hum$condition==4,]
d1.hum.vri <- d1.hum.dti[,c("choice_left","visonset_left","sex_diff",
                            "left_young","left_adult","left_elderly",
                            "right_young","right_adult","right_elderly")]
m.pooled.vri <- map2stan(
  alist(
    choice_left ~ dbinom( 1 , p ),
    logit(p) <- lb + visonset_left*ob - sex_diff*vf +
      (right_young-left_young)*vy + (right_elderly-left_elderly)*ve,
    c(lb,ob,vf,vy,ve) ~ dnorm( 0, 3 )
  ),
  data=d1.hum.vri, warmup=1000 , iter=6000 , chains=1 , cores=1 )

precis(m.pooled.vri)

coeftab(m.pooled.dtt, m.pooled.dti, m.pooled.vrt, m.pooled.vri)
