
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


theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# ###############################################
# ### load & pre-process data ###################
# ###############################################

setwd("/DATEN/PhD/research/VRdriveBA/_data_preStudy/data_csv/")
preprocess = function(filename){
d1= read.csv(filename)

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
if( "subject.number" %in% colnames(d1)){
d1$sn_idx =d1$subject.number#<- rethinking::coerce_index(d1$subject.number)
}else{
d1$sn_idx =d1$X..subject.number#<- rethinking::coerce_index(d1$subject.number)
}

# filter out trials with empty lane
d1 <- d1[d1$obstacle.left<7,]
d1 <- d1[d1$obstacle.right<7,]

return(d1)
}
d1 = preprocess("trialmatrix.csv")
# transform condition into modality / abstract dummy coding
d1$modality <- ifelse( d1$condition %in% c(1,2) , 0 , 1 ) # 0 = desktop, 1 = vr
d1$abstraction <- ifelse( d1$condition %in% c(1,3) , 0 , 1 ) # 0 = slow, 1 = fast

d2 = preprocess("trialmatrix_2.csv")
d2$abstraction <- ifelse( d2$condition %in% c(1,2) , 0 , 1 ) # 0 = text, 1 = image
d2$speed <- ifelse( d2$condition %in% c(1,3) , 1 , 0 )       # 0 = slow, 1 = fast

d1 = d1%>%mutate(study=1)
d2 = d2%>%mutate(study=2)
d = bind_rows(d1,d2)
#----
# Let's go! summarize the data


x = d%>%group_by(study,sn_idx,abstraction,speed,modality)%>%summarise( 
  young   = sum((young_diff  ==1 & choice_left==0) | (young_diff  ==-1 & choice_left==1))/sum(young_diff  !=0),
  elderly = sum((elderly_diff==1 & choice_left==0) | (elderly_diff==-1 & choice_left==1))/sum(elderly_diff!=0),
  sex     = sum((sex_diff    ==1 & choice_left==0) | (sex_diff    ==-1 & choice_left==1))/sum(sex_diff    !=0)
  )%>%ungroup()
x = x%>%tidyr::gather(key=effect,value=percentage,young,elderly,sex)
x = x%>%mutate(study = factor(study),
               abstraction=factor(abstraction,levels=c(0,1),labels=c("text","image")),
               speed=factor(speed,levels=c(0,1),labels=c("slow","fast")),
               modality=factor(modality,levels=c(0,1),labels=c("desktop","vr")),
)
x = x%>%tidyr::gather(key=factor,value=group,modality,speed)
ggplot(x,aes(x=effect,
             y=percentage,
             group=interaction(group,study),
             shape=group,
             color=study))+
  stat_summary(position=position_dodge(width=0.5))+
  ggbeeswarm::geom_quasirandom(dodge.width = 0.5,alpha=0.1)+
  ylim(c(0,1))+facet_grid(.~abstraction)+
  theme_Publication()

