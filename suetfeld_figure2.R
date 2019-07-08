
library(tidyverse)
source('git_moraldm/publication_theme.R')

# ###############################################
# ### load & pre-process data ###################
# ###############################################


source('git_moraldm/suetfeld_preprocessing.R')
d1 = read.csv("trialmatrix.csv")
d1 = preprocessing(d1,study='1')

d2 = read.csv("trialmatrix_2.csv")
d2 = preprocessing(d2,study='2')


d1 = d1%>%mutate(study=1)
d2 = d2%>%mutate(study=2)
d = bind_rows(d1,d2)

#----
# Let's go! summarize the data



x = d%>%group_by(study,sn_idx,abstraction,speed,modality)%>%summarise( 
  young   = sum((young_diff  ==1 & choice_left==1) | (young_diff  ==-1 & choice_left==0))/sum(young_diff  !=0),
  elderly = sum((elderly_diff==1 & choice_left==1) | (elderly_diff==-1 & choice_left==0))/sum(elderly_diff!=0),
  sex     = sum((sex_diff    ==1 & choice_left==1) | (sex_diff    ==-1 & choice_left==0))/sum(sex_diff    !=0)
  )%>%ungroup()

x = x%>%tidyr::gather(key=effect,value=percentage,young,elderly,sex)

x = x%>%mutate(study = factor(study),
               abstraction=factor(abstraction,levels=c(0,1),labels=c("text","image")),
               speed=factor(speed,levels=c(0,1),labels=c("slow","fast")),
               modality=factor(modality,levels=c(0,1),labels=c("desktop","vr")),
)

x = x%>%tidyr::gather(key=factor,value=group,abstraction,modality,speed)

x$group = forcats::fct_relevel(x$group,"image","text","desktop","vr","slow","fast")
x = x%>% mutate(plt_position=interaction(study,group))

x$effect = forcats::fct_relevel(x$effect,"young","sex","elderly")


ggplot(x%>%subset(!is.na(group)),aes(x=group,
             y=percentage*100,
             group=study,
             color=study))+

  #ggbeeswarm::geom_quasirandom(dodge.width = 0.5,alpha=0.1)+
  geom_point(position=position_jitterdodge(dodge.width=0.7,jitter.width = 0.8),alpha=0.1)+
  geom_hline(yintercept=50,linetype='dashed')+
  stat_summary(position=position_dodge(width=0.7),fun.data = "mean_cl_boot",size=0.8,color='black')+
  facet_grid(.~effect)+
  ylim(c(0,100))+ggsci::scale_color_d3()+
  scale_shape_manual(values=c(15,16,17,18,4,13),guide=F)+
  theme_Publication()+ylab("")+
  theme(axis.text.x=element_text(angle=35,hjust=1))
