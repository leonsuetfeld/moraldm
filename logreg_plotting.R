library(cowplot)
library(ggplot2)
library(bayesplot)
library(brms)
library(dplyr)
load("stan_model_result")
 


theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
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

#-----

           
m1 = mres_full%>%fixef%>%data.frame%>%tibble::rownames_to_column(var="effect")%>%mutate(study=1)
m2 = mres_agesexsdsvg%>%fixef%>%data.frame%>%tibble::rownames_to_column(var="effect")%>%mutate(study=2)

m = rbind(m1,m2) 

label_effects = function(x){
  
  tmp = x %>% mutate(twoway = grepl(':',effect),
                     threeway = grepl(':.+:',effect))
  tmp = tmp %>% mutate(twoway = ifelse(threeway,FALSE,twoway))
  tmp = tmp %>% mutate(effecttype = ifelse(threeway,'3xinteraction',ifelse(twoway,'2xinteraction','maineffect')))
  grouplist = c("abstraction","modality","speed","sex","age_center","sds17_center","vg_center")
  
  k = 0
  for (g in grouplist){
      k = k + 1
      tmp_rename = tmp%>%mutate(effect = gsub("sex_diff","gender_diff",effect))
      tmp_rename = tmp_rename%>%mutate(group_tmp = ifelse(grepl(g,effect),k,NaN))
      tmp[!is.na(tmp_rename$group_tmp),'group'] = k
  }
  if(length(unique(tmp$group)) == 4){
    grouplist = grouplist[1:3]
  }
  tmp$group = factor(tmp$group,labels = grouplist)
  return(tmp)
}

ranef_to_df = function (mres){
r = ranef(mres)$sn_idx
for (s in seq(dim(r)[1])){
  tmp = t(fixef(mres))
  ix = dimnames(tmp)[[2]] %in% dimnames(r)[[3]]
  tmp = tmp[,ix]
  r[s,,] = r[s,,]+tmp
}
r = reshape2::melt(r) %>% magrittr::set_colnames(c('subject','column','effect','value'))
r = reshape2::dcast(r,subject+effect~column)
colnames(r)[c(5,6)] = c("Q2.5","Q97.5")
return(r)
}
r1 = mres_full%>%ranef_to_df%>%mutate(study=1)
r2 = mres_agesexsdsvg%>%ranef_to_df%>%mutate(study=2)

r = rbind(r1,r2)
m = m%>%label_effects
r = r%>%label_effects


#-----
plot_names = function(x,short=FALSE){
  rename_frame = list(
    "modality:abstraction"  = "± diff. left lane [abstr.+mod.]",
    "speed:abstraction"     = "± diff. left lane [speed+abstr.]",
    
    "abstraction"  = "± diff. left lane [abstr.]",
    "speed"        = "± diff. left lane [speed]",
    "modality"     = "± diff. left lane [mod.]",
    
    "Intercept"    = "adv. left lane",
    
    "visonset_left"= "omission bias",
    
    "sex_diff"     = "adv. female",
    "young_diff"   = "adv. young",
    "elderly_diff" = "adv. elderly",
    
    
    
    "sds17_center" = "SDS 17",
    "vg_center"    = "video game exp.",
    "age_center"   = "subject age",
    "sex"          = "subject gender",
    
  "visonset_left:abstraction"          = "± diff. omission bias [abstr.]",
  "visonset_left:modality"             = "± diff. omission bias [mod.]",
  "visonset_left:speed"                = "± diff. omission bias [speed]",
  "visonset_left:speed:abstraction"    = "± diff. omission bias [abstr.+speed]",
  "visonset_left:modality:abstraction" = "± diff. omission bias [abstr.+mod.]",
  "visonset_left:sds17_center"         = "± diff. omission bias [SDS17]",
  "visonset_left:age_center"           = "± diff. omission bias [subject age]",
  "visonset_left:vg_center"            = "± diff. omission bias [videogame]",
  "visonset_left:sex"                  = "± diff. omission bias [subject sex]",
  
  "sex_diff:abstraction"          = "± diff. female [abstr.]",
  "sex_diff:modality"             = "± diff. female [mod.]",
  "sex_diff:speed"                = "± diff. female [speed]",
  "sex_diff:speed:abstraction"    = "± diff. female [abstr.+speed]",
  "sex_diff:modality:abstraction" = "± diff. female [abstr.+mod.]",
  "sex_diff:sds17_center"         = "± diff. female [SDS17]",
  "sex_diff:age_center"           = "± diff. female [subject age]",
  "sex_diff:vg_center"            = "± diff. female [videogame]",
  "sex_diff:sex"                  = "± diff. female [subject sex]",
  
  
  "young_diff:abstraction"          = "± diff. young [abstr.]",
  "young_diff:modality"             = "± diff. young [mod.]",
  "young_diff:speed"                = "± diff. young [speed]",
  "young_diff:speed:abstraction"    = "± diff. young [abstr.+speed]",
  "young_diff:modality:abstraction" = "± diff. young [abstr.+mod.]",
  "young_diff:sds17_center"         = "± diff. young [SDS17]",
  "young_diff:age_center"           = "± diff. young [subject age]",
  "young_diff:vg_center"            = "± diff. young [videogame]",
  "young_diff:sex"                  = "± diff. young [subject sex]",
  
  
  
  "elderly_diff:abstraction"          = "± diff. elderly [abstr.]",
  "elderly_diff:modality"             = "± diff. elderly [mod.]",
  "elderly_diff:speed"                = "± diff. elderly [speed]",
  "elderly_diff:speed:abstraction"    = "± diff. elderly [abstr.+speed]",
  "elderly_diff:modality:abstraction" = "± diff. elderly [abstr.+mod.]",
  "elderly_diff:sds17_center"         = "± diff. elderly [SDS17]",
  "elderly_diff:age_center"           = "± diff. elderly [subject age]",
  "elderly_diff:vg_center"            = "± diff. elderly [videogame]",
  "elderly_diff:sex"                  = "± diff. elderly [subject sex]"
  

  

  )
  
  
  labels = unname(rename_frame )
if (short){
  labels = gsub(" \\[.*\\]","",labels)
}
x = x%>%mutate(effect=factor(effect,
  levels= names(rename_frame),
    labels = labels))
return(x)

}


logToLog10 = scales::trans_new('logToLog10',transform=function(x)log10(exp(x)),inverse=function(x)log(10^x))

plot_misc = list(geom_hline(yintercept = 0,linetype="dotted"),
                 xlab(''),
                 ggsci::scale_color_d3(guide=F),
                 #scale_y_continuous(trans=logToLog10),
                 scale_y_continuous(labels=function(x)10^x,breaks=seq(-6,6,3),name = "Odds Ratios",limits = c(-4.8,6.5),position="right"),
                 annotation_logticks(base=10,sides="r",scaled=T),
                 theme_Publication(),
                 theme(axis.text.x=element_text(angle=90,hjust=0)),
                 theme(axis.text.y=element_text(angle=90,hjust=0))
                 )



## main effects only
ggplot(m%>%subset(effect %in% c("Intercept","visonset_left","young_diff","sex_diff","elderly_diff"))%>%plot_names(),
       aes(y=Estimate %>%exp%>%log10,
           ymin=Q2.5 %>%exp%>%log10,
           ymax=Q97.5%>%exp%>%log10,x=effect,group=study,color=factor(study)))+
  ggbeeswarm::geom_quasirandom(data=r%>%subset(effect %in% c("Intercept","visonset_left","young_diff","sex_diff","elderly_diff"))%>%plot_names(),alpha=0.1,dodge.width=0.5)+
  geom_errorbar(position=position_dodge(width=0.5),width=0.1,color="black")+
  geom_point(position=position_dodge(width=0.5),color="black")+
  scale_color_discrete(guide=F)+
  plot_misc


## Interactions only
p1 = ggplot(m%>%subset(effecttype!="3xinteraction"&!is.na(group)&!(group %in% c("age_center","sds17_center","vg_center")))%>%plot_names(short=T),
       aes(y=Estimate %>%exp%>%log10,
           ymin=Q2.5 %>%exp%>%log10,
           ymax=Q97.5%>%exp%>%log10,x=effect,group=study,color=factor(study)))+
  ggbeeswarm::geom_quasirandom(data=r%>%subset(effecttype!="3xinteraction"&!is.na(group)&!(group %in% c("age_center","sds17_center","vg_center")))%>%plot_names(short=T),alpha=0.1,dodge.width=0.5)+
  geom_errorbar(position=position_dodge(width=0.5),width=0.1,color="black")+
  geom_point(position=position_dodge(width=0.5),color="black")+
  scale_color_discrete(guide=F)+
  facet_wrap(~group,ncol = 4,scales = "free_x")+
  plot_misc
p1


p2 = ggplot(m%>%subset(effecttype!="3xinteraction"&!is.na(group)&(group %in% c("age_center","sds17_center","vg_center")))%>%plot_names(short=T),
       aes(y=Estimate %>%exp%>%log10,
           ymin=Q2.5 %>%exp%>%log10,
           ymax=Q97.5%>%exp%>%log10,x=effect,group=study,color=factor(study)))+
  ggbeeswarm::geom_quasirandom(data=r%>%subset(effecttype!="3xinteraction"&!is.na(group)&(group %in% c("age_center","sds17_center","vg_center")))%>%plot_names(short=T),alpha=0.1,dodge.width=0.5)+
  geom_errorbar(position=position_dodge(width=0.5),width=0.1,color="black")+
  geom_point(position=position_dodge(width=0.5),color="black")+
  scale_color_discrete(guide=F)+
  facet_wrap(~group,ncol = 3,scales = "free")+
  plot_misc+
  scale_y_continuous(labels=function(x)(10^x),breaks=log10(c(0.33,0.5,1,2,3)),name = "Odds Ratios",limits = c(-0.5,0.5),position="right")

p2
scale = 3
ggsave("figure4a.pdf",p1,width=scale*4,height=scale*2)
ggsave("figure4b.pdf",p2,width=scale*3,height=scale*2)
