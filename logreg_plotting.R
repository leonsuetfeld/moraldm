library(cowplot)
library(ggplot2)
library(bayesplot)
library(brms)
library(dplyr)
#load("stan_model_result")
load("C:/Users/behinger/cloud/PhD/Supervision/Leon/stan_model_1.RData")
load("C:/Users/behinger/cloud/PhD/Supervision/Leon/stan_model_2.RData")

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
    "modality:abstraction"  = " left lane bias [abstr.+mod.]",
    "speed:abstraction"     = " left lane bias [speed+abstr.]",
    
    "abstraction"  = " left lane bias [abstr.]",
    "speed"        = " left lane bias [speed]",
    "modality"     = " left lane bias [mod.]",
    
    "Intercept"    = "adv. left lane",
    
    "visonset_left"= "omission bias",
    
    "sex_diff"     = "adv. female",
    "young_diff"   = "adv. young",
    "elderly_diff" = "adv. elderly",
    
    
    
    "sds17_center" = "SDS 17",
    "vg_center"    = "video game",
    "age_center"   = "subject age",
    "sex"          = "subject gender",
    
  "visonset_left:abstraction"          = " omission bias [abstr.]",
  "visonset_left:modality"             = " omission bias [mod.]",
  "visonset_left:speed"                = " omission bias [speed]",
  "visonset_left:speed:abstraction"    = " omission bias [abstr.+speed]",
  "visonset_left:modality:abstraction" = " omission bias [abstr.+mod.]",
  "visonset_left:sds17_center"         = " omission bias [SDS17]",
  "visonset_left:age_center"           = " omission bias [subject age]",
  "visonset_left:vg_center"            = " omission bias [videogame]",
  "visonset_left:sex"                  = " omission bias [subject sex]",
  
  "sex_diff:abstraction"          = " gender bias [abstr.]",
  "sex_diff:modality"             = " gender bias [mod.]",
  "sex_diff:speed"                = " gender bias [speed]",
  "sex_diff:speed:abstraction"    = " gender bias [abstr.+speed]",
  "sex_diff:modality:abstraction" = " gender bias [abstr.+mod.]",
  "sex_diff:sds17_center"         = " gender bias [SDS17]",
  "sex_diff:age_center"           = " gender bias [subject age]",
  "sex_diff:vg_center"            = " gender bias [videogame]",
  "sex_diff:sex"                  = " gender bias [subject sex]",
  
  
  "young_diff:abstraction"          = " young bias [abstr.]",
  "young_diff:modality"             = " young bias [mod.]",
  "young_diff:speed"                = " young bias [speed]",
  "young_diff:speed:abstraction"    = " young bias [abstr.+speed]",
  "young_diff:modality:abstraction" = " young bias [abstr.+mod.]",
  "young_diff:sds17_center"         = " young bias [SDS17]",
  "young_diff:age_center"           = " young bias [subject age]",
  "young_diff:vg_center"            = " young bias [videogame]",
  "young_diff:sex"                  = " young bias [subject sex]",
  
  
  
  "elderly_diff:abstraction"          = " elderly bias [abstr.]",
  "elderly_diff:modality"             = " elderly bias [mod.]",
  "elderly_diff:speed"                = " elderly bias [speed]",
  "elderly_diff:speed:abstraction"    = " elderly bias [abstr.+speed]",
  "elderly_diff:modality:abstraction" = " elderly bias [abstr.+mod.]",
  "elderly_diff:sds17_center"         = " elderly bias [SDS17]",
  "elderly_diff:age_center"           = " elderly bias [subject age]",
  "elderly_diff:vg_center"            = " elderly bias [videogame]",
  "elderly_diff:sex"                  = " elderly bias [subject sex]"
  

  

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
                 scale_y_continuous(labels=function(x)10^x,breaks=seq(-6,6,3),name = "Odds Ratios",limits = c(-4.7,6.0),position="right"),
                 #annotation_logticks(base=10,sides="r",scaled=T),
                 theme_Publication(),
                 theme(axis.text.x=element_text(angle=90,hjust=0)),
                 theme(axis.text.y=element_text(angle=90,hjust=0))
                 )



## main effects only
p_main = ggplot(m%>%subset(effect %in% c("Intercept","visonset_left","young_diff","sex_diff","elderly_diff"))%>%plot_names(),
       aes(y=Estimate %>%exp%>%log10,
           ymin=Q2.5 %>%exp%>%log10,
           ymax=Q97.5%>%exp%>%log10,x=effect,group=study,color=factor(study)))+
  ggbeeswarm::geom_quasirandom(data=r%>%subset(effect %in% c("Intercept","visonset_left","young_diff","sex_diff","elderly_diff"))%>%plot_names(),alpha=0.1,dodge.width=0.5)+
  geom_errorbar(position=position_dodge(width=0.5),width=0.1)+
  geom_point(position=position_dodge(width=0.5))+
  scale_color_discrete(guide=F)+
  plot_misc

## Interactions only
p1 = ggplot(m%>%subset(effecttype!="3xinteraction"&!is.na(group)&!(group %in% c("age_center","sds17_center","vg_center")))%>%plot_names(short=T),
       aes(y=Estimate %>%exp%>%log10,
           ymin=Q2.5 %>%exp%>%log10,
           ymax=Q97.5%>%exp%>%log10,x=effect,group=study,color=factor(study)))+
  ggbeeswarm::geom_quasirandom(data=r%>%subset(effecttype!="3xinteraction"&!is.na(group)&!(group %in% c("age_center","sds17_center","vg_center")))%>%plot_names(short=T),alpha=0.1,dodge.width=0.5)+
  geom_errorbar(position=position_dodge(width=0.5),width=0.1)+
  geom_point(position=position_dodge(width=0.5))+
  scale_color_discrete(guide=F)+
  facet_wrap(~group,ncol = 1,scales = "free")+
  plot_misc+
  scale_y_continuous(labels=function(x)10^x,breaks=seq(-6,6,3),name = "Odds Ratios",limits = c(-3.0,4.5),position="right")
  
p1


p2 = ggplot(m%>%subset(effecttype!="3xinteraction"&!is.na(group)&(group %in% c("age_center","sds17_center","vg_center")))%>%plot_names(short=T),
       aes(y=Estimate %>%exp%>%log10,
           ymin=Q2.5 %>%exp%>%log10,
           ymax=Q97.5%>%exp%>%log10,x=effect,group=study,color=factor(study)))+
  ggbeeswarm::geom_quasirandom(data=r%>%subset(effecttype!="3xinteraction"&!is.na(group)&(group %in% c("age_center","sds17_center","vg_center")))%>%plot_names(short=T),alpha=0.1,dodge.width=0.5)+
  geom_errorbar(position=position_dodge(width=0.5),width=0.1)+
  geom_point(position=position_dodge(width=0.5))+
  scale_color_discrete(guide=F)+
  facet_wrap(~group,ncol = 1,scales = "free")+
  plot_misc+
  scale_y_continuous(labels=function(x)(10^x),breaks=log10(c(0.33,0.5,1,2,3)),name = "Odds Ratios",limits = c(-0.5,0.5),position="right")

p2
scale = 3

ggsave("figure3.pdf",p_main,width=scale,height=scale*2)
ggsave("figure4a.pdf",p1,width=scale*0.75,height=scale*2*2)
ggsave("figure4b.pdf",p2,width=scale,height=scale*2*1.5)



