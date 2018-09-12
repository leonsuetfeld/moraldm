library(cowplot)
library(ggplot2)
library(bayesplot)
library(brms)
library(dplyr)
load("stan_model_result")
 
#-----

           
m1 = mres_full%>%fixef%>%data.frame%>%tibble::rownames_to_column(var="effect")%>%mutate(study=1)
m2 = mres_agesexsdsvg%>%fixef%>%data.frame%>%tibble::rownames_to_column(var="effect")%>%mutate(study=2)

m = rbind(m1,m2) 

ranef_to_df = function (mres){
r = ranef(mres)$sn_idx
for (s in seq(dim(r)[1])){
  r[s,,] = r[s,,]+t(fixef(mres))
}
r = reshape2::melt(r) %>% magrittr::set_colnames(c('subject','column','effect','value'))
r = reshape2::dcast(r,subject+effect~column)
colnames(r)[c(5,6)] = c("Q2.5","Q97.5")
return(r)
}
r1 = mres_full%>%ranef_to_df%>%mutate(study=1)
r2 = mres_full%>%ranef_to_df%>%mutate(study=2)

r = rbind(r1,r2)


#-----
plot_names = function(x){
  rename_frame = list("elderly_diff:modality:abstraction"="± diff. elderly [abstr.+mod.]",
  "young_diff:modality:abstraction"=             "± diff. young [abstr.+mod.]",
  "sex_diff:modality:abstraction"=           "± diff. male [abstr.+mod.]",
  "visonset_left:modality:abstraction"=         "± diff. left lane [abstr.+mod.]",
  "elderly_diff:modality"=       "± diff. elderly [mod.]",
  "young_diff:modality"=     "± diff. young [mod.]",
  "sex_diff:modality"  =   "± diff. female [mod.]",
  "visonset:abstraction"= "± diff. left lane [mod.]",
  "elderly_diff:abstraction"="± diff. elderly [abstr.]",
  "young_diff:abstraction"="± diff. young [abstr.]",
  "sex_diff:abstraction"="± diff. female [abstr.]",
  "abstraction"="± omission bias [abstr.]",
  "visonset_left:abstraction"= "± diff. left lane [abstr.]",
  "elderly_diff"  = "adv. elderly",
  "young_diff"="adv. young",
  "sex_diff"="adv. female",
  "Intercept"= "omission bias",
  "visonset_diff"="adv. left lane")
x = x%>%mutate(effect=factor(effect,
  levels= names(rename_frame),
    labels = unname(rename_frame)))
return(x)
  
}
ggplot(m%>%filter(effect%in%c('Intercept','visonset_left','sex_diff','young_diff','elderly_diff'))%>%plot_names(),aes(y=Estimate,ymin=Q2.5,ymax=Q97.5,x=effect,group=study),)+
  geom_errorbar(position=position_dodge(width=0.5),width=0.1)+
  geom_point(position=position_dodge(width=0.5))+
  scale_y_continuous(labels=exp,breaks = log(10^seq(-9,9,3)),name = "Odds Ratios",limits = c(-10,10))+
  geom_hline(yintercept = 0,linetype="dotted")+
  ggbeeswarm::geom_quasirandom(data=r%>%plot_names(),alpha=0.1)+
  xlab('')+
  coord_flip()

       