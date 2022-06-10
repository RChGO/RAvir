
#alpha diversity#################################################################################
prof = read.csv('data/profile/virome.vOTU.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[,map$ID]
shannon1 = diversity(t(prof))


prof = read.csv('data/profile/metaphlan.L7.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[,map$ID]
shannon2 = diversity(t(prof))

dat = as.data.frame(t(rbind(shannon1,shannon2)))
dat$ID = row.names(dat)
dat =merge(dat,map,by= 'ID')

dat$Group = factor(dat$Group,Or1)
ggplot(dat,aes(x=shannon1,y=shannon2,color=Group,fill=Group))+
  geom_point(shape=21,size=1)+
  geom_smooth(method = 'lm')+
  scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
  scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
  theme_bw()

lab=Or1[3]
cor.test(dat$shannon1[dat$Group==lab],dat$shannon2[dat$Group==lab],method = 'spearman')
Or2 = c("Control","Onset","Treated")
fig=c()
for(x in Or1){
  dd=dat[dat$Group==x,]
  dd$Type = factor(dd$Type,Or2)
  fig[[x]]=ggplot(dd,aes(x=shannon1,y=shannon2,color=Type,fill=Type))+
    geom_point(shape=21,size=1)+
    geom_smooth(method = 'lm')+
    labs(x='Shannon index (Virus)',y='Shannon index (Archaea/Bacteria)')+
    scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
    scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
    theme_bw()
  
}

cowplot::plot_grid(plotlist = fig,nrow=1)

lab=Or2[1]
for(x in Or1){
  aa=c()
  for(y in Or2){
    aa = c(aa,cor.test(dat$shannon1[dat$Group==x&dat$Type==y],
                       dat$shannon2[dat$Group==x&dat$Type==y],method = 'spearman')$estimate)
  }
  print(aa)
}