
#alpha########################################################################################
prof = read.csv('data/profile/virome.vOTU.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[,map$ID]
shannon = diversity(t(prof))
vOTU = apply(prof, 2, function(x){sum(x!=0)})
Tab = data.frame(ID=names(shannon),Group=map$Group,Type=map$Type,Shannon=as.numeric(shannon),vOTU=as.numeric(vOTU),stringsAsFactors = F)

Or1=c('dental','saliva','faeces')
Tab$Group = factor(Tab$Group,Or1)
Or2 = unique(Tab$Type)
aa= combn(Or2,2)
my_comparisons = c()
for(x in 1:ncol(aa)){
  my_comparisons = c(my_comparisons,list(aa[,x]))
}

vm = c('Shannon','vOTU')
fig= c()
for(i in Or1){
  for(j in vm){
    TT = Tab[Tab$Group==i,]
    fig[[paste(i,j,sep='_')]] = ggplot(
                        TT,aes_string(x='Type',
                        #y=mapRate,
                        y=j,
                        #y=vOTU,
                        fill="Type",color="Type"))+
      geom_boxplot(outlier.size = 1)+
      labs(x='',
           #y='Mapping reads (%)'
           #y='Shannon diversity'
           y='The number of observed vOTUs'
      )+
      #facet_grid(~Group,scales = 'free',switch = 'x',space = 'free')+
      guides(fill=F,color=F)+
      scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
      scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
      stat_compare_means(comparisons=my_comparisons,size=2,tip.length = 0.01)+
      theme_classic()+
      theme(
        axis.title = element_text(
          size=8
        ),
        axis.text = element_text(
          size=6,
          angle=45,
          hjust=1
        )
      )
  }
}


cowplot::plot_grid(plotlist = fig)



#composition at the family level#################################################################################################
prof = read.csv('data/profile/virome.family.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[order(-rowMeans(prof)),map$ID]
prof = prof[,order(prof[1,],prof[2,],prof[3,])]
row.names(prof) = sub('.*f__','',row.names(prof))

dat = melt(as.matrix(prof))
dat = merge(dat,map,by.x = 'Var2',by.y = 'ID',all.x = T)
dat$Var1 = factor(dat$Var1,row.names(prof))
dat$Var2 = factor(dat$Var2,colnames(prof))

Or1=c('dental','saliva','faeces')
fig = c()
for(x in Or1){
  dd = dat[dat$Group==x,]
  dd = aggregate(dd$value,dd[,c('Var1','Type')],mean)
  dd$Var1 = factor(dd$Var1,row.names(prof))
  
  fig[[x]] = ggplot(dd)+
    geom_bar(aes(x=Type,y=x,fill=Var1),position = 'stack',stat = 'identity',width = 0.7)+
    scale_fill_manual(values = mycolor)+
    scale_y_continuous(expand = c(0,0))+
    guides(fill=guide_legend(ncol=1))+
    theme_classic()+
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      #axis.text.x = element_blank(),
      #axis.ticks.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust=1
      ),
      legend.key.height = unit(0.015, "cm"),
      legend.text=element_text(
        size=6
      )
    )
}
cowplot::plot_grid(plotlist = fig)

#test###
prof = read.csv('data/profile/virome.family.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[order(-rowSums(prof)),map$ID]
tt = c()
Or1=c('dental','saliva','faeces')
for(x in Or1){
  m=map[map$Group==x,]
  p=prof[,m$ID]
  tt[[x]] = kw_test(p,m,'Type')
}
View(tt$dental)



#dbRDA#######################################################################################
library(vegan) 
library(ape)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)


inf='data/profile/virome.vOTU.prof'
inp='data/mapping.file/mapping.file'
inc= 'Type'

fig=c()
for(x in Or1){
  lab=x
  
  Prof <- read.csv(file=inf,sep="\t",row.names=1,check.names = F,stringsAsFactors = F)  
  prof = Prof
  map <- read.csv(file=inp,sep="\t",check.names = F,stringsAsFactors = F)
  map = map[map$Group==lab,]
  #map = map[map$GA!=1,]
  #map[map[,2]=='HC',inc]<-1
  #map[map[,2]=='RA',inc]<-2
  map<- map[order(map[,1]),]
  prof <- t(prof[,map[,1]])
  prof <- prof[,colSums(prof)!=0]
  prof2 <- sqrt(prof)
  phen <- data.frame(Group=map[,inc],row.names = map[,1],stringsAsFactors = F)
  phen2 <- data.frame(ID = map[,1],Group=map[,inc],stringsAsFactors = F)
  
  decorana(prof2)
  
  ord <- capscale(prof2~Group,phen,distance = 'bray',scale=T) 
  p<-plot(ord)
  site.scores <- as.data.frame(scores(p, "site")) 
  site.scores$lab <- row.names(site.scores)
  site.scores <- merge(site.scores,phen2,by.x = 'lab',by.y = 'ID',on='left')
  
  species.scores <- as.data.frame(scores(p, "species")) 
  sp_or<-apply(species.scores,1,function(x){x[1]^2+x[2]^2})
  species.scores <- head(species.scores[order(sp_or,decreasing = T),],15)
  species.scores$lab <- species(row.names(species.scores))
  
  #ptest = permutest(ord,permu=999)
  ef=envfit(ord,phen,permu=999)
  
  cutoff=3
  fig[[lab]]=ggplot(data = site.scores,
                    aes(x=CAP1,y=CAP2)) +
    geom_point(aes(color=Group,fill=Group),size=0.5,shape=21)+
    stat_ellipse(aes(fill=Group),size=1, geom="polygon", level=0.8, alpha=0.3) +
    scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
    scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
    geom_vline(xintercept = 0,size=0.2,alpha=0.5,linetype=2) +
    geom_hline(yintercept = 0,size=0.2,alpha=0.5,linetype=2) +
    theme_bw()+
    labs(title="", x="CAP1", y="CAP2")+
    theme(
      panel.grid = element_blank(),
      text = element_text(
        size=6
      )
    )+
    annotate("text", x=-Inf, y=Inf,label=paste("Permutation test:\nP =",round(ef$factors$pvals,4),
                                               "\nR2 =",round(ef$factors$r,4)),size=2.5,hjust=-0.2,vjust=1.2)
  
}

cowplot::plot_grid(plotlist = fig)
