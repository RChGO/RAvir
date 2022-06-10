
#rarefaction#####################################################################################
Stat  = read.csv('data/profile/alpha_raref.txt',sep='\t',stringsAsFactors = F,check.names = F)
map = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)

#mapf = map
mapf = map[map$Group=='dental',]
mapf$Group = mapf$Type

dat = merge(Stat,mapf,by='ID')
dat = aggregate(dat[,4:5],dat[,c(1,2,6)],mean)

m = aggregate(dat$Obs,dat[,c('Step','Group')],mean)
sd = aggregate(dat$Obs,dat[,c('Step','Group')],function(x){sd(x)/sqrt(length(x)-1)})
m$sd = sd$x
m = m[m$Step<4000000,]

#Or=c('dental','saliva','faeces')
#m$Group = factor(m$Group,Or)

f1 = ggplot(m,aes(x=Step,y=x,group=Group,color=Group))+
  geom_line()+
  scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
  scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
  geom_errorbar(aes(ymin=x-sd, ymax=x+sd), width=500, position=position_dodge(.9))+
  labs(
    #y='Shannon index',
    y='The number of Observed vOTUs',
       x='Sampling')+
  theme_bw()

cowplot::plot_grid(f1,f2,nrow = 1)

#reads map rate##########################################################################################
Stat  = read.csv('data/data_generation/00.data/all.contig.stat',sep='\t',stringsAsFactors = F,check.names = F)[,1:3]
prof = read.csv('data/data_generation/00.data/virome.reads_count',sep='\t',stringsAsFactors = F,check.names = F)
aa = colSums(prof[,-1])
aa = data.frame(sampleid=names(aa),ViromeReads=aa)

Stat = merge(Stat,aa,by='sampleid',all=T)
Stat[,'ViromeData(bp)'] = Stat$ViromeReads*100
Stat[,'ViromePercent(%)'] = Stat[,'ViromeData(bp)']/Stat[,'CleanData(bp)']*100

#summary(Stat$`ViromePercent(%)`[grep('^s',Stat$sampleid)])
#write.table(Stat,'ReadsMap.rate',sep='\t',row.names = F,quote=F)


#overlap of vOTU########################################################################################
prof = read.csv('data/profile/virome.vOTU.prof',sep='\t',stringsAsFactors = F,
                check.names = F,row.names = 1)
map = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = melt(as.matrix(prof))

prof = merge(prof,map,by.x='Var2',by.y = 'ID')
dat = unique(prof[prof$value!=0,c('Var1','Group')])
ls = c()
for(x in unique(dat$Group)){
  ls[[x]] = dat$Var1[dat$Group==x]
  cat(x,length(ls[[x]]),'\n')
}

length(intersect(ls$dental,ls$saliva))/length(ls$dental)
length(intersect(ls$dental,ls$saliva))/length(ls$saliva)
length(intersect(ls$dental,ls$faeces))/length(ls$dental)
length(intersect(ls$saliva,ls$faeces))/length(ls$saliva)

length(ls$Faeces)


library(VennDiagram)
venn.plot <- venn.diagram(ls,filename = NULL, col = "transparent",
                          fill = c("#f76f6f", "#9dbcd6", "#7faa81"), alpha = 0.5,
                          label.col = c("darkred", "white", "darkblue", "white",
                                        "white", "white", "darkgreen"),
                          cex = 1.5, cat.default.pos = "text",
                          cat.col = c("darkred", "darkblue", "darkgreen"),
                          cat.cex = 2,cat.dist = c(0.06, 0.06, 0.03), cat.pos = 0)

pdf(file="analysis/1.alpha_div/Venn.pdf")
grid.draw(venn.plot)
dev.off()


#alpha####################################################################################################
prof = read.csv('data/data_generation/00.data/ReadsMap.rate',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[map$ID,]
mapRate = prof$`ViromePercent(%)`
prof = read.csv('data/profile/virome.vOTU.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[,map$ID]
shannon = diversity(t(prof))
vOTU = apply(prof, 2, function(x){sum(x!=0)})
Tab = data.frame(ID=names(shannon),Group=map$Group,mapRate=mapRate,Shannon=as.numeric(shannon),vOTU=as.numeric(vOTU),stringsAsFactors = F)

Or=c('dental','saliva','faeces')
Tab$Group = factor(Tab$Group,Or)

aa= combn(Or,2)
my_comparisons = c()
for(x in 1:ncol(aa)){
  my_comparisons = c(my_comparisons,list(aa[,x]))
}


f2 = ggplot(Tab,aes(x=Group,
                    #y=mapRate,
                    #y=shannon,
                    y=vOTU,
                    fill=Group,color=Group))+
  geom_boxplot(outlier.size = 1)+
  labs(x='',
       #y='Mapping reads (%)'
       #y='Shannon diversity'
       y='The number of observed vOTUs'
       )+
  scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
  scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
  stat_compare_means(comparisons=my_comparisons,size=2)+
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
f2
cowplot::plot_grid(f1,f2,f3,nrow =2)  



#pca,Group#################################################################################################
library(fpc)
library(ade4)
library(ggrepel)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(vegan)

prof = read.csv('data/profile/virome.vOTU.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[,map$ID]

vector_fd<-function(tab,sam){
  ms<-max(sam[,1]^2+sam[,2]^2)
  mt<-max(tab[,1]^2+tab[,2]^2)
  sqrt(ms)/sqrt(mt)
}
#=================
incol<-'Group'
sl <- 1
xx <- 1.2
yy <- 1.2
inlift <- "Axis1"
inright <- "Axis2"
#=================
data1 = map
data1 = data1[order(data1[,1]),]
data1 = data1[!is.na(data1[,incol]),]
data1 = data.frame(row.names = data1[,1],Group=data1[,incol],check.names=F,stringsAsFactors=F)
dat <- prof
dat <- dat[,row.names(data1)]

dat <- dat[rowSums(dat)!=0,]
data <- sweep(dat, 2, apply(dat,2,sum), "/")
data <- t(sqrt(data))
data <- data[order(row.names(data)),]
data.dudi <- dudi.pca(data, center=TRUE, scale=F, scan=F, nf=10)
data2 <- data.dudi$li

row.names(data2) = row.names(data1)
classified_c = as.character(unique(data1[,1]))
data3 = merge(data2,data1,by = "row.names")
row.names(data3) = data3[,1]
data3 = data3[,-1]

adonis1<-adonis(data ~ data1[order(row.names(data1)),],permutations = 2000,method = "bray")

phenotype <- data1[,0+1]
f = classified_c
Type <- factor(phenotype,levels=f)

m = data.dudi$li
n = data.dudi$c1

lift_data = m[as.character(inlift)]
#row.names(lift_data) = gsub(pattern = "[.]",replacement = "-",row.names(lift_data))
row.names(lift_data) = row.names(data2)
right_data = m[as.character(inright)]
#row.names(right_data) = gsub(pattern = "[.]",replacement = "-",row.names(right_data))
row.names(right_data) = row.names(data2)
data.dudi$li = cbind(lift_data,right_data)
data.dudi$li = data.dudi$li[sort(row.names(data.dudi$li)),]
num1 = substring(as.character(inlift),5,6)
num2 = substring(as.character(inright),5,6)
num1_data = n[paste("CS",num1,sep = '')]
num2_data = n[paste("CS",num2,sep = '')]
data.dudi$c1 = cbind(num1_data,num2_data)

right_data_box= merge(data1,right_data,by = "row.names")[,-1]
lift_data_box = merge(data1,lift_data,by = "row.names")[,-1]


d <- data.dudi$li
eig <- ( data.dudi$eig / sum( data.dudi$eig ) )


bb <- head(data.dudi$c1[order(sqrt((data.dudi$c1[,1])^2+(data.dudi$c1[,2])^2),decreasing=T),],n=6L)[1:6,]
cutoff <- vector_fd(bb,d)*0.4
d2 <- data.frame(X=bb[1:dim(bb)[1],1]*cutoff, Y=bb[1:dim(bb)[1],2]*cutoff, LAB=rownames(bb)[1:dim(bb)[1]])
d2[[3]] <- species(as.character(d2[[3]]))


ggdata <- data.frame(d)
Or=c('dental','saliva','faeces')
Type = factor(Type,Or)
  


p<-ggplot(ggdata) +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = 0,linetype=2) +
  geom_hline(yintercept = 0,linetype=2) +
  geom_point(aes(x=d[,1], y=d[,2],fill=Type),color='black',shape=21,size=1) +
  stat_ellipse(aes(x=d[,1], y=d[,2], fill=Type),size=1, geom="polygon", level=0.8, alpha=0.3) +
  #geom_text_repel(data=d2, aes(x=X, y=Y, label=LAB),
  #                fontface="italic", size=3, check_overlap=TRUE) +
  #geom_segment(data=d2, aes(x=0, y=0, xend=X, yend=Y),
  #             arrow = arrow(length = unit(0.2, "cm")), size=0.8, alpha=0.5)+
  scale_shape_manual(values=1:20)+
  scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
  scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
  guides(color=guide_legend(colnames(data1)[0+1]),
         fill=guide_legend(colnames(data1)[0+1]),
         shape=guide_legend(colnames(data1)[0+1]) ) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position=c(0.9,0.9),
        text=element_text(
          size=8
        ),
        legend.title=element_text(
          size=10
        ),
        legend.text=element_text(
          size=8
        )
  )+xlim(min(d[,1], 0)*xx, max(d[,1])*xx)+ylim(min(d[,2], 0)*yy, max(d[,2])*yy)
#xlim(min(d[,1], 0)*3, max(d[,1])*3)+ylim(min(d[,2], 0)*3, max(d[,2])*3)
p<- ggplotGrob(p)
right_data_box$Group = factor(right_data_box$Group, levels=Or)
d <- ggplot(right_data_box)+geom_boxplot(aes(x = Group,y = right_data_box[,2],fill = Group,color=Group),width = 0.5,outlier.size = 0.8)+
  theme_bw()+theme(panel.grid =element_blank())+
  scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
  scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
  #scale_fill_wsj("colors6", "Group")
  # scale_fill_manual(values=brewer.pal(9,"Set1"),breaks =f)+
  guides(fill=FALSE,color=F)+theme(axis.text.x = element_blank())+
  theme(
    #plot.margin = unit(c(-0.2,-0.2,-0.2,-0.2),"mm"),
    axis.ticks = element_blank(),
        text=element_text(
          size=8
        ))+
  ylim(min(right_data_box[,2], 0)*yy, max(right_data_box[,2])*yy)+
  # ylim(min(right_data_box[,2], 0)*3, max(right_data_box[,2])*3)+
  xlab("")+ylab(paste("PC",num2," (",round(eig[as.numeric(num2)]*100,2),"%)",sep=""))
lift_data_box$Group = factor(lift_data_box$Group, levels=Or)
b<- ggplot(lift_data_box)+geom_boxplot(aes(x = Group,y = lift_data_box[,2],fill = Group,color=Group),width = 0.5,outlier.size = 0.8)+
  theme_bw()+theme(panel.grid =element_blank())+coord_flip()+
  scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
  scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
  guides(fill=FALSE,color=F)+theme(axis.text.y = element_blank())+
  theme(axis.ticks = element_blank(),
        text=element_text(
          size=8
        ))+
  #scale_fill_manual(values=brewer.pal(9,"Set1"))+
  ylim(min(lift_data_box[,2], 0)*xx, max(lift_data_box[,2])*xx)+
  xlab("")+ylab(paste("PC",num1," (",round(eig[as.numeric(num1)]*100,2),"%)",sep=""))
a<-ggplot()+theme_bw()+theme(panel.border = element_blank(),panel.grid =element_blank(),
                             axis.text = element_blank(),axis.title = element_blank(),
                             axis.ticks = element_blank())+
  annotate("text", x=1, y=40, label=paste("Adonis R2 =",round(adonis1[[1]][5][[1]][1],4),'\n',
                                          "P.value   =",round(adonis1[[1]][6][[1]][1],4)), size=3)
a <- ggplotGrob(a)
d <- ggplotGrob(d)
b <- ggplotGrob(b)


grid.arrange(d,p,a,b,ncol=2,widths=c(1.2,4),heights = c(4,1.2))


#科水平组成#################################################################################################
Or = c('dental','saliva','faeces')
prof = read.csv('data/profile/virome.family.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
prof = prof[order(-rowMeans(prof)),map$ID]
prof = prof[,order(prof[1,],prof[2,],prof[3,])]
row.names(prof) = sub('.*f__','',row.names(prof))

dat = melt(as.matrix(prof))
dat = merge(dat,map[,1:2],by.x = 'Var2',by.y = 'ID',all.x = T)
dat$Var1 = factor(dat$Var1,row.names(prof))
dat$Var2 = factor(dat$Var2,colnames(prof))
dat$Group = factor(dat$Group,Or)

ggplot(dat)+
  geom_bar(aes(x=Var2,y=value,fill=Var1),position = 'stack',stat = 'identity')+
  facet_grid(~Group,scales = 'free',space = 'free',switch = 'x')+
  scale_fill_manual(values = mycolor[c(1:3,10,5:9,4,11:30)])+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key.height = unit(0.02, "cm"),
    legend.text=element_text(
      size=6
    )
  )

dd = aggregate(dat$value,dat[,c('Var1','Group')],mean)
dd$Var1 = factor(dd$Var1,row.names(prof))
dd$Group = factor(dd$Group,Or)

ggplot(dd)+
  geom_bar(aes(x=Group,y=x,fill=Var1),position = 'stack',stat = 'identity',width = 0.7)+
  scale_fill_manual(values = mycolor[c(1:3,10,5:9,4,11:30)])+
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
