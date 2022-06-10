
#adonis (RA indices)############################################################################################################
prof = read.csv('analysis/2.beta_div/compare/all.adonis2',sep='\t',header = F,as.is = T)
prof$R2adj = RsquareAdj(prof$V3,n = prof$V5,m=1)
aa = do.call(rbind,lapply(prof$V1, function(x){unlist(strsplit(x,'[.]'))}))
colnames(aa) = c('Micro','level','Group','ID')
pp = prof[,3:6]
colnames(pp) = c('R2','pvalue','sample','R2adj')
pp =cbind(aa,pp)

map = read.csv('data/mapping.file/medicine.map',sep='\t',header = F,as.is = T) 
pp = merge(pp,map,by.x = 'ID',by.y = 'V1')
pp$level[pp$level=='species'] = 'species / vOTUs'
pp$level[pp$level=='vOTUs'] = 'species / vOTUs'
aa = split(pp,pp[,c('Micro','level','Group')],drop = T)
for(x in 1:length(aa)){
  aa[[x]]$q=p.adjust(aa[[x]]$pvalue,method = 'BH')
}
pp = do.call(rbind,aa)

#ppf = pp[pp$pvalue<0.05,]
#pp = pp[pp$ID%in%unique(ppf$ID),]

gg = unique(pp[,c('ID','V3')])
gg$V4 = as.numeric(factor(gg$V3,c('Demographic data','RA-associated indice','Drugs','Others')))
gg = gg[order(gg$V4),]


fig = c()
for(x in unique(gg$ID)){
  dat = pp[pp$ID==x,]
  dat$level  = factor(dat$level,c('species / vOTUs','family'))
  dat$Micro[dat$Micro=='vs'] = 'Virus'
  dat$Micro[dat$Micro=='bt'] = 'Archaea/Bacteria'
  dat$Micro  = factor(dat$Micro,c('Virus','Archaea/Bacteria'))
  dat$lab = ''
  dat$lab[dat$pvalue<0.05] = '+'
  dat$lab[dat$q<0.05] = '*'
  dat$labp = dat$R2adj+max(dat$R2adj)/6
  dat$labp[dat$R2adj<0] = dat$R2adj-max(dat$R2adj)/6
  dat$Group = factor(dat$Group,c('dental','saliva','faecal'))
  
  fig[[x]] = ggplot(dat)+
    geom_bar(aes(x=Group,y=R2adj,fill=Group,color=Group),position = 'dodge',stat = 'identity',width = 0.7)+
    geom_text(aes(x=Group,y=labp,label=lab),position = 'dodge',stat = 'identity')+
    facet_grid(level~Micro)+
    scale_color_manual(values = c("darkred", "darkgreen","darkblue"))+
    scale_fill_manual(values = c("#f76f6f", "#7faa81", "#9dbcd6"))+
    labs(x=x,y='adj-R square')+
    theme_bw()+
    guides(fill=F,color=F)+
    theme(
      panel.grid = element_blank(),
      panel.spacing = unit(0, "lines"),
      text = element_text(
        size = 5
      ),
      axis.text = element_text(
        angle = 45,
        hjust = 1,
        size = 5
      ),
      strip.text = element_text(
        size=5,
        margin = margin(0.08,0.08,0.08,0.08, "cm")
        )
    )
}

cowplot::plot_grid(plotlist = fig)

#network (indices and microbes)#########################################################################################
gg = c('CDAI','CRP','DAS28','Disease_activity','GH','RA_duration')
Or=c('dental','saliva','faeces')
p1 = read.csv('data/profile/virome.vOTU.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
p2 = read.csv('data/profile/metaphlan.L7.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
tax = read.csv('analysis/4.abundance.sig/vOTU_L7.rb_gt_0.0005.txt',sep='\t',stringsAsFactors = F,check.names = F)
prof = rbind(p1[,map$ID],p2[,map$ID])
prof = prof[tax$ID,]
mtax = data.frame(ID=c(row.names(p1),row.names(p2)),Level=rep(c('Virus_votu','bt_species'),c(nrow(p1),nrow(p2))),stringsAsFactors = F)

p1 = read.csv('data/profile/virome.family.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
p2 = read.csv('data/profile/metaphlan.L5.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
prof = rbind(prof,p1[rowMeans(p1)>0.0005,map$ID],p2[rowMeans(p2)>0.0005,map$ID])
mtax = rbind(mtax,data.frame(ID=c(row.names(p1),row.names(p2)),Level=rep(c('Virus_family','bt_family'),c(nrow(p1),nrow(p2))),stringsAsFactors = F))

med = read.csv('data/mapping.file/medicine.txt',sep='\t',row.names = 1,stringsAsFactors = F,check.names = F)
med = med[map$ID,gg]

dat = c()
for(x in Or){
  mm = map[map$Group==x,]
  prof_f = t(prof[,mm$ID])
  med_f = med[mm$ID,]
  aa = cor_v2(prof_f,med_f)
  aa = two_matr_rp(aa$r,aa$p)
  aa$lab = x
  aa$q = NA
  aa$q[!is.na(aa$pvalue)] = p.adjust(aa$pvalue[!is.na(aa$pvalue)],method = 'BH')
  aa$node_size = colMeans(prof_f)
  dat = rbind(dat,aa)
}
dat = dat[which(dat$pvalue<0.01),]
#dat = dat[which(dat$q<0.05),]
dat = merge(dat,tax,by='ID',all.x = T)
la = sub('.*[|;]f__','',dat$ID)
la = sub('[|;]g__.*','',la)
dat$Rename[is.na(dat$Rename)] = la[is.na(dat$Rename)]
dat$Family[is.na(dat$Family)] = la[is.na(dat$Family)]
dat = dat[,1:9]
dat$linetype=1
dat$linetype[dat$q<0.05]=2
dat$edge_size = abs(dat$value)
dat$edge_color = 1
dat$edge_color[dat$value<0] =2
dat = merge(dat,mtax,by = 'ID',all.x = T)
dat$ID = paste(dat$lab,dat$Rename,sep='_')
Edge = dat[,c('ID','variable','linetype','edge_color','edge_size')] 
AA=aggregate(Edge$edge_size,list(ID=Edge$ID),sum)
colnames(AA)[2]  = 'node_size2' 

Node = unique(dat[,c('ID','Rename','node_size','lab','Level')])
aa = data.frame(ID=unique(dat$variable),Rename=unique(dat$variable),node_size=0.7,lab='clinical',Level='clinical',stringsAsFactors = F)
Node = rbind(Node,aa)
Node = merge(Node,AA,by.x = 'ID',all.x = T)
Node$node_size2[is.na(Node$node_size2)] = 4

write.table(Edge,'analysis/2.beta_div/compare/adonis/edge3.txt',sep='\t',row.names = F,quote = F)
write.table(Node,'analysis/2.beta_div/compare/adonis/node3.txt',sep='\t',row.names = F,quote = F)






