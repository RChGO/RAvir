#high-abundant marker (randomForest)##########################################################################################################
cutoff = 0.0005
Prof = read.csv('data/profile/metaphlan.L7.prof',sep='\t',stringsAsFactors = F,check.names = F,row.names = 1)
Map = read.csv('data/profile/mapping.file',sep='\t',stringsAsFactors = F,check.names = F)
Prof = Prof[order(-rowSums(Prof)),Map$ID]
Prof = Prof[rowMeans(Prof)>cutoff,]


tt = c()
Or1=c('dental','saliva','faeces')
for(x in Or1){
  m=Map[Map$Group==x,]
  tt[[x]] = kw_test(Prof,m,'Type')
}



library(pROC)
inc = 'Type'
intype = 'Discrete'
Threads=1
lab='saliva'
gg=list(c('Control','Onset'),c('Control','Treated'),c('Onset','Treated'))

Random = function(prof,Sam,intype,y){
  library(randomForest)
  Sample = Sam[Sam$group==y,]
  Pred=c()
  if(intype=="Discrete"){
    lab = as.character(prof$Group)[1]
    for(x in 1:nrow(Sample)){
      train = prof[-Sample$order[x],]
      test = prof[Sample$order[x],]
      fit <- randomForest(Group~ .,data=train,ntree=2000,
                          proximity=TRUE,importance=T) 
      test_true <- test$Group
      test <- subset(test,select = -c(Group))
      test_pred <- predict(fit,type="prob",newdata = test)
      pred <- data.frame(ID=row.names(test),true=as.character(test_true),
                         pred=test_pred[,lab],stringsAsFactors = F)
      Pred = rbind(Pred,pred)
      #print(pred)
    }
  }else if(intype=='Continuous'){
    for(x in 1:nrow(Sample)){
      train = prof[-Sample$order[x],]
      test = prof[Sample$order[x],]
      fit <- randomForest(Group~ .,data=train,ntree=2000,
                          proximity=TRUE,importance=T) 
      test_true <- test$Group
      test <- subset(test,select = -c(Group))
      test_pred <- predict(fit,newdata = test)
      test_pred <- as.numeric(test_pred)
      pred <- data.frame(ID=row.names(test),true=test_true,pred=test_pred,stringsAsFactors = F)
      Pred = rbind(Pred,pred)
      #print(pred)
    }
  }
  return(Pred)
}
func = function(y){
  Pred = Random(prof,Sam,intype,y)
  return(Pred)
}

all = c();Fig=c();auc=c()
for(i in 1:length(gg)){
  map <- Map[Map$Group==lab & Map[,inc]%in%gg[[i]],]
  map = data.frame(ID=map[,1],Group=map[,inc],check.names = F,stringsAsFactors = F)
  
  prof = Prof[row.names(Prof)%in%tt[[lab]]$ID[which(tt[[lab]]$P_value<0.05)],map$ID]
  #prof = Prof[,map$ID]
  
  map = map[!is.na(map[,'Group']),]
  #map= map[c(1:10,1000:1010),]
  map = map[order(map$Group),]
  prof= as.data.frame(t(prof[,map$ID]))
  colnames(prof) = make.names(colnames(prof))
  
  
  if(intype=='Discrete'){
    map[,'Group']<- as.factor(map[,'Group'])
    prof[,'Group'] = map[,'Group']
    
  }else if(intype=='Continuous'){
    map[,'Group']<- as.numeric(map[,'Group'])
    prof[,'Group'] = map[,'Group']
  }
  
  if(Threads>1){
    Sam = data.frame(order=1:nrow(prof),group=as.numeric(cut(1:nrow(prof),breaks = Threads)))
    cl <- makeCluster(Threads) 
    clusterExport(cl, c("prof","Sam","intype","Random")) 
    results <- parLapply(cl,1:Threads,func) 
    all <- do.call('rbind',results) 
    #colnames(pred) <- c('True','Predict')
    stopCluster(cl) 
  }else{
    Sam = data.frame(order=1:nrow(prof),group=1)
    all[[i]] <- Random(prof,Sam,intype,1)
  }
  
  
  #roc=
  myroc <- roc(all[[i]]$true,all[[i]]$pred,percent=TRUE,print.auc=TRUE,auc.polygon=TRUE,auc.polygon.col=NA,plot = F)
  auc[[i]]=as.character(lapply(as.numeric(ci.auc(myroc)), function(x){round(x,2)}))
  #plot.roc(myroc)
  #text(20,20,paste('AUC: ',auc[2],'%','\n(95% CI: ',auc[1],'% ~ ',auc[3],'%)',sep = ''))
  
  
  plot1 <- plot.roc(myroc,col='red')
  Fig[[i]] <- data.frame(specificities=plot1$specificities,sensitivities=plot1$sensitivities)

}

color_transparent <- adjustcolor(c("darkred", "darkgreen","darkblue"), alpha.f = 0.8) 

#pdf(paste(lab,'.bt','.pdf',sep=''),width = 4,height = 4)
for(i in 1:length(gg)){
  if(i==1){plot.roc(myroc,col='white')}
  lines(Fig[[i]],col=color_transparent[i],lwd=3,pch=0.5)
  text(32,45-13*i,paste('AUC: ',auc[[i]][2],'%','\n(95% CI: ',auc[[i]][1],'% ~ ',auc[[i]][3],'%)',sep = ''),
       col = color_transparent[i],cex =0.833)
}
#dev.off()
