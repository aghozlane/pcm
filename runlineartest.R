#! /usr/bin/Rscript --vanilla

# classification liblinear with bootstrap only
library(LiblineaR)
library(ROCR)
library(pROC)

#==============================================================
# Arguments
#==============================================================
# Parse aguments
input=function(){
        x=commandArgs(trail=T)
        if(length(x) != 6){
                print("Usage : ./runlineartest.r prediction.csv universal_model.csv number_bootstrap reference_output.pdf ref_output.txt prediction_output.tsv")
                q("yes",status=1)
        }
        return(x)
}

#==============================================================
# Analysis function
#==============================================================
trainset.select=function(dataframe, r) {
   label=names(dataframe)
   label.length=length(label)
   r.label=round(label.length*r)
   label.sample=sample(label, r.label, replace=F)
   return(dataframe[,label.sample])
}


logistic.classification=function(dataframe, groups, r=NULL, nboot=100) {
   require(LiblineaR)
   result.matrix=matrix(nr=ncol(dataframe), nc=nboot)
   rownames(result.matrix)=names(dataframe)
   metadata=data.frame(row.names=names(dataframe), groups)
   if(is.null(r)) r=0.5
   for (i in 1:nboot){
      # select train set and test set
      xTrain=t(trainset.select(dataframe,r))
      yTrain=metadata[rownames(xTrain),"groups"]
      xTest=t(dataframe[,!(names(dataframe) %in% rownames(xTrain))])
      yTest=metadata[rownames(xTest),"groups"]
      # Center and scale data
      #s=scale(xTrain,attr(s,"scaled:center"),attr(s,"scaled:scale"))
      s=log10(xTrain+0.000001)
      # Logistic Regression
      ty=6 #L1 logistic regression
      # Tune the cost parameter of a logistic regression via a 10-fold cross-validation
      #try=c(10000, 3000, 1000,300, 100, 30, 10, 3, 1, 0.3, 0.1)
      #tryTypes=c(0:7)
      try=c(10000,  1000, 100,  10,  1, 0.1)
      res=c()
      #for(ty in tryTypes){
         for(co in try){
            acc=LiblineaR(data=s, target=yTrain,type=ty,cost=co,bias=TRUE,cross=10,verbose=FALSE)
            res=c(res,acc)
            #cat("type",ty,"Results for C=",co," : ",acc," accuracy.\n",sep="")
            #flush.console()
         }
      #}
      # Re-train a model with best C value.
      best=which.max(res)
      co=try[best]
      m=LiblineaR(data=s, target=yTrain, type=ty, cost=co, bias=TRUE,
                  verbose=FALSE)
      # Scale the test data
      #s2=scale(xTest,attr(s,"scaled:center"),attr(s,"scaled:scale"))
      s2=log10(xTest+0.000001)
      # Make prediction
      pred=predict(m,s2,proba=TRUE)
      p=pred$predictions
      result.matrix[rownames(xTest),i]=p
      cat(i,"\r")
      flush.console()
   }
   result.table=apply(result.matrix,1, function(x){x=factor(x, levels=levels(factor(groups))); table(x)})
   result.prob=t(prop.table(result.table,2))
   groups.est=apply(result.prob,1, function(x) names(which(x==max(x))[1]))
   groups.prob=apply(result.prob,1, function(x) x[which(x==max(x))[1]]   )
   result=data.frame(result.prob, groups.input=groups, groups.output=groups.est,
                     confidence=groups.prob)
   return(result)
}


# classification liblinear function with CV, bootstrap and model performance estimator

cost_param=function(data, labels, type, cost,co, cross, verbose=FALSE,
                    cost.try=c(10000,  1000, 100,  10,  1, 0.1)) {
    res=c()
    for(co in cost.try){
        acc=LiblineaR(data=data, target=labels, type=type, cost=co, bias=TRUE,
                      cross=10,verbose=FALSE)
        res=c(res,acc)
    }
    best=which.max(res)
    co=cost.try[best]
    return(co) #best cost
}


create_folds=function(data, nfolds) {
    perm = sample(1:nrow(data), nrow(data)) / nrow(data)
    id = rep(0, nrow(data))
    for (f in nfolds:1) {
        id[perm <= f/nfolds] = f
    }
    return(id)
}


logistic.CV=function(data, labels, nfolds, ty, cost.try, sparse=FALSE) {

    if(length(levels(as.factor(labels))) != 2) stop("length of labels is not equal 2")
    list=1:nfolds
    id = create_folds(data, nfolds)
    #actual.id=1:nrow(data)
    prediction = NULL
    probabilities = NULL
    actual = NULL
    weight = vector("list", nfolds)
    model = vector("list", nfolds)
    names = NULL
    for (i in 1:nfolds) {
        xTrain = subset(data, id %in% list[-i])
        if(!sparse) {
            xTrain = scale(xTrain,center=TRUE,scale=TRUE) #scale the data
        } else {
            xTrain = t(t(scale(xTrain, scale=FALSE))/(apply(xTrain,2,sd)+quantile(apply(xTrain,2,sd),0.1))) # add 10th sd
            xTrain[,attr(na.omit(t(xTrain)), "na.action")] = 0
        }
        yTrain = subset(labels, id %in% list[-i])
        xTest = subset(data, id %in% c(i))
        if(!sparse) {
            xTest = scale(xTest,center=TRUE,scale=TRUE) # independently scale the test set
        } else {
            xTest = t(t(scale(xTest, scale=FALSE))/(apply(xTest,2,sd)+quantile(apply(xTest,2,sd),0.1))) # add 10th sd
            xTest[,attr(na.omit(t(xTest)), "na.action")] = 0
        }
        yTest = subset(labels, id %in% c(i))
        #print(dim(xTest))
        # Tune the cost parameter of a logistic regression via a nested 10-fold cross-validation
        co=cost_param(data=xTrain, labels=yTrain, type=ty,cost=co, cross=10,verbose=FALSE, cost.try)
        #print(co)
        # Re-train a model with best C value.
        mymodel=LiblineaR(data=xTrain,target=yTrain,type=ty,cost=co,bias=TRUE,verbose=FALSE)
        weight[[i]]=mymodel$W
        model[[i]]=mymodel
        pred=predict(mymodel,xTest,proba=TRUE)
        p=pred$predictions
        prob=pred$probabilities[,1]

        prediction = c(prediction,p)
        probabilities = c(probabilities,prob)
        actual = c(actual, as.character(yTest))
        names=c(names,subset(rownames(data), id %in% c(i)))
   }
   pred=data.frame(prediction=prediction,actual=actual,probabilities=probabilities,
                   names=names)
   W=sapply(weight, rbind, simplify=TRUE)
   row.names(W)=colnames(weight[[1]])
   return(list(pred=pred, W=W, model=model))
}


logistic.CV.boot=function(data=dd, labels=labels, nfolds=10, ty=6, cost.try=c(100000, 10000,  1000, 100,  10,  1, 0.1,0.01,0.001), sparse=TRUE, nboot=10, p.bar=TRUE){
   res=vector("list",nboot)
   if(p.bar) {pb <- txtProgressBar(min = 0, max = nboot, style = 3)}
   for (b in 1:nboot) {
   #do logistic CV
      res[[b]]=logistic.CV(data=data, labels=labels, nfolds=10, ty=6, cost.try=c(100000, 10000,  1000, 100,  10,  1, 0.1,0.01,0.001), sparse=sparse)
      setTxtProgressBar(pb, b)
   }
   if(p.bar) {close(pb)}
   return(res)
}


model.performance=function(res, nboot=10){
   perf=vector("list",nboot)
   model.W=vector("list",nboot)
   for(b in 1:nboot) {
   # calculate performance
      #levels(res[[b]]$pred$actual) = c(1,-1)
      pre=prediction(prediction=res[[b]]$pred$probabilities+(rnorm(length(res[[b]]$pred$probabilities))/(10^9)), labels=res[[b]]$pred$actual) #add some noise to fix a bug
      per=performance(pre,"tpr","fpr")
      AUC=(performance(pre,"auc"))@y.values[[1]]
      fpr=per@"x.values"[[1]]
      tpr=per@"y.values"[[1]]
      perf[[b]]=list(AUC=AUC,fpr=fpr,tpr=tpr)
      model.W[[b]]=res[[b]]$W
   }
   AUC=unlist(lapply(perf, function(x){x$AUC}))
   tpr=sapply(perf, function(x){x$tpr})
   fpr=sapply(perf, function(x){x$fpr})
   w=apply(sapply(model.W, function(x) { apply(x, 1, mean) }), 1, mean)
   w=sapply(model.W, function(x) { apply(x, 1, mean) })
   w=t(w/sum(abs(w)))*nboot
   model.perf=list(AUC=AUC, tpr=tpr, fpr=fpr, features.w=w)
   return(model.perf)

}



model.ext.validation=function(res.ext, nboot=10){
   perf=vector("list",nboot)
   for(b in 1:nboot) {
   # calculate performance
      pre=prediction(prediction=res.ext[[b]]$pred$probabilities+(rnorm(length(res.ext[[b]]$pred$probabilities))/(10^9)), labels=res.ext[[b]]$pred$actual) #add some noise to fix a bug
      per=performance(pre,"tpr","fpr")
      AUC=(performance(pre,"auc"))@y.values[[1]]
      fpr=per@"x.values"[[1]]
      tpr=per@"y.values"[[1]]
      perf[[b]]=list(AUC=AUC,fpr=fpr,tpr=tpr)
   }
   AUC=unlist(lapply(perf, function(x){x$AUC}))
    tpr=sapply(perf, function(x){x$tpr})
    fpr=sapply(perf, function(x){x$fpr})
    model.perf=list(AUC=AUC, tpr=tpr, fpr=fpr)
    return(model.perf)

}



##########################
# Specific to bla data #
##########################

# Arguments input
x=input()
# Read data
#bla=read.csv2(x[1], header=T, row=1, )
data_candidates=read.csv(x[1], header=T, row=1, sep="\t")
# print(bla)
# Split data
#blaref=bla[which(bla$Type=="Ref" | bla$Type=="tneg"),]
blaref=read.csv2(x[2], header=T, row=1, sep=";")
# print(blaref)
lab=blaref$Type
lab=as.factor(as.character(lab))
print(levels(lab))
type_prot=c(levels(lab))
levels(lab)=c("-1","1")
## Herre
#blaref_dd=blaref[,c(3,6:12)]
#blaref_dd=blaref[,c(6:12,14,25)]
blaref_dd=blaref[,c(13:34)]
#blaref_dd=blaref[,c(5,7:13)]
#rownames(blaref_dd) = blaref[,1]
rownames(blaref_dd) = rownames(blaref)
# Run logistic analysis
nboot=as.integer(x[3])
res=logistic.CV.boot(data=blaref_dd, labels=lab, nboot=nboot)
model.perf=model.performance(res, nboot)
tpr=model.perf$tpr
fpr=model.perf$fpr
AUC=model.perf$AUC
w=model.perf$features.w[,-which(colnames(model.perf$features.w)=="Bias")]

# Prediction for the reference
result_prediction=matrix(0, dim(blaref)[1], 3)
## Herre
rownames(result_prediction) = rownames(blaref)
#rownames(result_prediction) = blaref[,1]
for(b in 1:nboot) {
        for (i in seq(1, dim(blaref_dd)[1])){
            a = which(rownames(result_prediction) == as.character(res[[b]]$pred$names[i]))
            if(res[[b]]$pred$prediction[i] == "1"){
                result_prediction[a, 2] = result_prediction[a, 2] + 1
            } else {
                result_prediction[a, 1] = result_prediction[a, 1] + 1
            }
        }
}
result_prediction[,1] = 100.0 * result_prediction[,1] / nboot
result_prediction[,2] = 100.0 * result_prediction[,2] / nboot



# Plot result
pdf(x[4])
    # Roc curve on of each bootstrap
    plot(apply(fpr,1,median), apply(tpr, 1, median),main="ROC Curve",xlab="False Positive Rate",ylab="True Positive Rate",sub=paste("AUC:",round(median(AUC),2)," bootstrap:",nboot), type="n")
    # #plot(apply(fpr,1,median), apply(tpr, 1, median),main="ROC Curve",ylab="True Positive Rate",xlab="False Positive Rate",sub=paste("AUC:",round(median(AUC),2)), type="n")
    lines(apply(fpr,1,function(x) {quantile(x,0.1)}), apply(tpr,1,function(x) {quantile(x,0.9)}), lty=2)
    lines(apply(fpr,1,median), apply(tpr,1,median), lwd=2)
    lines(apply(fpr,1,function(x) {quantile(x,0.9)}), apply(tpr,1,function(x) {quantile(x,0.1)}), lty=2)
    par(mar=c(5,10,3,3))
    boxplot(w, las=2, horizontal=TRUE, xlab="% total model weight", main="feature weight")
    # Roc curve on the mean results of each bootstrap
    roc_data=roc(response=as.character(blaref$Type), predictor=result_prediction[,1], plot=TRUE,print.auc=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, auc.polygon.col="gray",  print.thres=TRUE,reuse.auc=FALSE, ci = TRUE, boot.n=2000)
    temp=coords(roc_data, "b", ret=c("se","sp","t"), best.method="youden")
    #temp=coords(roc_data, "b", ret=c("se","sp","t"), best.method="closest.topleft",best.weights=c(0.1, 0.2)
    points(temp[2],temp[1],col="black",pch=16)
dev.off()

# Write reference result
result_prediction[,3] = as.character(blaref$Type)
colnames(result_prediction) = c(type_prot, "Actual")
write.table(result_prediction, file=x[5], sep="\t")

#TODO : Model interpretation and marker extraction (50% of models, logs mean of weight in total weight of model)

# Candidate
##Herre

#data_type=levels(bla[,1])
# Remove ref and tneg
#data_candidates = data_type[which(!(data_type %in% c("Ref", "tneg")))]
# for(dtyp in data_candidates){
#     ##Herre
#     #blacandidate=bla[bla$Type=="Candidate",]
#     blacandidate=bla[bla$Type==dtyp,]
#     ## Herre
#     #blacandidate_dd=blacandidate[,c(4,6:12)]
#     blacandidate_dd=blacandidate[,c(5,7:13)]
#blacandidate_dd = data_candidates[,c(3,6:12)]
#blacandidate_dd = data_candidates[,c(6:12,14,25)]
blacandidate_dd = data_candidates[,c(13:34)]
    #predict pour le model i
    res_predict=vector("list",nboot)
    result_prediction=matrix(0, dim(blacandidate_dd)[1], 3)
    for(b in 1:nboot) {
        res_predict[[b]]=lapply(res[[b]]$model,function(x) {predict(x, blacandidate_dd, prob=TRUE)})
        for (pred in res_predict[[b]]){
            #print(pred$probabilities)
            for (i in seq(1, dim(blacandidate_dd)[1])){
                if(pred$predictions[i] == "1"){
                    result_prediction[i, 2] = result_prediction[i, 2] + 1
                } else {
                    result_prediction[i, 1] = result_prediction[i, 1] + 1
                }
            }
        }
    }


    for (i in seq(1, dim(blacandidate_dd)[1])){
        result_prediction[i, 3] = prop.test(result_prediction[i, 1], 10 * nboot, p=0.5, alternative='greater')$p.value
        result_prediction[i, 1] = round(100.0 * result_prediction[i, 1] / (10.0 * nboot),1)
        result_prediction[i, 2] = round(100.0 * result_prediction[i, 2] / (10.0 * nboot),1)
    }
    
    #correction 
    result_prediction[i, 3] = p.adjust(result_prediction[i, 3],method="BH")
    #print(blacandidate)
    ##Herre
    rownames(result_prediction) = rownames(blacandidate_dd)
    #rownames(result_prediction) = blacandidate[,1]
    colnames(result_prediction) = c(type_prot, paste("Adjusted p.value",type_prot[1]))
    #print(cbind(blacandidate[,1], result_prediction))
    #write.table(result_prediction, file=paste(x[6],dtyp,".txt",sep="") , sep="\t")
    write.table(result_prediction, file=paste(x[6], sep="") , sep="\t")
# }

