#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[3] Proof-of-concept (Curtis data)  --- Oncotype DX
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.1] discovery --> validaiton 
rm(list=ls())

#clinical information for each sample
myinf1A = "D:/Gene-Signature-Predictability-Prediction/discovery_Patient_tumor_info.txt"
myinf1B = "D:/Gene-Signature-Predictability-Prediction/validation_Patient_tumor_info.txt"

#oncotype DX scores for each sample
myinf2 = "D:/Gene-Signature-Predictability-Prediction/Curtis_oncotypeDX_result.txt"

#define function that will help us with plotting graphs later on
resize.win <- function(Width=5, Height=5)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}

#------------------------------
#Read and format data
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
info.A <- read.table(myinf1A, header=T, sep="\t",  quote="")
info.B <- read.table(myinf1B, header=T, sep="\t",  quote="")
source = c(rep("discovery", nrow(info.A)), rep("validation", nrow(info.B)))
info = rbind(info.A, info.B)
info = cbind(info, source)
info = info[!is.na(info$last_follow_up_status),]
info = info[!is.na(info$T),]
info[,1] = as.character(info[,1])
xx = info[,1]
xx = gsub("-", ".", xx)
info[,1] = xx
row.names(info) = xx
xx = rep(0, nrow(info))
xx[info[, "last_follow_up_status"]=="d-d.s."] = 1
e.rfs = xx
t.rfs = as.numeric(info[, "T"])
info = cbind(t.rfs, e.rfs, info)

#------------------------------
#Find common samples between data and clinical information
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)

#convert lymph_nodes_positive to binary variable
xx = data$lymph_nodes_positive
xx = ifelse(xx>0, 1, 0)
data$lymph_nodes_positive = xx

#remove data columns that have more than 500 na's
xx = apply(is.na(data), 2, sum)
data = data[, xx<500]

#add columns describing whether the sample had received chemotherapy (CT), radiation therapy (RT), or hormone therapy (HT)
se = grep("CT", data$Treatment)
CT = rep("no", nrow(data))
CT[se] = "yes"
se = grep("RT", data$Treatment)
RT = rep("no", nrow(data))
RT[se] = "yes"
se = grep("HT", data$Treatment)
HT = rep("no", nrow(data))
HT[se] = "yes"
data = cbind(CT, RT, HT, data)
se = c("CT", "RT","HT", "mytf","t.rfs", "e.rfs", "age_at_diagnosis", "menopausal_status_inferred", "grade","size", "stage","lymph_nodes_positive", "NPI", "histological_type","ER_IHC_status", "cellularity", "Pam50Subtype", "ER.Expr", "Her2.Expr","PR.Expr", "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="discovery",]

#Fit a Cox Proportional Hazards model using the OncoType DX score to predict survival
mycox = coxph(Surv(t.rfs, e.rfs)~mytf, data) 
summary(mycox)$concordance
## C: 0.67888999

#combine results into xx to calculate concordance
xx = cbind(data$t.rfs, mycox$linear.predictors, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

#calculate concordance by selecting all valid pairs and counting the number of concordant pairs
cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
    if(xx[i,"event"]==0 & xx[j,"event"]==0) next
    x1 = xx[i,"time"]-xx[j,"time"]
    x2 = xx[i,"score"]-xx[j,"score"]
    
    if(xx[i,"event"]==1 & xx[j,"event"]==1)
    {
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x1*x2>0)  
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x1*x2<0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
      next
    }
    if(xx[i,"event"]==1 & xx[j,"event"]==0)
    {
      if(x1>0)  next
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x2<0) 
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x2>0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
      next
    }
    if(xx[i,"event"]==0 & xx[j,"event"]==1)
    {
      if(x1<0)  next
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x2>0) 
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x2<0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
    }
  }

#----------------------------------------------
cc #112436
dd #52671
mycox$concordance
#Output---------------------------------------
  concordant   discordant       tied.x       tied.y      tied.xy  concordance 
1.124360e+05 5.267100e+04 1.937000e+03 8.000000e+00 0.000000e+00 6.788900e-01 
std 
1.722481e-02
#---------------------------------------------
#Calculate accuracy and create a plot to show that the accuracy vary dramatically between samples
resize.win(5,5) #skip this line of code if you are not using windows, or if it gives a error
res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

plot(count, accuracy, pch=20)
abline(v=50, col=2)

#only keep the samples that have at least 20% of comparisons being effective
se = which(res[, "count"]>nrow(data)*0.2)
res = res[se,]


#----------------------------------------------
#Divide data into well predicted and poorly predicted group
library(randomForest)
data = raw.data
data = data[data$source=="discovery",]
comxx = intersect(row.names(data), row.names(res))
acc = res[comxx, "accuracy"]
thr = summary(mycox)$concordance[1] #choose the threshold to be 0.679, the overall concordance of the model, as stated in the manuscript
plot(density(acc))
abline(v=thr)
sum(acc>=thr)
class = ifelse(acc>=thr, "Y", "N")
data = cbind(class, data[comxx,])
data=rfImpute(class~., data=data)
data = cbind(acc, data)
#----------------------------------------------
# Predictability prediction model (PPM) on poor and well predicted samples
mod1 = randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold Cross Validation

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)


myres = NULL
for(k in 1:kk)
{
  cat("\r", k)
  if(k==kk)
  {
    se = ((k-1)*psiz+1):nrow(pdat)
  }else
  {
    se = ((k-1)*psiz+1):(k*psiz)	
  }
  ptr = pdat[-se,]	
  pte = pdat[se,]
  if(k==kk)
  {
    se = ((k-1)*nsiz+1):nrow(ndat)
  }else
  {
    se = ((k-1)*nsiz+1):(k*nsiz)	
  }
  ntr = ndat[-se,]	
  nte = ndat[se,]
  
  
  tr = rbind(ptr, ntr)
  te = rbind(pte, nte)
  fit <- randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=tr, ntree=10000)
  tmp = predict(fit, te[,-1], type="prob")
  tmp = data.frame(te[,"class"], tmp)
  myres = rbind(myres, tmp)
}


thr = (1:99)*0.01
yy =  xx =  rep(0, length(thr))
fdr = rep(0,99)
for(i in 1:length(thr))
{
  aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
  bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
  cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
  dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
  fdr[i] = aa/sum(myres[,2]>=thr[i])
  yy[i] = aa/(aa+bb)
  xx[i] = cc/(cc+dd)
}
xx = c(1, xx, 0)
yy = c(1, yy, 0)
tmp1 = tmp2 = rep(0,100)
for(i in 1:100)
{
  tmp1[i] = xx[i]-xx[i+1]
  tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)
##  AUC = 0.6582878

#----------------------------------------------
#------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[data$Pam50Subtype!="NC", ]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(mycox, newdata= data, na.action=na.exclude)
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx, ]
risk = risk[comxx]

##
xx = cbind(data$t.rfs, risk, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
    if(xx[i,"event"]==0 & xx[j,"event"]==0) next
    x1 = xx[i,"time"]-xx[j,"time"]
    x2 = xx[i,"score"]-xx[j,"score"]
    
    if(xx[i,"event"]==1 & xx[j,"event"]==1)
    {
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x1*x2>0)  
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x1*x2<0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
      next
    }
    if(xx[i,"event"]==1 & xx[j,"event"]==0)
    {
      if(x1>0)  next
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x2<0) 
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x2>0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
      next
    }
    if(xx[i,"event"]==0 & xx[j,"event"]==1)
    {
      if(x1<0)  next
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x2>0) 
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x2<0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
    }
  }

#----------------------------------------------
cc
dd

count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
cv.res = cbind(predictability, f.mat, count, accuracy)
cv.res = cv.res[order(cv.res[,"Y"], decreasing=T), ]
plot(cv.res[,"Y"], cv.res[, "accuracy"])
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs")
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs", method="s")
## positive correlation 


#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]

w.siz = 200
step = 20
cnum = ceiling((nrow(tmp)-w.siz)/step)
yy1 = rep(0, cnum)
for(k in 1:cnum)
{
  sta = (k-1)*step+1
  end = min(sta+w.siz-1, length(yy1))
  tmpxx = tmp[sta:end,]
  yy1[k] = (sum(tmpxx[, "concordance"]) + 0.5*sum(tmpxx[, "tie"]))/sum(tmpxx[, "count"])
}
plot(yy1, pch=20, type="b")


#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]			## check here
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
  tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
  sta = (k-1)*step+1
  end = min(sta+w.siz-1, length(cum.con))
  yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")
abline(h=cum.con[length(cum.con)], col=2)		## smoothed curves

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.2] validation  --> discovery
rm(list=ls())
myinf1A = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/discovery_Patient_tumor_info.txt"
myinf1B = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/validation_Patient_tumor_info.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/Curtis_oncotypeDX_result.txt"


data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info.A <- read.table(myinf1A, header=T, sep="\t",  quote="")
info.B <- read.table(myinf1B, header=T, sep="\t",  quote="")
source = c(rep("discovery", nrow(info.A)), rep("validation", nrow(info.B)))
info = rbind(info.A, info.B)
info = cbind(info, source)
info = info[!is.na(info$last_follow_up_status),]
info = info[!is.na(info$T),]
info[,1] = as.character(info[,1])
xx = info[,1]
xx = gsub("-", ".", xx)
info[,1] = xx
row.names(info) = xx
xx = rep(0, nrow(info))
xx[info[, "last_follow_up_status"]=="d-d.s."] = 1
e.rfs = xx
t.rfs = as.numeric(info[, "T"])
info = cbind(t.rfs, e.rfs, info)

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)

xx = data$lymph_nodes_positive
xx = ifelse(xx>0, 1, 0)
data$lymph_nodes_positive = xx
xx = apply(is.na(data), 2, sum)
data = data[, xx<500]
se = grep("CT", data$Treatment)
CT = rep("no", nrow(data))
CT[se] = "yes"
se = grep("RT", data$Treatment)
RT = rep("no", nrow(data))
RT[se] = "yes"
se = grep("HT", data$Treatment)
HT = rep("no", nrow(data))
HT[se] = "yes"
data = cbind(CT, RT, HT, data)
se = c("CT", "RT","HT", "mytf","t.rfs", "e.rfs", "age_at_diagnosis", "menopausal_status_inferred", "grade","size", "stage","lymph_nodes_positive", "NPI", "histological_type","ER_IHC_status", "cellularity", "Pam50Subtype", "ER.Expr", "Her2.Expr","PR.Expr", "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="validation",]

mycox = coxph(Surv(t.rfs, e.rfs)~mytf, data) 
summary(mycox)$concordance
## 0.64135400

##
xx = cbind(data$t.rfs, mycox$linear.predictors, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
    if(xx[i,"event"]==0 & xx[j,"event"]==0) next
    x1 = xx[i,"time"]-xx[j,"time"]
    x2 = xx[i,"score"]-xx[j,"score"]
    
    if(xx[i,"event"]==1 & xx[j,"event"]==1)
    {
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x1*x2>0)  
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x1*x2<0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
      next
    }
    if(xx[i,"event"]==1 & xx[j,"event"]==0)
    {
      if(x1>0)  next
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x2<0) 
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x2>0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
      next
    }
    if(xx[i,"event"]==0 & xx[j,"event"]==1)
    {
      if(x1<0)  next
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x2>0) 
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x2<0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
    }
  }

#----------------------------------------------
cc
dd
mycox$concordance

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

plot(count, accuracy, pch=20)			#  this figure indicates the accuracy vary dramatically between samples

se = which(res[, "count"]>nrow(data)*0.2)
res = res[se,]


#----------------------------------------------
data = raw.data
data = data[data$source=="validation",]
comxx = intersect(row.names(data), row.names(res))
acc = res[comxx, "accuracy"]
thr = summary(mycox)$concordance[1]
plot(density(acc))
abline(v=thr)
sum(acc>=thr)
class = ifelse(acc>=thr, "Y", "N")
data = cbind(class, data[comxx,])
data=rfImpute(class~., data=data)
data = cbind(acc, data)
#----------------------------------------------
mod1 = randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)
## 0.7373272

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)


myres = NULL
for(k in 1:kk)
{
  cat("\r", k)
  if(k==kk)
  {
    se = ((k-1)*psiz+1):nrow(pdat)
  }else
  {
    se = ((k-1)*psiz+1):(k*psiz)	
  }
  ptr = pdat[-se,]	
  pte = pdat[se,]
  if(k==kk)
  {
    se = ((k-1)*nsiz+1):nrow(ndat)
  }else
  {
    se = ((k-1)*nsiz+1):(k*nsiz)	
  }
  ntr = ndat[-se,]	
  nte = ndat[se,]
  
  
  tr = rbind(ptr, ntr)
  te = rbind(pte, nte)
  fit <- randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=tr, ntree=10000)
  tmp = predict(fit, te[,-1], type="prob")
  tmp = data.frame(te[,"class"], tmp)
  myres = rbind(myres, tmp)
}


thr = (1:99)*0.01
yy =  xx =  rep(0, length(thr))
fdr = rep(0,99)
for(i in 1:length(thr))
{
  aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
  bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
  cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
  dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
  fdr[i] = aa/sum(myres[,2]>=thr[i])
  yy[i] = aa/(aa+bb)
  xx[i] = cc/(cc+dd)
}
xx = c(1, xx, 0)
yy = c(1, yy, 0)
tmp1 = tmp2 = rep(0,100)
for(i in 1:100)
{
  tmp1[i] = xx[i]-xx[i+1]
  tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)
myauc
## 0.8066951

#----------------------------------------------
#------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[data$Pam50Subtype!="NC", ]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(mycox, newdata= data, na.action=na.exclude)
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx,]
risk = risk[comxx]

##
xx = cbind(data$t.rfs, risk, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
    if(xx[i,"event"]==0 & xx[j,"event"]==0) next
    x1 = xx[i,"time"]-xx[j,"time"]
    x2 = xx[i,"score"]-xx[j,"score"]
    
    if(xx[i,"event"]==1 & xx[j,"event"]==1)
    {
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x1*x2>0)  
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x1*x2<0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
      next
    }
    if(xx[i,"event"]==1 & xx[j,"event"]==0)
    {
      if(x1>0)  next
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x2<0) 
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x2>0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
      next
    }
    if(xx[i,"event"]==0 & xx[j,"event"]==1)
    {
      if(x1<0)  next
      if(x2==0) 
      {	
        tt=tt+1
        f.mat[i,"tie"] = f.mat[i,"tie"] +1
        f.mat[j,"tie"] = f.mat[j,"tie"] +1
      }
      if(x2>0) 
      {
        cc = cc+1
        f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
        f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
      }
      if(x2<0) 
      {
        dd = dd+1
        f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
        f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
      }
    }
  }

#----------------------------------------------
cc
dd

count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
cv.res = cbind(predictability, f.mat, count, accuracy)
cv.res = cv.res[order(cv.res[,"Y"], decreasing=T), ]
plot(cv.res[,"Y"], cv.res[, "accuracy"])
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs")
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs", method="s")
## 0.3481439

#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]

w.siz = 200
step = 20
cnum = ceiling((nrow(tmp)-w.siz)/step)
yy1 = rep(0, cnum)
for(k in 1:cnum)
{
  sta = (k-1)*step+1
  end = min(sta+w.siz-1, length(yy1))
  tmpxx = tmp[sta:end,]
  yy1[k] = (sum(tmpxx[, "concordance"]) + 0.5*sum(tmpxx[, "tie"]))/sum(tmpxx[, "count"])
}
plot(yy1, pch=20, type="b")

#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
  tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
  sta = (k-1)*step+1
  end = min(sta+w.siz-1, length(cum.con))
  yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")
abline(h=cum.con[length(cum.con)], col=2)
