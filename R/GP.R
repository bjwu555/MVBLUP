# The script was used to obtain the trait predictions of training set and test set using GBLUP. 
# Notes:
# dt: phenotype of training set
# A: kinship matrix of training set 
# A_train: kinship matrix of training set 
# A_test: kinship matrix 
# nut: the location of test set 
# num_test: the number of individuals in test set   
# trainop: phenotype of training set
# testop: phenotype of test set

PR<-function(dt,A,nut){
dt_M<-dt
dt_M[nut]<-NA
data <- data.frame(y=dt_M,gid=rownames(A))
ans <- kin.blup(data=data,geno="gid",pheno="y",K=A)
pl_te<-ans$g[nut]+mean(dt[-nut])
pcor_te<-cor(as.numeric(dt[nut]),pl_te)
Model=list(pcor_te=pcor_te,pl_te=pl_te) 
return(Model)
}

Pre_ER5<-function(dt,A_train,num_test){
A <-A_train
res1<-PR(dt,A,1:num_test)
res2<-PR(dt,A,(num_test+1):(num_test*2))
res3<-PR(dt,A,(num_test*2+1):(num_test*3))
res4<-PR(dt,A,(num_test*3+1):(num_test*4))
res5<-PR(dt,A,(num_test*4+1):nrow(A))
res_Pcor<-c(res1$pcor_te,res2$pcor_te,res3$pcor_te,res4$pcor_te,res5$pcor_te)
res_Pl<-c(res1$pl_te,res2$pl_te,res3$pl_te,res4$pl_te,res5$pl_te)
Pcor5<-sum(res_Pcor)/5
res=c(res_Pl,res_Pcor,Pcor5)
return(res)
}


T_ER<-function(trainop,testop,A_test){
A <- A_test
length_test<-nrow(A_test)-length(trainop)
sop<-c(rep(1,length_test),trainop)
sop[1:length_test]<-NA
data <- data.frame(y=sop,gid=rownames(A))
ans <- kin.blup(data=data,geno="gid",pheno="y",K=A)
pl_te<-ans$g[1:length_test]+mean(trainop)
pcor_te<-cor(as.numeric(testop),pl_te)
res=c(pl_te,pcor_te)
return(res)
}