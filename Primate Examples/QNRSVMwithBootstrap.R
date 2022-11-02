
# Package Dependencies ----------------------------------------------------

library("apex")
library("combinat")
library("ape")
library("phangorn")
library("ips")
library("phyclust")
library("e1071")
library("class")
library("insect")
library("caret")
library("beepr")
library("randomForest")
library("apex")


# Load QNR Package and Prep Sequence Data ---------------------------------
load("SVMmodelName.R")
source("C:/path/FunctionsParallel.R")


setwd("C:/path to data")
cwd_path <-getwd()


# Seq prep if your data is in.fasta format --------------------------------



## If your data is separated be gene: Run this to collect different sequence files
## otherwise load your sequence data (.fasta format as AllSeqData)

# files<-dir(cwd_path, pattern="aligned.fasta",full=TRUE)
# cleanSeqfiles<-dir(cwd_path, pattern="cleaned.fasta",full=TRUE)
# AllSeqData<-read.multiFASTA(cleanSeqfile)
# AllSeqData<-concatenate(AllSeqData)

## use if your data in in Paup.nex format
PrimateData<-read.nexus.data("primate.paup.nex")
PrimateData<-nexus2DNAbin(PrimateData)
AllSeqData<-as.matrix(PrimateData)


## Data naming for output files
DataType<-"Primates"
ModelType<-"9_1_22"
DataDescription<-"PrimateSampling"

##Choose subset of Samples to use, by tip label index not name If you with to run all Samples
## nothing to change if you wwant to dellete all gaps, and to run all subseted of your data
Samples<-1:nrow(AllSeqData)

SeqData<-AllSeqData[Samples,]
SeqData<-del.colgapsonly(SeqData,0.01)


# QNR-SVM Should run, but not need to modify. -----------------------------


gc()


quartetgroups<-combn(Samples,4)

### Estimate Network for Each Quartet Group

Networks<- matrix(nrow=ncol(quartetgroups),ncol=Nbootstrap)

#ProbVector<-matrix(nrow=ncol(quartetgroups),ncol=24)
#NetworksProb<- rep(NA,ncol(quartetgroups))
for (i in 1:ncol(quartetgroups)){
  print(i)
  quartet<-quartetgroups[,i]
  #  qdna<-SeqData[quartet,]
  qdna<-PrimateData[quartet,]
  qdna<-del.colgapsonly(qdna,0.1)
  BioVector<-Bio_Data_to_Vector(qdna)
  total<-sum(BioVector)
  BioVector<-BioVector/total
  BioQ<-GeneTree_dist_Qmaker(BioVector)
  BioInv<-ScoreAllInv_Sym(BioQ)

  BioInv<-t(as.matrix(BioInv))
  colnames(BioInv) <- {}
  Networks[i,1]<-predict(svm.fit.All, BioInv)
  for (j in 2:100){
    trials<-rep(NA,dim(qdna)[2])
    trials<-sample(1:256,size=dim(qdna)[2],replace = TRUE, prob=BioVector)

    bsSample<-1:256
    for (k in 1:256){
      bsSample[k]<-sum(trials==k)
    }
    bsSample<-bsSample/length(trials)
    BioQ<-GeneTree_dist_Qmaker(bsSample)
    BioInv<-ScoreAllInv_Sym(BioQ)

    BioInv<-t(as.matrix(BioInv))
    colnames(BioInv) <- {}
    Networks[i,j]<-predict(svm.fit.All, BioInv)
    rm(trials)
  }

  gc()
}

## BS Analaysis
BootData<-as.data.frame(cbind(t(quartetgroups),Networks))
BootData$SelectedNetwork<-NA
BootData$BSValue<-NA



for(i in 1:dim(BootData)[1]){
  BootData[i,"SelectedNetwork"]<-BootData[i,5]
  BootData[i,"BSValue"]<-sum(BootData[i,6:105]==BootData[i,"SelectedNetwork"])
}

outputfile<-cbind(BootData[,1:4],BootData$SelectedNetwork,BootData$BSValue)

# #output to CSV ----------------------------------------------------------


#  Can comment out depending on the data you wish to use.

write.csv(outputfile,paste(DataType,ModelType,DataDescription,"BS.csv",sep=""))
write.csv(BootData,paste(DataType,ModelType,DataDescription,"All_Data.csv",sep=""))
write.csv(Samples,paste(DataType,ModelType,DataDescription,"Samples.csv",sep=""))


