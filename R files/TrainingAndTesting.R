# set working directory to source file location

source("SamplingFunctions.R")

# required packages
library("ape")
library("phangorn")
library("ips")
library("phyclust")
library("e1071")
library("class")
library("caret")
library("beepr")
library("randomForest")
library("ggplot2")
library("reshape2")
library("foreach")
library("doParallel")

#############Training and Testing#############################

#Parameters

numCores <- detectCores()
registerDoParallel(numCores-4)

sites  <- 1000000

# The parameter "size" is used for the training set. 
#Will get this number*8 samples from each network structure
size <- 1200 

#Range for the mixing parameter
gam    <- c(.25, .75)

#Min branch length
mn <- .05

#Max branch length
mx <- .4

# The parameter "testsize" is used for the training set. 
#Will get this number*8 samples from each network structure

testsize <- 25 


##################Create Training Set#################################

# Sample Fourier coordinates

TreeSamples      <- SampleStructure(sites = sites, size = size, mx = mx, mn=mn
                                    , gamma_range = gam, structure = "Tree")
S3CycleSamples   <- SampleStructure(sites = sites, size = 2*ceiling(size), mx = mx, mn=mn
                                    , gamma_range = gam, structure = "S3Cycle")
S4CycleSamples   <- SampleStructure(sites = sites, size = 4*ceiling(size), mx = mx, mn=mn
                                    , gamma_range = gam, structure = "S4Cycle")
Two3CycleSamples <- SampleStructure(sites = sites, size = ceiling(size), mx = mx, mn=mn
                                    ,gamma_range = gam, structure = "Two3Cycle")


# Create SVM training set by substituting each row into the invariants and
# labeling each row by the structure of the level one network (1,..,24)


TreeSamples.SVM       <- cbind(t(apply(TreeSamples     , MARGIN = 1, ScoreAllInv_Sym)),
                               rep(PermuteLabel_Tree_vec     , dim(TreeSamples)[1]/24) )
S3CycleSamples.SVM    <- cbind(t(apply(S3CycleSamples  , MARGIN = 1, ScoreAllInv_Sym)),
                               rep(PermuteLabel_3Cycle_vec   , dim(S3CycleSamples)[1]/24))
S4CycleSamples.SVM    <- cbind(t(apply(S4CycleSamples  , MARGIN = 1, ScoreAllInv_Sym)),
                               rep(PermuteLabel_4Cycle_vec   , dim(S4CycleSamples)[1]/24))
Two3CycleSamples.SVM  <- cbind(t(apply(Two3CycleSamples, MARGIN = 1, ScoreAllInv_Sym)),
                               rep(PermuteLabel_Two3Cycle_vec, dim(Two3CycleSamples)[1]/24))




##################Fit SVM#################################



AllSamples.SVM <- rbind(
  TreeSamples.SVM
  , S3CycleSamples.SVM
  , S4CycleSamples.SVM
  , Two3CycleSamples.SVM)


svm.fit.All <- svm(V1127 ~ .
                  , data = AllSamples.SVM
                  , type = "C-classification"
                  , kernel = "linear"
                  , probability = TRUE
                  , cost = 10^{-1})


##################Create Test Set#################################

Test.Set <- NULL
Test.Set <- foreach (i = 1:testsize, .combine=rbind) %dopar% {
  L <-NULL
  for (j in 1:24){
    PS <- Tree_Sampler(lengths = NULL, mn = mn, mx = mx, sites = sites)
    L<-rbind(L,(c(GeneTree_dist_Qmaker(Prob_permuter(PS, Allperms[j , ]))
                  , PermuteLabel_Tree_vec[j])))
  }
  for (j in c(1:24,1:24)){
    PS <-S3Cycle_Sampler(lengths = NULL, mn = mn, mx = mx, sites = sites, gamma1 = runif(1, gam[1], gam[2]) )
    L <- rbind(L
               , (c(GeneTree_dist_Qmaker(Prob_permuter(PS, Allperms[j , ]))
                    , PermuteLabel_3Cycle_vec[j])))
  }
  for (j in c(1:24,1:24, 1:24, 1:24)){
    PS <-S4Cycle_Sampler(lengths = NULL, mn = mn, mx = mx, sites = sites
                         , gamma1 = runif(1, gam[1], gam[2]) )
    L <- rbind(L
               , (c(GeneTree_dist_Qmaker(Prob_permuter(PS, Allperms[j , ]))
                    , PermuteLabel_4Cycle_vec[j])))
  }
  for (j in 1:24){
    PS <-Two3Cycle_Sampler(lengths = NULL, mn = mn, mx = mx, sites = sites, gamma1 = runif(1, gam[1], gam[2]), gamma2 = runif(1,gam[1], gam[2]))
    L <- rbind(L
               , (c(GeneTree_dist_Qmaker(Prob_permuter(PS, Allperms[j , ]))
                    , PermuteLabel_Two3Cycle_vec[j])))
  }
  L
}


##################Classify Points From Test Set#################################


Predictions <- apply(Test.Set[,-16], MARGIN = 1,
                     function(r){
                       data_point_Inv <- ScoreAllInv_Sym(r)
                       data_point_Inv <- t(as.matrix(data_point_Inv))
                       colnames(data_point_Inv) <- {}
                       predict(svm.fit.All, data_point_Inv, probability = TRUE)
                     }
)


# Plot confusion matrix

percent <- function(in_table)
{
  return(sum(diag(in_table))/sum(in_table))
}

CF <- confusionMatrix(as.factor(Predictions),as.factor(Test.Set[,16]))
CFtable <- CF$table


melted_cormat <- melt(CFtable)
ggplot(data = melted_cormat, aes(x=Reference, y=Prediction, fill=value)) + 
  geom_tile() + scale_fill_gradient(low = "white", high = "deepskyblue4") + 
  ggtitle(paste("Accuracy: ", round(percent(CFtable), 4))) + 
  geom_text(aes(Reference, Prediction, label = value), color = "black", size = 4) + 
  scale_x_discrete(limits = c(1:24)) + 
  scale_y_discrete(limits = c(1:24)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Truth") +
  ylab("Predictions")




