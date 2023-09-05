#This file contains functions called when training and 
# and testing the model

#In this file we use the poset in Figure 3 as the reference poset


#Data_to_Vector: Function that converts a 4-taxon seqgen sequence
# to a length 4^n probability vector or tensor

Data_to_Vector <- function(h, structure = "vector"){

  x1 <- strsplit(substring(h[2], 11), "")[[1]]
  x2 <- strsplit(substring(h[3], 11), "")[[1]]
  x3 <- strsplit(substring(h[4], 11), "")[[1]]
  x4 <- strsplit(substring(h[5], 11), "")[[1]]

  l <- length(x1)

  Branch1 <- 1:length(l)
  Branch2 <- 1:length(l)
  Branch3 <- 1:length(l)
  Branch4 <- 1:length(l)

  for(i in 1:l){
    Branch1[i] <- x1[i]
    Branch2[i] <- x2[i]
    Branch3[i] <- x3[i]
    Branch4[i] <- x4[i]
  }

  Branch <- rbind(Branch1, Branch2, Branch3, Branch4)

  V <- matrix(0, dim(Branch)[1], dim(Branch)[2])

  for (i in 1:dim(V)[1]){
    for (j in 1:dim(V)[2]){
      if (Branch[i,j] == "A") V[i,j] <- 1
      else if (Branch[i,j] == "C") V[i,j] <- 2
      else if (Branch[i,j] == "G") V[i,j] <- 3
      else if (Branch[i,j] == "T") V[i,j] <- 4
    }
  }


  Data_Tensor <- array(0,rep(4,(length(h) - 1)))

  for (k in 1:dim(V)[2]){
    Data_Tensor[t(V[,k])] <-  Data_Tensor[t(V[ , k])] + 1
  }

  Ch <- t(gtools::permutations(4,4,repeats.allowed = T))

  Data_Vector <- integer(dim(Ch)[2])

  for (k in 1:length(Data_Vector)){
    Data_Vector[k] <-  Data_Tensor[t(matrix(Ch[ , k]))]
  }

  if (structure == "vector") return(Data_Vector)
  else if (structure == "tensor") return(Data_Tensor)
  else print("structure must equal `vector` or `tensor`")

}

# Construct the 15x256 matrix ``Transform" for converting from probability
# coordinates to Fourier coordinates

qcoords <-
  rbind(c(1, 1, 1, 1), c(1, 1, 2, 2), c(1, 2, 1, 2), c(1, 2, 2, 1),
        c(1, 2, 3, 4), c(2, 1, 1, 2), c(2, 1, 2, 1), c(2, 1, 3, 4),
        c(2, 2, 1, 1), c(2, 2, 2, 2), c(2, 3, 1, 4), c(2, 3, 2, 3),
        c(2, 3, 4, 1), c(2, 2, 3, 3), c(2, 3, 3, 2))

CH <- gtools::permutations(4, 4,repeats.allowed = TRUE)
Chi <- rbind(c(1,  1,  1,  1),
             c(1, -1,  1, -1),
             c(1, 1, -1,  -1),
             c(1,  -1, -1, 1))

Transform <- NULL
for (i in 1:15){
  Transform <-  rbind(Transform,
                      apply(CH, MARGIN = 1, function(x){
                          Chi[qcoords[i,1], x[1] ]*
                          Chi[qcoords[i,2], x[2] ]*
                          Chi[qcoords[i,3], x[3] ]*
                          Chi[qcoords[i,4], x[4] ]}))
}

# Classes: Vector showing the Jukes-Cantor equivalence classes represented by
# the rows of t(gtools::permutations(4,4,repeats.allowed = T))

Classes <-
  c(1, 2, 2, 2, 3, 4, 5, 5, 3, 5, 4, 5, 3, 5, 5, 4, 6, 7, 8, 8,
  9, 10, 11, 11, 12, 13, 14, 15, 12, 13, 15, 14, 6, 8, 7, 8, 12,
  14, 13, 15, 9, 11, 10, 11, 12, 15, 13, 14, 6, 8, 8, 7, 12, 14,
  15, 13, 12, 15, 14, 13, 9, 11, 11, 10, 10, 9, 11, 11, 7, 6, 8,
  8, 13, 12, 14, 15, 13, 12, 15, 14, 4, 3, 5, 5, 2, 1, 2, 2, 5,
  3, 4, 5, 5, 3, 5, 4, 14, 12, 13, 15, 8, 6, 7, 8, 11, 9, 10, 11,
  15, 12, 13, 14, 14, 12, 15, 13, 8, 6, 8, 7, 15, 12, 14, 13, 11,
  9, 11, 10, 10, 11, 9, 11, 13, 14, 12, 15, 7, 8, 6, 8, 13, 15,
  12, 14, 14, 13, 12, 15, 11, 10, 9, 11, 8, 7, 6, 8, 15, 13, 12,
  14, 4, 5, 3, 5, 5, 4, 3, 5, 2, 2, 1, 2, 5, 5, 3, 4, 14, 15, 12,
  13, 15, 14, 12, 13, 8, 8, 6, 7, 11, 11, 9, 10, 10, 11, 11, 9,
  13, 14, 15, 12, 13, 15, 14, 12, 7, 8, 8, 6, 14, 13, 15, 12, 11,
  10, 11, 9, 15, 13, 14, 12, 8, 7, 8, 6, 14, 15, 13, 12, 15, 14,
  13, 12, 11, 11, 10, 9, 8, 8, 7, 6, 4, 5, 5, 3, 5, 4, 5, 3, 5,
  5, 4, 3, 2, 2, 2, 1)

# GeneTree_dist_Qmaker: Function that converts a length 256 probability vector
# to Fourier coordinates (the entries of the probability vector are first
# averaged over Jukes-Cantor equivalence classes)

GeneTree_dist_Qmaker <- function(Prob_Vec){

  dim(Prob_Vec) <- c(1,256)

  for (i in 1:15){
    Prob_Vec[which(Classes == i)] <- mean(Prob_Vec[Classes == i])
  }

  Q <- Prob_Vec%*%t(Transform)

}


# Tree_Sampler: Function that returns the length 256 probability vector generated
# by a tree 12|34 with random branch lengths 

Tree_Sampler <- function(lengths = NULL, mn = 0.1, mx =0.2
                         , sites=100, gamma1 = "dummy",gamma2 = "dummy"){
  if (is.null(lengths) == T)
  {
    lengths <- c(runif(5, mn, mx))
  }
  s <- paste("((1:", lengths[5], ", 2:", lengths[4], "):", lengths[3]/2
             ,", (3:", lengths[1], ", 4:", lengths[2], "):", lengths[3]/2
             , ");", sep = "")

  anc <- sample(c(0, 1, 2, 3), sites,c(0.25, 0.25, 0.25, 0.25), replace=TRUE)
  Seq <- gen.seq.HKY(read.tree(text = s), pi = c(0.25, 0.25, 0.25, 0.25)
                     , kappa = 1, sites, anc.seq = anc)
  Seq <- Data_to_Vector(Seq)

  return(Seq/sum(Seq))

}


# S3Cycle_Sampler: Function that returns a length 256 probability vector generated
# by one of the three networks in the equivalence class of the 3-cycle 
# at far left in the reference poset (structure 4)


S3Cycle_Sampler <- function(lengths = NULL, mn = 0.1, mx =0.2
                            , sites = 100, gamma1=0.5,gamma2 = "dummy")
{
  if (is.null(lengths) == T)
  {
    lengths <- c(runif(6, mn, mx/2), runif(2, mn, mx))
  } 
  
  
  config = runif(1)
  
  if (config < .33) {
    sites1 <- floor(gamma1*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6])/2,", (2:", lengths[7], ", 1:"
                 , lengths[8], "):", (lengths[5]+lengths[6])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace=TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)
    
    sites2 <- sites - sites1
    Tr2 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2] + lengths[5]
                 , "):", (lengths[6])/2,", (2:", lengths[7], ", 1:", lengths[8]
                 , "):", (lengths[6])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3),sites2,c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) + Data_to_Vector(sort(Seq2))
  }
  else if (config < .66) {
    sites1 <- floor(gamma1*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6])/2,", (2:", lengths[7], ", 1:"
                 , lengths[8], "):", (lengths[5]+lengths[6])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace=TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)
    
    sites2 <- sites - sites1
    Tr2 <- paste("((3:", lengths[1], ", 4:", lengths[2] +lengths[3], "):", (lengths[4]+lengths[6])/2,", (2:", lengths[7], ", 1:", lengths[8], "):", (lengths[4]+lengths[6])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3),sites2,c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) + Data_to_Vector(sort(Seq2))
  }
  else {sites1 <- floor(gamma1*sites)
  Tr1 <- paste("((3:", lengths[1], ", 4:", lengths[2] +lengths[3], "):", (lengths[4]+lengths[6])/2,", (2:", lengths[7], ", 1:", lengths[8], "):", (lengths[4]+lengths[6])/2, ");", sep = "")
  anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace=TRUE)
  Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                      , kappa = 1, sites1, anc.seq = anc)

  sites2 <- sites - sites1
  Tr2 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):", (lengths[5]+lengths[6])/2,", (2:", lengths[7], ", 1:", lengths[8], "):", (lengths[5]+lengths[6])/2, ");", sep = "")
  anc <- sample(c(0,1,2,3),sites2,c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
  Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                      , kappa = 1, sites2, anc.seq = anc)

  Seq <- Data_to_Vector(sort(Seq1)) + Data_to_Vector(sort(Seq2))
  }
  return(Seq/sum(Seq))

}

# S4Cycle_Sampler: Function that returns a length 256 probability vector 
# generated by 4-cycle at far left in poset (structure 10)

S4Cycle_Sampler <- function(lengths = NULL, mn = 0.1, mx =0.2
                            , sites=100, gamma1=0.5, gamma2 = "dummy")
{

  if (is.null(lengths) == TRUE)
  {
    lengths <- c(runif(2, mn, mx/2)
                 , runif(1, mn, mx)
                 , runif(5, mn, mx/2))
  }
  sites1 <- floor(gamma1*sites)
  Tr1 <- paste("(((1:", lengths[8], ", 2:", lengths[6] + lengths[7], "):"
               , (lengths[5]),", 4:", lengths[3], "):"
               , lengths[2]/2, ",3:", lengths[1] + lengths[2]/2, ");", sep = "")
  
  
  anc <- sample(c(0,1,2,3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
  Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                      , kappa = 1, sites1, anc.seq = anc)

  sites2 <- sites - sites1
  Tr2 <- paste("((2:", lengths[4]+lengths[6], ", 3:", lengths[1], "):"
               , (lengths[2])/2,", (4:", lengths[3], ", 1:"
               , lengths[5]+lengths[8], "):", (lengths[2])/2, ");", sep = "")
  anc <- sample(c(0, 1, 2, 3), sites2, c(0.25, 0.25, 0.25, 0.25), replace=TRUE)
  Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                      , kappa = 1, sites2, anc.seq = anc)

  Seq <- Data_to_Vector(sort(Seq1)) + Data_to_Vector(sort(Seq2))
  return(Seq/sum(Seq))

}

# Two3Cycle_Sampler: Function that returns a length 256 probability vector
# generated by one of the eight networks in the equivalence class of 
# the double 3-cycle 12|34 (structure 22)

Two3Cycle_Sampler <- function(lengths = NULL, mn = 0.1, mx =0.2
                              , sites=100, gamma1=.5,gamma2=.5)
{
  if (is.null(lengths) == T)
  {
    lengths <- c(runif(3, mn, (1/2)*mx), runif(5, mn, (1/3)*mx), runif(3, mn, (1/2)*mx))
  }

 alpha <- c(gamma1*gamma2,gamma1*(1-gamma2),(1-gamma1)*gamma2,(1-gamma1)*(1-gamma2))
  
  config = runif(1)
  
  if (config < .125) {

  sites1 <- floor(alpha[1]*sites)
  Tr1 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
               , (lengths[5]+lengths[6])/2,", (2:", lengths[7]+lengths[10]
               , ", 1:", lengths[8]+lengths[11], "):"
               , (lengths[5]+lengths[6])/2, ");", sep = "")
  anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
  Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                      , kappa = 1, sites1, anc.seq = anc)

  sites2 <- floor(alpha[2]*sites)
  Tr2 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
               , (lengths[5]+lengths[6]+lengths[8])/2,", (2:"
               , lengths[9]+lengths[10], ", 1:", lengths[11], "):"
               , (lengths[5]+lengths[6]+lengths[8])/2, ");", sep = "")
  anc <- sample(c(0,1,2,3), sites2, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
  Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                      , kappa = 1, sites2, anc.seq = anc)

  sites3 <- floor(alpha[3]*sites)
  Tr3 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
               , "):", (lengths[6])/2,", (2:", lengths[7]+lengths[10], ", 1:"
               , lengths[8]+lengths[11], "):", (lengths[6])/2, ");", sep = "")
  anc <- sample(c(0,1,2,3), sites3, c(0.25, 0.25, 0.25, 0.25),replace = TRUE)
  Seq3 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                      , kappa = 1, sites3, anc.seq = anc)

  sites4 <- sites - sites1 - sites2 - sites3
  Tr4 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
               , "):", (lengths[6]+lengths[8])/2,", (2:", lengths[9]+lengths[10]
               , ", 1:", lengths[11], "):", (lengths[6]+lengths[8])/2, ");"
               , sep = "")
  anc <- sample(c(0, 1, 2, 3), sites4, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
  Seq4 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                      , kappa = 1, sites4, anc.seq = anc)

  Seq <- Data_to_Vector(sort(Seq1)) +
    Data_to_Vector(sort(Seq2)) +
    Data_to_Vector(sort(Seq3)) +
    Data_to_Vector(sort(Seq4)) }
  
  else if (config < .25) {
    
    
    sites1 <- floor(alpha[1]*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6] + lengths[7])/2,", (2:", lengths[10], ", 1:"
                 , lengths[9]+lengths[11], "):", (lengths[5] + lengths[6]+lengths[7])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)
    
    sites2 <- floor(alpha[2]*sites)
    Tr2 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6]+lengths[8])/2,", (2:"
                 , lengths[9]+lengths[10], ", 1:", lengths[11], "):"
                 , (lengths[5]+lengths[6]+lengths[8])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites2, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)
    
    sites3 <- floor(alpha[3]*sites)
    Tr3 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6]+lengths[7])/2,", (2:", lengths[10], ", 1:"
                 , lengths[9]+lengths[11], "):", (lengths[6]+lengths[7])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites3, c(0.25, 0.25, 0.25, 0.25),replace = TRUE)
    Seq3 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites3, anc.seq = anc)
    
    sites4 <- sites - sites1 - sites2 - sites3
    Tr4 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6]+lengths[8])/2,", (2:", lengths[9]+lengths[10]
                 , ", 1:", lengths[11], "):", (lengths[6]+lengths[8])/2, ");"
                 , sep = "")
    anc <- sample(c(0, 1, 2, 3), sites4, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq4 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites4, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) +
      Data_to_Vector(sort(Seq2)) +
      Data_to_Vector(sort(Seq3)) +
      Data_to_Vector(sort(Seq4)) }
  
  else if (config < .375) {
  
    sites1 <- floor(alpha[1]*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6]+lengths[7])/2,", (2:", lengths[10], ", 1:"
                 , lengths[9]+lengths[11], "):", (lengths[5] + lengths[6]+lengths[7])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)
 
    sites2 <- floor(alpha[2]*sites)
    Tr2 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6])/2,", (2:", lengths[7]+lengths[10]
                 , ", 1:", lengths[8]+lengths[11], "):"
                 , (lengths[5]+lengths[6])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites2, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)
 
    sites3 <- floor(alpha[3]*sites)
    Tr3 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6]+lengths[7])/2,", (2:", lengths[10], ", 1:"
                 , lengths[9]+lengths[11], "):", (lengths[6]+lengths[7])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites3, c(0.25, 0.25, 0.25, 0.25),replace = TRUE)
    Seq3 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites3, anc.seq = anc)
    
    sites4 <- sites - sites1 - sites2 - sites3

    Tr4 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6])/2,", (2:", lengths[7]+lengths[10], ", 1:"
                 , lengths[8]+lengths[11], "):", (lengths[6])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites4, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq4 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites4, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) +
      Data_to_Vector(sort(Seq2)) +
      Data_to_Vector(sort(Seq3)) +
      Data_to_Vector(sort(Seq4)) }
  
  else if (config < .5) {

    sites1 <- floor(alpha[1]*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6]+lengths[8])/2,", (2:", lengths[9]+lengths[10]
                 , ", 1:", lengths[11], "):", (lengths[6]+lengths[8])/2, ");"
                 , sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)

    sites2 <- floor(alpha[2]*sites)
    Tr2 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6])/2,", (2:", lengths[7]+lengths[10], ", 1:"
                 , lengths[8]+lengths[11], "):", (lengths[6])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites2, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)

    sites3 <- floor(alpha[3]*sites)
    Tr3 <- paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                 , "):", (lengths[4]+lengths[6]+lengths[8])/2,", (2:", lengths[9]+lengths[10]
                 , ", 1:", lengths[11], "):", (lengths[4]+lengths[6]+lengths[8])/2, ");"
                 , sep = "")
    anc <- sample(c(0,1,2,3), sites3, c(0.25, 0.25, 0.25, 0.25),replace = TRUE)
    Seq3 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites3, anc.seq = anc)
  
    sites4 <- sites - sites1 - sites2 - sites3
    Tr4 <- paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                 , "):", (lengths[4]+lengths[6])/2,", (2:", lengths[7]+lengths[10], ", 1:"
                 , lengths[8]+lengths[11], "):", (lengths[4]+lengths[6])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites4, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq4 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites4, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) +
      Data_to_Vector(sort(Seq2)) +
      Data_to_Vector(sort(Seq3)) +
      Data_to_Vector(sort(Seq4)) }
  
  else if (config < .625) {
    
    sites1 <- floor(alpha[1]*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6]+lengths[7])/2,", (2:", lengths[10], ", 1:"
                 , lengths[9]+lengths[11], "):", (lengths[6]+lengths[7])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)
    
    sites2 <- floor(alpha[2]*sites)
    Tr2 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6])/2,", (2:", lengths[7]+lengths[10], ", 1:"
                 , lengths[8]+lengths[11], "):", (lengths[6])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites2, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)
    
    sites3 <- floor(alpha[3]*sites)
    Tr3 <- paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                 , "):", (lengths[4]+lengths[6]+lengths[7])/2,", (2:", lengths[10]
                 , ", 1:", lengths[9] + lengths[11], "):", (lengths[4]+lengths[6]+lengths[7])/2, ");"
                 , sep = "")
    anc <- sample(c(0,1,2,3), sites3, c(0.25, 0.25, 0.25, 0.25),replace = TRUE)
    Seq3 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites3, anc.seq = anc)
    
    sites4 <- sites - sites1 - sites2 - sites3
    Tr4 <- paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                 , "):", (lengths[4]+lengths[6])/2,", (2:", lengths[7]+lengths[10], ", 1:"
                 , lengths[8]+lengths[11], "):", (lengths[4] +lengths[6])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites4, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq4 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites4, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) +
      Data_to_Vector(sort(Seq2)) +
      Data_to_Vector(sort(Seq3)) +
      Data_to_Vector(sort(Seq4)) }
  
  else if (config < .75) {
    
    sites1 <- floor(alpha[1]*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6]+lengths[7])/2,", (2:", lengths[10], ", 1:"
                 , lengths[9]+lengths[11], "):", (lengths[6]+lengths[7])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)
    
    sites2 <- floor(alpha[2]*sites)
    Tr2 <- paste("((3:", lengths[1]+lengths[4], ", 4:", lengths[2]+lengths[5]
                 , "):", (lengths[6]+lengths[8])/2,", (2:", lengths[9]+lengths[10]
                 , ", 1:", lengths[11], "):", (lengths[6]+lengths[8])/2, ");"
                 , sep = "")
    anc <- sample(c(0,1,2,3), sites2, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)
    
    sites3 <- floor(alpha[3]*sites)
    Tr3 <- paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                 , "):", (lengths[4]+lengths[6]+lengths[7])/2,", (2:", lengths[10]
                 , ", 1:", lengths[9] + lengths[11], "):", (lengths[4]+lengths[6]+lengths[7])/2, ");"
                 , sep = "")
    anc <- sample(c(0,1,2,3), sites3, c(0.25, 0.25, 0.25, 0.25),replace = TRUE)
    Seq3 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites3, anc.seq = anc)
    
    sites4 <- sites - sites1 - sites2 - sites3
    Tr4 <-  paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                  , "):", (lengths[4]+lengths[6]+lengths[8])/2,", (2:", lengths[9]+lengths[10]
                  , ", 1:", lengths[11], "):", (lengths[4]+lengths[6]+lengths[8])/2, ");"
                  , sep = "")
    anc <- sample(c(0, 1, 2, 3), sites4, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq4 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites4, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) +
      Data_to_Vector(sort(Seq2)) +
      Data_to_Vector(sort(Seq3)) +
      Data_to_Vector(sort(Seq4)) }
  
  else if (config < .875) {
    
    sites1 <- floor(alpha[1]*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6] + lengths[7])/2,", (2:", lengths[10], ", 1:"
                 , lengths[9]+lengths[11], "):", (lengths[5] + lengths[6]+lengths[7])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)
    
    sites2 <- floor(alpha[2]*sites)
    Tr2 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6])/2,", (2:", lengths[7]+lengths[10]
                 , ", 1:", lengths[8]+lengths[11], "):"
                 , (lengths[5]+lengths[6])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites2, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)
    
    sites3 <- floor(alpha[3]*sites)
    Tr3 <- paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                 , "):", (lengths[4]+lengths[6]+lengths[7])/2,", (2:", lengths[10]
                 , ", 1:", lengths[9] + lengths[11], "):", (lengths[4]+lengths[6]+lengths[7])/2, ");"
                 , sep = "")
    anc <- sample(c(0,1,2,3), sites3, c(0.25, 0.25, 0.25, 0.25),replace = TRUE)
    Seq3 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites3, anc.seq = anc)
    
    sites4 <- sites - sites1 - sites2 - sites3
    Tr4 <-  paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                  , "):", (lengths[4]+lengths[6])/2,", (2:", lengths[7]+lengths[10], ", 1:"
                  , lengths[8]+lengths[11], "):", (lengths[4] + lengths[6])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites4, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq4 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites4, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) +
      Data_to_Vector(sort(Seq2)) +
      Data_to_Vector(sort(Seq3)) +
      Data_to_Vector(sort(Seq4)) }
  
   else{
    
    sites1 <- floor(alpha[1]*sites)
    Tr1 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6])/2,", (2:", lengths[7]+lengths[10]
                 , ", 1:", lengths[8]+lengths[11], "):"
                 , (lengths[5]+lengths[6])/2, ");", sep = "")
    anc <- sample(c(0, 1, 2, 3), sites1, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq1 <- gen.seq.HKY(read.tree(text = Tr1), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites1, anc.seq = anc)
    
    sites2 <- floor(alpha[2]*sites)
    Tr2 <- paste("((3:", lengths[1]+lengths[3], ", 4:", lengths[2], "):"
                 , (lengths[5]+lengths[6]+lengths[8])/2,", (2:"
                 , lengths[9]+lengths[10], ", 1:", lengths[11], "):"
                 , (lengths[5]+lengths[6]+lengths[8])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites2, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq2 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites2, anc.seq = anc)
   
    sites3 <- floor(alpha[3]*sites)
    Tr3 <- paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                 , "):", (lengths[4]+lengths[6])/2,", (2:", lengths[7]+lengths[10], ", 1:"
                 , lengths[8]+lengths[11], "):", (lengths[4] + lengths[6])/2, ");", sep = "")
    anc <- sample(c(0,1,2,3), sites3, c(0.25, 0.25, 0.25, 0.25),replace = TRUE)
    Seq3 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites3, anc.seq = anc)
    
    sites4 <- sites - sites1 - sites2 - sites3
    Tr4 <- paste("((3:", lengths[1], ", 4:", lengths[2]+lengths[3]
                 , "):", (lengths[4]+lengths[6]+lengths[8])/2,", (2:", lengths[9]+lengths[10]
                 , ", 1:", lengths[11], "):", (lengths[4] + lengths[6]+lengths[8])/2, ");"
                 , sep = "")
    anc <- sample(c(0, 1, 2, 3), sites4, c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
    Seq4 <- gen.seq.HKY(read.tree(text = Tr2), pi = c(0.25, 0.25, 0.25, 0.25)
                        , kappa = 1, sites4, anc.seq = anc)
    
    Seq <- Data_to_Vector(sort(Seq1)) +
      Data_to_Vector(sort(Seq2)) +
      Data_to_Vector(sort(Seq3)) +
      Data_to_Vector(sort(Seq4)) }

  return(Seq/sum(Seq))

}


# Prob_permuter: Function that applies the permutation "permute" to the
# coordinates of a length 256 probabilty vector.

## Example: Applying the permutation (1,4,2,3) to a probability vector
## from a tree with leaf labels 1,2,3,4 returns a probability vector
## from a tree with leaf labels 1,3,4,2

Ch <- t(gtools::permutations(4,4,repeats.allowed = T))
Ch_permute <- t(gtools::permutations(4,4,repeats.allowed = T))

Prob_permuter <- function(P, permute){

  P_permute <- numeric(256)

  Ch_permute[ 1, ] <- Ch[ permute[1], ]
  Ch_permute[ 2, ] <- Ch[ permute[2], ]
  Ch_permute[ 3, ] <- Ch[ permute[3], ]
  Ch_permute[ 4, ] <- Ch[ permute[4], ]

  for (i in 1:256){
    P_permute[i] <- P[
                      which(apply(Ch, MARGIN = 2
                      , function(x){sum(Ch_permute[,i] == x)}) == 4)
                     ]
  }

  return(P_permute)

}


# PermuteLeaves: Function that takes a length 256 probability vector from a
# 4-leaf phylogenetic network and returns the 24x15 matrix of the
# Fourier coordinates for each permutation of the original structure.

Allperms   <- gtools::permutations(4,4)

PermuteLeaves <- function(P){

  TT <- NULL
  for (i in 1:24){
    TT <- rbind(TT, GeneTree_dist_Qmaker(Prob_permuter(P, Allperms[ i, ])))
  }

  return(TT)

}

# Vectors of structure labels after "PermuteLeaves" is applied to
# Tree_Sampler, S3Cycle_Sampler, S4Cycle_Sampler, Two3Cycle_Sampler

PermuteLabel_Tree_vec <- c(1, 1, 2, 2, 3, 3, 1, 1
                         , 3, 3, 2, 2, 2, 2, 3, 3
                         , 1, 1, 3, 3, 2, 2, 1, 1)

PermuteLabel_3Cycle_vec <- c(4, 4, 5, 5, 6, 6, 4, 4
                           , 7, 7, 8, 8, 5, 5, 7, 7
                           , 9, 9, 6, 6, 8, 8, 9, 9)

PermuteLabel_4Cycle_vec <- c(10, 12, 20, 13, 21, 11, 18, 16
                           , 20, 17, 21, 19, 18, 15, 10, 14
                           , 11, 19, 16, 15, 12, 14, 13, 17)

PermuteLabel_Two3Cycle_vec <- 21 + PermuteLabel_Tree_vec

# SampleStructure: Randomly samples from a tree, 3-cycle, 4-cycle, or double-triangle
# depending on if structure = "Tree", "S3Cycle", "S4Cycle", or "Two3Cycle."
# Output: A 24x15 matrix of the Fourier coordinates for each
# permutation of the original structure.

SampleStructure <- function(sites, size = 100, mx = .3, mn = .1
                            , gamma_range = c(.5, .5), structure = "Tree"){
  Sampler <- eval(parse(text = paste(structure, "_Sampler", sep = "")))
  gamma   <- runif(2, gamma_range[1], gamma_range[2])
  TT<-foreach (i = 1:size, .combine=rbind) %dopar% {
    PermuteLeaves(Sampler(mn = mn, mx = mx, sites = sites,
                                          gamma1 = gamma[1],gamma2=gamma[2]))
  }
  return(TT)
}


# ScoreAllInv_Sym: Takes a length 15 vector of Fourier coordinates and
# returns a length 1126 vector where each entry is obtained by substituting
# the Fourier coordinates into a polynomial invariant for a 4-leaf level-one
# network.

ScoreAllInv_Sym  <- function(vec)
{
  Q1 <- vec[1]
  Q2 <- vec[2]
  Q3 <- vec[3]
  Q4 <- vec[4]
  Q5 <- vec[5]
  Q6 <- vec[6]
  Q7 <- vec[7]
  Q8 <- vec[8]
  Q9 <- vec[9]
  Q10 <- vec[10]
  Q11 <- vec[11]
  Q12 <- vec[12]
  Q13 <- vec[13]
  Q14 <- vec[14]
  Q15 <- vec[15]
  Inv <- c(Q5*Q7*Q9*Q11^2+Q4*Q7*Q11^3-2*Q3*Q7*Q11^2*Q13+Q2*Q9*Q11^2*Q13+Q3*Q6*Q11*Q13^2+2*Q3*Q7*Q9*Q11*Q14-2*Q2*Q9^2*Q11*Q14-Q1*Q9*Q11*Q12*Q14-Q3*Q6*Q9*Q13*Q14-2*Q1*Q11^2*Q13*Q14+2*Q1*Q9*Q11*Q14^2-Q3*Q8*Q9^2*Q15-Q3*Q7*Q9*Q11*Q15+2*Q2*Q9^2*Q11*Q15+Q1*Q9*Q11*Q12* Q15-Q1*Q9*Q11*Q14*Q15,
           -2*Q5^2*Q8^2*Q11^2+2*Q5^2*Q6*Q8*Q11*Q12+2*Q2*Q5*Q8*Q11^2*Q12+Q3*Q5*Q8^2*Q11*Q14+Q2*Q5* Q8*Q11^2*Q14-Q3*Q5*Q6*Q8*Q12*Q14-3*Q2*Q5*Q6*Q11*Q12*Q14-Q2*Q3*Q8*Q11*Q12*Q14-Q2^2* Q11^2*Q12*Q14+Q5^2*Q6^2*Q14^2+Q3*Q5*Q6*Q8*Q14^2+Q2*Q3*Q6*Q12*Q14^2-Q5^2*Q6*Q8*Q11*Q15+ Q3*Q5*Q8^2*Q11*Q15-2*Q2*Q5*Q8*Q11^2*Q15+Q2^2*Q11^2*Q12*Q15-Q3^2*Q8^2*Q14*Q15+Q2^2*Q11^2*Q14*Q15-Q2*Q3*Q6*Q14^2*Q15+Q2*Q3*Q8*Q11* Q15^2,
           -Q5*Q9*Q11*Q12-2*Q5*Q11^2*Q13-Q3*Q11*Q13* Q14+Q5*Q9*Q11*Q15+Q4*Q11^2*Q15+Q3*Q11*Q13* Q15+Q3*Q9*Q15^2,
           Q5^2*Q7*Q9-Q4^2*Q8*Q11-Q2*Q4*Q9*Q12-2* Q3*Q5*Q7*Q13+Q2*Q4*Q11*Q13+2*Q1*Q5*Q12*Q13+Q2*Q3*Q13^2-2*Q1*Q5*Q13*Q14+Q3*Q4*Q7* Q15-Q1*Q4*Q12*Q15+Q1*Q4*Q14*Q15,
           Q3*Q4*Q6*Q9*Q12-Q1*Q4*Q11^2*Q12-Q3*Q5*Q6*Q9*Q13+Q3*Q4*Q6*Q11*Q13-2*Q1*Q5* Q11^2*Q13-Q3^2*Q6*Q13^2+Q1*Q5*Q9*Q11*Q14+Q1*Q4*Q11^2*Q14+3*Q1*Q3*Q11*Q13*Q14-Q1*Q3*Q11* Q13*Q15-Q1*Q3*Q9*Q14*Q15,
           -Q5^2*Q7*Q8*Q9-Q4^2*Q8^2*Q11+Q3*Q5*Q7*Q8*Q13-Q3*Q5*Q7^2*Q14+2*Q3*Q4*Q7*Q8*Q14+ Q2*Q5*Q7*Q9*Q14+Q1*Q4*Q8*Q12*Q14-Q1*Q5*Q8*Q13*Q14-Q1*Q5*Q7*Q14^2+Q2*Q5*Q7* Q9*Q15-Q2*Q3*Q7*Q13*Q15+Q1*Q5*Q8*Q13*Q15-Q1*Q5*Q7*Q14*Q15,
           -Q5^2*Q6*Q8*Q9+Q4*Q5*Q6*Q8*Q11+Q2*Q5*Q6*Q9*Q12-Q2*Q4*Q6*Q11*Q12+Q1*Q5*Q8*Q11*Q12-Q3^2*Q8^2*Q13-Q4*Q5*Q6^2*Q14+2*Q3* Q4*Q6*Q8*Q14+Q2*Q5*Q6*Q9*Q14-Q1*Q5*Q8*Q11*Q14-Q1*Q5*Q6*Q12*Q14-Q1*Q5*Q6* Q14^2+Q1*Q3*Q8*Q14*Q15,
           Q4*Q8*Q11-Q2*Q9*Q15,
           Q3*Q5*Q7+Q3*Q4*Q8-Q1*Q5*Q12-Q2*Q3* Q13,
           Q4*Q6*Q12-Q3*Q8*Q13+Q2*Q11*Q13-Q1* Q12*Q15,
           Q3*Q5*Q7-Q3*Q4*Q8-Q1*Q5*Q12+ Q2*Q3*Q13,
           Q4*Q7*Q11-Q1*Q13*Q15,
           -Q5*Q11*Q12*Q13+Q5*Q9*Q12*Q14+Q4*Q11*Q12* Q14+Q5*Q11*Q13*Q14-Q4*Q11*Q14^2-Q3*Q13*Q14^2-Q5*Q9*Q12*Q15-Q5*Q11*Q13*Q15+Q5*Q9* Q14*Q15+Q3*Q13*Q14*Q15,
           -Q5*Q7+Q4*Q8,
           Q5*Q7-Q4*Q8,
           Q2*Q4*Q5*Q6*Q11-Q2*Q3*Q4*Q8*Q11-2*Q1*Q5^2*Q8*Q11-Q2^2*Q4*Q11^2+Q1*Q5^2*Q6*Q12+ Q1*Q3*Q5*Q8*Q12+3*Q1*Q2*Q5*Q11*Q12+Q2* Q3*Q4*Q6*Q14-Q1*Q5^2*Q6*Q14-Q1*Q2*Q5* Q11*Q15-Q1*Q2*Q3*Q12*Q15,
           -Q5*Q8*Q9^2+Q4*Q7*Q11^2+Q3*Q8*Q9*Q13-2*Q4*Q6*Q11*Q13-2*Q1*Q11*Q12*Q13+Q3*Q6*Q13^2+ Q4*Q6*Q9*Q14+Q1*Q9*Q12*Q14-Q3*Q7*Q9*Q15+2*Q1*Q11*Q13*Q15-Q1*Q9*Q14*Q15,
           -Q3*Q5*Q7*Q9+Q3*Q4*Q8*Q9+Q2*Q5*Q9^ 2-Q2*Q4*Q9*Q11+Q1*Q5*Q9*Q12-2*Q1*Q5* Q11*Q13+Q1*Q3*Q13*Q14+Q1*Q5*Q9*Q15-Q1*Q3*Q13*Q15,
           Q5*Q8*Q9*Q11-Q5*Q6*Q9*Q12+2*Q4*Q6*Q11*Q12-2*Q3*Q7*Q11*Q12+2*Q1*Q11*Q12^2- Q2*Q11^2*Q13+Q3*Q7*Q11*Q14-Q3*Q6*Q13*Q14- Q1*Q11*Q12*Q15,
           -Q5^2*Q6*Q9^2+Q4^2*Q6*Q11^2+2*Q1*Q5*Q9*Q11* Q12+Q1*Q4*Q11^2*Q12+Q3^2*Q8*Q9*Q13-Q3*Q4*Q6*Q11*Q13-Q3^2*Q7*Q11*Q13-Q1*Q3*Q11* Q12*Q13+Q3^2*Q6*Q13^2+Q3*Q4*Q6*Q9*Q14-Q1*Q4*Q11^2*Q14-Q1*Q3*Q11*Q13*Q14-Q3*Q4*Q6* Q9*Q15-Q1*Q4*Q11^2*Q15-Q1*Q3*Q9*Q12*Q15+Q1*Q3*Q11*Q13*Q15+Q1*Q3*Q9*Q15^2,
           -Q6*Q12*Q13+Q7*Q11*Q15,
           Q6*Q12*Q13-Q7*Q11*Q15,
           -Q3^2*Q7+Q2*Q3*Q9-Q1*Q5*Q11+Q1*Q3*Q12,
           Q3*Q8^2*Q9-2*Q3*Q7*Q8*Q11+Q2*Q7*Q11^2-Q2*Q6*Q9*Q12+2*Q1*Q8*Q11*Q12-Q5*Q6^2*Q13+ Q2*Q6*Q11*Q13-2*Q1*Q8*Q11*Q14+Q3*Q6*Q7* Q15-Q1*Q6*Q12*Q15+Q1*Q6*Q14*Q15,
           -Q5^2*Q6*Q8*Q13+Q3*Q5*Q8^2*Q13-Q2^2*Q11*Q12*Q13+Q5^2*Q6*Q7*Q14-Q3*Q5*Q7*Q8*Q14+ Q2*Q3*Q8*Q13*Q14+Q5^2*Q6*Q7*Q15-Q3*Q5*Q7*Q8*Q15-Q2*Q5*Q6*Q13*Q15+Q2*Q3*Q8* Q13*Q15+Q2^2*Q11*Q13*Q15-Q2*Q3*Q7*Q14*Q15,
           Q3^2*Q5*Q7*Q8-Q3^2*Q4*Q8^2-Q1*Q5^2*Q6* Q12-Q1*Q3*Q5*Q8*Q12+Q1*Q2*Q5*Q11*Q12-Q2*Q3*Q5*Q6*Q13-Q2*Q3^2*Q8*Q13+Q1*Q5^2* Q6*Q14+2*Q1*Q3*Q5*Q8*Q14-Q1*Q2*Q3*Q12* Q14+Q2*Q3^2*Q7*Q15,
           -Q3*Q5^2*Q7^2+2*Q3*Q4*Q5*Q7*Q8-Q3*Q4^2*Q8^2-Q2*Q4^2*Q8*Q11+Q1*Q4*Q5*Q8*Q12+Q2* Q3*Q5*Q7*Q13-Q2*Q3*Q4*Q8*Q13-Q1*Q5^2* Q7*Q14+Q1*Q4*Q5*Q8*Q14+Q2*Q3*Q4*Q7* Q15-Q1*Q5^2*Q7*Q15,
           -Q8*Q9^2*Q12^2+Q7*Q11^2*Q12*Q13-Q6*Q11*Q12*Q13^2+Q7*Q9*Q11*Q12*Q14-Q7*Q11^2*Q13*Q14+Q6* Q11*Q13^2*Q14+Q8*Q9^2*Q12*Q15+Q7*Q9*Q11*Q12* Q15-Q7*Q11^2*Q13*Q15-Q6*Q9*Q12*Q13*Q15+Q6*Q11*Q13^2*Q15-Q7*Q9*Q11*Q14*Q15,
           Q2*Q3*Q7-Q1*Q5*Q8-Q2^2*Q9+Q1*Q2*Q14,
           -Q2*Q11*Q12+Q5*Q6*Q14,
           Q5*Q8*Q9*Q12+Q4*Q8*Q11*Q12+Q5*Q8*Q11* Q13-Q3*Q8*Q12*Q13-Q2*Q11*Q12*Q13+Q4*Q8* Q11*Q14-Q4*Q6*Q12*Q14+Q3*Q7*Q12*Q14-Q1* Q12^2*Q14-Q2*Q11*Q13*Q14+Q3*Q7*Q14^2-Q1*Q12*Q14^2-Q5*Q8*Q9*Q15-Q3*Q8*Q13*Q15+Q1* Q12*Q14*Q15+Q1*Q14^2*Q15,
           Q8*Q11-Q6*Q12,
           -Q5^2*Q7*Q8*Q11+Q4*Q5*Q8^2*Q11+Q5^2*Q6*Q7*Q12-Q4*Q5*Q6*Q8*Q12-Q2*Q5*Q7*Q11* Q12+Q2*Q4*Q8*Q11*Q12+Q2^2*Q11*Q12*Q13+Q5^2*Q6*Q7*Q14-Q4*Q5*Q6*Q8*Q14+Q2*Q4*Q8* Q11*Q14-Q2*Q4*Q6*Q12*Q14-Q2^2*Q11*Q13*Q15,
           -Q3*Q5*Q7*Q12+Q2*Q5*Q9*Q12-Q2*Q4*Q11* Q12+Q1*Q5*Q12^2-Q5^2*Q6*Q13-Q3*Q5*Q8* Q13+Q2*Q5*Q11*Q13+Q3*Q4*Q8*Q14-Q3*Q5*Q7*Q15+Q1*Q5*Q12*Q15+Q2*Q3*Q13*Q15-Q1*Q5*Q14*Q15,
           Q2^2*Q5*Q9*Q11-Q2^2*Q4*Q11^2+Q1*Q5^2*Q6*Q12+2*Q1*Q2*Q5*Q11*Q12-Q2*Q3*Q5* Q6*Q13-Q2^2*Q3*Q11*Q13-Q1*Q5^2*Q6*Q14+Q1*Q3*Q5*Q8*Q14-Q1*Q2*Q5*Q11*Q14-Q1*Q2*Q3*Q12*Q14+Q2^2*Q3*Q9*Q15,
           Q1*Q12^2-Q5*Q6*Q13+Q3*Q8*Q13-Q2*Q11*Q13,
           -Q1*Q12^2+Q4*Q6*Q14,
           Q3*Q7*Q14-Q1*Q15^2,
           -Q5*Q8*Q11*Q12+Q5*Q8*Q11*Q14+Q5*Q6*Q12* Q14+Q2*Q11*Q12*Q14-Q5*Q6*Q14^2-Q3*Q8*Q14^2-Q5*Q8*Q11*Q15-Q2*Q11*Q12*Q15+Q3*Q8* Q14*Q15+Q2*Q11*Q14*Q15,
           -Q5*Q8*Q12*Q13+Q4*Q8*Q12*Q14+Q5*Q8*Q13* Q14+Q2*Q12*Q13*Q14-Q5*Q7*Q14^2-Q4*Q8*Q14^2-Q5*Q8*Q13*Q15-Q2*Q12*Q13*Q15+Q5*Q7* Q14*Q15+Q2*Q13*Q14*Q15,
           Q2*Q4*Q11-Q1*Q5*Q14,
           Q5*Q6*Q7*Q9-Q3*Q7*Q8*Q9-Q4*Q6*Q7* Q11+Q3*Q7^2*Q11-2*Q1*Q8*Q11*Q13+Q1*Q6*Q12*Q13+Q1*Q7*Q11*Q14-Q1*Q6*Q13*Q14+Q1*Q7*Q11*Q15,
           Q5*Q7*Q9-Q1*Q12*Q13,
           -Q3*Q5*Q6*Q7+Q3^2*Q7*Q8-Q2*Q3*Q8*Q9-2*Q1*Q5*Q8*Q11+Q1*Q2*Q11*Q12+Q2*Q3* Q6*Q13+Q1*Q3*Q8*Q14+Q1*Q3*Q8*Q15-Q1*Q2*Q11*Q15,
           -Q5*Q7*Q8*Q12*Q13+Q4*Q8^2*Q12* Q13+2*Q5*Q8^2*Q13^2-2*Q5*Q7*Q8*Q13*Q14-Q4*Q8^ 2*Q13*Q14+Q2*Q7*Q12*Q13*Q14-Q2*Q8*Q13^2* Q14+Q5*Q7^2*Q14^2-Q2*Q7*Q12*Q13*Q15+Q2*Q7*Q13*Q14*Q15,
           -Q5*Q8*Q9-Q5*Q7*Q11+Q2*Q9* Q12+Q5*Q6*Q13-Q1*Q12*Q14+Q1*Q12*Q15,
           -Q8*Q11*Q13-Q7*Q11*Q14+Q8*Q9*Q15+Q7*Q11* Q15,
           Q5^2*Q7*Q8*Q9+Q2*Q4*Q8*Q9*Q12-Q1*Q5*Q8*Q12*Q13-Q2^2*Q9*Q12*Q13-Q2*Q5* Q7*Q9*Q14-Q1*Q4*Q8*Q12*Q14+Q1*Q5*Q8*Q13*Q14+Q1*Q2*Q12*Q13*Q14+Q1*Q5*Q7*Q14^2+ Q1*Q4*Q8*Q12*Q15-Q1*Q5*Q8*Q13*Q15-Q1*Q2*Q12*Q13*Q15,
           Q5*Q7*Q8^2*Q9+Q4*Q8^3*Q9-2*Q3*Q7^2*Q8*Q12+2*Q2*Q7*Q8*Q9*Q12+2*Q1*Q7*Q8*Q12^2+Q3*Q7*Q8^2*Q13-2*Q2*Q8^2*Q9*Q13- Q2*Q6*Q7*Q12*Q13-2*Q1*Q8^2*Q12*Q13+Q2*Q6* Q8*Q13^2-Q1*Q7*Q8*Q12*Q14+2*Q3*Q7^2*Q8* Q15-Q2*Q7*Q8*Q9*Q15-Q2*Q7^2*Q11*Q15-Q1*Q7*Q8*Q12*Q15+Q1*Q7*Q8*Q14*Q15,
           Q5*Q7*Q9-Q4*Q8*Q9+Q2*Q9*Q13-Q1*Q13* Q14,
           -Q5*Q7*Q9+Q4*Q8*Q9+Q2*Q9*Q13-Q1* Q13*Q14,
           Q8*Q13-Q7*Q14,
           -Q4*Q8*Q14+Q2*Q13*Q15,
           Q4*Q8*Q14-Q2*Q13*Q15,
           Q3*Q5*Q7*Q9-Q3*Q4*Q8*Q9-Q1*Q5*Q9* Q12+Q3*Q4*Q6*Q13-Q3^2*Q7*Q13-2*Q1*Q5*Q11*Q13+Q1*Q3*Q12*Q13+Q1*Q5*Q9*Q15+Q1*Q4*Q11*Q15,
           Q8*Q9*Q13-Q7*Q11*Q13-Q6*Q13^2+Q7*Q9*Q15,
           -Q8*Q9*Q13+Q7*Q11*Q13+Q6*Q13^2- Q7*Q9*Q15,
           Q8*Q9*Q13-Q7*Q11*Q13+Q6*Q13^2- Q7*Q9*Q15,
           -Q8*Q9*Q13+Q7*Q11*Q13-Q6*Q13^2+ Q7*Q9*Q15,
           Q8*Q9*Q13+Q7*Q11*Q13-Q6*Q13^ 2-Q7*Q9*Q15,
           Q4^2*Q5*Q7*Q11-Q4^3*Q8*Q11-Q3*Q4^2*Q8*Q13+Q2*Q4^2*Q11*Q13-Q1*Q4*Q5* Q12*Q13+Q2*Q3*Q4*Q13^2-Q1*Q5^2*Q13^2+Q1*Q4^2*Q12*Q14,
           -Q8^2*Q9*Q12*Q13+Q6*Q8*Q12*Q13^2+ Q7*Q8*Q9*Q12*Q14+Q8^2*Q9*Q13*Q14-Q6*Q8* Q13^2*Q14-Q7^2*Q11*Q14^2-Q7*Q8*Q9*Q12*Q15-Q8^2*Q9*Q13*Q15+Q6*Q8*Q13^2*Q15+Q7*Q8*Q9* Q14*Q15+Q7^2*Q11*Q14*Q15-Q6*Q7*Q13*Q14*Q15,
           Q5*Q8*Q9-Q5*Q7*Q11-Q5*Q6*Q13+Q1*Q14^2,
           -Q5*Q6*Q7*Q9+Q3*Q7*Q8*Q9+Q4*Q6*Q7*Q11-Q3*Q7^2*Q11-Q1*Q8*Q9*Q12+Q1*Q7* Q11*Q12-2*Q1*Q8*Q11*Q13+Q1*Q8*Q9*Q15+Q1*Q6*Q13*Q15,
           -Q5^2*Q6*Q7^2+Q3*Q5*Q7^2*Q8-Q2* Q5*Q7^2*Q11-Q2*Q4*Q7*Q8*Q11-Q1*Q5*Q7* Q8*Q12-Q1*Q4*Q8^2*Q12+Q1*Q2*Q8*Q12*Q13+2*Q1*Q5*Q7*Q8*Q14+Q1*Q4*Q8^2*Q14-Q1* Q2*Q7*Q12*Q14+Q2*Q3*Q7^2*Q15,
           Q4*Q5^2*Q6^2-Q3^2*Q4*Q8^2-Q2*Q4*Q5*Q6* Q11-Q2^2*Q5*Q9*Q11+Q2^2*Q4*Q11^2+Q2*Q3*Q4* Q6*Q12-Q1*Q5^2*Q6*Q12-Q1*Q2*Q5*Q11*Q12+Q2^2*Q3*Q11*Q13+Q1*Q5^2*Q6*Q14+2*Q1*Q3* Q5*Q8*Q14-Q1*Q2*Q5*Q11*Q14-Q2*Q3*Q4*Q6*Q15-Q1*Q5^2*Q6*Q15+Q1*Q2*Q5*Q11*Q15- Q1*Q2*Q3*Q14*Q15+Q1*Q2*Q3*Q15^2,
           Q5*Q7*Q8*Q11-Q4*Q8^2*Q11-Q2*Q8*Q11*Q13+Q2*Q6*Q12*Q13-Q5*Q6*Q7*Q14+Q4*Q6*Q8*Q14+Q1*Q8*Q12*Q14+Q4*Q6*Q8*Q15-Q3*Q7*Q8*Q15+Q2*Q7*Q11*Q15-Q1*Q8*Q14*Q15-Q1*Q8*Q15^2,
           -Q8*Q9+Q7*Q11,
           Q8*Q9-Q7*Q11,
           -Q4*Q11*Q12-Q5*Q11*Q13+Q4*Q11*Q14+Q3*Q13* Q14,
           Q2*Q5^2*Q6*Q7-Q1*Q5^2*Q8^2-Q2^2*Q5*Q7*Q11-Q1*Q2*Q5*Q8*Q12+Q2^2*Q5*Q6*Q13+ Q2^2*Q3*Q8*Q13-Q2^3*Q11*Q13+Q1*Q2^2*Q12*Q15,
           Q5*Q9*Q12+Q4*Q11*Q12-Q5*Q11*Q13-Q4*Q11*Q14,
           -Q4*Q8*Q12*Q13-2*Q5*Q8*Q13^2+Q5*Q7* Q13*Q14+Q4*Q8*Q13*Q14+Q2*Q13^2*Q14+Q4*Q7*Q14^2-Q5*Q7*Q13*Q15,
           Q4*Q5*Q7*Q9*Q11-Q4^2*Q7*Q11^2-Q3*Q4* Q8*Q9*Q13-Q1*Q5*Q11*Q13^2-Q1*Q4*Q9*Q12*Q14+2*Q1*Q4*Q11*Q13*Q14+Q1*Q3*Q13^2*Q14+Q1* Q4*Q9*Q12*Q15+Q1*Q5*Q9*Q13*Q15-Q1*Q3*Q13^2*Q15-Q1*Q4*Q9*Q14*Q15,
           -Q5^2*Q6*Q7^2+Q4^2*Q6*Q8^2+Q2*Q4*Q6*Q7*Q12-Q1*Q4*Q8^2*Q12-Q2*Q4*Q6*Q8*Q13-Q2^2* Q8*Q9*Q13+Q2^2*Q7*Q11*Q13-Q1*Q2*Q8*Q12* Q13+Q2^2*Q6*Q13^2+2*Q1*Q5*Q7*Q8*Q14+Q1* Q4*Q8^2*Q14-Q1*Q2*Q8*Q13*Q14-Q2*Q4*Q6*Q7*Q15-Q1*Q4*Q8^2*Q15+Q1*Q2*Q8*Q13*Q15- Q1*Q2*Q7*Q14*Q15+Q1*Q2*Q7*Q15^2,
           Q5^2*Q7*Q9-Q4*Q5*Q7*Q11+Q3*Q4*Q8*Q13+Q2*Q3*Q13^2-Q4^2*Q6*Q14+Q2*Q4*Q9*Q14- Q1*Q4*Q12*Q14-Q1*Q4*Q14^2+Q1*Q4*Q12*Q15-2*Q1*Q5*Q13*Q15+Q1*Q4*Q14*Q15,
           -Q5*Q8^2*Q9^2+2*Q4*Q6*Q8*Q9*Q12+Q4*Q6*Q8*Q11*Q13-Q4*Q6^2*Q12*Q13+Q3*Q6*Q7*Q12* Q13-Q1*Q8*Q11*Q12*Q13-Q1*Q6*Q12^2*Q13-Q3*Q6*Q8*Q13^2-Q4*Q6*Q7*Q11*Q14+Q3*Q6*Q7* Q13*Q14+Q1*Q8*Q11*Q13*Q14-Q1*Q6*Q12*Q13*Q14+Q1*Q8*Q9*Q12*Q15,
           Q5^2*Q7*Q9+Q4*Q5*Q7*Q11-Q3*Q4*Q8*Q13+Q2*Q3*Q13^2-Q4^2*Q6*Q14+Q2*Q4*Q9*Q14- Q1*Q4*Q12*Q14-Q1*Q4*Q14^2+Q1*Q4*Q12*Q15-2*Q1*Q5*Q13*Q15+Q1*Q4*Q14*Q15,
           -Q4*Q8*Q9*Q12+Q5*Q7*Q11*Q13+Q4*Q6*Q12* Q13-Q3*Q8*Q13^2-Q5*Q7*Q9*Q15-2*Q4*Q6* Q13*Q15+2*Q2*Q9*Q13*Q15-Q1*Q13*Q14*Q15+2*Q1*Q13*Q15^2,
           -Q4^2*Q8*Q9*Q11-Q4^2*Q7*Q11^2+Q4^2* Q6*Q9*Q12-Q3*Q4*Q8*Q9*Q13+Q4^2*Q6*Q11*Q13+2*Q1*Q4*Q11*Q13*Q14+Q1*Q3*Q13^2*Q14+Q1* Q5*Q9*Q13*Q15-Q1*Q4*Q11*Q13*Q15-Q1*Q3*Q13^2*Q15-Q1*Q4*Q9*Q14*Q15,
           Q4*Q8*Q12-Q5*Q8*Q13+Q2*Q12*Q13-Q2*Q13* Q15,
           Q3*Q5*Q7^2*Q13-Q3*Q4*Q7*Q8*Q13-Q1*Q5*Q7*Q12*Q13-Q2*Q3*Q7*Q13^2+2*Q1*Q5* Q8*Q13^2+Q1*Q2*Q12*Q13^2+Q3*Q4*Q7^2*Q14-Q1* Q4*Q7*Q12*Q14-Q1*Q4*Q8*Q13*Q15-Q1*Q2* Q13^2*Q15+Q1*Q4*Q7*Q14*Q15,
           -Q5^2*Q6*Q9-Q3*Q5*Q8*Q9+2*Q4*Q5*Q6*Q11-Q3*Q5*Q7*Q11-Q2*Q4*Q11^2+2*Q1*Q5*Q11* Q12-Q3*Q5*Q6*Q13+Q2*Q3*Q11*Q13+Q3*Q4*Q6*Q14-Q3^2*Q7*Q14+Q1*Q3*Q12*Q14-Q1*Q3* Q14*Q15,
           -Q3^2*Q4*Q8^2+Q2*Q3*Q4*Q8*Q11-Q1* Q5^2*Q8*Q11-Q1*Q5^2*Q6*Q12+Q1*Q2*Q5* Q11*Q12-Q2*Q3*Q5*Q6*Q13+Q1*Q5^2*Q6*Q14+2*Q1*Q3*Q5*Q8*Q14-Q1*Q2*Q3*Q12*Q14+Q1* Q2*Q3*Q12*Q15-Q1*Q2*Q3*Q14*Q15,
           -Q2*Q9*Q12+Q3*Q8*Q13,
           Q5*Q6*Q8*Q9+Q3*Q8^2*Q9+Q2*Q7*Q11^2- Q4*Q6^2*Q12+Q3*Q6*Q7*Q12-Q1*Q6*Q12^2-Q2*Q6*Q11*Q13-Q1*Q6*Q12*Q14-2*Q1*Q8*Q11*Q15+ Q1*Q6*Q12*Q15+Q1*Q6*Q14*Q15,
           -Q5*Q6*Q8*Q9+Q3*Q8^2*Q9+Q2*Q7*Q11^2-Q4*Q6^2*Q12+Q3*Q6*Q7*Q12-Q1*Q6*Q12^2+Q2* Q6*Q11*Q13-Q1*Q6*Q12*Q14-2*Q1*Q8*Q11*Q15+Q1*Q6*Q12*Q15+Q1*Q6*Q14*Q15,
           -Q4^2*Q8^2*Q9-Q4^2*Q7*Q8*Q11+Q4^2*Q6*Q8*Q13-Q2*Q4*Q7*Q11*Q13+2*Q1*Q4*Q8*Q12*Q13+ Q1*Q2*Q12*Q13^2+Q4^2*Q6*Q7*Q14-Q1*Q4*Q7* Q12*Q15+Q1*Q5*Q7*Q13*Q15-Q1*Q4*Q8*Q13* Q15-Q1*Q2*Q13^2*Q15,
           Q3*Q8*Q13-Q1*Q15^2,
           Q5^2*Q6*Q7-2*Q4*Q5*Q6*Q8+Q3*Q4*Q8^2+ Q2*Q5*Q7*Q11-2*Q1*Q5*Q8*Q12-Q2^2*Q11*Q13+Q2*Q4*Q6*Q14+Q1*Q2*Q12*Q14-Q2*Q3*Q7* Q15+2*Q1*Q5*Q8*Q15-Q1*Q2*Q14*Q15,
           -Q5*Q7*Q8*Q11-Q4*Q8^2*Q11+Q4*Q6*Q8*Q12+Q2*Q8*Q11*Q13-Q2*Q6*Q12*Q13+Q5*Q6*Q7*Q14+Q1*Q8*Q12*Q14+Q4*Q6*Q8*Q15-Q2*Q8*Q9*Q15+Q2*Q7*Q11*Q15-Q1*Q8*Q12*Q15-Q1*Q8*Q15^2,
           Q5*Q8*Q9-Q1*Q15^2,
           Q8*Q11-Q6*Q15,
           Q2*Q7*Q9-Q1*Q8*Q13-Q1*Q7*Q14+Q1*Q7* Q15,
           -Q8*Q11+Q6*Q15,
           Q5*Q6*Q8*Q9-Q3*Q8^2*Q9-Q5*Q6*Q7*Q11-Q4*Q6*Q8*Q11+2*Q3*Q7*Q8*Q11-Q2*Q7* Q11^2-Q2*Q6*Q11*Q13-Q4*Q6^2*Q14+Q3*Q6*Q7*Q14-Q1*Q6*Q12*Q14+2*Q1*Q8*Q11*Q15+Q1*Q6* Q14*Q15,
           -Q1*Q14^2+Q3*Q7*Q15,
           Q3*Q5^2*Q7^2-Q3*Q4^2*Q8^2-Q2*Q3*Q4*Q7* Q12-Q1*Q5^2*Q7*Q12+Q1*Q2*Q4*Q12^2-Q2*Q3* Q5*Q7*Q13-Q2^2*Q5*Q9*Q13+Q2^2*Q4*Q11*Q13+Q1*Q2*Q5*Q12*Q13+Q2^2*Q3*Q13^2+Q1*Q5^2* Q7*Q14+2*Q1*Q4*Q5*Q8*Q14-Q1*Q2*Q4*Q12*Q14-Q1*Q2*Q5*Q13*Q14+Q2*Q3*Q4*Q7*Q15- Q1*Q5^2*Q7*Q15-Q1*Q2*Q5*Q13*Q15,
           Q5*Q6*Q9*Q12-Q5*Q6*Q11*Q13+Q3*Q8*Q11* Q13-Q2*Q11^2*Q13-Q3*Q7*Q11*Q14+Q2*Q9*Q11* Q14+Q3*Q6*Q13*Q14-Q1*Q11*Q14^2-Q3*Q8*Q9*Q15+Q2*Q9*Q11*Q15+Q1*Q11*Q12*Q15-Q1*Q11*Q14*Q15,
           Q3*Q4*Q8*Q9-Q3*Q4*Q7*Q11+Q1*Q5*Q9*Q12+Q3^2*Q7*Q13-Q2*Q3*Q9*Q13-2* Q1*Q5*Q11*Q13+Q1*Q3*Q13*Q14-Q1*Q5*Q9*Q15+Q1*Q3*Q13*Q15,
           -Q5*Q7*Q9+Q1*Q13*Q15,
           Q7*Q11-Q6*Q13,
           -Q7*Q11+Q6*Q13,
           Q5*Q7*Q9-Q1*Q13*Q15,
           Q5*Q8*Q9*Q12+Q4*Q8*Q11*Q12+Q5*Q8*Q11* Q13-Q5*Q6*Q12*Q13-Q3*Q8*Q12*Q13-Q4*Q8* Q11*Q14-Q3*Q8*Q13*Q14+Q5*Q8*Q9*Q15+Q3* Q7*Q12*Q15-Q2*Q9*Q12*Q15-Q1*Q12^2*Q15-Q5* Q6*Q13*Q15+Q1*Q12*Q14*Q15+Q3*Q7*Q15^2-Q1*Q12*Q15^2+Q1*Q14*Q15^2,
           -Q5*Q7^2*Q9*Q11-Q4*Q7*Q8*Q9*Q11-Q4*Q7^2*Q11^2+Q4*Q6*Q7*Q9*Q12+Q4*Q6*Q8*Q9* Q13+2*Q4*Q6*Q7*Q11*Q13-Q4*Q6^2*Q13^2-Q1*Q6* Q12*Q13^2+Q1*Q7*Q11*Q13*Q14-Q1*Q6*Q13^2*Q14+Q1*Q7*Q11*Q13*Q15,
           Q4*Q7*Q8*Q11*Q12+Q4*Q8^2*Q11*Q13-Q4*Q6*Q8*Q12*Q13-Q2*Q8*Q11*Q13^2+Q2*Q6*Q12* Q13^2+Q5*Q7^2*Q11*Q14+Q4*Q7*Q8*Q11*Q14-Q4*Q6*Q7*Q12*Q14-Q4*Q6*Q8*Q13*Q14-Q2*Q7* Q11*Q13*Q14+Q2*Q6*Q13^2*Q14-Q5*Q7^2*Q11*Q15,
           Q4*Q8^2-Q2*Q7*Q15,
           Q3*Q4*Q6-Q1*Q5*Q11+Q1*Q3*Q14-Q1*Q3* Q15,
           -Q4*Q5*Q6*Q8*Q12+Q3*Q4*Q8^2*Q12+Q5^2*Q6*Q8*Q13-Q3*Q5*Q8^2*Q13+Q2*Q5*Q6* Q12*Q13-Q2*Q3*Q8*Q12*Q13+Q2^2*Q11*Q12*Q13-Q4*Q5*Q6*Q8*Q14+Q3*Q4*Q8^2*Q14-Q2*Q4* Q6*Q12*Q14+Q2*Q5*Q6*Q13*Q14-Q2^2*Q11*Q13*Q15,
           Q3*Q8^2*Q9+Q5*Q6*Q7*Q11+Q2*Q7*Q11^2- Q3*Q6*Q8*Q13-Q4*Q6^2*Q14+Q2*Q6*Q9*Q14-Q1*Q6*Q12*Q14-Q1*Q6*Q14^2-2*Q1*Q8*Q11* Q15+Q1*Q6*Q12*Q15+Q1*Q6*Q14*Q15,
           Q3*Q8^2*Q9-Q5*Q6*Q7*Q11+Q2*Q7*Q11^2+ Q3*Q6*Q8*Q13-Q4*Q6^2*Q14+Q2*Q6*Q9*Q14-Q1*Q6*Q12*Q14-Q1*Q6*Q14^2-2*Q1*Q8*Q11*Q15+ Q1*Q6*Q12*Q15+Q1*Q6*Q14*Q15,
           Q5^2*Q7*Q9+Q4*Q5*Q7*Q11-Q4^2*Q8*Q11-2* Q2*Q5*Q9*Q13-2*Q1*Q5*Q12*Q13+Q2*Q3*Q13^2-Q3*Q4*Q7*Q14+2*Q1*Q5*Q13*Q14+Q2*Q4*Q9* Q15+Q1*Q4*Q12*Q15-Q1*Q4*Q14*Q15,
           Q4^2*Q6*Q9*Q12-Q4*Q5*Q6*Q9*Q13+Q4^2* Q6*Q11*Q13-Q3*Q4*Q6*Q13^2+2*Q1*Q5*Q11*Q13^2+Q1*Q4*Q9*Q12*Q14-Q1*Q5*Q9*Q13*Q14-Q1* Q3*Q13^2*Q14-Q1*Q4*Q9*Q12*Q15-Q1*Q4*Q11*Q13*Q15+Q1*Q3*Q13^2*Q15,
           -Q5^2*Q6*Q7^2+2*Q4*Q5*Q6*Q7*Q8-Q4^2*Q6*Q8^2-Q2*Q5*Q7^2*Q11+Q2*Q4*Q6*Q7*Q12-Q1* Q4*Q8^2*Q12-Q2*Q5*Q6*Q7*Q13+Q2*Q4*Q6* Q8*Q13+Q1*Q5*Q7*Q8*Q14-Q1*Q4*Q8^2*Q14+Q1*Q5*Q7*Q8*Q15,
           Q2*Q3*Q13-Q1*Q5*Q14,
           -2*Q5^2*Q8^2*Q11^2+Q5^2*Q6*Q8*Q11*Q12+Q3*Q5* Q8^2*Q11*Q12+Q2*Q5*Q6*Q11*Q12^2+Q2^2*Q11^2*Q12^2+Q5^2*Q6*Q8*Q11*Q14-2*Q3*Q5*Q8^2*Q11*Q14-Q2*Q5*Q8*Q11^2*Q14-Q5^2*Q6^2*Q12*Q14+Q3^2* Q8^2*Q12*Q14-Q2*Q3*Q6*Q12^2*Q14+Q3*Q5*Q6* Q8*Q14^2+2*Q3*Q5*Q8^2*Q11*Q15+2*Q2*Q5*Q8* Q11^2*Q15-Q3*Q5*Q6*Q8*Q12*Q15-Q3^2*Q8^2*Q12* Q15-Q2*Q5*Q6*Q11*Q12*Q15-3*Q2*Q3*Q8*Q11*Q12*Q15+Q2*Q3*Q6*Q12^2*Q15+Q3^2*Q8^2*Q14* Q15,
           Q3*Q7*Q8*Q9*Q11-Q3*Q7^2*Q11^2-Q5*Q6^2 *Q9*Q13-Q3*Q6*Q8*Q9*Q13+2*Q3*Q6*Q7*Q11*Q13+Q1*Q6*Q11*Q12*Q13-Q3*Q6^2*Q13^2-Q1* Q7*Q11^2*Q14+Q1*Q6*Q11*Q13*Q14+Q3*Q6*Q7*Q9*Q15-Q1*Q7*Q11^2*Q15,
           -Q5*Q8*Q9^2+Q5*Q7*Q9*Q11+Q4*Q7*Q11^2-2*Q4*Q6*Q11*Q13-2*Q1*Q11*Q12*Q13+Q3*Q6*Q13^2+ Q4*Q6*Q9*Q14+Q1*Q9*Q12*Q14-Q3*Q7*Q9*Q15+2*Q1*Q11*Q13*Q15-Q1*Q9*Q14*Q15,
           Q2*Q5*Q8*Q9*Q11-Q2*Q5*Q7*Q11^2+2*Q2* Q5*Q6*Q9*Q12+Q2*Q3*Q7*Q11*Q12-Q1*Q5*Q8*Q11*Q12-Q2^2*Q9*Q11*Q12-Q1*Q2*Q11*Q12^2- Q5^2*Q6^2*Q13+Q1*Q5*Q6*Q12*Q14-Q2*Q3*Q8*Q9*Q15+Q2*Q3*Q7*Q11*Q15+Q1*Q5*Q8*Q11* Q15-Q1*Q2*Q11*Q12*Q15,
           Q5*Q6*Q9-Q3*Q8*Q9-Q4*Q6*Q11+Q3*Q6* Q13,
           2*Q8^2*Q9*Q11*Q12*Q13-2*Q8^2*Q11^2*Q13^2+2*Q6*Q8*Q11*Q12*Q13^2-Q8^2*Q9^2*Q12*Q14-Q7*Q8* Q9*Q11*Q12*Q14+Q8^2*Q9*Q11*Q13*Q14+Q7*Q8*Q11^2*Q13*Q14-3*Q6*Q8*Q9*Q12*Q13*Q14-Q6*Q7* Q11*Q12*Q13*Q14+Q6*Q7*Q9*Q12*Q14^2+Q6*Q7* Q11*Q13*Q14^2+Q6^2*Q13^2*Q14^2+Q8^2*Q9^2*Q12*Q15-2 *Q8^2*Q9*Q11*Q13*Q15+Q7*Q8*Q11^2*Q13*Q15-Q6*Q8*Q11*Q13^2*Q15+Q8^2*Q9^2*Q14*Q15-Q7^2*Q11^2 *Q14*Q15-Q6*Q7*Q9*Q14^2*Q15+Q7*Q8*Q9*Q11*Q15^2,
           Q5*Q6*Q13-Q3*Q7*Q15,
           -Q5*Q8^2*Q9+Q2*Q8*Q11*Q13-Q5*Q6*Q7*Q14+Q4*Q6*Q8*Q14-2*Q4*Q6*Q8*Q15+2*Q3*Q7* Q8*Q15-Q2*Q7*Q11*Q15-Q1*Q8*Q12*Q15+2*Q1*Q8*Q15^2,
           -Q5*Q7*Q11+Q4*Q8*Q11+Q3*Q7*Q12- Q2*Q9*Q12-Q1*Q12^2-Q5*Q6*Q13+Q3*Q8*Q13+ Q3*Q7*Q15-Q1*Q12*Q15+Q1*Q14*Q15,
           Q5*Q7*Q11+Q4*Q8*Q11+Q3*Q7*Q12-Q2*Q9* Q12-Q1*Q12^2-Q5*Q6*Q13-Q3*Q8*Q13+Q3*Q7* Q15-Q1*Q12*Q15+Q1*Q14*Q15,
           -2*Q8^2*Q9*Q11*Q12*Q13-Q7*Q8*Q11^2*Q12*Q13+Q6*Q8*Q9*Q12^2*Q13-2*Q8^2*Q11^2*Q13^2+Q6*Q8* Q11*Q12*Q13^2+Q8^2*Q9^2*Q12*Q14+Q8^2*Q9*Q11*Q13* Q14+Q6*Q8*Q11*Q13^2*Q14-Q6^2*Q12*Q13^2*Q14+Q7^2*Q11^2*Q14^2-Q6*Q7*Q9*Q12*Q14^2+Q6*Q7*Q11* Q13*Q14^2+Q8^2*Q9^2*Q12*Q15+2*Q8^2*Q9*Q11*Q13*Q15+2*Q7*Q8*Q11^2*Q13*Q15-Q8^2*Q9^2*Q14*Q15-3* Q7*Q8*Q9*Q11*Q14*Q15-Q6*Q8*Q9*Q13*Q14*Q15-Q6*Q7*Q11*Q13*Q14*Q15+Q6*Q7*Q9*Q14^2* Q15,
           Q5*Q7*Q11-Q4*Q8*Q11+Q3*Q7*Q12-Q2* Q9*Q12-Q1*Q12^2+Q5*Q6*Q13-Q3*Q8*Q13+Q3*Q7*Q15-Q1*Q12*Q15+Q1*Q14*Q15,
           -Q5*Q7*Q11-Q4*Q8*Q11+Q3*Q7*Q12-Q2*Q9* Q12-Q1*Q12^2+Q5*Q6*Q13+Q3*Q8*Q13+Q3*Q7*Q15-Q1*Q12*Q15+Q1*Q14*Q15,
           -Q5*Q8*Q11*Q12-Q5*Q8*Q11*Q14-Q5*Q6*Q12* Q14+Q5*Q8*Q11*Q15+Q5*Q6*Q12*Q15+Q2*Q11* Q12*Q15+Q5*Q6*Q14*Q15+Q3*Q8*Q14*Q15-Q3* Q8*Q15^2-Q2*Q11*Q15^2,
           Q2*Q3*Q9-Q1*Q5*Q11,
           Q4*Q6*Q7-Q1*Q8*Q13,
           -Q5*Q8*Q9-Q5*Q7*Q11+Q5*Q6*Q13+Q3*Q7* Q14-Q1*Q12*Q14+Q1*Q14*Q15,
           -Q5^2*Q6^2*Q7+Q4*Q5*Q6^2*Q8+Q2*Q4*Q6^2*Q12-Q2*Q5*Q6^2*Q13-Q2*Q3*Q6*Q8*Q13+2*Q1* Q5*Q6*Q8*Q14+Q1*Q3*Q8^2*Q14-Q1*Q5*Q6* Q8*Q15-Q1*Q3*Q8^2*Q15+Q1*Q2*Q8*Q11*Q15-Q1*Q2*Q6*Q14*Q15,
           -Q5*Q8*Q9^2+Q4*Q8*Q9*Q11+Q4*Q7*Q11^2-Q4*Q6*Q9*Q12-2*Q3*Q7*Q11*Q13+2*Q1*Q11*Q12* Q13+Q3*Q6*Q13^2+Q3*Q7*Q9*Q14-Q1*Q9*Q12*Q14-2*Q1*Q11*Q13*Q15+Q1*Q9*Q14*Q15,
           -Q8*Q11*Q13-Q6*Q12*Q13+Q7*Q11*Q14+Q6*Q13* Q14,
           -Q5^2*Q7^2*Q9-Q4*Q5*Q7*Q8*Q9-Q4* Q5*Q7^2*Q11+Q2*Q4*Q7*Q9*Q12+2*Q2*Q5*Q7*Q9*Q13+Q2*Q4*Q8*Q9*Q13-Q2^2*Q9*Q13^2-Q1* Q2*Q12*Q13^2+Q1*Q5*Q7*Q13*Q14+Q1*Q5*Q7* Q13*Q15-Q1*Q2*Q13^2*Q15,
           Q3*Q8^2*Q9*Q11+Q5*Q6*Q7*Q11^2+Q4*Q6* Q8*Q11^2-2*Q3*Q7*Q8*Q11^2+Q2*Q7*Q11^3+2*Q4*Q6^ 2*Q11*Q14-Q3*Q6*Q7*Q11*Q14+Q1*Q6*Q11*Q12*Q14-Q3*Q6^2*Q13*Q14-Q3*Q6*Q8*Q9*Q15-2* Q4*Q6^2*Q11*Q15+2*Q3*Q6*Q7*Q11*Q15-2*Q1*Q8* Q11^2*Q15-Q1*Q6*Q11*Q12*Q15-Q1*Q6*Q11*Q14* Q15+2*Q1*Q6*Q11*Q15^2,
           -Q4^2*Q8^2*Q11+2*Q2*Q4*Q8*Q9*Q12+Q2*Q5*Q8*Q9*Q13+Q2*Q3*Q7*Q12*Q13-Q1*Q5*Q8* Q12*Q13-Q2^2*Q9*Q12*Q13-Q1*Q2*Q12^2*Q13-Q2*Q3*Q8*Q13^2+Q1*Q4*Q8*Q12*Q14-Q2*Q5*Q7* Q9*Q15+Q2*Q3*Q7*Q13*Q15+Q1*Q5*Q8*Q13*Q15-Q1*Q2*Q12*Q13*Q15,
           -Q5^2*Q7^2*Q9+Q4^2*Q8^2*Q9+Q4^2*Q7*Q8*Q11+ Q2*Q4*Q7*Q9*Q12-Q4^2*Q6*Q8*Q13-Q2*Q4* Q8*Q9*Q13-Q1*Q4*Q8*Q12*Q13+Q2^2*Q9*Q13^2-Q1*Q2*Q12*Q13^2-Q2*Q4*Q7*Q9*Q14+Q1* Q4*Q8*Q13*Q14-Q1*Q2*Q13^2*Q14+Q1*Q4*Q7*Q14^2+2*Q1*Q5*Q7*Q13*Q15-Q1*Q4*Q8*Q13*Q15+ Q1*Q2*Q13^2*Q15-Q1*Q4*Q7*Q14*Q15,
           -Q4*Q7*Q8*Q9*Q11-Q4*Q7^2*Q11^2+Q4*Q6*Q7*Q9*Q12+Q4*Q6*Q7*Q11*Q13-2*Q1*Q8*Q11* Q13^2-Q1*Q6*Q12*Q13^2+Q1*Q8*Q9*Q13*Q14+3*Q1* Q7*Q11*Q13*Q14+Q1*Q6*Q13^2*Q14-Q1*Q7*Q11* Q13*Q15-Q1*Q7*Q9*Q14*Q15,
           Q10-Q14,
           -Q3*Q4*Q5*Q8*Q9-Q1*Q5^2*Q9*Q12+Q1*Q4*Q5*Q11*Q12+Q3^2*Q5*Q7*Q13-Q3^2*Q4*Q8* Q13-Q1*Q3*Q5*Q12*Q13-Q2*Q3^2*Q13^2+Q3^2*Q4* Q7*Q14+Q1*Q5^2*Q9*Q15-Q1*Q3*Q4*Q12*Q15+ 2*Q1*Q3*Q5*Q13*Q15,
           Q5*Q6*Q13-Q2*Q9*Q15,
           -Q4*Q8^2*Q11-2*Q3*Q7*Q8*Q12+2*Q2*Q8*Q9*Q12+2*Q1*Q8*Q12^2+Q5*Q6*Q8*Q13-Q2*Q6*Q12* Q13-Q1*Q8*Q12*Q14+Q3*Q7*Q8*Q15-Q2*Q7*Q11*Q15,
           -Q5*Q8*Q13-Q5*Q7*Q14+Q5*Q7*Q15+ Q2*Q13*Q15,
           Q8^2*Q9^2*Q12+Q7*Q8*Q9*Q11*Q12-2*Q8^2*Q9*Q11*Q13+3*Q6*Q8*Q9*Q12*Q13+Q6* Q7*Q11*Q12*Q13-2*Q6*Q8*Q11*Q13^2+Q6^2*Q12*Q13^ 2-Q6*Q7*Q9*Q12*Q14-Q6^2*Q13^2*Q14-Q8^2*Q9^2* Q15-Q6*Q7*Q9*Q12*Q15+Q6*Q7*Q9*Q14*Q15,
           Q5*Q8*Q9-Q5*Q7*Q11+Q4*Q6*Q12-Q1*Q12* Q15,
           Q5*Q8^2*Q9^2-Q5*Q7^2*Q11^2-Q5*Q6*Q7* Q9*Q12+Q5*Q6*Q8*Q9*Q13+Q5*Q6^2*Q13^2-Q5* Q6*Q7*Q9*Q14+Q1*Q8*Q11*Q13*Q14-Q1*Q6* Q12*Q13*Q14-Q5*Q6*Q7*Q9*Q15-Q1*Q8*Q9* Q12*Q15+Q1*Q8*Q11*Q13*Q15+Q1*Q7*Q11*Q14* Q15,
           -Q2*Q4*Q8*Q9+Q2*Q4*Q7*Q11-Q2*Q3* Q7*Q13-2*Q1*Q5*Q8*Q13+Q2^2*Q9*Q13+Q1* Q2*Q12*Q13+Q1*Q5*Q7*Q14-Q1*Q5*Q7*Q15+Q1*Q2*Q13*Q15,
           Q4*Q8^3*Q9+Q4*Q7*Q8^2*Q11+2*Q4*Q6*Q7*Q8*Q12-2*Q3*Q7^2*Q8*Q12+2*Q1*Q7*Q8*Q12^2-2*Q4*Q6*Q8^2*Q13+Q3*Q7*Q8^2*Q13- Q2*Q6*Q7*Q12*Q13-2*Q1*Q8^2*Q12*Q13+Q2*Q6* Q8*Q13^2-Q5*Q6*Q7^2*Q14-Q4*Q6*Q7*Q8* Q14+2*Q3*Q7^2*Q8*Q14-Q1*Q7*Q8*Q12*Q14-Q1*Q7*Q8*Q12*Q15+Q1*Q7*Q8*Q14*Q15,
           -Q5*Q6*Q7*Q8*Q9-Q5*Q6*Q7^2*Q11+Q3*Q7^2*Q8*Q11-Q2*Q7^2*Q11^2-Q1*Q8^2*Q9*Q12-Q1* Q7*Q8*Q11*Q12+Q1*Q6*Q8*Q12*Q13+Q3*Q6*Q7^2*Q14+Q1*Q8^2*Q9*Q15+2*Q1*Q7*Q8*Q11*Q15-Q1*Q6*Q7*Q12*Q15,
           -Q5^2*Q7*Q9-Q4*Q5*Q8*Q9-Q4*Q5*Q7*Q11-Q4*Q5*Q6*Q13+2*Q3*Q5*Q7*Q13+Q2*Q4* Q11*Q13-Q2*Q3*Q13^2-Q4^2*Q6*Q14+Q3*Q4*Q7*Q14-Q1*Q4*Q12*Q14+2*Q1*Q5*Q13*Q15+Q1*Q4* Q14*Q15,
           -Q4*Q5*Q6*Q9*Q12+Q2*Q4*Q9*Q11*Q12+Q4*Q5*Q6*Q11*Q13-Q2*Q4*Q11^2*Q13+Q1* Q5*Q11*Q12*Q13-Q3^2*Q8*Q13^2-Q4^2*Q6*Q11*Q14+Q2*Q4*Q9*Q11*Q14-Q1*Q4*Q11*Q12*Q14+2* Q3*Q4*Q6*Q13*Q14-Q1*Q5*Q11*Q13*Q14-Q1*Q4*Q11*Q14^2+Q1*Q3*Q13*Q14*Q15,
           Q5^2*Q6*Q7^2-Q4*Q5*Q6*Q7*Q8-Q3*Q5* Q7^2*Q8+Q4^2*Q6*Q8^2+Q2*Q5*Q7^2*Q11-Q1*Q5* Q7*Q8*Q12+Q1*Q4*Q8^2*Q12+2*Q1*Q2*Q8*Q12*Q13-Q2^2*Q6*Q13^2+Q2*Q4*Q6*Q7*Q14-Q1* Q5*Q7*Q8*Q14-Q1*Q4*Q8^2*Q14-Q2*Q4*Q6*Q7*Q15+Q1*Q5*Q7*Q8*Q15-Q1*Q4*Q8^2*Q15- Q1*Q2*Q7*Q12*Q15+Q1*Q2*Q7*Q15^2,
           Q2*Q5*Q7*Q11^2-Q2*Q3*Q7*Q11*Q12+Q1*Q5*Q8*Q11*Q12+Q1*Q2*Q11*Q12^2+Q3*Q5*Q6* Q7*Q14-Q3^2*Q7*Q8*Q14-Q1*Q5*Q8*Q11*Q14-Q1*Q5*Q6*Q12*Q14+Q1*Q3*Q8*Q12*Q14-Q1* Q5*Q8*Q11*Q15+Q1*Q5*Q6*Q14*Q15-Q1*Q3*Q8*Q14*Q15,
           -Q5*Q9*Q12+Q4*Q11*Q14,
           -Q5*Q8*Q9^2*Q12-Q4*Q8*Q11^2*Q13+Q3*Q8*Q11*Q13^2+Q4*Q7*Q11^2*Q14+Q3*Q8*Q9*Q13*Q14- Q3*Q7*Q11*Q13*Q14+Q5*Q8*Q9^2*Q15-Q4*Q8*Q9*Q11*Q15+Q4*Q7*Q11^2*Q15+Q3*Q8*Q9*Q13* Q15-Q3*Q7*Q11*Q13*Q15-Q3*Q7*Q9*Q14*Q15,
           Q5^3*Q6^2*Q7-2*Q3*Q5^2*Q6*Q7*Q8+Q3^2*Q5* Q7*Q8^2-Q2*Q3*Q5*Q8^2*Q9+Q2*Q5^2*Q6*Q7*Q11-Q1*Q5^2*Q6*Q8*Q12+2*Q2*Q3*Q5*Q6* Q7*Q14-Q1*Q5^2*Q6*Q8*Q14-Q2*Q3^2*Q7*Q8*Q14+Q1*Q3*Q5*Q8^2*Q14+Q2^2*Q3*Q8*Q9*Q14- Q1*Q2*Q5*Q8*Q11*Q14+Q1*Q2*Q5*Q6*Q12* Q14-Q1*Q2*Q3*Q8*Q14^2-Q2*Q3*Q5*Q6*Q7* Q15+Q1*Q3*Q5*Q8^2*Q15+Q2^2*Q3*Q8*Q9* Q15-Q2^2*Q3*Q7*Q11*Q15+Q1*Q2*Q5*Q8*Q11*Q15-Q1*Q2*Q3*Q8*Q14*Q15,
           Q5*Q7*Q11-Q2*Q11*Q13+Q4*Q6*Q14-Q1*Q14* Q15,
           Q5*Q8*Q11*Q13-Q5*Q6*Q12*Q13-Q2*Q11* Q12*Q13-Q4*Q8*Q11*Q14+Q5*Q6*Q13*Q14+Q3* Q8*Q13*Q14-Q2*Q11*Q13*Q14-Q4*Q8*Q11*Q15+ Q3*Q8*Q13*Q15-Q3*Q7*Q14*Q15+Q2*Q9*Q14* Q15+Q1*Q12*Q14*Q15-Q1*Q14^2*Q15+Q2*Q9*Q15^ 2+Q1*Q12*Q15^2-Q1*Q14*Q15^2,
           Q4*Q8*Q9*Q11+Q4*Q7*Q11^2+Q3*Q7*Q9*Q12-Q2*Q9^2*Q12-Q1*Q9*Q12^2-Q5*Q6*Q9*Q13+ Q3*Q6*Q13^2+Q1*Q9*Q12*Q14-2*Q1*Q11*Q13*Q14-Q1*Q9*Q12*Q15+Q1*Q9*Q14*Q15,
           -Q4*Q8*Q9*Q11+Q4*Q7*Q11^2+Q3*Q7*Q9*Q12-Q2*Q9^2*Q12-Q1*Q9*Q12^2+Q5*Q6*Q9*Q13+ Q3*Q6*Q13^2+Q1*Q9*Q12*Q14-2*Q1*Q11*Q13*Q14-Q1*Q9*Q12*Q15+Q1*Q9*Q14*Q15,
           Q4*Q7*Q11^2*Q13+Q3*Q8*Q9*Q13^2-2*Q3*Q7* Q11*Q13^2+Q2*Q9*Q11*Q13^2+Q3*Q6*Q13^3-Q4*Q7* Q9*Q11*Q14+2*Q3*Q7*Q9*Q13*Q14-2*Q2*Q9^2*Q13*Q14-Q1*Q9*Q12*Q13*Q14-2*Q1*Q11*Q13^2*Q14+2* Q1*Q9*Q13*Q14^2-Q5*Q7*Q9^2*Q15-Q3*Q7*Q9* Q13*Q15+2*Q2*Q9^2*Q13*Q15+Q1*Q9*Q12*Q13* Q15-Q1*Q9*Q13*Q14*Q15,
           Q5^2*Q6*Q9+Q3*Q4*Q8*Q11-2*Q2*Q5*Q9* Q11+Q2*Q4*Q11^2+Q2*Q3*Q9*Q12-Q3^2*Q8*Q13-Q3*Q4*Q6*Q14+2*Q1*Q5*Q11*Q14-Q1*Q3*Q12* Q14-2*Q1*Q5*Q11*Q15+Q1*Q3*Q12*Q15,
           Q5*Q8*Q9^2*Q12+Q4*Q8*Q9*Q11*Q12+Q4*Q8*Q11^2*Q13-Q3*Q8*Q9*Q12*Q13-Q4*Q6*Q11* Q12*Q13-Q3*Q8*Q11*Q13^2+Q3*Q6*Q12*Q13^2+Q4*Q8*Q9*Q11*Q14-Q4*Q6*Q9*Q12*Q14-Q4*Q6* Q11*Q13*Q14+Q3*Q6*Q13^2*Q14-Q5*Q8*Q9^2*Q15,
           Q2*Q9*Q12+Q5*Q6*Q13-Q3*Q8*Q13-Q1*Q12* Q14,
           Q3*Q8^3*Q9+Q4*Q6*Q8^2*Q11-2*Q3*Q7* Q8^2*Q11+Q2*Q7*Q8*Q11^2+Q3*Q6*Q8^2*Q13-Q5* Q6^2*Q7*Q14+2*Q4*Q6^2*Q8*Q14-Q3*Q6*Q7*Q8*Q14+Q1*Q6*Q8*Q12*Q14-2*Q4*Q6^2*Q8*Q15+2* Q3*Q6*Q7*Q8*Q15-Q2*Q6*Q7*Q11*Q15-2*Q1* Q8^2*Q11*Q15-Q1*Q6*Q8*Q12*Q15-Q1*Q6*Q8* Q14*Q15+2*Q1*Q6*Q8*Q15^2,
           Q3*Q7*Q9-Q2*Q9^2-Q1*Q11*Q13+Q1*Q9*Q14,
           -Q4*Q5*Q6^2+Q3*Q4*Q6*Q8+Q2*Q5*Q6* Q9-2*Q1*Q5*Q8*Q11-Q2*Q3*Q6*Q13+Q1*Q3*Q8*Q14+Q1*Q2*Q11*Q14+Q1*Q5*Q6*Q15-Q1* Q3*Q8*Q15,
           -Q3*Q4*Q8*Q9+Q3*Q4*Q7*Q11-Q1*Q4*Q11*Q12-Q3^2*Q7*Q13+Q2*Q3*Q9*Q13-2* Q1*Q5*Q11*Q13+Q1*Q3*Q12*Q13+Q1*Q5*Q9*Q14+Q1*Q4*Q11*Q14,
           Q2*Q11*Q13-Q3*Q7*Q14,
           -Q5*Q8*Q9-Q5*Q7*Q11+Q2*Q9*Q12+Q3*Q8* Q13+Q2*Q11*Q13-Q4*Q6*Q14+Q2*Q9*Q14-Q1* Q12*Q14-Q1*Q14^2+Q1*Q12*Q15,
           -Q5*Q8*Q9+Q5*Q7*Q11+Q2*Q9*Q12-Q3*Q8* Q13+Q2*Q11*Q13-Q4*Q6*Q14+Q2*Q9*Q14-Q1* Q12*Q14-Q1*Q14^2+Q1*Q12*Q15,
           Q5*Q8*Q9-Q5*Q7*Q11+Q2*Q9*Q12+Q3*Q8* Q13-Q2*Q11*Q13-Q4*Q6*Q14+Q2*Q9*Q14-Q1* Q12*Q14-Q1*Q14^2+Q1*Q12*Q15,
           Q5*Q8*Q9+Q5*Q7*Q11+Q2*Q9*Q12-Q3*Q8* Q13-Q2*Q11*Q13-Q4*Q6*Q14+Q2*Q9*Q14-Q1* Q12*Q14-Q1*Q14^2+Q1*Q12*Q15,
           -Q2*Q4*Q5*Q7*Q11+Q2*Q3*Q4*Q8*Q13-Q1* Q5^2*Q8*Q13-Q2^2*Q3*Q13^2-Q1*Q5^2*Q7*Q14+ Q1*Q4*Q5*Q8*Q14+Q1*Q2*Q4*Q12*Q14+Q1*Q5^2*Q7*Q15-Q1*Q2*Q4*Q12*Q15+2*Q1*Q2*Q5* Q13*Q15-Q1*Q2*Q4*Q14*Q15,
           Q5*Q11-Q3*Q14,
           -Q5^2*Q6*Q7+2*Q4*Q5*Q6*Q8-Q3*Q4*Q8^2-Q2*Q5*Q8*Q9-Q2*Q5*Q7*Q11+Q2*Q4*Q6* Q12-Q2^2*Q9*Q12-Q2*Q5*Q6*Q13+Q2*Q3*Q8*Q13+2*Q1*Q5*Q8*Q14+Q1*Q2*Q12*Q14-Q1*Q2* Q12*Q15,
           -2*Q5*Q8^2*Q11+Q5*Q6*Q8*Q12+Q3*Q8^ 2*Q12+Q2*Q8*Q11*Q12+Q2*Q6*Q12^2-Q5*Q6* Q8*Q14-Q2*Q8*Q11*Q15,
           -Q4*Q8*Q9^2*Q12+Q4*Q7*Q11^2*Q13-Q4*Q6*Q9*Q12*Q13+2*Q2*Q9^2*Q12*Q13+Q5*Q6*Q9*Q13^2-2* Q4*Q6*Q11*Q13^2+Q2*Q9*Q11*Q13^2+Q3*Q6*Q13^3-Q4*Q7*Q9*Q11*Q14+2*Q4*Q6*Q9*Q13*Q14-2* Q2*Q9^2*Q13*Q14-Q1*Q9*Q12*Q13*Q14-2*Q1*Q11* Q13^2*Q14+2*Q1*Q9*Q13*Q14^2+Q1*Q9*Q12*Q13*Q15-Q1*Q9*Q13*Q14*Q15,
           Q5^2*Q7^2*Q9-Q4^2*Q8^2*Q9+Q4*Q5*Q7^2*Q11- Q3*Q5*Q7^2*Q13-Q2*Q5*Q7*Q9*Q13-Q1*Q5* Q7*Q12*Q13+2*Q1*Q4*Q8*Q12*Q13+Q2^2*Q9*Q13^2+Q1*Q2*Q12*Q13^2-Q2*Q4*Q7*Q9*Q14-Q1* Q4*Q7*Q12*Q14+Q1*Q5*Q7*Q13*Q14-Q1*Q2*Q13^2*Q14+Q1*Q4*Q7*Q14^2+Q2*Q4*Q7*Q9*Q15- Q1*Q5*Q7*Q13*Q15-Q1*Q2*Q13^2*Q15,
           Q5^2*Q6*Q9*Q13-Q1*Q5*Q11*Q12*Q13-Q4^2* Q6*Q11*Q14-Q1*Q4*Q11*Q12*Q14+Q3*Q4*Q6*Q13*Q14-Q1*Q5*Q11*Q13*Q14+Q1*Q3*Q12*Q13*Q14-Q4*Q5*Q6*Q9*Q15+Q1*Q5*Q11*Q13*Q15+Q1*Q4*Q11*Q14*Q15-Q1*Q3*Q13*Q14*Q15+Q1*Q5*Q9*Q15^2,
           -Q5*Q8*Q9*Q12+Q5*Q8*Q11*Q13- Q5*Q6*Q12*Q13+Q5*Q7*Q11*Q14-Q2*Q11*Q13*Q14+Q4*Q6*Q14^2+Q1*Q12*Q14^2+Q5*Q8*Q9*Q15+ Q5*Q7*Q11*Q15-Q5*Q6*Q13*Q15-Q2*Q11*Q13*Q15+Q4*Q6*Q14*Q15-Q3*Q7*Q14*Q15+Q1*Q12* Q14*Q15-Q1*Q14^2*Q15-Q1*Q14*Q15^2,
           Q5^2*Q7^2*Q11+Q4*Q5*Q7*Q8*Q11+Q4^2*Q8^2* Q11-Q2*Q4*Q7*Q11*Q12+Q1*Q5*Q8*Q12*Q13-Q2^2*Q11*Q13^2-Q2*Q4*Q7*Q11*Q14-Q1*Q4*Q8* Q12*Q14-Q2*Q4*Q7*Q11*Q15+Q1*Q5*Q8*Q13*Q15+Q1*Q2*Q12*Q13*Q15-Q1*Q5*Q7*Q14*Q15,
           -Q3^2*Q8^2*Q9+Q3^2*Q7*Q8*Q11+Q1*Q5*Q6*Q11*Q12-Q1*Q3*Q8*Q11*Q12-Q1*Q2*Q11^2*Q12- Q3^2*Q6*Q8*Q13-Q2*Q3*Q6*Q11*Q13+Q3^2*Q6* Q7*Q14+2*Q1*Q3*Q8*Q11*Q15+Q1*Q2*Q11^2*Q15-Q1*Q3*Q6*Q12*Q15,
           Q5*Q8*Q9-Q5*Q7*Q11-Q5*Q6*Q13+Q3*Q7* Q15-Q1*Q12*Q15+Q1*Q14*Q15,
           Q2*Q7*Q8*Q9*Q11-Q2*Q7^2*Q11^2+Q2*Q6* Q7*Q9*Q12-Q1*Q8^2*Q9*Q12-Q2*Q6*Q7*Q11*Q13-2*Q1*Q8^2*Q11*Q13-Q1*Q7*Q8*Q11*Q14+Q1* Q8^2*Q9*Q15+3*Q1*Q7*Q8*Q11*Q15+Q1*Q6*Q8* Q13*Q15-Q1*Q6*Q7*Q14*Q15,
           -Q3*Q5*Q6*Q8*Q9-Q3^2*Q8^2*Q9+Q2*Q5*Q6*Q9*Q11+2*Q2*Q3*Q8*Q9*Q11-Q2^2*Q9*Q11^2+ Q2*Q3*Q6*Q9*Q12-Q1*Q2*Q11^2*Q12-Q3^2*Q6* Q8*Q13+Q1*Q3*Q8*Q11*Q14+Q1*Q3*Q8*Q11* Q15-Q1*Q2*Q11^2*Q15,
           Q2*Q4*Q6-Q1*Q5*Q8+Q1*Q2*Q12-Q1*Q2* Q15,
           Q5^2*Q9^2*Q12+3*Q4*Q5*Q9*Q11*Q12+Q4^2* Q11^2*Q12-2*Q5^2*Q9*Q11*Q13-2*Q4*Q5*Q11^2*Q13+Q3*Q5*Q9*Q12*Q13+Q3*Q4*Q11*Q12*Q13-Q4^2* Q11^2*Q14-Q3*Q4*Q9*Q12*Q14-Q5^2*Q9^2*Q15-Q3* Q4*Q9*Q12*Q15+Q3*Q4*Q9*Q14*Q15,
           -Q5*Q8*Q9*Q11-Q4*Q8*Q11^2+Q5*Q6*Q9*Q12+Q3*Q8*Q11*Q13+Q4*Q6*Q11*Q14+Q1*Q11*Q12*Q14-Q3*Q6*Q13*Q14+Q3*Q8*Q9*Q15+Q4*Q6*Q11*Q15-Q3*Q7*Q11*Q15-Q1*Q11*Q14*Q15-Q1*Q11*Q15^2,
           -Q4*Q11*Q14+Q5*Q9*Q15,
           Q4*Q11*Q14-Q5*Q9*Q15,
           -Q3*Q4*Q8*Q9*Q11+Q1*Q5*Q9*Q11*Q12-Q1* Q4*Q11^2*Q12-Q3^2*Q8*Q9*Q13+Q3^2*Q7*Q11* Q13-Q1*Q3*Q11*Q12*Q13-Q3^2*Q6*Q13^2+Q1*Q4*Q11^2*Q14-Q1*Q3*Q9*Q12*Q14+2*Q1*Q3*Q11*Q13* Q14+Q3^2*Q7*Q9*Q15,
           Q5^2*Q7*Q8*Q12-Q4*Q5*Q8^2*Q12-Q2^2*Q12^2* Q13-Q5^2*Q7*Q8*Q14+Q4*Q5*Q8^2*Q14+Q2*Q5* Q7*Q12*Q14-Q5^2*Q7*Q8*Q15+Q4*Q5*Q8^2*Q15+Q2*Q5*Q7*Q12*Q15-Q2*Q4*Q8*Q12*Q15+Q2^2*Q12*Q13*Q15-Q2*Q5*Q7*Q14*Q15,
           -Q3*Q5*Q7^2+Q3*Q4*Q7*Q8+Q2*Q5*Q7*Q9-Q2*Q4*Q7*Q11+Q1*Q5*Q7*Q12-Q1*Q4*Q8*Q12-2*Q1*Q5*Q8*Q13+Q1*Q4*Q8*Q14+Q1* Q2*Q13*Q14,
           -Q5^2*Q8*Q9+Q4*Q5*Q8*Q11+Q4*Q5*Q6*Q12-Q3*Q5*Q7*Q12-Q2*Q4*Q11*Q12+ Q1*Q5*Q12^2-Q3*Q5*Q8*Q13-Q3*Q5*Q7*Q14+Q3*Q4*Q8*Q14+Q1*Q5*Q12*Q14+Q2*Q3*Q13* Q15-Q1*Q5*Q14*Q15,
           Q5*Q6*Q12-Q3*Q8*Q14,
           Q5*Q7^2*Q8*Q9+Q4*Q7*Q8^2*Q9-Q5*Q7^3* Q11-Q4*Q7^2*Q8*Q11+Q5*Q6*Q7^2*Q13-Q1*Q8^2* Q13^2-Q1*Q7*Q8*Q13*Q15+Q1*Q7^2*Q14*Q15,
           Q4*Q6*Q12^2+Q5*Q8*Q11*Q13-Q3*Q8*Q12*Q13+Q2*Q11*Q12*Q13-Q5*Q7*Q11*Q14-Q4*Q8*Q11*Q14+Q1*Q12^2*Q14+Q5*Q7*Q11*Q15-Q4*Q8* Q11*Q15+Q4*Q6*Q12*Q15-Q2*Q9*Q12*Q15-Q1*Q12^2*Q15-Q3*Q8*Q13*Q15+Q2*Q11*Q13*Q15+Q1* Q12*Q14*Q15-Q1*Q12*Q15^2,
           Q4*Q5*Q6*Q7-Q3*Q5*Q7^2-Q2*Q4*Q7*Q11+Q1*Q5*Q7*Q12+Q2*Q3*Q7*Q13-2*Q1*Q5* Q8*Q13-Q1*Q2*Q12*Q13+Q1*Q4*Q8*Q15+Q1*Q2*Q13*Q15,
           Q5^2*Q6*Q9-2*Q4*Q5*Q6*Q11+Q2*Q4*Q11^2+Q3*Q4*Q6*Q12-Q3^2*Q8*Q13+Q2*Q3* Q11*Q13-2*Q1*Q5*Q11*Q14+Q1*Q3*Q12*Q14-Q2* Q3*Q9*Q15+2*Q1*Q5*Q11*Q15-Q1*Q3*Q12*Q15,
           -Q4*Q5*Q7*Q9*Q11-Q4^2*Q7*Q11^2+Q3*Q4*Q7*Q11*Q13-Q1*Q4*Q11*Q12*Q13-2*Q1*Q5*Q11* Q13^2-Q1*Q4*Q9*Q12*Q14+Q1*Q5*Q9*Q13*Q14+3*Q1*Q4*Q11*Q13*Q14+Q1*Q3*Q13^2*Q14+Q3*Q4* Q7*Q9*Q15-Q1*Q3*Q13^2*Q15,
           Q3*Q7^2*Q11^2*Q13-Q2*Q7*Q9*Q11^2*Q13+Q3* Q6*Q8*Q9*Q13^2-2*Q3*Q6*Q7*Q11*Q13^2-Q1*Q6* Q11*Q12*Q13^2+Q3*Q6^2*Q13^3-Q3*Q7^2*Q9*Q11*Q14+Q2*Q7*Q9^2*Q11*Q14+2*Q3*Q6*Q7*Q9*Q13* Q14-Q1*Q8*Q9*Q11*Q13*Q14+Q1*Q7*Q11^2*Q13* Q14+Q1*Q6*Q9*Q12*Q13*Q14-Q1*Q6*Q11*Q13^2* Q14-Q1*Q7*Q9*Q11*Q14^2-Q3*Q7*Q8*Q9^2* Q15+Q2*Q7*Q9^2*Q11*Q15-Q3*Q6*Q7*Q9*Q13*Q15+Q1*Q8*Q9*Q11*Q13*Q15+Q1*Q7*Q11^2*Q13* Q15-Q1*Q7*Q9*Q11*Q14*Q15,
           -Q5*Q7*Q9*Q11-Q4*Q8*Q9*Q11-Q4*Q7*Q11^ 2+Q5*Q6*Q9*Q13+2*Q3*Q7*Q11*Q13-Q2*Q9* Q11*Q13-Q3*Q6*Q13^2+2*Q1*Q11*Q13*Q14+Q3*Q7*Q9*Q15-Q2*Q9^2*Q15-Q1*Q9*Q12*Q15+Q1*Q9* Q14*Q15,
           -Q5^2*Q7^2*Q12+Q4^2*Q8^2*Q12-Q2*Q4*Q7*Q12^2+Q2*Q4*Q8*Q12*Q13+Q2^2*Q12*Q13^2+Q4* Q5*Q7*Q8*Q14-Q4^2*Q8^2*Q14-Q2*Q4*Q7*Q12* Q14+Q5^2*Q7^2*Q15-Q4*Q5*Q7*Q8*Q15-Q2^2* Q13^2*Q15+Q2*Q4*Q7*Q14*Q15,
           Q2*Q4*Q9^2*Q12+Q2*Q5*Q9^2*Q13-Q2*Q4* Q9*Q11*Q13-Q2*Q3*Q9*Q13^2+2*Q1*Q5*Q11*Q13^2-Q1*Q4*Q9*Q12*Q14-Q1*Q5*Q9*Q13*Q14+Q1* Q3*Q13^2*Q14+Q1*Q4*Q9*Q12*Q15-Q1*Q4*Q11*Q13*Q15-Q1*Q3*Q13^2*Q15,
           Q8*Q9*Q12^2-Q8*Q11*Q12*Q13+Q6*Q12^2*Q13+ Q7*Q11*Q12*Q14-Q8*Q11*Q13*Q14+Q7*Q11*Q12*Q15-Q8*Q11*Q13*Q15-Q7*Q11*Q14*Q15,
           Q4*Q5*Q6-Q3*Q4*Q8+Q2*Q4*Q11-Q1*Q5* Q15,
           Q4*Q5*Q6+Q3*Q4*Q8-Q2*Q4*Q11-Q1* Q5*Q15,
           Q4*Q8*Q9*Q12-Q5*Q8*Q9*Q13+Q5* Q7*Q11*Q13-Q5*Q6*Q13^2-Q4*Q7*Q11*Q14+Q4* Q6*Q13*Q14+Q1*Q12*Q13*Q14+Q5*Q7*Q9*Q15+ Q4*Q6*Q13*Q15-Q3*Q7*Q13*Q15-Q1*Q13*Q14* Q15-Q1*Q13*Q15^2,
           -Q5^2*Q7*Q9^2-Q4*Q5*Q8*Q9^2-Q4*Q5*Q7*Q9*Q11+2*Q3*Q5*Q7*Q9*Q13+Q3*Q4*Q7*Q11* Q13+Q1*Q5*Q9*Q12*Q13-Q3^2*Q7*Q13^2+Q3*Q4* Q7*Q9*Q14-Q1*Q3*Q13^2*Q14+Q1*Q5*Q9*Q13* Q15-Q1*Q3*Q13^2*Q15,
           -Q5^2*Q7*Q11+Q2*Q5*Q9*Q12-Q2*Q4*Q11*Q12+Q3*Q5*Q8*Q13+2*Q4*Q5*Q6*Q14-Q3*Q4* Q8*Q14-2*Q2*Q5*Q9*Q14+2*Q1*Q5*Q14^2-Q1*Q5*Q14*Q15,
           Q3*Q5*Q8^2*Q9+Q2*Q5*Q6*Q9*Q12- Q1*Q5*Q8*Q11*Q12-Q2^2*Q9*Q11*Q12-Q2*Q3* Q8*Q9*Q14+Q1*Q5*Q8*Q11*Q14-Q1*Q5*Q6* Q12*Q14+Q1*Q2*Q11*Q12*Q14+Q1*Q3*Q8*Q14^2- Q1*Q5*Q8*Q11*Q15+Q1*Q5*Q6*Q12*Q15-Q1* Q2*Q11*Q12*Q15,
           Q5^2*Q8*Q9^2-Q4^2*Q8*Q11^2-Q3*Q4*Q8*Q9* Q12+Q3*Q5*Q8*Q9*Q13+Q1*Q5*Q11*Q12*Q13+Q3^2*Q8*Q13^2-Q3*Q4*Q8*Q9*Q14+Q1*Q4*Q11* Q12*Q14+Q1*Q5*Q11*Q13*Q14-Q3*Q4*Q8*Q9*Q15-Q1*Q5*Q9*Q12*Q15-Q1*Q3*Q13*Q14*Q15,
           Q5*Q7*Q8*Q9*Q12+Q5*Q8^2*Q9*Q13-Q2*Q8*Q9*Q12*Q13-Q5*Q6*Q8*Q13^2+Q2*Q6*Q12* Q13^2-Q5*Q7^2*Q11*Q14+Q5*Q7*Q8*Q9*Q15+Q5*Q7^2*Q11*Q15-Q2*Q7*Q9*Q12*Q15-Q5*Q6*Q7* Q13*Q15-Q2*Q8*Q9*Q13*Q15+Q2*Q6*Q13^2*Q15,
           -Q5*Q8*Q9+Q4*Q8*Q11-Q3*Q8*Q13+Q1*Q15^ 2,
           -Q5*Q7*Q11-Q4*Q8*Q11+Q2*Q11*Q13+Q1* Q14^2,
           -Q2*Q3*Q8*Q9-Q2*Q4*Q6*Q11-2*Q1*Q5*Q8*Q11+Q2^2*Q9*Q11-Q1*Q5*Q6*Q12+Q1*Q2* Q11*Q12+Q2*Q3*Q6*Q13+Q1*Q5*Q6*Q14+Q1*Q2*Q11*Q15,
           Q3*Q8^2*Q9*Q12-Q2*Q8*Q9*Q11* Q12-Q3*Q8^2*Q11*Q13+Q2*Q8*Q11^2*Q13+Q5*Q6^2*Q12*Q13-Q3*Q6*Q8*Q12*Q13+Q2*Q6*Q11* Q12*Q13-Q5*Q6^2*Q13*Q14+Q3*Q8^2*Q9*Q15-Q2*Q8*Q9*Q11*Q15-Q2*Q6*Q9*Q12*Q15+Q2*Q6* Q11*Q13*Q15,
           Q3*Q6*Q7-Q1*Q8*Q11-Q1*Q6* Q12+Q1*Q6*Q14,
           -Q8^2*Q9-Q7*Q8*Q11+Q6*Q8*Q13+Q6*Q7*Q14,
           Q8^2*Q9+Q7*Q8*Q11-Q6*Q8*Q13-Q6*Q7*Q14,
           Q8^2*Q9-Q7*Q8*Q11+Q6*Q8*Q13-Q6*Q7*Q14,
           -Q8^2*Q9+Q7*Q8*Q11-Q6*Q8*Q13+Q6*Q7* Q14,
           -Q8^2*Q9+Q7*Q8*Q11+Q6*Q8*Q13-Q6*Q7* Q14,
           -Q5*Q7*Q8*Q9^2-Q4*Q8^2*Q9^2-Q5*Q6*Q7*Q9*Q13+Q2*Q8*Q9^2*Q13+2*Q1*Q8*Q9*Q12* Q13+Q1*Q6*Q12*Q13^2-Q1*Q7*Q9*Q12*Q14-Q1* Q8*Q9*Q13*Q14+Q1*Q7*Q11*Q13*Q14-Q1*Q6* Q13^2*Q14+Q2*Q7*Q9^2*Q15,
           -Q4*Q14^2+Q5*Q13*Q15,
           -Q5*Q8*Q9^3+Q5*Q7*Q9^2*Q11+Q5*Q6*Q9^2*Q13-Q3*Q8*Q9^2*Q13-Q1*Q9*Q11*Q12*Q13+Q3* Q6*Q9*Q13^2-Q1*Q11^2*Q13^2+Q1*Q9^2*Q12*Q15,
           Q4^2*Q8^2*Q12^2+Q5^2*Q7*Q8*Q12*Q13+Q2*Q4* Q8*Q12^2*Q13-2*Q5^2*Q8^2*Q13^2+Q2*Q5*Q8*Q12*Q13^2-Q5^2*Q7^2*Q12*Q14-3*Q4*Q5*Q7*Q8*Q12*Q14+ Q2*Q4*Q7*Q12^2*Q14+2*Q5^2*Q7*Q8*Q13*Q14+2*Q4*Q5*Q8^2*Q13*Q14-Q2*Q5*Q7*Q12*Q13*Q14- Q2*Q4*Q8*Q12*Q13*Q14+Q5^2*Q7^2*Q12*Q15-Q2* Q4*Q7*Q12^2*Q15-2*Q5^2*Q7*Q8*Q13*Q15-Q4*Q5*Q8^2*Q13*Q15+Q2*Q5*Q8*Q13^2*Q15-Q2^2*Q12* Q13^2*Q15+Q5^2*Q7^2*Q14*Q15+Q2*Q5*Q7*Q13*Q15^2,
           -Q5^2*Q7^2*Q12+Q4^2*Q8^2*Q12-Q2*Q4*Q7*Q12^2+Q2*Q4*Q8*Q12*Q13+Q2^2*Q12*Q13^2+Q5^2*Q7^2* Q14-Q4^2*Q8^2*Q14-Q2*Q5*Q7*Q13*Q14-Q2*Q4*Q7*Q12*Q15+Q2*Q5*Q7*Q13*Q15-Q2^2*Q13^2* Q15+Q2*Q4*Q7*Q14*Q15,
           Q4*Q7*Q8*Q9*Q11-Q4*Q7^2*Q11^2-Q5*Q6* Q7*Q9*Q13+Q1*Q8*Q9*Q12*Q13-Q1*Q8*Q11*Q13^2-Q1*Q6*Q12*Q13^2-Q1*Q7*Q9*Q12*Q14+2*Q1* Q7*Q11*Q13*Q14+Q1*Q6*Q13^2*Q14+Q1*Q7*Q9*Q12*Q15-Q1*Q7*Q9*Q14*Q15,
           Q4*Q5*Q7*Q8^2*Q9+Q4^2*Q8^3*Q9+2*Q2*Q4* Q7*Q8*Q9*Q12-2*Q2*Q4*Q8^2*Q9*Q13+Q2*Q3* Q7^2*Q12*Q13-Q1*Q5*Q7*Q8*Q12*Q13-Q1*Q4* Q8^2*Q12*Q13-Q2^2*Q7*Q9*Q12*Q13-Q1*Q2* Q7*Q12^2*Q13-Q2*Q3*Q7*Q8*Q13^2+Q2^2*Q8*Q9*Q13^2+Q1*Q2*Q8*Q12*Q13^2+Q1*Q4*Q7*Q8*Q12* Q14-Q1*Q4*Q8^2*Q13*Q14-Q2*Q5*Q7^2*Q9*Q15-Q2*Q4*Q7*Q8*Q9*Q15+Q2*Q3*Q7^2*Q13* Q15+Q1*Q5*Q7*Q8*Q13*Q15-Q1*Q2*Q7*Q12*Q13*Q15+Q1*Q2*Q8*Q13^2*Q15,
           -2*Q8^2*Q11*Q13-Q6*Q8*Q12*Q13+Q8^2*Q9*Q14+Q7*Q8*Q11*Q14+Q6*Q8*Q13*Q14+Q6*Q7*Q14^2- Q7*Q8*Q11*Q15,
           -Q4*Q5*Q7*Q9*Q11-Q4^2*Q8*Q9*Q11-Q4^2*Q7*Q11^2+Q3*Q5*Q7*Q9*Q13+2*Q3*Q4*Q7*Q11* Q13+Q1*Q4*Q11*Q12*Q13-Q3^2*Q7*Q13^2+Q1*Q4* Q11*Q13*Q14-Q1*Q3*Q13^2*Q14+Q3*Q4*Q7*Q9* Q15-Q1*Q3*Q13^2*Q15,
           Q8^2*Q9*Q11*Q12*Q13-Q7*Q8*Q11^2*Q12*Q13+ Q6*Q8*Q9*Q12^2*Q13-2*Q8^2*Q11^2*Q13^2-2*Q6*Q8*Q11*Q12*Q13^2+2*Q7*Q8*Q11^2*Q13*Q14+2*Q6*Q8*Q11* Q13^2*Q14+Q6^2*Q12*Q13^2*Q14-Q8^2*Q9^2*Q12*Q15+Q8^2*Q9*Q11*Q13*Q15+Q6*Q8*Q11*Q13^2*Q15+Q6^2* Q12*Q13^2*Q15-Q7*Q8*Q9*Q11*Q14*Q15-Q6*Q8* Q9*Q13*Q14*Q15-3*Q6*Q7*Q11*Q13*Q14*Q15-Q6^2*Q13^2*Q14*Q15+Q7*Q8*Q9*Q11*Q15^2+Q7^2*Q11^2* Q15^2-Q6*Q7*Q9*Q12*Q15^2+Q6*Q7*Q9*Q14*Q15^2,
           Q12-Q15,
           -Q12+Q15,
           Q5^2*Q6^2*Q9-Q3^2*Q8^2*Q9-Q4*Q5*Q6^2*Q11- Q2*Q5*Q6*Q9*Q11+Q2^2*Q9*Q11^2+Q2*Q3*Q6* Q9*Q12-Q1*Q5*Q6*Q11*Q12-Q1*Q2*Q11^2*Q12+Q3*Q5*Q6^2*Q13-Q2*Q3*Q6*Q9*Q14+Q1* Q5*Q6*Q11*Q14-Q1*Q2*Q11^2*Q14+Q1*Q3*Q6*Q14^2-Q1*Q5*Q6*Q11*Q15+2*Q1*Q3*Q8*Q11*Q15+ Q1*Q2*Q11^2*Q15-Q1*Q3*Q6*Q14*Q15,
           Q5^2*Q7^2*Q11-Q4^2*Q8^2*Q11-Q2*Q4*Q7*Q11* Q12+Q2*Q5*Q7*Q11*Q13+Q1*Q5*Q8*Q12*Q13+Q2^2*Q11*Q13^2-Q2*Q4*Q7*Q11*Q14+Q1*Q4*Q8* Q12*Q14+Q1*Q5*Q8*Q13*Q14-Q2*Q4*Q7*Q11*Q15-Q1*Q2*Q12*Q13*Q15-Q1*Q5*Q7*Q14*Q15,
           -Q5^2*Q7*Q9^2-Q4*Q5*Q8*Q9^2+Q2*Q4*Q9^2*Q12-Q3*Q4*Q8*Q9*Q13+Q2*Q5*Q9^2*Q13-Q1* Q5*Q9*Q13*Q14+Q1*Q4*Q11*Q13*Q14-Q1*Q3*Q13^2*Q14+2*Q1*Q5*Q9*Q13*Q15+Q1*Q3*Q13^2*Q15-Q1*Q4*Q9*Q14*Q15,
           -2*Q5^2*Q7*Q8*Q13-2*Q2*Q5*Q8*Q13^2-Q2^2*Q12* Q13^2-Q5^2*Q7^2*Q14+Q2*Q4*Q7*Q12*Q14+Q5^2*Q7^2*Q15+Q4*Q5*Q7*Q8*Q15-Q2*Q4*Q7*Q12* Q15+3*Q2*Q5*Q7*Q13*Q15+Q2*Q4*Q8*Q13*Q15+Q2^2*Q13^2*Q15-Q2*Q4*Q7*Q14*Q15,
           Q5^2*Q9^2*Q12-Q4^2*Q11^2*Q12-Q3*Q5*Q9*Q12* Q13-Q5^2*Q9^2*Q14+Q4^2*Q11^2*Q14+Q3*Q4*Q11*Q13*Q14+Q3^2*Q13^2*Q14-Q3*Q4*Q9*Q14^2+Q3*Q4* Q9*Q12*Q15+Q3*Q5*Q9*Q13*Q15-Q3^2*Q13^2*Q15- Q3*Q4*Q9*Q14*Q15,
           Q5^2*Q9*Q11*Q12+2*Q5^2*Q11^2*Q13-Q3*Q5*Q11* Q12*Q13-Q3*Q4*Q11*Q12*Q14-Q5^2*Q9*Q11*Q15-Q4*Q5*Q11^2*Q15+Q3*Q4*Q11*Q12*Q15-2*Q3*Q5* Q11*Q13*Q15+Q3*Q4*Q11*Q14*Q15+Q3^2*Q13*Q15^ 2,
           -Q5^2*Q6^2*Q9-Q3*Q5*Q6*Q8*Q9+2*Q2*Q5* Q6*Q9*Q11+Q2*Q3*Q8*Q9*Q11-Q2^2*Q9*Q11^2+ Q1*Q5*Q6*Q11*Q12-Q1*Q2*Q11^2*Q12-Q3*Q5*Q6^2*Q13+Q1*Q5*Q6*Q11*Q14+Q2*Q3*Q6* Q9*Q15-Q1*Q2*Q11^2*Q15,
           Q5*Q6*Q9+Q3*Q8*Q9-Q2*Q9*Q11-Q3*Q6* Q13,
           Q5^2*Q6*Q9-2*Q2*Q5*Q9*Q11+Q2*Q4* Q11^2+Q2*Q3*Q9*Q12+Q3*Q5*Q6*Q13-Q3^2*Q8*Q13-Q3*Q4*Q6*Q14+2*Q1*Q5*Q11*Q14-Q1*Q3* Q12*Q14-2*Q1*Q5*Q11*Q15+Q1*Q3*Q12*Q15,
           -Q5^2*Q6^2*Q7+2*Q3*Q5*Q6*Q7*Q8-Q3^2*Q7*Q8^2-Q2*Q5*Q6*Q7*Q11+Q2*Q3*Q7*Q8*Q11+ Q1*Q5*Q6*Q8*Q12-Q2*Q5*Q6^2*Q13+Q1*Q5* Q6*Q8*Q14-Q1*Q3*Q8^2*Q14+Q2*Q3*Q6*Q7* Q15-Q1*Q3*Q8^2*Q15,
           -Q5*Q6^2*Q7*Q11+Q2*Q6*Q7*Q11^2-Q1*Q8^2*Q11^2-Q1*Q6*Q8*Q11*Q12-Q5*Q6^3*Q13+Q3*Q6^2* Q8*Q13+Q2*Q6^2*Q11*Q13+Q1*Q6^2*Q12*Q14,
           -Q6*Q14^2+Q8*Q11*Q15,
           -Q4*Q5^2*Q6*Q7+Q4^2*Q5*Q6*Q8+Q2*Q4^2*Q6*Q12-Q2*Q4*Q5*Q6*Q13+2*Q1*Q5^2*Q8*Q13- Q1*Q5^2*Q7*Q14+Q1*Q2*Q4*Q12*Q14-Q1*Q2* Q5*Q13*Q14+Q1*Q5^2*Q7*Q15-Q1*Q4*Q5*Q8* Q15-Q1*Q2*Q4*Q12*Q15,
           Q4*Q11*Q12-Q3*Q13*Q15,
           -Q4*Q11*Q12+Q3*Q13*Q15,
           Q5*Q7*Q11-Q5*Q6*Q13+Q2*Q9*Q15-Q1*Q14* Q15,
           Q5^2*Q9-Q3*Q4*Q14,
           -Q5*Q6*Q12+Q2*Q11*Q15,
           Q5*Q7*Q11-Q1*Q14^2,
           Q5^2*Q6*Q7-2*Q4*Q5*Q6*Q8+Q3*Q4*Q8^2-2* Q1*Q5*Q8*Q12+Q2*Q3*Q8*Q13-Q2^2*Q11*Q13+Q2*Q4*Q6*Q14+Q1*Q2*Q12*Q14-Q2*Q3*Q7* Q15+2*Q1*Q5*Q8*Q15-Q1*Q2*Q14*Q15,
           Q5*Q8*Q9*Q11-Q5*Q7*Q11^2+Q5*Q6*Q9*Q12-Q4*Q6*Q11*Q12-Q5*Q6*Q11*Q13-Q1*Q11*Q12*Q14+Q3*Q6*Q13*Q14-Q3*Q8*Q9*Q15-Q4*Q6*Q11*Q15+Q2*Q9*Q11*Q15+Q1*Q11*Q12*Q15+Q1*Q11*Q15^2,
           -Q5*Q6*Q8*Q9-Q3*Q8^2*Q9-Q4* Q6*Q8*Q11+2*Q3*Q7*Q8*Q11-Q2*Q7*Q11^2-Q3* Q6*Q8*Q13+Q2*Q6*Q11*Q13-Q4*Q6^2*Q14+Q3* Q6*Q7*Q14-Q1*Q6*Q12*Q14+2*Q1*Q8*Q11*Q15+Q1*Q6*Q14*Q15,
           -Q4*Q8*Q11*Q12+Q5*Q8*Q11*Q13-Q2*Q11*Q12* Q13-Q5*Q8*Q9*Q14+Q5*Q7*Q11*Q14+Q4*Q6* Q14^2+Q1*Q12*Q14^2-Q5*Q8*Q9*Q15+Q5*Q7*Q11*Q15-Q4*Q8*Q11*Q15+Q2*Q11*Q13*Q15+Q4*Q6*Q14*Q15-Q3*Q7*Q14*Q15+Q1*Q12*Q14*Q15-Q1*Q14^2*Q15-Q1*Q14*Q15^2,
           Q2*Q11*Q13-Q4*Q6*Q14,
           -Q8^2*Q9^2*Q12-2*Q8^2*Q9*Q11*Q13-2*Q7*Q8*Q11^2* Q13-Q7^2*Q11^2*Q14+Q6*Q7*Q9*Q12*Q14+Q8^2* Q9^2*Q15+3*Q7*Q8*Q9*Q11*Q15+Q7^2*Q11^2*Q15-Q6* Q7*Q9*Q12*Q15+Q6*Q8*Q9*Q13*Q15+Q6*Q7* Q11*Q13*Q15-Q6*Q7*Q9*Q14*Q15,
           Q4*Q5*Q6*Q9*Q11-Q4^2*Q6*Q11^2+Q3*Q4* Q6*Q9*Q12-Q1*Q4*Q11^2*Q12-Q3*Q5*Q6*Q9*Q13-Q3^2*Q8*Q9*Q13+2*Q3*Q4*Q6*Q11*Q13-Q3^2* Q6*Q13^2-Q1*Q4*Q11^2*Q14+Q1*Q3*Q11*Q13*Q14+ Q1*Q3*Q11*Q13*Q15,
           Q4*Q8*Q9*Q12+Q5*Q7*Q11*Q13-Q4*Q8*Q11* Q13-Q2*Q11*Q13^2+Q4*Q7*Q11*Q14-Q4*Q6*Q13* Q14-Q1*Q12*Q13*Q14-Q5*Q7*Q9*Q15-Q4*Q6* Q13*Q15+Q3*Q7*Q13*Q15+Q1*Q13*Q14*Q15+Q1* Q13*Q15^2,
           -Q5^2*Q7*Q11-Q4*Q5*Q8*Q11+Q3*Q5*Q7*Q12-Q2*Q5*Q9*Q12+Q2*Q4*Q11*Q12- Q1*Q5*Q12^2+Q2*Q5*Q11*Q13+Q3*Q4*Q8*Q14+Q3*Q5*Q7*Q15-Q1*Q5*Q12*Q15-Q2*Q3*Q13* Q15+Q1*Q5*Q14*Q15,
           Q5*Q7^2*Q8*Q9-Q5*Q7^3*Q11+Q5*Q6*Q7^2* Q13-Q2*Q7^2*Q11*Q13+Q2*Q6*Q7*Q13^2-Q1*Q8^2* Q13^2-Q1*Q7*Q8*Q13*Q14+Q1*Q7^2*Q14*Q15,
           Q5*Q6*Q7*Q9-Q4*Q6*Q8*Q9+Q4*Q6^2*Q13-Q3*Q6*Q7*Q13-2*Q1*Q8*Q11*Q13+Q1*Q6* Q12*Q13-Q1*Q7*Q11*Q14+Q1*Q6*Q13*Q14+Q1*Q7*Q11*Q15,
           -Q3*Q5^2*Q6*Q7+Q3^2*Q5*Q7*Q8- Q2*Q3*Q5*Q7*Q11+2*Q1*Q5^2*Q8*Q11+Q1*Q5^2* Q6*Q12-Q1*Q3*Q5*Q8*Q12-Q1*Q5^2*Q6*Q14-Q1*Q2*Q5*Q11*Q14+Q2*Q3^2*Q7*Q15-Q1* Q2*Q3*Q12*Q15+Q1*Q2*Q3*Q14*Q15,
           Q3*Q13*Q14-Q5*Q9*Q15,
           Q4*Q5*Q8*Q11-2*Q3*Q5*Q7*Q12+2*Q2*Q5* Q9*Q12-Q2*Q4*Q11*Q12+2*Q1*Q5*Q12^2-Q5^2*Q6*Q13-Q1*Q5*Q12*Q14+Q3*Q5*Q7*Q15-Q2*Q3* Q13*Q15,
           -Q3*Q4*Q5*Q8*Q9-Q2*Q4^2*Q11^2+Q1* Q5^2*Q9*Q12+2*Q1*Q4*Q5*Q11*Q12+Q2*Q3* Q4*Q11*Q13-Q1*Q5^2*Q11*Q13-Q1*Q3*Q4*Q12*Q14-Q1*Q5^2*Q9*Q15-Q1*Q3*Q4*Q12*Q15+Q1* Q3*Q5*Q13*Q15+Q1*Q3*Q4*Q14*Q15,
           Q4*Q8^2*Q12*Q13-Q2*Q8*Q12*Q13^2-Q4*Q8^2* Q13*Q14+Q2*Q7*Q12*Q13*Q14+Q2*Q8*Q13^2*Q14-Q5*Q7^2*Q14^2+Q4*Q8^2*Q13*Q15-Q2*Q7*Q12*Q13* Q15-Q2*Q8*Q13^2*Q15+Q5*Q7^2*Q14*Q15-Q4*Q7*Q8*Q14*Q15+Q2*Q7*Q13*Q14*Q15,
           Q4*Q6*Q12-Q5*Q6*Q13-Q3*Q8*Q13+Q2*Q11* Q13+Q1*Q12*Q14-Q1*Q12*Q15,
           Q5*Q8*Q9+Q4*Q8*Q11-Q5*Q6*Q13-Q2*Q11* Q13-Q3*Q7*Q14+Q2*Q9*Q14-Q1*Q14^2+Q2*Q9* Q15+Q1*Q12*Q15-Q1*Q14*Q15,
           Q5*Q8*Q9-Q4*Q8*Q11+Q5*Q6*Q13-Q2*Q11* Q13-Q3*Q7*Q14+Q2*Q9*Q14-Q1*Q14^2+Q2*Q9* Q15+Q1*Q12*Q15-Q1*Q14*Q15,
           -Q5*Q8*Q9-Q4*Q8*Q11+Q5*Q6*Q13+Q2*Q11* Q13-Q3*Q7*Q14+Q2*Q9*Q14-Q1*Q14^2+Q2*Q9*Q15+Q1*Q12*Q15-Q1*Q14*Q15,
           -Q5*Q8*Q9+Q4*Q8*Q11-Q5*Q6*Q13+Q2*Q11* Q13-Q3*Q7*Q14+Q2*Q9*Q14-Q1*Q14^2+Q2*Q9*Q15+Q1*Q12*Q15-Q1*Q14*Q15,
           Q5^2*Q8*Q9*Q11-Q4*Q5*Q8*Q11^2+Q3*Q5* Q8*Q9*Q12-Q2*Q5*Q9*Q11*Q12+Q2*Q4*Q11^2*Q12-Q3^2*Q8*Q13*Q14+Q3*Q5*Q8*Q9*Q15-Q3* Q4*Q8*Q11*Q15-Q2*Q5*Q9*Q11*Q15+Q2*Q4*Q11^2*Q15-Q2*Q3*Q9*Q12*Q15+Q3^2*Q8*Q13*Q15,
           Q5^2*Q7*Q9^2+Q4*Q5*Q8*Q9^2-Q4^2*Q7*Q11^2- Q3*Q4*Q7*Q9*Q12+Q1*Q4*Q9*Q12^2-Q3*Q5* Q7*Q9*Q13-Q2*Q5*Q9^2*Q13+Q1*Q5*Q9*Q12* Q13+Q3^2*Q7*Q13^2-Q1*Q3*Q12*Q13^2-Q1*Q4* Q9*Q12*Q14-Q1*Q5*Q9*Q13*Q14+2*Q1*Q4*Q11*Q13*Q14+Q1*Q3*Q13^2*Q14+Q3*Q4*Q7*Q9*Q15- Q1*Q5*Q9*Q13*Q15-Q1*Q3*Q13^2*Q15,
           Q4*Q5^2*Q6^2*Q8-2*Q3*Q4*Q5*Q6*Q8^2+Q3^2* Q4*Q8^3-Q2*Q5^2*Q6*Q8*Q9+Q2*Q3*Q4*Q8^ 2*Q11-Q2*Q3*Q4*Q6*Q8*Q12+Q1*Q5^2*Q6* Q8*Q12+Q2^2*Q5*Q6*Q9*Q12-Q2^2*Q4*Q6*Q11*Q12+Q1*Q2*Q5*Q8*Q11*Q12-Q2*Q4*Q5*Q6^2* Q14+2*Q2*Q3*Q4*Q6*Q8*Q14+Q1*Q5^2*Q6*Q8* Q14-Q1*Q3*Q5*Q8^2*Q14+Q2^2*Q5*Q6*Q9*Q14-Q1*Q2*Q5*Q8*Q11*Q14-Q1*Q2*Q5*Q6*Q12*Q14-Q1*Q2*Q5*Q6*Q14^2-Q1*Q3*Q5*Q8^2* Q15+Q1*Q2*Q3*Q8*Q14*Q15,
           -Q5^2*Q7*Q11+Q4*Q5*Q6*Q12-Q2*Q4*Q11*Q12+Q3*Q5*Q8*Q13-2*Q4*Q5*Q6*Q15+2*Q2*Q5* Q9*Q15-Q2*Q3*Q13*Q15-Q1*Q5*Q14*Q15+2*Q1*Q5*Q15^2,
           Q2*Q4*Q9-Q1*Q5*Q13,
           -Q4*Q8*Q12^2+Q5*Q8*Q12*Q13-Q2*Q12^2*Q13+Q5*Q7*Q12*Q14+Q4*Q8*Q12*Q14-Q5*Q8*Q13* Q14+Q5*Q7*Q12*Q15-Q5*Q8*Q13*Q15+Q2*Q12*Q13*Q15-Q5*Q7*Q14*Q15,
           -Q4*Q8*Q14+Q5*Q7*Q15,
           Q6*Q8*Q9*Q12*Q13-Q6^2*Q12*Q13^2-Q8^2*Q9^2* Q14+Q7^2*Q11^2*Q14-Q6*Q7*Q9*Q12*Q14+Q6*Q7* Q11*Q13*Q14+Q6^2*Q13^2*Q14-Q6*Q7*Q9*Q14^2+Q8^2*Q9^2*Q15-Q7^2*Q11^2*Q15+Q6*Q7*Q9*Q12*Q15- Q6*Q8*Q9*Q13*Q15,
           -Q5*Q8+Q2*Q14,
           Q5*Q8-Q2*Q14,
           Q3*Q7*Q8*Q9*Q12-Q2*Q8*Q9^2*Q12-Q1*Q8*Q9*Q12^2-Q3*Q8^2*Q9*Q13+Q2*Q8*Q9*Q11* Q13+2*Q2*Q6*Q9*Q12*Q13-Q1*Q8*Q11*Q12*Q13-Q5*Q6^2*Q13^2+Q1*Q6*Q12*Q13*Q14+Q3*Q7*Q8* Q9*Q15-Q2*Q7*Q9*Q11*Q15-Q1*Q8*Q9*Q12*Q15+Q1*Q8*Q11*Q13*Q15,
           Q4*Q8*Q9*Q12+Q5*Q8*Q9*Q13-Q4*Q8*Q11* Q13+Q3*Q7*Q12*Q13-Q2*Q9*Q12*Q13-Q1*Q12^2* Q13-Q3*Q8*Q13^2+Q4*Q7*Q11*Q14-Q5*Q7*Q9*Q15+Q3*Q7*Q13*Q15-Q1*Q12*Q13*Q15+Q1*Q13*Q14*Q15,
           -Q5*Q7*Q11+Q4*Q8*Q11-Q2*Q11* Q13+Q3*Q7*Q14-Q1*Q12*Q14+Q1*Q14*Q15,
           Q3*Q8-Q2*Q11,
           -Q3*Q8+Q2*Q11,
           Q5^3*Q7^2*Q9+Q4*Q5^2*Q7*Q8*Q9-Q2*Q4* Q5*Q7*Q9*Q12-Q2*Q4^2*Q8*Q9*Q12-2*Q2*Q5^2* Q7*Q9*Q13+Q2*Q4^2*Q6*Q12*Q13+Q1*Q4*Q5* Q8*Q12*Q13-Q2*Q4*Q5*Q6*Q13^2+Q2^2*Q5* Q9*Q13^2+Q1*Q2*Q5*Q12*Q13^2-Q1*Q5^2*Q7*Q13*Q14+2*Q2*Q4*Q5*Q7*Q9*Q15+Q2*Q4^2*Q6*Q13* Q15-Q1*Q5^2*Q7*Q13*Q15-Q1*Q4*Q5*Q8*Q13* Q15-Q2^2*Q4*Q9*Q13*Q15-Q1*Q2*Q4*Q12*Q13*Q15+Q1*Q2*Q5*Q13^2*Q15+Q1*Q4*Q5*Q7* Q14*Q15-Q1*Q2*Q4*Q13*Q15^2,
           Q2*Q6*Q13-Q1*Q8*Q14,
           -Q5*Q6*Q7*Q8*Q9+Q1*Q8^2*Q9*Q12+Q2*Q6*Q7*Q11*Q13-Q1*Q8^2*Q11*Q13+2*Q1*Q6*Q8* Q12*Q13-Q2*Q6^2*Q13^2-Q1*Q6*Q7*Q12*Q14-Q1*Q8^2*Q9*Q15+Q1*Q7*Q8*Q11*Q15-Q1*Q6*Q7* Q12*Q15+Q1*Q6*Q7*Q14*Q15,
           -Q5^2*Q6*Q11*Q12+Q2*Q5*Q11^2*Q12-Q5^2*Q6*Q11*Q14+Q2*Q5*Q11^2*Q14-Q3*Q5*Q6*Q12*Q14+ Q5^2*Q6*Q11*Q15-Q2*Q5*Q11^2*Q15+Q3*Q5*Q6* Q12*Q15+Q3*Q5*Q6*Q14*Q15+Q3^2*Q8*Q14*Q15- Q2*Q3*Q11*Q14*Q15-Q3^2*Q8*Q15^2,
           Q8^2*Q9^2*Q12-Q7^2*Q11^2*Q12-Q6*Q7*Q9*Q12^2+ Q6*Q8*Q9*Q12*Q13+Q6^2*Q12*Q13^2-Q6*Q7*Q9* Q12*Q14+Q6*Q7*Q11*Q13*Q14-Q6^2*Q13^2*Q14- Q8^2*Q9^2*Q15+Q7^2*Q11^2*Q15-Q6*Q7*Q11*Q13*Q15+ Q6*Q7*Q9*Q14*Q15,
           -Q5^2*Q6^2*Q7+Q3*Q5*Q6*Q7*Q8-Q2*Q5*Q6*Q7*Q11-2*Q1*Q5*Q8^2*Q11-Q1*Q5*Q6*Q8* Q12+3*Q1*Q5*Q6*Q8*Q14+Q1*Q3*Q8^2*Q14+Q1* Q2*Q8*Q11*Q14-Q1*Q2*Q6*Q12*Q14+Q2*Q3* Q6*Q7*Q15-Q1*Q3*Q8^2*Q15,
           -Q5*Q7*Q11+Q4*Q6*Q12+Q2*Q11*Q13-Q1*Q12* Q15,
           2*Q5^2*Q8^2*Q11+Q3*Q5*Q8^2*Q14-Q2*Q5*Q8*Q11*Q14-Q2*Q5*Q6*Q12*Q14-Q5^2*Q6*Q8* Q15-Q3*Q5*Q8^2*Q15-2*Q2*Q5*Q8*Q11*Q15+Q2* Q5*Q6*Q12*Q15+Q2*Q5*Q6*Q14*Q15+Q2^2*Q11*Q15^2,
           Q2*Q6*Q9-Q1*Q8*Q11,
           Q4*Q8*Q12-Q5*Q8*Q13+Q2*Q12*Q13-Q4*Q8* Q14,
           -Q4*Q5*Q6*Q11*Q12+Q2*Q4*Q11^2*Q12+Q5^2*Q6*Q11*Q13-Q2*Q5*Q11^2*Q13+Q3*Q5*Q6* Q12*Q13-Q4*Q5*Q6*Q11*Q14+Q2*Q4*Q11^2*Q14-Q3*Q4*Q6*Q12*Q14+Q3*Q5*Q6*Q13*Q14+Q3^2* Q8*Q13*Q14-Q2*Q3*Q11*Q13*Q14-Q3^2*Q8*Q13*Q15,
           Q5*Q6^2*Q8*Q9+Q3*Q6*Q8^2*Q9+Q5*Q6^ 2*Q7*Q11-Q1*Q8^2*Q11^2-Q1*Q6*Q8*Q11*Q12- Q5*Q6^3*Q13-Q3*Q6^2*Q8*Q13+Q1*Q6^2*Q12*Q14,
           -Q2*Q5^2*Q7*Q9-Q2*Q4*Q5*Q8*Q9+Q2^2*Q4*Q9*Q12+2*Q1*Q5^2*Q8*Q13+Q2^2*Q5*Q9*Q13+ Q1*Q5^2*Q7*Q14-Q1*Q2*Q4*Q12*Q14-Q1*Q2* Q5*Q13*Q14-Q1*Q5^2*Q7*Q15-Q1*Q4*Q5*Q8* Q15+Q1*Q2*Q4*Q12*Q15,
           Q4*Q8*Q9*Q12-Q5*Q8*Q9*Q13+Q5*Q7*Q11* Q13-Q5*Q6*Q13^2-Q4*Q7*Q11*Q14+Q3*Q7*Q13* Q14-Q2*Q9*Q13*Q14+Q1*Q13*Q14^2+Q5*Q7*Q9*Q15-Q2*Q9*Q13*Q15-Q1*Q12*Q13*Q15+Q1*Q13*Q14*Q15,
           -Q4^2*Q6*Q8*Q12+Q3*Q4*Q7*Q8* Q12-Q1*Q4*Q8*Q12^2+Q4*Q5*Q6*Q8*Q13-Q3*Q4*Q8^2*Q13+2*Q2*Q4*Q6*Q12*Q13-Q1*Q5*Q8* Q12*Q13-Q2^2*Q11*Q13^2-Q4*Q5*Q6*Q7*Q14+Q3* Q4*Q7*Q8*Q14-Q1*Q4*Q8*Q12*Q14+Q1*Q5* Q8*Q13*Q14+Q1*Q2*Q12*Q13*Q15,
           Q2*Q3*Q7-Q1*Q5*Q8,
           Q5*Q6*Q8*Q9^2+Q3*Q8^2*Q9^2-Q3*Q7*Q8* Q9*Q11-Q2*Q8*Q9^2*Q11+Q3*Q7^2*Q11^2-Q3*Q6* Q7*Q9*Q12+Q1*Q8*Q9*Q11*Q12-Q1*Q7*Q11^2* Q12+Q1*Q6*Q9*Q12^2-Q3*Q6^2*Q13^2-Q1*Q8* Q9*Q11*Q14+Q1*Q7*Q11^2*Q14-Q1*Q6*Q9*Q12*Q14+2*Q1*Q6*Q11*Q13*Q14+Q3*Q6*Q7*Q9*Q15- Q1*Q8*Q9*Q11*Q15-Q1*Q7*Q11^2*Q15,
           -Q2*Q5*Q7*Q8*Q9-Q2*Q4*Q8^2*Q9-Q1*Q5*Q7*Q8*Q12-Q1*Q4*Q8^2*Q12+2*Q1*Q5*Q8^2* Q13+Q2^2*Q8*Q9*Q13+Q1*Q4*Q8^2*Q14-Q1*Q2* Q8*Q13*Q14+Q2^2*Q7*Q9*Q15+Q1*Q2*Q7*Q12* Q15-Q1*Q2*Q7*Q14*Q15,
           -Q8^2*Q9^2*Q12+Q6*Q8*Q9*Q12*Q13-Q7^2*Q11^2*Q14+Q6*Q7*Q9*Q12*Q14-Q6*Q8*Q9*Q13*Q14+ Q6^2*Q13^2*Q14+Q8^2*Q9^2*Q15+Q7*Q8*Q9*Q11*Q15+Q7^2*Q11^2*Q15-Q6*Q7*Q9*Q12*Q15-Q6^2*Q13^2* Q15-Q6*Q7*Q9*Q15^2,
           2*Q5^2*Q9*Q11*Q12*Q13+2*Q4*Q5*Q11^2*Q12*Q13-2* Q5^2*Q11^2*Q13^2-Q5^2*Q9^2*Q12*Q14-3*Q4*Q5*Q9*Q11*Q12*Q14+Q5^2*Q9*Q11*Q13*Q14-Q3*Q5*Q9* Q12*Q13*Q14-Q3*Q4*Q11*Q12*Q13*Q14+Q3*Q5*Q11*Q13^2*Q14+Q4^2*Q11^2*Q14^2+Q3*Q4*Q9*Q12* Q14^2+Q3*Q4*Q11*Q13*Q14^2+Q5^2*Q9^2*Q12*Q15-2*Q5^ 2*Q9*Q11*Q13*Q15-Q4*Q5*Q11^2*Q13*Q15+Q3* Q5*Q11*Q13^2*Q15+Q5^2*Q9^2*Q14*Q15-Q3^2*Q13^2*Q14* Q15-Q3*Q4*Q9*Q14^2*Q15+Q3*Q5*Q9*Q13*Q15^2,
           -Q5*Q6*Q13-Q3*Q8*Q13+Q2*Q11*Q13+Q1* Q14^2,
           -Q5*Q7*Q8*Q9-Q4*Q8^2*Q9-Q4*Q7*Q8* Q11+Q5*Q6*Q7*Q13-Q3*Q7*Q8*Q13+2*Q2*Q8* Q9*Q13+2*Q1*Q8*Q12*Q13-Q2*Q6*Q13^2-Q3* Q7^2*Q15+Q2*Q7*Q9*Q15+Q1*Q7*Q12*Q15-Q1*Q7*Q14*Q15,
           Q3*Q4*Q7-Q1*Q4*Q12-Q1*Q5*Q13+Q1*Q4*Q14,
           Q2*Q7*Q11-Q1*Q8*Q12,
           -Q9*Q12+Q11*Q13,
           Q3*Q8^2*Q9+Q5*Q6*Q7*Q11-2*Q2*Q8*Q9* Q11+Q2*Q7*Q11^2-2*Q1*Q8*Q11*Q12-Q5*Q6^2*Q13-Q3*Q6*Q7*Q14+2*Q1*Q8*Q11*Q14+Q2*Q6*Q9* Q15+Q1*Q6*Q12*Q15-Q1*Q6*Q14*Q15,
           -Q1*Q12^2+Q5*Q6*Q13,
           -Q5*Q8^2*Q9*Q11+Q5*Q7*Q8*Q11^2-Q5*Q6^2*Q12*Q13-Q5*Q6*Q8*Q9*Q14+Q3*Q8^2*Q9*Q14+ Q5*Q6*Q7*Q11*Q14-Q3*Q7*Q8*Q11*Q14+Q5*Q6^2*Q13*Q14+Q3*Q8^2*Q9*Q15+Q5*Q6*Q7*Q11* Q15-Q3*Q7*Q8*Q11*Q15-Q3*Q6*Q7*Q14*Q15,
           -2*Q5^2*Q8^2*Q11^2+2*Q5^2*Q6*Q8*Q11*Q12+2*Q2*Q5*Q8*Q11^2*Q12-2*Q5^2*Q6*Q8*Q11*Q14+Q3*Q5* Q8^2*Q11*Q14-Q2*Q5*Q8*Q11^2*Q14+Q5^2*Q6^2*Q12*Q14+Q3*Q5*Q6*Q8*Q14^2+Q5^2*Q6*Q8*Q11* Q15+Q3*Q5*Q8^2*Q11*Q15-Q5^2*Q6^2*Q12*Q15-Q3* Q5*Q6*Q8*Q12*Q15-3*Q2*Q5*Q6*Q11*Q12*Q15- Q2*Q3*Q8*Q11*Q12*Q15+Q5^2*Q6^2*Q14*Q15- Q3^2*Q8^2*Q14*Q15+Q2*Q3*Q8*Q11*Q15^2+Q2^2*Q11^2* Q15^2+Q2*Q3*Q6*Q12*Q15^2-Q2*Q3*Q6*Q14* Q15^2,
           -Q5^2*Q6*Q7^2+Q4*Q5*Q6*Q7*Q8+Q2*Q4* Q6*Q7*Q12-Q1*Q4*Q8^2*Q12-Q2*Q5*Q6*Q7* Q13-2*Q1*Q5*Q8^2*Q13+3*Q1*Q5*Q7*Q8*Q14+ Q1*Q4*Q8^2*Q14+Q1*Q2*Q8*Q13*Q14-Q1*Q5*Q7*Q8*Q15-Q1*Q2*Q7*Q14*Q15,
           -Q5*Q8*Q9-Q4*Q8*Q11+Q3*Q8*Q13+Q2*Q9* Q15+Q1*Q12*Q15-Q1*Q14*Q15,
           Q4*Q8*Q11*Q12+Q2*Q9*Q12^2+Q5*Q8*Q11*Q13-Q3*Q8*Q12*Q13+Q5*Q7*Q11*Q14+Q4*Q8*Q11*Q14-Q4*Q6*Q12*Q14+Q2*Q9*Q12*Q14-Q1*Q12^2*Q14-Q3*Q8*Q13*Q14-Q2*Q11*Q13*Q14-Q1* Q12*Q14^2-Q5*Q7*Q11*Q15+Q1*Q12^2*Q15-Q2*Q11*Q13*Q15+Q1*Q12*Q14*Q15,
           -Q4*Q8^2*Q9^2+Q5*Q7^2*Q9*Q11+Q4*Q7^2*Q11^2- Q4*Q6*Q7*Q11*Q13-Q3*Q7^2*Q11*Q13+2*Q1*Q8*Q9*Q12*Q13-Q1*Q7*Q11*Q12*Q13+Q4*Q6^2* Q13^2+Q1*Q6*Q12*Q13^2+Q4*Q6*Q7*Q9*Q14-Q1*Q7*Q11*Q13*Q14-Q1*Q6*Q13^2*Q14-Q4*Q6*Q7* Q9*Q15-Q1*Q7*Q9*Q12*Q15+Q1*Q7*Q11*Q13*Q15-Q1*Q6*Q13^2*Q15+Q1*Q7*Q9*Q15^2,
           -Q4*Q6^2+Q3*Q6*Q7-Q1*Q8*Q11+Q1*Q6*Q15,
           Q5*Q9*Q12-Q3*Q13*Q15,
           Q5^2*Q7*Q8*Q11-Q4*Q5*Q8^2*Q11-Q2^2*Q11* Q12*Q13-Q3*Q5*Q7*Q8*Q14+Q3*Q4*Q8^2*Q14+Q2*Q5*Q7*Q11*Q14-Q3*Q5*Q7*Q8*Q15+Q3* Q4*Q8^2*Q15+Q2*Q5*Q7*Q11*Q15-Q2*Q4*Q8*Q11*Q15+Q2^2*Q11*Q13*Q15-Q2*Q3*Q7*Q14*Q15,
           Q3*Q8^3*Q9^2-2*Q3*Q7*Q8^2*Q9*Q11-Q4*Q6* Q7*Q8*Q11^2+Q3*Q7^2*Q8*Q11^2-Q1*Q8^2*Q9*Q11*Q12+Q3*Q6*Q8^2*Q9*Q13-Q3*Q6*Q7*Q8* Q9*Q14+Q4*Q6^2*Q7*Q11*Q14+Q1*Q7*Q8*Q11^2*Q14-Q3*Q6^2*Q7*Q13*Q14+Q1*Q6*Q8*Q11*Q13* Q14+2*Q3*Q6*Q7*Q8*Q9*Q15+Q4*Q6^2*Q7*Q11* Q15-Q3*Q6*Q7^2*Q11*Q15-Q1*Q8^2*Q9*Q11*Q15+Q1*Q7*Q8*Q11^2*Q15+Q1*Q6*Q8*Q9*Q12* Q15-Q1*Q6*Q8*Q11*Q13*Q15-Q1*Q6*Q7*Q11*Q14*Q15-Q1*Q6*Q7*Q11*Q15^2,
           Q11*Q13*Q14-Q9*Q15^2,
           2*Q8*Q11^2*Q13^2-Q8*Q9*Q11*Q13*Q14+Q7*Q11^2* Q13*Q14-Q6*Q9*Q12*Q13*Q14-2*Q8*Q9*Q11*Q13*Q15-Q7*Q11^2*Q13*Q15+Q6*Q9*Q12*Q13*Q15-Q6* Q11*Q13^2*Q15+Q6*Q9*Q13*Q14*Q15+Q8*Q9^2*Q15^ 2,
           -Q5*Q8*Q9+Q3*Q8*Q13+Q4*Q6*Q14-Q1* Q14*Q15,
           -Q5^2*Q8*Q9^2+Q4^2*Q8*Q11^2-Q3*Q4*Q8* Q9*Q12+Q3*Q4*Q8*Q11*Q13+Q1*Q5*Q11*Q12* Q13+Q3^2*Q8*Q13^2-Q3*Q4*Q8*Q9*Q14-Q1*Q4*Q11*Q12*Q14-Q3*Q4*Q8*Q9*Q15+Q1*Q5*Q9*Q12*Q15+Q1*Q5*Q11*Q13*Q15-Q1*Q3*Q13*Q14*Q15,
           Q4*Q11-Q3*Q13,
           -Q4*Q11+Q3*Q13,
           Q5^2*Q6*Q7*Q11-Q3*Q5*Q6*Q7*Q12+Q1*Q5*Q8*Q11*Q12+Q1*Q5*Q6*Q12^2-Q1*Q5*Q8* Q11*Q14-Q3^2*Q7*Q8*Q15+Q2*Q3*Q7*Q11*Q15-Q1*Q5*Q8*Q11*Q15+Q1*Q3*Q8*Q12*Q15-Q1* Q2*Q11*Q12*Q15-Q1*Q3*Q8*Q14*Q15+Q1*Q2*Q11*Q14*Q15,
           -Q5*Q6+Q2*Q11,
           Q5*Q6-Q2*Q11,
           -Q5^2*Q7*Q11+Q4*Q5*Q8*Q11-Q4*Q5*Q6*Q12+Q3*Q5*Q7*Q12+Q2*Q4*Q11*Q12-Q1*Q5*Q12^2-Q2*Q5*Q11*Q13+Q3*Q5*Q7*Q14-Q3*Q4* Q8*Q14-Q1*Q5*Q12*Q14+Q2*Q3*Q13*Q15+Q1*Q5*Q14*Q15,
           -Q5*Q11*Q12*Q13-Q4*Q11*Q12*Q14- Q5*Q11*Q13*Q14+Q5*Q9*Q12*Q15+Q4*Q11*Q12*Q15+Q5*Q11*Q13*Q15+Q4*Q11*Q14*Q15+Q3*Q13*Q14*Q15-Q5*Q9*Q15^2-Q3*Q13*Q15^2,
           -Q3*Q4*Q8+Q1*Q5*Q14,
           Q3*Q4*Q8-Q1*Q5*Q14,
           Q4*Q8*Q9*Q12+Q5*Q7*Q11*Q13-Q4*Q8*Q11* Q13-Q2*Q11*Q13^2+Q4*Q7*Q11*Q14-Q3*Q7*Q13* Q14+Q2*Q9*Q13*Q14-Q1*Q13*Q14^2-Q5*Q7*Q9*Q15+Q2*Q9*Q13*Q15+Q1*Q12*Q13*Q15-Q1*Q13*Q14*Q15,
           -Q5*Q7*Q11-Q4*Q8*Q11+Q2*Q11* Q13+Q3*Q7*Q15-Q1*Q12*Q15+Q1*Q14*Q15,
           -Q5^2*Q6*Q9^2+Q4*Q5*Q6*Q9*Q11+3*Q1*Q5*Q9*Q11*Q12+Q1*Q4*Q11^2*Q12-Q3*Q5*Q6*Q9* Q13-2*Q1*Q5*Q11^2*Q13+Q1*Q3*Q11*Q12*Q13+Q3* Q4*Q6*Q9*Q14-Q1*Q4*Q11^2*Q14-Q1*Q5*Q9* Q11*Q15-Q1*Q3*Q9*Q12*Q15,
           -Q4*Q5*Q7*Q8*Q12+Q4^2*Q8^2*Q12-Q2^2*Q12*Q13^2-Q5^2*Q7^2*Q14+Q4*Q5*Q7*Q8*Q14+Q2*Q4* Q7*Q12*Q14+Q5^2*Q7^2*Q15-Q4^2*Q8^2*Q15+Q2*Q5*Q7*Q13*Q15+Q2^2*Q13^2*Q15-Q2*Q4*Q7*Q14* Q15-Q2*Q4*Q7*Q15^2,
           Q4*Q9*Q12^2+Q5*Q9*Q12*Q13+Q4*Q11*Q12*Q13-2*Q5*Q11*Q13^2+Q3*Q12*Q13^2-Q4*Q11*Q13*Q14- Q5*Q9*Q13*Q15,
           -Q4*Q6^2+Q2*Q6*Q9-Q1*Q8*Q11+Q1*Q6*Q15,
           Q2*Q4*Q11*Q12-Q5^2*Q6*Q13+Q3*Q5*Q8*Q13-Q2*Q5*Q11*Q13+Q4*Q5*Q6*Q14-Q3*Q4* Q8*Q14+Q1*Q5*Q12*Q14+Q4*Q5*Q6*Q15-Q3*Q5*Q7*Q15+Q2*Q3*Q13*Q15-Q1*Q5*Q14*Q15-Q1*Q5*Q15^2,
           Q5^2*Q7+Q4*Q5*Q8-Q2*Q4*Q12-Q2*Q5*Q13,
           -Q5^2*Q7-Q4*Q5*Q8+Q2*Q4*Q12+ Q2*Q5*Q13,
           -Q5^2*Q7+Q4*Q5*Q8+Q2*Q4*Q12- Q2*Q5*Q13,
           Q5^2*Q7-Q4*Q5*Q8-Q2*Q4*Q12+Q2*Q5*Q13,
           -Q5^2*Q7+Q4*Q5*Q8-Q2*Q4* Q12+Q2*Q5*Q13,
           Q5*Q8*Q11*Q12-Q5*Q6*Q12^2- Q2*Q11*Q12^2-Q5*Q8*Q11*Q14+Q5*Q6*Q12*Q14+Q3*Q8*Q12*Q14-Q5*Q8*Q11*Q15+Q3*Q8*Q12* Q15+Q2*Q11*Q12*Q15-Q3*Q8*Q14*Q15,
           Q5^2*Q6-Q2*Q3*Q15,
           Q5^2*Q7*Q9*Q11-Q3*Q5*Q7*Q9*Q12+Q1*Q5*Q9*Q12^2+Q1*Q5*Q11*Q12*Q13+Q3*Q4*Q7* Q11*Q14-Q1*Q4*Q11*Q12*Q14-Q3^2*Q7*Q13*Q14-Q1*Q5*Q11*Q13*Q14+Q1*Q3*Q12*Q13*Q14-Q1* Q5*Q11*Q13*Q15+Q1*Q4*Q11*Q14*Q15-Q1*Q3*Q13*Q14*Q15,
           Q2*Q12*Q13-Q5*Q7*Q15,
           Q3^2*Q5*Q8*Q9+Q3^2*Q4*Q8*Q11+Q2*Q3* Q4*Q11^2-Q1*Q5^2*Q11^2-Q3^3*Q8*Q13-Q2*Q3^2*Q11* Q13-Q1*Q3*Q5*Q11*Q14+Q1*Q3^2*Q14*Q15,
           -Q4^2*Q11*Q12^2-Q5^2*Q9*Q12*Q13+Q3*Q5*Q12*Q13^2-Q4*Q5*Q9*Q12*Q14+Q4^2*Q11*Q12*Q14+Q5^2* Q9*Q13*Q14+Q3*Q4*Q12*Q13*Q14-Q3*Q5*Q13^2* Q14+Q5^2*Q9*Q13*Q15+Q3*Q4*Q12*Q13*Q15-Q3*Q5*Q13^2*Q15-Q3*Q4*Q13*Q14*Q15,
           -Q4*Q8*Q12+Q2*Q13*Q15,
           -Q5*Q8*Q9+Q5*Q7*Q11-Q5*Q6*Q13+Q4*Q6* Q14+Q1*Q12*Q14-Q1*Q14*Q15,
           Q3*Q5*Q6*Q7*Q8-Q3^2*Q7*Q8^2-Q2*Q5* Q6*Q7*Q11+2*Q2*Q3*Q7*Q8*Q11-Q2^2*Q7*Q11^2+Q1*Q2*Q8*Q11*Q12-Q2^2*Q6*Q11*Q13+Q2*Q3* Q6*Q7*Q14-Q1*Q3*Q8^2*Q14-Q1*Q3*Q8^2*Q15+ Q1*Q2*Q8*Q11*Q15,
           -Q5*Q8*Q11-Q3*Q8*Q14+Q3*Q8*Q15+Q2*Q11* Q15,
           -Q5*Q8*Q11-Q2*Q11*Q12+Q3*Q8*Q15+Q2*Q11*Q15,
           -Q5*Q8*Q12*Q13-Q4*Q8*Q12*Q14- Q5*Q8*Q13*Q14+Q4*Q8*Q12*Q15+Q5*Q8*Q13*Q15+Q2*Q12*Q13*Q15+Q5*Q7*Q14*Q15+Q4*Q8*Q14*Q15-Q5*Q7*Q15^2-Q2*Q13*Q15^2,
           -Q4*Q6*Q12+Q3*Q8*Q13,
           Q5*Q8*Q9-Q5*Q6*Q13+Q3*Q7*Q15-Q1*Q12* Q15,
           -Q3*Q5^2*Q7*Q9-Q3*Q4*Q5*Q7*Q11+Q1*Q5^2*Q9*Q12+Q3^2*Q5*Q7*Q13+2*Q1*Q5^2*Q11* Q13-Q1*Q3*Q5*Q12*Q13+Q3^2*Q4*Q7*Q14-Q1* Q3*Q4*Q12*Q14-Q1*Q5^2*Q9*Q15-Q1*Q4*Q5* Q11*Q15+Q1*Q3*Q4*Q14*Q15,
           Q5*Q7*Q8*Q11+Q2*Q8*Q9*Q12-Q3*Q8^2*Q13-Q2*Q6*Q12*Q13-Q5*Q6*Q7*Q14+2*Q4*Q6* Q8*Q14-2*Q2*Q8*Q9*Q14+2*Q1*Q8*Q14^2-Q1*Q8*Q14*Q15,
           Q8*Q12*Q13-Q7*Q14^2,
           2*Q8^2*Q9*Q11*Q12*Q13-2*Q8^2*Q11^2*Q13^2+2*Q6*Q8*Q11*Q12*Q13^2-Q8^2*Q9*Q11*Q13*Q14+Q7*Q8* Q11^2*Q13*Q14-2*Q6*Q8*Q11*Q13^2*Q14+Q6^2*Q12*Q13^ 2*Q14+Q6*Q7*Q11*Q13*Q14^2-Q7*Q8*Q9*Q11* Q12*Q15+Q7*Q8*Q11^2*Q13*Q15-3*Q6*Q8*Q9*Q12*Q13*Q15-Q6*Q7*Q11*Q12*Q13*Q15+Q6*Q8*Q11* Q13^2*Q15-Q6^2*Q12*Q13^2*Q15-Q7^2*Q11^2*Q14*Q15+Q6^2*Q13^2*Q14*Q15+Q8^2*Q9^2*Q15^2+Q7*Q8*Q9*Q11* Q15^2+Q6*Q7*Q9*Q12*Q15^2-Q6*Q7*Q9*Q14*Q15^2,
           Q4*Q6*Q7-Q1*Q8*Q13+Q1*Q7*Q14-Q1*Q7*Q15,
           2*Q5^2*Q8*Q13^2-Q5^2*Q7*Q13*Q14-2*Q4*Q5*Q8*Q13*Q14+Q2*Q4*Q12*Q13*Q14-Q2*Q5*Q13^2*Q14+Q4^2*Q8*Q14^2+Q5^2*Q7*Q13*Q15-Q4*Q5* Q8*Q13*Q15-Q2*Q4*Q12*Q13*Q15+Q2*Q4*Q13*Q14*Q15,
           -Q5*Q7*Q8+Q4*Q8^2+Q2*Q8*Q13-Q2* Q7*Q15,
           -Q5*Q7*Q8-Q4*Q8^2+Q2*Q8*Q13+Q2*Q7*Q15,
           Q5*Q7*Q8+Q4*Q8^2-Q2*Q8*Q13-Q2*Q7*Q15,
           Q5*Q7*Q8-Q4*Q8^2+Q2*Q8*Q13-Q2*Q7*Q15,
           -Q2*Q7*Q11+Q1*Q8*Q15,
           Q2*Q7*Q11-Q1*Q8*Q15,
           Q5*Q7*Q8-Q4*Q8^2-Q2*Q8*Q13+Q2*Q7*Q15,
           Q3*Q7-Q1*Q12,
           -Q4*Q8^2*Q9+Q4*Q7*Q8*Q11-Q5*Q6*Q7*Q13-Q3*Q7*Q8*Q13+2*Q2*Q8*Q9*Q13-Q2*Q7* Q11*Q13+2*Q1*Q8*Q12*Q13-Q2*Q6*Q13^2-Q3*Q7^2*Q15+Q2*Q7*Q9*Q15+Q1*Q7*Q12*Q15-Q1*Q7* Q14*Q15,
           Q4*Q5*Q9*Q11*Q12-Q4^2*Q11^2*Q12-Q5^2 *Q9^2*Q14+Q4^2*Q11^2*Q14-Q3*Q4*Q9*Q12*Q14+Q3*Q4*Q11*Q13*Q14+Q3^2*Q13^2*Q14-Q3*Q4*Q9* Q14^2+Q5^2*Q9^2*Q15-Q4*Q5*Q9*Q11*Q15+Q3*Q4*Q9*Q12*Q15-Q3^2*Q13^2*Q15,
           Q8^2*Q9*Q12*Q13-Q7*Q8*Q11*Q12*Q13+2*Q8^2* Q11*Q13^2-Q6*Q7*Q12*Q13*Q14-Q8^2*Q9*Q13*Q15-2* Q7*Q8*Q11*Q13*Q15+Q6*Q7*Q12*Q13*Q15-Q6* Q8*Q13^2*Q15+Q6*Q7*Q13*Q14*Q15+Q7^2*Q11*Q15^2,
           Q5*Q7*Q8*Q11+Q4*Q6*Q8*Q12-Q3*Q8^2*Q13-Q2*Q6*Q12*Q13-2*Q4*Q6*Q8*Q15+2*Q2*Q8* Q9*Q15-Q2*Q7*Q11*Q15-Q1*Q8*Q14*Q15+2*Q1* Q8*Q15^2,
           Q5*Q7*Q8*Q9+Q4*Q8^2*Q9-Q5* Q7^2*Q11+Q4*Q6*Q7*Q12-2*Q4*Q6*Q8*Q13+Q2*Q6*Q13^2+Q1*Q7*Q12*Q14-2*Q1*Q8*Q13*Q14-Q2* Q7*Q9*Q15-Q1*Q7*Q12*Q15+2*Q1*Q8*Q13*Q15,
           -Q8*Q9*Q12+Q7*Q11*Q14,
           Q8*Q9*Q12-Q7*Q11*Q14,
           -Q5*Q8*Q9*Q12-Q4*Q8*Q11*Q12+Q5*Q8*Q11* Q13+Q3*Q8*Q13*Q14-Q2*Q11*Q13*Q14+Q4*Q6* Q14^2+Q1*Q12*Q14^2+Q5*Q8*Q9*Q15-Q4*Q8*Q11*Q15+Q3*Q8*Q13*Q15-Q2*Q11*Q13*Q15+Q4*Q6*Q14*Q15-Q3*Q7*Q14*Q15+Q1*Q12*Q14*Q15-Q1*Q14^2*Q15-Q1*Q14*Q15^2,
           Q5*Q6*Q13-Q3*Q8*Q13-Q2*Q11*Q13+Q3*Q7* Q14-Q1*Q12*Q14+Q1*Q14*Q15,
           -Q5*Q6*Q9*Q12+Q2*Q9*Q11*Q12+Q5*Q6*Q11* Q13-Q3*Q8*Q11*Q13-Q2*Q11^2*Q13-Q4*Q6*Q11*Q14+Q2*Q9*Q11*Q14-Q1*Q11*Q12*Q14+Q3*Q6*Q13*Q14-Q1*Q11*Q14^2+Q3*Q8*Q9*Q15+Q1* Q11*Q12*Q15,
           -Q5^2*Q6^2*Q9-Q3*Q5*Q6*Q8*Q9+Q2*Q5*Q6*Q9*Q11-2*Q1*Q5*Q8*Q11^2+3*Q1*Q5*Q6*Q11*Q12+Q1*Q3*Q8*Q11*Q12+Q1*Q2*Q11^2*Q12-Q1*Q5*Q6*Q11*Q14-Q1*Q3*Q6*Q12* Q14+Q2*Q3*Q6*Q9*Q15-Q1*Q2*Q11^2*Q15,
           -Q5^2*Q6*Q7+2*Q3*Q5*Q7*Q8-Q3*Q4*Q8^2-Q2*Q5*Q8*Q9-Q2*Q5*Q7*Q11+Q2*Q4*Q8* Q11-Q2*Q5*Q6*Q13+2*Q1*Q5*Q8*Q14+Q2*Q3*Q7*Q15-Q2^2*Q9*Q15-Q1*Q2*Q12*Q15+Q1*Q2* Q14*Q15,
           -Q8*Q9*Q11-Q7*Q11^2+Q6*Q9*Q12+Q6* Q11*Q13,
           Q8*Q9*Q11+Q7*Q11^2-Q6*Q9*Q12-Q6*Q11*Q13,
           Q8*Q9*Q11-Q7*Q11^2+Q6*Q9*Q12-Q6*Q11*Q13,
           -Q4*Q8*Q9*Q12+Q5*Q8*Q9*Q13+2* Q4*Q6*Q12*Q13-2*Q3*Q7*Q12*Q13+2*Q1*Q12^2*Q13- Q2*Q11*Q13^2-Q4*Q7*Q11*Q14+Q3*Q7*Q13*Q14-Q1*Q12*Q13*Q15,
           -Q4*Q5*Q6+Q3*Q4*Q8+Q2*Q4*Q11-Q2*Q3* Q13,
           -Q8*Q9*Q11+Q7*Q11^2-Q6*Q9*Q12+Q6* Q11*Q13,
           Q8*Q9*Q11-Q7*Q11^2-Q6*Q9*Q12+Q6* Q11*Q13,
           -Q2*Q8^2*Q9^2-Q5*Q6*Q7^2*Q11+2*Q2*Q7*Q8*Q9*Q11-Q2*Q7^2*Q11^2+Q2*Q6*Q7*Q9* Q12-Q1*Q8^2*Q9*Q12+Q2*Q6*Q8*Q9*Q13-Q2*Q6*Q7*Q11*Q13+Q1*Q7*Q8*Q11*Q14-Q1*Q8^2* Q9*Q15+Q1*Q7*Q8*Q11*Q15,
           Q7*Q11^2*Q12*Q13-Q6*Q11*Q12*Q13^2+Q7*Q11^2* Q13*Q14-Q6*Q9*Q12*Q13*Q14-Q6*Q11*Q13^2*Q14+Q8*Q9^2*Q12*Q15-Q7*Q9*Q11*Q12*Q15-Q7*Q11^2* Q13*Q15+Q6*Q9*Q12*Q13*Q15+Q6*Q11*Q13^2*Q15+ Q6*Q9*Q13*Q14*Q15-Q8*Q9^2*Q15^2,
           -Q4^2*Q6*Q8+Q3*Q4*Q7*Q8-Q2*Q4*Q7*Q11+Q1*Q5*Q7*Q12+Q2*Q4*Q6*Q13-2*Q1*Q5* Q8*Q13+Q1*Q2*Q12*Q13+Q1*Q4*Q8*Q15-Q1*Q2*Q13*Q15,
           -Q3*Q7*Q8^2*Q9+Q3*Q7^2*Q8*Q11+ Q1*Q8^2*Q9*Q12-Q1*Q7*Q8*Q11*Q12-Q3*Q6*Q7*Q8*Q13+2*Q1*Q8^2*Q11*Q13+Q3*Q6*Q7^2*Q14-Q1*Q6*Q7*Q12*Q14-Q1*Q8^2*Q9*Q15-Q1* Q6*Q8*Q13*Q15+Q1*Q6*Q7*Q14*Q15,
           Q5^2*Q6^2*Q12-Q3*Q5*Q6*Q8*Q12-Q2^2*Q11^2* Q12+Q3*Q5*Q6*Q8*Q14-Q3^2*Q8^2*Q14+Q2*Q3* Q6*Q12*Q14-Q5^2*Q6^2*Q15+Q3^2*Q8^2*Q15+Q2*Q3*Q8*Q11*Q15+Q2^2*Q11^2*Q15-Q2*Q3*Q6*Q14* Q15-Q2*Q3*Q6*Q15^2,
           Q4*Q5*Q7*Q8*Q12+Q4^2*Q8^2*Q12-2*Q4*Q5* Q8^2*Q13+Q2*Q5*Q7*Q12*Q13+3*Q2*Q4*Q8*Q12*Q13-2*Q2*Q5*Q8*Q13^2+Q2^2*Q12*Q13^2-Q4^2*Q8^2* Q14-Q2*Q4*Q7*Q12*Q14-Q2*Q4*Q7*Q12*Q15-Q2^2*Q13^2*Q15+Q2*Q4*Q7*Q14*Q15,
           -Q3*Q4^2*Q8^2-Q2*Q4*Q5*Q7*Q11+Q2*Q3*Q4*Q8*Q13-Q1*Q5^2*Q8*Q13+Q1*Q5^2*Q7*Q14+2* Q1*Q4*Q5*Q8*Q14-Q1*Q2*Q4*Q12*Q14-Q1*Q5^2*Q7*Q15+Q1*Q2*Q4*Q12*Q15+Q1*Q2*Q5* Q13*Q15-Q1*Q2*Q4*Q14*Q15,
           -2*Q5^2*Q6*Q8*Q11-2*Q3*Q5*Q8^2*Q11-Q5^2*Q6^2* Q12+Q5^2*Q6^2*Q14+3*Q3*Q5*Q6*Q8*Q14+Q3^2* Q8^2*Q14+Q2*Q5*Q6*Q11*Q14+Q2*Q3*Q8*Q11*Q14-Q2*Q3*Q6*Q12*Q14-Q3^2*Q8^2*Q15+Q2*Q3* Q6*Q12*Q15-Q2*Q3*Q6*Q14*Q15,
           -Q4*Q8^2*Q11+Q5*Q6*Q8*Q13-Q5*Q6*Q7*Q14+2*Q3*Q7*Q8*Q14-2*Q2*Q8*Q9*Q14-Q1*Q8* Q12*Q14+2*Q1*Q8*Q14^2+Q2*Q8*Q9*Q15-Q2*Q7*Q11*Q15,
           -Q5^2*Q7^2*Q9-Q4*Q5*Q7^2*Q11+Q3* Q5*Q7^2*Q13-Q2*Q4*Q7*Q11*Q13-Q1*Q5*Q7*Q12*Q13+Q1*Q4*Q8*Q12*Q13-Q1*Q2*Q12*Q13^2+ Q3*Q4*Q7^2*Q14-Q1*Q4*Q7*Q12*Q15+2*Q1*Q5*Q7*Q13*Q15+Q1*Q2*Q13^2*Q15,
           -Q2*Q5^2*Q9^2-Q3*Q4^2*Q8*Q11+2*Q2*Q4*Q5*Q9*Q11-Q2*Q4^2*Q11^2-Q1*Q5^2*Q9*Q12+Q1*Q4*Q5*Q11*Q12+Q2*Q3*Q5*Q9*Q13-Q2*Q3*Q4*Q11*Q13+Q1*Q4*Q5*Q11*Q14+Q2*Q3*Q4*Q9*Q15-Q1*Q5^2*Q9*Q15,
           -Q5^2*Q6-Q3*Q5*Q8+Q2*Q5*Q11+Q2*Q3*Q15,
           Q5^2*Q6-Q3*Q5*Q8+Q2*Q5*Q11-Q2*Q3*Q15,
           -Q5^2*Q6+Q3*Q5*Q8+Q2*Q5*Q11-Q2*Q3* Q15,
           -Q5^2*Q6+Q3*Q5*Q8-Q2*Q5*Q11+Q2*Q3* Q15,
           Q5^2*Q6+Q3*Q5*Q8-Q2*Q5*Q11-Q2*Q3*Q15,
           -Q5*Q6*Q7*Q9*Q11+Q3*Q6*Q8*Q9* Q13-Q1*Q8*Q11^2*Q13-Q3*Q6^2*Q13^2+Q1*Q7*Q11^2* Q14-Q1*Q6*Q9*Q12*Q14+2*Q1*Q6*Q11*Q13*Q14+ Q1*Q8*Q9*Q11*Q15-Q1*Q7*Q11^2*Q15+Q1*Q6*Q9*Q12*Q15-Q1*Q6*Q9*Q14*Q15,
           Q3*Q4*Q6-Q1*Q5*Q11,
           Q4*Q7*Q8*Q11^2+Q4*Q6*Q8*Q9*Q12-Q4*Q6^2*Q12*Q13-Q1*Q8*Q11*Q12*Q13+Q1*Q8*Q9* Q12*Q14-Q1*Q8*Q11*Q13*Q14-Q1*Q6*Q12*Q13*Q14-Q4*Q6*Q7*Q11*Q15-Q1*Q8*Q9*Q12*Q15+Q1*Q8*Q11*Q13*Q15+Q1*Q6*Q12*Q13*Q15+Q1*Q7*Q11*Q15^2,
           -Q5*Q7*Q8*Q9-Q4*Q8^2*Q9-Q4* Q7*Q8*Q11+2*Q4*Q6*Q8*Q13-Q3*Q7*Q8*Q13+Q2*Q7*Q11*Q13+2*Q1*Q8*Q12*Q13-Q2*Q6*Q13^2+ Q4*Q6*Q7*Q14-Q3*Q7^2*Q14+Q1*Q7*Q12*Q14- Q1*Q7*Q14*Q15,
           -Q5*Q7^2*Q11^2-Q4*Q6*Q8*Q9*Q12+Q4*Q6*Q8*Q11*Q13+Q2*Q6*Q9*Q12*Q13+Q1*Q8*Q11* Q12*Q13-Q2*Q6*Q11*Q13^2+2*Q4*Q6*Q7*Q11*Q14-Q4*Q6^2*Q13*Q14+Q2*Q6*Q9*Q13*Q14-Q1*Q8* Q11*Q13*Q14-Q1*Q6*Q12*Q13*Q14-Q1*Q6*Q13*Q14^2+Q1*Q7*Q11*Q14*Q15,
           Q2*Q5*Q6*Q11*Q12-Q2^2*Q11^2*Q12+Q5^2*Q6^2* Q14-Q3^2*Q8^2*Q14-Q2*Q5*Q6*Q11*Q14+Q2*Q3* Q6*Q12*Q14-Q5^2*Q6^2*Q15+Q3^2*Q8^2*Q15+Q2*Q3*Q8*Q11*Q15+Q2^2*Q11^2*Q15-Q2*Q3*Q6*Q12* Q15-Q2*Q3*Q6*Q15^2,
           Q5*Q9-Q3*Q13,
           -Q5*Q9+Q3*Q13,
           -Q5*Q6*Q11+Q3*Q8*Q11-Q2*Q11^2+Q3*Q6*Q14,
           -Q5*Q6*Q11+Q3*Q8*Q11+Q2*Q11^2-Q3*Q6* Q14,
           Q2*Q5*Q7*Q9-Q2*Q4*Q7*Q11+Q2*Q4* Q6*Q13-2*Q1*Q5*Q8*Q13-Q2^2*Q9*Q13-Q1*Q5*Q7*Q14+Q1*Q2*Q13*Q14+Q1*Q5*Q7*Q15+Q1* Q4*Q8*Q15,
           Q5*Q6*Q11-Q3*Q8*Q11-Q2*Q11^2+ Q3*Q6*Q14,
           Q8*Q9*Q12-Q8*Q11*Q13+Q6*Q12* Q13-Q8*Q9*Q15,
           Q5*Q6*Q11-Q3*Q8*Q11+Q2*Q11^2-Q3*Q6*Q14,
           Q5*Q6*Q11+Q3*Q8*Q11-Q2*Q11^2-Q3*Q6*Q14,
           -Q4^2*Q11^2*Q12-2*Q4*Q5*Q11^2*Q13-2*Q3*Q5* Q11*Q13^2+Q4*Q5*Q9*Q11*Q14+Q4^2*Q11^2*Q14-Q3*Q4*Q9*Q12*Q14+Q3*Q5*Q9*Q13*Q14+3*Q3*Q4* Q11*Q13*Q14+Q3^2*Q13^2*Q14+Q3*Q4*Q9*Q12*Q15- Q3^2*Q13^2*Q15-Q3*Q4*Q9*Q14*Q15,
           -Q5*Q8*Q9+Q5*Q6*Q13+Q3*Q7*Q14-Q1*Q12* Q14,
           -Q4*Q7*Q11+Q1*Q13*Q14,
           Q4*Q7*Q11-Q1*Q13*Q14,
           Q4^2*Q5*Q6*Q9*Q12-Q2*Q4^2*Q9*Q11*Q12- Q4*Q5^2*Q6*Q9*Q13+Q2*Q5^2*Q9^2*Q13-Q2*Q3*Q4*Q9*Q12*Q13+Q1*Q5^2*Q9*Q12*Q13+Q1*Q4* Q5*Q11*Q12*Q13-2*Q2*Q3*Q5*Q9*Q13^2+Q2*Q3* Q4*Q11*Q13^2+Q2*Q3^2*Q13^3-Q1*Q3*Q5*Q13^2* Q14+Q4^2*Q5*Q6*Q9*Q15-Q2*Q4*Q5*Q9^2*Q15- Q1*Q4*Q5*Q9*Q12*Q15+2*Q2*Q3*Q4*Q9*Q13* Q15+Q1*Q5^2*Q9*Q13*Q15-Q1*Q4*Q5*Q11*Q13*Q15-Q1*Q3*Q5*Q13^2*Q15+Q1*Q3*Q4*Q13* Q14*Q15-Q1*Q4*Q5*Q9*Q15^2,
           Q4*Q5*Q6^2-Q3*Q4*Q6*Q8-Q2*Q5*Q6*Q9-2*Q1*Q5*Q8*Q11+Q1*Q5*Q6*Q12-Q1*Q2* Q11*Q12+Q2*Q3*Q6*Q13+Q1*Q5*Q6*Q14+Q1*Q2*Q11*Q15,
           Q3*Q5*Q7*Q8*Q11-Q3*Q4*Q8^2*Q11-Q2^2*Q11^2*Q13-Q3*Q5*Q6*Q7*Q14+Q3*Q4* Q6*Q8*Q14+Q1*Q5*Q8*Q11*Q14+Q3*Q4*Q6*Q8*Q15-Q3^2*Q7*Q8*Q15+2*Q2*Q3*Q7*Q11*Q15- Q1*Q5*Q8*Q11*Q15+Q1*Q2*Q11*Q12*Q15-Q1*Q3*Q8*Q14*Q15-Q1*Q3*Q8*Q15^2,
           Q4*Q6*Q9+Q1*Q9*Q12-Q1*Q11*Q13-Q1*Q9* Q15,
           Q3*Q8^2*Q9*Q11+Q4*Q6*Q8*Q11^2-2*Q2* Q8*Q9*Q11^2+Q2*Q7*Q11^3-Q5*Q6^2*Q9*Q12+2*Q4* Q6^2*Q11*Q12-Q2*Q6*Q9*Q11*Q12+Q2*Q6*Q11^2*Q13+Q1*Q6*Q11*Q12*Q14-Q3*Q6*Q8*Q9* Q15-2*Q4*Q6^2*Q11*Q15+2*Q2*Q6*Q9*Q11*Q15-2*Q1* Q8*Q11^2*Q15-Q1*Q6*Q11*Q12*Q15-Q1*Q6*Q11* Q14*Q15+2*Q1*Q6*Q11*Q15^2,
           Q6*Q12*Q13-Q7*Q11*Q14,
           -Q4*Q8^2*Q12*Q13+Q2*Q8*Q12*Q13^2-Q4*Q7*Q8*Q12*Q14-Q4*Q8^2*Q13*Q14+Q2*Q8*Q13^2*Q14+ Q4*Q7*Q8*Q12*Q15+Q4*Q8^2*Q13*Q15-Q2*Q8*Q13^2*Q15+Q5*Q7^2*Q14*Q15+Q4*Q7*Q8*Q14*Q15- Q2*Q7*Q13*Q14*Q15-Q5*Q7^2*Q15^2,
           Q1*Q5*Q7*Q8*Q12+Q1*Q4*Q8^2*Q12-Q2*Q5*Q6*Q7*Q13+Q2*Q4*Q6*Q8*Q13-2*Q1*Q5* Q8^2*Q13+3*Q1*Q2*Q8*Q12*Q13-Q2^2*Q6*Q13^2+Q2* Q4*Q6*Q7*Q14-Q1*Q4*Q8^2*Q14-Q1*Q2*Q7* Q12*Q15-Q1*Q2*Q8*Q13*Q15,
           Q5^2*Q6*Q7-2*Q3*Q5*Q7*Q8+Q3*Q4*Q8^2+ Q2*Q4*Q8*Q11-Q2*Q4*Q6*Q12+2*Q1*Q5*Q8*Q12-Q2^2*Q11*Q13+Q2*Q3*Q7*Q14-Q1*Q2*Q12* Q14-2*Q1*Q5*Q8*Q15+Q1*Q2*Q14*Q15,
           -Q5^2*Q7*Q9+Q4*Q5*Q7*Q11-Q4^2*Q6*Q12+Q2*Q4*Q9*Q12-Q4*Q5*Q6*Q13-Q3*Q4*Q8* Q13+2*Q2*Q5*Q9*Q13-Q2*Q4*Q11*Q13-Q2*Q3*Q13^2-Q1*Q4*Q12*Q14+Q1*Q4*Q12*Q15+2*Q1*Q5* Q13*Q15,
           -Q3*Q7+Q2*Q9,
           Q3*Q7-Q2*Q9,
           -Q5^2*Q6^2*Q9+Q3*Q5*Q6*Q8*Q9-Q1*Q5*Q8*Q11^2+2*Q1*Q5*Q6*Q11*Q12+Q1*Q2*Q11^2*Q12- Q2*Q3*Q6*Q11*Q13-Q1*Q3*Q6*Q12*Q14+Q1*Q3*Q8*Q11*Q15-Q1*Q2*Q11^2*Q15-Q1*Q3*Q6* Q12*Q15+Q1*Q3*Q6*Q14*Q15,
           -Q2*Q12^2+Q5*Q8*Q14,
           Q5*Q8^2*Q9*Q11-Q5*Q7*Q8*Q11^2+Q5*Q6* Q8*Q9*Q12-Q5*Q6*Q7*Q11*Q12-Q2*Q8*Q9*Q11*Q12+Q2*Q7*Q11^2*Q12+Q5*Q6^2*Q12*Q13-Q5* Q6^2*Q13*Q14+Q5*Q6*Q8*Q9*Q15-Q2*Q8*Q9*Q11*Q15+Q2*Q7*Q11^2*Q15-Q2*Q6*Q9*Q12*Q15,
           Q5^2*Q7-Q2*Q4*Q12,
           -Q5*Q8*Q9*Q12+Q5*Q8*Q11*Q13+Q5*Q6*Q12* Q13-Q3*Q8*Q12*Q13+Q2*Q11*Q12*Q13-Q5*Q8* Q9*Q14-Q4*Q6*Q12*Q14+Q3*Q7*Q12*Q14-Q1* Q12^2*Q14+Q5*Q6*Q13*Q14+Q3*Q7*Q14^2-Q1*Q12*Q14^2-Q3*Q8*Q13*Q15-Q2*Q11*Q13*Q15+Q1* Q12*Q14*Q15+Q1*Q14^2*Q15,
           Q8*Q9*Q12-Q8*Q11*Q13+Q6*Q12*Q13-Q6*Q13* Q14,
           Q5^2*Q6*Q9-Q3*Q5*Q8*Q9+Q2*Q4*Q11^2-2*Q1*Q5*Q11*Q12+Q2*Q3*Q11*Q13+Q1*Q3* Q12*Q14+Q3*Q4*Q6*Q15-Q3^2*Q7*Q15+Q1*Q3*Q12*Q15-Q1*Q3*Q14*Q15-Q1*Q3*Q15^2,
           Q5^2*Q6*Q9+Q3*Q5*Q8*Q9+Q2*Q4*Q11^2-2* Q1*Q5*Q11*Q12-Q2*Q3*Q11*Q13+Q1*Q3*Q12*Q14+Q3*Q4*Q6*Q15-Q3^2*Q7*Q15+Q1*Q3*Q12* Q15-Q1*Q3*Q14*Q15-Q1*Q3*Q15^2,
           -Q5*Q6*Q13-Q3*Q8*Q13+Q2*Q11*Q13+Q3*Q7* Q15-Q1*Q12*Q15+Q1*Q14*Q15,
           Q4^2*Q5*Q8*Q9-Q4^3*Q8*Q11+Q3*Q4^2*Q8* Q13-Q2*Q4^2*Q11*Q13+Q2*Q3*Q4*Q13^2-Q1*Q5^2* Q13^2+Q1*Q4^2*Q12*Q14-Q1*Q4*Q5*Q13*Q14,
           Q2*Q9*Q12^2+Q5*Q8*Q11*Q13+Q5*Q6*Q12*Q13-Q3*Q8*Q12*Q13-Q5*Q8*Q9*Q14+Q5*Q7*Q11*Q14-Q4*Q6*Q12*Q14+Q2*Q9*Q12*Q14-Q1*Q12^2*Q14+Q5*Q6*Q13*Q14-Q3*Q8*Q13*Q14-Q1* Q12*Q14^2-Q5*Q8*Q9*Q15-Q5*Q7*Q11*Q15+Q1*Q12^2*Q15+Q1*Q12*Q14*Q15,
           -Q4*Q8+Q2*Q13,
           Q4*Q8-Q2*Q13,
           -Q6*Q13*Q14+Q8*Q9*Q15,
           Q6*Q13*Q14-Q8*Q9*Q15,
           Q3*Q8*Q9-Q1*Q11*Q12,
           Q5^2*Q6*Q9^2+Q3*Q5*Q8*Q9^2-Q4*Q5*Q6* Q9*Q11-Q2*Q5*Q9^2*Q11+Q4^2*Q6*Q11^2+Q3*Q4* Q6*Q9*Q12-Q1*Q5*Q9*Q11*Q12-Q1*Q4*Q11^2* Q12-Q3^2*Q6*Q13^2-Q1*Q5*Q9*Q11*Q14+Q1* Q4*Q11^2*Q14+2*Q1*Q3*Q11*Q13*Q14-Q3*Q4*Q6*Q9*Q15+Q1*Q5*Q9*Q11*Q15-Q1*Q4*Q11^2*Q15- Q1*Q3*Q9*Q14*Q15+Q1*Q3*Q9*Q15^2,
           Q5*Q8*Q9^2*Q12-Q5*Q7*Q9*Q11*Q12+Q4*Q7*Q11^2*Q12-Q5*Q7*Q11^2*Q13+Q5*Q6*Q9*Q12* Q13-Q4*Q6*Q11*Q12*Q13+Q5*Q6*Q11*Q13^2+Q4*Q7*Q11^2*Q14-Q4*Q6*Q9*Q12*Q14+Q5*Q6*Q9* Q13*Q14-Q4*Q6*Q11*Q13*Q14-Q5*Q8*Q9^2*Q15,
           -Q5*Q8*Q9*Q12+Q5*Q8*Q11*Q13-Q5*Q6*Q12* Q13-Q5*Q8*Q9*Q14+Q5*Q7*Q11*Q14-Q4*Q8* Q11*Q14+Q5*Q6*Q13*Q14+Q5*Q7*Q11*Q15-Q4* Q8*Q11*Q15-Q3*Q7*Q14*Q15+Q2*Q9*Q14*Q15+ Q1*Q12*Q14*Q15-Q1*Q14^2*Q15+Q2*Q9*Q15^2+Q1*Q12*Q15^2-Q1*Q14*Q15^2,
           -Q3*Q4*Q8-Q2*Q5*Q9+Q2*Q4*Q11+Q2*Q3* Q13,
           Q8^2*Q9^2*Q12^2+Q7*Q8*Q11^2*Q12*Q13+Q6*Q8*Q9*Q12^2*Q13-2*Q8^2*Q11^2*Q13^2+Q6*Q8*Q11* Q12*Q13^2+Q7^2*Q11^2*Q12*Q14-Q6*Q7*Q9*Q12^2*Q14- Q8^2*Q9*Q11*Q13*Q14-2*Q7*Q8*Q11^2*Q13*Q14+ Q6*Q8*Q11*Q13^2*Q14-Q6^2*Q12*Q13^2*Q14+Q6*Q7* Q11*Q13*Q14^2-3*Q7*Q8*Q9*Q11*Q12*Q15-Q7^2*Q11^2*Q12*Q15+Q6*Q7*Q9*Q12^2*Q15+2*Q8^2*Q9*Q11* Q13*Q15+2*Q7*Q8*Q11^2*Q13*Q15-Q6*Q8*Q9*Q12*Q13*Q15-Q6*Q7*Q11*Q12*Q13*Q15+Q7^2*Q11^2* Q14*Q15,
           -Q2*Q5^2*Q6*Q9-Q2*Q3*Q5*Q8*Q9+2* Q1*Q5^2*Q8*Q11+Q2^2*Q5*Q9*Q11-Q1*Q5^2*Q6*Q12-Q1*Q3*Q5*Q8*Q12+Q1*Q5^2*Q6*Q14- Q1*Q2*Q5*Q11*Q14+Q2^2*Q3*Q9*Q15+Q1*Q2*Q3*Q12*Q15-Q1*Q2*Q3*Q14*Q15,
           Q5*Q7*Q13-Q4*Q8*Q13+Q2*Q13^2-Q4*Q7*Q14,
           Q5*Q7*Q13+Q4*Q8*Q13-Q2*Q13^2-Q4*Q7*Q14,
           Q5*Q7*Q13-Q4*Q8*Q13-Q2*Q13^2+Q4*Q7*Q14,
           -Q5*Q7*Q13+Q4*Q8*Q13+Q2*Q13^2-Q4*Q7* Q14,
           -Q5*Q7*Q13+Q4*Q8*Q13-Q2*Q13^2+Q4*Q7* Q14,
           Q5^2*Q6*Q9*Q12-Q2*Q5*Q9*Q11*Q12-Q5^2*Q6*Q11*Q13+Q2*Q5*Q11^2*Q13+Q2*Q3*Q11* Q12*Q13-Q3^2*Q8*Q13*Q14+Q5^2*Q6*Q9*Q15-Q2*Q5*Q9*Q11*Q15-Q2*Q3*Q9*Q12*Q15-Q3*Q5* Q6*Q13*Q15+Q3^2*Q8*Q13*Q15+Q2*Q3*Q11*Q13*Q15,
           -Q3*Q5*Q6*Q8*Q9-Q3^2*Q8^2*Q9+Q2*Q3*Q8*Q9*Q11-2*Q1*Q5*Q8*Q11^2+Q2*Q3*Q6* Q9*Q12-Q1*Q2*Q11^2*Q12-Q1*Q3*Q8*Q11*Q14+Q1*Q5*Q6*Q11*Q15+3*Q1*Q3*Q8*Q11*Q15+Q1* Q2*Q11^2*Q15-Q1*Q3*Q6*Q14*Q15,
           Q2*Q4*Q8*Q12*Q13-Q2^2*Q12*Q13^2-Q5^2*Q7^2* Q14+Q4^2*Q8^2*Q14+Q2*Q4*Q7*Q12*Q14-Q2*Q4* Q8*Q13*Q14+Q5^2*Q7^2*Q15-Q4^2*Q8^2*Q15-Q2*Q4*Q7*Q12*Q15+Q2*Q5*Q7*Q13*Q15+Q2^2*Q13^2* Q15-Q2*Q4*Q7*Q15^2,
           2*Q5^2*Q8*Q11^2+Q5^2*Q6*Q11*Q12-Q3*Q5*Q8* Q11*Q12-Q5^2*Q6*Q11*Q14-2*Q3*Q5*Q8*Q11*Q14- Q2*Q5*Q11^2*Q14+Q2*Q3*Q11*Q12*Q14+Q3^2*Q8*Q14^2-Q2*Q3*Q11*Q12*Q15+Q2*Q3*Q11*Q14* Q15,
           Q5^2*Q7*Q12+Q4*Q5*Q8*Q12+Q2*Q4*Q12^2-2* Q5^2*Q8*Q13+Q2*Q5*Q12*Q13-Q4*Q5*Q8*Q14- Q2*Q5*Q13*Q15,
           Q5^2*Q8*Q9^2+Q4*Q5*Q8*Q9*Q11+Q4^2*Q8* Q11^2-Q3*Q4*Q8*Q9*Q12-Q3^2*Q8*Q13^2-Q3*Q4* Q8*Q9*Q14-Q1*Q4*Q11*Q12*Q14+Q1*Q5*Q11* Q13*Q14-Q3*Q4*Q8*Q9*Q15-Q1*Q5*Q9*Q12* Q15+Q1*Q5*Q11*Q13*Q15+Q1*Q3*Q13*Q14*Q15,
           -Q7*Q8*Q11*Q12*Q13+2*Q8^2*Q11*Q13^2+Q6*Q8*Q12*Q13^2+Q7*Q8*Q9*Q12*Q14-Q8^2*Q9*Q13*Q14-2* Q7*Q8*Q11*Q13*Q14-Q6*Q8*Q13^2*Q14+Q7^2*Q11*Q14^2-Q7*Q8*Q9*Q12*Q15+Q7*Q8*Q9*Q14* Q15,
           Q4*Q6*Q7-Q3*Q7^2+Q1*Q7*Q12-Q1*Q8* Q13,
           -Q5*Q8*Q11*Q12-Q5*Q8*Q11*Q14+Q2*Q11*Q12*Q14+Q5*Q6*Q14^2+Q3*Q8*Q14^2-Q5*Q8*Q11* Q15-Q2*Q11*Q12*Q15+Q2*Q11*Q14*Q15,
           Q4*Q8*Q11*Q13-Q5*Q6*Q13^2-Q4*Q7*Q11*Q14+2*Q3*Q7*Q13*Q14-2*Q2*Q9*Q13*Q14-Q1*Q12* Q13*Q14+2*Q1*Q13*Q14^2-Q5*Q7*Q9*Q15+Q2*Q9*Q13*Q15,
           Q4^2*Q8*Q12^2-Q5^2*Q7*Q12*Q13-2*Q4*Q5*Q8*Q12*Q13+2*Q5^2*Q8*Q13^2-Q2*Q5*Q12*Q13^2+ Q4*Q5*Q7*Q12*Q14+Q4*Q5*Q7*Q12*Q15-Q4* Q5*Q8*Q13*Q15+Q2*Q5*Q13^2*Q15-Q4*Q5*Q7* Q14*Q15,
           Q3*Q8*Q9+Q3*Q7*Q11-Q1*Q11*Q12- Q3*Q6*Q13,
           -Q3*Q8*Q9+Q3*Q7*Q11-Q1*Q11*Q12+Q3*Q6*Q13,
           -Q4*Q8*Q9*Q12-Q5*Q7*Q11*Q13+Q4*Q8*Q11* Q13+Q2*Q9*Q12*Q13-Q2*Q11*Q13^2+Q4*Q7*Q11*Q14-Q4*Q6*Q13*Q14+Q2*Q9*Q13*Q14-Q1*Q12*Q13*Q14-Q1*Q13*Q14^2+Q5*Q7*Q9*Q15+Q1* Q12*Q13*Q15,
           Q5^2*Q6*Q9*Q11-2*Q4*Q5*Q6*Q11^2+ Q3*Q5*Q7*Q11^2+Q3*Q4*Q8*Q11^2+Q2*Q4*Q11^ 3-Q3*Q5*Q6*Q9*Q12+2*Q3*Q4*Q6*Q11*Q12-2* Q3^2*Q7*Q11*Q12-2*Q1*Q5*Q11^2*Q12+2*Q1*Q3*Q11* Q12^2-Q3*Q4*Q6*Q11*Q14+2*Q3^2*Q7*Q11*Q14- Q1*Q3*Q11*Q12*Q14-Q3^2*Q6*Q13*Q14-Q1*Q3*Q11*Q12*Q15+Q1*Q3*Q11*Q14*Q15,
           Q5*Q8*Q9-Q3*Q7*Q14,
           -Q4*Q5*Q6*Q9+Q3*Q4*Q8*Q9+Q2*Q5*Q9^ 2+Q1*Q5*Q9*Q12-Q1*Q4*Q11*Q12-Q2*Q3*Q9*Q13-2*Q1*Q5*Q11*Q13+Q1*Q4*Q11*Q14+Q1* Q5*Q9*Q15,
           Q5*Q6*Q7*Q9-Q3*Q7*Q8*Q9+ Q2*Q8*Q9^2+Q1*Q8*Q9*Q12-Q2*Q6*Q9*Q13-2*Q1*Q8*Q11*Q13+Q1*Q7*Q11*Q14+Q1*Q8*Q9* Q15-Q1*Q7*Q11*Q15,
           Q2*Q6*Q7-Q1*Q8^2,
           -Q5*Q8*Q9-Q4*Q8*Q11+Q1*Q12^2+Q3*Q8*Q13,
           -Q8*Q11*Q12*Q13-Q8*Q11*Q13*Q14-Q6*Q12* Q13*Q14+Q8*Q9*Q12*Q15+Q8*Q11*Q13*Q15+Q6*Q12*Q13*Q15+Q7*Q11*Q14*Q15+Q6*Q13*Q14*Q15-Q8*Q9*Q15^2-Q7*Q11*Q15^2,
           Q5*Q7*Q9-Q4*Q8*Q9+Q4*Q7*Q11-Q3*Q7* Q13,
           -Q5^2*Q6*Q9+Q3*Q5*Q8*Q9+2*Q4*Q5*Q6*Q11-Q3*Q5*Q7*Q11-Q3*Q4*Q8*Q11-Q2* Q4*Q11^2+2*Q1*Q5*Q11*Q12-Q2*Q3*Q11*Q13+Q3*Q4*Q6*Q14-Q3^2*Q7*Q14+Q1*Q3*Q12*Q14-Q1* Q3*Q14*Q15,
           -Q5*Q9*Q12-Q5*Q11*Q13+Q5*Q9*Q15+Q3*Q13*Q15,
           -Q5*Q11*Q12*Q13+Q5*Q9*Q12*Q14-Q5*Q11*Q13* Q14+Q4*Q11*Q14^2+Q3*Q13*Q14^2-Q5*Q9*Q12* Q15-Q5*Q11*Q13*Q15+Q5*Q9*Q14*Q15,
           Q8*Q11-Q6*Q14,
           -Q5*Q8*Q9^3+Q5*Q7*Q9^2*Q11-Q4*Q8*Q9^2*Q11+Q4*Q7*Q9*Q11^2+Q5*Q6*Q9^2*Q13-Q1*Q11^2* Q13^2+Q1*Q9^2*Q12*Q15-Q1*Q9*Q11*Q13*Q15,
           Q5*Q11*Q12-Q3*Q14^2,
           Q5*Q6*Q9-Q1*Q11*Q12,
           -Q5*Q6*Q9+Q1*Q11*Q12,
           Q7*Q12-Q8*Q13,
           -Q5*Q8*Q13+Q5*Q7*Q14+Q4*Q8*Q14-Q5*Q7* Q15,
           Q5*Q7*Q9-Q1*Q13*Q14,
           -Q7*Q12+Q8*Q13,
           Q5^2*Q7*Q9*Q13+Q4*Q5*Q6*Q13^2-2*Q3*Q5* Q7*Q13^2+Q3*Q4*Q8*Q13^2+Q2*Q3*Q13^3-Q4^2*Q7* Q11*Q14+2*Q4^2*Q6*Q13*Q14-Q3*Q4*Q7*Q13* Q14+Q1*Q4*Q12*Q13*Q14-Q4*Q5*Q7*Q9*Q15-2*Q4^2*Q6*Q13*Q15+2*Q3*Q4*Q7*Q13*Q15-Q1*Q4* Q12*Q13*Q15-2*Q1*Q5*Q13^2*Q15-Q1*Q4*Q13*Q14* Q15+2*Q1*Q4*Q13*Q15^2,
           -Q3*Q8*Q9+Q1*Q11*Q15,
           Q2*Q4*Q11*Q12-Q5^2*Q6*Q13+Q3*Q5*Q8*Q13-Q2*Q5*Q11*Q13+Q3*Q5*Q7*Q14-Q3*Q4*Q8*Q14-Q2*Q5*Q9*Q14+Q1*Q5*Q14^2-Q2*Q5* Q9*Q15-Q1*Q5*Q12*Q15+Q2*Q3*Q13*Q15+Q1*Q5*Q14*Q15,
           Q3*Q8*Q9-Q1*Q11*Q15,
           -Q4*Q11*Q12+Q5*Q9*Q15,
           -Q4^2*Q6*Q8*Q12+Q5^2*Q6*Q7*Q13+Q2*Q4*Q6*Q12*Q13-Q1*Q5*Q8*Q12*Q13-Q1*Q4*Q8* Q12*Q14-Q1*Q5*Q8*Q13*Q14+Q1*Q2*Q12*Q13*Q14-Q4*Q5*Q6*Q7*Q15+Q1*Q4*Q8*Q12*Q15+Q1*Q5*Q8*Q13*Q15-Q1*Q2*Q12*Q13*Q15+Q1*Q5*Q7*Q15^2,
           Q5*Q6*Q7*Q9-Q4*Q6*Q8*Q9+Q2*Q8*Q9^2-Q2*Q7*Q9*Q11+Q1*Q8*Q9*Q12-2* Q1*Q8*Q11*Q13-Q1*Q6*Q12*Q13+Q1*Q6*Q13*Q14+Q1*Q8*Q9*Q15,
           -Q5^2*Q7*Q9-Q4*Q5*Q8*Q9-Q4*Q5*Q7*Q11-Q4^2*Q6*Q12+Q2*Q4*Q9*Q12-Q4*Q5*Q6* Q13+Q3*Q4*Q8*Q13+2*Q2*Q5*Q9*Q13-Q2*Q3*Q13^2-Q1*Q4*Q12*Q14+Q1*Q4*Q12*Q15+2*Q1*Q5* Q13*Q15,
           -Q4*Q7*Q8*Q9*Q12+Q4*Q8^2*Q9*Q13+2* Q4*Q6*Q7*Q12*Q13-2*Q3*Q7^2*Q12*Q13+2*Q1* Q7*Q12^2*Q13+Q5*Q6*Q7*Q13^2-2*Q4*Q6*Q8*Q13^2+ Q3*Q7*Q8*Q13^2-2*Q1*Q8*Q12*Q13^2+Q2*Q6* Q13^3-Q4*Q7^2*Q11*Q14-Q4*Q6*Q7*Q13*Q14+2*Q3* Q7^2*Q13*Q14-Q1*Q7*Q12*Q13*Q14-Q1*Q7*Q12* Q13*Q15+Q1*Q7*Q13*Q14*Q15,
           Q3*Q7*Q9-Q1*Q11*Q13,
           -Q5^2*Q8*Q9+2*Q4*Q5*Q6*Q12-2*Q3*Q5*Q7*Q12-Q2*Q4*Q11*Q12+2*Q1*Q5*Q12^2+Q2*Q5*Q11* Q13+Q3*Q5*Q7*Q14-Q3*Q4*Q8*Q14-Q1*Q5*Q12*Q15,
           Q5*Q6^2*Q8*Q9+Q5*Q6^2*Q7*Q11+Q2* Q6*Q7*Q11^2-Q1*Q8^2*Q11^2-Q5*Q6^3*Q13-Q2* Q6^2*Q11*Q13-Q1*Q6*Q8*Q11*Q14+Q1*Q6^2*Q12* Q14,
           -Q5^2*Q6*Q7*Q13+Q3*Q5*Q7*Q8*Q13-Q2^2*Q11*Q13^2+Q4*Q5*Q6*Q7*Q14-Q3*Q4*Q7* Q8*Q14+Q1*Q5*Q8*Q13*Q14+Q4*Q5*Q6*Q7*Q15-Q3*Q5*Q7^2*Q15+2*Q2*Q3*Q7*Q13*Q15-Q1* Q5*Q8*Q13*Q15+Q1*Q2*Q12*Q13*Q15-Q1*Q5*Q7*Q14*Q15-Q1*Q5*Q7*Q15^2,
           -Q4*Q5^2*Q6^2+Q3*Q4*Q5*Q6*Q8+2*Q2*Q4*Q5*Q6*Q11-Q2*Q3*Q4*Q8*Q11-Q2^2*Q4*Q11^2- Q1*Q5^2*Q6*Q12+Q1*Q2*Q5*Q11*Q12-Q2^2*Q3* Q11*Q13+Q2*Q3*Q4*Q6*Q14-Q1*Q5^2*Q6*Q14+Q1*Q2*Q5*Q11*Q15,
           Q5*Q8*Q9*Q13-Q2*Q11*Q13^2-Q4*Q7*Q11*Q14+Q4*Q6*Q13*Q14-Q5*Q7*Q9*Q15-2*Q4*Q6* Q13*Q15+2*Q3*Q7*Q13*Q15-Q1*Q12*Q13*Q15+2*Q1*Q13*Q15^2,
           Q4*Q5*Q6*Q9-Q3*Q4*Q8*Q9-Q4^2* Q6*Q11+Q3*Q4*Q7*Q11+Q1*Q5*Q9*Q12-2*Q1* Q5*Q11*Q13+Q1*Q3*Q12*Q13-Q1*Q5*Q9*Q15+ Q1*Q4*Q11*Q15,
           -Q5*Q11*Q12*Q13-Q4*Q11*Q12*Q14-Q5*Q11*Q13* Q14+Q4*Q11*Q12*Q15-Q5*Q11*Q13*Q15+Q4*Q11* Q14*Q15+Q5*Q9*Q15^2+Q3*Q13*Q15^2,
           -Q3*Q5*Q8^2*Q9+Q3*Q5*Q7*Q8*Q11-Q5^2*Q6^2*Q13+2*Q3*Q5*Q6*Q7*Q14-Q3^2*Q7*Q8*Q14+ Q2*Q3*Q8*Q9*Q14-Q1*Q5*Q8*Q11*Q14+Q1*Q5*Q6*Q12*Q14-Q1*Q3*Q8*Q14^2+Q2*Q3*Q8* Q9*Q15-Q2*Q3*Q7*Q11*Q15+Q1*Q5*Q8*Q11*Q15-Q1*Q3*Q8*Q14*Q15,
           -Q5*Q6*Q7*Q9+Q4*Q6*Q7*Q11-Q4*Q6^2*Q13+Q2*Q6*Q9*Q13-2*Q1*Q8*Q11*Q13+Q1*Q8* Q9*Q14+Q1*Q7*Q11*Q14-Q1*Q7*Q11*Q15+Q1*Q6*Q13*Q15,
           Q5*Q8*Q9-Q4*Q8*Q11+Q4*Q6*Q12-Q3*Q8*Q13+Q1*Q12*Q14-Q1*Q12*Q15,
           Q3*Q8^2*Q11*Q12-Q2*Q8*Q11^2*Q12-Q5*Q6^2* Q12^2-Q3*Q8^2*Q11*Q14+Q2*Q8*Q11^2*Q14+Q5*Q6^2* Q12*Q14+Q3*Q6*Q8*Q12*Q14-Q2*Q6*Q11*Q12* Q14-Q3*Q8^2*Q11*Q15+Q2*Q8*Q11^2*Q15+Q3*Q6*Q8*Q12*Q15-Q3*Q6*Q8*Q14*Q15,
           Q5*Q7*Q11-Q4*Q8*Q11-Q2*Q11*Q13+Q4*Q6* Q14+Q1*Q12*Q14-Q1*Q14*Q15,
           Q5*Q8*Q9*Q11-Q5*Q7*Q11^2+Q5*Q6*Q9*Q12+Q3*Q7*Q11*Q12-Q2*Q9*Q11*Q12-Q1*Q11*Q12^2-Q5*Q6*Q11*Q13+Q3*Q6*Q13*Q14-Q3*Q8* Q9*Q15+Q3*Q7*Q11*Q15-Q1*Q11*Q12*Q15+Q1*Q11*Q14*Q15,
           -2*Q5^2*Q8^2*Q11^2+Q5^2*Q6*Q8*Q11*Q12-Q3*Q5*Q8^2*Q11*Q12-2*Q2*Q5*Q8*Q11^2*Q12+ Q2*Q5*Q6*Q11*Q12^2+Q5^2*Q6*Q8*Q11*Q14+Q2* Q5*Q8*Q11^2*Q14-Q5^2*Q6^2*Q12*Q14+Q2^2*Q11^2*Q12*Q14+Q3*Q5*Q6*Q8*Q14^2+Q3^2*Q8^2*Q14^2- Q2*Q3*Q6*Q12*Q14^2+2*Q3*Q5*Q8^2*Q11*Q15+2*Q2* Q5*Q8*Q11^2*Q15+Q2^2*Q11^2*Q12*Q15-Q3*Q5*Q6*Q8*Q14*Q15-Q2*Q5*Q6*Q11*Q14*Q15-3*Q2* Q3*Q8*Q11*Q14*Q15-Q2^2*Q11^2*Q14*Q15+Q2*Q3*Q6*Q14^2*Q15,
           -Q2*Q3*Q7*Q12*Q13+Q1*Q5*Q8* Q12*Q13+Q1*Q2*Q12^2*Q13+Q2*Q3*Q8*Q13^2-Q3*Q5*Q7^2*Q14+Q3*Q4*Q7*Q8*Q14+Q1*Q5* Q7*Q12*Q14-Q1*Q4*Q8*Q12*Q14-Q1*Q5*Q8*Q13*Q14-Q1*Q5*Q8*Q13*Q15-Q1*Q5*Q7*Q14*Q15+Q1*Q4*Q8*Q14*Q15,
           -Q8*Q11*Q12*Q13+Q8*Q9*Q12*Q14-Q8*Q11*Q13* Q14+Q7*Q11*Q14^2+Q6*Q13*Q14^2-Q8*Q9*Q12* Q15-Q8*Q11*Q13*Q15+Q8*Q9*Q14*Q15,
           Q4^2*Q11*Q12^2-Q5^2*Q9*Q12*Q13-2*Q4*Q5*Q11* Q12*Q13+2*Q5^2*Q11*Q13^2-Q3*Q5*Q12*Q13^2+Q3*Q4*Q12*Q13*Q14+Q5^2*Q9*Q13*Q15-Q4*Q5*Q11* Q13*Q15+Q3*Q4*Q12*Q13*Q15-Q3*Q4*Q13*Q14*Q15,
           -Q5^2*Q9^2*Q12+Q4^2*Q11^2*Q12-Q3*Q4*Q11*Q12* Q13+Q3*Q4*Q9*Q12*Q14+Q3*Q4*Q11*Q13*Q14- Q3^2*Q13^2*Q14+Q5^2*Q9^2*Q15-Q4^2*Q11^2*Q15+Q3* Q5*Q9*Q13*Q15+Q3^2*Q13^2*Q15-Q3*Q4*Q9*Q14*Q15-Q3*Q4*Q9*Q15^2,
           Q1*Q5*Q8*Q11*Q12+Q5^2*Q6^2*Q13+Q3*Q5* Q6*Q8*Q13+Q3^2*Q8^2*Q13-Q2^2*Q11^2*Q13-Q2*Q3* Q6*Q12*Q13-Q1*Q5*Q6*Q12*Q14-Q2*Q3*Q6* Q13*Q14+Q1*Q5*Q8*Q11*Q15+Q1*Q2*Q11*Q12* Q15-Q2*Q3*Q6*Q13*Q15-Q1*Q3*Q8*Q14*Q15,
           -Q4*Q5*Q8^2*Q9-Q5^2*Q7^2*Q11+Q2*Q4*Q8*Q9*Q12+Q4*Q5*Q6*Q8*Q13-Q2*Q4*Q6*Q12* Q13+Q1*Q5*Q8*Q12*Q13+2*Q4*Q5*Q6*Q7*Q14-Q4^2*Q6*Q8*Q14+Q2*Q4*Q8*Q9*Q14-Q1*Q4* Q8*Q12*Q14-Q1*Q5*Q8*Q13*Q14-Q1*Q4*Q8*Q14^2+Q1*Q5*Q7*Q14*Q15,
           -Q5*Q6*Q9+Q3*Q8*Q9-Q3*Q7*Q11+Q3*Q6* Q13,
           -Q5*Q8*Q9^2*Q12+Q5*Q7*Q11^2*Q13-Q5* Q6*Q11*Q13^2+Q5*Q7*Q9*Q11*Q14-Q3*Q7*Q11*Q13*Q14+Q3*Q6*Q13^2*Q14+Q5*Q8*Q9^2*Q15+Q5* Q7*Q9*Q11*Q15-Q5*Q6*Q9*Q13*Q15-Q3*Q7*Q11*Q13*Q15+Q3*Q6*Q13^2*Q15-Q3*Q7*Q9*Q14* Q15,
           -Q4*Q8^2*Q9*Q11-Q5*Q7^2*Q11^2+Q4*Q6*Q8*Q9*Q12+Q2*Q8*Q9*Q11*Q13-Q2*Q6*Q9* Q12*Q13+Q1*Q8*Q11*Q12*Q13+Q4*Q6*Q8*Q9*Q15-Q2*Q8*Q9^2*Q15+2*Q2*Q7*Q9*Q11*Q15-Q1* Q8*Q9*Q12*Q15-Q1*Q8*Q11*Q13*Q15+Q1*Q7*Q11*Q14*Q15-Q1*Q8*Q9*Q15^2,
           Q4*Q5^2*Q7*Q9+Q4^2*Q5*Q8*Q9-Q4^2*Q5* Q7*Q11-Q4^3*Q8*Q11+Q3*Q4^2*Q8*Q13-Q1*Q4*Q5*Q12*Q13-Q1*Q5^2*Q13^2+Q1*Q4^2*Q12*Q14,
           Q2*Q5^2*Q6*Q7-Q1*Q5^2*Q8^2+Q2^2*Q5*Q7* Q11+Q2^2*Q4*Q8*Q11-Q2^2*Q5*Q6*Q13-Q2^3*Q11* Q13-Q1*Q2*Q5*Q8*Q15+Q1*Q2^2*Q12*Q15,
           Q5^2*Q6*Q7*Q8-2*Q4*Q5*Q6*Q8^2+Q3*Q4* Q8^3+Q2*Q5*Q8^2*Q9+Q2*Q4*Q8^2*Q11-Q2*Q4* Q6*Q8*Q12+2*Q2^2*Q8*Q9*Q12-Q2^2*Q6*Q12*Q13-Q2*Q5*Q6*Q7*Q14+2*Q2*Q4*Q6*Q8*Q14-2* Q1*Q5*Q8^2*Q14-2*Q2^2*Q8*Q9*Q14-Q1*Q2*Q8* Q12*Q14+2*Q1*Q2*Q8*Q14^2+Q1*Q2*Q8*Q12* Q15-Q1*Q2*Q8*Q14*Q15,
           Q3*Q6*Q13-Q1*Q11*Q15,
           Q5*Q6*Q9-Q1*Q11*Q15,
           Q2*Q6*Q9-Q1*Q8*Q11+Q1*Q6*Q12-Q1*Q6* Q14,
           -Q3*Q4*Q7*Q8*Q12+Q1*Q4*Q8*Q12^2+Q3*Q4*Q8^2*Q13+Q1*Q5*Q8*Q12*Q13-Q1*Q5* Q8*Q13*Q14-Q3*Q5*Q7^2*Q15+Q1*Q5*Q7*Q12*Q15+Q2*Q3*Q7*Q13*Q15-Q1*Q5*Q8*Q13*Q15- Q1*Q2*Q12*Q13*Q15-Q1*Q5*Q7*Q14*Q15+Q1*Q2*Q13*Q14*Q15,
           -Q4*Q5*Q7*Q8*Q9-Q4^2*Q8^2* Q9-Q4^2*Q7*Q8*Q11+Q2*Q5*Q7*Q9*Q13+2*Q2* Q4*Q8*Q9*Q13+Q1*Q4*Q8*Q12*Q13-Q2^2*Q9* Q13^2-Q1*Q2*Q12*Q13^2+Q1*Q4*Q8*Q13*Q14+ Q2*Q4*Q7*Q9*Q15-Q1*Q2*Q13^2*Q15,
           -Q4*Q8*Q11+Q2*Q11*Q13+Q3*Q7*Q15-Q1*Q12* Q15,
           -Q3*Q4*Q8*Q9*Q11+Q1*Q5*Q9*Q11*Q12-Q1*Q4*Q11^2*Q12+Q3*Q5*Q6*Q9*Q13-Q1* Q5*Q11^2*Q13-Q3^2*Q6*Q13^2+Q1*Q4*Q11^2*Q14-Q1* Q3*Q9*Q12*Q14+2*Q1*Q3*Q11*Q13*Q14+Q1*Q3* Q9*Q12*Q15-Q1*Q3*Q9*Q14*Q15,
           Q4*Q8*Q9*Q12-Q5*Q8*Q9*Q13+Q4*Q8*Q11* Q13-Q4*Q6*Q12*Q13+Q3*Q7*Q12*Q13-Q1*Q12^2* Q13-Q3*Q8*Q13^2-Q4*Q7*Q11*Q14+Q3*Q7*Q13*Q14-Q1*Q12*Q13*Q14+Q5*Q7*Q9*Q15+Q1*Q13*Q14*Q15,
           -Q8^2*Q9*Q11*Q12+Q7*Q8*Q11^2*Q12- Q6^2*Q12^2*Q13+Q8^2*Q9*Q11*Q14-Q7*Q8*Q11^2*Q14-Q6*Q8*Q9*Q12*Q14+Q6*Q7*Q11*Q12*Q14+Q6^2*Q12*Q13*Q14+Q8^2*Q9*Q11*Q15-Q7*Q8*Q11^2* Q15+Q6*Q7*Q11*Q12*Q15-Q6*Q7*Q11*Q14*Q15,
           Q5^2*Q7*Q8*Q12-Q4*Q5*Q8^2*Q12+Q5^2*Q7* Q8*Q14-Q4*Q5*Q8^2*Q14-Q2*Q4*Q8*Q12*Q14-Q5^2*Q7*Q8*Q15+Q4*Q5*Q8^2*Q15-Q2*Q5*Q7* Q12*Q15+Q2*Q4*Q8*Q12*Q15+Q2^2*Q12*Q13*Q15+ Q2*Q4*Q8*Q14*Q15-Q2^2*Q13*Q15^2,
           -Q5*Q8^2*Q9+2*Q4*Q6*Q8*Q12-2*Q3*Q7*Q8*Q12+2*Q1*Q8*Q12^2+Q2*Q8*Q11*Q13-Q2*Q6*Q12* Q13-Q5*Q6*Q7*Q14+Q3*Q7*Q8*Q14-Q1*Q8*Q12*Q15,
           Q3*Q8^2-Q2*Q6*Q12,
           -Q5*Q8*Q11*Q12+Q5*Q6*Q12^2+Q2*Q11*Q12^2-Q5*Q8*Q11*Q14+Q3*Q8*Q12*Q14-Q5*Q8*Q11* Q15+Q3*Q8*Q12*Q15-Q3*Q8*Q14*Q15,
           Q5*Q8*Q9+Q5*Q7*Q11-Q4*Q6*Q12+Q3*Q7* Q12-Q1*Q12^2-Q3*Q8*Q13-Q2*Q11*Q13+Q3*Q7* Q14-Q1*Q12*Q14+Q1*Q14*Q15,
           Q5*Q8*Q9-Q5*Q7*Q11-Q4*Q6*Q12+Q3*Q7* Q12-Q1*Q12^2+Q3*Q8*Q13-Q2*Q11*Q13+Q3*Q7* Q14-Q1*Q12*Q14+Q1*Q14*Q15,
           -Q5*Q8*Q9-Q5*Q7*Q11-Q4*Q6*Q12+Q3*Q7* Q12-Q1*Q12^2+Q3*Q8*Q13+Q2*Q11*Q13+Q3*Q7*Q14-Q1*Q12*Q14+Q1*Q14*Q15,
           -Q5*Q8*Q9+Q5*Q7*Q11-Q4*Q6*Q12+Q3*Q7* Q12-Q1*Q12^2-Q3*Q8*Q13+Q2*Q11*Q13+Q3*Q7*Q14-Q1*Q12*Q14+Q1*Q14*Q15,
           -Q4*Q9*Q12-Q5*Q9*Q13+Q4*Q11*Q13+Q3*Q13^ 2,
           Q4*Q9*Q12-Q5*Q9*Q13+Q4*Q11*Q13-Q3*Q13^2,
           -Q4*Q9*Q12+Q5*Q9*Q13+Q4*Q11*Q13-Q3* Q13^2,
           -Q4*Q9*Q12+Q5*Q9*Q13-Q4*Q11*Q13+Q3*Q13^2,
           Q4*Q9*Q12+Q5*Q9*Q13-Q4*Q11*Q13-Q3* Q13^2,
           -Q5^2*Q7*Q9^2+Q4*Q5*Q7*Q9*Q11-Q3*Q4*Q8*Q9*Q13-Q1*Q5*Q11*Q13^2+Q1*Q4*Q9* Q12*Q14+Q1*Q4*Q11*Q13*Q14-Q1*Q3*Q13^2*Q14-Q1*Q4*Q9*Q12*Q15+2*Q1*Q5*Q9*Q13*Q15+Q1* Q3*Q13^2*Q15-Q1*Q4*Q9*Q14*Q15,
           Q5*Q7*Q11-Q2*Q9*Q12,
           Q5*Q8*Q9-Q4*Q6*Q14,
           -Q5^2*Q7*Q11+Q4*Q5*Q8*Q11-Q2*Q5*Q9*Q12+Q2*Q4*Q11*Q12-Q2*Q5*Q11*Q13+Q4*Q5*Q6*Q14-Q3*Q4*Q8*Q14-Q2*Q5*Q9*Q14+Q1*Q5*Q12*Q14+Q1*Q5*Q14^2-Q1*Q5*Q12*Q15+Q2* Q3*Q13*Q15,
           -Q5^2*Q6^2*Q9+Q4*Q5*Q6^2*Q11+2*Q1* Q5*Q6*Q11*Q12+Q1*Q2*Q11^2*Q12-Q3*Q5*Q6^2*Q13-Q2*Q3*Q6*Q11*Q13+Q3*Q4*Q6^2*Q14- Q1*Q5*Q6*Q11*Q15+Q1*Q3*Q8*Q11*Q15-Q1*Q2*Q11^2*Q15-Q1*Q3*Q6*Q12*Q15,
           Q4*Q5^2*Q7*Q9-Q4^2*Q5*Q8*Q9+Q4^2*Q5* Q7*Q11-Q4^3*Q8*Q11+Q2*Q4^2*Q11*Q13-Q1*Q5^2* Q13^2+Q1*Q4^2*Q12*Q14-Q1*Q4*Q5*Q13*Q14,
           Q2*Q3*Q4*Q8^2-Q1*Q5^2*Q8^2-Q2^2*Q4*Q8* Q11+Q2^2*Q5*Q6*Q13+Q2^2*Q3*Q8*Q13-Q2^3*Q11* Q13-Q1*Q2*Q5*Q8*Q15+Q1*Q2^2*Q12*Q15,
           Q3*Q5^2*Q7^2-Q4^2*Q5*Q6*Q8-Q3*Q4*Q5* Q7*Q8+Q3*Q4^2*Q8^2+Q2*Q4^2*Q8*Q11-Q2*Q3* Q4*Q7*Q12-Q1*Q5^2*Q7*Q12+Q1*Q4*Q5*Q8* Q12+Q1*Q2*Q4*Q12^2-Q2^2*Q3*Q13^2+Q2*Q3* Q4*Q7*Q14-Q1*Q5^2*Q7*Q14-Q1*Q4*Q5*Q8*Q14+Q1*Q5^2*Q7*Q15-Q1*Q4*Q5*Q8*Q15-Q1* Q2*Q4*Q12*Q15+2*Q1*Q2*Q5*Q13*Q15,
           Q5^2*Q7*Q12*Q13-Q2*Q5*Q12*Q13^2-Q4*Q5* Q7*Q12*Q14+Q4^2*Q8*Q12*Q14-Q5^2*Q7*Q13*Q14+Q2*Q4*Q12*Q13*Q14+Q2*Q5*Q13^2*Q14-Q4^2*Q8* Q14^2+Q5^2*Q7*Q13*Q15-Q2*Q4*Q12*Q13*Q15-Q2* Q5*Q13^2*Q15+Q2*Q4*Q13*Q14*Q15,
           -Q1*Q12^2+Q2*Q11*Q13,
           Q4*Q8*Q12-Q5*Q7*Q15,
           -Q4*Q8*Q12+Q5*Q7*Q15,
           Q3*Q5*Q6*Q8^2*Q9+Q3^2*Q8^3*Q9-2*Q2*Q3* Q8^2*Q9*Q11-Q2*Q4*Q6*Q8*Q11^2+Q2^2*Q8*Q9* Q11^2-Q2*Q5*Q6^2*Q9*Q12-Q2*Q3*Q6*Q8* Q9*Q12+Q2*Q4*Q6^2*Q11*Q12+Q1*Q5*Q6*Q8*Q11*Q12+Q1*Q2*Q8*Q11^2*Q12-Q1*Q3*Q8^2*Q11* Q14+2*Q2*Q3*Q6*Q8*Q9*Q15+Q2*Q4*Q6^2*Q11* Q15-Q1*Q5*Q6*Q8*Q11*Q15-Q1*Q3*Q8^2*Q11* Q15-Q2^2*Q6*Q9*Q11*Q15+Q1*Q2*Q8*Q11^2* Q15-Q1*Q2*Q6*Q11*Q12*Q15+Q1*Q3*Q6*Q8*Q14*Q15-Q1*Q2*Q6*Q11*Q15^2,
           Q4*Q6*Q12-Q1*Q14^2,
           Q2*Q5*Q9-Q2*Q4*Q11+Q2*Q3*Q13-Q1*Q5* Q14,
           -Q5*Q13+Q4*Q15,
           Q2*Q5*Q9+Q2*Q4*Q11-Q2*Q3*Q13-Q1*Q5* Q14,
           Q5*Q13-Q4*Q15,
           Q2*Q4*Q8*Q11^2-Q1*Q5*Q8*Q11*Q12-Q4*Q5*Q6^2*Q14+Q3*Q4*Q6*Q8*Q14-Q1*Q5*Q8* Q11*Q14-Q1*Q5*Q6*Q12*Q14+Q1*Q3*Q8*Q12*Q14-Q2*Q4*Q6*Q11*Q15+Q1*Q5*Q8*Q11*Q15+Q1*Q5*Q6*Q14*Q15-Q1*Q3*Q8*Q14*Q15+Q1*Q2*Q11*Q15^2,
           Q4*Q8^2*Q9*Q12-Q4*Q8^2*Q11*Q13-Q2*Q8*Q9*Q12*Q13+Q2*Q7*Q11*Q12*Q13+Q2* Q8*Q11*Q13^2-Q5*Q7^2*Q11*Q14+Q4*Q8^2*Q9*Q15+ Q5*Q7^2*Q11*Q15-Q4*Q7*Q8*Q11*Q15-Q2*Q7*Q9*Q12*Q15-Q2*Q8*Q9*Q13*Q15+Q2*Q7*Q11*Q13*Q15,
           Q4*Q8*Q9-Q1*Q12*Q13,
           -Q4*Q8*Q9+Q1*Q12*Q13,
           Q4*Q5*Q8*Q9*Q12+Q4^2*Q8*Q11*Q12+Q5^2* Q8*Q9*Q13-Q3*Q4*Q8*Q12*Q13-Q2*Q5*Q9*Q12*Q13-Q3*Q5*Q8*Q13^2+Q2*Q3*Q12*Q13^2-Q4^2* Q8*Q11*Q14+Q4*Q5*Q8*Q9*Q15-Q2*Q4*Q9*Q12*Q15-Q2*Q5*Q9*Q13*Q15+Q2*Q3*Q13^2*Q15,
           -2*Q3*Q5*Q8^2*Q11-2*Q2*Q5*Q8*Q11^2-Q2^2*Q11^2* Q12-Q3^2*Q8^2*Q14+Q2*Q3*Q6*Q12*Q14+Q3* Q5*Q6*Q8*Q15+Q3^2*Q8^2*Q15+Q2*Q5*Q6*Q11*Q15+3*Q2*Q3*Q8*Q11*Q15+Q2^2*Q11^2*Q15-Q2*Q3* Q6*Q12*Q15-Q2*Q3*Q6*Q14*Q15,
           -Q4^2*Q8^2*Q12+Q2*Q4*Q8*Q12*Q13+Q5^2*Q7^2*Q14+Q4*Q5*Q7*Q8*Q14+Q4^2*Q8^2*Q14-Q2*Q4* Q7*Q12*Q14-Q2^2*Q13^2*Q14-Q2*Q4*Q7*Q14^2-Q5^2*Q7^2*Q15+Q2*Q4*Q7*Q12*Q15-Q2*Q4*Q8* Q13*Q15+Q2^2*Q13^2*Q15,
           Q4*Q6-Q3*Q7+Q2*Q9+Q1*Q12-Q1*Q14-Q1* Q15,
           -Q4*Q6+Q3*Q7+Q2*Q9-Q1*Q12-Q1*Q14+ Q1*Q15,
           Q4*Q6+Q3*Q7-Q2*Q9-Q1*Q12+Q1* Q14-Q1*Q15,
           -Q5*Q7*Q9*Q11+Q4*Q7*Q11^2+ Q3*Q8*Q9*Q13+Q3*Q6*Q13^2+Q1*Q9*Q12*Q14-2*Q1*Q11*Q13*Q14+Q4*Q6*Q9*Q15-Q2*Q9^2*Q15- Q1*Q9*Q12*Q15+Q1*Q9*Q14*Q15-Q1*Q9*Q15^2,
           Q5*Q7*Q9*Q11+Q4*Q7*Q11^2-Q3*Q8*Q9*Q13+Q3*Q6*Q13^2+Q1*Q9*Q12*Q14-2*Q1*Q11*Q13* Q14+Q4*Q6*Q9*Q15-Q2*Q9^2*Q15-Q1*Q9*Q12*Q15+Q1*Q9*Q14*Q15-Q1*Q9*Q15^2,
           Q4*Q8*Q11-Q2*Q11*Q13+Q3*Q7*Q14-Q1*Q12* Q14,
           Q2*Q11*Q12-Q3*Q8*Q15,
           -Q4*Q5*Q6^2+Q3*Q5*Q6*Q7+Q2*Q4*Q6*Q11-2*Q1*Q5*Q8*Q11+Q1*Q3*Q8*Q12+Q1*Q2* Q11*Q12-Q2*Q3*Q6*Q13+Q1*Q5*Q6*Q15-Q1*Q2*Q11*Q15,
           Q2*Q11^2-Q3*Q6*Q14,
           -Q5*Q8^2*Q9^2-Q4*Q7*Q8*Q11^2+Q3*Q7*Q8*Q11*Q13+Q4*Q6*Q7*Q11*Q14-Q3*Q6*Q7*Q13* Q14+Q1*Q8*Q11*Q13*Q14+2*Q3*Q7*Q8*Q9*Q15+Q4*Q6*Q7*Q11*Q15-Q3*Q7^2*Q11*Q15+Q1*Q8* Q9*Q12*Q15-Q1*Q8*Q11*Q13*Q15-Q1*Q7*Q11*Q14*Q15-Q1*Q7*Q11*Q15^2,
           -Q4*Q8*Q9*Q12-Q5*Q7*Q11*Q13+Q4*Q8*Q11* Q13+Q4*Q6*Q12*Q13-Q3*Q7*Q12*Q13+Q1*Q12^2 *Q13-Q2*Q11*Q13^2+Q4*Q7*Q11*Q14-Q3*Q7*Q13*Q14+Q1*Q12*Q13*Q14+Q5*Q7*Q9*Q15-Q1*Q13*Q14*Q15,
           Q3*Q7^2*Q8*Q9*Q12-Q2*Q7*Q8*Q9^2*Q12-Q1*Q7*Q8*Q9*Q12^2-Q3*Q7*Q8^2*Q9* Q13+Q2*Q8^2*Q9^2*Q13+2*Q2*Q6*Q7*Q9*Q12*Q13+Q1*Q8^2*Q9*Q12*Q13-Q1*Q7*Q8*Q11*Q12* Q13-2*Q2*Q6*Q8*Q9*Q13^2+Q2*Q6*Q7*Q11*Q13^2-Q1*Q6*Q8*Q12*Q13^2+Q2*Q6^2*Q13^3+Q1*Q6*Q7*Q12*Q13*Q14-Q1*Q6*Q8*Q13^2*Q14+Q3*Q7^2* Q8*Q9*Q15-Q2*Q7^2*Q9*Q11*Q15-Q1*Q7*Q8*Q9*Q12*Q15-Q2*Q6*Q7*Q9*Q13*Q15+Q1*Q8^2* Q9*Q13*Q15+Q1*Q7*Q8*Q11*Q13*Q15,
           Q4^2*Q5*Q6*Q8-Q3*Q4^2*Q8^2-Q2*Q4*Q5* Q7*Q11-Q2*Q4^2*Q8*Q11+Q2*Q4^2*Q6*Q12+Q1*Q5^2*Q7*Q14+2*Q1*Q4*Q5*Q8*Q14-Q1*Q5^2*Q7*Q15-Q1*Q4*Q5*Q8*Q15+Q1*Q2*Q5*Q13*Q15-Q1*Q2*Q4*Q14*Q15,
           -Q4^2*Q8^2*Q12-2*Q5^2*Q7*Q8*Q13-2*Q4*Q5*Q8^2* Q13+Q5^2*Q7^2*Q14+3*Q4*Q5*Q7*Q8*Q14+Q4^2* Q8^2*Q14-Q2*Q4*Q7*Q12*Q14+Q2*Q5*Q7*Q13*Q14+Q2*Q4*Q8*Q13*Q14-Q5^2*Q7^2*Q15+Q2*Q4* Q7*Q12*Q15-Q2*Q4*Q7*Q14*Q15,
           -Q5*Q6*Q13+Q3*Q8*Q13-Q2*Q11*Q13+Q4*Q6* Q14+Q1*Q12*Q14-Q1*Q14*Q15,
           -Q8*Q9*Q12*Q13-2*Q8*Q11*Q13^2-Q7*Q11*Q13*Q14+Q8*Q9*Q13*Q15+Q7*Q11*Q13*Q15+Q6*Q13^2* Q15+Q7*Q9*Q15^2,
           Q3*Q4*Q6-Q3^2*Q7-Q1*Q5*Q11+Q1*Q3*Q12,
           -Q5*Q6*Q7*Q8*Q9-Q2*Q7^2*Q11^2-Q1*Q8^2* Q9*Q12+Q2*Q6*Q7*Q11*Q13-Q1*Q8^2*Q11*Q13+ Q1*Q6*Q8*Q12*Q13+Q1*Q6*Q7*Q12*Q14+Q1* Q8^2*Q9*Q15+2*Q1*Q7*Q8*Q11*Q15-Q1*Q6*Q7*Q12*Q15-Q1*Q6*Q7*Q14*Q15,
           Q5*Q7*Q11-Q1*Q15^2,
           -Q4^2*Q6*Q7*Q8*Q12+Q3*Q4*Q7^2*Q8*Q12-Q1*Q4*Q7*Q8*Q12^2+Q4^2*Q6*Q8^2*Q13-Q3*Q4* Q7*Q8^2*Q13+2*Q2*Q4*Q6*Q7*Q12*Q13-Q1*Q5* Q7*Q8*Q12*Q13+Q1*Q4*Q8^2*Q12*Q13+Q2*Q5*Q6*Q7*Q13^2-2*Q2*Q4*Q6*Q8*Q13^2-Q1*Q2* Q8*Q12*Q13^2+Q2^2*Q6*Q13^3-Q4*Q5*Q6*Q7^2*Q14+Q3*Q4*Q7^2*Q8*Q14-Q1*Q4*Q7*Q8*Q12* Q14-Q2*Q4*Q6*Q7*Q13*Q14+Q1*Q5*Q7*Q8*Q13*Q14+Q1*Q4*Q8^2*Q13*Q14+Q1*Q2*Q7*Q12* Q13*Q15-Q1*Q2*Q8*Q13^2*Q15,
           -Q5*Q8*Q13-Q2*Q12*Q13+Q5*Q7*Q15+Q2*Q13* Q15,
           -Q4*Q8*Q11+Q3*Q8*Q13+Q2*Q9*Q15-Q1*Q14*Q15,
           -Q5*Q8*Q12*Q13-Q5*Q8*Q13*Q14+ Q2*Q12*Q13*Q14+Q5*Q7*Q14^2+Q4*Q8*Q14^2-Q5*Q8*Q13*Q15-Q2*Q12*Q13*Q15+Q2*Q13*Q14*Q15,
           -Q2*Q5*Q9^2*Q12+Q2*Q4*Q9*Q11*Q12-Q1*Q5*Q11*Q12*Q13+Q2*Q3*Q11*Q13^2+Q1*Q5*Q9* Q12*Q14-Q1*Q4*Q11*Q12*Q14-Q2*Q3*Q9*Q13*Q14+Q1*Q5*Q11*Q13*Q14+Q1*Q3*Q13*Q14^2-Q1* Q5*Q9*Q12*Q15+Q1*Q4*Q11*Q12*Q15-Q1*Q5*Q11*Q13*Q15,
           -2*Q5*Q8*Q11^2-Q5*Q6*Q11*Q12+Q5* Q6*Q11*Q14+Q3*Q8*Q11*Q14+Q2*Q11^2*Q14+Q3*Q6*Q14^2-Q3*Q8*Q11*Q15,
           -Q5*Q6*Q7*Q9+Q3*Q7*Q8*Q9-Q2*Q8*Q9^ 2+Q1*Q7*Q11*Q12+Q2*Q6*Q9*Q13-2*Q1*Q8* Q11*Q13+Q1*Q6*Q12*Q13+Q1*Q8*Q9*Q14-Q1*Q6*Q13*Q14,
           2*Q5*Q8^2*Q11^2-2*Q5*Q6*Q8*Q11*Q12- Q3*Q8^2*Q11*Q12-Q2*Q8*Q11^2*Q12+Q5*Q6^2*Q12^2+Q3*Q6*Q8*Q12*Q14-Q5*Q6*Q8*Q11*Q15+ Q2*Q8*Q11^2*Q15+Q3*Q6*Q8*Q12*Q15-Q3*Q6*Q8*Q14*Q15,
           -Q4*Q8*Q9+Q4*Q7*Q11+Q4*Q6* Q13-Q1*Q13*Q15,
           Q4*Q8*Q9-Q4*Q7*Q11+Q4*Q6*Q13-Q1*Q13* Q15,
           -Q5^2*Q6*Q9-Q3*Q5*Q7*Q11-Q3*Q4*Q8*Q11+2*Q2*Q5*Q9*Q11-Q2*Q4*Q11^2+2*Q1*Q5* Q11*Q12+Q3*Q5*Q6*Q13-Q2*Q3*Q11*Q13-Q3^2*Q7*Q15+Q2*Q3*Q9*Q15+Q1*Q3*Q12*Q15-Q1* Q3*Q14*Q15,
           -Q5^2*Q7*Q11-Q4*Q5*Q8*Q11-Q4* Q5*Q6*Q12+Q2*Q4*Q11*Q12+Q2*Q5*Q11*Q13+ Q3*Q4*Q8*Q14-Q1*Q5*Q12*Q14-Q4*Q5*Q6* Q15+Q2*Q5*Q9*Q15+Q1*Q5*Q12*Q15-Q2*Q3* Q13*Q15+Q1*Q5*Q15^2,
           Q4*Q8^2*Q9-Q5*Q7^2*Q11+Q4*Q7*Q8*Q11+ Q2*Q7*Q9*Q12-2*Q2*Q8*Q9*Q13+Q2*Q6*Q13^2-Q4*Q6*Q7*Q14-Q1*Q7*Q12*Q14+2*Q1*Q8*Q13* Q14+Q1*Q7*Q12*Q15-2*Q1*Q8*Q13*Q15,
           Q4*Q8*Q9*Q11-Q4*Q7*Q11^2-Q5*Q6*Q9*Q13-Q3*Q8*Q9*Q13+2*Q3*Q7*Q11*Q13-Q2*Q9* Q11*Q13-Q3*Q6*Q13^2+2*Q1*Q11*Q13*Q14+Q3*Q7*Q9*Q15-Q2*Q9^2*Q15-Q1*Q9*Q12*Q15+Q1*Q9* Q14*Q15,
           -Q9*Q12^2+Q11*Q13*Q14,
           -Q5*Q7*Q11^2-Q5*Q6*Q9*Q12+Q2*Q9*Q11*Q12+Q3*Q8*Q11*Q13+2*Q4*Q6*Q11*Q14-2*Q2*Q9* Q11*Q14-Q3*Q6*Q13*Q14+2*Q1*Q11*Q14^2-Q1*Q11*Q14*Q15,
           Q2*Q4*Q6-Q1*Q5*Q8-Q2^2*Q9+Q1* Q2*Q14,
           -Q5*Q8*Q9-Q5*Q7*Q11+Q5*Q6*Q13+ Q1*Q15^2,
           Q12-Q14,
           -Q12+Q14,
           Q5*Q8*Q9*Q12-Q5*Q7*Q11*Q12-Q4*Q8*Q11* Q12+Q5*Q8*Q11*Q13+Q5*Q6*Q12*Q13-Q5*Q7* Q11*Q14-Q5*Q6*Q13*Q14+Q5*Q8*Q9*Q15-Q4* Q8*Q11*Q15+Q3*Q7*Q12*Q15-Q2*Q9*Q12*Q15- Q1*Q12^2*Q15+Q1*Q12*Q14*Q15+Q3*Q7*Q15^2-Q1*Q12*Q15^2+Q1*Q14*Q15^2,
           -Q2*Q8*Q9*Q12+Q5*Q6*Q8*Q13-Q3*Q8^2*Q13-Q2*Q8*Q11*Q13+Q2*Q6*Q12*Q13-Q5*Q6*Q7*Q14+Q4*Q6*Q8*Q14-Q2*Q8*Q9*Q14+Q1*Q8*Q12*Q14+Q1*Q8*Q14^2+Q2*Q7*Q11*Q15-Q1* Q8*Q12*Q15,
           Q5*Q7*Q11-Q4*Q6*Q12,
           Q2*Q4*Q9+Q1*Q4*Q12-Q1*Q5*Q13-Q1*Q4* Q14,
           -Q4*Q8*Q11^2+Q5*Q6*Q11*Q13+2*Q3*Q7*Q11*Q14-2*Q2*Q9*Q11*Q14-Q1*Q11*Q12*Q14-Q3* Q6*Q13*Q14+2*Q1*Q11*Q14^2-Q3*Q8*Q9*Q15+Q2*Q9*Q11*Q15,
           2*Q5^2*Q11*Q13^2+Q4*Q5*Q9*Q12*Q14-Q5^2*Q9*Q13*Q14-2*Q4*Q5*Q11*Q13*Q14-Q3* Q5*Q13^2*Q14+Q4^2*Q11*Q14^2-Q4*Q5*Q9*Q12*Q15- Q4*Q5*Q11*Q13*Q15+Q3*Q5*Q13^2*Q15+Q4*Q5* Q9*Q14*Q15,
           Q5*Q7*Q8*Q9-Q4*Q8^2*Q9-Q5*Q6*Q7*Q13+2*Q4*Q6*Q8*Q13-Q3*Q7*Q8* Q13-Q2*Q7*Q11*Q13+2*Q1*Q8*Q12*Q13-Q2*Q6*Q13^2+Q4*Q6*Q7*Q14-Q3*Q7^2*Q14+Q1*Q7*Q12* Q14-Q1*Q7*Q14*Q15,
           -Q5*Q8*Q9+Q4*Q8*Q11+Q2*Q9*Q12-Q3*Q8* Q13-Q1*Q12*Q14+Q1*Q12*Q15,
           -Q5*Q6*Q8*Q9^2-Q3*Q8^2*Q9^2+2*Q3*Q7*Q8*Q9*Q11-Q3*Q7^2*Q11^2+Q1*Q8*Q9*Q11*Q12-Q3* Q6*Q8*Q9*Q13+Q3*Q6*Q7*Q11*Q13+Q3*Q6* Q7*Q9*Q14-Q1*Q7*Q11^2*Q14+Q1*Q8*Q9*Q11* Q15-Q1*Q7*Q11^2*Q15,
           -Q5*Q8^2*Q9-Q5*Q7*Q8*Q11+Q4*Q6*Q8*Q12-Q3*Q7*Q8*Q12+Q1*Q8*Q12^2+Q5*Q6*Q8* Q13-Q2*Q6*Q12*Q13+Q5*Q6*Q7*Q14-Q3*Q7*Q8*Q14+Q1*Q8*Q12*Q14+Q2*Q7*Q11*Q15-Q1*Q8*Q14*Q15,
           -Q5^2*Q7*Q8*Q12*Q13-2*Q4*Q5*Q8^2* Q12*Q13+Q2*Q4*Q8*Q12^2*Q13-2*Q5^2*Q8^2*Q13^2+Q2*Q5*Q8*Q12*Q13^2+Q4^2*Q8^2*Q12*Q14+2*Q5^2* Q7*Q8*Q13*Q14+2*Q4*Q5*Q8^2*Q13*Q14+Q4^2*Q8^2* Q12*Q15+Q4*Q5*Q8^2*Q13*Q15+Q2*Q5*Q8*Q13^ 2*Q15-Q2^2*Q12*Q13^2*Q15-3*Q4*Q5*Q7*Q8*Q14* Q15-Q4^2*Q8^2*Q14*Q15-Q2*Q5*Q7*Q13*Q14*Q15- Q2*Q4*Q8*Q13*Q14*Q15+Q5^2*Q7^2*Q15^2-Q2*Q4*Q7*Q12*Q15^2+Q2*Q5*Q7*Q13*Q15^2+Q2*Q4* Q7*Q14*Q15^2,
           Q5*Q8*Q9*Q12-Q5*Q7*Q11*Q12+ Q4*Q6*Q12^2+Q5*Q8*Q11*Q13-Q4*Q8*Q11*Q14+ Q1*Q12^2*Q14-Q3*Q8*Q13*Q14+Q5*Q8*Q9*Q15- Q5*Q7*Q11*Q15-Q4*Q8*Q11*Q15+Q4*Q6*Q12* Q15-Q2*Q9*Q12*Q15-Q1*Q12^2*Q15+Q3*Q8*Q13* Q15+Q1*Q12*Q14*Q15-Q1*Q12*Q15^2,
           -Q2*Q5*Q6*Q7*Q11+Q2*Q3*Q7*Q8*Q11-2*Q1*Q5*Q8^2*Q11-Q2^2*Q7*Q11^2-Q1*Q2*Q8*Q11* Q12+Q2*Q3*Q6*Q7*Q14-Q1*Q3*Q8^2*Q14+Q1*Q5*Q6*Q8*Q15+Q1*Q3*Q8^2*Q15+3*Q1*Q2*Q8* Q11*Q15-Q1*Q2*Q6*Q12*Q15,
           Q4*Q8*Q9-Q1*Q13*Q15,
           Q5^2*Q9^2*Q12+Q4*Q5*Q9*Q11*Q12+Q4^2*Q11^2* Q12-Q3*Q4*Q9*Q12^2-Q3^2*Q12*Q13^2-Q4^2*Q11^2*Q14-Q3*Q5*Q9*Q13*Q14+Q3^2*Q13^2*Q14-Q5^2*Q9^2* Q15-Q3*Q4*Q9*Q12*Q15+Q3*Q5*Q9*Q13*Q15+ Q3*Q4*Q9*Q14*Q15,
           Q5^2*Q6*Q7+Q3*Q4*Q8^2+Q2*Q4*Q8*Q11+ Q2*Q3*Q7*Q12-Q2^2*Q9*Q12-Q1*Q2*Q12^2-Q2*Q5*Q6*Q13-2*Q1*Q5*Q8*Q14+Q1*Q2*Q12*Q14- Q1*Q2*Q12*Q15+Q1*Q2*Q14*Q15,
           Q5*Q8*Q14-Q2*Q15^2,
           Q5^2*Q6*Q7+Q3*Q4*Q8^2-Q2*Q4*Q8*Q11+ Q2*Q3*Q7*Q12-Q2^2*Q9*Q12-Q1*Q2*Q12^2+Q2*Q5*Q6*Q13-2*Q1*Q5*Q8*Q14+Q1*Q2*Q12*Q14- Q1*Q2*Q12*Q15+Q1*Q2*Q14*Q15,
           Q4*Q11*Q12-Q3*Q13*Q14,
           -Q4*Q11^2*Q12*Q13+Q3*Q11*Q12*Q13^2-Q4*Q9*Q11*Q12*Q14-Q4*Q11^2*Q13*Q14+Q3*Q11*Q13^2*Q14+ Q5*Q9^2*Q12*Q15+Q4*Q9*Q11*Q12*Q15+Q4*Q11^2* Q13*Q15-Q3*Q9*Q12*Q13*Q15-Q3*Q11*Q13^2*Q15+ Q4*Q9*Q11*Q14*Q15-Q5*Q9^2*Q15^2,
           -Q5*Q6*Q14+Q3*Q8*Q15,
           Q5*Q7*Q11-Q4*Q8*Q11+Q2*Q9*Q15-Q1*Q14* Q15,
           Q4*Q7*Q9-Q1*Q13^2,
           -2*Q7*Q8*Q11^2*Q13-2*Q6*Q8*Q11*Q13^2-Q6^2*Q12* Q13^2+Q7*Q8*Q9*Q11*Q14+Q7^2*Q11^2*Q14-Q6*Q7*Q9*Q12*Q14+Q6*Q8*Q9*Q13*Q14+3*Q6*Q7* Q11*Q13*Q14+Q6^2*Q13^2*Q14-Q7^2*Q11^2*Q15+Q6*Q7* Q9*Q12*Q15-Q6*Q7*Q9*Q14*Q15,
           Q5*Q7*Q14-Q2*Q13*Q15,
           -Q5*Q8*Q9*Q11-Q5*Q7*Q11^2+Q5*Q6*Q9*Q12-Q2*Q9*Q11*Q12+Q5*Q6*Q11*Q13+Q4*Q6*Q11*Q14-Q2*Q9*Q11*Q14+Q1*Q11*Q12*Q14-Q3*Q6*Q13*Q14+Q1*Q11*Q14^2+Q3*Q8*Q9*Q15-Q1* Q11*Q12*Q15,
           -Q5*Q8*Q11+Q5*Q6*Q14+Q3*Q8*Q14-Q3*Q8*Q15,
           Q5^2*Q6*Q7-2*Q3*Q5*Q7*Q8+Q3*Q4*Q8^2-Q2*Q4*Q6*Q12+2*Q1*Q5*Q8*Q12+ Q2*Q5*Q6*Q13-Q2^2*Q11*Q13+Q2*Q3*Q7*Q14- Q1*Q2*Q12*Q14-2*Q1*Q5*Q8*Q15+Q1*Q2*Q14*Q15,
           -Q2*Q12*Q13+Q4*Q8*Q14,
           Q2*Q8*Q9-Q2*Q7*Q11+Q2*Q6*Q13-Q1*Q8* Q14,
           Q3*Q5^2*Q7^2*Q8-2*Q3*Q4*Q5*Q7*Q8^2+ Q3*Q4^2*Q8^3-Q2*Q5^2*Q7*Q8*Q9-Q1*Q4*Q5* Q8^2*Q12+Q2*Q3*Q4*Q8^2*Q13-Q2*Q3*Q5* Q7^2*Q14+2*Q2*Q3*Q4*Q7*Q8*Q14+Q1*Q5^2*Q7*Q8*Q14-Q1*Q4*Q5*Q8^2*Q14+Q2^2*Q5*Q7*Q9* Q14+Q1*Q2*Q4*Q8*Q12*Q14-Q1*Q2*Q5*Q8* Q13*Q14-Q1*Q2*Q5*Q7*Q14^2-Q2*Q3*Q4*Q7* Q8*Q15+Q1*Q5^2*Q7*Q8*Q15+Q2^2*Q5*Q7* Q9*Q15-Q2^2*Q3*Q7*Q13*Q15+Q1*Q2*Q5*Q8*Q13*Q15-Q1*Q2*Q5*Q7*Q14*Q15,
           Q2*Q8*Q9+Q2*Q7*Q11-Q2*Q6*Q13-Q1*Q8* Q14,
           Q5*Q7*Q9*Q11-Q4*Q7*Q11^2+Q4*Q6*Q9*Q12-Q2*Q9^2*Q12-Q5*Q6*Q9*Q13-Q3*Q8* Q9*Q13+2*Q4*Q6*Q11*Q13-Q2*Q9*Q11*Q13-Q3*Q6*Q13^2+Q1*Q9*Q12*Q14+2*Q1*Q11*Q13*Q14-Q1* Q9*Q12*Q15,
           Q5*Q7*Q11-Q4*Q8*Q11+Q1*Q12^2- Q2*Q11*Q13,
           Q2*Q5^2*Q9^2-Q4^2*Q5*Q6*Q11+ Q3*Q4^2*Q8*Q11-Q2*Q4*Q5*Q9*Q11+Q2*Q4^2*Q11^2+Q2*Q3*Q4*Q9*Q12-Q1*Q5^2*Q9*Q12-Q1* Q4*Q5*Q11*Q12-Q2*Q3^2*Q13^2-Q2*Q3*Q4*Q9* Q14-Q1*Q5^2*Q9*Q14+Q1*Q4*Q5*Q11*Q14+Q1*Q3*Q4*Q14^2+Q1*Q5^2*Q9*Q15-Q1*Q4*Q5* Q11*Q15+2*Q1*Q3*Q5*Q13*Q15-Q1*Q3*Q4*Q14*Q15,
           Q4*Q5*Q9*Q11*Q12^2+Q4^2*Q11^2*Q12^2+Q5^2*Q9*Q11*Q12*Q13-2*Q5^2*Q11^2*Q13^2+Q3*Q5*Q11*Q12* Q13^2-Q4*Q5*Q9*Q11*Q12*Q14+Q3*Q4*Q9*Q12^ 2*Q14+2*Q4*Q5*Q11^2*Q13*Q14-Q3*Q5*Q9*Q12* Q13*Q14-3*Q3*Q4*Q11*Q12*Q13*Q14+2*Q3*Q5*Q11*Q13^2*Q14-Q3^2*Q12*Q13^2*Q14-Q5^2*Q9^2*Q12*Q15- Q3*Q4*Q9*Q12^2*Q15+Q5^2*Q9*Q11*Q13*Q15-Q4*Q5*Q11^2*Q13*Q15-2*Q3*Q5*Q11*Q13^2*Q15+Q3^2* Q12*Q13^2*Q15+Q3^2*Q13^2*Q14*Q15+Q3*Q5*Q9*Q13*Q15^2,
           -2*Q5^2*Q6*Q8*Q11-2*Q2*Q5*Q8*Q11^2+Q5^2* Q6^2*Q12+Q3*Q5*Q6*Q8*Q12+3*Q2*Q5*Q6*Q11* Q12+Q2*Q3*Q8*Q11*Q12+Q2^2*Q11^2*Q12-Q5^2* Q6^2*Q14-Q2*Q3*Q6*Q12*Q14-Q2^2*Q11^2*Q15-Q2* Q3*Q6*Q12*Q15+Q2*Q3*Q6*Q14*Q15,
           Q2*Q4*Q11^2*Q13-Q1*Q5*Q11*Q12*Q13-Q2*Q4*Q9*Q11*Q14+Q1*Q5*Q11*Q13*Q14+Q1*Q4*Q11*Q14^2-Q2*Q5*Q9^2*Q15-Q1*Q5*Q9*Q12*Q15+ Q2*Q3*Q9*Q13*Q15-Q1*Q5*Q11*Q13*Q15+Q1*Q3*Q12*Q13*Q15+Q1*Q5*Q9*Q14*Q15-Q1*Q3* Q13*Q14*Q15,
           -Q5*Q8*Q9*Q11-Q4*Q8*Q11^2+Q5*Q6*Q9*Q12+Q3*Q8*Q11*Q13+Q3*Q7*Q11*Q14- Q2*Q9*Q11*Q14-Q3*Q6*Q13*Q14+Q1*Q11*Q14^2+Q3*Q8*Q9*Q15-Q2*Q9*Q11*Q15-Q1*Q11*Q12* Q15+Q1*Q11*Q14*Q15,
           Q5*Q7*Q8*Q9^2+Q4*Q8^2*Q9^2-Q4*Q7^2*Q11^2+ Q4*Q6*Q7*Q9*Q12-Q4*Q6*Q8*Q9*Q13-Q2*Q8*Q9^2*Q13-Q1*Q8*Q9*Q12*Q13+Q4*Q6^2*Q13^2 -Q1*Q6*Q12*Q13^2-Q1*Q8*Q9*Q13*Q14+2*Q1*Q7*Q11*Q13*Q14+Q1*Q6*Q13^2*Q14-Q4*Q6*Q7* Q9*Q15+Q1*Q8*Q9*Q13*Q15-Q1*Q6*Q13^2*Q15-Q1*Q7*Q9*Q14*Q15+Q1*Q7*Q9*Q15^2,
           -Q3*Q4*Q8*Q9-Q4^2*Q6*Q11+Q2*Q4*Q9*Q11+Q3*Q4*Q6*Q13-2*Q1*Q5*Q11*Q13+Q1*Q5* Q9*Q14+Q1*Q3*Q13*Q14+Q1*Q4*Q11*Q15-Q1*Q3*Q13*Q15,
           Q2*Q9*Q12-Q1*Q15^2,
           -Q4*Q6*Q8^2*Q9-Q4*Q6*Q7*Q8*Q11-Q1*Q8^2*Q9*Q12-Q1*Q7*Q8*Q11*Q12+Q4*Q6^2*Q8* Q13+2*Q1*Q8^2*Q11*Q13+Q4*Q6^2*Q7*Q14+Q1*Q6* Q7*Q12*Q14+Q1*Q8^2*Q9*Q15-Q1*Q6*Q8*Q13* Q15-Q1*Q6*Q7*Q14*Q15,
           -Q4^2*Q8*Q9*Q12+Q5^2*Q7*Q9*Q13+2*Q4^2*Q6*Q12*Q13-Q2*Q4*Q9*Q12*Q13+Q4*Q5*Q6*Q13^2-2* Q2*Q5*Q9*Q13^2+Q2*Q4*Q11*Q13^2+Q2*Q3*Q13^3+Q1*Q4*Q12*Q13*Q14-Q4*Q5*Q7*Q9*Q15-2* Q4^2*Q6*Q13*Q15+2*Q2*Q4*Q9*Q13*Q15-Q1*Q4*Q12*Q13*Q15-2*Q1*Q5*Q13^2*Q15-Q1*Q4*Q13*Q14* Q15+2*Q1*Q4*Q13*Q15^2,
           -Q5*Q7^2*Q9*Q11-Q4*Q7^2*Q11^2-Q5*Q6*Q7*Q9*Q13+Q3*Q7^2*Q11*Q13+Q1*Q8*Q9*Q12*Q13- Q1*Q7*Q11*Q12*Q13-Q1*Q6*Q12*Q13^2-Q1*Q7*Q9*Q12*Q14+2*Q1*Q7*Q11*Q13*Q14+Q1*Q6*Q13^2* Q14+Q3*Q7^2*Q9*Q15,
           -Q5^2*Q7*Q9^2+Q4^2*Q8*Q9*Q11+Q4^2*Q7*Q11^2- Q3*Q4*Q7*Q9*Q12+Q1*Q4*Q9*Q12^2-Q4^2*Q6*Q11*Q13-Q3*Q4*Q7*Q11*Q13+Q1*Q4*Q11*Q12*Q13+Q3^2*Q7*Q13^2-Q1*Q3*Q12*Q13^2+Q3*Q4* Q7*Q9*Q14-Q1*Q4*Q11*Q13*Q14-Q1*Q3*Q13^2*Q14-Q1*Q4*Q9*Q12*Q15+2*Q1*Q5*Q9*Q13*Q15- Q1*Q4*Q11*Q13*Q15+Q1*Q3*Q13^2*Q15,
           -Q5*Q8^2*Q9^2+Q5*Q7^2*Q11^2-Q5*Q6*Q7*Q9*Q12+Q5*Q6*Q7*Q11*Q13+Q1*Q8*Q11*Q12*Q13+ Q5*Q6^2*Q13^2-Q5*Q6*Q7*Q9*Q14-Q1*Q6*Q12* Q13*Q14-Q5*Q6*Q7*Q9*Q15+Q1*Q8*Q9*Q12* Q15+Q1*Q8*Q11*Q13*Q15-Q1*Q7*Q11*Q14*Q15,
           Q5*Q6*Q7-Q1*Q8*Q12,
           -Q4^2*Q6+Q3*Q4*Q7-Q1*Q5*Q13+Q1*Q4*Q15,
           Q5*Q8*Q9*Q12+Q4*Q6*Q12^2+Q5*Q8*Q11*Q13-Q3*Q8*Q12*Q13-Q5*Q7*Q11*Q14+Q1*Q12^2* Q14-Q5*Q6*Q13*Q14+Q5*Q8*Q9*Q15+Q5*Q7*Q11*Q15+Q4*Q6*Q12*Q15-Q2*Q9*Q12*Q15-Q1* Q12^2*Q15-Q5*Q6*Q13*Q15-Q3*Q8*Q13*Q15+Q1*Q12*Q14*Q15-Q1*Q12*Q15^2,
           Q5^2*Q9^2*Q12+Q4*Q5*Q9*Q11*Q12+Q4^2*Q11^2* Q12-Q3*Q4*Q9*Q12^2-Q3^2*Q12*Q13^2-Q4^2*Q11^2*Q14-Q3*Q4*Q9*Q12*Q14+Q3*Q4*Q11*Q13*Q14- Q5^2*Q9^2*Q15-Q3*Q4*Q11*Q13*Q15+Q3^2*Q13^2*Q15+ Q3*Q4*Q9*Q14*Q15,
           -Q4*Q8^2*Q9^2-Q4*Q7*Q8*Q9*Q11+Q4*Q6*Q8*Q9*Q13+3*Q1*Q8*Q9*Q12*Q13+Q1*Q7*Q11* Q12*Q13-2*Q1*Q8*Q11*Q13^2+Q1*Q6*Q12*Q13^2+Q4* Q6*Q7*Q9*Q14-Q1*Q6*Q13^2*Q14-Q1*Q7*Q9* Q12*Q15-Q1*Q8*Q9*Q13*Q15,
           2*Q5^2*Q8^2*Q13+Q5^2*Q7*Q8*Q14-Q2*Q4*Q8* Q12*Q14-Q2*Q5*Q8*Q13*Q14-Q5^2*Q7*Q8*Q15-Q4*Q5*Q8^2*Q15+Q2*Q4*Q8*Q12*Q15-2*Q2*Q5* Q8*Q13*Q15+Q2*Q4*Q8*Q14*Q15+Q2^2*Q13*Q15^ 2,
           -Q5*Q6*Q13+Q2*Q11*Q13+Q3*Q7*Q15-Q1* Q12*Q15,
           -Q5*Q6*Q7*Q9+Q4*Q6*Q8*Q9+Q1*Q8*Q9*Q12+Q1*Q7*Q11*Q12-Q4*Q6^2*Q13+Q3* Q6*Q7*Q13-2*Q1*Q8*Q11*Q13-Q1*Q8*Q9*Q15+Q1*Q6*Q13*Q15,
           Q4*Q7*Q8*Q9*Q11^2+Q4*Q7^2* Q11^3-Q4*Q6*Q8*Q9^2*Q12-Q4*Q6*Q7*Q9*Q11*Q12-2*Q4*Q6*Q7*Q11^2*Q13+Q2*Q6*Q9^2*Q12* Q13+Q1*Q8*Q9*Q11*Q12*Q13+Q4*Q6^2*Q11*Q13^2- Q2*Q6*Q9*Q11*Q13^2+Q1*Q6*Q11*Q12*Q13^2+2*Q4*Q6*Q7*Q9*Q11*Q14-Q4*Q6^2*Q9*Q13*Q14+ Q2*Q6*Q9^2*Q13*Q14-Q1*Q8*Q9*Q11*Q13*Q14-Q1*Q7*Q11^2*Q13*Q14-Q1*Q6*Q9*Q12*Q13*Q14+ Q1*Q6*Q11*Q13^2*Q14-Q1*Q6*Q9*Q13*Q14^2-Q1* Q7*Q11^2*Q13*Q15+Q1*Q7*Q9*Q11*Q14*Q15,
           -Q1*Q5*Q12+Q2*Q3*Q13,
           -Q5^2*Q8*Q9-Q4*Q5*Q8*Q11+Q2*Q4*Q11*Q12+Q3*Q5*Q8*Q13-Q4*Q5*Q6*Q14+Q3*Q4*Q8*Q14-Q1*Q5*Q12*Q14-Q4*Q5*Q6*Q15+Q3*Q5*Q7*Q15-Q2*Q3*Q13*Q15+Q1*Q5*Q14*Q15+Q1*Q5*Q15^2,
           -Q4^2*Q6+Q2*Q4*Q9-Q1*Q5*Q13+ Q1*Q4*Q15,
           Q5*Q8*Q9+Q4*Q8*Q11-Q5*Q6* Q13-Q2*Q11*Q13+Q4*Q6*Q14+Q1*Q12*Q14+Q4*Q6*Q15-Q3*Q7*Q15-Q1*Q14*Q15-Q1*Q15^2,
           Q5*Q8*Q9*Q11-Q2*Q11^2*Q13+Q4*Q6*Q11*Q14-Q3*Q6*Q13*Q14-Q3*Q8*Q9*Q15-2*Q4*Q6* Q11*Q15+2*Q3*Q7*Q11*Q15-Q1*Q11*Q12*Q15+2*Q1*Q11*Q15^2,
           -Q5*Q8*Q9-Q4*Q8*Q11+Q5*Q6*Q13+ Q2*Q11*Q13+Q4*Q6*Q14+Q1*Q12*Q14+Q4*Q6*Q15-Q3*Q7*Q15-Q1*Q14*Q15-Q1*Q15^2,
           -Q5*Q8*Q9+Q4*Q8*Q11-Q5*Q6*Q13+Q2*Q11* Q13+Q4*Q6*Q14+Q1*Q12*Q14+Q4*Q6*Q15-Q3* Q7*Q15-Q1*Q14*Q15-Q1*Q15^2,
           Q5*Q8*Q9-Q4*Q8*Q11+Q5*Q6*Q13-Q2*Q11* Q13+Q4*Q6*Q14+Q1*Q12*Q14+Q4*Q6*Q15-Q3* Q7*Q15-Q1*Q14*Q15-Q1*Q15^2,
           -Q5*Q8^2*Q9-Q5*Q7*Q8*Q11+Q2*Q8*Q9*Q12+Q5*Q6*Q8*Q13-Q2*Q6*Q12*Q13+Q5*Q6*Q7*Q14-Q4*Q6*Q8*Q14+Q2*Q8*Q9*Q14-Q1*Q8*Q12*Q14-Q1*Q8*Q14^2+Q2*Q7*Q11*Q15+Q1* Q8*Q12*Q15,
           -Q3*Q4*Q6*Q8^2*Q11+Q3^2*Q7*Q8^2* Q11+Q2*Q5*Q6*Q7*Q11^2-2*Q2*Q3*Q7*Q8*Q11^2+Q2^2*Q7*Q11^3-Q1*Q2*Q8*Q11^2*Q12-Q3*Q5* Q6^2*Q7*Q14+Q3*Q4*Q6^2*Q8*Q14-Q2*Q3*Q6* Q7*Q11*Q14+Q1*Q5*Q6*Q8*Q11*Q14+Q1*Q3* Q8^2*Q11*Q14+Q3*Q4*Q6^2*Q8*Q15-Q3^2*Q6*Q7*Q8*Q15+2*Q2*Q3*Q6*Q7*Q11*Q15-Q1*Q5* Q6*Q8*Q11*Q15+Q1*Q3*Q8^2*Q11*Q15-Q1*Q2*Q8*Q11^2*Q15+Q1*Q2*Q6*Q11*Q12*Q15-Q1*Q3* Q6*Q8*Q14*Q15-Q1*Q3*Q6*Q8*Q15^2,
           -Q5*Q8^2*Q9+Q5*Q7*Q8*Q11-Q5*Q6*Q8*Q13+Q2*Q6*Q12*Q13+Q5*Q6*Q7*Q14-Q3*Q7*Q8*Q14+Q2*Q8*Q9*Q14-Q1*Q8*Q14^2+Q2*Q8* Q9*Q15-Q2*Q7*Q11*Q15+Q1*Q8*Q12*Q15-Q1*Q8*Q14*Q15,
           -Q5*Q7*Q11-Q4*Q8*Q11+Q4*Q6* Q12+Q2*Q11*Q13+Q1*Q12*Q14-Q1*Q12*Q15,
           Q5^2*Q6*Q11*Q12-Q2*Q5*Q11^2*Q12-Q5^2*Q6* Q11*Q14+Q2*Q5*Q11^2*Q14+Q2*Q3*Q11*Q12*Q14-Q3^2*Q8*Q14^2+Q5^2*Q6*Q11*Q15-Q2*Q5*Q11^2*Q15-Q2*Q3*Q11*Q12*Q15-Q3*Q5*Q6*Q14*Q15+Q3^2*Q8*Q14*Q15+Q2*Q3*Q11*Q14*Q15,
           Q4*Q5*Q6*Q7*Q8-Q4^2*Q6*Q8^2-Q1*Q4* Q8^2*Q12-Q2*Q5*Q6*Q7*Q13+2*Q2*Q4*Q6*Q8*Q13-Q2^2*Q7*Q11*Q13+Q1*Q2*Q8*Q12*Q13-Q2^2* Q6*Q13^2+Q2*Q4*Q6*Q7*Q14-Q1*Q4*Q8^2*Q14+ Q1*Q2*Q8*Q13*Q15,
           Q3*Q5*Q7^2-Q3*Q4*Q7*Q8-Q2*Q5*Q7*Q9+Q2*Q4*Q7*Q11-2*Q1*Q5*Q8*Q13+Q1*Q2* Q12*Q13+Q1*Q5*Q7*Q14+Q1*Q5*Q7*Q15-Q1*Q2*Q13*Q15,
           Q5*Q7*Q8*Q9+Q4*Q8^2*Q9-Q2*Q7*Q11*Q13-2*Q1*Q8*Q12*Q13+Q2*Q6*Q13^2+Q1* Q7*Q12*Q14+Q4*Q6*Q7*Q15-Q3*Q7^2*Q15+Q1*Q7*Q12*Q15-Q1*Q7*Q14*Q15-Q1*Q7*Q15^2,
           Q5*Q8*Q9*Q11-Q4*Q8*Q11^2-Q5*Q6*Q9*Q12-Q3*Q7*Q11*Q12+Q2*Q9*Q11*Q12+Q1*Q11*Q12^2-Q3*Q8*Q11*Q13+Q3*Q6*Q13*Q14+Q3*Q8* Q9*Q15-Q3*Q7*Q11*Q15+Q1*Q11*Q12*Q15-Q1*Q11*Q14*Q15,
           -Q5*Q7*Q8*Q9+Q4*Q8^2*Q9+Q2* Q7*Q11*Q13-2*Q1*Q8*Q12*Q13+Q2*Q6*Q13^2+Q1*Q7*Q12*Q14+Q4*Q6*Q7*Q15-Q3*Q7^2*Q15+Q1* Q7*Q12*Q15-Q1*Q7*Q14*Q15-Q1*Q7*Q15^2,
           -Q5*Q9^2*Q12^2-Q4*Q11^2*Q12*Q13+Q3*Q11*Q12*Q13^2+Q4*Q11^2*Q13*Q14+Q3*Q9*Q12*Q13*Q14-Q3* Q11*Q13^2*Q14+Q5*Q9^2*Q12*Q15-Q4*Q9*Q11*Q12* Q15+Q4*Q11^2*Q13*Q15+Q3*Q9*Q12*Q13*Q15-Q3*Q11*Q13^2*Q15-Q3*Q9*Q13*Q14*Q15,
           -Q5^2*Q7^2*Q9-Q4*Q5*Q7*Q8*Q9+Q2*Q4*Q7*Q9*Q12+Q2*Q5*Q7*Q9*Q13-2*Q1*Q5*Q8* Q13^2-Q1*Q2*Q12*Q13^2-Q1*Q5*Q7*Q13*Q14+3*Q1* Q5*Q7*Q13*Q15+Q1*Q4*Q8*Q13*Q15+Q1*Q2* Q13^2*Q15-Q1*Q4*Q7*Q14*Q15,
           -Q4*Q6*Q8*Q9*Q11-Q4*Q6*Q7*Q11^2+Q4*Q6^2*Q9*Q12+Q4*Q6^2*Q11*Q13+2*Q1*Q8*Q11^2*Q13- Q1*Q8*Q9*Q11*Q14-Q1*Q7*Q11^2*Q14+Q1*Q6* Q9*Q12*Q14+Q1*Q7*Q11^2*Q15-Q1*Q6*Q9*Q12* Q15-Q1*Q6*Q11*Q13*Q15,
           -Q5*Q6*Q9+Q3*Q8*Q9+Q2*Q9*Q11-Q1*Q11* Q14,
           -Q5*Q7*Q9+Q4*Q8*Q9+Q4*Q7*Q11-Q4*Q6*Q13,
           Q5*Q6*Q9-Q3*Q8*Q9+Q2*Q9*Q11-Q1*Q11*Q14,
           -Q5*Q7^2*Q8*Q9+Q4*Q7*Q8^2* Q9-Q5*Q7^3*Q11+Q4*Q7^2*Q8*Q11+Q2*Q7^2*Q11* Q13-Q1*Q8^2*Q13^2-Q1*Q7*Q8*Q13*Q14+Q1* Q7^2*Q14*Q15,
           Q4*Q8*Q9*Q12+Q5*Q8*Q9*Q13- Q4*Q8*Q11*Q13-Q4*Q6*Q12*Q13-Q3*Q8*Q13^2+Q4*Q7*Q11*Q14-Q1*Q12*Q13*Q14-Q5*Q7*Q9* Q15-Q4*Q6*Q13*Q15+Q2*Q9*Q13*Q15+Q1*Q12*Q13*Q15+Q1*Q13*Q15^2,
           -Q1*Q12^2+Q2*Q9*Q15,
           Q5*Q8*Q9*Q11-Q4*Q8*Q11^2-Q5*Q6*Q9*Q12+Q4*Q6*Q11*Q12-Q3*Q8*Q11*Q13+Q1*Q11*Q12*Q14+Q3*Q6*Q13*Q14+Q3*Q8*Q9*Q15+Q4*Q6*Q11*Q15-Q2*Q9*Q11*Q15-Q1*Q11*Q12*Q15-Q1*Q11*Q15^2,
           2*Q4*Q5*Q8^2*Q12*Q13-2*Q5^2*Q8^2*Q13^2+2*Q2*Q5*Q8*Q12*Q13^2+Q4^2*Q8^2*Q12*Q14+Q5^2 *Q7*Q8*Q13*Q14-2*Q4*Q5*Q8^2*Q13*Q14-Q2*Q5*Q8*Q13^2*Q14+Q4*Q5*Q7*Q8*Q14^2-Q4*Q5* Q7*Q8*Q12*Q15-Q4^2*Q8^2*Q12*Q15+Q5^2*Q7*Q8* Q13*Q15+Q4*Q5*Q8^2*Q13*Q15-Q2*Q5*Q7*Q12*Q13*Q15-3*Q2*Q4*Q8*Q12*Q13*Q15-Q5^2*Q7^2* Q14*Q15+Q4^2*Q8^2*Q14*Q15+Q2*Q4*Q7*Q12*Q15^2+ Q2*Q5*Q7*Q13*Q15^2+Q2^2*Q13^2*Q15^2-Q2*Q4* Q7*Q14*Q15^2,
           -Q3^2*Q4*Q8*Q11+Q2*Q3*Q4*Q11^ 2-Q1*Q5^2*Q11^2+Q3^2*Q5*Q6*Q13-Q3^3*Q8*Q13+ Q2*Q3^2*Q11*Q13-Q1*Q3*Q5*Q11*Q15+Q1*Q3^2* Q14*Q15,
           Q5^2*Q9^2*Q12^2+Q4*Q5*Q9*Q11*Q12^2+Q4*Q5*Q11^2*Q12*Q13-2*Q5^2*Q11^2*Q13^2+Q3*Q5* Q11*Q12*Q13^2-Q4^2*Q11^2*Q12*Q14-Q3*Q4*Q9*Q12^2* Q14-Q5^2*Q9*Q11*Q13*Q14+Q4*Q5*Q11^2*Q13* Q14-2*Q3*Q5*Q11*Q13^2*Q14+Q3^2*Q12*Q13^2*Q14+Q3* Q4*Q11*Q13*Q14^2-Q4*Q5*Q9*Q11*Q12*Q15+Q3* Q4*Q9*Q12^2*Q15+2*Q5^2*Q9*Q11*Q13*Q15-3*Q3* Q5*Q9*Q12*Q13*Q15-Q3*Q4*Q11*Q12*Q13*Q15+2*Q3*Q5*Q11*Q13^2*Q15-Q3^2*Q12*Q13^2*Q15+Q3^2* Q13^2*Q14*Q15,
           Q11*Q13-Q9*Q15,
           Q5*Q6*Q7-Q1*Q8*Q15,
           -Q8*Q9+Q6*Q13,
           Q10+Q12-Q14-Q15,
           Q8*Q12*Q13-Q7*Q15^2,
           Q8*Q9-Q6*Q13,
           Q10-Q12+Q14-Q15,
           Q10-Q12-Q14+Q15,
           -2*Q5^2*Q8^2*Q11^2+Q3*Q5*Q8^2*Q11*Q12+Q2*Q5* Q8*Q11^2*Q12+Q5^2*Q6^2*Q12^2+Q2*Q5*Q6*Q11*Q12^2+2*Q5^2*Q6*Q8*Q11*Q14+2*Q3*Q5*Q8^2*Q11* Q14-3*Q3*Q5*Q6*Q8*Q12*Q14-Q3^2*Q8^2*Q12*Q14-Q2*Q5*Q6*Q11*Q12*Q14-Q2*Q3*Q8*Q11*Q12* Q14+Q2*Q3*Q6*Q12^2*Q14-Q5^2*Q6*Q8*Q11*Q15-2* Q3*Q5*Q8^2*Q11*Q15+Q2*Q5*Q8*Q11^2*Q15+ Q3^2*Q8^2*Q12*Q15-Q2^2*Q11^2*Q12*Q15-Q2*Q3*Q6* Q12^2*Q15+Q3^2*Q8^2*Q14*Q15+Q2*Q3*Q8*Q11*Q15^2,
           -Q5*Q7*Q11*Q12+Q4*Q8*Q11*Q12+Q2*Q9* Q12^2+Q5*Q8*Q11*Q13-Q5*Q8*Q9*Q14-Q5*Q7*Q11*Q14+Q4*Q8*Q11*Q14-Q4*Q6*Q12*Q14+Q2* Q9*Q12*Q14-Q1*Q12^2*Q14+Q3*Q8*Q13*Q14-Q1*Q12*Q14^2-Q5*Q8*Q9*Q15+Q1*Q12^2*Q15-Q3*Q8* Q13*Q15+Q1*Q12*Q14*Q15,
           -Q5*Q9*Q12+Q3*Q13*Q14,
           Q5*Q9*Q12-Q3*Q13*Q14,
           Q5*Q6*Q13-Q2*Q11*Q13+Q3*Q7*Q14-Q1*Q12* Q14,
           Q5*Q6*Q8*Q9+Q3*Q8^2*Q9-2*Q3*Q7* Q8*Q11+Q2*Q7*Q11^2-Q2*Q6*Q9*Q12+2*Q1*Q8*Q11*Q12-Q5*Q6^2*Q13-2*Q1*Q8*Q11*Q14+Q3*Q6* Q7*Q15-Q1*Q6*Q12*Q15+Q1*Q6*Q14*Q15,
           Q5*Q9*Q11+Q4*Q11^2-Q3*Q11*Q13-Q3*Q9*Q15,
           Q5*Q9*Q11-Q4*Q11^2+Q3*Q11*Q13-Q3*Q9*Q15,
           Q2*Q3*Q4-Q1*Q5^2,
           -Q5*Q9*Q11+Q4*Q11^2+Q3*Q11*Q13-Q3*Q9*Q15,
           -Q5*Q9*Q11-Q4*Q11^2+Q3*Q11*Q13+Q3*Q9* Q15,
           Q5*Q9*Q11-Q4*Q11^2-Q3*Q11*Q13+Q3*Q9* Q15,
           -Q5*Q6*Q7*Q9*Q11+Q4*Q6^2*Q9*Q12-Q5* Q6^2*Q9*Q13+Q4*Q6^2*Q11*Q13-Q3*Q6^2*Q13^2+Q1*Q7*Q11^2*Q14+2*Q1*Q6*Q11*Q13*Q14+Q1*Q8* Q9*Q11*Q15-Q1*Q7*Q11^2*Q15-Q1*Q6*Q11*Q13* Q15-Q1*Q6*Q9*Q14*Q15,
           -Q5^2*Q6*Q7^2-Q2*Q4*Q7*Q8*Q11-Q1*Q4*Q8^2*Q12+Q2*Q5*Q6*Q7*Q13-Q1*Q5*Q8^2*Q13+ Q1*Q2*Q8*Q12*Q13+2*Q1*Q5*Q7*Q8*Q14+Q1*Q4*Q8^2*Q14-Q1*Q2*Q7*Q12*Q14+Q1*Q2*Q7* Q12*Q15-Q1*Q2*Q7*Q14*Q15,
           Q4*Q5*Q9*Q11*Q12^2+Q5^2*Q9*Q11*Q12*Q13-2* Q4*Q5*Q11^2*Q12*Q13-2*Q5^2*Q11^2*Q13^2-Q3*Q5*Q11*Q12*Q13^2+Q4^2*Q11^2*Q12*Q14+2*Q4*Q5*Q11^2*Q13* Q14+2*Q3*Q5*Q11*Q13^2*Q14-Q5^2*Q9^2*Q12*Q15+Q4^2*Q11^2*Q12*Q15+Q5^2*Q9*Q11*Q13*Q15+Q4*Q5* Q11^2*Q13*Q15-Q4*Q5*Q9*Q11*Q14*Q15-Q4^2*Q11^2*Q14*Q15-Q3*Q5*Q9*Q13*Q14*Q15-3*Q3*Q4* Q11*Q13*Q14*Q15-Q3*Q4*Q9*Q12*Q15^2+Q3*Q5*Q9*Q13*Q15^2+Q3^2*Q13^2*Q15^2+Q3*Q4*Q9*Q14*Q15^2,
           Q5*Q6*Q9+Q4*Q6*Q11-Q3*Q6*Q13-Q1*Q11*Q15,
           Q3*Q5*Q7^2*Q8-Q3*Q4*Q7*Q8^2-Q1*Q5*Q7*Q8*Q12+Q1*Q4*Q8^2*Q12-Q2*Q3*Q7* Q8*Q13+2*Q1*Q5*Q8^2*Q13-Q1*Q4*Q8^2*Q14-Q1*Q2*Q8*Q13*Q14+Q2*Q3*Q7^2*Q15-Q1*Q2* Q7*Q12*Q15+Q1*Q2*Q7*Q14*Q15,
           -Q5*Q6*Q9+Q4*Q6*Q11+Q3*Q6*Q13-Q1*Q11* Q15,
           Q4*Q5*Q6*Q7-Q4^2*Q6*Q8+Q2*Q4*Q8*Q9-Q2*Q4*Q7*Q11-2*Q1*Q5*Q8*Q13+Q1* Q5*Q7*Q14+Q1*Q2*Q13*Q14-Q1*Q5*Q7*Q15+Q1*Q4*Q8*Q15,
           -Q4^2*Q8^2*Q12-Q2*Q5*Q7*Q12* Q13+Q2^2*Q12*Q13^2+Q5^2*Q7^2*Q14+Q4*Q5*Q7*Q8*Q14+Q4^2*Q8^2*Q14-Q2^2*Q13^2*Q14-Q2*Q4*Q7* Q14^2-Q5^2*Q7^2*Q15+Q2*Q4*Q7*Q12*Q15+Q2*Q5*Q7*Q13*Q15-Q2*Q4*Q7*Q14*Q15,
           Q2*Q3*Q7-Q1*Q5*Q8-Q1*Q2*Q12+Q1*Q2* Q15,
           Q4*Q5*Q6*Q9*Q12-Q4^2*Q6*Q11*Q12-Q1*Q5*Q11*Q12*Q13+Q3*Q5*Q6*Q13^2+Q1*Q5* Q9*Q12*Q14-Q1*Q4*Q11*Q12*Q14-Q1*Q5*Q11*Q13*Q14-Q1*Q5*Q9*Q12*Q15+Q1*Q4*Q11*Q12*Q15-Q3*Q4*Q6*Q13*Q15+Q1*Q5*Q11*Q13*Q15+Q1*Q3*Q13*Q15^2,
           Q3*Q7*Q8*Q11*Q13-Q2*Q7*Q11^2*Q13-Q5*Q6^2*Q13^2-Q3*Q7^2*Q11*Q14+Q2*Q7* Q9*Q11*Q14+2*Q3*Q6*Q7*Q13*Q14-Q1*Q8*Q11* Q13*Q14+Q1*Q6*Q12*Q13*Q14-Q1*Q7*Q11*Q14^2 -Q3*Q7*Q8*Q9*Q15+Q2*Q7*Q9*Q11*Q15+Q1* Q8*Q11*Q13*Q15-Q1*Q7*Q11*Q14*Q15,
           Q5^3*Q6^2*Q9+Q3*Q5^2*Q6*Q8*Q9-2*Q2*Q5^2* Q6*Q9*Q11-Q2*Q3*Q5*Q7*Q11^2+Q2^2*Q5*Q9* Q11^2+2*Q2*Q3*Q5*Q6*Q9*Q12-Q1*Q5^2*Q6* Q11*Q12+Q2*Q3^2*Q7*Q11*Q12-Q1*Q3*Q5*Q8*Q11*Q12-Q2^2*Q3*Q9*Q11*Q12+Q1*Q2*Q5*Q11^2* Q12-Q1*Q2*Q3*Q11*Q12^2-Q1*Q5^2*Q6*Q11*Q14+ Q1*Q3*Q5*Q6*Q12*Q14-Q2*Q3*Q5*Q6*Q9* Q15-Q2*Q3^2*Q8*Q9*Q15+Q2*Q3^2*Q7*Q11*Q15+Q1*Q3*Q5*Q8*Q11*Q15+Q1*Q2*Q5*Q11^2* Q15-Q1*Q2*Q3*Q11*Q12*Q15,
           Q2*Q4*Q11-Q1*Q5*Q12,
           -Q2*Q4*Q11+Q1*Q5*Q12,
           -Q5*Q6*Q13+Q3*Q8*Q13+Q2*Q9*Q15-Q1*Q14* Q15,
           -Q5*Q8*Q9+Q5*Q7*Q11-Q5*Q6*Q13+Q2*Q9*Q15+Q1*Q12*Q15-Q1*Q14*Q15,
           -Q5^2*Q6*Q8*Q12+Q3*Q5*Q8^2*Q12-Q5^2*Q6*Q8*Q14+Q3*Q5*Q8^2*Q14-Q2*Q5*Q6*Q12*Q14+ Q5^2*Q6*Q8*Q15-Q3*Q5*Q8^2*Q15+Q2*Q5*Q6* Q12*Q15-Q2*Q3*Q8*Q12*Q15+Q2^2*Q11*Q12*Q15+ Q2*Q5*Q6*Q14*Q15-Q2^2*Q11*Q15^2,
           -Q6*Q12*Q13+Q8*Q9*Q15,
           Q3*Q8^2*Q11*Q12-Q2*Q8*Q11^2*Q12-Q3*Q8^2* Q11*Q14+Q2*Q8*Q11^2*Q14+Q5*Q6^2*Q12*Q14-Q3*Q6*Q8*Q12*Q14+Q2*Q6*Q11*Q12*Q14-Q5*Q6^2* Q14^2+Q3*Q8^2*Q11*Q15-Q2*Q8*Q11^2*Q15-Q2*Q6*Q11*Q12*Q15+Q2*Q6*Q11*Q14*Q15,
           Q4*Q5*Q6*Q12-Q2*Q4*Q11*Q12-Q5^2*Q6*Q13-Q3*Q5*Q8*Q13+Q2*Q5*Q11*Q13+Q3*Q4*Q8*Q14+Q1*Q5*Q12*Q14+Q4*Q5*Q6*Q15-Q2*Q5*Q9*Q15-Q1*Q5*Q12*Q15+Q2*Q3*Q13*Q15-Q1*Q5*Q15^2,
           -2*Q5*Q8^2*Q13-Q2*Q8*Q12*Q13-Q5* Q7*Q8*Q14+Q5*Q7*Q8*Q15+Q4*Q8^2*Q15+Q2* Q8*Q13*Q15+Q2*Q7*Q15^2,
           -Q4*Q8*Q9*Q12+Q5*Q8*Q9*Q13-Q5*Q7*Q11* Q13-Q3*Q7*Q12*Q13+Q2*Q9*Q12*Q13+Q1*Q12^2 *Q13-Q5*Q6*Q13^2+Q4*Q7*Q11*Q14+Q5*Q7*Q9*Q15-Q3*Q7*Q13*Q15+Q1*Q12*Q13*Q15-Q1*Q13*Q14*Q15,
           Q5^3*Q6*Q9^2-2*Q4*Q5^2*Q6*Q9*Q11+ Q4^2*Q5*Q6*Q11^2-Q3*Q4*Q5*Q7*Q11^2+2*Q3*Q4*Q5*Q6*Q9*Q12-Q3*Q4^2*Q6*Q11*Q12+Q3^2* Q4*Q7*Q11*Q12-Q1*Q5^2*Q9*Q11*Q12+Q1*Q4* Q5*Q11^2*Q12-Q1*Q3*Q4*Q11*Q12^2+Q3*Q5^2*Q6*Q9*Q13-Q1*Q3*Q5*Q11*Q12*Q13-Q3*Q4* Q5*Q6*Q9*Q14+Q3^2*Q4*Q7*Q11*Q14+Q1*Q4*Q5*Q11^2*Q14-Q1*Q3*Q4*Q11*Q12*Q14-Q3^2*Q4* Q6*Q13*Q14+Q1*Q3*Q5*Q11*Q13*Q14-Q1*Q5^2*Q9*Q11*Q15+Q1*Q3*Q5*Q9*Q12*Q15,
           Q2*Q5*Q6*Q9+Q2*Q3*Q7*Q11-2*Q1*Q5*Q8*Q11-Q2^2*Q9*Q11+Q1*Q5*Q6*Q12+Q1*Q3* Q8*Q12-Q2*Q3*Q6*Q13-Q1*Q5*Q6*Q14+Q1*Q2*Q11*Q14,
           -Q2*Q12*Q13+Q5*Q7*Q14,
           Q2*Q12*Q13-Q5*Q7*Q14,
           Q3*Q5^2*Q6*Q9-Q3^2*Q5*Q8*Q9-Q1*Q5^2* Q11^2+Q3^2*Q5*Q6*Q13-Q3^3*Q8*Q13+Q2*Q3^2*Q11* Q13-Q1*Q3*Q5*Q11*Q14+Q1*Q3^2*Q14*Q15,
           Q4*Q5*Q6^2-Q3*Q5*Q6*Q7-Q2*Q4*Q6*Q11-2*Q1*Q5*Q8*Q11+Q1*Q5*Q6*Q12+Q2*Q3* Q6*Q13+Q1*Q5*Q6*Q14-Q1*Q3*Q8*Q14+Q1*Q3*Q8*Q15,
           -Q5*Q7*Q8*Q9^2-Q4*Q8^2*Q9^2-Q4* Q7*Q8*Q9*Q11+2*Q4*Q6*Q8*Q9*Q13+Q4*Q6* Q7*Q11*Q13+Q1*Q8*Q9*Q12*Q13-Q4*Q6^2*Q13^ 2-Q1*Q6*Q12*Q13^2+Q4*Q6*Q7*Q9*Q14-Q1* Q6*Q13^2*Q14+Q1*Q8*Q9*Q13*Q15,
           Q4*Q8*Q9*Q11^2+Q4*Q7*Q11^3-Q5*Q6*Q9^2* Q12-Q4*Q6*Q9*Q11*Q12+2*Q2*Q9^2*Q11*Q12-2*Q4* Q6*Q11^2*Q13+Q2*Q9*Q11^2*Q13+Q3*Q6*Q11*Q13^2+2*Q4*Q6*Q9*Q11*Q14-2*Q2*Q9^2*Q11*Q14-Q1* Q9*Q11*Q12*Q14-Q3*Q6*Q9*Q13*Q14-2*Q1*Q11^2* Q13*Q14+2*Q1*Q9*Q11*Q14^2+Q1*Q9*Q11*Q12* Q15-Q1*Q9*Q11*Q14*Q15,
           Q5*Q8*Q9-Q4*Q8*Q11+Q3*Q7*Q15-Q1*Q12* Q15,
           Q2*Q8^2*Q9^2+Q5*Q6*Q7^2*Q11-Q3*Q7^2* Q8*Q11-Q2*Q7*Q8*Q9*Q11+Q2*Q7^2*Q11^2+Q1* Q8^2*Q9*Q12-Q1*Q7*Q8*Q11*Q12+2*Q1*Q6*Q8*Q12*Q13-Q2*Q6^2*Q13^2-Q2*Q6*Q7*Q9*Q14- Q1*Q8^2*Q9*Q14+Q1*Q7*Q8*Q11*Q14-Q1*Q6*Q7*Q12*Q14+Q1*Q6*Q7*Q14^2+Q2*Q6*Q7*Q9* Q15-Q1*Q8^2*Q9*Q15-Q1*Q7*Q8*Q11*Q15,
           Q4*Q6*Q9-Q2*Q9^2-Q1*Q11*Q13+Q1*Q9*Q14,
           -Q5^2*Q6*Q7+2*Q3*Q5*Q7*Q8-Q3*Q4*Q8^2- Q2*Q5*Q8*Q9-Q2*Q4*Q8*Q11+Q2*Q5*Q6*Q13-Q2*Q3*Q8*Q13+2*Q1*Q5*Q8*Q14+Q2*Q3* Q7*Q15-Q2^2*Q9*Q15-Q1*Q2*Q12*Q15+Q1*Q2* Q14*Q15,
           Q1*Q8^2*Q9*Q12+Q1*Q7*Q8*Q11*Q12+Q2*Q6*Q8*Q9*Q13-Q2*Q6*Q7*Q11*Q13-2* Q1*Q8^2*Q11*Q13+3*Q1*Q6*Q8*Q12*Q13-Q2*Q6^2*Q13^2-Q1*Q6*Q7*Q12*Q14-Q1*Q6*Q8*Q13*Q14+ Q2*Q6*Q7*Q9*Q15-Q1*Q8^2*Q9*Q15,
           Q2*Q7*Q11-Q1*Q8*Q14,
           -Q5*Q7*Q8*Q11-Q4*Q8^2*Q11-Q3*Q7*Q8*Q12+Q2*Q8*Q9*Q12+Q1*Q8*Q12^2+Q2*Q8*Q11* Q13-Q2*Q6*Q12*Q13+Q5*Q6*Q7*Q14-Q3*Q7*Q8*Q15+Q2*Q7*Q11*Q15+Q1*Q8*Q12*Q15-Q1*Q8*Q14*Q15,
           Q7*Q11^2-Q6*Q9*Q12,
           Q4*Q8^2*Q9-Q5*Q7^2*Q11+Q4*Q6*Q7*Q12-2* Q4*Q6*Q8*Q13+Q2*Q7*Q11*Q13+Q2*Q6*Q13^2+Q1*Q7*Q12*Q14-2*Q1*Q8*Q13*Q14-Q2*Q7*Q9* Q15-Q1*Q7*Q12*Q15+2*Q1*Q8*Q13*Q15,
           Q5*Q6*Q9*Q12-Q5*Q6*Q11*Q13+Q3*Q8*Q11* Q13-Q2*Q11^2*Q13-Q4*Q6*Q11*Q14-Q1*Q11*Q12* Q14+Q3*Q6*Q13*Q14-Q3*Q8*Q9*Q15-Q4*Q6* Q11*Q15+Q3*Q7*Q11*Q15+Q1*Q11*Q14*Q15+Q1* Q11*Q15^2,
           -Q2*Q8*Q9^2*Q12+Q2*Q7*Q11^2*Q13+Q2*Q6*Q9*Q12*Q13-Q1*Q8*Q11*Q12*Q13-Q2* Q7*Q9*Q11*Q14+Q1*Q8*Q9*Q12*Q14+Q1*Q8*Q11*Q13*Q14-Q1*Q6*Q12*Q13*Q14+Q1*Q7*Q11* Q14^2-Q1*Q8*Q9*Q12*Q15-Q1*Q8*Q11*Q13*Q15+Q1*Q6*Q12*Q13*Q15,
           -Q5*Q6*Q8*Q9^2-Q3*Q8^2*Q9^2-Q5*Q6*Q7*Q9*Q11+Q2*Q8*Q9^2*Q11+Q2*Q6*Q9^2*Q12-Q1* Q8*Q9*Q11*Q14-Q1*Q7*Q11^2*Q14+Q1*Q6*Q11* Q13*Q14+2*Q1*Q8*Q9*Q11*Q15+Q1*Q7*Q11^2*Q15-Q1*Q6*Q9*Q14*Q15,
           Q4*Q8*Q9*Q12-Q5*Q8*Q9*Q13+Q4*Q8*Q11* Q13-Q2*Q9*Q12*Q13-Q3*Q8*Q13^2-Q4*Q7*Q11* Q14+Q4*Q6*Q13*Q14-Q2*Q9*Q13*Q14+Q1*Q12* Q13*Q14+Q1*Q13*Q14^2+Q5*Q7*Q9*Q15-Q1*Q12* Q13*Q15,
           -Q4*Q8*Q11^2-Q5*Q6*Q9*Q12-2*Q3* Q7*Q11*Q12+2*Q2*Q9*Q11*Q12+2*Q1*Q11*Q12^2+Q5* Q6*Q11*Q13-Q1*Q11*Q12*Q14-Q3*Q8*Q9*Q15+ Q3*Q7*Q11*Q15,
           Q5*Q6*Q7+Q3*Q7*Q8-Q2*Q7*Q11-Q1*Q8* Q12,
           Q5*Q6*Q7-Q3*Q7*Q8+Q2*Q7*Q11-Q2* Q6*Q13,
           -Q5*Q6*Q7+Q3*Q7*Q8+Q2*Q7*Q11- Q1*Q8*Q12,
           Q3*Q4*Q6*Q8-Q3^2*Q7*Q8+Q2*Q3*Q7*Q11-2*Q1*Q5*Q8*Q11+Q1*Q3*Q8* Q12-Q1*Q2*Q11*Q12-Q2*Q3*Q6*Q13+Q1*Q5*Q6*Q15+Q1*Q2*Q11*Q15,
           -Q5*Q7*Q11+Q2*Q9*Q12+Q5*Q6*Q13-Q1*Q12* Q14,
           -Q5^2*Q8*Q9^2+Q3*Q5*Q7*Q11*Q13-Q3* Q5*Q6*Q13^2-Q3*Q4*Q7*Q11*Q14+Q3*Q4*Q6*Q13*Q14+Q1*Q5*Q11*Q13*Q14+2*Q3*Q5*Q7*Q9* Q15+Q1*Q5*Q9*Q12*Q15+Q3*Q4*Q6*Q13*Q15-Q3^2*Q7*Q13*Q15-Q1*Q5*Q11*Q13*Q15-Q1*Q3* Q13*Q14*Q15-Q1*Q3*Q13*Q15^2,
           -Q5^2*Q9*Q12*Q13+Q3*Q5*Q12*Q13^2+Q4*Q5*Q9*Q12*Q14+Q4^2*Q11*Q12*Q14+Q5^2*Q9*Q13*Q14- Q3*Q4*Q12*Q13*Q14-Q3*Q5*Q13^2*Q14-Q4^2*Q11* Q14^2-Q4*Q5*Q9*Q12*Q15-Q5^2*Q9*Q13*Q15+Q3*Q5*Q13^2*Q15+Q4*Q5*Q9*Q14*Q15,
           -Q5*Q7*Q8*Q12*Q13+2*Q5*Q8^2*Q13^2+Q2*Q8*Q12*Q13^2-Q4*Q7*Q8*Q12*Q14+Q4*Q7*Q8*Q12* Q15-2*Q5*Q7*Q8*Q13*Q15-Q4*Q8^2*Q13*Q15-Q2* Q8*Q13^2*Q15+Q4*Q7*Q8*Q14*Q15+Q5*Q7^2*Q15^2,
           -Q5^2*Q7*Q8*Q12-Q4*Q5*Q8^2*Q12+2*Q5^2*Q8^2*Q13-2*Q2*Q5*Q8*Q12*Q13+Q2^2*Q12^2*Q13+Q4* Q5*Q8^2*Q14+Q2*Q5*Q7*Q12*Q14-Q2*Q5*Q8* Q13*Q14+Q2*Q5*Q7*Q12*Q15-Q2*Q5*Q7*Q14* Q15,
           -Q5^2*Q8*Q9+Q4*Q5*Q8*Q11+Q2*Q5*Q9*Q12-Q2*Q4*Q11*Q12-Q3*Q5*Q8*Q13-Q4*Q5*Q6*Q14+Q3*Q4*Q8*Q14+Q2*Q5*Q9*Q14-Q1*Q5*Q12*Q14-Q1*Q5*Q14^2+Q1*Q5*Q12*Q15+ Q2*Q3*Q13*Q15,
           Q8^2*Q9-Q6*Q7*Q14,
           Q3*Q8*Q13-Q2*Q11*Q13+Q4*Q6*Q14-Q1*Q14* Q15,
           Q5*Q8*Q9-Q4*Q8*Q11-Q3*Q8*Q13+Q1* Q14^2,
           -Q4^2*Q8*Q12^2+Q5^2*Q7*Q12*Q13-Q2*Q5*Q12*Q13^2+Q4*Q5*Q7*Q12*Q14+Q4^2*Q8*Q12*Q14- Q5^2*Q7*Q13*Q14-Q2*Q4*Q12*Q13*Q14+Q2*Q5* Q13^2*Q14+Q4*Q5*Q7*Q12*Q15-Q5^2*Q7*Q13* Q15+Q2*Q5*Q13^2*Q15-Q4*Q5*Q7*Q14*Q15,
           -Q5*Q7*Q11*Q12+Q4*Q6*Q12^2+Q5*Q8*Q11*Q13+Q2*Q11*Q12*Q13+Q1*Q12^2*Q14-Q5*Q6*Q13* Q14-Q3*Q8*Q13*Q14-Q5*Q7*Q11*Q15+Q4*Q6*Q12*Q15-Q2*Q9*Q12*Q15-Q1*Q12^2*Q15-Q5*Q6* Q13*Q15+Q3*Q8*Q13*Q15+Q2*Q11*Q13*Q15+Q1*Q12*Q14*Q15-Q1*Q12*Q15^2,
           -Q5^2*Q6*Q8*Q12+Q3*Q5*Q8^2*Q12-Q2^2*Q11*Q12^2+Q5^2*Q6*Q8*Q14-Q3*Q5*Q8^2*Q14+Q2*Q3* Q8*Q12*Q14+Q5^2*Q6*Q8*Q15-Q3*Q5*Q8^2*Q15-Q2*Q5*Q6*Q12*Q15+Q2*Q3*Q8*Q12*Q15+Q2^2*Q11*Q12*Q15-Q2*Q3*Q8*Q14*Q15,
           Q3*Q8^2*Q11*Q13-Q2*Q8*Q11^2*Q13-Q5*Q6^2* Q12*Q13-Q3*Q7*Q8*Q11*Q14+Q2*Q7*Q11^2*Q14+Q5*Q6^2*Q13*Q14+Q3*Q6*Q8*Q13*Q14-Q2*Q6* Q11*Q13*Q14-Q3*Q7*Q8*Q11*Q15+Q2*Q7*Q11^2*Q15+Q3*Q6*Q8*Q13*Q15-Q3*Q6*Q7*Q14*Q15,
           Q2*Q4*Q11-Q1*Q5*Q15,
           -Q5^2*Q6*Q9^2-Q3*Q5*Q8*Q9^2-Q3*Q4*Q8*Q9*Q11+Q2*Q5*Q9^2*Q11+2*Q1*Q5*Q9*Q11*Q12+ Q1*Q4*Q11^2*Q12-Q1*Q5*Q9*Q11*Q14-Q1*Q4* Q11^2*Q14-Q1*Q3*Q9*Q12*Q14+Q1*Q3*Q11*Q13* Q14+Q2*Q3*Q9^2*Q15,
           -Q8^2*Q9^2*Q12-Q6*Q7*Q11*Q12*Q13+Q6^2*Q12*Q13^2-Q7^2*Q11^2*Q14+Q6*Q7*Q9*Q12*Q14+Q6*Q7* Q11*Q13*Q14+Q8^2*Q9^2*Q15+Q7*Q8*Q9*Q11*Q15+Q7^2*Q11^2*Q15-Q6^2*Q13^2*Q15-Q6*Q7*Q9*Q14* Q15-Q6*Q7*Q9*Q15^2,
           -Q5*Q6*Q7-Q2*Q8*Q9+Q2*Q7*Q11+Q2*Q6* Q13,
           -Q5*Q7*Q11+Q4*Q8*Q11-Q2*Q11*Q13+Q1*Q15^2,
           -Q4*Q8*Q9*Q12+Q5*Q7*Q11*Q13+Q2* Q9*Q12*Q13-Q3*Q8*Q13^2-Q4*Q7*Q11*Q14+2*Q4*Q6*Q13*Q14-2*Q2*Q9*Q13*Q14+2*Q1*Q13*Q14^2-Q1*Q13*Q14*Q15,
           Q4*Q5*Q11^2*Q12+2*Q5^2*Q11^2*Q13- Q3*Q5*Q11*Q12*Q13-Q5^2*Q9*Q11*Q14-Q4*Q5* Q11^2*Q14+Q3*Q5*Q9*Q12*Q14-2*Q3*Q5*Q11* Q13*Q14+Q3^2*Q13*Q14^2-Q3*Q5*Q9*Q12*Q15+Q3*Q5*Q9*Q14*Q15,
           -Q3*Q5*Q6*Q7*Q11+Q3^2*Q7* Q8*Q11-Q2*Q3*Q7*Q11^2+2*Q1*Q5*Q8*Q11^2-Q1*Q3*Q8*Q11*Q12+Q1*Q2*Q11^2*Q12+Q3^2*Q6* Q7*Q14-Q1*Q3*Q6*Q12*Q14-Q1*Q5*Q6*Q11*Q15-Q1*Q2*Q11^2*Q15+Q1*Q3*Q6*Q14*Q15,
           Q2^2*Q8*Q9*Q11-Q2^2*Q7*Q11^2+Q2^2*Q6*Q9* Q12-Q2*Q3*Q6*Q8*Q13-Q2^2*Q6*Q11*Q13+Q1*Q5*Q6*Q8*Q14-Q1*Q3*Q8^2*Q14-Q1*Q2*Q8* Q11*Q14+Q1*Q3*Q8^2*Q15+2*Q1*Q2*Q8*Q11*Q15- Q1*Q2*Q6*Q14*Q15,
           -Q5^2*Q6*Q9^2-Q3*Q4*Q8*Q9*Q11+2*Q1*Q5*Q9*Q11*Q12+Q1*Q4*Q11^2*Q12+Q3*Q5*Q6*Q9* Q13-Q1*Q5*Q11^2*Q13-Q1*Q4*Q11^2*Q14-Q1*Q3* Q9*Q12*Q14+Q1*Q3*Q11*Q13*Q14-Q1*Q3*Q9* Q12*Q15+Q1*Q3*Q9*Q14*Q15,
           -Q8*Q11*Q13+Q7*Q11*Q14+Q6*Q13*Q14-Q7*Q11* Q15,
           -Q1*Q5*Q8*Q12*Q13+Q2*Q5*Q6*Q13^2+ Q4*Q5*Q6*Q7*Q14-Q4^2*Q6*Q8*Q14+Q1*Q5*Q7*Q12*Q14-Q1*Q4*Q8*Q12*Q14-Q1*Q5*Q8* Q13*Q14-Q2*Q4*Q6*Q13*Q15+Q1*Q5*Q8*Q13*Q15-Q1*Q5*Q7*Q14*Q15+Q1*Q4*Q8*Q14*Q15+Q1*Q2*Q13*Q15^2,
           Q1*Q5*Q8*Q11*Q12-Q5^2*Q6^2*Q13+Q3^2*Q8^2*Q13+Q2*Q3*Q8*Q11*Q13+Q2^2*Q11^2*Q13-Q2*Q3*Q6*Q12*Q13+Q1*Q5*Q8*Q11* Q14+Q1*Q5*Q6*Q12*Q14-Q2*Q3*Q6*Q13*Q14-Q1*Q2*Q11*Q12*Q15-Q2*Q3*Q6*Q13*Q15-Q1*Q3*Q8*Q14*Q15,
           -Q4*Q7*Q8*Q9*Q12+Q4*Q8^2* Q9*Q13-2*Q3*Q7^2*Q12*Q13+2*Q2*Q7*Q9*Q12*Q13+2* Q1*Q7*Q12^2*Q13+Q3*Q7*Q8*Q13^2-2*Q2*Q8* Q9*Q13^2+Q2*Q7*Q11*Q13^2-2*Q1*Q8*Q12*Q13^2+Q2* Q6*Q13^3-Q1*Q7*Q12*Q13*Q14-Q5*Q7^2*Q9* Q15+2*Q3*Q7^2*Q13*Q15-Q2*Q7*Q9*Q13*Q15-Q1*Q7*Q12*Q13*Q15+Q1*Q7*Q13*Q14*Q15,
           -Q5^2*Q6^2*Q7+Q2*Q5*Q6*Q7*Q11-Q1*Q5*Q8^2*Q11-Q2*Q3*Q6*Q8*Q13+2*Q1*Q5*Q6*Q8* Q14+Q1*Q3*Q8^2*Q14-Q1*Q2*Q6*Q12*Q14-Q1*Q3*Q8^2*Q15+Q1*Q2*Q8*Q11*Q15+Q1*Q2*Q6* Q12*Q15-Q1*Q2*Q6*Q14*Q15,
           Q3*Q4*Q8^2*Q11-Q4*Q5*Q6^2*Q12+Q2*Q4* Q6*Q11*Q12-Q1*Q5*Q8*Q11*Q12-Q1*Q5*Q8*Q11*Q14-Q1*Q5*Q6*Q12*Q14+Q1*Q2*Q11*Q12*Q14-Q3*Q4*Q6*Q8*Q15+Q1*Q5*Q8*Q11*Q15+Q1*Q5*Q6*Q12*Q15-Q1*Q2*Q11*Q12*Q15+Q1*Q3*Q8*Q15^2,
           Q8*Q13-Q7*Q15,
           Q4*Q8*Q11-Q1*Q14^2,
           Q2*Q8*Q9^2*Q11-Q2*Q7*Q9*Q11^2+Q2*Q6* Q9^2*Q12-Q2*Q6*Q9*Q11*Q13+2*Q1*Q8*Q11^2*Q13-Q1*Q8*Q9*Q11*Q14+Q1*Q7*Q11^2*Q14-Q1*Q6* Q9*Q12*Q14-Q1*Q7*Q11^2*Q15+Q1*Q6*Q9*Q12* Q15-Q1*Q6*Q11*Q13*Q15,
           Q4*Q5*Q8*Q11-Q5^2*Q6*Q13+2*Q3*Q5*Q7* Q14-Q3*Q4*Q8*Q14-2*Q2*Q5*Q9*Q14-Q1*Q5*Q12*Q14+2*Q1*Q5*Q14^2+Q2*Q5*Q9*Q15-Q2*Q3* Q13*Q15,
           -Q5*Q8*Q9+Q4*Q8*Q11+Q3*Q7*Q14-Q1*Q12*Q14,
           -Q4*Q6*Q8*Q12-Q5*Q6*Q8*Q13- Q3*Q8^2*Q13+Q2*Q8*Q11*Q13+Q2*Q6*Q12*Q13+ Q5*Q6*Q7*Q14-Q1*Q8*Q12*Q14-Q4*Q6*Q8* Q15+Q2*Q8*Q9*Q15-Q2*Q7*Q11*Q15+Q1*Q8* Q12*Q15+Q1*Q8*Q15^2,
           Q5^3*Q7*Q9+Q4*Q5^2*Q7*Q11+Q4*Q5^2*Q6* Q13-2*Q3*Q5^2*Q7*Q13+Q2*Q3*Q5*Q13^2+2*Q4^2*Q5*Q6*Q14-Q3*Q4*Q5*Q7*Q14-Q3*Q4^2*Q8* Q14+Q1*Q4*Q5*Q12*Q14-2*Q4^2*Q5*Q6*Q15+2*Q3*Q4*Q5*Q7*Q15-Q1*Q4*Q5*Q12*Q15-Q2*Q3* Q4*Q13*Q15-2*Q1*Q5^2*Q13*Q15-Q1*Q4*Q5*Q14* Q15+2*Q1*Q4*Q5*Q15^2,
           -Q5^2*Q9^2*Q12-2*Q5^2*Q9*Q11*Q13-2*Q3*Q5*Q11* Q13^2+Q3*Q4*Q9*Q12*Q14-Q3^2*Q13^2*Q14+Q5^2*Q9^2*Q15+Q4*Q5*Q9*Q11*Q15-Q3*Q4*Q9*Q12* Q15+3*Q3*Q5*Q9*Q13*Q15+Q3*Q4*Q11*Q13*Q15+Q3^2*Q13^2*Q15-Q3*Q4*Q9*Q14*Q15,
           Q2*Q4*Q5*Q9*Q11-Q2*Q4^2*Q11^2+Q1*Q5^2* Q9*Q12+3*Q1*Q4*Q5*Q11*Q12-Q2*Q3*Q4*Q11*Q13-2*Q1*Q5^2*Q11*Q13+Q1*Q3*Q5*Q12*Q13-Q1* Q4*Q5*Q11*Q14-Q1*Q3*Q4*Q12*Q14+Q2*Q3*Q4*Q9*Q15-Q1*Q5^2*Q9*Q15,
           -Q5*Q7*Q11^2-Q5*Q6*Q9*Q12+Q4*Q6*Q11*Q12+Q3*Q8*Q11*Q13-Q3*Q8*Q9*Q15-2*Q4*Q6* Q11*Q15+2*Q2*Q9*Q11*Q15-Q1*Q11*Q14*Q15+2*Q1*Q11*Q15^2,
           -Q5*Q8*Q11-Q5*Q6*Q12+Q5*Q6*Q14+ Q3*Q8*Q14,
           -Q1*Q8*Q11*Q12*Q13+Q2*Q6*Q11* Q13^2-Q2*Q6*Q9*Q13*Q14+Q1*Q8*Q11*Q13*Q14+ Q1*Q6*Q13*Q14^2-Q2*Q8*Q9^2*Q15+Q2*Q7* Q9*Q11*Q15-Q1*Q8*Q9*Q12*Q15+Q1*Q7*Q11*Q12*Q15-Q1*Q8*Q11*Q13*Q15+Q1*Q8*Q9*Q14*Q15-Q1*Q7*Q11*Q14*Q15,
           Q4*Q5*Q7*Q11^2-Q3*Q4*Q7*Q11*Q12+Q1*Q4*Q11*Q12^2+Q1*Q5*Q11*Q12*Q13-Q1*Q5*Q11* Q13*Q14+Q3*Q5*Q7*Q9*Q15-Q1*Q5*Q9*Q12*Q15-Q3^2*Q7*Q13*Q15-Q1*Q5*Q11*Q13*Q15+Q1* Q3*Q12*Q13*Q15+Q1*Q5*Q9*Q14*Q15-Q1*Q3*Q13*Q14*Q15,
           Q3*Q5*Q6*Q8*Q9-Q3^2*Q8^2*Q9-Q1*Q5*Q8*Q11^2+Q1*Q5*Q6*Q11*Q12-Q1*Q2* Q11^2*Q12-Q2*Q3*Q6*Q11*Q13+Q1*Q3*Q6*Q12*Q14+2*Q1*Q3*Q8*Q11*Q15+Q1*Q2*Q11^2*Q15-Q1* Q3*Q6*Q12*Q15-Q1*Q3*Q6*Q14*Q15,
           -Q2*Q5*Q6*Q9-Q2*Q3*Q7*Q11-2*Q1*Q5*Q8*Q11+Q2^2*Q9*Q11+Q1*Q2*Q11*Q12+Q2*Q3* Q6*Q13+Q1*Q3*Q8*Q14-Q1*Q3*Q8*Q15+Q1*Q2*Q11*Q15,
           Q2*Q5*Q9^2*Q11-Q2*Q4*Q9*Q11^2-Q1*Q4*Q11^2*Q12-Q2*Q3*Q9*Q11*Q13+2*Q1*Q5* Q11^2*Q13-Q1*Q3*Q11*Q12*Q13-Q1*Q5*Q9*Q11* Q14+Q1*Q4*Q11^2*Q14+Q2*Q3*Q9^2*Q15+Q1*Q3*Q9*Q12*Q15-Q1*Q3*Q9*Q14*Q15,
           -Q5*Q7*Q11+Q4*Q8*Q11+Q2*Q9*Q12-Q2*Q11* Q13-Q1*Q12*Q14+Q1*Q12*Q15,
           -Q5^2*Q8*Q9-Q4*Q5*Q8*Q11+Q2*Q4*Q11*Q12+Q3*Q5*Q8*Q13-Q3*Q5*Q7*Q14+Q3*Q4*Q8*Q14+Q2*Q5*Q9*Q14-Q1*Q5*Q14^2+Q2*Q5* Q9*Q15+Q1*Q5*Q12*Q15-Q2*Q3*Q13*Q15-Q1*Q5*Q14*Q15,
           -Q5^2*Q9*Q11*Q12+Q4*Q5*Q11^2*Q12+ Q5^2*Q9*Q11*Q14-Q4*Q5*Q11^2*Q14+Q3*Q5*Q9* Q12*Q14-Q3^2*Q13*Q14^2-Q5^2*Q9*Q11*Q15+Q4*Q5*Q11^2*Q15-Q3*Q5*Q9*Q12*Q15+Q3*Q5*Q9* Q14*Q15-Q3*Q4*Q11*Q14*Q15+Q3^2*Q13*Q14*Q15,
           -Q5^2*Q6^2*Q12+Q2*Q5*Q6*Q11*Q12+Q5^2*Q6^2*Q14+Q3*Q5*Q6*Q8*Q14+Q3^2*Q8^2*Q14-Q2^2*Q11^2*Q14-Q2*Q3*Q6*Q12*Q14-Q2*Q3*Q6*Q14^2- Q3^2*Q8^2*Q15-Q2*Q5*Q6*Q11*Q15+Q2^2*Q11^2*Q15+ Q2*Q3*Q6*Q12*Q15,
           -Q8^2*Q9*Q11*Q12-Q7*Q8*Q11^2*Q12+2*Q8^2*Q11^2* Q13-2*Q6*Q8*Q11*Q12*Q13+Q6^2*Q12^2*Q13+Q6*Q7*Q11*Q12*Q14+Q8^2*Q9*Q11*Q15+Q6*Q7*Q11* Q12*Q15-Q6*Q8*Q11*Q13*Q15-Q6*Q7*Q11*Q14*Q15,
           Q8*Q9*Q11*Q12+Q7*Q11^2*Q12+Q6*Q9*Q12^2-2* Q8*Q11^2*Q13+Q6*Q11*Q12*Q13-Q6*Q11*Q13*Q14- Q8*Q9*Q11*Q15,
           -Q4*Q6*Q8^2*Q9*Q11+Q2*Q8^2*Q9^2*Q11-2*Q2*Q7*Q8*Q9*Q11^2+Q2*Q7^2*Q11^3+Q4*Q6^2*Q8* Q9*Q12-Q2*Q6*Q7*Q9*Q11*Q12+Q1*Q8^2*Q9*Q11*Q12+Q2*Q6*Q7*Q11^2*Q13-Q2*Q6^2*Q9*Q12* Q13+Q1*Q6*Q8*Q11*Q12*Q13-Q1*Q7*Q8*Q11^2*Q14+Q4*Q6^2*Q8*Q9*Q15-Q2*Q6*Q8*Q9^2*Q15+2 *Q2*Q6*Q7*Q9*Q11*Q15+Q1*Q8^2*Q9*Q11*Q15-Q1*Q7*Q8*Q11^2*Q15-Q1*Q6*Q8*Q9*Q12* Q15-Q1*Q6*Q8*Q11*Q13*Q15+Q1*Q6*Q7*Q11*Q14*Q15-Q1*Q6*Q8*Q9*Q15^2,
           -Q5^2*Q7*Q8*Q12*Q13+Q4*Q5*Q8^2*Q12*Q13+Q2*Q4*Q8*Q12^2*Q13-2*Q5^2*Q8^2*Q13^2-2*Q2*Q5* Q8*Q12*Q13^2-Q4^2*Q8^2*Q12*Q14+Q4*Q5*Q8^2*Q13* Q14+Q2*Q5*Q8*Q13^2*Q14+Q2^2*Q12*Q13^2*Q14+Q5^2*Q7^2*Q14^2+Q4*Q5*Q7*Q8*Q14^2-Q2*Q4*Q7* Q12*Q14^2+2*Q5^2*Q7*Q8*Q13*Q15+2*Q2*Q5*Q8*Q13^2*Q15+Q2^2*Q12*Q13^2*Q15-Q4*Q5*Q7*Q8*Q14* Q15-3*Q2*Q5*Q7*Q13*Q14*Q15-Q2*Q4*Q8*Q13* Q14*Q15-Q2^2*Q13^2*Q14*Q15+Q2*Q4*Q7*Q14^2* Q15,
           -Q4*Q12^2+Q5*Q13*Q15,
           -Q5*Q8*Q9*Q12-Q5*Q7*Q11*Q12+Q4*Q8*Q11* Q12+Q5*Q8*Q11*Q13+Q2*Q11*Q12*Q13-Q5*Q8* Q9*Q14+Q4*Q8*Q11*Q14-Q4*Q6*Q12*Q14+Q3* Q7*Q12*Q14-Q1*Q12^2*Q14+Q3*Q7*Q14^2-Q1*Q12*Q14^2-Q5*Q7*Q11*Q15-Q2*Q11*Q13*Q15+Q1* Q12*Q14*Q15+Q1*Q14^2*Q15,
           Q4*Q8^2*Q9+Q4*Q7*Q8*Q11-Q5*Q6*Q7*Q13-2*Q1*Q8*Q12*Q13+Q2*Q6*Q13^2-Q3*Q7^2*Q14+ Q2*Q7*Q9*Q14+Q1*Q7*Q12*Q14-Q1*Q7*Q14^2+Q1*Q7*Q12*Q15-Q1*Q7*Q14*Q15,
           -Q2*Q4*Q5*Q7*Q11+Q2^2*Q4*Q9*Q12+Q2^2*Q5*Q9*Q13-Q2^2*Q4*Q11*Q13-Q2^2*Q3*Q13^2-Q1* Q5^2*Q7*Q14+Q1*Q4*Q5*Q8*Q14-Q1*Q2*Q5* Q13*Q14+Q1*Q5^2*Q7*Q15+2*Q1*Q2*Q5*Q13*Q15-Q1*Q2*Q4*Q14*Q15,
           Q4*Q8^2*Q9-Q4*Q7*Q8*Q11+Q5*Q6*Q7*Q13-2*Q1*Q8*Q12*Q13+Q2*Q6*Q13^2-Q3*Q7^2*Q14+ Q2*Q7*Q9*Q14+Q1*Q7*Q12*Q14-Q1*Q7*Q14^2+Q1*Q7*Q12*Q15-Q1*Q7*Q14*Q15,
           Q3*Q8^2*Q9-2*Q2*Q8*Q9*Q11+Q2*Q7*Q11^2-2* Q1*Q8*Q11*Q12-Q5*Q6^2*Q13+Q3*Q6*Q8*Q13-Q3*Q6*Q7*Q14+2*Q1*Q8*Q11*Q14+Q2*Q6*Q9* Q15+Q1*Q6*Q12*Q15-Q1*Q6*Q14*Q15,
           Q10-Q12,
           -Q5*Q8*Q9^2+Q4*Q7*Q11^2-Q4*Q6*Q9*Q12+Q5*Q6*Q9*Q13-2*Q3*Q7*Q11*Q13+2*Q1*Q11*Q12* Q13+Q3*Q6*Q13^2+Q3*Q7*Q9*Q14-Q1*Q9*Q12*Q14-2*Q1*Q11*Q13*Q15+Q1*Q9*Q14*Q15,
           Q4*Q5*Q7*Q9*Q11^2+Q4^2*Q7*Q11^3-2*Q3*Q4* Q7*Q11^2*Q13-Q1*Q4*Q11^2*Q12*Q13+Q3^2*Q7*Q11* Q13^2-Q2*Q3*Q9*Q11*Q13^2+2*Q3*Q4*Q7*Q9* Q11*Q14+Q1*Q4*Q9*Q11*Q12*Q14-Q3^2*Q7*Q9*Q13*Q14+Q2*Q3*Q9^2*Q13*Q14-Q1*Q5*Q9*Q11* Q13*Q14-Q1*Q4*Q11^2*Q13*Q14+Q1*Q3*Q11*Q13^2* Q14-Q1*Q3*Q9*Q13*Q14^2-Q3*Q5*Q7*Q9^2*Q15-Q3*Q4*Q7*Q9*Q11*Q15+Q2*Q3*Q9^2*Q13* Q15+Q1*Q5*Q9*Q11*Q13*Q15+Q1*Q3*Q11*Q13^2*Q15-Q1*Q3*Q9*Q13*Q14*Q15,
           Q5*Q11*Q12-Q3*Q15^2,
           Q5^3*Q7*Q9+Q4*Q5^2*Q8*Q9+2*Q4^2*Q5*Q6* Q12-Q2*Q4*Q5*Q9*Q12-Q2*Q4^2*Q11*Q12+Q4*Q5^2*Q6*Q13-2*Q2*Q5^2*Q9*Q13+Q2*Q3*Q5*Q13^2+Q1*Q4*Q5*Q12*Q14-2*Q4^2*Q5*Q6*Q15+2*Q2* Q4*Q5*Q9*Q15-Q1*Q4*Q5*Q12*Q15-Q2*Q3*Q4*Q13*Q15-2*Q1*Q5^2*Q13*Q15-Q1*Q4*Q5*Q14* Q15+2*Q1*Q4*Q5*Q15^2,
           Q2*Q3*Q4*Q8^2-Q1*Q5^2*Q8^2+Q2^2*Q5*Q7* Q11+Q2^2*Q4*Q8*Q11-Q1*Q2*Q5*Q8*Q12-Q2^2* Q3*Q8*Q13-Q2^3*Q11*Q13+Q1*Q2^2*Q12*Q15,
           Q5*Q7*Q8*Q11-Q4*Q8^2*Q11-Q2*Q8*Q11*Q13+Q2*Q6*Q12*Q13-Q5*Q6*Q7*Q14+Q3*Q7*Q8*Q14-Q2*Q8*Q9*Q14+Q1*Q8*Q14^2-Q2*Q8* Q9*Q15+Q2*Q7*Q11*Q15-Q1*Q8*Q12*Q15+Q1*Q8*Q14*Q15,
           Q4*Q5^2*Q6^2*Q11-Q3*Q5^2*Q6*Q7*Q11-2*Q2*Q4*Q5*Q6*Q11^2+Q2*Q3*Q4*Q8*Q11^2 +Q2^2*Q4*Q11^3-Q3*Q4*Q5*Q6^2*Q12+Q3^2*Q5*Q6*Q7*Q12+2*Q2*Q3*Q4*Q6*Q11*Q12+Q1*Q5^2* Q6*Q11*Q12-Q1*Q3*Q5*Q8*Q11*Q12-Q1*Q2* Q5*Q11^2*Q12-Q1*Q3*Q5*Q6*Q12^2+Q3^2*Q5*Q6*Q7*Q14-Q3^2*Q4*Q6*Q8*Q14-Q2*Q3*Q4* Q6*Q11*Q14+Q1*Q5^2*Q6*Q11*Q14+Q1*Q3*Q5*Q8*Q11*Q14-Q1*Q3*Q5*Q6*Q12*Q14-Q1*Q2* Q5*Q11^2*Q15+Q1*Q2*Q3*Q11*Q12*Q15,
           2*Q5^2*Q8^2*Q11-Q5^2*Q6*Q8*Q12-Q3*Q5*Q8^2* Q12-2*Q2*Q5*Q8*Q11*Q12+Q2^2*Q11*Q12^2+Q5^2*Q6*Q8*Q14-Q2*Q5*Q8*Q11*Q14+Q2*Q3*Q8* Q12*Q14+Q2*Q3*Q8*Q12*Q15-Q2*Q3*Q8*Q14*Q15,
           Q2*Q3*Q4*Q8*Q11-Q1*Q5^2*Q8*Q11-Q2^2*Q4*Q11^2+Q1*Q5^2*Q6*Q12+2*Q1*Q2*Q5*Q11*Q12- Q2*Q3*Q5*Q6*Q13-Q1*Q5^2*Q6*Q14+Q1*Q3* Q5*Q8*Q14-Q1*Q2*Q3*Q12*Q14-Q1*Q2*Q3* Q12*Q15+Q1*Q2*Q3*Q14*Q15,
           -Q6*Q12^2+Q8*Q11*Q15,
           -Q4*Q8*Q11*Q12+Q5*Q8*Q11*Q13+Q5*Q6*Q12* Q13-Q3*Q8*Q12*Q13+Q2*Q11*Q12*Q13-Q5*Q6* Q13*Q14-Q3*Q8*Q13*Q14-Q4*Q8*Q11*Q15+Q3* Q7*Q12*Q15-Q2*Q9*Q12*Q15-Q1*Q12^2*Q15+Q2* Q11*Q13*Q15+Q1*Q12*Q14*Q15+Q3*Q7*Q15^2-Q1*Q12*Q15^2+Q1*Q14*Q15^2,
           Q4*Q8^3*Q9^2+Q4*Q7*Q8^2*Q9*Q11+2*Q4*Q6* Q7*Q8*Q9*Q12-2*Q4*Q6*Q8^2*Q9*Q13-Q4*Q6^2* Q7*Q12*Q13+Q3*Q6*Q7^2*Q12*Q13-Q1*Q8^2* Q9*Q12*Q13-Q1*Q7*Q8*Q11*Q12*Q13-Q1*Q6*Q7*Q12^2*Q13+Q4*Q6^2*Q8*Q13^2-Q3*Q6*Q7*Q8* Q13^2+Q1*Q6*Q8*Q12*Q13^2-Q4*Q6*Q7*Q8*Q9* Q14-Q4*Q6*Q7^2*Q11*Q14+Q3*Q6*Q7^2*Q13*Q14+Q1*Q7*Q8*Q11*Q13*Q14-Q1*Q6*Q7*Q12*Q13*Q14+Q1*Q6*Q8*Q13^2*Q14+Q1*Q7*Q8*Q9* Q12*Q15-Q1*Q8^2*Q9*Q13*Q15,
           -Q5*Q8*Q9+Q4*Q8*Q11-Q3*Q8*Q13+Q3*Q7* Q14-Q1*Q12*Q14+Q1*Q14*Q15,
           -Q4*Q5*Q6*Q7*Q13+Q4^2*Q6*Q8*Q13-Q1*Q5*Q7*Q12*Q13-Q2*Q4*Q6*Q13^2+2*Q1*Q5*Q8* Q13^2-Q1*Q2*Q12*Q13^2+Q4^2*Q6*Q7*Q14+Q1*Q4* Q7*Q12*Q14-Q1*Q4*Q8*Q13*Q15+Q1*Q2*Q13^2* Q15-Q1*Q4*Q7*Q14*Q15,
           2*Q4*Q5*Q8^2*Q12*Q13-2*Q5^2*Q8^2*Q13^2+2*Q2*Q5*Q8*Q12*Q13^2-Q4*Q5*Q7*Q8*Q12*Q14+Q5^2* Q7*Q8*Q13*Q14-Q2*Q5*Q7*Q12*Q13*Q14-3*Q2*Q4*Q8*Q12*Q13*Q14+Q2*Q5*Q8*Q13^2*Q14-Q2^2* Q12*Q13^2*Q14+Q4*Q5*Q7*Q8*Q14^2+Q4^2*Q8^2*Q14^2+Q2*Q4*Q7*Q12*Q14^2+Q5^2*Q7*Q8*Q13*Q15- Q4*Q5*Q8^2*Q13*Q15-2*Q2*Q5*Q8*Q13^2*Q15+Q2^2*Q12*Q13^2*Q15-Q5^2*Q7^2*Q14*Q15+Q2^2*Q13^2*Q14*Q15-Q2*Q4*Q7*Q14^2*Q15+Q2*Q5*Q7*Q13* Q15^2,
           -Q3*Q5^2*Q7*Q9*Q11+Q2*Q5^2*Q9^2*Q11-2*Q2*Q4*Q5*Q9*Q11^2+Q2*Q4^2*Q11^3+Q3^2*Q5*Q7* Q9*Q12-Q2*Q3*Q5*Q9^2*Q12+2*Q2*Q3*Q4*Q9* Q11*Q12+Q1*Q5^2*Q9*Q11*Q12-Q1*Q4*Q5*Q11^2*Q12-Q1*Q3*Q5*Q9*Q12^2+Q2*Q3*Q4*Q11^2* Q13-Q1*Q3*Q5*Q11*Q12*Q13-Q1*Q4*Q5*Q11^2*Q14+Q1*Q3*Q4*Q11*Q12*Q14+Q3^2*Q5*Q7*Q9* Q15-Q2*Q3*Q4*Q9*Q11*Q15+Q1*Q5^2*Q9*Q11* Q15-Q1*Q3*Q5*Q9*Q12*Q15-Q2*Q3^2*Q9*Q13* Q15+Q1*Q3*Q5*Q11*Q13*Q15,
           Q4*Q5^2*Q6^2-Q3*Q4*Q5*Q6*Q8-Q3^2*Q5* Q7*Q8+Q3^2*Q4*Q8^2-Q2^2*Q4*Q11^2+Q1*Q5^2*Q6*Q12-Q1*Q3*Q5*Q8*Q12+2*Q1*Q2*Q5*Q11* Q12+Q2*Q3^2*Q8*Q13+Q2*Q3*Q4*Q6*Q14-Q1*Q5^2*Q6*Q14-Q1*Q3*Q5*Q8*Q14-Q2*Q3*Q4* Q6*Q15-Q1*Q5^2*Q6*Q15+Q1*Q3*Q5*Q8*Q15-Q1*Q2*Q3*Q12*Q15+Q1*Q2*Q3*Q15^2,
           -Q5*Q8*Q9*Q11-Q5*Q7*Q11^2+Q5*Q6*Q9*Q12-Q4*Q6*Q11*Q12+Q3*Q7*Q11*Q12-Q1*Q11*Q12^2+Q5*Q6*Q11*Q13+Q3*Q7*Q11*Q14-Q1*Q11* Q12*Q14-Q3*Q6*Q13*Q14+Q3*Q8*Q9*Q15+Q1*Q11*Q14*Q15,
           -Q5^2*Q9^2*Q12+Q4*Q5*Q9*Q11*Q12- Q4*Q5*Q9*Q11*Q14+Q4^2*Q11^2*Q14+Q3*Q4*Q9* Q12*Q14-Q3^2*Q13^2*Q14+Q5^2*Q9^2*Q15-Q4^2*Q11^2* Q15-Q3*Q4*Q9*Q12*Q15+Q3*Q5*Q9*Q13*Q15+ Q3^2*Q13^2*Q15-Q3*Q4*Q9*Q15^2,
           Q5*Q8^2*Q9^2+Q5*Q7*Q8*Q9*Q11+Q5*Q7^2* Q11^2-Q5*Q6*Q7*Q9*Q12+Q1*Q8*Q11*Q12*Q13-Q5*Q6^2*Q13^2-Q5*Q6*Q7*Q9*Q14+Q1*Q8*Q11* Q13*Q14+Q1*Q6*Q12*Q13*Q14-Q5*Q6*Q7*Q9*Q15-Q1*Q8*Q9*Q12*Q15-Q1*Q7*Q11*Q14*Q15,
           -Q2*Q3*Q13+Q1*Q5*Q15,
           Q2*Q3*Q13-Q1*Q5*Q15,
           2*Q8^2*Q11^2*Q13-Q8^2*Q9*Q11*Q14-Q7*Q8*Q11^2* Q14+Q6*Q8*Q9*Q12*Q14-2*Q6*Q8*Q11*Q13*Q14+ Q6^2*Q13*Q14^2+Q7*Q8*Q11^2*Q15-Q6*Q8*Q9*Q12*Q15-Q6*Q8*Q11*Q13*Q15+Q6*Q8*Q9*Q14* Q15,
           Q5*Q8*Q9+Q4*Q6*Q12-Q3*Q8*Q13-Q1* Q12*Q15,
           Q5*Q6*Q7-Q4*Q6*Q8-Q2*Q7*Q11+ Q2*Q6*Q13,
           -Q6*Q13*Q14+Q7*Q11*Q15,
           Q5*Q11-Q3*Q12,
           -Q5*Q11+Q3*Q12,
           2*Q5*Q11^2*Q13^2-Q4*Q9*Q11*Q12*Q14-Q5*Q9* Q11*Q13*Q14+Q3*Q11*Q13^2*Q14+Q4*Q9*Q11*Q12*Q15-2*Q5*Q9*Q11*Q13*Q15-Q4*Q11^2*Q13*Q15-Q3* Q11*Q13^2*Q15+Q4*Q9*Q11*Q14*Q15+Q5*Q9^2*Q15^ 2,
           Q3*Q8*Q13-Q1*Q14^2,
           Q4*Q5*Q6^2*Q8-Q3*Q4*Q6*Q8^2-Q2*Q4* Q6*Q8*Q11+2*Q1*Q5*Q8^2*Q11+Q2*Q4*Q6^2*Q12-Q1*Q3*Q8^2*Q14-Q1*Q2*Q8*Q11*Q14+Q1*Q2* Q6*Q12*Q14-Q1*Q5*Q6*Q8*Q15+Q1*Q3*Q8^2* Q15-Q1*Q2*Q6*Q12*Q15,
           -Q3*Q6*Q8*Q9*Q13+Q3*Q6*Q7*Q11*Q13-2*Q1*Q8*Q11^2*Q13-Q1*Q6*Q11*Q12*Q13-Q3*Q6^2* Q13^2+Q1*Q8*Q9*Q11*Q14+Q1*Q7*Q11^2*Q14-Q1*Q6*Q9*Q12*Q14+3*Q1*Q6*Q11*Q13*Q14+Q3*Q6* Q7*Q9*Q15-Q1*Q7*Q11^2*Q15,
           -Q4*Q9*Q12+Q3*Q13^2,
           -Q1*Q11*Q12+Q3*Q6*Q13,
           -Q5*Q6*Q7*Q9+Q4*Q6*Q8*Q9-Q2*Q8*Q9^ 2+Q2*Q7*Q9*Q11-2*Q1*Q8*Q11*Q13+Q1*Q8* Q9*Q14-Q1*Q7*Q11*Q14+Q1*Q7*Q11*Q15+Q1*Q6*Q13*Q15,
           Q5*Q6*Q13-Q3*Q8*Q13-Q2*Q11*Q13+Q1*Q15^2,
           -Q5*Q6*Q7+Q4*Q6*Q8+Q2*Q6* Q13-Q1*Q8*Q15,
           Q5*Q6*Q7+Q4*Q6*Q8-Q2* Q6*Q13-Q1*Q8*Q15,
           Q4*Q8*Q12^2-Q5*Q8*Q12*Q13+Q2*Q12^2*Q13+ Q5*Q7*Q12*Q14-Q5*Q8*Q13*Q14+Q5*Q7*Q12*Q15-Q5*Q8*Q13*Q15-Q5*Q7*Q14*Q15,
           Q3*Q8*Q9-Q1*Q11*Q14,
           -Q3*Q5*Q7*Q9*Q11-Q3*Q4*Q7*Q11^2+Q1*Q4*Q11^2*Q12+Q3^2*Q7*Q11*Q13+2*Q1*Q5*Q11^2*Q13- Q1*Q3*Q11*Q12*Q13-Q1*Q5*Q9*Q11*Q14-Q1*Q4*Q11^2*Q14+Q3^2*Q7*Q9*Q15-Q1*Q3*Q9*Q12* Q15+Q1*Q3*Q9*Q14*Q15,
           Q5^3*Q6*Q7-2*Q3*Q5^2*Q7*Q8+Q3*Q4*Q5* Q8^2+Q2*Q5^2*Q8*Q9+Q2*Q5^2*Q7*Q11+2*Q2*Q3* Q5*Q7*Q14-Q2*Q3*Q4*Q8*Q14-2*Q1*Q5^2* Q8*Q14-2*Q2^2*Q5*Q9*Q14-Q1*Q2*Q5*Q12*Q14+2*Q1*Q2*Q5*Q14^2-Q2*Q3*Q5*Q7*Q15+2*Q2^2*Q5*Q9*Q15+Q1*Q2*Q5*Q12*Q15-Q2^2*Q3*Q13* Q15-Q1*Q2*Q5*Q14*Q15,
           Q10-Q15,
           -Q4*Q8*Q9*Q12+Q5*Q8*Q9*Q13-Q5*Q7*Q11* Q13+Q4*Q6*Q12*Q13-Q5*Q6*Q13^2+Q4*Q7*Q11*Q14+Q1*Q12*Q13*Q14+Q5*Q7*Q9*Q15+Q4*Q6*Q13*Q15-Q2*Q9*Q13*Q15-Q1*Q12*Q13*Q15-Q1*Q13*Q15^2,
           Q5*Q7*Q11-Q4*Q8*Q11-Q2*Q11*Q13+Q2*Q9*Q15+Q1*Q12*Q15-Q1*Q14*Q15,
           Q2*Q4*Q6-Q1*Q5*Q8,
           -Q5^2*Q8*Q9^2-Q4*Q5*Q7*Q11^2+2*Q4*Q5*Q6*Q9*Q12-Q4^2*Q6*Q11*Q12+Q3*Q4*Q7*Q11*Q12- Q1*Q4*Q11*Q12^2+Q4*Q5*Q6*Q11*Q13-Q1*Q5* Q11*Q12*Q13+Q3*Q4*Q7*Q11*Q14-Q1*Q4*Q11* Q12*Q14-Q3*Q4*Q6*Q13*Q14+Q1*Q5*Q11*Q13* Q14+Q1*Q5*Q9*Q12*Q15,
           -Q5*Q11*Q13+Q4*Q11*Q14+Q3*Q13*Q14-Q3*Q13* Q15,
           -Q2*Q8^2*Q9^2+Q2*Q7*Q8*Q9*Q11-Q1* Q8^2*Q9*Q12-Q5*Q6^2*Q7*Q13+2*Q2*Q6*Q8*Q9* Q13-Q2*Q6*Q7*Q11*Q13+Q1*Q6*Q8*Q12*Q13- Q2*Q6^2*Q13^2+Q1*Q6*Q8*Q13*Q14+Q2*Q6*Q7*Q9*Q15-Q1*Q8^2*Q9*Q15,
           -Q5^2*Q6^2*Q7+Q3^2*Q7*Q8^2-Q2*Q3*Q7*Q8*Q11-Q2^2*Q8*Q9*Q11+Q2^2*Q7*Q11^2-Q2*Q3*Q6*Q7*Q12-Q1*Q3*Q8^2*Q12+Q1*Q2*Q8*Q11* Q12+Q1*Q2*Q6*Q12^2+Q2^2*Q6*Q11*Q13+2*Q1*Q5*Q6*Q8*Q14+Q1*Q3*Q8^2*Q14-Q1*Q2*Q8*Q11* Q14-Q1*Q2*Q6*Q12*Q14+Q2*Q3*Q6*Q7*Q15-Q1*Q3*Q8^2*Q15-Q1*Q2*Q8*Q11*Q15,
           Q4*Q5*Q6*Q9-Q3*Q4*Q8*Q9-Q2*Q5*Q9^2 +Q2*Q3*Q9*Q13-2*Q1*Q5*Q11*Q13+Q1*Q5*Q9*Q14-Q1*Q3*Q13*Q14+Q1*Q4*Q11*Q15+Q1*Q3*Q13*Q15,
           Q2*Q5^2*Q9^2-Q2*Q4^2*Q11^2+Q1*Q5^2* Q9*Q12+2*Q1*Q4*Q5*Q11*Q12-Q3^2*Q5*Q7* Q13+Q3^2*Q4*Q8*Q13-Q2*Q3*Q5*Q9*Q13-Q1*Q3*Q5*Q12*Q13+Q2*Q3^2*Q13^2-Q2*Q3*Q4*Q9* Q14-Q1*Q5^2*Q9*Q14-Q1*Q3*Q4*Q12*Q14+Q1*Q3*Q5*Q13*Q14+Q1*Q3*Q4*Q14^2+Q2*Q3*Q4* Q9*Q15-Q1*Q5^2*Q9*Q15-Q1*Q3*Q5*Q13*Q15,
           Q4*Q5*Q9*Q11*Q12^2-2*Q5^2*Q9*Q11*Q12*Q13+ Q4*Q5*Q11^2*Q12*Q13-2*Q5^2*Q11^2*Q13^2-Q3*Q5*Q11*Q12*Q13^2+Q5^2*Q9^2*Q12*Q14-Q4^2*Q11^2*Q12*Q14+ Q5^2*Q9*Q11*Q13*Q14+Q4*Q5*Q11^2*Q13*Q14-Q3*Q4*Q9*Q12*Q14^2+Q3*Q4*Q11*Q13*Q14^2+Q3^2* Q13^2*Q14^2+Q5^2*Q9^2*Q12*Q15+2*Q5^2*Q9*Q11*Q13*Q15+2*Q3*Q5*Q11*Q13^2*Q15-Q5^2*Q9^2*Q14*Q15-Q4*Q5*Q9*Q11*Q14*Q15-3*Q3*Q5*Q9*Q13*Q14* Q15-Q3*Q4*Q11*Q13*Q14*Q15+Q3*Q4*Q9*Q14^2*Q15,
           -Q4*Q8*Q11*Q12+Q5*Q8*Q11*Q13-Q2*Q11* Q12*Q13+Q5*Q7*Q11*Q14+Q4*Q8*Q11*Q14-Q5*Q6*Q13*Q14-Q2*Q11*Q13*Q14+Q5*Q7*Q11*Q15- Q5*Q6*Q13*Q15-Q3*Q7*Q14*Q15+Q2*Q9*Q14*Q15+Q1*Q12*Q14*Q15-Q1*Q14^2*Q15+Q2*Q9*Q15^2+ Q1*Q12*Q15^2-Q1*Q14*Q15^2,
           Q5*Q9*Q12+Q4*Q11*Q12-Q5*Q11*Q13-Q5*Q9* Q15,
           -Q5*Q6*Q8*Q9-Q3*Q8^2*Q9+Q5*Q6*Q7*Q11-Q4*Q6*Q8*Q11+2*Q2*Q8*Q9*Q11-Q2* Q7*Q11^2-Q4*Q6^2*Q12+Q2*Q6*Q9*Q12-Q3*Q6*Q8*Q13-Q1*Q6*Q12*Q14+2*Q1*Q8*Q11*Q15+Q1* Q6*Q12*Q15,
           Q5^2*Q6*Q7*Q8-2*Q3*Q5*Q7*Q8^2 +Q3*Q4*Q8^3+Q2*Q5*Q8^2*Q9+Q2*Q3*Q8^2*Q13-Q2*Q5*Q6*Q7*Q14+2*Q2*Q3*Q7*Q8*Q14-2* Q1*Q5*Q8^2*Q14-2*Q2^2*Q8*Q9*Q14-Q1*Q2*Q8*Q12*Q14+2*Q1*Q2*Q8*Q14^2-Q2*Q3*Q7*Q8* Q15+2*Q2^2*Q8*Q9*Q15-Q2^2*Q7*Q11*Q15+Q1*Q2* Q8*Q12*Q15-Q1*Q2*Q8*Q14*Q15,
           -Q5*Q7*Q9*Q11-Q4*Q8*Q9*Q11-Q4*Q7*Q11^ 2+Q4*Q6*Q9*Q12-Q2*Q9^2*Q12+Q3*Q8*Q9* Q13+2*Q4*Q6*Q11*Q13-Q2*Q9*Q11*Q13-Q3*Q6*Q13^2+Q1*Q9*Q12*Q14+2*Q1*Q11*Q13*Q14-Q1*Q9* Q12*Q15,
           Q7*Q11*Q14-Q8*Q9*Q15,
           -Q4*Q12+Q5*Q13,
           Q5*Q8*Q9-Q4*Q8*Q11-Q3*Q8*Q13+Q3*Q7* Q15-Q1*Q12*Q15+Q1*Q14*Q15,
           Q8^2*Q9*Q11*Q12*Q13+Q7*Q8*Q11^2*Q12*Q13+ Q6*Q8*Q9*Q12^2*Q13-2*Q8^2*Q11^2*Q13^2+Q6^2*Q12^2*Q13^2-Q7*Q8*Q9*Q11*Q12*Q14-Q7^2*Q11^2*Q12*Q14+ Q6*Q7*Q9*Q12^2*Q14+2*Q7*Q8*Q11^2*Q13*Q14-Q6*Q8*Q9*Q12*Q13*Q14-3*Q6*Q7*Q11*Q12*Q13* Q14+2*Q6*Q8*Q11*Q13^2*Q14-Q8^2*Q9^2*Q12*Q15+Q7^2* Q11^2*Q12*Q15-Q6*Q7*Q9*Q12^2*Q15+Q8^2*Q9* Q11*Q13*Q15-2*Q7*Q8*Q11^2*Q13*Q15-Q6*Q8*Q11*Q13^2*Q15+Q7^2*Q11^2*Q14*Q15+Q7*Q8*Q9*Q11*Q15^2,
           Q2*Q4*Q8*Q9-Q2*Q4*Q7*Q11+Q1*Q5*Q7*Q12+Q1*Q4*Q8*Q12+Q2*Q3*Q7*Q13-2*Q1* Q5*Q8*Q13-Q2^2*Q9*Q13-Q1*Q4*Q8*Q14+Q1*Q2*Q13*Q14,
           -Q4*Q5*Q6*Q7+Q3*Q5*Q7^2+Q2* Q4*Q7*Q11+Q1*Q4*Q8*Q12-Q2*Q3*Q7*Q13-2* Q1*Q5*Q8*Q13+Q1*Q5*Q7*Q14-Q1*Q4*Q8* Q14+Q1*Q5*Q7*Q15,
           Q5^2*Q6*Q9*Q11+Q3*Q5*Q7*Q11^2-2*Q2*Q5* Q9*Q11^2+Q2*Q4*Q11^3-Q3*Q5*Q6*Q9*Q12-2*Q3^2* Q7*Q11*Q12+2*Q2*Q3*Q9*Q11*Q12-2*Q1*Q5* Q11^2*Q12+2*Q1*Q3*Q11*Q12^2+Q2*Q3*Q11^2*Q13-Q1* Q3*Q11*Q12*Q14-Q3^2*Q8*Q9*Q15+2*Q3^2*Q7*Q11*Q15-Q2*Q3*Q9*Q11*Q15-Q1*Q3*Q11*Q12* Q15+Q1*Q3*Q11*Q14*Q15,
           Q3*Q5*Q7*Q9-Q3*Q4*Q8*Q9-Q2*Q5*Q9^2 +Q2*Q4*Q9*Q11+Q1*Q4*Q11*Q12-2*Q1*Q5*Q11*Q13+Q1*Q3*Q12*Q13+Q1*Q5*Q9*Q14-Q1*Q4*Q11*Q14,
           -2*Q5^2*Q8^2*Q11^2-2*Q5^2*Q6*Q8*Q11* Q12-Q3*Q5*Q8^2*Q11*Q12+Q2*Q5*Q8*Q11^2*Q12+Q2*Q5*Q6*Q11*Q12^2+2*Q5^2*Q6*Q8*Q11*Q14+2*Q3*Q5*Q8^2*Q11*Q14+Q5^2*Q6^2*Q12*Q14+Q5^2*Q6* Q8*Q11*Q15+Q2*Q5*Q8*Q11^2*Q15+Q5^2*Q6^2*Q12*Q15-Q2^2*Q11^2*Q12*Q15-Q5^2*Q6^2*Q14*Q15-3*Q3* Q5*Q6*Q8*Q14*Q15-Q2*Q5*Q6*Q11*Q14*Q15- Q2*Q3*Q8*Q11*Q14*Q15+Q3^2*Q8^2*Q15^2+Q2*Q3*Q8*Q11*Q15^2-Q2*Q3*Q6*Q12*Q15^2+Q2*Q3* Q6*Q14*Q15^2,
           Q5*Q11-Q3*Q15,
           Q5*Q6*Q13-Q1*Q14^2,
           -Q2*Q4*Q7*Q8*Q11+Q1*Q4*Q8^2*Q12+Q2^2*Q8*Q9*Q13-Q2^2*Q7*Q11*Q13+2*Q1*Q2*Q8*Q12* Q13-Q2^2*Q6*Q13^2+Q1*Q5*Q7*Q8*Q14-Q1*Q4* Q8^2*Q14-Q1*Q2*Q7*Q12*Q14-Q1*Q2*Q8*Q13* Q14+Q2^2*Q7*Q9*Q15,
           -Q5*Q6*Q7*Q8*Q9+Q1*Q8^2*Q9*Q12-Q5*Q6^2*Q7*Q13+Q4*Q6^2*Q8*Q13+2*Q1*Q6*Q8*Q12* Q13-Q2*Q6^2*Q13^2+Q4*Q6^2*Q7*Q14-Q1*Q8^2*Q9*Q15+Q1*Q7*Q8*Q11*Q15-Q1*Q6*Q7*Q12*Q15-Q1*Q6*Q8*Q13*Q15,
           -Q3*Q4*Q6*Q8+Q3^2*Q7*Q8-Q2*Q3*Q7*Q11-2*Q1*Q5*Q8*Q11+Q1*Q5*Q6*Q12+Q2*Q3* Q6*Q13-Q1*Q5*Q6*Q14+Q1*Q3*Q8*Q14+Q1*Q3*Q8*Q15,
           Q5^2*Q6^2*Q12-Q3^2*Q8^2*Q12+Q2*Q5* Q6*Q11*Q12+Q2^2*Q11^2*Q12-Q2*Q3*Q6*Q12^2-Q5^2*Q6^2*Q14+Q3^2*Q8^2*Q14-Q2*Q3*Q8*Q11*Q14+ Q2*Q3*Q8*Q11*Q15-Q2^2*Q11^2*Q15-Q2*Q3*Q6*Q12*Q15+Q2*Q3*Q6*Q14*Q15,
           Q5^2*Q6^2*Q7-Q4*Q5*Q6^2*Q8-Q3*Q5*Q6* Q7*Q8+Q3^2*Q7*Q8^2-Q2^2*Q7*Q11^2-Q2*Q3*Q6* Q7*Q12+Q1*Q5*Q6*Q8*Q12-Q1*Q3*Q8^2*Q12+Q1*Q2*Q6*Q12^2+Q2*Q5*Q6^2*Q13+Q2*Q3* Q6*Q7*Q14-Q1*Q5*Q6*Q8*Q14-Q1*Q3*Q8^2*Q14-Q1*Q5*Q6*Q8*Q15+Q1*Q3*Q8^2*Q15+2*Q1* Q2*Q8*Q11*Q15-Q1*Q2*Q6*Q12*Q15,
           Q2*Q11*Q13-Q1*Q15^2,
           -Q5^2*Q6*Q7+2*Q4*Q5*Q6*Q8-Q3*Q4*Q8^2-Q2*Q5*Q8*Q9+Q2*Q5*Q7*Q11-Q2*Q4*Q8* Q11+Q2*Q4*Q6*Q12-Q2^2*Q9*Q12-Q2*Q3*Q8*Q13+2*Q1*Q5*Q8*Q14+Q1*Q2*Q12*Q14-Q1*Q2* Q12*Q15,
           -Q8*Q11*Q12*Q13-Q8*Q11*Q13*Q14-Q6*Q12*Q13*Q14-Q8*Q11*Q13*Q15+Q6*Q12*Q13*Q15+ Q6*Q13*Q14*Q15+Q8*Q9*Q15^2+Q7*Q11*Q15^2,
           -Q3*Q5*Q7+Q3*Q4*Q8-Q2*Q4*Q11+Q2*Q3* Q13,
           -Q3*Q5*Q7*Q9+Q3*Q4*Q8*Q9+Q1*Q4*Q11*Q12-Q3*Q4*Q6*Q13+Q3^2*Q7*Q13-2*Q1* Q5*Q11*Q13-Q1*Q4*Q11*Q14+Q1*Q3*Q13*Q14+Q1*Q3*Q13*Q15,
           Q3*Q4*Q5*Q6*Q8-Q3^2*Q4*Q8^2-Q2*Q3*Q4*Q8*Q11-2*Q1*Q5^2*Q8*Q11+Q2* Q3*Q4*Q6*Q12-Q1*Q5^2*Q6*Q12+Q1*Q5^2*Q6*Q14+3*Q1*Q3*Q5*Q8*Q14+Q1*Q2*Q5*Q11* Q14-Q1*Q3*Q5*Q8*Q15-Q1*Q2*Q3*Q14*Q15,
           -Q8*Q11*Q12*Q13+Q8*Q9*Q12*Q14+Q8*Q11*Q13* Q14+Q6*Q12*Q13*Q14-Q7*Q11*Q14^2-Q6*Q13*Q14^2-Q8*Q9*Q12*Q15-Q8*Q11*Q13*Q15+Q8*Q9* Q14*Q15+Q7*Q11*Q14*Q15,
           Q2*Q5*Q8*Q9*Q11-Q2*Q4*Q8*Q11^2-Q2*Q5*Q6*Q9*Q12+Q2*Q4*Q6*Q11*Q12+Q1*Q5*Q8*Q11*Q12-Q3^2*Q8^2*Q13+2*Q2*Q3*Q8*Q9*Q15+ Q2*Q4*Q6*Q11*Q15-Q1*Q5*Q8*Q11*Q15-Q2^2*Q9*Q11*Q15-Q1*Q2*Q11*Q12*Q15+Q1*Q3*Q8* Q14*Q15-Q1*Q2*Q11*Q15^2,
           Q3*Q6*Q13-Q1*Q11*Q14,
           -Q3*Q6*Q13+Q1*Q11*Q14,
           2*Q5*Q8^2*Q11^2-2*Q5*Q6*Q8*Q11*Q14-Q3*Q8^2* Q11*Q14-Q2*Q8*Q11^2*Q14+Q2*Q6*Q11*Q12*Q14+ Q5*Q6^2*Q14^2-Q5*Q6*Q8*Q11*Q15+Q3*Q8^2*Q11*Q15-Q2*Q6*Q11*Q12*Q15+Q2*Q6*Q11*Q14* Q15,
           Q4*Q5*Q6^2*Q11-Q3*Q4*Q6*Q8*Q11-Q2* Q4*Q6*Q11^2+2*Q1*Q5*Q8*Q11^2-Q1*Q3*Q8*Q11* Q12-Q1*Q2*Q11^2*Q12+Q3*Q4*Q6^2*Q14+Q1*Q3*Q6*Q12*Q14-Q1*Q5*Q6*Q11*Q15+Q1*Q2*Q11^2*Q15-Q1*Q3*Q6*Q14*Q15,
           Q5*Q6*Q9-Q1*Q11*Q14,
           -Q3*Q4*Q5*Q8*Q9+Q4^2*Q5*Q6*Q11-Q3*Q4^2*Q8*Q11-Q2*Q4^2*Q11^2+Q1*Q5^2*Q9*Q12+2*Q1* Q4*Q5*Q11*Q12+Q3*Q4^2*Q6*Q14-Q1*Q5^2*Q9*Q15-Q1*Q4*Q5*Q11*Q15-Q1*Q3*Q4*Q12*Q15+Q1*Q3*Q5*Q13*Q15,
           Q4*Q8^2*Q9*Q11-Q1*Q8*Q11*Q12*Q13+Q4*Q6*Q7*Q11*Q14+Q1*Q7*Q11*Q12*Q14-Q4*Q6^2* Q13*Q14-Q1*Q8*Q11*Q13*Q14-Q1*Q6*Q12*Q13*Q14-Q4*Q6*Q8*Q9*Q15+Q1*Q8*Q11*Q13*Q15-Q1*Q7*Q11*Q14*Q15+Q1*Q6*Q13*Q14*Q15+Q1*Q8*Q9*Q15^2,
           2*Q5^2*Q8*Q11^2-Q3*Q5*Q8*Q11*Q12+ Q2*Q5*Q11^2*Q12-Q3*Q5*Q6*Q12*Q14-Q5^2*Q6*Q11*Q15-2*Q3*Q5*Q8*Q11*Q15-Q2*Q5*Q11^2* Q15+Q3*Q5*Q6*Q12*Q15+Q3*Q5*Q6*Q14*Q15+Q3^2*Q8*Q15^2,
           Q2*Q3*Q5*Q7*Q13-Q2*Q3*Q4*Q8*Q13-2*Q1*Q5^2*Q8*Q13-Q1*Q2*Q5*Q12*Q13- Q2^2*Q3*Q13^2+Q2*Q3*Q4*Q7*Q14-Q1*Q5^2*Q7*Q14+Q1*Q5^2*Q7*Q15+Q1*Q4*Q5*Q8*Q15- Q1*Q2*Q4*Q12*Q15+3*Q1*Q2*Q5*Q13*Q15,
           Q2*Q8^2*Q9^2-Q2*Q7^2*Q11^2+Q2*Q6*Q7*Q9* Q12-Q1*Q8^2*Q9*Q12+Q5*Q6^2*Q7*Q13-Q4*Q6^2* Q8*Q13-Q2*Q6*Q8*Q9*Q13-Q1*Q6*Q8*Q12* Q13+Q2*Q6^2*Q13^2-Q2*Q6*Q7*Q9*Q14-Q1*Q8^2*Q9*Q14+Q1*Q6*Q8*Q13*Q14+Q1*Q6*Q7* Q14^2+Q1*Q8^2*Q9*Q15+2*Q1*Q7*Q8*Q11*Q15-Q1*Q6*Q8*Q13*Q15-Q1*Q6*Q7*Q14*Q15,
           Q8*Q9^2*Q12^2-2*Q8*Q9*Q11*Q12*Q13-Q7*Q11^2* Q12*Q13+2*Q8*Q11^2*Q13^2-Q6*Q11*Q12*Q13^2+Q7*Q9*Q11*Q12*Q14-Q8*Q9*Q11*Q13*Q14+Q6*Q11*Q13^2*Q14+Q7*Q9*Q11*Q12*Q15-Q7*Q9*Q11*Q14* Q15,
           -Q8^2*Q9*Q11*Q12+Q7*Q8*Q11^2*Q12+Q8^2*Q9* Q11*Q14-Q7*Q8*Q11^2*Q14+Q6*Q8*Q9*Q12*Q14- Q6*Q7*Q11*Q12*Q14+Q6^2*Q12*Q13*Q14-Q6^2* Q13*Q14^2-Q8^2*Q9*Q11*Q15+Q7*Q8*Q11^2*Q15-Q6*Q8*Q9*Q12*Q15+Q6*Q8*Q9*Q14*Q15,
           Q3*Q6*Q9-Q1*Q11^2,
           -Q4*Q5*Q7*Q8*Q9-Q4^2*Q8^2*Q9+Q2*Q4*Q8*Q9*Q13+Q1*Q5*Q7*Q12*Q13+3*Q1*Q4*Q8* Q12*Q13-2*Q1*Q5*Q8*Q13^2+Q1*Q2*Q12*Q13^2-Q1* Q4*Q7*Q12*Q14-Q1*Q4*Q8*Q13*Q14+Q2*Q4* Q7*Q9*Q15-Q1*Q2*Q13^2*Q15,
           -Q5*Q8*Q11+Q5*Q6*Q12+Q2*Q11*Q12-Q5*Q6* Q14,
           Q5^3*Q7*Q9^2+Q4*Q5^2*Q7*Q9*Q11-2*Q3*Q5^2*Q7*Q9*Q13-Q1*Q5^2*Q9*Q12*Q13-Q3*Q4* Q5*Q6*Q13^2+Q3^2*Q5*Q7*Q13^2-Q3*Q4*Q5*Q7*Q9*Q14-Q3*Q4^2*Q7*Q11*Q14+Q3*Q4^2*Q6* Q13*Q14+Q1*Q4*Q5*Q11*Q13*Q14+Q1*Q3*Q5*Q13^2*Q14+2*Q3*Q4*Q5*Q7*Q9*Q15+Q1*Q4*Q5* Q9*Q12*Q15+Q3*Q4^2*Q6*Q13*Q15-Q3^2*Q4*Q7* Q13*Q15-Q1*Q5^2*Q9*Q13*Q15-Q1*Q4*Q5*Q11*Q13*Q15+Q1*Q3*Q5*Q13^2*Q15-Q1*Q3*Q4* Q13*Q14*Q15-Q1*Q3*Q4*Q13*Q15^2,
           Q5*Q6*Q12-Q3*Q8*Q15,
           -Q5*Q6*Q12+Q3*Q8*Q15,
           -Q4*Q8^2*Q9^2+Q4*Q7*Q8*Q9*Q11-Q5*Q6*Q7*Q9*Q13+2*Q1*Q8*Q9*Q12*Q13-Q1*Q8*Q11* Q13^2+Q1*Q6*Q12*Q13^2-Q1*Q7*Q9*Q12*Q14+Q1*Q7*Q11*Q13*Q14-Q1*Q6*Q13^2*Q14-Q1*Q7*Q9* Q12*Q15+Q1*Q7*Q9*Q14*Q15,
           Q2*Q3*Q4*Q9*Q12-Q1*Q5^2*Q9*Q12+Q2*Q3*Q5*Q9*Q13-Q2*Q3*Q4*Q11*Q13-2*Q1*Q5^2* Q11*Q13-Q2*Q3^2*Q13^2-Q1*Q3*Q5*Q13*Q14+Q1*Q5^2*Q9*Q15+Q1*Q4*Q5*Q11*Q15+3*Q1*Q3*Q5* Q13*Q15-Q1*Q3*Q4*Q14*Q15,
           -Q5^2*Q8*Q9*Q11+Q4*Q5*Q8*Q11^2+Q5^2*Q6*Q9*Q12-Q4*Q5*Q6*Q11*Q12+Q3*Q4*Q8*Q11* Q12+Q5^2*Q6*Q9*Q14-Q3*Q5*Q8*Q9*Q14-Q4*Q5*Q6*Q11*Q14+Q3*Q4*Q8*Q11*Q14-Q3*Q4* Q6*Q12*Q14+Q3^2*Q8*Q13*Q14-Q3^2*Q8*Q13*Q15,
           -Q5^2*Q9-Q4*Q5*Q11+Q3*Q5*Q13+Q3*Q4*Q14,
           -Q5^2*Q9+Q4*Q5*Q11+Q3*Q5*Q13-Q3*Q4* Q14,
           Q5^2*Q9-Q4*Q5*Q11+Q3*Q5*Q13-Q3*Q4* Q14,
           Q5^2*Q9+Q4*Q5*Q11-Q3*Q5*Q13-Q3*Q4* Q14,
           -Q5^2*Q9+Q4*Q5*Q11-Q3*Q5*Q13+Q3*Q4* Q14,
           Q4*Q8^2*Q9*Q12-Q5*Q8^2*Q9*Q13+Q5*Q6*Q7*Q12*Q13-Q4*Q6*Q8*Q12*Q13+Q5*Q6* Q8*Q13^2-Q5*Q7*Q8*Q9*Q14+Q4*Q8^2*Q9*Q14+Q5*Q7^2*Q11*Q14-Q4*Q6*Q7*Q12*Q14+Q5*Q6* Q7*Q13*Q14-Q4*Q6*Q8*Q13*Q14-Q5*Q7^2*Q11*Q15,
           -Q5*Q8*Q9^3+Q4*Q8*Q9^2*Q11-Q5*Q6*Q9^2*Q13+Q3*Q8*Q9^2*Q13+Q3*Q6*Q9*Q13^2-Q1* Q11^2*Q13^2+Q1*Q9^2*Q12*Q15-Q1*Q9*Q11*Q13*Q15,
           Q5^2*Q7*Q9+Q4*Q5*Q8*Q9-Q4^2*Q8*Q11- Q2*Q4*Q9*Q12-2*Q3*Q5*Q7*Q13+2*Q1*Q5*Q12*Q13+Q2*Q3*Q13^2-2*Q1*Q5*Q13*Q14+Q3*Q4*Q7* Q15-Q1*Q4*Q12*Q15+Q1*Q4*Q14*Q15,
           -Q5*Q8*Q9^3-Q5*Q7*Q9^2*Q11+Q4*Q8*Q9^2*Q11+Q4*Q7*Q9*Q11^2+Q3*Q8*Q9^2*Q13-Q1*Q9* Q11*Q12*Q13-Q1*Q11^2*Q13^2+Q1*Q9^2*Q12*Q15,
           Q3*Q7*Q9-Q1*Q9*Q12-Q1*Q11*Q13+Q1*Q9* Q15,
           Q2*Q8*Q9^2*Q13-Q2*Q7*Q9*Q11*Q13-Q1*Q7*Q11*Q12*Q13-Q2*Q6*Q9*Q13^2+2*Q1*Q8* Q11*Q13^2-Q1*Q6*Q12*Q13^2-Q1*Q8*Q9*Q13*Q14+Q1*Q6*Q13^2*Q14+Q2*Q7*Q9^2*Q15+Q1*Q7*Q9* Q12*Q15-Q1*Q7*Q9*Q14*Q15,
           -Q2*Q5*Q7*Q9+Q2*Q4*Q7*Q11-Q1*Q4*Q8* Q12-Q2*Q4*Q6*Q13-2*Q1*Q5*Q8*Q13+Q2^2* Q9*Q13+Q1*Q2*Q12*Q13+Q1*Q4*Q8*Q14+Q1*Q2*Q13*Q15,
           -Q5*Q8*Q9*Q12-Q4*Q8*Q11*Q12+ Q5*Q8*Q11*Q13-Q5*Q8*Q9*Q14+Q4*Q8*Q11*Q14-Q5*Q6*Q13*Q14+Q3*Q8*Q13*Q14-Q5*Q6*Q13*Q15+Q3*Q8*Q13*Q15-Q3*Q7*Q14*Q15+Q2*Q9*Q14*Q15+Q1*Q12*Q14*Q15-Q1*Q14^2*Q15+Q2* Q9*Q15^2+Q1*Q12*Q15^2-Q1*Q14*Q15^2,
           Q3*Q4*Q8*Q9+Q4^2*Q6*Q11-Q2*Q4*Q9*Q11-Q1*Q5*Q9*Q12+Q1*Q4*Q11*Q12-Q3*Q4*Q6*Q13-2*Q1*Q5*Q11*Q13+Q1*Q4*Q11*Q14+Q1* Q5*Q9*Q15,
           Q2*Q13^2-Q4*Q7*Q14,
           -Q4^2*Q8*Q11*Q12+Q5^2*Q7*Q11*Q13-Q2*Q5*Q11*Q13^2+Q4*Q5*Q7*Q11*Q14+Q4^2*Q8*Q11*Q14- Q3*Q5*Q7*Q13*Q14-Q2*Q4*Q11*Q13*Q14+Q2*Q3*Q13^2*Q14+Q4*Q5*Q7*Q11*Q15-Q3*Q5*Q7* Q13*Q15+Q2*Q3*Q13^2*Q15-Q3*Q4*Q7*Q14*Q15,
           Q2*Q3*Q9-Q1*Q5*Q11-Q1*Q3*Q14+Q1*Q3* Q15,
           Q4*Q8*Q12-Q5*Q7*Q14,
           Q5*Q8-Q2*Q12,
           -Q5^2*Q7*Q9+Q4*Q5*Q8*Q9-Q4*Q5*Q6*Q13+2*Q3*Q5*Q7*Q13-Q3*Q4*Q8*Q13-Q2*Q4* Q11*Q13-Q2*Q3*Q13^2-Q4^2*Q6*Q14+Q3*Q4*Q7*Q14-Q1*Q4*Q12*Q14+2*Q1*Q5*Q13*Q15+Q1*Q4* Q14*Q15,
           -2*Q5^2*Q8*Q11-Q2*Q5*Q11*Q12-Q3*Q5* Q8*Q14+Q5^2*Q6*Q15+Q3*Q5*Q8*Q15+Q2*Q5*Q11*Q15+Q2*Q3*Q15^2,
           -Q2*Q4*Q7*Q8*Q11+Q1*Q4*Q8^2*Q12+Q2*Q5*Q6*Q7*Q13-Q1*Q5*Q8^2*Q13+2*Q1*Q2*Q8* Q12*Q13-Q2^2*Q6*Q13^2+Q1*Q5*Q7*Q8*Q14-Q1*Q4*Q8^2*Q14-Q1*Q2*Q7*Q12*Q14-Q1*Q2*Q7* Q12*Q15+Q1*Q2*Q7*Q14*Q15,
           -Q5*Q6^2*Q8*Q9+Q3*Q6*Q8^2*Q9-Q1*Q8^2*Q11^2-Q5*Q6^3*Q13+Q3*Q6^2*Q8*Q13+Q2*Q6^2*Q11*Q13-Q1*Q6*Q8*Q11*Q14+Q1*Q6^2*Q12*Q14,
           -Q4*Q5*Q6*Q9^2*Q12+Q2*Q4*Q9^2*Q11*Q12+Q4^2*Q6*Q11^2*Q13-Q2*Q4*Q9*Q11^2*Q13-Q3*Q4* Q6*Q9*Q12*Q13+Q1*Q5*Q9*Q11*Q12*Q13+Q1*Q4*Q11^2*Q12*Q13+Q3*Q5*Q6*Q9*Q13^2-2*Q3*Q4*Q6*Q11*Q13^2+Q3^2*Q6*Q13^3-Q4^2*Q6*Q9*Q11* Q14+Q2*Q4*Q9^2*Q11*Q14-Q1*Q4*Q9*Q11*Q12* Q14+2*Q3*Q4*Q6*Q9*Q13*Q14-Q1*Q5*Q9*Q11* Q13*Q14+Q1*Q4*Q11^2*Q13*Q14-Q1*Q3*Q11*Q13^2*Q14-Q1*Q4*Q9*Q11*Q14^2-Q1*Q3*Q11*Q13^2* Q15+Q1*Q3*Q9*Q13*Q14*Q15,
           -Q5^2*Q7*Q9*Q11-Q4^2*Q8*Q11^2+Q3*Q5*Q7*Q9*Q12-Q2*Q5*Q9^2*Q12+2*Q2*Q4*Q9*Q11*Q12- Q1*Q5*Q9*Q12^2+Q2*Q5*Q9*Q11*Q13-Q1*Q5* Q11*Q12*Q13+Q1*Q4*Q11*Q12*Q14+Q3*Q5*Q7* Q9*Q15-Q1*Q5*Q9*Q12*Q15-Q2*Q3*Q9*Q13* Q15+Q1*Q5*Q11*Q13*Q15,
           Q5*Q13-Q4*Q14,
           Q5*Q9*Q12^2+Q4*Q11*Q12^2-Q5*Q11*Q12*Q13- Q5*Q11*Q13*Q14+Q3*Q12*Q13*Q14-Q5*Q11*Q13*Q15+Q3*Q12*Q13*Q15-Q3*Q13*Q14*Q15,
           Q8^2*Q9^2*Q12-Q7*Q8*Q9*Q11*Q12-Q6^2*Q12* Q13^2-Q8^2*Q9^2*Q14+Q7^2*Q11^2*Q14+Q6*Q7*Q11*Q13*Q14+Q6^2*Q13^2*Q14-Q6*Q7*Q9*Q14^2+Q7*Q8* Q9*Q11*Q15-Q7^2*Q11^2*Q15+Q6*Q7*Q9*Q12*Q15- Q6*Q7*Q9*Q14*Q15,
           -Q4^2*Q8*Q11*Q12-Q5^2*Q8*Q9*Q13+Q3*Q5*Q8*Q13^2+Q5^2*Q7*Q9*Q14-Q4*Q5*Q8*Q9*Q14+ Q4^2*Q8*Q11*Q14-Q3*Q5*Q7*Q13*Q14+Q3*Q4*Q8*Q13*Q14+Q5^2*Q7*Q9*Q15-Q3*Q5*Q7*Q13* Q15+Q3*Q4*Q8*Q13*Q15-Q3*Q4*Q7*Q14*Q15,
           Q1*Q8*Q12-Q2*Q6*Q13,
           -Q1*Q8*Q12+Q2*Q6*Q13,
           Q5^2*Q6*Q9-Q3*Q4*Q8*Q11+Q2*Q4*Q11^2-2* Q1*Q5*Q11*Q12+Q3*Q5*Q6*Q13-Q3^2*Q7*Q14+Q2*Q3*Q9*Q14+Q1*Q3*Q12*Q14-Q1*Q3*Q14^2+ Q1*Q3*Q12*Q15-Q1*Q3*Q14*Q15,
           Q5^2*Q6*Q9+Q3*Q4*Q8*Q11+Q2*Q4*Q11^2-2* Q1*Q5*Q11*Q12-Q3*Q5*Q6*Q13-Q3^2*Q7*Q14+Q2*Q3*Q9*Q14+Q1*Q3*Q12*Q14-Q1*Q3*Q14^2+ Q1*Q3*Q12*Q15-Q1*Q3*Q14*Q15,
           -Q5*Q7*Q11+Q4*Q8*Q11+Q2*Q9*Q12-Q1*Q12* Q14,
           -Q4*Q8*Q9*Q12+Q4*Q8*Q11*Q13-2*Q3* Q7*Q12*Q13+2*Q2*Q9*Q12*Q13+2*Q1*Q12^2*Q13-Q5*Q6*Q13^2-Q1*Q12*Q13*Q14-Q5*Q7*Q9*Q15+Q3* Q7*Q13*Q15,
           -Q5*Q6*Q13+Q3*Q8*Q13-Q2*Q11* Q13+Q2*Q9*Q15+Q1*Q12*Q15-Q1*Q14*Q15,
           -Q5*Q9+Q4*Q11,
           Q5*Q9-Q4*Q11,
           Q5^3*Q6*Q9+Q3*Q5^2*Q8*Q9+Q3*Q5^2*Q7* Q11-2*Q2*Q5^2*Q9*Q11+Q2*Q4*Q5*Q11^2-2*Q3^2*Q5*Q7*Q12+2*Q2*Q3*Q5*Q9*Q12-Q2*Q3*Q4* Q11*Q12-2*Q1*Q5^2*Q11*Q12+2*Q1*Q3*Q5*Q12^2-Q1* Q3*Q5*Q12*Q14+2*Q3^2*Q5*Q7*Q15-Q2*Q3*Q5*Q9*Q15-Q1*Q3*Q5*Q12*Q15-Q2*Q3^2*Q13* Q15+Q1*Q3*Q5*Q14*Q15,
           Q2*Q11*Q12-Q3*Q8*Q14,
           -Q2*Q11*Q12+Q3*Q8*Q14,
           Q5*Q6*Q7*Q9+Q3*Q7^2*Q11-Q2*Q7*Q9*Q11+Q1*Q8*Q9*Q12-Q3*Q6*Q7*Q13-2*Q1*Q8* Q11*Q13+Q1*Q7*Q11*Q14-Q1*Q8*Q9*Q15+Q1*Q7*Q11*Q15,
           Q3*Q4*Q5*Q7*Q8-Q3*Q4^2*Q8^2-Q1*Q4*Q5*Q8*Q12-Q2*Q3*Q4*Q8*Q13-2*Q1* Q5^2*Q8*Q13+Q1*Q5^2*Q7*Q14+3*Q1*Q4*Q5*Q8* Q14-Q1*Q2*Q4*Q12*Q14+Q1*Q2*Q5*Q13*Q14+ Q2*Q3*Q4*Q7*Q15-Q1*Q5^2*Q7*Q15,
           -Q5^2*Q7*Q9^2-Q4*Q5*Q7*Q9*Q11+Q3*Q5*Q7*Q9*Q13-Q1*Q5*Q9*Q12*Q13-2*Q1*Q5*Q11* Q13^2+Q3*Q4*Q7*Q9*Q14-Q1*Q3*Q13^2*Q14-Q1*Q4*Q9*Q12*Q15+3*Q1*Q5*Q9*Q13*Q15+Q1*Q4* Q11*Q13*Q15+Q1*Q3*Q13^2*Q15,
           Q2*Q7*Q9-Q1*Q8*Q13,
           Q4*Q8*Q11+Q2*Q9*Q12-Q3*Q8*Q13-Q1*Q12* Q14,
           Q5^2*Q7*Q9+Q4*Q5*Q8*Q9-Q4^2*Q6* Q12+Q3*Q4*Q7*Q12-Q1*Q4*Q12^2-Q2*Q4*Q11*Q13+Q2*Q3*Q13^2-Q1*Q4*Q12*Q14+Q1*Q4*Q12* Q15-2*Q1*Q5*Q13*Q15+Q1*Q4*Q14*Q15,
           Q5^2*Q7*Q9-Q4*Q5*Q8*Q9-Q4^2*Q6*Q12+ Q3*Q4*Q7*Q12-Q1*Q4*Q12^2+Q2*Q4*Q11*Q13+Q2*Q3*Q13^2-Q1*Q4*Q12*Q14+Q1*Q4*Q12*Q15-2* Q1*Q5*Q13*Q15+Q1*Q4*Q14*Q15,
           Q4*Q6-Q3*Q7,
           -Q4*Q6+Q3*Q7,
           -Q5*Q8*Q11+Q5*Q6*Q12+Q2*Q11*Q12-Q2*Q11* Q15,
           -Q5*Q8*Q9+Q5*Q7*Q11+Q4*Q6*Q14-Q1*Q14*Q15,
           -Q5*Q7*Q11*Q12+Q4*Q8*Q11*Q12+ Q5*Q8*Q11*Q13-Q5*Q6*Q12*Q13+Q2*Q11*Q12*Q13-Q5*Q7*Q11*Q14-Q4*Q8*Q11*Q14+Q3*Q7*Q12*Q15-Q2*Q9*Q12*Q15-Q1*Q12^2*Q15-Q5*Q6* Q13*Q15+Q2*Q11*Q13*Q15+Q1*Q12*Q14*Q15+Q3*Q7*Q15^2-Q1*Q12*Q15^2+Q1*Q14*Q15^2,
           -Q8*Q9*Q12^2+Q8*Q11*Q12*Q13-Q6*Q12^2*Q13+Q7*Q11*Q12*Q14-Q8*Q11*Q13*Q14+Q6*Q12*Q13* Q14+Q8*Q9*Q12*Q15+Q7*Q11*Q12*Q15-Q8*Q11*Q13*Q15-Q7*Q11*Q14*Q15,
           Q5^2*Q7*Q9*Q12-Q4*Q5*Q7*Q11*Q12+Q4^2* Q8*Q11*Q12-Q5^2*Q7*Q11*Q13-Q2*Q5*Q9*Q12*Q13+Q2*Q4*Q11*Q12*Q13+Q2*Q5*Q11*Q13^2-Q4^2* Q8*Q11*Q14+Q5^2*Q7*Q9*Q15-Q2*Q4*Q9*Q12*Q15-Q2*Q5*Q9*Q13*Q15+Q2*Q4*Q11*Q13*Q15,
           -Q5^2*Q6*Q9-Q3*Q5*Q8*Q9-Q3*Q5*Q7*Q11+Q3*Q4*Q8*Q11+2*Q2*Q5*Q9*Q11-Q2*Q4* Q11^2+2*Q1*Q5*Q11*Q12-Q3*Q5*Q6*Q13-Q3^2*Q7*Q15+Q2*Q3*Q9*Q15+Q1*Q3*Q12*Q15-Q1*Q3* Q14*Q15,
           -Q4*Q5*Q11*Q12-2*Q5^2*Q11*Q13+Q5^2*Q9*Q14+Q4*Q5*Q11*Q14+Q3*Q5*Q13*Q14+Q3*Q4*Q14^2-Q3*Q5*Q13*Q15,
           -Q2*Q5^2*Q9^2+Q2*Q4*Q5*Q9*Q11+Q2*Q3*Q4*Q9*Q12-Q1*Q5^2*Q9*Q12-Q3^2*Q4*Q8*Q13+2* Q2*Q3*Q5*Q9*Q13-Q2*Q3*Q4*Q11*Q13-Q2*Q3^2*Q13^2+Q1*Q3*Q5*Q13*Q14-Q1*Q5^2*Q9*Q15+Q1*Q3*Q5*Q13*Q15,
           Q4*Q6-Q1*Q15,
           Q4*Q6-Q2*Q9,
           -Q4*Q6+Q2*Q9,
           -Q3*Q7^2+Q2*Q7*Q9+Q1*Q7*Q12-Q1*Q8*Q13,
           -Q5*Q8*Q9+Q5*Q7*Q11+Q1*Q12^2-Q5*Q6* Q13,
           Q5^3*Q6*Q7^2-2*Q4*Q5^2*Q6*Q7*Q8+Q4^2* Q5*Q6*Q8^2-Q2*Q4*Q5*Q8^2*Q9-Q2*Q4*Q5*Q6*Q7*Q12+Q1*Q4*Q5*Q8^2*Q12+Q2^2*Q4*Q8* Q9*Q12+Q2*Q5^2*Q6*Q7*Q13-Q2^2*Q4*Q6*Q12* Q13+Q1*Q2*Q5*Q8*Q12*Q13+2*Q2*Q4*Q5*Q6*Q7*Q14-Q2*Q4^2*Q6*Q8*Q14-Q1*Q5^2*Q7* Q8*Q14+Q1*Q4*Q5*Q8^2*Q14+Q2^2*Q4*Q8*Q9*Q14-Q1*Q2*Q4*Q8*Q12*Q14-Q1*Q2*Q5*Q8* Q13*Q14-Q1*Q2*Q4*Q8*Q14^2-Q1*Q5^2*Q7*Q8* Q15+Q1*Q2*Q5*Q7*Q14*Q15,
           Q5*Q8-Q2*Q15,
           -Q5*Q7*Q11*Q12+Q2*Q9*Q12^2+Q5*Q8*Q11*Q13+Q5*Q6*Q12*Q13-Q5*Q7*Q11*Q14-Q4*Q6*Q12*Q14+Q2*Q9*Q12*Q14-Q1*Q12^2*Q14+Q5*Q6* Q13*Q14+Q3*Q8*Q13*Q14-Q2*Q11*Q13*Q14-Q1*Q12*Q14^2+Q1*Q12^2*Q15-Q3*Q8*Q13*Q15-Q2*Q11* Q13*Q15+Q1*Q12*Q14*Q15,
           -Q4^2*Q8*Q11^2+Q3*Q5*Q7*Q11*Q13-Q2*Q3*Q11*Q13^2+2*Q3*Q4*Q7*Q11*Q14+Q1*Q4*Q11*Q12* Q14-Q3^2*Q7*Q13*Q14+Q2*Q3*Q9*Q13*Q14-Q1*Q5*Q11*Q13*Q14-Q1*Q3*Q13*Q14^2-Q3*Q5*Q7* Q9*Q15+Q2*Q3*Q9*Q13*Q15+Q1*Q5*Q11*Q13*Q15-Q1*Q3*Q13*Q14*Q15,
           -Q8*Q9*Q12+Q6*Q13*Q14,
           Q5*Q7*Q9+Q4*Q8*Q9-Q4*Q7*Q11-Q2*Q9* Q13,
           Q5*Q6-Q3*Q8,
           -Q5*Q6+Q3*Q8,
           Q5*Q7-Q2*Q13,
           -Q5*Q7+Q2*Q13,
           -Q3*Q8^2*Q9^2-Q5*Q6*Q7*Q9*Q11+Q3*Q6*Q8*Q9*Q13-Q1*Q8*Q11^2*Q13-Q1*Q7*Q11^2*Q14+ Q1*Q6*Q9*Q12*Q14+Q1*Q6*Q11*Q13*Q14+2*Q1*Q8*Q9*Q11*Q15+Q1*Q7*Q11^2*Q15-Q1*Q6*Q9* Q12*Q15-Q1*Q6*Q9*Q14*Q15,
           -Q5*Q11*Q13-Q3*Q13*Q14+Q5*Q9*Q15+Q3*Q13* Q15,
           -Q3*Q7*Q8*Q9*Q13+Q3*Q7^2*Q11*Q13- Q1*Q7*Q11*Q12*Q13-Q3*Q6*Q7*Q13^2+2*Q1*Q8*Q11*Q13^2+Q1*Q6*Q12*Q13^2-Q1*Q8*Q9*Q13*Q14- Q1*Q6*Q13^2*Q14+Q3*Q7^2*Q9*Q15-Q1*Q7*Q9* Q12*Q15+Q1*Q7*Q9*Q14*Q15,
           Q4*Q8^2*Q9-Q5*Q7^2*Q11+Q2*Q7*Q9*Q12+ Q5*Q6*Q7*Q13-2*Q2*Q8*Q9*Q13+Q2*Q6*Q13^2-Q4*Q6*Q7*Q14-Q1*Q7*Q12*Q14+2*Q1*Q8*Q13* Q14+Q1*Q7*Q12*Q15-2*Q1*Q8*Q13*Q15,
           Q5*Q6*Q8^2*Q9+Q3*Q8^3*Q9+Q4*Q6*Q8^2* Q11-2*Q2*Q8^2*Q9*Q11+Q2*Q7*Q8*Q11^2+2*Q4*Q6^ 2*Q8*Q12-Q2*Q6*Q8*Q9*Q12-Q2*Q6^2*Q12* Q13+Q1*Q6*Q8*Q12*Q14-2*Q4*Q6^2*Q8*Q15+2*Q2*Q6*Q8*Q9*Q15-Q2*Q6*Q7*Q11*Q15-2*Q1*Q8^2* Q11*Q15-Q1*Q6*Q8*Q12*Q15-Q1*Q6*Q8*Q14* Q15+2*Q1*Q6*Q8*Q15^2,
           Q4*Q8*Q9-Q1*Q13*Q14,
           Q2*Q6*Q13-Q1*Q8*Q15,
           Q3*Q4*Q8-Q1*Q5*Q12,
           Q3*Q5^2*Q6*Q9+Q3^2*Q5*Q8*Q9+Q3^2*Q4* Q8*Q11-Q1*Q5^2*Q11^2-Q3^2*Q5*Q6*Q13-Q3^3*Q8* Q13-Q1*Q3*Q5*Q11*Q15+Q1*Q3^2*Q14*Q15,
           -Q5*Q9*Q12^2-Q4*Q11*Q12^2+Q5*Q11*Q12*Q13+Q4*Q11*Q12*Q14-Q5*Q11*Q13*Q14+Q3*Q12*Q13* Q14+Q5*Q9*Q12*Q15-Q5*Q11*Q13*Q15+Q3*Q12*Q13*Q15-Q3*Q13*Q14*Q15,
           -Q5*Q8^2*Q9+Q5*Q7*Q8*Q11-Q5*Q6*Q8*Q13+Q2*Q6*Q12*Q13+Q5*Q6*Q7*Q14-Q4*Q6*Q8*Q14-Q1*Q8*Q12*Q14-Q4*Q6*Q8*Q15+Q3*Q7*Q8*Q15-Q2*Q7*Q11*Q15+Q1*Q8*Q14*Q15+Q1*Q8*Q15^2,
           Q4*Q5*Q7*Q8*Q9-Q4^2*Q8^2*Q9-Q2*Q4*Q7*Q11*Q13+2*Q1*Q4*Q8*Q12*Q13-Q1* Q5*Q8*Q13^2+Q1*Q2*Q12*Q13^2-Q1*Q4*Q7*Q12* Q14-Q1*Q4*Q7*Q12*Q15+Q1*Q5*Q7*Q13*Q15- Q1*Q2*Q13^2*Q15+Q1*Q4*Q7*Q14*Q15,
           Q5^3*Q6*Q7-2*Q4*Q5^2*Q6*Q8+Q3*Q4*Q5* Q8^2+Q2*Q5^2*Q8*Q9-Q2*Q4*Q5*Q6*Q12+2*Q2^2* Q5*Q9*Q12-Q2^2*Q4*Q11*Q12+Q2*Q5^2*Q6* Q13+2*Q2*Q4*Q5*Q6*Q14-Q2*Q3*Q4*Q8*Q14-2*Q1*Q5^2*Q8*Q14-2*Q2^2*Q5*Q9*Q14-Q1*Q2*Q5* Q12*Q14+2*Q1*Q2*Q5*Q14^2+Q1*Q2*Q5*Q12*Q15-Q1*Q2*Q5*Q14*Q15,
           Q3*Q8*Q14-Q2*Q11*Q15,
           -Q3*Q7*Q8*Q9*Q12+Q1*Q8*Q9*Q12^2+Q3*Q8^2*Q9*Q13+Q1*Q8*Q11*Q12*Q13-Q3*Q7^2*Q11* Q14+Q1*Q7*Q11*Q12*Q14+Q3*Q6*Q7*Q13*Q14-Q1*Q8*Q11*Q13*Q14-Q1*Q6*Q12*Q13*Q14-Q1* Q8*Q11*Q13*Q15-Q1*Q7*Q11*Q14*Q15+Q1*Q6*Q13*Q14*Q15,
           -Q3*Q4*Q5*Q8*Q9-Q1*Q5^2*Q9* Q12+Q1*Q4*Q5*Q11*Q12+Q2*Q3*Q4*Q11*Q13-Q1*Q5^2*Q11*Q13-Q2*Q3^2*Q13^2+Q1*Q3*Q4*Q12* Q14+Q1*Q5^2*Q9*Q15-Q1*Q3*Q4*Q12*Q15+2*Q1* Q3*Q5*Q13*Q15-Q1*Q3*Q4*Q14*Q15,
           Q5^3*Q6*Q9-2*Q4*Q5^2*Q6*Q11+Q3*Q5^2*Q7* Q11+Q2*Q4*Q5*Q11^2+2*Q3*Q4*Q5*Q6*Q12-2*Q3^ 2*Q5*Q7*Q12-Q2*Q3*Q4*Q11*Q12-2*Q1*Q5^2* Q11*Q12+2*Q1*Q3*Q5*Q12^2+Q3*Q5^2*Q6*Q13-Q3* Q4*Q5*Q6*Q14+2*Q3^2*Q5*Q7*Q14-Q3^2*Q4*Q8*Q14-Q1*Q3*Q5*Q12*Q14-Q1*Q3*Q5*Q12* Q15+Q1*Q3*Q5*Q14*Q15,
           -Q5^2*Q7^2*Q11-Q2*Q4*Q8*Q9*Q12+Q2*Q5*Q8*Q9*Q13+Q2*Q4*Q6*Q12*Q13+Q1*Q5*Q8* Q12*Q13-Q2*Q5*Q6*Q13^2+2*Q2*Q5*Q7*Q9*Q15+Q2*Q4*Q6*Q13*Q15-Q1*Q5*Q8*Q13*Q15-Q2^2* Q9*Q13*Q15-Q1*Q2*Q12*Q13*Q15+Q1*Q5*Q7*Q14*Q15-Q1*Q2*Q13*Q15^2,
           Q4*Q8*Q11-Q1*Q12^2,
           Q5^2*Q9*Q11*Q12-Q4*Q5*Q11^2*Q12+Q5^2*Q9* Q11*Q14-Q4*Q5*Q11^2*Q14-Q3*Q4*Q11*Q12*Q14-Q5^2*Q9*Q11*Q15+Q4*Q5*Q11^2*Q15+Q3*Q4*Q11* Q12*Q15-Q3*Q5*Q9*Q14*Q15+Q3*Q4*Q11*Q14*Q15+Q3^2*Q13*Q14*Q15-Q3^2*Q13*Q15^2,
           -Q5*Q8*Q11*Q12-Q5*Q8*Q11*Q14-Q5*Q6*Q12* Q14-Q5*Q8*Q11*Q15+Q5*Q6*Q12*Q15+Q5*Q6* Q14*Q15+Q3*Q8*Q15^2+Q2*Q11*Q15^2,
           Q6*Q13^2-Q7*Q9*Q15,
           Q4*Q11^2-Q3*Q9*Q15,
           Q2*Q5*Q6*Q7*Q11-Q1*Q5*Q8^2*Q11-Q2^2* Q7*Q11^2-Q2*Q3*Q6*Q8*Q13+Q1*Q5*Q6*Q8*Q14-Q1*Q3*Q8^2*Q14+Q1*Q2*Q6*Q12*Q14+Q1* Q3*Q8^2*Q15+2*Q1*Q2*Q8*Q11*Q15-Q1*Q2*Q6*Q12*Q15-Q1*Q2*Q6*Q14*Q15,
           Q5*Q8*Q9-Q5*Q7*Q11+Q4*Q6*Q12-Q5*Q6* Q13+Q1*Q12*Q14-Q1*Q12*Q15,
           -Q3*Q8^2*Q9^2+Q3*Q7*Q8*Q9*Q11-Q1*Q8*Q9*Q11*Q12-Q3*Q6*Q8*Q9*Q13-2*Q1*Q8*Q11^2* Q13+Q3*Q6*Q7*Q9*Q14-Q1*Q7*Q11^2*Q14+3*Q1* Q8*Q9*Q11*Q15+Q1*Q7*Q11^2*Q15-Q1*Q6*Q9* Q12*Q15+Q1*Q6*Q11*Q13*Q15,
           -Q5^2*Q8*Q9+Q2*Q5*Q11*Q13+Q4*Q5*Q6*Q14-Q3*Q4*Q8*Q14-2*Q4*Q5*Q6*Q15+2*Q3*Q5* Q7*Q15-Q1*Q5*Q12*Q15-Q2*Q3*Q13*Q15+2*Q1*Q5*Q15^2,
           Q8^2*Q9^2*Q12-Q7^2*Q11^2*Q12-Q6*Q7*Q9*Q12^2+Q6*Q8*Q9*Q12*Q13+Q6^2*Q12*Q13^2-Q7* Q8*Q9*Q11*Q14+Q7^2*Q11^2*Q14-Q6^2*Q13^2*Q14-Q8^2*Q9^2*Q15+Q7*Q8*Q9*Q11*Q15-Q6*Q7*Q9* Q12*Q15+Q6*Q7*Q9*Q14*Q15,
           -Q4*Q5*Q6*Q9+Q3*Q4*Q8*Q9+Q4^2*Q6*Q11-Q3*Q4*Q7*Q11+Q1*Q4*Q11*Q12-2*Q1*Q5* Q11*Q13+Q1*Q4*Q11*Q14-Q1*Q3*Q13*Q14+Q1*Q3*Q13*Q15,
           Q5*Q6*Q7*Q9-Q4*Q6*Q7*Q11-Q1*Q8*Q9*Q12+Q4*Q6^2*Q13-Q2*Q6*Q9*Q13-2* Q1*Q8*Q11*Q13+Q1*Q6*Q12*Q13+Q1*Q6*Q13*Q14+Q1*Q8*Q9*Q15,
           -Q3*Q8^2*Q9-Q5*Q6*Q7*Q11-Q4*Q6*Q8*Q11+2*Q2*Q8*Q9*Q11-Q2*Q7*Q11^2-Q4*Q6^2*Q12+ Q2*Q6*Q9*Q12+Q3*Q6*Q8*Q13-Q2*Q6*Q11*Q13-Q1*Q6*Q12*Q14+2*Q1*Q8*Q11*Q15+Q1*Q6* Q12*Q15,
           Q5^2*Q7*Q9-Q4^2*Q8*Q11+Q3*Q4*Q8* Q13-2*Q2*Q5*Q9*Q13-2*Q1*Q5*Q12*Q13+Q2*Q3*Q13^2-Q3*Q4*Q7*Q14+2*Q1*Q5*Q13*Q14+Q2* Q4*Q9*Q15+Q1*Q4*Q12*Q15-Q1*Q4*Q14*Q15,
           -Q5^2*Q6^2*Q9+Q3^2*Q8^2*Q9-Q3^2*Q7*Q8*Q11- Q2*Q3*Q8*Q9*Q11+Q2^2*Q9*Q11^2+2*Q1*Q5*Q6*Q11*Q12-Q1*Q3*Q8*Q11*Q12+Q1*Q2*Q11^2* Q12+Q3^2*Q6*Q8*Q13-Q2*Q3*Q6*Q9*Q14+Q1*Q3*Q8*Q11*Q14-Q1*Q2*Q11^2*Q14-Q1*Q3*Q6* Q12*Q14+Q1*Q3*Q6*Q14^2+Q2*Q3*Q6*Q9*Q15-Q1*Q3*Q8*Q11*Q15-Q1*Q2*Q11^2*Q15,
           -Q3*Q6*Q7*Q12*Q13+Q1*Q8*Q11*Q12*Q13+Q1* Q6*Q12^2*Q13+Q3*Q6*Q8*Q13^2-Q1*Q8*Q11* Q13*Q14+Q3*Q7*Q8*Q9*Q15-Q3*Q7^2*Q11*Q15-Q1*Q8*Q9*Q12*Q15+Q1*Q7*Q11*Q12*Q15-Q1* Q8*Q11*Q13*Q15+Q1*Q8*Q9*Q14*Q15-Q1*Q7*Q11*Q14*Q15,
           -Q4*Q5^2*Q6*Q7*Q13+Q3*Q5^2*Q7^2* Q13-2*Q2*Q3*Q5*Q7*Q13^2+Q2*Q3*Q4*Q8*Q13^ 2-Q1*Q2*Q5*Q12*Q13^2+Q2^2*Q3*Q13^3+Q4^2*Q5* Q6*Q7*Q14-Q3*Q4^2*Q7*Q8*Q14-Q2*Q3*Q4* Q7*Q13*Q14+Q1*Q5^2*Q7*Q13*Q14+Q1*Q4*Q5* Q8*Q13*Q14+Q4^2*Q5*Q6*Q7*Q15-Q3*Q4*Q5*Q7^2*Q15+2*Q2*Q3*Q4*Q7*Q13*Q15+Q1*Q5^2* Q7*Q13*Q15-Q1*Q4*Q5*Q8*Q13*Q15+Q1*Q2*Q4*Q12*Q13*Q15-Q1*Q2*Q5*Q13^2*Q15-Q1*Q4* Q5*Q7*Q14*Q15-Q1*Q4*Q5*Q7*Q15^2,
           -Q5^2*Q7^2*Q9+Q4*Q5*Q7*Q8*Q9-Q2*Q4*Q7*Q11*Q13+Q1*Q4*Q8*Q12*Q13-Q1*Q5*Q8* Q13^2-Q1*Q2*Q12*Q13^2+Q1*Q4*Q7*Q12*Q14-Q1*Q4*Q7*Q12*Q15+2*Q1*Q5*Q7*Q13*Q15+Q1*Q2* Q13^2*Q15-Q1*Q4*Q7*Q14*Q15,
           2*Q5^2*Q9*Q11*Q12*Q13+2*Q4*Q5*Q11^2*Q12*Q13-2* Q5^2*Q11^2*Q13^2+Q4^2*Q11^2*Q12*Q14-Q5^2*Q9*Q11*Q13*Q14-2*Q4*Q5*Q11^2*Q13*Q14+Q3*Q5*Q11*Q13^2 *Q14+Q3*Q4*Q11*Q13*Q14^2-3*Q4*Q5*Q9*Q11*Q12*Q15-Q4^2*Q11^2*Q12*Q15+Q4*Q5*Q11^2*Q13*Q15- Q3*Q5*Q9*Q12*Q13*Q15-Q3*Q4*Q11*Q12*Q13*Q15+Q3*Q5*Q11*Q13^2*Q15+Q4^2*Q11^2*Q14*Q15-Q3^2*Q13^2*Q14*Q15+Q5^2*Q9^2*Q15^2+Q3*Q4*Q9*Q12* Q15^2+Q3*Q5*Q9*Q13*Q15^2-Q3*Q4*Q9*Q14*Q15^2,
           Q14-Q15,
           -Q14+Q15,
           -Q5*Q8*Q9-Q4*Q8*Q11+Q3*Q8*Q13+Q4*Q6* Q14+Q1*Q12*Q14-Q1*Q14*Q15,
           Q8^2*Q9*Q12*Q13-Q6*Q8*Q12*Q13^2+Q8^2*Q9* Q13*Q14-Q6*Q7*Q12*Q13*Q14-Q6*Q8*Q13^2*Q14-Q8^2*Q9*Q13*Q15+Q6*Q7*Q12*Q13*Q15+Q6*Q8* Q13^2*Q15-Q7*Q8*Q9*Q14*Q15+Q7^2*Q11*Q14*Q15+ Q6*Q7*Q13*Q14*Q15-Q7^2*Q11*Q15^2,
           -Q4*Q5*Q6*Q7+Q4^2*Q6*Q8-Q2*Q4*Q8*Q9+Q2*Q4*Q7*Q11+Q1*Q4*Q8*Q12-2*Q1*Q5* Q8*Q13-Q1*Q2*Q12*Q13+Q1*Q4*Q8*Q14+Q1*Q2*Q13*Q15,
           -Q5*Q8*Q12*Q13-Q4*Q8*Q12*Q14- Q5*Q8*Q13*Q14+Q4*Q8*Q12*Q15-Q5*Q8*Q13*Q15+Q4*Q8*Q14*Q15+Q5*Q7*Q15^2+Q2*Q13*Q15^2,
           -Q5^2*Q6*Q7*Q11+Q4*Q5*Q6*Q8*Q11-Q4*Q5*Q6^2*Q12+Q3*Q5*Q6*Q7*Q12+2*Q2*Q4*Q6* Q11*Q12-Q1*Q5*Q8*Q11*Q12-Q1*Q5*Q6*Q12^2-Q2^2*Q11^2*Q13+Q3*Q5*Q6*Q7*Q14-Q3*Q4*Q6* Q8*Q14+Q1*Q5*Q8*Q11*Q14-Q1*Q5*Q6*Q12*Q14+Q1*Q2*Q11*Q12*Q15,
           Q5^2*Q6*Q7+Q3*Q4*Q8^2+Q2*Q5*Q7*Q11- Q2*Q3*Q8*Q13-2*Q1*Q5*Q8*Q14+Q1*Q2*Q12*Q14+Q2*Q4*Q6*Q15-Q2^2*Q9*Q15-Q1*Q2*Q12* Q15+Q1*Q2*Q14*Q15-Q1*Q2*Q15^2,
           Q5^2*Q6*Q7+Q3*Q4*Q8^2-Q2*Q5*Q7*Q11+ Q2*Q3*Q8*Q13-2*Q1*Q5*Q8*Q14+Q1*Q2*Q12*Q14+Q2*Q4*Q6*Q15-Q2^2*Q9*Q15-Q1*Q2*Q12* Q15+Q1*Q2*Q14*Q15-Q1*Q2*Q15^2,
           Q5*Q9^2*Q12^2-2*Q5*Q9*Q11*Q12*Q13-Q4*Q11^2* Q12*Q13+2*Q5*Q11^2*Q13^2-Q3*Q11*Q12*Q13^2-Q5*Q9*Q11*Q13*Q14+Q4*Q11^2*Q13*Q14+Q3*Q9*Q12* Q13*Q14+Q3*Q9*Q12*Q13*Q15-Q3*Q9*Q13*Q14*Q15,
           Q3*Q4*Q8-Q1*Q5*Q15,
           -Q5^2*Q7^2*Q11+Q4^2*Q8^2*Q11-Q2*Q4*Q7*Q11*Q12+Q2*Q4*Q8*Q11*Q13+Q2^2*Q11*Q13^2-Q2*Q4* Q7*Q11*Q14-Q1*Q4*Q8*Q12*Q14+Q1*Q5*Q8* Q13*Q14-Q2*Q4*Q7*Q11*Q15+Q1*Q5*Q8*Q13* Q15-Q1*Q2*Q12*Q13*Q15+Q1*Q5*Q7*Q14*Q15,
           -Q4*Q11*Q14+Q3*Q13*Q15,
           Q4*Q6*Q9-Q1*Q11*Q13,
           -Q5*Q6*Q14+Q2*Q11*Q15,
           Q5*Q6*Q14-Q2*Q11*Q15,
           Q3*Q4*Q7-Q1*Q5*Q13,
           Q5^2*Q6*Q8*Q9-Q1*Q5*Q8*Q11*Q12-Q2*Q5*Q6*Q9*Q14+Q1*Q5*Q8*Q11*Q14+Q1*Q5*Q6*Q14^2+Q2*Q3*Q8*Q9*Q15-Q1*Q5*Q8*Q11* Q15-Q2^2*Q9*Q11*Q15+Q1*Q3*Q8*Q12*Q15-Q1*Q2*Q11*Q12*Q15-Q1*Q3*Q8*Q14*Q15+Q1*Q2* Q11*Q14*Q15,
           Q4^2*Q6*Q8-Q3*Q4*Q7*Q8+Q2* Q4*Q7*Q11+Q1*Q4*Q8*Q12-Q2*Q4*Q6*Q13-2*Q1*Q5*Q8*Q13-Q1*Q5*Q7*Q14+Q1*Q4*Q8* Q14+Q1*Q5*Q7*Q15,
           -Q8*Q9*Q12-Q8*Q11*Q13+Q8*Q9*Q15+Q7*Q11* Q15,
           Q5^2*Q6^2*Q12-Q3^2*Q8^2*Q12+Q2*Q5*Q6*Q11*Q12+Q2^2*Q11^2*Q12-Q2*Q3*Q6*Q12^2-Q5^2* Q6^2*Q14+Q3*Q5*Q6*Q8*Q14-Q2*Q3*Q6*Q12*Q14-Q3*Q5*Q6*Q8*Q15+Q3^2*Q8^2*Q15-Q2^2*Q11^2* Q15+Q2*Q3*Q6*Q14*Q15,
           Q5^2*Q6*Q9+Q3*Q5*Q8*Q9-2*Q4*Q5*Q6* Q11+Q2*Q4*Q11^2+Q3*Q4*Q6*Q12-Q3^2*Q8*Q13-2*Q1*Q5*Q11*Q14+Q1*Q3*Q12*Q14-Q2*Q3*Q9* Q15+2*Q1*Q5*Q11*Q15-Q1*Q3*Q12*Q15,
           -Q5*Q7*Q9+Q4*Q7*Q11+Q3*Q7*Q13-Q1*Q12* Q13,
           Q5*Q7*Q9-Q4*Q7*Q11+Q3*Q7*Q13-Q1* Q12*Q13,
           Q4*Q7*Q11-Q1*Q12*Q13,
           Q11*Q13-Q9*Q14,
           -Q11*Q13+Q9*Q14,
           -Q5*Q6*Q7+Q1*Q8*Q14,
           Q5*Q6*Q7-Q1*Q8*Q14,
           Q2*Q3*Q8*Q9+Q2*Q4*Q6*Q11-2*Q1*Q5*Q8*Q11-Q2^2*Q9*Q11-Q2*Q3*Q6*Q13-Q1*Q3* Q8*Q14+Q1*Q2*Q11*Q14+Q1*Q5*Q6*Q15+Q1*Q3*Q8*Q15,
           Q5^2*Q6^2*Q13-Q3^2*Q8^2*Q13+Q2*Q5* Q6*Q11*Q13+Q2^2*Q11^2*Q13-Q2*Q3*Q6*Q12*Q13+Q1*Q5*Q8*Q11*Q14-Q1*Q5*Q6*Q12*Q14-Q2*Q3*Q6*Q13*Q14+Q1*Q5*Q8*Q11*Q15-Q1*Q2*Q11*Q12*Q15-Q2*Q3*Q6*Q13*Q15+Q1*Q3*Q8*Q14*Q15,
           Q4*Q5*Q6*Q9*Q12-Q2*Q4*Q9*Q11*Q12-Q5^2*Q6*Q9*Q13+Q2*Q5*Q9*Q11*Q13+ Q1*Q5*Q11*Q12*Q13-Q3^2*Q8*Q13^2+Q4*Q5*Q6*Q9*Q15-Q2*Q5*Q9^2*Q15-Q1*Q5*Q9*Q12*Q15+2* Q2*Q3*Q9*Q13*Q15-Q1*Q5*Q11*Q13*Q15+Q1*Q3*Q13*Q14*Q15-Q1*Q5*Q9*Q15^2,
           Q3*Q5*Q6*Q7-Q3^2*Q7*Q8+Q2*Q3*Q8*Q9-2*Q1*Q5*Q8*Q11-Q1*Q5*Q6*Q12+Q1*Q3* Q8*Q12-Q2*Q3*Q6*Q13+Q1*Q5*Q6*Q14+Q1*Q2*Q11*Q14,
           Q5^2*Q7*Q8*Q12*Q13+Q4*Q5*Q8^2*Q12*Q13+Q2*Q4*Q8*Q12^2*Q13-2*Q5^2*Q8^2*Q13^2+ Q2^2*Q12^2*Q13^2+Q5^2*Q7^2*Q12*Q14-Q4^2*Q8^2*Q12*Q14-Q2*Q4*Q7*Q12^2*Q14-2*Q5^2*Q7*Q8*Q13*Q14+ Q4*Q5*Q8^2*Q13*Q14-Q2*Q5*Q8*Q13^2*Q14+Q4* Q5*Q7*Q8*Q14^2-Q5^2*Q7^2*Q12*Q15-Q4*Q5* Q7*Q8*Q12*Q15+Q2*Q4*Q7*Q12^2*Q15+2*Q5^2*Q7* Q8*Q13*Q15-3*Q2*Q5*Q7*Q12*Q13*Q15-Q2*Q4* Q8*Q12*Q13*Q15+2*Q2*Q5*Q8*Q13^2*Q15+Q5^2* Q7^2*Q14*Q15,
           -Q4*Q8*Q12-Q5*Q8*Q13+Q5*Q7*Q14+Q4*Q8*Q14,
           Q3*Q6*Q7-Q1*Q8*Q11,
           -Q5^2*Q6*Q9^2-Q3*Q5*Q8*Q9^2+2*Q4*Q5*Q6*Q9*Q11-Q4^2*Q6*Q11^2+Q1*Q5*Q9*Q11*Q12-Q1* Q4*Q11^2*Q12-Q3*Q5*Q6*Q9*Q13+Q3*Q4*Q6* Q11*Q13+Q3*Q4*Q6*Q9*Q14-Q1*Q4*Q11^2*Q14+Q1*Q5*Q9*Q11*Q15,
           Q5*Q8*Q11*Q13-Q5*Q6*Q12*Q13-Q2*Q11*Q12* Q13-Q5*Q8*Q9*Q14+Q3*Q8*Q13*Q14+Q4*Q6* Q14^2+Q1*Q12*Q14^2-Q5*Q8*Q9*Q15-Q5*Q6*Q13*Q15+Q3*Q8*Q13*Q15+Q2*Q11*Q13*Q15+Q4*Q6*Q14*Q15-Q3*Q7*Q14*Q15+Q1*Q12*Q14*Q15-Q1*Q14^2*Q15-Q1*Q14*Q15^2,
           Q8*Q9*Q12-Q7*Q11*Q15,
           -Q4*Q6*Q8*Q12+Q3*Q7*Q8*Q12-Q1*Q8*Q12^ 2+Q5*Q6*Q8*Q13-Q3*Q8^2*Q13-Q2*Q8*Q11* Q13+Q2*Q6*Q12*Q13-Q5*Q6*Q7*Q14+Q3*Q7*Q8*Q14-Q1*Q8*Q12*Q14+Q2*Q7*Q11*Q15+Q1*Q8*Q14*Q15,
           -Q3*Q5^2*Q7^2+Q3*Q4*Q5*Q7*Q8+2* Q2*Q3*Q5*Q7*Q13-Q2*Q3*Q4*Q8*Q13-Q2^2* Q4*Q11*Q13+Q1*Q2*Q5*Q12*Q13-Q2^2*Q3*Q13^ 2+Q2*Q3*Q4*Q7*Q14-Q1*Q5^2*Q7*Q14-Q1* Q5^2*Q7*Q15+Q1*Q2*Q5*Q13*Q15,
           Q5*Q8*Q9-Q1*Q12^2,
           -Q5*Q6*Q9*Q12+Q4*Q6*Q11*Q12-Q3*Q7*Q11* Q12+Q1*Q11*Q12^2+Q5*Q6*Q11*Q13-Q3*Q8*Q11*Q13-Q2*Q11^2*Q13-Q3*Q7*Q11*Q14+Q1*Q11* Q12*Q14+Q3*Q6*Q13*Q14+Q3*Q8*Q9*Q15-Q1*Q11*Q14*Q15,
           Q5*Q7*Q11-Q4*Q8*Q11+Q4*Q6*Q12+Q5*Q6*Q13-Q3*Q8*Q13+Q1*Q12*Q14+Q4*Q6*Q15-Q2*Q9*Q15-Q1*Q12*Q15-Q1*Q15^2,
           -Q5*Q7*Q11-Q4*Q8*Q11+Q4*Q6*Q12+Q5*Q6* Q13+Q3*Q8*Q13+Q1*Q12*Q14+Q4*Q6*Q15-Q2* Q9*Q15-Q1*Q12*Q15-Q1*Q15^2,
           -Q5*Q7*Q11+Q4*Q8*Q11+Q4*Q6*Q12-Q5*Q6* Q13+Q3*Q8*Q13+Q1*Q12*Q14+Q4*Q6*Q15-Q2* Q9*Q15-Q1*Q12*Q15-Q1*Q15^2,
           Q5*Q7*Q11+Q4*Q8*Q11+Q4*Q6*Q12-Q5*Q6* Q13-Q3*Q8*Q13+Q1*Q12*Q14+Q4*Q6*Q15-Q2* Q9*Q15-Q1*Q12*Q15-Q1*Q15^2,
           -Q2*Q5*Q6*Q8*Q9-Q2*Q3*Q8^2*Q9+2*Q1*Q5*Q8^2*Q11+Q2^2*Q8*Q9*Q11+Q2^2*Q6*Q9*Q12+ Q1*Q3*Q8^2*Q14-Q1*Q2*Q8*Q11*Q14-Q1*Q2* Q6*Q12*Q14-Q1*Q5*Q6*Q8*Q15-Q1*Q3*Q8^2* Q15+Q1*Q2*Q6*Q12*Q15,
           -Q4*Q5^2*Q6^2+2*Q3*Q4*Q5*Q6*Q8-Q3^2*Q4*Q8^2+Q2*Q4*Q5*Q6*Q11-Q2*Q3*Q4*Q8*Q11+ Q2*Q3*Q4*Q6*Q12-Q1*Q5^2*Q6*Q12-Q2*Q3^2* Q8*Q13-Q1*Q5^2*Q6*Q14+Q1*Q3*Q5*Q8*Q14+Q1*Q3*Q5*Q8*Q15,
           Q5*Q6*Q8-Q3*Q8^2-Q2*Q8*Q11+Q2*Q6*Q12,
           Q5*Q6*Q8+Q3*Q8^2-Q2*Q8*Q11-Q2*Q6*Q12,
           -Q5*Q6*Q8+Q3*Q8^2+Q2*Q8*Q11-Q2*Q6* Q12,
           -Q5*Q6*Q8-Q3*Q8^2+Q2*Q8*Q11+Q2*Q6* Q12,
           Q5*Q6*Q8-Q3*Q8^2+Q2*Q8*Q11-Q2*Q6*Q12,
           -Q5*Q6*Q7*Q9-Q3*Q7^2*Q11+Q2*Q7* Q9*Q11+Q1*Q7*Q11*Q12+Q3*Q6*Q7*Q13-2*Q1*Q8*Q11*Q13-Q1*Q6*Q12*Q13+Q1*Q8*Q9*Q14+ Q1*Q6*Q13*Q14,
           Q2*Q9*Q12+Q5*Q6*Q13-Q3* Q8*Q13-Q2*Q11*Q13-Q1*Q12*Q14+Q1*Q12*Q15,
           Q3*Q7*Q8*Q12-Q2*Q8*Q9*Q12-Q1*Q8*Q12^2 -Q5*Q6*Q8*Q13-Q3*Q8^2*Q13+Q2*Q8*Q11*Q13+Q2*Q6*Q12*Q13+Q5*Q6*Q7*Q14+Q3*Q7*Q8*Q15-Q2*Q7*Q11*Q15-Q1*Q8*Q12*Q15+Q1*Q8*Q14*Q15,
           Q5*Q8*Q9*Q12-Q5*Q7*Q11*Q12+Q5*Q8*Q11*Q13+Q5*Q6*Q12*Q13-Q2*Q11*Q12*Q13-Q4*Q6*Q12*Q14+Q3*Q7*Q12*Q14-Q1*Q12^2* Q14+Q5*Q6*Q13*Q14-Q2*Q11*Q13*Q14+Q3*Q7*Q14^2-Q1*Q12*Q14^2-Q5*Q8*Q9*Q15-Q5*Q7*Q11* Q15+Q1*Q12*Q14*Q15+Q1*Q14^2*Q15,
           Q2*Q9-Q1*Q14,
           -Q3*Q8^2*Q9^2+Q3*Q7^2*Q11^2-Q3*Q6*Q7*Q9*Q12-Q1*Q7*Q11^2*Q12+Q1*Q6*Q9*Q12^2+Q5*Q6^2 *Q9*Q13-Q4*Q6^2*Q11*Q13-Q3*Q6*Q7*Q11*Q13+Q1*Q6*Q11*Q12*Q13+Q3*Q6^2*Q13^2+Q3*Q6* Q7*Q9*Q14-Q1*Q7*Q11^2*Q14-Q1*Q6*Q11*Q13*Q14+2*Q1*Q8*Q9*Q11*Q15+Q1*Q7*Q11^2*Q15-Q1* Q6*Q9*Q12*Q15-Q1*Q6*Q11*Q13*Q15,
           -Q4*Q5^2*Q6*Q9+Q4^2*Q5*Q6*Q11-Q1*Q5^2*Q9*Q12-Q3*Q4*Q5*Q6*Q13+2*Q1*Q5^2*Q11*Q13- Q1*Q3*Q5*Q12*Q13+Q3*Q4^2*Q6*Q14+Q1*Q3* Q4*Q12*Q14+Q1*Q5^2*Q9*Q15-Q1*Q4*Q5*Q11* Q15-Q1*Q3*Q4*Q14*Q15,
           Q4*Q5*Q8^2*Q9-Q1*Q5*Q8*Q12*Q13-Q2*Q4*Q8*Q9*Q14+Q1*Q5*Q8*Q13*Q14+Q1*Q4*Q8*Q14^2+Q2*Q5*Q7*Q9*Q15+Q1*Q5*Q7*Q12* Q15-Q1*Q5*Q8*Q13*Q15-Q2^2*Q9*Q13*Q15-Q1*Q2*Q12*Q13*Q15-Q1*Q5*Q7*Q14*Q15+Q1*Q2* Q13*Q14*Q15,
           -Q5*Q7^3*Q11+Q4*Q7^2*Q8*Q11-Q5*Q6*Q7^2*Q13+Q2*Q7^2*Q11*Q13+Q2*Q6*Q7*Q13^2 -Q1*Q8^2*Q13^2-Q1*Q7*Q8*Q13*Q15+Q1*Q7^2*Q14*Q15,
           -Q5^2*Q6^2*Q12-Q2*Q3*Q8*Q11*Q12+Q2^2 *Q11^2*Q12+Q5^2*Q6^2*Q14+Q3*Q5*Q6*Q8*Q14+Q3^2*Q8^2*Q14-Q2^2*Q11^2*Q14-Q2*Q3*Q6*Q14^2- Q3^2*Q8^2*Q15+Q2*Q3*Q8*Q11*Q15+Q2*Q3*Q6*Q12*Q15-Q2*Q3*Q6*Q14*Q15,
           Q4*Q8*Q11-Q3*Q7*Q15)


#  Inv = unlist(lapply(Inv, abs))
  return(Inv)
}
