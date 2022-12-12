# We introduce two functions for implementing the CV1 and CV2 cross validation schemes as follows:
# 1. CV1 Cross Validation Scheme
CV1_Splitter<-function(DIR,DATAset,SEED,FOLD_TEST){
  ### Arguments
  # DIR: Full path to the directory where the synergy DATA folders locates
  # DATAset: Merck, Sanger, and etc.
  # SEED: see number for random number generation
  # FOLD_TEST: the fold that is considered as the test fold (a number between 1 and 5)
  ### Synergy data and unique drug pairs
  data<-read.csv(paste0(DIR,"/Data/",DATAset,".csv"))
  DD<-dim(data)[1]
  U<-unique(data[,c("Index1","Index2")])
  D<-dim(U)[1]
  ### Splitting the unique drug pairs
  require(caret)
  set.seed(SEED)
  S<- createFolds(1:D, k = 5, list = TRUE, returnTrain = FALSE)
  ### Assigning folds to the unique drug pairs 
  TEMP<-numeric(D)
  for (i in 1:5){TEMP[S[[i]]]<-i}
  ### Assigning folds to the samples in the drug synergy data 
  FOLDS<-numeric(DD)
  for (i in 1:D){FOLDS[intersect(which(data[,"Index1"]==U[i,1]),which(data[,"Index2"]==U[i,2]))]<-TEMP[i]}
  ### Generating the Train, Test and Validation Indices
  TEST<-sort(which(FOLDS==FOLD_TEST))
  ALL_TRAIN<-which(FOLDS!=FOLD_TEST)
  set.seed(SEED)
  VAL<-sort(sample(ALL_TRAIN,round(length(ALL_TRAIN)/5)))
  TRAIN<-sort(setdiff(ALL_TRAIN,VAL))
  ### Considering both orientations and applying the (-1) to correct the indices for the Python code
  TEST<-(c(TEST,(TEST+DD))-1)
  VAL<-(c(VAL,(VAL+DD))-1)
  TRAIN<-(c(TRAIN,(TRAIN+DD))-1)
  ### Writing the files, which will be used as inputs in the CCSynergy software
  write.csv(data.frame(Index=TEST),paste0(DIR,"/Data/Cross_Validation/Testing_index.csv"),col.names=TRUE,row.names=FALSE)
  write.csv(data.frame(Index=VAL),paste0(DIR,"/Data/Cross_Validation/Validation_index.csv"),col.names=TRUE,row.names=FALSE)
  write.csv(data.frame(Index=TRAIN),paste0(DIR,"/Data/Cross_Validation/Training_index.csv"),col.names=TRUE,row.names=FALSE)
}


# 2. CV2 Cross Validation Scheme
CV2_Splitter<-function(DIR,DATAset,SEED,FOLD_TEST,TINDEX){
  ### Arguments
  # DIR: Full path to the directory where the synergy DATA folders locates
  # DATAset: Merck, Sanger, and etc.
  # SEED: see number for random number generation
  # FOLD_TEST: the fold that is considered as the test fold (a number between 1 and 5)
  # TINDEX: the index of the given tissue
  ### Synergy data and unique drug pairs
  data<-read.csv(paste0(DIR,"/Data/",DATAset,".csv"))
  DD<-dim(data)[1]
  U<-unique(data[,c("Index1","Index2")])
  D<-dim(U)[1]
  ### Splitting the unique drug pairs
  require(caret)
  set.seed(SEED)
  S<- createFolds(1:D, k = 5, list = TRUE, returnTrain = FALSE)
  ### Assigning folds to the unique drug pairs 
  TEMP<-numeric(D)
  for (i in 1:5){TEMP[S[[i]]]<-i}
  ### Assigning folds to the samples in the drug synergy data 
  FOLDS<-numeric(DD)
  for (i in 1:D){FOLDS[intersect(which(data[,"Index1"]==U[i,1]),which(data[,"Index2"]==U[i,2]))]<-TEMP[i]}
  ### Generating the Train, Test and Validation Indices
  TEST<-sort(intersect(which(FOLDS==FOLD_TEST),which(data$Tindex==TINDEX)))
  ALL_TRAIN<-intersect(which(FOLDS!=FOLD_TEST),which(data$Tindex!=TINDEX))
  set.seed(SEED)
  VAL<-sort(sample(ALL_TRAIN,round(length(ALL_TRAIN)/5)))
  TRAIN<-sort(setdiff(ALL_TRAIN,VAL))
  ### Considering both orientations and the applying (-1) to correct the indices for the Python code
  TEST<-(c(TEST,(TEST+DD))-1)
  VAL<-(c(VAL,(VAL+DD))-1)
  TRAIN<-(c(TRAIN,(TRAIN+DD))-1)
  ### Writing the files, which will be used as inputs in the CCSynergy software
  write.csv(data.frame(Index=TEST),paste0(DIR,"/Data/Cross_Validation/Testing_index.csv"),col.names=TRUE,row.names=FALSE)
  write.csv(data.frame(Index=VAL),paste0(DIR,"/Data/Cross_Validation/Validation_index.csv"),col.names=TRUE,row.names=FALSE)
  write.csv(data.frame(Index=TRAIN),paste0(DIR,"/Data/Cross_Validation/Training_index.csv"),col.names=TRUE,row.names=FALSE)
}


### Example 
# In this example, the index files that are on the GitHub repository are generated as follows
# The cross validation scheme is CV1
# DIR: The working directory
# DATAset: Merck dataset is used
# SEED: we used 28
# FOLD_TEST: 1
CV1_Splitter("/home/User/CCSynergy","Merck",28,1)
