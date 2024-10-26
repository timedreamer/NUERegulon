# This is for the RANDOM NXTime genes. 

# Author: Ji Huang

library(here)
library(tictoc)
library(xgboost)
library(data.table)
library(optparse)

# source(here("src", "machine_learning", "001_load_dataset.R"))
tictoc::tic()

data <- readRDS(here("result", "machine_learning", "xgboost",
                     "input_data_zma_nxt_shoot.RDS"))

# Define the command line arguments
option_list <- list(
    make_option(c("--eta"), type = "numeric", default = 0.01, 
                help = "Learning rate for XGBoost", metavar = "number"),
    make_option(c("--r"), type = "integer", default = 100, 
                help = "Number of boosting rounds", metavar = "number"),
    make_option(c("--max_depth"), type = "integer", default = 3, 
                help = "Maximum depth of a tree", metavar = "number"),
    make_option(c("--subsample"), type = "numeric", default = 0.6, 
                help = "Subsample ratio of the training instances", metavar = "number"),
    make_option(c("--colsample"), type = "numeric", default = 0.6, 
                help = "Subsample ratio of columns when constructing each tree", metavar = "number")
)

# Parse command line arguments
args <- parse_args(OptionParser(option_list = option_list))

## Setting parameters

n=16 # number of genotypes
c=2  # number of samples per genotype 
# m=331
m=24 ## This is for random MYB-DIV1 genes

### About the run
k=100 # number of iteration
jj=k-1

### Hyperparameters
r=args$r # nrounds
colsample=args$colsample #0.25
eta=args$eta #0.1
subsample=args$subsample # 0.25 #1 
max_depth_user <- args$max_depth

num_parallel_tree=1

### About the output structure
### Need to reset the following otherwise you'll get 'out of subscript' error
y=0 
obs = matrix(nrow = k*n, ncol=2)
pred = matrix(nrow = k*n, ncol=2)
rmse = vector()
train.rmse = vector()
impt.out=vector('list',k*n)

for (i in 1:n){
    for (j in 0:jj){
        set.seed(j)
        gene.index <- sample(2:ncol(data),m)
        
        p=2*i
        test.index <-c(p-1,p)
        testing<-data[test.index,c(1,as.numeric(gene.index))]
        training<-data[-test.index,c(1,as.numeric(gene.index))]
        
        #convert data frame to data table
        setDT(training) 
        setDT(testing)
        
        #using one hard encoding 
        train.trait <- training$trait 
        test.trait <- testing$trait
        new_training <- model.matrix(~.+0,data = training[,-c("trait"),with=F]) 
        new_testing <- model.matrix(~.+0,data = testing[,-c("trait"),with=F])
        
        #preparing matrix 
        dtrain <- xgb.DMatrix(data = new_training,label = train.trait) 
        dtest <- xgb.DMatrix(data = new_testing,label=test.trait)
        watchlist <- list(train=dtrain, test=dtest)
        
        #Run XGBoost
        params <- list(booster = "gbtree", 
                       objective = "reg:linear", 
                       eta=eta, 
                       gamma=150, 
                       max_depth=max_depth_user, 
                       min_child_weight=1, 
                       subsample=subsample,  
                       eval_metric="rmse",
                       colsample_bytree=colsample,
                       verbosity =0) 
        
        set.seed(j)
        bst.val<-xgb.train( params = params, 
                            data = dtrain, 
                            nrounds = r,
                            nfold = 5, 
                            showsd = T, 
                            stratified = T, 
                            print_every_n = 10, 
                            early_stop_round = 5, 
                            watchlist = watchlist,
                            maximize = F,
                            verbose = F)
        
        y=y+1
        
        pred[y,1:2]<- predict(bst.val, dtest)
        obs[y,1:2]<-test.trait
        rmse[y]<- as.numeric(bst.val$evaluation_log[r,3])
        train.rmse[y]<-as.numeric(bst.val$evaluation_log[r,2])
        
        # extract important features
        importance_matrix <- xgb.importance(model = bst.val)
        impt.out[[y]]<-paste(importance_matrix$Feature,importance_matrix$Gain, sep = ",")
    }}


# Calculate COR for each model
COR=vector()

pred.mat = matrix(pred, nrow=(jj+1))
obs.mat = matrix(obs, nrow=(jj+1))
for (i in 1:n){
    O=c(obs.mat[,c(i,i+n)])
    P=c(pred.mat[,c(i,i+n)])
    COR[i]=cor(P,O)
}

COR
mean(COR) #0.6195525
cor.random<-COR
tictoc::toc()

save(data,obs,pred,rmse,train.rmse,impt.out,k,jj,n,r,c,subsample,colsample,eta,num_parallel_tree,COR,
     file=here("result", "machine_learning", "xgboost",
               "XGBoost.Zm-leaf-random331-TotalNUE-output.RData"))


load(file=here("result", "machine_learning", "xgboost",
                   "XGBoost.Zm-leaf-conservNxT-TotalNUE-output.RData"))

wilcox.test(COR, cor.random, paired=T, alternative= c("greater")) 


