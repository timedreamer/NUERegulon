# This script does the XGBoost for the conserved 353 genes. The script is based on Chia-Yi's original XGBoost pipeline.

# Author: Ji Huang

# 0. Prep -----------------------------------------------------------------

library(here)
library(tictoc)
library(xgboost)
library(data.table)
library(optparse)

# source(here("src", "machine_learning", "001_save_dataset.R"))
tictoc::tic()

data <- readRDS(here("result", "machine_learning", "xgboost",
             "input_data_both_nxt_shoot.RDS"))

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


# 1. Run XGboost ----------------------------------------------------------

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
        
        p=c*i
        test.index <-c(p-1,p)
        testing<-data[test.index,]
        training<-data[-test.index,]
        
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
        
        pred[y,1:c]<- predict(bst.val, dtest)
        obs[y,1:c]<-test.trait
        rmse[y]<- as.numeric(bst.val$evaluation_log[r,3])
        train.rmse[y]<-as.numeric(bst.val$evaluation_log[r,2])
        
        # extract important features
        importance_matrix <- xgb.importance(model = bst.val)
        impt.out[[y]]<-paste(importance_matrix$Feature,importance_matrix$Gain, sep = ",")
    }}

# Calculate COR for each model
pred.mat = matrix(pred, nrow=(jj+1))
obs.mat = matrix(obs, nrow=(jj+1))

COR=vector()

for (i in 1:n){
    O=c(obs.mat[,c(i,i+n)])
    P=c(pred.mat[,c(i,i+n)])
    COR[i]=cor(P,O)
}

mean(COR)   # 0.8257
tictoc::toc()


# 2. Save the result ------------------------------------------------------

genotype<-unlist(strsplit(rownames(data)[seq(2,by=2,32)],"_",fixed=T
))[seq(2,by=3,48)]

genotype = factor(genotype, levels =c("B73xIHP1","B73xILP1","B73xLH82","B73xMo17",
                                      "B73xMo18W","B73xOh7B","B73xPH207","B73xPHG47","B73xPHG84",
                                      "B73", "IHP1","ILP1","LH82","Mo17","PH207","PHG84"))


save(data,obs,pred,rmse,train.rmse,impt.out,k,jj,n,r,c,subsample,colsample,eta,num_parallel_tree,COR,
     file=here("result", "machine_learning", "xgboost",
         "XGBoost.Zm-leaf-conservNxT-TotalNUE-output.RData"))


# 3. Extract important genes  ----------------------------------------------------------------------

library(tidyverse)
library(data.table)

load(file=here("result", "machine_learning", "xgboost",
               "XGBoost.Zm-leaf-conservNxT-TotalNUE-output.RData"))

annotation <- read.csv(here("data", "external_data", "Cheng_NC",
                          "B73-AGPv4.gene-symbol-description-AtHomolog.202005.tsv"),
                     header = F,sep = "\t")
names(annotation)=c("Gene","Symbol","Description","Arabidopsis homolog")

s=16 # number of genotypes, n
t=100 # number of iteration, jj+1
weighted.impt=vector("list",s)

# Convert list to data frame

for (i in 1:s){
    m=(i-1)*t+1
    n=t*i
    
    d2 <- as.data.frame(do.call(rbind,flatten(impt.out[m:n]))) %>% 
        separate(.,V1, into =c("Gene","Importance"), sep=",")
    
    ## convert importance from character to numeric
    d2$Importance=as.numeric(d2$Importance)
    
    ## convert data frame to data table for easy calculation
    setDT(d2)
    d2[,sum(Importance),by=Gene]
    output=d2[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
    weighted.impt[[i]]=paste(output$Gene,output$SUM, sep = ",")
}

# Combine the impt
d3 <- as.data.frame(do.call(rbind,flatten(weighted.impt))) %>% 
    separate(.,V1, into =c("Gene","Importance"), sep=",")

d3$Importance=as.numeric(d3$Importance)

# Calculate the sum of importance scores
setDT(d3)
d3[,sum(Importance),by=Gene]

output=d3[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
nrow(output) #248

output2=merge(output,annotation,by="Gene",all.x = T)[order(-SUM)]

# Count the frequency of each feature
output.freq=d3[,.N,by=Gene][order(-N)]

table(output.freq$N==16)
table(output.freq$N>10)  

output3=merge(output2,output.freq,by="Gene",all.x = T)[order(-SUM)]

both_nxt_shoot <- readr::read_tsv(here("result", "network_comparison",
    "orthofinder2_group", "shoot",
    "both_nxt_shoot.tsv"))
gene_list <- unique(both_nxt_shoot$zma_gene)


SUMzero_genes <- data.frame(Gene = setdiff(gene_list, output3$Gene), SUM = 0)  |> 
left_join(annotation, by = "Gene") |> 
mutate(N = 0)

output4  <- bind_rows(output3, SUMzero_genes) 

# Output important genes
write.table(output4,row.names = F,quote = F,sep = "\t",col.names = F,
            file = here("result", "machine_learning", "xgboost",
                        "conservNxT_ZmTotalNUE.XGBoost-importantgene-Athomolog-frequency.tsv")
)

