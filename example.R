rm(list=ls())

source("functions.R")

## input GE as G, regulators as R 
n<-dim(G)[[1]]
pG<-dim(G)[[2]]
pR<-dim(R)[[2]]

### Step I: identify regulatory modules
alpha = 1 ## 1=lasso 0=ridge
nfolds= 5 
theta_hat0<-cv.multi_elnet(G,R,alpha,nfolds) 
## theta_update: estimation of regulatory relationship between GEs and regulators
theta_update<-theta_hat0$theta_esti

theta_list<-lay_list<-list()
S= 2
for(s in 1:S){ 
  theta_s<-theta_update
  lay_list[[s]]<-testbic(theta_update) ## sparse clustering 
  theta_update<-theta_list[[s]]<-update_bic(theta_s,lay_list[[s]]) ## subtract identified module
}


### Step II: integrate omics data based on identified modules
COV_list<-s_d_list<-s_w_list<-list()
idx_list<-idz_list<-list()
pc_list<-list()
s_cutoff<-c()
c=0.8
for(s in 1:S){
  slay<-lay_list[[s]]
  
  idG<-which(apply(slay,2,sum)!=0)
  idR<-which(apply(slay,1,sum)!=0)
  COV<-cov(cbind(G[,idG],R[,idR]))
  
  temp<-svd(COV)
  w<-temp$u
  d<-temp$d
  ind_d<-order(d,decreasing = T)
  s_cutoff[s]<-which.min(cumsum(sort(d/sum(d),decreasing = T))<c)
  
  index_pc<-ind_d[1:s_cutoff[s]] ##
  pc_list[[s]]<-cbind(cbind(G[,idG],R[,idR]))%*%w[,index_pc] 
}

X<- Reduce(cbind,pc_list)

lay_reduce<-Reduce("+",lay_list)
Z<-cbind(G[,which(apply(lay_reduce,2,sum)==0)],R[,which(apply(lay_reduce,1,sum)==0)])

### Step III: Joint model for hierarchical G-E interactions
main<-cbind(X,Z)
# input E for environment factors
inter<-list()
for(i in 1:dim(E)[2]){
  inter[[i]]<-apply(main,2,function(x) {x*E[,i]}) 
}

group_ind<-cbind(s_cutoff[which(s_cutoff>1)],1:length(s_cutoff[which(s_cutoff>1)]))
group_g<-unlist(apply(group_ind,1,function(x){rep(x[2],x[1])}))
group_s<-seq(from=max(group_ind[,2])+1,to=max(group_ind[,2])+dim(Z)[2])
group<-c(group_g,group_s)

p =dim(main)[2]
L =length(inter)

count_group<-as.vector(table(group))
groups_id<-which(count_group>1) 
single_id<-which(count_group==1)
ind_groups<-which(group%in%groups_id) 
ind_single<-which(group%in%single_id)
# input response
# set l1=lambda_1, l2=lambda_2
fit<-joint_model(E,main,inter,response,l1,l2)

