library(grpreg)
library(MASS)
library(refund)
library(fda)
library(igraph)
library(flexclust)
set.seed(1)
use_sample<-sample(1:1e5,5,replace = FALSE)

n<-100
N <- 200   
t <- seq(0,1,by=1/N)
J <- 2   
rangeval<-c(0,1)
mfunctions<-10
phi_list<-matrix(1,nrow = J+1,ncol = N+1)
for(j in 1:J){
  phi_list[j+1,] <- sqrt(2)*cos(pi*j*t)
}
set.seed(use_sample[1])
setsamplez1<-sample(1:1e5,mfunctions,replace = FALSE)
Z1<-list()
for(mset in 1:mfunctions){
  set.seed(setsamplez1[mset])
  Z1[[mset]]<-mvrnorm(n=n,mu=rep(0,(J+1)),Sigma = diag(rep(c(1:(J+1))^(-2))))
}

X<-matrix(0,nrow = n,ncol = mfunctions*(N+1))
X_list<-list()
for(mset in 1:mfunctions){
  for(i in 1:(n)){
    X[i,((mset-1)*(N+1)+1):(mset*(N+1))]<-t(phi_list)%*%(Z1[[mset]][i,])
  }
  X_list[[mset]]<-X[,((mset-1)*(N+1)+1):(mset*(N+1))]
}

beta<-matrix(0,N+1,mfunctions)
beta[,1]<-t(phi_list)%*%(c(1,1.1,1.2))
beta[,2]<-t(phi_list)%*%(c(1,1.1,1.2))
true_beta<-beta

yita<-rep(0,n)
for (mset in 1:mfunctions) {
  for(i in 1:(n)){
    ina<-rep(0,N)
    for(j in 1:N){
      ina[j]<-0.5*(X_list[[mset]][i,j]%*%beta[j,mset]+(X_list[[mset]][i,j+1]%*%beta[j+1,mset]))
    }
    yita[i]<-yita[i]+sum(ina)*1/N
  }
}

K<-2
mu1<-1
mu2<--1

set.seed(use_sample[2])
e<-rnorm(n,0,0.1)

set.seed(use_sample[3])
mu_group<-rbinom(n,1,0.5)

y<-as.matrix(mu_group*mu1+(1-mu_group)*mu2,n,1)+yita+e

sore_list<-list()
group<-c()
pca_list<-list()

for(mset in 1:mfunctions){
  nowdata<-X_list[[mset]]
  hgtbasis <- create.bspline.basis(c(0,1), length(t) + 4 - 2, 4, t)
  growfdPar <- fdPar(hgtbasis, Lfdobj=2, lambda=10^(-5) )
  datafd <- smooth.basis(t, t(nowdata), growfdPar)$fd
  daytemppcaobj <- pca.fd(datafd, nharm=10, growfdPar)
  varsum<-cumsum(daytemppcaobj$varprop)
  nhar<-min(which(varsum>0.95))
  group<-c(group,rep(mset,nhar))
  daytemppcaobj <- pca.fd(datafd, nharm=nhar, growfdPar)
  sore_list[[mset]]<-daytemppcaobj$scores
  pca_list[[mset]]<-predict(daytemppcaobj$harmonics,seq(0,1,by=1/200)) 
}


group<-factor(x=group)
p_group<-group
Sore<-matrix(0,nrow = n,ncol = length(group))
for(mset in 1:mfunctions){
  Sore[,which(group==mset)]<-sore_list[[mset]]
}



mysubgroup_function<-function(x,y,max_iter,RELTOL,theta,gamma_mcp,lambda_sub,lambda_group,p_group){
  n<-dim(x)[1]
  p<-dim(x)[2]
  
  ST_MCP<-function(deta,lambda,gamma,theta){
    if(abs(deta)<=(gamma*lambda)){
      t<-ST(deta,lambda/theta)/(1-1/(gamma*theta))
    }else{
      t<-deta
    }
    return(t)
  }
  
  ST<-function(t,lambda){
    return(sign(t)*max((abs(t)-lambda),0))
  }
  
  
  eifunction <- function(p,i) {
    e<-rep(0,p)
    e[i]<-1
    return(e)
  }
  deta<-matrix(0,n*(n-1)/2,n)
  flag<-1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      deta[flag,]<- eifunction(n,i)-eifunction(n,j)
      flag<-flag+1
    }
  }
  
  beta<-rep(0,p)
  mu<-rep(0,n)
  yitaij<-rep(0,n*(n-1)/2)
  vij<-rep(0,n*(n-1)/2)
  
  glpre<-cv.grpreg(x, y, p_group, penalty="grLasso",seed = 1)
  beta<-glpre$fit$beta[,glpre$min][-1]
  mu<-y-x%*%beta
  in_deta<-rep(0,n*(n-1)/2)
  flag<-1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      in_deta[flag]<-mu[i]-mu[j]
      yitaij[flag]<- in_deta[flag] 
      flag<-flag+1
    }
  }
   
  beta_new<-beta
  mu_new<-mu
  yitaij_new<-yitaij
  vij_new<-vij
  mysolve<- solve(theta*t(deta)%*%deta+diag(n))
  
  for(update_i in 1:max_iter){
    afit<-grpreg(x, y-mu, p_group, penalty="grMCP",lambda = lambda_group)
    beta_new<-afit$beta[-1]
    mu_new<-mysolve%*%(diag(n)%*%y-x%*%beta_new+theta*t(deta)%*%(yitaij-1/theta*vij))
    in_deta<-rep(0,length(vij))
    flag<-1
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        in_deta[flag]<-mu_new[i]-mu_new[j]+vij[flag]/theta
        yitaij_new[flag]<-ST_MCP(in_deta[flag],lambda_sub,gamma_mcp,theta)
        
        flag<-flag+1
      }
    }
    flag<-1
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        vij_new[flag]<-vij[flag]+theta*(mu_new[i]-mu_new[j]-yitaij_new[flag])
        flag<-flag+1
      }
    }
    r<-deta%*%mu_new-yitaij_new
    nowtol<-norm(r,type="2")
    nowtol
    if(nowtol<RELTOL) break
    beta<- beta_new
    mu<-mu_new
    yitaij<-yitaij_new
    vij<-vij_new
  }
  
  create_adjacency <- function(zeta, n) {
    connected_ix <- which(abs(zeta) <= 1e-5);
    index = t(combn(n,2));
    i <- index[connected_ix,1]
    j <- index[connected_ix,2]
    A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
    A[(j-1)*n + i] <- 1
    return(A)
  }
  Ad_final <- create_adjacency(yitaij, n);
  G_final <- graph.adjacency(Ad_final, mode = 'upper')
  class_group <- components(G_final);
  
  return(list(mu=mu,beta=beta,class_group=class_group))
}


sub_fit<-mysubgroup_function(Sore,y,max_iter=1e4,RELTOL=1e-3,theta=1,gamma_mcp=3,lambda_sub=0.15,lambda_group=0.05,p_group)




 


