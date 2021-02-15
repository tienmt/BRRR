
# random data generation
n = 100  # samples
l = 8  # response
p = 12  # predictors
r = 3   # true rank
sig2 = 1 # noise variance

# generate covariance matrix  
rho.x = 0.0  
S = matrix(rho.x, ncol = p, nrow = p);diag(S) = 1
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)

ystar = diag(p)
library(Matrix)


Iters = 200
h =  3/p/n
lam = 3

Est.mala = Pred.mala = nmse.mala = rank.mala = ac.mala =
  Est.lmc = Pred.lmc = nmse.lmc = rank.lmc = 
  Est.rrr = Pred.rrr = nmse.rrr = rank.rrr =
  Est.Gibb = Pred.Gibb = nmse.gibb = rank.gibb = c()
for(ss in 1:100){
  X = matrix(rnorm(n*p),nr=n,nc=p) %*% S.sqrt
  B = matrix(rnorm(p*r),nr=p,nc=r)%*%t(matrix(rnorm(l*r),nr=l,nc=r))# + rnorm(p*l)
  Y = X%*%B  + rnorm(n*l)
  
  XX = t(X)%*%X
  XY = t(X)%*%Y
  M =  solve(XX+0.1*diag(p),XY) 
  a = 0
  # MALA
  Bm = matrix(data=0,nr=p,nc=l)
  burnin = Iters -100
  for(s in 1:Iters){
    tam1 = glmnet::glmnet(M,ystar,
                          family = 'mgaussian',
                          alpha = 0,lambda = lam^2,
                          standardize.response = F,
                          intercept = F)$beta
    tam = M + h*(XY - XX%*%M)/sig2 + h*(p+l+2)*do.call(rbind, lapply(tam1,as.vector)) + sqrt(2*h)*matrix(rnorm(p*l),nr = p)
    
    logdet = determinant(lam^2*diag(p)+tam%*%t(tam))
    pro.tam = -sum((Y-X%*%tam)^2)/(2*sig2) - 0.5*(p+l+2)*logdet$modulus*logdet$sign
    logdet = determinant(lam^2*diag(p)+M%*%t(M))
    pro.M = -sum((Y-X%*%M)^2)/(2*sig2) - 0.5*(p+l+2)*logdet$modulus*logdet$sign
    
    tam2 = glmnet::glmnet(tam,ystar,
                          family = 'mgaussian',
                          alpha = 0,lambda = lam^2,
                          standardize.response = F,
                          intercept = F)$beta
    tran.m = -sum((M-tam-h*(XY - XX%*%tam)/sig2 - h*(p+l+2)*do.call(rbind, lapply(tam2,as.vector)) )^2)/(4*h)
    tran.tam = -sum((tam-M-h*(XY - XX%*%M)/sig2 - h*(p+l+2)*do.call(rbind, lapply(tam1,as.vector)) )^2)/(4*h)
    
    pro.trans = pro.tam+tran.m-pro.M-tran.tam
    if(log(runif(1)) < pro.trans){
      M = tam
      a = a+1
    } 
    if (s>burnin){
      Bm = Bm + M
    } 
  }
  Bm = Bm/(Iters-burnin)
  ac.mala[ss] = a/Iters
  Est.mala[ss] = mean((Bm - B )^2)
  rank.mala[ss] = Matrix::rankMatrix(Bm,.1)[1]
  nmse.mala[ss] = mean((Bm - B )^2)/mean(B^2) 

  xtest = matrix(rnorm(n*p),nr=n,nc=p) %*% S.sqrt
  ytest = xtest%*%B  + rnorm(n*l)
  Pred.mala[ss] = mean((xtest%*%Bm - ytest )^2)
  
  ### LMC ##########################
  M = solve(XX+0.1*diag(p),XY)
  burnin = Iters -100
  Mmean = matrix(data=0,nr=p,nc=l)
  for(s in 1:Iters){
    tam = glmnet::glmnet(M,ystar,
                         family = 'mgaussian',
                         alpha = 0,lambda = lam^2,
                         standardize.response = F,
                         intercept = F)$beta
    M = M + h*(XY - XX%*%M)/sig2 +
      h*(p+l+2)*do.call(rbind, lapply(tam,as.vector)) +
      sqrt(2*h)*matrix(rnorm(p*l),nr = p)
    if (s>burnin){
      Mmean = Mmean + M     #  print(s)
    } 
  }
  Mmean = Mmean/(Iters-burnin)
  Est.lmc[ss] = mean((Mmean - B )^2)
  rank.lmc[ss] = Matrix::rankMatrix(Mmean,.1)[1]
  Pred.lmc[ss] = mean((xtest%*%Mmean - ytest )^2)
  nmse.lmc[ss] = mean((Mmean - B )^2)/mean(B^2) 
  
  # CV.RRR from rrpack package
  fit.rrr <- rrpack::cv.rrr(Y, X, nfold = 3)
  Est.rrr[ss] = mean((coef(fit.rrr) - B )^2)
  Pred.rrr[ss] = mean((xtest%*%coef(fit.rrr) - ytest )^2)
  nmse.rrr[ss] = mean((coef(fit.rrr) - B )^2)/mean(B^2) 
  rank.rrr[ss] = Matrix::rankMatrix(coef(fit.rrr),.1)[1]
  
  #Gibbs sampler
  m=l
  sigma=1
  k = min(p,m)
  a = 1
  b = 1/1000
  lambda = 1/(2*sigma^2)
  Mstep = matrix(data=0,nr=p,nc=k)
  Nstep = matrix(data=0,nr=l,nc=k)
  gamma = rep(b/a,k)
  Bmean = matrix(data=0,nr=p,nc=l)
  L2 = 2*lambda
  
  Nmcmc = 200
  burnin = Nmcmc-100
  for (step in 1:Nmcmc){
    # update M
    NN = t(Nstep)%*%Nstep
    D  = L2*t(Nstep)%*%t(XY)
    for (i in 1:p) for (j in 1:k){
      Msteptrou = Mstep
      Msteptrou[i,j] = 0
      V  = (1/gamma[j])+L2*NN[j,j]*XX[i,i]
      Mstep[i,j] = rnorm(1,(D-L2*(NN%*%t(Msteptrou)%*%XX))[j,i]/sqrt(V),1) /sqrt(V)
    }
    # update N
    V = X%*%Mstep
    V = L2*t(V)%*%V+diag(1/gamma)
    E = eigen(V)
    E = E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors)
    D = L2*t(Mstep)%*%XY
    Nstep = t(E%*%(matrix(data=rnorm(m*k,0,1),nr=k,nc=m)+E%*%D))
    # update gamma
    for (j in 1:k) gamma[j] = 1/rgamma(1,a+(p+m)/2,b+(sum(Mstep[,j]^2)+sum(Nstep[,j]^2))/2)
    if (step>burnin){
      Bmean = Bmean + Mstep%*%t(Nstep)
    } 
  }
  Bmean = Bmean/(Nmcmc-burnin)
  Est.Gibb[ss] = mean((Bmean - B)^2)
  rank.gibb[ss] = Matrix::rankMatrix(Bmean,.1)[1]
  nmse.gibb[ss] = mean((Bmean - B )^2)/mean(B^2)
  Pred.Gibb[ss] = mean((xtest%*%Bmean - ytest )^2)
}
mala = cbind(Est.mala =Est.mala,
             Pred.mala =Pred.mala,
             nmse.mala =nmse.mala,
             rank.mala =rank.mala,
             ac.mala =ac.mala)
lmc = cbind(Est.lmc = Est.lmc,
            Pred.lmc = Pred.lmc,
            nmse.lmc = nmse.lmc,
            rank.lmc = rank.lmc)
cvrrr = cbind(Est.rrr = Est.rrr,
              Pred.rrr = Pred.rrr,
              nmse.rrr = nmse.rrr,
              rank.rrr = rank.rrr)
gibb = cbind(Est.Gibb = Est.Gibb,
             Pred.Gibb = Pred.Gibb,
             nmse.gibb = nmse.gibb,
             rank.gibb = rank.gibb)
save.image(file = '/data2/thetm/LangevinMCforBRRR/model1.rhox05.rda')



colMeans(lmc)
colMeans(mala)
colMeans(cvrrr)
colMeans(gibb)

apply(lmc, 2, sd)
apply(mala, 2, sd)
apply(cvrrr, 2, sd)
apply(gibb, 2, sd)



