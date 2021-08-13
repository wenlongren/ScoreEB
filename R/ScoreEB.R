#2010 EM_Bayes
ebayes_EM<-function(x,z,y)
{
  n<-nrow(z);k<-ncol(z)
  if(abs(min(eigen(crossprod(x,x))$values))<1e-6)
    b<-solve(crossprod(x,x)+diag(ncol(x))*0.01)%*%crossprod(x,y)
  else
    b<-solve(crossprod(x,x))%*%crossprod(x,y)
  v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
  u<-matrix(rep(0,k),k,1)
  v<-matrix(rep(0,k),k,1)
  s<-matrix(rep(0,k),k,1)
  for(i in 1:k)
  {
    zz<-z[,i]
    s[i]<-((crossprod(zz,zz)+1e-100)^(-1))*v0
    u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
    v[i]<-u[i]^2+s[i]
  }
  vv<-matrix(rep(0,n*n),n,n);
  for(i in 1:k)
  {
    zz<-z[,i]
    vv<-vv+tcrossprod(zz,zz)*v[i]
  }
  vv<-vv+diag(n)*v0
  iter<-0;err<-1000;iter_max<-150;err_max<-1e-6
  tau<-0;omega<-0
  while((iter<iter_max)&&(err>err_max))
  {
    iter<-iter+1
    v01<-v0
    v1<-v
    b1<-b
    vi<-solve(vv)
    xtv<-crossprod(x,vi)
    if(ncol(x)==1)
    {
      b<-((xtv%*%x)^(-1))*(xtv%*%y)
    }else
    {
      if(abs(min(eigen(xtv%*%x)$values))<1e-6){
        b<-solve((xtv%*%x)+diag(ncol(x))*0.01)%*%(xtv%*%y)
      }
      else{
        b<-solve(xtv%*%x)%*%(xtv%*%y)
      }
    }
    r<-y-x%*%b
    ss<-matrix(rep(0,n),n,1)
    for(i in 1:k)
    {
      zz<-z[,i]
      zztvi<-crossprod(zz,vi)
      u[i]<-v[i]*zztvi%*%r
      s[i]<-v[i]*(1-zztvi%*%zz*v[i])
      v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
      ss<-ss+zz*u[i]
    }
    v0<-as.numeric(crossprod(r,(r-ss))/n)
    vv<-matrix(rep(0,n*n),n,n)
    for(i in 1:k)
    {
      zz<-z[,i]
      vv<-vv+tcrossprod(zz,zz)*v[i]
    }
    vv<-vv+diag(n)*v0
    err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
    beta<-t(b)
    sigma2<-v0
  }
  return (u)
}

#Multinormal distribution density function
multinormal<-function(y,mean,sigma)
{
  pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
  return (pdf_value)
}
#LOD value test
likelihood<-function(xxn,xxx,yn,bbo)#xxn:fix matrix;xxx:gene matrix;yn:pheno matrix;bbo:gene effect from EM_bayes2010
{
  nq<-ncol(xxx)
  ns<-nrow(yn)
  at1<-0
  ww1<-as.matrix(which(abs(bbo)>1e-5))
  at1<-dim(ww1)[1]
  lod<-matrix(rep(0,nq),nq,1)
  if(at1>0.5)
    ad<-cbind(xxn,xxx[,ww1])
  else
    ad<-xxn
  #if(abs(det(crossprod(ad,ad)))<1e-6)
  if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6)
    bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
  else
    bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
  vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns);
  ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
  
  sub<-1:ncol(ad);
  if(at1>0.5)
  {
    for(i in 1:at1)
    {
      #ij<-which(sub!=sub[i+1])
      ij<-which(sub!=sub[i+ncol(xxn)])
      ad1<-ad[,ij]
      #if(abs(det(crossprod(ad1,ad1)))<1e-6)
      if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6)
        bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
      else
        bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn) 
      vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
      ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
      lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
    }
  }
  return (lod)
}

### Jacobi preconditioner, A=(1/M)*sigma.k2*GG'+sigma.e2*I, Ax=b, tol is maximum tolerance
### miter is the maxinum iteration number of times
PCG <- function(G,b,m.marker,sigma.k2,sigma.e2,tol,miter){
  tau <- c(sigma.k2,sigma.e2)
  ### k is the iteration number of times
  k <- 0
  x <- matrix(0,length(b),1)
  ### r = b - Ax0
  r <- b
  ### p = M^(-1)*r, here is the Jacobi preconditioner, diagonal of M = diag(A)
  M <- 1/(tau[1]*(1/m.marker)*rowSums(G^2) + tau[2])
  p <- M*r
  ### initial value of min.tol, here min.tol = norm(r,"2")
  min.tol <- sqrt(sum(p^2))  
  while((min.tol > tol) && (k < miter)){
    Ap <- tau[1]*(1/m.marker)*tcrossprod(G,crossprod(p,G)) + tau[2]*p
    alpha <- diag(crossprod(r,p))/diag(crossprod(p,Ap))
    x1 <- x + p*alpha
    r1 <- r - Ap*alpha
    min.tol <- sqrt(colSums(r1^2))
    if(min.tol < tol){
      x <- x1
      r <- r1
      p <- p1
      break
    } 
    p1 <- M*r1
    beta <- diag(crossprod(p1,r1))/diag(crossprod(p,r))
    p1 <- p1 + p*beta
    x <- x1
    r <- r1
    p <- p1
    k <- k + 1
  } 
  print(k)
  print(min.tol)
  ### return the solution x=A^(-1)*b
  return (x) 
}

ScoreEB <- function(genofile, phenofile, popfile = NULL, trait.num = 1, B.Moment = 20, tol.pcg = 1e-4, iter.pcg = 100, bin = 100, lod.cutoff = 3.0, dir_out)
{
  t.start <- proc.time()

  geno <- fread(genofile, header = TRUE)
  pheno <- fread(phenofile, header = TRUE)
  pheno <- as.matrix(pheno[,-1])
  n.geno <- dim(geno)[2]
  X <- t(geno[,5:n.geno])
  X.scale <- scale(X, center = TRUE, scale = TRUE)
  n.sample <- dim(X)[1]
  m.marker <- dim(X)[2]
  F.fix <- as.matrix(rep(1, n.sample))
  if(is.null(popfile) == FALSE){
    popstr <- fread(popfile, header = TRUE)
    popstr <- popstr[-1,]
    popstr <- popstr[,-1]
    F.fix <- cbind(F.fix, popstr)
  }
  
  result.total <- NULL
  for(jj in 1:trait.num){
    
    result <- NULL
    result.final <- NULL
    
    ##Center the phenotype, center and scale the genotype
    Y <- as.matrix(pheno[,jj])
    Y.center <- scale(Y, center = TRUE, scale = FALSE) 
    
    ##Moment estimation for variance component
    set.seed(10000)
    B <- B.Moment
    Zb <- matrix(0,n.sample,B)
    for(i in 1:B)
    {
      Zb[,i] <- rnorm(n.sample,0,1)
    }
    XX.Zb <- X.scale%*%(t(X.scale)%*%Zb)
    XX.Fix.Zb <- X.scale%*%(t(X.scale)%*%(F.fix%*%(solve(crossprod(F.fix, F.fix))%*%crossprod(F.fix, Zb))))
    Minus.Tmp <- XX.Zb - XX.Fix.Zb
    
    Zb.Minus <- matrix(0,B,1)
    Minus.Minus <- matrix(0,B,1)
    for(i in 1:B){
      Zb.Minus[i] <- crossprod(as.matrix(Zb[,i]),as.matrix(Minus.Tmp[,i]))
      Minus.Minus[i] <- crossprod(as.matrix(Minus.Tmp[,i]),as.matrix(Minus.Tmp[,i]))
    }
    
    Trace.KV <- (1/B)*(1/m.marker)*sum(Zb.Minus)
    Trace.KVKV <- (1/B)*(1/m.marker)^2*sum(Minus.Minus)
    N.C <- n.sample - dim(F.fix)[2]
    Coef.Var <- matrix(c(Trace.KVKV, Trace.KV, Trace.KV, N.C), 2, 2, byrow = TRUE)
    
    YY <- crossprod(Y.center, Y.center)
    Y.Fix.Y <- t(Y.center)%*%(F.fix%*%(solve(crossprod(F.fix,F.fix))%*%crossprod(F.fix, Y.center)))
    YVY <- as.numeric(YY) - as.numeric(Y.Fix.Y)
    
    XY <- crossprod(X.scale, Y.center)
    X.Fix.Y <- t(X.scale)%*%(F.fix%*%(solve(crossprod(F.fix,F.fix))%*%crossprod(F.fix, Y.center)))
    Minus.Y <- XY - X.Fix.Y
    Y.VKV.Y <- (1/m.marker)*as.numeric(crossprod(Minus.Y, Minus.Y))
    
    Predict.Vec <- as.matrix(c(Y.VKV.Y, YVY))
    var.com <- solve(Coef.Var)%*%Predict.Vec
    
    ##Variance Components
    sigma.k2 <- abs(var.com[1])
    sigma.e2 <- abs(var.com[2])
    
    ##T-score statistics
    if(n.sample <= 1000){
      M0 <- sigma.k2*tcrossprod(X, X)/m.marker + sigma.e2*diag(n.sample)
      M0Y <- solve(M0)%*%Y
      M0F <- solve(M0)%*%F.fix
    }else{
      M0Y <- PCG(X,Y,m.marker,sigma.k2,sigma.e2,tol.pcg,iter.pcg)
      M0F <- PCG(X,F.fix,m.marker,sigma.k2,sigma.e2,tol.pcg,iter.pcg) 
    }
    
    FtM0F <- solve(crossprod(F.fix, M0F))
    FtM0Y <- crossprod(F.fix, M0Y)
    right.part <- M0F%*%FtM0F%*%FtM0Y
    PY <- M0Y - right.part
    
    tmp <- crossprod(X, PY)
    t.score <- 0.5*tmp^2
    p.value <- pchisq(t.score, 1, lower.tail = FALSE)
    result <- cbind(as.matrix(jj, m.marker), as.matrix(1:m.marker), geno[,2:3], as.matrix(t.score), as.matrix(p.value))
    
    ##divide SNPs into bin
    if((m.marker%%bin)==0){
      group <- m.marker/bin
    }else{
      group <- floor(m.marker/bin) + 1
    }
    
    find.bin.max <- NULL
    for(i in 1:(group-1)){
      tmp <- (1+(i-1)*bin):(i*bin)
      max.score <- result[(i-1)*bin + max(which(result[tmp, 5] == max(result[tmp, 5]))),]
      find.bin.max <- rbind(find.bin.max, matrix(max.score,1,))
    }
    tmp.last <- (1+(group-1)*bin):m.marker
    max.score.last <- result[(group-1)*bin + max(which(result[tmp.last, 5] == max(result[tmp.last, 5]))),]
    find.bin.max <- rbind(find.bin.max, matrix(max.score.last,1,))
    
    find.bin.max <- find.bin.max[order(as.numeric(find.bin.max[,5]), decreasing = TRUE),]
    nrow.find.bin.max <- dim(find.bin.max)[1]
    nrow.select <- min(n.sample, nrow.find.bin.max)
    find.bin.max <- find.bin.max[1:nrow.select,]
    
    ##Empirical Bayes and likelihood ratio test
    geno.bayes <- X[,as.numeric(find.bin.max[,2])]
    b.bayes <- ebayes_EM(F.fix,geno.bayes,Y)
    lod <- likelihood(F.fix,geno.bayes,Y,b.bayes)
    
    result.final <- cbind(find.bin.max, as.matrix(b.bayes), as.matrix(lod))
    result.final <- result.final[which(result.final[,8]>=lod.cutoff),]
    result.final[,6] <- as.matrix(pchisq(as.numeric(result.final[,8])*4.605,1,lower.tail = FALSE)) 
    result.final <- cbind(result.final, result.final[,6])
    
    result.total <- rbind(result.total,result.final[,-6])
  }
  colnames(result.total) <- c("Trait", "Id", "Chr", "Pos", "Score", "Beta", "Lod", "Pvalue")

  t.end <- proc.time()
  t.use <- t.end - t.start

  time.use <- as.matrix(c(t.use[[1]], t.use[[2]], t.use[[3]]))
  rownames(time.use) <- c("User", "System", "Elapse")

  write.table(time.use, paste0(dir_out, "ScoreEB.time.csv"), sep = ",", quote = FALSE, row.names = TRUE, col.names = FALSE)
  write.table(result.total, paste0(dir_out, "ScoreEB.Result.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

  return (result.total)
}









