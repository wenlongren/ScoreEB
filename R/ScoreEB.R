ScoreEBtest <- function(genfile,famfile,bimfile,phenofile,IniPthreshold,LOD,repetition){
  XX.gen <- as.matrix(fread(file = genfile, sep = "\t",header = FALSE))
  XX.fam <- fread(file = famfile, sep = "\t",header = TRUE)
  XX.bim <- fread(file = bimfile, sep = "\t",header = TRUE)
  Ytotal <- fread(file = phenofile, sep = ",",header = TRUE)
  XX.gen <- t(XX.gen)
  x <- as.bed.matrix(XX.gen, XX.fam, XX.bim)
  standardize(x) <- "p"
  Ytotal <- as.matrix(Ytotal[,-1])
  Xfix <- matrix(1,dim(x)[1],1)
  KM <- (1/dim(XX.gen)[2])*GRM(x)
  nsam <- dim(XX.gen)[1]
  nmaker <- dim(XX.gen)[2]
  retainP <- 0.05/nmaker
  LODthred <- LOD
  iniP <- IniPthreshold
  resLODTotal2 <- numeric()
  repeatn <- repetition
  for(ii in 1:repeatn){
    resEach <- numeric()
    resEtotal <- numeric()
    resLodTot2 <- numeric()
    Yphe <- as.matrix(Ytotal[,ii])
    resScore <- association.test(x = x,Y = Yphe,X = Xfix, method = "lmm", response = "quantitative", test = "score", K = KM)
    initialChose <- resScore[which(resScore[,8]<=iniP),]
    if(dim(initialChose)[1] > 0){
      retainSNP <-  resScore[which(resScore[,8]<=retainP),]
      if(dim(retainSNP)[1] > 0){
        nintial <- dim(initialChose)[1]
        ntime <- as.matrix(rep(ii,nintial))
        initialChose <- cbind(ntime,initialChose)
        diffSNP <- setdiff(initialChose[,4],retainSNP[,3])
        initialChose <- initialChose[match(diffSNP,initialChose[,4]),]
        retainGen <- as.matrix(XX.gen[,retainSNP[,3]])
        bb1 <- ebayes_EM(Xfix,retainGen,as.matrix(Yphe))
        lod1 <- likelihood(Xfix,retainGen,as.matrix(Yphe),bb1)
        retainSave <- cbind(ii,retainSNP,lod1,abs(bb1))
        retainSave <- retainSave[which(retainSave[,10]>=LODthred),]
        if(dim(retainSave)[1] > 0){
          colnames(retainSave) <- c("ntime","chr","pos","id","A1","A2","freq","score","p","lod","lrtb")
        }
      
        choseGene <- XX.gen[,diffSNP]
        bb2 <- ebayes_EM(Xfix,choseGene,as.matrix(Yphe))
        lod2 <- likelihood(Xfix,choseGene,as.matrix(Yphe),bb2)
        resEach <- cbind(initialChose,lod2,abs(bb2))
        resEach <-resEach[which(resEach[,10]>=LODthred),]
        if(dim(resEach)[1] > 0){
          colnames(resEach) <- c("ntime","chr","pos","id","A1","A2","freq","score","p","lod","lrtb")
        }
        if((dim(retainSave)[1] > 0) && (dim(resEach)[1] > 0)){
          resEtotal <- rbind(retainSave,resEach)
        }else if((dim(retainSave)[1] > 0) && (dim(resEach)[1] == 0)){
          resEtotal <- retainSave
        }else if((dim(retainSave)[1] == 0) && (dim(resEach)[1] > 0)){
          resEtotal <- resEach
        }
      
        if(dim(resEtotal)[1] > 0){
          resLODTotal2 <- rbind(resLODTotal2,resEtotal)
        }
      
      }else{
        nintial <- dim(initialChose)[1]
        ntime <- as.matrix(rep(ii,nintial))
        initialChose <- cbind(ntime,initialChose)
        choseGene <- XX.gen[,initialChose[,4]]
        bb <- ebayes_EM(Xfix,choseGene,as.matrix(Yphe))
        lod <- likelihood(Xfix,choseGene,as.matrix(Yphe),bb)
        resEach <- cbind(initialChose,lod,abs(bb))
        resEach <-resEach[which(resEach[,10]>=LODthred),]
        if(dim(resEach)[1] > 0){
          colnames(resEach) <- c("ntime","chr","pos","id","A1","A2","freq","score","p","lod","lrtb")
          resEtotal <- rbind(resEtotal,resEach)
        }
        if(dim(resEtotal)[1] > 0){
          resLODTotal2 <- rbind(resLODTotal2,resEtotal)
        }
      
      }
    }
  }
  return (resLODTotal2)
}
