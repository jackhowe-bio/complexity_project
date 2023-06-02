#Just for the fixed effects

GenerateFixedEffectsTable<-function(model=NULL,responses=NULL,link=c("gaussian"),S2var=0,fixed_names=NULL,fixed_del="none",fixed_grp=NULL,fixed_diffdel="none",fixed_diffinc="all",fixed_diff_diffs =NULL,dec_PM=2,pvalues="include")
{ #Load packages
  pacman::p_load(MCMCglmm,coda,openxlsx)
  
  #rename model fixed effects
  model$Sol <- model$Sol[, which(colnames(model$Sol) %in% fixed_names)]
  if(is.null(fixed_names)) {
    colnames(model$Sol)<- colnames(model$Sol)
  } else  {colnames(model$Sol)<-fixed_names
  }
  
  #****************************************************
  #Fixed effects and Differences between levels
  #****************************************************
  #Main effects
  nF=model$Fixed$nfl
  fe1=paste(round(posterior.mode(model$Sol),dec_PM)," (",round(HPDinterval(model$Sol)[,1],dec_PM), ", ",round(HPDinterval(model$Sol)[,2],dec_PM),")",sep="")
  
  #P values using summary.MCMCglmm code
  #Pvalues = option to exclude 
  if(pvalues == "include") {
    fe1_p=pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,1:nF, drop = FALSE] > 0)/dim(model$Sol)[1], 1 - colSums(model$Sol[, 1:nF, drop = FALSE] > 0)/dim(model$Sol)[1]))*2
    fe1_p=round(as.numeric(fe1_p),3)
  } else  {
    if(pvalues == "exclude") {
      fe1_p=rep("-",length(fixed_names))
    } else  {
      fe1_p=pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,1:nF, drop = FALSE] > 0)/dim(model$Sol)[1], 1 - colSums(model$Sol[, 1:nF, drop = FALSE] > 0)/dim(model$Sol)[1]))*2
      fe1_p=round(as.numeric(fe1_p),3)
      fe1_p[pvalues]<-"-"
    }
  }
  
  fe1=data.frame(Fixed_Effects=colnames(model$Sol),Estimates=fe1, pMCMC=fe1_p)
  
  #Do any fixed effects need to be deleted
  if(any(fixed_del == "none")) {
    fe1=fe1
  } else  {
    fe1 = fe1 %>% dplyr::filter(Fixed_Effects %in% fixed_del == F) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  }
  
  #function for calculating differences between all columns of a matrix
  pairwise.diffs <- function(x, nF=1)
  {if(is.matrix(x) & nF>1) {
    #Differences
    # create column combination pairs
    prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
    col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
    
    #pairwise differences 
    result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
    # set colnames
    if(is.null(colnames(x)))
      colnames(x) <- 1:ncol(x)
    colnames(result) <- paste(colnames(x)[col.diffs[, 1]], " vs ", 
                              colnames(x)[col.diffs[, 2]], sep = "")
    ndiffs<-dim(result)[2]
    result<-as.mcmc(result)
    fe2=paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), ", ",round(HPDinterval(result)[,2],dec_PM),")",sep="")
    fe2_p=pmax(0.5/dim(result)[1], pmin(colSums(result[,1:ndiffs, drop = FALSE] > 0)/dim(result)[1], 1 - colSums(result[, 1:ndiffs, drop = FALSE] > 0)/dim(result)[1]))*2
    fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=round(as.numeric(fe2_p),3), check.names=FALSE)
    return(fe2)
  } 
  }
  
  #Differences between fixed effects
  fe2<-pairwise.diffs(model$Sol,nF=model$Fixed$nfl)
  
  #Estimate differences between specified groups of fixed effects
  if(is.null(fixed_grp)) {
    fixed=rbind(fe1,fe2)
  } else  {
    tmp_mat<-matrix(1,nrow=dim(model$Sol)[1],ncol=dim(model$Sol)[2])
    colnames(tmp_mat)<-fixed_grp
    combinations<-pairwise.diffs(x=tmp_mat,nF=2)
    combinations<-combinations %>% dplyr::select(Fixed_Effects) %>% separate(Fixed_Effects,c("comb1","comb2"),sep =" vs ")
    #dplyr::select desired combinations
    fe2<-cbind(fe2,combinations)
    fe2<-fe2 %>% dplyr::filter(as.numeric(comb1) == as.numeric(comb2)) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
    #Combine main effects and differences
    fixed=rbind(fe1,fe2)
  }
  
  #Do any comparisons need to be deleted
  if(any(fixed_diffdel == "none")) {
    fixed=fixed
  } else  {
    fixed = fixed %>% dplyr::filter(Fixed_Effects %in% fixed_diffdel == F) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  }
  
  #Should only specific effects be included
  if(any(fixed_diffinc == "all")) {
    fixed=fixed
  } else  {
    fixed = fixed %>% dplyr::filter(Fixed_Effects %in% c(fixed_names,fixed_diffinc) == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  }
  
  #Should all comparison be deleted
  if(any(fixed_diffdel == "all")) {
    fixed=fe1
  } else  {
    fixed = fixed
  }
  
  #Should any differences of differences be calculated
  if(is.null(fixed_diff_diffs)) {
    fixed=fixed
  } else  {
    #create matrix of diffs
    #function for calculating differences between all columns of a matrix
    pairwise.diffs.mat <- function(x)
    {if(is.matrix(x)) {
      # create column combination pairs
      prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
      col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
      
      #pairwise differences 
      result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
      # set colnames
      if(is.null(colnames(x)))
        colnames(x) <- 1:ncol(x)
      colnames(result) <- paste(colnames(x)[col.diffs[, 1]], " vs ", 
                                colnames(x)[col.diffs[, 2]], sep = "")
      ndiffs<-dim(result)[2]
      result<-as.mcmc(result)
      return(result)
    } 
    }
    
    #Create a matrix of differences
    diffs.mat<-pairwise.diffs.mat(model$Sol)
    
    #Calculate differences between differences
    pairwise.diffs2 <- function(x)
    {if(is.matrix(x)) {
      #Differences
      # create column combination pairs
      prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
      col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
      
      #pairwise differences 
      result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
      # set colnames
      if(is.null(colnames(x)))
        colnames(x) <- 1:ncol(x)
      colnames(result) <- paste(colnames(x)[col.diffs[, 1]], " - ", 
                                colnames(x)[col.diffs[, 2]], sep = "")
      ndiffs<-dim(result)[2]
      result<-as.mcmc(result)
      fe2=paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), ", ",round(HPDinterval(result)[,2],dec_PM),")",sep="")
      fe2_p=pmax(0.5/dim(result)[1], pmin(colSums(result[,1:ndiffs, drop = FALSE] > 0)/dim(result)[1], 1 - colSums(result[, 1:ndiffs, drop = FALSE] > 0)/dim(result)[1]))*2
      fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=round(as.numeric(fe2_p),3), check.names=FALSE)
      return(fe2)
    } 
    }
    diffs.diffs<-pairwise.diffs2(diffs.mat)
    #Select differences that are specified
    diffs.diffs = diffs.diffs %>% dplyr::filter(Fixed_Effects %in% fixed_diff_diffs == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
    #add results to fixed effects
    fixed = rbind(fixed,diffs.diffs)
  }
  
  #Output to excel
  #Fixed effects
  fixedeff <- fixed[!grepl(" vs ",fixed$Fixed_Effects),]
  fixedeff<-data.frame("Fixed Effects"=fixedeff$Fixed_Effects,"Posterior Mode (CI)"=fixedeff$Estimates,"pMCMC"=fixedeff$pMCMC,check.names=FALSE)
  #Fixed differences
  fixeddiff <- fixed[grepl(" vs ",fixed$Fixed_Effects),]
  fixeddiff<-data.frame("Fixed Effects Comparisons"=fixeddiff$Fixed_Effects,"Posterior Mode (CI)"=fixeddiff$Estimates,"pMCMC"=round(as.numeric(fixeddiff$pMCMC),3),check.names=FALSE)
  
  FixedEffectTables<- list(fixedeff,fixeddiff)
  names(FixedEffectTables)<- c('FixedEffects','FixedEffectComparison')
  return(FixedEffectTables)
}
