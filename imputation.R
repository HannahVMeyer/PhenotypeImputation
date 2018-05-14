imputeData <- function(data, fulldata, methods= c("phenix", "mice", "mvn"), 
                       kinship=NULL){
    cor_setups <- data.frame(pheno = factor(colnames(data)))
    # ii) using only phenotype correlations
    if (any("mvn" %in% methods)){
        ## a) via mvn
        imputeData_mvn <- mvn_impute(as.matrix(data), trace=TRUE)
        cor_mvn_data <- cor2matrices(imputeData_mvn$imp, fulldata)
        cor_setups$mvn <- cor_mvn_data$cor_r
    }
    if (any("mice" %in% methods)) {
        ### b) via mice
        imputeData_mice <- mice_impute(data=as.matrix(data), fulldata)
        cor_setups$mice_corrAll <- imputeData_mice$imp$cor_Corr0.0$cor_r 
        cor_setups$mice_corr0.1 <- imputeData_mice$imp$cor_Corr0.1$cor_r
        cor_setups$mice_corr0.2 <- imputeData_mice$imp$cor_Corr0.2$cor_r
        cor_setups$mice_corr0.3 <- imputeData_mice$imp$cor_Corr0.3$cor_r
    }
    # i) using genetics (phenix) 
    if (any("phenix" %in% methods)) {
        imputeData_phenix <- phenix::phenix(Y=as.matrix(data), 
                                            as.matrix(kinship))
        cor_phenix_data <- cor2matrices(imputeData_phenix$imp, fulldata)
        cor_setups$phenix <- cor_phenix_data$cor_r
    }
    
    cor_setups_melt <- melt(cor_setups, variable.name="setup")
    cor_setups_melt$type <- gsub("corr(.*)", "\\1", cor_setups_melt$setup)
    cor_setups_melt$method <- "non-genetic"
    cor_setups_melt$method[cor_setups_melt$type=="phenix"] <- "genetic"
    return(cor_setups_melt)
}


mice_impute <- function(data, fulldata) {
    ## a) create predictor matrices based in correlations of phenotypes: 
    ## if trait-trait corr > threshold, use as predictor
    ## if no correlation greater then threshold use all traits
    ## design: rows are targets, 0/1 in columns specify whether trait is used as 
    ## predictor or not
    
    corr_info0.0 <- quickpred(data, mincor=0.0, minpuc=0.2)
    corr_info0.1 <- quickpred(data, mincor=0.1, minpuc=0.2)
    corr_info0.2 <- quickpred(data, mincor=0.2, minpuc=0.2)
    corr_info0.3 <- quickpred(data,mincor=0.3, minpuc=0.2)
    
    ## b) impute with different predictor matrix schemes
    imputeData_Corr0.0 <- mice(data,m=20, predictorMatrix=corr_info0.1, 
                               maxit=30, meth='pmm', seed=500)
    imputeData_Corr0.1 <- mice(data, m=20, predictorMatrix=corr_info0.1, 
                               maxit=30,meth='pmm',seed=500) 
    imputeData_Corr0.2 <- mice(data, m=20, predictorMatrix=corr_info0.2, 
                               maxit=30,meth='pmm',seed=500) 
    imputeData_Corr0.3 <- mice(data, m=20, predictorMatrix=corr_info0.3, 
                               maxit=30,meth='pmm',seed=500) 
    
    imputeData <- list(
        complete_Corr0.0=combineImpute(imputeData_Corr0.0,npheno=ncol(data)),
        complete_Corr0.1=combineImpute(imputeData_Corr0.1,npheno=ncol(data)),
        complete_Corr0.2=combineImpute(imputeData_Corr0.2,npheno=ncol(data)),
        complete_Corr0.3=combineImpute(imputeData_Corr0.3,npheno=ncol(data))
    )
    correlation_impute <- list(
        cor_Corr0.0=cor2matrices(complete_Corr0.0, fulldata),
        cor_Corr0.1=cor2matrices(complete_Corr0.1, fulldata),
        cor_Corr0.2=cor2matrices(complete_Corr0.2, fulldata),
        cor_Corr0.3=cor2matrices(complete_Corr0.3, fulldata)
    )
    
    return(list(correlation_impute=correlation_impute,
                imp=imputeData))
}


mvn_impute  <- function( Y, reltol=1e-4, intercept=TRUE, maxit=1e2, 
                         trace=FALSE ){
    
    if( any( colMeans( is.na(Y) ) == 1 ) )
        stop( 'Completely blank columns of Y are not allowed' )
    
    N <- nrow(Y)
    P <- ncol(Y)    
    miss.indices	<- sapply( 1:N, function(n) which( is.na(Y[n,]) ) )
    nmiss		<- sapply( miss.indices, length )
    
    if( sum(nmiss) == 0 ) 
        return(list( Y=Y, Sigma=var(Y) ))
    
    loopers   <- which( nmiss > 0 & nmiss < P )
    emptyInds <- which( nmiss==P )
    fullInds  <- which( nmiss==0 )
    nempty    <- length(emptyInds)
    
    ######## initialization
    if( intercept ){
        mu		= colMeans( Y, na.rm=TRUE )
    } else {
        mu  <- rep( 0, P )
    }
    
    Yhat	= Y
    for(p in 1:P)	##### initialize missing values to column means
        Yhat[which( is.na(Yhat[,p])),p] = mu[p]
    
    Sigma   = .9*var(Yhat)+.1*diag(P) # identity causes weird singularity and var(Yhat) needn't be pd
    
    if( trace ){
        cat( sprintf("%10s, %12s\n",  "counter", '% imp change') )
        cat( sprintf("%10d, %12.2e\n", 0,         NA) )
    }
    
    for( it in 1:maxit ){
        Y.old   <- Yhat
        
        ####### E step
        Yhat		<- Y
        Sig.n	<- array( 0, dim=c(N,P,P) ) ## for M step
        for( n in loopers ){
            m   	    <- miss.indices[[n]]
            o   	    <- (1:P)[-m]
            sigma_o_o <- Sigma[o,o]
            diag(sigma_o_o) <- diag(sigma_o_o) + 10^(-4)
            Yhat[n,m]	    <- mu[m] + Sigma[m,o,drop=FALSE] %*% solve(sigma_o_o) %*% 
                t(Y[n,o,drop=FALSE]-mu[o])
            Sig.n[n,m,m]  <- phenix::schur( Sigma, m ) ### used in M
        }
        
        ####### E step cont'd
        if( nempty > 0 ){
            Yhat[emptyInds,]	<- matrix( mu, nempty, P, byrow=TRUE )
            for( n in emptyInds )
                Sig.n[n,,]	<- Sigma ### used in M
        }
        
        ##### determine convergence
        delta <- 100 * mean(( (Yhat-Y.old)[is.na(Y)] )^2 ) / mean( ( Yhat[is.na(Y)] )^2 )
        if( trace )
            cat( sprintf("%10d, %12.2e\n", it, delta ) )
        if ( delta/100 < reltol )
            break
        
        ####### M step
        if( intercept )
            mu	<- colMeans( Yhat )
        
        mu.mat    <- matrix( mu, N, P, byrow=TRUE )
        S.part.1  <- 1/N * t(Yhat - mu.mat)%*%(Yhat - mu.mat)
        Sigma     <- S.part.1 + 1/N * apply( Sig.n, 2:3, sum )
    }
    return(list(imp=Yhat, Sigma=Sigma))
}

combineImpute <- function(imputelist, npheno) { 
    mat <- complete(imputelist, 'r')
    nimpute <- ncol(mat)/npheno
    sets <- seq(1, ncol(mat), nimpute)
    imputemedian <- sapply(sets, function(s, mat) {
        apply(mat[, s:(s+nimpute-1)], 1, median)
    }, mat=mat)
    return(imputemedian)
}