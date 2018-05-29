imputeData <- function(data, methods= c("phenix", "mice", "mvn"), fulldata, 
                       imputable=NULL, cutoff=0.95, testing=FALSE, verbose=FALSE,
                       kinship=NULL){
    cor_setups <- data.frame(pheno = factor(colnames(data)))
    imp <- list()
    # i) using genetics (phenix) 
    if (any("phenix" %in% methods)) {
        if(verbose) message("Imputing phenotypes via phenix")
        if (!is.null(imputable)) {
            data_phenix <- as.matrix(data[,imputable$phenix$traits2keep])
        } else {
            data_phenix <- data
        }
        imp$phenix <- phenix_impute(Y=data_phenix, as.matrix(kinship))
        if (testing) {
            cor_phenix_data <- cor2matrices(imp$phenix$imp, fulldata)
            cor_setups$phenix <- cor_phenix_data$cor_r
        }
    }
    # ii) using only phenotype correlations
    if (any("mvn" %in% methods)){
        if (!is.null(imputable)) {
            data_mvn <- as.matrix(data[,imputable$mvn$traits2keep])
        } else {
            data_mvn <- data
        }
        if(verbose) message("Imputing phenotypes via mvn")
        imp$mvn <-  mvn_impute(data_mvn, trace=verbose)
        if (testing) {
            cor_mvn_data <- cor2matrices(imp$mvn$imp, fulldata)
            cor_setups$mvn <- cor_mvn_data$cor_r
        }
    }
    if (any("mice" %in% methods)) {
        if(verbose) message("Imputing phenotypes via mice")
        if (!is.null(imputable)) {
            predictors <- imputable$mice$predictors
        } else {
            predictors=NULL
        }
        imp$mice <- mice_impute(data=as.matrix(data), fulldata, 
                                       testing=testing, predictors=predictors)
        if (!is.null(imputable)) {
            imp$mice$imp  <- imp$mice$imp[,imputable$mice$traits2keep]
        }
        if (testing) {
            cor_setups$mice_corrAll <- 
                imp$mice$correlation_impute$cor_Corr0.0$cor_r 
            cor_setups$mice_corr0.1 <- 
                imp$mice$correlation_impute$cor_Corr0.1$cor_r
            cor_setups$mice_corr0.2 <- 
                imp$mice$correlation_impute$cor_Corr0.2$cor_r
            cor_setups$mice_corr0.3 <- 
                imp$mice$correlation_impute$cor_Corr0.3$cor_r
        }
    }
    if (any("mpmm" %in% methods)) {
        if(verbose) message("Imputing phenotypes via mpmm")
        imp$mpmm <- phenix::MPMM_impute(Y=as.matrix(data), 
                                     K=as.matrix(kinship))
        if (testing) {
            cor_mpmm_data <- cor2matrices(imp$mpmm$imp, fulldata)
            cor_setups$mpmm <- cor_mpmm_data$cor_r
        }
    }
    if (testing) {
        cor_setups_melt <- melt(cor_setups, variable.name="setup")
        cor_setups_melt$type <- gsub("corr(.*)", "\\1", cor_setups_melt$setup)
        cor_setups_melt$method <- "non-genetic"
        cor_setups_melt$method[cor_setups_melt$type %in% c("phenix", "mpmm")] <-
            "genetic"
        summaryImpute <- data.frame(median=reshape2::acast(cor_setups_melt, 
                                                           type ~ ., median),
                                     mean=reshape2::acast(cor_setups_melt,
                                                          type ~ ., mean),
                                     sd=reshape2::acast(cor_setups_melt,
                                                        type ~ ., sd)
        )
        colnames(summaryImpute) <- c("median", "mean", "sd")
        return(list(imp=imp, cor=cor_setups_melt, summary=summaryImpute))
    } else{
        return(list(imp=imp))
    }
}


mice_impute <- function(data, fulldata, seed=500, testing=FALSE, verbose=FALSE,
                        maxit=10, m=10, predictors=NULL) 
    {
    ## a) create predictor matrices based in correlations of phenotypes: 
    ## if trait-trait corr > threshold, use as predictor
    ## if no correlation greater then threshold use all traits
    ## design: rows are targets, 0/1 in columns specify whether trait is used as 
    ## predictor or not
    
    if(verbose) printFlag <- TRUE
    if(!verbose) printFlag <- FALSE
    
    if (is.null(predictors)) {
        if(verbose) message("Generating predictor matrices")
        corr_info <- list(corr_info0.0=quickpred(data, mincor=0.0, minpuc=0.2),
                          corr_info0.1=quickpred(data, mincor=0.1, minpuc=0.2),
                          corr_info0.2=quickpred(data, mincor=0.2, minpuc=0.2),
                          corr_info0.3=quickpred(data,mincor=0.3, minpuc=0.2))
        
        ## b) impute with different predictor matrix schemes
        if(verbose) message("Imputation with all traits as predictors")
    
        imputeData_Corr0.0 <- mice(data, m=m, 
                                   predictorMatrix=corr_info$corr_info0.0, 
                                   maxit=maxit, meth='pmm', seed=seed, 
                                   printFlag=printFlag)
        if(verbose) {
            message("Imputation only using predictors with corr > 0.1")
        }
        imputeData_Corr0.1 <- mice(data, m=m, 
                                   predictorMatrix=corr_info$corr_info0.1, 
                                   maxit=maxit,meth='pmm',seed=seed, 
                                   printFlag=printFlag) 
        if(verbose) {
            message("Imputation only using predictors with corr > 0.2")
        }
        imputeData_Corr0.2 <- mice(data, m=m, 
                                   predictorMatrix=corr_info$corr_info0.2, 
                                   maxit=maxit,meth='pmm',seed=seed, 
                                   printFlag=printFlag) 
        if(verbose) {
            message("Imputation only using predictors with corr > 0.3")
        }
        imputeData_Corr0.3 <- mice(data, m=m, 
                                   predictorMatrix=corr_info$corr_info0.3, 
                                   maxit=maxit, meth='pmm',seed=seed, 
                                   printFlag=printFlag) 
    } else {
        if(verbose) {
            message("Imputation using specified predictor matrix")
        }
        imputeData_predictors <- mice(data, m=m, 
                                   predictorMatrix=predictors, 
                                   maxit=maxit, meth='pmm',seed=seed, 
                                   printFlag=printFlag) 
    }  
    if(verbose) {
        message("Combine imputation results")
    }
    if (is.null(predictors)){
        imputed <- list(
            complete_Corr0.0=combineImpute(imputeData_Corr0.0,npheno=ncol(data)),
            complete_Corr0.1=combineImpute(imputeData_Corr0.1,npheno=ncol(data)),
            complete_Corr0.2=combineImpute(imputeData_Corr0.2,npheno=ncol(data)),
            complete_Corr0.3=combineImpute(imputeData_Corr0.3,npheno=ncol(data))
        )
        if(testing) {
            if(verbose) {
                message("Compute correlation between imputation results and
                         known data")
            }
            correlation_impute <- list(
                cor_Corr0.0=cor2matrices(imputed$complete_Corr0.0, fulldata),
                cor_Corr0.1=cor2matrices(imputed$complete_Corr0.1, fulldata),
                cor_Corr0.2=cor2matrices(imputed$complete_Corr0.2, fulldata),
                cor_Corr0.3=cor2matrices(imputed$complete_Corr0.3, fulldata)
            )
            return(list(imp=imputed, correlation_impute=correlation_impute, 
                        corr_info=corr_info))
        } else {
            return(list(imp=imputed))
        }
    } else {
        imputed <- combineImpute(imputeData_predictors,
                                              npheno=ncol(data))
        colnames(imputed) <- colnames(data)
        if(testing) {
            if(verbose) {
                message("Compute correlation between imputation results and
                        known data")
            }
            correlation_impute <- list(
                cor_predictors=cor2matrices(imputed$complete_predictors, 
                                            fulldata))
            return(list(imp=imputed, correlation_impute=correlation_impute))
        } else {
            return(list(imp=imputed))
        }
    }
}

imputableTraits <- function(m, imputed, cutoff) {
    if (grepl("mice", m)) {
        tmp <- selectPredictorsMICE(imputed, cutoff)
        return(list(traits2keep=tmp$traits2keep, predictors=tmp$predictors))
    } else {
        df <- dplyr::filter(imputed$cor, type==m) 
        return(list(traits2keep=df$value > cutoff, predictors=NULL))
    }
}

selectPredictorsMICE <- function(imputed, cutoff){
    tmp <- reshape2::acast(dplyr::filter(imputed$cor, grepl("mice", setup)), 
                           pheno ~ setup )
    corr_info <- imputed$imp$mice$corr_info
    traits2keep <- apply(tmp, 1, max) >= cutoff
    predictors <- colnames(tmp)[apply(tmp, 1, which.max)]
    
    predictorMatrix <- do.call(rbind, lapply(1:length(predictors), 
                                             function(x) {
        if (predictors[x] == "mice_corrAll") {
            tmp <- rep(1, length(predictors))
            tmp[x] <- 0
            return(tmp)
        }
        if (predictors[x] == "mice_corr0.1") {
            return(corr_info$corr_info0.1[x,])
        }
        if (predictors[x] == "mice_corr0.2") {
            return(corr_info$corr_info0.2[x,])
        }
        if (predictors[x] == "mice_corr0.3") {
            return(corr_info$corr_info0.3[x,])
        }
    }))
    colnames(predictorMatrix) <- predictors
    rownames(predictorMatrix) <- predictors
    return(list(traits2keep=traits2keep, predictors=predictorMatrix))
}





mvn_impute  <- function( Y, reltol=1e-4, intercept=TRUE, maxit=100, 
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
    
    Sigma   = .9*var(Yhat)+.1*diag(P)
    
    if( trace ){
        message( sprintf("%10s, %12s\n",  "counter", '% imp change') )
        message( sprintf("%10d, %12.2e\n", 0,         NA) )
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
            Sig.n[n,m,m]  <- phenix_schur( Sigma, m ) ### used in M
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

phenix_schur <- function (X, m) {
    if (min(m) < 0 | max(m) > ncol(X)) 
        stop("m must have entries in 1, 2, ..., ncol(X)")
    X_m_m <- X[-m, -m]
    diag(X_m_m) <- diag(X_m_m) + 10^(-4)
    X[m, m] - X[m, -m, drop = FALSE] %*% solve(X_m_m) %*% 
        t(X[m, -m, drop = FALSE])
}

phenix_impute <- function (Y, K, Q, lam_K, test = FALSE, test.frac = 0.05, 
                           test.cols = 1:ncol(Y), seed = 8473, 
                           quantnorm = FALSE, scale = TRUE, trim = FALSE, 
                           trim.sds = 4, ...) {
    set.seed(seed)
    intercept <- all(Y[, 1] == 1)
    if (test) {
        t0 <- proc.time()[3]
        obs <- which(!is.na(Y))
        if (any(test.cols < 1) | any(test.cols > ncol(Y))) 
            stop("test.cols must be an integer vector with entries in 1, 2, ..., ncol(Y)")
        if (test.frac <= 0 | test.frac >= 1) 
            stop("test.frac must be between 0 and 1")
        if (intercept & (1 %in% test.cols)) {
            if (!all(test.cols == 1:ncol(Y))) 
                warning("Testing intercept not allowed")
            test.cols <- test.cols[-which(test.cols == 1)]
        }
        obs <- intersect(obs, c(matrix(1:(nrow(Y) * ncol(Y)), 
                                       nrow(Y), ncol(Y))[, test.cols]))
        mask <- sample(obs, floor(test.frac * length(obs)))
        Y.input <- Y
        Y[mask] <- NA
    }
    if (trim) {
        outlier.mat <- matrix(0, nrow(Y), ncol(Y))
        for (p in 1:ncol(Y)) {
            y <- Y[, p]
            mu <- mean(y, na.rm = TRUE)
            sig <- sd(y, na.rm = TRUE)
            if (sig == 0 & p == 1) 
                next
            outliers_p <- which(y < mu - trim.sds * sig | y > 
                                    mu + trim.sds * sig)
            if (length(outliers_p) > 0) 
                outlier.mat[outliers_p, p] <- 1
            rm(outliers_p)
        }
        Y_outliers <- Y
        outliers <- which(outlier.mat == 1)
        if (length(outliers) > 0) 
            Y[outliers] <- NA
    }
    else {
        outliers <- NULL
    }
    if (any(colMeans(is.na(Y)) == 1)) 
        stop("Entirely blank phenotypes not allowed")
    if (scale) {
        phenmeans <- apply(Y, 2, mean, na.rm = TRUE)
        phensds <- apply(Y, 2, sd, na.rm = TRUE)
        Y0 <- Y
        if (intercept) {
            Y[, -1] <- ((Y - rep(1, nrow(Y)) %o% phenmeans)/
                            (rep(1,  nrow(Y)) %o% phensds))[, -1]
        }
        else {
            Y <- (Y - rep(1, nrow(Y)) %o% phenmeans)/(rep(1, 
                                                          nrow(Y)) %o% phensds)
        }
    }
    if (quantnorm) 
        Y <- quantnorm(Y)
    if (missing(K)) {
        if (missing(Q) | missing(lam_K)) 
            stop("Either K or its eigen-components, Q and lam_K, must be given")
        K <- Q %*% diag(lam_K) %*% t(Q)
    }
    if (nrow(K) != nrow(Y)) 
        stop("K and Y are non-conformable")
    if (missing(Q) | missing(lam_K)) {
        eig.K <- eigen(K, symmetric = TRUE)
        Q <- eig.K$vec
        lam_K <- eig.K$val
        rm(eig.K)
    }
    obs <- which(rowSums(!is.na(Y)) > 0)
    if (length(obs) == nrow(Y)) {
        eig.K <- list(vec = Q, val = lam_K)
    }
    else {
        eig.K <- eigen(K[obs, obs], symmetric = TRUE)
    }
    vb_out <- phenix_vb(Y = Y[obs, ], Q = eig.K$vec, lam_K = eig.K$val, 
                        ...)
    Yhat <- Y + NA
    Uhat <- Y + NA
    Yhat[obs, ] <- vb_out$Y
    Uhat[obs, ] <- vb_out$U
    if (length(obs) != nrow(Y)) {
        Yhat[-obs, ] <- K[-obs, obs] %*% eig.K$vec %*% diag(1/eig.K$val) %*% 
            t(eig.K$vec) %*% vb_out$U
        Uhat[-obs, ] <- Yhat[-obs, ]
    }
    if (scale) {
        if (intercept) {
            Yhat[, -1] <- (Yhat * (rep(1, nrow(Y)) %o% phensds) + 
                               rep(1, nrow(Y)) %o% phenmeans)[, -1]
            Uhat[, -1] <- (Uhat * (rep(1, nrow(Y)) %o% phensds) + 
                               rep(1, nrow(Y)) %o% phenmeans)[, -1]
        }
        else {
            Yhat <- Yhat * (rep(1, nrow(Y)) %o% phensds) + rep(1, 
                                                               nrow(Y)) %o% phenmeans
            Uhat <- Uhat * (rep(1, nrow(Y)) %o% phensds) + rep(1, 
                                                               nrow(Y)) %o% phenmeans
        }
    }
    if (test) {
        test_out <- list(time = proc.time()[3] - t0, 
                        cor_glob = cor(Yhat[mask], Y.input[mask]),
                        cors = sapply(test.cols, function(p) {
                            mask.p <- intersect(mask, (p - 1) * nrow(Y) + 1:nrow(Y))
                            if (length(mask.p) == 0) return(NA)
                            cor(Yhat[mask.p], Y.input[mask.p])
                         }), mse_glob = mean((Yhat[mask] - Y.input[mask])^2), 
                        mses = sapply(test.cols, function(p) {
                             mask.p <- intersect(mask, (p - 1) * nrow(Y) + 
                                                     1:nrow(Y))
                             if (length(mask.p) == 0) return(NA)
                             mean((Yhat[mask.p] - Y.input[mask.p])^2)
                         }))
        return(list(imp = Yhat, U = Uhat, S = vb_out$S, beta = vb_out$beta, 
                    h2 = vb_out$h2, Q = Q, lam_K = lam_K, vb_out = vb_out, 
                    B = vb_out$B, E = vb_out$E, outliers = outliers, 
                    test = test_out))
    }
    else {
        return(list(imp = Yhat, U = Uhat, S = vb_out$S, beta = vb_out$beta, 
                    h2 = vb_out$h2, Q = Q, lam_K = lam_K, vb_out = vb_out, 
                    B = vb_out$B, E = vb_out$E, outliers = outliers))
    }
}

phenix_vb <- function (Y, Q, lam_K, N = nrow(Y), P = ncol(Y), M = min(c(N, 
                                                           P)), tau = 0, e = P + 5, E.inv = solve(diag(P)/(e - P - 1)), 
          trace = 0, reltol = 1e-08, maxit = 1000, cutoff = 0.001, 
          ...) 
{
    if (any(lam_K < 0)) {
        warning("K has negative eigenvalues:")
        print(c("range:", range(lam_K[lam_K < 0])))
        print(c("number:", length(lam_K[lam_K < 0])))
        print("Setting negative eigenvalues to 0")
        lam_K[lam_K < 0] <- 0
    }
    N <- nrow(Y)
    P <- ncol(Y)
    eprime <- e + N
    fullInds <- which(rowMeans(is.na(Y)) == 0)
    Y.indic <- is.na(Y) + 0
    j <- 0
    done <- rep(0, N)
    if (length(fullInds) > 0) 
        done[fullInds] <- 1
    j2ns <- j2miss <- list()
    j2len <- numeric()
    while (mean(done) < 1) {
        j <- j + 1
        n <- which(done == 0)[1]
        miss.pattern <- Y.indic[n, ]
        row.indic <- apply(Y.indic, 1, function(y) mean((y + 
                                                             miss.pattern)%%2))
        chosen <- which(row.indic == 0)
        j2ns[[j]] <- chosen
        j2len[j] <- length(chosen)
        j2miss[[j]] <- which(miss.pattern == 1)
        done[chosen] <- 1
    }
    n.miss.types <- j
    if (trace > 0) 
        print("starting initialization...")
    if (all(Y[, 1] == 1)) {
        if (trace) 
            print("intercept found")
        out <- mvn_impute(Y[, -1], intercept = FALSE, ...)
        out$imp <- cbind(1, out$imp)
        out$Sigma <- cbind(c(1, rep(0, P - 1)), rbind(0, out$Sigma))
    }
    else {
        out <- mvn_impute(Y, intercept = FALSE, ...)
    }
    mu_Y <- out$imp
    Omega.inv <- out$Sigma
    diag(Omega.inv) <- diag(Omega.inv) + 10^(-4)
    Omega <- solve(Omega.inv)
    svd.Omega <- svd(Omega, LINPACK = TRUE)
    if (n.miss.types == 0) {
        Sigma_Y <- array(0, dim = c(2, P, P))
    }
    else {
        Sigma_Y <- array(0, dim = c(n.miss.types, P, P))
        for (j in 1:n.miss.types) {
            m <- j2miss[[j]]
            Sigma_Y[j, m, m] <- j2len[j] * solve(Omega[m, m])
        }
    }
    out <- svd(mu_Y, LINPACK = TRUE)
    mu_S <- (out$u %*% diag(sqrt(out$d)))[, 1:M]
    mu_beta <- (diag(sqrt(out$d)) %*% t(out$v))[1:M, ]
    V_beta <- V_S <- diag(M)
    svd.V_beta <- svd.V_S <- list(u = diag(M), d = rep(1, M))
    svd.Omega_beta <- svd.Omega
    counter <- 1
    lik.path <- zero.path <- numeric()
    lik.path[counter] <- lik <- phenix_likelihood(Y, Q, lam_K, 
                                                  N, P, M, tau, E.inv, mu_Y, Sigma_Y, Omega, mu_beta, mu_S, 
                                                  eprime, svd.Omega, svd.V_S, svd.V_beta, svd.Omega_beta, 
                                                  j2miss, j2len, n.miss.types)
    zero.path[counter] <- zeros <- sum(sapply(1:M, function(m) mean(abs(mu_S[, 
                                                                             m] %o% mu_beta[m, ]))) < cutoff)
    if (trace > 0) {
        print("initialization done")
        cat(sprintf("%10s, %10s, %10s, %10s", "counter", "loglik", 
                    "Zeroes", "delta l"), "\n")
        cat(sprintf("%10d, %10.3e, %10d, %10.3e \n", counter, 
                    lik, zeros, NA))
    }
    repeat {
        Omega_beta <- Omega
        svd.Omega_beta <- svd.Omega
        V_S <- t(mu_S) %*% mu_S + trp_KS(svd.keep = svd.V_beta, 
                                         lam.kill = 1/lam_K)
        svd.V_S <- svd(V_S, LINPACK = TRUE)
        mu_beta <- KS_dot_vec(svd.C = list(u = svd.Omega$u, d = tau/svd.Omega$d), 
                              svd.A = svd.V_S, B = t(mu_S) %*% mu_Y, inv = TRUE)
        V_beta <- mu_beta %*% Omega %*% t(mu_beta) + trp_KS(lam.kill = tau/svd.Omega$d, 
                                                            svd.keep = svd.V_S)
        svd.V_beta <- svd(V_beta, LINPACK = TRUE)
        mu_S <- KS_dot_vec(svd.C = svd.V_beta, svd.A = list(u = Q, 
                                                            d = 1/lam_K), B = mu_Y %*% Omega %*% t(mu_beta), 
                           inv = TRUE)
        V_S <- t(mu_S) %*% mu_S + trp_KS(svd.keep = svd.V_beta, 
                                         lam.kill = 1/lam_K)
        svd.V_S <- svd(V_S)
        mu_beta <- KS_dot_vec(svd.C = list(u = svd.Omega$u, d = tau/svd.Omega$d), 
                              svd.A = svd.V_S, B = t(mu_S) %*% mu_Y, inv = TRUE)
        mu_Y.d <- mu_Y - mu_S %*% mu_beta
        trp.s <- trp_KS(svd.keep = svd.V_beta, lam.kill = 1/lam_K)
        if (tau != 0) {
            Delta <- trp_KS(svd.keep = svd.Omega_beta, lam.kill = tau/svd.V_S$d)
        }
        else {
            Delta <- solve(Omega_beta) * M
        }
        Omega.inv <- 1/eprime * (t(mu_Y.d) %*% mu_Y.d + apply(Sigma_Y, 
                                                              2:3, sum) + t(mu_beta) %*% trp.s %*% mu_beta + Delta + 
                                     E.inv)
        svd.Omega.inv <- svd(Omega.inv, LINPACK = TRUE)
        svd.Omega <- list(u = svd.Omega.inv$u, d = 1/svd.Omega.inv$d)
        Omega <- svd.Omega$u %*% diag(svd.Omega$d) %*% t(svd.Omega$u)
        U <- mu_S %*% mu_beta
        mu_Y <- Y
        if (n.miss.types != 0) {
            Sigma_Y <- array(0, dim = c(n.miss.types, P, P))
            for (j in 1:n.miss.types) {
                m <- j2miss[[j]]
                o <- (1:P)[-m]
                ns <- j2ns[[j]]
                mu_Y[ns, m] <- U[ns, m] + (Y - U)[ns, o, drop = FALSE] %*% 
                    solve(Omega.inv[o, o]) %*% Omega.inv[o, m, 
                                                         drop = FALSE]
                Sigma_Y[j, m, m] <- j2len[j] * solve(Omega[m, 
                                                           m])
            }
        }
        counter <- counter + 1
        oldlik <- lik
        lik.path[counter] <- lik <- phenix_likelihood(Y, Q, lam_K, 
                                                      N, P, M, tau, E.inv, mu_Y, Sigma_Y, Omega, mu_beta, 
                                                      mu_S, eprime, svd.Omega, svd.V_S, svd.V_beta, svd.Omega_beta, 
                                                      j2miss, j2len, n.miss.types)
        zero.path[counter] <- zeros <- sum(sapply(1:M, function(m) mean(abs(mu_S[, 
                                                                                 m] %o% mu_beta[m, ]))) < cutoff)
        delta <- (lik - oldlik)/abs(lik)
        if (is.na(lik)) {
            save.image("/data/emu/not-backed-up/dahl/na_lik.Rdata")
            stop("na likelihood")
        }
        if (trace > 0) 
            cat(sprintf("%10d, %10.3e, %10d, %10.3e \n", counter, 
                        lik, zeros, delta))
        if (delta < -reltol) 
            if (abs(delta) < 1e-12) {
                warning(paste("KL decreased: delta ll =", delta))
            }
        else {
            stop(paste("KL decreased: delta ll =", delta))
        }
        if (delta < reltol) {
            break
        }
        else if (counter > maxit) {
            warning(paste("Failed Convergence: delta =", delta))
            break
        }
    }
    B <- t(mu_beta) %*% mu_beta
    E <- solve((eprime - P - 1)/eprime * Omega)
    h2 <- diag(B)/(diag(B) + diag(E))
    vb_pars <- list(S = mu_S, beta = mu_beta, Omega = Omega, 
                    Sigma_Y = Sigma_Y, V_beta = V_beta, V_S = V_S, eprime = eprime, 
                    tau = tau, M = M, Q = Q, E.inv = E.inv)
    return(list(Y = mu_Y, U = U, S = mu_S, beta = mu_beta, counter = counter, 
                lik.path = lik.path, zero.path = zero.path, lik = lik, 
                h2 = h2, B = B, E = E, vb_pars = vb_pars, tau = tau, 
                fitted.M = P - zeros))
}