fullPhenotypes <-function(data) {
    noNA_samples <- !apply(data, 1, function(x) any(is.na(x)))
    data_noNA <- pheno[noNA_samples,]
    return(data_noNA)
}

artificialMissingness <- function(data, fulldata, kinship, seed=3422) {
    set.seed(3422)
    NA_samples <- apply(data, 1, function(x) any(is.na(x)))
    data_NA <- data[which(NA_samples),]
    random_samples <- sample(1:nrow(data), nrow(fulldata))
    data_small <- data[random_samples,]
    data_addNA <- fulldata
    data_addNA[is.na(data_small)] <- NA
    kinship <- kinship[random_samples, random_samples]
    return(list(data_addNA=data_addNA, kinship=kinship))
}

correlationMissingness <- function(data) {
    data_missing <- as.data.frame(abs(is.na(data)))
    corrMiss <- t(sapply(1:ncol(data), function(y) {
        sapply(1:ncol(data), function(x) {
            rcorr(data_missing[!is.na(data[,y]),x], 
                  data[!is.na(data[,y]),y])$r[1,2]
        })
    }))
    colnames(corrMiss) <- colnames(pheno)
    rownames(corrMiss) <- colnames(pheno)
    return(corrMiss)
}

correlationPhenotypes <- function(data, type="spearman") {
    data_r <- rcorr(as.matrix(data), type=type)$r
    data_p <- rcorr(as.matrix(data), type=type)$P
    data_n <- diag(rcorr(as.matrix(data), type=type)$n)
    data_padjust <- apply(data_p, 1, p.adjust)
    return(list(r=data_r, p=data_p, n=data_n, padjust=data_padjust))
}