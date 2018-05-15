method2color <- function(vec_method) {
    sapply(vec_method, function(x) {
        switch(EXPR = x, "phenix"='#1b9e77', "mvn"='#d95f02', 
               "mice_0.3"='#7570b3',
               "mice_0.2"='#e7298a', 
               "mice_0.1"='#66a61e', "mice_All"='#e6ab02')
    })
}

plot_pattern_missingness <- function(data, color=c('#fc8d62','#8da0cb'), 
                                     name="missingness_pattern", 
                                     directory=".", savePdf=FALSE, 
                                     showPlot=TRUE, width=7, height=5,
                                     labels=FALSE) {
    if (labels) labels <- names(data)
    if (savePdf) {
        pdf(paste(directory, "/", name, ".pdf", sep=""), height=height, 
            width=width)
        aggr_plot <- capture.output(aggr(data, col=color, prop=TRUE, numbers=FALSE, 
                          sortVars=TRUE, labels=labels, cex.axis=.7, gap=3, 
                          ylab=c("Combinations"), combined=TRUE, border=NA,
                          bars=TRUE, cex.numbers=0.9, cex.lab=1, 
                          oma = c(8,4,4,2) + 0.1,
                          only.miss=FALSE, ylim=c(0, 0.42)))
        dev.off()
    }
    if (showPlot) {
        aggr_plot <- capture.output(aggr(data, col=color, prop=TRUE, numbers=FALSE, 
                          sortVars=TRUE, labels=labels, cex.axis=.7, cex.lab=1, 
                          gap=3, 
                          ylab=c("Combinations"), combined=TRUE, border=NA,
                          bars=TRUE, cex.numbers=0.9, oma = c(8,4,4,2) + 0.1,
                          only.miss=FALSE, ylim=c(0, 0.42)))
    }
}

plot_correlation_missingness <- function(data, 
                                         color=c('#f1a340','#f7f7f7','#998ec3'),
                                         directory=".", 
                                         name="correlation_missingness",
                                         savePdf=FALSE, showPlot=TRUE,
                                         width=7, height=7,
                                         addgrid.col=NA, labelsize=1
                                         ) {
    if (savePdf) {
        pdf(paste(directory, "/", name, ".pdf", sep=""), height=height, 
            width=width)
        par(xpd=TRUE)
        correlation_plot(data, tl.col='black', method ="color", 
                         col=colorRampPalette(col=color)(100),
                         cl.title = "Spearman correlation",
                         cl.lim = c(-max(abs(data), na.rm=TRUE),
                                    max(abs(data), na.rm=TRUE)), 
                         cl.cex=0.7, tl.cex=labelsize, tl.offset=0.2, 
                         cl.offset=0.2, 
                         cl.align.text='l', is.corr=FALSE, na.label="square", 
                         na.label.col="grey", mar = c(4,6,7,0), 
                         addgrid.col = addgrid.col)
        mtext(text="Phenotype", side=2, line=3, cex=1.1, adj=0.55)
        mtext(text="Missingness", side=1, line=2, cex=1.1, adj=0.5)
        dev.off()
    }
    if (showPlot) {
        par(xpd=TRUE)
        correlation_plot(data, tl.col='black', method ="color", 
                         col=colorRampPalette(col=color)(100),
                         cl.title = "Spearman correlation",
                         cl.lim = c(-max(abs(data), na.rm=TRUE),
                                    max(abs(data), na.rm=TRUE)), 
                         cl.cex=0.3, tl.cex=labelsize, tl.offset=0.2, 
                         cl.offset=0.2, 
                         cl.align.text='l', is.corr=FALSE, na.label="square", 
                         na.label.col="grey", mar = c(4,4,5,0), 
                         addgrid.col = addgrid.col)
        mtext(text="Phenotype", side=2, line=1, cex=0.8, adj=0.5)
        mtext(text="Missingness", side=1, line=2, cex=0.8, adj=0.5)
    }
}

plot_correlation_phenotypes <- function(data_r, data_p, 
                                        color=c('#f1a340','#f7f7f7','#998ec3'),
                                        directory=".", 
                                        name="correlation_pheno",
                                        savePdf=FALSE, showPlot=TRUE,
                                        width=7, height=7,
                                        addgrid.col=NA, labelsize=1) {
    if (savePdf) {
        pdf(paste(directory, "/", name, ".pdf", sep=""), 
            height=height, width=width)
        par(xpd=TRUE)
        correlation_plot(data_r, tl.col='black', method ="ellipse", 
               col=colorRampPalette(col=color)(100),
               order="hclust", insig="blank", p.mat=data_p, 
               tl.cex=labelsize, cl.cex=0.7, tl.offset=0.2, cl.offset=0.2, 
               mar = c(4,6,7,0),
               cl.align.text='l', addgrid.col = addgrid.col)
        dev.off() 
    }
    if (showPlot) {
        par(xpd=TRUE)
        correlation_plot(data_r, tl.col='black', method ="ellipse", 
                         col=colorRampPalette(col=color)(100),
                         order="hclust", insig="blank", p.mat=data_p, 
                         tl.cex=labelsize, cl.cex=0.7, tl.offset=0.2, 
                         cl.offset=0.2, 
                         mar = c(4,6,7,0),
                         cl.align.text='l', addgrid.col = addgrid.col)
    }
}
    
plot_overview_correlation_imputation <- function(dataframe, cutoff=0.95,
                                        directory= ".",
                                        name="correlation_imputation",
                                        savePdf=FALSE, showPlot=TRUE,
                                        width=8, height=5,
                                        text_size=12) {
    color <- method2color(unique(as.character(dataframe$type)))
    ymin <- min(dataframe$value) - 0.003
    
    p <- ggplot()
    p <- p + geom_boxplot(data=dataframe, aes(x=as.numeric(as.factor(type)), 
                                                    group=type,
                                                    y=value, fill=type)) + 
        scale_fill_manual(values=color, name='Predictors') +
        theme_bw() +
        scale_y_continuous(limits=c(ymin, 1.003), expand = c(0,0) ) +  
        theme(axis.title.y = element_text(size=text_size + 2),
              axis.title.x = element_text(size=text_size + 2),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(size=text_size),
              axis.text.y = element_text(size=text_size),
              legend.title = element_text(size=text_size + 2),
              legend.text = element_text(size=text_size),
              panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.minor.y=element_blank(),
              panel.grid.major.y=element_line(colour = 'grey', size = 0.5)) +
        labs(x="Method", y="Pearson Correlation") +
        geom_hline(yintercept=cutoff)
    if(length(unique(dataframe$type)) > 1 && any(grepl("phenix",
                                                       dataframe$type))) {
        nrTypes <- length(unique(dataframe$type))
        nongenetic <- nrTypes -1 
        p <- p + geom_vline(xintercept=(nrTypes - 0.5), linetype = "dashed") +
            scale_x_continuous(breaks = c(nongenetic/2 + 0.5, nrTypes), 
                               limits=c(0.5, nrTypes + 0.5),
                               labels = c("non-genetic", "genetic"), 
                               expand = c(0,0))
    }
    if (savePdf) {
        ggsave(plot=p, file=paste(directory, "/", name, ".pdf", 
                                  sep=""), height=height, 
               width=width, units="in")
    } 
    if (showPlot) print(p)
}

plot_individual_correlation_imputation <- function(dataframe, cutoff=0.95,
                                                 directory= ".",
                                                 name="correlation_imputation",
                                                 savePdf=FALSE,
                                                 width=8, height=5, 
                                                 text_size=12) {
    
    dataframe$pheno_numeric <- as.numeric(as.factor(dataframe$pheno))
    dataframe <- dataframe[order(dataframe$pheno_numeric), ]
    methods <- unique(dataframe$type)
    traits <- length(unique(dataframe$pheno))
    rect_left <- seq(0.5, (traits - 0.5), 2)
    
    ## i ) genetic via phenix
    corPlots <- lapply(methods, function(m) {
        df <- dplyr::filter(dataframe, type==m) 
        color <- method2color(m)
        Traits2Keep <- df$value > cutoff
        xaxis_color <- rep('black', nrow(df))
        xaxis_color[!Traits2Keep] <- 'darkred'
        
        ymin <- min(df$value) - 0.003
        rectangles <- data.frame(xmin = rect_left, 
                                 xmax = rect_left + 1, 
                                 ymin = ymin, 
                                 ymax = 1.003)
        
        p <- ggplot()
        p <- p + geom_rect(data=rectangles, 
                           aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                           fill='gray80', alpha=0.8) +
            geom_point(data=df, 
                       aes(x=pheno_numeric, y=value, 
                           color=type)) +
            scale_color_manual(values=color, guide=FALSE) +
            theme_bw() +
            scale_y_continuous(limits=c(ymin, 1.003), expand = c(0,0) ) +  
            scale_x_continuous(breaks = seq(1, traits, 1), 
                               limits=c(0.5, traits + 0.5),
                               minor_breaks=seq(0.5, traits + 0.5, 1),
                               labels = unique(df$pheno), 
                               expand = c(0,0)) +
            theme(axis.title.y = element_text(size=text_size + 2),
                  axis.title.x = element_text(size=text_size + 2),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, 
                                             size=text_size, 
                                             color=xaxis_color),
                  axis.text.y = element_text(size=text_size),
                  legend.title = element_text(size=text_size + 2),
                  legend.text = element_text(size=text_size),
                  plot.title = element_text(size=text_size + 2),
                  panel.grid.major.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.minor.y=element_blank(),
                  panel.grid.major.y=element_line(colour = 'grey',
                                                  size = 0.5)) +
            labs(x="Phenotype", y="Pearson Correlation", 
                 title=paste("Imputation with", m)) +
            geom_hline(yintercept=cutoff)
        if(savePdf) {
            ggsave(plot=p, file=paste(directory, "/", name, "_", m, ".pdf", 
                                      sep=""), 
                   height=height, width=width, units="in")
        }
        return(p)
    })
    return(corPlots)
}