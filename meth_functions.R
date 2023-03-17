
# generate a plot with a slope relating methylation beta scores with
# the clinical variable
effect_plot <- function(PROBE) {
  GENE <- res[which(res$Name==PROBE),"UCSC_RefGene_Name"]
  GENE <- paste0(unique(unlist(strsplit(GENE,";"))),collapse=";")
  HEADER=paste(GENE,PROBE)
  BETAS <- ilogit2(mvals[rownames(mvals) == PROBE,])
  plot(BETAS, ADOS, main = HEADER, ylab = "ADOS", xlab = "Beta Value")
  regl5 <- lm(ADOS ~ BETAS)
  SLOPE=signif(regl5$coefficients[2],3)
  abline(regl5)
  mycor <- cor.test(ADOS, BETAS, method = "pearson")
  mycor_stat <- signif(mycor$estimate,3)
  mycor_p <- signif(mycor$p.value,3)
  SUBHEAD <- paste("R=",mycor_stat," p=",mycor_p," slope=",SLOPE, sep="")
  mtext(SUBHEAD,cex=0.6)

}


# enrichment on genomic compartments, like exons, TSSs, etc
compartment_enrichment <- function(dma) {
  up <- subset(dma,logFC>0 & P.Value<1e-4)
  dn <- subset(dma,logFC<0 & P.Value<1e-4)
  all <- table(unique(dma)$Regulatory_Feature_Group)
  up <- table(unique(up)$Regulatory_Feature_Group)
  dn <- table(unique(dn)$Regulatory_Feature_Group)
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(up,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  rownames(xx)[1] <- "Intergenic"
  xx[,1] = NULL
  colnames(xx) <- c("all","up")
  xx[is.na(xx)] <- 0
  head(xx)
  x=xx$up
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$up)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$up)-x[2], sum(xx$all) - sum(xx$up) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat)
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  up_comp <- xx

  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(dn,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  rownames(xx)[1] <- "Intergenic"
  xx[,1] = NULL
  colnames(xx) <- c("all","dn")
  xx[is.na(xx)] <- 0
  x=xx$dn
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$dn)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$dn)-x[2], sum(xx$all) - sum(xx$dn) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  dn_comp <- xx
  list("up_comp"=up_comp,"dn_comp"=dn_comp)
}

# enrichment on genomic compartments, like exons, TSSs, etc
# but only for the top 1000 probes
compartment_enrichment2 <- function(dma) {
  all <- table(unique(dma)$Regulatory_Feature_Group)
  dma <- head(dma,1000)
  up <- subset(dma,logFC>0 )
  dn <- subset(dma,logFC<0 )
  up <- table(unique(up)$Regulatory_Feature_Group)
  dn <- table(unique(dn)$Regulatory_Feature_Group)
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(up,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  rownames(xx)[1] <- "Intergenic"
  xx[,1] = NULL
  colnames(xx) <- c("all","up")
  xx[is.na(xx)] <- 0
  head(xx)
  x=xx$up
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$up)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$up)-x[2], sum(xx$all) - sum(xx$up) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  up_comp <- xx

  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(dn,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  rownames(xx)[1] <- "Intergenic"
  xx[,1] = NULL
  colnames(xx) <- c("all","dn")
  xx[is.na(xx)] <- 0
  x=xx$dn
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$dn)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$dn)-x[2], sum(xx$all) - sum(xx$dn) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  dn_comp <- xx
  list("up_comp"=up_comp,"dn_comp"=dn_comp)
}

# enrichment analysis but for CpG island contexts incl shores, shelves, etc.
cgi_enrichment <- function(dma) {
  dma$Relation_to_Island <- gsub("N_","",dma$Relation_to_Island)
  dma$Relation_to_Island <- gsub("S_","",dma$Relation_to_Island)
  all <- table(unique(dma)$Relation_to_Island)
  up <- subset(dma,logFC>0 & P.Value < 1e-4)
  dn <- subset(dma,logFC<0 & P.Value < 1e-4)
  up <- table(unique(up)$Relation_to_Island)
  dn <- table(unique(dn)$Relation_to_Island)
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(up,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  xx[,1] = NULL
  colnames(xx) <- c("all","up")
  xx[is.na(xx)] <- 0
  head(xx)
  x=xx$up
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$up)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$up)-x[2], sum(xx$all) - sum(xx$up) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  up_comp <- xx
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(dn,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  xx[,1] = NULL
  colnames(xx) <- c("all","dn")
  xx[is.na(xx)] <- 0
  x=xx$dn
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$dn)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$dn)-x[2], sum(xx$all) - sum(xx$dn) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  dn_comp <- xx
  list("up_comp"=up_comp,"dn_comp"=dn_comp)
}

# enrichment analysis but for CpG island contexts incl shores, shelves, etc.
# Now just looking at the top 1000 probes
cgi_enrichment2 <- function(dma) {
  dma$Relation_to_Island <- gsub("N_","",dma$Relation_to_Island)
  dma$Relation_to_Island <- gsub("S_","",dma$Relation_to_Island)
  all <- table(unique(dma)$Relation_to_Island)
  dma <- head(dma,1000)
  up <- subset(dma,logFC>0 )
  dn <- subset(dma,logFC<0 )
  up <- table(unique(up)$Relation_to_Island)
  dn <- table(unique(dn)$Relation_to_Island)
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(up,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  xx[,1] = NULL
  colnames(xx) <- c("all","up")
  xx[is.na(xx)] <- 0
  head(xx)
  x=xx$up
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$up)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$up)-x[2], sum(xx$all) - sum(xx$up) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  up_comp <- xx
  xx=NULL
  xx <- merge(as.data.frame(all, row.names = 1),as.data.frame(dn,row.names = 1),by=0, all = TRUE)
  rownames(xx) <- xx[,1]
  xx[,1] = NULL
  colnames(xx) <- c("all","dn")
  xx[is.na(xx)] <- 0
  x=xx$dn
  m=xx$all
  n=sum(xx$all)-xx$all
  k=sum(xx$dn)
  xl <- apply(xx,1,function(x) {
    mat <- matrix(c(x[2],x[1]-x[2], sum(xx$dn)-x[2], sum(xx$all) - sum(xx$dn) -x [1] + x[2] ),2,2)
    mat
    fisher.test(mat) 
  })
  xx$OR <- unname(unlist(lapply(X=xl, FUN = function(x) {x$estimate})))
  xx$fisherPval <- unname(unlist(lapply(X=xl, FUN = function(x) {x$p.value})))
  xx$lowerCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[1]]})))
  xx$upperCI <- unname(unlist(lapply(X=xl, FUN = function(x) {x$conf.int[[2]]})))
  dn_comp <- xx
  list("up_comp"=up_comp,"dn_comp"=dn_comp)
}

# make a forest plot from the enrichment data
make_forest_plots <- function(comp) {

  comp_data <- 
    structure(list(
      "mean"  = comp$up_comp$OR , 
      "lower" = comp$up_comp$lowerCI ,
      "upper" = comp$up_comp$upperCI ,
      .Names = c("mean", "lower", "upper"), 
      row.names = c(NA, -11L), 
      class = "data.frame"))

  comp_data <- as.data.frame(comp_data[1:3],row.names = rownames(comp$up_comp) )

  forestplot(comp_data,title = "hypermethylated",
    labeltext = as.list(rownames(comp_data)),
    mean=mean,lower=lower,upper=upper)

comp_data <- 
  structure(list(
    "mean"  = comp$dn_comp$OR , 
    "lower" = comp$dn_comp$lowerCI ,
    "upper" = comp$dn_comp$upperCI ,
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame"))

comp_data <- as.data.frame(comp_data[1:3],row.names = rownames(comp$dn_comp) )

  forestplot(comp_data,title = "hypomethylated",
    labeltext = as.list(rownames(comp_data)),
    mean=mean,lower=lower,upper=upper)
}

# pmes is an enrichment technique
# get the median, mean and 1-sample t-test result for each gene
# pmea_res <- pmea(mval=mval,design=design,sets=head(sets,20),cores=detectCores()/2)
pmea <- function(mval,design,sets,cores=2) {
  fit <- lmFit(mval, design)
  fit <- eBayes(fit)
  top <- topTable(fit,coef=ncol(design),num=Inf, sort.by = "P")
  l <- mclapply(seq(1,length(sets)), function(i) {
    g <- names(sets[i])
    tstats <- top[rownames(top) %in% sets[[i]],"t"]
    myn <- length(tstats)
    mymean <- mean(tstats)
    mymedian <- median(tstats)
    if ( length(tstats) < 2 ) {
      pval=1
    } else {
      wtselfcont <- t.test(tstats)
      pval=wtselfcont$p.value
    }
    c("gene"=g,"nprobes"=myn,"mean"=mymean,"median"=mymedian,
      "P.Value"=pval)
  } , mc.cores=cores)

  df <- do.call(rbind, l)
  rownames(df) <- df[,1]
  df <- df[,-1]
  tmp <- apply(df,2,as.numeric)
  rownames(tmp) <- rownames(df)
  df <- as.data.frame(tmp)
  df$sig <- -log10(df[,4])
  df <- df[order(-df$sig),]
  df$FDR <- p.adjust(df$P.Value)
  out <- list("df"=df,"toptable"=top)
  return(out)
}

# Run the fry test for each gene. This is a more conservative test.
#fry_res <- run_fry(mval=mval,design=design,sets=sets,cores=cores)
run_fry <- function(mval,design,sets,cores=2) {
  split_sets <- split(sets, ceiling(seq_along(sets)/200))
  fry_l <- mclapply(split_sets,function(l) {
    fry(y=mval, index = l, design = design,
      contrast = ncol(design) )
  } , mc.cores=cores )
  fry_res <- do.call(rbind,fry_l)
  rownames(fry_res) <- sub("\\.","@",rownames(fry_res))
  rownames(fry_res) <- sapply(strsplit(rownames(fry_res),"@"),"[[",2)
  fry_res[is.na(fry_res$PValue),"PValue"] <- 1
  fry_res <- fry_res[order(fry_res$PValue),]
  fry_res$FDR <- p.adjust(fry_res$PValue,method="fdr")
  return(fry_res)
}

# main function to perform 1-sample t-test and fry test and merge the results.
# res <- main(mval,design,sets,cores=detectCores()/2)
main <- function(mval,design,sets,cores=2){
  pmea_res <- pmea(mval=mval,design=design,sets=sets,cores=cores)
  pmea_df <- pmea_res[[1]]
  limma_df <- pmea_res[[2]]
  fry_res <- run_fry(mval=mval,design=design,sets=sets,cores=cores)
  m <- merge(pmea_df,fry_res,by=0)
  rownames(m) <- m$Row.names
  m$Row.names = NULL
  m <- m[,c("nprobes","median","PValue","FDR.y")]
  colnames(m) <- c("nprobes","median","PValue","FDR")
  m <- m[order(m$PValue),]
  out <- list("res"=m,"limma_df"=limma_df)
  return(out)
}

# A chart for probe no bias
probe_bias <- function(res) {
  res$sig <- -log10(res$PValue)
  sig <- subset(res,FDR < 0.05)
  plot(res$nprobes,res$sig,log="x",
    pch=19,cex=0.6,
    xlab="no. probes",ylab="-log10(p-value)")
  points(sig$nprobes,sig$sig,col="red",pch=19,cex=0.62)
  SIG = nrow(sig)
  TOT = nrow(res)
  HEADER <- paste(TOT, "genes in total.", SIG, "with FDR<0.05.")
  mtext(HEADER)
}

# A nice volcano plot for the main pathway result
volcano_plot <- function(res) {
  res$sig <- -log10(res$PValue)
  sig <- subset(res,FDR < 0.05)
  plot(res$median,res$sig,
    pch=19,cex=0.6,
    xlab="median t-statistic",ylab="-log10(p-value)")
  points(sig$median,sig$sig,col="red",pch=19,cex=0.62)
  SIG = nrow(sig)
  UP = nrow(subset(sig,median>0))
  DN = nrow(subset(sig,median<0))
  TOT = nrow(res)
  HEADER <- paste(TOT, "genes in total.", SIG, "with FDR<0.05;",DN,"down,",UP,"up")
  mtext(HEADER)
}

# A boxplot of GMEA results
gmea_boxplot <- function(res,sets,n=50) {
  df <- res[[1]]
  limma_df <- res[[2]]
  # smallest pval
  par(mfrow=c(1,2))
  gs <- head(rownames(df),n)
  mysets <- sets[names(sets) %in% gs]
  tstats <- lapply(mysets, function(set) {
    limma_df[rownames(limma_df) %in% set,"t"]
  })
  tstats <- tstats[order(unlist(lapply(tstats,median)))]
  boxplot(tstats,horizontal=TRUE,las=1,
    main="smallest p-val",cex.axis=0.6,
    xlab="t-statistic")
  grid()
  # biggest effect size (median)
  sig <- subset(df,FDR < 0.05)
  gs <- head(rownames(sig[order(-abs(sig$median)),]),n)
  if ( length(gs) >2 ) {
    tstats <- lapply(gs, function(g) {
      df[which(df$genes==g),"tvals"]
    })
    names(tstats) <- gs
    tstats <- tstats[order(unlist(lapply(tstats,median)))]
    boxplot(tstats,horizontal=TRUE,las=1,
      main="biggest effect size(median)",cex.axis=0.6,
      xlab="t-statistic")
    grid()
  } else {
    plot(1)
    mtext("too few significant genes found")
  }
    par(mfrow=c(1,1))
}

# Aggregate mval to genes and examine pathway enrichment
# this actually works on top of main
gsmea <- function(mval,design,probesets,genesets,cores=2) {
  ag <- mclapply(probesets,function(ps) {
    mval <- mval[rownames(mval) %in% ps,,drop=FALSE]
    if (nrow(mval>1)) {
      o <- apply(mval,2,median)
    } else {
      o <- mval
    }
   o
  },mc.cores=cores)
  ag <- do.call(rbind,ag)
  res <- main(mval=ag,design=design,sets=genesets,cores=cores )
}



