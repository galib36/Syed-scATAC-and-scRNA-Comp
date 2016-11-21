library(scde)
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(genefilter)
library(M3Drop)
library(rtracklayer)
library(prabclus)
library(rGREAT)
library(ggplot2)
library(ComplexHeatmap)
library(Rtsne)

peakAccesibility <- function(mergedPeakFile, peakFolder, peakFilePattern, outputFolder, outputPeakFileName){
    if(length(list.files(path=peakFolder, pattern = paste0(peakFilePattern,'.bed'), full.names=TRUE))>0){
        files = list.files(path=peakFolder, pattern = paste0(peakFilePattern,'.bed'), full.names=TRUE)
        unlink(files)
    }
    formatingToBed(mergedPeakFile, peakFolder, peakFilePattern, outputFolder)
    peakFilePattern = paste0(peakFilePattern,'.bed')
    query = import(mergedPeakFile)
    queryDF <- data.frame(query)
    #queryDF <- queryDF[order(queryDF$name, decreasing=TRUE),]    
    
    totalOverlap <- data.frame(seqnames = queryDF$seqnames, start = queryDF$start, end = queryDF$end)    
    files <- list.files(path=peakFolder, pattern = peakFilePattern, full.names=TRUE)
    cellName <- list.files(path=peakFolder, pattern = peakFilePattern)
    for (i in 1:length(files)){
        subject = import(files[i])
        hits = findOverlaps(query, subject)
        hitsDF <- data.frame(hits)
        cellName[i] <- gsub(peakFilePattern, '', cellName[i])
        totalOverlap[hitsDF$queryHits, cellName[i]] <- 1
        totalOverlap[-hitsDF$queryHits, cellName[i]] <- 0
    }
    outputFile = paste0(outputFolder,outputPeakFileName)
    write.csv(totalOverlap, outputFile, row.names=FALSE)    
}    


formatingToBed <- function(mergedPeakFile, peakFolder, peakFilePattern, outputFolder){
    files <- list.files(path=peakFolder, pattern = peakFilePattern, full.names=TRUE)
    for (i in 1:length(files)){
        narrowPeak <- read.csv(files[i], header=FALSE, sep='\t')
        write.table(narrowPeak[,1:4], paste0(files[i],'.bed'), 
                    row.names=FALSE, sep='\t', col.names=FALSE, quote=FALSE)
    }
}


getJaccardDist <- function(cdBinary){
        
    if(colnames(cdBinary[,2:3])[1] == 'start' && colnames(cdBinary[,2:3])[2] == 'end'){
        SingleCell.Binary <- cdBinary[,4:(dim(cdBinary)[2])]
    }
    else
        SingleCell.Binary <- cdBinary
    
    
    SingleCell.Binary.Jaccard <- jaccard(as.matrix(SingleCell.Binary))
    
    return(SingleCell.Binary.Jaccard)
}


plotCluster <- function(cdBinary, k, groups=NULL, cell.names, ret.val=FALSE, text.label=FALSE){

    if(missing(k)){
        stop("ERROR: Number of groups \"k\" is missing")
    }
    if(missing(ret.val)){
        ret.val = FALSE
    }
    
    SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
    fit <- cmdscale(as.dist(SingleCell.Binary.Jaccard),eig=TRUE, k=k)
 
    if(is.null(groups)){
        df<-data.frame(x=fit$points[,1],y=fit$points[,2], Cell=colnames(SingleCell.Binary.Jaccard))
        p <- ggplot(df, aes_string(x="x",y ="y"))
    }
    else{
        df<-data.frame(x=fit$points[,1],y=fit$points[,2], Cell=colnames(SingleCell.Binary.Jaccard), groups=groups)
        p <- ggplot(df, aes_string(x="x",y ="y", color="groups"))
    }
    
    p <- p + ggtitle("MDS with Jaccard Distance Matrix") + theme(plot.title = element_text(size = 16, face = "bold"))
    p <- p + geom_point(size = 3)
    p <- p + xlab("Corrdinate 1") 
    p <- p + ylab("Corrdinate 2")
    #p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14),
    #legend.text = element_text(size = 14), legend.title = element_text(size = 14))
    if(text.label==TRUE){
        p<-p + geom_text(data=df,aes(label=colnames(SingleCell.Binary.Jaccard)),
                         alpha=0.6,size=3, vjust=1,hjust=0.7,angle=45)
    }        
    
    print(p)
    
    if(ret.val == TRUE)
        return(fit)
}


plotPCAJaccard <- function(cdBinary, k, groups=NULL, cell.names, ret.val=FALSE , text.label=FALSE, title=""){
    
    if(missing(k)){
        stop("ERROR: Number of groups \"k\" is missing")
    }
    if(missing(ret.val)){
        ret.val = FALSE
    }
    
    SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)
    FinalPCAData <- t(SingleCell.Binary.Jaccard)
    PCx=1
    PCy=2
    pcaPRComp <- prcomp(FinalPCAData)
    percentVar <- pcaPRComp$sdev^2/sum(pcaPRComp$sdev^2)
 
    if(is.null(groups)){
        df<-data.frame(PCX=pcaPRComp$x[,PCx],PCY=pcaPRComp$x[,PCy])
        p1<-ggplot(df, aes_string(x="PCX",y ="PCY"))        
    }
    else{
        df<-data.frame(PCX=pcaPRComp$x[,PCx],PCY=pcaPRComp$x[,PCy], groups=groups)
        p1<-ggplot(df, aes_string(x="PCX",y ="PCY", color="groups"))             
    }          
    if(title=="")
        p1<-p1+ggtitle("PCA with Jacard Matrix")
    else
        p1<-p1+ggtitle(title)
    p1<-p1+geom_point(size = 3)
    p1<-p1+xlab(paste("PC",PCx,": ", round(percentVar[PCx] * 100), "% variance"))
    p1<-p1+ylab(paste("PC",PCy,": ", round(percentVar[PCy] * 100), "% variance"))
    if(text.label == TRUE){
        p1<-p1+geom_text(data=df,aes(label=colnames(SingleCell.Binary.Jaccard)),
                       alpha=0.6,size=3, vjust=1,hjust=0.7,angle=45, color="black")
    }
    
    p1<-p1+scale_fill_hue(l=40)
    return(p1)
}

plotHeatmapJaccard <- function(cdBinary, ret.val=FALSE){
    
    SingleCell.Binary.Jaccard <- getJaccardDist(cdBinary)         
    p <- Heatmap(SingleCell.Binary.Jaccard,row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8), name='Jacard')
    print(p)
    if(ret.val == TRUE)
        return(hclust(as.dist(SingleCell.Binary.Jaccard)))
}





getGreatAnalysis <- function(Input, species, p_threshold){
    
    job = submitGreatJob(bed, version = "3.0", species=species)   
    tb = getEnrichmentTables(job)
    
    write.table(tb[[1]][tb[[1]]$Binom_Raw_PValue< p_threshold,],'Great_GO_Molecular_Functions')
    write.table(tb[[2]][tb[[2]]$Binom_Raw_PValue< p_threshold,],'Great_GO_Biological_Process')
    write.table(tb[[3]][tb[[3]]$Binom_Raw_PValue< p_threshold,],'Great_GO_Cellular_Component')
    
    print(job)
    
    res = plotRegionGeneAssociationGraphs(job)    
}


getDiffBind <- function(cdBinary, groups=NULL){

    if (length(levels(groups)) != 2) {
        stop(paste("ERROR: wrong number of levels in the grouping factor (", 
            paste(levels(groups), collapse = " "), "), but must be two.", 
            sep = ""))
    }
    if (is.null(groups)){
        stop("ERROR: groups factor is not provided")
    }
    
    
    SingleCell.Group1.CellNames <- names(groups[groups==levels(groups)[1]])
    SingleCell.Group2.CellNames <- names(groups[groups==levels(groups)[2]])
    
    SingleCell.Group1.Binary <- cdBinary[,SingleCell.Group1.CellNames]
    SingleCell.Group2.Binary <- cdBinary[,SingleCell.Group2.CellNames]

    SingleCell.Group1.ZeroCount <- rowSums(SingleCell.Group1.Binary==0)
    SingleCell.Group2.ZeroCount <- rowSums(SingleCell.Group2.Binary==0)
    SingleCell.Group1.OneCount <- rowSums(SingleCell.Group1.Binary==1)
    SingleCell.Group2.OneCount <- rowSums(SingleCell.Group2.Binary==1)
    
    SingleCell.Group1VsGroup2 <- data.frame(Chr = cdBinary[,1],
                                           Start = cdBinary[,2],
                                           end = cdBinary[,3])

    SingleCell.Group1VsGroup2$group1OneCounts <- SingleCell.Group1.OneCount
    SingleCell.Group1VsGroup2$group2OneCounts <- SingleCell.Group2.OneCount
    SingleCell.Group1VsGroup2$log2Fold <- log2(SingleCell.Group1.OneCount+1) - log2(SingleCell.Group2.OneCount+1)
    pvalue = vector(mode="numeric", length=length(SingleCell.Group2.OneCount))
    nGenes <- nrow(SingleCell.Group1.Binary)
    
    
    for(i in 1:length(SingleCell.Group2.OneCount))
    {
        contingencyTable <- matrix(c(SingleCell.Group1.ZeroCount[i], SingleCell.Group2.ZeroCount[i],
                                     SingleCell.Group1.OneCount[i], SingleCell.Group2.OneCount[i]), ncol=2)
        SingleCell.Group1VsGroup2$pvalue[i] <- fisher.test(contingencyTable)$p.value        
    }
    
    # Bonferroni correction
    SingleCell.Group1VsGroup2$p_adjust <- p.adjust(SingleCell.Group1VsGroup2$pvalue, method="bonferroni")

    write.csv(SingleCell.Group1VsGroup2[order(SingleCell.Group1VsGroup2$pvalue, decreasing=FALSE),], 
              paste0(levels(groups)[1],'_vs_',levels(groups)[2],'.csv'), row.names=FALSE)
    
    return(SingleCell.Group1VsGroup2)    
}


