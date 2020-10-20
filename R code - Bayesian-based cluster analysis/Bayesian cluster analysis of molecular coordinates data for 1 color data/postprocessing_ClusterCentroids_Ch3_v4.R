source("internal_v4.R")
foldernames=c("Cell 1")
makeplot=TRUE
superplot=FALSE
skeleton=FALSE
 
sapply(foldernames, function(expname){
    nexpname=expname
    if (skeleton){
        nexpname=paste("R_", expname, sep="")
        dir.create(file.path(nexpname))
        file.copy(file.path(paste(expname, "/config.txt", sep="")), paste(nexpname, "/config.txt", sep=""))
        file.copy(file.path(paste(expname, "/sim_params.txt", sep="")), paste(nexpname, "/sim_params.txt", sep=""))
    }

    r = readLines(con=file.path(paste(nexpname, "/config.txt", sep="")))
    get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
    as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
    model=get("model")
    {if (model=="Gaussian(prec)"){
        xlim = as.v(get("xlim"))
        ylim = as.v(get("ylim"))
        histbins = as.v(get("histbins"))
        histvalues = as.v(get("histvalues"))
        if (length(grep("pbackground",r))==0 | length(grep("alpha",r))==0){
            useplabel=FALSE; pb=NULL; alpha=NULL
        }
        else {
            useplabel=TRUE; 
            pb=as.numeric(get("pbackground"))
            alpha=as.numeric(get("alpha"))
        }
    }
     else {stop("Haven't implemented anything else!")}}
    
    all=list.files(expname)
    dirnames=all[file.info(file.path(paste(expname, "/", all, sep="")))$isdir]
    axes=FALSE;
    cex=1/(3*sqrt(length(dirnames)))
    if (makeplot & superplot) {
        ##These settings control image size and resolution
        png(file.path(paste(nexpname, "/together.png",sep="")), width=10, height=10, units="cm", res=1200)
        nrow=ceiling(sqrt(length(dirnames)))
        par(mfrow=c(nrow, nrow))
    }

    res=lapply(dirnames, function(dirname){
        foldername=file.path(paste(expname, "/", dirname, sep=""))
        nfoldername=file.path(paste(nexpname, "/", dirname, sep=""))
        if (skeleton){
            dir.create(nfoldername)  
            file.copy(file.path(paste(foldername, "/ClusterCentroids_Ch3_v1.txt", sep="")), file.path(paste(nfoldername, "/ClusterCentroids_Ch3_v1.txt", sep="")))
        }
        data=read.table(file.path(paste(nfoldername, "/ClusterCentroids_Ch3_v1.txt", sep="")))
        
        pts = data[,1:2]; sds = data[,3];
        if (skeleton){
            file.copy(file.path(paste(foldername, "/r_vs_thresh_ClusterCentroids_Ch3_v1.txt", sep="")), file.path(paste(nfoldername, "/r_vs_thresh_ClusterCentroids_Ch3_v1.txt", sep="")))
        }
        r = read.csv(file.path(paste(nfoldername, "/r_vs_thresh_ClusterCentroids_Ch3_v1.txt",sep="")), header=FALSE, sep="\t")
        
        m = as.matrix(r)
        cs=(m[1,])[-1]
        thr=(m[,1])[-1]
        m = as.matrix(m[2:length(m[,1]),2:length(m[1,])])
        which.maxm <- function(mat){
            indcol <- rep(1:ncol(mat), each=nrow(mat))[which.max(mat)] 
            indrow <- rep(1:nrow(mat), ncol(mat))[which.max(mat)]
            c(indrow, indcol)
        }
        best=which.maxm(m)
        bestcs=cs[best[2]]
        bestthr=thr[best[1]]
        bfile=file.path(paste(foldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
        nbfile=bfile
        if (skeleton){
            dir.create(paste(nfoldername, "/labels", sep=""))    
            nbfile=file.path(paste(nfoldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
            file.copy(bfile, nbfile)
        }
        labelsbest = strsplit(readLines(nbfile),",")[[1]]
        #################################################################################
        # Extra descriptors
        cluster_list=t(nClustersi(labelsbest))
        print(cluster_list)
        radii_list=c()
        nmol_list=c()
        loc_list=c()
        
        for (k in 1:length(cluster_list)){
          
          radii_list=rbind(radii_list, clusterRadiii(pts,labelsbest,cluster_list[k]))
          nmol_list=rbind(nmol_list, molsPerClusteri(labelsbest,cluster_list[k]))
          
          # print("clusterLocation(pts, labelsbest, cluster_list[k])")
          # print(clusterLocation(pts, labelsbest, cluster_list[k]))
          loc_list=rbind(loc_list, clusterLocation(pts, labelsbest, cluster_list[k]))
          
        }
        print(radii_list)
        print(nmol_list)
        print(loc_list)
        
        LIST=cbind(t(cluster_list), loc_list,radii_list, nmol_list)
        
        
        ##saving the file
        write.csv(LIST,file = paste(nfoldername, "/list_ClusterCentroids_Ch3_v1.csv", sep=""), sep=",")
        
        #################################################################################
        
        ##Some summaries
        wfile=file.path(paste(nfoldername, "/summary_ClusterCentroids_Ch3_v1.txt", sep=""))
        cat("The best: clusterscale", bestcs, " thresh", bestthr, "labels.txt\nNumber of clusters:", nClusters(labelsbest), "\nPercentage in clusters: ", percentageInCluster(labelsbest), "%\nMean number of molecules per cluster: ", nMolsPerCluster(labelsbest), "\nMean radius: ", mean(clusterRadii(pts, labelsbest)), sep="", file=wfile)
        ##browser()
        if (makeplot){
            if (!superplot){
                    pdf(file.path(paste(nfoldername, "/plot_ClusterCentroids_Ch3_v1.pdf", sep="")))
                    axes=TRUE
                    cex=1
                }           
            if (dim(data)[2]==4){
                labelstrue=sapply(as.numeric(data[,4]), function(n){if (n==0) paste(runif(1)) else {paste(n)}})

                par(pty="s")
                par(mfrow=c(1,2))
                par(mar=c(4,4,.5, .5))
                par(oma=c(1,1,1,1))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelstrue), sub="True labels", xlab="",ylab="")
                                        #X11()
                par(mar=c(4,4,.5, .5))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelsbest), sub="Estimated",xlab="",ylab="")
                if (!superplot) dev.off()
            }
            else {
                par(pty="s")
                par(mar=c(0,0,0,0))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelsbest), sub="Clustering",xlab="",ylab="", pch=16, cex=cex, axes=axes)
                box()
                if (!superplot) dev.off()
            }
        }
        ##browser()
        list(radii=clusterRadii(pts, labelsbest), nmols=molsPerCluster(labelsbest), nclusters=nClusters(labelsbest), pclustered=percentageInCluster(labelsbest), totalmols=length(labelsbest), reldensity=reldensity(pts, labelsbest, xlim, ylim))
    })

    if (makeplot & superplot) dev.off()
    nmols=c()
    for (i in 1:length(res)){
        nmols=c(nmols, res[[i]]$nmols)
    }

    ## for (i in 1:(length(res)/4)){
    ##   nmols=c(nmols, res[[4*(i-1)+2]])
    ## }
    h=hist(nmols, plot=FALSE)
    pdf(file.path(paste(nexpname, "/nmols_ClusterCentroids_Ch3_v1.pdf", sep="")))
    plot(h, xlab="Number of molecules", ylab="Number of clusters", main="")
    dev.off()
    f=file.path(paste(nexpname, "/nmols_ClusterCentroids_Ch3_v1.txt", sep="")); cat(nmols, file=f, sep=","); cat("\n", file=f, append=TRUE)
    
    radii=c()
    for (i in 1:length(res)){
        radii=c(radii, res[[i]]$radii)
    }

    ##   for (i in 1:(length(res)/4)){
    ##   radii=c(radii, res[[4*(i-1)+1]])
    ## }
    h=hist(radii, plot=FALSE)
    pdf(file.path(paste(nexpname, "/radii_ClusterCentroids_Ch3_v1.pdf", sep="")))
    plot(h, xlab="Cluster radius", ylab="Number of clusters", main="")
    dev.off()
    f=file.path(paste(nexpname, "/radii_ClusterCentroids_Ch3_v1.txt", sep="")); cat(radii, file=f, sep=","); cat("\n", file=f, append=TRUE)


    nclusters=c()
    for (i in 1:length(res)){
        nclusters=c(nclusters, res[[i]]$nclusters)
    }
    
    h=hist(nclusters, plot=FALSE)
    pdf(file.path(paste(nexpname, "/nclusters_ClusterCentroids_Ch3_v1.pdf", sep="")))
    plot(h, xlab="Number of clusters", ylab="Number of regions", main="")
    dev.off()
    f=file.path(paste(nexpname, "/nclusters_ClusterCentroids_Ch3_v1.txt", sep="")); cat(nclusters, file=f, sep=","); cat("\n", file=f, append=TRUE)

    pclustered=c()
    for (i in 1:length(res)){
        pclustered=c(pclustered, res[[i]]$pclustered)
    }
    ## for (i in 1:(length(res)/4)){
    ##   pclustered=c(pclustered, res[[4*(i-1)+4]])
    ## }
    
    h=hist(pclustered, plot=FALSE)
    pdf(file.path(paste(nexpname, "/pclustered_ClusterCentroids_Ch3_v1.pdf", sep="")))
    plot(h, xlab="Percentage clustered", ylab="Number of regions", main="")
    dev.off()
    f=file.path(paste(nexpname, "/pclustered_ClusterCentroids_Ch3_v1.txt", sep="")); cat(pclustered, file=f, sep=","); cat("\n", file=f, append=TRUE)   

    totalmols=c()
    for (i in 1:length(res)){
        totalmols=c(totalmols, res[[i]]$totalmols)
    }
    ## for (i in 1:(length(res)/4)){
    ##   pclustered=c(pclustered, res[[4*(i-1)+4]])
    ## }
    
    h=hist(totalmols, plot=FALSE)
    pdf(file.path(paste(nexpname, "/totalmols_ClusterCentroids_Ch3_v1.pdf", sep="")))
    plot(h, xlab="Total Mols per ROI", ylab="Number of regions", main="")
    dev.off()
    f=file.path(paste(nexpname, "/totalmols_ClusterCentroids_Ch3_v1.txt", sep="")); cat(totalmols, file=f, sep=","); cat("\n", file=f, append=TRUE)   


    reldensity=c()
    for (i in 1:length(res)){
        reldensity=c(reldensity, res[[i]]$reldensity)
    }
    
    h=hist(reldensity, plot=FALSE)
    pdf(file.path(paste(nexpname, "/reldensity_ClusterCentroids_Ch3_v1.pdf", sep="")))
    plot(h, xlab="Total Mols per ROI", ylab="Number of regions", main="")
    dev.off()
    f=file.path(paste(nexpname, "/reldensity_ClusterCentroids_Ch3_v1.txt", sep="")); cat(reldensity, file=f, sep=","); cat("\n", file=f, append=TRUE)

})
