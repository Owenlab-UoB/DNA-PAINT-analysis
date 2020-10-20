library(splancs)
library(igraph)
 
mcgaussprec <- function(pts, sds, xlim = c(0,1), ylim=c(0,1), psd=function(sd){0},minsd=0.1,maxsd=100, grid=100){
  N = dim(pts)[1]
  fsd <- Vectorize(function(sd){
    wts = 1/(sd^2 + sds^2); tildeN = sum(wts)
    mu = c(sum(wts * pts[,1])/tildeN, sum(wts*pts[,2])/tildeN)
    totdist = sum(c(wts*(pts[,1] - mu[1])^2, wts*(pts[,2]-mu[2])^2))    

    ##x-axis
    log(pnorm(sqrt(tildeN) * (xlim[2]-mu[1])) - pnorm(sqrt(tildeN) * (xlim[1] - mu[1])))+
      ##y-axis
      log(pnorm(sqrt(tildeN) * (ylim[2]-mu[2])) - pnorm(sqrt(tildeN) * (ylim[1] - mu[2])))+
        ##marginalised (with factor 2pi taken from standardisation above)
        -(N-1)*log(2 * pi)+sum(log(wts))-totdist/2+
          ##size of area
          -log(diff(xlim)*diff(ylim))+
            ##cst -- left in standardisation above
            -log(tildeN)+
              ##prior on sd
              psd(sd)
  })
  ##discrete prior:
  x = seq(minsd, maxsd, length=grid)[-1];
  values = fsd(x); dx=x[2]-x[1]; m = max(values); int=sum(exp(values-m))*dx;
  log(int)+m
}


mkcols <- function(labels){
  t=table(labels)
  cnames=names(t[t>1])  
  colors=sample(rainbow(length(cnames)))
  s=sapply(labels, function(l){
    i=which(names(t)==l);
    if (t[i]==1){"grey"}
    else {colors[which(cnames==l)]}
  })
  s  
}

toroid <- function(pts, xlim, ylim, range){
  xd=xlim[2]-xlim[1]; yd=ylim[2]-ylim[1]
  R=pts[pts[,1] >= (xlim[2]-range),,drop=FALSE]; Rshift=t(apply(R, 1, function(v) {v - c(xd, 0)}))
  L=pts[pts[,1] <= range,,drop=FALSE]; Lshift=t(apply(L, 1, function(v) {v + c(xd, 0)}))
  U=pts[pts[,2] >= ylim[2]-range,,drop=FALSE]; Ushift=t(apply(U, 1, function(v) {v - c(0, yd)}))
  D=pts[pts[,2] <= range,,drop=FALSE]; Dshift=t(apply(D, 1, function(v) {v + c(0, yd)}))

  LU=pts[(pts[,1] <= range) & (pts[,2] >= ylim[2]-range),,drop=FALSE]; LUshift=t(apply(LU, 1, function(v) {v + c(xd, -yd)}))
  RU=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] >= ylim[2]-range),,drop=FALSE]; RUshift=t(apply(RU, 1, function(v) {v + c(-xd, -yd)}))
  RD=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] <= range),,drop=FALSE]; RDshift=t(apply(RD, 1, function(v) {v + c(-xd, yd)}))
  LD=pts[(pts[,1] <= range) & (pts[,2] <= range),,drop=FALSE]; LDshift=t(apply(LD, 1, function(v) {v + c(xd, yd)}))
  if (length(Rshift)>0)
    pts=rbind(pts, Rshift)
  if (length(Lshift)>0)
    pts=rbind(pts, Lshift)
  if (length(Ushift)>0)
    pts=rbind(pts, Ushift)
  if (length(Dshift)>0)
    pts=rbind(pts, Dshift)
  if (length(LUshift)>0)
    pts=rbind(pts, LUshift)
  if (length(RUshift)>0)
    pts=rbind(pts, RUshift)
  if (length(RDshift)>0)
    pts=rbind(pts, RDshift)
  if (length(LDshift)>0)
    pts=rbind(pts, LDshift)
  pts
}

####################################################################################

get_max<-function(r,mat,av){
  num_mat=length(mat)
  L=mat[,3]
  Max=matrix()
  if (num_mat>3){
    first=0
    num_l=dim(mat)[1] 
    coord=cbind(mat[,1],mat[,2])
    dist_mat=as.matrix(dist(coord))
    for (i in 1:num_l){
      index=c(1:num_l)
      D=L[dist_mat[i,]<r]
      index=index[dist_mat[i,]<r]
      if(length(D)-1>(av)){
        ########### see if max
        MM=max(D)
        if(MM==L[i]){
          #the "plateau" case
          if (length(D[D[]==MM])>1){
            index=index[D[]==MM]
            if(i==max(index)){
              #print("The plateau case")
              if (first==0){
                Max=mat[i,1:3]
                first=1
              }else{
                Max=rbind(Max,mat[i,1:3])
              }
            }
          }else{
            if (first==0){
              Max=mat[i,1:3]
              first=1
            }else{
              Max=rbind(Max,mat[i,1:3])
            }
            
          }

        }
        }
    }
  }
  return(Max) 
}





get_tp_values<-function(r,Max_infot, mat){
  
  num_max=length(Max_infot)
  
  if (num_max==1){
    Max_with_TP=matrix()
    
  }
  
  if (num_max==3){
    Max_with_TP=c(0,0,0,0,0)
    Max_with_TP[1]=Max_infot[1] 
    Max_with_TP[2]=Max_infot[2] 
    Max_with_TP[3]=Max_infot[3]
    Max_with_TP[4]=0 
    Max_with_TP[5]=Max_infot[3] 
    
  }
  
  if (num_max>3){
    
    num_l=dim(mat)[1]
    num_c=dim(mat)[2]
    num_max=dim(Max_infot)[1]
    Max_info=matrix(0,num_max,4)
    Max_info[,1:3]=Max_infot
    
    for (h in 1:num_max){
      go=0
      for (g in 1:num_l){
        if (go==0){
          if (mat[g,1]==Max_info[h,1] && mat[g,2]==Max_info[h,2]){
            Max_info[h,4]=g
            go=1
          }
        }
      }
    }
    
    Max_with_TP=matrix()
    Max_with_TP=matrix(0,num_max,5)
    coord=matrix()
    coord=cbind(mat[,1],mat[,2])
    dist_mat=as.matrix(dist(coord))
    
    
    for (t in 1:num_max){
      
      
      mat_to_flood=matrix()
      labels_in_mat=1:num_l
      mat_to_flood=cbind(mat,labels_in_mat)
      num_max_sup=t-1
      saddle_point=0
      saddle_point_previous=0
      saddle_point_new=0
      flooding_coeff=Max_info[t,3]/2 
      flooding_previous=0
      flooding=0
      ##
      
      if (t==1){ # The Everest case
        Max_with_TP[t,1]=Max_info[t,1] 
        Max_with_TP[t,2]=Max_info[t,2] 
        Max_with_TP[t,3]=Max_info[t,3] 
        Max_with_TP[t,4]=Max_info[t,4] 
        Max_with_TP[t,5]=Max_info[t,3] 
        
      }else{
        
        finish=FALSE
        first_thing_first=0
        once=0
        
        while (finish==FALSE){
          
          mat_to_flood[,3]=mat[,3]-flooding
          C=which(mat_to_flood[,3]>0)
          
          if (length(C)>0){
            G=graph.adjacency(dist_mat[C,C]<r) 
            lab=clusters(G,"weak")
            init=num_l+1
            fin=2*num_l
            labels=init:fin
            labels[C]=lab$membership
            
            CONNECT=FALSE
            
            for (i in 1:num_max_sup){ 
              if (labels[Max_info[t,4]]==labels[Max_info[i,4]]){
                first_thing_first=1
                CONNECT=TRUE
                break
              }else{
                CONNECT=FALSE
              }
            }
            
            if (CONNECT==FALSE && first_thing_first==0){ 
              saddle_point=0 
              first_thing_first=1
              finish=TRUE
              break
            }
            if (CONNECT==FALSE && first_thing_first!=0){ 
              flooding=flooding_previous-flooding_coeff 
            }
            if (CONNECT==TRUE){ 
              flooding=flooding_previous+flooding_coeff
            }
            
            
            flooding_previous=flooding
            flooding_coeff=flooding_coeff/2
            
            if (flooding_coeff<0.5){ 
              if (CONNECT==TRUE){
                while (CONNECT==TRUE){
                  flooding=flooding+1
                  mat_to_flood[,3]=mat[,3]-flooding
                  C=which(mat_to_flood[,3]>0)
                  if (length(C)>1){
                    G=graph.adjacency(dist_mat[C,C]<r)
                    lab=clusters(G,"weak")
                    init=num_l+1
                    fin=2*num_l
                    labels=init:fin
                    labels[C]=lab$membership
                    clor=0
                    for (i in 1:num_max_sup){
                      if (clor==0){
                        if (labels[Max_info[t,4]]==labels[Max_info[i,4]]){
                          CONNECT=TRUE
                          clor=1
                          break
                        }else{
                          CONNECT=FALSE
                          break
                        }}
                      
                    }}else{
                      CONNECT=FALSE
                      break
                    }
                }
                saddle_point=flooding
                finish=TRUE
                break
              }
              
              if (CONNECT== FALSE){
                saddle_point=flooding
                finish=TRUE
                break
              }
              
            } 
            
          }
          
        } 
        Max_with_TP[t,1]=Max_info[t,1] 
        Max_with_TP[t,2]=Max_info[t,2] 
        Max_with_TP[t,3]=Max_info[t,3] 
        Max_with_TP[t,4]=Max_info[t,4] 
        Max_with_TP[t,5]=Max_info[t,3]-saddle_point 
        
      }
    } 
  } 
  
  return(Max_with_TP) 
}



above_th<-function(Max_TP_values,th){
  
  num_max=length(Max_TP_values)
  
  if (num_max==1){
    Max_above_th=matrix()
  }
  if (num_max==5){
    if (Max_TP_values[5]>th){
      Max_above_th=Max_TP_values
    }else{
      Max_above_th=matrix()
    }
  }
  if (num_max>5){
    num_l=dim(Max_TP_values)[1]
    index=0
    Max_above_th=matrix(0,1,5)
    if (num_l>1){
      for (i in 1:num_l){
        if (Max_TP_values[i,5]>th){
          index=index+1
          if (index==1){
            Max_above_th=Max_TP_values[i,]}else{
              Max_above_th=rbind(Max_above_th,Max_TP_values[i,])
            }
        }
      }
      if (index==0){
        Max_above_th=matrix()
      }}
  }
  return(Max_above_th)
}


get_clusters<-function(r, Max_above_th,  mat_croped, mat){
  
  num_l_croped=dim(mat_croped)[1]
  num_c_croped=dim(mat_croped)[2]
  num_l=dim(mat)[1]
  num_c=dim(mat)[2]
  num_max=length(Max_above_th)
  
  if (num_max==1){
    mat_with_clusters=matrix(0,num_l,3)
    mat_with_clusters[,1:2]=mat[,1:2]
    labels_in_mat=1:num_l
    mat_with_clusters=cbind(mat_with_clusters,labels_in_mat)
    
  }
  if (num_max==5){
    
    mat_with_clusters=matrix()
    mat_with_clusters=matrix(0,num_l,3)
    mat_with_clusters[,1:2]=mat[,1:2]
    mat_with_clusters_croped=matrix()
    mat_with_clusters_croped=matrix(0,num_l_croped,3)
    mat_with_clusters_croped[,1:2]=mat_croped[,1:2]
    coord=cbind(mat_croped[,1],mat_croped[,2])
    dist_mat=as.matrix(dist(coord))
    
    
    for (g in 1:num_l_croped){
      if (mat_croped[g,1]==Max_above_th[1] && mat_croped[g,2]==Max_above_th[2]){
        max_index=g
      }
    }
    
    C=which(mat_croped[,3]>0)
    G=graph.adjacency(dist_mat[C,C]<r)
    lab=clusters(G,"weak")
    init=num_l_croped+1
    fin=2*num_l_croped
    labels=init:fin
    labels[C]=lab$membership
    for (rr in 1:num_l_croped){
      if (labels[max_index]==labels[rr]){
        mat_with_clusters_croped[rr,3]=1
      }
    }
    
    
    for (k in 1:num_l_croped){
      yep=0
      for (l in 1:num_l){
        if (yep==0){
          if (mat_croped[k,1]==mat[l,1] && mat_croped[k,2]==mat[l,2]){
            mat_with_clusters[l,3]=mat_with_clusters_croped[k,3]
            yep=1
          } 
        }
      }
    }
    
    index=1
    for (i in 1:num_l){
      if (mat_with_clusters[i,3]==0){
        index=index+1
        mat_with_clusters[i,3]=index
      }
    }
    
  }
  
  if (num_max>5){
    
    num_max=dim(Max_above_th)[1]
    Max_info=matrix(0,num_max,4)
    Max_info[,1:3]=Max_above_th[,1:3]
    
    
    for (h in 1:num_max){
      go=0
      for (g in 1:num_l_croped){
        if (go==0){
          if (mat_croped[g,1]==Max_info[h,1] && mat_croped[g,2]==Max_info[h,2]){
            Max_info[h,4]=g
            go=1
          }
        }
      }
    }
    
    mat_with_clusters=matrix()
    mat_with_clusters=matrix(0,num_l,3)
    mat_with_clusters[,1:2]=mat[,1:2]
    mat_with_clusters_croped=matrix()
    mat_with_clusters_croped=matrix(0,num_l_croped,3)
    mat_with_clusters_croped[,1:2]=mat_croped[,1:2]
    
    coord=cbind(mat_croped[,1],mat_croped[,2])
    dist_mat=matrix()
    dist_mat=as.matrix(dist(coord))
    
    
    for (t in 1:num_max){
      t=num_max-t+1
      mat_to_flood=matrix()
      mat_to_flood=mat_croped[,1:3]
      flooding_previous=0
      flooding=Max_above_th[t,3]-Max_above_th[t,5]
      
      mat_to_flood[,3]=mat_croped[,3]-flooding
      C=which(mat_to_flood[,3]>0)
      
      if (length(C)>0){
        G=graph.adjacency(dist_mat[C,C]<r)
        lab=clusters(G,"weak")
        init=num_l_croped+1
        fin=2*num_l_croped
        labels=init:fin
        labels[C]=lab$membership
        CONNECT=FALSE
        for (rr in 1:num_l_croped){
          if (labels[Max_info[t,4]]==labels[rr]){
            if (mat_with_clusters_croped[rr,3]>0){
              CONNECT=TRUE
              cc=mat_with_clusters_croped[rr,3]
            }
          }}
        if (CONNECT==TRUE){
          
          n=which(Max_info[,4]==cc)
          
          mat_to_flood=matrix()
          mat_to_flood=mat_croped[,1:3]
          flooding_previous=0
          flooding=Max_above_th[n,3]-Max_above_th[n,5]# the TP value
          
          mat_to_flood[,3]=mat_croped[,3]-flooding
          C=which(mat_to_flood[,3]>0)
          
          if (length(C)>0){
            G=graph.adjacency(dist_mat[C,C]<r)
            lab=clusters(G,"weak")
            init=num_l_croped+1
            fin=2*num_l_croped
            labels=init:fin
            labels[C]=lab$membership}
        }
        
        for (rr in 1:num_l_croped){
          if (labels[Max_info[t,4]]==labels[rr]){
            mat_with_clusters_croped[rr,3]=Max_info[t,4]
          }}}}
    
    
    
    for (w in 1:num_max){
      count=0
      for (k in 1:num_l_croped){
        if  (mat_with_clusters_croped[k,3]==Max_info[w,4]){
          count=count+1
        }
      }
      if (count<5){
        for (k in 1:num_l_croped){
          if  (mat_with_clusters_croped[k,3]==Max_info[w,4]){
            mat_with_clusters_croped[k,3]=0
          }
        }
      }
    }
    
    
    for (k in 1:num_l_croped){
      yep=0
      for (l in 1:num_l){
        if (yep==0){
          if (mat_croped[k,1]==mat[l,1] && mat_croped[k,2]==mat[l,2]){
            mat_with_clusters[l,3]=mat_with_clusters_croped[k,3]
            yep=1
          } }}}
    
    
    index=max(Max_info[,4])
    for (i in 1:num_l){
      if (mat_with_clusters[i,3]==0){
        index=index+1
        mat_with_clusters[i,3]=index
      }
    }
    
  }
  return(mat_with_clusters) 
}

####################################################################################


#Kclust <- function(pts, sdsx=0, sdsy=0, sdsz=0, xlim, ylim, zlim, paramc, paramd, psd=NULL, minsd=NULL, maxsd=NULL, useplabel=TRUE, alpha=NULL, pb=.5, rseq=seq(10, 150, by=5), thseq=seq(0, 100, by=2.5), score=FALSE, rlabel=FALSE, report=TRUE){
  

Kclust <- function(pts, sds=0, xlim, ylim, psd=NULL, minsd=NULL, maxsd=NULL, useplabel=TRUE, alpha=NULL, pb=.5, rseq=seq(10, 200, by=5), thseq=seq(5, 500, by=5), score=FALSE, rlabel=TRUE, report=TRUE, clustermethod="K"){
    ##browser()
###################################################################################
  
  getting_in=0
  getting_in<<-getting_in
  N = dim(pts)[1]    
  if (N==1){
    rs=c()
    ths=c()
    for (r in rseq){
      for (th in thseq){
        rs=c(rs, r); ths=c(ths,th)
      }
    }
    labels=rep(1,length(rs)); dim(labels)=c(length(rs),1)
    return(list(scores=rep(0,length(rs)), scale=rs, thresh=ths, labels=labels))
  }
  tor=toroid(pts, xlim, ylim, max(rseq))
  D=as.matrix(dist(tor))
  D=D[1:N, 1:N]
  scores=c()
  retlabels=c()
  rs=c()
  ths=c()
  
  coord=cbind(pts[,1],pts[,2])
  dist_mat=matrix()
  num_l=dim(pts)[1]
  dist_mat=as.matrix(dist(coord))
  
  for (r in rseq){
    K=apply(D, 1, function(v){sum(v <= r)-1})
    L=sqrt((diff(xlim)+2*max(rseq))*(diff(ylim)+2*max(rseq))* K /(pi * (dim(tor)[1]-1)))
    
    NEW_R=TRUE
    
    num_l=dim(pts)[1]
    X=c()
    Y=c()
    
    X=runif(num_l, xlim[1], xlim[2])
    Y=runif(num_l, ylim[1], ylim[2])
    pts_sigma=matrix()
    pts_sigma=cbind(X,Y)
    
    
    coord=cbind(pts_sigma[,1],pts_sigma[,2])
    dist_mat=matrix()
    num_l=dim(pts_sigma)[1]
    dist_mat=as.matrix(dist(coord))
    Min=0
    
    for (i in 1:num_l){
      dist_mat_i=c()
      dist_mat_i=dist_mat[i, dist_mat[i,]!=0]
      Min=Min+min(dist_mat_i)
    }
    r_scanning=Min/num_l
    
    
    av_point_per_nm=num_l/(2000*2000)
    av_point_in_sphere=av_point_per_nm*(pi*r^2)
    
    
    
    tor_sigma=toroid(pts_sigma, xlim, ylim, max(rseq))
    D_sigma=as.matrix(dist(tor_sigma))
    D_sigma=D_sigma[1:N, 1:N]
    
    K_sigma=apply(D_sigma, 1, function(v){sum(v <= r)-1})
    L_sigma=sqrt((diff(xlim)+2*max(rseq))*(diff(ylim)+2*max(rseq))* K_sigma /(pi * (dim(tor_sigma)[1]-1)))
    L_mean=mean(L_sigma)
    print(L_mean)
    sigma_TBU=sd(L_sigma)
    print(sigma_TBU)
    
    
    if (NEW_R==TRUE){ 
      
      Lr=t(L)
      num_l=dim(pts)[1]
      mat=matrix(0,num_l,3)
      mat[,1]=pts[,1]
      mat[,2]=pts[,2]
     
      
      mat[,3]=Lr-(L_mean+sigma_TBU) # 95% confidence interval
      
      num_l=dim(mat)[1]
      print(dim(mat)[1])
      #print(mat[1:100,])
      #mat<<-mat
      index=c(1:num_l)
      index=index[mat[,3]>0]
      mat_croped=mat[index,]
      print(dim(mat_croped)[1])
      max_list=matrix()
      max_list=get_max(r,mat_croped,av_point_in_sphere) 
      print(dim(max_list)[1])
      
      NN=length(max_list)
      if (NN>3){
        
        num_max=dim(max_list)[1]
        
        max_list_ordered_coeff=order(max_list[,3],decreasing=TRUE)
        
        max_list_tempo=matrix(0,num_max,3)
        
        for (i in 1:num_max){
          for (j in 1:num_max){
            if (i==max_list_ordered_coeff[j]){
              max_list_tempo[j,]=max_list[i,]
            }
          }
        }
        max_list=c()
        max_list=max_list_tempo 
        num_max=dim(max_list)[1]
        max_list<<-max_list
        
      }
      if (NN==1){
        
        max_list=matrix()
      }
      
      max_tp_list=matrix()
     
      max_tp_list=get_tp_values(r_scanning,max_list, mat_croped)
      
      
      NEW_R=FALSE
    }
    
    
    NN=length(max_tp_list)
    if (NN>5){
      num_tp_max=dim(max_tp_list)[1]
      
      max_tp_list_ordered_coeff=order(max_tp_list[,5],decreasing=TRUE)
      
      max_tp_list_tempo=matrix(0,num_max,5)
      
      for (i in 1:num_tp_max){
        for (j in 1:num_tp_max){
          if (i==max_tp_list_ordered_coeff[j]){
            max_tp_list_tempo[j,]=max_tp_list[i,]
          }
        }
      }
      max_tp_list=c()
      max_tp_list=max_tp_list_tempo 
      num_tp_max=dim(max_tp_list)[1]
      max_tp_list<<-max_tp_list
      
    }
    if (NN==1){
      num_tp_max=0
      max_list=matrix()
    }
    if (NN==5){
      num_tp_max=1
      
      
    }
    
   
    if (num_tp_max>=1){    
      th_prec=-1
      for (p in 1:num_tp_max){
        q=(num_tp_max-p+1)

        if (num_tp_max==1){
          th=max_tp_list[5] 
        }
        else{
          th=max_tp_list[q,5] 
        }
        
        
        if (th!=th_prec && th>0){
          th_prec=th
          max_above_tpth_list=matrix()
          max_above_tpth_list=above_th(max_tp_list,th)
          mat_with_clusters_labels=matrix()
          
          mat_with_clusters_labels=get_clusters(r_scanning, max_above_tpth_list, mat_croped, mat)
          labels=matrix(0,num_l,1)
          labels=mat_with_clusters_labels[,3]
          s=0
          getting_in=getting_in+1
          getting_in<<-getting_in
          
          
          if (score){
            
            if (getting_in>1){
              
              getting_in=getting_in+1
              getting_in<<-getting_in
              
              
              s=scorewprec(labels=labels, pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
              if (s>scores){
                
                scores=s
                scores<<-scores
                ths_best=floor(th*10)
                rs_best=r
                ths_best<<-ths_best
                rs_best<<-rs_best
                
                retlabels_best=labels
                retlabels_best<<-retlabels_best
              }                 
              
            }else{
              
              s=scorewprec(labels=labels, pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
              print(s)
              getting_in=getting_in+1
              getting_in<<-getting_in
              scores=s
              scores<<-scores
              ths_best=floor(th)
              rs_best=r
              ths_best<<-ths_best
              rs_best<<-rs_best
              
              retlabels_best=labels
              retlabels_best<<-retlabels_best
            }
            
          }
          
          if (report){cat("Scale:", r, "Thr:", th, "Score: ", s, "\n")}

        }
      }
    }
    if (num_tp_max==0){
      
      labels=matrix(0,num_l,1)
      for (u in 1:num_l){
        labels[u]=u
      }
      s=0
      getting_in=getting_in+1
      getting_in<<-getting_in
      th=0
      if (score){
        if (getting_in>1){
          getting_in=getting_in+1
          getting_in<<-getting_in
          s=scorewprec(labels=labels, pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
          if (s>scores){
            scores=s
            scores<<-scores
            ths_best=floor(th)
            rs_best=r
            ths_best<<-ths_best
            rs_best<<-rs_best
          }                 
          
        }else{
          s=scorewprec(labels=labels, pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
          getting_in=getting_in+1
          getting_in<<-getting_in
          scores=s
          
          scores<<-scores
          ths_best=floor(th*10)
          rs_best=r
          ths_best<<-ths_best
          rs_best<<-rs_best
        }
        
      }
      
      if (report){cat("Scale:", r, "Thr:", th, "Score: ", s, "\n")}
      if (rlabel){
        if (getting_in>1){
          s=scorewprec(labels=labels, pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
          if (s>scores){
            retlabels_best=labels
            retlabels_best<<-retlabels_best
            
          }                 
          
        }else{
          retlabels_best=labels
          retlabels_best<<-retlabels_best
          
          
        }
        
        
      }
      
    }
  }
  
  
  print(retlabels_best)
  scores_best=scores
  scores_best<<-scores_best
  list(scores=scores_best, scale=rs_best, thresh=ths_best, labels=retlabels_best)
  print(list(scores=scores_best, scale=rs_best, thresh=ths_best, labels=retlabels_best))
###################################################################################  
  

}

writeRes <- function(res, rfile, labdir, bestonly=FALSE){
  scale=unique(res[["scale"]]); scale=scale[order(as.numeric(scale))]
  thresh = unique(res[["thresh"]]); thresh=thresh[order(as.numeric(thresh))]
  cat("0", scale, sep="\t", file=rfile); cat("\n", file=rfile, append=TRUE)
  for (line in thresh){
    scales=res[["scale"]][res[["thresh"]]==line]; o=order(scales); scales=scales[o]
    scores=res[["scores"]][res[["thresh"]]==line]; scores=scores[o]
    cat(line, "\t", sep="", file=rfile, append=TRUE); cat(scores, sep="\t", append=TRUE, file=rfile); cat("\n", file=rfile, append=TRUE)
  }
  dir.create(labdir)
  #print("tricky tricky")
  print(res[["labels"]])
  print(length(res[["labels"]]))
  
  if (bestonly) is=which.max(res[["scores"]])
  else 
    f=file.path(paste(labdir, "/clusterscale", res[["scale"]], "\ thresh", res[["thresh"]], "labels.txt", sep=""))
    cat(res[["labels"]], file=f, sep=","); cat("\n", file=f, append=TRUE)
 
}

nClusters <- function(labels){
  sum(table(labels)>1)
}

percentageInCluster <- function(labels){
  Nb=sum(table(labels)==1) ##That was lucky
  (length(labels)-Nb)/length(labels) * 100
}

molsPerCluster <- function(labels){
  ta=table(labels); ta[ta>1]
}

nMolsPerCluster <- function(labels){
  length(labels)*percentageInCluster(labels)/(100*nClusters(labels))
}

histnMols <- function(labels){
  ta=table(labels)[table(labels)>1]; h=hist(ta, plot=FALSE)
  plot(h, xlab="Number of molecules", ylab="Number of clusters", main="")
}

clusterRadii <- function(pts, labels){
  radii=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1) -1
    else {
        mean(c(sd(pts[v,1]), sd(pts[v,2])))
    }
    })
  radii[radii>=0]
}

convexHullAreas <- function(pts, labels){
  areas=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1) -1
    else {
      i<-chull(pts[v,1],pts[v,2])
      areapl(as.matrix(pts[v[i],]))
    }
    })
  areas[areas>=0]
}


reldensity <- function(pts, labels, xlim, ylim){
    rs=clusterRadii(pts, labels)
    tb=table(labels)
    nclustered = sum(tb[tb>=2])
    nb=length(labels)-nclustered
    areaclustered=sum(pi*rs^2)
    (nclustered/areaclustered)/(nb/(diff(xlim)*diff(ylim)-areaclustered))
}
##Old version with convex hull
## reldensity <- function(pts, labels, xlim, ylim){
##     tb=table(labels)
##     nclusteredgeq2 = sum(tb[tb>2])
##     nb=sum(tb[tb==1])
##     areaclustered=sum(convexHullAreas(pts, labels))
##     (nclusteredgeq2/areaclustered)/(nb/(diff(xlim)*diff(ylim)-areaclustered))
## }

plabel <- function(labels, alpha, pb){
  cnt <-tapply(1:length(labels), labels, length)
  cl =cnt[cnt!=1]
  B = length(labels)-sum(cl)
  Bcont = B*log(pb)+(1-B)*log(1-pb)
  ## Green 2001 p.357, Scand J Statist 28
  partcont=0
  if (length(cl) >0)
    partcont=length(cl)*log(alpha)+lgamma(alpha)+sum(lgamma(cl))-lgamma(alpha+sum(cl))
  Bcont+partcont
}

scorewprec <- function(labels, pts, sds, xlim, ylim, psd, minsd, maxsd, useplabel=TRUE, alpha=NULL, pb=.5){  
  s=sum(tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)>1) mcgaussprec(pts[v,], sds[v], xlim, ylim, psd=psd, minsd=minsd, maxsd=maxsd)
    else -log(diff(xlim)*diff(ylim))      
  }))  
  prlab=0
  if (useplabel){
    if (is.null(alpha)){
      cnt <-tapply(1:length(labels), labels, length)
      n =sum(cnt[cnt!=1])
      alpha=20
    }
    prlab=plabel(labels, alpha, pb)
  }
  s+prlab
}

#################################################################################################
## Additional functions for new descriptor list


nClustersi <- function(labels){
  
  s=length(labels)
  Nclus=c()
  M=0
  for (i in 1:s){
    if (i==1){
      M=as.numeric(labels[i])
    }else{
      if (M<as.numeric(labels[i])){
        M=as.numeric(labels[i])
      }
    }
  }
  
  
  for (j in 1:M){
    count=0
    for (i in 1:s){
      if (as.numeric(labels[i])==j) {
        count=count+1
      }
    }
    if (count>1){
      Nclus=c(Nclus, j)
    }
  }
  return(Nclus)

}

clusterRadiii <- function(pts, labels, i){
  id=i
  radii=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1) -1
    else {
      mean(c(sd(pts[v,1]), sd(pts[v,2])))
    }
  })

  radii[radii>=0]
  for (g in (1:length(radii))){
    l=as.numeric(row.names(radii)[g])
    if(l==id){
      RADII=radii[g]
    }    
  }

  return(RADII)
}


clusterLocation <- function(pts, labels, i){
  id=i
  locationx=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1){-1} 
    else {
      #print(pts[v,1])
      #print(mean(pts[v,1]))
      mean(pts[v,1])
    }
  })
  locationy=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1){-1} 
    else {
      #print(pts[v,1])
      #print(mean(pts[v,1]))
      mean(pts[v,2])
    }
  })
#print(dim(locationx))
#print(dim(locationy))
location=cbind(locationx,locationy)
#print(dim(location))
#print(location)
#print(dim(location)[1])
  
  for (g in (1:dim(location)[1])){
    l=as.numeric(row.names(location)[g])
   
    if(l==id){
      LOC=location[g,]
    }    
  }

  return(LOC)
}

molsPerClusteri <- function(labels,i){
  s=length(labels); 
  count=0
  n=as.numeric(i)
  for (h in (1:s)){
    a=as.numeric(labels[h])
    
    if (a==n){
      count=count+1
      
    }
  }
  
  return(count)
}

molsPerClusterii <- function(labels,i,D1,D2,D3){
  s=length(labels); 
  count_channel1=0
  count_channel2=0
  count_channel3=0
  count_channel4=0
  n=as.numeric(i)
  for (h in (1:s)){
    a=as.numeric(labels[h])
    
    if (h<=D1){
      if (a==n){
        count_channel1=count_channel1+1
        
      }
      
    }else{
      if (h<=D1+D2){
        if (a==n){
          count_channel2=count_channel2+1
          
        }
        
      }else{
        if (h<=D1+D2+D3){
          if (a==n){
            count_channel3=count_channel3+1
            
          }
          
        }else{
          #channel 4
          if (a==n){
            count_channel4=count_channel4+1
            
          }
        }
      }

      
    }

  }
  
  return(c(count_channel1,count_channel2,count_channel3,count_channel4))
}


#################################################################################################
