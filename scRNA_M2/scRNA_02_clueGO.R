library(data.table)
library(pheatmap)

enrich_IDcluster<-function(D,aa=0.35){
  colnames(D)<-c('identifier2cluster','Attribute')
  lst <- lapply(unique(D$identifier2cluster),function(x){unique(D$Attribute[D$identifier2cluster == x])})
  names(lst) <- unique(D$identifier2cluster)
  D <- stack(lst)[,2:1]
  D$ind <- as.character(D$ind)
  colnames(D)<-c('identifier2cluster','Attribute')
  D<-D[order(D$identifier2cluster),]
  Total_N<-nrow(D)
  D_unique<-D[!duplicated(D$identifier2cluster),]   #Extract unique gene name
  N<-nrow(D_unique)
  Kappa_matrix<-matrix(0,N+1, N+1) #create a matrix for kappa values
  for (i in 1:N){
    Kappa_matrix[i+1,1]=D_unique[i,1]
    Kappa_matrix[1,i+1]=D_unique[i,1]
  }
  
  Cat_unique<-D[!duplicated(D$Attribute),]
  M<-nrow(Cat_unique)
  paste("Total number of rows in your input file:", Total_N)
  paste ("Total number of unique identifiers:", N)
  paste ("Total number of unique attribute information:", M)
  
  ##########
  D <- as.data.table(D)
  dt <- data.table::dcast(D,identifier2cluster~Attribute,fun.aggregate=function(x){as.character(length(x))},value.var="Attribute")
  rownames(dt) <- dt$identifier2cluster
  dt <- as.matrix(dt[Kappa_matrix[2:(N+1),1],c("identifier2cluster",Cat_unique[,2])])
  dt1 <- t(data.frame(c("0",colnames(dt)[-1])))
  colnames(dt1) <- colnames(dt)
  dt <- rbind(dt1,dt)
  Gene_Cat_matrix <- dt
  
  dt <- as.data.table(t(combn(2:(N+1),2)))
  lst <- lapply(2:(N+1), function(x){
    D$Attribute[D$identifier2cluster==Kappa_matrix[x,1]]
  })
  dt <- dt[,.(kp={
    gene1=lst[[V1-1]]
    gene2=lst[[V2-1]]
    sum1 <- length(gene1)
    sum2 <- length(gene2)
    a <- length(intersect(gene1,gene2))
    b <- sum1-a; c <- sum2-a; d <- M+a-sum1-sum2
    kp <- ((a+d)*M-(a+b)*(a+c)-(c+d)*(b+d))/(M^2-(a+b)*(a+c)-(c+d)*(b+d))
  }),by=.(V1,V2)]
  dt <- rbind(dt,dt[,.(V1=V2,V2=V1,kp)],data.table(V1=2:(N+1),V2=2:(N+1),kp=1))
  dt <- dt[order(V1,V2)]
  mtx <- matrix(dt$kp,nrow = N,ncol = N)
  mtx <- cbind(Kappa_matrix[2:(N+1),1],mtx)
  mtx <- rbind(Kappa_matrix[1,],mtx)
  Kappa_matrix <- mtx
  #############
  write.table(Kappa_matrix,file="test2del.txt", sep="\t", col.names = F, row.names = F)
  
  #aa=0.35
  D <- as.data.frame(Kappa_matrix[-1,],stringsAsFactors=F)
  colnames(D)[-1] <- Kappa_matrix[1,-1]
  Row<-nrow(D)
  Col<-ncol(D)
  
  mtx <- apply(D[,-1], 1, as.numeric)
  colnames(mtx) <- colnames(D)[-1]
  mtx[mtx >= aa] <- 1
  mtx[mtx != 1] <- 0
  rownames(mtx) <- D[,1]
  idx <- hclust(dist(mtx),method="ward.D2")$order
  pheatmap(mtx[idx,idx],
           cluster_cols =F,
           cluster_rows = F,
           fontsize=8,
           fontsize_row=6,
           scale="none",
           show_colnames=T,
           show_rownames = T,
           cellheight = 8,cellwidth = 8,
           treeheight_row = 0,treeheight_col = 0,
           filename = "test2del.pdf")
  
  D <- cbind(D[,1],as.data.frame(mtx))
  colnames(D)[1] <- "X0"
  
  Group<-cbind(D,0)                 
  names = c(colnames(D), "Seed")
  colnames(Group)<- names
  
  for (i in 1:Row){                  
    S<-(sum(D[i,2:Col])-1)          
    if (S<3){
      Group[i,(Col+1)] <- 0
    }
    else{
      idx <- which(D[i,-1]==1)
      idx <- idx[idx != i]
      Temp <- D[idx,idx+1]
      a <- sum(as.matrix(Temp))
      
      if(2*a>=S*(S+1)){
        Group[i,(Col+1)]<-1
      } else{Group[i,(Col+1)]<-0}
    }
    #print (c(i,j))
  }
  #print("Complete selected seed genes")
  
  D <- Group
  Row <- nrow(D)
  Col <- ncol(D)
  Group <- D
  
  for(i in 1:(Row-1)){
    if(Group[i,Col]==1){
      for(j in (i+1):Row){
        if(Group[j,Col]==1 & Group[i,Col]==1){
          s1<-sum(Group[i,2:(Col-1)])
          s2<-sum(Group[j,2:(Col-1)])
          a <- sum(unlist(Group[i,2:(Col-1)])==1 & unlist(Group[j,2:(Col-1)])==1)
          if (2*a>s1){
            l <- which(unlist(Group[i,2:(Col-1)])==1 & unlist(Group[j,2:(Col-1)])==0)
            if (length(l)>0) {
              Group[j,l+1] <- 1
            }
            Group[i, Col]<-0          
          }
          else if(2*a>s2){
            m <- which(unlist(Group[i,2:(Col-1)])==0 & unlist(Group[j,2:(Col-1)])==1)
            if(length(m)>0) {
              Group[i,m+1] <- 1
            }
            Group[j, Col]<-0
          }
        }
      }
    }
  }
  
  names(Group)[Col] <- "Cluster"
  names(Group)[1] <- "Seed_gene"
  Group<-subset(Group,Cluster==1)
  rownames(Group)=Group[,1]
  Group<-as.data.frame(Group[,c(-1,-ncol(Group))])
  outlier=colnames(Group)[apply(Group,2,sum)==0]
  #print(outlier)
  result.list=list()
  result.list[["cluster"]]=Group
  result.list[["outlier"]]=outlier
  #Row<-nrow(Group)
  #Col<-ncol(Group)
  
  # for (i in 1:Row){
  #   
  #   Group[i,Col]<-i
  #   for (j in 2:(Col-1)){
  #     if (Group[i,j]==1){
  #       Group[i,j]<-colnames(Group)[j]
  #     }
  #   }
  # }
  return(result.list)
}




