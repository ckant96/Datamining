library(igraph)
library(Matrix)
library(R.matlab)
####
#population initialisation
popin<-function(grap,n){ 
  m<-matrix(nrow=n,ncol=ncol(grap))
  for(i in 1:n){
    m[i,]=ranpop(grap)
  }
  m
}
ranpop<-function(grap1)
{
  x<-c(1:ncol(grap1))
  for(i in 1:ncol(grap1)){
    y=sample(1:ncol(grap1),ncol(grap1),replace=F)
    for(j in y){
      if(grap1[i,j]==1){
        x[i]=j
        break
      }
    }
  }
  x
}
#######
comunity<-function(m){ 
  c<-matrix(nrow=nrow(m),ncol=ncol(m))
  for(i in 1:nrow(m)){
    k<-com(m[i,],ncol(m))
    
    for(j in 1:ncol(m)){
      c[i,j]<-k[j]
    }
  }
  c 
}
######
com<-function(a,n){
  x<-c(1:n)
  y<-c(1:n)
  ## y=a
  for(p in 1:n){
    y[p]=a[p]
  }
  d<-c(1:n)
  for(z in 1:n){
    d[z]<-0
  }
  l=0
  for(i in 1:n){
    if(i==1){
      d[i]=1
      l=1
    }
    else{
      if(i==2){
        if(x[1]==y[2]||y[1]==x[2]){
          d[2]=d[1]
        }
      }
      else{
        for(j in 1:(i-1) ){
          if(x[j]==y[i]||y[j]==x[i]){
            d[i]=d[j]
            break
          }
        }}
      if(d[i]==0){
        l=l+1
        d[i]=l
      }
    }
  }
  d
}
################
comm<-function(a,n){
  x<-c(1:n)
  y<-c(1:n)
  for(p in 1:n){
    y[p]<-a[p]
  }
  d<-c(1:n)
  for(z in 1:n){
    d[z]<-0
  }
  print(d[3])
  l=1
  for(i in 1:n){
    k<-0
    b<-0
    if(d[i]==0){
      print(d[i])
      while(b!=x[i]&&b!=y[i]){
        print(d[i])
        if(k==0){
          b<-y[i]
          d[i]<-l
          d[b]<-l
          b<-y[b]
          k<-k+1
        }
        else if(d[b]==0){
          d[b]<-l
          b<-y[b]}
      }
      l<-l+1
    }
  }
  d
}

############
com2<-function(a,n){
  x<-c(1:n)
  y<-c(1:n)
  for(p in 1:n){
    y[p]<-a[p]
  }
  d<-c(1:n)
  for(z in 1:n){
    d[z]<-0
  }
  b=0
  l=1
  for(i in 1:n){
    if(i==1){
      d[i]=l
      b=y[i]
      d[b]=l
      b=y[b]
    }
    else if(i!=1&&d[i]==0){
      if(i!=n){
        for(j in (i+1):n){
          if(x[i]==y[j]&&d[y[j]]!=0){
            d[i]=d[y[j]]
            break
          }
          else if(y[i]==x[j]&&d[x[j]]!=0){
            d[i]=d[x[j]]
            break
          }
        }
      }
    }
    if(i==1||d[i]==0){
      while(b!=y[i]||b!=x[i]){
        d[b]=d[i]
        b=y[b]
      }
    }
    if(d[i]==0){
      l=l+1
      d[i]=l
    }
  }
  d 
  
}
#############################################3
community_no<-function(a,n){
  visited<-c(1:n)
  com<-c(1:n)
  for(i in 1:n){
    visited[i]<-0
  }
  visited[1]<-1
  v<-which(q %in% 1)
  k<-1
  com[1]<-1
  for(i in v){
    com[i]<-1
  }
  m<-matrix(nrow=n,ncol=n)
  for(i in 1:n){
    if(visited[i]==0){
      v<-which(q %in% i)
    }
  }
}
#########################################################3
com_no<-function(chrom,col){
  m<-matrix(nrow=col,ncol=col)
  for( i in 1:col){
    for(j in 1:col){
      m[i,j]<-0
    }
  }
  for( i in 1:col){
    m[i,chrom[i]]<-1
    m[chrom[i],i]<-1
  }
  #print(m)
  com<-c(1:col)
  visited<-c(1:col)
  for(i in 1:col){
    com[i]<-0
    visited[i]<-0
  }
  k<-1
  for(i in 1:col){
    if(visited[i]==0){
      visited[i]<-k
      com[i]<-k
      for(j in 1:col){
        if(m[i,j]==1){
          if(com[i]==0){
            com[i]<-k
          }
          if(com[j]==0){
            com[j]<-k
          }
          for(l in 1:col){
            if(m[l,j]==1){
              if(com[l]==0){
                com[l]<-k
              }
              if(com[j]==0){
                com[j]<-k
              }
            }
          }
          visited[j]<-1
        }
      }
      k<-k+1
    }
  }
  com
}
####3####  function to get S matrices (takes input the grap and community array) ##########

##### functtion to append new elements in list #########3
###########################################3
get_smat<-function(grap,com){
  m<-max(com)
  counts<-table(com)
  s_mat<-list()
  for(i in 1:m){
    c<-counts[names(counts)==i]
    mat<-matrix(nrow=c,ncol=c)
    node<-list()
    k=1
    for(j in 1:length(com)){
      if(com[j]==i){
        node[k]=j
        k=k+1
      }
    }
    for(l in 1:c){
      for(m in 1:c){
        if(grap[ node[[l]],node[[m]] ]==1){
          mat[l,m]<-1
        }
        else{
          mat[l,m]<-0
        }
      }
    }
    s_mat[[i]]<-mat
  }
  s_mat
}
###################################################
########### comunity score ####################
comunity_score<-function(grap,com,r){
  s_mat<-get_smat(grap,com)
  qvalue<-c(1:length(s_mat))
  for( i in 1:length(s_mat)){
    s_i<- s_mat[[i]]
    row_means<-c(1:nrow(s_i))
    I<-nrow(s_i )
    for( j in 1:I){
      row_means[j]<-sum(s_i[j,])
      row_means<-(row_means)/I
    }
    vs<-sum(s_i)
    #print("VS")
    #print(vs)
    for( j in 1:I){
      row_means[j]<-row_means[j]^r
    }
    ms<-sum(row_means)
    #print("MS")
    #print(ms)
    ms<-row_means/I
    qvalue[i]<-ms*vs
  }
  sum(qvalue)
}
###############################
cs<-function(g,com,r){
  r<-nrow(com)
  cscore<-c(1:r)
  for(i in 1:r ){
    cscore[i]<-comunity_score(g,com[i,],r)
  }
  cscore
}

###################################################
#####modularity
mvalue<-function(g){
  l=0
  ##for(i in 1:nrow(g)){
  ###  for(j in 1:ncol(g)){
  ##   l=l+g[i,j]
  ## }
  ##}
  l=sum(g)/2
  l
}
kvalue<-function(g){
  l<-c(1:nrow(g))
  for(i in 1:nrow(g)){
    l[i]<-0
  }
  for(i in 1:n){
    for(j in 1:n ){
      l[i]=l[i]+g[i,j]
    }
  }
  l
}
delta<-function(a,i,j){
  l=0
  if(a[i]==a[j])
    l=1
  l
}
qvalue<-function(g,m,k,c){
  l=0
  for(i in 1:nrow(g)){
    for(j in 1:nrow(g)){
      l=l+((g[i,j]-((k[i]*k[j])/(2*m)))*delta(c,i,j))
    }
  }
  l/(2*m)  
}
modularity<-function(g,c){
  n<-nrow(g)
  q<-c(1:nrow(c))
  m<-mvalue(g)
  k<-kvalue(g)
  for(i in 1:nrow(c)){
    q[i]<-qvalue(g,m,k,c)
  }
  q
}
#############################################################
#######    selection #############
totalfit<-function(q,n){
  l=0
  for(i in 1:n){
    l=l+q[i]
  }
  l
}
cumprob<-function(q,t,n){
  z<-c(1:n)
  p<-c(1:n)
  for(i in 1:n){
    p[i]<-(q[i]/t)
  }
  l=0
  for(i in 1:n){
    z[i]=l+p[i]
    l=z[i]
  }
  z
}
selection<-function(c,q){
  n<-nrow(c)
  mm<-max(q)
  for(i in 1:length(q)){
    q[i]<-mm-q[i]+0.1
  }
  tfit<-totalfit(q,n)
  s<-cumprob(q,tfit,n)
  mm<-max(q)
  rrr<-runif(1,0,1)
  l=n
  for(i in 1:n){
    if(s[i] >= rrr){
      l=i
      break
    }
  }
  l
}
################################################
########## mutation ################
mutation<-function(g,m){
  n<-ncol(g)
  r<-sample(1:n,1)
  k<-sample(1:n,n,replace = F)
  for(i in 1:n){
    if(m[i]!=r&&g[r,k[i]]==1){
      m[i]=k[i]
      break
    }
  }
  m 
}
###############################################
#################### multipoint mutation ####################
multipmutation<-function(chromo,nocol,grap){
  newchromo<-chromo
  randchoice=sample(1:nocol,1,replace=F)
  #print("no of points to be mutetd")
  #print(randchoice)
  randindex=sample(1:nocol,randchoice,replace=F)
  #print("points that are muteted")
  #print(randindex)
  for(i in randindex){
    randnodes=sample(1:nocol,nocol,replace=F)
    breakcount=0
    for(j in randnodes ){
      if(grap[i,j]==1&&newchromo[i]!=j){
        newchromo[i]=j
        breakcount=1
        break
      }
    }
  }
  newchromo
}
########### crossover #####################
crossover<-function(a,b,n){
  f<-c(1:n)
  d<-sample(0:1,n,replace=T)
  s=0
  for(i in d){
    if(i==1){
      s=s+1
      f[s]=a[s]
    }
    else{
      s=s+1
      f[s]=b[s]
    }
  }
  f
}

#############################################################33
############## elitism #########################
elitism<-function(q){
  k=max(q)
  l=which(q %in% k)
  l[1]
}
##########################################################################
fitness<-function(q,cs){
  fit<-c(1:length(q))
  for(i in 1:length(q)){
    fit[i]<-((1-q[i])+(10/(1+cs[i])))
  }
  fit
}
############## ga function ####################

ga<-function(g,cr,mr,er,r,pop,itr){
  m<-popin(g,pop)
  c<-comunity(m)
  com_sc<-cs(g,c,r)
  q<-modularity(g,c)
  fit<-fitness(q,com_sc)
  newm=matrix(nrow=pop,ncol=ncol(m))
  q<-fit
  for(i in 1:itr){
    s=0
    while(s<pop){
      re<-runif(1,0,1)
      rm<-runif(1,0,1)
      rc<-runif(1,0,1) 
      ##elitism####
      if(re<er&&s<pop){
        choice=runif(1,0,1)
        if(choice<0.5){
          h<-elitism(q)
          s<-(s+1)
          newm[s,]=m[h,]
        }
        else{
          m1<-max(q)
          f<-which(q %in% m1)
          m2<-min(q)
          l<-which(q %in% m2)
          m[l[1],]<-m[f[1],]
        }
      }
      ##### mutation ######
      if(rm<=mr&&s<pop){
        ##### multip for choosing single or mutipoint mutation############
        multip<-sample(0:1,1)
        
        if(multip==0){
          t<-selection(c,q)
          #print("gene 1")
          #print(m[t,])
          mutate<-mutation(g, m[t,])
          #print("offspring")
          #print(mutate)
          s<-(s+1)
          for(i in 1:ncol(m)){
            newm[s,i]<-mutate[i]
          }
        }
        else{
          t<-selection(c,q)
          mutate<-multipmutation( m[t,],ncol(g),g)
          #print("gene 1")
          #print(m[t,])
          #print("offspring")
          #print(mutate)
          s<-(s+1)
          for(i in 1:ncol(m)){
            newm[s,i]<-mutate[i]
          }
        }
      }
      
      #######crossover#########
      if(rc<=cr&&s<pop){
        a<-0
        b<-0
        while(a==b){
          a<-selection(c,q)
          b<-selection(c,q)
        }
        crosschild<-crossover(m[a,],m[b,],ncol(m))
        s<-(s+1)
        for(i in 1:ncol(m)){
          newm[s,i]<-crosschild[i]
        }
      }
      
    }
    c1<-comunity(newm)
    c<-c1
    m<-newm
    css<-cs(g,c,r)
    qval<-modularity(g,c)
    q<-fitness(qval,css)
  }
  o=max(q)
  mod<-modularity(g,c)
  print(max(c))
  print(c)
  print(mod)
  o
}

