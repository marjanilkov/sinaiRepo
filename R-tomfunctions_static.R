#****************************************************************************************
#
#                       Gene co-Expression Network Analysis Package
#         
#
#        Objective: 
#                   Automatic construction of gene co-expression networks &
#                   Decomposition of a network into tightly co-expressed modules
#
#        Input:     DNA microarray data file
#
#        Output:    Topological overlap matrix (TOM) and gene co-expressed modules
#
#        Authors:   Bin Zhang
#
#
#        Contact:   binzhang.ucla@gmail.com
#
#
#        Date:      Aug. 18, 2004
#
#        
#        REV            DATE            BY           DESCRIPTION
#
#        
#****************************************************************************************



#****************************************************************************************
#
#      APCLUSTER Affinity Propagation Clustering (Frey/Dueck, Science 2007)
#
#         -- transformed by Bin from the original Matlab code, 08/29,2007
#
# [idx,netsim,dpsim,expref]=APCLUSTER(s,p) clusters data, using a set 
# of real-valued pairwise data point similarities as input. Clusters 
# are each represented by a cluster center data point (the "exemplar"). 
# The method is iterative and searches for clusters so as to maximize 
# an objective function, called net similarity.
# 
# For N data points, there are potentially N^2-N pairwise similarities; 
# this can be input as an N-by-N matrix 's', where s(i,k) is the 
# similarity of point i to point k (s(i,k) needn? equal s(k,i)).  In 
# fact, only a smaller number of relevant similarities are needed; if 
# only M similarity values are known (M < N^2-N) they can be input as 
# an M-by-3 matrix with each row being an (i,j,s(i,j)) triple.
# 
# APCLUSTER automatically determines the number of clusters based on 
# the input preference 'p', a real-valued N-vector. p(i) indicates the 
# preference that data point i be chosen as an exemplar. Often a good 
# choice is to set all preferences to median(s); the number of clusters 
# identified can be adjusted by changing this value accordingly. If 'p' 
# is a scalar, APCLUSTER assumes all preferences are that shared value.
# 
# The clustering solution is returned in idx. idx(j) is the index of 
# the exemplar for data point j; idx(j)==j indicates data point j 
# is itself an exemplar. The sum of the similarities of the data points to 
# their exemplars is returned as dpsim, the sum of the preferences of 
# the identified exemplars is returned in expref and the net similarity 
# objective function returned is their sum, i.e. netsim=dpsim+expref.
# 
# 	[ ... ]=apcluster(s,p,'NAME',VALUE,...) allows you to specify 
# 	  optional parameter name/value pairs as follows:
# 
#   'maxits'     maximum number of iterations (default: 1000)
#   'convits'    if the estimated exemplars stay fixed for convits 
#          iterations, APCLUSTER terminates early (default: 100)
#   'dampfact'   update equation damping level in [0.5, 1).  Higher 
#        values correspond to heavy damping, which may be needed 
#        if oscillations occur. (default: 0.9)
#   'plot'       (no value needed) Plots netsim after each iteration
#   'details'    (no value needed) Outputs iteration-by-iteration 
#      details (greater memory requirements)
#   'nonoise'    (no value needed) APCLUSTER adds a small amount of 
#      noise to 's' to prevent degenerate cases; this disables that.
# 
#****************************************************************************************
#
#   function [idx,netsim,dpsim,expref]=apcluster(s,p,varargin);
#
#s=1-dist1; p=rep(1,no.genes); maxits=1000; convits=100;
#       lam=0.5; plt=T; details=F; nonoise=F;

# turn list into matrix
#    mres     = lapply(termGenes, fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=totalbackground)
#    datMtrxB = data.frame(do.call(rbind, mres))
#    datMtrxB = t(datMtrxB)

#****************************************************************************
# notice that the different convention in R and Matlab
# row-based operation is 1 in R but is 2 in matlab
# vice versa for the column-based operation
#
# [tmp c]=max(S(:,I),[],2); find max for each row
#
#  Matlab  <===> R
#
# tmp = apply(s, 1, max)
# c   = apply(s, 1, which.max)
#
#****************************************************************************

repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

# collapse a list of matrices into a single matrix 
#final = data.frame(do.call(rbind, cnvProbes))

# first element is index
findIJmin=function(mvect){
   val = mvect[-1];       
   minv=ifelse(val>mvect[1], mvect[1], val)
   return (minv)
}

# modMtrx = lapply(exprAll, pfnParsing)
# modMtrx = data.frame(do.call(rbind, modMtrx))
#
apcluster = function(s,p,maxits=1000, convits=100, nochange_cut=15,
       lam=0.9, plt=T, details=F, nonoise=F, symmetric=T, debug=F, 
       overlappingModule=F, minModuleSize=10, useNumberAsLabel=F, startlabel=0)
{

    #if nargin==0, % display demo
    #print('Affinity Propagation (APCLUSTER) sample code\n\n');
    #print('N=100; x=rand(N,2); % Create N, 2-D data points\n');
    #print('M=N*N-N; s=zeros(M,3); % Make ALL N^2-N similarities\n');
    #print('j=1;\n');
    #print('for i=1:N\n');
    #print('  for k=[1:i-1,i+1:N]\n');
    #print('    s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);\n');
    #print('    j=j+1;\n');
    #print('  end;\n');
    #print('end;\n');
    #print('p=median(s(:,3)); #Set preference to median similarity\n');
    #print('[idx,netsim,dpsim,expref]=apcluster(s,p,plot);\n');
    #print('print(''Number of clusters: %%d\\n,length(unique(idx)));\n');
    #print('print(''Fitness (net similarity): %%f\\n,netsim);\n');
    #print('figure; % Make a figures showing the data and the clusters\n');
    #print('for i=unique(idx)''\n');
    #print('  ii=find(idx==i); h=plot(x(ii,1),x(ii,2),''o''); hold on;\n');
    #print('  col=rand(1,3); set(h,''Color'',col,''MarkerFaceColor'',col);\n');
    #print('  xi1=x(i,1)*ones(size(ii)); xi2=x(i,2)*ones(size(ii)); \n');
    #print('  line([x(ii,1),xi1]'',[x(ii,2),xi2]'',''Color'',col);\n');
    #print('end;\n');
    #print('axis equal tight;\n\n');
    #   return;
    #end;
    start = proc.time()[3];

    i=1;

    if (lam>0.9) {
       print('\n*** Warning: Large damping factor in use. Turn on plotting to monitor the net similarity. The algorithm will change decisions slowly, so consider using a larger value of convits');
    }

    # Check that standard arguments are consistent in size

    # Construct similarity matrix
    S=s
    N=dim(s)[1]
    idxN = c(1:N)

    # In case user did not remove degeneracies from the input similarities,
    # avoid degenerate solutions by adding a small amount of noise to the
    # input similarities
    #if (!nonoise) {
    #    rns=randn('state'); randn('state',0);
    #    S=S+(eps*S+realmin(class(s))*100).*rand(N,N);
    #    randn('state',rns);
    #}

    if (!nonoise) {
        rns=runif(N*N, min=0, max=1);
        rnm=matrix(rns, N, N);
        S=S+(2e-16*S)*rnm;
        
    }

    # Place preferences on the diagonal of S
    #if length(p)==1 for i=1:N S(i,i)=p; end;
    #else for i=1:N S(i,i)=p(i); end;
    #end;

    diag(S) <- p

    # Allocate space for messages, etc
    dS=diag(S); A=matrix(0, N,N); R=matrix(0, N,N); t=1;

    if (plt) {
        netsim=rep(0,maxits+1);
    }

    if (details) {
        idx   =matrix(0, N, maxits+1);
        netsim=rep(0,maxits+1); 
        dpsim =rep(0,maxits+1); 
        expref=rep(0,maxits+1); 
    }

    # Execute parallel affinity propagation updates
    e=matrix(0, N,convits); dn=F; i=0;

    #if symmetric, ST=S; else ST=S'; end; # saves memory if it's symmetric
    ST = S
    I0 = sample(c(1:N), sample(c(1:N), 1) ) # random exemplars
    KI0= length(I0)
    nochange_cnt = 0
    while (!dn ) {
        i=i+1; 

        # Compute responsibilities
        Rold=R;
        AS=A+S; 
        #[Y,I]=max(AS,[],2);, matlab code
        Y = apply(AS, 1, max)
        I = apply(AS, 1, which.max)        
        #for (i in c(1:N) ){ AS[i,I[i]]=-Inf; }
        IidxN = cbind(idxN, I)
        AS[IidxN] = -Inf;

        Y2= apply(AS, 1, max)
        I2= apply(AS, 1, which.max)

        #R=S-repmat(Y,[1,N]);
        tmp=repmat(Y,1,N);
        R=S-tmp;
        rm(tmp)

        #for (i in c(1:N) ) { R[i,I[i]]=S[i,I[i]]-Y2[i];}
        R[IidxN ]= S[IidxN ] - Y2[idxN]
 
        R=(1-lam)*R+lam*Rold; # Dampen responsibilities
        
        # Compute availabilities
        Aold=A;        
        Rp=ifelse(R>0, R, 0); #Rp=max(R,0);
        #for (k in c(1:N)) { Rp[k,k]=R[k,k];}
        diag(Rp) <- diag(R)

        #A=repmat(sum(Rp,1),[N,1])-Rp;
        psum = apply(Rp, 2, sum)
        tmp=repmat(rbind(psum),N,1);
        A=tmp-Rp;
        rm(tmp)

        dA=diag(A); 
        A=ifelse(A<0,A, 0); #A=min(A,0);
        #for (k in c(1:N)) { A[k,k]=dA[k];}
        diag(A) <- dA

        A=(1-lam)*A+lam*Aold; # Dampen availabilities

        # Check for convergence
        AB = diag(A)+diag(R);
        E  = AB>0;
        e[,(i-1)%%convits+1]=E; #e(:,mod(i-1,convits)+1)=E;
        K  = sum(E);
        #K

        if(K>0) {
            I=c(1:length(E))[E]
            if(length(I0)!=K) {
               nochange_cnt = 0; # exemplar changes
            } else{
               if (sum(I0!=I)>0 ){
                  nochange_cnt = 0; # exemplar changes
                }else {
                  nochange_cnt = nochange_cnt + 1    
                }
            }
            I0=I
        }

        if (i>=convits) {
            #se=sum(e,2); sum for each row
            se=apply(e, 1, sum)
            unconverged=(sum((se==convits)+(se==0)) !=N);
            if ( (!unconverged & (K>0)) |(i==maxits) ) {
               dn=T;
            }
        }

        # Handle plotting and storage of details, if requested
        if (plt | details ) {
            if (K==0) {
                tmpnetsim=NA; tmpdpsim=NA; tmpexpref=NA; tmpidx=NA;
            } else {
                #I=find(E); [tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c);

                I=c(1:length(E))[E];# [tmp c]=max(S[,I],[],2); 
                Si  = S[,I]
                tmp = apply(Si, 1, max)
                c   = apply(Si, 1, which.max)

                c[I]=c(1:K); tmpidx=I[c];
                #tmpnetsim=sum(S((tmpidx-1)*N+[1:N]'));
                tmpnetsim =sum( S[(tmpidx-1)*N+c(1:N) ] )

                tmpexpref=sum(dS[I]); tmpdpsim=tmpnetsim-tmpexpref;
            }
        }

        if (details) {
            netsim[i]=tmpnetsim; dpsim[i]=tmpdpsim; expref[i]=tmpexpref;
            idx[,i]=tmpidx;
        }

        #if plt & toc>1, tic;
        if (plt){
            netsim[i]=tmpnetsim;
            tmp  = c(1:i); 
            #tmpi= find(~isnan(netsim(1:i)));
            tsel = !is.na(netsim[1:i])
            if(sum(tsel)>0) {
               tmpi = tmp[ tsel ];
               plot(tmp[tmpi],netsim[tmpi], xlab="# Iterations",
                    ylab="Fitness (net similarity) of quantized intermediate solution");
            }
        }

        if(debug){
           print(round(A,3))
           print(round(R,3))
           print(round(A+R,3))
           print(round(A+R+s,3))
        }

        if(nochange_cnt>=nochange_cut){
           dn=T
           unconverged = F
           break
        }

    } #while

    #I=find(diag(A+R)>0); 
    ARSel = diag(A+R)>0;
    I     = c(1:length(ARSel) )[ARSel]
    K     = length(I); # Identify exemplars

    if (K>0) {
        # assign the module labels
        #[tmp c]=max(S(:,I),[],2); c(I)=1:K; idx=I(c); % Assignments
        Si = S[,I]
        if(K>1) {
          tmp = apply(Si, 1, max)
          c   = apply(Si, 1, which.max)
        }else{
          tmp = Si
          c   = rep(1, length(Si) )
        }

        c[I]=c(1:K); # for exemplars themselves
        tmpidx=I[c];

        tmpnetsim=sum(S[(tmpidx-1)*N+c(1:N) ]); 
        tmpexpref=sum(dS[I]);
        

    }else {
        tmpidx=rep(NA, N); tmpnetsim=NA; tmpexpref=NA;
    }

    if (details) {
        netsim[i+1]=tmpnetsim; netsim=netsim[1:(i+1)];
        dpsim[i+1]=tmpnetsim-tmpexpref; dpsim=dpsim[1:(i+1)];
        expref[i+1]=tmpexpref; expref=expref[1:(i+1)];
        idx[,i+1]=tmpidx; idx=idx[,1:(i+1)];
    }else{
        netsim=tmpnetsim; dpsim=tmpnetsim-tmpexpref;
        expref=tmpexpref; idx=tmpidx;
    }

    if (plt | details) {
        print(paste("Number of exemplars identified: ", K, " (for ", N, " data points)") );
        print(paste("Net similarity: ",tmpnetsim) );
        print(paste("  Similarities of data points to exemplars: ",dpsim[length(dpsim)], "\n") );
        print(paste("  Preferences of selected exemplars: ",tmpexpref, "\n"));
        print( paste("Number of iterations: ",i, "\n\n") );
            #print('Elapsed time: %f sec\n',etime(clock,start));
        totaltime  <- proc.time()[3]-start
        print( paste("Elapsed time: ", round(totaltime,3), " sec; ", i, " iterations") );
    }

    if (unconverged) {
        print("\n*** Warning: Algorithm did not converge. Activate plotting\n");
        print("    so that you can monitor the net similarity. Consider\n");
        print("    increasing maxits and convits, and, if oscillations occur\n");
        print("    also increasing dampfact.\n\n");
    }

    totaltime  <- proc.time()[3]-start
    print( paste(K, " exemplars, Elapsed time: ", round(totaltime,3), " sec; ", i, " iterations") );


    # overlapping modules
    #
    ARS = A[,I]+R[,I]+s[,I]
    for(i in c(1:K) ){ # handle exemplars, belong to themselves
       ARS[I(i),i] = 1
    }    
    ARS = ifelse(ARS >0, 1, 0)
    multiModule <- apply(ARS, 1, sum)
    
    # hand grey genes
    #
    idxG    = ifelse(multiModule ==0, 0, idx)
    greySel = multiModule ==0
    nogreys = sum(multiModule) # no of grey nodes


    # assign module colors
    #
    modulecolor = assignModuleColor(as.character(idxG), minsize1=minModuleSize-1,  
                                       anameallmodules=F, auseblackwhite=F, 
                                       useNumberAsLabel=useNumberAsLabel, 
                                       startlabel=startlabel)
    # locate the exemplar for each module
    mc      = as.character(modulecolor)
    mcOrder = order(idxG)
    mcD     = mc[mcOrder]
    idxD    = idxG[mcOrder]
    explIdx = findBoundary(mcD)    
    me      = cbind(mcD[explIdx], idxD[explIdx])
    colnames(me) <- c("module", "exemplar")
    meSel = as.integer(me[,2]) >0 # remove (grey) exemplar '0'
    me = me[meSel,]

    # non-overlapping modules
    #
    if(!overlappingModule) {
       retlist = list(modulecolor, idxG, as.integer(me[,2]), me) 
       names(retlist) <- c("mcolor", "mexemplar", "exemplar", "exemplar4module")
       return (retlist)
    }
    
    # ARS doesn't have the grey module, so we need add it as first column
    #
    if (nogreys>0) {
        greyCol = greySel + 0
        ARSwg   = cbind(greyCol, ARS)
    }else{
        ARSwg   = ARS
    }

    colnames(ARSwg) <- mcD[explIdx]
   
    retlist = list(modulecolor, idx, idxD[explIdx], me, ARSwg)    
    names(retlist) <- c("mcolor", "mexemplar", "exemplar", 
                        "exemplar4module", "overlapmodule")
    return (retlist)

    #return (tmpidx) #(idx)
}


########################### basic version ########################################
#
#
apclusterBasic = function(s,p,maxits=1000, convits=100,
       lam=0.9, plt=T, details=F, nonoise=F, symmetric=T, debug=F,
       overlappingModule=F, minModuleSize=10, useNumberAsLabel=F, startlabel=0)
{

    #if nargin==0, % display demo
    #print('Affinity Propagation (APCLUSTER) sample code\n\n');
    #print('N=100; x=rand(N,2); % Create N, 2-D data points\n');
    #print('M=N*N-N; s=zeros(M,3); % Make ALL N^2-N similarities\n');
    #print('j=1;\n');
    #print('for i=1:N\n');
    #print('  for k=[1:i-1,i+1:N]\n');
    #print('    s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);\n');
    #print('    j=j+1;\n');
    #print('  end;\n');
    #print('end;\n');
    #print('p=median(s(:,3)); #Set preference to median similarity\n');
    #print('[idx,netsim,dpsim,expref]=apcluster(s,p,plot);\n');
    #print('print(''Number of clusters: %%d\\n,length(unique(idx)));\n');
    #print('print(''Fitness (net similarity): %%f\\n,netsim);\n');
    #print('figure; % Make a figures showing the data and the clusters\n');
    #print('for i=unique(idx)''\n');
    #print('  ii=find(idx==i); h=plot(x(ii,1),x(ii,2),''o''); hold on;\n');
    #print('  col=rand(1,3); set(h,''Color'',col,''MarkerFaceColor'',col);\n');
    #print('  xi1=x(i,1)*ones(size(ii)); xi2=x(i,2)*ones(size(ii)); \n');
    #print('  line([x(ii,1),xi1]'',[x(ii,2),xi2]'',''Color'',col);\n');
    #print('end;\n');
    #print('axis equal tight;\n\n');
    #   return;
    #end;
    start = proc.time()[3];

    i=1;

    if (lam>0.9) {
       print('\n*** Warning: Large damping factor in use. Turn on plotting to monitor the net similarity. The algorithm will change decisions slowly, so consider using a larger value of convits');
    }

    # Check that standard arguments are consistent in size

    # Construct similarity matrix
    S=s
    N=dim(s)[1]

    # In case user did not remove degeneracies from the input similarities,
    # avoid degenerate solutions by adding a small amount of noise to the
    # input similarities
    #if (!nonoise) {
    #    rns=randn('state'); randn('state',0);
    #    S=S+(eps*S+realmin(class(s))*100).*rand(N,N);
    #    randn('state',rns);
    #}

    if (!nonoise) {
        rns=runif(N*N, min=0, max=1);
        rnm=matrix(rns, N, N);
        S=S+(2e-16*S)*rnm;
        
    }

    # Place preferences on the diagonal of S
    #if length(p)==1 for i=1:N S(i,i)=p; end;
    #else for i=1:N S(i,i)=p(i); end;
    #end;

    diag(S) <- p

    # Allocate space for messages, etc
    dS=diag(S); A=matrix(0, N,N); R=matrix(0, N,N); t=1;

    if (plt) {
        netsim=rep(0,maxits+1);
    }

    if (details) {
        idx   =matrix(0, N, maxits+1);
        netsim=rep(0,maxits+1); 
        dpsim =rep(0,maxits+1); 
        expref=rep(0,maxits+1); 
    }

    # Execute parallel affinity propagation updates
    e=matrix(0, N,convits); dn=F; i=0;

    #if symmetric, ST=S; else ST=S'; end; # saves memory if it's symmetric
    ST=S
    while (!dn ) {
        i=i+1; 

        # Compute responsibilities
        A=t(A); R=t(R);
        for (ii in c(1:N) ) {
            old = R[,ii];
            AS = A[,ii] + ST[,ii]; 
            mi = maxValueIndex(AS) # [Y,I]=max(AS);
            Y  = mi[1]
            I  = mi[2] 
            AS[I]=-Inf;
            
            mi = maxValueIndex(AS) # [Y2,I2]=max(AS);
            Y2 = mi[1]
            I2 = mi[2] 

            R[, ii]=ST[, ii]-Y;
            R[I,ii]=ST[I,ii]-Y2;
            R[, ii]=(1-lam)*R[,ii]+lam*old; # Damping
        }
        A=t(A); R=t(R);
        
        # Compute availabilities
        for (jj in c(1:N) ) {
            old = A[,jj];

            #Rp = max(R[,jj],0); 
            # positve responsibilities candidate exemplar jj receives
            #
            Rp  = ifelse(R[,jj]<0, 0, R[,jj])

            Rp[jj] = R[jj,jj];
            A[,jj] = sum(Rp)-Rp;
            dA     = A[jj,jj]; 

            #A[,jj] = min(A(:,jj),0); 
            A[,jj]  = ifelse(A[,jj]>0, 0, A[,jj])

            A[jj,jj] = dA;
            A[,jj]   = (1-lam)*A[,jj] + lam*old; # Damping
        }
        
        # Check for convergence
        AB = diag(A)+diag(R);
        E  = AB>0;
        e[,(i-1)%%convits+1]=E; #e(:,mod(i-1,convits)+1)=E;
        K  = sum(E);
        #K

        if (i>=convits) {
            #se=sum(e,2); sum for each row
            se=apply(e, 1, sum)
            unconverged=(sum((se==convits)+(se==0)) !=N);
            if ( (!unconverged & (K>0)) |(i==maxits) ) {
               dn=T;
            }
        }

        # Handle plotting and storage of details, if requested
        if (plt | details ) {
            if (K==0) {
                tmpnetsim=NA; tmpdpsim=NA; tmpexpref=NA; tmpidx=NA;
            } else {
                #I=find(E); [tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c);

                I=c(1:length(E))[E];# [tmp c]=max(S[,I],[],2); 
                Si  = S[,I]
                tmp = apply(Si, 1, max)
                c   = apply(Si, 1, which.max)

                c[I]=c(1:K); tmpidx=I[c];
                #tmpnetsim=sum(S((tmpidx-1)*N+[1:N]'));
                tmpnetsim =sum( S[(tmpidx-1)*N+c(1:N) ] )

                tmpexpref=sum(dS[I]); tmpdpsim=tmpnetsim-tmpexpref;
            }
        }

        if (details) {
            netsim[i]=tmpnetsim; dpsim[i]=tmpdpsim; expref[i]=tmpexpref;
            idx[,i]=tmpidx;
        }

        #if plt & toc>1, tic;
        if (plt){
            netsim[i]=tmpnetsim;
            tmp  = c(1:i); 
            #tmpi= find(~isnan(netsim(1:i)));
            tsel = !is.na(netsim[1:i])
            if(sum(tsel)>0) {
               tmpi = tmp[ tsel ];
               plot(tmp[tmpi],netsim[tmpi], xlab="# Iterations",
                    ylab="Fitness (net similarity) of quantized intermediate solution");
            }
        }

       if(debug){
         print(round(A,3))
         print(round(R,3))
         print(round(A+R,3))
       }

    } #while

    #I=find(diag(A+R)>0); 
    ARSel = diag(A+R)>0;
    I     = c(1:length(ARSel) )[ARSel]
    K     = length(I); # Identify exemplars
    if (K>0) {
        #
        # notice that the different convention in R and Matlab
        # row-based operation is 1 in R but is 2 in matlab
        # vice versa for the column-based operation
        #
        # [tmp c]=max(S(:,I),[],2); find max for each row
        # 
        # Si is a similariy matrix with columns as exemplars and rows as genes
        #  now we need assign each row element with a cluster, one of exemplars
        #
        Si = S[,I]
        if(K>1) {
          tmp = apply(Si, 1, max)       # max similarity with exemplars
          c   = apply(Si, 1, which.max) # indices of exemplars in I
        }else{
          tmp = Si
          c   = rep(1, length(Si) )
        }
    
        c[I]=c(1:K); # handle exemplars themselves, the cluster of an exemplar is itself
        
        # R.1) Refine the final set of exemplars and clusters and return results
        for (k in c(1:K) ) {
            #ii=find(c==k); 
            #[y j]=max(sum(S(ii,ii),1)); # column based operation

            # ii are the indices of the nodes assigned to the k-th exemplar
            # S[ii,ii] is the module subnetwork initially labeled as k-th exemplar,
            #   but will be updated with the node of maximum connectivity in the subnetwork
            #
            #>ii
            # 1a 1b 1c 1d 
            # 1  2  3  4 
            # > S[ii,ii]
            #    1a 1b 1c 1d
            # 1a  9  1  1  1
            # 1b  1  1  0  0
            # 1c  1  0  1  0
            # 1d  1  0  0  1

            ii = which(c==k);
            if(length(ii)>1) {
               Ssum = apply(S[ii,ii], 2, sum);
               j    = which.max(Ssum)
            }else{
               Ssum = sum(S[ii,ii])
               j    = c(1)
            }                

            # update the module label, not necessarily k-th exemplar
            I[k] = ii[ j[1] ];
        }
        
        # R.2) reassign the module labels
        #
        #[tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c);
        Si = S[,I]
        if(K>1) {
          tmp = apply(Si, 1, max)
          c   = apply(Si, 1, which.max)
        }else{
          tmp = Si
          c   = rep(1, length(Si) )
        }

        c[I]=c(1:K); # for exemplars themselves
        tmpidx=I[c];

        tmpnetsim=sum(S[(tmpidx-1)*N+c(1:N) ]); 
        tmpexpref=sum(dS[I]);

    }else {
        tmpidx=rep(NA, N); tmpnetsim=NA; tmpexpref=NA;
    }

    if (details) {
        netsim[i+1]=tmpnetsim; netsim=netsim[1:(i+1)];
        dpsim[i+1]=tmpnetsim-tmpexpref; dpsim=dpsim[1:(i+1)];
        expref[i+1]=tmpexpref; expref=expref[1:(i+1)];
        idx[,i+1]=tmpidx; idx=idx[,1:(i+1)];
    }else{
        netsim=tmpnetsim; dpsim=tmpnetsim-tmpexpref;
        expref=tmpexpref; idx=tmpidx;
    }

    if (plt | details) {
        print(paste("Number of exemplars identified: ", K, " (for ", N, " data points)") );
        print(paste("Net similarity: ",tmpnetsim) );
        print(paste("  Similarities of data points to exemplars: ",dpsim[length(dpsim)], "\n") );
        print(paste("  Preferences of selected exemplars: ",tmpexpref, "\n"));
        print( paste("Number of iterations: ",i, "\n\n") );
            #print('Elapsed time: %f sec\n',etime(clock,start));
        totaltime  <- proc.time()[3]-start
        print( paste("Elapsed time: ", round(totaltime,3), " sec; ", i, " iterations") );
    }

    if (unconverged) {
        print("\n*** Warning: Algorithm did not converge. Activate plotting\n");
        print("    so that you can monitor the net similarity. Consider\n");
        print("    increasing maxits and convits, and, if oscillations occur\n");
        print("    also increasing dampfact.\n\n");
    }

    totaltime  <- proc.time()[3]-start
    print( paste(K, " exemplars, Elapsed time: ", round(totaltime,3), " sec; ", i, " iterations") );



   if(!overlappingModule) {
       modulecolor = assignModuleColor(idx, minsize1=minModuleSize-1,  
                                       anameallmodules=F, auseblackwhite=F, 
                                       useNumberAsLabel=useNumberAsLabel, 
                                       startlabel=startlabel)
       # locate the exemplar for each module
       mc = as.character(modulecolor)
       mcOrder = order(mc)
       mcD = mc[mcOrder]
       idxD= idx[mcOrder]
       explIdx = findBoundary(mcD)
       me = cbind(mcD[idxD], idxD[idxD])
       colnames(me) <- c("module", "exemplar")

       return (list(modulecolor, idx, me) )
   }


    #return (idx)
    #return (tmpidx) #(idx)

}


#------------------------------------------------------------------------------------------------------------------
#Function: cluster idnetification based on dynamic programming
#
#
#
#Input:
#
#    adjmatrix           ~ adjacency matrix
#    hierclust           ~ hierarchical clustering structure
#    minModuleSize       ~  min size of module
#
#    Choices of clustering coefficient (L= # of links): 
#       1) L/n,   favour bigger modules 
#       2) L/n^2, favour small clusters as small modules are likely to have high CC
#       3) L/(n*1.5+const), no bias to big or small modules
#
#       4) foldchange_cc is used to terminate the search when the module coherence is already small
#
#       5) apower is used for defining module efficiency CC*n^apower, so smaller apower favour small modules
#           and so bigger apower favour big modules
#
#Output:  module colors for every node
# 
#------------------------------------------------------------------------------------------------------------------
identifyModuleByDynamicProgram = function(adjmatrixSorted, hierclust, minModuleSize=5, npower=1.5, nconst=50, maxModules=169, heicutoff=0.99, useNumberAsLabel=F, startlabel=1, foldchange_cc=100.0, heistep = 0.02)
{
 
  nodes = dim(adjmatrixSorted)[1]
  cc    = matrix(0, nodes, nodes)

  # find break points to avoid false modules
  #
  #staticClu  = cutTreeStatic(hiercluster=hierclust, heightcutoff=heicutoff, minsize1=0)

  # make sure that there is at least two clusters after the static cut
  curheicutoff=heicutoff  
  no.bounds   = 0
  while(no.bounds <=1 & curheicutoff>0){
    staticClu  = cutTreeStatic(hiercluster=hierclust, heightcutoff=curheicutoff, minsize1=0)
    ostaticClu = staticClu[ hierclust$order] #ordered
    bounds     = findBoundary(ostaticClu)
    bounds     =c(bounds, nodes+1)
    no.bounds  = length(bounds)
    curheicutoff = curheicutoff - heistep
  }
  if (curheicutoff<=0 | no.bounds <=1){     
     return(rep("grey", nodes))
  }

  breakpoints= rep(0, nodes)
  for(i in c(1:(no.bounds-1) )){
     # index of elements belong to the current cluster i
     # as bounds[i+1] is the start point of next cluster, 
     # bounds[i+1]-1 is the end point of the cluster i
     idx = c(bounds[i]:(bounds[i+1]-1))
     breakpoints[idx] = bounds[i+1]-1
  }
  breakpoints

  #assign colors for modules
  module.assign = rep(0, nodes)
  module.cnt=1

  # 1) number of links of a network with nodes from 1 to j, use an accumulated ways to save time
  i= 1
  for (j in c((i+1):nodes) ){
     cc[i,j] = cc[i, j-1] + sum(adjmatrixSorted[j, c(i:j)] )*2
  }

  # 2) number of links of a network with nodes from  i to j
  for (i in c(2:nodes) ){
    n = cumsum( adjmatrixSorted[i-1, ] )
    n = n-n[i-1]
    
    #if (i+1>nodes){
    #   break
    #}
    if (i==nodes){
       idx = c(nodes)
    }else{
       idx = c((i+1):nodes)
    }

    cc[i, idx ] = cc[i-1, idx] - n[idx]*2
  }

  # 2.1) enforce breakpoints
  for (i in c(1:(nodes-1)) ){
      if(breakpoints[i] <nodes){
          bidx         = c((breakpoints[i]+1):nodes)
      }else{
          bidx         = c(breakpoints[i]:nodes)
          break
      }
      cc[i, bidx ] = 0 
  }

  # 3) normalization by the total number of nodes in each network
  for (i in c(1:nodes) ){
    n       = rep(1, nodes)
    idx     = c(i:nodes)
    n[idx]  = idx-i+1
    n=ifelse(n==0, 1, n)

    #apower = n**npower + nconst
    #apower = n + nconst
    #cc[i, ] = cc[i,]/(apower)
    #cc[i, ] = cc[i,]/n
    #cc[i,]  = ifelse(n<minModuleSize, 0, cc[i,]/n)
    #cc[i, ] = cc[i,]/(apower)    
    #cc[i, ] = cc[i,]/(n-1)/((n)^0.5)

    cc[i, ] = cc[i,]/(n-1)/((n)^(1-npower))
    cc[i,]  = ifelse(n<minModuleSize, 0, cc[i,])
  }

  if(F) {
  # module identification
  cc1d = -as.numeric(cc) # convert a matrix into a vector
  maxIdx     = order(cc1d)
  Mc = cc1d[ maxIdx[1] ]
  # location of the max in the original matrix
  Ii  = maxIdx[1] %% nodes  #Ii as row of max
  if(Ii==0){#last column
    Ii=nodes
  }
  i = as.integer(maxIdx[1]/nodes) #i as column of max
  fi = maxIdx[1]/nodes
  if (fi>i){
      i=i+1
  }
  }

  # MCI[1, ]: max value and MCI[2, ] is index
  #
  mcI = apply(cc, 2, maxValueIndex)
  mc  = mcI[1, ]
  I   = mcI[2, ]
  Mci = maxValueIndex(mc)
  Mc  = Mci[1]
  i   = Mci[2]
  Ii  = I[i]

  #maxMc=max(mc)
  maxMc=mean(cc)
  cutMc= maxMc/foldchange_cc

  mresult = NULL
  #while(Mc > 0 & module.cnt <maxModules){ #original definition
  while(Mc > cutMc & module.cnt <maxModules){
      
     nres = c(Ii, i, Mc, i-Ii)
     mresult = rbind(mresult, nres)

     mstr = paste("max_cc=", as.character(Mc), " size=", as.character(i-Ii+1),"\n",sep="")
     cat(mstr)

     #assign modules
     module.size = i-Ii+1
     if(module.size >=minModuleSize){
        #assign module lable
        module.assign[c(Ii:i)] = rep(module.cnt, module.size)
        module.cnt = module.cnt + 1
     }

     mstr = paste("max_cc=", as.character(Mc), " size=", as.character(i-Ii+1),"\n",sep="")
     cat(mstr)

     # reset the nodes already in the module
     cc[1:i, c(Ii:nodes) ] = 0 #original

     #cc[1:i, c(Ii:nodes) ] = 0
     #cc[Ii:i, c(Ii:nodes) ] = 0
     #cc[1:Ii, c(Ii:i) ] = 0

     if(F) {
     cc1d  = -as.numeric(cc)
     maxIdx= order(cc1d)
     Mc    = cc1d[ maxIdx[1] ]
     # location of the max in the original matrix
     Ii  = maxIdx[1] %% nodes  #Ii as row of max
     if(Ii==0){#last column
       Ii=nodes
     }

     i = as.integer(maxIdx[1]/nodes) #i as column of max
     fi = maxIdx[1]/nodes
     if (fi>i){
        i=i+1
     }
     }
 
     # MCI[1, ]: max value and MCI[2, ] is index
     #
     mcI = apply(cc, 2, maxValueIndex)
     mc  = mcI[1, ]
     I   = mcI[2, ]
     Mci = maxValueIndex(mc)
     Mc  = Mci[1]
     i   = Mci[2]
     Ii  = I[i]
  }

  #mcol     = dim(mresult)[2]
  #selected = mresult[,mcol]>=minModuleSize
  #sel.results = mresult[selected, ]

  colcode.reduced.order = assignModuleColor(module.assign, minsize1=minModuleSize-1,  
               anameallmodules=F, auseblackwhite=F, useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)

  #re-order to the normal one with sequential signleton index
  recov.order     = order( hierclust$order)
  colcode.reduced = colcode.reduced.order[recov.order]

  rm(cc)
  collect_garbage()

  return (colcode.reduced)

}

maxValueIndex = function(mvect){
   index  = c(1:length(mvect))
   maxval = max(mvect, na.rm=T)
   sel    = maxval==mvect
   sel    = ifelse(is.na(sel), F, sel)
   mindex = index[sel]
   return ( c(maxval, mindex[1]) )
}
#a=cbind(c(1:3), c(5,6,7), c(9,10,8))
#b=apply(a, 1, maxValueIndex)

# a= 1 1 2 2 3 3 3
# findBoundary(a) => [1] 1 3 5
#
findBoundary=function(vect){
   no.elements = length(vect)
   esel = vect[no.elements]==vect[1]
   esel = ifelse(is.na(esel), F, esel)

   #if(esel){
   #  return (1);
   #}

   shifted= c(vect[no.elements], vect[1:(no.elements-1)] )
   sel = shifted != vect

   # handle NA
   sel = ifelse(is.na(sel), F, sel)

   sel[1] = TRUE

   return ( c(1:no.elements)[sel] )
}

findSame=function(vect){
   no.elements = length(vect)
   if(vect[no.elements]==vect[1]){
     return (T);
   }

   shifted= c(vect[2:no.elements], vect[1])
   sel = shifted == vect

   return (sel)
}


findSameByDistance = function(vect, distance=0){
   no.elements = length(vect)
   if(vect[no.elements]==vect[1]){
     return (1);
   }

   shifted= c(vect[2:no.elements], vect[1])
   sel = abs(shifted- vect)<=distance
   return (sel)
}


# assume vect is sorted, we try to find the indices of the signulars
#
findSingularsIdx=function(vect){
  nelems    = length(vect)
  startIdx  = findBoundary(vect)
  no.members= c(startIdx[-1],nelems+1) - startIdx  
  singleIdx  = startIdx[no.members==1]
  return(singleIdx)
}

# for unsorted vectors
findUniqueIdx=function(vect){

   no.elements = length(vect)
   if(no.elements<=1){return (1)}

   vectIdxed   = cbind(vect, c(1:no.elements) )
   mo = order(vectIdxed[,1])
   merged = vectIdxed[mo,]

   bounds = findBoundary(merged[,1])

   uniqueIdx = as.integer(merged[,2])[bounds]

   # keep original order
   return ( sort(uniqueIdx) )
}

findUniqueIndex=function(vect){
 res = findUniqueIdx(vect)
 return (res)
}

# for unsorted vectors
findUniqueIdx2=function(vect){
   no.elements = length(vect)
   uniqueIdx      = c(1)
   for (i in c(2:no.elements) ){
       found = is.element(vect[i], vect[c(1:(i-1))])
       if(!found){
           uniqueIdx= c(uniqueIdx, i)  
       }
   }
   return ( uniqueIdx)
}


getElementIndex = function(vect, element){
  
   no.elements = length(vect)
   msel = vect == element
   if (sum(msel)>0) {
      return( c(1:no.elements)[msel] ) 
   } else{return (-1); }

}



#------------------------------------------------------------------------------------------------------------------
#Function: cut hierarchical clusering tree, base don internal structure of ordered dendrogram 
#          use the information of dendrogram height with the mean height of a module to determine the cut points,
#             if this doesn't split  the module, try the differences with (max_height+mean_height)/2
#          the whole process can be iterated until number of clusters become stable (if deepSplit is set as TRUE)
#Input:
#    hierclust           ~  hierarchical clustering object
#    deepSplit           ~  perform iterative searching of sub clusters if set to be TRUE, otherwise, do
#    minModuleSize       ~  min size of module
#    minAttachModuleSize ~  min size of module to be attached, as a major module with larger size
#                                startpoint of current module 
#Output:  module colors for every node
# 
#------------------------------------------------------------------------------------------------------------------

cutreeDynamic = function(hierclust, maxTreeHeight=1, deepSplit=TRUE, minModuleSize=50, 
                         minAttachModuleSize=50, maxMeanHeight=0.95, nameallmodules=FALSE, 
                         useblackwhite=FALSE, useNumberAsLabel=F, startlabel=1, 
                         minTailRunlength=8, fimgPCsHier=NULL, 
                         use_PCA=TRUE, pc_heightcutoff=0.3, data=NULL, coexppRes=NULL, pv_cut=0.05, use_Absolute_PC_Corrl=TRUE)
{
    #dim(hierclust$merge)
    #length(hierclust$height)
    #length(hierclust$order)
    #hierclust$merge[1:10, ]
    #round(hierclust$height[1:20], 2)
    #orderHei   = hierclust$height[hierclust$order]
    #orderMerge = hierclust$merge[, ]
    #plclust(hierclust, labels=F, xlab="",ylab="",main="",sub="")

    if(maxTreeHeight >=1){
      #staticCutCluster = cutTreeStatic(hiercluster=hierclust, heightcutoff=0.99, minsize1=minModuleSize)
      staticCutCluster = cutTreeStatic_new(hiercluster=hierclust, heightcutoff=0.99, minsize1=minModuleSize)

    }else{
      #staticCutCluster = cutTreeStatic(hiercluster=hierclust, heightcutoff=maxTreeHeight, minsize1=minModuleSize)
      staticCutCluster = cutTreeStatic_new(hiercluster=hierclust, heightcutoff=maxTreeHeight, minsize1=minModuleSize)
    }      

    #get tree height for every singleton
    #node_index   tree_height
    demdroHeiAll= rbind( cbind(hierclust$merge[,1], hierclust$height), cbind(hierclust$merge[,2], hierclust$height) )

    #singletons will stand at the front of the list
    myorder = order(demdroHeiAll[,1])

    #get # of singletons
    no.singletons = length(hierclust$order)

    #> demdroHei.sort[1:10,]
    #       [,1]      [,2]
    #[1,] -3000 0.7943783
    #[2,] -2999 0.7863851
    demdroHeiAll.sort = demdroHeiAll[myorder, ]
    demdroHei.sort    = demdroHeiAll.sort[c(1:no.singletons), ]

    #finally, we got tree height for each of the singleton inorder of 1 to no.singletons
    #     [,1]      [,2]
    #[1,]    1 0.8389184
    #[2,]    2 0.8772433
    #[3,]    3 0.8308482
    demdroHei      = demdroHei.sort[seq(no.singletons, 1, by=-1), ]
    demdroHei[,1]  = -demdroHei[,1]

    # combine with prelimilary cluster-cutoff results
    demdroHei  = cbind(demdroHei, as.integer(staticCutCluster))

    # re-order the order based on the dendrogram order hierclust$order
    demdroHei.order = demdroHei[hierclust$order, ]

    # get start and end posiiton of every cluster
    # [1,]  173  315
    # [2,]  676  793
    static.clupos = locateCluster(demdroHei.order[, 3])
    static.no     = dim(static.clupos)[1]
    static.sizes  = static.clupos[,2] - static.clupos[,1] + 1
  
    #split individual cluster if there are sub clusters embedded
    if(F){
        clusterDemdroHei=mydemdroHei.order;
        cminModuleSize       = minModuleSize;
        cminAttachModuleSize = minAttachModuleSize;
    }

        clupos = NULL
        for (i in c(1:static.no)){
           mydemdroHei.order = demdroHei.order[ c(static.clupos[i,1]:static.clupos[i,2]), ] #index to [1, clusterSize]
           mydemdroHei.order[, 1] = mydemdroHei.order[, 1] - static.clupos[i, 1] + 1

           #cat("Cycle ", as.character(mcycle), "cluster (", static.clupos[i,1], static.clupos[i,2], ")\n")
           #cat("i=", as.character(i), "\n")
           #iclupos = processInvididualCluster(mydemdroHei.order, cminModuleSize=minModuleSize, cminAttachModuleSize=minAttachModuleSize)

           # new approach
           iclupos1 = processInvididualCluster_new(clusterDemdroHei=mydemdroHei.order, cminModuleSize=minModuleSize, minTailRunlength=minTailRunlength, nstep=8, pv_cut=pv_cut)
           iclupos= attachMinorToMajor(clusterDemdroHei=mydemdroHei.order, cluposI=iclupos1, cminAttachModuleSize=minAttachModuleSize, maxMeanHeight=maxMeanHeight)

           #iclupos2 = attachMinorToMajor(clusterDemdroHei=mydemdroHei.order, cluposI=iclupos, cminAttachModuleSize=minAttachModuleSize, maxMeanHeight=maxMeanHeight)
           #iclupos  = removeLooseModules(clusterDemdroHei=mydemdroHei.order, cluposI=iclupos2, maxMeanHeight=maxMeanHeight)

           if(is.null(iclupos)) {next}

           iclupos[,1] = iclupos[,1] + static.clupos[i, 1] -1 #recover the original index
           iclupos[,2] = iclupos[,2] + static.clupos[i, 1] -1

           clupos  = rbind(clupos, iclupos) #put in the final output buffer
        }

        #clupos  = pruneModules(clusterDemdroHei=demdroHei.order, cluposI=clupos, maxMeanHeight=maxMeanHeight)

    final.cnt = dim(clupos)[1]

    # try to make use of all the colors
    msizes = clupos[, 2] - clupos[, 1] + 1
    allcolors = setdiff(colors(), "grey")
      
    if(length(allcolors)>=final.cnt){minModuleSize2=min(msizes)
    } else {
      minModuleSize2=-sort(-msizes)[length(allcolors)]
    }

    #assign colors for modules
    module.assign = rep(0, no.singletons)
    module.cnt=1

    for (i in c(1:final.cnt))
    {
       sdx = clupos[i, 1] #module start point
       edx = clupos[i, 2] #module end point

       module.size = edx - sdx +1 

       #if(module.size <minModuleSize2){
       #  next
       #}

       #assign module lable
       module.assign[sdx:edx] = rep(module.cnt, module.size)

       #update module label for the next module
       module.cnt = module.cnt + 1
    }
    #clupos  = pruneModules(clusterDemdroHei=demdroHei.order, cluposI=clupos, maxMeanHeight=maxMeanHeight)

    #colcode.reduced.order = assignModuleColor(module.assign, minsize1=minModuleSize2,  anameallmodules=nameallmodules, 
    #                         auseblackwhite=useblackwhite, useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)
    #tapply(demdroHei.order[1:no.singletons, 2], colcode.reduced.order, mean)

    # no restriction on size
    #
    colcode.reduced.order = assignModuleColor(module.assign, minsize1=1,  anameallmodules=TRUE, 
                             auseblackwhite=FALSE, useNumberAsLabel=TRUE, startlabel=1)

    #re-order to the normal one with sequential signleton index
    recov.order     = order( demdroHei.order[,1])
    mcolcode = colcode.reduced.order[recov.order]

    plotDendrogramModuleLabels(mdendro=hierclust, modulecolors=mcolcode, save2file=NULL, plotLabels=FALSE)

    mcolcode2=IterativeModuleMergeByPC(hierclust=hierclust, modules=as.character(mcolcode), 
      data=data, coexppRes=coexppRes, use_PCA=use_PCA, fimg=NULL, useNumberAsLabel=TRUE, 
      startlabel=1, heightcutoff=pc_heightcutoff, useAbsoluteValue=use_Absolute_PC_Corrl)

    if(useNumberAsLabel) {
      retcolors = mcolcode2
    } else {
      retcolors = reassignModuleNames(mcolcode2, minmodulesize=minModuleSize, anameallmodules=nameallmodules, 
                    auseblackwhite=auseblackwhite, useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)
    }  

    plotDendrogramModuleLabels(mdendro=hierclust, modulecolors=retcolors, save2file=NULL, plotLabels=FALSE)


    # hierarchical clustering of module-module correlations
    #
    moduleClusterByPC_plot(modules=retcolors, data=data, coexppRes=coexppRes, use_PCA=use_PCA, fimg=fimgPCsHier) 

    if(1==2){
        aveheight = averageSequence(demdroHei.order[,2], 2)
        procHei = demdroHei.order[,2]-mean(demdroHei.order[,2])

        par(mfrow=c(3,1), mar=c(0,0,0,0) )
        plot(hierclust, labels=F, xlab="",ylab="",main="",sub="",axes = F)
        #barplot(demdroHei.order[,2]-min(demdroHei.order[,2]), 
        #        col= "black", space=0,
        #        border=F,main="", axes = F, axisnames = F)
        #barplot(aveheight-mean(aveheight), 
        barplot(procHei, 
                col= "black", space=0,
                border=F,main="", axes = F, axisnames = F)
        barplot(height=rep(1, length(colcode.reduced)), 
                col= as.character(colcode.reduced[hierclust$order]), space=0,
                border=F,main="", axes = F, axisnames = F)
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    }

    retcolors 
}


#input is the cluster demdrogram of an individual cluster, we want to find its embbeded subclusters
#execution order: mean-height ==> (mean+max)/2 ==> (mean+min)/2
#useMean: =0 ~ use mean-height   as calibation line
#         =1 ~ use (mean+max)/2  as calibation line to detect relatively a small cluster sitting on the head of a bigger one,
#                      so mean-height is too low to detect the two modules.
#         =-1~ use (mean+min)/2  as calibation line to detect relatively a small cluster sitting on the tail of a bigger one,
#                      so mean-height & (mean+max)/2 are too high to detect the two modules

#processInvididualCluster = function(clusterDemdroHei, cminModuleSize=50, cminAttachModuleSize=100, minTailRunlength=12, useMean=0){
processInvididualCluster = function(clusterDemdroHei, cminModuleSize=50, cminAttachModuleSize=100, minTailRunlength=12, useMean=0){
    #for debug: use all genes
    #clusterDemdroHei =demdroHei.order

    no.cnodes = dim(clusterDemdroHei)[1]
    
    cmaxhei   = max(clusterDemdroHei[, 2])
    cminhei   = min(clusterDemdroHei[, 2])
    
    cmeanhei  = mean(clusterDemdroHei[, 2])
    cmidhei = (cmeanhei + cmaxhei)/2.0
    cdwnhei = (cmeanhei + cminhei)/2.0

    #if(mean(clusterDemdroHei[, 2]) < 0.2){
    #   new.clupos = cbind( 1, no.cnodes)
    #   return (new.clupos)
    #}

    if (useMean==1){
        comphei = cmidhei
    }else if (useMean==-1){
        comphei = cdwnhei
    }else{ #normal case
        comphei = cmeanhei
    }
        
    # compute height diffrence with mean height
    heidiff       = clusterDemdroHei[,2] - comphei
    heidiff.shift = shiftSequence(heidiff, -1)

    # get cut positions
    # detect the end point of a cluster, whose height should be less than meanhei 
    #  and the node behind it is the start point of the next cluster which has a height above meanhei
    cuts.bool = (heidiff<0) & (heidiff.shift > 0)
    cuts.bool[1]         = TRUE
    cuts.bool[no.cnodes] = TRUE
    
    if(sum(cuts.bool)==2){
          if (useMean==0){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=1)
          }else if(useMean==1){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=-1)          
          }else{
             new.clupos = rbind(c(1, no.cnodes))
          }
          return (new.clupos)
    }

    #a good candidate cluster-end point should have significant # of ahead nodes with head < meanHei
    cutindex =c(1:no.cnodes)[cuts.bool]
    no.cutps = length(cutindex)
    runlens  = rep(999, no.cutps)
    cuts.bool2 = cuts.bool
    for(i in c(2:(no.cutps-1)) ){
       seq = c( (cutindex[i-1]+1):cutindex[i] )
       runlens[i] = runlengthSign(heidiff[seq], leftOrright=-1, mysign=-1)

       if( (runlens[i]<minTailRunlength) & (runlens[i]<cminModuleSize) ){
          #cat("run length=", runlens[i], "\n")
          cuts.bool2[ cutindex[i] ] = FALSE
       }
    }

    #attach SMALL cluster to the left-side BIG cluster if the small one has smaller mean height
    cuts.bool3=cuts.bool2
    if(sum(cuts.bool2) > 3) {
       curj = 2
       while (1==1){
           cutindex2 =c(1:no.cnodes)[cuts.bool2]
           no.clus = length(cutindex2) -1
           if (curj>no.clus){
              break
           }
           pre.sdx = cutindex2[ curj-1 ]+1 #previous module start point
           pre.edx = cutindex2[ curj ] #previous module end   point
           pre.module.size = pre.edx - pre.sdx +1 
           pre.module.hei  = mean(clusterDemdroHei[c(pre.sdx:pre.edx) , 2])
         
           cur.sdx = cutindex2[ curj ]+1 #previous module start point
           cur.edx = cutindex2[ curj+1 ] #previous module end   point
           cur.module.size = cur.edx - cur.sdx +1 
           cur.module.hei  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])

           #merge to the leftside major module, don't change the current index "curj"
           #if( (pre.module.size >minAttachModuleSize)&(cur.module.hei<pre.module.hei)&(cur.module.size<minAttachModuleSize) ){
           if( (cur.module.hei<pre.module.hei)&(cur.module.size<cminAttachModuleSize) ){
                cuts.bool2[ cutindex2[curj] ] = FALSE
           }else{ #consider next cluster
                curj = curj + 1
           }
       }#while
   }#if 

   cutindex2 =c(1:no.cnodes)[cuts.bool2]
   no.cutps = length(cutindex2)

    #we don't want to lose the small cluster at the tail, attch it to the previous big cluster
    #cat("Lclu= ", cutindex2[no.cutps]-cutindex2[no.cutps-1]+1, "\n")
    if(no.cutps > 2){
      if( (cutindex2[no.cutps] - cutindex2[no.cutps-1]+1) < cminModuleSize ){
        cuts.bool2[ cutindex2[no.cutps-1] ] =FALSE  
      }
    }

    if(1==2){
        myseqnce = c(2300:3000)
        cutdisp = ifelse(cuts.bool2==T, "red","grey" )
        #re-order to the normal one with sequential signleton index
        par(mfrow=c(3,1), mar=c(0,0,0,0) )
        plot(hierclust, labels=F, xlab="",ylab="",main="",sub="",axes = F)
        barplot(heidiff[myseqnce],
                col= "black", space=0,
                border=F,main="", axes = F, axisnames = F)
        barplot(height=rep(1, length(cutdisp[myseqnce])), 
                col= as.character(cutdisp[myseqnce]), space=0,
                border=F,main="", axes = F, axisnames = F)
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    }

   cutindex2  = c(1:no.cnodes)[cuts.bool2]
   cutindex2[1]=cutindex2[1]-1 #the first 
   no.cutps2  = length(cutindex2)

   if(no.cutps2 > 2){
     new.clupos = cbind( cutindex2[c(1:(no.cutps2-1))]+1, cutindex2[c(2:no.cutps2)] )
   }else{
     new.clupos = cbind( 1, no.cnodes)
   }

   if ( dim(new.clupos)[1] == 1 ){   
          if (useMean==0){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=1)
          }else if(useMean==1){
             new.clupos=processInvididualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=-1)          
          }   
   }
   new.clupos
}


processInvididualCluster_new = function(clusterDemdroHei, cminModuleSize=50, minTailRunlength=5, nstep=8,  pv_cut=0.05){
    #for debug: use all genes
    #clusterDemdroHei =demdroHei.order

    no.cnodes = dim(clusterDemdroHei)[1]

    #pv_cut2= pv_cut/no.cnodes

    if(no.cnodes<cminModuleSize) {return (rbind(c(1, no.cnodes)) )}

    cmeanhei  = mean(clusterDemdroHei[, 2]) ; cSDhei = sd(clusterDemdroHei[, 2])
    cmaxhei   = max(clusterDemdroHei[, 2]) #cmeanhei - cSDhei
    cminhei   = min(clusterDemdroHei[, 2]) # +  cSDhei

    hei.shift = shiftSequence(clusterDemdroHei[,2], -1)
    heiNgbDiff= clusterDemdroHei[,2]-hei.shift
    meanDiff  = mean(abs(heiNgbDiff)); sdDiff=sd(abs(heiNgbDiff)); maxDiff= max(abs(heiNgbDiff))
    heidiff_pv= pnorm(q=abs(heiNgbDiff),mean=meanDiff, sd=sdDiff, lower.tail=FALSE, log.p=FALSE)

    heidiff_pv.shift = shiftSequence(heidiff_pv, 1)

    if(min(heidiff_pv[-c(1,no.cnodes)]) >pv_cut) {return (rbind(c(1, no.cnodes)) ) }

    heidiff_pv[1]=0; heidiff_pv[no.cnodes]=0; 
    heidiff_pv.shift[1]=0; heidiff_pv.shift[no.cnodes]=0; 

    #midx = which.max(abs(heiNgbDiff)[-1]) + 1 # don't consider the first point
    midx = which.min(heidiff_pv[-c(1,no.cnodes)]) + 1 # don't consider the first point
    cmaxhei = clusterDemdroHei[midx+1, 2]
    cminhei = clusterDemdroHei[midx, 2]

    # find the best cut to derive most clusters
    #
    if(cmaxhei > cminhei) {
      delta     = (cmaxhei - cminhei)/nstep
      cuts      = seq(from=cminhei+delta, to=cmaxhei-delta, by=delta)
    } else {
      cuts = cmeanhei
    }

    no.cuts   = length(cuts)
    max_ii = -1; max_clus = 0; #max_clus=no.cnodes;

    for(ii in c(1:no.cuts) ) {
      # compute height diffrence with mean height
      heidiff       = clusterDemdroHei[,2] - cuts[ii]
      heidiff.shift = shiftSequence(heidiff, -1)

      # get cut positions
      # detect the end point of a cluster, whose height should be less than meanhei 
      #  and the node behind it is the start point of the next cluster which has a height above meanhei
      cuts.bool = (heidiff<0) & (heidiff.shift > 0)
      #cuts.bool = (heidiff<0) & (heidiff_pv<pv_cut) & ( heidiff.shift > 0) #& (heidiff_pv.shift<pv_cut)

      if(sum(cuts.bool)>max_clus){ max_clus=sum(cuts.bool); max_ii=ii}
      #if((sum(cuts.bool)<max_clus) & (sum(cuts.bool)>0) ){ max_clus=sum(cuts.bool); max_ii=ii}

      #print(paste(ii, length(clusterDemdroHei[,2]), cuts[ii], sum(cuts.bool)))
    }

    if(max_ii==-1){
      new.clupos = cbind( 1, no.cnodes)
      return (new.clupos)
    }

    heidiff       = clusterDemdroHei[,2] - cuts[max_ii]
    heidiff.shift = shiftSequence(heidiff, -1)
    cuts.bool = (heidiff<0) & (heidiff.shift > 0)
    cuts.bool[1]         = TRUE
    cuts.bool[no.cnodes] = TRUE

    if(sum(cuts.bool)<=2){
      new.clupos = cbind( 1, no.cnodes)
      return (new.clupos)
    }

    #a good candidate cluster-end point should have significant # of ahead nodes with head < meanHei
    cutindex =c(1:no.cnodes)[cuts.bool]
    no.cutps = length(cutindex)
    runlens  = rep(999, no.cutps)
    cuts.bool2 = cuts.bool

    #if(F) {
     for(i in c(2:(no.cutps-1)) ){
       seq = c( (cutindex[i-1]+1):cutindex[i] )
       runlens[i] = runlengthSign(heidiff[seq], leftOrright=-1, mysign=-1)

       #if( (runlens[i]<minTailRunlength) & (runlens[i]<cminModuleSize) ){
       if(runlens[i]<minTailRunlength){
          #cat("run length=", runlens[i], "\n")
          cuts.bool2[ cutindex[i] ] = FALSE
       }
     } 
     if(sum(cuts.bool2)<=2){
       new.clupos = cbind( 1, no.cnodes)
       return (new.clupos)
     }
   #}

   cutindex2  = c(1:no.cnodes)[cuts.bool2]
   cutindex2[1]=cutindex2[1]-1 #the first 
   no.cutps2  = length(cutindex2)

   new.clupos = cbind( cutindex2[c(1:(no.cutps2-1))]+1, cutindex2[c(2:no.cutps2)] )
   no.new.clupos = nrow(new.clupos)

   clupos = NULL
   for (i in c(1:no.new.clupos)){
       mydemdroHei = clusterDemdroHei[ c(new.clupos[i,1]:new.clupos[i,2]), ] #index to [1, clusterSize]
       mydemdroHei[, 1] = mydemdroHei[, 1] - new.clupos[i, 1] + 1

       #cat("Cycle ", as.character(mcycle), "cluster (", static.clupos[i,1], static.clupos[i,2], ")\n")
       #cat("i=", as.character(i), "\n")

       iclupos = processInvididualCluster_new(mydemdroHei,cminModuleSize = cminModuleSize, minTailRunlength=minTailRunlength, nstep=nstep,  pv_cut= pv_cut)

       iclupos[,1] = iclupos[,1] + new.clupos[i, 1] -1 #recover the original index
       iclupos[,2] = iclupos[,2] + new.clupos[i, 1] -1

       clupos  = rbind(clupos, iclupos) #put in the final output buffer
   }

   return (clupos)
}





#attach SMALL cluster to the left-side BIG cluster if the small one has smaller mean height
#
attachMinorToMajor = function(clusterDemdroHei, cluposI, cminAttachModuleSize=50, maxMeanHeight=0.90){

    no.cnodes = dim(clusterDemdroHei)[1]
    no.clu.org= nrow(cluposI)
    
    if(no.clu.org ==1) { return(cluposI)}
    clupos = cluposI

  while(TRUE) {

     maxheis    = rep(0, no.clu.org)
     meanheis   = rep(0, no.clu.org)
     newcluIDs  = rep(0, no.clu.org); newcluIDs[1] = 1
     merged.clu = rbind(clupos[1,])

     # 
     maxheis2   = rep(0, no.clu.org)
     meanheis2  = rep(0, no.clu.org)

       cnt  = 1
       for(curj in c(2:no.clu.org) ) {

           #pre.sdx = merged.clu[cnt, 1]#previous module start point
           #pre.edx = merged.clu[cnt, 2]#previous module end   point

           pre.sdx = clupos[curj-1, 1]#previous module start point
           pre.edx = clupos[curj-1, 2]#previous module end   point
           pre.module.size = pre.edx - pre.sdx +1 
           pre.module.hei  = mean(clusterDemdroHei[c(pre.sdx:pre.edx) , 2])
         
           cur.sdx = clupos[curj, 1] #previous module start point
           cur.edx = clupos[curj, 2] #previous module end   point
           cur.module.size = cur.edx - cur.sdx +1 
           cur.module.hei  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])

           meanheis[curj-1]=pre.module.hei; meanheis[curj]=cur.module.hei; 
           maxheis[curj-1] =max(clusterDemdroHei[c(pre.sdx:pre.edx) , 2]); 
           maxheis[curj]   =max(clusterDemdroHei[c(cur.sdx:cur.edx) , 2]); 

           if(curj==2) {
              maxheis2[1]  = maxheis[1]
              meanheis2[1] = meanheis[1]
           }

           # nearest merge point
           sel = maxheis[curj]<=maxheis[c(1:(curj-1))]
           if(sum(sel)>0){
              sidxS = c(1:(curj-1))[sel];  sidx = sidxS[length(sidxS) ]
              mgClu = newcluIDs[sidx]
           } else {
              mgClu = cnt
           }

           #merge to the leftside major module, don't change the current index "curj"
           #if( (pre.module.size >minAttachModuleSize)&(cur.module.hei<pre.module.hei)&(cur.module.size<minAttachModuleSize) ){
           tmpClu = c(merged.clu[cnt, 1], clupos[curj, 2])
           tmp.module.hei  = mean(clusterDemdroHei[c(tmpClu[1]:tmpClu[2]) , 2])

           cond1 = (cur.module.hei<pre.module.hei)&(cur.module.size<cminAttachModuleSize)
           cond2 = #maxheis[curj]>maxheis[curj-1] & (mgClu==cnt) &(cur.module.size<cminAttachModuleSize)
           cond2 = (cur.module.hei<pre.module.hei) & (pre.module.size<cminAttachModuleSize) #&(cur.module.size>pre.module.size)
           condM1 = (cur.module.hei<maxMeanHeight) == (tmp.module.hei<maxMeanHeight) # make sure the combination won't lead to the removal of the current good module
           condM2 = (tmp.module.hei<maxMeanHeight) # make sure the combination won't lead to the removal of the current good module
           condN  = (maxheis[curj]<maxheis[curj-1])

           if( (cond1 | cond2) & (condM1 | condM2) ){
                merged.clu[cnt, ] = c(merged.clu[cnt, 1], clupos[curj, 2])
                maxheis2[cnt] = max(clusterDemdroHei[merged.clu[cnt, 1]:merged.clu[cnt, 2] ])

           }else{ #consider next cluster
                cnt = cnt + 1

                merged.clu=rbind( merged.clu, clupos[curj, ])
                maxheis2[cnt] = maxheis[curj]
           }
           newcluIDs[curj] = cnt
       }
       
      no.clu = nrow(merged.clu)
      if( (no.clu==no.clu.org) | (no.clu==1)) {break}

      # prepare for the next run
      clupos = merged.clu
      no.clu.org = no.clu

   } # while
    
    #we don't want to lose the small cluster at the tail, attch it to the previous big cluster
    #cat("Lclu= ", cutindex2[no.cutps]-cutindex2[no.cutps-1]+1, "\n")
    if(F) {
    no.cutps = nrow(merged.clu)
    dist2end = no.cnodes-merged.clu[no.cutps, 2]
    cur.sdx = merged.clu[no.cutps, 2] #previous module start point
    cur.edx = no.cnodes #previous module end   point
    cur.module.size = cur.edx - cur.sdx +1 
    cur.module.hei  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])
    cur.module.heiM = max(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])
    rst.module.heiM = max(clusterDemdroHei[c(1:(cur.sdx-1)) , 2])

    if( (sum(cur.module.hei<meanheis)>0) & (cur.module.heiM<rst.module.heiM) ){
       merged.clu[no.cutps, 2]= no.cnodes 
    }

    }

   merged.clu
}


removeLooseModules = function(clusterDemdroHei, cluposI, maxMeanHeight=0.90){

    no.cnodes = dim(clusterDemdroHei)[1]
    no.clu.org= nrow(cluposI)
    
    if(no.clu.org ==1) { return(cluposI)}

    meanheis = rep(0, no.clu.org)
    for(curj in c(1:no.clu.org) ) {
      cur.sdx = cluposI[curj, 1] #previous module start point
      cur.edx = cluposI[curj, 2] #previous module end   point
      cur.module.size = cur.edx - cur.sdx +1 
      meanheis[curj]  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])
    }
    
    sel = meanheis < maxMeanHeight
    if(sum(sel)==0) {return (NULL) } 
    return (rbind(cluposI[sel, ]))
}

pruneModules = function(clusterDemdroHei, cluposI, maxMeanHeight=0.96){

    no.cnodes = dim(clusterDemdroHei)[1]
    no.clu.org= nrow(cluposI)

    allcolors = setdiff(colors(), "grey")
    if(no.clu.org < length(allcolors)) {return (cluposI) }

    if(no.clu.org ==1) { return(cluposI) }

    meanheis = rep(0, no.clu.org)
    for(curj in c(1:no.clu.org) ) {
      cur.sdx = cluposI[curj, 1] #previous module start point
      cur.edx = cluposI[curj, 2] #previous module end   point
      cur.module.size = cur.edx - cur.sdx +1 
      meanheis[curj]  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])
    }
    heiOrder = order(meanheis)
    sel = heiOrder <= length(allcolors)
    sel = sel & (meanheis < maxMeanHeight) 

    if(sum(sel)==0) {return (NULL) } 
    return (rbind(cluposI[sel, ]))
}


module_denrogram_stats= function(hierclust, modules){
    demdroHeiAll= rbind( cbind(hierclust$merge[,1], hierclust$height), cbind(hierclust$merge[,2], hierclust$height) )

    #singletons will stand at the front of the list
    myorder = order(demdroHeiAll[,1])

    #get # of singletons
    no.singletons = length(hierclust$order)

    #> demdroHei.sort[1:10,]
    #       [,1]      [,2]
    #[1,] -3000 0.7943783
    #[2,] -2999 0.7863851
    demdroHeiAll.sort = demdroHeiAll[myorder, ]
    demdroHei.sort    = demdroHeiAll.sort[c(1:no.singletons), ]

    #finally, we got tree height for each of the singleton inorder of 1 to no.singletons
    #     [,1]      [,2]
    #[1,]    1 0.8389184
    #[2,]    2 0.8772433
    #[3,]    3 0.8308482
    demdroHei      = demdroHei.sort[seq(no.singletons, 1, by=-1), ]
    demdroHei[,1]  = -demdroHei[,1]

    heiMean = tapply(demdroHei[,2], INDEX=modules, FUN=mean, na.rm=TRUE)
    heiSD   = tapply(demdroHei[,2], INDEX=modules, FUN=sd, na.rm=TRUE)

    res = rbind(heiMean, heiSD)
    colnames(res) = names(heiMean)
    rownames(res) = c("mean_DendroHei", "sd_DendroHei")

    return (res)
}


findClustersSignificant=function(mysequence, modulecolor)
{
  
   modnames= names( table(modulecolor) )
   mysize     = length(modulecolor)
   validseq     = rep(TRUE, mysize)
   for (each in modnames ){
      mybool = (modulecolor==each)
      mymodulesig = mysequence[mybool]
      mydiff      = abs(mymodulesig - mymodulesig[1])
      if(sum(mydiff)==0){
         validseq = ifelse(mybool==TRUE, FALSE, validseq)
      }
   }
  validseq
}


#leftOrright >0 : running length (with same sign) to right, otherwise to the left
#mysign = -1: negative value, mysign = -1: positive value
runlengthSign = function(mysequence, leftOrright=-1, mysign=-1){
    seqlen = length(mysequence)   
    if(leftOrright<0){
        pseq = rev(mysequence)
    }else{
        pseq = mysequence
    }

    if(mysign<0){ #see where the first POSITIVE number occurs
       nonezero.bool = (pseq > 0)
    }else{ #see where the first NEGATIVE number occur
       nonezero.bool = (pseq < 0)
    }
    if( sum(nonezero.bool) > 0){
      runlength = min( c(1:seqlen)[nonezero.bool] ) - 1
    }else{
      runlength = 0
    }
}


#delta >0 : shift to right, otherwise to the left
shiftSequence = function(mysequence, delta){
    seqlen = length(mysequence)
    if(delta>0){
        finalseq=c(mysequence[1:delta], mysequence[1:(seqlen-delta)])
    }else{
        posdelta = -delta 
        finalseq=c(mysequence[(posdelta+1):seqlen], mysequence[(seqlen-posdelta+1):seqlen])
    }
    finalseq
}


#no of neighbors behind and before the point used for average
averageSequence=function(mysequence, noneighbors){
   sumseq = mysequence
   for(i in c(1:noneighbors)){
        iseq = shiftSequence(mysequence, i)
        sumseq =    sumseq + iseq
        iseq = shiftSequence(mysequence, -i)
        sumseq =    sumseq + iseq
   }
   sumseq = sumseq/(1+2*noneighbors)
}

#delta >0 : shift to right, otherwise to the left
shiftSequence = function(mysequence, delta){
    seqlen = length(mysequence)
    if(delta>0){
        finalseq=c(mysequence[1:delta], mysequence[1:(seqlen-delta)])
    }else{
        posdelta = -delta 
        finalseq=c(mysequence[(posdelta+1):seqlen], mysequence[(seqlen-posdelta+1):seqlen])
    }
    finalseq
}

# compute frequence in bins for a given vector
# return a matrix with 1st column as the interval middles and 2nd as the frequences
#
freqenceInBins = function(myvect, minval=0, maxval=NULL, binsize=20000, cisvect=NULL, elementnames=NULL){
    
    maxpos        = max(myvect)
    minpos        = min(myvect)

    no.elements = length(myvect)
    
    if (!is.null(minval) ){
        minpos    = minval
    }
    
    if (!is.null(maxval) ){
        maxpos    = maxval
    }

    if(is.null(cisvect) ){
        mycisvect = rep(F, no.elements)
    }else{
        mycisvect = ifelse(cisvect==1, T, F)
        mycisvect = ifelse(is.na(mycisvect), F, mycisvect)
    }

    if(is.null( elementnames) ){
        myelementnames = as.character( c(1:no.elements) )
    } else{
        myelementnames = elementnames
    }
    
    noBins        = as.integer( (maxpos-minpos)/binsize) + 1
    # intervals
    intervals     = c(0:noBins) * binsize
    
    # middle position of intervals
    intervalsMidd = minpos + c(1:noBins) * binsize - binsize/2

    no.intervals  = length(intervalsMidd)

    freqs         = rep(0,  no.intervals)
    binElements   = rep("", no.intervals)

    binCisFreqs   = rep(0,  no.intervals)
    binCisElements= rep("", no.intervals)

    halfbinsize   = binsize/2

    # count number of elements in each bin
    for (i in c(1:no.intervals)) {
        iselL    = myvect>= intervalsMidd[i] - binsize/2
        iselR    = myvect < intervalsMidd[i] + binsize/2
        isel     = iselL & iselR
        freqs[i] = sum(isel)

        if (freqs[i] == 0){
           next
        }

        # get element names in the bin
        binElements[i] = concatenate(elementnames[isel],'; ')
        
        # look at the cisQTL genes        
        icisSel           = mycisvect & isel
        binCisFreqs[i]    = sum( icisSel )
        if (sum(icisSel) >0) {
            binCisElements[i] = concatenate(elementnames[icisSel],'; ')
        }

    }
    
    retNumerical = cbind(intervalsMidd, freqs)
    retCharacter = cbind(binElements, binCisFreqs, binCisElements)

    colnames(retNumerical) <- c("interval(middle)", "QTL count") 
    colnames(retCharacter) <- c("QTL genes", "cis-QTL count", "cis-QTL genes") 

    return ( list(retNumerical,retCharacter) )
}


#find the middle of each cluster and label the middle position with the corrsponding color
getDisplayColorSequence=function(colordered){
 mylen = length(colordered)
 colordered2 = c(colordered[1], colordered[1:(mylen-1)] )
 colordiff   = (colordered != colordered2)
 colordiff[1] = TRUE
 colordiff[mylen] = TRUE
 #mydispcolor = ifelse(colordiff==TRUE, colordered, "")
 mydispcolor = rep("", mylen)
 mytrueseq = c(1:mylen)[colordiff]
 for (i in c(1:(length(mytrueseq)-1)) ){
    midi =  (mytrueseq[i] + mytrueseq[i+1])/2
    mydispcolor[midi] = colordered[midi]
 }
 fdispcolor = ifelse(mydispcolor=="grey", "", mydispcolor)
 fdispcolor
}

# find the start and end indice of each color
#
getDisplayColorPos=function(colordered){

 mylen       = length(colordered)
 myidx       = c(1:mylen)
 colordered2 = c(colordered[1], colordered[1:(mylen-1)] )
 colordiff   = (colordered != colordered2)
 colordiff[1]= TRUE

 xColor       = colordered[colordiff]
 xcolStartIdx = myidx[colordiff]
 xcolEndIdx   = c(xcolStartIdx[-1]-1, mylen)

 # excluding grey
 selNoGrey = xColor != "grey" 

 zColor = xColor[selNoGrey]
 zStart = xcolStartIdx[selNoGrey]
 zEnd   = xcolEndIdx[selNoGrey]

 zcolpos= data.frame(rbind(zStart, zEnd ))
 colnames(zcolpos) <- zColor

 zcolpos
}



cutTreeStatic_new = function(hiercluster, heightcutoff=0.99, minsize1=50){
    #for debug: use all genes
    #clusterDemdroHei =demdroHei.order

    demdroHeiAll= rbind( cbind(hiercluster$merge[,1], hiercluster$height), cbind(hiercluster$merge[,2], 
                         hiercluster$height) )

    #singletons will stand at the front of the list
    myorder = order(demdroHeiAll[,1])

    #get # of singletons
    no.singletons = length(hiercluster$order)

    #> demdroHei.sort[1:10,]
    #       [,1]      [,2]
    #[1,] -3000 0.7943783
    #[2,] -2999 0.7863851
    demdroHeiAll.sort = demdroHeiAll[myorder, ]
    demdroHei.sort    = demdroHeiAll.sort[c(1:no.singletons), ]

    #finally, we got tree height for each of the singleton inorder of 1 to no.singletons
    #     [,1]      [,2]
    #[1,]    1 0.8389184
    #[2,]    2 0.8772433
    #[3,]    3 0.8308482
    demdroHei      = demdroHei.sort[seq(no.singletons, 1, by=-1), ]
    demdroHei[,1]  = -demdroHei[,1]

    demdroHei.order = demdroHei[hiercluster$order, ]

    no.cnodes = dim(demdroHei)[1]
            
    # compute height diffrence with mean height
    heidiff       = demdroHei.order[,2] - heightcutoff
    heidiff.shift.s = shiftSequence(heidiff, 1)
    heidiff.shift.e = shiftSequence(heidiff, -1)

    #heidiff2      = demdroHei.order[,2] - shiftSequence(heightcutoff, -1)
    #mean2= abs(heidiff2); sd2 = sd(heidiff2)

    # get cut positions
    # detect the end point of a cluster, whose height should be less than meanhei 
    #  and the node behind it is the start point of the next cluster which has a height above meanhei
    cuts.bool.s = (heidiff<0) & (heidiff.shift.s > 0)
    cuts.bool.e = (heidiff<0) & (heidiff.shift.e > 0)

    # process boundary
    if(demdroHei.order[no.cnodes, 2]<heightcutoff) {
       cuts.bool.e[no.cnodes] = TRUE
    }
    if(demdroHei.order[1, 2]<heightcutoff) {
       cuts.bool.s[1] = TRUE
    }

    #a good candidate cluster-end point should have significant # of ahead nodes with head < meanHei
    mindex.s =c(1:no.cnodes)[cuts.bool.s]
    mindex.e =c(1:no.cnodes)[cuts.bool.e]

    no.cutps = length(mindex.s)-1
    msizes   = mindex.e - mindex.s + 1

    msel = msizes>= minsize1
    no.cutps = sum(msel)

    if(no.cutps >= 2){
       new.clupos = cbind(mindex.s, mindex.e)[msel, ]
    }else{
       #new.clupos = cbind( 1, no.cnodes)
       res=cutTreeStatic_new(hiercluster=hiercluster, heightcutoff=heightcutoff-0.05, minsize1=minsize1)
       return (res)
    }
    #new.clupos

    colorhelp = rep(-1, no.cnodes)
    for (i in c(1:no.cutps)) {
       idx = c(new.clupos[i,1]:new.clupos[i,2]) 
       colorhelp[idx] = i
    }

    #restore the orginal order
    od = order(demdroHei.order[, 1])
    colorhelp2 = colorhelp[od]
    return(colorhelp2)
}


#use height cutoff to remove
cutTreeStatic = function(hiercluster,heightcutoff=0.99, minsize1=50) {

    # here we define modules by using a height cut-off for the branches
    labelpred= cutree2(hiercluster,h=heightcutoff)
    sort1=-sort(-table(labelpred))
    sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 > minsize1
    no.modules=sum(modulebranch)

    colorhelp = rep(-1, length(labelpred) )
    if ( no.modules==0){
        print("No mudule detected\n")
    }
    else{
        for (i in c(1:no.modules)) {
            colorhelp=ifelse(labelpred==modulename[i],i ,colorhelp)
        }
    }
    colorhelp
}




cutree2 = function (tree, k = NULL, h = NULL)
{

    if (is.null(n1 <- nrow(tree$merge)) || n1 < 1)
        stop("invalid 'tree' (merge component)")
    n <- n1 + 1
    if (is.null(k) && is.null(h))
        stop("either 'k' or 'h' must be specified")
    if (is.null(k)) {
        ### tree$height must be sorted
        temp= tree$height-c(0, tree$height[-length(tree$height)])

        sorted=T

        #if(sum(sign(temp[-1])<0)>0) sorted=F

        if (sorted==F)
            stop("the 'height' component of 'tree' is not sorted\n(increasingly); consider applying as.hclust() first")

        k <- integer(length(h))

        k <- n + 1 - apply(outer(c(tree$height, Inf), h, ">"), 2, which.max)

        if (getOption("verbose"))
            cat("cutree(): k(h) = ", k, "\n")

    } else {
        k <- as.integer(k)
        if (min(k) < 1 || max(k) > n)
            stop(gettextf("elements of 'k' must be between 1 and %d",n), domain = NA)
    }

    #ans <- .Call("R_cutree", tree$merge, k, PACKAGE = "stats")
    ans <- cutree(tree, k=k, h=h)

    if (length(k) == 1) {
        ans <- as.vector(ans)
        names(ans) <- tree$labels
    } else {
        colnames(ans) <- if (!is.null(h)) h else k

        rownames(ans) <- tree$labels
    }

    return(ans)

}


#locate the start/end positions of each cluster in the ordered cluster label sequence 
#where "-1" indicating no cluster
#3-1 -1 1 1 1 1 2 2 2
#3 3 -1-1 1 1 1 1 2 2 2  (shift)
#---------------------------------
#0-4  0 2 0 0 0 1 0 0 0   (difference)
#       *     * @
locateCluster = function(clusterlabels)
{
 no.nodes = length(clusterlabels)
 clusterlabels.shift = c(clusterlabels[1], c(clusterlabels[1:(no.nodes-1)]) )
 
 #a non-zero point is the start point of a cluster and it previous point is the end point of the previous cluster
 label.diff = abs(clusterlabels - clusterlabels.shift)

 #process the first and last positions as start/end points if they belong to a cluster instead of no cluster "-1"
 if(clusterlabels[1]       >0) {label.diff[1]=1} 
 if(clusterlabels[no.nodes]>0) {label.diff[no.nodes]=1} 

 flagpoints.bool = label.diff > 0
 flagpoints = c(1:no.nodes)[flagpoints.bool]
 no.points  = length(flagpoints)

 myclupos=NULL
 for(i in c(1:(no.points-1)) ){
   idx = flagpoints[i]
   if(clusterlabels[idx]>0){
      if(flagpoints[i+1]==no.nodes) {#boundary effect
         myclupos = rbind(myclupos, c(idx, flagpoints[i+1]) )
         break
      }else{
         myclupos = rbind(myclupos, c(idx, flagpoints[i+1]-1) )
      }
   }
 }
 myclupos
}


#row-wise reorder a matrix based on "orderedList", default key column in the
#matrix is the first column
orderMergedMatrix = function(disorderMatrix, orderedList, keyCol=1){
  no.samples = length(orderedList)
  cnt = 1
  seqc  =c(1:no.samples)
  rightOrder = rep(0, no.samples)
  orderedMatrix = NULL
  for(each in orderedList){
    whichTrue = ( as.character(each)==as.character(disorderMatrix[, keyCol]) )
    idx = sum(whichTrue * seqc)
    #cat(as.character(each), " : ", as.character(disorderMatrix[idx,keyCol]),as.character(idx), "\n" )
    rightOrder[idx] = cnt
    cnt = cnt + 1
    orderedMatrix = rbind(orderedMatrix, disorderMatrix[idx,])
  }
  colnames(orderedMatrix) <- colnames(disorderMatrix)
  orderedMatrix
}

#-------------------------------------------------------------------------
#Function: convert the upper
#Input:   
#         mvector     ~ adjacency vector with first column is 
#                       the ROW index (i) of the vector in the original matrix A=[a(i,j)]
#         mcutoff     ~ color codes (modules) for genes
#         mgraphfname ~ file for storing output tripples of (i, j, a(i,j)), j>i
#         mfudgefactor~ fudge factor to spread (<1) or compress (>1) a network
#         symmetric   ~ the original adj matrix is symmetric or not
#                        if symmetric, then we just need consider the pairs [i,(i+1):no.nodes],
#                       otherwise, we need consider [i,(i+1):no.nodes]
#
#Output:  tripples of (i, j, a(i,j)), j>i
#-------------------------------------------------------------------------
vectorToPairs=function(mvector, mcutoff=NA, mgraphfname,mfudgefactor=1, symmetric=T, use_signif=F, separator=" ")
{
#mvector = kdatEdge[12, ]
#mcutoff =cutoff
#mgraphfname = graphfname
   i = as.integer(mvector[1])

   #actual data vector
   datv = as.numeric(as.vector(mvector[-c(1)]))
   mno.nodes = length(datv)

   #node index 
   if(symmetric) {
     index <- c((i+1):mno.nodes)   
   }else{
     index <- setdiff(c(1:mno.nodes), i)
   }

   av = datv[index]
   if (is.na(mcutoff)) {
      bv = rep(T, mno.nodes)
   } else{
      bv = av >= mcutoff
   }

   bv = ifelse(is.na(bv), 0, bv)

   no.selected = sum(bv)

   ret=c(0,0,0)
   if(no.selected>0){
     c1=rep(i, no.selected)
     c2=index[bv]
     c3=as.numeric(av[bv])^mfudgefactor
     if(use_signif) {
        c3 = signif(c3, 4)
     }
     ret = as.matrix(cbind(c1,c2,c3))
     write.table(ret, mgraphfname, append=T, sep=separator, quote=F, col.names=F, row.names=F)
   }

   #cat(i," :",as.character(no.selected),  "\n")  
   
   return(ret)
}

adjmatrixToPairs = function(adj, fout, tau=NA, return_pairs=F) {
       nnodes = dim(adj)[1]

       write.table(rbind(colnames(adj)), fout, sep="\t", quote=F, col.names=F, row.names=F)
              
       #redo generating edge frame by considering only nodes with at least one connection   
       jdatEdge <- cbind(c(1:nnodes),adj)
       jdatEdge <- jdatEdge[-nnodes, ]   
       edgeframe = apply(jdatEdge, 1, vectorToPairs, mcutoff= tau, mgraphfname=fout, use_signif=T, separator="\t")

       if(return_pairs){
           allMatrix <- read.delim(fout, sep="\t", header=F, skip=1)
           allMatrix <- as.matrix(allMatrix)
           mynames= colnames(adj)
           pairs = cbind(mynames[as.integer(allMatrix[,1])], mynames[as.integer(allMatrix[,2])])
           return(pairs)
       }
}

adjmatrixToPairsNew = function(adj, symmetric=TRUE, diagonal=FALSE, tau=NA, returnIndexToo=FALSE) {
  nnodes = dim(adj)[1]
  
  idxpair = getMatrixIndex(dim(adj), symmetric=symmetric, diagonal=diagonal)
  if(is.na(tau)) {
    sel = rep(TRUE, nrow(idxpair))
  } else {
    sel = adj[idxpair] >= tau
  }
  rnames = rownames(adj); cnames = colnames(adj)

  idxpair = idxpair[sel, ]    
    
  pairs = cbind(rnames[idxpair[,1]], cnames[idxpair[,2]], adj[idxpair])

  if( returnIndexToo){
    res = as.list(rep(NA,2))
    res[[1]] = pairs
    res[[2]] = idxpair
    return (res)
  } else {
    return(pairs)
  }
}


# get all indices of a matrix
# for symetric matrix size is a number otherwise size=(row, col)
getMatrixIndex = function(size, symmetric=T, diagonal=F)
{
   allidx = NULL

   if(symmetric){
      for(i in c(1:(size[1]-1)) ) {
         iv = cbind(i, (i+1):size[1])
         allidx = rbind(allidx, iv)
      }
      if(diagonal) {allidx = rbind(allidx, cbind(1:size[1], 1:size[1]) )}

   } else {
      for(i in c(1:(size[1])) ) {
         iv = cbind(i, 1:(size[2]))
         allidx = rbind(allidx, iv)
      }
   }

   return (allidx)
}

getMatrixIndexAll = function(ncol, nrow)
{
    idxpair = NULL
    for(i in c(1:ncol) ) {
       iv = cbind(1:nrow, i)
       idxpair = rbind(idxpair, iv)
    }
    return(idxpair)
}


# link type
# Undirected = 0
# Directed = 3
# Inhibitory = 1
# Activation = 2
#
vectorToPairsTGINAV=function(mvector,mgraphfname, mcutoff=0.5, directed=F, linkcolor="grey")
{
#mvector = kdatEdge[12, ]
#mcutoff =cutoff
#mgraphfname = graphfname
   i = as.integer(mvector[1])

   #print(i)

   #actual data vector
   datv = as.numeric(as.vector(mvector[-c(1)]))
   mno.nodes = length(datv)

   if(mno.nodes==i & (!directed) ){
     return ("")
   }

   if (directed){
     dstr = "3"
     index <- c(1:mno.nodes)
   }else{
     #node index 
     dstr = "0"
     index <- c((i+1):mno.nodes) 
   }

   av = datv[index] 
   bv = av >= mcutoff
   no.selected = sum(bv)
   
   ret=c(0,0,0)
   if(no.selected>0){
     c1=rep(i, no.selected)
     c2=index[bv]
     c3=as.numeric(av[bv])
     
     clinkcolor = linkcolor 
     if(length(linkcolor)==1){
        clinkcolor = rep(linkcolor, no.selected)
     }

     #ret = as.matrix(cbind(c1,c2,c3))
     ret = paste('   <link sourceNode="', c1, '"', 
                 ' targetNode="', c2, '"',
                 ' displayValue=""',
                 ' customValue="UNDIRECTED"', 
	         ' linkColor="', clinkcolor, '"',
                 ' linkWeight="0" linkStyle="0" ',
                 ' type="', dstr,'"', 
                 ' sourceNetwork="x"', 
       	         ' originalNetwork="x" originalId="" labelFont="Microsoft Sans Serif, 8pt" />', sep="")
        
     write.table(cbind(ret), mgraphfname, append=T, sep=" ", quote=F, col.names=F, row.names=F)
   }
   #cat(i," :",as.character(no.selected),  "\n")  
   #ret
}


if(F) {
colcode.reduced3= colcode.reduced
mdendro= h1row; genecluster = colcode.reduced;
pccluster= pccolcode;
pcnames = names(pcs[[1]]);
MaxGap = 10;
}


mergeClusterByPCA=function(mdendro, genecluster, pccluster, pcnames, MaxGap=0,
               anameallmodules=T, auseblackwhite=FALSE, useNumberAsLabel=FALSE, startlabel=1){

    geneclusize    = table(genecluster)
    finalgenecolor = as.character(genecluster)
    pcclunames = names(table(pccluster))

    ngenes = length(genecluster)

    # get ordered gene cluster assignment
    orderedCluRaw = as.character(genecluster[mdendro$order])
    fdispcolor = getDisplayColorSequence( orderedCluRaw)
    orderedgeneclus= fdispcolor[fdispcolor!=""]
    
    # get the start and end indices for each module excluding the grey module
    #
    geneModulePos = getDisplayColorPos(orderedCluRaw)

    #we need include grey module, but use different position to differentiate it from the ordinary clusters
    # so that the recombination will not consider the grey module
    orderedgeneclus = c(orderedgeneclus, "grey")

    geneclupos     = c(1:(length(orderedgeneclus)-1), 99999)
    names(geneclupos) <- orderedgeneclus

    # get ordered cluster-sizes
    orderedGeneclusize  = NULL
    for (each in orderedgeneclus){
       orderedGeneclusize  = c(orderedGeneclusize, geneclusize[[each]])
    }
    names(orderedGeneclusize)  <- orderedgeneclus

    for (each in pcclunames ){
       if (each =="grey"){
         next
       }
       
       # PCs belongs to the current cluster 
       sel.pcs      = as.character(pccluster) == each
       if (sum(sel.pcs)<=1){
          next
       }
       sel.pcnames  = pcnames[sel.pcs]
       no.selpcs    = length(sel.pcnames) 

       #get the PCs' sequential positions in the original dendrogram
       epcmGeneIndex = NULL
       for (eachmerged in sel.pcnames){
           if (eachmerged =="grey"){ #don't include the grey genes
              next
           }
           epos = geneModulePos[eachmerged]
           epcmGeneIndex = c(epcmGeneIndex, epos[1,1]:epos[2,1])
       }
       epcmGeneIndex = sort(epcmGeneIndex)
       no.egenes = length(epcmGeneIndex)

       # get the neighboring segments of the PCs, so we only merge the neighboring PCs
       #> order.selPCpos: 30 31 32 36 37 38 39 40
       #> shift.selPCpos: 29 30 31 32 36 37 38 39
       #>  diff.selPCpos:  1  1  1  4  1  1  1  1
       #> allow gaps
       shift.epcmGeneIndex  = c(epcmGeneIndex[1]-100, epcmGeneIndex[c(1:(no.egenes-1))])
       diff.egenes     = epcmGeneIndex - shift.epcmGeneIndex
       bool.startpos   = (diff.egenes > MaxGap + 1)

       startpos        = epcmGeneIndex [bool.startpos]
       endpos          = c( epcmGeneIndex[bool.startpos[-1]], epcmGeneIndex[no.egenes])
       nosegments      = length(endpos)
       
       #we choose the cluster with maximal size as the module assigment for this SEGMENT of clusters
       # here we update only the genes within each segment
       for (iseg in c(1:nosegments) ){
           seg = c(startpos[iseg]:endpos[iseg])
           eiallmodules = names(table(orderedCluRaw [seg]))
           eiallmodulesNogrey =setdiff(eiallmodules, "grey")
           for (x in  eiallmodules){
               orderedCluRaw[seg]=ifelse(orderedCluRaw[seg]==x, eiallmodulesNogrey[1], orderedCluRaw[seg]);
           }

           cat("merged the modules: ",  eiallmodules, "\n")
       }#for (iseg 
      
    }
   
    # now we need sort the module assignment with new names by using an existing function
    # to do so, we need convert the current colors into numerals which are the input of this 
    # existing function
    #retcolors = reassignModuleNames(finalgenecolor)
    if(F){
       a=sample(1:10, 5)
       ao=order(a)
       b=a[ao]#sorted
       c=b[ c(1:5)[order(ao)] ]
    }

    backOrder = c(1:ngenes)[order(mdendro$order)]
    orderedCluRawX = orderedCluRaw[backOrder]
    
    retcolors = reassignModuleNames(orderedCluRawX, minmodulesize=0, anameallmodules=anameallmodules, 
                    auseblackwhite=auseblackwhite, useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)

    retcolors
}

# now we need sort the module assignment with new names by using an existing function
# to do so, we need convert the current colors into numerals which are the input of this 
# existing function
reassignModuleNames=function(mycolorvect, minmodulesize=0, anameallmodules=TRUE, auseblackwhite=FALSE, useNumberAsLabel=FALSE, startlabel=0){
    fgenecolor = as.character(mycolorvect)
    ztable = table(fgenecolor)
    zfinal = rep(0, length(fgenecolor))
    iclu   = 1 
    for (each in names(ztable)){
       if (each=="grey")
          next
       if (ztable[[each]] < minmodulesize)
          next

       zfinal = ifelse(fgenecolor==each, iclu, zfinal)
       iclu=iclu+1
    }

    retcolors=assignModuleColor(labelpred=zfinal, minsize1=1, anameallmodules=anameallmodules, 
                    auseblackwhite=auseblackwhite, useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)
    retcolors
}

# data: columns as genes and rows as samples
# modulemembership: module assignment
# coexppRes: use coexpp object to compute kin
# use_PCA: using Principal component analysis
# fimg: output dendrogram
#
IterativeModuleMergeByPC = function(hierclust, modules, data=NULL, coexppRes=NULL, use_PCA=TRUE, 
              fimg=NULL, useNumberAsLabel=FALSE, startlabel=1, heightcutoff=0.3, useAbsoluteValue=TRUE) 
{

  tbcolor = table(modules)
  pre_nomodules= length(names(tbcolor))

  if(pre_nomodules<=4 ) {
     return (modules)
  }

  modulemembership = modules
  while(pre_nomodules>4){

    tbcolor = table(modulemembership)
    pre_nomodules= length(names(tbcolor))

    # based on PC
    #
    if(use_PCA) {
       pcs = ModulePrinComps(datexpr=as.matrix(data), couleur=as.character(modulemembership) )
    } else {
       #k.all       k.in     k.out     k.diff    k.in.normed  k.all.normed
       #6.128548 3.94784336 2.1807048  1.7671385 0.197392168  0.013294031
       #
       modStat<-intraModularStatistics(coexppRes,modules=as.character(modulemembership),stats=c("connectivity"))

       # take the profile of the most connected gene in each module as PC
       #    
        pcs=ModuleHubProfiles_NoAdjM(data  =as.matrix(data), kin=modStat[,2], 
                couleur=as.character(modulemembership), min_modulesize=10)
    }

    if(useAbsoluteValue) {
      distCorPCs = 1-abs(cor(pcs[[1]],use="p"))
    } else {
      distCorPCs = (1-cor(pcs[[1]],use="p"))/2
    }

    distCorPCs = ifelse(is.na(distCorPCs), 0, distCorPCs)
    h1PCs      = hclust(as.dist(distCorPCs),method="a")

    #displace PC hierarchical dendrogram on screen
    #par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
    #plot(h1PCs, xlab="",ylab="",main="",sub="")
    #par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    #merge a minor cluster to a major cluster, needed in some cases
    #modulemembership = merge2Clusters(modulemembership, mainclusterColor="pink", minorclusterColor="salmon")
    #---- show/save the hierarchical clustering dendrogram and the cluster color bars in the same image
    #plotDendrogramModuleLabels(mdendro=h1row, modulecolors=modulemembership, save2file=NULL, plotLabels=FALSE)
    #plotDendrogramModuleLabels(mdendro=h1row, modulecolors=modulemembership, save2file=imgHierClustModule, plotLabels=FALSE)
    #---- display/save Dendrogram, Cluster Color-bars, and Color Names together
    #plotDendrogramModuleLabels(mdendro=h1row, modulecolors=modulemembership, save2file=NULL, plotLabels=TRUE)
    #plotDendrogramModuleLabels(mdendro=h1row, modulecolors=modulemembership, save2file=imgHierClustModuleColorlabel, plotLabels=TRUE)

    #plot(h1PCs, xlab="",ylab="",main="",sub="")

    # clusters of PCs
    pccolcode = moduleDetectLabel(hiercluster=h1PCs, heightcutoff, minsize1=1, nameallmodules=FALSE, 
                        useblackwhite=FALSE, useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)

    #check whether only "grey" module is detected
    #
    pcmodules   = names(table(pccolcode))
    if( (length(pcmodules)==1) & (pcmodules[1]=="grey") ){
        break
    }

    # modulemembership2 =  modulemembership
    modulemembership = mergeClusterByPCA(mdendro=hierclust, genecluster = modulemembership,
                                        pccluster= pccolcode,
                                        pcnames = names(pcs[[1]]),MaxGap=50,
                                        anameallmodules=FALSE, auseblackwhite=FALSE, 
                                        useNumberAsLabel=useNumberAsLabel, startlabel=startlabel)
    tbcolor = table(modulemembership)
    cur_nomodules= length(names(tbcolor))
    if( pre_nomodules==cur_nomodules  ){
       break
    }
    pre_nomodules= cur_nomodules
  }

  #hiercutinforPC=paste("PC heicutoff=", as.character(heightcutoff))

  if(!is.null(fimg)) {
    #save PC hierarchical dendrogram on screen
    openImgDev(fimg, iwidth = 1200, iheight =800)
    par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
    plot(h1PCs, xlab="",ylab="",main="",sub="")
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off()
  }

  return (modulemembership)
}

moduleClusterByPC_plot = function(modules, data=NULL, coexppRes=NULL, use_PCA=TRUE, fimg=NULL) 
{
    # based on PC
    #
    if(use_PCA) {
       pcs = ModulePrinComps(datexpr=as.matrix(data), couleur=as.character(modules) )
    } else {
       #k.all       k.in     k.out     k.diff    k.in.normed  k.all.normed
       #6.128548 3.94784336 2.1807048  1.7671385 0.197392168  0.013294031
       #
       modStat<-intraModularStatistics(coexppRes,modules=as.character(modules),stats=c("connectivity"))

       # take the profile of the most connected gene in each module as PC
       #    
        pcs=ModuleHubProfiles_NoAdjM(data  =as.matrix(data), kin=modStat[,2], 
                couleur=as.character(modules), min_modulesize=10)
    }

    distCorPCs = 1-abs(cor(pcs[[1]],use="p"))
    distCorPCs = ifelse(is.na(distCorPCs), 0, distCorPCs)
    h1PCs      = hclust(as.dist(distCorPCs),method="a")

    if(!is.null(fimg)) {
      #save PC hierarchical dendrogram on screen
      openImgDev(fimg, iwidth = 1200, iheight =800)
      par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)
      plot(h1PCs, xlab="",ylab="",main="",sub="")
      par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
      dev.off()
    }
}



#orderMergedMatrix = function(disorderMatrix, orderedList, keyCol=1){
#  orderM = order(disorderMatrix[keyCol, ])
#  orderL = order(orderedList)
#  orderL2=order(orderL)  
#  orderM2= orderM[orderL2]
#  disorderMatrix[orderM2, ]
#}

#-------------------------------------------------------------------------
#Function: compute whithin-module cluster coefficient for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module cluster coefficients
#-------------------------------------------------------------------------
computeModuleCC = function(adjMatrix, colorcodeC, weighted=T)
{
   modnames= names( table(colorcodeC) )

   no.genes     = dim(adjMatrix)[1]
   clustercoeff = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   for (each in modnames ){
      if(as.character(each)=="grey" ){
         next
      }
      whichmod = each==colorcodeC
      if( sum(whichmod)>1){
         icc = computeClusterCoefficient(adjMatrix[whichmod,whichmod], weighted)
      }else{
         icc = 1
      }

      #indices of genes in the current module
      idxmod = idxseq * whichmod
      
      #put the cc's to the
      clustercoeff[idxmod] = icc
   }
   clustercoeff
}

#------------------------------------------------------------------------
#Function: compute whithin-module conformity for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module conformity
#-------------------------------------------------------------------------
computeModuleConformity = function(adjMatrix, colorcodeC)
{
   modnames= names( table(colorcodeC) )

   no.genes     = dim(adjMatrix)[1]
   conformity = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   for (each in modnames ){
      if(as.character(each)=="grey" ){
         next
      }
      whichmod = each==colorcodeC

      module.size = sum(whichmod)
      if (module.size==1){
        next
      }

      iconform = SDADJ1(adjMatrix[whichmod,whichmod])     

      #indices of genes in the current module
      idxmod = idxseq * whichmod
      
      #put the conformity to the global list
      conformity[idxmod] = iconform
   }
   conformity
}


#-------------------------------------------------------------------------
#Function: compute total number of connections for each gene excluding 
#          genes in the grey module
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#
#Output:  the number of connections of each gene to the gens 
#          in non-grey modules
#-------------------------------------------------------------------------
computeTotalLinksToNongreygenes = function(adjMatrix, colorcodeC, isAdjacency=TRUE, normalized=FALSE, usegreymodule=F)
{
   modnames= names( table(colorcodeC) )
   no.genes     = dim(adjMatrix)[1]
   links        = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   nongrey      = colorcodeC != "grey"

   total.nonegrey = sum(nongrey)

   totallinks <- apply(adjMatrix[, nongrey], 1,sum, na.rm=TRUE)
   if(!isAdjacency){
      totallinks <- (total.nonegrey - totallinks)
   }

   #normalize against the module size
   if(normalized==TRUE){
      totallinks = totallinks /total.nonegrey
   }      
    
  #put the links's to the buffer      
  #links[nongrey] = totallinks
  
  totallinks
}


#-------------------------------------------------------------------------
#Function: compute whithin-module number of connections for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module number of connections of each gene
#-------------------------------------------------------------------------
computeModuleLinks = function(adjMatrix, colorcodeC, isAdjacency=TRUE, normalized=FALSE, usegreymodule=F)
{
   modnames= names( table(colorcodeC) )

   no.genes     = dim(adjMatrix)[1]
   links        = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   for (each in modnames ){
      if((usegreymodule==F) & (as.character(each)=="grey" ) ){
         next
      }

      whichmod    = each==colorcodeC
      module.size = sum(whichmod)

      if (module.size==1){
        next
      }

      modk <- apply(adjMatrix[whichmod,whichmod],2,sum, na.rm=TRUE) 
      if(!isAdjacency){
          modk <- (module.size -modk)
      }

      #normalize against the module size
      if(normalized==TRUE){
         modk = modk/module.size
      }

      #indices of genes in the current module
      idxmod = idxseq * whichmod
      
      #put the links's to the buffer      
      links[idxmod] = modk
   }
   links
}


#---------------------------------------------------------------------------------
#Function: compute whithin-module number of connections for each row-based element
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  a matrix of cross-module numbers of connections of each of row-based elements
#         columns as modules and rows as elements of interest
#
#Notice that here it is different from previous function, here adjMatrix may be
#     asymmetric, columns are color-coded
#
#----------------------------------------------------------------------------------
computeMultiModuleConnectivity = function(adjMatrix, colorcodeC)
{
   modnames   = names( table(colorcodeC) )
   kin.matrix = NULL
   for (each in modnames ){
      whichmod    = each==colorcodeC
      if( sum(whichmod) >1){
         modk <- apply(adjMatrix[, whichmod], 1, sum, na.rm=TRUE)
      }else{
         modk <- as.integer(adjMatrix[, whichmod])
      }
      kin.matrix = cbind(kin.matrix, modk)
   }
   colnames(kin.matrix) <- modnames
   kin.matrix
}

#process markers location into a sequence in order of chromsomes 1 to 21
processChromPosSequence = function(avector, colorcodeC)
{
  #initilaization
  no.elements= length(avector)
  outpos = avector
  idxseq =c(1:no.elements)

   modnames= names( table(colorcodeC) )
   currMax=0
   for (each in modnames ){      
      whichmod = each==colorcodeC
      #indices of elements in the current module
      idxmod = idxseq * whichmod
      outpos[idxmod] = outpos[idxmod] + currMax 
      
      #last marker pos in the current Xmosome
      imax   = max( avector[whichmod] )
      currMax=  currMax + imax + 1
   }
   
   #tried to find the Xmosoome boundaries
   shifted = c(0, colorcodeC[1:(no.elements-2)], 0)
   boundary= (shifted != colorcodeC)
   boundary.colors=ifelse(boundary==TRUE, colorcodeC, "")
   boundary.colors[ length(boundary.colors) ] = "."

   list(outpos, boundary.colors)
}

is.constant=function(vect, threshold=10^-50){
  diff = vect-vect[1]
  dsum= sum(abs(diff)<threshold, na.rm=TRUE)
  if(is.na(dsum)) {return (TRUE)}
  if(dsum==length(vect)){ return (TRUE)}
  return (FALSE)
}

#perform ttest for two groups (expression values across samples)
#two group are represented by two boolean vectors with the same size as exprvector
#
# if usePseudoSD=T, then, for the group with 1 sample, we use an addiitonal pseudo sample
#  to ensure pvalue is computable, in this way we can know that the standard deviation of 
#   the group with more than one sample
# 
#
#ttestVector (exprvector=c(1:10), tgroupA=c(T, rep(F,9)), tgroupB=c(rep(F,8), T,T), usePseudoSD=F)
ttestVector = function(exprvector, tgroupA, tgroupB, usePseudoSD=F){

    #all the elements of one group are NA
    if( (sum(is.na(exprvector[tgroupA]))==sum(tgroupA)) ||  (sum(is.na(exprvector[tgroupB]))==sum(tgroupB)) )
      return (1)

    #too few elements in one group to be used to perform ttest
    if (!usePseudoSD) {
       if( (sum(!is.na(exprvector[tgroupA]))<2) ||  (sum(!is.na(exprvector[tgroupB]))<2) )
           return (1)
    }

    # identical
    diff= abs(exprvector-exprvector[1])
    dsum= sum(diff, na.rm=TRUE)
    if(is.na(dsum)) {return (1)}
    if(dsum==0){ return (1)}

    myexprA=exprvector[tgroupA]
    myexprB=exprvector[tgroupB]

    # the two vectors are constant and have the same value
    constant = is.constant(c(myexprA, myexprB))
    if(constant){ return (1)}

    # the two vectors are constant but have different values
    constantA = is.constant(myexprA)
    constantB = is.constant(myexprB)
    if(constantA & constantB){ return (0)}

    delta=5
    #too small group-size
    if ( (sum(tgroupA)<2)  || (sum(tgroupB)<2) ){

       if(!usePseudoSD) #no need for further computation of pvalue
           return (1)

       if ( sum(tgroupA)==1)
          myexprA =c(myexprA-delta, myexprA+delta)
       
       if ( sum(tgroupB)==1)
          myexprB =c(myexprB-delta, myexprB+delta)
    }

    #print(c(myexprA,myexprB, sum(diff)))
    
    tt    = t.test(myexprA,myexprB)
    tt$p.value

    #print(tt$p.value)
}

ttest2Vectors = function(exprvectorA, exprvectorB, usePseudoSD=F){

    #all the elements of one group are NA
    if( (sum(!is.na(exprvectorA))==0) ||  (sum(!is.na(exprvectorB))==0) )
      return (1)

    #too few elements in one group to be used to perform ttest
    if (!usePseudoSD) {
       if( (sum(!is.na(exprvectorA))<2) ||  (sum(!is.na(exprvectorB))<2) )
           return (1)
    }

    # identical
    diff= abs(c(exprvectorA,exprvectorB)-exprvectorA[1])
    dsum= sum(diff, na.rm=TRUE)
    if(is.na(dsum)) {return (1)}
    if(dsum==0){ return (1)}
    if(dsum<10^-50){ return (1)}

    delta=5
    #too small group-size
    if ( (length(exprvectorA)<2)  || (length(exprvectorB)<2) ){

       if(!usePseudoSD) #no need for further computation of pvalue
           return (1)

       if ( length(exprvectorA)==1)
          exprvectorA=c(exprvectorA-delta, exprvectorA+delta)
       
       if (length(exprvectorB)==1)
          exprvectorB=c(exprvectorB-delta, exprvectorB+delta)
    }

    #print(c(myexprA,myexprB, sum(diff)))
    
    tt    = t.test(exprvectorA,exprvectorB)
    tt$p.value

    #print(tt$p.value)
}

ttest_VectIndexed_Matrix = function(exprvectorA, exprMtrxB, usePseudoSD=F){

    res1=ttest2Vectors(exprvectorA[-1], exprMtrxB[exprvectorA[1],], usePseudoSD=usePseudoSD)
    return(res1)
}

ttest2Matrixes = function(exprMtrxA, exprMtrxB, usePseudoSD=F){

    iexprMtrxA = cbind(c(1:nrow(exprMtrxA)), exprMtrxA)
    
    res=apply(iexprMtrxA, 1, ttest_VectIndexed_Matrix, exprMtrxB, usePseudoSD=usePseudoSD)

    return (res)
}

compute_robust_residuals1<-function(m,er){
   lm.obj<-tryCatch(rlm(m~as.factor(er),na.action="na.exclude",init="lts", psi=psi.huber),
                     lm(m~as.factor(er)) ) 
   resid.data<-residuals(lm.obj)

}
#allMatrix.adjust <-apply(allMatrix,1,compute.robust.residuals1,er=xer)

compute_robust_residuals5<-function(m,fr1, fr2, fr3, nr1, nr2){
   lm.obj <- tryCatch(rlm(m~as.factor(fr1)+as.factor(fr2)+ as.factor(fr3)+ nr1+nr2,na.action="na.exclude",init="lts", psi=psi.huber),                     lm(m~as.factor(fr1)+as.factor(fr2)+ as.factor(fr3)+ nr1+nr2) ) 

   resid.data<-residuals(lm.obj)
}

compute_residuals5<-function(m,fr1, fr2, fr3, nr1, nr2, KeepMin=FALSE){
   lm.obj <- lm(m~as.factor(fr1)+as.factor(fr2)+ as.factor(fr3)+ nr1+nr2)
   resid.data <- residuals(lm.obj)
   if(KeepMin) {
     resid.data2 = resid.data - min(resid.data)  + min(m)
     return (resid.data2)   
   } else {
     return (resid.data)
   }
}


compute_residuals<-function(m,er){
   lm.obj<-lm(m ~ ., data=er, na.action="na.exclude")
   resid.data<-residuals(lm.obj)
}


coxSurvModel = function(exprvector, survector, minvar=0.000001){
    exprvector.num = as.numeric(exprvector)
    mpvalue=1
    #for PM/MM truncated cases, a few changes across samples, cox model will crash
    levels = names( table(exprvector.num) )

    if( (var(exprvector.num) >minvar) & (length(levels) >2) ) {
      cox1=coxph(survector ~ exprvector.num, na.action="na.omit")
      mpvalue=1-pchisq(cox1$score, 1)
    }
    mpvalue
}

coxSurvModelComplex = function(exprvector, survector, minvar=0.000001, outimage=NULL, unit="days", imgsize=400,  mtitle="", uselogrank=F, showPval=TRUE){
    exprvector.num = as.numeric(exprvector)
    mpvalue = 1
    sde     =-1
    lower.95=-1
    upper.95=-1
    hr      =-1
    kmpvalue= 1

    #print(exprvector)

    #for PM/MM truncated cases, a few changes across samples, cox model will crash
    levels = names( table(exprvector.num) )

    if( (var(exprvector.num, na.rm =T) >=minvar) & (length(levels) >2) ) {

       cox1=coxph(survector ~ exprvector.num, na.action="na.omit")
       mpvalue=1-pchisq(cox1$score, 1) #logrank test pvalue

       groups=I(exprvector.num > median(exprvector.num, na.rm =T) )

       if( length(unique(groups)) >1 ) {

       cox2 = coxph(survector ~ groups)
       sde= sqrt(cox2$var)

       lower.95=exp(cox2$coef-1.96*sde)
       upper.95=exp(cox2$coef+1.96*sde)
       hr      =exp(cox2$coef)

       # for KM plot 
       groupsnum = as.numeric(groups+0)
       fit <- survfit(survector ~ groupsnum)
       diff = survdiff(survector ~ groupsnum)
       kmpvalue= 1-pchisq(diff$chisq,1)

       if(!is.null(outimage) ) {

           if(uselogrank) {
             mtitlep = paste("P<", signif(mpvalue, 2), sep="")
           } else{
             mtitlep = paste("P<", signif(kmpvalue, 2), sep="")
           }

           openImgDev(imgname=outimage, iwidth = imgsize, iheight = imgsize, ipointsize = 15)
           par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
           plot(fit, xlab=paste("Survival Time (", unit, ")",sep=""), 
                     ylab="Survival Probability", col=c("blue", "red"))#, main=mtitle)
           title(mtitle, cex.main=1)
           #tx=0.75*max(fit[[2]])+ min(fit[[2]])*0.25; ty=0.95
           tx=0.6*max(fit[[2]])+ min(fit[[2]])*0.25; ty=0.1
           
           if(showPval) {
              text(tx, ty, labels=mtitlep, cex=1.2)
           }

           par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
           dev.off()
       }

      } # unique
    }
    
    outphr = c(mpvalue, hr, lower.95, upper.95, sde, kmpvalue)
}

coxSurvModelComplex_LeaveOneOut = function(exprvector, survector, minvar=0.000001, outimage=NULL, unit="days", imgsize=400,  mtitle=""){
    exprvector.num = as.numeric(exprvector)
    mpvalue = 1
    sde     =-1
    lower.95=-1
    upper.95=-1
    hr      =-1
    kmpvalue= 1

    #print(exprvector)

    #for PM/MM truncated cases, a few changes across samples, cox model will crash
    levels = names( table(exprvector.num) )

    if( (var(exprvector.num, na.rm =T) >=minvar) & (length(levels) >2) ) {

       cox1=coxph(survector ~ exprvector.num, na.action="na.omit")
       mpvalue=1-pchisq(cox1$score, 1)

       #ihr = predict(cox1, exprvector.num[1])

       groups=I(exprvector.num > median(exprvector.num, na.rm =T) )

       if(FALSE){ # search the best cutoff
           minp = 1; mincut = 0; allres=NULL
           mrank= rank(exprvector.num)
           mpers= c(2:8)/10
           mcuts= length(mrank)*mpers
           for(kkk in c(1:length(mcuts)) ) {
               mcut=mcuts[kkk]
               mper=mpers[kkk]
               groups=I(mrank>mcut)
               cox2 = coxph(survector ~ groups)
               sde= sqrt(cox2$var)

               lower.95=exp(cox2$coef-1.96*sde)
               upper.95=exp(cox2$coef+1.96*sde)
               hr      =exp(cox2$coef)

               # for KM plot 
               groupsnum = as.numeric(groups+0)
               fit <- survfit(survector ~ groupsnum)
               diff = survdiff(survector ~ groupsnum)
               kmpvalue= 1-pchisq(diff$chisq,1)
               if(kmpvalue<minp){minp=kmpvalue; mincut=mcut}
               allres= rbind(allres, c(mper, mcut, kmpvalue))
           }
           colnames(allres) = c("TopExpr", "TopNumber", "logrankP")

       }

       cox2 = coxph(survector ~ groups)
       sde= sqrt(cox2$var)

       lower.95=exp(cox2$coef-1.96*sde)
       upper.95=exp(cox2$coef+1.96*sde)
       hr      =exp(cox2$coef)

       # for KM plot 
       groupsnum = as.numeric(groups+0)
       fit <- survfit(survector ~ groupsnum)
       diff = survdiff(survector ~ groupsnum)
       kmpvalue= 1-pchisq(diff$chisq,1)

       if(!is.null(outimage) ) {

           mtitlep = paste("P=", signif(kmpvalue, 2), sep="")

           openImgDev(imgname=outimage, iwidth = imgsize, iheight = imgsize, ipointsize = 11)
           par(mfrow=c(1,1), mar=c(5, 4, 2, 2) + 0.1)
           plot(fit, xlab=paste("Survival Time (", unit, ")",sep=""), 
                     ylab="Survival Probability", col=c("blue", "red"), main=mtitle)
           tx=0.75*max(fit[[2]])+ min(fit[[2]])*0.25; ty=0.25
           text(tx, ty, labels=mtitlep, cex = 2.5)
           par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
           dev.off()
       }
    }

    outphr = c(mpvalue, hr, lower.95, upper.95, sde, kmpvalue)
}

returnNULL = function(){return (NULL)}

# colnames of exprMatrix have to be "x1, x2, x3, x4"
coxSurvModelComplex_multiVariables = function(exprMatrix, survector, minvar=1, outimage=NULL, unit="days", imgsize=400, mtitle="", uselogrank=F, showPval=TRUE){
    exprvector.num = as.numeric(as.matrix(exprMatrix))
    mpvalue = 1
    sde     =-1
    lower.95=-1
    upper.95=-1
    hr      =-1
    kmpvalue= 1

    # no missing value is allowed in any case
    naMtrx = is.na(exprMatrix)
    nas = apply(naMtrx, 1, sum)
    selNoNAs = nas <1

    #print(exprvector)

    #for PM/MM truncated cases, a few changes across samples, cox model will crash
    levels = names( table(exprvector.num))

    if( (var(exprvector.num, na.rm =T) >=minvar) & (length(levels) >2) ) {

       #cox1=coxph(survector[selNoNAs,] ~ ., data=exprMatrix[selNoNAs, ], na.action="na.omit")
       cox1=tryCatch(coxph(survector[selNoNAs,] ~ ., data=exprMatrix[selNoNAs, ], na.action="na.omit"),
                      error=function(e) returnNULL() )

       if(is.null(cox1)){return (c(1, 0, 0, 0, 0, 1) )} # a singular fit. 

       mpvalue=1-pchisq(cox1$score, dim(exprMatrix)[2])

       groups=I(cox1$linear.predictors > median(cox1$linear.predictors, na.rm =T) )

       cox2 = coxph(survector[selNoNAs,] ~ groups)
       sde= sqrt(cox2$var)

       lower.95=exp(cox2$coef-1.96*sde)
       upper.95=exp(cox2$coef+1.96*sde)
       hr      =exp(cox2$coef)

       # for KM plot
       groupsnum = as.numeric(groups+0)
       fit <- survfit(survector[selNoNAs,] ~ groupsnum)
       diff = survdiff(survector[selNoNAs,] ~ groupsnum)
       kmpvalue= 1-pchisq(diff$chisq,1)

       if(!is.null(outimage) ) {
           if(uselogrank) {
             mtitlep = paste("P<", signif(mpvalue, 2), sep="")
           } else{
             mtitlep = paste("P<", signif(kmpvalue, 2), sep="")
           }

           openImgDev(imgname=outimage, iwidth = imgsize, iheight = imgsize, ipointsize = 15)
           par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
           plot(fit, xlab=paste("Survival Time (", unit, ")",sep=""), 
                     ylab="Survival Probability", col=c("blue", "red"), main=mtitle)
           #tx=0.75*max(fit[[2]])+ min(fit[[2]])*0.25; ty=0.95
           tx=0.6*max(fit[[2]])+ min(fit[[2]])*0.25; ty=0.1

           if(showPval) {
             text(tx, ty, labels=mtitlep, cex=1.2)
           }

           par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
           dev.off()
       }

    }

    outphr = c(mpvalue, hr, lower.95, upper.95, sde, kmpvalue)
}


# colnames of exprMatrix have to be "x1, x2, x3, x4"
corBYlm = function(trait, exprMatrix){

  if(F){
    a= sample(1:100,20, replace=T)
    b= sample(1:100,20, replace=T)
    #b[1]=NA
    lm1 = lm(a ~b)
  }

  if(sum(!is.na(trait))<2){return (c(0,1)) }

  lm1=lm(trait ~ ., data=exprMatrix, na.action="na.omit")
  slm1 = summary(lm1)
  pvalue= 1-pf(q=slm1$fstatistic[1], df1=slm1$fstatistic[2], 
                 df2=slm1$fstatistic[3], lower.tail = TRUE, log.p = FALSE)

  res = c( (slm1$r.squared)^0.5, pvalue)

  res[1]=ifelse(is.na(res[1]), 0, res[1])
  res[2]=ifelse(is.na(res[2]), 1, res[2])
  
  return(res)
}

# colnames of exprMatrix have to be "x1, x2, x3, x4"
logit = function(trait, exprMatrix){

  if(F){
    a= sample(1:100,20, replace=T)
    b= sample(1:100,20, replace=T)
    #b[1]=NA
    lm1 = lm(a ~b)
  }

  if(sum(!is.na(trait))<2){return (c(0,1)) }

  #lm1=lm(trait ~ ., data=exprMatrix, na.action="na.omit")
  #slm1 = summary(lm1)
  #pvalue= 1-pf(q=slm1$fstatistic[1], df1=slm1$fstatistic[2], 
  #               df2=slm1$fstatistic[3], lower.tail = TRUE, log.p = FALSE)
  #res = c( (slm1$r.squared)^0.5, pvalue)

  mylogit<- glm(trait~ ., data=exprMatrix, family=binomial(link="logit"), na.action=na.pass)

  smylogit= summary(mylogit)
  coeff_pval = smylogit$coefficients[,4]

  # You can also exponentiate the coefficients and interpret them as odds-ratios. 
  #exp(mylogit$coefficients)

  pvalue = dchisq(mylogit$null.deviance-mylogit$deviance, mylogit$df.null-mylogit$df.residual)
  pvalue = ifelse(pvalue >1, 1, pvalue)
  res   = list(mylogit$aic, pvalue, exp(mylogit$coefficients), coeff_pval)

  #res[[1]]=ifelse(is.na(res[[1]]), 10000, res[[1]])
  #res[[2]]=ifelse(is.na(res[[2]]), 1, res[[2]])
  names(res) <- c("AIC", "PV", "OddsRatio", "CoeffPV")
  return(res)
}


computeTOM= function(adjMatrix) {
   no.singletons <- dim(adjMatrix)[1]
   nolinks.reduced  <- apply(adjMatrix, 2, sum, na.rm=TRUE)
   #Let?s calculate topological overlap matrix
   numTOM= adjMatrix %*% adjMatrix + adjMatrix
   dist1 = matrix(NA, no.singletons, no.singletons)
   diag(dist1) <- 0
   for (i in 1:(no.singletons-1) ){
      for (j in (i+1):no.singletons){
          denomTOMij = min(nolinks.reduced[i], nolinks.reduced[j]) + 1 - adjMatrix[i,j]          
          dist1[i,j] = 1- numTOM[i,j]/denomTOMij
          dist1[j,i] = dist1[i,j]
      }
   }
  dist1
}

computeLinksInNeighbors <- function(x, imatrix)
{
  y= x %*% imatrix %*% x
  y
}

computeClusterCoefficient = function(adjMatrix, weighted=F) {

        no.genes <- dim(adjMatrix)[1]
        nolinksNeighbors <- c(rep(-666,no.genes))
        total.edge <- c(rep(-666,no.genes))

        #for (i in 1:no.genes){
        #     nolinksNeighbors[i] <-  adjMatrix[i,] %*% adjMatrix %*% adjMatrix[,i]
        #     #total.edge[i] <-  adjMatrix[i,] %*% Pmax %*% adjMatrix[,i]
        #}
        nolinksNeighbors <- apply(adjMatrix, 1, computeLinksInNeighbors, imatrix=adjMatrix)

        plainsum  <- apply(adjMatrix, 1, sum)
        if(weighted) {
           squaresum <- apply(adjMatrix^2, 1, sum)
           total.edge = plainsum^2 - squaresum
        }else{ # for unweighted network, 1^2 = 1
           total.edge = plainsum^2 - plainsum
        }

        # in case of single node, this will not affect the CC computation
        #
        total.edge = ifelse(total.edge==0, 1, total.edge)

        cluster.coef = nolinksNeighbors/total.edge
        cluster.coef = ifelse(total.edge>0,cluster.coef, 0) 

        cluster.coef
}

pajekColorcode2 = function(bincolorcode)
{
   colorcodeC=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow","grey")
   colorcodeN=c("28",        "4",   "14",   "1",     "21",   "3",  "13",   "5",   "17",     "8",     "24",         "30", "22",   "0",    "18",           "31",       "33",     "15",         "16",         "38" )

   rcolors=bincolorcode
   clevels = length(colorcodeC)
   for (i in c(1:clevels) ){
      whichcolor = rcolors==colorcodeC[i]
      rcolors = ifelse(whichcolor, colorcodeN[i], rcolors)
   }
   rcolors
}

pajekColorcode = function(bincolorcode)
{
   colorcodeC=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow","grey")
   colorcodeN=c("28",        "4",   "14",   "1",     "21",   "3",  "13",   "5",   "17",     "8",     "24",         "30", "22",   "0",    "18",           "31",       "33",     "15",         "16",         "38" )

   colorcodeCN= cbind(colorcodeC, colorcodeN)
   merged = mergeTwoMatricesByKeepAllPrimary(cbind(bincolorcode), colorcodeCN, 
                     missinglabel=NA, keepAllPrimary=T, keepPrimaryOrder=T)
   merged = as.matrix(merged)

   rcolors = merged[,2]
   rcolors = ifelse( is.na(rcolors), "38", rcolors)
   rcolors
}



#------- Pajakificator -----
#
Pajakify_adj = function (nodenames, nodemodules, adjmatrix, ofname="tmp", cutoff=0.5, remove_singular=F) {

    graphfname    = paste(fname, "_paj.net", sep='')
    clustfname    = paste(fname, "_paj.clu", sep='')
    pairfname     = paste(fname, ".pair",    sep='')
    tmpfname      = "tmp.pjk"

    no.nodes <- length(nodenames)

    #-----------------------------------------------------------------------------------
    #first column is the ROW index (i) of the vector in the original matrix A=[a(i,j)]
    #AND remove the last row, as all its links have been processed by the previous nodes
    #
    idatEdge <- cbind(c(1:no.nodes), as.matrix(adjmatrix))
    idatEdge <- idatEdge[-no.nodes, ]

    #---------------------- START from here when try to use new cutoff -----------------
    #---------------------- START from here when try to use new cutoff -----------------
    #
    close(file(graphfname, "w+") )
    close(file(clustfname, "w+") )
    close(file(pairfname,  "w+") )

    appendStringToFile(pairfname, "nodei nodej strength module")

    #cluster assignment
    #clusters = as.character(nodemodules)
    clusters  = pajekColorcode(as.character(nodemodules))
    
    if (remove_singular==FALSE){

       #------------------------------- Cluster File(*.clu)------------------------------------

       clu.keywords = paste("*vertices ", as.character(no.nodes), sep="")
       #write out cluster file
       finalclusters <- c(clu.keywords, clusters)
       write.table(finalclusters, clustfname, sep="\t", quote=FALSE, col.names=F, row.names=F)

       #------------------------------- Network File (*.net) ------------------------------------
       #write out edge files: Part 1: node names
       #*Vertices  36
       #  1 "ASPM"
       #  2 "BUB1B"

       vtitle=paste("*vertices ", as.character(no.nodes),sep="")
       appendStringToFile(graphfname, vtitle)

       quotenodenames = paste("\"", nodenames, "\"", sep="")
       nodelabelMatrix = cbind(as.character(c(1:no.nodes)), quotenodenames)
       write.table(nodelabelMatrix, graphfname, append=T,sep=" ", quote=F, col.names=F, row.names=F)

       appendStringToFile(graphfname, "*Edges ")

       edgeframe = apply(idatEdge , 1, vectorToPairs, mcutoff=cutoff, mgraphfname=graphfname)

    }else{

       #save the output to a temporary file, not the result file
       edgeframe = apply(adjmatrix, 1, vectorToPairs, mcutoff=cutoff, mgraphfname=tmpfname)

       #read in edge files in format of "nodei, nodej, strength"   
       edgeMatrix <- read.delim(tmpfname,sep=" ", header=T)
       
       #count number of links
       nodesInNet = c(as.integer(edgeMatrix[,1]),  as.integer(edgeMatrix[,2]))
       mhist = hist(nodesInNet, breaks=c(0:no.nodes) )

       #get index of nodes with connectivity >0
       nodeIndex =   c(1:no.nodes)[mhist$counts>0]

       no.nodes.sel = length(nodeIndex)

       #------------------------------- Cluster File(*.clu)------------------------------------
       #cluster assignment

       clu.keywords = paste("*vertices ", as.character(no.nodes.sel), sep="")
       #write out cluster file
       finalclusters <- c(clu.keywords, clusters[nodeIndex])
       write.table(finalclusters, clustfname, sep="\t", quote=FALSE, col.names=F, row.names=F)

       #------------------------------- Network File (*.net) ------------------------------------
       #write out edge files: Part 1: node names
       #*Vertices  36
       #  1 "ASPM"
       #  2 "BUB1B"

       #*Vertices  36
       vtitle=paste("*vertices ", as.character(no.nodes.sel),sep="")
       appendStringToFile(graphfname, vtitle)

       #  1 "ASPM"
       #  2 "BUB1B"
       quotenodenames = paste("\"", nodenames[nodeIndex], "\"", sep="")
       nodelabelMatrix = cbind(as.character(c(1:no.nodes.sel)), quotenodenames)
       write.table(nodelabelMatrix, graphfname, append=T, sep=" ", quote=F, col.names=F, row.names=F)

       appendStringToFile(graphfname, "*Edges ")

       #redo generating edge frame by considering only nodes with at least one connection   
       jdatEdge <- cbind(c(1:no.nodes.sel), as.matrix(adjmatrix[nodeIndex, nodeIndex]))
       jdatEdge <- jdatEdge[-no.nodes.sel, ]   
       edgeframe = apply(jdatEdge, 1, vectorToPairs, mcutoff=cutoff, mgraphfname=graphfname)   
    }

}


#<?xml version="1.0" standalone="yes" ?>
#<nav>
#<network name="CONVERTED-6-n-3-e-exampleInputFile" bgColor="White" thresholdType="0" thresholdValue="0" aggregateNodes="false" aggregateLinks="false" labelsInside="false" resolve="false"/>
#<nodes>
#  <node id="1" displayValue="18" fillColor="grey" borderColor="red" borderWidth="4" labelColor="Black"
#	labelFont="Microsoft Sans Serif, 8pt" model="" label="ASD" shape="Circle"
#	pinned="false" species="mouse" symbol="ASD" type="0"
#	centerX="0" centerY="0" width="1" height="1" scaleFactorValue="2"/>
#</nodes>
#<links>
#  <link sourceNode="1" targetNode="2" displayValue="" customValue="DIRECTED" 
#	linkColor="blue" linkWeight="0" linkStyle="0" type="3" sourceNetwork="x"
#	originalNetwork="x" originalId="" labelFont="Microsoft Sans Serif, 8pt" />
#</links>
#<notes />
#</nav>
#
TGINAV_adj = function (nodenames, adjmatrix, directed=F, 
                       nodemodules="gray", nodeShapes="circle", nlabelColor="black",
                       linkcolor="grey", 
                       ofname="tmp", cutoff=0.5, remove_singular=F,
                       nspecies    ="human", notes="", fontsize=12) {


    close(file(ofname, "w+") )

    xmlHeader = '<?xml version="1.0" standalone="yes" ?>'
    titlestr  = paste('<network name="', ofname, '"',
                      ' bgColor="White" thresholdType="0"',
                      ' thresholdValue="', cutoff,'"',
                      ' aggregateNodes="false" aggregateLinks="false" labelsInside="false"',
                      ' resolve="false"/>', sep="")

    no.nodes <- length(nodenames)

    # 1.--------- header ------------------
    appendStringToFile(ofname, xmlHeader)
    appendStringToFile(ofname, "<nav>")
    appendStringToFile(ofname, titlestr)

    # 2.--------- nodes -------------------
    #
    if(length(nodemodules)==1){
      nboardcolor = rep(nodemodules, no.nodes)
      ncolor      = rep(nodemodules, no.nodes)
    } else{
      nboardcolor = nodemodules
      ncolor      = nodemodules
    }

    nboardwidth = rep(4,      no.nodes)
    #nlabelColor = rep("Black",no.nodes) 

    appendStringToFile(ofname, "<nodes>")
    #for (k in c(1:no.nodes)) {
      kstr = paste('   <node id="', c(1:no.nodes), '"',
                   ' displayValue="18"',
                   ' fillColor="',   ncolor, '"',
                   ' borderColor="', nboardcolor, '"',
                   ' borderWidth="4"',
                   ' labelColor="', nlabelColor, '"',
	           #' labelFont="Microsoft Sans Serif, 8pt" model=""',
	           ' labelFont="Microsoft Sans Serif, ', fontsize, 'pt" model=""',
                   ' label="', nodenames, '"',
                   ' shape="', nodeShapes, '"',
	           ' pinned="false"', 
                   ' species="', nspecies, '"',
                   ' symbol="', nodenames, '"',
                   ' type="0"',
	           ' centerX="0" centerY="0" width="1" height="1" scaleFactorValue="2"/>', sep="")
      #appendStringToFile(ofname, kstr)

      write.table(cbind(kstr), ofname, quote=FALSE,sep='\t',col.names=F,row.names=FALSE,append=T) 

    #}

    appendStringToFile(ofname, "</nodes>")

    # 3.--------- links -------------------
    #
    jdatEdge <- cbind(c(1:no.nodes), as.matrix(adjmatrix))

    appendStringToFile(ofname, "<links>")
    edgeframe = apply(jdatEdge , 1, vectorToPairsTGINAV, mcutoff=cutoff, 
                                    mgraphfname=ofname, directed=directed, linkcolor=linkcolor)
    appendStringToFile(ofname, "</links>")
    mynotes = paste('<notes component="', notes, '"/>', sep="")
    appendStringToFile(ofname,  mynotes)   
    appendStringToFile(ofname, "</nav>")
}


TGINAV_pair = function(netpairs, directed=F, filename="tmp", species="human", nodecolor=NULL, nodeshape=NULL) {

    nodenames       = union(netpairs[,1], netpairs[,2])    
    uniquenames     = sort(nodenames)
    no.uniquenames  = length(uniquenames)

    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )

    if(is.null(nodecolor)){
       mnodecolor="gray"
    } else{
       merged = merge(uniquenames, nodecolor, by.x=1, by.y=1, all.x=T)
       merged = as.matrix(merged)
       
       mnodecolor = merged[,2]
       mnodecolor = ifelse(is.na(mnodecolor), "gray", mnodecolor)
    }

    if(is.null(nodeshape)){
       mnodeshape="circle"
    } else{
       merged = merge(uniquenames, nodeshape, by.x=1, by.y=1, all.x=T)
       merged = as.matrix(merged)
       
       mnodeshape = merged[,2]
       mnodeshape = ifelse(is.na(mnodeshape), "circle", mnodeshape)
    }


    # 3. adjacency matrix
    adjmatrix = makeAjacencyMatrix(inetmatrix=netpairs, coding=c(0,1),
                                     matrixsize=no.uniquenames, directed=directed, 
                                     myname2idxMatrix = name2idxMatrix)

    TGINAV_adj(uniquenames, adjmatrix, directed=directed, 
               nodemodules=mnodecolor, nodeShapes=mnodeshape, nlabelColor="black",
               linkcolor="gray", 
               ofname=filename, cutoff=0.5, remove_singular=F,
               nspecies    =species, notes="")

}

BN2linkpairs = function(fnameBN, seperator="->") {

    allMatrixA <- read.delim(fnameBN,sep="\t", header=F)
    allMatrixA <- as.matrix(allMatrixA)
    dim(allMatrixA)

    no.links =  dim(allMatrixA)[1]
    links = NULL
    for(j in c(1:no.links)) {
      jsplitted = splitString(allMatrixA[j,], seperator)
      if(length(jsplitted)==1){
        next
      }

      jsplitted2= splitString(jsplitted[2], " ")
      
      links = rbind(links, c(jsplitted[1], jsplitted2[1]) )
    }
    return (links)
}


statistics_edgelinks = function(netmatrix, filename, keyword="", usedirect=F)
{
  no.links = dim(netmatrix)[1]
  mnetmatrix = as.matrix(netmatrix)
  if (usedirect){
    nodenames = c(as.character(mnetmatrix[,1]), as.character(mnetmatrix[,2]))
  }else{
    #nodenames = as.character(mnetmatrix[,1])
    nodenames = c(as.character(mnetmatrix[,1]), as.character(mnetmatrix[,2]))
  }

  # number of links for each unique gene
  ntab = table(nodenames)
  
  #no of unique genes
  no.uniquegenes = length(ntab)

  # average links per gene
  avglinks=  no.links/no.uniquegenes
 
  # output table to a file
  tabfname = paste(filename, "_nolinks.xls",sep="")
  ordertab = order(-ntab )
  write.table(ntab[ordertab], tabfname, sep="\t",quote=FALSE, col.names=F, row.names=T)
  
  # output unique nodeIDs
  idfname = paste(filename, "_nodeIDs.txt",sep="")
  write.table(as.matrix(names(ntab)), idfname, sep="\t",quote=FALSE, col.names=F, row.names=F)

  # draw histogram, compute frequency of no of links
  ntabtab = table(as.numeric(ntab))

  linkfreq = as.numeric(ntabtab)
  linkidx  = as.numeric( names(ntabtab) )

  maxlinks = max(linkidx)
  mylabel = paste(keyword, "(max no of links=", as.character(maxlinks), ")",sep="" )

  outimg = paste(filename, "_nolinksHistogram.png",sep="")
  openImgDev(outimg, iwidth = 600, iheight = 600)
  barplot(linkfreq, names.arg= as.character(linkidx), 
     xlab="number of links of a node", ylab="frequency", main= mylabel)

  #plot(linkidx, linkfreq,xlab="number of links of a node", ylab="frequency", main= keyword, type="h")
  #histogram(as.numeric(ntab), br=max(linkidx)-1, 
  #         xlab="number of links of a node", ylab="frequency", main= keyword)
  dev.off()

  # return
  c(no.links, no.uniquegenes, avglinks)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- abs(cor(x, y))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex * r)
}

cor1=function(x) {
 if(dim(x)[1]==1){
    omatrix = matrix(0, dim(x)[2], dim(x)[2] )
    colnames(omatrix) <- colnames(x)
    rownames(omatrix) <- colnames(x)
    return (omatrix)
 }

 col=cor(x,use="p", method="pearson")
 col=ifelse(is.na(col), 0, col)
 signif(col,2)
}

corSpearman=function(x) {
 #print(dim(x))
 if(dim(x)[1]==1){
    omatrix = matrix(0, dim(x)[2], dim(x)[2])
    colnames(omatrix) <- colnames(x)
    rownames(omatrix) <- colnames(x)
    return (omatrix)
 }

 col=cor(x,use="p", method ="spearman")
 col=ifelse(is.na(col), 0, col)
 signif(col,2)
}


# this function computes the standard error
stderror <- function(x){ sqrt( var(x)/length(x) ) }

# Error bars for barplot
# written by: Uli Flenker
# Institute of Biochemistry
# German Sports University
# Cologne Carl-Diem-Weg 6
# 50933 Cologne / Germany
# Phone 0049/0221/4982-493 
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)

err.bp<-function(daten,error,two.side=F){
 if(!is.numeric(daten)) {
      stop("All arguments must be numeric")
 }
 if(is.vector(daten)){ 
    xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
 }else{
    if (is.matrix(daten)){
      xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
dim=c(1,length(daten))))+0:(length(daten)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") 
    }
 }
 MW<-0.25*(max(xval)/length(xval)) 
 ERR1<-daten+error 
 ERR2<-daten-error
 for(i in 1:length(daten)){
    segments(xval[i],daten[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
    if(two.side){
      segments(xval[i],daten[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    } 
 } 
} 

#works for binary outcomes
krusktest=function( x ) {
   len1=dim(x)[[2]]-1
   out1=rep(666, len1);

   #if outcomes are the same value, no need to do the test
   totalOnes=sum(x[,1])
   if (totalOnes==0 | totalOnes == dim(x)[[1]]){
      for (i in c(1:len1) ) {out1[i]= 1.0}
   } else{
      for (i in c(1:len1) ) {out1[i]= signif( kruskal.test(x[,i+1], x[,1] )$p.value ,2) }
   }

   data.frame( variable=names(x)[-1] , kruskP=out1)
}

# access the pvalues of correlation between continuous vectors
#
CORtest=function( x ) {
   len1=dim(x)[[2]]-1
   no.observations= dim(x)[[1]]
   out1=rep(666, len1);
   for (i in c(1:len1) ) {
       if(no.observations>2){
         out1[i]= signif( cor.test(x[,i+1], x[,1] ,method="p",use="p")$p.value ,2) 
       }else{
         out1[i]= 1
       }
   }
   data.frame( variable=names(x)[-1] , cor.Pvalue=out1)
}

xcortest<-function(x,y){
   res<-cor.test(x, y)
   return( c(res$p.value, res$estimate) )
}

cor_multitraits <- function(x, ys){
    no.ys = dim(ys)[2]
    mres  = NULL
    for (j in c(1:no.ys)){
        jres = xcortest(x, ys[,j])
        mres = rbind(mres, jres)
    }
    return (mres)
}

# the following r<-->p can be verified by cor.test()
#
# r is correlation, n is the sample size
# t test is one-sided but F-test it two-sided
#
# t = -1.579, df = 8, p-value = 0.153
# alternative hypothesis: true correlation is not equal to 0 
#
# so, pvalue is the sum of the areas under the two tails  
#
# but t-test is only one-sided
#
computePvalueByR=function(r, n){
   t   =r*((n-2)/(1-r^2))^0.5
   prob=pt(t, df=n-2, lower.tail=F)
   return (2*prob)
}

compute_R_byPvalue=function(pvalues, n){
   zr = NULL
   for( pvalue in pvalues) {
   t=qt(pvalue/2,  df=n-2, lower.tail=F)
   r=t/(t^2+n-2)^0.5
   zr = c(zr, r)
   }

   return (zr)
}


corpval_multitraits <- function(x, ys){
    no.ys = dim(ys)[2]
    mres  = NULL
    for (j in c(1:no.ys)){
        sel = !(is.na(x) | is.na(ys[,j]) )
        if(sum(sel)>2) {
           jres = cor.test(x, ys[,j])
           mres = c(mres, signif(jres$p.value,3))
        } else {
           mres = c(mres, 1)
        }
    }
    return (mres)
}


# pvalue = 0.00021 ==>"2.1e-4"
pvalue2scientificrep=function(mypvalue){
   if(mypvalue==0){
      return("<e-22")
   }
   if(mypvalue>=1 | mypvalue<0){
      return(as.character(mypvalue))
   }

   base = 10
   cnt  = 0
   while(T){
     cnt     = cnt + 1
     mybase  = base^cnt
     foldval = mypvalue*mybase
     if(foldval>1.0){
        retvalue = paste(as.character(foldval), "e-", as.character(cnt),sep="")
        return (retvalue)      
     }
  }  
}

simpleCorTest=function(x,y){
 signif( cor.test(x,y,method="p",use="p")$p.value ,2) 
}

# no of rows of amatrix is the same as the length of myvect
corTest4multivects=function(myvect, amatrix){
 pvals = apply(amatrix, 2, simpleCorTest, y=myvect)
 #cat(pvals[], "\n")
 as.numeric(pvals)
}

# compute correlation coefficients (spearman, pearson), pvalues of the columns 
corRhoPvalSpearmanPearson = function (datMatrix) {

  rho=cor(datMatrix, method="pearson", use="complete.obs")
  pval=apply(datMatrix, 2, corTest4multivects, datMatrix)

  #datMatrixRank = apply(datMatrix, 2, rank)
  rhoR=cor(datMatrix, method="spearman", use="complete.obs")
  #pvalR=apply(datMatrixRank, 2, corTest4multivects, xdatMatrixRank)

  midx = getMatrixIndex(size=dim(rho), symmetric=TRUE, diagonal=FALSE)
  id1  = colnames(datMatrix)[midx[,1]]
  corMatrix = cbind(colnames(datMatrix)[midx[,1]], colnames(datMatrix)[midx[,2]],
              signif(rho[midx],3),signif(rhoR[midx],3), signif(pval[midx],3))
  colnames(corMatrix) = c("TraitA", "TraitB", "rho_pearson", "rho_spearman", "pvalue")

  return (corMatrix)
 
}


# compute correlation coefficients (spearman, pearson), pvalues of the columns 
corRhoPvalSpearmanPearson_TwoMatrices = function (datMatrix, datMatrix2) {

  rho=cor(datMatrix, datMatrix2, method="pearson", use="complete.obs")
  pval=apply(datMatrix, 2, corTest4multivects, datMatrix2)
  pval=t(pval)

  datMatrixRank = apply(datMatrix, 2, rank)
  datMatrixRank = matrix(as.integer(datMatrixRank), nrow=nrow(datMatrixRank))
  #datMatrixRankT = t(datMatrixRank)

  datMatrixRank2 = apply(datMatrix2, 2, rank)
  datMatrixRank2 = matrix(as.integer(datMatrixRank2), nrow=nrow(datMatrixRank2))
  #datMatrixRankT2 = t(datMatrixRank2)

  rhoR=cor(datMatrixRank, datMatrixRank2, method="pearson", use="complete.obs")
  pvalR=apply(datMatrixRank, 2, corTest4multivects, datMatrixRank2)
  pvalR=t(pvalR)

  midx = getMatrixIndex(size=dim(rho), symmetric=FALSE, diagonal=TRUE)
  id1  = colnames(datMatrix)[midx[,1]]
  corMatrix = cbind(colnames(datMatrix)[midx[,1]], colnames(datMatrix2)[midx[,2]],
              signif(rho[midx],3),signif(pval[midx],3), signif(rhoR[midx],3), signif(pvalR[midx],3))
  colnames(corMatrix) = c("TraitA", "TraitB", "rho_pearson", "pvalue_pearson", "rho_spearman", "pvalue_spearman")

  return (corMatrix)
 
}




## get P values of correlation coefficient matrix rmat, degree of
# freedom dfr (= no. samples used to calc each r - 2)
#
cor_pearson_pvalue <- function(rmat, dfr, is.sym = T)
{
   pmat <- matrix(0, nrow=nrow(rmat), ncol=ncol(rmat))
   if(is.sym == T)
   {
       above <- upper.tri(rmat)
       r2 <- rmat[above]^2
       Fstat <- r2 * dfr / (1 - r2)
       pmat[above] <- 1 - pf(Fstat, 1, dfr)
       pmat <- pmat + t(pmat)

   }
   else
   {
       r2 <- rmat^2
       Fstat <- r2 * dfr / (1 - r2)
       pmat <- 1 - pf(Fstat, 1, dfr)
   }

   return(pmat)
}


# get index of elements in a given vector
#
getMatchedIndex=function(cvector, subvect){
 subindex=NULL

 if(is.null(cvector)){return (NULL)} 

 fullindex = c(1:length(cvector) )
  for (each in subvect){
      sel=cvector== each
      if(sum(sel,na.rm=T)>0) {
         sel =ifelse(is.na(sel), F, sel)
         idx = fullindex[sel]
         subindex=c(subindex, idx)
      }
  }

  return (subindex)
}

# assume that subvect is a subset of cvector
getMatchedIndexFast=function(cvector, subvect){
  fullindex = c(1:length(cvector) )
  orgIdx    = cbind(cvector, fullindex)

  index2    = c(1:length(subvect))
  subIdex   = cbind(subvect, index2)

  merged    = merge(subIdex, orgIdx, by.x=1, by.y=1, all.x=T)
  if( dim(merged)[1]==0 ) {return(NULL)}

  merged    = as.matrix(merged)

  if(dim(merged)[1]>1){
    od        = order(as.integer(merged[,2]))  # restore the original order of subvect
    merged    = merged[od, ]
  }
  
  outIndex  = as.integer(merged[,3])

  return (outIndex)
}


equalSum = function(vect1, vect2){
   return( sum(vect1!=vect2))
}

# assume that subvect is a subset of cvector
getMatchedIndexFast_MatrixVect=function(cmatrix, cvect){
  
  msums = apply(cmatrix,1, equalSum, vect2=cvect)
  idx = c(1:nrow(cmatrix))[msums==0]

  if(length(idx)==0){
     return (-1)
  }
  return (idx)
}





# assume that find indices of the common components in the two sets
#
getMatchedIndex2way=function(cvector, dvect){
  fullindex = c(1:length(cvector) )
  orgIdx    = cbind(cvector, fullindex)

  index2    = c(1:length(dvect))
  subIdex   = cbind(dvect, index2)

  merged    = merge(orgIdx, subIdex, by.x=1, by.y=1, all=F)
  merged    = as.matrix(merged)
  
  outIndex  = cbind(as.integer(merged[,2]), as.integer(merged[,3]) )
  mo = order(outIndex[,1])

  outIndex = outIndex[mo,]

  return (outIndex)
}





factorization = function(mvect){
  mlevels   = names(table(mvect) )
  no.levels = length(mlevels)
  levelIdx  = c(1:no.levels)
  newvect = mvect
  for (i in levelIdx){
     newvect = ifelse(newvect==mlevels[i], i, newvect)
  }
  return ( as.integer(newvect))
}


# to split "abc|123", use sep="\\|", "abc.123" use "\\."
splitString =function(mystring, separator="; "){
  splitted = NULL
  for (each in mystring){
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     splitted =c(splitted, a)
  }
  #a=unlist( strsplit(mystring, separator) )
  return(splitted )
}

getTopNChars = function(mystrs, N=4){
    res = NULL
    for(str in mystrs){
       chars = splitString(str, "")
       cstr = concatenate(chars[c(1:N)], "")
       res  = c(res, cstr)
    }
    return (res)
}

splitStringsAsLists =function(mystring, separator="; "){
  nelemts = length(mystring)
  splitted = as.list(rep(NA, nelemts))
  for (i in c(1:nelemts) ){
     each = mystring[i]
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     splitted[[i]]=a
  }
  #a=unlist( strsplit(mystring, separator) )
  return(splitted )
}

# "HOXA9K"                  "RTHB6 /// LOC650428"     "TAKEDA_NUP8_HOXA9_8D_DN"
#   ===>
#
# "HOXA9K" "RTHB6"     "TAKEDA_NUP8_HOXA9_8D_DN"
# "HOXA9K" "LOC650428" "TAKEDA_NUP8_HOXA9_8D_DN"
#
splitMiddleElementNames = function(vect, elementIdx)
{
    #print(vect)

    igenes = splitString(vect[elementIdx], " /// ")
    igenes = setdiff(igenes, c("T", "L", "I") )
    igenes = replaceString(igenes, " ///", "")
    igenes = replaceString(igenes, " //", "")
    if(length(igenes)==1){
        newvect = vect; newvect[elementIdx]=igenes
        return (newvect)
    }

    # matrix construction row based
    #
    newMatrix = t(matrix(rep(vect, length(igenes)), nrow=length(vect)))
    newMatrix[,elementIdx] = igenes
    
    return (newMatrix)
}

elementInSplitString = function(vect, sep=" ", element=c("breast"))
{
    #print(vect)
    igenes = replaceString(vect, ", ", " ")
    igenes = replaceString(igenes, ",", " ")
    igenes = toupper(splitString(igenes, sep))
    
    return (is.element(element, igenes))
}


getStringLength = function(mystring, separator="; "){
  nelemts = length(mystring)
  mylen = rep(0, nelemts)
  for (i in c(1:nelemts) ){
     each = mystring[i]
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     mylen[i]=length(a)
  }
  return(mylen)
}


strINstr = function(substr, istring){
   splitted = splitString(istring, substr)
   if(length(splitted)==1){return (F);
   }else{return(T);}
}

strINstrings = function(substr, strings){
   ns  = length(strings)
   res = rep(F, ns)
   for(i in c(1:ns) ){
      istring = strings[i]
      if(istring==""){next}
      splitted = splitString(istring, substr)
      if(length(splitted)>1){res[i]=T;}
   }
   return (res)
}

union_length =function(mvect)
{ 
  return (length(union(mvect,NULL)))
}

concatenate=function(myvect, mysep="", do_union=FALSE, do_sort=FALSE)
{
  noitems = length(myvect)
  if (noitems==0){
    return ("")
  }else if (noitems==1){
    return (as.character(myvect) )
  }
  
  if(do_union) {
     myvect2=unique(as.character(myvect))
  } else{
     myvect2 = myvect
  }
  if(do_sort){myvect2=sort(myvect2)}
  concatenated <- paste(myvect2, sep="", collapse=mysep)

  return (concatenated)

  #tmpfn = "tmp.txt"
  #write.table(t(as.character(myvect)),tmpfn,sep=mysep,quote=FALSE, col.names=F, row.names=FALSE)
  #concatenated <- read.delim(tmpfn, sep="!", header=F)
  #return (as.character(as.matrix(concatenated) ))
}

concatenateSelectedElelments=function(boolvect, mnames, sep="-"){
    return(concatenate(mnames[boolvect], mysep=sep))
}


concatenateSelectedElelmentNames = function(vect, vnames, selected_value, mysep="_")
{
   sel = vect==selected_value
   if(sum(sel)==1){return(vnames[sel]) }

   res = concatenate(vnames[sel], mysep=mysep)   
   res
}


concatenateSelectedElelments2vect = function(boolvect, boolvect2, mnames, sep="-"){
    return(concatenate(mnames[boolvect&boolvect2], mysep=sep))
}

concatenateSelectedElelments_vectMatrix =function(boolvect, boolMatrix, mnames, sep="-"){
    res = apply(boolMatrix, 2, concatenateSelectedElelments2vect, boolvect2=boolvect, mnames=mnames, sep=sep)
    return (res)
}



concatenateSelectedElelmentNamesValues = function(vect, vnames, pcut=0.05, mysep="; ")
{
   od  = order(vect)
   vect2 = vect[od]
   vnames2 = vnames[od]

   sel = vect2<pcut
   combs = paste(vnames2, ": p=", vect2, sep="")

   if(sum(sel)==0){return ("") }

   if(sum(sel)==1){return(combs[sel]) }

   res = concatenate(combs[sel], mysep=mysep)   
   res
}

concatenateSelectedElelmentsValues_vectMatrix =function(valueMatrix, mnames, pcut = 0.05, sep="; "){
    #mnames = colnames(valueMatrix)

    res = apply(valueMatrix, 1, concatenateSelectedElelmentNamesValues, vnames=mnames, pcut =pcut, mysep=sep)
    return (res)
}




# correlation based ---------------------------------------------------------------------------------

concatenateSelectedElelmentNamesValues_ABS = function(vect, vnames, pcut=0.2, abs=TRUE, mysep="; ")
{
   od  = order(vect)
   vect2 = vect[od]
   vnames2 = vnames[od]

   if(abs) {
     sel = abs(vect2) >= pcut
   } else {
     sel = vect2 >= pcut
   }
   combs = paste(vnames2, ": r=", vect2, sep="")

   if(sum(sel)==0){return ("") }

   if(sum(sel)==1){return(combs[sel]) }

   res = concatenate(combs[sel], mysep=mysep)   
   res
}
concatenateSelectedElelmentsValues_vectMatrix_ABS =function(valueMatrix, mnames, pcut = 0.2, abs=TRUE, sep="; "){
    #mnames = colnames(valueMatrix)

    res = apply(valueMatrix, 1, concatenateSelectedElelmentNamesValues_ABS, vnames=mnames, pcut =pcut, abs=abs, mysep=sep)
    return (res)
}


concatenateSelectedElelmentNamesValues_ABS_module = function(vect, vnames, pcut=0.2, abs=TRUE)
{
   od  = order(vect)
   vect2 = vect[od]
   vnames2 = vnames[od]

   if(abs) {
     sel = abs(vect2) >= pcut
   } else {
     sel = vect2 >= pcut
   }
   combs = paste(vnames2, vect2, sep="\t")

   if(sum(sel)==0){return (" \t ") }

   if(sum(sel)==1){return(combs[sel]) }

   res = concatenate(combs[sel], mysep="\n")   
   res
}
concatenateSelectedElelmentsValues_vectMatrix_ABS_module =function(valueMatrix, mnames, pcut = 0.2, abs=TRUE){
    #mnames = colnames(valueMatrix)
    res = apply(valueMatrix, 1, concatenateSelectedElelmentNamesValues_ABS_module, vnames=mnames, pcut =pcut, abs=abs)
    return (res)
}



concatenateSelectedElelmentNamesValues_Neg_module = function(vect, vnames, pcut=0.2)
{
   od  = order(vect)
   vect2 = vect[od]
   vnames2 = vnames[od]

   sel = vect2 <= pcut

   combs = paste(vnames2, signif(vect2, 3), sep="\t")

   if(sum(sel)==0){return ("") }

   if(sum(sel)==1){return(combs[sel]) }

   res = concatenate(combs[sel], mysep="\n")  
   print(res)
 
   res
}
concatenateSelectedElelmentsValues_vectMatrix_Neg_module =function(valueMatrix, mnames, pcut = 0.2){
    #mnames = colnames(valueMatrix)
    res = apply(valueMatrix, 1, concatenateSelectedElelmentNamesValues_Neg_module, vnames=mnames, pcut =pcut)
    return (res)
}



########################################################################################################

# list files with certain keywords embedded inside
#
dirfiles=function(path, pattern, extension="."){
   mylist=dir(path, pattern=extension)
   matches=NULL
   for (each in mylist){
       splitted = splitString(each, pattern)
       if (length(splitted) > 1){
           matches = c(matches, each)
       }
   }

   return (matches)
}


#get the filename without extension
#
getFileExtension=function(fullfname){
    splitted=unlist( strsplit(fullfname, "\\.") )
    
    if( length(splitted) >1){
      return (splitted[length(splitted)])
    } else{
      return ("")
    }
}

#get the filename without extension
getFileName=function(fullfname){
    ext=getFileExtension(fullfname)
    if(ext ==""){
       return (fullfname)
    }
    extd = paste(".", ext, sep="")
    splitted=splitString(fullfname, extd)

    splitted[1]
}

#get the filename without extension
getFileNames=function(fullfnames){

  final = NULL
  for(ef in fullfnames) {
     fn = getFileName(ef)
     final = c(final, fn)
  }
  return (final)
}


#get the filename without extension
getFileNameOld=function(fullfname){
    splitted=unlist( strsplit(fullfname, "\\.") )
    mfname=""
    for ( i in c(1:(length(splitted)-1)) ){
       mfname =paste(mfname, splitted[i], sep="")
    }
    mfname
}



# get second part: 31357-31351 ==> 31351 
#
getSecondPart=function(fullfnames, sep="-", whichpart=-1, retBlank_ifNoMatch=F){

   n.elements = length(fullfnames)
   p1 = rep(T, n.elements);
   if(sep!="") {
     p1= strINstrings(sep, fullfnames)
   }
  
  ret = NULL
  for(each in fullfnames) {
    splitted=unlist( strsplit(each, sep) )
    if (whichpart==-1) {
       reti=splitted[ length(splitted) ]
    } else {
       if ( whichpart > length(splitted) ) {
          reti= each
       } else {
          reti=splitted[whichpart]
       }
    }

    ret = c(ret, reti)
  }

  if(retBlank_ifNoMatch) {
    ret[!p1] = ""
  }

  ret
}

getSecondPartFast=function(fullfnames, sep="-", whichpart=-1, retBlank_ifNoMatch=F){
  ret=sapply(fullfnames, getSecondPart, sep=sep, whichpart=whichpart, retBlank_ifNoMatch=retBlank_ifNoMatch)
  as.character(ret)
}

getInBetween=function(fullfnames, presep="-", postsep="", skipfirst=F, skiplast=F){
   n.elements = length(fullfnames)
   p1 = rep(T, n.elements); p2 = rep(T, n.elements); 
   if(presep!="") {
     p1= strINstrings(presep, fullfnames)
   }
   if(postsep!="") {
     p2= strINstrings(postsep, fullfnames)
   }

   splitted = getSecondPart(fullfnames,   presep, 2)
   if(skipfirst){splitted=splitted[-1]}

   splitted2= getSecondPart(splitted, postsep, 1); nl = length(splitted2)
   if(skiplast){splitted2=splitted2[-nl]}

   sel = !(p1 & p2)
   splitted2[sel] =""

   return(splitted2)
}

getUrlTopField = function(urls) {
   xUrls = replaceString(urls, "http:\\/\\/", "")
   xUrls = replaceString(xUrls, "www.", "")
   xUrls = getSecondPart(xUrls, "\\/",1)
   return(xUrls)
}


removeLastPart=function(fullfnames, sep="-"){

  ret = NULL
  for(each in fullfnames) {
    splitted=unlist( strsplit(each, sep) )
    reti=splitted[ -length(splitted) ]
    newstr = concatenate(reti, sep)
    ret = c(ret, newstr)
  }

  ret
}


getAllParts=function(fullfnames, sep="-", max_len=NULL, ret_len=F){

  splitted=unlist( strsplit(fullfnames[1], sep) )
  nn  = length(fullfnames)
  nf  = length(splitted)
  if(!is.null(max_len)){nf=max(nf, max_len)}

  ret = matrix("", nn, nf)
  lens= rep(0, nn)
  for(i in c(1:nn) ) {
    each = fullfnames[i]
    splitted=unlist( strsplit(each, sep) )
    ino     = length(splitted)
    if(ino >=nf) {
       ret[i,] = splitted[1:nf]
    }else{
       ret[i,] = c(splitted, rep("",nf-ino ))
    }

    lens[i] = ino
  }

  if ( ret_len ){
     return (lens)
  } else {
     return (ret) 
  }
  
}


getAllPartsVect = function(fullfname, sep="-"){

  splitted=unlist( strsplit(fullfname, sep) )
  return (splitted) 
   
}



#get the filename without extension
getOntologyNameFromPath=function(fullfname){
    fn=getFileName(fullfname)
    splitted=unlist( strsplit(fn, "\\_") )
    splitted[length(splitted) ]
}


#get the filename without extension and path information
getFileNameNopath=function(fullfname){
   res = NULL
   for(each in fullfname) {
    myfilename = getFileName(each)
    splitted=unlist( strsplit(myfilename, "/") )
     res= c(res, splitted[length(splitted) ])
   }
   return (res)
}

#get the filename without path information
getFileFullNameNopath=function(fullfnames){
   res = NULL
   for(each in fullfnames) {
     splitted=unlist( strsplit(each, "/") )
     res= c(res, splitted[length(splitted) ])
   }
   return (res)
}


# get a particulr field of splitted string
getAFieldBySplit=function(strvect,separator="_", fieldId=1)
{
  splitted = NULL
  for (each in strvect){
     a=unlist( strsplit(each, separator) )
     splitted =c(splitted, a[fieldId])
  }
  return(splitted )
}

maxNumberCharacters = function(mstrings){
   maxlen = 0
   for(intStr in mstrings) {
     digits    = unlist(strsplit(intStr,''))# get single digits for the given number
     ndigits   = length(digits)
     if(ndigits>maxlen){maxlen=ndigits}
   }
   maxlen
}

patchNcharsVect = function(intStrs, patched="0", number=3, prefix=TRUE)
{
  res=NULL
  for(ec in intStrs) {
    ires= patchNchars(ec, patched=patched, number=number, prefix=prefix)
    res = c(res, ires)
  }
  return (res)
}

patchNchars=function(intStr, patched="0", number=3, prefix=TRUE)
{

   #print(intStr)

   zeros     = rep(patched, number)

   digits    = unlist(strsplit(intStr,''))# get single digits for the given number
   ndigits   = length(digits)
   
   # put back digits in the zeros
   for (i in c(1:ndigits) ) {
      zidx       = number-i+1
      if(prefix) {
         zeros[zidx]= digits[ndigits-i+1]
      } else {
         zeros[i]= digits[i]
      }
   }
   finalstr=""
   for (each in zeros){
      finalstr=paste(finalstr, each, sep="")
   }
   return(finalstr)
}


replaceChars = function(stringvect, oldchar, newchar){
   nelems=length(stringvect)

   final = rep("", nelems)
   for (i in c(1:nelems) ){
       each = stringvect[i]
       final[i] = replaceCharsCore(each, oldchar, newchar)
   }
   return(final)
}


replaceCharsCore=function(fullfname, oldchar, newchar){

    if( (fullfname=="") | is.na(fullfname) ){
        return (fullfname)
    }

    #splitted=strsplit(fullfname,oldchar)
    #if(length(splitted)==1){return(fullname)}
    #paste()

    splitted=strsplit(fullfname,'')
    
    i= length( splitted[[1]] )
    for (j in c(1:i) ){
	if (splitted[[1]][j] == oldchar){
           splitted[[1]][j] = newchar
        }
    }
    mfname =''
    for (j in c(1:i) )
       mfname = paste(mfname,splitted[[1]][j],sep='')
    mfname
}


replaceString=function(fullfnames, oldstr, newstr){

  no.files = length(fullfnames)
  res = NULL
  for(each in fullfnames) {
    #print(paste(i, "/", no.files, ":", each) )

    each2 = paste(each, oldstr, sep="")
    splitted = splitString(each2, oldstr)

    neweach = concatenate(splitted, newstr)

    if(F) {
    neweach  =""
    for (is in splitted) {
       neweach= paste(neweach, newstr, sep=is)
    }
    }

    #oldeach  = paste(pathnet, each,   sep="")
    #neweach  = paste(pathnet, newstr, splitted[2], sep="")
    #a=file.rename(from=oldeach, to=neweach)
    #print(a)

    res = c(res, neweach)
  }
  return (res)  
}

replaceStringRecursive=function(fullfnames, oldstr, newstr){
   res=NULL
   for(each in fullfnames){
      old = each; new2=""
      while(T) {
          new2 = replaceString(old, oldstr, newstr)
          if(new2==old){break}
          old = new2
      }
      res = c(res, new2)
   }
   return (res)
}


#get the filename without extension
reverseString=function(mystring){
    splitted=strsplit(mystring,'')
    i= length( splitted[[1]] )
    mfname =''
    for (j in c(i:1) )
       mfname = paste(mfname,splitted[[1]][j],sep='')
    mfname
}

splitAndAverage = function(str, sep="," ){
    pt = getAllParts(str, sep=sep)
    pt = ifelse(pt=="X", NA, pt)
    pt = ifelse(pt=="x", NA, pt)
    pt = ifelse(pt=="" , NA, pt)

    mn = mean(as.numeric(pt), na.rm=TRUE)
    return (mn)
}

splitAndAverageVect = function(strvect, sep=",") {
    nvect = lapply(strvect, splitAndAverage, sep=sep)
    return (unlist(nvect))
}


#get the filename without extension
DNApair=c("A","C","G","T")
names(DNApair) <- c("C","A","T","G")
complementDNAString=function(mystring, mappingTable=DNApair){
    splitted=strsplit(mystring,'')
    i= length( splitted[[1]] )
    mfname =''
    for (j in c(i:1) ){
       jchar = splitted[[1]][j]
       cchar = as.character(mappingTable[jchar]) #get complement character
       mfname = paste(mfname, cchar,sep='')
    }
    mfname
}


#get the filename without extension
reverseString4apply=function(mystringvect){
    splitted=strsplit(as.character(mystringvect[1]),'')
    i= length( splitted[[1]] )
    mfname =''
    for (j in c(i:1) )
       mfname = paste(mfname,splitted[[1]][j],sep='')
    mfname
}


appendStringToFile=function(fname, mstring, newline=T){
    fp <- file(fname, "a")
    if(newline){
     cat(mstring, "\n", file=fp, sep="")
    }else{
     cat(mstring, file=fp)
    }
    close(fp)    
}

appendListToFile=function(fname, mlist, listTitle=""){
  fp <- file(fname, "a")
  if (length(listTitle)>0){
       cat(as.character(listTitle), "\n", file=fp)
  }  
      no.fields=length(mlist)
      #write column title
      for (z in 1:no.fields ){
         cat(as.character(names(mlist[z])),"\t", file=fp)
      }      
      cat("\n", file=fp)

      for (z in 1:no.fields ){
         itab=mlist[z]
         cat(as.character(itab[[1]]),"\t", file=fp)
      }      
      cat("\n", file=fp)
      close(fp)
 }


appendTableToFile=function(fname, mtable, tableTitle="", myappend=T){
   if ( is.null(mtable) ){
     return
   }

   if(myappend==T){
     fp <- file(fname, "a")    
   }else{
     fp <- file(fname, "w")
   }
   if (tableTitle != "" ){
     cat(as.character(tableTitle), "\n", file=fp)    
   }

  if ( (!is.na( dim(mtable)[1])) & is.na( dim(mtable)[2]) ) {#only one row in the table
      #write column title
      coltitles=names(mtable)
      for (z in 1:length(coltitles) ){
         cat(as.character(coltitles[z]),"\t", file=fp)
      }
      cat("\n", file=fp)

      for (i in 1:(dim(mtable)[1]) ){
          cat(as.character(mtable[i]), "\t", file=fp)
      }
      cat("\n", file=fp)

   }else{ # normal table
       cat(" \t", file=fp)
       #write column title
       coltitles=colnames(mtable)
       for (z in 1:length(coltitles) ){
          cat(as.character(coltitles[z]),"\t", file=fp)
       }
       cat("\n", file=fp)

       rowsname = rownames(mtable)
       for (i in 1:(dim(mtable)[1]) ){
           cat(as.character(rowsname[i]), "\t", file=fp)
           for(j in 1:(dim(mtable)[2])){
              cat(as.character(mtable[i, j]), "\t", file=fp)
            }
           cat("\n", file=fp)
       }   
   }
   cat("\n", file=fp)
   close(fp)
}

appendMultiTablesToFile=function(fname, multitables){
    #table of tables
    if(is.na( dim(multitables)[2]) ) {
        titles=names(multitables)
        for (i in 1:(dim(multitables)[1]) ){
            if ( is.null(multitables[[i]]) )
                 next
            appendTableToFile(fname, multitables[[i]], as.character(titles[i]))
         }
    }else{#single table
      appendTableToFile(fname, multitables)
    }
}

#iunits: "px", "cm", "mm"
#
# Letter Size 8.5 x 11 inches
# images widths: 8.3cm ~ 3.26in, 12.35cm ~ 4.86in, 17.35cm ~ 6.83in
#        heights: max 23.35cm ~ 9.19in
#
openImgDev=function(imgname, iwidth = 1024, iheight=1024, ipointsize=12, iunits="px", ires=72, icompression="lzw")
{
  imgtype = getFileExtension(imgname)
  
  if (imgtype=="ps"){
     postscript(imgname,width=iwidth, height=iheight, pointsize=ipointsize)
  }else if (imgtype=="png"){
     png(imgname, width=iwidth, height=iheight, pointsize=ipointsize, units=iunits, res=ires)
  }else if ( (imgtype=="jpg") | (imgtype=="jpeg") ){
     jpeg(imgname, width=iwidth, height=iheight, pointsize=ipointsize,units=iunits, res=ires,quality =100)
  }else if ( (imgtype=="tif")|(imgtype=="tiff") ){
     tiff(imgname, width=iwidth, height=iheight, pointsize=ipointsize,units=iunits, res=ires,compression=icompression)
  }else if ( (imgtype=="pdf") | (imgtype=="PDF") ){
     pdf(imgname, width=iwidth, height=iheight, pointsize=ipointsize)
     #return
  }else{
     png(imgname, width=iwidth, height=iheight, pointsize=ipointsize, units=iunits, res=ires)
  }
  trellis.device(new = FALSE, col = TRUE) 
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: compute sigmoid functions, single/vector inputs
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#compute sigmoid function
sigmoid <- function(x,alpha,tau)
{
  y=1.0+exp( -alpha*(x-tau) )
  y=1.0/y
  y
}

#compute sigmoid values for a matrix of inputs
sigmoidMatrix <- function(xmatrix,alpha,tau0)
{
 ey = exp( -alpha*(xmatrix - tau0) )
 #ey = 3^( -alpha*(xmatrix - tau0) )
 y= 1.0/(1.0+ey)
 rm(ey)
 y
}



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: correlation cutoff Computation 
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#input is a matrix with Pearson correlations

#iteratively call gc() untill there is no change in memory hit
#i.e., all freed memory has been released to the system
#Usage: 1. immediately call this function after you call a function or
#       2. rm()
collect_garbage=function(){
    while (gc()[2,4] != gc()[2,4]){}
}

dichotcut = function(corrlMatrix,noOFsamples, mlogFname, samplingSize=3000, RsquaredCut=0.8, 
mincutval=0.2, maxcutval=0.95, cutstep=0.01) {
        orgSize   = dim(corrlMatrix)[1] 
        #cat('samplingSize', samplingSize)

        # perform sampling 	
        subsetsize=min(samplingSize, orgSize)
        subset1=sample(c(1:orgSize),subsetsize,replace=F)
        cor1 <- corrlMatrix[subset1, subset1]

        # perform sampling
        sampled = T	
        if(samplingSize<=0 |samplingSize >= orgSize){
            sampled = F
            no.genes= orgSize
        } else{
            subsetsize=min(samplingSize, orgSize)
            subset1   =sample(c(1:orgSize),subsetsize,replace=F)
            cor1 <- corrlMatrix[subset1, subset1]
            no.genes   <- dim(cor1)[[2]]
        }

        cutvector=seq(mincutval, maxcutval, by=cutstep)

        colname1=c("Cut","p-value", "Adj R^2","Truncated Adj R^2", "slope","mean(k)","median(k)","max(k)")

        datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
        names(datout)=colname1   

        for (i in c(1:length(cutvector) ) ){
             cut1=cutvector[i]

             if(sampled) {
                 dichotCorhelp  <-  I(abs(cor1)>cut1 )+0.0
             }else{
                 dichotCorhelp  <-  I(abs(corrlMatrix)>cut1 )+0.0
             }

             diag(dichotCorhelp)<- 0
             nolinkshelp <- apply(dichotCorhelp,2,sum, na.rm=TRUE) 

             print( paste(cut1, ": ", no.genes, "genes, ", sum(nolinkshelp), " links") )

             str_cutoff=paste("tau=", as.character(cut1), sep="")
             fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             pvalue  = 2*(1-pt(sqrt(noOFsamples-1)*cut1/sqrt(1-cut1^2),noOFsamples-1))
             datout[i,]=signif(c( cut1, pvalue, fitness[1], fitness[2], fitness[3], mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)

             if(1==2){# oldway to do it
             # let's check whether there is a scale free topology
             no.breaks=15
             cut2=cut(nolinkshelp,no.breaks)

             freq1=tapply(nolinkshelp,cut2,length)/length(nolinkshelp)
             binned.k=tapply(nolinkshelp, cut2,mean)

             #remove NAs
             noNAs = !(is.na(freq1) | is.na(binned.k))
             freq.noNA= freq1[noNAs]
             k.noNA  = binned.k[noNAs]

             #remove Zeros
             noZeros  = !(freq.noNA==0 | k.noNA==0) 
             freq.log = as.numeric(log10(freq.noNA[noZeros]))
             k.log    = as.numeric(log10(k.noNA[noZeros]))

             lm1=lm( freq.log ~ k.log, na.action=na.omit) 
             lm2=lm( freq.log ~ k.log +I(10^freq.log) );

             pvalue=2*(1-pt(sqrt(noOFsamples-1)*cut1/sqrt(1-cut1^2),noOFsamples-1))
             datout[i,]=signif(c( cut1, pvalue, summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared, lm1$coef[[2]],mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)
             }
        }

        #output data
        #datout[length(cutvector)+1, ] = c( "Selected Cutoff = ", cutvector[indcut][[1]],"","","","","")
        print(datout);
        write.table(datout, mlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

        corrlcutoff=NA
        msqrcut = RsquaredCut
        while(is.na(corrlcutoff)){
             ind1   = datout[,3] > msqrcut
             ind1   = ifelse(is.na(ind1), F, ind1)

             indcut = NA
             indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
             corrlcutoff = cutvector[indcut][[1]]
             msqrcut     = msqrcut - 0.05
        }

        # consider the truncated scalefree index if the scalefree index is not well satisfied
        if (msqrcut+0.05<0.75){
            truncatedSel = (datout[,4] >= 0.9) & (datout[,6]<=60) & (datout[,3] >=0.5)
            if (sum(truncatedSel)>0){
                 colcutoffs = datout[truncatedSel, 1]
                 corrlcutoff= colcutoffs[1]
            }
        }

        corrlcutoff
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: correlation POWER Computation 
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#input is a matrix with Pearson correlations
powercut = function(corrlMatrix,noOFsamples, mlogFname, samplingSize=3000, RsquaredCut=0.8, 
mincutval=1,maxcutval=12, bystep=0.5, plotimg=T, maxMedianLinks=60, minMedianLinks=30) {
        orgSize   = dim(corrlMatrix)[1] 
        #cat('samplingSize', samplingSize)

        # perform sampling
        sampled = T	
        if(samplingSize<=0 |samplingSize >= orgSize){
            sampled = F
            no.genes= orgSize
        } else{
            subsetsize=min(samplingSize, orgSize)
            subset1   =sample(c(1:orgSize),subsetsize,replace=F)
            cor1 <- corrlMatrix[subset1, subset1]
            no.genes   <- dim(cor1)[[2]]
        }

        cutvector=seq(mincutval, maxcutval, by=bystep)

        colname1 =c("Cut","p-value", "Adj R^2","Truncated Adj R^2", "slope","mean(k)","median(k)","max(k)")

        datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
        names(datout)=colname1   

        for (i in c(1:length(cutvector) ) ){
             cut1=cutvector[i]
             if(sampled) {
                 dichotCorhelp      <-  abs(cor1)^cut1
             }else{
                 dichotCorhelp      <-  abs(corrlMatrix)^cut1
             }

             diag(dichotCorhelp)<- 0
             nolinkshelp <- apply(dichotCorhelp,2,sum, na.rm=TRUE) 

             print( paste(cut1, ": ", no.genes, "genes, ", sum(nolinkshelp), " links") )

             str_cutoff=paste("beta=", as.character(cut1), sep="")
             #fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             if(plotimg) {
               fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             }else{
               fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="sclfree.png", outputFitness=TRUE)
             }

             pvalue  = -1

             datout[i,]=signif(c( cut1, pvalue, fitness[1], fitness[2], fitness[3], mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)

             if(1==2){# oldway to do it
             # let's check whether there is a scale free topology
             no.breaks=15
             cut2=cut((nolinkshelp+1),no.breaks)

             freq1=tapply(nolinkshelp,cut2,length)
             binned.k=tapply(nolinkshelp+1,cut2,mean)

             #remove NAs
             noNAs = !(is.na(freq1) | is.na(binned.k))
             freq.noNA= freq1[noNAs]
             k.noNA  = binned.k[noNAs]

             #remove Zeros
             noZeros  = !(freq.noNA==0 | k.noNA==0) 
             freq.log = as.numeric(log10(freq.noNA[noZeros]))
             k.log    = as.numeric(log10(k.noNA[noZeros]))

             lm1=lm( freq.log ~ k.log, na.action=na.omit) 
             lm2=lm( freq.log ~ k.log +I(10^freq.log) );

             #pvalue=2*(1-pt(sqrt(noOFsamples-1)*cut1/sqrt(1-cut1^2),noOFsamples-1))
             pvalue=-1
             datout[i,]=signif(c( cut1, pvalue, summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared, lm1$coef[[2]],mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)
             }
        }

        print(datout);

        #in case that the largest R^2 is not bigger than RsquaredCut,
        #we have to reduce it gradually
        corrlcutoff=NA
        msqrcut = RsquaredCut
        xcnt = 0
        while(is.na(corrlcutoff) & (xcnt<length(cutvector) ) ){
             ind1   = (datout[,3] > msqrcut) & (datout[,7]<=maxMedianLinks)  & (datout[,7]>=minMedianLinks)
             if(sum(ind1)>0){
             indcut = NA
             indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
             corrlcutoff = cutvector[indcut][[1]]
             }
             msqrcut     = msqrcut - 0.05;
             xcnt = xcnt + 1;
             
        }

        # consider the truncated scalefree index if the scalefree index is not well satisfied
        if (F) {
        if (msqrcut+0.05<0.75){
            truncatedSel = (datout[,4] >= 0.9) & (datout[,7]<=maxMedianLinks) & (datout[,3] >=0.5)
            if (sum(truncatedSel)>0){
                 colcutoffs = datout[truncatedSel, 1]
                 corrlcutoff= colcutoffs[1]
            }
        }}

        if(is.na(corrlcutoff)){ # just based on medianLinks
            ind1 = datout[,7]<=maxMedianLinks
            if ( sum(ind1) >0 ) {
               indcut=min(c(1:length(ind1))[ind1])
               corrlcutoff = cutvector[indcut][[1]]
            } else{# so max link > maxMedianLinks
               sel=(datout[,5]<0); mR2=ifelse(sel, datout[,3], 0)
               corrlcutoff = datout[which.max(mR2), 1]
            }
        }

       #output data
       #datout[length(cutvector)+1, ] = c( "Selected Cutoff = ", cutvector[indcut][[1]],"","","","","")

       write.table(datout, mlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

       corrlcutoff
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: correlation Sigmoid Computation, search for alpha
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#input is a matrix with Pearson correlations
sigmoidcut = function(corrlMatrix,noOFsamples, mlogFname, samplingSize=3000, RsquaredCut=0.8, tau=0.5, out=0) {
        orgSize   = dim(corrlMatrix)[1] 
        #cat('samplingSize', samplingSize)

        mincutval=1
        maxcutval=12

        # perform sampling
        sampled = T	
        if(samplingSize<=0 |samplingSize >= orgSize){
            sampled = F
            no.genes= orgSize
        } else{
            subsetsize=min(samplingSize, orgSize)
            subset1   =sample(c(1:orgSize),subsetsize,replace=F)
            cor1 <- corrlMatrix[subset1, subset1]
            no.genes   <- dim(cor1)[[2]]
        }


        cutvector=seq(mincutval, maxcutval, by=0.2)

        colname1=c("tau", "alpha","p-value", "Adj R^2","Truncated Adj R^2", "slope","mean(k)","median(k)","max(k)")

        datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
        names(datout)=colname1   

        for (i in c(1:length(cutvector) ) ){
             cut1=cutvector[i]
             dichotCorhelp      <-  sigmoidMatrix( abs(cor1), cut1, tau)
             collect_garbage()

             if(sampled) {
                 dichotCorhelp   <-  sigmoidMatrix( abs(cor1), cut1, tau)
             }else{
                 dichotCorhelp   <-  sigmoidMatrix( abs(corrlMatrix), cut1, tau)
             }


             diag(dichotCorhelp)<- 0
             nolinkshelp <- apply(dichotCorhelp,2,sum, na.rm=TRUE) 

             str_cutoff=paste("beta=", as.character(cut1), sep="")
             fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             pvalue  = -1
             datout[i,]=signif(c(tau, cut1, pvalue, fitness[1], fitness[2], fitness[3], mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)



             if(1==2){# oldway to do it
             # let's check whether there is a scale free topology
             no.breaks=15
             cut2=cut((nolinkshelp+1),no.breaks)

             freq1=tapply(nolinkshelp,cut2,length)
             binned.k=tapply(nolinkshelp+1,cut2,mean)

             #remove NAs
             noNAs = !(is.na(freq1) | is.na(binned.k))
             freq.noNA= freq1[noNAs]
             k.noNA  = binned.k[noNAs]

             #remove Zeros
             noZeros  = !(freq.noNA==0 | k.noNA==0) 
             freq.log = as.numeric(log10(freq.noNA[noZeros]))
             k.log    = as.numeric(log10(k.noNA[noZeros]))

             lm1=lm( freq.log ~ k.log, na.action=na.omit) 
             lm2=lm( freq.log ~ k.log +I(10^freq.log) );

             #pvalue=2*(1-pt(sqrt(noOFsamples-1)*cut1/sqrt(1-cut1^2),noOFsamples-1))
             pvalue=-1
             datout[i,]=signif(c(tau, cut1, pvalue, summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared, lm1$coef[[2]],mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)
             }

             rm(dichotCorhelp)
             collect_garbage()
        }

        #in case that the largest R^2 is not bigger than RsquaredCut,
        #we have to reduce it gradually
        corrlcutoff=NA
        msqrcut = RsquaredCut
        while(is.na(corrlcutoff)){
             ind1   = datout[,4] > msqrcut
             indcut = NA
             indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
             corrlcutoff = cutvector[indcut][[1]]
             msqrcut     = msqrcut - 0.05
        }

       #output data
       #datout[length(cutvector)+1, ] = c( "Selected Cutoff = ", cutvector[indcut][[1]],"","","","","")
       print(datout);

       if (out==0){
          returnval = corrlcutoff
          write.table(datout, mlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
       }
       else
          returnval = datout
       #return returnval
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##        Function: correlation Sigmoid Computation, search for alpha and tau
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigmoidcutTauAlpha = function(mcorrlMatrix,mnoOFsamples, mmlogFname, msamplingSize=3600, mRsquaredCut=0.8) 
{
   mintau = 0.0
   maxtau = 1.0

   tauvector=seq(mintau, maxtau, by = 0.02)
   
   allData = c()
   for (i in tauvector ){
      idata = sigmoidcut(corrlMatrix=mcorrlMatrix, noOFsamples=mnoOFsamples, mlogFname=mmlogFname, 
                         samplingSize=msamplingSize, RsquaredCut=mRsquaredCut, tau=i, out=1) 
      allData = rbind(allData, idata)
   }
   write.table(allData, mmlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
   
   corrlcutoff=NA
   msqrcut = mRsquaredCut
   while(is.na(corrlcutoff[1])){
     ind1   = allData[,4] > msqrcut
     indcut = NA
     indcut = ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)

     corrlcutoff = allData[ indcut[[1]], c(1:2)]
     msqrcut     = msqrcut - 0.05
   }

   #a 2-d vector with [1]: tau and [2]: alpha
   corrlcutoff
}


#------------------------------------------------------------------------------------------------------------------
#Function: get most connected genes based on a scale free network derived from a power adjacency function
#
#Input:
#    minputfname    ~  gene expresion matrix
#    geneinforCols  ~  number of columns as gene information
#    R2Cut          ~  scale free fitness:ouput most connected genes from a the network with scale free fitness>R2Cut OR
#    tructR2Cut     ~  scale free fitness > tructR2Cut
#    topN           ~  top N most connected genes will be output 
#    mincutval      ~  min power value
#    maxcutval      ~  max power value
#
#Output:  1) scale free plot 2) logfile 3) array data with most connected genes
#
#------------------------------------------------------------------------------------------------------------------
getMostConnectedGenesByPowerAdj = function(minputfname, geneinforCols, R2Cut=0.85, tructR2Cut=0.90, topN=3600, mincutval=4, maxcutval=16, maxMeanLinks=30)
{
       fname       =getFileName(minputfname)
       mlogFname   =paste(fname, "_logo.txt",       sep='')
       imgScaleFree=paste(fname, "_imgScaleFree.png",sep='')

       #------- STEP 0: read in gene information, expression data
       allMatrix <- read.delim(minputfname,sep="\t", header=T)
       dim(allMatrix)

       genesInfor <- allMatrix[,c(1:geneinforCols)]
       rowTitles=names(allMatrix)

       #These are the expression values
       datExpr <- t(allMatrix[,-c(1:geneinforCols)])
       no.samples <- dim(datExpr)[1]
       dim(datExpr)

       corhelp <- cor(datExpr, use = "pairwise.complete.obs") 
       no.genes <- dim(corhelp)[[2]]
       dim(corhelp)

       diag(corhelp) <- 0

        cutvector=seq(mincutval, maxcutval, by=1)

        colname1=c("Cut","p-value", "Adj R^2","Truncated Adj R^2", "slope","mean(k)","median(k)","max(k)")
        datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
        names(datout)=colname1

        for (i in c(1:length(cutvector) ) ){
             cut1=cutvector[i]
             
             #get no links for each gene
             nolinkshelp   = rep(0.0, no.genes)
             for (j in c(1:no.genes)){
                sel =  as.numeric(abs(corhelp[j,])) ^ cut1
                nolinkshelp[j]  = sum(sel)
             }

             str_cutoff= paste("beat=", as.character(cut1), sep="")
             fitness = ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="", outputFitness=TRUE)
             pvalue=-1
             datout[i,]=signif(c( cut1, pvalue, fitness[1], fitness[2], fitness[3], mean(nolinkshelp), median(nolinkshelp), max(nolinkshelp) ), 3)

             if( mean(nolinkshelp)<maxMeanLinks & (fitness[1] >= R2Cut || fitness[2]>tructR2Cut) ){
                break
             }
        }
       print(datout);
       write.table(datout[c(1:i),], mlogFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

       str_cutoff=paste("beta=",as.character(cut1), sep="")
       suminfor=ScaleFreePlot(nolinkshelp,no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile="")
       ScaleFreePlot(nolinkshelp, no.breaks=15, mtitle=str_cutoff,truncated1=TRUE, tofile=imgScaleFree)

       #*----------------------------- choose the top N genes ----------------------------------------
       orderLink = order(-nolinkshelp)
       nolinks.ordered = nolinkshelp[orderLink]
     
       restFname = paste(fname, "_p", as.character(topN), ".xls",       sep='')
       finalMatrix = allMatrix[orderLink[1:topN], ]
       write.table(finalMatrix, restFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE) 
       str_cutoff
}



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                          Function: automatically detect & label modules
##
##
##   nameallmodules=FALSE: label modules with all possible colors
##                 =TRUE:  when # of modules exceeds length(colorcode), we use false color 
##                          names to label the reamining modules
##
##   useblackwhite=FALSE: label as normal
##                =TRUE:  label extra modules by black and white alternatively
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#The input is an hclust object.
moduleDetectLabel = function(hiercluster,heightcutoff=0.5,minsize1=20, nameallmodules=FALSE, useblackwhite=FALSE, useNumberAsLabel=F, startlabel=1) {

    # here we define modules by using a height cut-off for the branches
    labelpred= cutree2(hiercluster,h=heightcutoff)
    sort1    =-sort(-table(labelpred))
    #sort1
    modulename   = as.numeric(names(sort1))
    modulebranch = sort1 > minsize1
    no.modules   = sum(modulebranch)

    if (useNumberAsLabel){
       # now make cluster label
       #
       colorcode=NULL
       for (i in c(startlabel:(startlabel+no.modules-1)) ){
          ipad = patchZeros(i)
          colorcode= c(colorcode, ipad)
       }
    }else{
    # now we assume that there are fewer than 10 modules
    colorcodex=c("turquoise",    "blue",     "brown",   "yellow",      "green",      "red",     "black",
                "pink",         "magenta",  "purple",  "greenyellow", "tan",        "salmon",  "cyan", 
                "midnightblue", "lightcyan","grey60",  "lightgreen",  "lightyellow","coral",   "sienna",
                "gold",         "peru",     "wheat",   "chocolate",   "seashell",   "khaki",   "bisque",
                "forestgreen",  "navy",     "plum",    "mediumblue",  "violet",     "hotpink",
                "thistle",      "orchid",   "maroon",  "violetred",   "firebrick",  "honeydew","chartreuse",
                "deeppink",     "darkcyan", "beige",   "snow",        "burlywood",  "goldenrod",
                "brown2",       "red2",     "gold2",   "yellow2",     "green2",     "cyan2",    "blue2",
                "brown3",       "red3",     "gold3",   "yellow3",     "green3",     "cyan3",    "blue3",
                "brown4",       "red4",     "gold4",   "yellow4",     "green4",     "cyan4",    "blue4",
                "gray1","gray2","gray3","gray4","gray5","gray6","gray7","gray8","gray9","gray10",
                "gray11","gray12","gray13","gray14","gray15","gray16","gray17","gray18","gray19","gray20",
                "gray21","gray22","gray23","gray24","gray25","gray26","gray27","gray28","gray29","gray30",
                "gray31","gray32","gray33","gray34","gray35","gray36","gray37","gray38","gray39","gray40",
                "gray41","gray42","gray43","gray44","gray45","gray46","gray47","gray48","gray49","gray50",
                "gray51","gray52","gray53","gray54","gray55","gray56","gray57","gray58","gray59","gray60",
                "gray61","gray62","gray63","gray64","gray65","gray66","gray67","gray68","gray69","gray70",
                "gray71","gray72","gray73","gray74","gray75","gray76","gray77","gray78","gray79","gray80",
                "gray81","gray82","gray83","gray84","gray85","gray86","gray87","gray88","gray89","gray90",
                "gray91","gray92","gray93","gray94","gray95","gray96","gray97","gray98","gray99")
    allcolors = colors()
    colorcoder=setdiff(allcolors, c(colorcodex, "grey") )
    colorcode = c(colorcodex, colorcoder)
    }

    # ?grey" means not in any module;
    colorhelp=rep("grey",length(labelpred))
    if ( no.modules==0){
        print("No mudule detected\n")
    }
    else{
        if ( no.modules > length(colorcode)  ){
            print( paste("Too many modules \n", as.character(no.modules)) )
        }

        if ( (nameallmodules==FALSE) || (no.modules <=length(colorcode)) ){
            labeledModules = min(no.modules, length(colorcode) )
            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)
            }
            if(!useNumberAsLabel){
               colorhelp=factor(colorhelp,levels=c(colorcode[1:labeledModules],"grey"))
            }
        }else{#nameallmodules==TRUE and no.modules >length(colorcode)
            maxcolors=length(colorcode)
            labeledModules = no.modules
            extracolors=NULL
            blackwhite=c("black", "white")
            for(i in c((maxcolors+1):no.modules)){
              if(useblackwhite==FALSE){
                  icolor=paste("module", as.character(i), sep="")
              }else{#use balck white alternatively represent extra colors, for display only
                  icolor=blackwhite[1+(i%%2) ]
              }
              extracolors=c(extracolors, icolor)
            }
            allcolorcode=c(colorcode, extracolors)
            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],allcolorcode[i],colorhelp)
            }
            
            if(!useNumberAsLabel){
               colorhelp=factor(colorhelp,levels=c(allcolorcode[1:labeledModules],"grey"))
            }
        }
    }
    colorhelp
}

# make label by patching zeros
patchZeros = function(intnumber, digits=4){
  strnumber = as.character(intnumber)
  splitted  = splitString(strnumber,"")
  mydigits  = length(splitted) 

  patchedDigits = digits - mydigits
  if (patchedDigits<=0){
     return(strnumber)
  }
  myzeros=c("0","00","000","0000", "00000","000000","0000000", "00000000", "000000000")

  patched = paste(myzeros[patchedDigits], strnumber, sep="")

  return (patched)
}




#"0" is for grey module
assignModuleColor = function(labelpred, minsize1=50, anameallmodules=FALSE, auseblackwhite=FALSE, useNumberAsLabel=F, startlabel=0) {
    # here we define modules by using a height cut-off for the branches
    #labelpred= cutree2(hiercluster,h=heightcutoff)
    #cat(labelpred)

    #"0", grey module doesn't participate color assignment, directly assigned as "grey"
    labelpredNoZero = labelpred[ labelpred >0 ]
    sort1=-sort(-table(labelpredNoZero))
    sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 > minsize1
    no.modules=sum(modulebranch)

    if (useNumberAsLabel){
       # now make cluster label
       #
       colorcode=NULL
       for (i in c(startlabel:(startlabel+no.modules-1)) ){
          ipad = patchZeros(i)
          colorcode= c(colorcode, ipad)
       }
    }else{
    # now we assume that there are fewer than 10 modules
    #colorcode=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow")
    # now we assume that there are fewer than 10 modules
    colorcodex=c("turquoise",    "blue",     "brown",   "yellow",      "green",      "red",     "black",
                "pink",         "magenta",  "purple",  "greenyellow", "tan",        "salmon",  "cyan", 
                "midnightblue", "lightcyan","grey60",  "lightgreen",  "lightyellow","coral",   "sienna",
                "gold",         "peru",     "wheat",   "chocolate",   "seashell",   "khaki",   "bisque",
                "forestgreen",  "navy",     "plum",    "mediumblue",  "violet",     "hotpink",
                "thistle",      "orchid",   "maroon",  "violetred",   "firebrick",  "honeydew","chartreuse",
                "deeppink",     "darkcyan", "beige",   "snow",        "burlywood",  "goldenrod",
                "brown2",       "red2",     "gold2",   "yellow2",     "green2",     "cyan2",    "blue2",
                "brown3",       "red3",     "gold3",   "yellow3",     "green3",     "cyan3",    "blue3",
                "brown4",       "red4",     "gold4",   "yellow4",     "green4",     "cyan4",    "blue4",
                "gray1","gray2","gray3","gray4","gray5","gray6","gray7","gray8","gray9","gray10",
                "gray11","gray12","gray13","gray14","gray15","gray16","gray17","gray18","gray19","gray20",
                "gray21","gray22","gray23","gray24","gray25","gray26","gray27","gray28","gray29","gray30",
                "gray31","gray32","gray33","gray34","gray35","gray36","gray37","gray38","gray39","gray40",
                "gray41","gray42","gray43","gray44","gray45","gray46","gray47","gray48","gray49","gray50",
                "gray51","gray52","gray53","gray54","gray55","gray56","gray57","gray58","gray59","gray60",
                "gray61","gray62","gray63","gray64","gray65","gray66","gray67","gray68","gray69","gray70",
                "gray71","gray72","gray73","gray74","gray75","gray76","gray77","gray78","gray79","gray80",
                "gray81","gray82","gray83","gray84","gray85","gray86","gray87","gray88","gray89","gray90",
                "gray91","gray92","gray93","gray94","gray95","gray96","gray97","gray98","gray99")
    allcolors = colors()
    colorcoder=setdiff(allcolors, c(colorcodex,"grey"))
    colorcode = c(colorcodex, colorcoder)

    }

    #"grey" means not in any module;
    colorhelp=rep("grey",length(labelpred))
    if ( no.modules==0){
        print("No mudule detected\n")
    } else{
        if ( no.modules > length(colorcode)  ){
            print( paste("Too many modules \n", as.character(no.modules)) )
        }

        if ( (anameallmodules==FALSE) | (no.modules <=length(colorcode)) ){
            labeledModules = min(no.modules, length(colorcode) )
            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)
            }
            if(!useNumberAsLabel){
               colorhelp=factor(colorhelp,levels=c(colorcode[1:labeledModules],"grey"))
            }

        }else{#nameallmodules==TRUE and no.modules >length(colorcode)
            maxcolors=length(colorcode)
            labeledModules = no.modules
            extracolors=NULL
            blackwhite=c("red", "black")
            for(i in c((maxcolors+1):no.modules)){
              if(auseblackwhite==FALSE){
                  icolor=paste("module", as.character(i), sep="")
              }else{#use balck white alternatively represent extra colors, for display only
                  #here we use the ordered label to avoid put the same color for two neighboring clusters
                  icolor=blackwhite[1+(as.integer(modulename[i])%%2) ]
              }
              extracolors=c(extracolors, icolor)
            }

            #combine the true-color code and the extra colorcode into a uniform colorcode for 
            #color assignment
            allcolorcode=c(colorcode, extracolors)

            for (i in c(1:labeledModules)) {
               colorhelp=ifelse(labelpred==modulename[i],allcolorcode[i],colorhelp)
            }

            if(!useNumberAsLabel){
               colorhelp=factor(colorhelp,levels=c(allcolorcode[1:labeledModules],"grey"))
            }
        }
    }

    colorhelp
}


# label each grey node as a single module until all colors used up
#
turnGreyNodesIntoModules = function(incolorcode) {

    # now we assume that there are fewer than 10 modules
    #colorcode=c("turquoise","blue","brown","yellow","green","red","black","pink","magenta","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow")
    # now we assume that there are fewer than 10 modules
    colorcodex=c("turquoise",    "blue",     "brown",   "yellow",      "green",      "red",     "black",
                "pink",         "magenta",  "purple",  "greenyellow", "tan",        "salmon",  "cyan", 
                "midnightblue", "lightcyan","grey60",  "lightgreen",  "lightyellow","coral",   "sienna",
                "gold",         "peru",     "wheat",   "chocolate",   "seashell",   "khaki",   "bisque",
                "forestgreen",  "navy",     "plum",    "mediumblue",  "violet",     "hotpink",
                "thistle",      "orchid",   "maroon",  "violetred",   "firebrick",  "honeydew","chartreuse",
                "deeppink",     "darkcyan", "beige",   "snow",        "burlywood",  "goldenrod",
                "brown2",       "red2",     "gold2",   "yellow2",     "green2",     "cyan2",    "blue2",
                "brown3",       "red3",     "gold3",   "yellow3",     "green3",     "cyan3",    "blue3",
                "brown4",       "red4",     "gold4",   "yellow4",     "green4",     "cyan4",    "blue4",
                "gray1","gray2","gray3","gray4","gray5","gray6","gray7","gray8","gray9","gray10",
                "gray11","gray12","gray13","gray14","gray15","gray16","gray17","gray18","gray19","gray20",
                "gray21","gray22","gray23","gray24","gray25","gray26","gray27","gray28","gray29","gray30",
                "gray31","gray32","gray33","gray34","gray35","gray36","gray37","gray38","gray39","gray40",
                "gray41","gray42","gray43","gray44","gray45","gray46","gray47","gray48","gray49","gray50",
                "gray51","gray52","gray53","gray54","gray55","gray56","gray57","gray58","gray59","gray60",
                "gray61","gray62","gray63","gray64","gray65","gray66","gray67","gray68","gray69","gray70",
                "gray71","gray72","gray73","gray74","gray75","gray76","gray77","gray78","gray79","gray80",
                "gray81","gray82","gray83","gray84","gray85","gray86","gray87","gray88","gray89","gray90",
                "gray91","gray92","gray93","gray94","gray95","gray96","gray97","gray98","gray99")
    allcolors = colors()
    #colorcoder=setdiff(allcolors, colorcodex)
    colorcoder=setdiff(allcolors, c(colorcodex,"grey"))
    colorcode = c(colorcodex, colorcoder)

    maxcolors = length(colorcode)
    colorcodeIndex = cbind(colorcode, c(1:maxcolors)) 
    
    cincolorcode   = as.character(incolorcode)

    # find the maximum index of the last color used
    #
    merged = merge(incolorcode, colorcodeIndex, by.x=1, by.y=1, all=F)
    maxColIdx = max(as.integer(as.character(merged[,2])))


    # turn each grey into a new color
    #
    no.items = length(cincolorcode)
    selgrey  = cincolorcode=="grey"
    if(sum(selgrey)==0 | maxColIdx>maxcolors){
       return (incolorcode)
    }

    greyIdx  = c(1:no.items)[selgrey]
    no.greys = length(greyIdx)
    
    curcolorIdx = maxColIdx + 1
    for (g in c(1:no.greys)){
       cincolorcode[ greyIdx[g] ] = colorcode[curcolorIdx]
       curcolorIdx = curcolorIdx + 1
       if(curcolorIdx>maxcolors){
           break;
       }
    }

    colorhelp=factor(cincolorcode)

    return (colorhelp)
}



#merge the minor cluster into the major cluster
merge2Clusters = function(mycolorcode, mainclusterColor, minorclusterColor){
  mycolorcode2 = mycolorcode
  for(mc in minorclusterColor) {
     mycolorcode2 = ifelse(as.character(mycolorcode2)==mc, mainclusterColor, as.character(mycolorcode2) )
  }

  fcolorcode   =factor(mycolorcode2)
  fcolorcode
}


mergeMultiClusters = function(mycolorcode, clustersForMergeList, mdendro){
  mycolorcode2 = mycolorcode
  no.lists = length( clustersForMergeList)
  for(x in c(1:no.lists)) {
     xclusters = clustersForMergeList[[x]]
     mainclusterColor = xclusters[1]
     minorclusterColor= xclusters[-1]
     for(mc in minorclusterColor) {
        mycolorcode2 = ifelse(as.character(mycolorcode2)==mc, mainclusterColor, as.character(mycolorcode2) )
     } 
  }
  tbAll = table(mycolorcode); minsize = min(tbAll)
 
  ngenes       = length(mycolorcode2)
  orderMtrx    = cbind(c(1:ngenes), mdendro$order); colnames(orderMtrx)=c("DenfrogramOrder", "OrignalOrder")
  colorOrdered = as.character(mycolorcode2[mdendro$order])
   
  fdispcolor     = getDisplayColorSequence(colorOrdered)
  labeledColors  = fdispcolor[fdispcolor!=""]
  labeledColXcoor= c(1:length(fdispcolor))[fdispcolor!=""]
  tb             = table(labeledColors)
  if(sum(tb>1)==0) {
    brokenColors = names(tb)
    mycolorcodeX = mycolorcode2
  } else {
    replColors = names(tb)[tb>1]
    for(ec in replColors) {
       idx = c(1:ngenes)[colorOrdered ==ec];
       colorOrdered[min(idx):max(idx)] = ec
    }
    od = order(orderMtrx[,2])
    orderMtrx2= orderMtrx[od, ]
    mycolorcodeX = colorOrdered[ orderMtrx2[,1] ]
  }

  mycolorcode3 = reassignModuleNames(mycolorcodeX, minmodulesize=minsize-1, anameallmodules=FALSE, auseblackwhite=FALSE,
useNumberAsLabel=FALSE, startlabel=0)

  fcolorcode   =factor(mycolorcode3)
  fcolorcode
}
        


showOrderedModuleColors = function(mycolorcode, mdendro){
   colorOrdered = as.character(mycolorcode[mdendro$order])
  fdispcolor     = getDisplayColorSequence(colorOrdered)
  labeledColors  = fdispcolor[fdispcolor!=""]
  labeledColors
} 
 


# unnamed modules are represented by like "module12", thus these module names include keyword "module"
# > unlist(strsplit("module21","module"))
# [1] ""   "21"
# > unlist(strsplit("yellow","module"))
# [1] "yellow"
assignPseudoColors = function(mdendro, incolorcode){

  #locate ordered locations of clusters in the dendrogram
  fdispcolor = getDisplayColorSequence(as.character(incolorcode[mdendro$order]))
  labeledColors= fdispcolor[fdispcolor!=""]

  mycolorcode= as.character(incolorcode)

  pseudocolors=c("white","black")
  cnt.unnamedcolor = 0
  for ( each in labeledColors){
    is.unnamedcolor = unlist(strsplit(each,"module"))	
    if (length(is.unnamedcolor) == 1){ #the module is labeled with true color, and keep it
        next
    }
    
    # the module is labeled with true color, so we need replace it with 
    if (cnt.unnamedcolor==1){
        cnt.unnamedcolor =2
    }else{
        cnt.unnamedcolor =1
    }

    ipseduocol = pseudocolors[cnt.unnamedcolor]

    #cat(each, " to ", ipseduocol, "\n")

    mycolorcode = ifelse(mycolorcode==each, ipseduocol, mycolorcode)
  }

  fcolorcode   =factor(mycolorcode)
  fcolorcode
}


exchange2Clusters = function(mycolorcode, colorA, colorB){
  colorC   = "great"
  fcolcode = mycolorcode
  fcolcode = merge2Clusters(fcolcode, mainclusterColor=colorC, minorclusterColor=colorA) #A->C
  fcolcode = merge2Clusters(fcolcode, mainclusterColor=colorA, minorclusterColor=colorB) #B->A
  fcolcode = merge2Clusters(fcolcode, mainclusterColor=colorB, minorclusterColor=colorC) #C->B
  fcolorcode   =factor(fcolcode)
  fcolorcode
}

# return the name of the max element from a set of element names
# 
getNameOfMaxElement=function(mytable, selectednames){
   maxsize=0
   maxname=""
   for (each in selectednames){
       if (mytable[[each]] >maxsize){
           maxsize = mytable[[each]]
           maxname = each
       }
   }
   return (maxname)
}


# cut dendrogram to get the required number of clusters
#
moduleDetectByFixedClusterno = function(ihcluster, nocluster=2) {
   maxhei   = max(ihcluster$height)
   minhei   = min(ihcluster$height)
   curheightcutoff = minhei
   curno = 0
   cnt=0
   while(curno != nocluster){
       cnt = cnt + 1
       if (cnt > 6)
          break
       curheightcutoff = (curheightcutoff + maxhei)/2.0
       colcode  = moduleDetectLabel(hiercluster=ihcluster, curheightcutoff, minsize1=1)
       colcode  = as.character(colcode)
       # we need as.character(colcode) to deal with the case of grey is zero
       #table(colcolors)
       #turquoise      blue      grey 
       #12             6         0 
       colr = as.integer(table(colcode))
       curno    = sum(colr>0)
   }
   colcode   
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##                           Get Most Varying Genes  
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get nvars most varying genes, consider only genes with less than 20% missing values
## default use of standard deviation
## other option: coefficient of variation (CV), i.e., the standard deviation divided by mean
getMostVaryingGenes=function(myinputfname, nvars, myHeadCol=8, missingRate=0.2, mysep="\t", vartype=0){
    
    mfname      =getFileName(myinputfname)
    mvgOutFname =paste(mfname, "_", nvars, "mvgenes.xls", sep='')

    allMatrix <- read.delim(myinputfname,sep=mysep, header=T)
    dims=dim(allMatrix)

    genesInfor <- allMatrix[,c(1:myHeadCol)]
    rowTitles=names(allMatrix)
    rowTitles=rowTitles[1:dims[2]]

    dat1 <-t(allMatrix[,-c(1:myHeadCol)])

    ## --------- compute genes with missing values
    naMatrix <- is.na(dat1)
    nasumVector  <- apply(naMatrix,2,sum, na.rm=TRUE)
    naprobVector <- nasumVector/dim(dat1)[[1]]

    dat1 <- dat1[,naprobVector<missingRate]
    dim(dat1)
    if(myHeadCol>1){
      genesInfor_NAclean <- genesInfor[naprobVector<missingRate,]
    }else{#geneinfor is a vector
      genesInfor_NAclean <- genesInfor[naprobVector<missingRate]
    }

    var1 <- apply(dat1,2,sd,   na.rm=TRUE)
    means <- apply(dat1,2,mean, na.rm=TRUE)

    if (vartype !=0 ){
       var1 <- var1/means        
    }

    rk1  <- rank(-var1)
    dat2 <- dat1[,rk1<nvars+1]

    noGenes   <- dim(dat2)[2]
    noSamples <- dim(dat2)[1]

    if(myHeadCol>1){
      filteredGenesInfor <- genesInfor_NAclean[rk1<nvars+1, ]
    }else{
      filteredGenesInfor <- genesInfor_NAclean[rk1<nvars+1]
    }

    newCol = length(rowTitles)
    yfile <- file(mvgOutFname, "w+")
    #first row
    for (z in 1:newCol){
       cat(as.character(rowTitles[z]), file=yfile)
      if (z==newCol){
         cat("\n", file=yfile)
      }
      else{
         cat("\t", file=yfile)
      }
    }

    for (z in 1:noGenes){
      if(myHeadCol>1){
          for (i in 1:dim(filteredGenesInfor)[2] ){
             cat(as.character(filteredGenesInfor[z,i]), file=yfile)
             cat("\t", file=yfile)
          }
      }else{
          cat(as.character(filteredGenesInfor[z]), file=yfile)
          cat("\t", file=yfile)
      }

      for (i in 1:noSamples){
         cat(as.character(dat2[i,z]), file=yfile)
         if (i==noSamples){
            cat("\n", file=yfile)
         }
         else{
            cat("\t", file=yfile)		
         }
      }
    }
    close(yfile)

    mvgOutFname
}

# trim the values away from 4 standard deviations
Ztransform = function(datvect, cutoff=4)
{
 mymean = mean(datvect, na.rm=TRUE)
 mysd = sd(datvect, na.rm=TRUE)
 zscores = (datvect - mymean)/mysd
 if(!is.na(cutoff) ) { 
   zscores2 = ifelse(zscores > cutoff, cutoff, ifelse(zscores < -cutoff, -cutoff, zscores) )
   return (cbind(zscores2))
 } else {
   return (cbind(zscores));
 }
}


#----------------------------------------------------------------------------------------------------
#write a table with row and colnames into Excel or Latex output
#Input
#     myMatrix: table 
#     fname:    output filename
#     latex:    output Latex(T) or Excel file(F)
#     myappend: append to the file or not(T/F)
#     usesignif: use short repsentation of floating numbers (T/F)
#----------------------------------------------------------------------------------------------------
writeTableIntoLatex=function(myMatrix, fname, myappend=F, latex=T, mycol.names=T, myrow.names=T, usesignif=F){
   
   mycolnames = colnames(myMatrix)
   myrownames = rownames(myMatrix)

   #put colnames at the first row and rownames as the 1st column, leave a blank for the cell [1,1]
   latexMatrix1 = rbind(mycolnames, as.matrix(myMatrix) )
   latexMatrix  = cbind( c(" ", myrownames), as.matrix(latexMatrix1))
   no.rows=dim(latexMatrix)[1]
   no.cols=dim(latexMatrix)[2]

   if(latex==FALSE){
      if(usesignif==TRUE){
         for (i in c(2:no.rows) ){
           for (j in c(2:no.cols) ){
              latexMatrix[i, j] = signif(as.numeric(latexMatrix[i,j]),2)
           } 
         }
      }
      write.table(latexMatrix,fname, append=myappend, sep="\t", quote=FALSE, col.names=F, row.names=F)
      return (1)
   }

   for (i in c(1:no.rows) ){
     #add & to each cell except those in the last column
     for (j in c(1:(no.cols-1)) ){
         if((usesignif==TRUE)&(i>1)&(j>1) ){
           latexMatrix[i, j] = paste(signif(as.numeric(latexMatrix[i,j]),2), "&", sep="")
         }else{
           latexMatrix[i, j] = paste(latexMatrix[i,j], "&", sep="")
         }
     }

     #\\\\\\hline to each row
     if((usesignif==TRUE)&(i>1) ){
         latexMatrix[i, no.cols] = paste(signif(as.numeric(latexMatrix[i, no.cols]),2), "\\\\\\hline", sep="")
     }else{
         latexMatrix[i, no.cols] = paste(latexMatrix[i, no.cols], "\\\\\\hline", sep="")
     }
   }

   write.table(latexMatrix,fname, append=myappend, sep="\t", quote=FALSE, col.names=F, row.names=F)
}


compareTwoModulesOldOld = function(gifA, gifB, moduleNameA, moduleNameB, uniqueIdCol=1, moduleColA, moduleColB, totalGenes, removeDuplicate=TRUE)
{
  if(! removeDuplicate) {
    restrictA = gifA[, moduleColA]== moduleNameA 
    restrictB = gifB[, moduleColB]== moduleNameB
    noA= sum(restrictA)
    noB= sum(restrictB)

    moduleSetA = gifA[restrictA, ]
    moduleSetB = gifB[restrictB, ]

    ovlp = intersect(moduleSetA[,uniqueIdCol], moduleSetB[,uniqueIdCol])
    #mergeMatrix = merge(moduleSetA, moduleSetB, by.x=uniqueIdCol, by.y=uniqueIdCol, sort=F,all=FALSE)
    #intersectNo = dim(mergeMatrix)[1]
  } else {
    
    selA = gifA[, moduleColA]== moduleNameA
    restrictA = unique(gifA[selA, uniqueIdCol])

    selB = gifB[, moduleColB]== moduleNameB 
    restrictB = unique(gifB[selB, uniqueIdCol])
    
    noA= length(restrictA)
    noB= length(restrictB)

    ovlp = intersect(restrictA, restrictB)
  }

  intersectNo = length(ovlp)
  #phyper(89,702,4000-702,280, lower.tail=F)
  pval = phyper(intersectNo-1, noA, totalGenes-noA, noB, lower.tail=F)

  ret=c(intersectNo, pval, noA, noB)
  ret
}

compareTwoModulesOld = function(gifA, gifB, moduleNameA, moduleNameB, uniqueIdCol=1, moduleColA, moduleColB, totalGenes)
{
  restrictA = gifA[, moduleColA]== moduleNameA 
  restrictB = gifB[, moduleColB]== moduleNameB
  noA= sum(restrictA)
  noB= sum(restrictB)

  moduleSetA = gifA[restrictA, ]
  moduleSetB = gifB[restrictB, ]

  mergeMatrix = merge(moduleSetA, moduleSetB, by.x=uniqueIdCol, by.y=uniqueIdCol, sort=F,all=FALSE)
  intersectNo = dim(mergeMatrix)[1]

  #phyper(89,702,4000-702,280, lower.tail=F)
  pval = phyper(intersectNo-1, noA, totalGenes-noA, noB, lower.tail=F)

  ret=c(intersectNo, pval, noA, noB)
  ret
}

overlapTwoLists = function(listA, listB, totalGenes, retentries=F)
{
  noA = length(listA)
  noB = length(listB)
  overlap = intersect(listA, listB); #merge(listA, listB, by.x=1, by.y=1, sort=F,all=FALSE)
  intersectNo = length(overlap)
  
  if (intersectNo >0){
    #phyper(89,702,4000-702,280, lower.tail=F)
    pval = phyper(intersectNo-1, noA, totalGenes-noA, noB, lower.tail=F)
    commongenes = concatenate(sort(overlap),";")
  }else{
    pval = 1
    commongenes = ""
  }

  if (retentries){
    ret=c(intersectNo, pval, commongenes)
  }else{
    ret=c(intersectNo, pval)
  }
  ret
}



#************************************
#* Given d and r, we compute s, s= d/(sin(pi/2 - atan(r))
#  /|
#/  |  /
#   |/
#  /|\ d  /
#/ s|  \/
#   | /
#   /)r
# /-------------
loessfitbounds=function(y, x){
     myseq = as.numeric(seq(min(x),to=max(x), length=20 )) # length=length(x)) )
     length(myseq)
     new <- data.frame( myseq )

     corrl= cor(x,y)
     lm1=lm(y ~ x)
     r= lm1$coefficients[2]
     d=1.0
     if (r >= 5*pi/180){
       s= d/( sin(pi/2 - atan(r)) )
     }else{
       s= d
     }

     loe1=loess(y ~ x, span=1, degree=1)
     predloe1=predict(loe1, new, se=T)
     myfit= predloe1$fit

     #return x_fitted, y_fitted, mean(se)
     #list(myseq, myfit, mean(predloe1$se))
     list(myseq, myfit, predloe1$se)
     #list(myseq, myfit, s)
}

precision=function(x, digits=3){
   if(abs(x)>1){
       xx=round(x, digits)
   } else{
       xx=signif(x, digits)
   }
   return(xx)
}

# textPosY as the portion of the range Y to be shifted for text display
#
#
scatterPlot=function(aa.pvector, bb.pvector, myxlab="", myylab="", imgfilename=NULL, 
                     pointcolor=NULL, textlabels=NULL, textcolor=NULL, shortnote="", 
                     useloess=F, splinefit=F,powerfit=F, lmline=TRUE, imgwidth=400, imgheight=400, 
                     icex=0.6, pointsize = 12, bcolor="blue",
                     textPosY=0, ret_title=T, cormethod="pearson", plotSquare=F, pretitle="")
{   
   lm1=lm(bb.pvector ~ aa.pvector)

   a.min= min(aa.pvector); a.max=max(aa.pvector)

   xv = seq(a.min, a.max, (a.max-a.min)/length(aa.pvector))

   if(splinefit) {
       fit.spl <- smooth.spline(aa.pvector, bb.pvector, df=6.4)   
   }

   if(useloess) { #return x_fitted, y_fitted, mean(se)
      fitted=loessfitbounds(y=bb.pvector, x=aa.pvector)
   }

   no.points = length(aa.pvector)

   nonnas = !(is.na(aa.pvector) | is.na(bb.pvector)) 

   if( length(nonnas) >2 ) {

    if(cormethod=="spearman") {
      corrl= cor(rank(aa.pvector), rank(bb.pvector))
      corrl.p = signif( cor.test(rank(aa.pvector), rank(bb.pvector), method="p",use="p")$p.value ,2)
    } else{
      corrl= cor(aa.pvector, bb.pvector)
      corrl.p = signif( cor.test(aa.pvector, bb.pvector, method="p",use="p")$p.value ,2)
    }

    if(is.na(corrl.p)) {
       corrp.string = "=NA"
       mytitle=paste(shortnote, "r=", as.character(signif(corrl, 2))," (p", corrp.string, ")", sep="")
    } else if (corrl.p>0){
       corrp.string = as.character(format(corrl.p, scientific=TRUE))
       mytitle=paste(shortnote, "r=", as.character(signif(corrl, 2))," (p=", corrp.string, ")", sep="")
    }else{
       corrp.string = " < e-22"
       mytitle=paste(shortnote, "r=", as.character(signif(corrl, 2))," (p", corrp.string, ")", sep="")
    }

   } else {
       mytitle=shortnote
   }

   pmodel  = ""
   if(powerfit){ # log(y) ~const + power*x
     bb.log = log10(bb.pvector)
     lm3    = lm(bb.log ~ aa.pvector);
     slm3   = summary(lm3)
     const  = slm3$coeff[1,1]
     power  = slm3$coeff[2,1]
     bb.xv  = 10^(const + xv*power)
     pmodel = paste("y=", precision(10^const,3), "*10^(",precision(power,3), "x)",sep="")
     mytitle= paste(mytitle, ", model fitting: ", pmodel, sep="")
   }

   #plot(aa.pvector, bb.pvector, main=mytitle,xlab=myxlab, ylab=myylab )

   mse=2

   #ostring=paste(as.character(signif(corrl,2)), "\t", as.character(corrl.p),"\n" )
   #cat(ostring, file=fp)

   if (is.null(pointcolor) ){
        pcolors=rep("black", no.points)
   }else{
        pcolors=pointcolor
   }
   if (is.null(textcolor) ){
        tcolors=rep("black", no.points)
   }else{
        tcolors=textcolor
   }

   if(plotSquare) {
      allvect = c(aa.pvector, bb.pvector)
      mxlim = c(min(allvect), max(allvect)); mylim = mxlim;
      mytitle=pretitle
   } else {
      mxlim  = c(min(aa.pvector), max(aa.pvector)); mylim = c(min(bb.pvector), max(bb.pvector));
      mytitle= paste(pretitle, mytitle, sep=" ")
   }


   if(!is.null(imgfilename)){
     openImgDev(imgfilename, iwidth =imgwidth, iheight =imgheight, ipointsize = pointsize)
     #openImgDev(imgfilename)
     par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)

     plot(aa.pvector, bb.pvector, main=mytitle,xlab=myxlab, ylab=myylab, col=pcolors, xlim=mxlim, ylim=mylim)

     if(useloess){
       lines(x=fitted[[1]], y=fitted[[2]], col=bcolor)
       lines(x=fitted[[1]], y=fitted[[2]]+mse*fitted[[3]], col=bcolor,lty=2)
       lines(x=fitted[[1]], y=fitted[[2]]-mse*fitted[[3]], col=bcolor,lty=2)
     }else{
       if(plotSquare) {
          lines(mxlim, mylim, col=bcolor)
       } else{
          if(lmline) {
            abline(lm1, col=bcolor)
          }
       }
     }

     if(splinefit) {
       lines(pp <- predict(fit.spl, xv), col = "red")
     }

     if(powerfit){ # log(y) ~const + power*x
        lines(xv, bb.xv, col=2)
     }

     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
     if (!is.null(textlabels) ){
        adj=( max(bb.pvector)-min(bb.pvector))*textPosY
        text(aa.pvector, bb.pvector+adj, labels = textlabels, col=tcolors, cex=icex)
     }
     dev.off()
   }

     #par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
     plot(aa.pvector, bb.pvector, main=mytitle,xlab=myxlab, ylab=myylab, col=pcolors, xlim=mxlim, ylim=mylim)
     

     #par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
     if (!is.null(textlabels) ){
        #yshift= rep((max(bb.pvector)-min(bb.pvector))*textPosY, length(bb.pvector) )
        adj=( max(bb.pvector)-min(bb.pvector))*textPosY
        text(aa.pvector, bb.pvector+adj, labels = textlabels, col=tcolors, cex=icex)
     }

     if(useloess){
       lines(x=fitted[[1]], y=fitted[[2]], col=bcolor)
       lines(x=fitted[[1]], y=fitted[[2]]+mse*fitted[[3]], col=bcolor,lty=2)
       lines(x=fitted[[1]], y=fitted[[2]]-mse*fitted[[3]], col=bcolor,lty=2)
     }else{
       if(plotSquare) {
          lines(mxlim, mylim, col=bcolor)
       } else{
          if(lmline) {
             abline(lm1, col=bcolor)
          }
       }

     }

     if(splinefit) {
       lines(pp <- predict(fit.spl, xv), col = "red")
     }

     if(powerfit){ # log(y) ~const + power*x
       lines(xv, bb.xv, col=2)
     }

   if(ret_title) {
      return (mytitle)
   } else{
      return ( c(corrl, corrl.p) )
   }
}

# here we identify outliers with difference of fitted value and true value is 
#   equal to or bigger than mnoSDsAsOutlier*sd(fitted differences)
identify_linearfit_outliers=function(yy, xx, mnoSDsAsOutlier=2){
   lm1=lm(yy~xx)
   #plot(a, c)
   #abline(lm1)
   fitted= as.numeric(lm1$fitted.values)
   diffs = abs(fitted-c)
   mean.diff =  mean(diffs)
   sd.diff   =  sd(diffs)
   threshold = mean.diff + sd.diff*mnoSDsAsOutlier

   goods  = (diffs < threshold) # points with fitted value less than mean + n*standard deviation
   goods
}

# we consider an option of removing ouliers
singlePlotComplex = function(outputImage=NULL, rawValA, rawValB, plotTitle="", myylab="",myxlab="", dataSet="", userawdata=F, numPoints=20, removeOutliers=F, noSDsAsOutlier=2){

   if(!userawdata){
	#in general, we will have the thing being tested be rawValA, and K be rawValB
	cut_factor = cut(rawValB,numPoints)
	# Now give me the mean proportion of essentiality in each chunk of the essentiality vector 
	ValA = tapply(rawValA, cut_factor, mean)
	# Now give me the mean proportion of k in each chunk of the essentiality vector 
	ValB = tapply(rawValB, cut_factor, mean)
        
        ValA = ValA[!is.na(ValA)]
        ValB = ValB[!is.na(ValB)]

   }else{
        ValA = rawValA
        ValB = rawValB
   }


   if(removeOutliers){
      good.points = identify_linearfit_outliers(ValA, ValB, mnoSDsAsOutlier=noSDsAsOutlier)
      #print(length(good.points) )
      #print(sum(good.points) )

      ValA= ValA[good.points]
      ValB= ValB[good.points]
   }
        #print(ValA)
        #print(ValB)


	#Now we need to just make a line fit
	line =lm(ValA ~ ValB)

	corr_p   = signif( cor(ValB,ValA),2)	#pearson 
        corr_p.p = signif( cor.test(ValB,ValA, method="p",use="p")$p.value ,2)

	corr_s   = signif(cor(ValB,ValA, method = "s"),2)	#spearman 
        corr_s.p = signif( cor.test(ValB,ValA, method="s",use="p")$p.value ,2)

	titleinfor = paste(as.character(dataSet),",",as.character(plotTitle), ",", "p_r=", as.character(corr_p), ",p_p=", as.character(corr_p.p)
	,"(s_r=", as.character(corr_s),",s_p=", as.character(corr_s.p),")",sep="" )

	# do a plot() of these two vectors against each other. 
	if (!is.null(outputImage) )
	    openImgDev(outputImage)

        plot(ValB,ValA,main=titleinfor, xlab=myxlab, ylab=myylab)
     	abline(line, col="red")

	# and then draw the line
	if (!is.null(outputImage) )
             dev.off()
        #titleinfor
    
    return ( c(corr_p, corr_p.p, corr_s, corr_s.p) )
}



#save the image into file (Dendrogram, Cluster Color-bars, and Color Names together)
plotDendrogramModuleLabels=function(mdendro, modulecolors, save2file=NULL, plotLabels=FALSE, secondColorbar=NULL, secondLabel=NULL, 
      gramcolor="black", plotSampleNames=F, pwid=1024,phei=600, acex=0.6)
{
   fdispcolor = getDisplayColorSequence(as.character(modulecolors[mdendro$order]))
   labeledColors= fdispcolor[fdispcolor!=""]
   labeledColXcoor= c(1:length(fdispcolor))[fdispcolor!=""]

   pseudocolor = assignPseudoColors(mdendro=mdendro, incolorcode=modulecolors)

   if(!is.null(save2file)){ 
      openImgDev(save2file, iwidth = pwid, iheight = phei)
   }

   if(!is.null(secondColorbar) & plotLabels) {
     par(mfrow=c(3,1), mar=c(5,0,0,0), cex=acex)
   } else if(plotLabels){ #show color labels
     par(mfrow=c(2,1), mar=c(10, 0, 0, 0), cex=acex)
   } else if(!is.null(secondColorbar) ) {
     par(mfrow=c(3,1), mar=c(0,0,0,0) )
   } else{
     par(mfrow=c(2,1), mar=c(0,0,0,0) )
   } 

   if (plotSampleNames){
     plot(mdendro, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor)
   }else{
     plot(mdendro, labels=F, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor)
   }

   #plot(c(1:length(modulecolors)), rep(1, length(modulecolors)), ylim=c(0,1),xlab="", main="",
   #          ylab="",  type = "h",col=as.character(modulecolors[mdendro$order]), las=2,
   #          lwd=1, axes=F, frame=F)

   plot(c(1:length(modulecolors)), rep(1, length(pseudocolor)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=as.character(pseudocolor[mdendro$order]), las=2,
             lwd=8, axes=F, frame=F)


   if(plotLabels){ #show color labels
      axis(1, at =labeledColXcoor, labels = labeledColors, las=2)
   }
   #barplot(height=rep(1, length(modulecolors)), names.arg=fdispcolor, 
   #       col= as.character(modulecolors[mdendro$order]), space=0, las=2,
   #      border=F,main="", axes = T, axisnames = T)#, cex.lab=0.5)

   if(!is.null(secondColorbar) ) {
     colorlabel2 = as.character(secondColorbar[mdendro$order])
     plot(c(1:length(secondColorbar)), rep(1, length(secondColorbar)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=colorlabel2, las=2,
             lwd=8, axes=F, frame=F)     

     if(!is.null(secondLabel)) {
        labeledColXcoor2 = c(1:length(colorlabel2))[colorlabel2!=""]
        lablels2         = secondLabel[colorlabel2 !=""]
        axis(1, at =labeledColXcoor2, labels =lablels2 , las=2)
     }

   }

   par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

   if(!is.null(save2file)){ 
      dev.off()
   }
}



#save the image into file (Dendrogram, Cluster Color-bars, and Color Names together)
plotDendrogramModuleLabelsHubs=function(mdendro, modulecolors, save2file=NULL, plotLabels=FALSE, secondColorbar=NULL, secondLabel=NULL, 
               gramcolor="black", plotSampleNames=F, pwid=1024,phei=800)
{
   fdispcolor = getDisplayColorSequence(as.character(modulecolors[mdendro$order]))
   labeledColors= fdispcolor[fdispcolor!=""]
   labeledColXcoor= c(1:length(fdispcolor))[fdispcolor!=""]

   pseudocolor = assignPseudoColors(mdendro=mdendro, incolorcode=modulecolors)

   if(!is.null(save2file)){ 
      openImgDev(save2file, iwidth = pwid, iheight = phei)
   }

   if(!is.null(secondColorbar) & plotLabels) {
     par(mfrow=c(1,1), mar=c(10,0,0,0), cex=0.8)
   } else if(plotLabels){ #show color labels
     par(mfrow=c(2,1), mar=c(10, 0, 0, 0), cex=0.6)
   } else if(!is.null(secondColorbar) ) {
     par(mfrow=c(3,1), mar=c(0,0,0,0) )
   } else{
     par(mfrow=c(2,1), mar=c(0,0,0,0) )
   } 


   if (plotSampleNames){
     if(!is.null(secondLabel)) {
       plot(mdendro, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor, labels=secondLabel)
       #plclust(mdendro, hang=2, xlab="",ylab="",main="",sub="",axes = T, labels=secondLabel)
     }else{
       plot(mdendro, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor)
     }
   }else{
     plot(mdendro, labels=F, xlab="",ylab="",main="",sub="",axes = T, col=gramcolor)
   }


   if(!is.null(secondColorbar) ) {
     colorlabel2 = as.character(secondColorbar[mdendro$order])
     secondLabel2 = secondLabel[mdendro$order]
     #plot(c(1:length(secondColorbar)), rep(1, length(secondColorbar)), ylim=c(0,1),xlab="", main="",
     #        ylab="",  type = "h",col=colorlabel2, las=2,
     #        lwd=1, axes=F, frame=F)     

     if(!is.null(secondLabel)) {
        labeledColXcoor2 = c(1:length(colorlabel2))[colorlabel2!=""]
        lablels2         = secondLabel2[colorlabel2 !=""]
        axis(1, at =labeledColXcoor2, labels =lablels2 , las=2, 
               col = "violet", col.axis="dark violet", lwd = 1)
     }

   }

   if(1==2){
   #plot(c(1:length(modulecolors)), rep(1, length(modulecolors)), ylim=c(0,1),xlab="", main="",
   #          ylab="",  type = "h",col=as.character(modulecolors[mdendro$order]), las=2,
   #          lwd=1, axes=F, frame=F)

   plot(c(1:length(modulecolors)), rep(1, length(pseudocolor)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=as.character(pseudocolor[mdendro$order]), las=2,
             lwd=1, axes=F, frame=F)


   if(plotLabels){ #show color labels
      axis(1, at =labeledColXcoor, labels = labeledColors, las=2)
   }
   #barplot(height=rep(1, length(modulecolors)), names.arg=fdispcolor, 
   #       col= as.character(modulecolors[mdendro$order]), space=0, las=2,
   #      border=F,main="", axes = T, axisnames = T)#, cex.lab=0.5)

   }

   par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

   if(!is.null(save2file)){ 
      dev.off()
   }
}


#save the image into file (Dendrogram, Cluster Color-bars, and Color Names together)
plotDendrogramModuleMitosis=function(mdendro, modulecolors, addmitocolor, save2file=NULL)
{
   if(!is.null(save2file)){ 
      openImgDev(save2file)
   }
   par(mfrow=c(3,1), mar=c(0,0,0,0) )
   plot(mdendro, labels=F, xlab="",ylab="",main="",sub="",axes = T)
   plot(c(1:length(modulecolors)), rep(1, length(modulecolors)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=as.character(modulecolors[mdendro$order]), las=2,
             lwd=1, axes=F, frame=F)
   plot(c(1:length(addmitocolor)), rep(1, length(addmitocolor)), ylim=c(0,1),xlab="", main="",
             ylab="",  type = "h",col=as.character(addmitocolor[mdendro$order]), las=2,
             lwd=1, axes=F, frame=F)
   par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
   if(!is.null(save2file)){ 
      dev.off()
   }
}


#---------------------------------------------------------------------------------------------------------
#ModulePrincComps, which computes the first principal component of each module. 
#Also it provides the variance explained for the first 5 PCs of a module.
#It is based on the singular value decomposition.
#It takes as input datExpr  for which the rows are samples and the columns are genes.
#Here is how you would use it
#PC1=ModulePrinComps2[datExpr,color1]
#Then PC1[[1]] provides a data frame with the first PC of each module
#PC1[[2]] provides a data frame where each column contains the percentages of 
#the variance explained by the first 10 PCs by module.
#---------------------------------------------------------------------------------------------------------
# the output is a list and each element is a matrix whose number of columns equals the number of modules
# the first element of the list contains the 1st principal components of all modules, and 
# the second element of the list contains the 2nd principal components of all modules, etc...
# the last element is the expained variations of top "no.pcs" PCs in each modules
#---------------------------------------------------------------------------------------------------------

ModulePrinComps_alternative = function(datexpr,couleur, min_modulesize=10, no.pcs=10) {

  allmean = mean(datexpr, na.rm=T)

  modlevels= names(table(couleur)) #levels(factor(couleur))

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels

  for(i in c(1:length(modlevels)) ){
 
    print(paste(i, modlevels[i]) )
    modulename    = modlevels[i]
    restrict1= as.character(couleur)== modulename
    restrict1= ifelse(is.na(restrict1), FALSE, restrict1)

    if(modlevels[i]=="grey"){
       next
    }

    # in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])
    dim(datModule)

    # check whether some samples have missing rate (missing value in column) >80%
    naDat      = is.na(datModule)
    nasBySample= apply(naDat, 2, sum)
    nacutoff   = 0.6*dim(datModule)[1]
    selSamples = nasBySample>=nacutoff
    if(sum(selSamples)>0){
        print("patch samples with missing rate >=0.8")
        naSampleIdx = c(1:length(nasBySample))[selSamples]
        for(x in naSampleIdx){
            #print(paste("Sample Idx=", x) )
            datModule[,x]=ifelse(is.na( datModule[,x]), allmean, datModule[,x])
        }
    }
 
    if(sum(restrict1, na.rm=TRUE)<min_modulesize){
       listPCs[[1] ][,i] = datModule[1,]
       next
    }

    if(F) {
    imputed=impute.knn(as.matrix(datModule)) #$data for R version < 2.0
    datModule2 =imputed$data #=imputed for R < 2.9
    datModule2=t(scale(t(datModule2)))

    # use the hub's profuile as PC if svd doesn't converge
    #svd1=svd(datModule2)
    svd1=tryCatch(svd(datModule), error=function(e) ModuleHubProfileSingle(datModule), no_pcs=no.pcs)
    }

    tdatModule <- prep(t(datModule), scale='UV', center=T)
    pcaRes=tryCatch(pca(tdatModule, method='ppca', nPcs=no.pcs), error=function(e) ModuleHubProfileSingle(datModule), no_pcs=no.pcs)
    pcs <- scores(pcaRes)

    #signif(cor(svd1$v[,1:10], pcs ), 3)

    mtitle=paste("PCs of ", modulename," module", sep="")
    no.samples = dim(datModule)[2]

    actualpcs=min(dim(svd1$v)[2], no.pcs)

    #cat(modulename, as.character(i), as.character(actualpcs), "\n")

    #explained variation (%)
    #listPCs[[ no.pcs+1] ][,i]= (svd1$d[1:no.pcs])^2/sum(svd1$d^2)
    tmpMtrx = summary(pcaRes)
    listPCs[[ no.pcs+1] ][,i]= tmpMtrx[1,]

    # this is the first principal component, the ith column of the j-th element in the list
    for (j in c(1:actualpcs) ){
         #print(j)        
         pcj= pcs[,j]; #svd1$v[,j]

         # detect NAs
         jnas = is.na(pcj)
         if (sum(jnas)>0){
             break;
         }

         #signhj=sign(sum(cor(pcj,  t(datModule2))))
         signhj=sign(sum(cor(pcj,  t(tdatModule))))
         #print(is.na(signhj))

         if( !(signhj == 0 | is.na(signhj)) ) pcj=signhj* pcj
         listPCs[[j] ][,i] = pcj
    }

  }

  listPCs
}



ModulePrinComps = function(datexpr,couleur, min_modulesize=10) {

  no.pcs=10

  allmean = mean(datexpr, na.rm=T)

  modlevels= names(table(couleur)) #levels(factor(couleur))

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels

  for(i in c(1:length(modlevels)) ){
 
    #print(paste(i, modlevels[i]) )
    modulename    = modlevels[i]
    restrict1= as.character(couleur)== modulename
    restrict1= ifelse(is.na(restrict1), FALSE, restrict1)

    if(modlevels[i]=="grey"){
       next
    }

    # in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])
    dim(datModule)

    # check whether some samples have missing rate (missing value in column) >80%
    naDat      = is.na(datModule)
    nasBySample= apply(naDat, 2, sum)
    nacutoff   = 0.6*dim(datModule)[1]
    selSamples = nasBySample>=nacutoff
    if(sum(selSamples)>0){
        print("patch samples with missing rate >=0.8")
        naSampleIdx = c(1:length(nasBySample))[selSamples]
        for(x in naSampleIdx){
            #print(paste("Sample Idx=", x) )
            datModule[,x]=ifelse(is.na( datModule[,x]), allmean, datModule[,x])
        }
    }

 
    if(sum(restrict1)<min_modulesize){
       listPCs[[1] ][,i] = datModule[1,]
       next
    }

    imputed=impute.knn(as.matrix(datModule)) #$data for R version < 2.0
    datModule =imputed$data #=imputed for R < 2.9

    datModule=t(scale(t(datModule)))

    # use the hub's profuile as PC if svd doesn't converge
    #svd1=svd(datModule)
    svd1=tryCatch(svd(datModule), error=function(e) ModuleHubProfileSingle(datModule), no_pcs=no.pcs)

    mtitle=paste("PCs of ", modulename," module", sep="")

    no.samples = dim(datModule)[2]

    actualpcs=min(dim(svd1$v)[2], no.pcs)

    #cat(modulename, as.character(i), as.character(actualpcs), "\n")

    #explained variation (%)
    listPCs[[ no.pcs+1] ][,i]= (svd1$d[1:no.pcs])^2/sum(svd1$d^2)

    # this is the first principal component, the ith column of the j-th element in the list
    for (j in c(1:actualpcs) ){
         #print(j)        
         pcj=svd1$v[,j]

         # detect NAs
         jnas = is.na(pcj)
         if (sum(jnas)>0){
             break;
         }

         signhj=sign(sum(cor(pcj,  t(datModule))))

         if( !(signhj == 0 | is.na(signhj)) )  pcj=signhj* pcj
         listPCs[[j] ][,i] = pcj
    }

  }

  listPCs
}


ModuleHubProfiles = function(datexpr, adjmatrix, couleur, min_modulesize=10, whichmodule=NULL) {

  no.pcs=10

  modlevels= names(table(couleur)) #levels(factor(couleur))

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels


    # take the profile of the most connected gene in each module as PC
    #
    colorI = as.character(colcode.reduced)
    kin    = computeModuleLinks(adjmatrix, couleur)

    # find the hub in each module 
    #
    kinIdxed= cbind(kin, c(1:length(couleur)), couleur)
    orderK  = order(-kin)
    kinIdxed= kinIdxed[orderK, ]
    orderK  = order(kinIdxed[,3])
    kinIdxed= kinIdxed[orderK, ]
    
    hubIdx    = rep(0, length(modlevels) )
    for(z in c(1:length(modlevels)) ) {        
        isel      = modlevels[z] == kinIdxed[,3]
        ikinIdxed = kinIdxed[isel, ]

        # extract hubs' profiles
        #
        listPCs[[1] ][,z] = datexpr[,as.integer(ikinIdxed [1,2])]
        hubIdx[z] = ikinIdxed [1,2]
    }

    # extract hubs' profiles
    #
    #listPCs[[1]] = t(datexpr[,as.integer(hubIdx)])
    names(listPCs[[1]]) <- modlevels

    if(!is.null(whichmodule) ){
       wsel = modlevels==whichmodule
       widx = c(1:length(modlevels))[widx]
       return(listPCs[[1] ][,wsel])
    }

    return(listPCs)
}

ModuleHubProfiles_NoAdjM = function(datexpr, kin, couleur, min_modulesize=10, whichmodule=NULL) {

  no.pcs=10

  modlevels= names(table(couleur)) #levels(factor(couleur))

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels


    # take the profile of the most connected gene in each module as PC
    #
    colorI = as.character(colcode.reduced)
    #kin    = computeModuleLinks(adjmatrix, couleur)

    # find the hub in each module 
    #
    kinIdxed= cbind(kin, c(1:length(couleur)), couleur)
    orderK  = order(-kin)
    kinIdxed= kinIdxed[orderK, ]
    orderK  = order(kinIdxed[,3])
    kinIdxed= kinIdxed[orderK, ]
    
    hubIdx    = rep(0, length(modlevels) )
    for(z in c(1:length(modlevels)) ) {        
        isel      = modlevels[z] == kinIdxed[,3]
        ikinIdxed = kinIdxed[isel, ]

        # extract hubs' profiles
        #
        listPCs[[1] ][,z] = datexpr[,as.integer(ikinIdxed [1,2])]
        hubIdx[z] = ikinIdxed [1,2]
    }

    # extract hubs' profiles
    #
    #listPCs[[1]] = t(datexpr[,as.integer(hubIdx)])
    names(listPCs[[1]]) <- modlevels

    if(!is.null(whichmodule) ){
       wsel = modlevels==whichmodule
       widx = c(1:length(modlevels))[widx]
       return(listPCs[[1] ][,wsel])
    }

    return(listPCs)
}

ModuleHubProfileSingle = function(datexpr, no_pcs=10) {

    ngenes = dim(datexpr)[1]

    print("PC by ModuleHubProfileSingle .................. ")

    corrlmatrix = cor( t(datexpr), use = "pairwise.complete.obs")
    corrlmatrix = abs(corrlmatrix)

    diag(corrlmatrix)<- 0
    kin <- apply(corrlmatrix,2,sum, na.rm=TRUE) 

    orderK  = order(-kin)

    # select no.pcs genes which are ebv
    step = as.integer(ngenes/no_pcs)
    if(step<1){step=1; }

    selIdx= seq(from=1,to=ngenes, by=step)

    # mimic values from SVD
    #
    v = t(datexpr[selIdx[1:no_pcs], ])
    d = rep(0, no_pcs)

    ret   = NULL
    ret$v = v
    ret$d = d

    return(ret)
}


plotPrinComps=function(pcVariances,fieldname="", save2file=FALSE)
{
  modulename = names(pcVariances)
  imgname = paste(fieldname, "__PCA-", modulename, ".png", sep="")

  if (save2file)
    openImgDev(imgname)
 
  no.dispPSs = min(8, length( as.matrix(pcVariances)) )
  
  dispnames=rep("Comp.", no.dispPSs )
  for (i in c(1:no.dispPSs) ){
     dispnames[i] =paste("Comp.", as.character(i), sep="")
  }

  mytitle = paste("PCs of ", modulename, " module",sep="")

  par(mfrow=c(1,1), mar=c(4, 6, 4, 2) + 0.1, cex=1)
  barplot(as.numeric(as.matrix(pcVariances))[c(1:no.dispPSs)], 
          names.arg=dispnames,ylim=c(0,1), col=heat.colors(no.dispPSs),
          ylab="Variances",main=mytitle )
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

  if (save2file)
     dev.off()

}

# module PCs and prejection using prcomp()
#
ModulePC_Project = function(datexpr, couleur, var =0.8, maxPCs=3, min_modulesize=2, testdata=NULL, impute_by_knn=T) {

  no.pcs=10

  allmean = mean(datexpr, na.rm=T)

  modlevels= names(table(couleur)) #levels(factor(couleur))

  pcMatrix = NULL; testdata_proj = NULL; pcnames = NULL
  for(i in c(1:length(modlevels)) ){
 
    #print(paste(i, modlevels[i]) )
    modulename    = modlevels[i]
    restrict1= as.character(couleur)== modulename
    restrict1= ifelse(is.na(restrict1), FALSE, restrict1)

    if(modlevels[i]=="grey"){
       next
    }

    # in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])
    dim(datModule)

    # check whether some samples have missing rate (missing value in column) >80%
    naDat      = is.na(datModule)
    nasBySample= apply(naDat, 2, sum)
    nacutoff   = 0.6*dim(datModule)[1]
    selSamples = nasBySample>=nacutoff
    if(sum(selSamples)>0){
        print("patch samples with missing rate >=0.8")
        naSampleIdx = c(1:length(nasBySample))[selSamples]
        for(x in naSampleIdx){
            #print(paste("Sample Idx=", x) )
            datModule[,x]=ifelse(is.na( datModule[,x]), allmean, datModule[,x])
        }
    }

 
    if(sum(restrict1)<min_modulesize){
       pcMatrix = cbind(pcMatrix, t(datModule)); npcs=sum(restrict1)
       pcnames  = c(pcnames, paste(modulename, "_PC", c(1:npcs), sep=""))
       testdata_proj = c(testdata_proj, testModule)
       next
    }

    if(impute_by_knn) {
      imputed=impute.knn(as.matrix(datModule)) #$data for R version < 2.0
      datModule =imputed$data #=imputed for R < 2.9
      rm(imputed);
    } else{
      naSel = is.na(datModule)
      datModule[naSel] = mean(datModule, na.rm=T)
    }

    #datModule=t(scale(t(datModule)))
    # use the hub's profuile as PC if svd doesn't converge
    #svd1=svd(datModule)
    #svd1=tryCatch(svd(datModule), error=function(e) ModuleHubProfileSingle(datModule), no_pcs=no.pcs)
    
    pals   = prcomp(t(datModule ))

    # choose the top N PCs which can explain at least certain variation
    #
    expVar = cumsum(pals$sdev^2/sum(pals$sdev^2))
    sel = c(1:length(expVar))[expVar >= var]
    npcs= ifelse(sel[1]>maxPCs, maxPCs, sel[1])

    # rows as samples and cols as PCs
    pcMatrix = cbind(pcMatrix, pals$x[,c(1:npcs)])
    pcnames  = c(pcnames, paste(modulename, "_PC", c(1:npcs),sep=""))

    #print(paste("explained var=", concatenate(signif(expVar[1:npcs],3), " ")))

    if(!is.null(testdata)) {
       testModule = testdata[restrict1]
       projected  = (testModule-pals$center)%*%(pals$rotation)
       testdata_proj = c(testdata_proj, projected[c(1:npcs)]) # take only the part mapped to top PCs
    }
  }

  colnames(pcMatrix) = pcnames
  if(!is.null(testdata)) {testdata_proj = rbind(testdata_proj); colnames(testdata_proj) <- pcnames }

  ret = as.list(NA, 2)
  ret[[1]] = pcMatrix
  ret[[2]] = testdata_proj 

  rm(testdata_proj)
  rm(pcMatrix)
  rm(datModule)

  return(ret)

}



#take the PCs_list, Traits Matrix, moduleIndex (which module)  as input and output correlation between each of
# top N PCs and each of 20 traits, into a figure and output the table into the logo file
#columns of mydatTraits are traits and rows are mice
ComputePlotCorrelOfPCsAndTraits_InModule=function(pcs_list, moduleIndex, mydatTraits, TopNpcs = 4, fieldname=NULL, logfile=NULL){

        topNpcs = min(TopNpcs, length(pcs_list)-1 )

        #module name
        modulename=colnames(pcs_list[[1]])[moduleIndex]

        corrlPCtraits = NULL
        pcnames = NULL
        for (i in c(1:topNpcs)){
           corrlvectA=abs( cor(pcs_list[[i]][,moduleIndex], mydatTraits, use = "pairwise.complete.obs") )
           corrlvect =ifelse(is.na(corrlvectA), 0, corrlvectA)
           corrlPCtraits=rbind(corrlPCtraits, corrlvect)
           pcnames =c(pcnames, paste("PC.", as.character(i), sep=""))
        }
        rownames(corrlPCtraits) <- pcnames

        mycolors=heat.colors( dim(corrlPCtraits)[1])
        mytitle =paste(modulename, " module:", " correlations of traits and the top ", 
                        as.character(topNpcs), " PC.s", sep="")

        if( !is.null(fieldname)){
          pcTraitsfname=paste(fieldname, "_", modulename, "_PCs_corrlTraits.png", sep='')         
          openImgDev(pcTraitsfname,iwidth = 1278, iheight = 768)
        }
        par(mfrow=c(1,1), mar=c(12, 6, 4, 2) + 0.1, cex=1)
        mp  <- barplot(corrlPCtraits)
        tot <- colMeans(corrlPCtraits)
        text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
        barplot(corrlPCtraits, beside = TRUE, ylab="correlation", las=2, 
             col =mycolors, 
             main=mytitle,
             legend = rownames(corrlPCtraits) )
        abline(h=0,  col="black")
        abline(h=0.3, col="grey")
        #abline(h=-0.3, col="grey")
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
        if( !is.null(fieldname)){
           dev.off()
        }
        
        if( !is.null(logfile))
           appendTableToFile(logfile, round(corrlPCtraits,3), tableTitle=mytitle, myappend=T)

        corrlPCtraits
}


ComputePlotPCsQTL_InModule=function(pcs_list, moduleIndex, TopNpcs = 4, fieldname=NULL, logfile=NULL,minLOD=2.0){
#ComputePlotPCsQTL_InModule=function(pcs_list, moduleIndex, mydatGeno, TopNpcs = 4, fieldname=NULL, logfile=NULL){       

        topNpcs = min(TopNpcs, length(pcs_list)-1 )

        #module name
        modulename=colnames(pcs_list[[1]])[moduleIndex]

        qtl_list = as.list( rep(NULL,topNpcs+1) )
        pcnames = NULL
        for (i in c(1:topNpcs)){
           datGeno$pheno[,2] <- as.numeric( pcs_list[[i]][,moduleIndex] )           
           qtl_list[[i]] = scanone(datGeno, pheno.col=2, method="hk")
           pcnames =c(pcnames, paste("eQTL of PC.", as.character(i), sep=""))
        }

        mycolors=heat.colors(topNpcs)
        mytitle =paste(modulename, " module:", "eQTL of the top ", 
                        as.character(topNpcs), " PC.s",sep="")


        if( !is.null(fieldname)){
          pcQTLfname=paste(fieldname, "_", modulename, "_PCs_QTL.png", sep='')
          openImgDev(pcQTLfname,iwidth = 1278, iheight = 1400)
        }

        par(mfrow=c(topNpcs,1), mar=c(4, 5, 4, 2) + 0.1, cex=0.8)
        for (i in c(1:topNpcs)){            
            mytitle = paste(modulename, " module:", pcnames[i], sep="")
            plot(qtl_list[[i]], col = "green", main=mytitle)
            abline(h=minLOD, col="grey")
        }
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

        if( !is.null(fieldname)){
           dev.off()
           plot(qtl_list[[1]], col = "green", main=pcnames[1])
        }
        qtl_list
}

#~~~~~~~~~~~~~~~~~~ rank-sum based QTL detection ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
RankSumSimple= function(snpvect, exprvect, method="kruskal"){
   sel = !is.na(exprvect) & !is.na(snpvect)
   groups = union(snpvect[sel], NULL)
   
   if(sum(sel)<=2 | length(groups)<=1) {return(1)}

   dat = cbind(value=exprvect[sel], grp=snpvect[sel])
   
   if(method=="kruskal") {
       res = kruskal.test(value~grp, data=dat)
   }else{
       res = wilcox.test(value~grp, data=dat)
   }

   return (res$p.value)
}

# ret=apply(datSnp, 1, RankSumSimple, exprvect=datExpr[2,], method="kruskal")
#
RankSumComplex= function(exprvect, snpDat, method="kruskal"){
   res=apply(snpDat, 1, RankSumSimple, exprvect=exprvect, method=method)
   return ( cbind(res) )
}


# columns are samples for both inputs, method could be wilcox or kruskal
#  but wilcox can only work for 2 classes while kruskal for multiple classes
#
# res = apply(datExpr[c(1:2), ], 1, RankSumComplex, snpDat=datSnp, method="kruskal")
#
QTL_by_RankSum = function(expDat, snpDat, method="kruskal"){
  
  # return result with rows as SNP and cols as genes 
  ps = apply(expDat, 1, RankSumComplex, snpDat, method=method)
  
  return( list(ps=t(ps)) )
}



# return the interval of the maximum valley (blank area in the profile)
#--
#  \                 /\             |\
#   ----------------/  \------------| \----
#      max valley
#
findMaxValleyInProfile =function(profile, xpos, ratio=0.75){   

   nopoints = length(profile)
   index    = c(1:nopoints)

   sel   = profile<Inf
   nonInfprofile = profile[sel]
   pmin = min(nonInfprofile )
   pmax = max(nonInfprofile )
   yrange = pmax-pmin

   # explore the best cut
   newrat = ratio
   for(nr in seq(ratio, 0.95, 0.05) ) {
     threshB= pmin + yrange*nr
     pminusB= profile  - threshB
     selB   = pminusB >0
     transitpoints = c(1, index[selB], nopoints); nlen = length(transitpoints)
     shiftpos = c(transitpoints[2:nlen], transitpoints[nlen])
     diff     = shiftpos - transitpoints
     if(max(diff)>2){break;}
   }

   # transition points are those points that are just exactly above the threshold
   #

   # get the transition points' physical position on x-axis
   #
   transitxpos   = c(xpos[1], xpos[selB])
   notransits    = length(transitxpos)
   shifttransxpos= c(transitxpos[c(2:notransits)], xpos[nopoints])
   distance      = as.vector(shifttransxpos-transitxpos)
   maxDidx       = which.max(abs(distance)) #left side of max valley 

   maxValleyLeftxPos= transitxpos[maxDidx]

   # give a little space between
   #
   if(distance[maxDidx]>0){
      maxValleyLeft = maxValleyLeftxPos + distance[maxDidx]/8
   
   } else{ # descendant order of the index
      maxValleyLeft = maxValleyLeftxPos + distance[maxDidx] - distance[maxDidx]/8
   }

   return (maxValleyLeft)   
}

#######################################################################################################
#                    Plot multiple eQTL profiles in one figure
#
#  eQTL_profiles: row as markers and column as genes
#
#    type: 1-character string giving the type of plot desired.  The
#          following values are possible, for details, see 'plot': '"p"'
#          for points, '"l"' for lines, '"o"' for overplotted points and
#          lines, '"b"', '"c"') for (empty if '"c"') points joined by
#          lines, '"s"' and '"S"' for stair steps and '"h"' for
#          histogram-like vertical lines.  Finally, '"n"' does not
#          produce any points or lines.

plotMultieQTLProfiles=function(eQTL_profiles, smarkerpos, smarkerXlabel, sboundaries, filename=NULL, minLOD=2.5, singleView=F, 
                ylab="LOD score", plottype='l', colorcode=NULL, mcex=1, mywid = 1278, plotLegend=TRUE, maintitle=""){
        no.profiles  = dim(eQTL_profiles)[2]
        no.pmarkers   = dim(eQTL_profiles)[1]
        profilenames <- colnames(eQTL_profiles)

        if(is.null(colorcode)) {
           colorcodeC=c("red","black","blue","green","magenta","yellow","brown","pink","purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow","grey")
        }else{
           colorcodeC=colorcode
        }
        # make x labels by concatenating chromosome labels
            bool.chroms=smarkerXlabel != ""
            chroms=smarkerXlabel[bool.chroms]
            no.chroms=length(chroms)-1 # last one is "."
            if (no.chroms>=16){
               xlabel="genome "
            }else{
               xlabel="markers on chromosome "
               for (i in c(1:no.chroms) ){
                    if (i==1){
                      if (no.chroms==1)
                         xlabel=paste(xlabel, chroms[i], sep="") 
                      else
                         xlabel=paste(xlabel, chroms[i], ",", sep="")
                    }else if (i==no.chroms){
                      xlabel=paste(xlabel, chroms[i], sep="")
                    }else{
                      xlabel=paste(xlabel, chroms[i], ",", sep="")
                    }
               }
            }

    showaxis=T
    if (no.chroms>=2)
      showaxis=F

    if(singleView){
        myhei = 400
    }else{
        myhei = 300*no.profiles
    }

    if( !is.null(filename)){
          openImgDev(filename,iwidth = mywid, iheight = myhei)
    }

    if(singleView){
            #itype=c(1:no.profiles)
            itype=rep(1, no.profiles)
            icolor=colorcodeC[c(1:no.profiles)]

            # simple way to figure out where to put the legend, x position is tricky
            maxprofiles = apply(eQTL_profiles, 1, max)
            if(1==2){
            maxprforders = order(-maxprofiles)
            peakposIdx = maxprforders[1]     # peak index
            peakpos = smarkerpos[peakposIdx] # peak position
            leftmean = mean(maxprofiles[c(1:(peakposIdx-1))]) # mean max-profiles on the left of Peak value
            rightmean= mean(maxprofiles[c(peakposIdx:no.pmarkers)]) # mean max-profiles on the right of Peak value

            minpos=min(smarkerpos)
            maxpos=max(smarkerpos)
            adjpos= maxpos-(maxpos-minpos)/4
            if ( abs(peakpos-mean(smarkerpos)) < (maxpos-minpos)/8 ){
                # in this case, the peak is in the middle, so we need consider 
                # the mean profiles in both sides of the peak
                legendx = ifelse(leftmean<rightmean, minpos, adjpos)
            }else{
               legendx = ifelse(peakpos-minpos>maxpos-peakpos, minpos, adjpos)
            }
            }

            legendx = findMaxValleyInProfile(maxprofiles, xpos=smarkerpos)

            if( !is.null(filename)){
              par(mfrow=c(1, 1), mar=c(5, 5, 3, 2) + 0.1, cex=mcex)
            }
            matplot(x=smarkerpos, y=eQTL_profiles, type=plottype, lty=itype, col = icolor,
                xlab=xlabel, 
                ylab=ylab, main=maintitle, axes=showaxis)

            if(plotLegend) { 
              if ( max(maxprofiles) < Inf){
                 legend(legendx, max(eQTL_profiles), legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)
              }else{
                 legend(legendx, max(maxprofiles[maxprofiles < Inf]), legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)
              }
            }

            if(!showaxis){
               axis(1, at =smarkerpos[sboundaries], labels = smarkerXlabel[sboundaries], las=1.2)
               axis(2,las=2)
            }
            abline(h=minLOD, col="grey")

    }else{ # multiple sub plots

        par(mfrow=c(no.profiles, 1), mar=c(4, 5, 4, 2) + 0.1, cex=mcex)
        for (i in c(1:no.profiles)){            
            #plot(eQTL_profiles[,i], col = "green", main=profilenames[i])

            plot(smarkerpos, eQTL_profiles[,i], type=plottype, col = colorcodeC[1],
                xlab=xlabel, 
                ylab=ylab, main=profilenames[i], axes=showaxis)
            if(!showaxis){
               axis(1, at =smarkerpos[sboundaries], labels = smarkerXlabel[sboundaries], las=1)
               axis(2,las=2)
            }
            abline(h=minLOD, col="grey")
        }
     }

    if( !is.null(filename)){
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
        dev.off()
    }
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ A simplfied version of multi curves in one plot $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

plotMultiProfiles=function(eQTL_profiles, xindex, xlabel="",ylabel="", mainlabel="", 
                           filename=NULL, plottype="l", ltype=1, myhei = 400, mywid = 500, mypointsize=12, myunits="px", myres = 72, mycolors=NULL,
                           mycex=0.9, mycexText=1, deltaLegendX=0, deltaLegendY=0, curvecolor=NULL, show_legend=T,  showaxis=TRUE,
                           horizontal_line=NULL, display_string=NULL, dispYp=0, dispXp=4, linewidth=1, pointcolor=NULL, ylim=NULL)
{
    no.profiles  = dim(eQTL_profiles)[2]
    no.pmarkers  = dim(eQTL_profiles)[1]
    profilenames = colnames(eQTL_profiles)
 
    if(!is.null(mycolors)) {
      colorcodeC=mycolors
    }else {
      colorcodeC=c("red", "black", "blue","green","grey","yellow","pink","magenta", "purple", "greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", "lightyellow", "grey")
    }

    #itype=c(1:no.profiles)
    itype=rep(1, no.profiles)
    if(is.null(curvecolor)) {
       icolor=colorcodeC[c(1:no.profiles)]
    } else{
       icolor=curvecolor
    }

    if(no.profiles==1 & !is.null(pointcolor)) {
       icolor= pointcolor
    }

    # simple way to figure out where to put the legend, x position is tricky
    maxprofiles = apply(eQTL_profiles, 1, max)
    if(1==2){
     maxprforders = order(-maxprofiles)
     peakposIdx = maxprforders[1]     # peak index
     peakpos = xindex[peakposIdx] # peak position
     leftmean = mean(maxprofiles[c(1:(peakposIdx-1))]) # mean max-profiles on the left of Peak value
     rightmean= mean(maxprofiles[c(peakposIdx:no.pmarkers)]) # mean max-profiles on the right of Peak value

     minpos=min(xindex)
     maxpos=max(xindex)
     adjpos= maxpos-(maxpos-minpos)/4
     if ( abs(peakpos-mean(xindex)) < (maxpos-minpos)/8 ){
         # in this case, the peak is in the middle, so we need consider 
         # the mean profiles in both sides of the peak
         legendx = ifelse(leftmean<rightmean, minpos, adjpos)
     }else{
         legendx = ifelse(peakpos-minpos>maxpos-peakpos, minpos, adjpos)
     }
    }

     maxprforders = order(-maxprofiles)
     peakposIdx = maxprforders[1]     # peak index
     peakpos = xindex[peakposIdx] # peak position

     legendx = findMaxValleyInProfile(profile=maxprofiles, xpos=xindex) 

    if( !is.null(filename)){
          openImgDev(filename,iwidth = mywid, iheight = myhei, ipointsize=mypointsize, iunits=myunits, ires=myres)
          par(mfrow=c(1, 1), mar=c(5, 5, 3, 2) + 0.1, cex=mycex)
    }
     #colnames(eQTL_profiles)<- c("l")
     
     #matplot(x=xindex, y=eQTL_profiles, type=ltype, lty=itype, col = icolor,
     #matpoints(x=xindex, y=eQTL_profiles, type=ltype, lty=itype, col = icolor,

     if(ltype[1]=="b") {
         matplot(x=xindex, y=eQTL_profiles, type=plottype, lty=ltype, col = icolor, lwd=linewidth,
                xlab=xlabel, ylab=ylabel, main=mainlabel, axes=showaxis, pch="*", cex=mycex, ylim=ylim)
     }else{
         matplot(x=xindex, y=eQTL_profiles, type=plottype, lty=ltype, col = icolor, lwd=linewidth,
                xlab=xlabel, ylab=ylabel, main=mainlabel, axes=showaxis, pch="*", cex=mycex, ylim=ylim)
     }

     if(no.profiles==1 & !is.null(pointcolor)) {
        idx = findBoundary(pointcolor); no.ii= length(idx)
        #print(idx)
        idx2= c(idx, length(pointcolor))
        for(ii in c(1:no.ii) ){
            iisel = c( idx2[ii]: idx2[ii+1])
            lines(cbind(xindex[iisel], eQTL_profiles[iisel,1]),col=pointcolor[idx2[ii]], lwd=linewidth)
        }
     }


     if(!is.null(horizontal_line)){
        abline(h=horizontal_line, col="gray", lty="dotted")
     }

     #  legend(legendx, max(eQTL_profiles), legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)      
     if ( max(maxprofiles) < Inf){
        legdx =legendx+ deltaLegendX; legdy= max(eQTL_profiles)+ deltaLegendY
     }else{
        legdx =legendx+ deltaLegendX; legdy= max(maxprofiles[maxprofiles < Inf])+ deltaLegendY
     }    
     if(no.profiles>1 & show_legend) {
        legend(legdx, legdy, legend=profilenames,lty=itype,lwd=2,col=icolor, ncol =1, cex=0.8,pt.cex=1)
     }

     # pos: Values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified coordinates.
     if(!is.null(display_string)){
         #text(x=legdx+deltaLegendX*2,y=legdy-2*deltaLegendY, labels= display_string, col="black", pos=1)
         xt = min(xindex) + (max(xindex)-min(xindex))/dispXp
         if(abs(xt-peakpos) < 0.1*(max(xindex)-min(xindex))){
           xt = min(xindex) + (max(xindex)-min(xindex))/6  
         }
         text(x=xt,y=legdy-2*deltaLegendY+dispYp, labels= display_string, col="black", pos=1, cex=mycexText)
     }


     if( !is.null(filename)){
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
        dev.off()
     }
}



PlotMultipleTimeSeries=function(datedposts, legend_ncol=1, maintitle="", fimg=NULL, 
                         imgwid=1000, imghei=500, mcex=1.3) {

  # get all dates, posts= overall posts on each date
  #
  tbdates = table(datedposts[,1]); no.points = length(tbdates); xlabel=names(tbdates); posts = as.integer(tbdates)
  maxv = max(as.integer(tbdates)); maxv = as.integer((maxv/100+0.5))*100
  if(max(posts)>maxv) {
    maxv=maxv+50
  }
  step = as.integer(maxv/100+0.5)*10
  if(step>50) {
     step = as.integer(step/100+0.5)*100
  }
  yseq = seq(0, maxv, by=step)

  xcategories = names(table(datedposts[,2])); no.categ = length(xcategories)

  # find number of posts for each individual catgeory, blog, microblog and forum
  final = cbind(names(tbdates), as.integer(tbdates))
  if(no.categ>1) {
     xtabs = tapply(datedposts[,1], as.factor(datedposts[,2]), FUN=table)
     for(xcg in names(xtabs)){
       xtmp = unlist(xtabs[xcg]); xnames = getSecondPart(names(xtmp), "\\.", 2)
       xres = cbind(xnames, as.integer(xtmp) )
       final= merge(final, xres, by.x=1, by.y=1, all=T)
     }
     final = as.matrix(final)
     colnames(final) <- c("date", "all", names(xtabs))
     final = ifelse(is.na(final), 0, final)
     no.points = dim(final)[1]

     finalDat  = matrix(as.integer(final[,-1]), nrow=no.points)
     colnames(finalDat) <- c("all", names(xtabs))

  } else{
     colnames(final) <- c("date", xcategories)
     finalDat  = cbind(as.integer(final[,2]))
     colnames(finalDat) <- xcategories
  }
  fcols = dim(finalDat)[2]
  xindex = c(1:no.points);

  colors=c("red", "blue", "green", "brown")

  openImgDev(fimg, iwidth =imgwid, iheight = imghei)
  par(mfrow=c(1,1), mar=c(8, 5, 3, 0) + 0.1, cex=mcex)#, srt=90)

  matplot(x=repmat(cbind(xindex),1,fcols), y=finalDat, type='l', lty="solid", col = colors, lwd=2,
                xlab="", ylab="", main=maintitle, axes=F, pch="*", frame.plot=F, lend="square")

  axis(1, at =xindex, labels =xlabel, las=2, 
               col = "brown", col.axis="black", lwd = 1, tick = T, line =1)
  axis(2, at =yseq, labels =yseq, las=2, 
               col = "brown", col.axis="black", lwd = 1, tick = T, line =1)
  for(j in yseq){
    abline(h=j, col="gray", lty="dotted")
  }
  #text(x=xindex, y=posts, label="*", cex=1.8, col="brown")
  #delt = (max(posts)-min(posts))/100
  hcol="salmon";deltx = 0.2; delty=0.2
  hcol="salmon";deltx = length(xindex)/200; delty=(max(posts)-min(posts))/100

  for(x in c(1:fcols)){
    for(i in c(1:length(xindex)) ) {
       rect(xleft=xindex[i]-deltx, ybottom=finalDat[i,x]-delty,xright=xindex[i]+deltx,
             ytop=finalDat[i,x]+delty, col="white", border=colors[x], lwd=1, lty=1)
    }
  }

  if(fcols>1) {
     ltitle = colnames(finalDat)
     ltitle = replaceString(ltitle, "microblogs", "mblog")  
     ltitle = replaceString(ltitle, "blogs", "blog")  
     ltitle = replaceString(ltitle, "forums", "forum")  

     legendx = findMaxValleyInProfile(c(posts,max(posts)), xpos=c(1:(no.points+1)), ratio=0.5)
     legendy = max(posts)
     legend(legendx-1, legendy, legend=ltitle, box.lty=1,box.lwd=1,box.col="red", ncol=legend_ncol, 
                fill=colors[1:fcols], cex=0.8,pt.cex=0.8)
  }

  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
}


# input is a globale genotype data structure "datGeno" and 
#         a gene's expression vector across all samples
# output is a vector of the QTL LOD scores on the whole genome
# iqtl[1:10,]
#       chr       pos        lod       AA       AB       BB    sigma
#p45178   1  0.205459 0.04779396 4.255317 4.111057 4.560112 5.634021
#p45787   1  1.872961 0.04778489 4.255144 4.111150 4.560071 5.634021
#p44584   1  6.692455 0.04311972 4.166914 4.149415 4.560018 5.634316
#
qtlOnGenome = function(exprvector){
    datGeno$pheno[,2] <- as.numeric(exprvector)
    iqtl = scanone(datGeno, pheno.col=2, method="hk")
    round(as.numeric(iqtl$lod),2)
}



# plot correlation of each trait pf PCs across modules, different from pervious plot which is based on
# module(all correlations are plot in one figure for each module)
# list_corrlModulePCsTraits[[1]], the first module; list_corrlModulePCsTraits[[2]], 2nd module
#        WeightG      AbFat    OtherFat    TotalFat X100xfat.weight      Trigly
#PC.1 -0.3852039  0.2714077 -0.01301723  0.18511831      0.37228074 -0.17832184
#PC.2 -0.2425847 -0.1467176 -0.09857543 -0.15055786     -0.03087848 -0.08506316
#PC.3  0.3856605 -0.4090452  0.01609531 -0.26830352     -0.45573164  0.32946309
#PC.4 -0.2344006  0.1694582 -0.11794072  0.07000631      0.15279737 -0.06884160

plotCorrlPCsTraitsCrossModule=function(myCorrlModulePCsTraits_list,modulenames, whichPC=1, fieldname=NULL){

        no.modules = length(myCorrlModulePCsTraits_list)
        no.pcs   = dim(myCorrlModulePCsTraits_list[[1]])[1]
        no.traits=dim(myCorrlModulePCsTraits_list[[1]])[2]
        subfigures=4 # number of traits in one figure
        no.figures= as.integer(no.traits/4) + (no.traits%%4>0)

        traitnames = colnames(myCorrlModulePCsTraits_list[[2]])
        pcnames    = rownames(myCorrlModulePCsTraits_list[[1]])

        for ( ifig in c(1:no.figures) ){
            start= (ifig-1)*subfigures + 1   #first trait
            end= ifig*subfigures             #last trait
            if(end>no.traits)
                end=no.traits
            
            actualsubfigures = end - start + 1

            #here we re-organize the correlation data into a matrix with rows as traits and 
            # columns as modules
            corrlMatrix = NULL
            mytitles =  NULL
            for (j in c(start:end) ){ #traits as rows
                rvect=NULL
                for(m in c(1:no.modules) ) #module as columns
                   rvect=c(rvect, (myCorrlModulePCsTraits_list[[m]])[whichPC, j] )
                corrlMatrix = rbind(corrlMatrix, rvect)
                mytitle = paste("r(", traitnames[j], " , ", pcnames[whichPC], ")", sep="")
                mytitles = c(mytitles, mytitle)    
            }

            if( !is.null(fieldname)){
              imgfname=paste(fieldname, "_corrl_PC", as.character(whichPC), 
                             "_Traits", as.character(start), "to", as.character(end), ".png", sep='')
              openImgDev(imgfname,iwidth = 1278, iheight = 1400)
            }
            par(mfrow=c(actualsubfigures,1), mar=c(6, 5, 4, 2) + 0.1, cex=0.8)

            #display corelations of PC with trait in each module
            for (j in c(1:actualsubfigures) ){               
               barplot(as.numeric(corrlMatrix[j,]), ylab="correlation",                       
                 col = modulenames, main=mytitles[j], names.arg=modulenames, las=2)
                 abline(h=0.3, col="grey")
                 abline(h=0, col="black")
                 abline(h=-0.3, col="grey")
            }
            par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
            if( !is.null(fieldname)){
              dev.off()
            }
      }
}



# plot Gene Significance across modules,  a single subfigure contains the means of a trait in
# all modules, a figure will have 4 sub figures
#> genesig_means[1:5,1:5]
#         WeightG     AbFat   OtherFat   TotalFat X100xfat.weight
#black 0.33819315 0.2371843 0.04326499 0.16219943       0.3253637
#blue  0.51582538 0.1813136 0.19555568 0.12881716       0.2365684
#brown 0.62463205 0.1791032 0.24054805 0.11426613       0.2576446
#cyan  0.26215453 0.1969778 0.11479211 0.10464036       0.1801414
#green 0.08332405 0.0964498 0.09881134 0.07656996       0.1036234

# subfigures=4: # number of traits in one figure

plotMeanGeneSignifCrossModule=function(means_matrix, sderrs_matrix, fieldname=NULL, 
signifiance_level=0.3, signifiance_label="", imgwid=1278, imghei=1400, xspace=6,
ylabel="mean(|r|)", mcex=0.8, xlab_cex=1, xaxisLabel_rotate=NULL, usr_delta=0,  xlabel_orientation=0,
ipointsize=12, iunits="px", ires=300, icompression="lzw", subfigures=4){

        no.modules = dim(means_matrix)[1]
        no.traits  = dim(means_matrix)[2]       
        no.figures= as.integer(no.traits/subfigures) + (no.traits%%subfigures>0)

        traitnames = colnames(means_matrix)
        modulenames= rownames(means_matrix)

        allcolors = colors()
        overlaps  = intersect(allcolors, modulenames)
        if(length(overlaps) < length(modulenames) ) {
            colornames = rep("blue", length(modulenames))
        }else{
            colornames = modulenames
        }

       interv = 0.5; space4grp = interv+1
       xcoord = seq(1, space4grp*no.modules, space4grp)

        for ( ifig in c(1:no.figures) ){
              
            start= (ifig-1)*subfigures + 1   #first trait
            end= ifig*subfigures             #last trait
            if(end>no.traits)
                end=no.traits
            
            actualsubfigures = end - start + 1            

            if( !is.null(fieldname)){
              imgfname=paste(fieldname, "_genesignif",
                             "_Traits", as.character(start), "to", as.character(end), ".png", sep='')
              openImgDev(imgfname,iwidth = imgwid, iheight = imghei, ipointsize=ipointsize, iunits=iunits, ires=ires, icompression=icompression)
              
            }
            par(mfrow=c(actualsubfigures,1), mar=c(xspace, 5, 3, 1) + 0.1, cex=mcex)

            #display trait-based gene significance in each module
            for (j in c(start:end) ){
               # find the appropriate position for putting text legend
               #
               jabove = means_matrix[,j]>min(signifiance_level)
               if(sum(jabove)==0){
                  jlpos = 0
               } else{
                  jabove[1]=T
                  jabove[no.modules]=T
                  jleftbound = c(1:no.modules)[jabove]
                  jrightbound =c(jleftbound[-1], no.modules)
                  jdistw = which.max(jrightbound - jleftbound)
                  jdist  = max(jrightbound - jleftbound)
                  jlpos  = jleftbound[jdistw] + 0.2*jleftbound[jdistw] + jdist/4
                  
               }

                 barplot(as.numeric(means_matrix[,j]), ylab=ylabel, axisnames=F,  beside = TRUE, space=interv,                
                      col = colornames, main=traitnames[j], legend=FALSE)#names.arg=modulenames, las=2, xpd=F)

                 #xcoord =  c(1:length(means_matrix[,j]))
                 if(is.null(xaxisLabel_rotate) ) {
                     axis(1, at =xcoord, labels =modulenames, las=xlabel_orientation, 
                          col = "black", col.axis="black", lwd = 0, tick = T, line =1)
                     abline(h=0, col="black")
                     err.bp(as.vector(means_matrix[,j]), as.vector(sderrs_matrix[,j]), two.side=T)

                 } else{
                      axis(1, at =xcoord, labels =rep("", length(means_matrix[,j])), las=xlabel_orientation, 
                           col = "black", col.axis="black", lwd = 0, tick = T, line =1)
                      text(xcoord, par("usr")[3] - usr_delta, srt = xaxisLabel_rotate, adj = 1, labels = modulenames, 
                             xpd = TRUE, cex=xlab_cex)
                      par(srt=0) 
                 }

                 for (x in c(1:length(signifiance_level)) ) {
                     abline(h=signifiance_level[x], col="grey60")
                     text(jlpos,signifiance_level[x], signifiance_label[x], col = "black", adj = c(0, -.2) )
                 }

            }
            par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
            if( !is.null(fieldname))
              dev.off()
      }
}


#********************************************************************************************************
# the corrlmatrix is the correlations between a few network properties (as row) 
#   and a trait (traitname) across modules (column)
# 
#
plotCorrelOfNetiesAndTrait_AcrossModules=function(mycorrlmatrix, mytraitname, fieldname=NULL){

        no_neties   = dim(mycorrlmatrix)
        #corrlmatrix = abs(mycorrlmatrix)
        corrlmatrix = mycorrlmatrix

        traitname = replaceChars(mytraitname, ".", "")

        #mycolors=heat.colors( dim(corrlmatrix)[1])
        mycolors=c("lightgreen","blue", "red", "yellow")
        mytitle =paste("correlations between gene network properties and gene significance based on ", traitname, sep="")

        if( !is.null(fieldname)){
          netiesTraitfname=paste(fieldname, "_genesignif_Neties_", traitname, ".png", sep='')         
          openImgDev(netiesTraitfname,iwidth = 1278, iheight = 768)
        }
        par(mfrow=c(1,1), mar=c(11, 6, 4, 2) + 0.1, cex=1)
        mp  <- barplot(corrlmatrix)
        tot <- colMeans(corrlmatrix)
        text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
        barplot(corrlmatrix, beside = TRUE, ylab="correlation", las=2, 
             col =mycolors, 
             main=mytitle,
             legend = rownames(corrlmatrix) )
        abline(h=0, col="black")
        abline(h=0.3, col="grey")
        abline(h=-0.3, col="grey")
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
        if( !is.null(fieldname)){
           dev.off()
        }
}



#count LOD scores on genome for all the genes in a module
LODscoreCOUNT=function(TRAITS, LODscoreCut=2.5, title1="eQTL profile"){
   QTLlist=list(NULL)

   maxLOD = 0

   datGeno$pheno=data.frame( datGeno$pheno[,1], TRAITS)
   no.traits= dim(datGeno$pheno)[[2]]

   QTLlist=scanone(datGeno, pheno.col=2, method="hk")

   maxLODmarker=data.frame(matrix(-666, nrow=no.traits-1, ncol= length(QTLlist$pos)  +1) )
   chromo=QTLlist$chr
   names(maxLODmarker)[[1]]="Trait"

   names(maxLODmarker)[-1]=paste(names(maxLODmarker)[-1],"chr",chromo,sep=".") 
   maxLODmarker[,1]=names(TRAITS)
   for (i in c(1:no.traits) ) {
      datGeno$pheno[,2] <- as.numeric( pcs_list[[i]][,moduleIndex] )           
      QTLlist=scanone(datGeno, pheno.col=i, method="hk");
      maxLODmarker[i-1,-1]=QTLlist$lod
   }

   count1=apply( maxLODmarker[,-1]> LODscoreCut, 2,sum)
   myylab=paste("no. of LODs above ", as.character(LODscoreCut), sep="")

   plot(count1,col=as.character(chromo), main=title1,xlab="marker number", ylab=myylab)
   list(maxLODmarker =maxLODmarker, chromo=chromo)

}

#here, we attach the information in the minor Matrix to the major matrix, if some row keys 
#in the major matrix don't match any in the minor matrix, we use missing label for those keys
# keepAllPrimary=T: return all primary row-elements
# keepAllPrimary=F: return all MISSED primary row-elements
# keepAllPrimary=NULL: return the common elements
mergeTwoMatricesByKeepAllPrimary=function(primaryMatrix, minorMatrix, missinglabel="", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
{
  no.promarycols <- dim(primaryMatrix)[2]
  no.mustbegenes <- dim(primaryMatrix)[1]

  # we add in one more column to indicate which genes are mustbeincluded after being merged with mcg
  keyword="mustbeused"
  mustbeGenesMatrix = cbind(primaryMatrix, c(1:no.mustbegenes), rep(keyword, no.mustbegenes) )

  if (is.null(colnames(primaryMatrix)) ){
    colnames(mustbeGenesMatrix) <- c( c(1:no.promarycols), "primorder", keyword)
  }else{
    colnames(mustbeGenesMatrix) <- c( colnames(primaryMatrix), "primorder", keyword)
  }
  dim(mustbeGenesMatrix)

  if(is.null(keepAllPrimary) ){ #normal merge: to have the common elements
    myMatrix = merge(mustbeGenesMatrix, minorMatrix, by.x=1, by.y=1,all.x=F,sort=F,all=F)
  }else{
    myMatrix = merge(mustbeGenesMatrix, minorMatrix, by.x=1, by.y=1,all.x=T,sort=F,all=T)
  }

  dim(myMatrix)
  nocols.mymatrix <- dim(myMatrix)[2]

  #the mustbeused genes which are not included in minor have NAs in the column $mustbeused
  #so we can use this information to figure out which mustbeused genes missing in minorMatrix
  myMatrix[,nocols.mymatrix] = ifelse( is.na(myMatrix[,nocols.mymatrix]), missinglabel, as.character(myMatrix[,nocols.mymatrix]) )

  orders = order( as.numeric(as.matrix(myMatrix[, no.promarycols+1])))
  if (keepPrimaryOrder)
      myMatrix = myMatrix[orders,]

  if (is.null(keepAllPrimary) ){
     selected = rep(T, dim(myMatrix)[1])
  }else{
     if (keepAllPrimary)
       selected = !(is.na(myMatrix[, no.promarycols+2]))
     else #return the row-elements in minor which are missed in primary
       selected = is.na(myMatrix[, no.promarycols+2])
  }

  sum(selected)

  #keep the primary matrix and remove the mustbeused column
  myMatrix[selected, -c(no.promarycols+1, no.promarycols+2)]
}



# Here we return the boolean vactor which indicates whether each element in the primary list is
# in theminor list, the boolean vector is in the same order as the primary list
#
findListInPrimarysPosition=function(primaryMatrix, minorMatrix)
{
  #no.mustbegenes <- dim(minorMatrix)[1]
  # we add in one more column to indicate which genes are mustbeincluded after being merged with mcg
  #keyword="mustbeused"
  #mustbeGenesMatrix = cbind(minorMatrix, rep(keyword, no.mustbegenes) )
  #colnames(mustbeGenesMatrix) <- c( colnames(minorMatrix), keyword)
  #dim(mustbeGenesMatrix)
  #It is always dangerous to use merge if you need keep the selected.bool.vector in the same order 
  # as primaryMatrix, whatever, the order will be changed after merging
  #
  #myMatrix = merge(mustbeGenesMatrix, primaryMatrix, by.x=1, by.y=1,all.y=T, sort=F,all=F)
  #dim(myMatrix)
  #nocols.mymatrix <- dim(myMatrix)[2]
  #selected = ifelse(is.na(myMatrix$mustbe), FALSE, TRUE)
  #sum(selected)

  minorlist  = as.character( as.matrix(minorMatrix[,1]) )
  primarylist= as.character( as.matrix(primaryMatrix[,1]) )

  no.primary = dim(primaryMatrix)[1]
  selected=rep(FALSE, no.primary)
  primaryseq = c(1:no.primary)

  for (each in minorlist){
      isel = (each==primarylist)
      if (sum(isel) >0){
        selected[ primaryseq[isel] ] = TRUE
	#cat(sum(isel),primaryseq[isel] , "\n")
      }
  }
  
  selected
}

# re-organize the QTL in the gene order 
#
# turn gene name into gene symbols
# 
# 8	chr8 1.3e+08	7	709	25	11466	6.00E-04	NM_018412;NM_024525;NM_024604;NM_015990;Contig2399_RC;NM_022051;NM_021205
#
# called by C:\Zhang\WntSignaling\HypoxiaNew_Oct20\R-decompose_qtl-tables.R
#
expandQTLenrich_by_genes = function(fqtlenrich, modulecolname="LocusBin", 
                                 genecolname="Probe.Set", geneinfo, signat_source, fout)
{
    allMatrix<- read.delim(fqtlenrich, sep="\t", header=T)
    allMatrix<- as.matrix(allMatrix)
    dim(allMatrix)

    no.rows = dim(allMatrix)[1]

    acolnames = colnames(allMatrix)
    aindices  = getMatchedIndex(cvector=acolnames, 
                                subvect=c(modulecolname, genecolname) )
    midx      = aindices[1]
    gidx      = aindices[2]

    # to shorten the gene information matrix
    #
    allgenes  =splitString(allMatrix[, gidx], ";")

    # remove the duplicates
    #
    uniques   = names(table(allgenes))
    geneinfo2 = merge( cbind(uniques),geneinfo, by.x=1, by.y=1, all=F)
    
    for(i in c(1:no.rows)){
        imodule  = allMatrix[i, midx]
        iallgenes= allMatrix[i, gidx]
        igenes  =splitString(iallgenes, ";")

        # remove the duplicates
        #
        iuniques= names(table(igenes))
        imerged = merge( cbind(iuniques),geneinfo2, by.x=1, by.y=1, all=F)
        imerged = as.matrix(imerged)
        irows   = dim(imerged)[1]

        final   = cbind(imerged[,2], imerged[,1], 
                        rep(imodule,irows), rep(signat_source, irows) )

        write.table(final, fout, sep="\t",quote=FALSE, col.names=F, 
                    row.names=FALSE, append=T)
    }

    rm(geneinfo2)        
}



#--------------------------------------- Steve -----------------------------------------------------------
#--------------------------------------- Steve -----------------------------------------------------------

#The function ScaleFreePlot1 creates a plot for checking scale free topology
ScaleFreePlot = function(kk,no.breaks=15, mtitle="",truncated1=TRUE, tofile="", outputFitness=FALSE){
    
    cut1=cut(kk,no.breaks)
    binned.k=tapply(kk,cut1,mean)
    freq1=tapply(kk,cut1,length)/length(kk)

    #remove NAs
    noNAs = !(is.na(freq1) | is.na(binned.k))
    freq.noNA= freq1[noNAs]
    k.noNA  = binned.k[noNAs]

    #remove Zeros
    noZeros  = !(freq.noNA==0 | k.noNA==0) 
    freq.log = as.numeric(log10(freq.noNA[noZeros]))
    k.log    = as.numeric(log10(k.noNA[noZeros]))

    if(tofile!=""){
       openImgDev(tofile, iwidth = 700, iheight = 700, ipointsize = 10)
    }

    par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1, cex=1.8)

    plot(k.log, freq.log, xlab="log10(k)",ylab="log10(p(k))" )
    lm1=lm( freq.log ~ k.log, na.action=na.omit) 
    lines(k.log,predict(lm1),col=1)
    if (truncated1==TRUE) { 
       lm2=lm(freq.log ~ k.log +I(10^k.log) );
       lines(k.log,predict(lm2),col=2);

       if (mtitle !=""){
           ititle=paste(mtitle, ", scale R^2=",
                          as.character(round(summary(lm1)$adj.r.squared,2)) , 
                          ", trunc. R^2=",
                          as.character(round(summary(lm2)$adj.r.squared,2)),
                          ", slope=", as.character(round(lm1$coef[[2]], 2)),
                          sep="")
       }else{
           ititle=paste("scale R^2=",
                          as.character(round(summary(lm1)$adj.r.squared,2)) , 
                          ", trunc. R^2=",
                          as.character(round(summary(lm2)$adj.r.squared,2)),
                          ", slope=", as.character(round(lm1$coef[[2]], 2)),
                          sep="")
       }

    }else{
       if (mtitle !=""){
            ititle = paste(mtitle, ", scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)), 
                     ", slope=", as.character(round(lm1$coef[[2]], 2)), sep="")
       }else{
            ititle = paste("scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)), 
                     ", slope=", as.character(round(lm1$coef[[2]], 2)), sep="")
       }
    }
    title(ititle)

    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    if(tofile !=""){
       dev.off()
    }
        
    if(outputFitness==TRUE){
      if(truncated1){
        myoutput = c(summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared, lm1$coef[[2]])
      }else{
        myoutput = c(summary(lm1)$adj.r.squared, lm1$coef[[2]])
      }
    }else{
      myoutput=ititle
    }

    myoutput
} # end of function

# The function VariableRelationship visualizes the relationships between the columns of a data frame.
# It is an MDS plot on the basis of the Spearman correlations between the variables
# the relationship between variables
VariableRelationshipPlot = function(dath, fname, colorcode,myimageType="png")
{
   globalfname = paste(fname, "_allmodules.",myimageType, sep="")
   openImgDev(globalfname)
   mds=cmdscale(1-abs(cor(dath, use="p",method="s")),k=2)
   plot(mds,type="n" ,main="Across Modules")
   text(mds,label=names(dath))
   dev.off()

   colorlevels = names(table(colorcode))
   for (icolor in  colorlevels){
       if (icolor=="grey"){
          next
       }
       sel = as.character(colorcode)==as.character(icolor)
       modulefname = paste(fname, "_", as.character(icolor), ".",myimageType,    sep="")
       mtitle      = paste("Inside ",  as.character(icolor), " Module", sep="")
       openImgDev(modulefname)
       plot(varclus(as.matrix(dath[sel,])), main=mtitle, xlab=mtitle)
       dev.off()
   }
}

#System        /GeneCategory      /ListHits/ListTotal/PopulationHits/PopulationTotal/EASE score
#BiologicalProc/mitotic cell cycle/54      /166      /305           /10642          /7.23e-042
# a=rbind( c(305-54, 10642-251-112-54), c(54, 112) )
# a
# signif(fisher.test(a)$p.value,2)
#[1] 4.1e-43
# signif(chisq.test(a)$p.value,2)
#
# pop_list_hits: in order of PopulationTotal, PopulationHits, List Total, List Hits
#             Category           Non-Category
#Population     a[1,1]             a[1,2]
#List           a[2,1]             a[2,2]
#
#
#************ example, Fisher exact Test & Hypergeometric test ********************
#
# 196 of 673 genes overlap the 948 genes, total geens 24000
#a=rbind( c(673-196, 24000-(673-196)-(948-196)-196), 
#          c(196,     948-196))
#
#  signif(fisher.test(a)$p.value,2)
# hypergeometric test: 
#    phyper(196, 673, 24000-673, 948, lower.tail=F)
#
# By chance you expect 196/23756*948 = 8 overlap
#
# if population hits < minPopulationHits, then the test doesn't continue
#
#input vector:=[population, population_hit, list, list_hit]
fisherTest = function(poplisthits, minPopulationHits=5)
{
   if ( (poplisthits[2]<minPopulationHits) | (poplisthits[4]==0) )
      return (1)

   #old way, maybe wrong
   #a=rbind( c(poplisthits[2]-poplisthits[4], poplisthits[1]-poplisthits[2]-poplisthits[3]+poplisthits[4]), 
   #         c(poplisthits[4],           poplisthits[3]-poplisthits[4] ) )
   #cat(as.character(poplisthits[1]), as.character(poplisthits[2]),as.character(poplisthits[3]),as.character(poplisthits[4]),"\n")
   #signif(fisher.test(a)$p.value,3)
   q= poplisthits[4] #list hit
   m= poplisthits[2] #population hit
   n= poplisthits[1]-poplisthits[2] #population non-hit
   k= poplisthits[3] #list total

   #myp=1-phyper(q-1, m, n, k)
   myp=phyper(q-1, m, n, k, lower.tail=F)
   signif(myp,3)
}

# 11 out of 14 overlap 38, total 2588
fisherExactTest=function(overlap, poolA, poolB, total){
   a=rbind( c(poolA-overlap, total- (poolA-overlap)-(poolB-overlap)-overlap),
            c(overlap,       poolB-overlap)
     )
   signif(fisher.test(a)$p.value,2)
}


# compute FET for one module against N categories 
#
# the first vector has two elements (population, population_hit) and overlaps(N)
#   listsize is the sizes of the N categories
#
fisherTestVect = function(PopHitListhits, listsizes){
   N = length(listsizes)
   minPopulationHits=1;
   res = rep(1, N)
   for (i in c(1:N) ) {
     myvect = c(PopHitListhits[c(1:2)], listsizes[i], PopHitListhits[i+2])
     pv = fisherTest(myvect, minPopulationHits=minPopulationHits)
     res[i] = pv
   }
   
   return (rbind(res))
}

multiplyVect = function(myvector){
  res = 1
  for (i in myvector) {
    res = res*i
  }
  return (res)
}

forcenumericSum=function(myvector){
   sum( as.numeric(as.matrix(myvector)), na.rm=T)
}


tapplySumFunc=function(myvector, mycolor){
   restab=tapply(as.numeric(as.matrix(myvector)), mycolor, sum)
   as.matrix(restab)
}

#imput is a boolean vector and another vector whose components are to be selected
# here we use a trick of temporary file
applySelectFunc=function(selvector0or1, myvector, removeduplicates=F){
   selvector= as.integer(as.matrix(selvector0or1)) & TRUE
   if(sum(selvector)<=0)
      return (" ")

   tmpfn="tmp.txt"

   selcomp = as.character(myvector[selvector])
   if (removeduplicates){
     selgenes = sort( as.character(myvector[selvector]) )
     no.s =length(selgenes)
     selgenes.shift = c("0000000", selgenes[c(1:(no.s-1))])
     unique.bool = (selgenes.shift != selgenes)
     selcomp  = selgenes[unique.bool]
   }

   write.table( t(as.matrix(selcomp)), tmpfn, sep=";",quote=FALSE, col.names=F, row.names=FALSE)
   mystring <- read.delim(tmpfn,sep="\t", header=F)
   as.character(as.matrix(mystring))
}



# The function fisherPvector allows one to compute Fisher exact p-values
# Thus it allows one to carry out an EASE analysis
# Output: a table of Fisher?s exact p-values
# Input: annotation1 is a vector of gene annotations
# Input: couleur (French for color) denotes the module color of each gene
# Only those gene functions (levels of annotation1) that occur a certain mininum number of times
# (parameter= minNumberAnnotation) in the data set will be considered.  
fisherPvectorOld=function(couleur, annotation1, minNumberAnnotation=10) {

    levelsannotation1 =levels(annotation1)
    levelscouleur     =levels(couleur)
    no.couleur        =length(levelscouleur)
    restannotation1   =table(annotation1)>minNumberAnnotation
    no.annotation     =sum( restannotation1)

    datout=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )

    names(datout)        =levelscouleur
    restlevelsannotation1= levelsannotation1[restannotation1]
    row.names(datout)    = restlevelsannotation1

    for (i in c(1:no.annotation) ) {
      for (j in c(1:no.couleur) ){
        tab1=table( annotation1==restlevelsannotation1[i], couleur==levelscouleur[j])
        datout[i,j]=signif(fisher.test(tab1)$p.value,2)
    }}
    datout

} # end of function fisherPvector


#outformat="pvalue":     output pvalue table
#outformat="ratio":      output tables of ratios of two proportions
#outformat="both":       output pvalue(proportion) table
fisherPvector =function(couleur, annotation1, outformat="pvalue", minNumberAnnotation=10) {

    levelsannotation1 =levels(annotation1)
    levelscouleur     =levels(couleur)
    no.couleur        =length(levelscouleur)
    restannotation1   =table(annotation1)>minNumberAnnotation
    no.annotation     =sum( restannotation1)

    datout=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )

    names(datout)        =levelscouleur
    restlevelsannotation1= levelsannotation1[restannotation1]
    row.names(datout)    = restlevelsannotation1

    annotation1.chars = as.character(annotation1)

    for (i in c(1:no.annotation) ) {
      for (j in c(1:no.couleur) ){

        bool.annoted= (annotation1.chars==restlevelsannotation1[i])
        bool.colorj = (couleur==levelscouleur[j])

        # some genes are not annoted, so we need remove those from our proportion computation
        bool.annoted.noNAs = bool.annoted[ !is.na(bool.annoted) ]
        bool.colorj.noNas  = bool.colorj[  !is.na(bool.annoted) ]

        tab1=table(bool.annoted.noNAs, bool.colorj.noNas)
        pvalue = signif(fisher.test(tab1)$p.value,2)

        #                    rest module
        #bool.annoted.noNAs FALSE TRUE
        # nofunc       FALSE 3257   292
        #   func       TRUE    29     1
        portion.rest   = tab1[2,1]/(tab1[1,1]+tab1[2,1]) # portion of 
        portion.colorj = tab1[2,2]/(tab1[1,2]+tab1[2,2]) #=length(bool.colorj.noNas)

        #proportion
        ratioOFportions   = portion.colorj/portion.rest
        
        if (outformat=="pvalue"){
            datout[i,j]=pvalue
        }else if( outformat=="ratio" ){
            datout[i,j]=signif(ratioOFportions,2)
        }else{
            myout = paste(as.character(pvalue),"(", 
                          as.character(signif(portion.colorj,2)),",",
                          as.character(signif(portion.rest,2)),  ")", sep="")
            datout[i,j]= myout
        }

    }}
    datout

} # end of function fisherPvector





# this function computes the standard error
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }

choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {out1[i] <- 0}
    else {out1[i] <- choose(n[i], k)}}
  out1
}


####  Function pamPsClNo computes prediction strength and returns the estimated number of clusters  
#### based on PAM clustering
                                        #k=6  at most k clusters
                                        #v=2 fold cross validation
                                        #c=5  no. of cross validations
pamPsClNo <- function(original.data,
                      klimit=5,  # max no. of clusters
                      cvFold=2,  # how many fold cross validation
                      cvCount=5, # how many cross validations
                      m=2,       # number of points
                      double1=FALSE, diss=FALSE,
                      cut1=1){
  if (double1) {original.data <- rbind(original.data,original.data) }
  clustlearn <- list(NULL);
  clusttest  <- list(NULL);
  nData <- nrow(original.data);
  cps1 <- matrix(1,nrow=klimit,ncol=cvCount)
  criterion1 <- matrix(1,nrow=klimit,ncol=cvCount) 
  if (diss) {
    alldist <- as.matrix(original.data)
  } else {
    alldist <- as.matrix(dist(original.data))
  }
  alldist <- signif(alldist,10)
  ## for each cross-validation set
  for (cc in 1:cvCount) {
    ## two matrices used to store prediction strength calculated by
    ## four kinds of metrics
    ps1 <- matrix(1, nrow=klimit, ncol=cvFold)

    ## utility vector indicating split of data set into cvFold sets 
    rest1 <- nData-as.integer(nData/cvFold)*cvFold
    sam1 <-  sample(c(rep(c(1:cvFold), as.integer(nData/cvFold)),
                      sample(1:cvFold, rest1)))

    ## for each possible number of clusters,
    for (kk in 2:klimit){
      ## cvFold fold splitting for cross validation
      for (vv in 1:cvFold){

        ## indices of test and training sets
        test.index <- c(1:nData)[sam1==vv]
        learn.index <- c(1:nData)[sam1!=vv]
        no.test <- length(test.index)
        no.learn <- length(learn.index)

        if (no.test <= kk || no.learn <= kk) {
          ## clustering too few points into too many clusters
          ps1[kk,vv] <- 0
          next
        }
        ## distances between points in test and training sets
        test.dist <- alldist[test.index, test.index]
        learn.dist <- alldist[learn.index, learn.index]

        ## perform clusterings on test and training sets
        ##print(paste("Performing pam with kk=", kk, "no.test", no.test, "no.learn", no.learn))
        clustlearn[[kk]] <- pam(as.dist(learn.dist), kk, diss=T)
        clusttest[[kk]] <-  pam(as.dist(test.dist), kk, diss=T)
        
        ## this assigns a cluster to each test set observation, based on the
        ## clustering of the training set.
        Cl <- rep(0, no.test)
        d <- rep(10000000, kk) #difference matrix for assigning Cl
        
        ## determine which medoid each test set point is closest to/i.e.,
        ## which cluster it belongs to
        index1 <-  clustlearn[[kk]]$medoids # length is kk
        for (i in 1:no.test){
          for (j in 1:kk){
            d[j] <- alldist[index1[[j]], test.index[i]]
            ## note: this assumes that the medoids are in the original dataset
          }
          mincluster <- c(1:kk)[rank(d) == min(rank(d))]
          if (length(mincluster) == 1) {
            Cl[i] <- mincluster
          }
          else if (length(mincluster)>1) {
            Cl[i] <- sample(mincluster, 1)
          }
        }  # end for for over i in 1:no.test
        
        ## now we compute how often m samples are co-clustered
        tab0 <- table(Cl, clusttest[[kk]]$clustering)
        pshelp <- rep(10000000, kk)
        for (l in 1:kk){
          ## marginals
          tab1=tab0
          tab1[tab0<cut1]=0
          nkl <- sum(tab1[,l])
          if (nkl < m)  {
            pshelp[l] <- 0
          } else { 
            pshelp[l] <- sum(choosenew( tab1[,l], m ) )/ choosenew(nkl,m)     
          }
        } # end of  for l in 1:kk
        ps1[kk,vv] <- min(pshelp)              # Min
      } # end of vv in 1:cvFold
    } # end of kk in 2:klimit
    cps1[,cc] <- apply(ps1,1,mean)
    ## gives max over vv for cvFold=2
    criterion1[,cc] <-  cps1[, cc] + apply(ps1, 1, stderr1)
  } # end of for cc in 1:cvCount
  psse <- signif(apply(criterion1, 1, mean), 3)
  ##psse <- signif(apply(cps1, 1, mean), 3)
  psse
  ##print(psse)
  ##thres <- 0.8-(m-2)*0.05 
  ##for(nn in 1:klimit) { 
  ##  if(psse[nn]>=thres) {clusterNo <- nn}
  ##}
  ##clusterNo
} # end of function pamPsClNo


kpredictPS2=function(ps,span1=0.75,title1="prediction strength") {
    kk=c(1:length(ps))
    kout=1
    if (length(kk)>1) { 
        lm1=loess(ps~ kk,span=span1,degree=2)    
        rank1=rank(round(resid( lm1 ),digits=10))
        kouthelp=  kk[   rank1==max(rank1) ]
        kestimate=kouthelp[length(kouthelp)]
        plot(ps,main=title1,xlab="no. of clusters k",ylab="PSE",ylim=c(0,1))
        lines(kk,predict(lm1))
        abline(v=kestimate,col="red")
        abline(h=.8)
        points(kk,resid(lm1) ,pch="*",col="blue")
        kestimate
    }
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##     Function: Make a plot, run a couple tests and output everything to proper files 
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

scatterCorrelTester=function(rawValA, rawValB, binaryvar=FALSE, PlotsTitle ="", myxlab="",myylab="",outputImage="tmp.png", outputFile=tests_Fname, numPoints=20 ){

        if( (sd(rawValA, na.rm = T)==0) | (sd(rawValB, na.rm = T)==0) ){
            titleinfor = paste(PlotsTitle, "error in scatterCorrelTester, sd()=0", paste="")
            appendListToFile(outputFile, titleinfor) 
            return
        }

	#in general, we will have the thing being tested be rawValA, and K be rawValB
	cut_factor = cut(rawValB,numPoints)

	# Now give me the mean proportion of essentiality in each chunk of the essentiality vector 
	mean_rawValA = tapply(rawValA, cut_factor, mean)

	# Now give me the mean proportion of k in each chunk of the essentiality vector 
	mean_rawValB = tapply(rawValB, cut_factor, mean)

	#now we need a p-value for this:
        corrl   = signif(cor(rawValB,rawValA),2)
        #if (binaryvar==FALSE){
           corrl.p = signif( cor.test(rawValB,rawValA, method="p",use="p")$p.value ,2)
        #}else{
	#   kruskal = kruskal.test(rawValB,rawValA)
        #   corrl.p = signif(kruskal$p.value, 2)
        #}
        titleinfor = paste(PlotsTitle, ": r=", as.character(corrl),",p=", as.character(corrl.p), paste="")

	#Now we need to just make a line fit
	line =lm(mean_rawValA ~ mean_rawValB)
	sum1 = summary(line)
	sum1

	#Now put in a title for this part of the data:
	appendListToFile(outputFile, titleinfor)
	appendListToFile(outputFile, sum1)

	# do a plot() of these two vectors against each other. 
	openImgDev(outputImage)
        par(mfrow=c(1,1), mar=c(5, 5, 4, 2) + 0.1)
	plot(mean_rawValB,mean_rawValA, main=titleinfor, xlab=myxlab, ylab=myylab)
	# and then draw the line
	abline(line, col="red")
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
	dev.off()
}

conservationVectorCleaner = function(rawVal, nasVal) {

	# I need a way to always filter out the NAs from the mean_e_val vector...
	no.nas = nasVal
	# variable log_mean_e_val must be defined for use in plotting conservation plots
	mean_e_valNoNA = rawVal[no.nas]
	# 2nd I need to replace the 0 values with extremely small numbers (smaller than you see so e-200)
	mean_e_valNoNANorZero = ifelse(mean_e_valNoNA==0, 10^(-200), mean_e_valNoNA)
	# 3rd I need to transform the e value since its range is just HUGE
	log(mean_e_valNoNANorZero)
	#R just returns the last "thing" listed in a function
}


# the following function computes the Rand index between 2 clusterings
# assesses how similar to clusterings are
choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {
        out1[i] <- 0
    }else {
        out1[i] <- choose(n[i], k)
    }
  }
  out1
}

Rand <- function(itab,adjust=T) {
  no.rows <- nrow(itab)
  no.cols <- ncol(itab)
    tab  =  matrix(as.numeric(as.matrix(itab)), no.rows, no.cols)

     a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
     n <- nrow(tab)
     
     for (i in 1:n) {
          for(j in 1:n) {
               a <- a+choosenew(tab[i,j],2)
               nj <- sum(tab[,j])
               c <- c+choosenew(nj,2)
          }
          ni <- sum(tab[i,])
          b <- b+choosenew(ni,2)
          nn <- nn+ni
     }

     if(adjust==T) {
          d <- choosenew(nn,2)
          adrand <- (a-(b*c/n)/d)/(0.5*(b+c/n)-(b*c/n)/d)
          adrand
     } else {
          b <- b-a
          c <- c/n-a
          d <- choosenew(nn,2)-a-b-c
          rand <- (a+d)/(a+b+c+d)
          rand
     }
}

combinBy2=function(a)
{
  a*(a-1)*0.5
}

#based on the paper: http://faculty.washington.edu/kayee/pca/supp.pdf
RandIndex= function(tab)
{
  no.rows <- nrow(tab)
  no.cols <- ncol(tab)
  
  A  =  matrix(as.numeric(as.matrix(tab)), no.rows, no.cols)
  ni <- apply(A, 1, sum) #sum over columns for each row
  nj <- apply(A, 2, sum) #sum over rows for each column
  n  =  sum(ni)

  n.comb2   = combinBy2(n)
  ni.comb2  = combinBy2(ni)
  nj.comb2  = combinBy2(nj)
  nij.comb2 = combinBy2(A)

  sum.ni.comb2 = sum(ni.comb2)
  sum.nj.comb2 = sum(nj.comb2)

  sum.nij = sum( apply(nij.comb2, 1, sum) )
  
  a= sum.nij - sum.ni.comb2*sum.nj.comb2/n.comb2
  b= 0.5*(sum.ni.comb2+sum.nj.comb2) - sum.ni.comb2*sum.nj.comb2/n.comb2
  rand=a/b
  abs(rand)
}


# This function computes the GTOMm DISSIMILARITY
#
# Input:
# - adjmat1, a symmetric adjacency matrix with binary entries
# - m, the order of GTOM
#
# Output:
# - The GTOMm dissimilarity matrix
#
# Andy M. Yip and Steve Horvath
# January 2005

#if(exists("GTOMmdist1")) rm(GTOMmdist1);
GTOMmdist1 = function(adjmat1,m=1){
    if (m!=round(abs(m))){
        stop("m must be a positive integer!!!", call.=TRUE);}
    if (any(adjmat1!=0 & adjmat1!=1)){
        stop("The adjacency matrix must be binary!!!", call.=TRUE);}

    B <- adjmat1;
    if (m>=2) {
        for (i in 2:m) {
            diag(B) <- diag(B) + 1;
            B = B %*% adjmat1;}}   # number of paths with length at most m connecting each pair
    B <- (B>0);                    # m-step reachability matrix
    diag(B) <- 0;                  # exclude each node being its own neighbor
    B <- B %*% B;                  # number of common k-step neighbors that each pair of nodes share

    Nk <- diag(B);                 # number of common k-step neighbors that each node possesses
    B <- B +adjmat1;
    diag(B) <- 1;
    denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
    diag(denomTOM) <- 1;
    1 - B/denomTOM                 # turn the GTOM matrix into a dissimilarity
}

#ModulePrincComps, which computes the first principal component of each module. 
#Also it provides the variance explained for the first 5 PCs of a module.
#It is based on the singular value decomposition.
#It takes as input datExpr  for which the rows are samples and the columns are genes.
#Here is how you would use it
#PC1=ModulePrinComps2[datExpr,color1]
#Then PC1[[1]] provides a data frame with the first PC of each module
#PC1[[2]] provides a data frame where each column contains the percentages of 
#the variance explained by the first 10 PCs by module.

ModulePrinComps.old = function(datexpr,couleur) {
  no.pcs=10

  modlevels=levels(factor(couleur))

  PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
  PrinComps2=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
  varexplained= data.frame(matrix(666,nrow=no.pcs,ncol= length(modlevels)))

  colnames(PrinComps)    <- modlevels
  colnames(PrinComps2)   <- modlevels
  colnames(varexplained) <- modlevels

  for(i in c(1:length(modlevels)) ){
    #print(i)   
    modulename    = modlevels[i]
    restrict1= as.character(couleur)== modulename
    # in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])
    datModule=impute.knn(as.matrix(datModule))$data

    datModule=t(scale(t(datModule)))
    svd1=svd(datModule)
    mtitle=paste("PCs of ", modulename," module", sep="")
    varexplained[,i]= (svd1$d[1:no.pcs])^2/sum(svd1$d^2)

    # this is the first principal component
    pc1=svd1$v[,1]
    signh1=sign(sum(cor(pc1,  t(datModule))))
    if (signh1 != 0)  pc1=signh1* pc1

    # this is the second principal component
    pc2=svd1$v[,2]
    signh2=sign(sum(cor(pc2,  t(datModule))))
    if (signh2 != 0)  pc2=signh2* pc2

    PrinComps[,i]= pc1
    PrinComps2[,i]= pc2
  }

  list(PrinComps, PrinComps2, varexplained)
}

#============================================================
# adj is the input adjacency matrix. It should be symmetric and all its
#diagonal elements should be 1. Do NOT assign zeros to the diagonal entries
#as before!
#
# kpred=T will output the predicted within-module connectivity.
#
# Umat=T will output the U matrix in the Spectral Decomposition. 
# CAUTION: U is an N*N matrix which will eat up a lot of memory and CPU when N is big
#  (such as N=3000).
#
# D=T will output the eigenvalues of the matrix adj.
#
# In default, the output includes Module Size, Module Cohesiveness, Weight,
#   Module Conformity and Within-module connectivity.

SDADJ1 = function(inadj, kPred=F, Umat=F, D=F) 
{
    # Check if adj is a valid adjacency matrix:  square matrix, positive
    #entries, symmetric and diagonal elements being 1.
    adj=inadj
    diag(adj) <- 1

    n=dim(adj)
    tol=0.000000001

    if ( n[1] != n[2])               stop("The adjacency matrix is not a square matrix!")
    if ( sum(is.na(adj))>0 )         stop("There are missing values in the adjacency matrix!")
    if ( sum(adj<0)>0 )              stop("There are negative entries in the adjacency matrix!")
    if ( max(abs(adj-t(adj))) > tol) stop("The adjacency matrix is not symmetric!")
    if ( max(abs(diag(adj)-1))> tol) warning("The diagonal elements are not all one!")

    # Calculate kWithin, MCH, MCF, and/or kPred, Umat, D
    n=n[1]
    sd=eigen(adj) # Spectral Decomposition

    kWithin=apply(adj, 2, sum)                # Within Module Connectivities
    MCohesiveness=sd$values[1]/sum(sd$values) # Module Cohesiveness

    MConformity=abs(sd$vectors[,1])*sqrt(n) # Module Conformity

    Weight= sqrt(MCohesiveness*sum(kWithin))

    output=list(MSize=n,MCohesiveness=MCohesiveness,Weight=Weight,MConformity=MConformity,kWithin=kWithin)

    if (kPred) {
      kPred= Weight* MConformity
      output$kPred=kPred
    }

    if (Umat) output$Umat=sd$vectors
    if (D)    output$D=sd$values
    #output

    MConformity
}
#=======================================================================================================================
#=======================================================================================================================
#=======================================================================================================================

# whichway="row" : cluster by "column","row", or "rowcolumn"
# clusterno.row  : <=0 no row-wise cluster detection, >0 identify only the specified number of clusters
# clusterno.col  : <=0 no col-wise cluster detection, >0 identify only the specified number of clusters

#************** Hierarchical Clustering  ********************************************

## whichway = "row"/"column"/"rowcolumn"
#methods:  "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
hiercluster2way=function(datMatrix, similarity="cor",  whichway="row",clusterno.row=2, clusterno.col=3, imagename=NULL,
                         showlabels=FALSE, clustmethod="average", scale="row",maxlabels= 100, 
                         iwidth = 1024, iheight = 1024, ipointsize = 12,  margin=c(8, 4, 4, 2), cex=0.9, mcexrow=1, ires=300,
                         heatmapcolor=rgcolors.func(50)){
  mno.cols = dim(datMatrix)[2]
  mno.rows = dim(datMatrix)[1]
  #mcex=0.9
  cormatrix = NULL

  iint = as.integer(mno.cols/50) + 1
  mcex = cex*(0.9^iint)

  iint = as.integer(mno.rows/50) + 1
  mcexRow = cex*(0.98^iint)

  if (!is.null(imagename) ){ #get file names for column and row dendrograms
       generalname    = getFileName(imagename)
       extensionname  = getFileExtension(imagename)
       myname.row = paste(generalname, "_dendrogramRow.", extensionname, sep="")
       myname.col = paste(generalname, "_dendrogramColumn.", extensionname, sep="")
  }
  
  if( (whichway=="column") | (whichway=="rowcolumn") ){
     #--------- compute correlation coefficient 
     if(similarity =="binary") {
        cormatrix <- similarity_binary(datMatrix)
     } else{
        cormatrix <- cor(datMatrix, use = "pairwise.complete.obs") 
     }

     cormatrix <- ifelse(is.na(cormatrix), 0, cormatrix)
     #--------- hierarchical clustering
     distmatrix = 1-cormatrix
     colnames(distmatrix) = colnames(datMatrix)
     rownames(distmatrix) = colnames(datMatrix)

     hierclu.col <- hclust(as.dist(distmatrix), method= clustmethod) 

     # we expect two or three clusters of samples from  group comparison
     if (clusterno.col>0){
       colcolors = moduleDetectByFixedClusterno(ihcluster=hierclu.col, nocluster=clusterno.col)
     }else{
       colcolors   = NA#rep("white", mno.cols)
     }

     if (!is.null(imagename) ) {
         if(mno.cols<50){
           openImgDev(myname.col)
           par(mfrow=c(1,1), cex=1.1)
         } else{
           openImgDev(myname.col, iwidth =iwidth, iheight = iheight, ipointsize = ipointsize, ires=ires)
           par(mfrow=c(1,1), cex=mcex)
         }
     }
     if ((showlabels) && (mno.cols<= maxlabels)){
        plot(hierclu.col, xlab="",ylab="",main="",sub="")
     }else{
        plot(hierclu.col, labels=F, xlab="",ylab="",main="",sub="")
     }
     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
     if (!is.null(imagename) )
          dev.off()

  }else {
     hierclu.col = NA
     colcolors   = NA#rep("white", mno.cols)
  }

  if(whichway=="row" | (whichway=="rowcolumn") ) {

     if(similarity =="binary") {
        cormatrix <- similarity_binary(t(datMatrix))
     } else{
        cormatrix <- cor(t(datMatrix), use = "pairwise.complete.obs")
     }

     cormatrix <- ifelse(is.na(cormatrix), 0, cormatrix)

     distmatrix = 1-cormatrix

     colnames(distmatrix) = rownames(datMatrix)
     rownames(distmatrix) = rownames(datMatrix)

     hierclu.row <- hclust(as.dist(distmatrix), method= clustmethod) 

     # we expect exactly two  clusters of genes from  group comparison
     if (clusterno.row>0){
       rowcolors = moduleDetectByFixedClusterno(ihcluster=hierclu.row, nocluster=clusterno.row)
     }else{
       rowcolors   = NA#rep("white", mno.rows)
     }
   
    if (!is.null(imagename) ) {
         if(mno.rows<50){
            openImgDev(myname.row)
            par(mfrow=c(1,1), cex=1.1)
         }else{
            openImgDev(myname.row, iwidth =iwidth, iheight = iheight, ipointsize = ipointsize, ires=ires)
            par(mfrow=c(1,1), cex=mcex)
         }
    }

     if ( (showlabels) && (mno.rows<=maxlabels)){
        plot(hierclu.row, xlab="",ylab="",main="",sub="")
     }else{
        plot(hierclu.row, labels=F, xlab="",ylab="",main="",sub="")
     }
     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    if (!is.null(imagename) ) 
          dev.off()
  }else {
     hierclu.row = NA
     rowcolors   = NA#rep("white", mno.rows)
  }

  #heatmap(as.matrix(datMatrix), Rowv=as.dendrogram(hierclu.row), Colv=as.dendrogram(hierclu.col), 
  #           RowSideColors=rowcolors, ColSideColors=colcolors,
  #           scale="row", revC=F, xlab="", col = rgcolors.func(50))

  #cat(colcolors, "\n")

  plotTwowayDendrogramsOnExprArrays(exprarray=datMatrix, hiercol=hierclu.col, hierrow=hierclu.row, 
                                              columnmodules=colcolors, rowmodules=rowcolors, imagename=imagename, 
                                              scale_default=scale, 
                                              iwidth =iwidth, iheight = iheight, ipointsize = ipointsize,
                                              margin=margin, cex0=cex,  cexrow=  mcexrow, heatmapcolor=heatmapcolor)

  if(!is.null(cormatrix)){
     rm(cormatrix)
     rm(distmatrix)
     collect_garbage()
  }

  list(rowcolors, colcolors)
}


plotTwowayDendrogramsOnExprArrays = function(exprarray, hiercol=NA, hierrow=NA, 
                             columnmodules=NA, rowmodules=NA, imagename=NA, 
                             scale_default="row", iwidth = 1024, iheight = 1024, 
                             ipointsize = 12, margin=c(10, 4, 4, 5), cex0=1, cexrow=0.8, 
                             ires=300, heatmapcolor=rgcolors.func(50))
{
  max.elements = 2500
  no.cols = dim(exprarray)[2]
  no.rows = dim(exprarray)[1]

  iint = as.integer(no.cols/50) + 1
  mcex = cex0*(0.9^iint)

  if(is.null(cexrow)){
  iint = as.integer(no.rows/50) + 1
  if (no.rows<50) {
    mcexrow = 1.1
  }else{
    mcexrow = cex0*(0.9^iint)
  }
  } else{
    mcexrow = cexrow
  }


  if (!is.na(imagename)){
     openImgDev(imagename, iwidth =iwidth, iheight =iheight, ipointsize =ipointsize)
     par(mfrow=c(1,1), mar=margin + 0.1, cex=cex0)
  }

  orderedMatrix = as.matrix(exprarray)
  if( no.rows<=max.elements && !is.na(hierrow) ){ #normal case
     hierrow4disp  = as.dendrogram(hierrow) #have to use as.dendrogram(..)
     rowmodulesZ = as.character(rowmodules)
  }else if( no.rows>max.elements && !is.na(hierrow) ){
     #too many genes,so we order the arrays, but don't draw the dendrogram
     orderedMatrix = orderedMatrix[hierrow$order,]
     hierrow4disp  = NA
     rowmodulesZ   = as.character(rowmodules[hierrow$order])
  }else{
     hierrow4disp  = NA
     rowmodulesZ = rowmodules
  }

  if( (no.cols <= max.elements) && (!is.na(hiercol)) ){
     hiercol4disp=as.dendrogram(hiercol)
     columnmodulesZ = as.character(columnmodules)
  }else if( (no.cols >max.elements) && (!is.na(hiercol)) ){
     orderedMatrix = orderedMatrix[, hiercol$order]
     hiercol4disp  = NA
     columnmodulesZ = as.character(columnmodules[hiercol$order])
  }else{
     hiercol4disp  = NA
     columnmodulesZ = as.character(columnmodules)
  }

  #mycexCol = 0.2 + 1/log10(nc)

  
  if( !(is.na(columnmodules)) & !(is.na(rowmodules)) ) {
  heatmap(orderedMatrix, Rowv=hierrow4disp, Colv=hiercol4disp, 
             RowSideColors=rowmodulesZ, ColSideColors=as.character(columnmodulesZ),
             margins = margin[c(1,4)], cexCol = mcex,cexRow=mcexrow,
             scale=scale_default, revC=F, xlab="", col = heatmapcolor)

  }else if( is.na(columnmodules) & !(is.na(rowmodules)) ) {
  heatmap(orderedMatrix, Rowv=hierrow4disp, Colv=hiercol4disp, 
             RowSideColors=rowmodulesZ,
             margins = margin[c(1,4)],cexCol = mcex,cexRow=mcexrow,
             scale=scale_default, revC=F, xlab="", col = heatmapcolor)

  }else if( is.na(rowmodules) & !(is.na(columnmodules)) ) {
  heatmap(orderedMatrix, Rowv=hierrow4disp, Colv=hiercol4disp, 
             ColSideColors=as.character(columnmodulesZ),
             margins =margin[c(1,4)],cexCol = mcex,cexRow=mcexrow,
             scale=scale_default, revC=F, xlab="", col = heatmapcolor)

  }else{
  heatmap(orderedMatrix, Rowv=hierrow4disp, Colv=hiercol4disp, 
             margins =margin[c(1,4)],cexCol = mcex,cexRow=mcexrow,
             scale=scale_default, revC=F, xlab="", col = heatmapcolor)

  }

  if (!is.na(imagename)){
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off()
  }

}

##################################### EASE ANALYSIS ###############################################################
#

EASEOntologyEnrichmentAnalysis=function(genesInfor, ontologyfn, fname, maxSignifLevel=0.05,myMinPopulationHits=3, OntologyType="Ontology", background=0)
{
    #ontologyname = getOntologyNameFromPath(ontologyfn)

    cols.geneinfor <- dim(genesInfor)[2]

    modulecolor = as.character(genesInfor$module)
    ctable = table(modulecolor )
    #ctable

    #*-------------------------------------------------------------------------------------
    #* STEP 1: read in ontology information, columns as ontology category
    #*
    #*FORMAT:
    #*  System	             GeneCategory	                     PopulationHits	PopulationTotal   LocusLinkNumbers
    #* GOMolecularFunction   ubiquitin C-terminal hydrolase activity	39	        11229             23123; 10234; 33342; 
    #*
    #*
    #*-------------------------------------------------------------------------------------
    ontologyMatrix.all <- read.delim(ontologyfn, sep="\t", header=T)
    ontologyMatrix.all <- as.matrix(ontologyMatrix.all)

    #get unique annotated genes
    allgogenes = NULL
    for (each in as.character(ontologyMatrix.all[, 5]) ){
         #each_lowercase= tolower(each)
         each_lowercase= toupper(each)
         llids=unlist( strsplit( as.character(each_lowercase), "; ") )
         allgogenes = c(allgogenes, llids)
    }
    length(allgogenes )
    #use old trick to get unique LLIDs based on shift operations
    allgogenes.sorted = sort(allgogenes )
    allgogenes.shift  = c(-1, allgogenes.sorted[1:(length(allgogenes)-1)] )
    uniques = (allgogenes.sorted != allgogenes.shift)

    annotatedgenes = allgogenes.sorted[uniques]
    no.annotatedgenes = length(annotatedgenes)
    length(annotatedgenes)
    no.categories = dim(ontologyMatrix.all)[1]

    # population/background
    #
    if(background>0) {
       ontologyMatrix.all[,4]=as.character(background)
    }else{
       ontologyMatrix.all[,4]=as.character(no.annotatedgenes)
    }

    modulenames      = names(ctable)
    no.modules       = length(modulenames)

    #eventualColnames=c("System","Gene Category", "List Hits", "List Total", 
    #                   "Population Hits", "Population Total", 
    eventualColnames=c("System","GeneCategory", "ModuleOverlap", "GeneSetCount", 
                       "PopulationOverlap", "Population", 
                       "pvalue_corrected", "Fisher_pvalue", "Probe Set", "Gene Symbol") 

    #process module by module
    for (eachmodule in modulenames){

       cat("Module ", eachmodule, "\n")
       modulesel    = modulecolor==eachmodule
       modulesize   = sum(modulesel)
       modgeneinfor = genesInfor[modulesel, ]
       #modgeneinfor[,1] = toupper(as.character(modgeneinfor[,1]))

       uidx = findUniqueIdx(modgeneinfor[,1])
       modgeneinfor = modgeneinfor[uidx, ]
       dim(modgeneinfor)

       #find module genes annotated 
       
       mergedmodule = merge(annotatedgenes, modgeneinfor, by.x=1, by.y=1, sort=F,all=FALSE)

       #cat(mergedmodule,"\n")

       moduletotals = dim(mergedmodule )[1]

       #write column names
       outfilename  = paste(fname, "_", OntologyType, "_",  eachmodule, ".xls",  sep='')
       #write.table(t(as.matrix(eventualColnames)), outfilename, sep="\t", quote=FALSE, col.names=F, row.names=FALSE, append=F)
       
       if(background>0) { # use actual module size and given background
          a=apply(ontologyMatrix.all, 1, EASE4module, 
             modulegeneInfor=modgeneinfor, ModuleListTotals=modulesize, resultfname=outfilename, 
             maxSignifLevel=maxSignifLevel, myMinPopHits=myMinPopulationHits, no.goTerms=no.categories)

       }else{ # use actual matched genes in module and all annotated genes as background

          a=apply(ontologyMatrix.all, 1, EASE4module, 
             modulegeneInfor=modgeneinfor, ModuleListTotals=moduletotals, resultfname=outfilename, 
             maxSignifLevel=maxSignifLevel, myMinPopHits=myMinPopulationHits, no.goTerms=no.categories)
       }

       # sort the table by pvalue in ascendant order
       mymoduleGOMatrix <- read.delim(outfilename, sep="\t", header=T)
       pvalue.order     <- order(as.numeric(as.matrix(mymoduleGOMatrix$Fisher_pvalue)) )
       
       if( dim(mymoduleGOMatrix)[1]==1) {
          mytitle <- colnames(mymoduleGOMatrix)
          mymoduleGOMatrix <- as.matrix(mymoduleGOMatrix )
          mymoduleGOMatrix <- rbind(mymoduleGOMatrix)
          mymoduleGOMatrix[,4] = rep(as.character(moduletotals), dim(mymoduleGOMatrix)[1])
          fm = mymoduleGOMatrix
          fm2= cbind( rep(modulesize, dim(fm)[1]), fm)
          colnames(fm2) = c("ModuleSize", colnames(mymoduleGOMatrix))
       }else{
          mymoduleGOMatrix <- as.matrix(mymoduleGOMatrix )
          mymoduleGOMatrix[,4] = rep(as.character(moduletotals), dim(mymoduleGOMatrix)[1])
          fm = mymoduleGOMatrix[pvalue.order,]
          #print(fm)
          fm2= cbind( rep(modulesize, dim(fm)[1]), fm)
          colnames(fm2) = c("ModuleSize", colnames(fm))
       }

       write.table(fm2, outfilename, sep="\t", quote=FALSE, col.names=T, row.names=FALSE, append=F)
    }
}

EASE4module=function(OneGOcategoryInfor, modulegeneInfor, ModuleListTotals, resultfname, 
                     maxSignifLevel=0.05, myMinPopHits=3, no.goTerms=1)
{

    #*-------------------------------------------------------------------------------------
    #*OneGOcategoryInfor FORMAT:
    #*
    #*  System	             GeneCategory	                     PopulationHits	PopulationTotal   LocusLinkNumbers
    #* GOMolecularFunction   ubiquitin C-terminal hydrolase activity	39	        11229             23123; 10234; 33342; 
    #*
    #*
    #*-------------------------------------------------------------------------------------
    go.fields = length(OneGOcategoryInfor)

    # ++++++++++ number of annotated genes for this category ++++++++++
    System           <- as.character(as.matrix(OneGOcategoryInfor[1]))
    GeneCategory     <- as.character(as.matrix(OneGOcategoryInfor[2]))
    PopulationHits   <- as.integer(as.matrix(OneGOcategoryInfor[3]))
    population_total <- as.integer(as.matrix(OneGOcategoryInfor[4]))

    #cat(System, GeneCategory, "\n")

    goLLIDsA = as.character(as.matrix(OneGOcategoryInfor[5]))
    #goLLIDs = tolower(goLLIDsA)
    goLLIDs  = toupper(goLLIDsA)

    ontologyGenes = unlist( strsplit(goLLIDs, "; ") )

    #case 1. blank
    if( length(ontologyGenes )<1 ){
        print(OneGOcategoryInfor)
        cat("no genes gtt splitted at all\n")
        return (0)
    }

    #*-------------------------------------------------------------------------------------
    #* STEP 2: EASE ANALYSIS of Ontology 
    #*         by merge module genes with annotated genes
    #*         
    #*        
    #* 
    #*-------------------------------------------------------------------------------------
    # Here, we avoid the merge of the enormous matrix, instead we make a matrix with geneIds and 
    # their corresponding ROW indices in the big matrix, then perform merge operation
    mergedAnnotGeneInfor  = merge(modulegeneInfor, ontologyGenes, by.x=1, by.y=1, sort=F,all=FALSE)

    dim(mergedAnnotGeneInfor)

    # ++++++++++ number of annotated genes for each category ++++++++++ 
    ModuleListHits <- dim(mergedAnnotGeneInfor)[1]

    if ( (ModuleListHits ==0) | (PopulationHits<myMinPopHits) )
       return (-1)

    #eventualColnames=c("System","Gene Category", "List Hits", "List Total", 
    #                   "Population Hits", "Population Total", 
    #                   "Pvalue_FisherExactTest","Gene Symbol","Unique ID") 

    # make a big matrix including all number required for computing pvalue
    fishervector = c(population_total, PopulationHits,
                     ModuleListTotals, ModuleListHits)

    #cat(population_total, PopulationHits, ModuleListTotals, ModuleListHits, "\n")

    mypvalue     = fisherTest(fishervector, minPopulationHits=myMinPopHits)

    # The corrected p-values represent the Bonferroni-corrected p-values 
    #    (nominal p-value multiplied by the number of GO categories searched)
    #
    mypvalue.corrected = mypvalue*no.goTerms

    if ( is.na(mypvalue) )
       return (-2)

    if(mypvalue > maxSignifLevel)
       return (-3)

    if(mypvalue.corrected>1) { mypvalue.corrected=1}

    # Concatenate selected probset names, eg, "1367659_s_at;1367659_s_at;1367659_s_at;1367659_s_at"
    # genes are row-wise, singificant function categories are column-wise
    #
    sel.genes       = rep(T, ModuleListHits)
    signif.probsets = applySelectFunc(sel.genes, as.character(mergedAnnotGeneInfor[,1]) )
    signif.symbols  = applySelectFunc(sel.genes, as.character(mergedAnnotGeneInfor[,2]) )
    signif.mypvalues= as.character(signif(mypvalue, 2))
    signif.pvalcorrect= as.character(signif(mypvalue.corrected, 2))

    sel_rslt = c(System,GeneCategory,    ModuleListHits, ModuleListTotals,
                       PopulationHits,   population_total, 
                       signif.pvalcorrect, signif.mypvalues, 
                       signif.symbols,   signif.probsets)

    write.table( t(as.matrix(sel_rslt)), resultfname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE, append=T)
 
    return (1)
}


###################################################################################################################

OntologyEnrichmentAnalysis=function(genesInfor, ontologyfn, fname, maxSignifLevel=0.05,myMinPopulationHits=3,OntologyType="Ontology")
{
    ontologyname = getOntologyNameFromPath(ontologyfn)

    cols.geneinfor <- dim(genesInfor)[2]

    color1 = genesInfor$module
    ctable = table(color1)
    #ctable

    #*-------------------------------------------------------------------------------------
    #* STEP 1: read in ontology information, columns as ontology category
    #*
    #*
    #*-------------------------------------------------------------------------------------
    ontologyMatrix.all <- read.delim(ontologyfn, sep="\t", header=T)
    ontologyGenes      <- as.character(ontologyMatrix.all[,2])
    dim(ontologyMatrix.all)

    # total number of categories
    no.categories       <- ( dim(ontologyMatrix.all)[1] - 1 )

    # ++++++++++ number of annotated genes for each category ++++++++++ 
    # population/background
    #
    if(background>0) {
       population_total =  background
    }else{
       population_total =  dim(ontologyMatrix.all)[1]
    }

    PopulationHits   <- apply(ontologyMatrix.all[, c(2:(no.categories+1))], 2, sum ) #vector


    #*-------------------------------------------------------------------------------------
    #* STEP 2: EASE ANALYSIS of Ontology 
    #*           by merge module genes with annotated genes
    #*         
    #*        
    #* 
    #*-------------------------------------------------------------------------------------
    # Here, we avoid the merge of the enormous matrix, instead we make a matrix with geneIds and 
    # their corresponding ROW indices in the big matrix, then perform merge operation
    #
    #mergeMatrix.all = merge(genesInfor, ontologyMatrix.all, by.x=1, by.y=1, sort=F,all=FALSE)
    #dim(mergeMatrix.all)
    no.ontogenes = dim(ontologyMatrix.all)[1]
    tmpMatrix    = cbind(as.matrix(ontologyMatrix.all[,1]), c(1:no.ontogenes) )
    colnames(tmpMatrix) <- c("gene_id", "index")

    tmpmergeMatrix = merge(genesInfor, tmpMatrix, by.x=1, by.y=1, sort=F,all=FALSE)
    dim(tmpmergeMatrix)

    selectedIndices = as.integer(as.matrix(tmpmergeMatrix$index))

    mergeMatrix.all = cbind(as.matrix(tmpmergeMatrix), as.matrix(ontologyMatrix.all[selectedIndices, -c(1)]) )
    mergeMatrix.all = as.data.frame(mergeMatrix.all)
    dim(mergeMatrix.all)

    cat("  finished merge, sizes (rows, cols): ", as.character(dim(mergeMatrix.all)[1]), as.character(dim(mergeMatrix.all)[2]), "\n")

    ordered.genes      = order(as.character(as.matrix(mergeMatrix.all$module) ))
    orderedMergedMatrix <- mergeMatrix.all[ordered.genes,]

    #mergedAnnotGeneInfor <- orderedMergedMatrix[,  c(1, gene_symbolcol,cols.geneinfor, cols.geneinfor+1)]
    mergedAnnotGeneInfor <- orderedMergedMatrix[,  c(1:(cols.geneinfor+1))]
    mergedAnnotMatrix    <- orderedMergedMatrix[, -c(1:(cols.geneinfor+1))]
    dim(mergedAnnotMatrix)

    cols.mergedinfo = dim(mergedAnnotGeneInfor)[2]

    # ++++++++++ number of annotated genes for each category ++++++++++ 
    OverallListTotal = dim(mergedAnnotMatrix)[1]
    OverallListHits <- apply(mergedAnnotMatrix, 2, forcenumericSum) #vector across all categorries

    color1.merged    = as.character(mergedAnnotGeneInfor$module)
    ModuleListTotals = table(color1.merged)
    modulenames      = names(ModuleListTotals)
    no.modules       = length(modulenames)

    # ++++++++++ 
    # we count how many genes in each module have certain function category
    ModuleListHits = apply(mergedAnnotMatrix, 2, tapplySumFunc, color1.merged)

    if (!is.matrix(ModuleListHits)){
       ModuleListHits = t(as.matrix(ModuleListHits))       
    }

    #print(modulenames)
    #print(ModuleListHits)
    #print(dim(ModuleListHits))

    rownames(ModuleListHits) <- modulenames

    eventualColnames=c("System","Gene Category", "List Hits", "List Total", 
                       "Population Hits", "Population Total", 
                       "Pvalue_FisherExactTest","Gene Symbol","Unique ID") 

    for (imod in c(1:no.modules) ){
       
       cat(modulenames[imod],"\n")

       # make a big matrix including all number required for computing pvalue
       fisherMatrix = rbind(rep(population_total, no.categories), 
                                PopulationHits,
                                rep(ModuleListTotals[imod], no.categories),
                                ModuleListHits[imod, ])

       mypvalues     = apply(fisherMatrix, 2, fisherTest, minPopulationHits=myMinPopulationHits)
       mypvalues.sel = mypvalues <= maxSignifLevel #for column selection

       if (sum(mypvalues.sel) ==0)
           next

       # for row selection
       mymodule      = modulenames[imod]
       mod.sel       = as.character(mergedAnnotGeneInfor[,cols.geneinfor])==mymodule

       #annotated genes as a significant group for each function category
       geneSelMatrix = ( (mergedAnnotMatrix==1) * mod.sel)[, mypvalues.sel]
       geneSelMatrix = as.matrix(geneSelMatrix )

       # Concatenate selected probset names, eg, "1367659_s_at;1367659_s_at;1367659_s_at;1367659_s_at"
       # genes are row-wise, singificant function categories are column-wise
       #
       signif.probsets = apply(geneSelMatrix, 2, applySelectFunc, as.character(mergedAnnotGeneInfor[,1]))
       signif.symbols  = apply(geneSelMatrix, 2, applySelectFunc,  as.character(mergedAnnotGeneInfor[,2]))

       signif.ModuleListHits= ModuleListHits[imod, mypvalues.sel]
       signif.PopulationHits=PopulationHits[mypvalues.sel]
       signif.mypvalues     = mypvalues[mypvalues.sel]

       finaltestmatrix=NULL
       categorynames = names(signif.mypvalues)
       for (jc in c(1:length(categorynames)) ){
          sel_rslt = c(ontologyname,categorynames[jc], signif.ModuleListHits[jc], ModuleListTotals[imod],
                       signif.PopulationHits[jc], population_total, 
                       signif.mypvalues[jc], signif.symbols[jc], signif.probsets[jc])
          finaltestmatrix = rbind(finaltestmatrix, sel_rslt)
       }

       colnames(finaltestmatrix ) <- eventualColnames
       orderpvalue = order(signif.mypvalues)

       #zfinaltestmatrix = as.matrix(rbind(eventualColnames, finaltestmatrix))
       if (length(orderpvalue)>1 ){
         zfinaltestmatrix=finaltestmatrix[orderpvalue, ]
       }else{
         zfinaltestmatrix=t(finaltestmatrix[orderpvalue, ])
       }

       testoutputfname = paste(fname, "_", OntologyType, "_",  modulenames[imod], ".xls",  sep='')
       write.table(zfinaltestmatrix, testoutputfname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

}

##------------------------------------------------------------------------------------------------------------------
## inputfname contains: gene information, expression data, module information (last column) 
## identifier_col     : unique gene identifier column
##                      =1 for  normal case (probeset)
##                      =3 for bxh (locus link number)
##
## gene_symbolcol     : we assume the first col as probset and this column gives the real gene names
## ontologyfnlist     : lists of GeneFisher ontology files 
## maxSignifLevel     : report only the categories with FisherExactTest Pvalue < maxSignifLevel
## useEASEannotation  : use Affymetrix or EASE annotation file
## outdir             : =="", put the results under the same directory as the input, otherwise under the new directory
##
##  OntologyType="Ontology"/"TF"/"QTL"
##
OntologyAnalysisDull=function(inputfname, identifier_col=1, gene_symbolcol=1, ontologyfnlist, signifLevel=0.05, minPopHits=3, 
                              useEASEannotation=T, useAllModules=F,OntologyType="Ontology", outdir="", ctopNrows=2, background=0,
                              removeIndividualFiles=TRUE)
{

    #* STEP 0: read in gene information, expression data, module information

    allMatrix <- read.delim(inputfname,sep="\t", header=T)
    #attach(allMatrix)
    dim(allMatrix)

    if(outdir=="") {
       mfname        = getFileName(inputfname)
    } else{
       mfname        = getFileNameNopath(inputfname)
       mfname        = paste(outdir, "/", mfname, sep="")
    }

    no.cols  = dim(allMatrix)[2]

    if (!useAllModules){ # consider modules
       inforCols = c(identifier_col, gene_symbolcol, no.cols)
       mgenesInfor <- allMatrix[,inforCols]
       colnames(mgenesInfor) <- c("ID", "GeneSymbol", "module")

    }else{ # consider all gene in the list
       inforCols = c(identifier_col, gene_symbolcol)
       mgenesInforA <- allMatrix[,inforCols]
       allinfo = rep("all", dim(allMatrix)[1]) # make a false module with all genes
       allinfo=factor(allinfo,levels=allinfo)
       mgenesInfor <- cbind(mgenesInforA, allinfo)
       colnames(mgenesInfor) <- c( colnames(mgenesInforA), "module")
    }

    #print(mgenesInfor[1:5,])

    #unify upper/lower cases
    mgenesInfor[, 1] = toupper( as.character(mgenesInfor[, 1]))

    rowTitles=names(allMatrix)
    dim(mgenesInfor)

    color1 = mgenesInfor$module
    ctable = table(color1)
    ctable
    
    if(OntologyType=="QTL"){
        top2colnames=c("Chromosome","LocusBin")
    }else if (OntologyType=="TF"){
        top2colnames=c("TranscriptionFactor","TF")
    } else if (OntologyType=="knockout"){
        top2colnames=c("knockout","knockout")
    } else{
        top2colnames=c("System","Gene Category")
    }
    
    #eventualColnames=c(top2colnames, "List Hits", "List Total", 
    #                   "Population Hits", "Population Total",
    eventualColnames=c(top2colnames, "ModuleOverlap", "GeneSetCount", 
                       "PopulationOverlap", "Population",
    #                   "Pvalue_FisherExactTest","Gene Symbol","Unique ID")
                       "pvalue_corrected", "Fisher_pvalue", "Probe Set", "Gene Symbol") 


    #delete exiting reuslts
    firstfile = NULL
    outfiles  =NULL
    for (each in names(ctable) ){
       testoutputfname = paste(mfname, "_",OntologyType, "_",  each, ".xls",  sep='')
       write.table(t(as.matrix(eventualColnames)), testoutputfname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
       if ( is.null(firstfile) ){
           firstfile = testoutputfname 
       }
       outfiles = c(outfiles, testoutputfname)
    }

    for (each in ontologyfnlist){
       cat(each, "\n")
       if (useEASEannotation){
           EASEOntologyEnrichmentAnalysis(genesInfor=mgenesInfor, ontologyfn=each, 
                            fname=mfname, maxSignifLevel=signifLevel, myMinPopulationHits=minPopHits,
                            OntologyType=OntologyType, background=background)
       }else{
           OntologyEnrichmentAnalysis(genesInfor=mgenesInfor, ontologyfn=each, 
                                      fname=mfname, maxSignifLevel=signifLevel, myMinPopulationHits=minPopHits,
                                      OntologyType=OntologyType, background=background)
       }
    }

    # combine all
    topMatrix=NULL
    coutputfname = paste(mfname, "_",OntologyType, ".xls",  sep='')
    for(j in c(1:length(outfiles) ) ) {
       efo = outfiles[j]
       print (efo)
       ematrix= retrieveTopNrows(efo, topNrows=ctopNrows,  OntologyType=OntologyType)
       if (!is.null(ematrix) ){
          topMatrix=rbind(topMatrix, ematrix)
       }
      
    }
    #write.table(topMatrix[,-c(2)], coutputfname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

    if(removeIndividualFiles){file.remove(outfiles)}

    # use shortnames
    systems= as.character(topMatrix[,3])
    topMatrix[,3]=ifelse(systems=="GO Biological Process", "BP",
                     ifelse(systems=="GO Cellular Component","CC", 
                     ifelse(systems=="GO Molecular Function","MF", systems) ))

    msizes = as.integer(topMatrix[,2])
    szorder= order(-msizes)

    write.table(topMatrix[szorder,], coutputfname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)

    # return the first Ontology output file
    return (firstfile )
}

# get LLIds of the last column
#* GOMolecularFunction   ubiquitin C-terminal hydrolase activity	39	        11229             23123; 10234; 33342; 
unlistLLI = function(easeline){
   llids=unlist( strsplit( as.character(easeline[5]), "; ") )
   llids
}


###########################################################################################################
#
#  for small gene set: in output file, genes as column and markers as row
#
###########################################################################################################

computeQTLmatrix_4smallgeneset = function(inputfname, headCol=8, genoCSVfname, outdir="")
{
    #----------------------------- 1. READ expression -------------------------------------
    allMatrix <- read.delim(inputfname,sep="\t", header=T)
    attach(allMatrix)
    dim(allMatrix)

    qtlFileDir ="./"

    no.cols <- dim(allMatrix)[2]
    rowTitles <- rownames(allMatrix)
    colTitles <- colnames(allMatrix)

    fname      = getFileName(inputfname)
    eqtlFname  = paste(outdir,fname, "_eQTLmatrix.xls", sep='')
    countFname  = paste(outdir,fname, "_eQTLcounter.xls", sep='')

    #These are the expression values, the last three cols are module infomation
    datExpr <- t(allMatrix[, -c(1:headCol)])
    dim(datExpr)

    genesInfor  <- allMatrix[,c(1:headCol)]
    genenames   <- allMatrix[,1]

    #----------------------------- 2. READ geno markers -------------------------------------
    #datGenoTable  <- read.cross("csv", dir=qtlFileDir, file=genoCSVfname)
    #datGeno        = calc.genoprob(datGenoTable, step=0, error=0.01)

    genoMatrix <- read.delim(genoCSVfname,sep=",", header=F)
    dim(genoMatrix)
    markerInfor <- t(genoMatrix[c(1:3),-c(1,2)])
    dim(markerInfor)

    colnames(markerInfor) <- c("SNP","chr","pos")

    #----------------------------- 3. Compute eQTL matrix -----------------------------------
    #
    # eQTL for all genes, rows as SNPs and column as genes
    eQTLmatrix = apply(datExpr, 2, qtlOnGenome) 

    colnames( eQTLmatrix ) <- as.character(genesInfor[,1])
    eQTLmatrix.4save = cbind(markerInfor, eQTLmatrix)

    #here we threshold eQTL matrix to get the total number of genes with eQTL at each marker
    thresholds = seq(2, 6, 0.5)
    countermatrix = NULL
    for (ithresh in thresholds){
	bool.eQTLmatrix = eQTLmatrix >= ithresh
        icount=apply(bool.eQTLmatrix, 1, sum, na.rm =T)
        countermatrix = cbind(countermatrix, icount)
    }
    colnames(countermatrix) <- as.character(thresholds)
    countmatrix.4save = cbind(markerInfor, countermatrix)
    write.table(countmatrix.4save,countFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE) 

    write.table(eQTLmatrix.4save, eqtlFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE) 
    dim(eQTLmatrix)

    rm(eQTLmatrix)
    rm(allMatrix)
}









































###########################################################################################################
#
#  for huge gene set: markers as column and genes as row
#
###########################################################################################################

computeQTLmatrix_4hugegeneset = function(inputfname, headCol=8, genoCSVfname, outdir="")
{
    #----------------------------- 1. READ expression -------------------------------------
    allMatrix <- read.delim(inputfname,sep="\t", header=T)
    attach(allMatrix)
    dim(allMatrix)

    qtlFileDir ="./"

    no.cols <- dim(allMatrix)[2]
    rowTitles <- rownames(allMatrix)
    colTitles <- colnames(allMatrix)

    fname      = getFileName(inputfname)
    eqtlFname  = paste(outdir,fname, "_eQTLmatrix.xls", sep='')
    countFname  = paste(outdir,fname, "_eQTLcounter.xls", sep='')

    #These are the expression values, the last three cols are module infomation
    datExpr <- t(allMatrix[, -c(1:headCol)])
    dim(datExpr)

    no.genes = dim(datExpr)[2]

    genesInfor  <- allMatrix[,c(1:headCol)]
    genenames   <- allMatrix[,1]

    #----------------------------- 2. READ geno markers -------------------------------------
    #datGenoTable  <- read.cross("csv", dir=qtlFileDir, file=genoCSVfname)
    #datGeno        = calc.genoprob(datGenoTable, step=0, error=0.01)

    genoMatrix <- read.delim(genoCSVfname,sep=",", header=F)
    dim(genoMatrix)
    markerInfor <- t(genoMatrix[c(1:3),-c(1,2)])
    dim(markerInfor)

    no.markers = dim(markerInfor)[1]
    colnames(markerInfor) <- c("SNP","chr","pos")

    #----------------------------- 3. Compute eQTL matrix -----------------------------------
    #
    # eQTL for all genes, rows as SNPs and column as genes

    #here we threshold eQTL matrix to get the total number of genes with eQTL at each marker
    thresholds    = seq(2, 6, 0.5)
    no.thresh     = length(thresholds)
    countermatrix = matrix(0, no.markers, no.thresh)
    colnames(countermatrix) <- as.character(thresholds)
    
    # decompose the genes into several small datasets
    intvls    = seq(1, no.genes, 50)
    intvls    = c(intvls, no.genes)
    no.intvls = length(intvls)
    for (i in c(1:(no.intvls-1)) ){
       if (i < no.intvls-1){
          indices        = c(intvls[i]:(intvls[i+1]-1) )
       }else{
          indices        = c(intvls[i]:(intvls[i+1]) )
       }

       # returned matrix with genes as cols and markers as rows
       eQTLmatrix     = apply(datExpr[, indices], 2, qtlOnGenome)
       
       if (i==1){
          #row names for the first subset
          combeQTLmatrix = t(cbind(markerInfor, eQTLmatrix) ) # marker as the top rows
          
          #row names for the first subset
          irownames = c(colnames(markerInfor), as.character(genesInfor[indices,1]) )
          eQTLmatrix.4save = cbind(irownames, combeQTLmatrix)
          write.table(eQTLmatrix.4save, eqtlFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE) 

       }else {
          irownames = c(as.character(genesInfor[indices,1]) )
          eQTLmatrix.4save = cbind(irownames, t(eQTLmatrix))
          write.table(eQTLmatrix.4save, eqtlFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T) 
       }


       for ( k in c(1:no.thresh) ){
          kthresh = thresholds[k]
    	  bool.eQTLmatrix = eQTLmatrix >= kthresh
          icount=apply(bool.eQTLmatrix, 1, sum, na.rm =T)
          countermatrix[,k] = countermatrix[,k] + icount
       }

    }

    countmatrix.4save = cbind(markerInfor, countermatrix)
    write.table(countmatrix.4save,countFname, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)
}

############################### NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN #######################################

# input: "literature-network_unique_tommodules_Ontology_blue.xls"
# output: "blue"
getModuleNameFromEASEresfile=function(easeRESfile, OntologyPatt="_Ontology_"){
  a=unlist( strsplit(easeRESfile, OntologyPatt) )
  if (length(a) <=1)
       return ("")
  b=unlist( strsplit(a[2], ".xls") )
    if (length(b) <=0)
       return ("")
  return (b)
}

retrieveTopNrows=function(inputfname, topNrows=3, OntologyType="Ontology"){

  ontopatt = paste("_", OntologyType, "_",sep="")
  
  modulename = getModuleNameFromEASEresfile(inputfname, OntologyPatt=ontopatt)
  if (modulename=="")
     return (NULL)

  allMatrix <- read.delim(inputfname,sep="\t", header=T)
  dim(allMatrix)
  nrows = dim(allMatrix)[1]
  if (nrows<=0)
     return (NULL)
  actualrows = min(topNrows, nrows)

  extracolumn = rep(modulename, actualrows)
  
  outmatrix   = cbind(extracolumn, allMatrix[c(1:actualrows),])
  newcolnames = c("module", colnames(allMatrix) )

  colnames(outmatrix)  <- newcolnames

  return (outmatrix)
}


writeHugeTable = function(mymatrix, outfname, colnames=F, myappend=F, step=3000, use_signif=NA)
{

  if (colnames){
     write.table(rbind(colnames(mymatrix)), outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=myappend)
  } else{
    if (!myappend) {
      yfile <- file(outfname, "w+")
      close(yfile)
    }
  }


  # write into file in many times for saving memory
  norows=dim(mymatrix)[1]
  
  intervals = seq(from=1, to=norows, by=step )

  # test if the last one equal the no.rows, otherwise, append no.rows to the end
  no.intrervals= length(intervals)
  if (intervals[no.intrervals] != norows){
    intervals = c(intervals, norows)
  }
  no.intrervals= length(intervals)
  intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

  total_matched=0
  for (i in c(2:no.intrervals) ){
   sidx= intervals[i-1]
   eidx= intervals[i]-1
   #cat(sidx,"\n")
   #cat(eidx,"\n")
   irows= c(sidx:eidx)
   imatrix=mymatrix[irows, ]

   if(is.na(use_signif)) {
       write.table(imatrix, outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
   }else{
       write.table(signif(imatrix,use_signif), outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
   }

   #write.table(imatrix, outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
  }

}

# restrict certain cols and rows for the output
# if we need save the transposed matrix, we need save the rows as whole, 
#   therefore, the intervals are based on the columns
#
writeHugeTableSELs = function(mymatrix, outfname, rowSel, colSel, colnames=F, myappend=F, step=3000, transpose=T)
{

  if (colnames){
    if(transpose){
       write.table(t(as.matrix( rownames(mymatrix) )), outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=myappend)
    }else{
       write.table(t(as.matrix( colnames(mymatrix) )), outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=myappend)
    }  
  } else{
    yfile <- file(outfname, "w+")
    close(yfile)
  }

  # write into file in many times for saving memory
  if(transpose) {
    norows=dim(mymatrix)[2]
    rowSelIdx = c(1:norows)[colSel]
  }else{
    norows=dim(mymatrix)[1]
    rowSelIdx = c(1:norows)[rowSel]
  }

  intervals = seq(from=1, to=norows, by=step )

  # test if the last one equal the no.rows, otherwise, append no.rows to the end
  no.intrervals= length(intervals)
  if (intervals[no.intrervals] != norows){
    intervals = c(intervals, norows)
  }
  no.intrervals= length(intervals)
  intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

  total_matched=0
  for (i in c(2:no.intrervals) ){
   sidx= intervals[i-1]
   eidx= intervals[i]-1
   #cat(sidx,"\n")
   #cat(eidx,"\n")
   irows= c(sidx:eidx)
   irows= intersect(rowSelIdx, irows)   

   if(transpose){
      imatrix=t(mymatrix[rowSel, irows])
   }else{
      imatrix=mymatrix[irows, colSel]
   }

   write.table(imatrix, outfname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

  }

}




#################### Histogram -based Functions ###########################################


# Here we compute the histogram of a list of correlation coefficients [0,1]
# but we customize the output as the frequency in fixed bins specified by
#  the input parameter "nobreaks", eg, 20 intervals in the range of [0,1]
#
histogram_customize = function(corvect,nobreaks=20,xmin=NULL,xmax=NULL,vectorname="", 
             xlabel="", imgfilename=NULL, showplot=T)
{
  # draw histogram, compute frequency of no of links
  #numvect = as.numeric(corvect)
  maxsize = 1000000

  nfold = nobreaks/10
  if(length(corvect) <= maxsize ){
      numvect = as.numeric(corvect)+0.05/nfold
  }else{
      tmpvect= sample(corvect, maxsize, replace=F)
      numvect = as.numeric(tmpvect)+0.05/nfold
  }

  ntab = table(round(numvect*nfold, 1) )
  linkfreq = as.numeric(ntab)
  linkidx  = as.numeric( names(ntab) )/nfold
  linkidx

  no.items = length(numvect)

  histMatrix = cbind(as.character(linkidx), as.character(linkfreq))

  if(is.null(xmin)){
     mxmin=min(numvect,na.rm=T)
  }else{
     mxmin=xmin
  }
  if( is.null(xmax) ){
     mxmax=max(numvect,na.rm=T)
  }else{
     mxmax=xmax
  }

  completeIntervals  = seq(from=mxmin, to=mxmax, by=1/nobreaks)


  completeHistMatrix = merge(as.matrix(as.character(completeIntervals)), as.matrix(histMatrix), by.x=1, by.y=1, all.x=T)
  completeHistMatrix <- as.matrix(completeHistMatrix)
  #completeHistMatrix
  completeHistMatrix = ifelse(is.na(completeHistMatrix),0, completeHistMatrix)
  #completeHistMatrix
  completeHistMatrix2 <- matrix(as.numeric(completeHistMatrix), ncol=2)
  #completeHistMatrix2

  completeHistMatrix2[,2] = completeHistMatrix2[,2] /no.items

  #mylabel = vectorname
  #outimg = paste(filename, "-HIST", vectorname, ".png", sep="")

     index=completeHistMatrix2[,1]
     freq =completeHistMatrix2[,2]

     maxval = round(max(corvect,na.rm=T),2)
     meanval= round(mean(corvect,na.rm=T), 2)

     mylabel = paste(vectorname, ": max=", as.character(maxval), 
                           ", mean=", as.character(meanval), sep="" )

  if (!is.null(imgfilename) ){
     openImgDev(imgfilename, iwidth = 600, iheight = 600)
  }

  if (showplot) {
     barplot(freq, names.arg= as.character(index),
       xlab=xlabel, ylab="frequency", main= mylabel, col="blue")
  }

    #plot(c(0,1), c(0,1), xlab=xlabel,ylab="frequency",col='white' )
    #lines(binned.x, linkfreq, col="black");

    #plot(linkidx, linkfreq,xlab="number of links of a node", ylab="frequency", main= keyword, type="h")
    #histogram(as.numeric(corvect), br=seq(from=0,to=1, by=0.05), col="blue",
    #         xlab=xlabel, ylab="frequency", main= mylabel)

  if (!is.null(imgfilename) ){
    dev.off()
  }

  return (completeHistMatrix2)
}

# bar plot of count and percent of integer vector
#
# if the no of unique values >100, the program automatically compute the binned distribution
#  otherwise, no bins are used for distribution computation
#
histogram_4integers =function(ivector,fkeyname=NULL,keyword="",hxlabel="k",no.cuts= 15, imgwid=600, imghei=600){

    maxsize = max(ivector)
    minsize = min(ivector)
    if (length(ivector) >1) {
       szrange = c(minsize : maxsize)
    }else{
       szrange = c((ivector-1):(ivector+1) )
    }

    uniques = names(table(ivector))
    if (length(uniques) < 100) {
        no.cuts = length(szrange)
        binned.k= szrange
    } else{
        intervalSize = (maxsize -minsize)/no.cuts
        binned.k     = rep(minsize-intervalSize/2, no.cuts)
        for (kk in c(1:no.cuts) ) {
           binned.k[kk] = as.integer( (binned.k[kk] + kk*intervalSize) )
        }
    }
    
    cut1 = cut(ivector, no.cuts )
    frequence= tapply(ivector,cut1,length)
    frequence= as.numeric(frequence)
    frequence= ifelse( is.na(frequence),0, frequence )
    percent  = frequence/length(ivector)

    imgCount  =paste( fkeyname, "-Count.png",sep='')  #png only
    imgPercent=paste( fkeyname, "-Percent.png",sep='')  #png only

    # save statistics into log file
    #
    counts = c("count",   length(ivector), frequence)
    freqs  = c("percent(%)", 100,         signif(percent,3)*100)
    title  = c(keyword, "total", as.character(binned.k) )
    logMatrix = rbind(title, counts, freqs)

    if(is.null(fkeyname)){
        return (logMatrix)
    }

    maintitle =""
    if (keyword !=""){
        #maintitle = paste("distribution of number of ", keyword, sep="") 
        maintitle = paste("distribution of ", keyword, sep="") 
    }

    openImgDev(imgCount,  iwidth = imgwid, iheight = imghei)
    barplot(frequence, names.arg= as.character(binned.k),  col="green",
     ylab="count", xlab="k", main= maintitle)
    dev.off()

    openImgDev(imgPercent,  iwidth = 600, iheight = 600)
    #histogram(kcliquesSize, type="percent",xlab=hxlabel)
    barplot(percent*100, names.arg= as.character(binned.k), col="green",
     ylab="percent", xlab=hxlabel, main= maintitle )
    dev.off()

    return(logMatrix)
}



#
# wil.cox test is very slow for huge dataset, so here we use sampling method to solve this dilemma
# if a vector is too big (>maxsize), then we need sample the dataset into a set with length
#   < maxsize
#
wilcoxtest4Hugedata = function(numericA, numericB, maxsize = 50000)
{
   #maxsize = 50000
   if ( length(numericA)>maxsize ){
       snumericA = sample(numericA, maxsize, replace=F)
   }else{
       snumericA = numericA
   }
   if ( length(numericB)>maxsize ){
        snumericB = sample(numericB, maxsize, replace=F)
   }else{
        snumericB = numericB
   }

   wtest = wilcox.test(snumericA, snumericB)
   wtest$p.value
}


# perform wil.cox ranking test on the column based vectors
# return a table of pvalue table
wilcoxtest4sets = function(datalist, usewilcox=T)
{
  no.vects = length(datalist)
  datout    = matrix('-',nrow=no.vects,ncol=no.vects)
  sets      = names(datalist)
  for (i in c(1:(no.vects-1) ) ){
   for (j in c((i+1):(no.vects) ) ){
       if(usewilcox){
          ijpval      = wilcoxtest4Hugedata(datalist[[i]], datalist[[j]])
       }else{
          tt    = t.test(datalist[[i]], datalist[[j]])
          ijpval=tt$p.value
       }
       datout[i,j] = signif(ijpval, 3)
   }
  }

  # add a column for the set names
  fdatoutA = cbind(sets, datout)
  #print(fdatout)
  #print(sets)
  colnames(fdatoutA) <- c("Wilcoxon test", sets)

  fdatout = rbind(c("Wilcoxon.pvalue", sets), fdatoutA)

  fdatout
}


# transform single point based histograms into step-based representation
#  
transformHisto2Steps=function(mindex, mhist){
  npoints = length(mindex)
  mstep   = mindex[2] - mindex[1]
  hstep   = mstep/2
  newindex=NULL
  newhist =NULL
  for (i in c(1:npoints) ){
       idx = mindex[i] # current value
       iy  = mhist[i]
       newindex = c(newindex, idx-hstep)     
       newindex = c(newindex, idx+hstep)
       newhist  = c(newhist, iy)
       newhist  = c(newhist, iy)
  }
  newpoints = length(newindex)

  #process start and end which are now out of the original limit
  # because of operations + and - hstep
  newindex[1]=mindex[1]
  newindex[newpoints]=mindex[npoints]

  retMatrix = cbind(newindex, newhist)  

  return (retMatrix )
}

# plot multi-histograms in one figure based on the step-style plot
# mindex:     the mean values of intervals
# histograms: the histograms, all were computed in the same intervals
#
# dispXp: intervals for displaying string
#
plotMultiHistsograms= function(histograms,mxindex,xlabel="",ylabel="",maintitle="",filename, 
                               phei=400, pwid=500, ppointsize=12, punits="px", pres = 72, pcolors=NULL,
                               stepshape=T, mycex=1, mycexText=1, show_legend=T, horizontal_line=NULL,
                               display_string=NULL, dispYp=0, dispXp=4,showaxis=TRUE, ylim=NULL )
{
 no.hists = dim(histograms)[2]

 if(stepshape){
  stepHists = NULL
 
  for (i in c(1:no.hists) ){
    istephist=transformHisto2Steps(mxindex, histograms[,i])
    stepHists = cbind(stepHists, istephist[,2])
  }
  stepindex = istephist[,1]
  stepHists <- as.matrix(stepHists)
  colnames(stepHists) <- colnames(histograms)

  plotMultiProfiles(eQTL_profiles=stepHists, xindex=stepindex, 
            xlabel=xlabel,ylabel=ylabel,mainlabel=maintitle, 
            plottype="l", ltype=rep(1,no.hists), filename=filename, mycex=mycex, showaxis=showaxis,
            show_legend=show_legend, horizontal_line=horizontal_line, display_string=display_string, 
            dispYp=dispYp, dispXp=dispXp, ylim=ylim,
            myhei = phei, mywid =pwid, mypointsize=ppointsize, myunits=punits, myres =pres, mycolors=pcolors)
 }else{
  plotMultiProfiles(eQTL_profiles=histograms, xindex=mxindex,
            xlabel=xlabel,ylabel=ylabel, mainlabel=maintitle, 
            plottype="l", ltype=rep(1,no.hists), filename=filename, mycex=mycex, showaxis=showaxis,
            show_legend=show_legend, horizontal_line=horizontal_line, display_string=display_string, 
            dispYp=dispYp, dispXp=dispXp, ylim=ylim,
            myhei = phei, mywid =pwid, mypointsize=ppointsize, myunits=punits, myres =pres, mycolors=pcolors)
 }
 
}

# inputs:
# networkfname ~ network link file where each entry is a pair of genes, separated a separator
# idxMatrix    ~ a matrix where the 1st column includes the gene names and the 2nd column are the 
#                coresponding index in the correlation matrix
#
# The output is a 1x2 matrix, where each row, containing two integer numbers, 
#         corersponds to a link in the network file
#
matchNetworkLinks2Corrmatrixindex = function(networkfname, idxMatrix)
{
  bbMatrixAll <- read.delim(networkfname, sep="\t", header=F)
  bbMatrixAll <- as.matrix(bbMatrixAll)
  #dim(bbMatrixAll)

  bbMatrix <- bbMatrixAll[,c(1,2)]
  indexMatrixInt = matchNetworkLinks2Corrmatrixindex_Subfunc(netlinks=bbMatrix, 
                 idxMatrix=idxMatrix)

  return (indexMatrixInt)
}

matchNetworkLinks2Corrmatrixindex_Subfunc = function(netlinks, idxMatrix) {
  leftMatched = merge(idxMatrix, netlinks, by.x=1, by.y=1,all=F)
  #dim(leftMatched)

  bothMatched = merge(idxMatrix, leftMatched, by.x=1, by.y=3,all=F)
  #dim(bothMatched)

  indexMatrix  = bothMatched[,c(2,4)]
  indexMatrix  <- as.matrix(indexMatrix)

  indexMatrixS = cbind(as.matrix(bothMatched[,4]), as.matrix(bothMatched[,2]) )
  indexMatrixT = data.frame(indexMatrixS)

  indexMatrixInt = cbind(as.integer(indexMatrix[,1]), as.integer(indexMatrix[,2]))

  return (indexMatrixInt)
}


# turn one row of an adj matrix into pairwise representation
#
adj2pair = function(adjrow, ithresh=0.8, genenames, genenamesRow, outfile=NULL)
{  
   if(is.null(genenamesRow)){
      genenames2= genenames
   } else {
      genenames2= genenamesRow
   }

   rlen = length(adjrow)
   sidx = adjrow[rlen] #source node
   if(sidx >=(rlen-1) ){
       return (NULL)
   }   

   #consider only the upper part of the matrix above diagonal
   mask = c(rep(F, sidx), rep(T, rlen-sidx))

   sel = adjrow >= ithresh # destination nodes
   sel = (sel & mask)

   sel[rlen] = F # index should be reset   

   if(sum(sel)>0){
     selgenenames <- genenames[sel]
     pairnet=paste(genenames2[sidx], selgenenames, sep="\t")

     if(!is.null(outfile) ) {
        write.table(as.matrix(pairnet), outfile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
     }

     return ( cbind(genenames2[sidx], selgenenames) )
   }
   return (NULL) 
}


adj2pairSimple = function(adjrow, sidx, genenames, genenamesRow=NULL, outfile=NULL, directed =F, ouputweight=F, symmetric=TRUE)
{
   if(is.null(genenamesRow)){
      genenames2= genenames
   } else {
      genenames2= genenamesRow
   }
  
   if(ouputweight){
      sel  = rep(T, length(adjrow)) # destination nodes
   } else{
      sel  = adjrow ==1 # destination nodes
   }

   # for undirected network, we consider only the upper diaganal part
   #  the trick is to mask out the lower part of the adjrow
   #
   if( (!directed) & symmetric) {
     maskoutIdx = c(1:sidx)
     sel[maskoutIdx ] = FALSE
   }

   if(sum(sel)>0){
     selgenenames <- genenames[sel]
     if (ouputweight) {
         pairnet=paste(genenames2[sidx], selgenenames, adjrow[sel], sep="\t")
     } else{
         pairnet=paste(genenames2[sidx], selgenenames, sep="\t")
     }
   
     if(!is.null(outfile)) {
        write.table(as.matrix(pairnet), outfile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
     }

     return ( pairnet )
   }
   return (NULL) 
}


adjmatrix2linkpairs = function(adjmatrix, genenames, genenamesRow=NULL, outfile=NULL, 
                               directed = FALSE, ouputweight=FALSE, symmetric=TRUE)
{
   #clean up
   if(!is.null(outfile)) {
     yfile <- file(outfile, "w+")
     close(yfile)
   }

   nrows = dim(adjmatrix)[1]
   netMatrix = NULL
   for(i in c(1:nrows) ){
       ires = adj2pairSimple(adjrow=adjmatrix[i,], sidx=i, genenames=genenames, 
           genenamesRow=genenamesRow, outfile=outfile, directed=directed, ouputweight=ouputweight,
           symmetric=symmetric)
       netMatrix = c(netMatrix, ires)
   }

   if(!is.null(outfile)) {
     netMatrix <- read.delim(outfile, sep="\t", header=F)
     netMatrix  = as.matrix(netMatrix)
   }

   return(netMatrix)
}


# convert multiple matrices, with the first column as IDs
#
general_matrices_to_pairs = function (ymtrxNumb, ymtrxPval, ymtrxFreq) {

  otitleNumb <- colnames(ymtrxNumb)
  otitlePval <- colnames(ymtrxPval)
  otitleFreq <- colnames(ymtrxFreq)

  colnames(ymtrxNumb) <- paste("Number",      otitleNumb)
  colnames(ymtrxFreq) <- paste("Fold_Change", otitleFreq)
  colnames(ymtrxPval) <- paste("Pvalue",      otitlePval)

  ymtrxNumb[1:3,]
  ymtrxFreq[1:3,]
  ymtrxPval[1:3,]

  dim(ymtrxNumb)
  dim(ymtrxFreq)
  dim(ymtrxPval)

  if(FALSE) {
  # merge
  final = merge(ymtrxNumb, ymtrxFreq, by.x=1, by.y=1, all=T)
  final = merge(final,    ymtrxPval, by.x=1, by.y=1, all=T)
  final = as.matrix(final)
  colnames(final) <- c("Module", colnames(final)[-1])

  # matrix form
  #write.table(final, foutymtrx, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
  }

  # pairwise form
  yymtrxNumb = as.matrix(ymtrxNumb)[,-1]
  yymtrxFreq = as.matrix(ymtrxFreq)[,-1]
  yymtrxPval = as.matrix(ymtrxPval)[,-1]
  dim(yymtrxNumb)
  dim(yymtrxFreq)
  dim(yymtrxPval)

  if(is.null(dim(yymtrxNumb)) ) {
     index     = c(1: length(yymtrxNumb) )
  } else {
     index     = getMatrixIndex(dim(yymtrxNumb), symmetric=F, diagonal=T)
  }

  yrownames = as.character(as.matrix(ymtrxFreq)[,1])
  ycolnames = getSecondPart(colnames(ymtrxFreq)[-1], " ", 2)

  datMatrix = cbind(yymtrxNumb[index], yymtrxPval[index], yymtrxFreq[index])

  if(is.null(dim(yymtrxNumb)) ) {
    final2    = cbind(yrownames[index], rep(ycolnames, length(index)), datMatrix)     
  } else {
    final2    = cbind(yrownames[index[,1]], ycolnames[index[,2]], datMatrix)
  }
 
  colnames(final2) <- c("ID", "signature", "Number","pvalue", "fold_change")

  return (final2)
}



# convert multiple matrices
#
gomatrices_to_pairs = function (ymtrxNumb, ymtrxPval, ymtrxFreq) {

  otitleNumb <- colnames(ymtrxNumb)
  otitlePval <- colnames(ymtrxPval)
  otitleFreq <- colnames(ymtrxFreq)

  colnames(ymtrxNumb) <- paste("Number",      otitleNumb)
  colnames(ymtrxFreq) <- paste("Fold_Change", otitleFreq)
  colnames(ymtrxPval) <- paste("Pvalue",      otitlePval)

  ymtrxNumb[1:3,]
  ymtrxFreq[1:3,]
  ymtrxPval[1:3,]

  dim(ymtrxNumb)
  dim(ymtrxFreq)
  dim(ymtrxPval)

  # merge
  final = merge(ymtrxNumb, ymtrxFreq, by.x=1, by.y=1, all=T)
  final = merge(final,    ymtrxPval, by.x=1, by.y=1, all=T)
  final = as.matrix(final)
  colnames(final) <- c("Module", colnames(final)[-1])

  # matrix form
  #write.table(final, foutymtrx, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

  # pairwise form
  if(nrow(ymtrxPval)>2 & ncol(ymtrxPval)>2) {
    yymtrxNumb = as.matrix(ymtrxNumb)[-1,-c(1,2)]
    yymtrxFreq = as.matrix(ymtrxFreq)[-1,-1]
    yymtrxPval = as.matrix(ymtrxPval)[-1,-1]
  } else if(nrow(ymtrxPval)==2) {
    yymtrxNumb = rbind(as.matrix(ymtrxNumb)[-1,-c(1,2)])
    yymtrxFreq = rbind(as.matrix(ymtrxFreq)[-1,-1])
    yymtrxPval = rbind(as.matrix(ymtrxPval)[-1,-1])
  } else {
    yymtrxNumb = cbind(as.matrix(ymtrxNumb)[-1,-c(1,2)])
    yymtrxFreq = cbind(as.matrix(ymtrxFreq)[-1,-1])
    yymtrxPval = cbind(as.matrix(ymtrxPval)[-1,-1])
  }
  dim(yymtrxNumb)
  dim(yymtrxFreq)
  dim(yymtrxPval)

  modulesizes   = ymtrxNumb[-1,2]
  singaturesize = ymtrxNumb[1,-c(1,2)]

  if(is.null(dim(yymtrxNumb)) ) {
     index     = c(1: length(yymtrxNumb) )
  } else {
     index     = getMatrixIndex(dim(yymtrxNumb), symmetric=F, diagonal=T)
  }

  yrownames = as.character(as.matrix(ymtrxFreq)[-1,1])
  ycolnames = getSecondPart(colnames(ymtrxFreq)[-1], " ", 2)

  datMatrix = cbind(yymtrxNumb[index], yymtrxFreq[index], yymtrxPval[index])

  if(is.null(dim(yymtrxNumb)) ) {
    final2    = cbind(yrownames[index], modulesizes[index], 
                    rep(ycolnames, length(index)), rep(singaturesize,length(index)), datMatrix)     
  } else {
    final2    = cbind(yrownames[index[,1]], modulesizes[index[,1]], 
                    ycolnames[index[,2]], singaturesize[index[,2]], datMatrix)
  }
 
  sel = !is.na(datMatrix[,3])  
  final2 = final2[sel, ] 
  colnames(final2) <- c("module", "modulesize", "signature", "signature_size", "overlap","fold_change", "pvalue" )

  return (final2)
}

# merge three GO files
#
merge_GO_PvalFreqNumb = function(prefixes, indir="", fname="", specialHBTRC=F)
{

nocases = length(prefixes)

flag = c("restricted network B modules to network A", "restricted network B to network A")

# output
#foutMtrx = paste(indir, fname, "_matrix.xls", sep="")
foutPair = paste(indir, fname, "_pair.xls",   sep="")
foutPair2= paste(indir, fname, "_pair-restricted.xls",   sep="")

xfinal = NULL
xfinal2= NULL

for (i in c(1:nocases) ) {

  fmemb = paste(indir, prefixes[i], ".xls",           sep="")
  fpval = paste(indir, prefixes[i], "_pvalue.xls",    sep="")
  ffreq = paste(indir, prefixes[i], "_frequence.xls", sep="")
  fnumb = paste(indir, prefixes[i], "_number.xls",    sep="")

  print( paste(i, prefixes[i]) )

  mtrxPval = getSecondTable(fpval, firstTable=T, endflag=flag, 
               blankline_before=1, blankline_after=1)

  mtrxFreq = getSecondTable(ffreq, firstTable=T, endflag=flag, 
               blankline_before=1, blankline_after=1)

  mtrxNumb <- read.delim(fnumb,sep="\t", header=T)
  mtrxNumbR<- as.matrix(mtrxNumb)[-1,] # restricted
  mtrxNumb <- as.matrix(mtrxNumb)[-2,]
  mtrxNumb[1,1]="ALL"
  mtrxNumbR[1,1]="ALL"

  mtrxFreq <- as.matrix(mtrxFreq)
  mtrxPval <- as.matrix(mtrxPval)

  mtrxMemb <- read.delim(fmemb,sep="\t", header=T)
  mtrxMemb <- as.matrix(mtrxMemb)
  mtrxMemb[1:3,]

  ncc= dim(mtrxMemb)[2]
  membBYmod = tapply(X=mtrxMemb[,1], INDEX=mtrxMemb[,ncc], FUN=concatenate, mysep=",")

  #[1,] "CR"  "A"  "N"  "all" 
  #[2,] "CR"  "A"  "N"  "down"
  #
  ncols = dim(mtrxNumb)[2]

  if(specialHBTRC) {
    parts = getAllParts( colnames(mtrxNumb)[-c(1:2)], "\\.")
    sel=parts[,4]!="all"; 
    selIdx1   = c(1,2, c(3:ncols)[sel])
    mtrxNumb  = mtrxNumb[, selIdx1]
    mtrxNumbR = mtrxNumbR[, selIdx1]
    ncols = dim(mtrxPval)[2]
    selIdx2   = c(1, c(2:ncols)[sel])
    mtrxPval  = mtrxPval[,selIdx2]
    mtrxFreq  = mtrxFreq[,selIdx2]
  }

  #ymtrxNumb=mtrxNumb; ymtrxPval=mtrxPval; ymtrxFreq=mtrxFreq 
  tmp=gomatrices_to_pairs(ymtrxNumb=mtrxNumb, ymtrxPval=mtrxPval, ymtrxFreq=mtrxFreq )
  rnTmp = paste(tmp[,1], tmp[,3], sep=".")
  tmp2  = cbind(tmp, membBYmod[rnTmp])
  colnames(tmp2) <- c(colnames(tmp), "members")
  xfinal = rbind(xfinal, tmp2)
  

  # 2nd table
  #
  mtrxPval2 = getSecondTable(fpval, firstTable=FALSE, endflag=flag, 
               blankline_before=1, blankline_after=1)

  mtrxFreq2 = getSecondTable(ffreq, firstTable=FALSE, endflag=flag, 
               blankline_before=1, blankline_after=1)
  mtrxFreq2 <- as.matrix(mtrxFreq2)
  mtrxPval2 <- as.matrix(mtrxPval2)

  if(specialHBTRC) {
     mtrxPval2  = mtrxPval2[,selIdx2]
     mtrxFreq2  = mtrxFreq2[,selIdx2]
  }

  nsize=dim(mtrxPval2)
  mtrxPval3  = rbind(rep(1,nsize[2]), mtrxPval2); colnames(mtrxPval3) <- colnames(mtrxPval2)
  tmp2  = gomatrices_to_pairs(ymtrxNumb=mtrxNumbR, ymtrxPval=mtrxPval3, ymtrxFreq=mtrxFreq2)

  rnTmp = paste(tmp2[,1], tmp2[,3], sep=".")
  tmp3  = cbind(tmp2, membBYmod[rnTmp])
  colnames(tmp3) <- c(colnames(tmp2), "members")

  xfinal2 = rbind(xfinal2, tmp3) 

}

ncols = dim(xfinal2)[2]-1
od = order(as.numeric(xfinal2[,ncols]))
xfinal2 = xfinal2[od, ] 

od = order(as.numeric(xfinal[,ncols]))
xfinal = xfinal[od, ] 

# matrix form
write.table(xfinal, foutPair, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
write.table(xfinal2, foutPair2, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

}








# merge three GO files
#
merge_GO_PvalFreqNumb_Matrix = function(prefixes, indir="", fname="", specialHBTRC=F, rowchoices=NULL)
{

nocases = length(prefixes)
flag = "restricted network B modules to network A"

for (i in c(1:nocases) ) {

  fmemb = paste(indir, prefixes[i], ".xls",    sep="")
  fpval = paste(indir, prefixes[i], "_pvalue.xls",    sep="")
  ffreq = paste(indir, prefixes[i], "_frequence.xls", sep="")
  fnumb = paste(indir, prefixes[i], "_number.xls",    sep="")

  print( paste(i, prefixes[i]) )

  mtrxPval = getSecondTable(fpval, firstTable=T, endflag=flag, 
               blankline_before=1, blankline_after=1)

  mtrxFreq = getSecondTable(ffreq, firstTable=T, endflag=flag, 
               blankline_before=1, blankline_after=1)

  mtrxNumb <- read.delim(fnumb,sep="\t", header=T)
  mtrxNumbR<- as.matrix(mtrxNumb)[-1,] # restricted
  mtrxNumb <- as.matrix(mtrxNumb)[-2,]
  mtrxNumb[1,1]="ALL"
  mtrxNumbR[1,1]="ALL"

  mtrxFreq <- as.matrix(mtrxFreq)
  mtrxPval <- as.matrix(mtrxPval)

  mtrxMemb <- read.delim(fmemb,sep="\t", header=T)
  mtrxMemb <- as.matrix(mtrxMemb)
  mtrxMemb[1:3,]

  ncc= dim(mtrxMemb)[2]
  membBYmod = tapply(X=mtrxMemb[,1], INDEX=mtrxMemb[,ncc], FUN=concatenate, mysep=",")

  #[1,] "CR"  "A"  "N"  "all" 
  #[2,] "CR"  "A"  "N"  "down"
  #
  ncols = dim(mtrxNumb)[2]

  if(specialHBTRC) {
    parts = getAllParts( colnames(mtrxNumb)[-c(1:2)], "\\.")
    sel=parts[,4]!="all"; 
    selIdx1   = c(1,2, c(3:ncols)[sel])
    mtrxNumb  = mtrxNumb[, selIdx1]
    mtrxNumbR = mtrxNumbR[, selIdx1]
    ncols = dim(mtrxPval)[2]
    selIdx2   = c(1, c(2:ncols)[sel])
    mtrxPval  = mtrxPval[,selIdx2]
    mtrxFreq  = mtrxFreq[,selIdx2]
  }

  # transpose the matrix
  mtrxPval2 = t(mtrxPval); mtrxFreq2 = t(mtrxFreq); mtrxNumb2 = t(mtrxNumb[,-2])
  mtrxPval2[1:3, ]; mtrxFreq2[1:3, ]; mtrxNumb2[1:3, ]

  totalsizesRow = as.integer(mtrxNumb[,2]) 
  if( !is.null(rowchoices) ) {
    idx = getMatchedIndex(mtrxPval[,1], rowchoices)
    mtrxPval2 = mtrxPval2[,idx]; mtrxFreq2 = mtrxFreq2[,idx]; mtrxNumb2 = mtrxNumb2[,idx];
    totalsizesRow = as.integer(mtrxNumb[,2])[idx]
  }
  mtrxNumb2[1,] <- paste(mtrxNumb2[1,], " (", totalsizesRow, ")", sep="") 

  mtrxNumb2[1,] <- paste("Overlap",         mtrxNumb2[1,])
  mtrxFreq2[1,] <- paste("Fold Enrichment", mtrxFreq2[1,])
  mtrxPval2[1,] <- paste("FET p-value",     mtrxPval2[1,])

  nc = ncol(mtrxFreq2);   nr = nrow(mtrxFreq2)
  mtrxFreq3 = NULL
  for(jj in c(1:nc)){
     if(mtrxFreq2[1,jj] == "Fold Enrichment ALL_freq"){
        ivect <- c(mtrxFreq2[1,jj], signif(as.numeric(mtrxFreq2[2:nr,jj]), 3))
     } else{
        ivect <- c(mtrxFreq2[1,jj], round(as.numeric(mtrxFreq2[2:nr,jj]), 3))
     }
     mtrxFreq3 = cbind(mtrxFreq3, ivect)
  }
  #print(mtrxFreq2)

  final= cbind(mtrxNumb2, mtrxFreq3, mtrxPval2)

  final2 = cbind(rownames(mtrxNumb2), final)

  # output
  foutMtrx = paste(indir, prefixes[i], "_3in1-matrix.xls",    sep="")

  # matrix form
  write.table(final2, foutMtrx, sep="\t",quote=FALSE, col.names=FALSE, row.names=FALSE)  
}

}





#--------------------- GO Annotation ----------------------------------
#

if(F) {
gifA=aaMatrix
gifB=bbMatrix
moduleNameA=modlevelsA[i]
moduleNameB=modlevelsB[j]
uniqueIdCol=1
moduleColA=imoduleColA
moduleColB=imoduleColB
totalGenes=totalbackground 
removeDuplicate=removeDuplicate

gifA=aaMatrix
gifB=ccMatrix
moduleNameA=modlevelsA[i]
moduleNameB=modlevelsB[j]
uniqueIdCol=1
moduleColA=imoduleColA
moduleColB=imoduleColB
totalGenes=totalbackground 
removeDuplicate=removeDuplicate

gifA=AAX
gifB=bbMatrix
moduleNameA="ALL"
moduleNameB=modlevelsB[j]
uniqueIdCol=1
moduleColA=2
moduleColB=imoduleColB
totalGenes=totalbackground 
removeDuplicate=removeDuplicate

}

compareTwoModules = function(gifA, gifB, moduleNameA, moduleNameB, uniqueIdCol=1, moduleColA, moduleColB, totalGenes, removeDuplicate=F)
{
  restrictA = gifA[, moduleColA]== moduleNameA 
  restrictB = gifB[, moduleColB]== moduleNameB
  
  restrictA = ifelse(is.na(gifA[,uniqueIdCol]), F, restrictA)
  restrictB = ifelse(is.na(gifB[,uniqueIdCol]), F, restrictB)

  noA= sum(restrictA)
  noB= sum(restrictB)
  moduleSetA = rbind(gifA[restrictA, ])
  moduleSetB = rbind(gifB[restrictB, ])

  if (noA==0 | noB==0){
     ret  = c(0, 1, noA, noB)
     # here we include merged matrix
     finalret= list(ret, NULL)
     return (finalret)
  }
  
  if (removeDuplicate) {
    orderA    = order(moduleSetA[,uniqueIdCol])
    moduleSetA=       rbind(moduleSetA[orderA,])
    boundA    = findBoundary(moduleSetA[,uniqueIdCol])
    moduleSetA= rbind(moduleSetA[boundA,])

    orderB = order(moduleSetB[,uniqueIdCol])
    moduleSetB=    rbind(moduleSetB[orderB,])
    boundB    = findBoundary(moduleSetB[,uniqueIdCol])
    moduleSetB=    rbind(moduleSetB[boundB,])

    noA= dim(moduleSetA)[1]
    noB= dim(moduleSetB)[1]
  }

  mergeMatrix = merge(moduleSetA, moduleSetB, by.x=uniqueIdCol, by.y=uniqueIdCol, sort=F,all=FALSE)
  intersectNo = dim(mergeMatrix)[1]

  # new module assignment
  combmodulename=paste(moduleNameA, ".", moduleNameB, sep="")
  newmergeMatrix = cbind(mergeMatrix, rep(combmodulename, intersectNo) )

  # phyper(89,702,4000-702,280, lower.tail=F)
  if(intersectNo>0) {
    pval = phyper(intersectNo-1, noA, totalGenes-noA, noB, lower.tail=F)

  }else{
    pval = 1
  }

  ret  = c(intersectNo, pval, noA, noB)

  # here we include merged matrix
  finalret= list(ret, newmergeMatrix)
  finalret
}

# restrict_networkB_to_A ==T: take networkB as background when mapping modulesB to network A
# restrict_networkB_to_A ==F: consider modulesB as independent, so map individual modules to network A
#
moduleBasedIntersection=function(fnameA, fnameB,outputDir="", keywords="", uniqueIdCol=c(1,1), geneInforCols=8, itotalGenes = 3600, signifpvalue=0.001, latex=F, removeDuplicate=F, restrict_networkB_to_A=T)
{
    #fkeyB       =getFileName(fnameB)
    #outfname    paste(outputDir, "intersect.csv", sep='')

    # get unique labels for each input files
    keys= unlist(strsplit(keywords, "-vs-"))
    fkeyA       =keys[1]
    fkeyB       =keys[2]

    #outfname     = paste(outputDir,  "intersect_", keywords, ".xls", sep='')
    outfname     = paste(outputDir,  "intersect_", keywords, ".xls", sep='')

    nofname     = paste(outputDir,  "intersect_", keywords, "_number.xls", sep='')
    pvalfname   = paste(outputDir,  "intersect_", keywords, "_pvalue.xls", sep='')
    freqfname   = paste(outputDir,  "intersect_", keywords, "_frequence.xls", sep='')


    latexnofname  = paste(outputDir, "intersect_", keywords, "_number.tex", sep='')
    latexpvalfname= paste(outputDir, "intersect_", keywords, "_pvalue.tex", sep='')

    aaMatrixAll <- read.delim(fnameA, sep="\t", header=T)
    dim(aaMatrixAll)
    aaMatrixAll <- as.matrix(aaMatrixAll)

    bbMatrixAll <- read.delim(fnameB, sep="\t", header=T)
    dim(bbMatrixAll)
    bbMatrixAll <- as.matrix(bbMatrixAll)

    aaMatrixAll[, uniqueIdCol[1]] = toupper(aaMatrixAll[, uniqueIdCol[1]])
    bbMatrixAll[, uniqueIdCol[2]] = toupper(bbMatrixAll[, uniqueIdCol[2]])

    no.colsA = dim(aaMatrixAll)[2]
    no.colsB = dim(bbMatrixAll)[2]

    # include gene infor in the aaMatrix but not bbMatrix, in this way the merged matrix
    #  will not have duplicated columns
    aaMatrix <- aaMatrixAll[, c(uniqueIdCol[1], 1:geneInforCols, no.colsA)]
    colnames(aaMatrix) <- c( colnames(aaMatrixAll)[c(uniqueIdCol[1], 1:geneInforCols)], fkeyA)

    bbMatrix <- bbMatrixAll[, c(uniqueIdCol[2], no.colsB)]
    colnames(bbMatrix) <- c( colnames(bbMatrixAll)[1], fkeyB)


    #--------------  total overlap --------------------
    AA=union(as.character(as.matrix(aaMatrix[, 1])), NULL)
    BB=union(as.character(as.matrix(bbMatrix[, 1])), NULL)
    #totaloverlap    = merge(AA,BB, by.x=1, by.y=1, sort=F, all=F)
    #no.totalOverlap = dim(totaloverlap)[1]
    totaloverlap    = intersect(AA, BB)
    no.totalOverlap = length(totaloverlap)

    totalbackground = length(union(AA, BB))
    if (!is.na(itotalGenes)){
        totalbackground = itotalGenes
    }else{
        totalbackground = length(AA)
    }
    
    # global enrichment test
    #
    pval = phyper(no.totalOverlap-1, length(AA), itotalGenes-length(AA), length(BB), lower.tail=F)
    globalEnrich = c(no.totalOverlap , pval, length(AA), length(BB))
    
    # matched to aaMatrix (mnore conservative test)
    #
    ccMatrix <- merge(bbMatrix, cbind(totaloverlap), by.x=1, by.y=1, all=F)
    ccMatrix <- as.matrix(ccMatrix)

    #-------------- module-based overlap --------------------
    #
    imoduleColA=dim(aaMatrix)[2]
    imoduleColB=dim(bbMatrix)[2]

    modlevelsA= names(table(aaMatrix[,imoduleColA]) )
    modlevelsB= names(table(bbMatrix[,imoduleColB]) )
    
    #modlevelsA= levels(aaMatrix[,imoduleColA]) #levels(aaMatrix$module)
    #modlevelsA
    #modlevelsB= levels(bbMatrix[,imoduleColB]) #levels(bbMatrix$module)
    #modlevelsB

    no.modulesA = length(modlevelsA)
    no.modulesB = length(modlevelsB)

    latexchline2 = rep("\\\\\\hline", no.modulesA + 2)
    latexchline  = rep("\\\\\\hline", no.modulesA + 1)

    internoMatrix = matrix(NA, no.modulesA+1, no.modulesB)
    pvalMatrix    = matrix(NA, no.modulesA+1, no.modulesB)
    pvalMatrix_Rst= matrix(NA, no.modulesA,   no.modulesB) # restricted test (to AA)

    freqMatrix     = matrix(0, no.modulesA+1, no.modulesB)
    freqMatrix_Rst = matrix(0, no.modulesA+1,   no.modulesB)

    genenoInModuelsA =c()
    genenoInModuelsB =c()
    allmergedmatrix=NULL

    # overlap between AA and individual BB modules
    #
    #overlap_bbmod2AA = rep(0, no.modulesB+ 1); overlap_bbmod2AA[1]=no.totalOverlap;
    AAX = cbind(AA, rep("ALL", length(AA)) )
    for (j in c(1:no.modulesB) ) {
        comrslt = compareTwoModules(AAX, bbMatrix, "ALL", modlevelsB[j], uniqueIdCol=1, 
                     moduleColA=2,moduleColB=2, totalGenes=totalbackground, removeDuplicate=removeDuplicate)
        interrslt =comrslt[[1]];
        internoMatrix[1,j] = interrslt[1]
        pvalMatrix [1,j]   = interrslt[2]
                
        if(restrict_networkB_to_A) {
           freqMatrix[1, j]    = length(BB)/totalbackground # frequency of this moduleB in the whole population
           freqMatrix_Rst[1, j]= no.totalOverlap/length(AA)  # frequency of this moduleB in the network A
        }else{
           freqMatrix_Rst[1, j]= interrslt[1]/length(AA)  # frequency of this moduleB in the network A
           freqMatrix[1, j]    = interrslt[4]/totalbackground # frequency of this moduleB in the whole population
        }
    }
    
    for (i in c(1:no.modulesA) ) {
      for (j in c(1:no.modulesB) ) {
        
        comrslt = compareTwoModules(aaMatrix, bbMatrix, modlevelsA[i], modlevelsB[j], uniqueIdCol=1, 
                        moduleColA=imoduleColA,moduleColB=imoduleColB, totalGenes=totalbackground, removeDuplicate=removeDuplicate)

        comrsltRst= compareTwoModules(gifA=aaMatrix, gifB=ccMatrix, moduleNameA=modlevelsA[i], 
                        moduleNameB=modlevelsB[j], uniqueIdCol=1, 
                        moduleColA=imoduleColA,moduleColB=imoduleColB, 
                        totalGenes=length(AA), removeDuplicate=removeDuplicate)

        interrslt =comrslt[[1]]; interrslt_Rst =comrsltRst[[1]]
        internoMatrix[i+1,j] = interrslt[1]
        pvalMatrix [i+1,j]   = interrslt[2]
        pvalMatrix_Rst[i,j]= interrslt_Rst[2]

        freqMatrix[i+1, j] =(interrslt[1]/interrslt[3])/freqMatrix[1, j];
        freqMatrix_Rst[i+1, j] =(interrslt_Rst[1]/interrslt_Rst[3])/freqMatrix_Rst[1, j];

        if (i== 1) {
           genenoInModuelsB =c(genenoInModuelsB, interrslt[4])
        }

        #if ( (interrslt[2] <= signifpvalue) & (modlevelsA[i]!="grey") & (modlevelsB[j]!="grey") ){
            allmergedmatrix=rbind(allmergedmatrix, comrslt[[2]])
        #}
      }
      genenoInModuelsA =c(genenoInModuelsA, interrslt[3])
    }

    #add in gene no in second column, names first column
    xinternoMatrix = cbind(c(length(AA),as.character(genenoInModuelsA)), internoMatrix)
    yinternoMatrix = cbind(c("ALL", as.character(modlevelsA)),       xinternoMatrix)

    #add in gene no in second row, names first row
    yinternoMatrix = rbind( c("total", as.character(no.totalOverlap), as.character(genenoInModuelsB)), yinternoMatrix)
    zinternoMatrix = rbind( c("overlap", "total", as.character(modlevelsB)), yinternoMatrix)

    write.table(zinternoMatrix, nofname,   sep="\t", quote=FALSE, col.names=F, row.names=F)

    if(latex) {
    zinternoMatrixLatex = zinternoMatrix
    for (i in c(1: (no.modulesA+2) ))
         zinternoMatrixLatex[i, no.modulesB+2] = paste(zinternoMatrixLatex[i, no.modulesB+2], "\\\\\\hline", sep="")
    write.table(zinternoMatrixLatex, latexnofname,   sep=" &", quote=FALSE, col.names=F, row.names=F)
    }

    #xpvalMatrix = cbind(as.character(genenoInModuelsA), pvalMatrix)
    #ypvalMatrix = cbind(as.character(modlevelsA),       xpvalMatrix)
    #ypvalMatrix = rbind( c("total", as.character(no.totalOverlap), as.character(genenoInModuelsB)), ypvalMatrix)
    #zpvalMatrix = rbind( c("overlap", "total", as.character(modlevelsB)), ypvalMatrix)

    # --------------------- pvalue --------------------------------
    #
    sigMatrix = pvalMatrix < signifpvalue
    sigMatrix_Rst = pvalMatrix_Rst < signifpvalue

    threshpMatrix = ifelse(sigMatrix, signif(pvalMatrix,3), "")
    threshpMatrix_Rst= ifelse(sigMatrix_Rst, signif(pvalMatrix_Rst,3), "")

    #zpvalMatrix = cbind(as.character(modlevelsA), signif(pvalMatrix,3))
    zpvalMatrix = cbind(c("ALL", as.character(modlevelsA)), threshpMatrix)
    zpvalMatrix = rbind( c("p.value", as.character(modlevelsB)), zpvalMatrix)
    write.table(zpvalMatrix,    pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F)

    if(restrict_networkB_to_A){
      appendStringToFile(pvalfname, "\nrestricted network B to network A\n")
    }else{appendStringToFile(pvalfname, "\nrestricted network B modules to network A\n")}

    zpvalMatrix = cbind(as.character(modlevelsA), threshpMatrix_Rst)
    zpvalMatrix = rbind( c("p.value", as.character(modlevelsB)), zpvalMatrix)
    write.table(zpvalMatrix, pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)


    # --------------------- frequency/fold-change --------------------------------
    #
    
    zpvalMatrix = cbind(c("ALL_freq", as.character(modlevelsA)),  freqMatrix)
    zpvalMatrix = rbind( c("freq/fold_change", as.character(modlevelsB)), zpvalMatrix)
    write.table(zpvalMatrix,    freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F)

    if(restrict_networkB_to_A){
      appendStringToFile(freqfname, "\nrestricted network B to network A\n")
    }else{appendStringToFile(freqfname, "\nrestricted network B modules to network A\n")}

    zpvalMatrix = cbind(c("ALL_freq", as.character(modlevelsA)), freqMatrix_Rst)
    zpvalMatrix = rbind( c("freq/fold_change", as.character(modlevelsB)), zpvalMatrix)
    write.table(zpvalMatrix, freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)


    if(latex) {
    zpvalMatrixLatex = zpvalMatrix
    for (i in c(1: (no.modulesA+1)) )
         zpvalMatrixLatex[i, no.modulesB+1] = paste(zpvalMatrixLatex[i, no.modulesB+1], "\\\\\\hline", sep="")
    write.table(zpvalMatrixLatex,    latexpvalfname, sep=" &", quote=FALSE, col.names=F, row.names=F)
    }

    if ( !is.null(allmergedmatrix)){
      colsmerged=dim(allmergedmatrix)[2]
      colnames(allmergedmatrix) <- c( colnames(allmergedmatrix)[c(1:(colsmerged-1))], keywords)
      write.table(allmergedmatrix, outfname, sep="\t", quote=FALSE, col.names=T, row.names=F)
      retname= outfname
    }else{
      retname=NULL
    }

    fn = paste("intersect_", keywords, sep="")
    #prefixes=fn; indir=outputDir; fname=fn
    merge_GO_PvalFreqNumb(prefixes=fn, indir=outputDir, fname=fn)

    retname
}


# to prevent genes that are not in allnames, we use intersection of input names and all the names 
#
fillInVector = function(mynames, indexvector, emptyvector, allnames)
{
   mynamesFound = intersect(mynames, allnames)
   iidx = indexvector[mynames]
   ovect= emptyvector  
   ovect[iidx] = TRUE
   return (ovect)
}

# restrict_networkB_to_A ==T: take networkB as background when mapping modulesB to network A
# restrict_networkB_to_A ==F: consider modulesB as independent, so map individual modules to network A
#
moduleBasedIntersectionMatrix=function(fnameA, fnameB,outputDir="", keywords="", uniqueIdCol=c(1,1), geneInforCols=8, genesymbolIdx=NULL, itotalGenes = 3600, signifpvalue=0.001, latex=F, removeDuplicate=F, restrict_networkB_to_A=T)
{
    #fkeyB       =getFileName(fnameB)
    #outfname    paste(outputDir, "intersect.csv", sep='')

    # get unique labels for each input files
    keys= unlist(strsplit(keywords, "-vs-"))
    fkeyA       =keys[1]
    fkeyB       =keys[2]

    #outfname     = paste(outputDir,  "intersect_", keywords, ".xls", sep='')
    outfname     = paste(outputDir,  "intersect_", keywords, ".xls", sep='')

    nofname     = paste(outputDir,  "intersect_", keywords, "_number.xls", sep='')
    pvalfname   = paste(outputDir,  "intersect_", keywords, "_pvalue.xls", sep='')
    freqfname   = paste(outputDir,  "intersect_", keywords, "_frequence.xls", sep='')

    pairfname   = paste(outputDir,  "intersect_", keywords, "_pair.xls", sep='')
    pairfnameR  = paste(outputDir,  "intersect_", keywords, "_pair-restricted.xls", sep='')


    latexnofname  = paste(outputDir, "intersect_", keywords, "_number.tex", sep='')
    latexpvalfname= paste(outputDir, "intersect_", keywords, "_pvalue.tex", sep='')

    aaMatrixAll <- read.delim(fnameA, sep="\t", header=T)
    dim(aaMatrixAll)
    aaMatrixAll <- as.matrix(aaMatrixAll)

    bbMatrixAll <- read.delim(fnameB, sep="\t", header=T)
    dim(bbMatrixAll)
    bbMatrixAll <- as.matrix(bbMatrixAll)

    aaMatrixAll[, uniqueIdCol[1]] = toupper(aaMatrixAll[, uniqueIdCol[1]])
    bbMatrixAll[, uniqueIdCol[2]] = toupper(bbMatrixAll[, uniqueIdCol[2]])

    no.colsA = dim(aaMatrixAll)[2]
    no.colsB = dim(bbMatrixAll)[2]

    # include gene infor in the aaMatrix but not bbMatrix, in this way the merged matrix
    #  will not have duplicated columns
    #
    aaMatrix <- aaMatrixAll[, c(uniqueIdCol[1], 1:geneInforCols, no.colsA)]
    colnames(aaMatrix) <- c( colnames(aaMatrixAll)[c(uniqueIdCol[1], 1:geneInforCols)], fkeyA)

    bbMatrix <- bbMatrixAll[, c(uniqueIdCol[2], no.colsB)]
    colnames(bbMatrix) <- c( colnames(bbMatrixAll)[1], fkeyB)

    imoduleColA=dim(aaMatrix)[2]
    imoduleColB=dim(bbMatrix)[2]

    # remove duplicates
    strA = paste(aaMatrix[,1], aaMatrix[,imoduleColA])
    idxA = findUniqueIndex(strA)
    if(length(idxA)<nrow(aaMatrix)) {
       aaMatrix = aaMatrix[idxA, ]
    }

    strB = paste(bbMatrix[,1], bbMatrix[,imoduleColB])
    idxB = findUniqueIndex(strB)
    if(length(idxB)<nrow(bbMatrix)) {
      bbMatrix = bbMatrix[idxB, ]
    }

    #--------------  total overlap --------------------
    AA=unique(as.character(aaMatrix[, 1]))
    BB=unique(as.character(bbMatrix[, 1]))
    #totaloverlap    = merge(AA,BB, by.x=1, by.y=1, sort=F, all=F)
    #no.totalOverlap = dim(totaloverlap)[1]
    totaloverlap    = intersect(AA, BB)
    no.totalOverlap = length(totaloverlap)

    totalbackground = length(union(AA, BB))
    if (!is.na(itotalGenes)){
        allnodes = sort(union(AA, BB))
        totalbackground = itotalGenes
    }else{
        totalbackground = length(AA)
        allnodes = sort(AA)
    }
    no.allnodes = length(allnodes)

    # global enrichment test
    #
    pval = phyper(no.totalOverlap-1, length(AA), itotalGenes-length(AA), length(BB), lower.tail=F)
    globalEnrich = c(no.totalOverlap , pval, length(AA), length(BB))
    
    # matched to aaMatrix (more conservative test)
    #
    ccMatrix <- merge(bbMatrix, cbind(totaloverlap), by.x=1, by.y=1, all=F)
    ccMatrix <- as.matrix(ccMatrix)

    #----------------- module-based overlap -------------------------
    #
    # restrict B to A
    restricted     = setElementInSet_Fast(bbMatrix[,1], allnodes)
    bbMatrixRstrct = bbMatrix[restricted, ]

    modlevelsA= names(table(aaMatrix[,imoduleColA]) )
    modlevelsB= names(table(bbMatrix[,imoduleColB]) )
    no.modulesA = length(modlevelsA)
    no.modulesB = length(modlevelsB)

    #-------------------make probe-categories matrix -----------------
    #
    # TRUE means a probe belongs to the corresponding category
    #
    allnodesIndex = c(1:no.allnodes); names(allnodesIndex) = allnodes
    ncA = ncol(aaMatrix); ncB = ncol(bbMatrix)

    empty    = rep(FALSE, no.allnodes)
    mres     = tapply(aaMatrix[, 1], INDEX=aaMatrix[, ncA], fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=allnodes)
    datMtrxA = data.frame(do.call(rbind, mres))
    datMtrxA = t(datMtrxA)
    dim(datMtrxA)

    mres     = tapply(bbMatrix[, 1], INDEX=bbMatrix[, ncB], fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=allnodes)
    datMtrxB = data.frame(do.call(rbind, mres))
    datMtrxB = t(datMtrxB)
    dim(datMtrxB)

    # find out module sizes in A, B, B restricted to A 
    sizes.modulesA = apply(datMtrxA, 2, sum)
    sizes.modulesBr= apply(datMtrxB, 2, sum)
    sizes.modulesB = table(bbMatrix[,imoduleColB])

    #****************** overlap **************************************
    ovlpMtrx = t(datMtrxA) %*% datMtrxB
    
    # match column names and rownames to the A module names and B module names
    midx = getMatchedIndex(names(sizes.modulesA), rownames(ovlpMtrx))
    sizes.modulesA = sizes.modulesA[midx]    
    midx = getMatchedIndex(names(sizes.modulesB), colnames(ovlpMtrx))
    sizes.modulesB = sizes.modulesB[midx]
    midx = getMatchedIndex(names(sizes.modulesBr), colnames(ovlpMtrx))
    sizes.modulesBr = sizes.modulesBr[midx]

    # prepare FET test vectors
    idxy = getMatrixIndex(size=dim(ovlpMtrx), symmetric=FALSE, diagonal=TRUE)
    no.tests = nrow(idxy)
    fetMtrx  = cbind(rep(totalbackground, no.tests), sizes.modulesA[idxy[,1]], sizes.modulesB[idxy[,2]],  ovlpMtrx[idxy])
    fetMtrxR = cbind(rep(no.allnodes,     no.tests), sizes.modulesA[idxy[,1]], sizes.modulesBr[idxy[,2]], ovlpMtrx[idxy])

    foldEnrich = (fetMtrx[,4]/fetMtrx[,3])*(fetMtrx[,1]/fetMtrx[,2])
    foldEnrichR= (fetMtrxR[,4]/fetMtrxR[,3])*(fetMtrxR[,1]/fetMtrxR[,2])

    # turn vector into matrix
    foldEnrichMtrx = matrix(foldEnrich, nrow=no.modulesA)
    foldEnrichMtrxR= matrix(foldEnrichR, nrow=no.modulesA)

    #************ fisher's test *********************************
    #
    fetP     = apply(fetMtrx,  1, fisherTest, minPopulationHits=0)
    fetPR    = apply(fetMtrxR, 1, fisherTest, minPopulationHits=0)

    # turn vector into matrix, be careful of the the order of vector into matrix
    fetPMtrx = t(matrix(fetP,  nrow=no.modulesB))
    fetPMtrxR= t(matrix(fetPR, nrow=no.modulesB))

    # output pair
    correction    = length(setdiff(names(sizes.modulesA), "grey")) *length(setdiff(names(sizes.modulesB), "grey")) 
    fetPcorrected =  fetP*correction; fetPcorrected = ifelse(fetPcorrected >1, 1, fetPcorrected)
    fetPcorrectedR= fetPR*correction; fetPcorrectedR= ifelse(fetPcorrectedR>1, 1, fetPcorrectedR)
    od = order(fetP)

    ocompares = paste(names(sizes.modulesA)[idxy[,1]], names(sizes.modulesB)[idxy[,2]])
    opairMtrx = cbind(ocompares, names(sizes.modulesA)[idxy[,1]], names(sizes.modulesB)[idxy[,2]],  fetMtrx,  foldEnrich,  fetP,  fetPcorrected )[od,]
    opairMtrxR= cbind(ocompares, names(sizes.modulesA)[idxy[,1]], names(sizes.modulesBr)[idxy[,2]], fetMtrxR, foldEnrichR, fetPR, fetPcorrectedR)[od,]
    tmpStrs = paste("module_size (", c(fkeyA, fkeyB), ")", sep="")
    colnames(opairMtrx) = c("Comparison", fkeyA, fkeyB,"population", tmpStrs, "Overlap", "Fold_Enrichment", "FET_P", "FET_P_corrected")
    colnames(opairMtrxR)= c("Comparison", fkeyA, fkeyB,"population", tmpStrs, "Overlap", "Fold_Enrichment", "FET_P", "FET_P_corrected")

    # -------- global FET test ----------------------
    globalFetMtrx = cbind(rep(totalbackground, no.modulesB), rep(length(AA)), 
                      rep(length(BB), no.modulesB),  sizes.modulesB)
    fetPglobal    = apply(globalFetMtrx, 1, fisherTest, minPopulationHits=0)
    foldEglobal   = (globalFetMtrx[,4]/globalFetMtrx[,3])*(globalFetMtrx[,1]/globalFetMtrx[,2])
    freqGlobal    = globalFetMtrx[,3]/totalbackground
    freqGlobalR   = globalFetMtrx[,3]/no.allnodes

    #-------------------------------------------------------------------------------------
    #------------ output AxB matrix results ----------------------------------------------
    #
    oNumberMtrx = rbind(sizes.modulesB, sizes.modulesBr, ovlpMtrx)
    oNumberMtrx = cbind(c("total", "All", names(sizes.modulesA)),
                        c(totalbackground, no.allnodes, sizes.modulesA), oNumberMtrx)
    colnames(oNumberMtrx) = c("overlap", "total", names(sizes.modulesB))

    #--- FET fold enrichment
    oFoldMtrx = rbind(freqGlobal, foldEnrichMtrx)
    oFoldMtrx = cbind(c("ALL_freq", names(sizes.modulesA)),oFoldMtrx)
    colnames(oFoldMtrx) = c("freq/fold_enrichment", names(sizes.modulesB))

    oFoldMtrxR= rbind(freqGlobalR, foldEnrichMtrxR)
    oFoldMtrxR= cbind(c("ALL_freq", names(sizes.modulesA)),oFoldMtrxR)
    colnames(oFoldMtrxR) = c("freq/fold_enrichment", names(sizes.modulesB))

    #--- FET p
    oFetPMtrx = rbind(fetPglobal, fetPMtrx)
    oFetPMtrx = cbind(c("All", names(sizes.modulesA)),oFetPMtrx)
    colnames(oFetPMtrx) = c("pvalue", names(sizes.modulesB))

    oFetPMtrxR= fetPMtrxR
    oFetPMtrxR= cbind(names(sizes.modulesA),oFetPMtrxR)
    colnames(oFetPMtrxR) = c("pvalue", names(sizes.modulesBr))

    #-------------------------------------------------------------------------------------
    #------------ save Matrix form results into files ------------------------------------
    #
    write.table(rbind(colnames(oNumberMtrx)), nofname,   sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(oNumberMtrx, nofname,   sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

    write.table(rbind(colnames(oFetPMtrx)), pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(oFetPMtrx, pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)
    if(restrict_networkB_to_A){
      appendStringToFile(pvalfname, "\nrestricted network B to network A\n")
    }else{appendStringToFile(pvalfname, "\nrestricted network B modules to network A\n")}
    write.table(rbind(colnames(oFetPMtrx)), pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)
    write.table(oFetPMtrxR, pvalfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)

    write.table(rbind(colnames(oFoldMtrx)), freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(oFoldMtrx, freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)
    if(restrict_networkB_to_A){
      appendStringToFile(pvalfname, "\nrestricted network B to network A\n")
    }else{appendStringToFile(freqfname, "\nrestricted network B modules to network A\n")}
    write.table(rbind(colnames(oFoldMtrx)), freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)
    write.table(oFoldMtrxR, freqfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)

    #-------------------------------------------------------------------------------------
    #------------ find overlap genes for each group pair          ------------------------
    #
    #print("Merging two matrices")
    selA = aaMatrix[,ncol(aaMatrix)] != "grey"
    selB = bbMatrix[,ncol(bbMatrix)] != "grey"

    abMtrx = merge(aaMatrix[selA, ], bbMatrix[selB,c(1, ncol(bbMatrix))], by.x=1, by.y=1, all=FALSE)
    abMtrx = as.matrix(abMtrx); ncx = ncol(abMtrx)
    comparison = paste(abMtrx[, ncx-1], abMtrx[, ncx])
    abMtrx2= cbind(abMtrx, comparison)
    write.table(abMtrx2, outfname, sep="\t", quote=FALSE, col.names=TRUE, row.names=F)

    mergedMembers = tapply(abMtrx[,1], INDEX=comparison, concatenate, mysep=",", do_union=TRUE, do_sort=TRUE)

    selGS = !is.null(genesymbolIdx)
    if(selGS){selGS = genesymbolIdx!=uniqueIdCol[1] }

    if(selGS){
       mergedMembersGS= tapply(abMtrx[,genesymbolIdx+1], INDEX=comparison, concatenate, mysep=",", do_union=TRUE, do_sort=TRUE)
       mergedMembersMtrx = cbind(names(mergedMembers), mergedMembers, mergedMembersGS)
       colnames(mergedMembersMtrx) = c("Comparison", "SharedProbes", "SharedGenes")
    } else {
       mergedMembersMtrx = cbind(names(mergedMembers), mergedMembers)
       colnames(mergedMembersMtrx) = c("Comparison", "SharedMembers")
    }
    #mergedMembersModNames = getAllParts(names(mergedMembers), " ")
    #indexnamesA = c(1:length(sizes.modulesA)); names(indexnamesA) = names(sizes.modulesA)
    #indexnamesB = c(1:length(sizes.modulesB)); names(indexnamesB) = names(sizes.modulesB)
    #mergedModNamesIdxy = cbind(indexnamesA[mergedMembersModNames[,1]], indexnamesB[mergedMembersModNames[,2]])
    
    opairMtrx2 = merge(opairMtrx,  mergedMembersMtrx, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
    opairMtrx2R= merge(opairMtrxR, mergedMembersMtrx, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)

    pidx = ncol(opairMtrx) - (ncol(mergedMembersMtrx)-1) -1
    tmpP = as.numeric(opairMtrx[, pidx]); od = order(tmpP)

    write.table(rbind(colnames(opairMtrx2)), pairfname, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(opairMtrx2, pairfname, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)

    write.table(rbind(colnames(opairMtrx2R)), pairfnameR, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(opairMtrx2R, pairfnameR, sep="\t", quote=FALSE, col.names=F, row.names=F, append=T)

}


# Test the modules in the files for enrichment of gene ontology terms
moduleBasedIntersectionMatrix_GeneOntology = function(fnames, fnameB, outputDir="", uniqueIdCol=1, ctopNrows=10, signifLevel=0.05, 
   removeDuplicate=TRUE, removeGrey=TRUE)

{
    ########################### GO Matrix ##################################################
    #
    #System	GeneCategory	PopulationHits	PopulationTotal	GeneSymbol
    #GO Biological Process	cell-cell signaling	933	16773	NAALAD2; NAALADL1; ...
    #GO Biological Process	protein biosynthesis	912	16773	PIGK; FARSLB; NR1H3; ...
    #
    bbMatrixAll <- read.delim(fnameB, sep="\t", header=T)
    bbMatrixAll <- as.matrix(bbMatrixAll)
    dim(bbMatrixAll)

    print("Build up GO term membership matrix ....")

    no.colsB = ncol(bbMatrixAll)
    newIDs      <- paste(bbMatrixAll[, 1], bbMatrixAll[, 2], sep=";;;")
    termGenes   <- tapply(bbMatrixAll[, no.colsB], INDEX=newIDs, FUN=getAllPartsVect, sep="; ")
    termGenes   <- lapply(termGenes, unique)
    termGenesLen<- unlist(lapply(termGenes, length))

    totalbackground = unique(unlist(termGenes))
    no.allnodes     = length(totalbackground)
    allnodesIndex = c(1:no.allnodes); names(allnodesIndex) = totalbackground 

    empty    = rep(FALSE, no.allnodes)

    mres     = lapply(termGenes, fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=totalbackground)
    datMtrxB = data.frame(do.call(rbind, mres))
    datMtrxB = t(datMtrxB)
    dim(datMtrxB)

    modlevelsB= names(termGenes)
    no.modulesB = length(modlevelsB)
    sizes.modulesB = apply(datMtrxB, 2, sum)

    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    #
    for(fnameA in fnames) {

    print(fnameA)

    fname         = getFileNameNopath(fnameA)
    outfnameTop   = paste(outputDir,  fname, "_Ontology.xls", sep='')
    outfnameTop10 = paste(outputDir,  fname, "_OntologyTop10.xls", sep='')
    outfnameTop10m= paste(outputDir,  fname, "_OntologyTop10-sortByModule.xls", sep='')

    aaMatrixAll <- read.delim(fnameA, sep="\t", header=T)
    aaMatrixAll <- as.matrix(aaMatrixAll)
    dim(aaMatrixAll)
    modProbeSizes <- tapply(aaMatrixAll[,1], INDEX=aaMatrixAll[,ncol(aaMatrixAll)], FUN=length)

    if(removeGrey) { # remove grey module
        gsel = aaMatrixAll[, ncol(aaMatrixAll)] != "grey"
        aaMatrixAll = aaMatrixAll[gsel,]
    }

    #aaMatrixAll[, uniqueIdCol[1]] = toupper(aaMatrixAll[, uniqueIdCol[1]])

    no.colsA = dim(aaMatrixAll)[2]

    # include gene infor in the aaMatrix but not bbMatrix, in this way the merged matrix
    #  will not have duplicated columns
    #
    aaMatrix <- aaMatrixAll[, c(uniqueIdCol[1], no.colsA)]
    colnames(aaMatrix) <- colnames(aaMatrixAll)[c(uniqueIdCol,no.colsA)]


    #------ remove duplicates in each module -------------------
    #
    imoduleColA=dim(aaMatrix)[2]
    strA = paste(aaMatrix[,1], aaMatrix[,imoduleColA])
    idxA = findUniqueIndex(strA)
    if(length(idxA)<nrow(aaMatrix)) {
       aaMatrix = aaMatrix[idxA, ]
    }
    
    #----------------- module-based overlap -------------------------
    #
    # restrict A to B (go terms)
    #
    restricted     = setElementInSet_Fast(aaMatrix[,1], totalbackground)
    aaMatrixRstrct = aaMatrix[restricted, ]

    modlevelsA= sort(unique(aaMatrixRstrct[,imoduleColA]) )
    no.modulesA = length(modlevelsA)

    #-------------------make probe-categories matrix -----------------
    #
    # TRUE means a probe belongs to the corresponding category
    #
    ncA = ncol(aaMatrixRstrct); 

    print("Build up module membership matrix ....")
    mres     = tapply(aaMatrixRstrct[, 1], INDEX=aaMatrixRstrct[, ncA], fillInVector, indexvector=allnodesIndex, emptyvector=empty, allnames=totalbackground )
    datMtrxA = data.frame(do.call(rbind, mres))
    datMtrxA = t(datMtrxA)
    dim(datMtrxA)

    # find out module sizes in A, B, B restricted to A 
    sizes.modulesA = apply(datMtrxA, 2, sum)
    sum(termGenesLen!= sizes.modulesB)

    #****************** overlap **************************************
    print("Perform FET test ....")

    ovlpMtrx = t(datMtrxA) %*% datMtrxB
    
    # match column names and rownames to the A module names and B module names
    midx           = getMatchedIndexFast(names(sizes.modulesA), rownames(ovlpMtrx))
    sizes.modulesA = sizes.modulesA[midx]    
    midx           = getMatchedIndexFast(names(sizes.modulesB), colnames(ovlpMtrx))
    sizes.modulesB = sizes.modulesB[midx]

    # prepare FET test vectors
    idxy       = getMatrixIndex(size=dim(ovlpMtrx), symmetric=FALSE, diagonal=TRUE)
    no.tests   = nrow(idxy)
    fetMtrx    = cbind(rep(no.allnodes, no.tests), sizes.modulesA[idxy[,1]], sizes.modulesB[idxy[,2]],  ovlpMtrx[idxy])
    foldEnrich = (fetMtrx[,4]/fetMtrx[,3])*(fetMtrx[,1]/fetMtrx[,2])

    # turn vector into matrix
    foldEnrichMtrx = matrix(foldEnrich, nrow=no.modulesA)

    #*************** fisher's test *********************************
    #
    fetP     = apply(fetMtrx, 1, fisherTest, minPopulationHits=0)
    od       = order(fetP)

    # turn vector into matrix, be careful of the the order of vector into matrix
    fetPMtrx = t(matrix(fetP,  nrow=no.modulesB))

    # output pair
    correction    = no.modulesA * no.modulesB
    fetPcorrected = fetP*correction
    fetPcorrected = ifelse(fetPcorrected >1, 1, fetPcorrected)
    fetPcorrected = signif(fetPcorrected, 2)

    #-------------------------------------------------------------------------------------
    #------------ find overlap genes for each group pair          ------------------------
    #

    # to pull out only the overlaps from the significant tests
    selMtrx = fetPMtrx < signifLevel
    selMod  = (apply(selMtrx, 1, sum))>=1;
    selGot  = (apply(selMtrx, 2, sum))>=1;

    print(sprintf("Get overlap genes for %s modules and %s terms ....", sum(selMod), sum(selGot) ) )

    #tmp=concatenateSelectedElelments_vectMatrix(datMtrxA[,178], datMtrxB, mnames=totalbackground, sep=",")
    #mergedMembersGS = apply(datMtrxA, 2, concatenateSelectedElelments_vectMatrix, boolMatrix=datMtrxB, mnames=totalbackground, sep=",")

    mergedMembersGS = apply(datMtrxA[,selMod], 2, concatenateSelectedElelments_vectMatrix, boolMatrix=datMtrxB[,selGot], mnames=totalbackground, sep=",")
    mergedMembersGS = t(mergedMembersGS )

    idxy2      = getMatrixIndex(size=dim(mergedMembersGS), symmetric=FALSE, diagonal=TRUE)
    membersVect= mergedMembersGS[idxy2]
    names.membersVect = paste(rownames(mergedMembersGS)[idxy2[,1]], colnames(mergedMembersGS)[idxy2[,2]])

    # coimbine results from all the pair-wise tests
    xmodnames  = names(sizes.modulesA)[idxy[,1]]
    xmodSizesP = modProbeSizes[xmodnames]
    goSrcTerms = getAllParts(names(sizes.modulesB), ";;;")
    xgoSrcTerms= goSrcTerms[idxy[,2],]
    xgoTermSzs = cbind(sizes.modulesB[idxy[,2]], rep(no.allnodes, no.tests) )

    # matched the significant pairs to the all apirs
    names.alltests = paste(xmodnames, names(sizes.modulesB)[idxy[,2]])
    midx = getMatchedIndexFast(names.alltests, names.membersVect)
    sum(is.na(midx))
    membersVect2 = rep("", no.tests)
    membersVect2[midx] = membersVect

    ocompares = cbind(xmodnames,xmodSizesP, xgoSrcTerms)
    opairMtrx = cbind(ocompares, ovlpMtrx[idxy], sizes.modulesA[idxy[,1]], xgoTermSzs, fetP,  fetPcorrected, foldEnrich,  membersVect2)
    colnames(opairMtrx) = c("module", "ModuleSize", "System", "Gene.Category", "ModuleOverlap", "GeneSetCount", 
                            "PopulationOverlap", "Population", "FET_P", "Corrected_P", "Fold_Enrichment", "Genes")
    opairMtrx = opairMtrx[od,]

    pidx = getMatchedIndex(colnames(opairMtrx), "FET_P")
    tmpP = as.numeric(opairMtrx[, pidx]);     
    xsel = tmpP<signifLevel
    opairMtrx2 = opairMtrx[xsel,]
    write.table(rbind(colnames(opairMtrx2)), outfnameTop, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(opairMtrx2,                  outfnameTop, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

    modIndices = tapply(c(1:nrow(opairMtrx2)), INDEX=opairMtrx2[,1], return_original, topN=ctopNrows)
    modIndicesV= unlist(modIndices)
    opairMtrxTop = opairMtrx2[modIndicesV, ]

    write.table(rbind(colnames(opairMtrxTop)), outfnameTop10m, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(opairMtrxTop,                  outfnameTop10m, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

    tmpP = as.numeric(opairMtrxTop[, pidx]);
    od   = order(tmpP)
    write.table(rbind(colnames(opairMtrxTop)), outfnameTop10m, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(opairMtrxTop[od,],                  outfnameTop10m, sep="\t", quote=FALSE, col.names=F, row.names=F, append=TRUE)

    }

}











# here the saving process is decomposed into many small steps, each of which saves only
# a small matrix into the file so that it doesn't need a lot of memory for the whole operation
# but with a little bit more time
#
saveHugeMatrix =function(hugematrix, titlerow=NULL, outfilename,  use1minus=F){

    #title row
    if (!is.null(titlerow)){
       write.table(t(as.matrix(titlerow)), outfilename, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    }

    # write into file in many times for saving memory
    no.rows=dim(hugematrix)[1]
    step   =100
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=hugematrix[irows, ]

       if (use1minus){
         write.table(1-imatrix, outfilename, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
       }else{
         write.table(imatrix, outfilename, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
       }
    }

}

#******************************************************************************************
#********************* FOR converting pair link list into adjacency matrix ****************
#******************************************************************************************
#

############################## links pairs to adj maTRIX #########################
# merge to get indices of two genes forming a link
# netmatrix, mappingMatrix are bot Nx2 matrix
# netmatrix, each row contains source and dest nodes
# mappingMatrix, each rwo represents the mapping between node name and index
#
getNetworknodeIdxByNames = function(netmatrix, mappingMatrix)
{
  # node on the left side
  mergedleft = merge(netmatrix, mappingMatrix, by.x=1, by.y=1, all=F)

  # node on the right side
  mergedright = merge(mergedleft, mappingMatrix, by.x=2, by.y=1, all=F)
  mergedright = as.matrix(mergedright)
  indices2D = cbind(as.integer(mergedright[,3]), as.integer(mergedright[,4]) )
  indices2D
}

# note that name2idxMatrix is a global variable
# coding == NULL: use the third column in inetmatrix
#
makeAjacencyMatrix = function(inetmatrix, coding, matrixsize, directed=F, myname2idxMatrix)
{
  # initialization, no links between any two nodes
  adjMatrix = matrix(0, ncol=matrixsize, nrow=matrixsize)

  # find nodes of each edge, and replace the corresponding element with the edge code
  edgeIndices=getNetworknodeIdxByNames(netmatrix=inetmatrix[,c(1,2)], mappingMatrix=myname2idxMatrix)
  if( !is.null(coding) ) {
     adjMatrix[edgeIndices] = coding[2]
  }else{
     adjMatrix[edgeIndices] = inetmatrix[,3]
  }

  if(!directed){ # for undirected graph, the ajacency matrix is symmetric
    transposedIndices = edgeIndices[,c(2,1)]

    if( !is.null(coding) ) {
       adjMatrix[transposedIndices] = coding[2]
    }else{
       adjMatrix[transposedIndices] = inetmatrix[,3]
    }
  }

  adjMatrix
}


# return adj lists for all nodes
# if returnUniqueNames==T, return only the node names
#
# notice that a node without a link is marked as -1
#
makeAjacencyLists=function(inetmatrix, returnUniqueNames=F, directed=F)
{
    edgesInNet = dim(inetmatrix)[1]

    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(inetmatrix[,1]) )
    allnodenames = c(allnodenames, as.character(inetmatrix[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    if (returnUniqueNames){
       return ( list(uniquenames, no.uniquenames) )
    }

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )


    # initialization, no links between any two nodes
    #
    adjlists = as.list( rep(-1, no.uniquenames) )

    # find nodes of each edge, and replace the corresponding element with the edge code
    edgeIndices = getNetworknodeIdxByNames(netmatrix=inetmatrix[,c(1,2)], mappingMatrix=name2idxMatrix)
    no.edges    = dim(edgeIndices)[1]  

    for (k in c(1:no.edges) ){
       adjlists[[ edgeIndices[k,1] ]] = c(adjlists[[ edgeIndices[k,1] ]], edgeIndices[k,2])
       if (!directed){
          adjlists[[ edgeIndices[k,2] ]] = c(adjlists[[ edgeIndices[k,2] ]], edgeIndices[k,1])
       }
    }

    # remove redundant elements & -1
    #
    for (i in c(1:no.uniquenames)){
       if ( length(adjlists[[i]]) >1){
          adjlists[[i]]=setdiff(adjlists[[i]], -1)
       }
    }

    return (adjlists)
}


 
# input link pairs are node indices, so we construct 
# adjlist for the nodes from 1 to the maximum index 
#
makeAjacencyListsFromLinksIndex = function(linksindex, excludedNodes=NULL)
{
    no.edges = dim(linksindex)[1]
    maxIndex = max( max(linksindex[,1]), max(linksindex[,2]) )

    # initialization, no links between any two nodes
    #
    adjlists = as.list( rep(-1, maxIndex) )

    for (k in c(1:no.edges) ){
       adjlists[[ linksindex[k,1] ]] = c(adjlists[[ linksindex[k,1] ]], linksindex[k,2])
       adjlists[[ linksindex[k,2] ]] = c(adjlists[[ linksindex[k,2] ]], linksindex[k,1])
    }

    # remove redundant -1
    #
    for (i in c(1:maxIndex )){
       if ( length(adjlists[[i]]) >1){
          adjlists[[i]]=setdiff(adjlists[[i]], -1)
       }
    }

    return (adjlists)
}

# make adj lists for a subset of nodes
#
makeAjacencyListsFromSubsetnodesindex = function(orgAdjLists, subsetNodes=NULL)
{
    if(is.null(subsetNodes)) {
       return(orgAdjLists)
    }

    no.orgnodes = length(orgAdjLists)
    
    # initialization, no links between any two nodes
    #
    adjlists = as.list( rep(-1, no.orgnodes) )
    selNodes = rep(F, no.orgnodes)
    selNodes = T

    for (i in  subsetNodes){
       ioverlap = intersect(orgAdjLists[[i]],subsetNodes)
       if( length(ioverlap)>0) {
          adjlists[[ i ]] = ioverlap
       }
    }

    return (adjlists)
}


# convert dichotmize matrix into link pairs
#

dichoMatrix_to_linkpairs=function(dichotCor, genenames, pairExprnet)
{
  no.genes = dim(dichotCor)[1]
  
  for (i in 1:(no.genes) ){
     noleft     = no.genes - i
     iselpool   = c(rep(F, i), rep(T, noleft) )
     iconnected = (dichotCor[i,] & iselpool)

     if(sum(iconnected )==0){
        next
     }

     ipairs     = paste(genenames[i], "\t", genenames[iconnected],sep="") 
  
     if(i==1){
        write.table( cbind(ipairs), pairExprnet, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
     }else{
        write.table( cbind(ipairs), pairExprnet, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
     }
  }
}


anyMatrix_to_linkpairs=function(anyMatrix, irownames, icolnames, pairExprnet)
{
  no.genes   = dim(anyMatrix)[1]
  no.samples = dim(anyMatrix)[2]
  
  for (i in 1:(no.genes) ){
  
     iconnected = (anyMatrix[i,] ==1)

     if(sum(iconnected )==0){
        next
     }

     ipairs     = paste(irownames[i], "\t", icolnames[iconnected],sep="") 
  
     if(i==1){
        write.table( cbind(ipairs), pairExprnet, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
     }else{
        write.table( cbind(ipairs), pairExprnet, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
     }
  }
}



#********************************** INPUT **********************************
# each gene pair in "netpairMatrix" represents a link with two nodes
# for directed networks, the first node (element) is the source and 
#    the second is the destination
#
# codepair=c(0,1):  #[1] for no connection, [2] for connection; ==NULL, use confidence value from the third column of netpairMatrix
# directed= T: #TRUE    #FALSE # directed network or not

convert_links2adjmatrix= function(netpairMatrix,directed= T, codepair=c(0,1), keywords="", fadjname=NULL) {

    #keywords = getFileName(fname)
    #fadjname = paste(outputDir, keywords, ".adj", sep="")

    edgesInNet = dim(netpairMatrix)[1]

    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(netpairMatrix[,1]) )
    allnodenames = c(allnodenames, as.character(netpairMatrix[,2]) )
     
    #nametable = table(allnodenames)
    #length(nametable)
    #uniquenames     = names(nametable)
    
    uniquenames     = union(allnodenames , NULL)
    uniquenames     = sort(uniquenames)

    no.uniquenames  = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )


    adjmatrix = makeAjacencyMatrix(inetmatrix=netpairMatrix, coding=codepair,
                                     matrixsize=no.uniquenames, directed=directed, 
                                     myname2idxMatrix = name2idxMatrix)
    colnames(adjmatrix) <- uniquenames
    rownames(adjmatrix) <- uniquenames


    #edgeIndices=getNetworknodeIdxByNames(netmatrix=netpairMatrix[,c(1,2)], mappingMatrix=name2idxMatrix)

    if(is.null(fadjname)){
      return (adjmatrix)
    }


    #---------------------------- output adj matrix ------------------------------------------
    #title row
    write.table(t(as.matrix(uniquenames)), fadjname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)

    # write into file in many times for saving memory
    no.rows=dim(adjmatrix)[1]
    step=200
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=adjmatrix[irows, ]

       write.table(imatrix, fadjname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

    #---------------------------- output statistics --------------------------------------------
    # matched percentage
    #
    linkspernode = as.integer(100*edgesInNet/no.uniquenames)/100
    matchMatrix  = cbind(keywords, edgesInNet, no.uniquenames, linkspernode)
    fmatrix      = rbind( c("network", "# of edges", "# of nodes", "edges per node"), matchMatrix)

    #write.table(fmatrix, flogname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    fmatrix
}

# sna format: 
#*vertices 760		
#nodename	color	shape
#YAL005C	white	50
#YAL034C	white	50
#......
#YAL005C	YAL034C	YJL141C	YAL054C
#0	0	0	0
#0	0	1	0


# normColor/normShape for default nodes
# highColor/highShape for highlighted nodes (signature)
# koColor/koShape for knockout
#
#
makeSNA_from_netpairs=function(netpairs, highlightNodes=NULL, 
           knockouts=NULL,
           koColor  ="green",koShape  ="4",
           normColor="grey", highColor="red", 
           normShape="50",   highShape="50", directed= T, snafile="tmp.sna")
{
    edgesInNet = dim(netpairs)[1]

    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(netpairs[,1]) )
    allnodenames = c(allnodenames, as.character(netpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = sort(names(nametable))
    no.uniquenames  = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    
    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix
    verticesMatrix = NULL
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames) )
    } else{
      #shape
      if (length(highShape)==1){#use the same shape for all nodes
          nodesHighShape = cbind(highlightNodes,
                                 rep(highShape,length(highlightNodes)) )
      }else{
          nodesHighShape = cbind(highlightNodes,highShape)
      }

      # color
      if (length(highColor)==1){#use the same shape for all nodes
          nodesHighColor = cbind(highlightNodes,
                                 rep(highColor,length(highlightNodes)) )
      }else{
          nodesHighColor = cbind(highlightNodes,highColor)
      }

      #align highlighted nodes to the network index
      #
      colored = mergeTwoMatricesByKeepAllPrimary(cbind(uniquenames),nodesHighColor, 
                                            missinglabel=normColor,
                                            keepAllPrimary=T, keepPrimaryOrder=T)
      verticesMatrix = mergeTwoMatricesByKeepAllPrimary(colored,nodesHighShape, 
                                            missinglabel=normShape,
                                            keepAllPrimary=T, keepPrimaryOrder=T)
    }
    verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape")

    # change the settings of knockout
    #
    if( !is.null(knockouts) ){
       for(eko in knockouts) {
          if(!is.element(eko, uniquenames)){
             next
          }
          koIdx = c(1:no.uniquenames)[eko == uniquenames]
          verticesMatrix[koIdx,2]= koColor
          verticesMatrix[koIdx,3]= koShape
       }
    }

    # 3. adjacency matrix
    adjmatrix = makeAjacencyMatrix(inetmatrix=netpairs, coding=c(0,1),
                                     matrixsize=no.uniquenames, directed=directed, 
                                     myname2idxMatrix = name2idxMatrix)
    
    #edgeIndices=getNetworknodeIdxByNames(netmatrix=netpairMatrix[,c(1,2)], mappingMatrix=name2idxMatrix)

   
    #---------------------------- output adj matrix ------------------------------------------
    #
    write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(verticesMatrix, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    # title row
    write.table(t(as.matrix(uniquenames)), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # write into file in many times for saving memory
    no.rows=dim(adjmatrix)[1]
    step=200
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=adjmatrix[irows, ]

       write.table(imatrix, snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SNA from matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# nodenames are sort and in the same order as in the col/row names
# adjmatrix define the color of edges or simply connected or not
#
makeSNA_from_Adjmatrix=function(adjmatrix, nodenames, highlightNodes=NULL, 
           edgecolorlevels,
           normColor="grey", highColor="red", 
           normShape="50",   highShape="50", directed= T, snafile="tmp.sna")
{
    uniquenames     = nodenames
    no.uniquenames  = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    
    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix
    verticesMatrix = NULL
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames) )
    } else{
      #shape
      if (length(highShape)==1){#use the same shape for all nodes
          nodesHighShape = cbind(highlightNodes,
                                 rep(highShape,length(highlightNodes)) )
      }else{
          nodesHighShape = cbind(highlightNodes,highShape)
      }

      # color
      if (length(highColor)==1){#use the same shape for all nodes
          nodesHighColor = cbind(highlightNodes,
                                 rep(highColor,length(highlightNodes)) )
      }else{
          nodesHighColor = cbind(highlightNodes,highColor)
      }

      #align highlighted nodes to the network index
      #
      colored = mergeTwoMatricesByKeepAllPrimary(cbind(uniquenames),nodesHighColor, 
                                            missinglabel=normColor,
                                            keepAllPrimary=T, keepPrimaryOrder=T)
      verticesMatrix = mergeTwoMatricesByKeepAllPrimary(colored,nodesHighShape, 
                                            missinglabel=normShape,
                                            keepAllPrimary=T, keepPrimaryOrder=T)
    }
    verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape")
    
   
    #---------------------------- output adj matrix ------------------------------------------
    #
    write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(verticesMatrix, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    #edge color
    write.table(t(as.matrix(c("edge_color_levels", edgecolorlevels) )), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # title row
    write.table(t(as.matrix(uniquenames)), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # write into file in many times for saving memory
    no.rows=dim(adjmatrix)[1]
    step=200
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=adjmatrix[irows, ]

       write.table(imatrix, snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SNP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# the third column of netpairsWtype is the link type 
#  
# highlightNodes can be a vector of genes or a list of gene sets with
#    highColor and highShape for each gene set, i.e.,
#
#  highlightNodes[[i]] <-  highColor[i] and highShape[i]
#
makeSNP =function(netpairsWtype, highlightNodes=NULL,
           edgecolorlevels,
           normColor="grey", highColor="red", 
           normShape="50",   highShape="50", 
           normNodeSize ="1",    highNodeSize="2",
           normFontSize ="1",    highFontSize="2",
           directed= T, legendtable=NA, snafile="tmp.sna")
{

    xfname = getFileName(snafile)
    fcys  = paste(xfname, "_cys.txt", sep="")
    fcysn = paste(xfname, "_cys-nodes.txt", sep="")
    write.table(netpairsWtype, fcys, sep="\t", quote=FALSE, col.names=T, row.names=FALSE)

    pcols = dim(netpairsWtype)[2]

    # consider both columns
    uniquenames    = union(as.character(netpairsWtype[,1]), as.character(netpairsWtype[,2]) )
    uniquenames    = sort(uniquenames)
    no.uniquenames = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    

    # 0. make link index: A1 A2 T & A I
    # A1 A2 T & A I ==> A1 A2 T I1
    #
    leftIdx = merge(netpairsWtype, name2idxMatrix, by.x=1, by.y=1, all=F) 
    
    #  A1 A2 T I1 & A I ==> A2 A1 T I1 I2: I1 and I2 are the indices of A1 and A2 respectively
    #
    allIdx  = merge(leftIdx, name2idxMatrix, by.x=2, by.y=1, all=F)

    no.pairs = dim(allIdx)[1]

    if(pcols==2){
      no.elevels = length(edgecolorlevels)
      greyIdx    = c(1:no.elevels) [edgecolorlevels=="grey"]
      linksIdxMatrix = cbind( allIdx[,c(3,4)], rep(greyIdx,no.pairs) )
    } else{
      linksIdxMatrix = allIdx[,c(4,5,3)]
    }


    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix 
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames),
                               rep(normNodeSize, no.uniquenames),
                               rep(normFontSize, no.uniquenames))
    } else{
      verticesMatrix = matrix("", no.uniquenames, 5)
      verticesMatrix[,1] = uniquenames
      verticesMatrix[,2] = rep(normColor,no.uniquenames)
      verticesMatrix[,3] = rep(normShape,no.uniquenames)
      verticesMatrix[,4] = rep(normNodeSize, no.uniquenames)
      verticesMatrix[,5] = rep(normFontSize, no.uniquenames)

      xcolor= rep(normColor,no.uniquenames);    names(xcolor) <- uniquenames
      xshape= rep(normShape,no.uniquenames);    names(xshape) <- uniquenames
      xNsize= rep(normNodeSize,no.uniquenames); names(xNsize) <- uniquenames
      xSsize= rep(normFontSize,no.uniquenames); names(xSsize) <- uniquenames

      # set color and shape for highlighted nodes
      #
      if ( !is.list(highlightNodes) ){
          highlightNodes2 = intersect(highlightNodes, uniquenames)
          xcolor[ highlightNodes2] = highColor[1]
          xshape[ highlightNodes2] = highShape[1]
          xNsize[ highlightNodes2] = highNodeSize[1]
          xSsize[ highlightNodes2] = highFontSize[1]
      }else{
         no.highnodeSets = length(highlightNodes)
         for(il in c(1:no.highnodeSets) ) {
             highlightNodes2 = intersect(highlightNodes[[il]], uniquenames)
             if(length(highlightNodes2)==0){next};
             xcolor[highlightNodes2] = highColor[il]
             xshape[highlightNodes2] = highShape[il]
             xNsize[highlightNodes2] = highNodeSize[il]
             xSsize[highlightNodes2] = highFontSize[il]
         }
      }
      verticesMatrix[,2] = as.character(xcolor)
      verticesMatrix[,3] = as.character(xshape)
      verticesMatrix[,4] = as.character(xNsize)
      verticesMatrix[,5] = as.character(xSsize)
    }

    #verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape", "size", "font_size")

    write.table(verticesMatrix, fcysn, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

    #**************************************************************************
    #
    # 3. output indexed netpairs
    #
    # Legend
    if ( !is.na(legendtable) ) {
      mhead = paste("Legend", dim(legendtable)[1], sep=" ")
      write.table(as.matrix(mhead), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
      write.table(legendtable, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)
      # vertex
      write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    } else{
      # vertex
      write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    }

    write.table(verticesMatrix, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    #edge color
    write.table(t(as.matrix(c("edge_color_levels", edgecolorlevels) )), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # link pairs based index   
    #
    write.table(rbind(c("src","dst", "type")), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    write.table(linksIdxMatrix, snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    return (c(fcys, fcysn))
}


if(F){  # change A at positions on (row,col) pairs with value on the third column
  a=matrix(1,3,3)
  a[1,2]=3
  a[3,3]=4
  a[2,1]=2
  a[1,3]=-1
  a[2,3]=-2
  a[3,3]=-3
  b[a[,c(1,2)]] = a[,3]
}

# configure nodes's shapes and colors for cytoscape visualization
#
# allnodes: all nodes in the network for visusalization
# signature: signature genes
# keydriverMatrix: results from KDA (first column is geneID and last column is key-driver status  
#                   with 1 for global driver and 0 for local drivers)
#

configure_node_visualization = function(allnodes, signature, kdaMatrix, bNodeSz=40, bFontSz=12) {

   # SIG--signature; NSIG--not signature; GKD--Global KeyDriver; LKD--Local KeyDriver; NKD--Not KeyDriver
   #
   xcategories = c("SIG_GKD", "SIG_LKD", "SIG_NKD", "NSIG_GKD", "NSIG_LKD", "NSIG_NKD"); xcat2=c("NSIG_GKD", "NSIG_LKD", "NSIG_NKD")
   xcolors     = c("red",     "blue",    "lightgreen",  "red",      "blue",     "gray");  names(xcolors)<- xcategories
   xshapes     = c("square",  "square",  "circle",   "circle",    "circle",   "circle");names(xshapes)<- xcategories
   xsizes      = c(3*bNodeSz, 2*bNodeSz, bNodeSz,   3*bNodeSz,  2*bNodeSz,  bNodeSz);     names(xsizes) <- xcategories
   xfontsz     = c(3*bFontSz, 2*bFontSz, bFontSz,   3*bFontSz,  2*bFontSz,  bFontSz);     names(xfontsz)<- xcategories

   no.nodes = length(allnodes)

   # legend table 
   legendtb = cbind(xcategories, xshapes, xcolors, xcolors, xsizes, xfontsz)
   colnames(legendtb) <- c("label", "shape", "color", "border", "node_size", "font_size")

   sigInNet = intersect(allnodes, signature)
   sig_status = rep("NSIG", no.nodes); names(sig_status) <- allnodes; sig_status[sigInNet]="SIG"
   kdr_status = rep("NKD",  no.nodes); names(kdr_status) <- allnodes; 

   nf.cols = dim(kdaMatrix)[2]; nf.rows = dim(kdaMatrix)[1]
   keydrvNames = NULL
   if(nf.rows>0) {
     keydrv  = as.integer(kdaMatrix[,nf.cols])

     # global driver
     keysel  = c(1:nf.rows)[keydrv==1]; 
     keydrvNames = kdaMatrix[keysel,1];
     kdr_status[keydrvNames] = "GKD"

     # local driver
     if(sum(keydrv==0)>0) {
        keysel  = c(1:nf.rows)[keydrv==0];
        keydrvNames = kdaMatrix[keysel,1];
        kdr_status[keydrvNames] = "LKD"
     }

     # combined signature-keydriver status
     #
     sigkdr_status=paste(sig_status, kdr_status, sep="_")
     hnList = tapply(allnodes, sigkdr_status, list) # make a list for each category
     sigkdr_names = names(hnList)

     isNonSig = intersect(xcat2, sigkdr_names) # if all nodes are signatures, we use only circle for display
     if(length(isNonSig)==0){
       xshapes     = c("circle",   "circle",   "circle",  "circle",   "circle",   "circle");names(xshapes)<- xcategories
     }

     # set up actual visualization properties
     yHighColor = xcolors[sigkdr_names]
     yHighShape = xshapes[sigkdr_names]
     yHighSize  = xsizes[sigkdr_names]
     yHighFontSZ= xfontsz[sigkdr_names]

   } else {
     hnList     = list(sigInNet) # highlight only signature
     yHighColor = c("brown")
     yHighShape = c("circle")
     yHighSize  = c("1")
     yHighFontSZ= c("1")
   }

   return( list(hnList, cbind(yHighColor, yHighShape, yHighSize, yHighFontSZ), legendtb) )
}


# configure nodes's shapes and colors for cytoscape visualization
#
# allnodes: all nodes in the network for visusalization
# signature: signature genes
# keydriverMatrix: results from KDA (first column is geneID and last column is key-driver status  
#                   with 1 for global driver and 0 for local drivers)
#
# so we use different shapes for types of 
#
# genesets is a list of gene sets, each with specific graph shapes predefined 
#             and graph properties of GKD, LKA, Validation are the same for each set 
#
# zG: zoom-in for the global driver; zL: zoom-in for the local driver
#
configure_node_visualization_complex = function(genesets, shapes, kdaMatrices, validated, bNodeSz=40, bFontSz=12, zG=4,zL=2) {

   # VA--validated; NVA--Not Validated; GKD--Global KeyDriver; LKD--Local KeyDriver; NKD--Not KeyDriver
   #
   xcategories = c("GKD_VA",   "LKD_VA",  "NKD_VA", "GKD_NVA", "LKD_NVA", "NKD_NVA"); 
   xcolors     = c("green",    "green",   "gray",   "yellow",  "yellow",  "gray");    names(xcolors)<- xcategories
   xsizes      = c(zG*bNodeSz, zL*bNodeSz, bNodeSz,  zG*bNodeSz,  zL*bNodeSz,  bNodeSz);   names(xsizes) <- xcategories
   xfontsz     = c(zG*bFontSz, zL*bFontSz, bFontSz,  zG*bFontSz,  zL*bFontSz,  bFontSz);   names(xfontsz)<- xcategories
   #xshapes     = c("square",  "square",  "circle",   "circle",    "circle",   "circle");names(xshapes)<- xcategories

   no.categ = length(xcategories) 
   no.sets  = length(genesets)

   # legend table 
   legendtb   = NULL;
   hnList     = list(rep(NA,no.categ*no.sets));
   yHighColor = NULL;
   yHighShape = NULL;
   yHighSize  = NULL;
   yHighFontSZ= NULL;
 
   cnt = 0
   setnames=names(genesets)
   for(i in c(1:no.sets) ) {

      ilegendtb = cbind(paste(setnames[i],xcategories), rep(shapes[i], no.categ), xcolors, xcolors, xsizes, xfontsz)
      legendtb  = rbind(legendtb, ilegendtb)

      allnodes = genesets[[i]]
      no.nodes = length(allnodes)

      kdaMatrix = kdaMatrices[[i]]
      nf.cols   = dim(kdaMatrix)[2]; nf.rows = dim(kdaMatrix)[1]

      merged     = merge(cbind(genesets[[i]]), kdaMatrix[,c(1,nf.cols)], by.x=1,by.y=1, all.x=T)
      merged     = as.matrix(merged)
      merged[,2] = ifelse( is.na(merged[,2]), "-1", merged[,2])

      ivalid   = setElementInSet_Fast(merged[,1], validated)
      kdstatus = ifelse(merged[,2]=="1", "GKD", ifelse(merged[,2]=="0", "LKD","NKD") )
      vastatus = ifelse(ivalid, "VA", "NVA")

      istatus  = paste(kdstatus, vastatus, sep="_")
      ishape   = rep(shapes[i], no.nodes)

      ustatus  = union(istatus, NULL)
      for (js in ustatus) {
          cnt = cnt + 1
          jsel = istatus == js
          # set up actual visualization properties
          hnList[[cnt]]  = merged[jsel,1]
          yHighColor = c(yHighColor, xcolors[js])
          yHighShape = c(yHighShape, shapes[i])
          yHighSize  = c(yHighSize,  xsizes[js])
          yHighFontSZ= c(yHighFontSZ,xfontsz[js])
      }
   }

   colnames(legendtb) <- c("label", "shape", "color", "border", "node_size", "font_size")

   return( list(hnList[c(1:cnt)], cbind(yHighColor, yHighShape, yHighSize, yHighFontSZ), legendtb) )
}



configure_node_visualization_complex_flattable = function(genesets, shapes, kdaMatrices, validated, bNodeSz=40, bFontSz=12) {

   # VA--validated; NVA--Not Validated; GKD--Global KeyDriver; LKD--Local KeyDriver; NKD--Not KeyDriver
   #
   xcategories = c("GKD_VA", "LKD_VA", "NKD_VA", "GKD_NVA", "LKD_NVA", "NKD_NVA"); 
   xcolors     = c("green",  "green",  "gray",   "yellow",  "yellow",  "gray");    names(xcolors)<- xcategories
   xsizes      = c(2*bNodeSz, bNodeSz, bNodeSz,  2*bNodeSz,  bNodeSz,  bNodeSz);   names(xsizes) <- xcategories
   xfontsz     = c(2*bFontSz, bFontSz, bFontSz,  2*bFontSz,  bFontSz,  bFontSz);   names(xfontsz)<- xcategories
   #xshapes     = c("square",  "square",  "circle",   "circle",    "circle",   "circle");names(xshapes)<- xcategories

   no.categ = length(xcategories) 

   # legend table 
   legendtb = NULL;
   hnList   = NULL;
   yHighColor = NULL;
   yHighShape = NULL;
   yHighSize  = NULL;
   yHighFontSZ= NULL;
 
   no.sets = length(genesets)
   setnames=names(genesets)
   for(i in c(1:no.sets) ) {

      ilegendtb = cbind(paste(setnames[i],xcategories), rep(shapes[i], no.categ), xcolors, xcolors, xsizes, xfontsz)
      legendtb  = rbind(legendtb, ilegendtb)

      allnodes = genesets[[i]]
      no.nodes = length(allnodes)

      kdaMatrix = kdaMatrices[[i]]
      nf.cols = dim(kdaMatrix)[2]; nf.rows = dim(kdaMatrix)[1]

      merged = merge(cbind(genesets[[i]]), kdaMatrix[,c(1,nf.cols)], by.x=1,by.y=1, all.x=T)
      merged = as.matrix(merged)
      merged[,2] = ifelse( is.na(merged[,2]), "-1", merged[,2])

      ivalid = setElementInSet_Fast(merged[,1], validated)
      kdstatus = ifelse(merged[,2]=="1", "GKD", ifelse(merged[,2]=="0", "LKD","NKD") )
      vastatus = ifelse(ivalid, "VA", "NVA")

      istatus  = paste(kdstatus, vastatus, sep="_")
      ishape   = rep(shapes[i], no.nodes)

      # set up actual visualization properties
      hnList     = c(hnList, merged[,1])
      yHighColor = c(yHighColor, xcolors[istatus])
      yHighShape = c(yHighShape, ishape)
      yHighSize  = c(yHighSize,  xsizes[istatus])
      yHighFontSZ= c(yHighFontSZ,xfontsz[istatus])
   }

   colnames(legendtb) <- c("label", "shape", "color", "border", "node_size", "font_size")

   return( list(hnList, cbind(yHighColor, yHighShape, yHighSize, yHighFontSZ), legendtb) )
}




#
#**********************************  END  **********************************




# Here we assume that both input tables don't have the title row
#
# transform each pairof gene ids in each link in a pairlist file into
#  the coressponding gene symbols, if no match, we use the original IDs
#  so node link will be missed in the out put
#
mapPairlinkIDS2genesymbol = function(fpairlink, inpath="",fmapping, myoutdir="")
{

    # --------- 1. read in network (pair list) ---------------------------------
    #
    fullname = paste(inpath, fpairlink, sep="")
    keywords = getFileName(fpairlink)

    ofname = paste(myoutdir, keywords, "-gsymb.pair", sep="")

    aaMatrix <- read.delim(fullname, sep="\t", header=F)
    aaMatrix =  as.matrix(aaMatrix)
    dim(aaMatrix)

    # --------- 2. read in mapping table (MMT-id to Gene Symbol) ---------------
    #
    mapMatrix <- read.delim(fmapping, sep="\t", header=colheader)
    mapMatrix =  as.matrix(mapMatrix)

    # --------- 3. merge the two t?wice (left node and right node) --------------
    #[A1,A2] & [B1,B2] = [A1,A2,B2L]

    mergeleft = merge(aaMatrix, mapMatrix, by.x=1, by.y=1, all.x=1)
    mergeleft <- as.matrix(mergeleft)

    # make up NA, ehich mean nod gene sysmbol matched, so we use the MMT-id
    #
    mergeleft[,3] = ifelse(is.na(mergeleft[,3]), mergeleft[,1], mergeleft[,3])
    dim(mergeleft)


    #[A1,A2,B2] & [B1, B2] = [A2,A1,B2L, B2R]
    #
    mergeright = merge(mergeleft, as.matrix(mapMatrix), by.x=2, by.y=1, all.x=1)
    dim(mergeright)

    mergeright <- as.matrix(mergeright)

    # make up NA, ehich mean nod gene sysmbol matched, so we use the MMT-id
    mergeright[,4] = ifelse(is.na(mergeright[,4]), mergeright[,1], mergeright[,4])

    dim(mergeright)

    #---------------------------- output adj matrix ------------------------------------------
    #title row
    
    write.table(mergeright[,c(3,4)], ofname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)

    return (ofname)
}

# here network can have additional information for each link
#
mapNetworkIDS2genesymbol = function(fnetwork, inpath="",colheader=F, fmapping, myoutdir="")
{

    # --------- 1. read in network (pair list) ---------------------------------
    #
    fullname = paste(inpath, fnetwork, sep="")

    keywords = getFileName(fnetwork)
    extension= getFileExtension(fnetwork)

    ofname = paste(myoutdir, keywords, "-gsymbol.", extension, sep="")

    aaMatrix <- read.delim(fullname, sep="\t", header=colheader)
    aaMatrix =  as.matrix(aaMatrix)
    dim(aaMatrix)

    nocols = dim(aaMatrix)[2]

    # --------- 2. read in mapping table (MMT-id to Gene Symbol) ---------------
    #
    mapMatrix <- read.delim(fmapping, sep="\t", header=F)
    mapMatrix =  as.matrix(mapMatrix)

    # --------- 3. merge the two t?wice (left node and right node) --------------

    colnames( mapMatrix) <- c("dst", "dstGeneSymbol")
    mergeright = merge( mapMatrix, aaMatrix, by.x=1, by.y=2, all.y=1)
    mergeright <- as.matrix(mergeright)

    # make up NA, ehich mean nod gene sysmbol matched, so we use the MMT-id
    #
    mergeright[,2] = ifelse(is.na(mergeright[,2]), mergeright[,1], mergeright[,2])
    dim(mergeright)


    #[A1,A2,B2] & [B1, B2] = [A2,A1,B2L, B2R]
    #
    colnames( mapMatrix) <- c("src", "srcGeneSymbol")
    mergeleft = merge(mapMatrix, mergeright, by.x=1, by.y=3, all.y=1)
    dim(mergeleft)

    mergeleft <- as.matrix(mergeleft)
    # make up NA, ehich mean nod gene sysmbol matched, so we use the MMT-id
    mergeleft[,2] = ifelse(is.na(mergeleft[,2]), mergeleft[,1], mergeleft[,2])

    dim(mergeleft)

    #title row
    
    write.table(mergeleft, ofname, sep="\t",quote=FALSE, col.names=colheader, row.names=FALSE, append=F)

    return (ofname)
}



# perform union and overlap networks between two networks in form of 
#  pair lists and encode the union and overlap networks
#
merge2networksBYpair = function(fnetlist, directeds, keywords, outputDir="") 
{
    codepairs= cbind( c(-1, 1), c(-3,3) )

    flogname = paste(outputDir, keywords, "_log.xls", sep="")
    foutcomb = paste(outputDir, keywords, "_union.adj", sep="")
    foutshare= paste(outputDir, keywords, "_overlap.adj", sep="")
    fshareids= paste(outputDir, keywords, "_overlap-nodeIDs.txt", sep="")
    foutsharpajek= paste(outputDir, keywords, "_overlap4pajek.xls", sep="")

    totalfiles = length(fnetlist)

    nodesInNets=rep(0, totalfiles+2) #plus the union and overlap networks
    edgesInNets=rep(0, totalfiles+2)
    namesofNets=rep("", totalfiles+2)

    allnodenames = NULL
    indexrange  = c(0)
    for (each in fnetlist){
      
      aaMatrix <- read.delim(each, sep="\t", header=F)
      aaMatrix = as.matrix(aaMatrix)

      # consider both columns
      allnodenames = c(allnodenames, as.character(aaMatrix[,1]) )
      allnodenames = c(allnodenames, as.character(aaMatrix[,2]) )

      indexrange  = c(indexrange, length(allnodenames) )

    }

    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    #analyze how many genes overlapped in the two networks
    mergedGmatrix=NULL
    for (i in c(1:totalfiles) ) {
       sidx = indexrange[i] + 1
       eidx = indexrange[i+1]
       
       inametable   = table(allnodenames[c(sidx:eidx) ])
       iuniquenames = as.matrix( names(inametable) )
       if(is.null(mergedGmatrix)){
           mergedGmatrix = iuniquenames
       }else{
           mergedGmatrix = merge(mergedGmatrix, iuniquenames, by.x=1, by.y=1, all=F)
       }
    }


    # initialize matrix with the size equal to no.uniquenames
    # for the first network, edge is represented by 1 and non-edge by -1
    # for the 2nd network, edge is represented by 3 and non-edge by -3
    # so element in the ajacency matrix of the combined network:
    #   = 4: link exists in both networks 
    #   =-2: link exists in 1st network but not 2nd network
    #   = 2: link exists in 2nd network but not 1st network
    #   =-4: link doesn't exist in any of the two networks
    #
    cname2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )


    combinedMatrix = NULL
    for (i in c(1:totalfiles) ){
      each     = fnetlist[i]
      fname    = getFileName(each)
      outfname = paste(outputDir, fname, sep='')

      aaMatrix <- read.delim(each, sep="\t", header=F)
      dim(aaMatrix)

      iadj = makeAjacencyMatrix(inetmatrix=aaMatrix, coding=codepairs[,i], 
                                matrixsize=no.uniquenames, 
                                directed=directeds[i], 
                                name2idxMatrix=cname2idxMatrix)

      if (is.null(combinedMatrix)){
        combinedMatrix = iadj
      }else{
        combinedMatrix = combinedMatrix + iadj
      }

      # we output the summary information about the current network
      aaMatrix = as.matrix(aaMatrix)
      inodenames = c(as.character(aaMatrix[,2]), as.character(aaMatrix[,1])) 
      inametable = table(inodenames)
      iuniquenames     = names(inametable)
      ino.uniquenames  = length(iuniquenames)

      nodesInNets[i] = ino.uniquenames
      edgesInNets[i] = dim(aaMatrix)[1]
      namesofNets[i] = fname

      # output the node Ids for the current input networks
      #  ifids    = paste(outputDir, fname, "-nodeIDs.txt", sep='')
      ifids    = paste(fname, "-nodeIDs.txt", sep='')
      write.table(as.matrix(iuniquenames), ifids, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)

    }


    # statistics about the union network
    nodesInNets[totalfiles + 1] = no.uniquenames
    edgesInNets[totalfiles + 1] = sum( (combinedMatrix != -4) )
    namesofNets[totalfiles + 1] = "union network"

    # change the coding for simplicity
    # element in the final ajacency matrix of the combined network:
    #   = 4: ( 4)link exists in both networks 
    #   = 1: (-2) link exists in 1st network but not 2nd network
    #   = 2: ( 2) link exists in 2nd network but not 1st network
    #   = 0: (-4) link doesn't exist in any of the two networks
    #
    #need a lot of memory, so we put it in the later part
    #
    #combinedMatrix = ifelse(combinedMatrix==-4, 0, 
    #                        ifelse(combinedMatrix==-2, 1, combinedMatrix) )


    #---------------------------- output combined adj matrix ------------------------------------------

    #title row
    write.table(t(as.matrix(uniquenames)), foutcomb, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)

    # write into file in many times for saving memory
    no.rows=dim(combinedMatrix)[1]
    step=200
    intervals = seq(from=1, to=no.rows, by=step )

    # test if the last one equal the no.rows, otherwise, append no.rows to the end
    no.intrervals= length(intervals)
    if (intervals[no.intrervals] != no.rows){
      intervals = c(intervals, no.rows)
    }
    no.intrervals= length(intervals)
    intervals[no.intrervals] = intervals[no.intrervals] + 1 # a little trick here

    total_matched=0
    for (i in c(2:no.intrervals) ){
       sidx= intervals[i-1]
       eidx= intervals[i]-1
       #cat(sidx,"\n")
       #cat(eidx,"\n")
       irows= c(sidx:eidx)
       imatrix=combinedMatrix[irows, ]
       imatrix=ifelse(imatrix==-4, 0, 
                            ifelse(imatrix==-2, 1, imatrix) )

       write.table(imatrix, foutcomb, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

       #matched edges
       bmatrix= imatrix==4
       
       total_matched = total_matched + sum(bmatrix)
    }


    #---------------------------- output shared adj matrix ----------------------------------------
    bool.matrix = combinedMatrix==4

    rows.sum    = apply(bool.matrix, 1, sum)
    cols.sum    = apply(bool.matrix, 2, sum)
    shared.bool = (rows.sum>0) | (cols.sum>0)

    sharedNames    = uniquenames[shared.bool]
    no.sharednodes = length(sharedNames)
    sharedMatrix   = combinedMatrix[shared.bool, shared.bool]
    #sum(sharedMatrix==4)

    sharedMatrix = ifelse(sharedMatrix==4, 1, 0)
    #sum(sharedMatrix==1)

    #title row
    write.table(t(as.matrix(sharedNames)), foutshare, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(sharedMatrix, foutshare, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # node IDs
    write.table(as.matrix(sharedNames), fshareids, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)


    #-------------------------- shared matrix for Pajekificator ---------------------------------
    if(1==2){
    colors=rep(13, no.sharednodes)
    cc    =rep(0.5,no.sharednodes)
    pajekMatrix = cbind(sharedNames, cc, colors, as.matrix(sharedMatrix) )

    titlerow=c("pajek","cc","colorCode", sharedNames)

    #title row
    write.table(t(as.matrix(titlerow)), foutsharpajek, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    write.table(pajekMatrix, foutsharpajek, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

    # statistics about the shared network
    nodesInNets[totalfiles + 2] = no.sharednodes
    edgesInNets[totalfiles + 2] = sum( sharedMatrix )
    namesofNets[totalfiles + 2] = "overlap network"

    #---------------------------- output statistics --------------------------------------------
    # matched percentage
    #
    matchedpercents = signif(total_matched/edgesInNets, 4)
    matchedpercents = matchedpercents*100
    matcheds        = rep(total_matched, totalfiles+2)
    matchMatrix     = cbind(namesofNets, nodesInNets, edgesInNets, matcheds, matchedpercents)
    fmatrix = rbind( c("network", "# of nodes", "# of edges", 
                     "# of matched edges", "matched percentage (%)") ,
                     matchMatrix)

    write.table(fmatrix, flogname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)

    sharedgenes = dim(mergedGmatrix)[1]
    mystring = paste("\nnumber of common nodes in two networks =", sharedgenes)
    appendStringToFile(flogname, mystring)


    #mystring = paste("\nnumber of matched edges=", total_matched)
    #appendStringToFile(flogname, mystring)
    #cat("total matched edges are ", total_matched, "\n")
    fmatrix
}




############################# Causality test #######################################
#


#################################################################################
#
#                           Josh's algorithm
#
#
#Causality Test -- Intersection-Union Test; Likelihood Ratio Tests

#Example call: caustest(dat,"snp1","g1","g2")

#Input
#dat: a dataframe with locus data and expression data for gene1 and gene2
#The following 3 arguments are the colnames of the locus, gene1, and gene2.

#Output
# The output is two p-values, the first is for the test that g1 is causal for g2, 
#  and the second is for the test that g2 is causal for g1.

#Two p-values. The first is the p-value for the test:
#H0: independence, or g2 is causal for g1
#H1: gene1 is causal for gene2

#The second is the p-value for the test:
#H0: independence, or g1 is causal for g2
#H1: gene2 is causal for gene1


caustest = function(dat,loc,g1,g2){
	pvalues = rep(NA,2)
	dat1 = dat[,c(loc,g1,g2)]
	pvalues[1] = causp(dat1)
	dat1 = dat[,c(loc,g2,g1)]
	pvalues[2] = causp(dat1)
	pvalues
}


causp4 = function(dat1, no.bootstrap=50){
	#remove missing values
	for(j in 1:dim(dat1)[2]){
		dat1 = dat1[!is.na(dat1[,j]),]
	}

	names(dat1) = c("loc","r","c")
	llevels = as.integer(levels(as.factor(dat1$loc)))
	dfloc = length(llevels) - 1

	pvec = rep(NA,4)
	#fit0 = lm(c ~ 1)
	#fit1 = lm(c ~ as.factor(loc))
	#fit2 = lm(r ~ c)
	#fit3 = lm(c ~ r)
	#fit4 = lm(r ~ c + as.factor(loc))
	#fit5 = lm(c ~ r + as.factor(loc))

	fit0 = lm(c ~ 1,data=dat1)
	fit1 = lm(c ~ as.factor(loc),data=dat1)
	fit2 = lm(r ~ c,data=dat1)
	fit3 = lm(c ~ r,data=dat1)
	fit4 = lm(r ~ c + as.factor(loc),data=dat1)
	fit5 = lm(c ~ r + as.factor(loc),data=dat1)

	l = 2*(logLik(fit1)[1] - logLik(fit0)[1])
	pvec[1] = pchisq(l,df=dfloc,ncp=0,lower.tail=F)

	l = 2*(logLik(fit4)[1] - logLik(fit2)[1])
	pvec[2] = pchisq(l,df=dfloc,ncp=0,lower.tail=F)

	l = 2*(logLik(fit5)[1] - logLik(fit1)[1])
	pvec[3] = pchisq(l,df=1,ncp=0,lower.tail=F)

	l_ = 2*(logLik(fit5)[1] - logLik(fit3)[1])

	#randomize  within genotype to estimate ncp
	datalist = vector("list",(dfloc+1))
	for(i in 1:(dfloc + 1)){
		datalist[[i]] = dat1[dat1$loc == llevels[i],]
	}
	lvecr = rep(NA,no.bootstrap)
	rclist = vector("list",(dfloc+1))
	rdatlist = vector("list",(dfloc+1))

	for(i in 1:no.bootstrap){
		for(ii in 1:(dfloc + 1)){
			rclist[[ii]] = datalist[[ii]][order(runif(dim(datalist[[ii]])[1])),"c"]
			rdatlist[[ii]] = as.data.frame(cbind(datalist[[ii]][,c("r","loc")],rclist[[ii]]))
			names(rdatlist[[ii]])[3] = "c"
		}

		if(dfloc == 1) rdat = rbind(rdatlist[[1]],rdatlist[[2]],deparse.level=F)

		if(dfloc == 2) rdat = rbind(rdatlist[[1]],rdatlist[[2]],rdatlist[[3]],deparse.level=F)
		names(rdat) = c("r_","loc_","c_")

		fit0r = lm(c_ ~ 1,data=rdat)
		fit3r = lm(c_ ~ r_,data=rdat)
		fit5r = lm(c_ ~ r_ + as.factor(loc_),data=rdat)
		lvecr[i] = 2*(logLik(fit5r)[1] - logLik(fit3r)[1])
	}

	lvecr = lvecr[!is.na(lvecr)]
	el = mean(lvecr,na.rm=F)
	ncp = el - dfloc
	if(ncp < 0) ncp = 0
	pvec[4] = pchisq(l_,dfloc,ncp=ncp)
	pval = max(pvec)
        return(pval)
}


causp = function(dat1, no.bootstrap=50){
	#remove missing values
	dat_f = dat1
	for(j in 1:dim(dat_f)[2]){
		dat_f = dat_f[!is.na(dat_f[,j]),]
	}

   names(dat_f) = c("L","G","T")
   llevels = as.integer(levels(as.factor(dat_f$L)))
   dfL = length(llevels) - 1

   pvec = rep(NA,4)

   if(dfL == 2){
	dat_f$L1 = ifelse(dat_f$L == 1,1,0)
	dat_f$L2 = ifelse(dat_f$L == 2,1,0)

	fit0 = lm(T ~ 1,data=dat_f)
	fit1 = lm(T ~ L1 + L2,data=dat_f)
	fit2 = lm(G ~ T,data=dat_f)
	fit3 = lm(T ~ G,data=dat_f)
	fit4 = lm(G ~ T + L1 + L2,data=dat_f)
	fit5 = lm(T ~ G + L1 + L2,data=dat_f)

	pvec[1] = anova(fit0,fit1)$"Pr(>F)"[2]
	pvec[2] = anova(fit2,fit4)$"Pr(>F)"[2]
	pvec[3] = anova(fit1,fit5)$"Pr(>F)"[2]

	f_ = anova(fit3,fit5)$F[2]

	fit1G = lm(G ~ L1 + L2,data=dat_f)

	alg = summary(fit1G)$coefficients["(Intercept)",1]
	blg1 = summary(fit1G)$coefficients["L1",1]
	blg2 = summary(fit1G)$coefficients["L2",1]

	dat_f$eG = resid(fit1G)       

	ss = dim(dat_f)[1]
	lvecr = rep(NA,no.bootstrap)
	fvecr = rep(NA,no.bootstrap)

	for(i in 1:no.bootstrap){
		nni <- trunc(1 + ss*runif(ss, 0, 1)) ;
		dat_f$G_ = alg + blg1*dat_f$L1 + blg2*dat_f$L2 + dat_f$eG[nni]

		fit_0 = lm(T ~ G_,data=dat_f)
		fit_1 = lm(T ~ G_ + L1 + L2,data=dat_f)
		lvecr[i] = 2*(logLik(fit_1)[1] - logLik(fit_0)[1])
		fvecr[i] = anova(fit_0,fit_1)$F[2]
	}
   }#End dfL == 2
   if(dfL == 1){

	dat_f$L1 = ifelse(dat_f$L == 1,1,0)

	fit0 = lm(T ~ 1,data=dat_f)
	fit1 = lm(T ~ L1,data=dat_f)
	fit2 = lm(G ~ T,data=dat_f)
	fit3 = lm(T ~ G,data=dat_f)
	fit4 = lm(G ~ T + L1,data=dat_f)
	fit5 = lm(T ~ G + L1,data=dat_f)

	pvec[1] = anova(fit0,fit1)$"Pr(>F)"[2]
	pvec[2] = anova(fit2,fit4)$"Pr(>F)"[2]
	pvec[3] = anova(fit1,fit5)$"Pr(>F)"[2]

	f_ = anova(fit3,fit5)$F[2]

	fit1G = lm(G ~ L1,data=dat_f)

	alt = summary(fit1)$coefficients["(Intercept)",1]
	blt1 = summary(fit1)$coefficients["L1",1]

	alg = summary(fit1G)$coefficients["(Intercept)",1]
	blg1 = summary(fit1G)$coefficients["L1",1]

	dat_f$eG = resid(fit1G)       

	ss = dim(dat_f)[1]
	lvecr = rep(NA,no.bootstrap)
	fvecr = rep(NA,no.bootstrap)

	for(i in 1:no.bootstrap){
		nni <- trunc(1 + ss*runif(ss, 0, 1)) ;
		dat_f$G_ = alg + blg1*dat_f$L1 + dat_f$eG[nni]

		fit_0 = lm(T ~ G_,data=dat_f)
		fit_1 = lm(T ~ G_ + L1,data=dat_f)
		fvecr[i] = anova(fit_0,fit_1)$F[2]
	}
   } #End dfL == 1

   #####F Method
   fvecr = fvecr[!is.na(fvecr)]
   df1 = anova(fit3,fit5)$Df[2]
   df2 = anova(fit3,fit5)$Res.Df[2]
   fncp = mean(fvecr,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
   if(fncp < 0) fncp = 0

   ######### Transform F to normal
   npvals = pf(fvecr,df1,df2,ncp=fncp,lower.tail=TRUE)
   nfvecr = qnorm(npvals)

   npf = pf(f_,df1,df2,ncp=fncp,lower.tail=TRUE) #Transform observed F
   zf = qnorm(npf)
   pvec[4] = pnorm(zf,mean=0,sd=sd(nfvecr))
	pval = max(pvec)
	return(pval)              
}

# return (3,-2) means that too many missing values
#
causalityTest_Josh=function(LL, GG, TT, Bootstrap=F,no.bootstrap=50, pvaluecut=0.05){

  # here we have to preprocess the data to amke sure that no missing values there
  sel = (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  L=LL[sel]
  G=GG[sel]
  T=TT[sel]

  selno = sum(sel)
  # no enough samples
  if(selno < length(LL)/2){
    ret = c(3, -2)
    return(ret) 
  }

  # causality test
  # 
  datframe = data.frame( cbind(L, G, T) )
  pv1 <- causp(datframe, no.bootstrap=no.bootstrap);

  datframe2 = data.frame( cbind(L, T, G))
  pv2 <- causp(datframe2, no.bootstrap=no.bootstrap);


  if(is.na(pv1) ) {pv1=1}
  if(is.na(pv2) ) {pv2=1}

  result <- getRst(pv1, pv2, pvaluecut, pvaluecut);

  ret = c(signif(pv1,3), signif(pv2,3), result, 1)
  
  return (ret)
}




#################################################################################
#
#                           Eric's algorithm
#
#
# examine two models 1) G->T 2) T->G; the assumption here is that G and T are both
#   controlled by L, eg, 
# for each model, we use the following rationale to 
# for the 1st model, if lm(residual(G->T) ~ L) is significant (p2T), then we REJECT the 
#  model G->T
#
# pvalue for testing (T~G)~L, if the pval is significant, reject G->T
causalpval = function(L,H, G,T){
            # G --> T
            lm1T <- lm( T ~ G ) ;
            #lm1Ta <- anova(lm1T) ;

            #lm1Ts <- summary(lm1T) ; 

            resT <- resid(lm1T) ;

            lm2T <- lm( resT ~ L + H ) ;
            lm2Ts <- summary(lm2T) ; 
            #rsqrT <- lm2Ts$r.squared ; ntmp <- length(lm2Ts$residuals) ;
            if (!is.null(lm2Ts$fstatistic) ) {
               p2T <- 1 - pf(lm2Ts$fstatistic[1], lm2Ts$fstatistic[2], lm2Ts$fstatistic[3])
            }else {return (1)}

      return(p2T)   
}

getPv <- function(L, H, G, T) {
            p2T <- causalpval(L, H, G, T)
            p2G <- causalpval(L, H, T, G)

   return(c(p2G, p2T))
}

getPvOrg <- function(L, H, G, T) {
            # G --> T
            lm1T <- lm( T ~ G ) ;
            #lm1Ta <- anova(lm1T) ;
            resT <- resid(lm1T) ;
            lm2T <- lm( resT ~ L + H ) ;
            lm2Ts <- summary(lm2T) ; 
            #rsqrT <- lm2Ts$r.squared ; ntmp <- length(lm2Ts$residuals) ;
            p2T <- 1 - pf(lm2Ts$fstatistic[1], lm2Ts$fstatistic[2], lm2Ts$fstatistic[3])

            # T --> G
            lm1G <- lm( G ~ T ) ;
            #lm1Ga <- anova(lm1G) ;
            resG <- resid(lm1G) ;
            lm2G <- lm( resG ~ L + H ) ;
            lm2Gs <- summary(lm2G) ; 
            #rsqrG <- lm2Gs$r.squared ; #ntmp <- length(lm2Gs$residuals) ;
            p2G <- 1 - pf(lm2Gs$fstatistic[1], lm2Gs$fstatistic[2], lm2Gs$fstatistic[3])
   return(c(p2G, p2T))
}

getRst <- function(p2G, p2T, pvth1, pvth2) {
            rst <- 0
            if (p2T > pvth1 & p2G < pvth2) rst <- 1 ; #accept G->T & reject T->G
            if (p2T < pvth2 & p2G > pvth1) rst <- 2 ; #accept T->G & reject G->T
            if (p2T < pvth2 & p2G < pvth2) rst <- 3 ; #reject both G->T & T->G, G & T are independent, but both dependent on L
   return(rst)
}

getRstCmplx <- function(p2GandT, p2T, pvth1, pvth2, pfoldchange=10) {
            
            p2G=p2GandT[1]
            p2T=p2GandT[2]
            
            rst <- 0

            pGT = p2G/p2T
            pTG = p2T/p2G

            if (p2T > pvth1 & p2G < pvth2) rst <- 1 ; #accept G->T & reject T->G
            if (p2T < pvth2 & p2G > pvth1) rst <- 2 ; #accept T->G & reject G->T

            #if (p2T > pvth1 & p2G < pvth2 & pTG>=pfoldchange ) rst <- 1 ; #accept G->T & reject T->G
            #if (p2T < pvth2 & p2G > pvth1 & pGT>=pfoldchange ) rst <- 2 ; #accept T->G & reject G->T

            #reject both G->T & T->G, G & T are independent, but both dependent on L
            #if (p2T < pvth2/pfoldchange & p2G < pvth2/pfoldchange) { 
            if (p2T < pvth2 & p2G < pvth2) { 
                rst <- 3
            }

   return(rst)
}

getRstCmplx2 <- function(p2GandT, p2T, pvth1, pvth2, pfoldchange=100) {
            
            p2G=p2GandT[1]
            p2T=p2GandT[2]
            
            rst <- 0

            pGT = p2G/p2T
            pTG = p2T/p2G

            if (p2T > pvth1 & p2G < pvth2) rst <- 1 ; #accept G->T & reject T->G
            if (p2T < pvth2 & p2G > pvth1) rst <- 2 ; #accept T->G & reject G->T

            #if (p2T > pvth1 & p2G < pvth2 & pTG>=pfoldchange ) rst <- 1 ; #accept G->T & reject T->G
            #if (p2T < pvth2 & p2G > pvth1 & pGT>=pfoldchange ) rst <- 2 ; #accept T->G & reject G->T

            #reject both G->T & T->G, G & T are independent, but both dependent on L
            if (p2T < pvth2 & p2G < pvth2) { 
                rst <- 3
                if (pTG>=pfoldchange ) rst <- 1 ; #accept G->T & reject T->G
                if (pGT>=pfoldchange ) rst <- 2 ; #accept T->G & reject G->T
            }

   return(rst)
}




bootstrap <- function(L, H, G, T, callrst, id = 1, nkkk = 100) {
  # id = 1 to get reliability score for Eric's test
  #      2 for Eric's test + AIC
  #      3 for Eric's test + AIC + BIC
  # callrst is test rst with fields of eric, or aic or bic

    ntmp <- length(T) ; 
    sst1 <- NULL ; sst2 <- NULL ; sst3 <- NULL ;
    for (kkk in 1:nkkk) {
      nni <- trunc(1 + ntmp*runif(ntmp, 0, 1)) ; 
      tmp <- getPv(L[nni], H[nni], G[nni], T[nni]) ; 
      sst1 <- c(sst1, getRst(tmp[1], tmp[2], .05, .05)) 
     }

    relia <- NULL
    relia$eric <- sum(sst1 == callrst$eric)/nkkk
    if (id >= 2) relia$aic <- sum(sst2 == callrst$aic)/nkkk
    if (id >= 3) relia$bic <- sum(sst3 == callrst$bic)/nkkk
    return(relia)
}

bootstrap_multiLoci <- function(L, H, G, T, callrst, id = 1, nkkk = 100) {
  # id = 1 to get reliability score for Eric's test
  #      2 for Eric's test + AIC
  #      3 for Eric's test + AIC + BIC
  # callrst is test rst with fields of eric, or aic or bic

    ntmp <- length(T) ; 
    sst1 <- NULL ; sst2 <- NULL ; sst3 <- NULL ;
    for (kkk in 1:nkkk) {
      nni <- trunc(1 + ntmp*runif(ntmp, 0, 1)) ; 
      tmp <- getPv(L[nni, ], H[nni, ], G[nni], T[nni]) ; 
      sst1 <- c(sst1, getRst(tmp[1], tmp[2], .05, .05)) 
     }

    relia <- NULL
    relia$eric <- sum(sst1 == callrst$eric)/nkkk
    if (id >= 2) relia$aic <- sum(sst2 == callrst$aic)/nkkk
    if (id >= 3) relia$bic <- sum(sst3 == callrst$bic)/nkkk
    return(relia)
}


strlookup <- function(str1, str2) {
  q <- vector("numeric", length(str1))
  for (i in 1:length(str1)) {
    qq <- which(str1[i] == str2) ;
    if (length(qq) > 0) q[i] <- qq[1] 
  }
  return(q)
}

# return (3,-2) means that too many missing values
#
causalityTest=function(LL, HH, GG, TT, Bootstrap=F, no.bootstrap=10,  pvaluecut=0.05){

  # here we have to preprocess the data to amke sure that no missing values there
  sel = (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  L=LL[sel]
  H=HH[sel]
  G=GG[sel]
  T=TT[sel]

  selno = sum(sel)
  # no enough samples
  if(selno < length(LL)/2){
    ret = c(3, -2)
    return(ret) 
  }

  # causality test
  pv <- getPv(L, H, G, T) ; 
  callrst <- NULL;
  callrst$eric <- getRst(pv[1], pv[2], pvaluecut, pvaluecut);
 
  #bootstrap
  if (Bootstrap){
     bs <- bootstrap(L, H, G, T, callrst, 1, no.bootstrap);
     ret=c(signif(pv[1],3), signif(pv[2],3), callrst$eric, bs$eric)
  }else{
     ret = c(signif(pv[1],3), signif(pv[2],3), callrst$eric, -1)
  }
  
  return (ret)
}


########################################################################################
# 
################### causality test by multiple loci     ################################
#
#
# return (3,-2) means that too many missing values
# now LL and GG are matrix, columns as markers and rows as genes
#
if(F){
  i = 1  
  #  test causality btw all QTL genes
  #  get ID # for gene A & B and marker
    cidx = as.integer(neighborpairs[i,1])
    tidx = as.integer(neighborpairs[i,2])
    markx= as.integer(splitString(neighborpairs[i,3],";"))

   LL=t(genoMatrix[markx,]);  HH=t(hgenoMatrix[markx,]);
  GG=as.numeric(exprMatrix[cidx,]);    TT=as.numeric(exprMatrix[tidx,]);
  Bootstrap=T; no.bootstrap=100;
  pvaluecut=0.05
}

causalityTest_MultiLoci=function(LL, HH, GG, TT, Bootstrap=F, no.bootstrap=10,  pvaluecut=0.05){

  # here we have to preprocess the data to make sure that no missing values there
  naL   = is.na(LL)
  naLsum= apply(naL, 1, sum)

  naH   = is.na(HH)
  naHsum= apply(naH, 1, sum)

  sel = (!is.na(GG)) & (!is.na(TT)) & (naLsum==0) & (naHsum==0)
  
  L=LL[sel,]
  H=HH[sel,]
  G=GG[sel]
  T=TT[sel]

  selno = sum(sel)
  # no enough samples
  if(selno < length(GG)/2){
    ret = c(3, -2)
    return(ret) 
  }

  # causality test
  pv <- getPv(L, H, G, T) ; 
  callrst <- NULL;
  callrst$eric <- getRst(pv[1], pv[2], pvaluecut, pvaluecut);
 
  #bootstrap
  if (Bootstrap){
     bs <- bootstrap_multiLoci(L, H, G, T, callrst, 1, no.bootstrap);
     ret=c(signif(pv[1],3), signif(pv[2],3), callrst$eric, bs$eric)
  }else{
     ret = c(signif(pv[1],3), signif(pv[2],3), callrst$eric, -1)
  }
  
  return (ret)
}


icausalityTest_multiLoci=function(LL, HH, GG, TT, Bootstrap=F, no.bootstrap=100,  pvTtestCut=0.05, pvNocallCut=0.05){

  # here we have to preprocess the data to amke sure that no missing values there
  # here we have to preprocess the data to make sure that no missing values there
  naL   = is.na(LL)
  naLsum= apply(naL, 1, sum)

  naH   = is.na(HH)
  naHsum= apply(naH, 1, sum)

  sel = (!is.na(GG)) & (!is.na(TT)) & (naLsum==0) & (naHsum==0)

  L=LL[sel,]
  H=HH[sel,]
  G=GG[sel]
  T=TT[sel]

  selno = sum(sel)
  # no enough samples
  if(selno < length(GG)/2){
    ret = c(3, -2)
    return(ret) 
  }

  # causality test
  pv <- getPv(L, H, G, T) ; 
  #callrst <- NULL;
  #callrst$eric <- getRst(pv[1], pv[2], pvaluecut, pvaluecut);
  
  xx = xbootstrap_multiLoci(L, H, G, T, nkkk = no.bootstrap)
  xx = rbind(c(pv[1],pv[2]), xx)
  colnames(xx) <- c("p2G", "p2T")
  
  ttest   = t.test(xx[,1], xx[,2])
  tpval   = ttest$p.value
  medpval = c(median(xx[,1],na.rm=T), median(xx[,2],na.rm=T) )

  callrst <- NULL;
  callrst$eric <- igetRst(medpval[1], medpval[2], tpval, pvTtestCut, pvNocallCut);
  
  ret=c(signif(pv[1],3), signif(pv[2],3), callrst$eric, 1-tpval)
  
  return (ret)
}


########################################################################################
# 
################### causality test by interaction terms ################################
#
if(F){

  library(lme4)
 
  i = 1
  
  #  test causality btw all QTL genes
  #  get ID # for gene A & B and marker
  cidx = neighborpairs$GeneIdxA[i]
  tidx = neighborpairs$GeneIdxB[i]
  markx= neighborpairs$QTLmarkIdx[i]

  LL=as.numeric(genoMatrix[markx,]);   HH=as.numeric(hgenoMatrix[markx,]);
                                GG=as.numeric(exprMatrix[cidx,]);    TT=as.numeric(exprMatrix[tidx,]);
                                Bootstrap=T; no.bootstrap=100;
                                pvaluecut=0.05

  sel = (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  L=LL[sel]
  H=HH[sel]
  G=GG[sel]
  T=TT[sel]

  lmerT2G = lmer(G ~ T + L+T*L + (1|L), method="ML")
  lmerG2T = lmer(T ~ G + L+G*L+ (1|L), method="ML")
  anova(lmerT2G, lmerG2T)

  lmerT2G = lm(G ~ T + L+T*L )
  lmerG2T = lm(T ~ G + L+G*L)

  summary(lmerG2T)
  summary(lmerT2G)

  anova(lmerT2G, lmerG2T)
  anova(lmerG2T, lmerT2G)


  xx = xbootstrap(L, H, G, T, nkkk = 100)
  colnames(xx) <- c("T2G", "G2T")
  
  t.test(xx[,1], xx[,2])
  means = apply(xx, 2, median, na.rm=T)
  means

  yy = -log10(xx)

  pval = wilcoxtest4Hugedata(xx[,1],xx[,2])
  pval
  pval = wilcoxtest4Hugedata(yy[,1],yy[,2])
  pval

  means = apply(xx, 2, median, na.rm=T)
  means
  pval

  ss1 = hist(yy[,1])
  ss2 = hist(yy[,2])
  br = ss1$breaks + (ss1$breaks[2]-ss1$breaks[1])/2
  lines(ss1$breaks, ss1$density, color="blue")
  lines(ss2$breaks, ss2$density, col="red")

  LMp2T <- iAcausalpval(L, G, T)
  LMp2G <- iAcausalpval(L, T, G)

  p2T <- icausalpval(L, G, T)
  p2G <- icausalpval(L, T, G)
  p2T
  p2G

  lm1Ti <- lm( T ~ G + L + G*L) ;
  lm1Gi <- lm( G ~ T + L + T*L) ;
}


icausalityTest=function(LL, HH, GG, TT, Bootstrap=F, no.bootstrap=100,  pvTtestCut=0.05, pvNocallCut=0.05){

  # here we have to preprocess the data to amke sure that no missing values there
  sel = (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  L=LL[sel]
  H=HH[sel]
  G=GG[sel]
  T=TT[sel]

  selno = sum(sel)
  # no enough samples
  if(selno < length(LL)/2){
    ret = c(1,1,3, -2); #c(3, -2)
    return(ret) 
  }

  # causality test
  pv <- getPv(L, H, G, T) ; 
  #callrst <- NULL;
  #callrst$eric <- getRst(pv[1], pv[2], pvaluecut, pvaluecut);

  xx = xbootstrap(L, H, G, T, nkkk = no.bootstrap)
  xx = rbind(c(pv[1],pv[2]), xx)
  colnames(xx) <- c("p2G", "p2T")
  
  #ttest   = t.test(xx[,1], xx[,2])

  ttest   = tryCatch(t.test(xx[,1], xx[,2]), error=function(e) (return(NULL)) )

  if (!is.null(ttest) ) {
     tpval   = ttest$p.value
     if (is.na(ttest$p.value)) tpval =1;

  }else{tpval= 1}

  medpval = c(median(xx[,1],na.rm=T), median(xx[,2],na.rm=T) )
  
  callrst <- NULL;
  callrst$eric <- igetRst(medpval[1], medpval[2], tpval, pvTtestCut, pvNocallCut);
  
  #ret=c(signif(pv[1],3), signif(pv[2],3), callrst$eric, 1-tpval)
  ret=c(signif(medpval[1],3), signif(medpval[2],3), callrst$eric, 1-tpval)
  
  return (ret)
}

igetRst <- function(p2G, p2T, testpv, pvthrshold, pvnocall) {
            if (p2G>=pvnocall & p2T>= pvnocall){return(0);}

            rst <- 3 # by default, #reject both G->T & T->G, G & T are independent, but both dependent on L

            if (p2T > p2G & testpv<pvthrshold) rst <- 1 ; #accept G->T & reject T->G
            if (p2T < p2G & testpv<pvthrshold) rst <- 2 ; #accept T->G & reject G->T

            return(rst)
}



# pvalue for testing (T~G)~L, if the pval is significant, reject G->T
# Rationale: if T is better explained by G + G*L, ( T ~ G + G*L) ~ L will be less significant
#
icausalpvalA = function(L,G,T){
            # G --> T
            lm1T <- lm( T ~ G ) ;
            #lm1Ta <- anova(lm1T) ;
            resT <- resid(lm1T) ;
            lm2T <- lm( resT ~ L ) ;
            lm2Ts <- summary(lm2T) ; 
            #rsqrT <- lm2Ts$r.squared ; ntmp <- length(lm2Ts$residuals) ;
            p2T <- 1 - pf(lm2Ts$fstatistic[1], lm2Ts$fstatistic[2], lm2Ts$fstatistic[3])

            # G+L*G --> T
            lm1T <- lm( T ~ G + G*L) ;
            #lm1Ta <- anova(lm1T) ;
            resT <- resid(lm1T) ;
            lm2T <- lm( resT ~ L) ;
            lm2Ts <- summary(lm2T) ; 
            #rsqrT <- lm2Ts$r.squared ; ntmp <- length(lm2Ts$residuals) ;
            p2Ti <- 1 - pf(lm2Ts$fstatistic[1], lm2Ts$fstatistic[2], lm2Ts$fstatistic[3])

      return( min(p2T, p2Ti), max(p2T, p2Ti))   
}

# if p2T is sigificant, it means that T cannot be explained by G and it sinteraction with L
# otherwise, T can be explained by G and it sinteraction with L
#
icausalpval = function(L,G,T, pcut=0.05){

            lm1Ti <- lm( T ~ G) ;
            resT <- resid(lm1Ti) ;

            lm2T <- lm( resT ~ L + L*G) ;

            #lm2T <- lm( resT ~ L)
            lm2Ts <- summary(lm2T) ; 
            #rsqrT <- lm2Ts$r.squared ; ntmp <- length(lm2Ts$residuals) ;
            p2T <- 1 - pf(lm2Ts$fstatistic[1], lm2Ts$fstatistic[2], lm2Ts$fstatistic[3])

            return(p2T)
}


igetPv <- function(L, H, G, T) {
   p2T <- icausalpval(L, G, T)
   p2G <- icausalpval(L, T, G)

   return(c(p2G, p2T))
}


icausalpval_old = function(L,G,T, pcut=0.05){
            # G --> T
            lm1T <- lm( T ~ G ) ;
            #lm1Ta <- anova(lm1T) ;
            lm1Ti <- lm( T ~ G + G*L) ;

            resT <- resid(lm1T) ;

            # compare the two models
            avMC = anova(lm1T,lm1Ti)
            if (!is.na(avMC[2,6])){
               slm1Ti = summary(lm1Ti)
               # interactive model and the coefficient of the interactive term should be both significant
               if(avMC[2,6]<pcut & slm1Ti$coeff[3,4]<pcut){
                   resT <- resid(lm1Ti) ;
               }
            }
            
            lm2T <- lm( resT ~ L ) ;
            lm2Ts <- summary(lm2T) ; 
            #rsqrT <- lm2Ts$r.squared ; ntmp <- length(lm2Ts$residuals) ;
            p2T <- 1 - pf(lm2Ts$fstatistic[1], lm2Ts$fstatistic[2], lm2Ts$fstatistic[3])

            return(p2T)
}


ibootstrap <- function(L, H, G, T, callrst, id = 1, nkkk = 100) {
  # id = 1 to get reliability score for Eric's test
  #      2 for Eric's test + AIC
  #      3 for Eric's test + AIC + BIC
  # callrst is test rst with fields of eric, or aic or bic

    ntmp <- length(T) ; 
    sst1 <- NULL ; sst2 <- NULL ; sst3 <- NULL ;
    for (kkk in 1:nkkk) {
      nni <- trunc(1 + ntmp*runif(ntmp, 0, 1)) ; 
      tmp <- igetPv(L[nni], H[nni], G[nni], T[nni]) ; 
      sst1 <- c(sst1, getRst(tmp[1], tmp[2], .05, .05)) 
     }

    relia <- NULL
    relia$eric <- sum(sst1 == callrst$eric)/nkkk
    if (id >= 2) relia$aic <- sum(sst2 == callrst$aic)/nkkk
    if (id >= 3) relia$bic <- sum(sst3 == callrst$bic)/nkkk
    return(relia)
}


xbootstrap <- function(L, H, G, T, nkkk = 100) {
  # id = 1 to get reliability score for Eric's test
  #      2 for Eric's test + AIC
  #      3 for Eric's test + AIC + BIC
  # callrst is test rst with fields of eric, or aic or bic

  ntmp <- length(T) ; 
  sst1 <- NULL ;
  for (kkk in 1:nkkk) {
      nni <- trunc(1 + ntmp*runif(ntmp, 0, 1)) ; 
      tmp <- getPv(L[nni], H[nni], G[nni], T[nni]) ; 
      sst1 <- rbind(sst1, tmp)
  }

  return(sst1)
}


xbootstrap_multiLoci <- function(L, H, G, T, nkkk = 100) {
  # id = 1 to get reliability score for Eric's test
  #      2 for Eric's test + AIC
  #      3 for Eric's test + AIC + BIC
  # callrst is test rst with fields of eric, or aic or bic

  ntmp <- length(T) ; 
  sst1 <- NULL ;
  for (kkk in 1:nkkk) {
      nni <- trunc(1 + ntmp*runif(ntmp, 0, 1)) ; 
      tmp <- getPv(L[nni,], H[nni,], G[nni], T[nni]) ; 
      sst1 <- rbind(sst1, tmp)
  }

  return(sst1)
}



##########################################################################
# identify QTL for one gene at one chromosome
# lod score should be ranked from small position to large
# th1 = 1, th2 = 2 for default
# rst has 3 fields: LOD score, max position and index for each QTL
#
# th1 is used to define the bell shape boundary
# th2 is used to define the min QTL LOD
#
identifyQTLorg <- function(lod, pos, th1 = 1, th2 = 2) {

  q1 <- lod > th1 ;
  q2 <- lod <= th1 ;

  tmp <- c(0, q2[1:(length(q2)-1)]) ;
  tmp1 <- which(tmp + q1 == 2) ;
  tmp <- c(q2[2:length(q2)], 0) ;
  tmp2 <- which(tmp + q1 == 2) ;

  if (lod[1] > th1) tmp1 <- c(1, tmp1)
  if (lod[length(lod)] > th1) tmp2 <- c(tmp2, length(lod))

  rst <- NULL;
  rst$lod <- NULL ; rst$pos <- NULL ; rst$idx <- NULL ;
  rst$idx1 <- NULL ; rst$idx2 <- NULL ;
  if (length(tmp1) == 0) return(rst)
  for (i in 1:length(tmp1)) {
    tmp <- max(lod[tmp1[i]:tmp2[i]]) ;
    I   <- which.max(lod[tmp1[i]:tmp2[i]]) ;
    if (tmp > th2) {

        j <- tmp1[i] + I - 1;
        j1 <- max(c(1, j-1)) ;
        j2 <- min(c(j+1, length(pos)))

        if (sum(lod[j1:j2] > th1) < length(j1:j2)) next

        rst$lod  <- c(rst$lod, tmp) ;
        rst$pos  <- c(rst$pos, pos[j]) ;
        rst$idx  <- c(rst$idx, j) ;
        rst$idx1 <- c(rst$idx1, tmp1[i]) ;
        rst$idx2 <- c(rst$idx2, tmp2[i]) ;
    }
  }
  return(rst)
}

# ##############################################################################
# 
# returned results are a matrix of QTL info and the columns are in the order of
# 
# QTLlod QTLchrom  QTLpos    QTLmarkeridx   QTLmarkerLeftBoundIdx     QTLmarkerRightBoundIdx
#

identifyQTL <- function(lod, pos, chrom, lodth1 = 1, lodth2 = 2) {

  q1 <- lod > lodth1 ;
  q2 <- lod <= lodth1 ;

  tmp <- c(0, q2[1:(length(q2)-1)]) ;
  tmp1 <- which(tmp + q1 == 2) ;
  tmp <- c(q2[2:length(q2)], 0) ;
  tmp2 <- which(tmp + q1 == 2) ;

  if (lod[1] > lodth1) tmp1 <- c(1, tmp1)
  if (lod[length(lod)] > lodth1) tmp2 <- c(tmp2, length(lod))

  rst <- NULL;

  # return if no LOD bigger than the thresholds
  if (length(tmp1) == 0) return(rst)
  peaks = (lod > lodth2)
  if ( sum(peaks) == 0) return(rst)

  for (i in 1:length(tmp1)) {
    tmp <- max(lod[tmp1[i]:tmp2[i]]) ;
    I   <- which.max(lod[tmp1[i]:tmp2[i]]) ;
    if (tmp > lodth2) {

        j  <- tmp1[i] + I - 1;
        j1 <- max(c(1, j-1)) ;
        j2 <- min(c(j+1, length(pos)))

        if (sum(lod[j1:j2] > lodth1) < length(j1:j2)) next

        # notice that here cbind generate a matrix instead of a vector
        #
        #            lod  chrom  pos    idx   idx1     idx2
        ires = cbind(tmp, chrom, pos[j], j,  tmp1[i],  tmp2[i])
        rst  = rbind(rst, ires)
    }
  }

  return(rst)
}

############################################################################################
# 
# search for genome-wide QTLs of a given gene
#
# returned results are a matrix of QTL info and the columns are in the order of
# 
# QTLlod QTLchrom QTLpos QTLmarkIdx QTLmarkLeftIdx QTLmarkRightIdx geneChrom genePosBeg genePosEnd
#

identifyQTLGenome <- function(LOD, POS, CHROM, geneChrom=NA, genePosBeg=NA, genePosEnd=NA,Llod=1,Hlod=2,maxcisdist=20000) {
  allres <- NULL;
  
  # how many chromosomes here
  clevels = levels(factor(CHROM))
  chroms  = as.integer(clevels)

  #check LOD by chromosome
  #
  index = c( 1:length(LOD) )# marker index
  for (ec in chroms){
     sel     = (CHROM == ec)
     selnona = sel & (!is.na(sel) )

     marksIdx = index[selnona]

     # returned results are a matrix of QTL info and the columns are in the order of
     # 
     # QTLlod QTLchrom  QTLpos    QTLmarkeridx   QTLmarkerLeftBoundIdx     QTLmarkerRightBoundIdx

     ecRes = identifyQTL(lod=LOD[selnona], pos=POS[selnona], chrom=ec, lodth1 = Llod, lodth2 = Hlod)

     # get absolute Index of markers
     for (mi in c(4:6) ){
         ecRes[,mi] = marksIdx [ ecRes[,mi] ]
     }

     allres= rbind (allres, ecRes)
  }

  if (is.null(allres) ){
     return (NULL)
  }

  noQTLs = dim(allres)[1]

  extracols = cbind( rep(geneChrom,noQTLs), rep(genePosBeg,noQTLs),  rep(genePosEnd,noQTLs) )

  finalres = cbind(allres, extracols)

  return (finalres)
}

# compare two intervals to see if they significantly overlap each other
#
qtlOverlap = function(intervalA, intervalB, minoverlap=0.6)
{
  # 1) there are no overlap at all
  if ( (intervalA[2] <= intervalB[1]) | (intervalA[1]>=intervalB[2]) ){
     return (F)
  }

  # 2) one is inside the other
  #
  if ( (intervalA[1]>=intervalB[1]) & (intervalA[2]<=intervalB[2]) ){
     return (T)
  }
  if ( (intervalB[1]>=intervalA[1]) & (intervalB[2]<=intervalA[2]) ){
     return (T)
  }
  
  #3) there are some overlap, then look at the percentage of overlap
  oS = max(intervalA[1], intervalB[1])
  oE = min(intervalA[2], intervalB[2])
  od = oE-oS

  ovlpA = od/(intervalA[2]- intervalA[1])
  ovlpB = od/(intervalB[2]- intervalB[1])

  if ( (ovlpA>minoverlap) | (ovlpB>minoverlap) ){
      return (T)#(c(T, ovlpA, ovlpB ) )
  }else{
      return (F)#(c(F, ovlpA, ovlpB ) )
  }

}

#
#
############################# END of Causality Test #######################################

# plot SNA format
#
# so, the adj matrix in SNA can be either "the link colors" or binary values 
#  representing connection/disconnection, we need differentiate the two situations
#

plotNetwork = function(input, directed, fimg="tmp.png", disphubs=3, plotfigure=T, nodenamecolor="purple", labelpos=0){

    #filename = getFileName(input)
    #fimg     = paste(filename, ".png",sep="")

    linkcolors = c("white",  "gray", "blue", "red")# "Green", "Gray"

    # 1) get the number of nodes in the network 
    #
    verticesno <- read.delim(input, sep="\t", header=F, skip=0, nrows=1)
    splitted = splitString(as.character(as.matrix(verticesno)), separator=" ")
    nonodes = as.integer(splitted[2])
    nonodes

    # 2) get the node information
    #
    vertices <- read.delim(input, sep="\t", header=T, skip=1, nrows=nonodes)
    dim(vertices)
    colnames(vertices)

    vertexName  = as.character(as.matrix(vertices$nodename) )
    vertexColor = as.character(as.matrix(vertices$color) )
    vertexShape = as.integer(as.matrix(vertices$shape) )
    
    vertexColor = ifelse(vertexColor =="white", "gray", vertexColor)
    vertexColor = ifelse(vertexColor =="White", "gray", vertexColor)

    labColor    = rep(nodenamecolor, nonodes)

    # 3) get the link information
    #

    firstrow <- read.delim(input, sep="\t", header=F, skip=nonodes+2, nrows=1)
    firstrow <- as.character(as.matrix(firstrow))

    # so, the adj matrix shows index of the link colors instead of binary values for 
    #  connection/disconnection
    #
    if (firstrow[1]== "edge_color_levels"){

        edge_color_levels= firstrow[-1] #remove flag of "edge_color_levels"

        edges <- read.delim(input, sep="\t", header=T, skip=nonodes+3) #**
        dim(edges)

        enum = as.numeric(as.matrix(edges))
        edges = matrix(enum, nonodes, nonodes)
        
        edgeColor = matrix(edge_color_levels[enum+1], nonodes, nonodes)
        #edges = ifelse(edges ==0, 0, 1)

    } else{
        edges <- read.delim(input, sep="\t", header=T, skip=nonodes+2)
        dim(edges)
        
        edges = matrix(as.numeric(as.matrix(edges)), nonodes, nonodes)

        # consider python vector index starts from '0'
        edgeColor = matrix("white", nonodes, nonodes)
        for (i in c(2:length(linkcolors)) ){
          edgeColor = ifelse( edges==(i-1), linkcolors[i], edgeColor)
        }
    }

    dim(edgeColor)

    # find hubs
    edges = ifelse(edges>0, 1, edges)#force to be adjacency matrix
    inlinks   = apply(edges, 2, sum) #degree(edges,cmode="indegree")
    outlinks  = apply(edges, 1, sum) #degree(edges,cmode="outdegree")

    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }
    hubidx    = order(-totallinks)

    top2hubs  = NULL
    top2links = NULL
    actualhubs= min(disphubs, length(inlinks))
    for (j in c(1:actualhubs)){
        top2hubs  = c(top2hubs,  vertexName[hubidx[j]])
        top2links = c(top2links, totallinks[hubidx[j]])
    }


      #make title
      #concatenate hubs and no of links
      nohubs = length(top2hubs)
      hubs   = ""
      links  = ""
      for (k in c(1:nohubs) ){
          if (k==1){
              hubs = top2hubs[k]
              links= top2links[k]
          }else{
              hubs = paste(hubs,  ", ", top2hubs[k], sep="")
              links= paste(links, ", ", top2links[k], sep="")
          }
      }
      hubsNnets = c(hubs, links)

      if (!plotfigure) {
         rm(edges)
         collect_garbage()
         return (hubsNnets)
      }

      if (disphubs>0){
         ititle = paste("Top hubs: ",hubs,sep="")
      }else{
         ititle = ""
      }

    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 3600){
      imgsize = 3600
    }

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(imgsize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes-5)


    objcnst     = 0.08
    objinter    = (0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes-5)
    objscale    = objscale/(1+vercnst^0.5)

    openImgDev(fimg,iwidth = imgsize, iheight = imgsize, ipointsize = 12)
    if (disphubs>0){
      par(mfrow=c(1,1), mar=c(0,1,2,1) )
    }else{
      par(mfrow=c(1,1), mar=c(0,1,0,1) )
    }
    gplot(   edges,                  diag=T,
             gmode= gmod,
             pad=2,                label.pad=0.1,
             vertex.cex=4/vercexscale, vertex.sides=vertexShape,
             vertex.col=vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=edgeColor,
             edge.lty  = 0,          edge.lwd     = 1.3,
             displaylabels=TRUE,     label.bg="white",label.col=labColor,
             label=vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=1/labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/vercexscale,
             arrowhead.cex = 4/vercexscale,     
             object.scale  = objscale, main=ititle)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off()

    rm(edges)
    collect_garbage()

    return (hubsNnets)
}


plotMultiNetworks = function(inputs, directeds, fimg="tmp_multinets.png", disphubs=3, plotfigure=T, labelpos=0){

    if (length(inputs)==1){
       plotNetwork(inputs, directeds)
       return
    }

    nnets    = length(inputs)
    imgsize  = 1000
    nets     = as.list(rep(0, nnets) )

    hubsNnets = NULL

    for (i in c(1:nnets) ) {
      nets[[i]] = plotNetworkCore(inputs[i], directeds[i], idisphubs=disphubs, rimagesize = imgsize)

      #concatenate hubs and no of links
      nohubs = length(nets[[i]]$top2hubs)
      hubs   = ""
      links  = ""
      for (k in c(1:nohubs) ){
          if (k==1){
              hubs = nets[[i]]$top2hubs[k]
              links= nets[[i]]$top2links[k]
          }else{
              hubs = paste(hubs,  ", ", nets[[i]]$top2hubs[k], sep="")
              links= paste(links, ", ", nets[[i]]$top2links[k], sep="")
          }
      }
      hubsNnets = c(hubsNnets, hubs, links)
    }

    if (!plotfigure){
       rm(nets)
       #collect_garbage()
       return (hubsNnets)
    }

    openImgDev(fimg,iwidth = nnets*imgsize, iheight = imgsize, ipointsize = 12)

    if (disphubs>0){
      par(mfrow=c(1,1), mar=c(0,1,2,1) )
    }else{
      par(mfrow=c(1,1), mar=c(0,1,0,1) )
    }

    #plotNetworkCore(inputs[i], directeds[i], imagesize = imgsize)

    for (i in c(1:nnets) ) {
      #make title
      if (disphubs>0){
         ititle = paste("Top hubs: ",hubsNnets[1+(i-1)*2],sep="")
      }else{
         ititle = ""
      }

      gplot( nets[[i]]$edges,                  diag=T,
             gmode=nets[[i]]$gmod,  #mode="target",#"kamadakawai",
             pad=2,                label.pad=0.1,
             vertex.cex=4/nets[[i]]$vercexscale, vertex.sides=nets[[i]]$vertexShape,
             vertex.col=nets[[i]]$vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=nets[[i]]$edgeColor,
             edge.lty  = 0,          edge.lwd     = 1,
             displaylabels=TRUE,     label.bg="white",label.col=nets[[i]]$labColor,
             label=nets[[i]]$vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=0.8/nets[[i]]$labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/nets[[i]]$vercexscale,
             arrowhead.cex = 2/nets[[i]]$vercexscale,     
             object.scale  = nets[[i]]$objscale,main=ititle) 
      #text(0, 150, ititle, cex =nets[[i]]$labcexscale)
    }
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    dev.off()

    rm(nets)
    collect_garbage()

    return (hubsNnets)
}


plotNetworkCore = function(input, directed, idisphubs=3, rimagesize=1000){

    linkcolors = c("white",  "gray", "blue", "red")# "Green", "Gray"

    # 1) get the number of nodes in the network 
    #
    verticesno <- read.delim(input, sep="\t", header=F, skip=0, nrows=1)
    splitted = splitString(as.character(as.matrix(verticesno)), separator=" ")
    nonodes = as.integer(splitted[2])
    nonodes

    # 2) get the node information
    #
    vertices <- read.delim(input, sep="\t", header=T, skip=1, nrows=nonodes)
    dim(vertices)
    colnames(vertices)

    vertexName  = as.character(as.matrix(vertices$nodename) )
    vertexColor = as.character(as.matrix(vertices$color) )
    vertexShape = as.integer(as.matrix(vertices$shape) )

    vertexColor = ifelse(vertexColor =="white", "gray", vertexColor)
    vertexColor = ifelse(vertexColor =="White", "gray", vertexColor)

    labColor    = rep("purple", nonodes)

    # 3) get the link information
    #

    edges <- read.delim(input, sep="\t", header=T, skip=nonodes+2)
    dim(edges)

    edges = matrix(as.numeric(as.matrix(edges)), nonodes, nonodes)

    # consider python vector index starts from '0'
    edgeColor = matrix("white", nonodes, nonodes)
    for (i in c(2:length(linkcolors)) ){
      edgeColor = ifelse( edges==(i-1), linkcolors[i], edgeColor)
    }
    dim(edgeColor)

    # find hubs
    edges = ifelse(edges>0, 1, edges)#force to be adjacency matrix
    inlinks   = apply(edges, 2, sum) #degree(edges,cmode="indegree")
    outlinks  = apply(edges, 1, sum) #degree(edges,cmode="outdegree")
    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }
    hubidx    = order(-totallinks)

    top2hubs  = NULL
    top2links = NULL
    actualhubs= min(idisphubs, length(inlinks))
    for (j in c(1:actualhubs)){
        top2hubs  = c(top2hubs,  vertexName[hubidx[j]])
        top2links = c(top2links, totallinks[hubidx[j]])
    }


    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 3600){
      imgsize = 3600
    }

    pr = rimagesize/imgsize

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(rimagesize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes*pr-5)

    objcnst     = 0.08
    objinter    = pr*(0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes*pr-5)
    objscale    = objscale/(1+vercnst^0.5)

    netpara=list(edges=edges,
               gmod=gmod,
               vercexscale =vercexscale,
               vertexShape =vertexShape,
               vertexColor=vertexColor,
               edgeColor   =edgeColor,
               labColor    =labColor,
               vertexName  =vertexName,
               labcexscale =labcexscale,
               objscale    =objscale,
               top2hubs    =top2hubs,
               top2links   =top2links)
               
   return (netpara)
}

#------------------------ plot from a matrix with pairs -------------------------------------
#
# Here we provide the gene information matrix and specify the gene symbol column which to be
#  shown in the network
#
# Also output in/out/total links information as "*_inoutlinks.xls"
#

plotNetwork_inPairs_DisplayGeneSymbol = function(linkpairs, directed, geneinfo=NULL, genesymbolCol=2,
                   nodecolor="blue", nodeshape=30, nodenamecolor="purple", labelpos=1,
                   rankhubby="totallinks",
                   disphubs=3, plotfigure=T, fimg="tmp.png", 
                   saveDegrees=T, colorNodes=NULL, shapeNodes=NULL, shownodename=T, placemode = "fruchtermanreingold", default_title=""){

    #filename = getFileName(input)
    #fimg     = paste(filename, ".png",sep="")

    linkcolors = c("white",  "gray", "blue", "red")# "Green", "Gray"

    netrows <- dim(linkpairs)[1]
    netcols <- dim(linkpairs)[2]

    #------------------------------------------------------------------------------------
    #                        make network adjacency matrix
    #

    # find unique node names 
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(linkpairs[,1]) )
    allnodenames = c(allnodenames, as.character(linkpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    # find corrsponding gene symbols
    if (!is.null(geneinfo)) {
       unmatrix <- cbind(uniquenames)
       colnames(unmatrix) <- c("mmtid")

       ordergeneinfo  = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=unmatrix, 
                                         minorMatrix=geneinfo, missinglabel="", 
                                         keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
       ordergeneinfo = as.matrix(ordergeneinfo)

       # use MMTid if there is no symbol found
       ordergeneinfo[,genesymbolCol] = ifelse(is.na(ordergeneinfo[,genesymbolCol]),
                                              ordergeneinfo[,1],
                                              ordergeneinfo[,genesymbolCol])

       uniquesymbols = as.character(ordergeneinfo[,genesymbolCol])
    } else{
       ordergeneinfo = cbind(uniquenames)
       uniquesymbols = uniquenames
    }

    # look at the third column which is the confidence value 
    # if missing, add one more column 
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    if (netcols ==3) {
       mycoding =NULL
    } else{
       mycoding = c(0, 1)
    }    

    adjmatrix = makeAjacencyMatrix(inetmatrix=linkpairs, 
                                   coding=mycoding,
                                   matrixsize=no.uniquenames, 
                                   myname2idxMatrix= name2idxMatrix,
                                   directed=directed)
    
    # -------------- find hubs --------------------------
    #
    inlinks   = apply(adjmatrix, 2, sum) #degree(adjmatrix,cmode="indegree")
    outlinks  = apply(adjmatrix, 1, sum) #degree(adjmatrix,cmode="outdegree")

    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }

    # rank by different ways
    if (rankhubby=="inlinks"){
       hubidx    = order(-inlinks)
    }else if (rankhubby=="outlinks"){
       hubidx    = order(-outlinks)
    }else{
       hubidx    = order(-totallinks)
    }

    #-- make title using top hubs ---------------------
    #
    #concatenate hubs and no of links
    hubnames = uniquesymbols[hubidx]
    if(disphubs>0){
       tophubs  = concatenate(hubnames[c(1:disphubs)],", ")

      if (rankhubby=="outlinks"){
         ititle = paste("Top causal hubs: ", tophubs,sep="")
      } else{
        ititle = paste("Top hubs: ", tophubs,sep="")
      }
    }else{
        ititle = default_title
    }

    # ++++++++++  save the link matrix +++++++++++++++++++++++
    #
    linkmatrix = cbind(ordergeneinfo, inlinks, outlinks, totallinks)
    colnames(linkmatrix) <-c( colnames(ordergeneinfo), "inlinks", "outlinks", "total_links")
    linkmatrix = linkmatrix [hubidx,]
    if (!is.null(fimg) & saveDegrees ) {
       myfn = getFileName(fimg)
       flink= paste(myfn, "_inoutlinks.xls",sep="")
       write.table(linkmatrix, flink, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=F)
    }


    # 1) get the number of nodes in the network 
    #
    nonodes = no.uniquenames
    nonodes

    # 2) get the node information
    #
    if ( shownodename){
       vertexName  = uniquesymbols
    } else{
       vertexName  = rep("",nonodes)
    }

    if(length(nodecolor)>1){
         vertexColor = nodecolor
    }else{
         vertexColor         = rep(nodecolor, nonodes)
         if (is.null(colorNodes)) {
            if (disphubs>0){
               vertexColor[hubidx[c(1:disphubs)] ] = "red" #highlight hubs
            }#else{
             #  vertexColor[hubidx[1] ] = "red" #highlight hubs
             #}
          }else{
            colormerged = merge(name2idxMatrix, colorNodes, by.x=1,by.y=1,all=F)
            colormerged = as.matrix(colormerged)
            coloridx    = as.integer(colormerged[,2])
            vertexColor[coloridx] = "red"
          }
    }

    if(length(nodeshape)>1){
         vertexShape = nodeshape
    } else {
         vertexShape = rep(nodeshape, nonodes)
         if (!is.null(shapeNodes)) {
            shapemerged = merge(name2idxMatrix, shapeNodes, by.x=1,by.y=1,all=F)
            shapemerged = as.matrix(shapemerged)
            shapeidx    = as.integer(shapemerged[,2])
            vertexShape [shapeidx] = 3
         }
    }   

    if ( length(nodenamecolor)>1){
       labColor    = nodenamecolor
    }else{
       labColor    = rep(nodenamecolor, nonodes)
    }

    # 3) get the link information
    #
    edges = matrix(as.numeric(as.matrix(adjmatrix)), nonodes, nonodes)

    # consider python vector index starts from '0'
    if (no.uniquenames < 4000) {# to save memory
       edgeColor = matrix("white", nonodes, nonodes)
       for (i in c(2:length(linkcolors)) ){
         edgeColor = ifelse( edges==(i-1), linkcolors[i], edgeColor)
       }
       dim(edgeColor)
    } else{
       edgeColor = "gray"
    }

    rm(adjmatrix)
    collect_garbage()

    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 2000){ #3600
      imgsize = 2000
    }

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(imgsize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes-5)


    objcnst     = 0.08
    objinter    = (0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes-5)
    objscale    = objscale/(1+vercnst^0.5)


    if (!is.null(fimg)) {
       openImgDev(fimg,iwidth = imgsize, iheight = imgsize, ipointsize = 12)
    }
    par(mfrow=c(1,1), mar=c(0,1,2,1) )
    gplot(   edges,                  diag=T,
             gmode= gmod,
             pad=2,                label.pad=0.1,
             vertex.cex=3/vercexscale, vertex.sides=vertexShape,
             vertex.col=vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=edgeColor,
             edge.lty  = 0,          edge.lwd     = 1.3,
             displaylabels=TRUE,     label.bg="white",label.col=labColor,
             label=vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=1/labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/vercexscale,
             arrowhead.cex = 3/vercexscale,     
             object.scale  = objscale, main=ititle, mode = placemode)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    if (!is.null(fimg)) {
       dev.off()
    }

    rm(edges)
    collect_garbage()

    return (linkmatrix )
}


#********************************************************************************************
#------------------------ plot network in form of matrix -------------------------------------
#
# 1) different from network inpairs where we can only show nodes with at least one link
#    here we can show isolated nodes
#
# 2) geneinfo contains the columns in the order of: nodeIndex, node name, node size
#
# 3) here the link color is proportional to the value
#
# 3) the rest of output is the same as  plotNetwork_inPairs_DisplayGeneSymbol 
#

plotNetwork_inMatrix_DisplayGeneSymbol = function(adjmatrix, directed, geneinfo=NULL, genesymbolCol=2,
                   nodecolor="blue", nodeshape=30, nodenamecolor="purple", labelpos=1,
                   rankhubby="totallinks",
                   disphubs=3, plotfigure=T, fimg="tmp.png", 
                   saveDegrees=T, colorNodes=NULL, shownodename=T, scaleVertex=F,
                   placemode = "fruchtermanreingold", default_title=""){

    linkcolors = c("white",  "gray", "blue", "red")# "Green", "Gray"

    no.uniquenames= dim(adjmatrix)[1]

    imagetype = getFileExtension(fimg)    

    # scale up the size of vertex
    #
    nodesize      = rep(1, no.uniquenames)
    if (!is.null(geneinfo) ){
      if (scaleVertex & (dim(geneinfo)[2]==3) ) { 
         nodesize      = as.integer(geneinfo[,3])
         nodesize      = as.integer(5*nodesize/max(nodesize))
         nodesize      = ifelse(nodesize <1, 1, nodesize)
      }
      ordergeneinfo = geneinfo
      uniquesymbols = geneinfo[,2]
      uniquenames   = geneinfo[,2]
    } else{      
      uniquesymbols = colnames(adjmatrix)
      uniquenames   = colnames(adjmatrix)
      ordergeneinfo = cbind(uniquesymbols , uniquesymbols)
    }

    # -------------- find hubs --------------------------
    #
    inlinks   = apply(adjmatrix, 2, sum) #degree(adjmatrix,cmode="indegree")
    outlinks  = apply(adjmatrix, 1, sum) #degree(adjmatrix,cmode="outdegree")

    if (directed ){
       totallinks= inlinks + outlinks
    }else{
       totallinks= inlinks
    }

    # rank by different ways
    if (rankhubby=="inlinks"){
       hubidx    = order(-inlinks)
    }else if (rankhubby=="outlinks"){
       hubidx    = order(-outlinks)
    }else{
       hubidx    = order(-totallinks)
    }

    #-- make title using top hubs ---------------------
    #
    #concatenate hubs and no of links
    hubnames = uniquesymbols[hubidx]
    if(disphubs>0){
       tophubs  = concatenate(hubnames[c(1:disphubs)],", ")

      if (rankhubby=="outlinks"){
         ititle = paste("Top causal hubs: ", tophubs,sep="")
      } else{
        ititle = paste("Top hubs: ", tophubs,sep="")
      }
    }else{
        ititle = default_title
    }

    # ++++++++++  save the link matrix +++++++++++++++++++++++
    #
    linkmatrix = cbind(ordergeneinfo, inlinks, outlinks, totallinks)
    colnames(linkmatrix) <-c( colnames(ordergeneinfo), "inlinks", "outlinks", "total_links")
    linkmatrix = linkmatrix [hubidx,]
    if (!is.null(fimg) & saveDegrees ) {
       myfn = getFileName(fimg)
       flink= paste(myfn, "_inoutlinks.xls",sep="")
       write.table(linkmatrix, flink, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=F)
    }


    # 1) get the number of nodes in the network 
    #
    nonodes = no.uniquenames
    nonodes

    # 2) get the node information
    #
    if ( shownodename){
       vertexName  = uniquesymbols
    } else{
       vertexName  = rep("",nonodes)
    }

    if(length(nodecolor)>1){
         vertexColor = nodecolor
    }else{
         vertexColor         = rep(nodecolor, nonodes)
         if (is.null(colorNodes)) {
            if (disphubs>0){
               vertexColor[hubidx[c(1:disphubs)] ] = "red" #highlight hubs
            }#else{
             #  vertexColor[hubidx[1] ] = "red" #highlight hubs
             #}
          }else{
            colormerged = merge(name2idxMatrix, colorNodes, by.x=1,by.y=1,all=F)
            colormerged = as.matrix(colormerged)
            coloridx    = as.integer(colormerged[,2])
            vertexColor[coloridx] = "red"
          }
    }

    if(length(nodeshape)>1){
         vertexShape = nodeshape
    } else {
         vertexShape = rep(nodeshape, nonodes)
    }   

    if ( length(nodenamecolor)>1){
       labColor    = nodenamecolor
    }else{
       labColor    = rep(nodenamecolor, nonodes)
    }

    # 3) get the link information
    #
    edges = matrix(as.numeric(as.matrix(adjmatrix)), nonodes, nonodes)

    # consider python vector index starts from '0'
    if (no.uniquenames < 4000) {# to save memory
       edgeColor = matrix("white", nonodes, nonodes)
       for (i in c(2:length(linkcolors)) ){
         edgeColor = ifelse( edges==(i-1), linkcolors[i], edgeColor)
       }
       dim(edgeColor)
    } else{
       edgeColor = "gray"
    }

    rm(adjmatrix)
    collect_garbage()

    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = imgcnst + szintercept * (nonodes-5)
    if(imgsize > 2000){ #3600
      imgsize = 2000
    }

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(imgsize + 3600)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes-5)

    objcnst     = 0.08
    objinter    = (0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes-5)
    objscale    = objscale/(1+vercnst^0.5)

    if (!is.null(fimg)) {
       openImgDev(fimg,iwidth = imgsize, iheight = imgsize, ipointsize = 12)
    }

    if(imagetype =="pdf"){
        labcexscale = labcexscale*2
    }

    par(mfrow=c(1,1), mar=c(0,1,2,1) )
    gplot(   edges,                  diag=T,
             gmode= gmod,
             pad=2,                label.pad=0.1,
             vertex.cex=nodesize*3/vercexscale, vertex.sides=vertexShape,
             vertex.col=vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=edgeColor,
             edge.lty  = 0,          edge.lwd     = 1.3,
             displaylabels=TRUE,     label.bg="white",label.col=labColor,
             label=vertexName,       label.pos=labelpos,     boxed.labels=F,
             label.cex=1/labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = 1/vercexscale,
             arrowhead.cex = 3/vercexscale,     
             object.scale  = objscale, main=ititle, mode = placemode)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    if (!is.null(fimg)) {
       dev.off()
    }

    rm(edges)
    collect_garbage()

    return (linkmatrix )
}

#---------------------- plot SNP ----------------------------------------
#displayfirstfield="\\." display the first filed of the node names
#nodecolor_nolable: don't display the nodes in some colors to enhance the readability
#
# nodenamecolor = NA
#
plotNetworkSNP = function(input, directed="directed", fimg="tmp.png", disphubs=3, 
                          nodenamecolor="purple", plotfigure=T, labelpos=0, 
                          figuretitle="",
                          vertexscaleCoeff=1, displayfirstfield=NA, 
                          legend_gsize=0.08, legend_tcex=1.2, legend_vertical=T,
                          nodecolor_nolable=NULL, plot2D=TRUE, 
                          zoom_view=1, xyz_scales=c(1.2,2,1.2),
                          theta = 0, phi =0, fov = 60, lighteffect=F){
    # 0) get extension
    ext = getFileExtension(fimg)

    legendtable=NA; hasLegend=F

    # 1) get the number of nodes in the network 
    #
    verticesno <- read.delim(input, sep="\t", header=F, skip=0, nrows=1)
    splitted = splitString(as.character(as.matrix(verticesno)), separator=" ")

    if(splitted[1]=="Legend"){
       #     label shape       color border
       #[1,] "PFC" "circle"    NA    "gray"
       #
       no.lengends = as.integer(splitted[2])
       legendtable=read.delim(input, sep="\t", header=T, skip=1, nrows=no.lengends)
       legendtable=as.matrix(legendtable)
       skiplegend =no.lengends + 2
       hasLegend=T

       verticesno <- read.delim(input, sep="\t", header=F, skip=skiplegend, nrows=1)
       splitted = splitString(as.character(as.matrix(verticesno)), separator=" ")
       nonodes = as.integer(splitted[2])

    } else{
       skiplegend = 0
       nonodes = as.integer(splitted[2])
       nonodes
    }

    # 2) get the node information
    #
    vertices <- read.delim(input, sep="\t", header=T, skip=1+skiplegend, nrows=nonodes)
    dim(vertices)
    colnames(vertices)

    vertexName  = as.character(as.matrix(vertices$nodename) )
    vertexColor = as.character(as.matrix(vertices$color) )
    #vertexShape = as.integer(as.matrix(vertices$shape) )
    vertexShape = as.character(as.matrix(vertices$shape) )
    
    vertexColor = ifelse(vertexColor =="white", "gray", vertexColor)
    vertexColor = ifelse(vertexColor =="White", "gray", vertexColor)

    labColor    = rep(nodenamecolor, nonodes)

    # 3) get the link information
    #
    firstrow <- read.delim(input, sep="\t", header=F, skip=nonodes+2+skiplegend, nrows=1)
    firstrow <- as.character(as.matrix(firstrow))

    # so, the adj matrix shows index of the link colors instead of binary values for 
    #  connection/disconnection
    #
    edge_color_levels= firstrow[-1] #remove flag of "edge_color_levels"

    linksIdx <- read.delim(input, sep="\t", header=T, skip=nonodes+3+skiplegend) #**
    dim(linksIdx)
    li.cols = dim(linksIdx)[2]
    li.rows = dim(linksIdx)[1]
    linksIdx= matrix(as.integer(as.matrix(linksIdx )), li.rows, li.cols)

    # assign 1 for each link
    #
    edges = matrix(0, nonodes, nonodes)
    edges[ linksIdx[,c(1:2)] ] = 1

    edgeColor = matrix("white", nonodes, nonodes)
    #edges = ifelse(edges ==0, 0, 1)
    edgeColor[ linksIdx[,c(1:2)] ] =  edge_color_levels[ linksIdx[,3] ]

    #dim(edgeColor)

    # find hubs
    #edges = ifelse(edges>0, 1, edges)#force to be adjacency matrix

    inlinks   = apply(edges, 2, sum) #degree(edges,cmode="indegree")
    outlinks  = apply(edges, 1, sum) #degree(edges,cmode="outdegree")
    if (directed=="directed" ){
       totallinks= outlinks
    }else{
       totallinks= inlinks + outlinks
    }
    hubidx    = order(-totallinks)

    top2hubs  = NULL
    top2links = NULL
    actualhubs= min(disphubs, length(inlinks))
    for (j in c(1:actualhubs)){
        top2hubs  = c(top2hubs,  vertexName[hubidx[j]])
        top2links = c(top2links, totallinks[hubidx[j]])
    }

      #make title
      #concatenate hubs and no of links
      nohubs = length(top2hubs)
      hubs   = ""
      links  = ""
      for (k in c(1:nohubs) ){
          if (k==1){
              hubs = top2hubs[k]
              links= top2links[k]
          }else{
              hubs = paste(hubs,  ", ", top2hubs[k], sep="")
              links= paste(links, ", ", top2links[k], sep="")
          }
      }
      hubsNnets = c(hubs, links)

      if (!plotfigure) {
         rm(edges)
         collect_garbage()
         return (hubsNnets)
      }

      if (disphubs>0){
         ititle = paste("Top hubs: ",hubs,sep="")
      }else{         
         ititle = figuretitle
      }

    # 4) define image size & scale factor CEX

    #vercexscale = 15*nonodes/760
    #labcexscale = 1.6

    imgcnst     = 400
    szintercept = (3000-imgcnst)/(1000-5)
    imgsize     = as.integer((imgcnst + szintercept * (nonodes-5))/10)*10
    if(imgsize > 3600){
      imgsize = 3600
    }
    #imgsize = 400+as.integer(1500*sigmoid(nonodes, 0.01, 500)/10)*10
    imgsize = 200+2000*sigmoid(nonodes,0.0015, 1000)

    gmod = ifelse(directed, "digraph", "graph")

    # play with the CEX to get right scale
    labcexscale = 2*3600/(imgsize + 3600)
    labcexscale = 0.5*sigmoid(nonodes, -0.01, 500)
    labcexscale = 0.8*sigmoid(nonodes, -0.01, 500) + 2.4*sigmoid(nonodes, -0.01, 1200)*(nonodes/1200)

    vercnst      = 5
    verintercept = (20-vercnst)/(1000-5)
    vercexscale  = vercnst + verintercept*(nonodes^0.5-5)
    vercexscale  = 2+4.5*sigmoid(nonodes, -0.01, 500)
    verarrowsize = 0.25+0.5*sigmoid(nonodes, -0.01, 500)
    verlabeldist = 0.2*sigmoid(nonodes, -0.01, 500) + 0.1*sigmoid(nonodes, -0.01, 1200)*(nonodes/1200)

    vercexscale  =  6*sigmoid(nonodes, -0.002, 800)
    verarrowsize =  1*sigmoid(nonodes, -0.002, 800)
    labcexscale  =  1.1*sigmoid(nonodes, -0.002, 1200)
    verlabeldist = 0.25*sigmoid(nonodes,-0.002, 800)

    #x=c(30:2000)
    #plot(x, sigmoid(x, -0.002, 800))

    
    #delta = 0.5*sigmoid(nonodes, 0.005, 1400)  # original position may change by size
    delta=0
    if(nonodes>1000) {delta=0.05}
    if(ext=="pdf"){
       vercexscale  =   4*sigmoid(nonodes, -0.002, 800)
       verarrowsize = 0.4*sigmoid(nonodes, -0.002, 800)
       labcexscale  = 0.5*sigmoid(nonodes, -0.002, 1000)
       verlabeldist = 0.18*sigmoid(nonodes,-0.002, 800)

       delta = 0
    }

    objcnst     = 0.08
    objinter    = (0.02-objcnst)/(1000-5)
    objscale    = objcnst + objinter*(nonodes-5)
    objscale    = objscale/(1+vercnst^0.5)
    objscale    = 0.01+0.01*sigmoid(nonodes, -0.01, 500)

    edgewid = 1.5
    if(nonodes <200){
      edgewid = 2
    } else if (nonodes >= 3000){
      edgewid = 0.8
    }
    edgewid = 0.02+2*sigmoid(nonodes, -0.01, 500)

    # do not show the names
    vertexNamex = vertexName
    if (nodenamecolor=="white" | is.na(nodenamecolor)){
       vertexNamex = rep("", length(vertexName))
    }

    if (is.na(nodenamecolor)){
        nodenamecolorx = vertexColor
        verlabeldist = 0
    } else {
        nodenamecolorx = nodenamecolor
    }

    colnames(edges) <- vertexName
    ig = graph.adjacency(adjmatrix=edges, mode=directed, weighted=NULL, diag=TRUE, #c("directed", "undirected", "max","min", "upper", "lower", "plus")
        add.colnames=NULL, add.rownames=NA)

    xvertexName = vertexName 
    if(!is.na(displayfirstfield)){
        xvertexName = getSecondPart(vertexName, displayfirstfield, 1)
    }
    if ( !is.null(nodecolor_nolable) ) {
       for(ec in nodecolor_nolable){
          xvertexName = ifelse(vertexColor==ec, "", xvertexName)
       }
       labcexscale  = 1+sigmoid(nonodes, 0.002, 800)
       if(ext=="pdf"){labcexscale  = 0.8+1.1*sigmoid(nonodes, 0.002, 800)}
    }

    #----------------------------------------------------------------------------------------
    if(plot2D) {
      openImgDev(fimg,iwidth = imgsize, iheight = imgsize, ipointsize = 12*vertexscaleCoeff)
      if (ititle !=""){
        par(mfrow=c(1,1), mar=c(0,1,2,1) )
      }else{
        par(mfrow=c(1,1), mar=c(0,1,0,1) )
      }

      plot(ig, vertex.label=xvertexName, vertex.shape=vertexShape, 
               vertex.size=vercexscale, vertex.size2=vercexscale*0.6,
               vertex.frame.color=vertexColor, vertex.color= vertexColor,
               vertex.label.cex=labcexscale*vertexscaleCoeff,   vertex.label.dist=verlabeldist, 
               vertex.label.color=nodenamecolorx,     vertex.label.degree=-pi/2,
               edge.arrow.size=verarrowsize, layout=layout.fruchterman.reingold)

        if( hasLegend) {
          plotLegends(legends, xpos=-1.05+delta, ypos=-1.05+delta, sizex=2/6, graphsize=legend_gsize, 
         #plotLegends(legends, xpos=-1.12+delta, ypos=-1.14+delta, sizex=2/6, graphsize=legend_gsize, 
                   interval=0.1, tcex=legend_tcex+delta*3, vertical=legend_vertical)
        }

      par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
      dev.off()

    } else {

      scales=c(1,1,1); ys=1 # scale up Y axis

      rgl.open();
      par3d(windowRect=c(0,0,imgsize,imgsize), cex=labcexscale, no.readonly=F, scale=scales )

      #limpia la ventana antes de dibujar el gr?fico
      rgl.clear(type=c("shapes", "lights"))
 
      rgl.bg(color='white')
      rgl.viewpoint(theta = theta, phi =phi, fov = fov, interactive = F, zoom=1/ys)

      # set up the lighting source
      if(lighteffect){
         lid = rgl.light( theta = 0, phi =0, viewpoint.rel = T, ambient = "#FFFFFF", diffuse = "#FFFFFF", specular = "#FFFFFF")
      } else{
         lid = rgl.light( theta = 0, phi =0, viewpoint.rel = F, ambient = "#111111", diffuse = "#FFFFFF", specular = "#111111")
      }

      rglplot(ig, vertex.label=xvertexName, vertex.shape=vertexShape, 
               vertex.size=vercexscale, vertex.size2=vercexscale*0.6,
               vertex.frame.color=vertexColor, vertex.color= vertexColor,
               vertex.label.cex=labcexscale*vertexscaleCoeff,   vertex.label.dist=verlabeldist, 
               vertex.label.color=nodenamecolorx,     vertex.label.degree=-pi/2,
               edge.arrow.size=verarrowsize, layout=layout.fruchterman.reingold)

     ext = getFileExtension(fimg)
     if(ext=="png") {
       rgl.snapshot(fimg)
     } else {
       rgl.postscript(fimg,"pdf",drawText=T)
     }

     rgl.pop(type="lights", lid)
     rgl.close()

  }


    if(F){
    gplot(   edges,                  diag=T,
             gmode= gmod,
             pad=2,                label.pad=0.1,
             #vertex.cex=vertexscaleCoeff/vercexscale, vertex.sides=vertexShape,
             vertex.cex=vercexscale, vertex.sides=vertexShape,
             vertex.col=vertexColor,   vertex.border= 0,
             vertex.rot=0,             vertex.lty   = 1,
             edge.col=edgeColor,
             edge.lty  = 0,          edge.lwd     = edgewid,
             displaylabels=T,     label.bg="white",label.col=labColor,
             label=vertexNamex,       label.pos=labelpos,     boxed.labels=F,
             label.cex=labcexscale,label.lty=0,     label.lwd = 0, 
             loop.cex      = vercexscale,
             arrowhead.cex = 0.8*vercexscale,     
             object.scale  = objscale, main=ititle)
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off()

    }

    rm(edges)
    collect_garbage()

    return (hubsNnets)
}

#1. legends
#
#    label shape       color  border
#    "PFC" "circle"    "gray" "gray"
#    "CR"  "square"    "gray" "gray"
#    "VC"  "rectangle" "gray" "gray"
#
#2. sizex: size of each legend including graph & text
#3. graphsize: graph size only
#4. interval: interval between legends in Y axis
#
plotLegends = function(legends, xpos=-1, ypos=-1, sizex=0.3, 
                       graphsize=0.08, interval=1/10, 
                       tcex=1.2, vertical=T)
{

    no.lengends = dim(legends)[1]

    startY = ypos; startYt = ypos+0.04;

    gspace=sizex/2;   tspace=sizex/2; 
    lgdsz=graphsize; lgdszSQR=lgdsz*0.8;

    if(vertical) {
       curX=xpos;
    } else {
       curX=-sizex*no.lengends/2; # horizontal legend
    }

    for(ix in c(1:no.lengends) ) {
          if (legends[ix, 2] =="circle"){
             circle(center=c(curX+lgdsz/2,startY+lgdsz/2), radius=lgdsz/2, lwd=2, lty=1, col=legends[ix,3], border=legends[ix,4]); 
             text(x=curX+2.3*lgdsz, y=startYt, labels=legends[ix,1], cex=tcex, col="black")

          } else if (legends[ix, 2] =="square"){
             rect(xleft=curX, ybottom=startY, xright=curX+lgdszSQR,ytop=startY+lgdszSQR, col=legends[ix,3], border=legends[ix,4], lwd=2, lty=1)
             text(x=curX+2*lgdsz, y=startYt, labels=legends[ix,1], cex=tcex, col="black")

          } else if (legends[ix, 2] =="rectangle"){
             rect(xleft=curX, ybottom=startY+lgdsz/5, xright=curX+lgdsz, 
                  ytop=startY+lgdsz/5+lgdsz/2, col=legends[ix,3], border=legends[ix,4], lwd=2, lty=1)
             text(x=curX+2*lgdsz, y=startYt, labels=legends[ix,1], cex=tcex, col="black")
          }

       if(vertical) {
           startY= startY + interval; startYt = startYt + interval
       } else {
           curX = curX + sizex
       }
   }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end of plot networks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### read in Roy's combined expression, trait and genotype data #########
#
ReadBinary<-function(bin.fn)
{
    #
    # read previously created bin file; process; save it as text
    #
    now<-proc.time()[3]
    bin.names<-load(bin.fn)
    t.readTime<-proc.time()[3]-now
    #
    # verify data structure name
    #
    tName<-c("all.data")
    if (bin.names[1]!=tName)
    {
        stop(paste("bin data struct ", tName , " not found.\n\t(",
                    bin.names[1],")", sep=""))
    }
    return(all.data)
}

# get Batch ID from cross, tissue, and sex
#
getBatchId = function(sqlserver, dbname, icross, itissue, isex, flag="CorrRegress", shortenBatchName=F){

  # make a string as a combination of icross, itissue, isex, flag
  ilabel = paste( icross, " ", itissue, " ", isex, " ", flag, sep="")  
  ilabel = tolower(ilabel)

  TBbatch   = "Batch"
  econdition= NULL
  selectedfields = "Batch_name, Batch_id"
  batchTable <- accessWholeTableFromDB(sqlserver, dbname, TBbatch,
                 condition=econdition, outfields = selectedfields)

  # no such a table
  if( is.null(batchTable )) {
     return (-1)
  }

  batchTable <- as.matrix(batchTable)
  no.batches <- dim(batchTable)[1]

  # no.words in label
  lwords = splitString(ilabel,sep=" ")
  no.words= length(lwords)

  # to lower
  for (i in c(1:no.batches) ) {
     #batchTable[i,1] = tolower(batchTable[i,1])

     # here we use only the first "no.words" in the batch name
     # since for mci_bxa, it is hard to specify the long batch name
     #
     # for instance, ilabel="mci_bxa Liver All W10 CorrRegress"
     #
     # mci_bxa Liver  Male W10 CorrRegress Express Adjusted
     # mci_bxa Liver All W10 CorrRegress Adjusted Data
     if (shortenBatchName) {
        iwords = splitString(batchTable[i,1],sep=" ")
        no.min = min(no.words, length(iwords))
        iNwords= concatenate(iwords[1:no.min], " ")
        itmp = tolower(iNwords)
     } else{
        itmp = tolower(batchTable[i,1])
     }

     batchTable[i,1] = replaceChars(itmp, "  ", " ")
     
  }

  merged = merge(batchTable, ilabel, by.x=1,by.y=1, all=F)
  merged = as.matrix(merged)

  # no match  
  if ( dim(merged)[1]==0) {
      return (-1)
  }

  # we take the largest Batch_ID if there are multiple matches
  #
  no.matches = dim(merged)[1]

  return (as.integer(merged[no.matches,2]))
}


# direct means the table is directly uner the connectDB specified
#
accessWholeTableFromDB = function(servername, dbname, tablename, condition=NULL, outfields = "*"){
   channel <- odbcConnect(servername)
   if (is.null(condition)){
     query   <- paste("select ", outfields, " from ", dbname, ".dbo.[", tablename, "]", sep="")
   }else{
     query   <- paste("select ", outfields, " from ", dbname, ".dbo.[", tablename, "] where ", condition, sep="")
   }
   
   tb      <- sqlQuery(channel, query)
   odbcClose(channel)
   return (tb)
}


findMarkersInQTLBycM = function (eQTLonegene, markerinfoMx, markerindex)
{ 
  # look at markers on the same chromosome as the given QTL gene
  #
  #selchrom = eQTLonegene$chrom == markerinfoMx$marker_chrom
  #dist     = abs(eQTLonegene$qtl_max_pos - markerinfoMx$marker_pos[selchrom])
  selchrom  = (eQTLonegene[1] == markerinfoMx$marker_chrom)
  selchrom  = ifelse(is.na(selchrom), F, selchrom)

  if ( sum(selchrom,na.rm=T) == 0 ) {
     return ( c(-1,-1,-1) )
  }

  dist     = abs(eQTLonegene[2] - markerinfoMx$marker_pos[selchrom])

  selindex = markerindex[selchrom]

  # find the nearest marker
  #
  mIdx     = which.min(dist)
  #selindex[mIdx]

  myret = c(selindex[mIdx], (markerinfoMx$marker_pos[selchrom])[mIdx], 
                            (markerinfoMx$Base_pair[selchrom])[mIdx] )

  return ( myret)

}


############# network properties ###############################################

# examine the total number of links, # of nodes, average links per node, scale free fitting R^2 & slope
#
# the input is the first two columns in the network file
# 
# if removeRedantLinks==T, then, we first remove redundant links and then do subsequent analysis
#
getNetwork_Statistics_Scalefree_GeneLink = function(networkinput, colheader=F, directed=F, removeRedantLinks=F, 
                 downstream_layers=NULL, outputDir="", msep="\t", returnDegree=F, plotimg=T){

    imagetype   ="png"
    fname       =getFileName(networkinput)
    extname     =getFileExtension(networkinput)

    # use different tag, if from non-redundant network
    #
    if(removeRedantLinks) {
      fname = paste(fname, "_uniq", sep="")
    }

    logFname    =paste(outputDir, fname, "-log.xls",       sep='')
    imgScaleFree=paste(outputDir, fname, "-imgScaleFree.",imagetype,sep='')

    codepair= c(0,1)  #[1] for no connection, [2] for connection

    #0) read in pairwise network

    aaMatrix   <- read.delim(networkinput, sep=msep, header=colheader)
    aaMatrix   <- as.matrix(aaMatrix)    
    dim(aaMatrix)    
    
    #1) remove self-links & redundant links and save into corresponding files
    #
    sel = aaMatrix[,1] == aaMatrix[,2]
    aaMatrix = aaMatrix[!sel, ]

    if(removeRedantLinks) {
      aaMatrix = removeDuplicatedLinks(linkpairs=aaMatrix[,c(1:2)], directed=directed)

      uniquefnet  = paste(fname, ".pair",sep="")
      colnames(aaMatrix) = c("node.x", "node.y")

      write.table(aaMatrix, uniquefnet, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    
    }


    dim(aaMatrix)    
    
    edgesInNet = dim(aaMatrix)[1]
    
    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(aaMatrix[,1]) )
    allnodenames = c(allnodenames, as.character(aaMatrix[,2]) )
    
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    #name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    #if (no.uniquenames<= 4000){
    if (1==2){
      adjmatrix = makeAjacencyMatrix(inetmatrix=aaMatrix, coding=codepair,
                                     matrixsize=no.uniquenames, directed=directed, 
                                     myname2idxMatrix=name2idxMatrix)
      edgeIndices=getNetworknodeIdxByNames(netmatrix=aaMatrix, mappingMatrix=name2idxMatrix)

      # find hubs
      inlinks   = apply(adjmatrix, 2, sum) #degree(adjmatrix,cmode="indegree")
      outlinks  = apply(adjmatrix, 1, sum) #degree(adjmatrix,cmode="outdegree")

      if (directed){
        totallinks= inlinks + outlinks
      }else{
        totallinks= (inlinks + outlinks)/2 #because of matrix symmetry for undirectecd adj matrix
      }
    }

    totallinks= as.integer(nametable) # no of links for each node
    totalmatrix = cbind(names(nametable),  totallinks)

    if( !is.null(downstream_layers) ) {
       dns=rep(0, no.uniquenames)
       for(i in c(1:no.uniquenames) ) {
          if(i%%50 ==0)print(i)
         igenes=downStreamGenes(netpairs=aaMatrix, seednodes=uniquenames[i], N=downstream_layers, directed=directed)
         dns[i] = length(igenes)
       }
    }

      if(directed){
        # outlines
        dnodenames = as.character(aaMatrix[,1])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        outlinks = as.integer(as.matrix(iolinks[,2])) 


        # inlines
        dnodenames = as.character(aaMatrix[,2])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        inlinks = as.integer(as.matrix(iolinks[,2])) 

      }else{
        inlinks   = totallinks
        outlinks  = totallinks
      }

    hubidx    = order(-totallinks)


    # let's check whether there is a scale free topology
    suminfor   =ScaleFreePlot(totallinks,no.breaks=40, mtitle="",truncated1=T, tofile="")
    if(plotimg) {
       scaleResult=ScaleFreePlot(totallinks,no.breaks=40, mtitle="",truncated1=T, tofile=imgScaleFree, outputFitness=T)
    } else{
       scaleResult=ScaleFreePlot(totallinks,no.breaks=40, mtitle="",truncated1=T, tofile="", outputFitness=T)
    }

    # ------------------ output log file -----------
    #appendStringToFile(logFname, suminfor)

    # output network statistics 
    #
    mtitle = c("# of links", "# of genes", "links/gene", "scale R^2", "trunc. R^2", "slope")
    
    linkpernode = 2*edgesInNet/no.uniquenames

    finalMatrix = rbind( c(edgesInNet,no.uniquenames, linkpernode, scaleResult) )
    finalMatrix = round(finalMatrix, 2)
    
    colnames(finalMatrix) <- mtitle
    write.table(finalMatrix, logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE,append=F)


    # output in/out links for each gene
    #

    if( !is.null(downstream_layers) ) {
       linksMatrix            = cbind(uniquenames, inlinks, outlinks, totallinks, dns)
       colnames(linksMatrix) <- c("gene", "inlinks", "outlinks", "totallinks", "downstream_links")
    } else{
       linksMatrix            = cbind(uniquenames, inlinks, outlinks, totallinks)
       colnames(linksMatrix) <- c("gene", "inlinks", "outlinks", "totallinks")
    }

    write.table(linksMatrix[hubidx, ], logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    #rm(linksMatrix)
    rm(inlinks)
    rm(outlinks)
    rm(totallinks)
    rm(aaMatrix)
    collect_garbage()

    if(returnDegree) {
       return(linksMatrix[hubidx,])
    } else{
       return (finalMatrix)
    }
}


#------------------------------- centrality by pairs ---------------------------------------------
#
# to perform centrality analysis on large scale network, we avoid adjacency matrix
# the output is ordered based on total links
#
degree_ByLinkPairs = function(linkpairs, directed=F, cleangarbage=F){

    codepair= c(0,1)  #[1] for no connection, [2] for connection

    edgesInNet = dim(linkpairs)[1]
    
    # consider both columns
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(linkpairs[,1]) )
    allnodenames = c(allnodenames, as.character(linkpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)
    
    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    totallinks  = as.integer(nametable) # no of links for each node
    totalmatrix = cbind(names(nametable),  totallinks)

    if(directed){
        # outlines
        dnodenames = as.character(linkpairs[,1])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        outlinks = as.integer(as.matrix(iolinks[,2])) 


        # inlines
        dnodenames = as.character(linkpairs[,2])
        dnametable = table(dnodenames)
        duniquenames= names(dnametable)
        dmatrix     = cbind(names(dnametable), as.integer(dnametable ) )
        colnames(dmatrix) <- c("node", "links")

        iolinks = mergeTwoMatricesByKeepAllPrimary(primaryMatrix=cbind(uniquenames), minorMatrix=dmatrix, 
                                                  missinglabel="0", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
        inlinks = as.integer(as.matrix(iolinks[,2])) 

    }else{
        inlinks   = totallinks
        outlinks = totallinks
    }

    #hubidx    = order(-totallinks)

    # output in/out links for each gene
    #
    linksMatrix            = cbind(inlinks, outlinks, totallinks)
    colnames(linksMatrix) <- c("inlinks", "outlinks", "totallinks")
    rownames(linksMatrix) <- uniquenames

    rm(inlinks)
    rm(outlinks)
    rm(totallinks)

    if(cleangarbage) {
        collect_garbage()
    }

    return ( data.frame(linksMatrix) )
}


centralities=function(adjmatrix, gmode="graph", weighted=F)
{
   # 1. degree
   outdegree = apply(adjmatrix, 1, sum)    
   indegree  = apply(adjmatrix, 2, sum)
   if(gmode=="digraph") {
      totaldegree = outdegree + indegree
   }else{
      totaldegree = outdegree
   }

   #2. betweenness(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE,
   #         tmaxdev=FALSE, cmode="directed", geodist.precomp=NULL, 
   #         rescale=FALSE)
   betweenness = betweenness(adjmatrix, gmod=gmode) 

   #3. closeness(graph, v=V(graph), mode = "all") #"in", "out", "all"
   #
   closeness = closeness(adjmatrix, gmod=gmode, cmode="undirected") 
   
   #4. Clustering coefficient
   clusteringcoef = computeClusterCoefficient(adjmatrix, weighted=weighted)
   
   centrls = data.frame(cbind(indegree, outdegree, 
                              totaldegree, betweenness, closeness, clusteringcoef) )
   return(centrls)
}


# pair a given integer with each element in a vector
#
pairI_to_neighbors = function(si, ivect, removeMinus1=T){
  ino    = length(ivect)
  ipairs = cbind( rep(si,ino), ivect)
  isel   = (ipairs[,1] !=-1) & (ipairs[,2] !=-1)

  if(removeMinus1){
     if(sum(isel)==0){
        return (NULL)
     }else{
        return (ipairs[isel, ]) 
     }
  } else{
     return (ipairs)     
  }
}

pairs_To_Index = function(linkpairs, nodeIndex){
  mergedleft=mergeTwoMatricesByKeepAllPrimary(linkpairs, minorMatrix=nodeIndex, 
                                 missinglabel="", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
  mergedright=mergeTwoMatricesByKeepAllPrimary(mergedleft[,c(2,1,3)], minorMatrix=nodeIndex, 
                                  missinglabel="", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
  linksIndex = mergedright[,c(3,4)]
  linksIndex = as.matrix(linksIndex)
  ino.links  = dim(linksIndex)[1]
  linksIndex = matrix(as.integer(linksIndex), ino.links,2)
  colnames(linksIndex ) =c("src","dst")
  return (linksIndex)
}

# remove duplications including reversed duplications, and self links
#
removeDuplicatedLinks= function(linkpairs, directed=F){

    if(dim(linkpairs)[1]==1){
       return (linkpairs)
    }

    links = paste(linkpairs[,1], linkpairs[,2],sep="\t")

    # 1. remove duplications 
    #
    cleanedlinkMatrix = union(links, NULL)
    length(cleanedlinkMatrix)

    #ofname ="tmp.txt"
    #write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)

    # 2. remove inversed duplications
    #
    #linkMatrix <- read.delim(ofname, sep="\t", header=F) # not good for batch operation, ie, many of the same jobs running
    linkMatrix  <- getAllParts(cleanedlinkMatrix, "\t")
    dim(linkMatrix)
    linkMatrix  = as.matrix(linkMatrix)

    if(directed){
       return(linkMatrix)
    }

    if (dim(linkMatrix)[1]==1) {return(linkMatrix);}

    #  first, remove self-interactions are also removed
    #
    selSelfLinks      = linkMatrix[,1]==linkMatrix[,2]
    linkMatrix        = linkMatrix[!selSelfLinks,]
    cleanedlinkMatrix = cleanedlinkMatrix[!selSelfLinks]

    reversedLinks = paste(linkMatrix[,2], linkMatrix[,1],sep="\t")

    no.links = length(reversedLinks)
    reversedLinksWithIndex = cbind(reversedLinks, c(1:no.links) )

    merged = merge(cleanedlinkMatrix, reversedLinksWithIndex, by.x=1,by.y=1,all=F)
    dim(merged)

    if (dim(merged)[1]>0) {
        merged      = as.matrix(merged)
        removedCols = as.integer(merged[,2])

        # cosntruct non-duplicated interactions
        #
        dupLinks    = cleanedlinkMatrix[removedCols]
        dupLinksRev = reversedLinksWithIndex[removedCols]
        uniques     = NULL
        for(i in c(1:length(dupLinks)) )
        {
           found    = is.element(dupLinks[i], uniques)
           foundRev = is.element(dupLinksRev[i], uniques)
           combined = found | foundRev
           if (!combined){
               uniques=c(uniques,dupLinks[i])
           }
        }
        length(uniques)
        xlinkMatrix  = c(cleanedlinkMatrix[-removedCols], uniques)
        #write.table( cleanedlinkMatrix[-removedCols], ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        #write.table( uniques, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE,append=T)
    }else{
        #write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        xlinkMatrix  = cleanedlinkMatrix
    }

    #linkMatrix <- read.delim(ofname, sep="\t", header=F)
    #dim(linkMatrix)
    #linkMatrix  = as.matrix(linkMatrix)

    linkMatrix = getAllParts(xlinkMatrix, "\t")

    return(linkMatrix)
}


# get the network with the links between all nodes in specified by subnet
# keep the original order
# the third column can be index, so the linkpair index is also returned
getSubnetwork_LinkPairs = function(linkpairs, subnetNodes){
   mergeright = merge(linkpairs, subnetNodes, by.x=2, by.y=1, all=F)
   if( dim(mergeright)[1] == 0) {
     return (NULL)
   }

   mergeleft  = merge(mergeright, subnetNodes, by.x=2, by.y=1, all=F)

   #mergeright = merge(linkpairs, subnetNodes, by.x=1, by.y=1, all=F)
   #mergeleft2 = merge(mergeright, subnetNodes, by.x=2, by.y=1, all=F)

   if( dim(mergeleft)[1] == 0) {
     return (NULL)
   }

   return (as.matrix(mergeleft) )   
}


findNLayerNeighbors_LinkPairs = function(linkpairs, subnetNodes, nlayers=1, directed=F){

  #linkpairs=linkpairs; subnetNodes=ineighbors; nlayers=nlayers-1; directed=directed

   if(nlayers<=0) {return (NULL)}

   merged = merge(linkpairs, subnetNodes, by.x=1, by.y=1, all=F)
   merged <- as.matrix(merged)

   if(!directed) { # undirected networks
     mergeleft  = merge(linkpairs, subnetNodes, by.x=2, by.y=1, all=F)
     mergeleft  <- as.matrix(mergeleft)
     mergeleft  <- mergeleft[,c(2,1)] # keep the original link direction
     merged     <- rbind(merged, mergeleft)
   }

   if ( dim(merged)[1] ==0){return (NULL) }

   if( (dim(merged)[1]==1) & (merged[1,1]==merged[1,2]) ){return(merged)}

   dim1 = dim(merged)[1]
   if (dim1==0) { # no links
         return (NULL)
   }else if (is.null(dim1)){ # only one link
       merged = rbind(merged)
   }

   ineighbors =union(merged[,1], merged[,2])
   
   if (nlayers==1){
      res=getSubnetwork_LinkPairs(linkpairs, subnetNodes=ineighbors)
      return (res)
   }

   # stop earlier if no chnage
   #
   common = intersect(ineighbors, subnetNodes)
   if (length(common)==length(ineighbors)) {
      res=getSubnetwork_LinkPairs(linkpairs, subnetNodes=ineighbors)
      return (res)
   }
   
   ret=findNLayerNeighbors_LinkPairs(linkpairs,ineighbors, nlayers-1, directed)

   return (ret)
}


findNLayerNeighbors_LinkPairsOLd = function(linkpairs, subnetNodes, nlayers=1, directed=F){

  #linkpairs=linkpairs; subnetNodes=ineighbors; nlayers=nlayers-1; directed=directed

   if(nlayers<=0) {return (NULL)}

   merged = merge(linkpairs, subnetNodes, by.x=1, by.y=1, all=F)
   merged <- as.matrix(merged)

   if(!directed) { # undirected networks
     mergeleft  = merge(linkpairs, subnetNodes, by.x=2, by.y=1, all=F)
     mergeleft  <- as.matrix(mergeleft)
     mergeleft  <- mergeleft[,c(2,1)] # keep the original link direction
     merged     <- rbind(merged, mergeleft)
   }

   if (dim(merged)[1] ==0){return (NULL) }

   if( (dim(merged)[1]) ==1) {
       if(merged[1,1]==merged[1,2]) {return(merged)}
   }

   if(!directed) { # undirected networks
     mergeleft  = merge(linkpairs, subnetNodes, by.x=2, by.y=1, all=F)
     if (dim(mergeleft)[1] >0) {
       mergeleft  <- as.matrix(mergeleft)
       mergeleft  <- mergeleft[,c(2,1)] # keep the original link direction
       merged     <- rbind(merged, mergeleft)
     }
   }

   dim1 = dim(merged)[1]
   if (dim1==0) { # no links
         return (NULL)
   }else if (is.null(dim1)){ # only one link
       merged = rbind(merged)
   }

   merged=removeDuplicatedLinks(merged, directed)
   
   if (nlayers==1){
      return (merged)
   }

   # get nodes   
   ineighbors =union(merged[,1], merged[,2])

   # stop earlier if no chnage
   #
   common = intersect(ineighbors, subnetNodes)
   if (length(common)==length(ineighbors)) {
      return (merged)
   }
   
   ret=findNLayerNeighbors_LinkPairs(linkpairs,ineighbors, nlayers-1, directed)

   return (ret)
}


# the downstream genes of a given gene
#  dsnodes     -- current nodes in the downstream (gloabl variable for easy update)
#  startidx    -- start idx (search only those starting from here
#  adjlists    -- global adjlist
#
downStreamGenes=function(netpairs, seednodes, N= 100, directed=T, return_size=FALSE)
{
   prenodes = seednodes
   cnt = N
   while(T) {
      retlinks = findNLayerNeighbors_LinkPairs(linkpairs=netpairs, subnetNodes=prenodes, 
                      nlayers=1, directed=directed)

      if(is.null(retlinks)){return (NULL); }

      curnodes = union(retlinks[,1],retlinks[,2]) 
      pcdiff   = setdiff(curnodes, prenodes)
      prenodes = curnodes
      
      if(length(pcdiff)==0){break}

      cnt= cnt-1
      if (cnt==0) {break}
   }

   if(is.null(retlinks)){

      if(return_size) {
        return (0)
      } else{
        return (NULL)
      }

   } else {

      if(return_size) {
        return( length(setdiff(curnodes, seednodes) ) )
      } else{
        return(curnodes)
      }
   }
}


downStreamGenes0=function(netpairs, seednodes)
{
   retlinks = findNLayerNeighbors_LinkPairs(linkpairs=netpairs, subnetNodes=seednodes, 
                      nlayers=100, directed=T)

   if(is.null(retlinks)){return (NULL)}
   else {
      return( union(retlinks[,1],retlinks[,2]) )
   }
}


mapPairlinkIDS2genesymbol2=function(fpairlink, inpath, fmapping, myoutdir)
{
   infullname = paste(inpath, fpairlink, sep="")
   linkMatrix <- read.delim(infullname,sep="\t", header=F)
   dim(linkMatrix)
   
   transMatrix <- read.delim(fmapping,sep="\t", header=F)
   dim(transMatrix)

   # src dst srcgs   
   leftmatrix = merge(linkMatrix, transMatrix, by.x=1, by.y=1, all=F)

   # src dst srcgs   & yid gs = dst src srcgs dstgs
   rightmatrix = merge(leftmatrix, transMatrix, by.x=2, by.y=1, all=F)

   fname     =getFileName(fpairlink)
   ofname    =paste(myoutdir, fname, "_gs.pair", sep='')

   write.table( rightmatrix[,c(3,4)], ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
       
   return (rightmatrix[,c(3,4)])
}

# Here we use hierarchical clustering to idenitify Connected Components
#
find_ConnectedComponents_linkpairs = function(linkpairs, minModuleSize=10, cclabelstart=1){

    # find unique node names 
    allnodenames = NULL
    allnodenames = c(allnodenames, as.character(linkpairs[,1]) )
    allnodenames = c(allnodenames, as.character(linkpairs[,2]) )
     
    nametable = table(allnodenames)
    length(nametable)

    uniquenames     = names(nametable)
    no.uniquenames  = length(uniquenames)

    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )

    #*-------------------------------------------------------------------------------------
    #* prepare matrix
    #*   
    adjmatrix = makeAjacencyMatrix(inetmatrix=linkpairs[,c(1,2)], 
                                   coding=c(0,1),
                                   matrixsize = no.uniquenames, 
                                   myname2idxMatrix= name2idxMatrix,
                                   directed=F)

    colnames(adjmatrix) <- uniquenames
    rownames(adjmatrix) <- uniquenames

    #"single" ensures that if a node is connected to at least one memeber in a candidate cluster
    # then this node will have a disctance 0 to the cluster, i.e, the concept of connect component
    #
    h1row <- hclust(as.dist(1-adjmatrix),method="single")

    #collect_garbage()
    #plot(h1row, labels=F, xlab="",ylab="",main="",sub="")
    #plot(h1row, xlab="",ylab="",main="",sub="")

    # ----------------- output Hierarchical Clustering image ------------------------------
    #plot(h1row, labels=F, xlab="",ylab="",main="",sub="")

    #*-------------------------------------------------------------------------------------
    #* detect and label modules in TOM based on the hierarchical clustering dendrogram
    #*              
    myheightcutoff    = 0.99
    colcode.reduced   = moduleDetectLabel(hiercluster=h1row, myheightcutoff, 
                                          minsize1=minModuleSize, 
                                          useNumberAsLabel=T, startlabel=cclabelstart)
    table(colcode.reduced)
    #plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=FALSE)

    cc= cbind(uniquenames, colcode.reduced)
    #colnames(cc) <- c("node", "ConnectedComponent")

    return(cc)
}


# here we return index
#
find_ConnectedComponents_adjmatrix = function(adjmatrix, minModuleSize=10, cclabelstart=1){

    dims = dim(adjmatrix)[1]

    # single item
    if (is.null(dims)){ # single item
          ipad = patchZeros(1)
          return ( rbind(c(1,ipad) ) )
    }

    # assign unique IDs
    no.uniquenames  = dims
    uniquenames     = c(1:no.uniquenames)

    colnames(adjmatrix) <- as.character(uniquenames)
    rownames(adjmatrix) <- as.character(uniquenames)

    #"single" ensures that if a node is connected to at least one memeber in a candidate cluster
    # then this node will have a disctance 0 to the cluster, i.e, the concept of connect component
    #
    h1row <- hclust(as.dist(1-adjmatrix),method="single")

    #collect_garbage()
    #plot(h1row, labels=F, xlab="",ylab="",main="",sub="")
    #plot(h1row, xlab="",ylab="",main="",sub="")

    # ----------------- output Hierarchical Clustering image ------------------------------
    #plot(h1row, labels=F, xlab="",ylab="",main="",sub="")

    #*-------------------------------------------------------------------------------------
    #* detect and label modules in TOM based on the hierarchical clustering dendrogram
    #*              
    myheightcutoff    = 0.99
    colcode.reduced   = moduleDetectLabel(hiercluster=h1row, myheightcutoff, 
                                          minsize1=minModuleSize, 
                                          useNumberAsLabel=T, startlabel=cclabelstart)
    table(colcode.reduced)
    #plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=FALSE)

    cc= cbind(uniquenames, colcode.reduced)
    #colnames(cc) <- c("node", "ConnectedComponent")

    return(cc)
}










##################################### K-CORES k-cores #################################################
#
# Recursively remove all vertices of degree less than K, UNTIL all vertices in the remaining# graph have at least  degree K
#
# return a matrix with the 1st column as node names
#                      the 2nd column as k-shell
#                      the rest columns as connected components in each core
#                               missing values (NA) means that these nodes are not in the cores
#                               -1 means that the nodes are isolated in the cores
#
# if fkey is not null, then it saves the node-kcore matrix and the link-kcore matrix into files
#
#
find_kcores = function (linkpairs, min_core=1, minCCsize=1, returnNodeCores=F, fkey=NULL) {
    kcoresCC      = NULL
    kshell        = NULL

    newlinkpairs  = linkpairs
    degrees       = degree_ByLinkPairs (linkpairs, directed=F)

    i = min_core
    while (T) {

       # i-core, recursively remove nodes of degree less than i
       #  until all nodes in the remaining graph have at least degree i
       #
       while ( T) {
           # get node names
           currnodenames = rownames(degrees)

           icoreSel  = degrees$totallinks >= i
           icoreSize = sum(icoreSel)

           # i-shell (actually i-1 shell)
           ishellSel  = !icoreSel
           ishellSize = sum(ishellSel, na.rm=T)
           if(ishellSize >0){
              ishell     = cbind(currnodenames[ishellSel],  rep(i-1,ishellSize) )
              kshell     = rbind(kshell, ishell)
           }

           # no more nodes with i-core
           if ( icoreSize <=1 ){
             break
           }

           # update adj matrix & network properties
           #
           newlinkpairs = getSubnetwork_LinkPairs(linkpairs=newlinkpairs, 
                                                  subnetNodes=currnodenames[icoreSel])
           if (is.null(newlinkpairs)) {
               icoreSize = 0 #no i-core found
               break
           }

           degrees      = degree_ByLinkPairs (newlinkpairs, directed=F)

           if (min(degrees$totallinks) >= i){
              break
           }
       }

       # no more nodes with i-core
       if ( icoreSize ==0){
         break
       }

       mystr=paste("Shell=", as.character(i), " core=", as.character(min(degrees$totallinks)), "\n",sep="")
       print(mystr)
       #if(i>=8){
       #   print(newlinkpairs )
       #   print(degrees)
       #}

       # finally, we got i-core nodes
       #
       # update node names in i-core
       currnodenames = rownames(degrees)

       icoreSize = length(currnodenames)
       icore     = cbind(currnodenames,  rep(i, icoreSize) )

       # get connected components for the current i-core
       #
       icorelinkpairs = getSubnetwork_LinkPairs(newlinkpairs, currnodenames)

       # nodes in the current icore have no links
       # so we assign a "grey" CC label
       #
       if( dim(icorelinkpairs)[1]==0){
         icclabel = rep("-1", icoreSize)
         icoreCC  = cbind(icore, icclabel)
         colnames(icoreCC) <- c("node", paste("core",as.character(i),sep=""),
                                        paste("cc_core",as.character(i),sep="") )

       }else{ # dectect cc, return is the (nodename, cclabel)
         icclabel = find_ConnectedComponents_linkpairs(linkpairs=icorelinkpairs, 
                                          minModuleSize=minCCsize,
                                          cclabelstart=1)
         icoreCC = merge(icore, icclabel, by.x=1, by.y=1, all.x=T)
         icoreCC = as.matrix(icoreCC)
         icoreCC = ifelse(is.na(icoreCC), "-1", icoreCC ) # -1: not in any CC
         icoreCC = ifelse(icoreCC=="grey", -1,  icoreCC )

         colnames(icoreCC) <- c("node", paste("core",as.character(i),sep=""),
                                        paste("cc_core",as.character(i),sep="") )
       }

       # so, CC's in i-core are put as a column, the missing vALUE indicates 
       #  these nodes are not in the current i-core
       #
       if (is.null(kcoresCC)){
          kcoresCC = icoreCC[,c(1,3)]
       } else{
          kcoresCC = merge(kcoresCC, icoreCC[,c(1,3)], by.x=1, by.y=1, all=T)
       }

       #print(paste(as.character(i), " coreness"))
       #print(icoreCC)

       # update adj matrix & network properties
       #
       #newlinkpairs = getSubnetwork_LinkPairs(linkpairs=newlinkpairs, subnetNodes=currnodenames[icoreSel])
       #degrees      = degree_ByLinkPairs (newlinkpairs, directed=F)

       i=i+1
    }

    colnames(kshell) <- c("node", "k_shell")
    kshell_coresCC = merge(kshell, kcoresCC, by.x=1, by.y=1, all=T)

    # coreness of each node
    # node	k_shell	cc_core1	cc_core2	cc_core3
    # 1a	1	1		NA		NA
    #
    if(!is.null(fkey)){
       fnodecores = paste(fkey, ".xls", sep="")
       write.table(kshell_coresCC, fnodecores, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
    }

    # core based link pairs
    #  a link belongs to a CC in a core iff both nodes are in the CC of the core
    #
    ncols      = dim(kshell_coresCC )[2]
    totalcores = ncols -2
    fnodecores = paste(fkey, "-linkpairs.xls", sep="")
    mergeleft  = merge(linkpairs, kshell_coresCC, by.x=2, by.y=1)
    mergeright = merge(mergeleft, kshell_coresCC, by.x=2, by.y=1)
    mergeright = as.matrix(mergeright )

    # src   dst  k_shell  cc_core1 cc_core2  cc_core3 k_shell cc_core1	cc_core2 cc_core3
    # 1     2    3                                     
    # link's shell is the min of the shells of the two nodes
    #
    newCCinCores = NULL
    shellmatrix = cbind( as.integer(mergeright[,3]), as.integer(mergeright[,totalcores+4]) )
    linkshell   = apply(shellmatrix, 1, min)
    newCCinCores = cbind(mergeright[,c(1:2)], linkshell)

    # look at whether two nodes of each link are in the same CC of the same core
    #
    for (i in c(1:totalcores) ){
        icorelinkflag = mergeright[,i+3]==mergeright[,totalcores+4+i]
        newCCinCores  = cbind(newCCinCores, mergeright[icorelinkflag, i+3] )
    }

    colnames(newCCinCores) <- c("src","dst", "kshell", colnames(kshell_coresCC)[3:(ncols)] )

    if(!is.null(fkey)){
        flinkcores = paste(fkey, "-cc-pairs.xls", sep="")
        write.table(newCCinCores, flinkcores, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
     }

    if(returnNodeCores){
        return(kshell_coresCC)
    }else{
        return (newCCinCores)
    }
}


#
#**********************************  END of K-cores  **********************************





##################################### K-clique #################################################
#
# Recursively remove all vertices of degree less than K, UNTIL all vertices in the remaining
# graph have at least  degree K
#
# return a matrix with the 1st column as node names
#                      the 2nd column as k-shell
#                      the rest columns as connected components in each core
#                               missing values (NA) means that these nodes are not in the cores
#                               -1 means that the nodes are isolated in the cores
#
# if fkey is not null, then it saves the node-kclique matrix and the link-kclique matrix into files
#
#

# check if setA is a subset of setB
setsub = function(setA,setB){
    if (length(setA) > length(setB)){
        return (F)
    }
    setAB = union(setA,setB)
    equal = setequal(setAB,setB)
    return (equal)
}

# return a vector of boolean flags indicating whether each element in A is in B
setElementInSet = function(setA,setB){
    found = rep(F, length(setA))
    for (i in c(1:length(setA)) ){
       idiff    = setdiff(setA[i], setB)
       ifound   = (length(idiff)==0)
       found [i]= ifound
    }
    return (found)
}

setElementInSet_Fast = function(setA,setB){
    found = rep(FALSE, length(setA))
    axed = cbind(setA, c(1:length(setA) ))
    merged=merge(axed, cbind(setB), by.x=1, by.y=1, all=F)
    if (dim(merged)[1]==0) {return(found) } 
    merged=as.matrix(merged)
   
    found[as.integer(merged[,2])] = TRUE

    return (found)
}


# check if a set is a subset of any set in a set list
#
setInSets = function (setC, setlist){
   nosets = length(setlist)
   for (i in c(1:nosets) ){
       isSub = setsub(setC, setlist[[i]])
       if (isSub){ #setC is a subset of setlist[[i]]
          return (T)
       }
   }
   return (F)
}

setInSets_ret_vect = function (setC, setlist){
   nosets = length(setlist)
   res    = rep(F, nosets)
   for (i in c(1:nosets) ){
       isSub = setsub(setC, setlist[[i]])
        res[i] = isSub
   }
   return (res)
}


# globalCliquesMatrix is a mtrix with cols as gene index
#   and rows as cliques (sets)
# So, we first get a matrix with only the cols given by setC 
#  then, compute the sum of each row
# If one sum value is equal to the number of elements in setC,
#  then, setC is a subset of setInSets_matrix, otherwise, not
#
setInSets_matrix = function (setC) {
   rsum  = apply(globalCliquesMatrix[,setC], 1, sum)
   found = length(setC)==rsum
   asum  = sum(found)
   return (asum>=1) 
}


# ------ global variables ------------

#~~~~~~~ input ~~~~~~~~~~~~~~~~~~~~
# globalAdjLists
# globalLinkpairsIndex
# globalFullLinked
# globalAdjSize

#~~~~~~~ output ~~~~~~~~~~~~~~~~~~~
# globalCliques
# globalCliqueSize

kcore_subnets = function(linksindex, nodesindex=NULL, kcore=3, name2integer=T) {

    if(!is.null(nodesindex)) {
        mright  = merge(linksindex, nodesindex, by.x=2, by.y=1, all=F)
        if( dim(mright)[1]==0){
            return (list(NULL, NULL, F) )
        }

        orgLinks= merge(mright,     nodesindex, by.x=2, by.y=1, all=F)
        if( dim(orgLinks)[1]==0){
            return (list(NULL, NULL, F) )
        }
    }else{
        orgLinks= linksindex
    }

 
    while(T) {
      no.orglinks   = dim(orgLinks)[1]
      degrees       = degree_ByLinkPairs (orgLinks, directed=F)
      dim(degrees)
      selected    = degrees$totallinks >= kcore
      no.selnodes = sum(selected)

      if (no.selnodes==0) {
          return (list(NULL, NULL, F) )
      }

      selectedNodes    = rownames(degrees)[selected]
      selected_nolinks = degrees$totallinks[selected]/no.selnodes

      # =========== important ======================
      if(name2integer){
          selectedNodes = sort( as.integer(selectedNodes) )
      }

      if (no.selnodes < kcore+1){ # no enough nodes, so no enough links for each node
          return (list(NULL, NULL, F))
      }

      selectedLinks = getSubnetwork_LinkPairs(orgLinks, subnetNodes=selectedNodes)
      if (is.null(selectedLinks)) {
          return (list(NULL, NULL, F))
      }

      selectedLinks = as.matrix(selectedLinks)
      no.newlinks   = dim(selectedLinks)[1]

      if (no.newlinks==0) {
          return (list(NULL, NULL, F))
      }

      orgLinks = selectedLinks

      if(no.newlinks==no.orglinks){
          break
      }

    }

    # check if these nodes are fully linked
    #
    fulllinks  = no.selnodes*(no.selnodes -1)/2
    fulllinked = no.newlinks == fulllinks

    rets = list(selectedNodes, selectedLinks, fulllinked)

    return (rets)
}

#v=775
#ret=kcore_subnets(selectedLinksIndex, c(v, globalAdjLists[[v]]), kcore=5)

subnet_is_clique = function(nodes, links){
   maxlinks = nodes*(nodes-1)/2
   ratio    = links/maxlinks
   isclique = ratio==1
   return (isclique)
}

####################################################################################
################ find all cores and return the core nodes ##########################
#
#  nodeMembershipOnly==T
#       return the T/F corenodesMatrix
#  nodeMembershipOnly==F  
#       return( list(allcores, coresnodes, coreslinks, oressizes, iscliques) )
#
#  returnLinkID=T: this means that the third column of linksindex is IDs of links
#
KCores_ALL = function(linksindex, nodesindex=NULL, minkcore=3, name2integer=T, nodeMembershipOnly=T, minkcoreOnly=F, returnLinkID=F) {

    if(!is.null(nodesindex)) {
        mright  = merge(linksindex, nodesindex, by.x=2, by.y=1, all=F)
        if( dim(mright)[1]==0){
            #return (list(NULL, NULL, F) )
            return (list(NULL, NULL, NULL, NULL, NULL) )
        }

        orgLinks= merge(mright,     nodesindex, by.x=2, by.y=1, all=F)
        if( dim(orgLinks)[1]==0){
            #return (list(NULL, NULL, F) )
            return (list(NULL, NULL, NULL, NULL, NULL) )
        }
    }else{
        orgLinks= linksindex
    }

    degrees      = degree_ByLinkPairs (orgLinks, directed=F)

    kcore        = minkcore
    reachMaxcore = F
    maxindex     = max(degrees$totallinks) #max(linksindex)

    #nodeIdxSequence = c(1:maxindex)
    corenodesMatrix = NULL
    allcores        = NULL

    if(!nodeMembershipOnly) {
        kcoreMemberNodes = as.list(rep(0,maxindex))
        kcoreMemberLinks = as.list(rep(0,maxindex))
    }
    kcoreFullLinked  = rep(0,maxindex)
    kcoresSizes      = rep(0,maxindex)

    corecnt= 0
    while(T) { #all cores

        while(T) { #kcore
          no.orglinks   = dim(orgLinks)[1]
          degrees       = degree_ByLinkPairs (orgLinks, directed=F)
          dim(degrees)
          selected    = degrees$totallinks >= kcore
          no.selnodes = sum(selected)

          if (no.selnodes==0) {
              reachMaxcore=T
              break
          }

          selectedNodes    = rownames(degrees)[selected]
          selected_nolinks = degrees$totallinks[selected]/no.selnodes

          # =========== important ======================
          if(name2integer){
              selectedNodes = sort( as.integer(selectedNodes) )
          }

          if (no.selnodes < kcore+1){ # no enough nodes, so no enough links for each node
              reachMaxcore=T
              break
          }

          selectedLinks = getSubnetwork_LinkPairs(orgLinks, subnetNodes=selectedNodes)
          if (is.null(selectedLinks)) {
              reachMaxcore=T
              break
          }

          selectedLinks = as.matrix(selectedLinks)
          no.newlinks   = dim(selectedLinks)[1]

          if (no.newlinks==0) {
              reachMaxcore=T
              break
          }

          orgLinks = selectedLinks

          # record the nodes in the core
          if(no.newlinks==no.orglinks){
              corecnt   = corecnt + 1
              allcores  = c(allcores, kcore)

              if(nodeMembershipOnly) {
                 coreMembership      = rep(F, maxindex)
                 coreMembership[selectedNodes ] = T
                 corenodesMatrix = cbind(corenodesMatrix, coreMembership)
              }else{
                 kcoreMemberNodes[[corecnt]] = selectedNodes
                 if(returnLinkID) {
                     kcoreMemberLinks[[corecnt]] = as.integer(selectedLinks[,3])
                 }else{
                     kcoreMemberLinks[[corecnt]] = selectedLinks
                 }
              }      

              # clique or not
              kcoreFullLinked[corecnt] = subnet_is_clique(length(selectedNodes), dim(selectedLinks)[1])
              kcoresSizes[corecnt]     = length(selectedNodes)

              break
          }

        } #kcore


        if (reachMaxcore){
          break
        }

        # run only on the minkcore
        #
        if( minkcoreOnly){
           break
        }

        kcore = kcore + 1

    }

    if(nodeMembershipOnly) {
       colnames(corenodesMatrix) <- allcores
       return (corenodesMatrix)   
    } else{
       # NULL is OK, list(NULL) ==> [[1]] NULL
       #
       if(!is.null(allcores) ){
           length(kcoreMemberNodes) <- corecnt
           length(kcoreMemberLinks) <- corecnt
           length(kcoresSizes)      <- corecnt
           length(kcoreFullLinked)  <- corecnt
           finalList = list(corenesses= allcores, 
                        coresnodes= kcoreMemberNodes, 
                        coreslinks= kcoreMemberLinks,
                        coressizes= kcoresSizes, 
                        iscliques = kcoreFullLinked)
       } else{
           finalList = list(corenesses= NULL,
                        coresnodes= NULL, 
                        coreslinks= NULL,
                        coressizes= NULL, 
                        iscliques = NULL)
       }

       return (finalList)
    }
}

#####################################################################################
#  kcores with individual nodes's neighbors being cored
#
#  inputs are global variables:
#      selectedLinksIndex
#      globalAdjLists
#      no.nodes
#
# Notice that if we want k-core, then for the neighborhood of a node we need only
#       (k-1)-core. Keep in mind that
#        K-1 is recorded in global_nodeneighbor_cores_sizes
#
# *Notice that to assignm value to individual elements in a global variable, 
#   you have to use <<- instead of "=", otherwise, there is no change, check the
#   following example
#
if(F){
 # global variable
 a=c(1:10)

 gvar=function(){
     b<<-as.list(rep(10, 5))
 }
 cf = function(x){
     b[[3]] <<- x
 } #good

 cf2 = function(x){a[[6]]=x}
 #cf(10), cf2(11), cf(20)
}

KCores_ALL_in_neighbors = function(minkcore=3, returnLinkID=F) {

    # global variables
    #
    # initialization for each node
    # 1) coreness node/link members
    # 2) no of nodes in each coreness
    # 3) min/max coreness
    #
    global_nodenghr_cores       <<- as.list( rep(NA, no.nodes) )
    global_nodenghr_minmaxcore  <<- matrix(0, no.nodes, 2)

    # estimated max coreness
    #    
    avglinks      = no.links/no.nodes
    maxcore_rough = as.integer(4*avglinks)
    
    # for each core of each node
    #
    #tmplist = as.list( rep(NA, maxcore_rough) )
    #global_nodenghr_cores  = ifelse(is.na(global_nodeneighbor_cores), tmplist)
    #for (i in c(1:no.nodes) ){
    #    global_nodenghr_cores[[i]]      = tmplist
    #}
    #source("C:/Zhang/BinBlast/R-tomfunctions.R")
    #start.time <- proc.time()[3]

    print(paste("find all cores for the neighbors of each of ", no.nodes, " nodes", sep="") )
    cnt=0

    # count how many links in neighbor's cores
    #
    global_nolinks_nghbr_cores <<- rep(0, maxcore_rough)
    
    for (v in c(1:no.nodes)){
       print(as.character(v))
       
       if(returnLinkID) {
         vCores = KCores_ALL(linksindex=linkpairsIndexed, 
                   nodesindex=globalAdjLists[[v]], minkcore=minkcore-1, name2integer=T, 
                   nodeMembershipOnly=F, returnLinkID=returnLinkID)
       }else{
         vCores = KCores_ALL(linksindex=selectedLinksIndex, 
                   nodesindex=globalAdjLists[[v]], minkcore=minkcore-1, name2integer=T, 
                   nodeMembershipOnly=F, returnLinkID=returnLinkID)
       }


       if(!is.null(vCores$corenesses) ){
          #print (v)
          cnt= cnt + 1
          global_nodenghr_cores[[v]]     <<- vCores
          global_nodenghr_minmaxcore[v,] <<- c(min(vCores$corenesses), max(vCores$corenesses) )
          
          # for each core, we count how many links there
          #
          vcoresizes = vCores$coressizes
          length(vcoresizes) = maxcore_rough
          vcoresizes = ifelse(is.na(vcoresizes), 0, vcoresizes)
          global_nolinks_nghbr_cores <<- global_nolinks_nghbr_cores + vcoresizes
       }
    
    }

    # count of links is duplicated
    global_nolinks_nghbr_cores <<- global_nolinks_nghbr_cores

    #preprocessing_Time <- proc.time()[3]-start.time
    #preprocessing_Time
    #cnt

    collect_garbage()

    return (cnt)
}

# remove one node in the neighbors
#
#
# global_nodenghr_cores[[]]$corenesses
# global_nodenghr_cores[[]]$coresnodes
# global_nodenghr_cores[[]]$coreslinks
# global_nodenghr_cores[[]]$coressizes
# global_nodenghr_cores[[]]$iscliques
#
# global_nolinks_nghbr_cores
# global_nodenghr_minmaxcore
#

update_nodenghr_cores=function(nodenghr_cores, removegene=NULL, coreidx, use_link_id=F)
{

           curcoreness = nodenghr_cores$corenesses[coreidx]

           #print(curcoreness)
           #print(nodenghr_cores$coresnodes[[coreidx]])

           finalList   = list(corenesses= NULL,
                        coresnodes= NULL, 
                        coreslinks= NULL,
                        coressizes= NULL, 
                        iscliques = NULL)

           if(curcoreness==0){
               return (NULL)
           }

           # the node to be removed is not a neighbor of the current node
           #
           isneighbor= is.element(removegene, nodenghr_cores$coresnodes[[coreidx]])
           if(!isneighbor){
               return (NULL)
           }

           # ii) remove the current node and its links
           #
           # ii.1) remove it from LinkPairs
           #
           if(use_link_id){
             # get actual links and link IDs from the node's link IDs
             #
             ilinks     = linkpairsIndexed[ nodenghr_cores$coreslinks[[coreidx]], ]
             selLeft    = ilinks[,1] != removegene
             selRight   = ilinks[,2] != removegene             
           } else{
             selLeft    = nodenghr_cores$coreslinks[[coreidx]][,1] != removegene
             selRight   = nodenghr_cores$coreslinks[[coreidx]][,2] != removegene
           }

           sel        = selLeft & selRight
           no.remains = sum(sel)
           
           if(no.remains==0){
               return (finalList)
           }

          # redo the core analysis on the updated neighbors but only look at the same coreness
          #
          if(use_link_id){
            vCores = KCores_ALL(linksindex=rbind(ilinks[sel, ]), 
                   nodesindex=NULL, minkcore=curcoreness, name2integer=T, 
                   nodeMembershipOnly=F, minkcoreOnly=T, returnLinkID=use_link_id)
          } else{
            vCores = KCores_ALL(linksindex=rbind(nodenghr_cores$coreslinks[[coreidx]][sel,]), 
                   nodesindex=NULL, minkcore=curcoreness, name2integer=T, 
                   nodeMembershipOnly=F, minkcoreOnly=T, returnLinkID=use_link_id)
          }
          return (vCores)
}


find_maxcoreness=function(orgLinks, returnCoreSize=F){

    if(is.null(orgLinks)) {
         return (0)
    }

    newlinks= orgLinks
    degreesOrg = degree_ByLinkPairs (newlinks, directed=F)

    # find max core
    icore = max( degreesOrg$totallinks)
    no.orglinks = dim(newlinks)[1]
    while(T) {
          selected    = degreesOrg$totallinks >= icore
          nosels      = sum(selected)

          # icore is too big, then there is not enough members
          #  an icore requres at least (icore + 1) nodes
          #
          if (nosels < icore+1) {
               icore = icore - 1
               next
          }else{
               break
          }
    }

    if (icore ==1 ){
       return (icore)
    }

    while (T) { # choose k-core

        # from the very start
        newlinks = orgLinks
        found    = F
        firstrun = T

        while(T) {
          if (firstrun) {
              degrees = degreesOrg # no need to compute links again
              firstrun = F
          } else{
              degrees = degree_ByLinkPairs (newlinks, directed=F)
          }

          no.orglinks = dim(newlinks)[1]
          selected    = degrees$totallinks >= icore
          nosels      = sum(selected)

          # if no nodes was selected or only one node was selected, stop
          if (nosels <=1 ){
             break
          }

          selectedNodes = rownames(degrees)[selected]
          selectedLinks = getSubnetwork_LinkPairs(orgLinks[,c(1,2)],subnetNodes=selectedNodes)
          if (is.null(selectedLinks)) {
              break
          }

          selectedLinks = as.matrix(selectedLinks)
          no.newlinks   = dim(selectedLinks)[1]

          newlinks = selectedLinks

          if(no.newlinks==no.orglinks){
             found=T
             break
          }

        } #while
        #dim(selectedLinks)

        if(found){
           break
        }
        #if (icore <=2){
        #    return (icore)
        #}

        icore = icore - 1

   } #while

   if ( !returnCoreSize) {
      return (icore)    
   } else{
      return ( c(icore,length(selectedNodes), no.newlinks) )
   }
}


#********************************************************************************
#
# * notice that union(list(c(1,2),c(1,3)),list(c(1,2)) ) is different from
#               union(list(c(1,2),c(1,3)),list(c(1:2)) ) 
# the first yields right answer
# also          union(list(c(1,2),c(1,3)),list(c(1:2)) ) gives right answer
# so, internally in R, c(1,2) is different from c(1:2)
#
#--------------------------------------------------------------------------------
#
# initially, setA is empty and setB includes all neighbors of curnode
#  Kclique is the number of nodes needed from setB in order to form a k-clique
#   which includes curnode
#
# global_nodenghr_cores[[]]$corenesses
# global_nodenghr_cores[[]]$coresnodes
# global_nodenghr_cores[[]]$coreslinks
# global_nodenghr_cores[[]]$coressizes
# global_nodenghr_cores[[]]$iscliques
#

recursive_Cliques_nodenghr = function(curnode, setA, setB, neighborlinks, Kclique=3, coreidx, use_link_id=F){

       #print (curnode)
       # at the beginning, global_nodenghr_cores[[v]]$corenesses[coreidx] is not NULL, but it may be NULL after
       # several runs when remove the nodes which have been visited and have core(s)
       #
       # NULL means that the current node don't have enough neighbors (>=Kclique)
       #  so it cannot be in a clique of size = Kclique
       #

       if( is.null(global_nodenghr_cores[[curnode]]$corenesses[coreidx]) ){
           return (0)
       }

       #print(paste("curnode ", curnode, " coresness:"))
       #print(global_nodenghr_cores[[curnode]]$corenesses)

       # this node doesn't have a neighborhood including k-core with k < coreidx-th
       #
       if( length(global_nodenghr_cores[[curnode]]$corenesses)<coreidx){
           return (0)
       }

       if( global_nodenghr_cores[[curnode]]$corenesses[coreidx] ==0 ){
           return (0)
       }
       
       # move up v's neighbors to setA
       #
       setAA      = union(setA, curnode)
       size.setAA = length(setAA)

       #print("Set A, , currnode, B: ")
       #print(setA)
       #print(curnode)
       #print(setB)

       if (is.null(setB) ){
          setBB = global_nodenghr_cores[[curnode]]$coresnodes[[coreidx]]
       }else{
          setBB = setdiff(setB, curnode) # remove v from setB
          #update nodes
          setBB = intersect(setB, global_nodenghr_cores[[curnode]]$coresnodes[[coreidx]] )
       }

       size.setBB = length(setBB)

       # check if it is already a clique
       #
       if (size.setBB <= 1 ){#setBB is null or has only one node, setAA+setBB already form a clique
          isclique = TRUE
       }else{
          # update links
          newneighborlinks = getSubnetwork_LinkPairs(neighborlinks, setBB)
          if(is.null(newneighborlinks)){
              isclique = -1
          } else{
              isclique = subnet_is_clique(nodes=length(setBB), links=dim(newneighborlinks)[1])
          }
       }

       if(isclique==T){
           #print("Early: ")
           #print(sort(union(setAA, setBB) ))
           if( size.setBB == Kclique-size.setAA){
               currclique  = sort(union(setAA, setBB) )
               currcliqueL = list(currclique)

               #set global variable
               #
               #found   = setInSets(currclique, globalCliques)
               found   = setInSets_matrix(currclique)

               #print("isclique==T")
               #print(currclique)
               if(!found){
                  globalCliques   <<- c(globalCliques,  currcliqueL)

                  xcliq = globalNodeZeroIndex
                  xcliq[currclique] = 1
                  globalCliquesMatrix   <<- rbind(globalCliquesMatrix, xcliq)

                  return (1)
               }
               return (1)
           }
           return (0)
       }

       # the remaining nodes in setBB don't have any link with each other
       # so we stop the recursive process
       #
       if(isclique==-1){
           if( size.setAA == Kclique){
               currclique = sort(setAA)
               currcliqueL = list(currclique)
               #set global variable
               #
               #found   = setInSets(currclique, globalCliques)
               found   = setInSets_matrix(currclique)

               #print("isclique==-1")
               #print(currclique)

               if(!found){
                  globalCliques   <<- c(globalCliques,  currcliqueL)

                  xcliq = globalNodeZeroIndex
                  xcliq[currclique] = 1
                  globalCliquesMatrix   <<- rbind(globalCliquesMatrix, xcliq)

                  return (1)
               }
               return (1)
           }else if(size.setAA == Kclique-1) {

               for(ek in setBB) {
                   currclique = sort(c(setAA, ek) )
                   currcliqueL = list(currclique)
                   #set global variable
                   #
                   #found   = setInSets(currclique, globalCliques)
                   found   = setInSets_matrix(currclique)

                   if(!found){
                      globalCliques   <<- c(globalCliques,  currcliqueL)

                      xcliq = globalNodeZeroIndex
                      xcliq[currclique] = 1
                      globalCliquesMatrix   <<- rbind(globalCliquesMatrix, xcliq)
                   }
               }
               
               return (1)
           }

           return (0)
       }

       #print( paste("SetAA=", concatenate(as.character( sort(setAA)),","), sep="") )
       #print( paste("SetBB=", concatenate(as.character( sort(setBB)),","), sep="") )       
       # bigger clique, not interest
       # 
       if (length(setAA) > Kclique){
          return (0)
       }
       
       # do they have enough neighbors
       setAB = sort(union(setAA,setBB))
       if (length(setAB)<Kclique ) {
          return (0)
       }

       # check if seA+setB is already identified as a cliuqe
       if (!is.null(globalCliques)) {
             #print( paste("Gloabl CLIQUE=", concatenate(as.character(cliquesFound[[1]]),","), 
             #       sep="") )
             #print( paste("SetAB=", concatenate(as.character(sort(setAB)),","), sep="") )
             #print( paste("SetAA=", concatenate(as.character(sort(setAA)),","), sep="") )
             #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )

             #found   = setInSets(setAB, globalCliques)
             found   = setInSets_matrix(setAB)

             if (found){
                #print("Already found (first)")
                return(1)
             }
       }
      
       #---------------------------------------------------------------------------------
       # continue the remaining combinations
       #
       # now we need bkclique=Kclique-length(setAA) nodes from setBB to form a k-clique 
       #  together with setAA, thus we need (bkclique-1) core from sedtBB
       #
       bkclique = Kclique-size.setAA

       nlinks = dim(newneighborlinks)[1]
       nnodes = length(setBB)
       ncc    = 2*nlinks/nnodes/(nnodes-1)

       if(F){
       #if(ncc <=tuo_coreCC  & nnodes >= tuo_coreSZ) {
       ret=kcore_subnets(linksindex=newneighborlinks, nodesindex=setBB,
                               kcore=bkclique-1 )
       setBB               = ret[[1]]
       newneighborlinks    = ret[[2]]
       neighborsfulllinked = ret[[3]]

       # if all remaining neighbors are fully linked, and they form a bkclique-clique,
       #   we simply return the union of AA+BB
       #  or we reject this clique because of the bigger/smaller size
       #
       if(neighborsfulllinked) {
                #print( paste("SetBBcore=", concatenate(as.character(sort(setBB)),","), sep="") )

                if ( length(setBB)==bkclique) {
                    sortedAB   = sort(union(setAA, setBB) )
                    currclique = list(sortedAB)
                    #print( paste("CLIQUEfull=", concatenate(as.character(currclique),","), sep="") )

                    #found    = setInSets(sortedAB, globalCliques)
                    found    = setInSets(sortedAB)

                    if(!found) {
                        #currCliques          = union(cliquesFound,  currclique)
                        #options(globalCliques= currCliques)
                        globalCliques       <<- c(globalCliques, currclique)

                        xcliq = globalNodeZeroIndex
                        xcliq[sortedAB] = 1
                        globalCliquesMatrix <<- rbind(globalCliquesMatrix, xcliq)

                        return (1)
                    }
                }
                #currclique = list(sort(union(setAA, setBB) ) )
                #print( paste("CLIQUEfull -- =", concatenate(as.character(currclique),","), sep="") )
                return (1)
       }     

       # the neighbors don't form a (bkclique-2)-core, so curnode and its neighbors cannot
       #  form a Kclique-clique
       if(is.null(setBB)) {
          #print( "no cores exisit to form a cilque")
          return (0)
       }

      }



       # try different combinations of the neighbors
       # as we need bkclique nodes from setBB, so the combinations with too few nodes (<bkclique )
       # will not be cosnsidered
       #
       total_cliques = 0
       bbsize        = length(setBB)

       # index of the last node to be moved into setA, so it still has enough nodes behind it 
       #  to possibly form a clique
       #
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # to avoid finding the same clique multiple times, the nodes have to be transfered 
       #  in ascendant/descendant order
       #
       # choose the last element to be fed into AA, if there are too few elements behind an element
       #  it won't be transfered into A
       #
       #  |        bbsize      |
       #  **********************
       #             | bkclique|
       #
       #             ^
       #             |
       #            endIdx

          endIdx        = (bbsize-bkclique) + 1
          if(endIdx<=0) {
              return(0)
          }

          #endIdx = bbsize -1

          for (m in c(1:endIdx) ){
             inode = setBB[m]
             irest = c( (m+1):bbsize)
             count = recursive_Cliques_nodenghr(curnode =inode, setA=setAA, setB=setBB[irest], 
                                    neighborlinks=newneighborlinks, 
                                    Kclique      =Kclique, coreidx=coreidx)
             total_cliques = total_cliques + count
          }


       return (total_cliques)
}

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------


recursive_Cliques_Super=function(curnode, setA, setB, adjlists, neighborlinks, Kclique=3){

       #cliquesFound = getOption("globalCliques")

       #since the current node don't have enough neighbors (>=Kclique)
       # so it cannot be in a clique of size = Kclique
       #

       #print (curnode)
       if( length(adjlists[[curnode]]) < Kclique-1){
           #print( paste("N(", as.character(curnode), ")", as.character(length(adjlists[[curnode]])), sep="") )
           return (0)
       }

       # move up v's neighbors to setA
       #
       setAA      = union(setA, curnode)
       size.setAA = length(setAA)

       if (is.null(setB)){
          setBB = adjlists[[curnode]]
       } else{
          setBB = setdiff(setB, curnode) # remove v from setB
          setBB = intersect(setB, adjlists[[curnode]])
       }

       #print( paste("SetAA=", concatenate(as.character( sort(setAA)),","), sep="") )
       #print( paste("SetBB=", concatenate(as.character( sort(setBB)),","), sep="") )
       
       # bigger clique, not interest
       # 
       if (length(setAA) > Kclique){
          return (0)
       }

       # do they have enough neighbors
          setAB = sort(union(c(setAA,setBB), NULL))
          #setAB = union(setAA,setBB) 

          if (length(setAB)<Kclique ) {
            return (0)
          }

          # check if seA+setB is already identified as a cliuqe
          if (!is.null(globalCliques)) {
             #print( paste("Gloabl CLIQUE=", concatenate(as.character(cliquesFound[[1]]),","), 
             #       sep="") )
             #print( paste("SetAB=", concatenate(as.character(sort(setAB)),","), sep="") )
             #print( paste("SetAA=", concatenate(as.character(sort(setAA)),","), sep="") )
             #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )

             found   = setInSets(setAB, globalCliques)
             if (found){
                #print("Already found (first)")
                return(0)
             }
          }

       if  (length(setBB)==1 ){
          currclique  = sort(union(setAA, setBB) )
          if( length(currclique) == Kclique ){
             currcliqueL          = list(currclique)

             #set global variable
             #
             #globalCliques       <<- union(globalCliques,  currcliqueL)
             globalCliques       <<- c(globalCliques,  currcliqueL)

             #options(globalCliques= currCliques)

             return (1)
          }
          return (0)

       }else if (length(setBB)==0 ){
          # record
          if( size.setAA == Kclique ){
             #newglobalCliques = union(globalCliques, list(sort(setAA)) )
             #print(list(sort(setAA)))
             #print("--:")
             
             currclique           = sort(setAA)
             currcliqueL          = list(currclique)

             #set global variable
             #
             #currCliques          = union(globalCliques,  currcliqueL)
             #options(globalCliques= currCliques)
             #globalCliques       <<- union(globalCliques,  currcliqueL)
             globalCliques       <<- c(globalCliques,  currcliqueL)

             return (1)
          }
          return (0)
       }
      
       #---------------------------------------------------------------------------------
       # continue the remaining combinations
       #
       # now we need bkclique=Kclique-length(setAA) nodes from setBB to form a k-clique 
       #  together with setAA, thus we need (bkclique-1) core from sedtBB
       #
       bkclique = Kclique-size.setAA

       # try different combinations of the neighbors
       # as we need bkclique nodes from setBB, so the combinations with too few nodes (<bkclique )
       # will not be cosnsidered
       #
       total_cliques = 0
       bbsize        = length(setBB)

       # index of the last node to be moved into setA, so it still has enough nodes behind it 
       #  to possibly form a clique
       #
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # to avoid finding the same clique multiple times, the nodes have to be transfered 
       #  in ascendant/descendant order
       #
       # choose the last element to be fed into AA, if there are too few elements behind an element
       #  it won't be transfered into A
       #
       #  |        bbsize      |
       #  **********************
       #             | bkclique|
       #
       #             ^
       #             |
       #            endIdx

          endIdx        = (bbsize-bkclique) + 1
          if(endIdx<=0) {
              return(0)
          }

       # special process
       #if (bkclique==2 & Kclique<=4){
       if (bkclique==2 & Kclique<= 5){
           blinks = getSubnetwork_LinkPairs(linkpairs=neighborlinks, subnetNodes=setBB)
           if( is.null(blinks) ){
               return (0)
           }

           no.blinks = dim(blinks)[1]
           qcliques  = cbind(blinks, rep(curnode,no.blinks))
           qcliquesO = apply(qcliques, 1, sort) #notice that the array is transposed after sort

           #print(paste("row,col=", dim(qcliques)[1],  dim(qcliques)[2],sep=" " ))
           #print(paste("row,col (O)=", dim(qcliquesO)[1],  dim(qcliquesO)[2],sep=" " ))

           #cliquesFound = NULL
           for (q in no.blinks){
             #cliquesFound = union(cliquesFound, list(qcliquesO[, q]) )
             found         = setInSets(list(qcliquesO[,q]), globalCliques)
             if(!found) {
                #print( paste("CLIQUEnormal2=", concatenate(as.character(currclique),","), sep="") )
                #cliquesFound = c(cliquesFound, list(qcliquesO[,q]) )
                
                globalCliques <<- c(globalCliques,  list(qcliquesO[,q]))
                
             }
           }
           #
           #globalCliques       <<- union(globalCliques,  cliquesFound)
           #options(globalCliques= cliquesFound)

       }else{
          for (m in c(1:endIdx) ){
             inode = setBB[m]
             irest = c( (m+1):bbsize)
             count = recursive_Cliques_Super(curnode =inode, setA=setAA, setB=setBB[irest], 
                                    adjlists=adjlists,
                                    neighborlinks=neighborlinks, 
                                    Kclique=Kclique)
             total_cliques = total_cliques + count
          }

       }

       return (total_cliques)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

recursive_Cliques=function(curnode, setA, setB, adjlists, neighborlinks, Kclique=3){

       #cliquesFound = getOption("globalCliques")
       #adjlists     = getOption("")

       #sytr= paste("A=",  concatenate(as.character(setA),","), 
       #            " B=", concatenate(as.character(setB),",") )
       #print(sytr)

       #since the current node don't have enough neighbors (>=Kclique)
       # so it cannot be in a clique of size = Kclique
       #

       #print (curnode)

       if( length(adjlists[[curnode]]) < Kclique-1){
           #print( paste("N(", as.character(curnode), ")", as.character(length(adjlists[[curnode]])), sep="") )
           return (0)
       }

       # move up v's neighbors to setA
       #
       setAA = union(setA, curnode)
       size.setAA = length(setAA)

       if (is.null(setB)){
          setBB = adjlists[[curnode]]
       } else{
          setBB = setdiff(setB, curnode) # remove v from setB
          setBB = intersect(setB, adjlists[[curnode]])
       }

       #print( paste("SetAA=", concatenate(as.character( sort(setAA)),","), sep="") )
       #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )

       # bigger clique, not interest
       #
       if (length(setAA) > Kclique){
          return (0)
       }

       # do they have enough neighbors
       if (!is.null(setA) | !is.null(setB)) {
          setAB = sort(union(c(setAA,setBB), NULL))
          if (length(setAB)<Kclique ) {
            return (0)
          }
          # check if seA+setB is already identified as a cliuqe
          if (!is.null(cliquesFound)) {
             #print( paste("Gloabl CLIQUE=", concatenate(as.character(cliquesFound[[1]]),","), 
             #       sep="") )
             #print( paste("SetAB=", concatenate(as.character(sort(setAB)),","), sep="") )
             #print( paste("SetAA=", concatenate(as.character(sort(setAA)),","), sep="") )
             #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )

             found   = setInSets(setAB, globalCliques)

             if (found){
                #print("Already found (first)")
                return(0)
             }
          }
       }
       
       if (length(setBB)==0 ){
          # record
          if( size.setAA == Kclique ){
             #newglobalCliques = union(globalCliques, list(sort(setAA)) )
             #print(list(sort(setAA)))
             #print("--:")
             
             sortedAA = sort(setAA)
             found    = setInSets(sortedAA, globalCliques)
             if(!found) {
               currclique = list(sortedAA)
               #print( paste("CLIQUEnormal1=", concatenate(as.character(currclique),","), sep="") )
               #currCliques          = union(cliquesFound,  currclique)
               #options(globalCliques= currCliques)

               globalCliques       <<- union(globalCliques,  currclique )

             }
             return (1)

          }
          return (0)

       } else if  (length(setBB)==1 ){
          currclique  = sort(union(setAA, setBB) )
          currcliqueL = list(currclique)

          if( length(currclique) == Kclique ){
             found    = setInSets(currclique, globalCliques)
             if(!found) {
                #print( paste("CLIQUEnormal2=", concatenate(as.character(currclique),","), sep="") )
                #currCliques          = union(cliquesFound,  currcliqueL)
                #options(globalCliques= currCliques)

                globalCliques       <<- union(globalCliques, currcliqueL)

                return (1)
             }
             return (1)
          }
          return (0)

       }

       #---------------------------------------------------------------------------------
       # continue the remaining combinations
       #
       # now we need bkclique=Kclique-length(setAA) nodes from setBB to form a k-clique 
       #  together with setAA, thus we need (bkclique-1) core from sedtBB
       #
       bkclique = Kclique-size.setAA
       ret=kcore_subnets(linksindex=neighborlinks, nodesindex=setBB,
                               kcore=bkclique-1 )

       setBB               = ret[[1]]
       neighborlinks2      = ret[[2]]
       neighborsfulllinked = ret[[3]]

       # if all remaining neighbors are fully linked, and they form a bkclique-clique,
       #   we simply return the union of AA+BB
       #  or we reject this clique because of the bigger/smaller size
       #
       if(neighborsfulllinked) {
                #print( paste("SetBBcore=", concatenate(as.character(sort(setBB)),","), sep="") )

                if ( length(setBB)==bkclique) {
                    sortedAB   = sort(union(setAA, setBB) )
                    currclique = list(sortedAB)
                    #print( paste("CLIQUEfull=", concatenate(as.character(currclique),","), sep="") )

                    found    = setInSets(sortedAB, cliquesFound)
                    if(!found) {
                        #currCliques          = union(cliquesFound,  currclique)
                        #options(globalCliques= currCliques)
                        globalCliques       <<- union(globalCliques, currclique)

                        return (1)
                    }
                }
                #currclique = list(sort(union(setAA, setBB) ) )
                #print( paste("CLIQUEfull -- =", concatenate(as.character(currclique),","), sep="") )
                return (1)
       }              

       # the neighbors don't form a (bkclique-2)-core, so curnode and its neighbors cannot
       #  form a Kclique-clique
       if(is.null(setBB)) {
          #print( "no cores exisit to form a cilque")
          return (0)
       }

       # try different combinations of the neighbors
       # as we need bkclique nodes from setBB, so the combinations with too few nodes (<bkclique )
       # will not be cosnsidered
       #
       currCliques   = NULL
       total_cliques = 0
       bbsize        = length(setBB)

       # index of the last node to be moved into setA, so it still has enough nodes behind it 
       #  to possibly form a clique
       #
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # to avoid finding the same clique multiple times, the nodes have to be transfered 
       #  in ascendant/descendant order
       #
       # choose the last element to be fed into AA, if there are too few elements behind an element
       #  it won't be transfered into A
       #
       #  |        bbsize       |
       #  **********************
       #             | bkclique |
       #
       #             ^
       #             |
       #            endIdx

       if(1==1) {
          endIdx        = (bbsize-bkclique) + 1
          if(endIdx<=0) {
              return(0)
          }

          for (m in c(1:endIdx) ){
             inode = setBB[m]
             irest = c( (m+1):bbsize)
             count = recursive_Cliques(curnode =inode, setA=setAA, setB=setBB[irest], 
                                    adjlists=adjlists,
                                    neighborlinks=neighborlinks2, 
                                    Kclique=Kclique)
             total_cliques = total_cliques + count
          }

       } else{
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          # try all possible combinations
          for (each in setBB){
             count = recursive_Cliques(curnode =each, setA=setAA, setB=setBB, 
                                    adjlists=adjlists,
                                    neighborlinks=neighborlinks2, 
                                    Kclique=Kclique)
             total_cliques = total_cliques + count
          }

       }
       return (total_cliques)
}



###############################################################################################################
############################### Start of normal clique identification #########################################
#

recursive_CliquesNormal=function(curnode, setA, setB, adjlists, neighborlinks, Kclique=3){

       #cliquesFound = getOption("globalCliques")

       #since the current node don't have enough neighbors (>=Kclique)
       # so it cannot be in a clique of size = Kclique
       #

       #print (curnode)
       if( length(adjlists[[curnode]]) < Kclique-1){
           #print( paste("N(", as.character(curnode), ")", as.character(length(adjlists[[curnode]])), sep="") )
           return (0)
       }

       # move up v's neighbors to setA
       #
       setAA = union(setA, curnode)
       size.setAA = length(setAA)

       if (is.null(setB)){
          setBB = adjlists[[curnode]]
       } else{
          setBB = setdiff(setB, curnode) # remove v from setB
          setBB = intersect(setB, adjlists[[curnode]])
       }

       #print( paste("SetAA=", concatenate(as.character( sort(setAA)),","), sep="") )
       #print( paste("SetBB=", concatenate(as.character( sort(setBB)),","), sep="") )

       # bigger clique, not interest
       # 
       if (length(setAA) > Kclique){
          return (0)
       }

       # do they have enough neighbors
          setAB = sort(union(c(setAA,setBB), NULL))
          if (length(setAB)<Kclique ) {
            return (0)
          }
          # check if seA+setB is already identified as a cliuqe
          if (!is.null(globalCliques)) {
             #print( paste("Gloabl CLIQUE=", concatenate(as.character(cliquesFound[[1]]),","), 
             #       sep="") )
             #print( paste("SetAB=", concatenate(as.character(sort(setAB)),","), sep="") )
             #print( paste("SetAA=", concatenate(as.character(sort(setAA)),","), sep="") )
             #print( paste("SetBB=", concatenate(as.character(sort(setBB)),","), sep="") )
             found   = setInSets(setAB, globalCliques)
             if (found){
                #print("Already found (first)")
                return(0)
             }
          }

       if (length(setBB)==0 ){
          # record
          if( size.setAA == Kclique ){
             #newglobalCliques = union(globalCliques, list(sort(setAA)) )
             #print(list(sort(setAA)))
             #print("--:")
             
             currclique           = sort(setAA)
             currcliqueL          = list(currclique)
             found   = setInSets(currclique, globalCliques)
             if(!found) {
                 globalCliques        <<- union(globalCliques,  currcliqueL)
                 #options(globalCliques= currCliques)
             }
             return (1)
          }
          return (0)

       } else if  (length(setBB)==1 ){
          currclique  = sort(union(setAA, setBB) )
          if( length(currclique) == Kclique ){
             currcliqueL          = list(currclique)
             found   = setInSets(currclique, globalCliques)
             if(!found) {
                globalCliques        <<- union(globalCliques,  currcliqueL)
                #currCliques          = union(cliquesFound,  currcliqueL)
                #options(globalCliques= currCliques)
             }
             return (1)
          }
          return (0)
       }
       
       #---------------------------------------------------------------------------------
       # continue the remaining combinations
       #
       # now we need bkclique=Kclique-length(setAA) nodes from setBB to form a k-clique 
       #  together with setAA, thus we need (bkclique-1) core from sedtBB
       #
       bkclique = Kclique-size.setAA

       # try different combinations of the neighbors
       # as we need bkclique nodes from setBB, so the combinations with too few nodes (<bkclique )
       # will not be cosnsidered
       #
       currCliques   = NULL
       total_cliques = 0
       bbsize        = length(setBB)

       # index of the last node to be moved into setA, so it still has enough nodes behind it 
       #  to possibly form a clique
       #
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # to avoid finding the same clique multiple times, the nodes have to be transfered 
       #  in ascendant/descendant order
       #
       # choose the last element to be fed into AA, if there are too few elements behind an element
       #  it won't be transfered into A
       #
       #  |        bbsize      |
       #  **********************
       #             | bkclique|
       #
       #             ^
       #             |
       #            endIdx

          endIdx        = (bbsize-bkclique) + 1
          if(endIdx<=0) {
              return(0)
          }

          for (m in c(1:endIdx) ){
             inode = setBB[m]
             irest = c( (m+1):bbsize)
             count = recursive_CliquesNormal(curnode =inode, setA=setAA, setB=setBB[irest], 
                                    adjlists=adjlists,
                                    neighborlinks=neighborlinks, 
                                    Kclique=Kclique)
             total_cliques = total_cliques + count
          }

       return (total_cliques)
}


# used global variables
#      selectedNodesName
#      selectedLinksIndex
#      degrees
#      no.nodes
#      no.links


############################################################################################################
#
#                                  Super-fast clique identification
# used global variables
#      selectedNodesName
#      selectedLinksIndex
#      degrees
#      no.nodes
#      no.links
#
#  global_nodenghr_cores: each element corresponding to a node is a list of following items 
#     for the node's neighbors:
#           finalList = list(corenesses= allcores, 
#                        coresnodes= kcoreMemberNodes, 
#                        coreslinks= kcoreMemberLinks,
#                        coressizes= kcoresSizes, 
#                        iscliques = kcoreFullLinked)
#    
#  global_nodenghr_minmaxcore: no.nodes x 2 matrix with (min coreness, max coreness)
#
#
# IF MemoyForSpeed= T, then the core identification will return all nodes and links
#    in each core so that there is no need to perform merge operation to get link members
#    in each core
#

find_kcliques = function (min_clique=3, returnNodeCliques=F, fkey=NULL, printout=F,
                          transMatrix=NULL, MemoyForSpeed=T, cliqueOnly=T, use_link_id=F) {
    DEBUG=F

    now<-proc.time()[3]

    # make output files
    #
    imgKclique         = paste(fkey, "_imgKcliques_histogram",       sep='')
    imgKcommunitySize  = paste(fkey, "_imgKcommunitySize_histogram", sep='')
    imgKcommunityCount = paste(fkey, "_imgKcommunityCount_histogram",sep='')

    logFname       = paste(fkey, "_log.xls",       sep='')
    fcliquesOnly   = paste(fkey, ".xls",           sep='') #cliques
    fcommunityOnly = paste(fkey, "_kcommunity-nodeMembership.xls",sep='') #k-community
    fcommunityLinks= paste(fkey, "_kcommunity-linkmembership.xls",sep='') #k-community

    linkpairsIndexed <<- cbind(selectedLinksIndex, c(1:no.links) )
    colnames(linkpairsIndexed) <<-c("src","dst","id")

    nulcoresngbr = list(corenesses= NULL,coresnodes= NULL, 
                        coreslinks= NULL, coressizes= NULL, iscliques = NULL)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # find all cores in each node's neighbors
    #
    # k-core is defined as  number of links for each node is >= k
    # so a k+1 clique is contained in (k-1)-core
    #
    nodes_with_cores = KCores_ALL_in_neighbors(minkcore=min_clique-1, returnLinkID=use_link_id)

    allcores.readTime <- proc.time()[3]-now

    # as the coreness is based on node's neighbor, so we need add the node in
    #  
    maxcoreness  = max(global_nodenghr_minmaxcore[,2]) + 1

    coreindex    = c( (min_clique-1):maxcoreness)
    cliqueindex  = coreindex   + 1 
    max_clique   = maxcoreness + 1

    maxcliqueVect= rep(max_clique, no.nodes)
    no.cores     = length(coreindex)

    allSelnodesIndex = c(1:no.nodes)

    #-----------------------------  1) finding all possible cliques ----------------------------
    #
    # explore all possible combinations by a recursive function
    #
    print("1. finding all possible cliques")

    iclqFindTimes = NULL
    for (i in c(no.cores:1) ) {

        iclq= cliqueindex[i]
        iclq.now<-proc.time()[3]

        # here coreness (global_nodenghr_minmaxcore) is based on a node's neighbor, 
        # so including the node itself leads to  (global_nodenghr_minmaxcore + 1)-core
        #
        iclqNodesSel   = global_nodenghr_minmaxcore[,2] + 1 >= iclq-1

        iclqNodesIndex = allSelnodesIndex[iclqNodesSel]
        no.nodesINiclq = length(iclqNodesIndex)

        if(no.nodesINiclq < iclq){
           # recording time
           iclq.readTime <- proc.time()[3]-iclq.now
           iclqFindTimes =rbind(iclqFindTimes, c(iclq, no.nodesINiclq, 
                            global_nolinks_nghbr_cores[i],iclq.readTime) )
           next
        }

        # iterate each node to identify all iclq-cliques including this node
        #
        jjindex = c(1:no.nodesINiclq)
        for ( j in jjindex ) {

           cliques_count = 0
           v = iclqNodesIndex[j]
           
           mystr = paste("Clique-", as.character(iclq), ": ", as.character(j), "/", 
                             as.character(no.nodesINiclq), sep="")
           print(mystr)


           # global_nodenghr_cores[[]]$coresnesses
           # global_nodenghr_cores[[]]$coresnodes
           # global_nodenghr_cores[[]]$coreslinks
           # global_nodenghr_cores[[]]$coressizes
           # global_nodenghr_cores[[]]$iscliques
           #neighbors    = global_nodenghr_cores[[v]]$coresnodes[[i]]

          # first call, whether its neighbors form a clique already has been checked
          #
          # when the current node has Kclique-1 neighbors, there are two scenarios:
          #  1) they form a (Kclique-1)-clique
          #  2) they don't form a (Kclique-1)-clique
          #
          if(global_nodenghr_cores[[v]]$coressizes[i] == iclq-1) {
             # 1)
             if(global_nodenghr_cores[[v]]$iscliques[i]){
                 currclique  = sort(c(v, global_nodenghr_cores[[v]]$coresnodes[[i]]) )
                 currcliqueL = list(currclique)

                 #set global variable
                 #
                 #globalCliques  <<- union(globalCliques,  currcliqueL)
                 #found   = setInSets(currclique, globalCliques)
                 found   = setInSets_matrix(currclique)

                 if(!found){
                    globalCliques <<- c(globalCliques,  currcliqueL)

                    xcliq = globalNodeZeroIndex
                    xcliq[currclique] = 1
                    globalCliquesMatrix   <<- rbind(globalCliquesMatrix, xcliq)

                 }
                 cliques_count = 1
             }
             cliques_count = 1

          } else if(global_nodenghr_cores[[v]]$coressizes[i] < iclq-1) {
             # if the currrent node does't have enough neighbors, it should be removed
             #  from the neighbors of each other node
             # Notice that global_nodenghr_cores[[v]]$coressizes is dynamically changed
             #
             cliques_count = 1
          } else{

             if(use_link_id) {
                    ijlinks = selectedLinksIndex[ global_nodenghr_cores[[v]]$coreslinks[[i]], ]
                    cliques_count = recursive_Cliques_nodenghr(curnode=v, setA=NULL, setB=NULL, 
                              neighborlinks= ijlinks,
                              Kclique=iclq, coreidx=i, use_link_id=use_link_id)#-1)
             }else{
                    cliques_count = recursive_Cliques_nodenghr(curnode=v, setA=NULL, setB=NULL, 
                              neighborlinks= global_nodenghr_cores[[v]]$coreslinks[[i]],
                              Kclique=iclq, coreidx=i, use_link_id=use_link_id)#-1)
             }

          }

           #cliques_count

           #-------------------- important -----------------------
           # i) if no cliques found, keep the node and its links
           #
           if (cliques_count ==0){
               next
           }


           # remove the current node and its links
           #
           # 1) reset the i-th core of the current node v
           #
           if(DEBUG){
           global_nodenghr_cores[[v]]$corenesses[i]    <- 0
           global_nodenghr_cores[[v]]$coresnodes[[i]]  <- NA
           global_nodenghr_cores[[v]]$coreslinks[[i]]  <- NA
           global_nodenghr_cores[[v]]$coressizes[i]    <- 0
           global_nodenghr_cores[[v]]$iscliques[i]     <- F
           }else{
           global_nodenghr_cores[[v]]$corenesses[i]    <<- 0
           global_nodenghr_cores[[v]]$coresnodes[[i]]  <<- NA
           global_nodenghr_cores[[v]]$coreslinks[[i]]  <<- NA
           global_nodenghr_cores[[v]]$coressizes[i]    <<- 0
           global_nodenghr_cores[[v]]$iscliques[i]     <<- F
           }

           # 2) remove the current node from LinkPairs and AdjList of
           #     the i-th core of each other node (don't include the current node v)
           #    so that the same cliques including v will not be found again
           #
           for ( m in jjindex[-j] ) {
                qq = iclqNodesIndex[m]
                isneighbor = is.element(v, global_nodenghr_cores[[qq]]$coresnodes[[i]])
                if(!isneighbor){
                    next
                }

                icorenew =update_nodenghr_cores(nodenghr_cores=global_nodenghr_cores[[qq]], 
                                                removegene=v, coreidx=i, use_link_id=use_link_id)

                # qq's neighbors don't form a core or v is not in qq's neighbors
                if(is.null(icorenew)){
                   next
                }


                # elements of List can be null, but not for vector
                #
                # IMPORTANT: update list entries with NA instead of NULL,
                #
                if (is.null(icorenew$corenesses)){
                   if(DEBUG){
                   global_nodenghr_cores[[qq]]$corenesses[i]   <- 0
                   global_nodenghr_cores[[qq]]$coressizes[i]   <- 0
                   global_nodenghr_cores[[qq]]$iscliques[i]    <- F
                   global_nodenghr_cores[[qq]]$coresnodes[[i]] <- NA
                   global_nodenghr_cores[[qq]]$coreslinks[[i]] <- NA
                   }else{
                   global_nodenghr_cores[[qq]]$corenesses[i]   <<- 0
                   global_nodenghr_cores[[qq]]$coressizes[i]   <<- 0
                   global_nodenghr_cores[[qq]]$iscliques[i]    <<- F
                   global_nodenghr_cores[[qq]]$coresnodes[[i]] <<- NA
                   global_nodenghr_cores[[qq]]$coreslinks[[i]] <<- NA
                   }

                }else{
                   if(DEBUG){
                   global_nodenghr_cores[[qq]]$corenesses[i]   <- icorenew$corenesses[1]
                   global_nodenghr_cores[[qq]]$coressizes[i]   <- icorenew$coressizes[1]
                   global_nodenghr_cores[[qq]]$iscliques[i]    <- icorenew$iscliques[1]
                   global_nodenghr_cores[[qq]]$coresnodes[[i]] <- icorenew$coresnodes[[1]]
                   global_nodenghr_cores[[qq]]$coreslinks[[i]] <- icorenew$coreslinks[[1]]
                   }else{
                   global_nodenghr_cores[[qq]]$corenesses[i]   <<- icorenew$corenesses[1]
                   global_nodenghr_cores[[qq]]$coressizes[i]   <<- icorenew$coressizes[1]
                   global_nodenghr_cores[[qq]]$iscliques[i]    <<- icorenew$iscliques[1]
                   global_nodenghr_cores[[qq]]$coresnodes[[i]] <<- icorenew$coresnodes[[1]]
                   global_nodenghr_cores[[qq]]$coreslinks[[i]] <<- icorenew$coreslinks[[1]]
                   }
                }
           }#qq


        }#for (j in c(

       # recording time
       iclq.readTime <- proc.time()[3]-iclq.now
       iclqFindTimes =rbind(iclqFindTimes, c(iclq, no.nodesINiclq, 
                            global_nolinks_nghbr_cores[i],iclq.readTime) )

    }#for (iclq in

    clique.readTime <- proc.time()[3]-now

    # save execution time
    #
    colnames(iclqFindTimes ) <- c("k-clique", "(k-1)-core nodes",  "(k-1)-core links", "time")
    write.table(iclqFindTimes, logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

    minutes = as.character( round(allcores.readTime/60.0, 3) )
    mystr = paste("\nTime for identifying all cores = ", as.character(allcores.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("\nTime for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    #----------------------------- 2) plot distribution of cliques ----------------------------
    #
    #globalCliques = getOption("globalCliques")

    no.kcliques= length(globalCliques)
    kcliquesSize = rep(0, no.kcliques)
    for(i in c(1:no.kcliques )){
       kcliquesSize[i] = length(globalCliques[[i]])
    }
    maxk=max(kcliquesSize)
    mink=min(kcliquesSize)

    cfmatrix   = histogram_4integers(ivector=kcliquesSize, fkeyname=imgKclique, keyword="k-cliques")
    write.table(cfmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)    


    #-------------------------------- 3) print out k-cliques -----------------------------------
    #
    print("3. Output all possible cliques")

    write.table(rbind(c("k-clique","member")), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    for (i in c(1:no.kcliques) ){
        icliques = paste(as.character( length(globalCliques[[i]]) ), "\t", sep="")
        for (jnode in globalCliques[[i]]){
            icliques = paste(icliques, selectedNodesName[jnode],", ", sep="")
        }
        #print(icliques)
        write.table(rbind(icliques), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }
    
    # return if clique only
    if(cliqueOnly) {
        return (clique.readTime)
    }


    # ======================== 4) prepare clique adjacency matrix ================================
    #
    print("4. prepare clique adjacency matrix")

    cliqueAdjMatrix = matrix(0,no.kcliques, no.kcliques)
    diag(cliqueAdjMatrix) <- kcliquesSize

    for (i in c(1:(no.kcliques-1)) ){
       for (j in c((i+1):no.kcliques) ){
          ijoverlap            = intersect(globalCliques[[i]], globalCliques[[j]])
          cliqueAdjMatrix[i,j] = length(ijoverlap)
          cliqueAdjMatrix[j,i] = length(ijoverlap)
       }
    }

    #++++++++++++++++++++++++ initialization ++++++++++++++++++++++++++++++++++++
    #
    kcommunities       = cbind(selectedNodesName)
    #kcommunitiesLinks = cbind(linkpairs)
    leftnodeNames      = selectedNodesName[ selectedLinksIndex[,1] ]
    rightnodeNames     = selectedNodesName[ selectedLinksIndex[,2] ]
    kcommunitiesLinks  = cbind(leftnodeNames, rightnodeNames)

    kcommunityTitle          = c("gene")
    kcommunitySize           = NULL
    kcommunityComponentCount = NULL
    kcommunityComponentCount2= NULL

    # ======================== 5) merge cliques into communities ================================
    #
    for (kc in c(maxk:mink) ){

         kcname = paste(as.character(kc), "-communities", sep="")

         print(paste("5... identify ", kcname) )

         # ~~~~~~~~~~~~  clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~
         #

         # <1> overlap size
         kcCliqueAdjMatrix = ifelse(cliqueAdjMatrix == kc-1, 1, 0) #off diagonal

         # <2> clique size
         selDiag           = kcliquesSize == kc

         #nselIndex         = c(1:no.kcliques)[!selDiag]
         #kcCliqueAdjMatrix[nselIndex,] = 0
         #kcCliqueAdjMatrix[,nselIndex] = 0
         #diag(kcCliqueAdjMatrix) = selDiag

         # <3> find the sub matrix with connected kcliques
         #
         diag(kcCliqueAdjMatrix) <- 0
         overlapvect = apply(kcCliqueAdjMatrix, 1, sum)

         selkcs      = overlapvect >=1
         selkcs      = (selkcs & selDiag) | selDiag # add in those isolated cliques

         if ( sum(selkcs)==0 ){
             next
         }

         #
         #
         # ~~~~~~~~~~~~  END of clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~

         kcommunityMatrix = kcCliqueAdjMatrix[selkcs, selkcs]
         kcindex          = c(1:no.kcliques)[selkcs]
         

         # ~~~~~~~~~~~  find connected components in k-cliques ~~~~~~~~~~~~~~~~~~~~
         #
         kc_CClabel = find_ConnectedComponents_adjmatrix(adjmatrix=kcommunityMatrix, 
                                  minModuleSize=0, cclabelstart=1)

         # merge cliques in each CC
         CCsize        = table(kc_CClabel[,2])
         CCs           = names(CCsize)
         no.components = length(CCs)

         kcommunityComponentCount = c(kcommunityComponentCount, rep(kc, no.components) )
         kcommunityComponentCount2= c(kcommunityComponentCount2,no.components)

         CCsizes = NULL
         for (cc in CCs){ # cc is like "0001","0002"

            # make title
            kcommCCgeneIndices = NULL
            kcommCCname        = paste(kcname, "-cc", cc,sep="")
            kcommunityTitle    = c(kcommunityTitle, kcommCCname)

            # get members of the current connected components in k-community
            #
            icc  = as.integer(cc)
            isel = cc==kc_CClabel[,2]
            cliquesIdx = kcindex[isel]
            mymembers = NULL
            for (m in cliquesIdx){
                 mymembers = union(mymembers, globalCliques[[m]] )
            }
            mymembers  = as.integer(mymembers) # index
            mvect      = rep(0, no.nodes)
            mvect[mymembers] = 1
 
            kcommunities          = cbind(kcommunities, mvect)         
            CCsizes               = c( CCsizes, length(mymembers) )

            kcommunitySize        = c(kcommunitySize, length(mymembers) )

            # network view
            imgKcommunity = paste(fkey, "_", kcommCCname, ".png", sep="")        
            kcLinksIndex  = getSubnetwork_LinkPairs(linkpairsIndexed,subnetNodes= mymembers)
            kcLinksIndex  = as.matrix(kcLinksIndex)
            kcLinks       = kcLinksIndex[,c(1:2)]
            kcPairsIndex  = as.integer(kcLinksIndex[,3])

            # from index to actual node name
            kcLinksByName = cbind(selectedNodesName[kcLinks[,1]], selectedNodesName[kcLinks[,2]])

            plotNetwork_inPairs_DisplayGeneSymbol(linkpairs=kcLinksByName, 
                           directed=F, geneinfo=transMatrix, genesymbolCol=2,
                           nodecolor="blue", nodeshape=30, nodenamecolor="purple",
                           labelpos=1, rankhubby="totallinks",
                           disphubs=0, plotfigure=T, 
                           fimg=imgKcommunity, saveDegrees=F)


            # output k-community links membership
            #
            mlinkvect = rep(0, no.links)
            mlinkvect[kcPairsIndex] = 1
            kcommunitiesLinks = cbind(kcommunitiesLinks, mlinkvect)
         }

    }

    # ======================== 6) output node/link membership & statistics ==============================
    #
    print("6. plot histogram of number of components in k-community")

    newtitle      = paste("numbers of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunityComponentCount, 
                                  fkeyname=imgKcommunityCount, 
                                  keyword=newtitle, hxlabel="k")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    newtitle      = paste("sizes of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunitySize, 
                                  fkeyname=imgKcommunitySize, 
                                  keyword=newtitle, hxlabel="community size")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("Time for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)", sep="")
    appendStringToFile(logFname, mystr, newline=F)

    # ======================== 7) output node/link membership & statistics ==============================
    #
    print("7. output node/link membership & statistics")

    # save the number of connected components in each k-community
    #
    #countMatrix = rbind(as.character(maxk:mink), kcommunityComponentCount2)
    #countMatrix = cbind(c("k-community","components"),countMatrix)
    #appendStringToFile(logFname, "\n", newline=T)
    #write.table(countMatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    # community-component NODE membership ony
    write.table(rbind(kcommunityTitle), fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunities, fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)
    
    # community-component LINKS membership ony
    kclinktitle = c("src", "dst", kcommunityTitle[-1])
    write.table(rbind(kclinktitle),fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunitiesLinks, fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

find_kcliquesNormal = function (min_clique=3, returnNodeCliques=F, fkey=NULL, printout=F,
                          transMatrix=NULL, cliqueOnly=T) {

    now<-proc.time()[3]

    # make output files
    #
    imgKclique         = paste(fkey, "_imgKcliques_histogram",sep='')
    imgKcommunitySize  = paste(fkey, "_imgKcommunitySize_histogram",sep='')
    imgKcommunityCount = paste(fkey, "_imgKcommunityCount_histogram",sep='')

    logFname       = paste(fkey, "_log.xls",       sep='')
    fcliquesOnly   = paste(fkey, ".xls",           sep='') #cliques
    fcommunityOnly = paste(fkey, "_kcommunity-nodeMembership.xls",sep='') #k-community
    fcommunityLinks= paste(fkey, "_kcommunity-linkmembership.xls",sep='') #k-community

    linkpairsIndexed = cbind(selectedLinksIndex, c(1:no.links) )

    # find the maximum coreness, which will be used for inferring max(k)-clique
    #
    #maxcoreness  = find_maxcoreness(selectedLinksIndex)
    #max_clique   = maxcoreness + 1

    # max clique
    #
    degrees       = degree_ByLinkPairs (linkpairsIndexed, directed=F)
    maxdegree     = max(degrees$totallinks)

    maxcliqueVect = rep(maxdegree, no.nodes)    
    max_clique    = maxdegree

    #globalAdjLists  = makeAjacencyListsFromLinksIndex(selectedLinksIndex)

    adjlist.readTime <- proc.time()[3]-now


    #-----------------------------  1) finding all possible cliques ----------------------------
    #
    # explore all possible combinations by a recursive function
    #
    print("1. finding all possible cliques")

    iclqFindTimes = NULL
    for (iclq in c(max_clique:min_clique) ) {
        
        iclq.now<-proc.time()[3]

        # for normal clique identification, we use all links & nodes
        #
        iclqNodesIndex     = c(1:no.nodes) #selected nodes (index)
        iclqLinkpairsIndex = linkpairsIndexed       #selected links (indices)

        # the index of the lists correspond to the original node index
        #
        #iclqAdjlists   = makeAjacencyListsFromLinksIndex(linkpairsIndexed)
        iclqAdjlists   = globalAdjLists

        no.linksINiclq     = dim(linkpairsIndexed)[1]
        no.nodesINiclq     = no.nodes        # iterate each node to identify all iclq-cliques including this node

        #
        for ( j in c(1:no.nodesINiclq) ) {
           cliques_count = 0
           v = iclqNodesIndex[j]
           
           if ( length(iclqAdjlists[[v]])  >= iclq-1 ){

               mystr = paste("Clique-", as.character(iclq), ": ", as.character(j), "/", 
                             as.character(no.nodesINiclq), sep="")
               print(mystr)


               cliques_count = recursive_CliquesNormal(curnode=v, setA=NULL, setB=NULL,
                                         adjlists     =iclqAdjlists,
                                         neighborlinks=iclqLinkpairsIndex, Kclique=iclq)#-1)
                  
              #cliquesFoundAfter = getOption("globalCliques")

           }#if ( length(iclqAdjlists[[v]]) >= iclq-1

           #-------------------- important -----------------------
           # i) if no cliques found, keep the node and its links
           #
           if (cliques_count ==0){
               next
           }

           # ii) remove the current node and its links
           #
           # ii.1) remove it from LinkPairs
           #
           selLeft    = iclqLinkpairsIndex[,1] != v
           selRight   = iclqLinkpairsIndex[,2] != v
           sel        = selLeft & selRight
           no.remains = sum(sel)

           # no enough links left, stop the current clique finding
           #
           if ( no.remains < iclq-1 ) {
              break
           } else {
              if(no.remains ==1) {
                 iclqLinkpairsIndex = rbind(iclqLinkpairsIndex[sel,])
              } else{
                 iclqLinkpairsIndex = iclqLinkpairsIndex[sel,]
              }
           }

           # ii.2) remove it from AdjList
           #
           iclqAdjlists[[v]] = -1           
           for ( m in c(1:no.nodesINiclq) ) {
                qq = iclqNodesIndex[m]
                mneighbors = setdiff(iclqAdjlists[[qq]], v)
                if (length(mneighbors)==0) {
                     iclqAdjlists[[qq]] = -1
                }else{
                     iclqAdjlists[[qq]] = mneighbors
                }
            }

            #iclqAdjlists
            #iclqLinkpairsIndex 

        }#for (j in c(

       # recording time
       iclq.readTime <- proc.time()[3]-iclq.now
       iclqFindTimes =rbind(iclqFindTimes, c(iclq, no.nodesINiclq, no.linksINiclq, iclq.readTime) )

    }#for (iclq in

    clique.readTime <- proc.time()[3]-now

    # save execution time
    #
    colnames(iclqFindTimes ) <- c("k-clique", "(k-1)-core nodes", "(k-1)-core links", "time")
    write.table(iclqFindTimes, logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

    minutes = as.character( round(adjlist.readTime /60.0, 3) )
    mystr = paste("\nTime for identifying all cores = ", as.character(adjlist.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("\nTime for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    #----------------------------- 2) plot distribution of cliques ----------------------------
    #
    #globalCliques = getOption("globalCliques")
    no.kcliques= length(globalCliques)
    kcliquesSize = rep(0, no.kcliques)
    for(i in c(1:no.kcliques )){
       kcliquesSize[i] = length(globalCliques[[i]])
    }
    maxk=max(kcliquesSize)
    mink=min(kcliquesSize)

    cfmatrix   = histogram_4integers(ivector=kcliquesSize, fkeyname=imgKclique, keyword="k-cliques")
    write.table(cfmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)


    #-------------------------------- 3) print out k-cliques -----------------------------------
    #
    print("3. Output all possible cliques")

    write.table(rbind(c("k-clique","member")), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    for (i in c(1:no.kcliques) ){
        icliques = paste(as.character( length(globalCliques[[i]]) ), "\t", sep="")
        for (jnode in globalCliques[[i]]){
            icliques = paste(icliques, selectedNodesName[jnode],", ", sep="")
        }
        #print(icliques)
        write.table(rbind(icliques), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }

    # return if clique only
    if(cliqueOnly) {
        return (1)
    }
        
    # ======================== 4) prepare clique adjacency matrix ================================
    #
    print("4. prepare clique adjacency matrix")

    cliqueAdjMatrix = matrix(0,no.kcliques, no.kcliques)
    diag(cliqueAdjMatrix) <- kcliquesSize

    for (i in c(1:(no.kcliques-1)) ){
       for (j in c((i+1):no.kcliques) ){
          ijoverlap            = intersect(globalCliques[[i]], globalCliques[[j]])
          cliqueAdjMatrix[i,j] = length(ijoverlap)
          cliqueAdjMatrix[j,i] = length(ijoverlap)
       }
    }

    #++++++++++++++++++++++++ initialization ++++++++++++++++++++++++++++++++++++
    #
    kcommunities       = cbind(selectedNodesName)
    #kcommunitiesLinks = cbind(linkpairs)
    leftnodeNames      = selectedNodesName[ selectedLinksIndex[,1] ]
    rightnodeNames     = selectedNodesName[ selectedLinksIndex[,2] ]
    kcommunitiesLinks  = cbind(leftnodeNames, rightnodeNames)

    kcommunityTitle          = c("gene")
    kcommunitySize           = NULL
    kcommunityComponentCount = NULL
    kcommunityComponentCount2= NULL

    # ======================== 5) merge cliques into communities ================================
    #
    for (kc in c(maxk:mink) ){

         kcname = paste(as.character(kc), "-communities", sep="")

         print(paste("5... identify ", kcname) )

         # ~~~~~~~~~~~~  clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~
         #

         # <1> overlap size
         kcCliqueAdjMatrix = ifelse(cliqueAdjMatrix == kc-1, 1, 0) #off diagonal

         # <2> clique size
         selDiag           = kcliquesSize == kc

         #nselIndex         = c(1:no.kcliques)[!selDiag]
         #kcCliqueAdjMatrix[nselIndex,] = 0
         #kcCliqueAdjMatrix[,nselIndex] = 0
         #diag(kcCliqueAdjMatrix) = selDiag

         # <3> find the sub matrix with connected kcliques
         #
         diag(kcCliqueAdjMatrix) <- 0
         overlapvect = apply(kcCliqueAdjMatrix, 1, sum)

         selkcs      = overlapvect >=1
         selkcs      = (selkcs & selDiag) | selDiag # add in those isolated cliques

         if ( sum(selkcs)==0 ){
             next
         }

         #
         #
         # ~~~~~~~~~~~~  END of clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~

         kcommunityMatrix = kcCliqueAdjMatrix[selkcs, selkcs]
         kcindex          = c(1:no.kcliques)[selkcs]
         

         # ~~~~~~~~~~~  find connected components in k-cliques ~~~~~~~~~~~~~~~~~~~~
         #
         kc_CClabel = find_ConnectedComponents_adjmatrix(adjmatrix=kcommunityMatrix, 
                                  minModuleSize=0, cclabelstart=1)

         # merge cliques in each CC
         CCsize        = table(kc_CClabel[,2])
         CCs           = names(CCsize)
         no.components = length(CCs)

         kcommunityComponentCount = c(kcommunityComponentCount, rep(kc, no.components) )
         kcommunityComponentCount2= c(kcommunityComponentCount2,no.components)

         CCsizes = NULL
         for (cc in CCs){ # cc is like "0001","0002"

            # make title
            kcommCCgeneIndices = NULL
            kcommCCname        = paste(kcname, "-cc", cc,sep="")
            kcommunityTitle    = c(kcommunityTitle, kcommCCname)

            # get members of the current connected components in k-community
            #
            icc  = as.integer(cc)
            isel = cc==kc_CClabel[,2]
            cliquesIdx = kcindex[isel]
            mymembers = NULL
            for (m in cliquesIdx){
                 mymembers = union(mymembers, globalCliques[[m]] )
            }
            mymembers  = as.integer(mymembers) # index
            mvect      = rep(0, no.nodes)
            mvect[mymembers] = 1
 
            kcommunities          = cbind(kcommunities, mvect)         
            CCsizes               = c( CCsizes, length(mymembers) )

            kcommunitySize        = c(kcommunitySize, length(mymembers) )

            # network view
            imgKcommunity = paste(fkey, "_", kcommCCname, ".png", sep="")        
            kcLinksIndex  = getSubnetwork_LinkPairs(linkpairsIndexed,subnetNodes= mymembers)
            kcLinksIndex  = as.matrix(kcLinksIndex)
            kcLinks       = kcLinksIndex[,c(1:2)]
            kcPairsIndex  = as.integer(kcLinksIndex[,3])

            # from index to actual node name
            kcLinksByName = cbind(selectedNodesName[kcLinks[,1]], selectedNodesName[kcLinks[,2]])

            plotNetwork_inPairs_DisplayGeneSymbol(linkpairs=kcLinksByName, 
                           directed=F, geneinfo=transMatrix, genesymbolCol=2,
                           nodecolor="blue", nodeshape=30, nodenamecolor="purple",
                           labelpos=1, rankhubby="totallinks",
                           disphubs=0, plotfigure=T, 
                           fimg=imgKcommunity, saveDegrees=F)


            # output k-community links membership
            #
            mlinkvect = rep(0, no.links)
            mlinkvect[kcPairsIndex] = 1
            kcommunitiesLinks = cbind(kcommunitiesLinks, mlinkvect)
         }

    }

    # ======================== 6) output node/link membership & statistics ==============================
    #
    print("6. plot histogram of number of components in k-community")

    newtitle      = paste("numbers of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunityComponentCount, 
                                  fkeyname=imgKcommunityCount, 
                                  keyword=newtitle, hxlabel="k")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    newtitle      = paste("sizes of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunitySize, 
                                  fkeyname=imgKcommunitySize, 
                                  keyword=newtitle, hxlabel="community size")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("Time for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)", sep="")
    appendStringToFile(logFname, mystr, newline=F)

    # ======================== 7) output node/link membership & statistics ==============================
    #
    print("7. output node/link membership & statistics")

    # save the number of connected components in each k-community
    #
    #countMatrix = rbind(as.character(maxk:mink), kcommunityComponentCount2)
    #countMatrix = cbind(c("k-community","components"),countMatrix)
    #appendStringToFile(logFname, "\n", newline=T)
    #write.table(countMatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    # community-component NODE membership ony
    write.table(rbind(kcommunityTitle), fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunities, fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)
    
    # community-component LINKS membership ony
    kclinktitle = c("src", "dst", kcommunityTitle[-1])
    write.table(rbind(kclinktitle),fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunitiesLinks, fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    return (globalCliques)
}

find_kcliques_use_globalCORES = function (min_clique=3, returnNodeCliques=F, fkey=NULL, printout=F,
                          transMatrix=NULL, MemoyForSpeed=T, cliqueOnly=T) {

    now<-proc.time()[3]

    # make output files
    #
    imgKclique         = paste(fkey, "_imgKcliques_histogram",sep='')
    imgKcommunitySize  = paste(fkey, "_imgKcommunitySize_histogram",sep='')
    imgKcommunityCount = paste(fkey, "_imgKcommunityCount_histogram",sep='')

    logFname       = paste(fkey, "_log.xls",       sep='')
    fcliquesOnly   = paste(fkey, ".xls",           sep='') #cliques
    fcommunityOnly = paste(fkey, "_kcommunity-nodeMembership.xls",sep='') #k-community
    fcommunityLinks= paste(fkey, "_kcommunity-linkmembership.xls",sep='') #k-community

    if(1==2) {
    newlinkpairs  = linkpairs
    degrees       = degree_ByLinkPairs (newlinkpairs, directed=F)
    rets          = makeAjacencyLists(inetmatrix=linkpairs, returnUniqueNames=T)
    #adjlists  = makeAjacencyLists(inetmatrix=linkpairs, returnUniqueNames=F)
    nodenames = rets[[1]]
    no.nodes  = rets[[2]]
    # put link index into the 3rd column, for k-community links output
    #
    no.links  = dim(newlinkpairs)[1]
    }

    linkpairsIndexed = cbind(selectedLinksIndex, c(1:no.links) )

    # find the maximum coreness, which will be used for inferring max(k)-clique
    #
    #maxcoreness  = find_maxcoreness(selectedLinksIndex)
    CORES = KCores_ALL(linksindex=selectedLinksIndex, 
                       nodesindex=NULL, minkcore=min_clique-1, name2integer=T, nodeMembershipOnly=F)#!(MemoyForSpeed) )
    #globalAdjLists  = makeAjacencyListsFromLinksIndex(selectedLinksIndex)

    allcores.readTime <- proc.time()[3]-now

    print(paste("time for searching cores: ", allcores.readTime, sep="") )

    if(MemoyForSpeed){
       coreindex = CORES[[1]]
    } else{
       coreindex = as.integer(colnames(CORES))
    }

    cliqueindex  = coreindex  + 1 
    max_coreness = max(coreindex )
    max_clique   = max(cliqueindex) + 1
    maxcliqueVect= rep(max_clique, no.nodes)
    no.cores     = length(coreindex)
    
    allSelnodesIndex=c(1:no.nodes)

    #-----------------------------  1) finding all possible cliques ----------------------------
    #
    # explore all possible combinations by a recursive function
    #
    print("1. finding all possible cliques")

    iclqFindTimes = NULL
    for (i in c(no.cores:1) ) {

        iclq= cliqueindex[i]
        iclq.now<-proc.time()[3]

        # remove the nodes with the total links less than min_clique
        #
        #rets = kcore_subnets(linksindex=selectedLinksIndex, 
        #                     nodesindex=NULL, kcore=iclq-1)
        #
        #if ( !is.null(rets[[1]]) ){
        #     iclqNodesIndex     = sort(rets[[1]]) #selected nodes (index)
        #     iclqLinkpairsIndex = rets[[2]]       #selected links (indices)
        #}else{
        #     next
        #}

        # nodes in (iclq-1)-core
        if (i==1) {
           iclqNodesIndex     = c(1:no.nodes)
           iclqLinkpairsIndex = selectedLinksIndex
        } else{
          if(MemoyForSpeed) {
             iclqNodesIndex     = (CORES[[2]])[[i]]
             iclqLinkpairsIndex = (CORES[[3]])[[i]]
          }else{
             iclqNodesIndex     = allSelnodesIndex[ CORES[,i] ]
             merged1= merge(iclqNodesIndex,selectedLinksIndex,by.x=1,by.y=1,all=F)
             iclqLinkpairsIndex = merge(iclqNodesIndex,merged1,by.x=1,by.y=2,all=F)
          }
        }

        # the index of the lists correspond to the original node index
        #
        #iclqAdjlists   = makeAjacencyListsFromLinksIndex(iclqLinkpairsIndex)
        if (i==1) {
           iclqAdjlists   = globalAdjLists
        } else{
           iclqAdjlists   = makeAjacencyListsFromSubsetnodesindex(
                            orgAdjLists=globalAdjLists, 
                            subsetNodes=iclqNodesIndex)
        }

        no.nodesINiclq = length(iclqNodesIndex)
        no.linksINiclq = dim(linkpairsIndexed)[1]

        # iterate each node to identify all iclq-cliques including this node
        #
        for ( j in c(1:no.nodesINiclq) ) {
           cliques_count = 0
           v = iclqNodesIndex[j]
           
           #if ( length(iclqAdjlists[[v]]) >= iclq-1 & maxcliqueVect[v]>=iclq ){
           if ( length(iclqAdjlists[[v]]) >= iclq-1){

               mystr = paste("Clique-", as.character(iclq), ": ", as.character(j), "/", 
                             as.character(no.nodesINiclq), sep="")
               print(mystr)

               #cliques_count = recursive_CliquesNormal(curnode=v, setA=NULL, setB=NULL,
               cliques_count = recursive_Cliques_Super(curnode=v, setA=NULL, setB=NULL,
                                         adjlists     = iclqAdjlists,
                                         neighborlinks=iclqLinkpairsIndex, Kclique=iclq)#-1)

           }#if ( length(iclqAdjlists[[v]]) >= iclq-1
           cliques_count

           #-------------------- important -----------------------
           # i) if no cliques found, keep the node and its links
           #
           if (cliques_count ==0){
               next
           }

           # ii) remove the current node and its links
           #
           # ii.1) remove it from LinkPairs
           #
           selLeft    = iclqLinkpairsIndex[,1] != v
           selRight   = iclqLinkpairsIndex[,2] != v
           sel        = selLeft & selRight
           no.remains = sum(sel)

           # no enough links left, stop the current clique finding
           #
           if ( no.remains < iclq-1 ) {
              break
           } else {
              if(no.remains ==1) {
                 iclqLinkpairsIndex = rbind(iclqLinkpairsIndex[sel,])
              } else{
                 iclqLinkpairsIndex = iclqLinkpairsIndex[sel,]
              }
           }

           # ii.2) remove it from AdjList
           #
           iclqAdjlists[[v]] = -1           
           for ( m in c(1:no.nodesINiclq) ) {
                qq = iclqNodesIndex[m]
                mneighbors = setdiff(iclqAdjlists[[qq]], v)
                if (length(mneighbors)==0) {
                     iclqAdjlists[[qq]] = -1
                }else{
                     iclqAdjlists[[qq]] = mneighbors
                }
            }

            #iclqAdjlists
            #iclqLinkpairsIndex 

        }#for (j in c(

       # recording time
       iclq.readTime <- proc.time()[3]-iclq.now
       iclqFindTimes =rbind(iclqFindTimes, c(iclq, no.nodesINiclq, length(iclqLinkpairsIndex), iclq.readTime) )

    }#for (iclq in

    clique.readTime <- proc.time()[3]-now

    # save execution time
    #
    colnames(iclqFindTimes ) <- c("k-clique", "(k-1)-core nodes", "(k-1)-core links", "time")
    write.table(iclqFindTimes, logFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

    minutes = as.character( round(allcores.readTime/60.0, 3) )
    mystr = paste("\nTime for identifying all cores = ", as.character(allcores.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("\nTime for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)\n", sep="")
    print(mystr)
    appendStringToFile(logFname, mystr, newline=T)

    #----------------------------- 2) plot distribution of cliques ----------------------------
    #
    #globalCliques = getOption("globalCliques")

    no.kcliques= length(globalCliques)
    kcliquesSize = rep(0, no.kcliques)
    for(i in c(1:no.kcliques )){
       kcliquesSize[i] = length(globalCliques[[i]])
    }
    maxk=max(kcliquesSize)
    mink=min(kcliquesSize)

    cfmatrix   = histogram_4integers(ivector=kcliquesSize, fkeyname=imgKclique, keyword="k-cliques")
    write.table(cfmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)    


    #-------------------------------- 3) print out k-cliques -----------------------------------
    #
    print("3. Output all possible cliques")

    write.table(rbind(c("k-clique","member")), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    for (i in c(1:no.kcliques) ){
        icliques = paste(as.character( length(globalCliques[[i]]) ), "\t", sep="")
        for (jnode in globalCliques[[i]]){
            icliques = paste(icliques, selectedNodesName[jnode],", ", sep="")
        }
        #print(icliques)
        write.table(rbind(icliques), fcliquesOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    }
    
    # return if clique only
    if(cliqueOnly) {
        return (globalCliques)
    }


    # ======================== 4) prepare clique adjacency matrix ================================
    #
    print("4. prepare clique adjacency matrix")

    cliqueAdjMatrix = matrix(0,no.kcliques, no.kcliques)
    diag(cliqueAdjMatrix) <- kcliquesSize

    for (i in c(1:(no.kcliques-1)) ){
       for (j in c((i+1):no.kcliques) ){
          ijoverlap            = intersect(globalCliques[[i]], globalCliques[[j]])
          cliqueAdjMatrix[i,j] = length(ijoverlap)
          cliqueAdjMatrix[j,i] = length(ijoverlap)
       }
    }

    #++++++++++++++++++++++++ initialization ++++++++++++++++++++++++++++++++++++
    #
    kcommunities       = cbind(selectedNodesName)
    #kcommunitiesLinks = cbind(linkpairs)
    leftnodeNames      = selectedNodesName[ selectedLinksIndex[,1] ]
    rightnodeNames     = selectedNodesName[ selectedLinksIndex[,2] ]
    kcommunitiesLinks  = cbind(leftnodeNames, rightnodeNames)

    kcommunityTitle          = c("gene")
    kcommunitySize           = NULL
    kcommunityComponentCount = NULL
    kcommunityComponentCount2= NULL

    # ======================== 5) merge cliques into communities ================================
    #
    for (kc in c(maxk:mink) ){

         kcname = paste(as.character(kc), "-communities", sep="")

         print(paste("5... identify ", kcname) )

         # ~~~~~~~~~~~~  clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~
         #

         # <1> overlap size
         kcCliqueAdjMatrix = ifelse(cliqueAdjMatrix == kc-1, 1, 0) #off diagonal

         # <2> clique size
         selDiag           = kcliquesSize == kc

         #nselIndex         = c(1:no.kcliques)[!selDiag]
         #kcCliqueAdjMatrix[nselIndex,] = 0
         #kcCliqueAdjMatrix[,nselIndex] = 0
         #diag(kcCliqueAdjMatrix) = selDiag

         # <3> find the sub matrix with connected kcliques
         #
         diag(kcCliqueAdjMatrix) <- 0
         overlapvect = apply(kcCliqueAdjMatrix, 1, sum)

         selkcs      = overlapvect >=1
         selkcs      = (selkcs & selDiag) | selDiag # add in those isolated cliques

         if ( sum(selkcs)==0 ){
             next
         }

         #
         #
         # ~~~~~~~~~~~~  END of clique overlap  definition ~~~~~~~~~~~~~~~~~~~~~~~~

         kcommunityMatrix = kcCliqueAdjMatrix[selkcs, selkcs]
         kcindex          = c(1:no.kcliques)[selkcs]
         

         # ~~~~~~~~~~~  find connected components in k-cliques ~~~~~~~~~~~~~~~~~~~~
         #
         kc_CClabel = find_ConnectedComponents_adjmatrix(adjmatrix=kcommunityMatrix, 
                                  minModuleSize=0, cclabelstart=1)

         # merge cliques in each CC
         CCsize        = table(kc_CClabel[,2])
         CCs           = names(CCsize)
         no.components = length(CCs)

         kcommunityComponentCount = c(kcommunityComponentCount, rep(kc, no.components) )
         kcommunityComponentCount2= c(kcommunityComponentCount2,no.components)

         CCsizes = NULL
         for (cc in CCs){ # cc is like "0001","0002"

            # make title
            kcommCCgeneIndices = NULL
            kcommCCname        = paste(kcname, "-cc", cc,sep="")
            kcommunityTitle    = c(kcommunityTitle, kcommCCname)

            # get members of the current connected components in k-community
            #
            icc  = as.integer(cc)
            isel = cc==kc_CClabel[,2]
            cliquesIdx = kcindex[isel]
            mymembers = NULL
            for (m in cliquesIdx){
                 mymembers = union(mymembers, globalCliques[[m]] )
            }
            mymembers  = as.integer(mymembers) # index
            mvect      = rep(0, no.nodes)
            mvect[mymembers] = 1
 
            kcommunities          = cbind(kcommunities, mvect)         
            CCsizes               = c( CCsizes, length(mymembers) )

            kcommunitySize        = c(kcommunitySize, length(mymembers) )

            # network view
            imgKcommunity = paste(fkey, "_", kcommCCname, ".png", sep="")        
            kcLinksIndex  = getSubnetwork_LinkPairs(linkpairsIndexed,subnetNodes= mymembers)
            kcLinksIndex  = as.matrix(kcLinksIndex)
            kcLinks       = kcLinksIndex[,c(1:2)]
            kcPairsIndex  = as.integer(kcLinksIndex[,3])

            # from index to actual node name
            kcLinksByName = cbind(selectedNodesName[kcLinks[,1]], selectedNodesName[kcLinks[,2]])

            plotNetwork_inPairs_DisplayGeneSymbol(linkpairs=kcLinksByName, 
                           directed=F, geneinfo=transMatrix, genesymbolCol=2,
                           nodecolor="blue", nodeshape=30, nodenamecolor="purple",
                           labelpos=1, rankhubby="totallinks",
                           disphubs=0, plotfigure=T, 
                           fimg=imgKcommunity, saveDegrees=F)


            # output k-community links membership
            #
            mlinkvect = rep(0, no.links)
            mlinkvect[kcPairsIndex] = 1
            kcommunitiesLinks = cbind(kcommunitiesLinks, mlinkvect)
         }

    }

    # ======================== 6) output node/link membership & statistics ==============================
    #
    print("6. plot histogram of number of components in k-community")

    newtitle      = paste("numbers of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunityComponentCount, 
                                  fkeyname=imgKcommunityCount, 
                                  keyword=newtitle, hxlabel="k")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    newtitle      = paste("sizes of k-clique-communities", sep="")
    kcommCCmatrix = histogram_4integers(ivector=kcommunitySize, 
                                  fkeyname=imgKcommunitySize, 
                                  keyword=newtitle, hxlabel="community size")
    appendStringToFile(logFname, "\n", newline=F)
    write.table(kcommCCmatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    minutes = as.character( round(clique.readTime/60.0, 3) )
    mystr = paste("Time for searching cliques = ", as.character(clique.readTime), " seconds", 
                  "(", as.character(minutes), " minutes)", sep="")
    appendStringToFile(logFname, mystr, newline=F)

    # ======================== 7) output node/link membership & statistics ==============================
    #
    print("7. output node/link membership & statistics")

    # save the number of connected components in each k-community
    #
    #countMatrix = rbind(as.character(maxk:mink), kcommunityComponentCount2)
    #countMatrix = cbind(c("k-community","components"),countMatrix)
    #appendStringToFile(logFname, "\n", newline=T)
    #write.table(countMatrix, logFname, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    # community-component NODE membership ony
    write.table(rbind(kcommunityTitle), fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunities, fcommunityOnly, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)
    
    # community-component LINKS membership ony
    kclinktitle = c("src", "dst", kcommunityTitle[-1])
    write.table(rbind(kclinktitle),fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(kcommunitiesLinks, fcommunityLinks, sep="\t",quote=FALSE, col.names=F, row.names=FALSE,append=T)

    return (globalCliques)
}



############################### end of normal clique identification #########################################
#############################################################################################################











binaryVectorSimilarity = function(n1,n2, n12, measure="Correlation"){
   s11=n12/2
   s00=n12/2
   s10=n1-n12
   s01=n2-n12
   if (measure=="Correlation") {
      alpha      = sqrt((s10+s11)*(s01+s00)*(s11+s01)*(s00+s10))
      similarity = (s11*s00-s10*s01)/alpha
      similarity = (similarity +1)/2.0

   }else if (measure=="Tanimoto" | measure=="Russell-Rao"| measure=="Sokal-Michener"|measure=="Jaccard-Needham") {
      similarity = n12/(n1+n2-n12)

   }else if (measure=="Dice") {
      s11=n12
      similarity = 2*s11/(2*s11+s01+s10)

   }else if (measure=="Yule") {
      similarity = (s11*s00-s10*s01)/(s11*s00+s10*s01)
      similarity = (similarity +1)/2.0

   }else if (measure=="Rogers-Tanimoto" | measure=="Kulzinsky") {
      s11=n12
      similarity = s11/(s11+2*s10+2*s01)

   } else if (measure=="Sigmoid") {
      x = (s11*s00-s10*s01)/(s11*s00+s10*s01)
      x = (x +1)/2.0
      similarity = 1/(1+exp(-5*(x-0.5)) )

   } else if (measure=="SigmoidDice") {
      s11=n12
      x = 2*s11/(2*s11+s01+s10)
      similarity = 1/(1+exp(-5*(x-0.5)) )

   }else{ # Tanimoto
      similarity = n12/(n1+n2-n12)
   }
   
   return(similarity)
}

options(a=c(1:10))
test_globalvariable = function()
{
  a=getOption("a")
  print(a)
  a=c(1:3)
  print (a)
  a=getOption("a")
  print(a)

  a=getOption("a")  
  print (a)
}



#
#**********************************  END of k-cliques **********************************




##################################### FUNCTION ##########################
#
# Either by CDF of histgram count
#
# 1) get PDF by either histogram or scdf function
# 2) for each input data point, look at the nearest index in PDF and get 
#    its probability (P)
# 3) for the given P, check qnorm to get the inverse normal transformed value
#
# toy example
#    xx=c(1,1,2,2,3,3,4,4,4,4,4.5)
#    xx2=ecdf(xx)
#    cbind(xx2$x,xx2$y)
#

inverseNormalTransform = function(input, which="cdf", breakpoints=100) {

    x = round(input, 5)
    N = length(x)

    ######################### BY CDF #################################
    #
    if(which =="cdf"){
        #ecdf(x, datadensity='density')
        #ecdf(x, datadensity='hist')
        xecdf     = ecdf(x, xlab="test", datadensity='density', pl=F)

        xecdf.x   = xecdf$x
        xecdf.y   = xecdf$y
        xecdf.np  = length(xecdf.y)

    } else {

        xhist  = hist(x, br = breakpoints, plot=F)
        #xhist2 = cbind(xhist$mids, xhist$counts, xhist$density)
        #xhist2

        xecdf.y = cumsum(xhist$counts)/N
        xecdf.x = xhist$mids
        xecdf.np  = length(xecdf.y)
    }

    # handle the boundary
    #
    xecdf.y[1]= (xecdf.y[2] + xecdf.y[1])/2
    xecdf.y[xecdf.np]= (xecdf.y[xecdf.np-1] + xecdf.y[xecdf.np])/2

    xmean = mean(input, na.rm=T)
    xsd   = sd(input,   na.rm=T)

    invernorm = qnorm(xecdf.y, mean=xmean, sd=xsd)
    xhist3    = cbind(xecdf.x, xecdf.y, invernorm)
    xhist3 


    # 2. inverse normal distribution
    #
    xp = rep(0, N)
    for(i in c(1:N) ){
      xp[i] = invernorm[ which.min(abs(x[i]-xecdf.x) ) ]
    }

    #xop = cbind(x, xp)
    #xop

    return (xp)
}

#
###################### END of Inverse Normal Transform ##############################




# ************************* Useful Functions ******************************************
#
#
adjusted_cortest<-function(x1, tt1, gender1, age1){
   
   nnaTT = !(is.na(tt1))
   nnaXX = !(is.na(x1))
   nna   = nnaTT & nnaXX
   x=x1[nna]
   tt=tt1[nna]

   gender=gender1[nna]
   age=age1[nna]

   res<-cor.test(resid(lm(tt~gender+age)),resid(lm(as.numeric(as.vector(x))~gender+age)))
   return( c(res$p.value, res$estimate) )
}

adjusted_cor_multitraits <- function(x, ys, gender, age){
    no.ys = dim(ys)[2]
    mres  = NULL
    for (j in c(1:no.ys)){
        jres = adjusted_cortest(x, ys[,j], gender, age)
        mres = c(mres, jres)
    }
    return (mres)
}

logSTD = function(vect){
    vect2= ifelse(vect==0, 10^-10,vect)
    return (sd(-log10(vect2),na.rm=T))
}

logMean = function(vect){
    vect2= ifelse(vect==0, 10^-10,vect)
    return (mean(-log10(vect2),na.rm=T))
}


################### get the indices of min/max's for the subclasses ###################
#
#  Notice that here iclass can be in random order, i.e., not necessarily in descendant 
#   or ascendant order
#
#ivect  = abs(b[1:30,5] - 390000)
#iclass = b[1:30,1]

findIndexOfMinMaxInSubclass <- function(ivect, iclass, minimum=T)
{
  # original index
  iidx   = c(1:length(ivect))

  # reorder all the vectors
  iorder = order(iclass)
  xclass = iclass[iorder]
  xvect  = ivect[iorder]
  xidx   = iidx[iorder]

  #
  if(minimum) {
    midx = tapply(xvect, xclass, which.min)
  }else{
    midx = tapply(xvect, xclass, which.max)
  }

  # find the index of the 1st member in each subclass along the xvect
  #
  startIdx = findUniqueIdx(xclass)
  
  # now we try to get the absolute index of the min/max element of each class in xvect
  #
  sIdxMatrix = cbind(xclass[startIdx], startIdx)
  mIdxMatrix = cbind(names(midx), midx)

  merged = merge(sIdxMatrix, mIdxMatrix, by.x=1, by.y=1, all=F)
  merged = as.matrix(merged)

  absIdx = as.integer(merged[,2]) + as.integer(merged[,3])-1

  # recover the original index
  #
  return (xidx[absIdx])
  
}

complete_pairs_self = function(grpA){
   ino = length(grpA)
   mres= NULL
   for(i in c(1:(ino-1)) ){
      ielem = paste(grpA[i], grpA[(i+1):ino], sep="\t")
      mres  = c(mres, ielem)
   }
   return(mres)
}


# --- TPR FPR False Positive Rates
# If do_balance=T, nFP will be divided by the ratio of no of negatives and that of positves to
#   make sure that no of positives is the same as no of negatives
#
compute_tpr_fpr_precision = function(positives, negatives, do_balance=False, maxv=NULL, minv=NULL, ncuts = 10, xstep=NULL) {

    if(is.null(maxv)){
       scores=c(positves, negatives)
       minc = min(scores)
       maxc = max(scores)
    }else{
       minc = minv
       maxc = maxv
    }

    # -----------------------------------------------------
    #
    no.total=length(positives)
    if(!is.null(xstep) ) {
       step = xstep
       ncuts2 = (maxc-minc)/step       
    } else{
       step = (maxc - minc)/ncuts
       ncuts2 = ncuts
    }

    cuts = NULL
    for ( i  in c(1:(ncuts2+1)) ) {
      ival = minc + (i-1)*step  
      cuts = c(cuts, ival)
    }
    cuts = sort(c(cuts, 0))
    ns   = length(cuts)

    final = matrix(0, ns, 4)
    for ( i  in c(1:ns) ) {
      ival = cuts[i]
      sel  = positives > ival

      final[i, 1] = ival

      nTP = sum(positives>ival); nFP= sum(negatives>ival); nTN= sum(negatives<=ival)
      final[i, 2] = nTP/no.total  # True Positive Rate, recall
      if( do_balance ){
          nFP = nFP/(length(negatives)/length(positives))
          nTN = nTN/(length(negatives)/length(positives))
      }

      final[i, 3] = nFP/(nFP+nTN)  # False Positive Rate, equivelant to nFP/length(negatives)
      if(nTP+nFP>0) {
        final[i, 4] = nTP/(nTP+nFP) # Precision (balanced)
      } else{ final[i, 4]=0}
    }

    colnames(final) = c("threshold", "TPR", "FPR", "Precision")
    final
}

plot_ROC = function(TPR, FPR, curvenames="", imgname=NULL, legendx=0.5, legendy=0.4,
                    xlab="False Positive Rate", ylab="True Positive Rate") {
     xcolors = c("blue", "green", "red", "brown", "magenta", "purple", "cyan")
     xpoints= seq(0, 1, 0.1)

     no.curves = dim(TPR)[2]
     icolor    = xcolors[c(1:no.curves)]
     if(!is.null(imgname) ){
        openImgDev(imgname,iwidth =600, iheight = 600)
     }

     par(mfrow=c(1, 1), mar=c(3, 3, 1, 1) + 0.1, cex=1.5, mgp=c(1,0.5,0) )

        matplot(x=FPR, y=TPR, 
                type="l", lty=1, lwd=2, col = icolor,
                xlab=xlab, ylab=ylab, main="", axes=F, pch="*", xlim=c(0,1), ylim=c(0,1))#xlim=c(0.015,0.98), ylim=c(0.015,0.98) )

        if(xlab=="False Positive Rate") {
          lines(x=c(0,1), y=c(0,1), col="black")
        }

        #abline(h=0,lty = 1, lwd = 1, col="black")
        for (i in xpoints ){
            lines(x=c(i,i), y=c(0,1), col="grey60", lty=3)
            lines(x=c(0,1), y=c(i,i), col="grey60", lty=3)
        }

        # boundaries
        lines(x=c(0,0), y=c(0,1), col="black")
        lines(x=c(0,1), y=c(0,0), col="black")
        lines(x=c(1,0), y=c(1,1), col="black")
        lines(x=c(1,1), y=c(1,0), col="black")

        axis(1, at =xpoints, labels = xpoints, pos=0, las=0, outer=F, tick = F)
        axis(2, at =xpoints, labels = xpoints, pos=0, las=1, outer=F, tick = F)

        legend(legendx, legendy, legend=curvenames,
                 lty=1,lwd=2,col=icolor, ncol =1, cex=1,pt.cex=1)

     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1, mgp=c(3,1,0))

     if(!is.null(imgname) ){
         dev.off()
     }
}

# ----------------------------- Permutation FDR ---------------------------

permuteVect = function(myvect){
  pidx = sample(c(1:length(myvect)), length(myvect), replace=F)
  return(myvect[pidx])
}

correlationFDR = function(myMatrix, byCol=T, permutes=1, local_fdr=F, cormethod="spearman") {
  if(byCol) {
     orgCor  = cor(myMatrix,    use = "pairwise.complete.obs", method=cormethod)
  } else{
     orgCor  = cor(t(myMatrix), use = "pairwise.complete.obs", method=cormethod)
  }
  orgCor  = abs(orgCor)
  diag(orgCor) <- 0

  cuts = seq(0, 1, 0.1); no.cuts= length(cuts)
  FDR  = matrix(0, no.cuts, 4); FDR[,1] = cuts

  localMatrix = matrix(0, dim(orgCor)[1], dim(orgCor)[2])
  print("permutation (random correlations ...")
  for (p in c(1:permutes) ) {

    if( (p%%100)==0 ) { print( paste("Permute ", p) ) }

    if(byCol) {
      pMatrix = apply(myMatrix, 2, permuteVect)
    }else{
      pMatrix = apply(myMatrix, 1, permuteVect)
    }

    permCor = cor(pMatrix, use = "pairwise.complete.obs", method=cormethod)
    permCor = abs(permCor)
    diag(permCor) = 0

    localMatrix = localMatrix + (permCor > orgCor)
    
    for(i in c(1:no.cuts) ){
        FDR[i,2] = sum(orgCor>=cuts[i],  na.rm=T)/2 + FDR[i,2]
        FDR[i,3] = sum(permCor>=cuts[i], na.rm=T)/2 + FDR[i,3]
    }
  }
  
  FDR[,2] = FDR[,2]/permutes
  FDR[,3] = FDR[,3]/permutes
  FDR[,4] = ifelse(FDR[,2]>0, FDR[,3]/FDR[,2], FDR[,4])

  localfdr = localMatrix/permutes

  res = list(FDR, localfdr)

  return (res)
  
}

# compute local FDR - single vector
# step=-round(log10(min(realvect)))+1
#
computeLocalFDR = function(realvect, perMatrix, minv, maxv, step=10^(-10), is_pvalue=T)
{
  ngenes = length(realvect)
  nperms = dim(perMatrix)[2]

  # pval cutoff
  if(is_pvalue) {
    delta= step; maxv=-log10(step)
    pc = 10^(seq(from=1, to=maxv, by=1))*step; 
    pc = sort(c(pc, pc*0.2, pc*0.5), decreasing = TRUE)
  } else {
    pc = seq(minv, maxv, by=step)
  }
  ncuts = length(pc)

  # FDR by pvcalue cuts
  #
  fdrtab  = NULL
  fdrlocal=NULL
  for (i in c(1:ncuts) ) {
    isel = realvect < pc[i]
    total= isel+0;
    
    if(sum(total)==0){pc=pc[1:(i-1)]; break}

    #print(lodcuts[i])

    iselCb = rep(FALSE, ngenes)
    fdr    = rep(0, ngenes)
    fdrAvg = 0

    jp   = repmat(realvect, 1, nperms)
    fpMatrix = perMatrix < pc[i]
    fdr  = apply(fpMatrix, 1, sum)/nperms
    fdr = ifelse(total==0, NA, fdr) # handle the case without true signal under the cut

    # pseudon global FDR by focusing on the genes with eQTL
    #
    fdrAvg = mean(fdr, na.rm=T)

    ifdr   = c(pc[i], fdrAvg, sum(total));
    fdrtab = rbind(fdrtab, ifdr)

    #--------- true local fdr ----------------
    fdrlocal = cbind(fdrlocal, fdr)
  }
  colnames(fdrtab) <- c("pvalue", "FDR", "Trues")
  colnames(fdrlocal) <- paste("P<", pc,sep="")

  # FDR by comparing trues and permuted
  jp   = repmat(realvect, 1, nperms)
  fpMatrix = perMatrix < jp
  fdrTF  = apply(fpMatrix, 1, sum)/nperms

  fdrlocal2 = cbind(fdrTF, fdrlocal)
  colnames(fdrlocal2) <- c("FDR(no_cutoff)", colnames(fdrlocal) )

  res = as.list(c(1:2))
  res[[1]] = fdrtab
  res[[2]] = fdrlocal2

  return(res)
}


# compute local FDR - Matrix
computeSimpleFDR_Matrix = function(realMatrix, perMatrixList, is_pvalue=T, use_absolute_value=FALSE)
{
  nrows = dim(realMatrix)[1]; ncols = dim(realMatrix)[2]; 
  nperms = length(perMatrixList)

    fdr    = matrix(0, nrows, ncols)
    for(j in c(1:nperms) ){
       if(is_pvalue) {
          jsel= perMatrixList[[j]] < realMatrix
       } else{
          if(use_absolute_value) {
             jsel= abs(perMatrixList[[j]]) > abs(realMatrix)
          } else{
             jsel= perMatrixList[[j]] > realMatrix
          }
       }

       fdr = fdr + jsel
    }
    fdr = fdr/nperms

    return(fdr)
}

# compute local FDR - Matrix
computeGlobalFDR_Matrix = function(realMatrix, perMatrixList, is_pvalue=T)
{
  nrows = dim(realMatrix)[1]; ncols = dim(realMatrix)[2]; 
  nperms = length(perMatrixList)

  fdr    = matrix(0, nrows, ncols)
  perdat = NULL
  for(j in c(1:nperms) ){
     perdat = c(perdat, as.numeric(perMatrixList[[j]]) )
  }
  perdat= abs(perdat)
  pmean = mean(perdat); psd=sd(perdat)

  if(is_pvalue) {
     fdr = pnorm(abs(realMatrix), mean=pmean, sd=psd, lower.tail=TRUE)
  } else {
     fdr = pnorm(abs(realMatrix), mean=pmean, sd=psd, lower.tail=FALSE)
  }

  return(fdr)
}

# a matrix with three cols: MDC random_MDC_mean random_MDC_sd
#
MDC_FDR = function(mdcvect){
  if(mdcvect[1] < 1) {
     fdr = sum(mdcvect[1]>mdcvect[-1])/length(mdcvect[-1])
  } else {
     fdr = sum(mdcvect[1]<mdcvect[-1])/length(mdcvect[-1])
  }
  return (fdr)
}

# a matrix with three cols: MDC random_MDC_mean random_MDC_sd
#
MDC_FDR_by_normal_distr = function(mdcvect){
  if(mdcvect[1] < 1) {
     fdr = pnorm(mdcvect[1], mean=mdcvect[2], sd=mdcvect[3], lower.tail=TRUE)
  } else {
     fdr = pnorm(mdcvect[1], mean=mdcvect[2], sd=mdcvect[3], lower.tail=FALSE)
  }
  return (fdr)
}



############################# ACE Functions ##################################
#-------------------------------------------------------
# remove the first nrm char in each string of the oldList
# if any string has length less than nrm, return the whole 
# original list
rmFirstchar <- function(oldList,nrm=1){

   no.ele  = length(oldList)
   newList = oldList
   vlength = nchar(oldList)
   isIllegal = ifelse(vlength>nrm,0,1)
   if (sum(isIllegal)==0){
      for (i in c(1:no.ele)){
         newList[i]=substr(oldList[i],nrm+1,vlength[i])
      }
   }
   return(newList)
}
#--------------------------------------------------------
tscores <- function(DataEx,class.vec,class.value=c(0,1)){

   no.g    = dim(DataEx)[1]
   tscores = matrix(0,no.g,2)
   group0  = which(class.vec==class.value[1])
   group1  = which(class.vec==class.value[2])

   for (i in c(1:no.g)){
      ttmp = t.test(DataEx[i,group0],DataEx[i,group1])
      tscores[i,]=c(ttmp$stat,ttmp$p.value)
   }

   return(tscores)
}
# ---------------------------------------------------------
getWmatrix <- function(Coord, sigma1=100){
   
   no.locus = length(Coord)
   D        = repmat(Coord,1,no.locus)
   DIS      = D-t(D)
   power1   = DIS^2/sigma1
   Wmatrix  = exp(-power1)
   return(Wmatrix)   #Wmatrix is symmetric
}
# --------------------------------------------------------
# get weight vector based on specified location
getWvector <- function(Coord, specLoc,sigma1=100){
   
   no.locus = length(Coord)
   DIS      = Coord - specLoc
   power1   = DIS^2/sigma1
   Wvector  = exp(-power1)   
   return(Wvector)   
}

# ---------------------------------------------------------
getWmatrixbyIdx <- function(Coord, sigma1=100){
   
   no.locus = length(Coord)
   ord_     = order(Coord)
   D        = repmat(ord_,1,no.locus)
   DIS      = D-t(D)
   power1   = DIS^2/sigma1
   Wmatrix  = exp(-power1)
   return(Wmatrix)   #Wmatrix is symmetric
}
# ---------------------------------------------------------
getNS <- function(Wmatrix,score_,noperm=1000){
   
   no.locus = length(score_)
   ns.obs   = Wmatrix%*%score_

   ns.pv    = numeric(no.locus)
   ns.rand  = matrix(0,no.locus,noperm)
   for (j in c(1:noperm)){
      score.r = sample(score_,replace=F)
      ns.rand[,j] = Wmatrix%*%score.r 
   }
   for (i in c(1:no.locus)){
      if (ns.obs[i]>0){
         ns.pv[i]=length(which(ns.rand[i,]>=ns.obs[i]))
      }else{
         ns.pv[i]=length(which(ns.rand[i,]<=ns.obs[i]))
      }
   }
   ns.pv =ns.pv/noperm

   rm(list="ns.rand") #,"tmp1")
   collect_garbage()

   return(cbind(ns.obs,ns.pv))
}

# ------------------------------------------------------------
# collapse the multiple probes
getUniqueSet <- function(ExpM,GID,ChrV,LocV,byID=1,addID=1){

   IDsys     = byID  #genesymbol
   probes.no = dim(ExpM)[1]
   infoA     = cbind(GID[,IDsys],c(1:probes.no))

   uniqID    = unique(GID[,IDsys])
   uniqID.no = length(uniqID)
   infoU     = cbind(uniqID,c(1:uniqID.no))

   collapInf = merge(infoU,infoA,by.x=1,by.y=1,all=F)
   collapInf = as.matrix(collapInf)

   samples.no = dim(ExpM)[2]
   Ginfo  = cbind(GID[,byID],GID[,addID],ChrV,LocV)
   Ginfo.n= Ginfo[1:uniqID.no,]
   ExpM.n = matrix(0,uniqID.no,samples.no)
   IdxTab = matrix(as.integer(collapInf[,2:3]),ncol=2)  # old_list = filnal_list[IdxTab[,2]]
   IA = IdxTab[!duplicated(IdxTab[,1]),1]  #final_list=uniqID[IA]
   
   kk = 1
   for (i in c(1:uniqID.no)){
      if(i%%500==0) print(paste("Collapse ", i) )
      ind_=NULL
      while (IdxTab[kk,1]==IA[i]){
         ind_ = c(ind_,IdxTab[kk,2])
         kk   = kk+1
         if (kk >probes.no){
            break
         }             
      }
     
      if (length(ind_)>1){
         ExpM.n[i,] = as.numeric(colMeans(ExpM[ind_,],na.rm=TRUE))
      }else{
         ExpM.n[i,] = as.numeric(ExpM[ind_,])
      }
      Ginfo.n[i,] = Ginfo[ind_[1],]
      #ChrV.n[i] = ChrV[ind_[1]]
      #LocV.n[i,]= LocV[ind_[1],]
      #GID.n[i]  = GID[ind_[1],]
      #GIDadd.n[i] = GID[ind_[1],addID]
   }
   #GID.n=cbind(GID.n,GIDadd.n)
   rm(list=c("ExpM","GID","ChrV","LocV"))
   collect_garbage()

   return(list(ExpData=ExpM.n,GInfo=Ginfo.n))
}
#--------------------------------------------------------------------
getSegmentsChr <- function(NSinf,cut.n=5){

   nogenes = length(NSinf)
   marked  = NULL
   count.p = 1
   segID   = 1
   init    = NSinf[1]

   for (i in c(2:nogenes)){
      prod = NSinf[i]*init
      if (prod>0){
         count.p = count.p+1
         if (count.p==cut.n){
            stIdx = i-cut.n+1
            endIdx= i
         }else{
            if (count.p>cut.n){
                endIdx = i
            }
            #endIdx=ifelse(count.p>cut.n,i,0)
         }
      }else{
         init = NSinf[i]
         if (count.p>=cut.n){
             marked=rbind(marked,c(segID,stIdx,endIdx))
             segID = segID+1
         }             
         count.p = 1
      }
   }
   if (count.p>=cut.n){  #end of Chr
      marked=rbind(marked,c(segID,stIdx,endIdx))
   }
   return(marked)
}

# ----------------------------------------
getZscoreMatrix <- function(ExpData){

   nsamples = dim(ExpData)[2]
   ngenes   = dim(ExpData)[1]

   Zmatrix  = matrix(0,ngenes,nsamples)
   for (j in c(1:nsamples)){
      tmpExp = ExpData[,-j]
      tmpMean= rowMeans(tmpExp)
      tmpSDT = sd(t(tmpExp))
      Zmatrix[,j]= (ExpData[,j]-tmpMean)/(tmpSDT^2)
   }
   return(Zmatrix)
}



readPartialTable = function(finput, sep="\t", startstr=NULL, endstr=NULL, startDelta=0)
{
  if( is.null(startstr)  ){
     mtitle = read.delim(finput, sep=sep, header=F, nrow=1)
     mtitle = as.matrix(mtitle)[1,]     
     mt = read.delim(finput, sep=sep, header=T)
     colnames(mt) <- mtitle
     return (mt)
  }
  custr =""; startpos = 0;
  while(is.na(custr) | custr != startstr) {
      mt = read.delim(finput, sep=sep, header=F, skip=startpos, nrow=1)
      #print(mt)
      custr=as.character(as.matrix(mt)[1,1])
      startpos = startpos + 1
  }
  startpos = startpos+startDelta

  if( is.null(endstr)  ){
     mtitle = read.delim(finput, sep=sep, header=F, skip=startpos, nrow=1)
     mtitle = as.matrix(mtitle)[1,]     
     mt = read.delim(finput, sep=sep, header=T, skip=startpos)
     colnames(mt) <- mtitle
  } else{
     endpos = startpos
     while(is.na(custr) | (custr != endstr) ) {
         mt = read.delim(finput, sep=sep, header=F, skip=endpos, nrow=1)
         custr=as.character(as.matrix(mt)[1,1])
         endpos = endpos + 1
     }
     mtitle = read.delim(finput, sep=sep, header=F, skip=startpos, nrow=1)
     mtitle = as.matrix(mtitle)[1,]     

     mt = read.delim(finput, sep=sep, header=T, skip=startpos, nrow=endpos-startpos-2 )
     colnames(mt) <- mtitle
  }

  return (mt)
}

##################################### wavelet ###############################
# library(wavelets)
#
wavelet_transform = function(X, ncut=3, nl=4, wavelet= "d10")
{
  #nl=4; wavelet= "d10";#"d10", "bl20"; #boundary="periodic"
  mra.out <- mra(X, filter=wavelet, n.levels=nl, boundary="reflection")
  x1 = mra.out@S[[ncut]]
  return(x1[1:length(X)])
}


################################# keydriver analsyis #################################
#
# Inputs: 
# 1. linkpairs: input network in form of Nx2 matrix
# 2. signature: input signature for which key drivers need be identified
# 3. directed:  input network is directed or un-directed
# 4. nlayers:   how many layers you want to expand to find a subnetwork centering around the signature
# 5. min_downstreamnodes: minimum number of neighbors a key driver should have
# 6. FET_pvalue_cut: Fisher exact test p-value cutoff
# 7. boost_hubs: should we boost the nodes with many directed links to be key drivers?
# 8. dynamic_search: should we search different layers for the best enrichment?
# 9. bonferroni_correction: should we make Bonferroni correction of Pvalue
#
# Outputs:
#  1. X_Y_keydriver.xls: keydrivers for a signature Y in the network X
#  2. X_Y_cys.txt: the expanded subnetwork with key drivers highlighted
#  3. X_Y_cys-nodes.txt: the display properties (shape, color etc) in the expanded subnetwork
#  4. X_KDx_combined.xls: a combination of all the tables "X_Y_keydriver.xls"
#
#
keydriver_in_subnetwork = function(linkpairs, signature, directed=T, nlayers=6, 
                                   min_downstreamnodes=-1, FET_pvalue_cut=0.05, 
                                   boost_hubs=T, dynamic_search=T, bonferroni_correction=T) {

   allnodes = union(linkpairs[,1], linkpairs[,2])
   no.subnetsize = length(allnodes)

   # whole network nodes as the signature
   network_as_signature = length(setdiff(allnodes, signature)) ==0

   overlapped    = intersect(allnodes, signature)
   no.overlapped = length(overlapped) # within the subnetwork

   keydrivers= NULL
   kdMatrix  = NULL
   kdIndex   = NULL # indices of keydrivers in dsnodes_list

   dsnodes_list = as.list(rep(0,no.subnetsize)); no.dsnodes = rep(0, no.subnetsize)
   cnt = 1

   intv = as.integer(no.subnetsize/10)
   print("find downstream genes")

   # set up searching range
   if(dynamic_search) { # dynamic search for optimal layer
      layers_4search = c(1:nlayers)
   } else{  # fixed layer for optimal layer
      layers_4search = c(nlayers)
   }
   # if the network itself is the signature, no need for dynamic search
   if(network_as_signature){layers_4search = c(nlayers)}

   for(i in c(1:no.subnetsize) ) {

     if(i%%intv==0){ print(paste(i, "/", no.subnetsize)) }

     # initialization
     minpv=2;min_nohits=0;min_noidn=0;min_layer=0;min_dn=0; min_fc=0
     for(y in layers_4search) {
         #netpairs=linkpairs; seednodes=allnodes[i]; N=nlayers; directed=directed
         idn = downStreamGenes(netpairs=linkpairs, seednodes=allnodes[i], N=y, directed=directed)
         idn = setdiff(idn, allnodes[i])
         no.idn=length(idn);
      
         if(!network_as_signature){# do enrichment test for only expanded subnetwork
            hits    = intersect(idn, overlapped)
            no.hits = length(hits)
            foldchg = (no.hits/no.idn)/(no.overlapped/no.subnetsize)

            if(no.hits==0){next}
            pv = phyper(no.hits-1, no.idn, no.subnetsize-no.idn, no.overlapped, lower.tail=F)
            if(pv<minpv){
               minpv=pv;min_nohits=no.hits;min_noidn=no.idn;min_layer=y;min_fc=foldchg
            }
         } else{ # for non-expanded subnetwork
            no.hits = no.idn
            minpv=0;min_nohits=no.idn;min_noidn=no.idn;min_layer=y;min_fc=1
         }  
     } #y
    
     # record the down stream genes for the biggest layer
     if(no.idn>0) {
        dsnodes_list[[i]] = idn
        no.dsnodes[i]     = no.idn
     }

     res= c(min_nohits,min_noidn,no.overlapped,length(signature),no.subnetsize,min_fc,min_layer,minpv,minpv*no.subnetsize)
     kdMatrix  = rbind(kdMatrix, res)
     #print(res)
   }

   mymincut = min_downstreamnodes
   if (min_downstreamnodes<=0){
      mymincut = mean(no.dsnodes) + sd(no.dsnodes)
   }
   cutmatrix = c( mean(no.dsnodes), sd(no.dsnodes), mymincut)

   # pick up key drivers by pvalue and no. of downstream genes
   ncols = dim(kdMatrix)[2]

   if(bonferroni_correction) { # use corrected pvalue
      kdSel = (kdMatrix[,ncols] < FET_pvalue_cut) & (kdMatrix[,2]>= mymincut)
   } else{
      kdSel = (kdMatrix[,ncols-1] < FET_pvalue_cut) & (kdMatrix[,2]>= mymincut)
   }

   if( sum(kdSel)==0){return (NULL)}

   keydrivers= allnodes[kdSel]
   kdIndex   = c(1:no.subnetsize)[kdSel]
   n.drivers = length(keydrivers)

   #******************* local driver or not **************************************
   #
   # check whether a driver is in the downstream of other drivers
   keydrv = rep(0, no.subnetsize)
   #if (!network_as_signature) {
   for ( i in c(1:n.drivers) ) {

       # Note that kdIndex[i] is the index of ith keydriver in kdMatrix  
       # restrict to only candidate drivers 
       iselA  = (kdMatrix[,2] > kdMatrix[ kdIndex[i],2]) & kdSel
       isel   = c(1:no.subnetsize)[iselA]

       if ( sum(isel)>0) {
          if(directed) {
             ilocal= setInSets(setC=allnodes[ kdIndex[i] ],     setlist=dsnodes_list[isel])
          } else{
             ilocal= setInSets(setC=dsnodes_list[[ kdIndex[i] ]], setlist=dsnodes_list[isel])
          }
          keydrv[ kdIndex[i] ] = !ilocal + 0
       } else{
          keydrv[ kdIndex[i] ] = T
       }
   }
   #}

   # promote genes with many direct links to be key drivers
   #
   #              inlinks outlinks totallinks
   #0610031J06Rik       2        0          2
   #1110001J03Rik       0        1          1
   #
   if(boost_hubs) {
  
     if(!network_as_signature){
        # for expanded network, restrict the boosted nodes to the key driver candidates
        kdSelB = rep(F, no.subnetsize); kdSelB[kdIndex]=T;
     } else{
        # for non-expanded network, consider all the nodes in the subnetwork
        kdSelB = rep(T, no.subnetsize);
     }

     mydegree  = degree_ByLinkPairs(linkpairs=linkpairs, directed=directed, cleangarbage=F)
     if(directed) {
        directSel = mydegree[,2]> mean(mydegree[,2]) + 2*sd(mydegree[,2])
        cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), mymincut,
                       mean(mydegree[,2]),sd(mydegree[,2]), mean(mydegree[,2]) + 2*sd(mydegree[,2]) ))
     }else{
        directSel = mydegree[,3]> mean(mydegree[,3]) + 2*sd(mydegree[,3])
        cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), mymincut,
                       mean(mydegree[,3]),sd(mydegree[,3]), mean(mydegree[,3]) + 2*sd(mydegree[,3])))
     }
     directSel = directSel & kdSelB

     directeHub  = rownames(mydegree)[directSel]
     isDirectHub = setElementInSet(allnodes, directeHub)

     keydrv[isDirectHub] = T
     kdSel = kdSel | isDirectHub
     colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", "cut_downstream",
                           "mean_degree", "sd_degree", "cut_degree")

   } else{
     cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), mymincut, "F"))
     colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", "cut_downstream", "boost_directhubs")
   }

   ##
   # in this case, signature is the network nodes themselves, so pvalue will be 0 for all nodes
   # so the driver will be the ones with most downsttream genes
   #
   is_signature = rep(0, no.subnetsize); names(is_signature) <- allnodes
   is_signature[overlapped]  = 1
   
   fkd = cbind(allnodes, is_signature, kdMatrix, keydrv)[kdSel,];

   if(sum(kdSel) >1 ) {
       nf.cols = dim(fkd)[2]
       if(network_as_signature){
          mo = order(-as.integer(fkd[,2]))
       } else{
          mo = order(as.numeric(fkd[,9]))
       }

       fkd = fkd[mo, ]
       # put key driver on the top
       mo  = order( -as.integer(fkd[,nf.cols]) )
       fkd = fkd[mo, ]
   } else{
       fkd = rbind(fkd)
   }

   colnames(fkd) <- c("keydrivers", "is_signature", "hits", "downstream", "signature_in_network", "signature", 
                      "network_size", "fold_change", "optimal_layer", "pvalue", "pvalue_corrected", "keydriver")

   ret = as.list(c(1:2))
   ret[[1]] = fkd
   ret[[2]] = cutmatrix

   return (ret)
}


########################################################################################################################
#
# identify the optimal neighborhood of a set of given genes such that the neighborhood is most enriched for a signatures
#

search_best_neighborhood_KDA = function(linkpairs, candidates, signature, directed=T, nlayers=6, FET_pvalue_cut=0.05, 
                                   dynamic_search=T, bonferroni_correction=T) {

   allnodes = union(linkpairs[,1], linkpairs[,2])
   no.subnetsize = length(allnodes)

   no.candidates = length(candidates)

   overlapped    = intersect(allnodes, signature)
   no.overlapped = length(overlapped) # within the subnetwork

   intv = as.integer(no.candidates/10)
   dsnodes_list = as.list(rep(0,no.candidates)); no.dsnodes = rep(0, no.candidates)

   # whole network nodes as the signature
   network_as_signature = length(setdiff(allnodes, signature)) ==0

   keydrivers= NULL
   kdMatrix  = NULL
   kdIndex   = NULL # indices of keydrivers in dsnodes_list
   
   # set up searching range
   if(dynamic_search) { # dynamic search for optimal layer
      layers_4search = c(1:nlayers)
   } else{  # fixed layer for optimal layer
      layers_4search = c(nlayers)
   }
   # if the network itself is the signature, no need for dynamic search
   if(network_as_signature){layers_4search = c(nlayers)}

   for(i in c(1:no.candidates) ) {

     # initialization
     minpv=1;min_nohits=0;min_noidn=0;min_layer=0;min_dn=0; min_fc=0
     for(y in layers_4search) {
         #netpairs=linkpairs; seednodes=candidates[i]; N=nlayers; directed=directed
         idn = downStreamGenes(netpairs=linkpairs, seednodes=candidates[i], N=y, directed=directed)
         idn = setdiff(idn, candidates[i])
         no.idn=length(idn);
      
         if(!network_as_signature){# do enrichment test for only expanded subnetwork
            hits    = intersect(idn, overlapped)
            no.hits = length(hits)
            foldchg = (no.hits/no.idn)/(no.overlapped/no.subnetsize)

            if(no.hits==0){next}
            pv = phyper(no.hits-1, no.idn, no.subnetsize-no.idn, no.overlapped, lower.tail=F)
            if(pv<minpv){
               minpv=pv;min_nohits=no.hits;min_noidn=no.idn;min_layer=y;min_fc=foldchg
            }
         } else{ # for non-expanded subnetwork
            no.hits = no.idn
            minpv=0;min_nohits=no.idn;min_noidn=no.idn;min_layer=y;min_fc=1
         }  
     } #y
     
     # record the down stream genes for the biggest layer
     if(no.idn>0) {
        dsnodes_list[[i]] = idn
        no.dsnodes[i]     = no.idn
     }

     res= c(min_nohits,min_noidn,no.overlapped,length(signature),no.subnetsize,min_fc,min_layer,minpv,minpv*no.subnetsize)
     kdMatrix  = rbind(kdMatrix, res)
     #print(res)
   }

   # pick up key drivers by pvalue and no. of downstream genes
   ncols = dim(kdMatrix)[2]

   if(bonferroni_correction) { # use corrected pvalue
      kdSel = (kdMatrix[,ncols] < FET_pvalue_cut)
   } else{
      kdSel = (kdMatrix[,ncols-1] < FET_pvalue_cut)
   }

   if( sum(kdSel)==0){return (NULL)}

   keydrivers= candidates[kdSel]
   kdIndex   = c(1:no.candidates)[kdSel]
   n.drivers = length(keydrivers)

   #******************* local driver or not **************************************
   #
   # check whether a driver is in the downstream of other drivers
   keydrv = rep(0, no.subnetsize)
   for ( i in c(1:n.drivers) ) {

       # Note that kdIndex[i] is the index of ith keydriver in kdMatrix  
       # restrict to only candidate drivers 
       iselA  = (kdMatrix[,2] > kdMatrix[ kdIndex[i],2]) & kdSel
       isel   = c(1:no.subnetsize)[iselA]

       if ( sum(isel)>0) {
          if(directed) {
             ilocal= setInSets(setC=allnodes[ kdIndex[i] ],     setlist=dsnodes_list[isel])
          } else{
             ilocal= setInSets(setC=dsnodes_list[[ kdIndex[i] ]], setlist=dsnodes_list[isel])
          }
          keydrv[ kdIndex[i] ] = !ilocal + 0
       } else{
          keydrv[ kdIndex[i] ] = T
       }
   }
   
   fkd = cbind(candidates, kdMatrix, keydrv)[kdSel,];

   if(sum(kdSel) >1 ) {
       nf.cols = dim(fkd)[2]
       mo = order(as.numeric(fkd[,nf.cols-1]))

       fkd = fkd[mo, ]
       # put key driver on the top
       mo  = order( -as.integer(fkd[,nf.cols]) )
       fkd = fkd[mo, ]
   } else{
       fkd = rbind(fkd)
   }

   colnames(fkd) <- c("keydrivers", "hits", "downstream", "signature_in_network", "signature", 
                      "network_size", "fold_change", "optimal_layer", "pvalue", "pvalue_corrected", "keydriver")

   return (fkd)
}




#####
concateNonNAs = function(ivect, title){
    nonna = !is.na(ivect)
    if(sum(nonna)==0){
       iret=c(0,"")
    } else {
       istr=concatenate(title[nonna], "/")
       iret= c(sum(nonna), istr)
    }
    return( rbind(iret))
}


# endflag="": end of the first table and start of the 2nd tabke
# blankline_before: a blank line before endflag?
# blankline_after=0: a blank line before endflag?
# skiplines: skip the first line

getSecondTable = function(inputfname, firstTable=T, endflag="", blankline_before=0, blankline_after=0, skiplines=0)
{
        iline <- read.delim(inputfname,sep=";", header=F)
        iline <- as.matrix(iline)[,1]

        sel = rep(FALSE, length(iline))
        for(flag in endflag) {
           isel = iline==flag
           sel = sel | isel
        }

     if(sum(sel)==0){print("I Couldn't find flags\n"); return (LL)}

     lastrow <- c(1:length(iline))[sel]    

     allMatrix <- read.delim(inputfname,sep="\t", header=T, nrow=lastrow-blankline_before-skiplines, skip=skiplines)
     if (firstTable){
        return(allMatrix)
     }
     
     allMatrix <- read.delim(inputfname,sep="\t", header=T, skip=lastrow + 1 + blankline_after)
     return (allMatrix)
}


getSecondTableOld = function(inputfname, firstTable=T, endflag="", blankline_before=0, blankline_after=0, skiplines=0)
{
     lastrow = 1
     while (1==1) {
        iline <- read.delim(inputfname,sep="\t", header=F, skip=lastrow, nrow=1)
        iline <- as.character(as.matrix(iline))
        if (length(iline)==1) {
           if (is.element(iline,endflag)) {break;}
        }
        if(!is.na(iline[1]) ) {
           if (is.element(iline,endflag)) {break;}
        }

        lastrow <- lastrow + 1    
     }

     allMatrix <- read.delim(inputfname,sep="\t", header=T, nrow=lastrow-blankline_before-skiplines, skip=skiplines)
     if (firstTable){
        return(allMatrix)
     }
     
     allMatrix <- read.delim(inputfname,sep="\t", header=T, skip=lastrow + 1 + blankline_after)
     return (allMatrix)
}


getMultiTables = function(inputfname, endflag="")
{
     endflag_index = NULL
     lastrow = 1
     while (1==1) {
        iline <- tryCatch(read.delim(inputfname,sep="\t", header=F, skip=lastrow, blank.lines.skip =F, nrow=1),
                          error=function(e) (return(NULL)) )
        if (is.null(iline) ) {endflag_index = c(endflag_index, lastrow+1);break;}
        #print(iline)
        iline <- as.character(as.matrix(iline))
        if (length(iline)==1) {
           if (is.element(iline,endflag)) {endflag_index = c(endflag_index, lastrow+1);}
        }
        if (is.na(iline[1])) { # sometimes, it is NA
           endflag_index = c(endflag_index, lastrow+1);
        }
        lastrow <- lastrow + 1
     }

     no.tabs = length(endflag_index)
     tab_skips = c(0, endflag_index);

     all_tables= as.list(rep(NA, no.tabs))
     for( i in c(1:no.tabs) ){ 
         ititle <- read.delim(inputfname,sep="\t", header=F, blank.lines.skip =F,
                                 nrow=1, skip=tab_skips[i])
         ititle <- as.matrix(ititle)[1,]

         allMatrix <- read.delim(inputfname,sep="\t", header=T, blank.lines.skip =F,
                                 nrow=tab_skips[i+1]-1-1-tab_skips[i], skip=tab_skips[i])
         allMatrix <- as.matrix(allMatrix)
         colnames(allMatrix) <- ititle

         all_tables[[i]]= allMatrix
     }
     all_tables

     return (all_tables)
}


################################### vendiagram #########################################

##############################################################################
#### R script to
#### 	Produce Venn Diagrams with 1 to 5 groups
####		an extension on the code from the limma package
####		
#### Written By: Matt Settles
####				Postdoctoral Research Associate
####				Washington State University
####
##############################################################################
####
#### Change Log: 
####	Feb 8, 2008: 
####		formalized code
####    Dec 23, 2008:
####	    added mixed type to vennCounts
##############################################################################
####
####	Usage:
####	source("http://bioinfo-mite.crb.wsu.edu/Rcode/Venn.R")
####	can change colors now on 4 way plot
##############################################################################

####################################
## Function ellipse
## Add an elipse to the current plot
"ellipse" <- 
function (center, radius, rotate, 
    segments = 360, add = FALSE, xlab = "", ylab = "", las = par("las"), 
    col = palette()[2], lwd = 2, lty = 1, ...) 
{
	# x' = x cos? + y sin?
	# y' = y cos? - x sin?
    if (!(is.vector(center) && 2 == length(center))) 
        stop("center must be a vector of length 2")
    if (!(is.vector(radius) && 2 == length(radius))) 
        stop("radius must be a vector of length 2")
	
	angles <- (0:segments) * 2 * pi/segments  
	rotate <- rotate*pi/180
	ellipse <- cbind(radius[1] * cos(angles), 
					 radius[2] * sin(angles))
	if(rotate != 0)
		ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate),
						  ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
	ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
    if (add) 
        lines(ellipse, col = col, lwd = lwd, lty = lty, ...)
    else    plot(ellipse, type = "l", xlim = c(-4, 4), ylim = c(-4, 4),
			xlab = "", ylab = "", axes = FALSE, col = col, lwd = lwd,
			lty = lty, ...)
}

###################################
## Function vennCounts
## Produce venn table object
"vennCounts" <- 
function (x, include = "both") 
{
    x <- as.matrix(x)
    include <- match.arg(include, c("both", "up", "down","mixed"))
    x <- sign(switch(include, both = abs(x), up = x > 0, down = x < 0, mixed = x))
    nprobes <- nrow(x)
    ncontrasts <- ncol(x)
    names <- colnames(x)
    if (is.null(names)) 
        names <- paste("Group", 1:ncontrasts)
    noutcomes <- 2^ncontrasts
    if ( include == "mixed" ) noutcomes <- 3^ncontrasts
    outcomes <- matrix(0, noutcomes, ncontrasts)
    colnames(outcomes) <- names
    for (j in 1:ncontrasts) {
    	if( include == "mixed"){
    		outcomes[, j] <- rep(-1:1, times = 3^(j - 1), each = 3^(ncontrasts - j))
    	} else {
    	    outcomes[, j] <- rep(0:1, times = 2^(j - 1), each = 2^(ncontrasts - j))
    	}	
    }
    xlist <- list()
    for (i in 1:ncontrasts) {
    	if( include == "mixed"){
	    	xlist[[i]] <- factor(x[, ncontrasts - i + 1], levels = c(-1,0, 1))
		} else {
	       	xlist[[i]] <- factor(x[, ncontrasts - i + 1], levels = c(0, 1))
	    }    
    }
    counts <- as.vector(table(xlist))
    structure(cbind(outcomes, Counts = counts), class = "vennCounts")
}

"vennDiagram" <-
function (object, include = "both", names, mar = rep(1, 4), cex = 1.5, 
    lwd = 2, circle.col, counts.col, show.include, ...) 
{
    if (!is(object, "vennCounts")) {
        if (length(include) > 2) 
            stop("Cannot plot Venn diagram for more than 2 sets of counts")
        if (length(include) == 2) 
            object.2 <- vennCounts(object, include = include[2])
        object <- vennCounts(object, include = include[1])
    } else if (length(include == 2)) {
        include <- include[1]
    }

    nsets <- ncol(object) - 1
    if (nsets > 4) 
        stop("Can't plot Venn diagram for more than 4 sets")
    if (missing(names)) 
        names <- colnames(object)[1:nsets]
    counts <- object[, "Counts"]
    #if (length(include) == 2)   
    if (length(include) == 2 && nsets < 4) 
        counts.2 <- object.2[, "Counts"]
    #Setup colors
	 
    if (missing(circle.col) & nsets == 4) {
        circle.col <- c("red","blue","orange","green")
    } else if (missing(circle.col)) { 
        circle.col <- par("col")
    }

    if (length(circle.col) < nsets) 
        circle.col <- rep(circle.col, length.out = nsets)
    if (missing(counts.col) & nsets == 4) 
        counts.col <- c("red","blue","orange","green")
	 else
	 	  counts.col <- par("col") 
    if (length(counts.col) < length(include) & nsets < 4) 
        counts.col <- rep(counts.col, length.out = length(include))
    else if (length(counts.col) < nsets) 
        counts.col <- rep(counts.col, length.out = nsets)
    
    if (missing(show.include)) 
        show.include <- as.logical(length(include) - 1)
    xcentres <- list(0, c(-1, 1), c(-1, 1, 0), c(-0.2,0.2,-1.05,1.05))[[nsets]]
    ycentres <- list(0, c(0, 0), c(1/sqrt(3), 1/sqrt(3), -2/sqrt(3)),c(.20,.20,-0.35,-0.35))[[nsets]]
    centers <- cbind(xcentres,ycentres)
	r1 <- c(1.5, 1.5, 1.5, 1.5)[nsets]
	r2 <- c(1.5, 1.5, 1.5, 2.7)[nsets]
	radius <- c(r1,r2)
	rotate <- list(0, c(0,0), c(0,0,0), c(-45,45,-45,45))[[nsets]]
	
    xtext <- list(-1.2, c(-1.2, 1.2), c(-1.2, 1.2, 0),c(-3.2,3.2,-3.2,3.2))[[nsets]]
    ytext <- list(1.8, c(1.8, 1.8), c(2.4, 2.4, -3),c(3.2,3.2,-3.2,-3.2))[[nsets]]
    old.par <- par(mar = mar)
    on.exit(par(old.par))
    #plot(x = 0, y = 0, type = "n", xlim = c(-4.0, 4.0), ylim = c(-4.0, 
    #    4.0), xlab = "", ylab = "", axes = FALSE, ...)
    xxx = 3.0; xxy=3.0
    plot(x = 0, y = 0, type = "n", xlim = c(-xxx, xxx), ylim = c(-xxy, 
        xxy), xlab = "", ylab = "", axes = FALSE, ...)

    for (circle in 1:nsets) {
		ellipse(centers[circle,],radius,rotate[circle],add=TRUE,
				, lwd = lwd, col = circle.col[circle])
        text(xtext[circle], ytext[circle], names[circle], 
			#pos=ifelse(circle != as.integer(circle/2)*2,4,2),
			offset = 0.5, col = circle.col[circle], cex = cex)
    }
#	rect(-3.5, -3.5, 3.5, 3.5)
 
	switch(nsets, {
        #rect(-3, -2.5, 3, 2.5)
        printing <- function(counts, cex, adj, col, leg) {
            col <- col[1]
			text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(0, 0, counts[2], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.3, -2.1, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        #rect(-3, -2.5, 3, 2.5)
        printing <- function(counts, cex, adj, col, leg) {
            col <- col[1]
			text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
            text(-1.5, 0.1, counts[3], cex = cex, col = col, 
                adj = adj)
            text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.3, -2.1, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        #rect(-3, -3.5, 3, 3.3)
        printing <- function(counts, cex, adj, col, leg) {
			col <- col[1]
            #text(2.5, -3, counts[1], cex = cex, col = col, adj = adj)
            text(0, -1.7, counts[2], cex = cex, col = col, adj = adj)
            text(1.5, 1, counts[3], cex = cex, col = col, adj = adj)
            text(0.75, -0.35, counts[4], cex = cex, col = col, 
                adj = adj)
            text(-1.5, 1, counts[5], cex = cex, col = col, adj = adj)
            text(-0.75, -0.35, counts[6], cex = cex, col = col, 
                adj = adj)
            text(0, 0.9, counts[7], cex = cex, col = col, adj = adj)
            text(0, 0, counts[8], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.5, -3, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        #rect(-3.5, -3.5, 3.5, 3.5)
        printing <- function(counts, cex, adj, col, leg) {
            text(0, -3, 		counts[1], cex = cex, col = "black", adj = adj)
            text(2.5, 0, 		counts[2], cex = cex, col = col[4], adj = adj)
			 lines(c(2.25,2.75),c(-0.2,-0.2),col=col[4])
			 
            text(-2.5, 0, 		counts[3], cex = cex, col = col[3], adj = adj)
             lines(c(-2.75,-2.25),c(-0.2,-0.2),col=col[3])
			 
			text(0, -2.0, 		counts[4], cex = cex, col = "black", adj = adj)
             lines(c(-0.25,0.25),c(-2.2,-2.2),col=col[3])
			 lines(c(-0.25,0.25),c(-2.25,-2.25),col=col[4])
			
			text(1.3, 2.1, 		counts[5], cex = cex, col = col[2], adj = adj)
             lines(c(1.05,1.55),c(1.9,1.9),col=col[2])
            
			text(1.7, 1.2, 		counts[6], cex = cex, col = "black", adj = adj)
             lines(c(1.45,1.95),c(1.0,1.0),col=col[2])
             lines(c(1.45,1.95),c(0.95,0.95),col=col[4])
            
			text(-1.6, -1.1, 	counts[7], cex = cex, col = "black", adj = adj)
             lines(c(-1.85,-1.35),c(-1.3,-1.3),col=col[2])
			 lines(c(-1.85,-1.35),c(-1.35,-1.35),col=col[3])
			 
			text(-0.8, -1.55, 	counts[8], cex = cex, col = "black", adj = adj)
			 lines(c(-0.55,-1.05),c(-1.75,-1.75),col=col[2])
			 lines(c(-0.55,-1.05),c(-1.8,-1.8),col=col[3])
			 lines(c(-0.55,-1.05),c(-1.85,-1.85),col=col[4])

			text(-1.3, 2.1, 	counts[9], cex = cex, col = col[1], adj = adj)
			 lines(c(-1.55,-1.05),c(1.9,1.9),col=col[1])
            
			text(1.6, -1.1, 	counts[10], cex = cex, col = "black", adj = adj)
             lines(c(1.85,1.35),c(-1.3,-1.3),col=col[1])
			 lines(c(1.85,1.35),c(-1.35,-1.35),col=col[4])

			text(-1.7, 1.2, 	counts[11], cex = cex, col = "black", adj = adj)
			 lines(c(-1.45,-1.95),c(1.0,1.0),col=col[1])
             lines(c(-1.45,-1.95),c(0.95,0.95),col=col[3])
            
			text(0.8, -1.55, 	counts[12], cex = cex, col = "black", adj = adj)
			 lines(c(0.55,1.05),c(-1.75,-1.75),col=col[1])
			 lines(c(0.55,1.05),c(-1.8,-1.8),col=col[3])
			 lines(c(0.55,1.05),c(-1.85,-1.85),col=col[4])
			 
			text(0, 1.6, 		counts[13], cex = cex, col = "black", adj = adj)
             lines(c(-0.25,0.25),c(1.4,1.4),col=col[1])
			 lines(c(-0.25,0.25),c(1.35,1.35),col=col[2])
			 
			text(0.9, 0.5, 		counts[14], cex = cex, col = "black", adj = adj)
             lines(c(1.15,0.65),c(0.3,0.3),col=col[1])
			 lines(c(1.15,0.65),c(0.25,0.25),col=col[2])
			 lines(c(1.15,0.65),c(0.2,0.2),col=col[4])
			 
			text(-0.9, 0.5, 	counts[15], cex = cex, col = "black", adj = adj)
             lines(c(-1.15,-0.65),c(0.3,0.3),col=col[1])
			 lines(c(-1.15,-0.65),c(0.25,0.25),col=col[2])
			 lines(c(-1.15,-0.65),c(0.2,0.2),col=col[3])
			 
			text(0, -0.5, 		counts[16], cex = cex, col = "black", adj = adj)
			 lines(c(-0.25,0.25),c(-0.7,-0.7),col=col[1])
			 lines(c(-0.25,0.25),c(-0.75,-0.75),col=col[2])
			 lines(c(-0.25,0.25),c(-0.8,-0.8),col=col[3])
			 lines(c(-0.25,0.25),c(-0.85,-0.85),col=col[4])			  
      } 
	} )
    adj <- c(0.5, 0.5)
    if (length(include) == 2 & nsets < 4) 
        adj <- c(0.5, 0)
    printing(counts, cex, adj, counts.col, include[1])
    if (length(include) == 2 & nsets < 4) 
        printing(counts.2, cex, c(0.5, 1), counts.col[2], include[2])
    invisible()
}


# get representative pathway for each module, avoid names too long
#
#     module      ModuleSize System                       Gene.Category                               
#  "turquoise" "2537"     "Panther Biological Process" "Pre-mRNA processing"                       
#  "turquoise" "2537"     "BP"                         "negative regulation of histone acetylation"
#  "turquoise" "2537"     "BP"                         "mRNA processing"    

# find the row with minimum value at vcol and return the retcol of the row
getMinLen = function(iMatrix, vcol, retcol){
   nrows = dim(iMatrix)[1]
   if (nrows ==1){
      return (as.integer(iMatrix[1,retcol]))
   }

   ivalues = as.integer(iMatrix[,vcol])
   imin = which( ivalues==min(ivalues))
   return ( as.integer(iMatrix[imin, retcol]) )
}

# return the respresentative pathway for each module
# module module_size pathway 
#
getRepresentativePathwayModule = function(goannot, topn=3, shortenname=T) 
{
  pathw_words <- getAllParts(fullfnames=goannot[,4], sep=" ", ret_len=T)
  allPathw    <- cbind(goannot[,c(1,2,4)], pathw_words) # module name, size, pathway
  dim(allPathw)
  no.pws = dim(allPathw)[1]

  iPathw = NULL
  xpw    = allPathw
  for ( i in c(1:topn) ) {
     psel   <- findUniqueIdx(xpw[,1])
     iPathw <- rbind(iPathw, xpw[psel,])
     xpw    <- xpw[-psel,]
  }
  iPathw

  mo = order(as.integer(iPathw[,4]))
  iPathw = iPathw[mo, ]

  mo = order(iPathw[,1])
  iPathw = iPathw[mo, ]

  uidx = findUniqueIdx(iPathw[,1])
  res  = iPathw[uidx ,-4]

  parts= getSecondPart(res[,3], " of ", 2)
  parts= getSecondPart(parts, " in ", 1)
  parts = ifelse(parts=="skeletal development", "extracellular matrix", parts)
  parts = ifelse(parts=="cell cycle phase", "cell cycle", parts) 

  res[,3] = parts

  return (res)
}

# return tthe top N most enriched pathways for each module
#
getTopNPathwayModule = function(goannot, topn=3) 
{
  n.rows = dim(goannot)[1]
  n.cols = dim(goannot)[2]

  iPathw = NULL
  xpw    = cbind(goannot, cbind(c(1:n.rows)) )
  for ( i in c(1:topn) ) {
     psel   <- findUniqueIdx(xpw[,1])
     iPathw <- rbind(iPathw, xpw[psel,])
     xpw    <- xpw[-psel,]
  }
  iPathw
  
  mo = order(as.integer(iPathw[,n.cols+1]))
  iPathw = iPathw[mo, ]

  return (iPathw[,-c(n.cols+1)])
}

#col the color for filling the polygon. The default, NA, is to leave polygons unfilled, unless density is specified. (For back-compatibility, NULL is equivalent to NA.) If density is specified with a positive value this gives the color of the shading lines. 
#border the color to draw the border. The default, NULL, means to use par("fg"). Use border = NA to omit borders. 
# rect(xleft, ybottom, xright, ytop,col="yellow", border="yellow")
#
circle = function(center, radius, col=NA, border=NA, lty=2, lwd=2){
    xc = seq(from=-radius, to=radius, by=radius/300)
    yc = (radius^2-xc^2)^0.5
    xx = c(xc, rev(xc))+center[1]; yy=c(yc, -yc) + center[2] # combine upper semi-circle and lower one
    polygon(x=xx,y=yy, col=col, lty=lty, lwd=lwd, border=border)
}

# bootstrap the orgMtrix and fixMatrix, then compute the correlations  
#
# orgMtrix, fixMatrix: genes as rows and samples as columns
#
correlation_bootstrap = function(orgMtrix, fixMatrix, no.bootstraps=20, sampling_rate=0.90, selected_entries = NULL)
{
   avg_rs = rep(0, no.bootstraps)
   no.arrays = dim(orgMtrix)[2]; aidx = c(1:no.arrays); no.samples = as.integer(sampling_rate*no.arrays)
   for (x in c(1:no.bootstraps)){
       idx = sample(aidx, no.samples, replace=F)
       xcor= cor(t(orgMtrix[,idx]), t(fixMatrix[,idx]), use="p")
       avg_rs[x]= mean(abs(xcor[selected_entries]))
   }
   return(avg_rs)
}


# ----------------------------------barplot3d_original<-----------------------------------------------
#
# library(rgl)

barplot3d_original<-function(height,names.arg=NULL,col=NULL,main=NULL, theta = 0, phi =0, fov = 60){

if(is.table(height)){
if(is.null(names.arg)){
names.arg<-names(height)
}
frec<-height[]
}else{
frec<-height
}


m<-length(frec)

#distancia entre las barras
dist<-c(0.2)

E<-matrix(rep(0,m*24),ncol=3,nrow=8*m)

#Base de la primera barra
E[1,]<-c(0,0,0)
E[2,]<-c(0,0,1)
E[3,]<-c(1,0,1)
E[4,]<-c(1,0,0)

#Altura de la primera barra
E[5,]<-c(0,frec[1],0)
E[6,]<-c(0,frec[1],1)
E[7,]<-c(1,frec[1],1)
E[8,]<-c(1,frec[1],0)


if(m>1){
j<-1;
for(i in seq(1,8*(m-1)-7,8)){

#Bases del resto de barras
E[i+8,]<-E[1,]+j*c(1+dist,0,0)
E[i+9,]<-E[2,]+j*c(1+dist,0,0)
E[i+10,]<-E[3,]+j*c(1+dist,0,0)
E[i+11,]<-E[4,]+j*c(1+dist,0,0)

#Alturas del resto de barras
E[i+12,]<-c(E[i+8,1],frec[j+1],E[i+8,3])
E[i+13,]<-c(E[i+9,1],frec[j+1],E[i+9,3])
E[i+14,]<-c(E[i+10,1],frec[j+1],E[i+10,3])
E[i+15,]<-c(E[i+11,1],frec[j+1],E[i+11,3])

j<-j+1

}
}

#rgl.open();
#limpia la ventana antes de dibujar el gr?fico
rgl.clear(type="shapes")

rgl.bg(color='white')
rgl.viewpoint(theta = theta, phi =phi, fov = fov, zoom = 1)

if(is.null(col)){
color<-rainbow(m)
}else{
color<-rep(col,m)
}


n<-1

for(t in seq(1,8*m-7,8)){

#pinta la base y la tapa
rgl.quads(E[t:(t+3),],col=color[n])
rgl.quads(E[(t+4):(t+7),],col=color[n])

#pinta la cara izquierda y la derecha
rgl.quads(rbind(E[t,],E[t+1,],E[t+5,],E[t+4,]),col=color[n])
rgl.quads(rbind(E[t+3,],E[t+2,],E[t+6,],E[t+7,]),col=color[n])

#pinta la cara frontal y la trasera
rgl.quads(rbind(E[t+1,],E[t+2,],E[t+6,],E[t+5,]),col=color[n])
rgl.quads(rbind(E[t,],E[t+3,],E[t+7,],E[t+4,]),col=color[n])

n<-n+1
}

axis3d(edge='y+',col="black",pos=c(-0.8,0,0.5))
if(m>5){
axis3d(edge='y+',col="black",pos=c(m*(1+dist),0,0.5))
}
title3d(main=main,col="black")

if(!is.null(names.arg)){
#Creando la posici?n de las etiquetas
txt<-matrix(rep(0,3*m),ncol=3,nrow=m)

txt[1,]<-c(0.5,0,1.5)
rgl.texts(txt[1,],text=names.arg[1],col="black")

for(i in seq(2,m)){
txt[i,]<-txt[1,]+(i-1)*c(1+dist,0,0)
rgl.texts(txt[i,],text=names.arg[i],col="black")
}

}

}

# ----------------------------------- barplot3d -------------------------------------------------------
#
# display each column as a group with rows as bars
#
#
# datMatrix: columns as groups and rows as group properties, 
#            row and col names, which will be displayed, should be valid
#
# firstColumnAsLegend: first column is a title or not, If yes, the program internally does some transformation
# rowcol: colors for the row elements
# legendxDelta: shift legend x coord
# show_value: show bar height values
# zoom_view: zoom in factor of the whole image
# xyz_scales: scales of x, y and z axes
# legend_xscale: scale the width of legend box
# cex: font size
# xcex: for labels along x axis
#
barplot3D_by_group = function(datMatrix, firstColumnAsLegend=F, rowcol=c("red", "green", "blue", "yellow"), 
                  title="", imgfile=NULL,legendxDelta=0, legendyDelta=0, show_value=T, 
                  imgwid=800, imghei=400, zoom_view=1, xyz_scales=c(1.2,2,1.2), legend_xscale=1, cex=1.0, xcex=1.1,
                  theta = 0, phi =0, fov = 60, lighteffect=F, Ymax=NULL) {

   if (firstColumnAsLegend) {
     tmp2 = datMatrix[-1,]; colnames(tmp2)=datMatrix[1,]
  
     xmatrixMj =  matrix(as.integer(tmp2), ncol=dim(tmp2)[2])
     rownames(xmatrixMj) <- rownames(tmp2)
     colnames(xmatrixMj) <- colnames(tmp2)
   }else {
     xmatrixMj = datMatrix
   }

   ncols  = dim(xmatrixMj)[2]; nrows = dim(xmatrixMj)[1]
   groups = colnames(xmatrixMj) # group names
   subs   = rownames(xmatrixMj) # sub group names

   myspace=0.25

   # find xcoord for labeling 
   interv = 1; space4grp = (interv+nrows)
   xcoord =  seq(0, space4grp*ncols, space4grp)
   gene_max = apply(xmatrixMj, 2, max)
   which_min= as.integer(which.min(gene_max))
   ymax = (as.integer(max(xmatrixMj)/20)+1)*20; #max(ymeans);
   ymax = max(xmatrixMj); ymin = 0
   legendx = (which_min-1)*space4grp+legendxDelta;

   group_xcoord = seq(nrows/2,((nrows+1)*ncols)*(1+myspace), (nrows+1)*(1+myspace) )

  # add in space between groups
  datMatrix2 = rbind(xmatrixMj, rep(0, ncols))

  if(is.null(Ymax)) {
     vmax = as.integer((max(xmatrixMj)/10 + 0.99))*10
  } else{
     vmax = Ymax
  }

  # turn into 1D vector for drawing
  frec_org<- as.numeric(datMatrix2)

  frec<- 10*as.numeric(datMatrix2)/vmax
  m   <-length(frec)

   legendx = findMaxValleyInProfile(c(frec_org,max(frec_org)), xpos=c(1:(length(frec_org)+1)), ratio=0.5)
   legendx =(legendx + legendxDelta) * (1+myspace)
   legendy = 0.98*max(frec) + legendyDelta

  # set up color
  rowcol2= c(rowcol[1:nrows], "white")
  color = rep(rowcol2, ncols)

  #distancia entre las barras
  dist<-c(myspace)

  E<-matrix(rep(0,m*24),ncol=3,nrow=8*m)

 #Base de la primera barra
 #    Y
 #    ^
 #    |
 #    5-----8
 #   /|    /|
 #  / |   / |
 # 6-----7  |
 # |  |  |  |
 # |  1--|--4->X
 # | /   | /
 # |/    |/
 # 2-----3
 #/
 #Z

  E[1,]<-c(0,0,0)
  E[2,]<-c(0,0,1)
  E[3,]<-c(1,0,1)
  E[4,]<-c(1,0,0)

  #Altura de la primera   
  E[5,]<-c(0,frec[1],0)
  E[6,]<-c(0,frec[1],1)
  E[7,]<-c(1,frec[1],1)
  E[8,]<-c(1,frec[1],0)


 if(m>1){
  j<-1;
  for(i in seq(1,8*(m-1)-7,8)){

    #Bases del resto de barras
    E[i+8,]<-E[1,]+j*c(1+dist,0,0)
    E[i+9,]<-E[2,]+j*c(1+dist,0,0)
    E[i+10,]<-E[3,]+j*c(1+dist,0,0)
    E[i+11,]<-E[4,]+j*c(1+dist,0,0)

    #Alturas del resto de barras
    E[i+12,]<-c(E[i+8,1],frec[j+1],E[i+8,3])
    E[i+13,]<-c(E[i+9,1],frec[j+1],E[i+9,3])
    E[i+14,]<-c(E[i+10,1],frec[j+1],E[i+10,3])
    E[i+15,]<-c(E[i+11,1],frec[j+1],E[i+11,3])
    j<-j+1
  }
 }

 scales=c(1,1,1); ys=1 # scale up Y axis
 max_x = max(E[,1]); xyratio = max_x/10
 if(max_x>20 & is.na(zoom_view)) {ys=2
 } else {ys=zoom_view}

 rgl.open();
 par3d(windowRect=c(0,0,imgwid,imghei), cex=cex, no.readonly=F, scale=xyz_scales )

 #limpia la ventana antes de dibujar el gr?fico
 rgl.clear(type=c("shapes", "lights"))
 
 rgl.bg(color='white')
 rgl.viewpoint(theta = theta, phi =phi, fov = fov, interactive = F, zoom=1/ys)

 # set up the lighting source
 if(lighteffect){
    lid = rgl.light( theta = 0, phi =0, viewpoint.rel = T, ambient = "#FFFFFF", diffuse = "#FFFFFF", specular = "#FFFFFF")
 } else{
   lid = rgl.light( theta = 0, phi =0, viewpoint.rel = F, ambient = "#111111", diffuse = "#FFFFFF", specular = "#111111")
 }

 # bbox=c(0,imgwid,0,imghei, 0, 200),

 #color<-rainbow(m)
 n<-1

 for(t in seq(1,8*m-7,8)){

   if (n%%(nrows+1)==0){n<-n+1; next} # skip the space between groups

   #pinta la base y la tapa
   rgl.quads(E[t:(t+3),],col=color[n])
   rgl.quads(E[(t+4):(t+7),],col=color[n])

   #pinta la cara izquierda y la derecha
   rgl.quads(rbind(E[t,],E[t+1,],E[t+5,],E[t+4,]),col=color[n])
   rgl.quads(rbind(E[t+3,],E[t+2,],E[t+6,],E[t+7,]),col=color[n])

   #pinta la cara frontal y la trasera
   rgl.quads(rbind(E[t+1,],E[t+2,],E[t+6,],E[t+5,]),col=color[n])
   rgl.quads(rbind(E[t,],E[t+3,],E[t+7,],E[t+4,]),col=color[n])

   if(show_value) {
     text3d(E[t+5,]+c(0.5,0.5,0), text=frec_org[n], family="mono", font=2, cex=0.9, col="gray")
   }

   n<-n+1
 }

 axiscol = "gray2"
 axis3d(edge='y-',col="gray", pos=c(-0.8,0,0.5), nticks=6, labels=as.integer(seq(0,10,2)*vmax/10))
 #axis3d(edge='y+',col="black",pos=c(-0.8,0,0.5))
 #if(m>5){ axis3d(edge='y+',col="black",pos=c(m*(1+dist),0,0.5))}
 title3d(main=title,col="black")

 # find multiple lines in a name
 ml = getStringLength(mystring=groups[1], separator="\n")
 if(ml>1) {
   parts = getAllParts(groups, "\n")
 }
 # plot group names along x axis
 for(j in c(1:ncols) ) {
    #mtext3d(text=groups[j], edge='x-', at = NULL, pos = c(0,-0.8, 0.5))
    #rgl.texts(c(group_xcoord[j], -0.8, 0.5), text=groups[j], cex=2)
    if(ml>1) {
       text3d(c(group_xcoord[j], -0.5, 0.5), text=parts[j,1], family="mono", font=2, cex=xcex, col=axiscol, adj=c(0.5,0.5))
       text3d(c(group_xcoord[j], -1.1, 0.5), text=parts[j,2], family="mono", font=2, cex=xcex, col=axiscol, adj=c(0.5,0.5))
    } else {
       text3d(c(group_xcoord[j], -0.8, 0.5), text=groups[j],  family="mono", font=2, cex=xcex, col=axiscol, adj=c(0.5,0.5))
    }
 }

 # ------ legend ------------------
 legend3D(texts=subs, x=legendx, y=legendy, color=rowcol[1:nrows], cex=1.1, delt = 0.5, scale=legend_xscale) 

 ext = getFileExtension(imgfile)
 if(ext=="png") {
   rgl.snapshot(imgfile)
 } else {
   rgl.postscript(imgfile,"pdf",drawText=T)
 }

 rgl.pop(type="lights", lid)
 rgl.close()

}

 # ------ draw 3D legend ------------------
 #
 #    Y
 #    ^
 #    |
 #    *-----*
 #   /|    /|
 #  / |   / |
 # 4-----3  |
 # |  |  |  |
 # |  *--|--*->X
 # | /   | /
 # |/    |/
 # 1-----2
 #/
 #Z

legend3D = function(texts, x=0, y=0, color="blue", cex=1, delt = 0.5, scale_y=1) {

  no.items = length(texts)
  LE <- matrix(0,4,3)
  for(j in c(1:no.items)) {
    jadj = delt*1.5*(j-1) # adjustment on Y axis, downwards
    LE[1,]<-c(x+0,   y-delt-jadj,1)
    LE[2,]<-c(x+delt*scale_y,y-delt-jadj,1)
    LE[3,]<-c(x+delt*scale_y,y-jadj,     1)
    LE[4,]<-c(x+0,   y-jadj,     1)
    rgl.quads(LE, col=color[j])
    txtE = (LE[2,] + LE[3,])/2 + c(delt/2,-delt/4,0)
    text3d(txtE, text=texts[j], family="mono", font=2, cex=cex, adj=c(0,0), col="black")  
  }

  sublengths = getStringLength(mystring=texts, separator=""); maxlen = max(sublengths)
  OLE <- matrix(0,4,3); recthei = delt*1.5*(no.items-1) + delt*1.5; rectwid = (delt*3+maxlen*0.5)*scale_y
  OLE[1,]<-c(x-delt,   y-recthei, 1);  OLE[2,]<- OLE[1,]+c(rectwid, 0, 0)
  OLE[4,]<-c(x-delt,   y+delt/2,  1); OLE[3,]<- OLE[4,]+c(rectwid, 0, 0)
  rgl.lines(OLE, col="red"); rgl.lines(OLE[c(1,4,2,3),], col="red")
}


# -------------------- barplot groups (2D) -------------------------------
#
barplot_by_group = function(datMatrix, maskMatrix=NULL, firstColumnAsLegend=F, maskAsSD=FALSE,
                  rowcol=c("red", "green", "blue", "grey"), vlines=NULL,
                  title="", ylab="", imgfile=NULL, imgwid=800, imghei=400, 
                  ipointsize=12, iunits="px", ires=72, icompression="lzw",
                  mcex=1.6,legendcols=1, legendxDelta=0, legendyDelta=0, legendCex=0.8, display_exponent=F,
                  show_value=T, showLegend=T, show_xaxisLabel=T,
                  xmargin=4, ymargin=4, xlabel_orientation=0, ValueXdelta=0, 
                  ValueYdelta=NULL, valuecex=0.8, yMAX=NULL, yMIN=NULL, ymax_step=20,
                  headlabel_deltX=0.5, xcoord_delta=0) {

   if (firstColumnAsLegend) {
     tmp2 = datMatrix[-1,]; colnames(tmp2)=datMatrix[1,]
     xmatrixMj =  matrix(as.integer(tmp2), ncol=dim(tmp2)[2])
     rownames(xmatrixMj) <- rownames(tmp2)
     colnames(xmatrixMj) <- colnames(tmp2)

     if(!is.null(maskMatrix)){
        tmp2 = maskMatrix[-1,]; colnames(tmp2)=maskMatrix[1,]
        xmaskMatrix =  matrix(as.integer(tmp2), ncol=dim(tmp2)[2])
        rownames(xmaskMatrix) <- rownames(tmp2)
        colnames(xmaskMatrix) <- colnames(tmp2)
     }

   }else {
     xmatrixMj   = datMatrix
     xmaskMatrix = maskMatrix
   }

   # turn into 1D vector for drawing
   ncols      = dim(xmatrixMj)[2]; nrows = dim(xmatrixMj)[1]
   datMatrix2 = rbind(xmatrixMj, rep(0, ncols))
   frec_org   = as.numeric(datMatrix2); vmax = max(frec_org)
   no.bars    = length(frec_org)
   bar_x = c(1:no.bars); bar_y = ifelse(frec_org==0, "", frec_org)

   if(!is.null(maskMatrix)){
      xmaskMatrix2 = rbind(xmaskMatrix , rep("", ncols))
      mask_org   = as.character(xmaskMatrix2);
      bar_y = paste(mask_org, "\n", bar_y, sep="")
   }

   tab_names = rownames(xmatrixMj)
   no.groupsX = dim(xmatrixMj)[2]; no.tabs=dim(xmatrixMj)[1]

   if(no.tabs>1) {
      interv = 1; space4grp = (interv+no.tabs)
      xcoord =  seq(2.5, space4grp*no.groupsX, space4grp) + no.groupsX%%2-0.5 + xcoord_delta
      bar_x  = bar_x +headlabel_deltX
   } else{
      interv = 1; space4grp = interv+1
      xcoord =  seq(1.5, space4grp*no.groupsX, space4grp)+no.groupsX%%2-0.5  + xcoord_delta
      bar_x  = bar_x -headlabel_deltX
   }
 
   gene_max = apply(xmatrixMj, 2, max)
   which_min= as.integer(which.min(gene_max))

   if( is.null(yMAX) ) {
      ymax = (as.integer(max(xmatrixMj)/200)+1)*200; #max(ymeans);
      ymax = (as.integer(max(xmatrixMj)/ymax_step)+1)*ymax_step; #max(ymeans);
   } else{
      ymax = yMAX
   }

   if(is.null(yMIN) ) {
     ymin=0 #min(xmatrix) #-0.2
   } else {ymin=yMIN}

   legendx = (which_min-1)*space4grp #findMaxValleyInProfile(meanPermodule[,1], xpos=c(1:no.modules))
   #legendx = findMaxValleyInProfile(as.vector(xmatrixMj), xpos=c(1:no.modules))

   legendx = findMaxValleyInProfile(c(frec_org,max(frec_org)), xpos=c(1:(length(frec_org)+1)), ratio=0.5) +legendxDelta

   #legendx = (which_min-1)*space4grp+1+legendxDelta
   legendy = ymax-0.02*ymax + legendyDelta

   mycolor = rowcol[c(1:no.tabs)]

   if(!is.null(imgfile)) {
     openImgDev(imgfile, iwidth =imgwid, iheight = imghei, ipointsize=ipointsize, iunits=iunits, ires=ires, icompression=icompression)
     par(mfrow=c(1,1), mar=c(xmargin, ymargin, 1, 0) + 0.1, cex=mcex,las=1)#, srt=90)
   }

   barplot(xmatrixMj, beside = TRUE, #space=interv,
        col = mycolor, axisnames=F, ylab=ylab, ylim = c(ymin, ymax),
        legend =F)#shortnames) #legend =F
   title(main = title, font.main = 1)

   if(!is.null(maskMatrix) & maskAsSD) {
      err.bp(xmatrixMj, maskMatrix, two.side=TRUE) 
   }

   if(show_xaxisLabel) {
       axis(1, at =xcoord, labels =colnames(xmatrixMj), las=xlabel_orientation, 
               col = "black", col.axis="black", lwd = 0, tick = T, line =1)
   }

   if((nrows>1) & showLegend) {
     legend(legendx+1, legendy, legend=tab_names, box.lty=1,box.lwd=1,box.col=mycolor, ncol =legendcols, 
          fill=mycolor, cex=legendCex,pt.cex=legendCex)
   }

   ValueYdelta2 = ValueYdelta
   if(is.null(ValueYdelta)){
      ValueYdelta2 = ymax/20
   }

   if(show_value) {
      if(display_exponent) {
          text(x=bar_x, y=frec_org+ValueYdelta2, labels=bar_y, cex=valuecex, col="black")
      } else {
          text(x=bar_x, y=frec_org+ValueYdelta2, labels=bar_y, cex=valuecex, col="black")
      }
   }

   if(!is.null(vlines)){
        for(xx in vlines){abline(h=xx,col="gray", lty=1)}
   }

   if(!is.null(imgfile)) {
     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
     dev.off()
   }

}

# -------------------- barplot groups horizon (2D) -------------------------------
#
#
barplot_by_group_Horizon = function(datMatrix, firstColumnAsLegend=F, 
                  rowcol=c("red", "green", "blue", "grey"), vlines=NULL,
                  title="", ylab="", imgfile=NULL, imgwid=800, imghei=400, 
                  ipointsize=12, iunits="px", ires=72, icompression="lzw",
                  mcex=1.6,legendcols=1, legendxDelta=0, show_value=T, axiscolor="greenyellow", 
                  xmargin=4, ymargin=6, xlabel_orientation=0, ValueYdelta=0.6, yMAX=NULL, plotAxes=TRUE) {

   if (firstColumnAsLegend) {
     tmp2 = datMatrix[-1,]; colnames(tmp2)=datMatrix[1,]
  
     xmatrixMj =  matrix(as.integer(tmp2), ncol=dim(tmp2)[2])
     rownames(xmatrixMj) <- rownames(tmp2)
     colnames(xmatrixMj) <- colnames(tmp2)
   }else {
     xmatrixMj = datMatrix
   }

   # turn into 1D vector for drawing
   ncols      = dim(xmatrixMj)[2]; nrows = dim(xmatrixMj)[1]
   datMatrix2 = rbind(xmatrixMj, rep(0, ncols))
   frec_org   = as.numeric(datMatrix2); vmax = max(frec_org)
   no.bars    = length(frec_org)
   bar_x = c(1:no.bars); bar_y = ifelse(frec_org==0, "", frec_org)

   tab_names = rownames(xmatrixMj)
   no.groupsX = dim(xmatrixMj)[2]; no.tabs=dim(xmatrixMj)[1]

   if(no.tabs>1) {
      interv = 1; space4grp = (interv+no.tabs)
      xcoord =  seq(2.5, space4grp*no.groupsX, space4grp)
   } else{
      interv = 1; space4grp = interv+1
      xcoord =  seq(1.5, space4grp*no.groupsX, space4grp)
      bar_x  = bar_x -0.5
   }
 
   gene_max = apply(xmatrixMj, 2, max)
   which_min= as.integer(which.min(gene_max))

   if(is.null(yMAX) ) {
      ymax = (as.integer(max(xmatrixMj)/200)+1)*200; #max(ymeans);
      ymax = (as.integer(max(xmatrixMj)/20)+1)*20; #max(ymeans);
   } else{
      ymax = yMAX
   }

   ymin=0 #min(xmatrix) #-0.2
   legendx = (which_min-1)*space4grp #findMaxValleyInProfile(meanPermodule[,1], xpos=c(1:no.modules))
   #legendx = findMaxValleyInProfile(as.vector(xmatrixMj), xpos=c(1:no.modules))

   legendx = findMaxValleyInProfile(c(frec_org,max(frec_org)), xpos=c(1:(length(frec_org)+1)), ratio=0.5)

   legendx = (which_min-1)*space4grp+1+legendxDelta
   legendy = ymax-0.02*ymax

   mycolor = rowcol[c(1:no.tabs)]

   openImgDev(imgfile, iwidth =imgwid, iheight = imghei, ipointsize=ipointsize, iunits=iunits, ires=ires, icompression=icompression)
 
   par(mfrow=c(1,1), mar=c(xmargin, ymargin, 1, 0) + 0.1, cex=mcex,las=2)#, srt=90)
   barplot(xmatrixMj, beside = TRUE, #space=interv,
        col = mycolor, border=mycolor, axisnames=F, xlab=ylab, xlim = c(ymin, ymax),
        legend =F, horiz=TRUE, axes=plotAxes)#shortnames) #legend =F
   title(main = title, font.main = 1)
   #err.bp(t(ymeans), t(ysds), two.side=F) 
   axis(2, at =xcoord, labels =colnames(xmatrixMj), las=xlabel_orientation, 
               col ="black", col.axis="black", lwd = 0, tick = T, line =1)

   if(nrows>1) {
     legend(legendx+1, legendy, legend=tab_names, box.lty=1,box.lwd=1,box.col=mycolor, ncol =legendcols, 
          fill=mycolor, cex=0.8,pt.cex=0.8)
   }

   if(show_value) {
      text(y=bar_x+ValueYdelta, x=frec_org+vmax/25, labels=bar_y, cex=0.8, col="gray", srt = 90)

   }

   if(!is.null(vlines)){
        for(xx in vlines){abline(h=xx,col="gray", lty=3)}
   }

   par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
   dev.off()

}


# xaxisLabel_rotate = 45 # slant x-labels
#
barplot_by_oneGroup = function(datMatrix, maskMatrix=NULL, 
                  rowcol=c("red", "green", "blue", "grey"), vlines=NULL,
                  title="", ylab="", imgfile=NULL, imgwid=800, imghei=400, 
                  mcex=1.6, xlab_cex=1, xaxisLabel_rotate=NULL, usr_delta=2, 
                  legendcols=1, legendxDelta=0, legendCex=0.8, 
                  ipointsize=12, iunits="px", ires=72, icompression="lzw",
                  show_value=T, showLegend=T, show_xaxisLabel=T, show_mask=T,
                  xmargin=4, ymargin=4, zmargin=1, xlabel_orientation=0, ValueXdelta=0, 
                  ValueYdelta=NULL, valuecex=0.8,value_rotate=0, yMAX=NULL, yMIN=0) {

   xmatrixMj   = datMatrix
   xmaskMatrix = maskMatrix
 
   if(!is.null(maskMatrix)){
      bar_y = paste(datMatrix[,1], maskMatrix[,1] , sep="")
   }

   tab_names  = rownames(xmatrixMj)
   no.groupsX = dim(xmatrixMj)[1]; no.tabs=dim(xmatrixMj)[2]

      interv = 1; space4grp = interv+1
      xcoord =  0.5+ seq(1, space4grp*no.groupsX, space4grp)
      
   gene_max = max(xmatrixMj)
   which_min= as.integer(which.min(gene_max))

   if( is.null(yMAX) ) {
      ymax = (as.integer(max(xmatrixMj)/200)+1)*200; #max(ymeans);
      ymax = (as.integer(max(xmatrixMj)/20)+1)*20; #max(ymeans);
   } else{
      ymax = yMAX
   }

   if(is.null(yMIN) ) {ymin=min(xmatrixMj)
   } else {ymin=yMIN}

   mycolor = rowcol[c(1:no.groupsX)]

   if(!is.null(imgfile)) {
     openImgDev(imgfile, iwidth =imgwid, iheight = imghei, ipointsize=ipointsize, iunits=iunits, ires=ires, icompression=icompression)
     par(mfrow=c(1,1), mar=c(xmargin, ymargin, zmargin, 0) + 0.1, cex=mcex,las=1)#, srt=90)
   }

   barplot(datMatrix[,1], beside = TRUE, space=interv,
        col = mycolor, axisnames=F, ylab=ylab, ylim = c(ymin, ymax),
        legend =F)#shortnames) #legend =F

   if(title!="") {
       title(main = title, font.main = 1)
   }

   #err.bp(t(ymeans), t(ysds), two.side=F) 
   if(show_xaxisLabel) {
       if(is.null(xaxisLabel_rotate) ) {
          axis(1, at =xcoord, labels =rownames(xmatrixMj), las=xlabel_orientation, 
               col = "black", col.axis="black", lwd = 0, tick = T, line =1)
       } else{
          axis(1, at =xcoord, labels =rep("", nrow(datMatrix)), las=xlabel_orientation, 
               col = "black", col.axis="black", lwd = 0, tick = T, line =1)
          text(xcoord, par("usr")[3] - usr_delta, srt = xaxisLabel_rotate, adj = 1, labels = rownames(xmatrixMj), 
          xpd = TRUE, cex=xlab_cex)
          par(srt=0) 
       }
   }

   ValueYdelta2 = ValueYdelta
   if(is.null(ValueYdelta)){
      ValueYdelta2 = ymax/20
   }

   if(show_value) {
      text(x=xcoord+ValueXdelta, y=datMatrix[,1]+ValueYdelta2, srt=value_rotate,adj = 1,labels=bar_y, cex=valuecex, xpd = TRUE, col="black")
   }

   if(show_mask) {
      #text(x=xcoord, y=datMatrix[,1]+ValueYdelta2, srt=value_rotate,adj = 1,labels=maskMatrix[,1], xpd = TRUE, cex=valuecex, col="black")
      text(x=xcoord+ValueXdelta, y=datMatrix[,1]+ValueYdelta2, srt=value_rotate,adj = 1,labels=maskMatrix[,1], xpd = TRUE, cex=valuecex, col="black")
   }


   if(!is.null(vlines)){
        for(xx in vlines){abline(h=xx,col="gray", lty=3)}
   }

   if(!is.null(imgfile)) {
     par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1, srt=0)
     dev.off()
   }

}


read_weka_features = function(fweka) {

   nskips = 2

   title <- NULL
   while(T)  {
      allMatrix <- read.delim(fweka,sep=" ", header=F, nrow=1, 
                          skip=nskips)
      allstr    <- as.matrix(allMatrix)[1,]

      #print(allstr)

     if(allstr[1]=="@ATTRIBUTE"){
        splitted  = splitString(allstr[2], "\t")
        title     = c(title, splitted[1]) 
     }else if ( allstr[1]=="@DATA"){
        nskips = nskips + 1
        break
     }
     nskips = nskips + 1
   }
   nskips = nskips + 1 # one blank lines

   #*-------------------------------------------------------------------------------------
   #* STEP 1. read feature values (column as feature, row as sample)
   #
   allMatrix <- read.delim(fweka,sep=",", header=F, skip=nskips)
   allMatrix <- as.matrix(allMatrix)
   allMatrix <- ifelse(allMatrix=="?", NA, allMatrix)

   colnames(allMatrix) <- title

   return (allMatrix)

}

SetTextContrastColor <- function(color)
{
  ifelse( mean(col2rgb(color)) > 127, "black", "white")
}


copyFiles = function(fromDir, toDir, keylabels=NULL, overwrite=FALSE, recursive=FALSE) {

   allfiles  = NULL
   if(!is.null(keylabels)) {
      for (ek in keylabels) {
        allfiles = union(allfiles, list.files(fromDir, ek) )
      }
   } else{
      allfiles = list.files(fromDir, "")
   }

   srcfiles = paste(fromDir, allfiles, sep="")
   dstfiles = paste(toDir,   allfiles, sep="")

   res  = file.copy(from=srcfiles, to=dstfiles, overwrite = overwrite, recursive = recursive)

   res1 = cbind(allfiles,res)

   #print(res1)
   return(res1)
}

renameFiles = function(oldfiles, oldlabel="", newlabel="") {

   
   for (each in oldfiles) {

      dstfile = replaceString(each, oldstr=oldlabel, newstr=newlabel)
      res  = file.rename(from=each, to=dstfile)      
      print( paste(each, "-->", dstfile, sep="") )
      print(res)
   }
   #print(res1)
   #return(res1)
}

# get all pair-wise combinations
#
getPairs = function(memIdxOrg)
{

  memIdx = unique(memIdxOrg)
  #print(memIdx)
  
  nmem = length(memIdx)
  if(nmem==1) {return (NULL)}

  index= NULL
  for(i in c(1:(nmem-1)) ){
     index= rbind(index, cbind(rep(i,nmem-i), (i+1):nmem)   )
  }

  indx = cbind(memIdx[index[,1]], memIdx[index[,2]])

  #print(indx)
  
  return (indx)
}


upper.tri.index=function(imatrix){
  nr = dim(imatrix)[1]
  nc = dim(imatrix)[2]
  last = ifelse(nr<nc, nr, nc-1)
  index= NULL
  for(i in c(1:last) ){
     index= rbind(index, cbind(rep(i,nc-i), (i+1):nc)   )
  }
  return (index)
}

lower.tri.index=function(imatrix){
  index=upper.tri.index(imatrix)
  return (index[,c(2,1)])
}


#
# remove leading and trailing white spaces
#  line can be a vector
#
trim <- function(line) gsub("(^ +)|( +$)", "", line)  


## "Mixed Case" Capitalizing - toupper( every first letter of a word ) :
simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}
#simpleCap("the quick red fox jumps over the lazy brown dog")
## ->  [1] "The Quick Red Fox Jumps Over The Lazy Brown Dog"

## and the better, more sophisticated version:
capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s,1,1)),
                  {s <- substring(s,2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
#capwords


# log ratio pavlue in Agilent platform
#
#rgSigErr = c(rSig, gSig, rError, gError)
#
computePValueEAgilentrrorModel= function(rgSigErr){

    xdv  = abs(rgSigErr[1]-rgSigErr[2])/(rgSigErr[3]^2+rgSigErr[4]^2)^0.5
    pval = pnorm(xdv/2^0.5,mean=0,sd=1, lower.tail=FALSE)
    return (pval)
}

#Input "ABHD2"    "IL6ST" "0.254035"    "KLF13" "0.433281"
#Output: ABHD2  IL6ST
#        ABHD2  KLF13
#
read_aracne_res = function(reslineSTR)
{
   resline = getAllParts(reslineSTR, "\t")
   sel = resline !=""
   resl= resline[sel]
   #print(resline)
   #print("yyy\n")
   idx = seq(from=2, to=sum(sel), by=2)
   #inet = cbind(rep(resl[1], length(idx)), resl[idx])
   inet = paste(resl[1], resl[idx], sep="\t")
   return (inet)
}




#----------------------------------------------------------------------------------------------
#
# Two Y-axes representing two different scales
#
#    pch=19: solid circle,
#    pch=20: bullet (smaller solid circle, 2/3 the size of 19),
#    pch=21: filled circle,
#    pch=22: filled square,
#    pch=23: filled diamond,
#    pch=24: filled triangle point-up,
#    pch=25: filled triangle point down. 

PlotTimeSeriesTwoYaxes = function(data, xlabel_ypos=0, xtitle="", yindex, y1title, y1labels, y2title, y2labels=NULL, y2PosDelta=0, linecolors, 
                           maintitle="", xlabel_rotate=45, xcex=1, plot_legend=TRUE, legendCex=0.6*xcex, xlabelCex=NULL, xpch=21) {

  npoints = dim(data)[1]
  xindex = c(1:npoints)

   
  if(npoints<=20){cexX = xcex
  } else {cexX=xcex*(20/npoints)^0.5}

  if(!is.null(xlabelCex)){
     cexX = xlabelCex
  }

  xindex = c(1:npoints);
  matplot(x=xindex, y=data, type='l', lty="solid", col = linecolors, lwd=1,
                xlab=xtitle, ylab=y1title, main=maintitle, axes=F, pch=xpch, lend="square")#, frame.plot=TRUE)

  par(srt=xlabel_rotate, cex=cexX)
  axis(1, at =xindex, labels =rep("", npoints), las=2, pos=-0.25,
               col = "black", col.axis="black", lwd = 1, tick = T, line =1, cex=cexX)

  text(xindex+0.1, par("usr")[3] - 2-xlabel_ypos, srt =xlabel_rotate, adj = 1, labels = rownames(data), xpd = TRUE, cex=cexX)
  par(srt=0, cex=xcex) 

  axis(2, at =yindex, labels =y1labels, las=2, pos=0.75,
               col = linecolors[2], col.axis=linecolors[2], lwd = 1, tick = T, line =1)

  if(!is.null(y2labels)) {
  axis(4, at =yindex, labels =y2labels, las=2, pos=npoints+y2PosDelta,
               col = linecolors[1], col.axis=linecolors[1], lwd = 1, tick = T, line =1)
  }

  points(x=xindex, y=data[,1], pch=21, col=linecolors[1], bg=linecolors[1])
  points(x=xindex, y=data[,2], pch=21, col=linecolors[2], bg=linecolors[2])

  mtext(text=y2title, line=2, side = 4)

  #for(j in yseq){
  #  abline(h=j, col="gray", lty="dotted")
  #}
  #text(x=xindex, y=posts, label="*", cex=1.8, col="brown")

  if(F) {
  #delt = (max(posts)-min(posts))/10
  deltx = 0.4/npoints; delty=step/10
  for(ii in c(1:length(xindex)) ) {
    for (jj in c(1:nvariabs) ) {
      ijval = fulMtrx2[ii,jj]
      if(ijval==0) {next}
      rect(xleft=xindex[ii]-deltx, ybottom=ijval-delty,xright=xindex[ii]+deltx,ytop=ijval+delty, col=colcode[jj], border=colcode[jj], lwd=1, lty=1)
    }
  }
  }

 if(plot_legend) {
   mprof = apply(data, 1, max)
   mprofMiddle = ifelse(mprof>2*max(mprof)/3, 2*max(mprof)/3, mprof)

   legendx = findMaxValleyInProfile(mprofMiddle, xindex)
   legend(legendx, max(data), legend=colnames(data),lty=1,lwd=2, col=linecolors, ncol =1, cex=legendCex,pt.cex=0.8)
 }

}

if(F) {
    #To rotate axis labels (using base graphics), you need to use text(), rather than mtext(), as the latter does not support par("srt"). 

     ## Increase bottom margin to make room for rotated labels
     par(mar = c(7, 4, 4, 2) + 0.1)
     ## Create plot with no x axis and no x axis label
     plot(1 : 8, xaxt = "n",  xlab = "")
     ## Set up x axis with tick marks alone
     axis(1, labels = FALSE)
     ## Create some text labels
     labels <- paste("Label", 1:8, sep = " ")
     ## Plot x axis labels at default tick marks
     text(1:8, par("usr")[3] - 0.25, srt = 45, adj = 1,
          labels = labels, xpd = TRUE)
     ## Plot x axis label at line 6 (of 7)
     mtext(1, text = "X Axis Label", line = 6)

    #When plotting the x axis labels, we use srt = 45 for text rotation angle, adj 
}


compute.robust.residuals1<-function(m,batch, add_mean=T){
	lm.obj<-lm(m~as.factor(batch),na.action="na.exclude",singular.ok = T)
	resid.data<-residuals(lm.obj)
        if (add_mean) {resid.data=resid.data -mean(resid.data, na.rm=T)  + mean(m, na.rm=T)}
        return (resid.data)
}

# Input
#    1. exprDat is a matrix with rows as genes and columns as samples
#    2. samplebatch is the batch IDs for all the samples
# Output 
#    expression data with the batch effect corrected

correct_batch <- function(exprDat, samplebatch) {
    exprDat2 <-apply(exprDat,1,compute.robust.residuals1,
                     batch= samplebatch, add_mean=T)
    exprDat2 <- t(exprDat2)

    return (exprDat2)
}


########################### Group TTest #######################################

# global variables: datExpr, mask
# here groupA and groupB are the column indices
groupsTtest=function(groupA, groupB, identDiff=T, minpvalue=0.05, foldchange=1, mindifference=0){

           foldchangeInv=1/foldchange

           no.samplesA = length(groupA)
           no.samplesB = length(groupB)

           exprA = datExpr[, groupA]
           exprB = datExpr[, groupB]
           dim(exprA)
           dim(exprB)

           #mask information
           maskA = mask[, groupA]
           maskB = mask[, groupB]

           if (no.samplesA > 1){
              presentnoA = apply(maskA, 1, sum, na.rm=TRUE)
              #by gene (row)
              meanA = apply(exprA, 1, mean, na.rm=TRUE)
              #sdA  = apply(exprA, 1, sd, na.rm=TRUE)
           }else{ #only one sample in the group
              presentnoA = as.integer(maskA+0)
              meanA      = as.numeric(exprA)
           }

           if (no.samplesB > 1){
              presentnoB = apply(maskB, 1, sum, na.rm=TRUE)
              meanB      = apply(exprB, 1, mean, na.rm=TRUE)
              #sdB  = apply(exprB, 1, sd, na.rm=TRUE)
           }else{ #only one sample in the group
              presentnoB = as.integer(maskB+0)
              meanB      = as.numeric(exprB)
           }

           folds    = meanA/meanB
           meandiff = meanA-meanB

           #ttest pvalues
           pvalues  = apply(datExpr, 1, ttestVector, tgroupA=groupA, tgroupB=groupB)
           pvalues  = ifelse( (presentnoA==0) & (presentnoB==0), 1, pvalues)
           qval     = qvalue(pvalues)
           qvalues  = qval$qvalues

          if(identDiff){
           #select genes by ttest p-value and fold-change
           selectedByPval = qvalues      <= minpvalue
           selectedByFold = ( abs(folds) >= foldchange | abs(folds) <= foldchangeInv )
           selectedByDiff = abs(meandiff)>=  mindifference

          }else{ #similarly expressed genes
           selectedByPval = qvalues      >= minpvalue
           selectedByFold = ( abs(folds) <= foldchange & abs(folds) >= foldchangeInv )
           selectedByDiff = abs(meandiff) <=  mindifference
          }

           #we need change NA into False, as True & NA = TRUE
           selectedByPval = ifelse(is.na(selectedByPval), FALSE, selectedByPval)
           selectedByFold = ifelse(is.na(selectedByFold), FALSE, selectedByFold)
           selectedByDiff = ifelse(is.na(selectedByDiff), FALSE, selectedByDiff)

          if(identDiff){
           selected    = selectedByPval & selectedByFold & selectedByDiff
          }else{
           selected    = selectedByPval & selectedByDiff
          }

          #no.selected = sum(selected[ !is.na(selected) ])
          #selected10 = ifelse(selected, 1, 0)

      res = cbind(selected, meanA, meanB, folds, meandiff, pvalues, qvalues)
      return (res)

      #return (selected)

}


tt = function(vect){
   nsize = length(vect)
   res=t.test(vect[1:(nsize/2)], vect[(nsize/2+1):nsize])
   res$p.value
}

cortest = function(v1, v2){return (cor.test(v1,v2)$p.value) }
corpV2M = function(vect, mtrx){
   res=apply(mtrx, 1, cortest, v2=vect)
   return (rbind(res))
}


#corpval_pc2trait  = apply(pc1st, 2, corpval_multitraits, traitDat)
corpM2M = function(mtrx1, mtrx2){
   res=apply(mtrx1, 2, corpV2M, mtrx=t(mtrx2))
   return (res)
}


closest_integer_bound = function(val){
  logv = log10(val)
  ii   = as.integer(logv)
  idec = val/(10^ii) # 2.50
  iremaind = idec - as.integer(idec)

  ix = ifelse(iremaind<0.5, 0.5, 1)
  return ( (as.integer(idec)+ix)*(10^ii) )
}


# returns the content of gc() (space report) and proc.time() (elapsed time report)
# as a single vector
# cum is a cumulative data frame of results.  If not null the answer will be added
# as a new row and the entire frame returned.  Thus one may call cum<-spaceTime(cum)
# to build up a data frame of results
spaceTime<-function(label, cum=NULL) {
	space <-gc()
	time <- proc.time()
	spaceAndTime <- t(as.matrix(c(space[1,], space[2,], time)))
	rownames(spaceAndTime)<-label
	if (is.null(cum)) {
		return(spaceAndTime)
	} else {
		return(rbind(cum, spaceAndTime))
	}
}

spaceTimePlot<-function(spaceTimeStats) {
	elapsed<-spaceTimeStats[,"elapsed"]
	delta<-elapsed[2:length(elapsed)]-elapsed[1:length(elapsed)-1]
	opar <- par(mfrow=c(3,1))
	pie(delta, col=rainbow(length(delta)), mfrow=c(2,1), new=T)
	title(paste("Breakdown of ", secToHMS(elapsed[length(elapsed)]-elapsed[1])))
	
	mem<-spaceTimeStats[,c(2,4,6)]+spaceTimeStats[,c(8,10,12)]
	colnames(mem)<-c("used (mb)", "GC trigger (mb)", "max used (mb)")
	
	plotMemCol(mem, 1, isNew=F)
	plotMemCol(mem, 3, isNew=F)
}

secToHMS<-function(sec) {
	hours = trunc(sec/3600)
	remainder = sec-hours*3600
	minutes = trunc(remainder/60)
	seconds = round(remainder-minutes*60)
	paste(paste(hours, "h", sep=""), paste(minutes, "m", sep=""), paste(seconds, "s", sep=""), sep=":")
}


plotMemCol<-function(mem, col, isNew) {
	# plot 'used (mb)'
	plot(1:(dim(mem)[1]), mem[,col], "s", xlab="step", ylab="mb", new=isNew) # draw lines
	text(1:(dim(mem)[1]), mem[,col], rownames(mem), cex=.8, pos=4, srt=-30)#label
	title(colnames(mem)[col])
}


# mean of subclasses
avgvect = function(mvect, mclass){
   avg = tapply(mvect, INDEX=mclass, FUN=mean, na.rm = TRUE);
   return (rbind(avg))
}

# compute similarity between binary vectors
#
similarity_binary_vv = function(v1, v2, normalize=TRUE){
   if(normalize) {
     vv = sum(v1*v2)/length(v1)
   } else {
     vv = sum(v1*v2) 
   }
   return (vv)
}
similarity_binary_vm = function(v, mtrx, normalize=TRUE){
   vm = apply(mtrx, 2, similarity_binary_vv, v2=v, normalize=normalize)
   return (rbind(vm))
}
similarity_binary = function(mtrx, normalize=TRUE){
   vm = apply(mtrx, 2, similarity_binary_vm, mtrx= mtrx, normalize=normalize)
   return (vm)
}

# multiply all columns
multiply_cols = function(mtrx) {
   ncols = ncol(mtrx)
   res = mtrx[,1]
   for(i in c(1:ncols)) {
     res = res * mtrx[,i]
   }
   return (res)
}

fillInClass = function(mytissues, alltissues){
   myvect = rep(0, length(alltissues)); names(myvect) = alltissues
   myvect[mytissues] = 1
   return (myvect)
}

# collapse a list to vect 
list2vect = function(mlist){
   res = NULL
   for(el in mlist){
      res = c(res, el)
   }
   return (res)
}

list2MatrixModules = function(mlist){
   res = NULL
   mnames = names(mlist)
   for(xx in c(1:length(mlist)) ){
      el  = mlist[[xx]]
      itm = cbind(el, rep(mnames[xx], length(el)) )
      res = rbind(res, itm)
   }

   return (res)
}


# single module
module_conserve= function(mdat, permutes=10, method="pearson"){
   mcor = cor(mdat,  use="complete.obs", method=mehod)
   mlinks = sum(abs(mcor))

   mvect = as.numeric(mdat); n.pts = length(mvect); ndx = c(1:n.pts)
   for(i in c(1:permutes)){
     idat =apply(datExpr,  1, permuteVect)

   }

   
}

pnorm_vect= function(q_mean_sd, lower.tail=FALSE){
  pval = pnorm(q=q_mean_sd[1], mean = q_mean_sd[2], sd = q_mean_sd[3], lower.tail = lower.tail, log.p = FALSE)
  return(pval)
}

modules_conservation= function(datExpr, modulelabel, no.perms=10, method="pearson", beta=1){

  modulenames = unique(modulelabel); modulenames=setdiff(modulenames, "grey")
  modulenames = sort(modulenames)
  no.modules  = length(modulenames)

  kins = matrix(0, nrow=no.modules, ncol=no.perms+1)
  
  for( y in c(1:(no.perms+1)) ) {
    if(y!=(no.perms+1) ){ # permuated values
      XdatExpr =apply(datExpr,  2, permuteVect)
    } else { # true value
      XdatExpr =datExpr
    }

    for ( x in c(1:no.modules) ) {

      #print(paste(x, "/", no.modules, ":", modulenames[x], sep="") )

      xsel = modulelabel== modulenames[x]

      corhelp <- cor(XdatExpr[, xsel],  use = "pairwise.complete.obs", method=method)
      diag(corhelp) <- 0
      corhelp = abs(corhelp)^beta
      #links   <- apply(abs(corhelp), 1, sum, na.rm=TRUE)

      kins[x, y] = mean(corhelp, na.rm=T)
      rm(corhelp)
      collect_garbage()
   }
   rm(XdatExpr)
   collect_garbage()
  }

  means=apply(kins[,-c(no.perms+1)], 1, mean, na.rm=TRUE)
  sds  =apply(kins[,-c(no.perms+1)], 1, sd,   na.rm=TRUE) 
  kinDat= cbind(kins[,no.perms+1], means, sds)
  pvals=apply(kinDat, 1, pnorm_vect, lower.tail=FALSE)

  kinDat2= cbind(pvals, kinDat)
  colnames(kinDat2) = c("Module_Connectivity_Pval", "Module_Connectivity", "Random_Connectivity", "Random_Connectivity_SD")
  rownames(kinDat2) = modulenames 

  return(kinDat2)
}


# text mining functions

# geneLink=http://flybase.org/.bin/fbidq.html? protLink=http://flybase.org/.bin/fbidq.html?
# find gene id and 
#http://may2009.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000197102
# read URL
getUrlPage = function(murl) {
  zz <- url(murl)
  lines = readLines(zz)
  close(zz)
  return (lines)
}

#x='<gene id=14086 geneId=ENSG00000187609 protId=ENSP00000340474/>'
paresGeneID = function(xs){
   xs2 = sub("<gene id=", "", xs)
   xs2 = sub("geneId=", "", xs2)
   xs2 = sub("protId=", "", xs2)
   xs2 = sub("/>", "", xs2)
   return (xs2)
}

# parse orthologGroup
#
#    <orthologGroup id="1">
#      <score id="bit" value="6795"/>
#      <geneRef id="1">
#        <score id="inparalog" value="1"/>
#        <score id="bootstrap" value="1.00"/>
#      </geneRef>
#      <geneRef id="2">
#        <score id="inparalog" value="1"/>
#        <score id="bootstrap" value="1.00"/>
#      </geneRef>
#    </orthologGroup>
#
parseOrtholog = function(xs, scorecut=0.5){
   #x= '<orthologGroup id="1">'
   xs2 = gsub('"', "", xs)
   m=regexpr("([^=]+)([0-9]+)",xs2[1])
   groupid=parse_general(xs2[1], patt='([^=]+)?([0-9]+)') 
   groupid=groupid[1]

   xidx = c(1:length(xs))
   sel = xs2 == "</geneRef>"
   refIdxE = xidx[sel]
   npt = sum(sel)

   if(npt>1) {
     refIdxS=c(3, refIdxE[-npt]+1)
     refINDX=cbind(refIdxS, refIdxE)
   } else{
     refINDX = rbind(c(3, refIdxE))
   }

   ires = NULL
   for(ii in c(1:npt)) {
      ixs = xs2[refINDX[ii,1]:refINDX[ii,2]]
      iref = parse_general(ixs[1], patt='([^=]+)?([0-9]+)') 
      isco = parse_general(ixs[2], patt='([^=]+)?([0-9]+)') 
      ires = rbind(ires, as.integer(c(iref[1], isco[1])))
   }

   # cut off 
   sel = ires[,2]>=scorecut
   nsel = sum(sel)
   ires2= ires[sel, ]
   if(nsel <=1 ) {return (NULL)}

   pairs=NULL
   for(ii in c(1:(nsel-1)) ) {
     for(jj in c((ii+1):nsel) ) {
        pairs = rbind(pairs, c(ires2[ii,1],ires2[jj,1]) )
     }
   }

   pairs2 = cbind(rep(as.integer(groupid), nrow(pairs)), pairs)
   colnames(pairs2) = c("GroupID", "GeneRef_A", "GeneRef_B")

   return (pairs2)
}





#murl= "http://may2009.archive.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000348965"
#murl= "http://flybase.org/reports/FBpp0073215.html"
findGeneSymbol_Fly = function(protein) {
   murl = paste("http://flybase.org/reports/", protein, ".html", sep="")
   postword = "</title>"
   page = getUrlPage(murl)
   page = trim(page)
   sel  = page!=""; page = page[sel]
   sel2 = page==postword
   gs =""
   if(sum(sel2)>0) {
     idx  = c(1:length(page))[sel2]; gs = page[idx-1]
   }
   return(gs)
}

#murl= "http://may2009.archive.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=ENSP00000348965"
findGeneSymbol_Hsapiens = function(protein) {
   murl = paste("http://may2009.archive.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;p=", 
                protein, sep="")
   preword = "<head>"
   page = getUrlPage(murl)
   page = trim(page)
   sel  = page!=""
   page = page[sel]

   sel2 = page==preword
   gs =""
   if(sum(sel2)>0) {
     idx  = c(1:length(page))[sel2];
     title= (strsplit(page[idx+1], " "))[[1]]
     sel3 = title == "Transcript:"
     if(sum(sel3)>0) { 
        idx2 = c(1:length(title))[sel3]
        gs   = title[idx2+1]
     }
   }

   return(gs)
}



parse_general <- function(x, patt) {
    m <- regexec(patt, x)
    #print(m)
    parts <- do.call(rbind,lapply(regmatches(x, m), `[`) )
    parts
}


# check whether a gene is close to or include any SNPs
# geneCoord is a vect of start and end 
# SNPx is a matrix of 2 columns
#
snpINgene = function(geneCoord, SNPx, distcut=1000000){
   sel1 = (SNPx >= geneCoord[1]) & (SNPx <= geneCoord[2])
   sel2 = (SNPx <= geneCoord[1]) & (SNPx >= geneCoord[2])

   diffS= abs(geneCoord[1] - SNPx)
   diffE= abs(geneCoord[2] - SNPx)
   sel3 = (diffS<=distcut) | (diffS<=distcut)
   sel = sel1 | sel2 | sel3

   if(sum(sel)>0) {
     #print(names(SNPx)[sel])
     return (1) 
   } else {
     return (0) 
   }
}

# heatmaps from co-expp
#
.heatmap2 <- function(m, colors, ...) {
        heatmap(
                m,
                Rowv=NA, Colv=NA, scale="none", revC=TRUE, symm=TRUE, labRow="", labCol="",             
                ColSideColors=as.character(colors),
                #RowSideColors=as.character(rev(colors)),
                RowSideColors=as.character(colors),
                ...
        )
}

plotTOMHeatmap2<- function(coexppClusters, geneModules, samplingThreshold=1000, enhance=3, ...) {
        if (length(coexppClusters) > samplingThreshold) {
                warning("Too many observations to construct heatmap directly. Sampling to tractable number of observations.")           
                samples <- clusters(coexppClusters)$order[sort(sample(length(coexppClusters), samplingThreshold))]
                .heatmap2(
                        (sampleMatrix(coexppClusters, samples, kind="tom"))^enhance,
                        geneModules[samples],
                        ...
                )
        } else {
                samples <- clusters(coexppClusters)$order
                .heatmap2(
                        (tom(coexppClusters)[samples, samples])^enhance, 
                        geneModules[samples],
                        ...
                )
        }
}


#' @rdname plotTOMHeatmap
plotAdjHeatmap2<- function(coexppClusters, geneModules, samplingThreshold=1000, enhance=3, ...) {
        if (length(coexppClusters) > samplingThreshold) {
                warning("Too many observations to construct heatmap directly. Sampling to tractable number of observations.")           
                samples <- clusters(coexppClusters)$order[sort(sample(length(coexppClusters), samplingThreshold))]
                .heatmap2(
                        (1-sampleMatrix(coexppClusters, samples, kind="adj"))^enhance,
                        geneModules[samples],
                        ...
                )
        } else {
                samples <- clusters(coexppClusters)$order
                .heatmap2(
                        (1-adj(coexppClusters)[samples, samples])^enhance, 
                        geneModules[samples],
                        ...
                )
        }
}



rcompute.robust.residuals4<-function(m,x1,x2,x3, x4, add_mean=TRUE){
        #print(m[1:5])
	lm.obj<-lm(m~x1+x2+x3+x4,
                      na.action="na.exclude",singular.ok = TRUE)
	resid.data<-residuals(lm.obj) 
        if (add_mean) {resid.data=resid.data -mean(resid.data, na.rm=TRUE)  + mean(m, na.rm=TRUE)}
        return (resid.data)
}
rcompute.robust.residuals3<-function(m,x1,x2,x3, add_mean=TRUE){
        #print(m[1:5])
	lm.obj<-lm(m~x1+x2+x3,
                      na.action="na.exclude",singular.ok = TRUE)
	resid.data<-residuals(lm.obj) 
        if (add_mean) {resid.data=resid.data -mean(resid.data, na.rm=TRUE)  + mean(m, na.rm=TRUE)}
        return (resid.data)
}


#"pearson" (default), "kendall", or "spearman",
ComputeCorPvalFDR = function(dat, cormethod="pearson", no_permutes=100) {

  #*-------------------------------------------------------------------------------------
  #* STEP 1: compute correlation coefficient 
  #*
  #*
  #*
  #*-------------------------------------------------------------------------------------
  corhelp <- cor(dat, use = "pairwise.complete.obs", method=cormethod)
  corhelp <- ifelse(is.na(corhelp), 0, corhelp)

  corpval <- apply(dat, 2, corTest4multivects, amatrix=dat) 
  corpval <- ifelse(is.na(corpval), 1, corpval)
  corpval <- t(corpval)

  #corhelp<- abs(corhelp)
  no.genes <- dim(dat)[2]
  no.arrays<- dim(dat)[1]

  #*-------------------------------------------------------------------------------------
  #* STEP 2: compute FDR 
  #*
  #*
  #*
  #*-------------------------------------------------------------------------------------
  ABScorhelp = abs(corhelp)

  FDR=matrix(0, no.genes, no.genes); myindex = c(1:no.arrays)
  for(x in c(1:no_permutes) ){
      
     if(x%%10 ==0){ print(paste("Permutation  ", x)) }

     rmatrix = apply(dat, 2, permuteVect)

     randCorhelp <- cor(rmatrix, use = "pairwise.complete.obs", method=cormethod) 
     randCorhelp <- ifelse(is.na(randCorhelp), 0, randCorhelp)
     FDR <- FDR + (abs(randCorhelp)>ABScorhelp)
  }
  FDR = FDR/no_permutes

  res= as.list(rep(NA,3))
  res[[1]]= corhelp 
  res[[2]]= corpval 
  res[[3]]= FDR 
  
  return (res)  
}

vectorDivideByVector = function(vect1, vect2){

  return (cbind(vect1/vect2)) 

}


vectorMultiplySumByVector = function(vect1, vect2){

  return (sum(vect1*vect2)) 

}

vectorElementMultiply = function(vect){
   res = 1
   for(v in vect){
    res = res*v
   }
   return(res)
}


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


 #module	modulesize	signature	signature_size	overlap	fold_change	pvalue
 #grey		20161		grey		16860		15750	1.0391129	4.36E-178
 #yellow	232		blue		526		117	21.50126196	1.04E-129
 #turquoise	382		red		192		101	30.88222622	1.22E-129
 #
network_module_conservation = function(intersectPairFiles, no.networks, minModuleSize=0) {

 # gene annotation based on HSG
 #
 xmtrx = matrix(0, ncol=no.networks+1, nrow=no.networks)
 xmtrx2= matrix(0, ncol=no.networks+1, nrow=no.networks)

 conservedModules = as.list(rep(NA, no.networks))

 nfiles = length(intersectPairFiles)
 for(i in c(1:nfiles)) {

 fmodo = intersectPairFiles[i]
 print(paste(i, fmodo))

 #module	modulesize	signature	signature_size	overlap	fold_change	pvalue
 #grey		20161		grey		16860		15750	1.0391129	4.36E-178
 #yellow	232		blue		526		117	21.50126196	1.04E-129
 #turquoise	382		red		192		101	30.88222622	1.22E-129
 #
 modMtrx  = read.delim(fmodo, header=TRUE)
 modMtrx  = as.matrix(modMtrx)

 msizeA = as.integer(modMtrx[,2])
 msizeB = as.integer(modMtrx[,4])

 sel = (modMtrx[,1] != "grey") & (modMtrx[,3] != "grey") & (msizeA >=minModuleSize) & (msizeB >=minModuleSize)

 modMtrx  = modMtrx[sel, ]
 dim(modMtrx)

 #     modulesize signature_size overlap    pvalue
 #	232            526     	 117 	    1.04e-129
 #	382            192       101        1.22e-129
 #
 keycols = c("modulesize", "signature_size","overlap","pvalue")
 idxs    = getMatchedIndex(colnames(modMtrx), keycols)
 dat     = modMtrx[, idxs]
 dat[1:3,1:3]
 dat     = matrix(as.numeric(dat), nrow=nrow(dat))
 dat[1:3,1:3]
 colnames(dat) = keycols 

 xmods   = unique(modMtrx[,1]); no.x = length(xmods)
 ymods   = unique(modMtrx[,3]); no.y = length(ymods)

 xpercent= dat[,3]/dat[,1]
 ypercent= dat[,3]/dat[,2]

 # cutoffs
 pcut    = 0.05/(no.x*no.y); percut  = 0.20

 # x
 sel  = (dat[,4] < pcut) & (xpercent>=percut)
 cons.xmods = modMtrx[sel, 1]
 n.cons.xmods= length(unique(cons.xmods))
 per.cons.xmods = n.cons.xmods/no.x
 xres = paste(n.cons.xmods, " (", 100*signif(per.cons.xmods,3), ")", sep="")

 # y
 sel  = (dat[,4] < pcut) & (ypercent>=percut)
 cons.ymods = modMtrx[sel, 3]
 n.cons.ymods= length(unique(cons.ymods))
 per.cons.ymods = n.cons.ymods/no.y
 yres = paste(n.cons.ymods, " (", 100*signif(per.cons.ymods,3), ")", sep="")
 
 xi = kindexs[compLables[i, 1]]
 yi = kindexs[compLables[i, 2]]

 xmtrx[xi,xi] = no.x; xmtrx[yi,yi] = no.y
 xmtrx[xi,yi] = xres; xmtrx[yi,xi] = yres;

 xmtrx2[xi,xi] = no.x; xmtrx2[yi,yi] = no.y
 xmtrx2[xi,yi] = concatenate(unique(cons.xmods), ", "); xmtrx2[yi,xi] = concatenate(unique(cons.ymods), ", ");
 
 if(is.na(conservedModules[[xi]])) {conservedModules[[xi]]=unique(cons.xmods)
 } else { conservedModules[[xi]]=intersect(conservedModules[[xi]], cons.xmods)}

 if(is.na(conservedModules[[yi]])) {conservedModules[[yi]]=unique(cons.ymods)
 } else { conservedModules[[yi]]=intersect(conservedModules[[yi]], cons.ymods)}

 }

 # conservation for one network across all the other networks
 for(k in c(1:no.networks) ) {
    if(length(conservedModules[[k]])==0){next}
    if(is.na(conservedModules[[k]])){next}
    no.conserved = length(conservedModules[[k]])
    per.conserved= no.conserved/as.integer(xmtrx[k,k])
    res = paste(no.conserved , " (", 100*signif(per.conserved,3), ")", sep="")
    xmtrx[k,no.networks+1] = res
    xmtrx2[k,no.networks+1]= concatenate(conservedModules[[k]], ", ")
 }

 res = as.list(rep(NA,2))
 res[[1]] = xmtrx
 res[[2]] = xmtrx2

 return (res)
}

#R script to calculate the probability of overlap shared by 3 sets
#Author: Minghui Wang, Bin Zhang
#
#Background:
# From an urn of N balls, sample without replacement three independent sets A, B, and C, which share M balls.
# It is assumed that the balls are unbiased, i.e., each ball is equally likely to be sampled into each set.
#
#Mathematics reference:
#http://stats.stackexchange.com/questions/52004/overlapping-sets-3-significance-test
# Overlapping sets (3) - significance test
#
fisherTest_3sets_vect = function(mvect)
{
        A=mvect[1]; B=mvect[2]; C=mvect[3]; M=mvect[4]; N=mvect[5];
	if(A>N | B>N | C>N | A<M | B<M | C<M) stop('Invalid input\n')
	sum(sapply(M:min(A,B,C),function(m) dmcom(A, B, C, m, N)))
}

fisherTest_3sets = function(A, B, C, M, N)
{
	if(A>N | B>N | C>N | A<M | B<M | C<M) stop('Invalid input\n')
	sum(sapply(M:min(A,B,C),function(m) dmcom(A, B, C, m, N)))
}

dmcom=function(A, B, C, m, N){
	log_den=log_comb(N,B)+log_comb(N,C)
	log_num=sapply(m:A,function(l) {
		log_comb(A,l)+log_comb(N-A, B-l)+log_comb(l,m)+log_comb(N-l,C-m)
	})
	sum(exp(log_num-log_den),na.rm=T)
}
log_comb = function(n, x) {
	if(x==0 | n==x) return(0)
	if(x<0 | x>n) return(NA)
	sum(log((x+1):n)) - sum(log(1:(n-x)))
}




perm = function(n, x) {
  return(factorial(n) / factorial(n-x))
}

comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}

#http://stats.stackexchange.com/questions/52004/overlapping-sets-3-significance-test
# Overlapping sets (3) - significance test
#
# a, b, and c share l elements. there are n elements in total
fisherTest_3sets_org = function (a,b,c, l, n)
{
  sums = 0
  for(i in c(1:l)) {
     isum = comb(a,l)*comb(n-a, b-l)*comb(l,m)*comb(n-l,c-m)   
     sums = isum + sums
  }
  pr=sums/comb(n,b)/comb(n,c)
  return (pr)
}

ScoreMatricByRankOrder = function(nummtrx)
{
  rankMtrx  = apply(nummtrx, 2, rank)
  rankMtrx2 = nrow(rankMtrx) + 1-rankMtrx
  nr        = nrow(rankMtrx2)
  scores    = rep(1, nr)
  for(i in c(1:ncol(rankMtrx2)) ){
     scores =   scores * rankMtrx2[,i]/nr
  }
  return (scores)
}

# Concatenate results for all signatures defined by Column 2

# Comparison     AA    module_size..DE. Overlap Fold_Enrichment    FET_P FET_P_corrected
# purple         AA-UP 328              35      6.550477 8.18e-19     4.24542e-16   
# gray47         AA-UP 328              8       30.372822 3.92e-11     2.03448e-08 
#
OlvpPairsToMatrix = function(pairMtrx) {

  signats   = unique(pairMtrx[,2]); 
  deAll     = NULL

  modules = unique(pairMtrx[,1])
  nmods   = length(modules)
  for(es in signats) {
     sel  = pairMtrx[,2]==es
     emtrx= pairMtrx[sel,-c(2)]
     colnames(emtrx)=paste(es, c("Module", "Size", "Overlap", "Fold_Enrichment", "FET_P", "Corrected_P"), sep=" ")
     if(is.null(deAll)) {
       deAll = emtrx 
     } else {
       deAll = merge(deAll, emtrx, by.x=1, by.y=1, all=TRUE) 
     }
  }

  deAll
}

#module	trait	correlation	Cor.p	FDR
#aliceblue	WHRatio	-0.403	0.197	0.002795367
#aliceblue	BMI	0.0887	0.189	0.6254645
#
corrlPairsToMatrix = function(pairMtrx) {

  signats   = unique(pairMtrx[,2]); 
  deAll     = NULL

  modules = unique(pairMtrx[,1])
  nmods   = length(modules)
  for(es in signats) {
     sel  = pairMtrx[,2]==es
     emtrx= pairMtrx[sel,-c(2)]
     colnames(emtrx)=paste(es, colnames(pairMtrx)[-2], sep=" ")
     if(is.null(deAll)) {
       deAll = emtrx 
     } else {
       deAll = merge(deAll, emtrx, by.x=1, by.y=1, all=TRUE) 
     }
  }

  deAll
}

# perform prediction of exponential relationship using logorithm transform and linear regression model
#
#
exponent_lm_predict = function(x, y, newx){
  xx = log10(x); yy = log10(y)
  fit <- lm(yy ~ xx)
  newy = predict(fit, data.frame(xx=(log10(newx))) )
  return (10^(newy))
}

# calculate coeficients, FDR and p values of correlations between two matrix, 

#method = c("pearson", "kendall", "spearman"
correlation_pvalue_fdr_cuts=function(mDataA, mDataB, selfcorrelation=TRUE, geneinfoA=NULL, geneinfoB=NULL, 
                                method="spearman", no.permutes=100,pvlcuts=0.05, fdrcuts=c(0.05,0.1), fname=NULL) {

  if(!is.null(geneinfoA)){
     infoA = geneinfoA
     if(is.null(colnames(infoA))){colnames(infoA)=paste("Coln", c(1:ncol(infoA)), sep="")}
  } else {
     infoA = cbind(rownames(mDataA));
     colnames(infoA) =c("Gene")
  }
  
  if(!is.null(geneinfoB)){
     infoB = geneinfoB
     if(is.null(colnames(infoB))){colnames(infoB)=paste("Coln", c(1:ncol(infoB)), sep="")}
  } else {
     infoB = cbind(rownames(mDataB))
     colnames(infoB) =c("Trait")
  }

  if(method=="spearman") {
    mDataRankA = apply(mDataA, 1, rank)
    mDataRankA = t(mDataRankA)

    mDataRankB = apply(mDataB, 1, rank)
    mDataRankB = mDataRankB
    colnames(mDataRankA) = colnames(mDataA) 
    rownames(mDataRankB) = colnames(mDataB) 
  } else {
    mDataRankA =   mDataA
    mDataRankB = t(mDataB)
  }

  traitnames = infoB[,1]
  no.traits  = length(traitnames)

  #------------ p value -------------------------
  genesignifMatrixPval=apply(mDataRankA, 1, corpval_multitraits, mDataRankB)
  dim(genesignifMatrixPval)
  genesignifMatrixPval = t(genesignifMatrixPval)

  #------------ correlation -------------------------
  geneCorMatrix =cor(t(mDataRankA), mDataRankB, method=method)
  dim(geneCorMatrix)

  #------------------------------Permutation--------------------------------------
  corMtrxRandoms = as.list(rep(NA, no.permutes))
  for(jj in c(1:no.permutes)) {
    if(jj%%10==0) {print(paste("random ..... ", jj))}
    mDataRankB2 = apply(mDataRankB, 2, permuteVect)
    corMtrxRandoms[[jj]]=cor(t(mDataRankA), mDataRankB2, method=method)
  }

  #------------------------------ FDR --------------------------------------
  fdrs = computeSimpleFDR_Matrix(realMatrix=geneCorMatrix, perMatrixList=corMtrxRandoms, is_pvalue=FALSE, use_absolute_value=TRUE)

  geneCorMatrix        = ifelse(is.na(geneCorMatrix),        0, geneCorMatrix)
  genesignifMatrixPval = ifelse(is.na(genesignifMatrixPval), 1, genesignifMatrixPval)
  fdrs                 = ifelse(is.na(fdrs),                 1, fdrs)

  # self correlation, we need remove the diagnal
  if(selfcorrelation){diag(fdrs)=1}

  #----------------------------- define significant correlations ------------
  #fdrstrs = getSecondPart(as.character(fdrcuts*100), "\\.")
  fdrstrs = as.integer(fdrcuts*1000)

  if(length(pvlcuts) < length(fdrcuts) ) {
    pvlcuts2 = rep(pvlcuts, length(fdrcuts) )
  } else {
    pvlcuts2 = pvlcuts
  }

  cntMtrx = NULL
  for (ic in c(1:length(fdrcuts))) {
    fdrcut= fdrcuts[ic]
    sel   = (fdrs <= fdrcut) & (genesignifMatrixPval <pvlcuts2[ic])
    sum(sel)/ncol(fdrs)
    modMtrx = NULL
    netMtrx = NULL

    for(ii in c(1:no.traits)) {
       isel  = sel[, ii]
       iselp = isel & (geneCorMatrix[,ii]>0)
       iseln = isel & (geneCorMatrix[,ii]<0)

       if(sum(isel)==0) {next}

       istr = paste(traitnames[ii], ".corAll", sep="")
       istrp= paste(traitnames[ii], ".corPos", sep="")
       istrn= paste(traitnames[ii], ".corNeg", sep="")
 
       if(sum(isel)>1) {
         ares = cbind(infoA[isel,  ], geneCorMatrix[isel, ii], genesignifMatrixPval[isel,  ii], fdrs[isel,  ii], rep(istr,  sum(isel)) )
       } else {
         ares = rbind(c(infoA[isel,  ], geneCorMatrix[isel, ii], genesignifMatrixPval[isel,  ii], fdrs[isel,  ii], rep(istr,  sum(isel)) ))
       }

       pres = NULL; nres=NULL
       if(sum(iselp)>1) {
         pres = as.matrix(cbind(infoA[iselp, ], geneCorMatrix[iselp,ii], genesignifMatrixPval[iselp, ii], fdrs[iselp, ii], rep(istrp, sum(iselp)) ))
       } else if (sum(iselp)==1) {
         pres = rbind(c(infoA[iselp, ], geneCorMatrix[iselp,ii], genesignifMatrixPval[iselp, ii], fdrs[iselp, ii], rep(istrp, sum(iselp)) ))
       }

       if(sum(iseln)>1) {
           nres = cbind(infoA[iseln, ], geneCorMatrix[iseln,ii], genesignifMatrixPval[iseln, ii], fdrs[iseln, ii], rep(istrn, sum(iseln)) )
       } else if (sum(iseln)==1) {
           nres = rbind(c(infoA[iseln, ], geneCorMatrix[iseln,ii], genesignifMatrixPval[iseln, ii], fdrs[iseln, ii], rep(istrn, sum(iseln)) ))
       }

       ires    = rbind(ares,pres,nres)

       icounts = c(traitnames[ii], fdrcut, sum(isel), sum(iselp), sum(iseln))
       cntMtrx = rbind(cntMtrx, icounts)
       modMtrx = rbind(modMtrx, ires)
    }
    colnames(modMtrx) = c(colnames(infoA), "r", "p", "FDR", "module")
    colnames(cntMtrx) = c("trait", "FDR", "total", "positive", "negative")
    
    if(!is.null(fname)) {
      fxlsSignif = paste(fname, "_p-fdr", fdrstrs[ic], "-AsModules.xls",  sep="");
      write.table(modMtrx, fxlsSignif, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE) 

      fxlsSignifNet= paste(fname, "_p-fdr", fdrstrs[ic], "-cys.txt",  sep="");
      gnames= modMtrx[, 1]
      tnames= getSecondPart(modMtrx[, ncol(modMtrx)], "\\.", 1)
      signs = getSecondPart(modMtrx[, ncol(modMtrx)], "\\.", 2)
      selp  = signs=="corPos"; seln = signs=="corNeg"; color=ifelse(selp, "red", "green")
      netMtrx = cbind(gnames, tnames, modMtrx[,-c(1:ncol(infoA), ncol(modMtrx))], color)
      write.table(netMtrx, fxlsSignifNet, sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE)       
    }

  }

  if(!is.null(fname)) {
    fxlsGS     = paste(fname, "_cor.xls",  sep="");
    fxlsGSfdr  = paste(fname, "_fdr.xls",  sep="");
    fxlsGSpval = paste(fname, "_pvl.xls",  sep="");
    fcount     = paste(fname, "_log.xls",  sep="");

    write.table(rbind( c(colnames(infoA), infoB[,1]) ), 
              fxlsGSpval, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    #write.table(cbind(geneInfo, genesignifMatrix), fxlsGSpval, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T) 
    writeHugeTable(cbind(infoA, genesignifMatrixPval), fxlsGSpval, colnames=F, myappend=TRUE)

    write.table(rbind( c(colnames(infoA), infoB[,1]) ),    fxlsGS, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(cbind(infoA, geneCorMatrix), fxlsGS, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T) 

    write.table(rbind( c(colnames(infoA), traitnames) ),    fxlsGSfdr, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)
    write.table(cbind(infoA, fdrs), fxlsGSfdr, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T) 

    write.table(cntMtrx, fcount,  sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE) 
  }

  NodeB = getSecondPart(modMtrx[, ncol(modMtrx)], "\\.", 1); 
  pp = getSecondPart(modMtrx[, ncol(modMtrx)], "\\.", 2); 
  psel = pp=="corAll"
  pairMtrx = modMtrx[psel,]
  pairMtrx[, ncol(modMtrx)] = NodeB[psel]
  return (list(pairMtrx, geneCorMatrix, genesignifMatrixPval, fdrs, modMtrx) )
}

#------------------------- RCircos functions ---------------------------------------------------------------------
#
#to line-up multiple variables into RCircos histogram format
module_var_hist=function(Data,var.col){
	#Data is a data.frame with at least 4 columns: Chromosome, chromStart, chromEnd and plot variable(s)
	#on input, each row of Data is a unique chromosome (or module) with chromStart being 1 and chromEnd the chromosome (or module) size
	#var.col specifies the column(s) containing the data to be plotted
	#
	n=length(var.col)
	if(n <= 1) return(Data)
	Data=Data[order(Data$Chromosome),]
	HistData=apply(Data,1,function(dat){
		Len=dat['chromEnd']/n
		Start=round(1+Len*c(0:(n-1)),0)
		End=round(Len*c(1:n),0)
		data.frame(dat['Chromosome'],Start,End,unlist(dat[var.col]),1:n,stringsAsFactors=FALSE,row.names=NULL)
	})
	HistData=do.call(rbind,HistData)
	colnames(HistData)=c('Chromosome','chromStart','chromEnd','V1','PlotColor')
	HistData
}
#to make module link coordinates in ordered RCircos link/Ribbon format
order_module_links=function(Data,module.size){
	#Data is a data.frame with at least 6 columns: chromA, chromStartA, chromEndA, chromB, chromStartB and chromEndB
	Data=Data[order(Data$chromA,Data$chromB),]
	Data=lapply(split(Data,Data$chromA),function(dat){
		n=nrow(dat)
		Len=module.size[dat$chromA[1]]/n
		dat$chromStartA=pmax(round(1+Len*c(n:1)-Len/2,0),1)
		dat$chromEndA=pmax(round(1+Len*c(n:1)-Len/4,0),1)
		dat
	})
	Data=do.call(rbind,Data)
	Data=lapply(split(Data,Data$chromB),function(dat){
		n=nrow(dat)
		Len=module.size[dat$chromB[1]]/n
		dat$chromStartB=pmax(round(1+Len*c(n:1)-Len/2,0),1)
		dat$chromEndB=pmax(round(1+Len*c(n:1)-Len/4,0),1)
		dat
	})
	Data=do.call(rbind,Data)
	Data[order(Data$chromA,Data$chromB),]
}
#convert upper triangular of a square matrix to 3-column data.frame
UpperTriangular2c3=function(Mat){
	NR=nrow(Mat)
	NC=ncol(Mat)
	if(NR!=NC) stop('Input must be a square matrix\n')
	rowid=unlist(sapply(1:NC,function(i) seq(1,i)))
	colid=unlist(sapply(1:NC,function(i) rep(i,i)))
	data=Mat[unlist(sapply(1:NC,function(i) (i-1)*NC+seq(1,i)))]
	data.frame(rowid,colid,data)
}
#create a 2D matrix given row/col ids
c3toMatrix=function(dat,row.col=1,col.col=2,val.col=3,default.value=0){
	#dat is a data.frame or matrix with each row containing a record in the 2D matrix;
	#row.col, col.col and val.col specify the row, column and value for each record of the 2D matrix
	M=matrix(default.value,nrow=length(unique(dat[,row.col])),ncol=length(unique(dat[,col.col])))
	rownames(M)=unique(dat[,row.col])
	colnames(M)=unique(dat[,col.col])
	for(i in 1L:nrow(dat)) M[dat[i,row.col],dat[i,col.col]]=dat[i,val.col]
	M
}

#Modified from RCircos
RCircos.Set.Core.ComponentsW=function (cyto.info, chr.exclude = NULL, tracks.inside, tracks.outside) 
{
    cyto.band.data <- RCircos.Validate.Cyto.InfoW(cyto.info, chr.exclude)
    RCircos.Initialize.Parameters(tracks.inside, tracks.outside)
    RCircos.Set.Cytoband.dataW(cyto.band.data)
    RCircos.Set.Base.Plot.Positions()
    cat("\nRCircos.Core.Components initialized.\n")
}
RCircos.Validate.Cyto.InfoW=function (cyto.info, chr.exclude) 
{
	if (ncol(cyto.info) < 5) {
		stop(paste("Cytoband data must have columns for chromosome,", 
			"chromStart, chromEnd, band, Stain\n", "Current cyto.info columns:", 
			ncol(cyto.info)))
	}
	colnames(cyto.info)[1:5] <- c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain")
	cyto.info$Chromosome <- as.character(cyto.info$Chromosome)
	rows <- grep("chr", cyto.info$Chromosome)
	if (length(rows) < length(cyto.info$Chromosome)) {
		index <- 1:length(cyto.info$Chromosome)
		# index <- index[-rows]
		if(length(rows)>0) index <- index[-rows] #added
		cyto.info$Chromosome[index] <- paste0("chr", cyto.info$Chromosome[index])
	}
	rows <- grep("_", cyto.info$Chromosome)
	if (length(rows) > 0) {
		cyto.info <- cyto.info[-rows, ]
	}
	chromosomes <- unique(cyto.info$Chromosome)
	if (length(chr.exclude) > 0) {
		num.chroms <- grep("chr", chr.exclude)
		if (length(num.chroms) != length(chr.exclude)) {
			stop("chromosomes in chr.exclude must have prefix of chr")
		}
		ignore <- chromosomes %in% chr.exclude
		if (sum(ignore) > 0) {
			chromosomes <- chromosomes[ignore == F]
		}
	}
	chr.names <- chromosomes[grep("[0-9]", chromosomes)]
	chr.names <- chr.names[order(as.numeric(sub("chr", "", chr.names)))]
	if (length(chr.names) < length(chromosomes)) {
		other.chr <- chromosomes[-grep("[0-9]", chromosomes)]
		other.chr <- other.chr[order(other.chr)]
		chr.names <- c(chr.names, other.chr)
	}
	chromosomes <- chr.names
	new.cyto.info <- cyto.info[which(cyto.info$Chromosome == chromosomes[1]), ]
	for (a.chr in 1:length(chromosomes)) {
		the.row <- which(cyto.info$Chromosome == chromosomes[a.chr])
		the.info <- cyto.info[the.row, ]
		the.start <- as.numeric(the.info$ChromStart)
		the.end <- as.numeric(the.info$ChromEnd)
		the.info <- the.info[order(the.start), ]
		if (the.start[1] > 1) {
			stop("Cytoband start position cannot be greater than 1.")
		}
		if(length(the.row)>1){ #added
			for (a.row in 2:nrow(the.info)) {
				if (the.start[a.row] < the.end[(a.row - 1)]) {
					stop("Cytoband start position cannot be less than previous end position.")
				}
			}
		} #added
		if (a.chr == 1) {
			new.cyto.info <- the.info
		} else {
			new.cyto.info <- rbind(new.cyto.info, the.info)
		}
	}
	return(new.cyto.info)
}
RCircos.Set.Cytoband.dataW=function (cyto.band.info)  #changed
{
	# stain2color <- as.character(cyto.band.info$Stain)
	# band.color <- rep(colors()[652], length(stain2color))
	# stains <- c("gneg", "acen", "stalk", "gvar", "gpos", "gpos100", 
		# "gpos75", "gpos66", "gpos50", "gpos33", "gpos25")
	# color.index <- c(1, 552, 615, 418, 24, 24, 193, 203, 213, 223, 233)
	# for (a.stain in 1:length(stains)) {
		# bands <- which(stain2color == stains[a.stain])
		# if (length(bands) > 0) {
			# band.color[bands] <- colors()[color.index[a.stain]]
		# }
	# }
	# cyto.band.info["BandColor"] <- band.color
	cyto.band.info$BandColor <- cyto.band.info$Stain #added
    # chrom.color <- c(552, 574, 645, 498, 450, 81, 26, 584, 524, 
    #     472, 32, 57, 615, 635, 547, 254, 100, 72, 630, 589, 8, 
    #     95, 568, 52)
    # chrom2color <- as.character(cyto.band.info$Chromosome)
    # chromosomes <- unique(chrom2color)
    # num.chrom <- length(chromosomes)
    # num.color <- length(chrom.color)
    # if (num.chrom > num.color) {
        # recycle.time <- floor(num.chrom/num.color)
        # if (recycle.time > 1) {
            # chrom.color <- rep(chrom.color, recycle.time)
        # }
        # remains <- num.chrom%%num.color
        # if (remains > 0) {
            # chrom.color <- c(chrom.color, chrom.color[1:remains])
        # }
    # }
    # for (a.chr in 1:length(chromosomes)) {
        # rows <- which(chrom2color == chromosomes[a.chr])
        # if (length(rows) > 0) {
            # chrom2color[rows] <- colors()[chrom.color[a.chr]]
        # }
    # }
    # cyto.band.info["ChrColor"] <- chrom2color
    plot.par <- RCircos.Get.Plot.Parameters()
    cyto.band.info$ChromStart <- as.numeric(cyto.band.info$ChromStart)
    cyto.band.info$ChromEnd <- as.numeric(cyto.band.info$ChromEnd)
    band.len <- cyto.band.info$ChromEnd - cyto.band.info$ChromStart
    cyto.band.info["Length"] <- band.len
	cyto.band.info["Unit"] <- band.len/plot.par$base.per.unit #added
	cyto.band.info["Location"] <- cumsum(cyto.band.info$Unit) #added
    if (plot.par$chrom.paddings > 0) {
        chroms <- unique(cyto.band.info$Chromosome)
        chroms <- chroms[(chroms == chroms[1]) == F]
        num.pad <- plot.par$chrom.paddings
        for (a.chr in 1:length(chroms)) {
            index <- grep(paste(chroms[a.chr], "$", sep = ""), cyto.band.info$Chromosome)
            cyto.band.info$Location[index] <- num.pad + cyto.band.info$Location[index]
            num.pad <- num.pad + plot.par$chrom.paddings
        }
    }
    RCircosEnvironment <- NULL
    RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
    RCircosEnvironment[["RCircos.Cytoband"]] <- cyto.band.info
}
RCircos.Chromosome.Ideogram.PlotW=function (plot.chromosome.id=TRUE,plot.chromosome.color.bar=TRUE,chr.id.replace=NULL, textcolors=NULL) 
{
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
	RCircos.Pos <- RCircos.Get.Plot.Positions()
	RCircos.Par <- RCircos.Get.Plot.Parameters()
	right.side <- nrow(RCircos.Pos)/2
	outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width
	inner.location <- RCircos.Par$chr.ideog.pos
	chroms <- unique(RCircos.Cyto$Chromosome)

        if(is.null(textcolors)) {
          textcolors2=rep("black", length(chroms))
        } else {
          textcolors2= textcolors
        }

	for (a.chr in 1:length(chroms)) {
		the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr],]
		start <- the.chr$Location[1] - the.chr$Unit[1] + 1
		end <- the.chr$Location[nrow(the.chr)]
		mid <- round((end - start + 1)/2, digits = 0) + start
		chr.color <- the.chr$ChrColor[nrow(the.chr)]
		pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
			RCircos.Pos[end:start, 1] * inner.location)
		pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
			RCircos.Pos[end:start, 2] * inner.location)
		if(plot.chromosome.id){ #added
			if(!is.null(chr.id.replace)) chr.name=chr.id.replace[a.chr]
			else chr.name <- sub(pattern = "chr", replacement = "", chroms[a.chr])
			text(RCircos.Pos[mid, 1] * RCircos.Par$chr.name.pos, 
				RCircos.Pos[mid, 2] * RCircos.Par$chr.name.pos, label = chr.name, 
				srt = RCircos.Pos$degree[mid], col=textcolors2[a.chr])
		} #added

		if(!plot.chromosome.color.bar){ next }

		polygon(pos.x, pos.y)

		chr.bar.inner=RCircos.Par$highlight.pos  #added
		chr.bar.outer=chr.bar.inner+RCircos.Par$chrom.width/2  #added
		pos.x <- c(RCircos.Pos[start:end, 1] * chr.bar.inner,  #added
			RCircos.Pos[end:start, 1] * chr.bar.outer)  #added
		pos.y <- c(RCircos.Pos[start:end, 2] * chr.bar.inner,  #added
			RCircos.Pos[end:start, 2] * chr.bar.outer)  #added
		polygon(pos.x, pos.y,col = chr.color)  #added
	}

	if(plot.chromosome.color.bar) {
        print("Plot bars")

	for (a.band in 1:nrow(RCircos.Cyto)) {
		a.color <- RCircos.Cyto$BandColor[a.band]
		if (a.color == "white") {
			next
		}
		start <- RCircos.Cyto$Location[a.band] - RCircos.Cyto$Unit[a.band] + 1
		end <- RCircos.Cyto$Location[a.band]
		pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
			RCircos.Pos[end:start, 1] * inner.location)
		pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
			RCircos.Pos[end:start, 2] * inner.location)
		polygon(pos.x, pos.y, col = a.color, border = NA)
	}
       }
}

RCircos.Reset.Plot.ParametersW.single=function (new.params)
{
    RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
	RCircosEnvironment[["RCircos.PlotPar"]] <- new.params
}
RCircos.Reset.Plot.ParametersW=function (new.params) #changed
{
    point.chr <- new.params$point.type
    text.color <- new.params$text.color
    heatmap.color <- new.params$heatmap.color
    hist.color <- new.params$hist.color
    line.color <- new.params$line.color
    scatter.color <- new.params$scatter.color
    tile.color <- new.params$tile.color
    bg.color <- new.params$track.background
    grid.color <- new.params$grid.line.color
    params <- unlist(new.params)
    params <- params[-which(params == point.chr)]
    params <- params[-which(params == text.color)]
    params <- params[-which(params == heatmap.color)]
    params <- params[-which(params == hist.color)]
    params <- params[-which(params == line.color)]
    params <- params[-which(params == scatter.color)]
    params <- params[-which(params == tile.color)]
    params <- params[-which(params == bg.color)]
    params <- params[-which(params == grid.color)]
    params <- as.numeric(params)
    if (sum(is.na(params)) > 0) {
        stop("Plot parameters except of point.type must be numeric.")
    }
    if (sum(params < 0) > 0) {
        stop("Plot parameters cannot have negative values")
    }
    old.params <- RCircos.Get.Plot.Parameters()
    cyto.band.data <- RCircos.Get.Plot.Ideogram()
    cyto.band.data$Unit <- NULL
    cyto.band.data$Location <- NULL
    if (new.params$chrom.paddings != 0) {
        padding.const <- 3000/1e+06
        genome.lenth <- sum(cyto.band.data$Length)
        total.units <- genome.lenth/new.params$base.per.unit
        the.padding <- round(padding.const * total.units, digits = 0)
        if (new.params$chrom.paddings > the.padding) {
            cat(paste("\nNote: chrom.padding", new.params$chrom.paddings, 
                "is too big,", "and was reset to", the.padding, 
                "\n"))
            new.params$chrom.paddings <- the.padding
        }
    }
    if (new.params$radius.len < 1) {
        cat("\nNote: radius.len needs be at least 1.0\n\n")
        new.params$radius.len <- 1
    }
    new.params$chr.ideog.pos <- new.params$radius.len + 0.1
    new.params$highlight.pos <- new.params$radius.len + 0.25
    new.params$chr.name.pos <- new.params$radius.len + 0.4
    new.params$highlight.width = round(new.params$radius.len, digits = 0)
    new.params$track.in.start <- new.params$radius.len + 0.05
    new.params$track.out.start <- new.params$radius.len + 0.5
    differ <- old.params$plot.radius - old.params$radius.len
    new.params$plot.radius <- new.params$radius.len + differ
    RCircosEnvironment <- NULL
    RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
    RCircosEnvironment[["RCircos.PlotPar"]] <- NULL
    RCircosEnvironment[["RCircos.Cytoband"]] <- NULL
    RCircosEnvironment[["RCircos.Base.Position"]] <- NULL
    RCircosEnvironment[["RCircos.PlotPar"]] <- new.params
    RCircos.Set.Cytoband.dataW(cyto.band.data) #changed here
    RCircos.Set.Base.Plot.Positions()
}
RCircos.Set.Base.Plot.PositionsW=function ()  # no change
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    total.points <- RCircos.Cyto$Location[nrow(RCircos.Cyto)] + 
        RCircos.Par$chrom.paddings
    interval <- 2 * pi/total.points
    base.val <- seq(0, 2 * pi, interval)
    cor.x <- sin(base.val)
    cor.y <- cos(base.val)
    degree <- rep(0, length(base.val))
    mid <- round((length(base.val) - 1)/2, digits = 0) + 1
    for (pt in 1:mid) {
        degree[pt] <- 90 - (base.val[pt] * 180/pi)
    }
    for (pt in (mid + 1):length(base.val)) {
        degree[pt] <- 270 - (base.val[pt] * 180/pi)
    }
    plot.postions <- data.frame(cor.x, cor.y, degree)
    RCircosEnvironment <- NULL
    RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
    RCircosEnvironment[["RCircos.Base.Position"]] <- plot.postions
}
RCircos.Validate.Genomic.DataW=function (genomic.data, plot.type = c("plot", "link")) #no change
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    plot.type <- tolower(plot.type)
    if (plot.type == "plot") {
        chrom.col <- 1
    }
    else if (plot.type == "link") {
        chrom.col <- c(1, 4)
    }
    else {
        stop("Plot type must be \"plot\" or \"line\"")
    }
    for (a.col in 1:length(chrom.col)) {
        the.col <- chrom.col[a.col]
        genomic.data[, the.col] <- as.character(genomic.data[, the.col])
        for (a.row in 1:nrow(genomic.data)) {
            if (length(grep("chr", genomic.data[a.row, the.col])) == 0) {
                genomic.data[a.row, the.col] <- paste("chr", 
                  genomic.data[a.row, the.col], sep = "")
            }
        }
        cyto.chroms <- unique(as.character(RCircos.Cyto$Chromosome))
        data.chroms <- unique(as.character(genomic.data[, the.col]))
        if (sum(data.chroms %in% cyto.chroms) < length(data.chroms)) {
            cat(paste("Some chromosomes are in genomic data only", 
                "and have been removed.\n\n"))
            all.chroms <- as.character(genomic.data[, the.col])
            genomic.data <- genomic.data[all.chroms %in% cyto.chroms,]
        }
        data.chroms <- unique(as.character(genomic.data[, the.col]))
        if (min(genomic.data[, the.col + 1]) < 0) {
            stop("Error! chromStart position less than 0.")
        }
        if (min(genomic.data[, the.col + 2]) < 0) {
            stop("Error! chromEnd position less than 0.")
        }
        for (a.chr in 1:length(data.chroms)) {
            the.chr <- data.chroms[a.chr]
            in.data <- genomic.data[genomic.data[, the.col] == the.chr, ]
            cyto.data <- RCircos.Cyto[grep(the.chr, RCircos.Cyto$Chromosome), ]
            if (max(in.data[, the.col + 1]) > max(cyto.data[, 3]) | max(in.data[, the.col + 2]) > max(cyto.data[, 3])) {
                cat(paste(the.chr, max(in.data[, 2]), max(in.data[, 3]), "\n"))
                stop("Error! Location is outside of chromosome length.")
            }
        }
        for (a.row in 1:nrow(genomic.data)) {
            if (genomic.data[a.row, the.col + 1] > genomic.data[a.row,the.col + 2]) {
                cat("chromStart greater than chromEnd.\n")
                stop(paste("Row:", a.row, genomic.data[a.row,2], genomic.data[a.row, 3]))
            }
        }
    }
    return(genomic.data)
}
RCircos.Histogram.PlotW=function (hist.data, data.col, track.num, side) 
{
	RCircos.Pos <- RCircos.Get.Plot.Positions()
	RCircos.Par <- RCircos.Get.Plot.Parameters()
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
	hist.data <- RCircos.Get.Plot.DataW(hist.data, "plot")
	hist.locations <- as.numeric(hist.data[, ncol(hist.data)])
	start <- hist.locations - RCircos.Par$hist.width
	end <- hist.locations + RCircos.Par$hist.width
	data.chroms <- as.character(hist.data[, 1])
	chromosomes <- unique(data.chroms)
	cyto.chroms <- as.character(RCircos.Cyto$Chromosome)
	for (a.chr in 1:length(chromosomes)) {
		cyto.rows <- which(cyto.chroms == chromosomes[a.chr])
		locations <- as.numeric(RCircos.Cyto$Location[cyto.rows])
		chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]]
		chr.end <- max(locations)
		data.rows <- which(data.chroms == chromosomes[a.chr])
		start[data.rows[start[data.rows] < chr.start]] <- chr.start
		end[data.rows[end[data.rows] > chr.end]] <- chr.end
	}
	locations <- RCircos.Track.Positions(side, track.num)
	out.pos <- locations[1]
	in.pos <- locations[2]
	hist.colors <- RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color)
	num.subtrack <- RCircos.Par$sub.tracks
	# RCircos.Track.Outline(out.pos, in.pos, RCircos.Par$sub.tracks)
	max.hist.height=max(hist.data[, data.col]) #added
	for (a.point in 1:nrow(hist.data)) {
		# hist.height <- hist.data[a.point, data.col]
		hist.height <- hist.data[a.point, data.col]/max.hist.height #added
		the.start <- start[a.point]
		the.end <- end[a.point]
		height <- in.pos + RCircos.Par$track.height * hist.height
		polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * height, 
			RCircos.Pos[the.end:the.start, 1] * in.pos)
		polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * height, 
			RCircos.Pos[the.end:the.start, 2] * in.pos)
		polygon(polygon.x, polygon.y, col = hist.colors[a.point], border = NA)
	}
}
#plot histogram at exact start and end positions
#multiple tracks of related histograms can be plotted at the same time
RCircos.Histogram2.PlotW=function (hist.data, data.cols, track.nums, sides, max.hist.height=NULL, plot.track.bg=FALSE, gradient.color.func=NULL, fixed.height=FALSE, gradient.color.func.numColor=100, height_delta=NULL) 
{
	nt=length(track.nums)
	if(length(data.cols) != nt | nt != length(sides)) stop('data.cols, track.nums and sides must be of equal length\n')
	RCircos.Pos <- RCircos.Get.Plot.Positions()
	RCircos.Par <- RCircos.Get.Plot.Parameters()
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
	hist.data <- RCircos.Get.Plot.Data.BothEnd(hist.data, "plot")
	start <- hist.data$Location1
	end <- hist.data$Location2
	data.chroms <- as.character(hist.data[, 1])
	chromosomes <- unique(data.chroms)
	cyto.chroms <- as.character(RCircos.Cyto$Chromosome)
	for (a.chr in 1:length(chromosomes)) {
		cyto.rows <- which(cyto.chroms == chromosomes[a.chr])
		locations <- as.numeric(RCircos.Cyto$Location[cyto.rows])
		chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]]
		chr.end <- max(locations)
		data.rows <- which(data.chroms == chromosomes[a.chr])
		start[data.rows[start[data.rows] < chr.start]] <- chr.start
		end[data.rows[end[data.rows] > chr.end]] <- chr.end
	}
	if(is.null(max.hist.height)) max.hist.height=max(hist.data[, data.cols]) #added
	hist.data[, data.cols]=hist.data[, data.cols]/max.hist.height  #convert to a range between 0 and 1
	for(it in 1:nt){
		side=sides[it]
		track.num=track.nums[it]
		data.col=data.cols[it]
		locations <- RCircos.Track.Positions(side, track.num)
		out.pos <- locations[1]
		in.pos <- locations[2]
		if(is.null(gradient.color.func)){
			hist.colors <- RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color)
		}else{
			hist.colors <- gradient.color.func(gradient.color.func.numColor)[pmax(1,ceiling(gradient.color.func.numColor*hist.data[, data.col]))]
		}
		num.subtrack <- RCircos.Par$sub.tracks
		if(plot.track.bg) RCircos.Track.Outline(out.pos, in.pos, RCircos.Par$sub.tracks)
		for (a.point in 1:nrow(hist.data)) {
			hist.height <- ifelse(fixed.height,1,hist.data[a.point, data.col])
			if(hist.height==0) next
			the.start <- start[a.point]+1
			the.end <- end[a.point]
                        if(!is.null(height_delta)) {
       			height <- in.pos + RCircos.Par$track.height * hist.height + RCircos.Par$track.height*height_delta[a.point, it]
                        } else {
                        height <- in.pos + RCircos.Par$track.height * hist.height + RCircos.Par$track.height
                        }

			polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * height, 
				RCircos.Pos[the.end:the.start, 1] * in.pos)
			polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * height, 
				RCircos.Pos[the.end:the.start, 2] * in.pos)
			polygon(polygon.x, polygon.y, col = hist.colors[a.point], border = NA)
		}
	}
}

Rcircos.GetValues_for_MultipleRings=function(mydat, no.rings, start_from_zero=TRUE) {
        if(start_from_zero) {
           sdat = mydat-min(mydat)
        } else {
           sdat = mydat
        }
        srange       = max(sdat)-min(sdat); sintervals = srange/no.rings
        scoreMtrx    = matrix(0, nrow=length(sdat), ncol=no.rings)
        heiDeltaMtrx = matrix(0, nrow=length(sdat), ncol=no.rings)
        ssum         = 0
        for(i in c(1:no.rings)){
           jthreshold = i*sintervals
           if(i==1) {
              scoreMtrx[,i] = ifelse(sdat >= i*sintervals, sintervals+min(sdat), sdat)
           } else {
              scoreMtrx[,i] = ifelse(sdat-ssum >= sintervals, sintervals, sdat-ssum)
              scoreMtrx[,i] = ifelse(scoreMtrx[,i-1] <=0, 0, scoreMtrx[,i])
              scoreMtrx[,i] = ifelse(scoreMtrx[,i] <=0, 0, scoreMtrx[,i])
              heiDeltaMtrx[,i-1] = ifelse( scoreMtrx[,i] >0, 0.2, 0)
           }
           ssum = scoreMtrx[,i] + ssum
        }

      return (list(ringValueMtrix=scoreMtrx, ringGapMatrix=heiDeltaMtrx) )
}


#get both start and end points
RCircos.Get.Plot.Data.BothEnd=function (genomic.data, plot.type)
{
	genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, plot.type)
	genomic.data["Location1"]=unsplit(lapply(split(1:nrow(genomic.data),genomic.data[,1]),function(i){
		RCircos.Data.PointW(genomic.data[i[1], 1], genomic.data[i, 2])
	}),genomic.data[,1])
	genomic.data["Location2"]=unsplit(lapply(split(1:nrow(genomic.data),genomic.data[,1]),function(i){
		RCircos.Data.PointW(genomic.data[i[1], 1], genomic.data[i, 3])
	}),genomic.data[,1])
	genomic.data <- genomic.data[order(genomic.data$Location1), ]
	return(genomic.data)
}
RCircos.Get.Plot.DataW=function (genomic.data, plot.type) # changed
{
	genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, plot.type)
	dat=round((genomic.data[,2]+genomic.data[,3])/2,digits=0) #added
	genomic.data["Location"] <- unsplit(lapply(split(1:nrow(genomic.data),genomic.data[,1]),function(i){ #added
		RCircos.Data.PointW(genomic.data[i[1],1],dat[i]) #added
		}),genomic.data[,1]) #added
	genomic.data <- genomic.data[order(genomic.data$Location), ]
	return(genomic.data)
}
RCircos.Link.PlotW=function (link.data, track.num, by.chromosome = FALSE) 
{
	RCircos.Pos <- RCircos.Get.Plot.Positions()
	RCircos.Par <- RCircos.Get.Plot.Parameters()
	link.data <- RCircos.Validate.Genomic.Data(link.data, plot.type = "link")
	one.track <- RCircos.Par$track.height + RCircos.Par$track.padding
	start <- RCircos.Par$track.in.start - (track.num - 1) * one.track
	base.positions <- RCircos.Pos * start
	data.points=matrix(0,nrow=nrow(link.data),ncol=2)
	data.points[,1]=unsplit(lapply(split(1:nrow(link.data),link.data[,1]),function(i){ #added
		RCircos.Data.PointW(link.data[i[1],1],link.data[i,2]) #added
	}),link.data[,1]) #added
	data.points[,2]=unsplit(lapply(split(1:nrow(link.data),link.data[,4]),function(i){ #added
		RCircos.Data.PointW(link.data[i[1],4],link.data[i,5]) #added
	}),link.data[,4]) #added
	link.colors <- RCircos.Get.Link.Colors(link.data, by.chromosome)
	sw=data.points[, 1]>data.points[, 2]
	data.points[sw,c(1,2)]=data.points[sw,c(2,1)]
	P0 <- as.matrix(base.positions[data.points[, 1], 1:2])
	P2 <- as.matrix(base.positions[data.points[, 2], 1:2])
	bc.point.num = 1000;tx=seq(0, 1, 1/bc.point.num);mtx2=(1 - tx)^2;tx2=tx^2
	for (a.link in 1:nrow(data.points)) {
		links <- RCircos.Link.LineW(P0[a.link,], P2[a.link,], mtx2=mtx2, tx2=tx2)
		lines(links$pos.x, links$pos.y, type = "l", col = link.colors[a.link])
	}
}
RCircos.Link.LineW=function (P0, P2, mtx2, tx2) 
{
	#bc.point.num <- 1000
	#t <- seq(0, 1, 1/bc.point.num)
	#link.x <- (1 - t)^2 * P0[1] + t^2 * P2[1]
	#link.y <- (1 - t)^2 * P0[2] + t^2 * P2[2]
	link.x <- mtx2 * P0[1] + tx2 * P2[1]
	link.y <- mtx2 * P0[2] + tx2 * P2[2]
	return(list(pos.x = link.x, pos.y = link.y))
}
RCircos.Data.PointW=function (chromosome, starts) # changed
{

	the.point <- 0
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
	chrom.rows <- grep(paste("^", chromosome, "$", sep = ""), RCircos.Cyto$Chromosome)
	sapply(starts,function(start){ #added
		the.row <- which(RCircos.Cyto$ChromStart[chrom.rows] <= start & RCircos.Cyto$ChromEnd[chrom.rows] >= start)[1]
		band.length <- RCircos.Cyto$Length[chrom.rows[the.row]]
		band.units <- RCircos.Cyto$Unit[chrom.rows[the.row]]
		band.location <- RCircos.Cyto$Location[chrom.rows[the.row]]
		the.bases <- start - RCircos.Cyto$ChromStart[chrom.rows[the.row]]
		the.units <- the.bases/band.length * band.units
		the.point <- band.location - band.units + the.units
		return(round(the.point, digits = 0))
	}) #added
}
RCircos.Ribbon.PlotW=function (ribbon.data, track.num, by.chromosome = FALSE, twist = FALSE) #no change
{
	RCircos.Pos <- RCircos.Get.Plot.Positions()
	RCircos.Par <- RCircos.Get.Plot.Parameters()
	ribbon.data <- RCircos.Validate.Genomic.Data(ribbon.data, plot.type = "link")
	one.track <- RCircos.Par$track.height + RCircos.Par$track.padding
	track.out <- RCircos.Par$track.in.start - (track.num - 1) * one.track
	base.positions <- RCircos.Pos * track.out
	data.points <- matrix(rep(0, nrow(ribbon.data) * 4), ncol = 4)
	for (a.link in 1:nrow(ribbon.data)) {
		data.points[a.link, 1] <- RCircos.Data.Point(ribbon.data[a.link,1], ribbon.data[a.link, 2])
		data.points[a.link, 2] <- RCircos.Data.Point(ribbon.data[a.link,1], ribbon.data[a.link, 3])
		data.points[a.link, 3] <- RCircos.Data.Point(ribbon.data[a.link,4], ribbon.data[a.link, 5])
		data.points[a.link, 4] <- RCircos.Data.Point(ribbon.data[a.link,4], ribbon.data[a.link, 6])
		if (data.points[a.link, 1] == 0 || data.points[a.link,2] == 0 || data.points[a.link, 3] == 0 || data.points[a.link,4] == 0) {
			stop("Error in chromosome locations ...")
		}
	}
	ribbon.colors <- RCircos.Get.Link.ColorsW(ribbon.data, by.chromosome)
	for (a.ribbon in 1:nrow(ribbon.data)) {
		start.one <- data.points[a.ribbon, 1]
		end.one <- data.points[a.ribbon, 2]
		if (twist == FALSE) {
			start.two <- data.points[a.ribbon, 3]
			end.two <- data.points[a.ribbon, 4]
		}
		else {
			start.two <- data.points[a.ribbon, 4]
			end.two <- data.points[a.ribbon, 3]
		}
		P0 <- as.numeric(base.positions[end.one, ])
		P2 <- as.numeric(base.positions[start.two, ])
		line.one <- RCircos.Link.Line(P0, P2)
		P0 <- as.numeric(base.positions[end.two, ])
		P2 <- as.numeric(base.positions[start.one, ])
		line.two <- RCircos.Link.Line(P0, P2)
		polygon.x <- c(base.positions[start.one:end.one, 1], 
			line.one$pos.x, base.positions[start.two:end.two,1], line.two$pos.x)
		polygon.y <- c(base.positions[start.one:end.one, 2], 
			line.one$pos.y, base.positions[start.two:end.two,2], line.two$pos.y)
		polygon(polygon.x, polygon.y, border = NA, col = ribbon.colors[a.ribbon])
	}
}
RCircos.Get.Link.ColorsW=function (link.data, by.chromosome) #no change
{
	red.color <- rgb(1, 0, 0, alpha = 0.5)
	blue.color <- rgb(0, 0, 1, alpha = 0.5)
	link.colors <- rep(blue.color, nrow(link.data))
	if (by.chromosome == TRUE) {
		for (a.row in 1:nrow(link.data)) {
			if (link.data[a.row, 1] == link.data[a.row, 4]) {
				link.colors[a.row] <- red.color
			}
		}
	}
	else {
		color.col <- grep("PlotColor", colnames(link.data))
		if (length(color.col == 1)) {
			the.color <- as.character(link.data[, color.col])
			for (a.row in 1:length(the.color)) {
				rgb.val <- as.vector(col2rgb(the.color[a.row]))/255
				link.colors[a.row] <- rgb(red = rgb.val[1], green = rgb.val[2], 
				blue = rgb.val[3], alpha = 0.5)
			}
		}
		else {
			for (a.row in 1:nrow(link.data)) {
				rgb.val <- as.vector(col2rgb(a.row + 1))/255
				link.colors[a.row] <- rgb(red = rgb.val[1], green = rgb.val[2], 
				blue = rgb.val[3], alpha = 0.5)
			}
		}
	}
	return(link.colors)
}


# correct dates names as regular gene names

correct_GeneDateNames = function(mygenes){
   mygenes1= correct_GeneDateNamesSub(allnodes=mygenes, gformat="xdot")
   mygenes2= correct_GeneDateNamesSub(allnodes=mygenes1,gformat="dash")
   return(mygenes2)
}


#xdot X1.Dec        X1.Mar        X1.Sep
#dash 1-Dec         1-Mar         1-Sep
correct_GeneDateNamesSub = function(allnodes, gformat="xdot", months=c("Dec", "Mar", "Sep")) {

  if(gformat=="xdot") {
    p2       = getSecondPart(fullfnames=allnodes, sep="\\.", 2)
  } else {
    p2       = getSecondPart(fullfnames=allnodes, sep="-", 2)
  }
  ismonth  = setElementInSet_Fast(p2, months)
  sum(ismonth)
  if(sum(ismonth)==0) {return(allnodes) }

  dategenes= allnodes[ismonth]
  dategenes
  restgenes= setdiff(allnodes, dategenes)

  if(gformat=="xdot") {
    parts    = getAllParts(fullfnames=dategenes, sep="\\.")
    p1       = replaceChars(parts[,1], "X", "")
  } else {
    parts    = getAllParts(fullfnames=dategenes, sep="-")
    p1       = parts[,1]
  }

  dategenes2 = paste(parts[,2], p1, sep="")
  newnamesAll= c(dategenes2, restgenes)
  names(newnamesAll) = c(dategenes, restgenes)

  retnames = as.character(newnamesAll[allnodes])

  return(retnames)
}

#*************************************************************************************************************************************************************
######################################### plot complex dendrogram, trait heatmap and exprression heatmap #####################################################
#
averageByGroup = function(vect, group){ 
  res = tapply(X=vect, INDEX=group, FUN=mean, na.rm=TRUE) 
  res2 = cbind(res)	
  rownames(res2) = names(res)
  return (res2)
}


# direction: vertical/horizontal
# color.palette =colorRampPalette(c("green", "black", "red"))
# mycolors = color.palette(nsteps)
#
plotColorIntensityLegend = function(vect=NULL, mseq.in=NULL, colorsvect=NULL,color.palette, no.intervals=50, display.intervals=5, 
                   foutimage=NULL, imgWid=200, imgHei=600, image.resolution=300, direction="vertical", cex.legend=2, margin=3.5) 
{
  # intensity legend plot
  if(is.null(mseq.in)) {  
     exprmin=min(vect); exprmax=max(vect)
     mseq  = seq(from=exprmin, to=exprmax, by= (exprmax-exprmin)/no.intervals)
     mseq  = ifelse(abs(mseq)<0.001, 0, mseq)
     if(exprmax>1) {mseq  =  round(mseq,1)} else {mseq  = signif(mseq,2) }
  } else{
     mseq = mseq.in
  }

  nsteps= length(mseq)

   midx = seq(from=1,to=length(mseq), by=display.intervals); midx = union(midx, length(mseq))
   mseqstr = rep("", length(mseq)); 
   mseqstr[midx] = mseq[midx]
   color.palette =colorRampPalette(c("green", "black", "red"))
   if(is.null(colorsvect)) {
      mycolors = color.palette(nsteps)
   } else {
      mycolors = colorsvect
   }

   intensity = repmat(mseq, 1, nsteps);   dim(intensity)
   colnames(intensity) = rep("", length(mseq)); rownames(intensity) = mseqstr
   #heatmap.plus(intensity, Rowv=NA, Colv=NA, 
   #          margins = c(5,0),cexCol = 1,
   #          scale="none", revC=F, xlab="", col = mycolors)#rgcolors.func(50))

   if(!is.null(foutimage)) {

   openImgDev(foutimage, iwidth =imgWid, iheight = imgHei, ires=image.resolution)

   if(direction=="vertical") {
     par(mfrow=c(1,1), mar=c(0, 0, 0, 4) + 0.1, cex=2)
     heatmap.plus(intensity, Rowv=NA, Colv=NA, 
              cexCol =0.1, cexRow =cex.legend, margins = c(0.3, margin),
             scale="none", revC=F, xlab="", col = mycolors) #heat.colors(50)) #rgcolors.func(50))
   } else {
     par(mfrow=c(1,1), mar=c(3.5, 0, 0, 0) + 0.1, cex=2)
     intensity = t(intensity)
     heatmap.plus(intensity, Rowv=NA, Colv=NA, 
           cexCol =cex.legend, cexRow =0.1, margins = c( margin,0.3),
           scale="none", revC=F, ylab="", col = mycolors) #heat.colors(50)) #rgcolors.func(50))
   }

   par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
   dev.off()

   }

   return (cbind(mseq, mseqstr, mycolors))
}

set_panel_size <- function(p=NULL, g=ggplotGrob(p), width=unit(3, "cm"), height=unit(3, "cm")){
  panel_index_w<- g$layout$l[g$layout$name=="panel"]
  panel_index_h<- g$layout$t[g$layout$name=="panel"]
  g$widths[[panel_index_w]] <- width
  g$heights[[panel_index_h]] <- height
  class(g) <- c("fixed", class(g), "ggplot")
  g
}

set_panel_height <- function(g, height.ratio=0.5){
  panel_index_h<- g$layout$t[g$layout$name=="panel"]
  curhei=g$heights[[panel_index_h]]
  curheiInt = as.numeric(replaceString(as.character(curhei), "null", ""))
  if(!is.na(curheiInt)) {curheiInt =1}
  g$heights[[panel_index_h]] <-  unit(height.ratio*curheiInt, "cm")

  class(g) <- c("fixed", class(g), "ggplot")
  g
}

set_panel_width_height <- function(g, width.ratio=0.5, height.ratio=0.5){
  panel_index_h<- g$layout$t[g$layout$name=="panel"]
  curhei=g$heights[[panel_index_h]]
  curheiInt = as.numeric(replaceString(as.character(curhei), "null", ""))
  if(is.na(curheiInt)) {curheiInt =1}
  g$heights[[panel_index_h]] <-  unit(height.ratio*curheiInt, "cm")


  panel_index_w<- g$layout$l[g$layout$name=="panel"]
  curwid=g$widths[[panel_index_w]]
  curwidInt = as.numeric(replaceString(as.character(curwid), "null", ""))
  if(is.na(curwidInt)) {curwidInt =1}
  g$widths[[panel_index_w]] <-  unit(width.ratio*curwidInt, "cm")

  class(g) <- c("fixed", class(g), "ggplot")
  g
}





unit.pmax2 <- function (...) 
{
    select.i <- function(unit, i) {
        unit[i, top = FALSE]
    }
    x <- list(...)
    numargs <- length(x)
    if (numargs == 0L) 
        stop("no arguments where at least one expected")
    maxlength <- 0L
    for (i in seq_len(numargs)) if (length(x[[i]]) > maxlength) 
        maxlength <- length(x[[i]])        
    ## result <- max(unit.list.from.list(lapply(x, select.i, 1L)))
    UL <- grid:::unit.list.from.list(lapply(x, select.i, 1L))                 ##
    result <- if(all(sapply(UL, attr, "unit")=="null")) {                     ##
                  UL[which.max(UL)]} else {max(UL)}                           ##
    if (maxlength > 1L) 
        for (i in 2L:maxlength) {
            ## result <- unit.c(result, max(unit.list.from.list(lapply(x, 
            ##             select.i, i))))
            UL <- grid:::unit.list.from.list(lapply(x, select.i, i))          ##
            temp <- if(all(sapply(UL, attr, "unit")=="null")) {               ##
                        UL[which.max(UL)]} else {max(UL)}                     ##
            result <- unit.c(result, temp)                                    ##
        }
    result
}

create.grid.layout.rows = function(nrows, width=1, relative.heights=rep(1,nrows) ) {

  w2 <- list(unit(width,"null"), unit(width,"null")); class(w2) <-  c("unit.list", "unit")
  w2[1] <- unit.pmax2(unit(width,"null"), unit(width,"null"))

  h.lists= NULL
  for(i in c(1:nrows)) {
    h <- unit(relative.heights[i], "null"); 
    h.lists= c(h.lists, h)
  }

  ## w2[[1]] <- unit.pmax(unit(1,"null"), unit(1,"null"))  ## For comparison
  g.layout <- grid.layout(nrows, 1, widths =rep(w2,nrows), heights =h.lists,respect = TRUE)

  return (g.layout)
}

create.grid.layout.matrix = function(nrows, ncols, relative.widths=rep(1,ncols), relative.heights=rep(1,nrows) ) {

  w.lists= NULL
  for(i in c(1:ncols)) {
     #w2 <- list(unit(relative.widths[i],"null"), unit(relative.widths[i],"null")); class(w2) <-  c("unit.list", "unit")
     #w2[1] <- unit.pmax2(unit(relative.widths[i],"null"), unit(relative.widths[i],"null"))
     w2= unit(relative.widths[i],"cm")
     w.lists= c(w.lists, w2)
  }

  h.lists= NULL
  for(i in c(1:nrows)) {
    h <- unit(relative.heights[i], "null"); 
    h.lists= c(h.lists, h)
  }

  ## w2[[1]] <- unit.pmax(unit(1,"null"), unit(1,"null"))  ## For comparison
  #g.layout <- grid.layout(nrows, ncols, widths =w.lists, heights =h.lists,respect = TRUE)
  g.layout <- grid.layout(nrows, ncols, widths =w.lists, heights =h.lists,respect = TRUE)

  return (g.layout)
}





create.grid.layout.matrix.org = function(nrows, ncols, relative.widths=rep(1,ncols), relative.heights=rep(1,nrows) ) {

  w.lists= NULL
  for(i in c(1:ncols)) {
     #w2 <- list(unit(relative.widths[i],"null"), unit(relative.widths[i],"null")); class(w2) <-  c("unit.list", "unit")
     #w2[1] <- unit.pmax2(unit(relative.widths[i],"null"), unit(relative.widths[i],"null"))
     w2= unit(relative.widths[i],"null")
     w.lists= c(w.lists, w2)
  }

  h.lists= NULL
  for(i in c(1:nrows)) {
    h <- unit(relative.heights[i], "null"); 
    h.lists= c(h.lists, h)
  }

  ## w2[[1]] <- unit.pmax(unit(1,"null"), unit(1,"null"))  ## For comparison
  g.layout <- grid.layout(nrows, ncols, widths =w.lists, heights =h.lists,respect = TRUE)

  return (g.layout)
}




## plot dendrogram with heatmap and labels
# columns of data.mat are genes and rows are samples
# columns of subtype.matrix are traits and rows are samples
#
plot.Dendrogram.Traits.Expression <- function(subtype.matrix,HC.output,data.mat,plotfname,sort.genes = TRUE, 
                     method = dist.type,linkage = link.type, imgWid=2000, imgHei=3000, imgresolution=300, 
                     traitFontSize=10, geneFontSize=7, hcluFontSize=8, legendFontSize=7)
{

 #traitFontSize=10; geneFontSize=7; hcluFontSize=8; imgWid=2000; imgHei=3000

 # load ordered label
 ordered.label <- HC.output$dendro$labels[HC.output$dendro$order]
 
 ############# plot labels according to dendrogram ticks.
 # get breaks for module ticks
 #ii <- ceiling(sapply(split(1:nrow(subtype.matrix),subtype.matrix[ordered.label,ncol(subtype.matrix)]),median))
 #break.vec <- ordered.label[ii];names(break.vec) <- [ii]

 # make color labels using geom_tile
 subtype.df <- melt(subtype.matrix)
 
 # make dendrogram
 traitColors = unique(as.matrix(subtype.df)[,3]); names(traitColors) = traitColors

 dendro.label = ggplot(data=subtype.df, aes(x=X1, y=X2, fill=value)) + geom_tile() +
 scale_fill_manual(limits = as.character(traitColors), breaks=as.character(traitColors), values =  traitColors) + 
 scale_x_discrete(limits = ordered.label)  +
 scale_y_discrete(limits = rev(colnames(subtype.matrix))) + theme_bw() +
 #theme(axis.ticks = element_blank(),axis.text.x = element_text(angle = -90, hjust = 0,size = 0, color="white"),axis.text.y = element_text(size = traitFontSize),
 theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size = traitFontSize),
 axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none") 

 #png("temp.png",width = 1000);print(dendro.label);dev.off()
 ########### make dendrogram object
 dendro.obj <- ggdendrogram(data = HC.output$dendro,segments = TRUE,labels = FALSE,leaf_labels = FALSE,rotate = FALSE,theme_dendro = TRUE)
 dendro.obj <- dendro.obj + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(angle =0, color="black", size = hcluFontSize))

 ########### make heatmap object
 data.mat2 = scale(data.mat[ordered.label, ], center = T,scale = F)
 emean = mean(data.mat2, na.rm=TRUE); esd = sd(as.numeric(data.mat2), na.rm=TRUE)
 ecutUp= emean + 4*esd; ecutDn= emean - 4*esd
 data.mat2=ifelse(data.mat2 <ecutDn, ecutDn,  data.mat2)
 data.mat2=ifelse(data.mat2 >ecutUp, ecutUp,  data.mat2)

 #------------------ plot legend of expression band ---------------------------------
  exprmin = min(data.mat2); exprmax = max(data.mat2)

  yfname = getFileName(plotfname)
  fintensityV= paste(yfname, "_ExprHeatLegendV.png", sep="")
  fintensityH= paste(yfname, "_ExprHeatLegendH.png", sep="")

  colorshema=plotColorIntensityLegend(vect=data.mat2, mseq.in=NULL, colorsvect=NULL, color.palette=colorRampPalette(c("green", "black", "red")), 
                   no.intervals=50, display.intervals=5, foutimage=NULL, #foutimage=fintensityV, 
                   imgWid=200, imgHei=900, image.resolution=150, direction="vertical", cex.legend=1.5, margin=3.5) 

  #colorshema=plotColorIntensityLegend(vect=data.mat2, mseq.in=NULL, colorsvect=NULL, color.palette=colorRampPalette(c("green", "black", "red")), 
  #                 no.intervals=50, display.intervals=5, foutimage=fintensityH, 
  #                 imgWid=900, imgHei=200, image.resolution=150, direction="horizontal", cex.legend=1.5, margin=3.5) 


  #----------------------------------- Heatmap of expressoin band ----------------------------------------------

  heatmap.df <- melt(data.mat2)
 
  if (sort.genes){
    gdist = get.Distance(t(data.mat2),method =method)
    gene.dendro <- hclust(as.dist(gdist), linkage)
    gene.limits <- gene.dendro$labels[gene.dendro$order];
  }

  dendro.heat <- ggplot(data=heatmap.df, aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient2(low="green", high="red",mid = "black") + ylab("") + 
  scale_x_discrete(limits = ordered.label) + scale_y_discrete(limits = gene.limits) + 
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = geneFontSize),
  axis.text.y = element_text(size = geneFontSize, color="black"),legend.position  = "none")#legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") #legend.position = "none")


  #----------------------------------- Heatmap of expressoin band legend ----------------------------------------------
  ncolors=nrow(colorshema)
  X1= paste("x",c(1:nrow(colorshema)), sep=""); X2=rep("R1", nrow(colorshema)); X3=rep("R2", nrow(colorshema));
  heatmapLeg.df = data.frame(cbind(X1, X2, colorshema[,3])); colnames(heatmapLeg.df) =c("X1", "X2", "value")
  heatmapLegText.df = data.frame(cbind(X1, X2, colorshema[,2])); colnames(heatmapLegText.df) =c("Y1", "Y2", "value2")
 
  #heatmapLeg.df <- melt(leg.mat)
  legColors = unique(as.matrix(heatmapLeg.df)[,3]); names(legColors) = legColors

  legend.heat <- ggplot(data=heatmapLeg.df, aes(x=X1, y=X2, fill=value)) + geom_tile() + ylab("") + 
  scale_fill_manual(limits = as.character(legColors), breaks=as.character(legColors), values =  legColors) + 
  scale_x_discrete(limits = X1)  +
  #scale_y_discrete(limits = unique(heatmapLeg.df[,2])) + theme_bw() +
  theme(axis.ticks = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),
  axis.text.x = element_blank(), axis.text.y=element_blank(),legend.position  = "none")

  #----------------------------------- Heatmap of expressoin band legend TEXT ----------------------------------------------
  heatmapLeg.df2 = data.frame(cbind(X1, X2, rep("white",ncolors))); colnames(heatmapLeg.df2) =c("X1", "X2", "value")
  heatmapLegText.df = data.frame(cbind(X1, X2, colorshema[,2])); colnames(heatmapLegText.df) =c("Y1", "Y2", "value2")
 
  #heatmapLeg.df <- melt(leg.mat)
  legColors2 = unique(as.matrix(heatmapLeg.df2)[,3]); names(legColors2) = legColors2

  legendText.heat <- ggplot(data=heatmapLeg.df2, aes(x=X1, y=X2, fill=value)) + geom_tile() + ylab("") + 
  scale_fill_manual(limits = as.character(legColors2), breaks=as.character(legColors2), values =  legColors2) + 
  scale_x_discrete(limits = X1)  +
  #scale_y_discrete(limits = unique(heatmapLeg.df[,2])) + theme_bw() +
  theme(axis.ticks = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),
  axis.text.x = element_blank(), axis.text.y=element_blank(),legend.position  = "none") + 
  geom_text(data=heatmapLegText.df, aes(x=Y1, y=Y2, label=value2), 
            colour="black", inherit.aes=FALSE, parse=FALSE, size = legendFontSize, angle=-90)


  #----------------------------------- Heatmap of trait band legend ----------------------------------------------
  ncolors=nrow(colorshema)
  X1= paste("x",c(1:nrow(colorshema)), sep=""); X2=rep("R1", nrow(colorshema)); X3=rep("R2", nrow(colorshema));
  heatmapTrait.df = data.frame(cbind(X1, X2, colorshema[,3])); colnames(heatmapTrait.df) =c("X1", "X2", "value")
  heatmapTraitText.df = data.frame(cbind(X1, X2, colorshema[,2])); colnames(heatmapTraitText.df) =c("Y1", "Y2", "value2")
 
  #heatmapLeg.df <- melt(leg.mat)
  traitColorsLeg = unique(as.matrix(heatmapTrait.df)[,3]); names(traitColorsLeg) = traitColorsLeg

  trait.heat <- ggplot(data=heatmapTrait.df, aes(x=X2, y=X1, fill=value)) + geom_tile() + ylab("") + 
  scale_fill_manual(limits = as.character(traitColorsLeg), breaks=as.character(traitColorsLeg), values = traitColorsLeg) + 
  #scale_x_discrete(limits = unique(as.matrix(heatmapTrait.df))[,2])  +
  scale_y_discrete(limits = unique(heatmapTrait.df[,1])) + theme_bw() +
  theme(axis.ticks = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),
  axis.text.x = element_blank(), axis.text.y=element_blank(),legend.position  = "none")


  ####----------------- combine plots --------------------------------------------------------
  # Get the widths
  #
  gA <- ggplot_gtable(ggplot_build(dendro.obj))
  gB <- ggplot_gtable(ggplot_build(dendro.label))
  gC <- ggplot_gtable(ggplot_build(dendro.heat))
  gD <- ggplot_gtable(ggplot_build(legend.heat))
  gE <- ggplot_gtable(ggplot_build(legendText.heat))
  gF <- ggplot_gtable(ggplot_build(trait.heat))

  maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3],gC$widths[2:3], gD$widths[2:3], gE$widths[2:3], gF$widths[2:3])
  #maxHeight= unit.pmax(gA$heights[2:3], gB$heights[2:3],gC$heights[2:3])

  # Set the widths
  gA$widths[2:3] <- maxWidth
  gB$widths[2:3] <- maxWidth
  gC$widths[2:3] <- maxWidth
  gD$widths[2:3] <- maxWidth
  gE$widths[2:3] <- maxWidth
  gF$widths[2:3] <- maxWidth

  resratio = 300/imgresolution

  # Set the heightss
  if(FALSE) {
  gA=set_panel_height(gA, height.ratio=3.8)
  gB=set_panel_height(gB, height.ratio=4.8)
  gC=set_panel_height(gC, height.ratio=16)
  gD=set_panel_height(gD, height.ratio=1)
  g.layout=create.grid.layout.rows(nrows=3, width=2, relative.heights=c(0.4, 0.4,2) )

  } else {

  gA=set_panel_height(gA, height.ratio=3.5*resratio)
  gB=set_panel_height(gB, height.ratio=4.5*resratio)
  gC=set_panel_height(gC, height.ratio=16*resratio)
  gD=set_panel_height(gD, height.ratio=0.8*resratio)
  gE=set_panel_height(gE, height.ratio=0.8*resratio)

  if(FALSE) {
  nprop = 0.5
  gA=set_panel_height(gA, height.ratio=3.5)
  gB=set_panel_height(gB, height.ratio=4.5)
  gC=set_panel_height(gC, height.ratio=16*nprop)
  gD=set_panel_height(gD, height.ratio=0.8*nprop)
  gE=set_panel_height(gE, height.ratio=0.8*nprop)
  }

  g.layout=create.grid.layout.rows(nrows=5, width=2, relative.heights=c(0.4, 0.4,2,0.05,0.05) )

  }

  if(FALSE) {
  widr = 14
  gA=set_panel_width_height(gA, width.ratio=widr, height.ratio=3.5)
  gB=set_panel_width_height(gB, width.ratio=widr, height.ratio=4.5)
  gC=set_panel_width_height(gC, width.ratio=widr, height.ratio=16)
  gD=set_panel_width_height(gD, width.ratio=widr, height.ratio=0.8)
  gE=set_panel_width_height(gE, width.ratio=widr, height.ratio=0.8)
  gF=set_panel_width_height(gF, width.ratio=1,  height.ratio=6)

  g.layout=create.grid.layout.matrix(nrows=5, ncols=2, relative.widths=c(1,0.1), relative.heights=c(0.4, 0.4,2,0.05,0.05) )
  grid.show.layout(g.layout)

  openImgDev(imgname=plotfname, iwidth = imgWid, iheight=imgHei, ipointsize=1, iunits="px", ires=300, icompression="lzw")
  grid.newpage()
  fg <- frameGrob(layout=g.layout)
  fg <- placeGrob(fg, gA, row=1, col=1)
  fg <- placeGrob(fg, gB, row=2, col=1)
  fg <- placeGrob(fg, gC, row=3, col=1)
  fg <- placeGrob(fg, gD, row=4, col=1)
  fg <- placeGrob(fg, gE, row=5, col=1)
  fg <- placeGrob(fg, gF, row=1, col=2)
  fg <- placeGrob(fg, gF, row=2, col=2)
  fg <- placeGrob(fg, gF, row=3, col=2)
  fg <- placeGrob(fg, gF, row=4, col=2)

  grid.draw(fg)
  dev.off()


  }

  grid.show.layout(g.layout)

  openImgDev(imgname=plotfname, iwidth = imgWid, iheight=imgHei, ipointsize=1, iunits="px", ires=imgresolution, icompression="lzw")
  # Arrange the four charts
  #grid.arrange(gA, gB,gC, nrow=3)

  grid.newpage()
  fg <- frameGrob(layout=g.layout)
  fg <- placeGrob(fg, gA, row=1)
  fg <- placeGrob(fg, gB, row=2)
  fg <- placeGrob(fg, gC, row=3)
  fg <- placeGrob(fg, gD, row=4)
  fg <- placeGrob(fg, gE, row=5)

  grid.draw(fg)

  dev.off()

 
 return(0)
}


##################### functions to obtain color maps
#col.names <- c("rosybrown","mediumblue","magenta","yellow","red","green","darkslategray1","goldenrod","hotpink","thistle","purple","orange","blue","sienna1","tan3","turquoise","aliceblue","lightsteelblue","ivory","lightcoral")

### functions to get gradient colors
get.rampColors <- function(low,high,class.names,is.ordered = F)
{
 if (!is.ordered) class.names <- sort(class.names)
 
 color.names <- colorRampPalette(c(low,high))(length(class.names))
 
 output <- data.frame(class.name = class.names,color.name = color.names);
 return(output)
}

### functions to get random colors
get.discreteColors <- function(class.names,col.palette = col.names)
{
 output <- data.frame(class.name = class.names,color.name = col.palette[1:length(class.names)])
 return(output)
}

### function to get continuous colors
get.continuousColors_old <- function(low = "blue",middle="black", high = "red",value.vec,n.cut = 10)
{
 # color ramp
 color.vector <- colorRampPalette(c(low,middle,high))(n.cut+1)
 
 # get ranges
 min.val <- min(value.vec,na.rm = T)
 max.val <- max(value.vec,na.rm = T)
 
 dval <- (max.val - min.val)/n.cut

 ival <- min.val + seq(0,dval * n.cut,dval);
 fval <- ival + dval; fval[length(fval)] = fval[length(fval)] + 0.01*dval;

 ival[1] <- ival[1] - 0.01*dval;

 # output color map
 output = data.frame(start = ival,end = fval,color.name = color.vector)
 
 return(output)
}

# return value = c(start value, end value, color, display value)
#
#   start   end color.name display
#1     -2  0.02    #0000FF       0
#2      0  2.02    #0000CC       2
#3      2  4.02    #000099       4
#
get.continuousColors <- function(low = "blue",middle="black", high = "red",value.vec,n.cut = 10)
{
 # color ramp
 color.vector <- colorRampPalette(c(low,middle,high))(n.cut+1)
 
 # get ranges
 min.val <- min(value.vec,na.rm = T)
 max.val <- max(value.vec,na.rm = T) 
 dval <- (max.val - min.val)/n.cut

 mval <- min.val + seq(0,dval * n.cut,dval);

 ival <- min.val + seq(0,dval * n.cut,dval)-dval;
 fval <- ival + dval; fval[length(fval)] = fval[length(fval)] + 0.01*dval;
 #ival[1] <- ival[1] + 0.01*dval;
 fval <- fval + 0.01*dval;

 # output color map
 output = data.frame(start = ival, end = fval,color.name = color.vector, display=mval)
 
 return(output)
}


find.colors = function(x,tbl, missing.color= "#FFFFFF") {

   if(is.na(x)){return(missing.color) }

   mcol = as.character(tbl$color.name[which(tbl$start <= x & tbl$end > x)])
   return (mcol)
}

### map vectors to colors
map.colors <- function(data.vector,color.map, missing.color= "#FFFFFF")
{
 # data.vector = can be factor, numeric, or character
 # color.map = data.frame object that gives color map. For numeric, it is 3 column data frame where first and second are range values, and third is color. For character or factor, it is two column data frame, where first column are the names, and second column is the color. 
 
 if (!is.data.frame(color.map)) stop("color.map variable is not a data frame.");
 
 output <- NULL
 # handle numeric case
 if (is.numeric(data.vector))
 {
  #cat("Data is numeric\n");
  if (ncol(color.map) != 3) stop("color.map does not have three columns: not an appropriate format to map continuous variable.");

  color.mapped <- sapply(data.vector, find.colors,tbl = color.map, missing.color= missing.color)
  color.mapped <- unlist(color.mapped)

  output <- data.frame(id = names(data.vector),data = data.vector,color.name = color.mapped)
 }
 
 # handle character or factor case. 
 if (is.character(data.vector) | is.factor(data.vector))
 {
  # check if all levels are present in the color.map object
  if (is.character(data.vector))
  {
   level.vec <- unique(data.vector[!is.na(data.vector)])
  }else{level.vec <- levels(data.vector);nm <- names(data.vector);data.vector <- as.character(data.vector);names(data.vector) <- nm;rm(nm)}
  
  if (!all(level.vec %in% color.map[[1]])) 
  {
   cat("warning: not all levels are present in the color map.\n")
   data.vector[which(!(data.vector %in% as.character(color.map[[1]])))] <- NA;
  }
  
  color.vec <- as.character(color.map$color.name);names(color.vec) <- as.character(color.map[[1]]);
  output <- data.frame(id = names(data.vector),data = data.vector,color.name = color.vec[data.vector]);  
 }
 
 return(output)
}

### assign colors automatically

match.colors <- function(data.vector,col.map.defined=NULL, is.discrete = F,low = "blue",high = "red",col.palette = col.names,n.cut = 10)
{
 output <- NULL
 col.map <- NULL
 
 if(is.null(col.map.defined)) {
 ####### get color map
 # handle factor/character vectors
 if (is.character(data.vector) | is.factor(data.vector))
 {
  if (is.factor(data.vector)) {nm <- names(data.vector);data.vector <- as.character(data.vector);names(data.vector) <- nm;}
  
  # get color map
  if (is.discrete)
  {
   col.map <- get.discreteColors(class.names = unique(data.vector[which(!is.na(data.vector))]),col.palette = col.palette)
  }else{
   col.map <- get.rampColors(low = low,high = high,class.names = unique(data.vector[which(!is.na(data.vector))]),is.ordered = F)
  }
 }
 # handle numeric vectors
 if (is.numeric(data.vector))
 {
  col.map <- get.continuousColors(low = low,high = high,value.vec = data.vector,n.cut = n.cut)
 }
 } else {
   col.map = col.map.defined
 }
 
 ####### assign colors to each data value
 if (!is.null(col.map))
 {
  output <- list(color.df = map.colors(data.vector,color.map = col.map),color.map = col.map)
 }
 
 return(output) 
}

################# Misc functions: manipulating many tables.

combine.table <- function(abl,bbl)
{
 common.id <- union(as.character(abl[[1]]),as.character(bbl[[1]]))
 abl.align <- do.call(cbind,lapply(abl[2:ncol(abl)],function(x,y,z) {
                                   names(x) <- z;
                                   vec <- rep(NA,length(y));names(vec) <- y;
								   if (is.numeric(x)) {vec[names(x)] <- x;}else{
								   vec[names(x)] <- as.character(x);}
								   
								   return(vec)
                          },y = common.id,z = as.character(abl[[1]])))
 bbl.align <- do.call(cbind,lapply(bbl[2:ncol(bbl)],function(x,y,z) {
                                   names(x) <- z;
                                   vec <- rep(NA,length(y));names(vec) <- y;
								   if (is.numeric(x)) {vec[names(x)] <- x;}else{
								   vec[names(x)] <- as.character(x);}
								   
								   return(vec)
                          },y = common.id,z = as.character(bbl[[1]])))
 out <- cbind.data.frame(data.frame(id = common.id),abl.align,bbl.align)
 return(out) 
}


coerce.manyTables <- function(table.lst)
{
 if (length(table.lst) > 1)
 {
  out.table <- table.lst[[1]]
  for (i in 2:length(table.lst))
  {
   out.table <- combine.table(out.table,table.lst[[i]])
  }
 }else{
  out.table <- table.lst[[1]];
 }
 
 return(out.table)
}

require(entropy)
KL.symmetric.div <- function(x)
{
 num.bin <- floor(1+log2(nrow(x)))
 binned.x <- apply(x,2,function(x,n) discretize(x,numBins = n),n = num.bin)
 ij <- do.call(rbind,lapply(1:(ncol(binned.x)-1),function(i,n) cbind(rep(i,n-i),(i+1):n),n = ncol(binned.x)))
 
 KL.div <- apply(ij,1,function(xy,m) {
                               x <- m[,xy[1]];y <- m[,xy[2]];
							   if (any(x == 0) | any(y == 0)) {i <- which(x > 0 & y > 0);x <- x[i];y <- y[i];}
							   out <- (KL.empirical(x,y) + KL.empirical(y,x))/2
							   return(out)
							   },m = binned.x)

 D <- matrix(0,ncol(x),ncol(x));rownames(D) <- colnames(x);colnames(D) <- colnames(x);
 for (i in 1:nrow(ij))
 {
  D[ij[i,1],ij[i,2]] <- KL.div[i]
 }
 D <- D + t(D)
 return(D)
}

get.Distance <- function(x,method = c("corDist","corInvert","euclidean","minkowski","KL.sym"),p = 2)
{
 D <- switch(method,corDist = as.dist(sqrt(2*(1-cor(t(x),use = "na.or.complete")))),
                    corInvert = as.dist(1-cor(t(x)/2,use = "na.or.complete")),
					euclidean = dist(x,"euclidean"),
					minkowski = dist(x,"minkowski",p = p),
					KL.sym = as.dist(KL.symmetric.div(t(x))))
 return(D)
}


### this function takes a dist object, and perform HC. 
hierarchy.cluster <- function(x,k,method = "euclidean",linkage = "complete",do.scale = F,dendro.only = F) 
{
 # is data scaled? 
 if (do.scale) x <- t(scale(t(x)))
 # compute distance
 dist.mat <- get.Distance(x,method)
 
 # compute dendrogram
 dendro <- hclust(dist.mat,method = linkage)
 
 # output
 if (dendro.only) # if no clustering results are required
 {
  return(dendro)
 }else{
  # output clustering results.
  vec <- cutree(dendro,k=k)
  output <- list(cluster = vec,dendrogram = dendro,dist.mat = dist.mat);
  return(output)
 }

 
}

### filter out top 90% variance of data matrix
filter.TopVariance <- function(data.mat,total.var = 0.9,balance.dim = T)
{
 ## filter out genes for higher variance. rows are samples and columns are genes in data.mat.
 if (is.null(colnames(data.mat))) colnames(data.mat) <- as.character(1:ncol(data.mat))
 
 colVar <- apply(data.mat,2,function(x) var(x,na.rm = T))
 vec <-sort(colVar,decreasing = T)
 vec <- cumsum(vec)/sum(vec);
 
 id <- names(vec)[which(vec <= total.var)]
 
 min.col <- ceiling(sqrt(nrow(data.mat)));
 
 if (nrow(data.mat) < 20 & balance.dim)
 {
  cat("Number of samples in the data.mat are too few, and not filtering for any genes...\n")
  min.col <- nrow(data.mat)
  balance.dim = FALSE;
  id <- names(vec);
 }
 
 if ((length(id) < min.col) & balance.dim)
 {
  cat("Filtering for variance drastically reduced number of genes under sqrt(nrow(data.mat)). Recovering some genes to ncol = ceiling(sqrt(nrow(data.mat)))...\n")
  id <- names(vec)[1:min(c(ceiling(sqrt(nrow(data.mat))),ncol(data.mat)))];
  cat(paste("Retaining ",signif(vec[min(c(ceiling(sqrt(nrow(data.mat))),ncol(data.mat)))],2)," variance...\n",sep = ""))
 }
 
 output <- data.mat[,id]
 return(output)
}
### perform clustering using Gap statistics to determine optimal k

estimate.k <- function(x,plot.Gap = NULL,method = "euclidean",linkage = "complete",do.scale = F,k.max = 10, no.bootstraps=50)
{
 gap.output <- clusGap(x = x,FUNcluster = hierarchy.cluster, method = method,linkage = linkage,do.scale = do.scale,K.max = min(nrow(x),k.max), B = no.bootstraps, verbose = interactive())
 k.optimal <- which(gap.output$Tab[-1,"gap"] == max(gap.output$Tab[-1,"gap"],na.rm = T)) + 1;
 
 if (!is.null(plot.Gap))
 {
  png(plot.Gap)
  plot(gap.output)
  dev.off()
 }
 
 return(k.optimal)
}

### all-in-one function to perform clustering seamlessly from input matrix, to clusters. 

do.HC <- function(x,method = "euclidean",linkage = "complete",k.max = 10,plot.Gap = NULL,do.scale = F, no.bootstraps=50)
{
 #x = x; plot.Gap = plot.Gap; method = method;linkage = linkage; do.scale = do.scale; k.max = k.max
 optimal.k <- estimate.k(x = x,plot.Gap = plot.Gap,method = method,linkage = linkage,do.scale = do.scale,k.max = k.max,no.bootstraps=no.bootstraps)
 HC.output <- hierarchy.cluster(x = x,k = optimal.k,method = method,linkage = linkage,do.scale = do.scale,dendro.only = F) 
 return(HC.output) 
}

### plot function
# produce subtype tabke

sort.subtypeMatrix <- function(subtype.table,HC.output,color.code = col.names)
{
 ###### convert subtype table to subtype matrix
 subtype.table <- subtype.table[as.character(subtype.table[[1]]) %in% HC.output$dendro$labels,]
 ordered.label <- HC.output$dendro$labels[HC.output$dendro$order]

 subtype.matrix <- as.matrix(subtype.table[,2:ncol(subtype.table)]);rownames(subtype.matrix) <- as.character(subtype.table[[1]]);
 subtype.matrix <- subtype.matrix[rownames(subtype.matrix) %in% ordered.label,]
 leftout <- setdiff(ordered.label,rownames(subtype.matrix))
 leftout.matrix <- matrix(NA,nrow = length(leftout),ncol = ncol(subtype.matrix));rownames(leftout.matrix) <- leftout
 subtype.matrix <- rbind(subtype.matrix,leftout.matrix)
 subtype.matrix <- subtype.matrix[ordered.label,];

 cls <- rep(0,nrow(subtype.matrix));names(cls) <- rownames(subtype.matrix);
 cls[HC.output$dendro$label] <- color.code[HC.output$cluster]
 subtype.matrix <- cbind(subtype.matrix,cls);colnames(subtype.matrix)[ncol(subtype.matrix)] <- "modules"
  
 return(subtype.matrix)
}

is.element.in.matrix = function(element, mmtrx){
   nrows = nrow(mmtrx)
   overlap=rep(FALSE, nrows)
   for(i in c(1:nrows)) {
     overlap[i] = is.element(element, mmtrx[i,])
   }
   overlap
}

return_original=function(vect, topN=10){
  if(length(vect) <topN){ return(vect)
  } else {
    return(vect[1:topN])
  }
}


which.min= function(vect){
 return(which(vect==min(vect, na.rm=TRUE))[1])
}

min.value= function(vect){
 return(min(vect, na.rm=TRUE))
}

multiply = function(vect){
  mres=1
  for(ev in vect){
   mres=mres*ev
  }
  return(mres)
}

greater_than = function(vect1, vect2, absolute=FALSE){
  if(absolute) {
    return (abs(vect1)>abs(vect2))
  } else {
    return (vect1>vect2)
  }
}

greater_than_absolute= function(vect1, vect2){
    return (abs(vect1)>abs(vect2))
}

# if(! any(Ban[,1] %in% rownames(dataMatrix) & Ban[,2] %in% rownames(dataMatrix))) Ban=NULL
#
# TFtargetMatrix=cbind(gene, TF)
#
findPriorSubset = function(TFtargetMatrix, allTFs, geneset){
   TFinSet = intersect(allTFs, geneset)
   if(length(TFinSet)==0){return(NULL)}
   
   mypriors = NULL
   for(eTF in TFinSet){
      etargets = TFtargetMatrix[TFtargetMatrix[,2]==eTF, 1]
      egenes   = intersect(etargets, geneset)
      
      mypriors = rbind(mypriors, cbind(rep(eTF, length(egenes)), egenes) )
   }
   return (mypriors)
}


# mutual information
get.mim <- function(dataset,estimator = "pearson")
{
    mim <- cor(dataset, method = estimator, use = "complete.obs")^2
    diag(mim) <- 0
    maxi <- 0.999999
    mim[which(mim > maxi)] <- maxi
    mim <- -0.5 * log(1 - mim)

    mim[which(mim > 1)] <- 1   

    return(mim)
}


# mutual information
# library(mixtools)
get.mim.from.cor <- function(cormatrix)
{
    mim <- cormatrix^2
    diag(mim) <- 0
    maxi <- 0.999999
    mim[which(mim > maxi)] <- maxi
    mim <- -0.5 * log(1 - mim)

    mim[which(mim > 1)] <- 1

    return(mim)
}

sumSubSet = function(indices, valuesvect){
        return (sum(valuesvect[indices]))
}

splineFitting = function( datvect){
  idat = datvect
  mixmdl = normalmixEM(idat)
  mydensity = density(idat, n=length(idat))
  return(mydensity)
}

pfnParsing = function( pfnCluster, mysep="\t"){
   parts = (getAllParts(pfnCluster, mysep))[1,]
   tmp = cbind(parts[-1], rep(parts[1], length(parts)-1 ))
   return (tmp)
}

varyLengthParsing = function( pfnCluster, mysep="\t"){
   parts = (getAllParts(pfnCluster, mysep))[1,]
   #tmp = cbind(parts[-1], rep(parts[1], length(parts)-1 ))
   return (parts)
}


# LRP10(65),ITPKB(31),DDR1(25),TRIM56(19),CHST6(17),TENC1(16)
#
reformatMegenaModuleKeyDrivers = function(kdstr){

  #myname = names(kdstr)
  #print(paste("Name:", myname))

  if(kdstr=="()"){ return (rbind(c("", ""))) }

  tmp  = getAllParts(kdstr, ",")[1, ]
  tmp2 = getAllParts(tmp, "\\(")
  tmp2[,2] = getSecondPart(tmp2[,2], "\\)", 1)

  #if(nrow(tmp2)==1) {return (tmp2[1,])}
  #mymodule = tmp2)

  return (tmp2)
}


# reformat MEGENA output:the first elelment in mynames is the module name and the rest shows Yes or "", indicating module membership
#
reformatMegenaModules = function(selvector, mynames, selstring = "YES"){
  sel = selvector==selstring
  if(sum(sel)==0){return (NULL)}
  if(sum(sel)==1){return (rbind(c(mynames[sel], mynames[1])))}

  selnames = mynames[sel]
  mymodule = cbind(selnames, rep(selvector[1], length(selnames)))
  return (mymodule)
}

# comp1_11	TIMM10|P62072

reformatMegenaModulesList = function(memvector){
 
  memvector2 = memvector[ memvector!="" ]
  members = getAllParts(memvector2[-1], "\\|")
   
  no.mem = nrow(members)
  if(no.mem==0){return (NULL)}
  if(no.mem==1){return (rbind(c(members, memvector2)))}

  mymodule = cbind(members, rep(memvector2[1], no.mem))
  return (mymodule)
}

fillVectByValue = function(intvalue, nsize, bgval=0, fgval=1){
   myvect = rep(bgval, nsize)
   if(intvalue>0){
      myvect[c(1:intvalue)] = fgval
   }
   return (rbind(myvect))
}

fillInMatrixByValues = function (intvalues, bgval=0, fgval=1){
  res = lapply(intvalues, FUN=fillVectByValue, nsize=length(intvalues), bgval=bgval, fgval=fgval)
  resmtrx = data.frame(do.call(rbind, res))
  return(resmtrx)
}
