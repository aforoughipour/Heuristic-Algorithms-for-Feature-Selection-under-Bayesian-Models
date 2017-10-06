function my_det_block=GSGjp(sf0s,sf1s,rfi,u0,TF,n0,n1,t1,t2)

check=1;
U=[u0]; %%%the good block we are trying to find
cfs=1;  %%%current feature set size, i.e. |U|

v0=diag(sf0s)';
v1=diag(sf1s)';

while (check>0)
    
    sv=-inf*ones(1,TF);
    
    R0inv=inv(sf0s(U,U));
    R1inv=inv(sf1s(U,U));
    
    found=find(rfi>0);
    
    w0=sf0s(U,found);
    w1=sf1s(U,found);
    
    k0=R0inv*w0;
    k1=R1inv*w1;
    
    if cfs>1
        t0s=sum(w0.*k0);
        t1s=sum(w1.*k1);
        
    else
        
        t0s=(w0.*k0);
        t1s=(w1.*k1);
        
    end
    
    qu0=-0.5*(cfs+1)*log(n0)+ logGammaMD(0.5*n0,cfs+1) - 0.5*(cfs+1)*n0*log(n0-1);
    qt0=-0.5*cfs*log(n0)+ logGammaMD(0.5*n0,cfs) - 0.5*cfs*n0*log(n0-1);
    qs0=-0.5*log(n0)+ logGammaMD(0.5*n0,1) - 0.5*n0*log(n0-1);
    
    qu1=-0.5*(cfs+1)*log(n1)+ logGammaMD(0.5*n1,cfs+1) - 0.5*(cfs+1)*n1*log(n1-1);
    qt1=-0.5*cfs*log(n1)+ logGammaMD(0.5*n1,cfs) - 0.5*cfs*n1*log(n1-1);
    qs1=-0.5*log(n1)+ logGammaMD(0.5*n1,1) - 0.5*n1*log(n1-1);
    
    lq=(qu0-qt0-qs0)+(qu1-qt1-qs1);
    
    sv(found)=-0.5*( n0*log(v0(found)-t0s) +   n1*log(v1(found)-t1s)  -  n0*log(v0(found))  - n1*log(v1(found)) );
    
    [mval place]=max(sv);
    
    
    if (-log(t1)+lq+mval-t2*cfs*log(n0+n1)>0)
        U=[U place];
        rfi(place)=0;
        cfs=cfs+1;
        
        if cfs>min([n0 n1])-3   %%%%to be able to detect a feature of a block we need the estimated covariance matrices to be full rank. This is to guarantee this indeed holds
            check=0;
        end
        
        
    else
        check=0;
    end
    
    
    
end

my_det_block=U;
