function my_det_block=GSGpp(sf0s,sf1s,rfi,u0,TF,n0,n1,t1,t2,prior)

sf0=prior.sf0;
sf1=prior.sf1;
%sft=prior.sft;
vf0=prior.vf0;
vf1=prior.vf1;
%vft=prior.vft;
kfe0=prior.kfe0;
kfe1=prior.kfe1;

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
    
    qu0=0.5*cfs*(cfs+3+kfe0)*log(sf0)- logGammaMD(0.5*(cfs+3+kfe0),cfs+1)+0.5*(cfs+1)*log(vf0/(vf0+n0))+ logGammaMD(0.5*(cfs+3+n0+kfe0),cfs+1) - 0.5*(cfs+1)*(cfs+3+n0+kfe0)*log(n0-1);
    qt0=0.5*cfs*(cfs+2+kfe0)*log(sf0)- logGammaMD(0.5*(cfs+2+kfe0),cfs)+0.5*cfs*log(vf0/(vf0+n0))+ logGammaMD(0.5*(cfs+2+n0+kfe0),cfs) - 0.5*cfs*(cfs+2+n0+kfe0)*log(n0-1);
    qs0=0.5*(3+kfe0)*log(sf0)- logGammaMD(0.5*(3+kfe0),1)+0.5*log(vf0/(vf0+n0))+ logGammaMD(0.5*(3+n0+kfe0),1) - 0.5*(3+n0+kfe0)*log(n0-1);
    
    qu1=0.5*cfs*(cfs+3+kfe1)*log(sf1)- logGammaMD(0.5*(cfs+3+kfe1),cfs+1)+0.5*(cfs+1)*log(vf1/(vf1+n1))+ logGammaMD(0.5*(cfs+3+n1+kfe1),cfs+1) - 0.5*(cfs+1)*(cfs+3+n1+kfe1)*log(n1-1);
    qt1=0.5*cfs*(cfs+2+kfe1)*log(sf1)- logGammaMD(0.5*(cfs+2+kfe1),cfs)+0.5*cfs*log(vf1/(vf1+n1))+ logGammaMD(0.5*(cfs+2+n1+kfe1),cfs) - 0.5*cfs*(cfs+2+n1+kfe1)*log(n1-1);
    qs0=0.5*(3+kfe1)*log(sf1)- logGammaMD(0.5*(3+kfe1),1)+0.5*log(vf1/(vf1+n1))+ logGammaMD(0.5*(3+n1+kfe1),1) - 0.5*(3+n1+kfe1)*log(n0-1);
    
    lq=(qu0-qt0-qs0)+(qu1-qt1-qs1);
    
    sv(found)=-0.5*( (n0+cfs+3+kfe0)*log(v0(found)-t0s) +    (n1+cfs+3+kfe1)*log(v1(found)-t1s)  -  (n0+cfs+2+kfe0)*log(v0(found))  - (n1+cfs+2+kfe1)*log(v1(found)) );
    
    [mval place]=max(sv);
    
    
    if (-log(t1)+lq+mval-t2*cfs*log(n0+n1)>0)
        U=[U place];
        rfi(place)=0;
        cfs=cfs+1;
        
        
        %%%can be used to avoid picking too large blocks
        %if cfs>n0-3
        %    check=0;
        %end
        
        
    else
        check=0;
    end
    
    
    
end

my_det_block=U;
