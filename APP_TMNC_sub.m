function app_post_mat=APP_TMNC_sub(x0,x1,pfm,prior)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%initilize
TF=size(x0,2);


n0=size(x0,1); %%%sample size in class 0
n1=size(x1,1); %%%sample size in class 1
nt=n0+n1; %%%%total sample size

xt=[x0;x1];




if nargin<4   %%%%if no prior structure is given, useJeffreys non-informative prior
    
    sf0=0;
    sf1=0;
    sft=0;
    
    kf0=0;
    kf1=0;
    kft=0;
    
    mf0=0;
    mf1=0;
    mft=0;
    
    vf0=0;
    vf1=0;
    vft=0;
    
    lq0=-log(n0)+ logGammaMD(0.5*n0,2);
    lq1=-log(n1)+ logGammaMD(0.5*n1,2);
    lqt=-log(nt) + logGammaMD(0.5*nt,2);
    
    lq=lq0+lq1-lqt;
    
else   %%%%if we do have a prior, use it!
    
    sf0=prior.sf0;
    sf1=prior.sf1;
    sft=prior.sft;
    
    kf0=prior.kf0;
    kf1=prior.kf1;
    kft=prior.kft;
    
    mf0=prior.mf0;
    mf1=prior.mf1;
    mft=prior.mft;
    
    vf0=prior.vf0;
    vf1=prior.vf1;
    vft=prior.vft;
    
    lq0=kf0*log(sf0)+log(vf0)- logGammaMD(0.5*kf0,2)-log(vf0+n0)+ logGammaMD(0.5*(kf0+n0),2);
    lq1=kf1*log(sf1)+log(vf1)- logGammaMD(0.5*kf1,2)-log(vf1+n1)+ logGammaMD(0.5*(kf1+n1),2);
    lqt=kf0*log(sft)+log(vft)- logGammaMD(0.5*kft,2)-log(vft+nt)+ logGammaMD(0.5*(kft+nt),2);
    
    
    
    lq=lq0+lq1-lqt;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%find moments

%%%find means
muf0=mean(x0)';
muf1=mean(x1)';
muft=mean(xt)';

%%find variances
cf0=cov(x0);
cf1=cov(x1);
cft=cov(xt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%find posterior for each feature


%%%the coefficients used for computing updated s value
mfq0=(n0*vf0)./(vf0+n0);
mfq1=(n1*vf1)./(vf1+n1);
mfqt=(nt*vft)./(vft+nt);


%%%updated kappa
kf0s=kf0+n0;
kf1s=kf1+n1;
kfts=kft+nt;

%%%%updated s
sf0s=sf0*eye(TF)+(n0-1)*cf0+mfq0*(muf0-mf0)*(muf0-mf0)';
sf1s=sf1*eye(TF)+(n1-1)*cf1+mfq1*(muf1-mf1)*(muf1-mf1)';
sfts=sft*eye(TF)+(nt-1)*cft+mfqt*(muft-mft)*(muft-mft)';


% % %%%normalize using sample size since it helps with numerics
%
% sf0s=sf0s/(n0-1);
% sf1s=sf1s/(n1-1);
% sfts=sfts/(nt-1);

%%%%%%%%%%%%%%%%%
%%%find determinants

sf0si=repmat(diag(sf0s),1,TF);
sf1si=repmat(diag(sf1s),1,TF);
sftsi=repmat(diag(sfts),1,TF);

dm0=sf0si.*sf0si'-sf0s.^2;
dm1=sf1si.*sf1si'-sf1s.^2;
dmt=sftsi.*sftsi'-sfts.^2;



%%%%%find approximate posterior
log_app_post=lq+log(pfm)-0.5*( kf0s*log(dm0)+kf1s*log(dm1)-kfts*log(dmt) );

for i=1:TF
    log_app_post(i,i)=-inf;
end

%app_post_mat=exp( log_app_post-max(log_app_post(:))+600     );

app_post_mat=exp( log_app_post    );

