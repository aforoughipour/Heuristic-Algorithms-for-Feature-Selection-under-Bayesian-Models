function feat_rank=TwoMNC_Robust_gaussian(x0,x1,pfm,prior)

%x0 and x1 are the observed sample in classes 0 and 1 respectively. each row is a sample point and each column is a feature (gene), i.e., each observation is a row.
%%pfm is a matrix that the value in row i and column j is the prior of set{i,j}. Unless good prior information is availbale use an all ones matrix. 
%%if no prior is inputted the code assumes Jeffreys noninformative prior
%%if a prior is being used it has be of the following form:
%%prior is a structure with the following fields:
%%We assume each set of size 2 has scale matrix S^G_y=s^f_y identity(size G) in class y and S^G=s^f identity(size G) where
%%prior.sf0 is s^f_0, prior.sf1 is s^f_1 and prior.sft is s^f.
%%prior.kf0, prior.kf1, and prior.kft are the degrees of freedom (kappa)
%%prior.mf0, prior.mf0, and prior.mft are the average mean of features, i.e., m^f_0, m^f_1, and m^f, respectively.
%%vf0, vf1 and vft are nu^f_0, ny^f_1, and nu^f in the model

%%feat_rank is a ranking of features. Given feat_rank the used can pick the
%%top D ones as the good features.

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


%%%normalize using sample size since it helps with numerics

sf0s=sf0s/(n0-1);
sf1s=sf1s/(n1-1);
sfts=sfts/(nt-1);

%%%%%%%%%%%%%%%%%
%%%find determinants

sf0si=repmat(diag(sf0s),1,TF);
sf1si=repmat(diag(sf1s),1,TF);
sftsi=repmat(diag(sfts),1,TF);

dm0=sf0si.*sf0si'-sf0s.^2;
dm1=sf1si.*sf1si'-sf1s.^2;
dmt=sftsi.*sftsi'-sfts.^2;


%%%%%find approximate posterior
log_app_post=log(pfm)-0.5*( kf0s*log(dm0)+kf1s*log(dm1)-kfts*log(dmt) );

for i=1:TF
    log_app_post(i,i)=-inf;
end

un_post=exp( log_app_post-max(log_app_post(:))+600     );

un_post_f=sum(un_post);
[feat_un_post, feat_rank]=sort(un_post_f,'descend');
