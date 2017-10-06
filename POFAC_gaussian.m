function [feat_rank, feat_beta]=POFAC_gaussian(x0,x1,pfm,pf,prior)

%x0 and x1 are the observed sample in classes 0 and 1 respectively. each row is a sample point and each column is a feature (gene), i.e., each observation is a row.
%%pfm is a matrix that the value in row i and column j is the prior of set{i,j}. Unless good prior information is availbale use an all ones matrix. 
%%pf is the prior probability vector that feature f is a good feature. If pf is inputted as a number then the same value is used for all features.
%%if no prior is inputted the code assumes Jeffreys noninformative prior
%%if a prior is being used it has be of the following form:
%%prior is a structure with the following fields:
%%We assume each set of size 2 has scale matrix S^G_y=s^f_y identity(size G) in class y and S^G=s^f identity(size G) where
%%prior.sf0 is s^f_0, prior.sf1 is s^f_1 and prior.sft is s^f.
%%prior.kf0, prior.kf1, and prior.kft are the degrees of freedom (kappa)
%%prior.mf0, prior.mf0, and prior.mft are the average mean of features, i.e., m^f_0, m^f_1, and m^f, respectively.
%%vf0, vf1 and vft are nu^f_0, ny^f_1, and nu^f in the model

%%note when using the prior for sets of size 1 the code automatically adjusts kappa. So the kappa inserted in the code for sets of size 2 is enough


%%%feat_rank is a ranking of features based on beta, and feat_beta is the
%%%beta of ranked features

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%initilize
TF=size(x0,2);

if nargin<5   %%%%if no prior structure is given, useJeffreys non-informative prior
    
    app_post_mat=APP_TMNC_sub(x0,x1,pfm);
    un_post_f=igib_sub(x0,x1,pf);
    
else   %%%%if we do have a prior, use it!
    
    app_post_mat=APP_TMNC_sub(x0,x1,pfm,prior);
    priorp=prior;
    priorp.kf0=priorp.kf0-1;
    priorp.kf1=priorp.kf1-1;
    priorp.kft=priorp.kft-1;
    un_post_f=igib_sub(x0,x1,pf,priorp);
    
end

beta_mat=app_post_mat./repmat(un_post_f(:),1,TF);

beta_f=sum(beta_mat)/(TF-1);


[feat_beta, feat_rank]=sort(beta_f,'descend');
