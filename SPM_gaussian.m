function good_blocks=SPM_gaussian(x0,x1,pfm,pf,t1,t2,T4,prior)

%x0 and x1 are the observed sample in classes 0 and 1 respectively. Each row is a sample point and each column is a feature (gene), i.e., each observation is a row.
%%pfm is a matrix that the value in row i and column j is the prior of set{i,j}.
%%pf is the prior probability feature f is a good feature. If pf is a number then the same value is used for all features, but if it is a vector it is treated as the vector of prior for each feature.
%If no good prior information is available set pf=0.5 and pfm to all ones matrix
%%if no prior is inputted the code assumes Jeffreys non-informative prior
%%t1, t2, and T4 are thresholds. typically one can use t=1/pf, t2=1 or 2,
%%and T4 large enough to report a reasonable number of blocks
%%if a prior is being used it has be of the following form:
%%prior is a structure with the following fields
%%We assume each set of size 2 has scale matrix S^G_y=s^f_y identity(size G) in class y and S^G=s^f identity(size G) where
%%prior.sf0 is s^f_0, prior.sf1 is s^f_1 and prior.sft is s^f.
%%prior.kf0, prior.kf1, and prior.kft are the degrees of freedom (kappa)
%%prior.mf0, prior.mf0, and prior.mft are the average mean of features, i.e., m^f_0, m^f_1, and m^f, respectively.
%%vf0, vf1 and vft are nu^f_0, ny^f_1, and nu^f in the model

%%Note when using GSG, we assume that for a set of size |G|,
%%in class 0 kappa is |G|+2+kf0e
%%in class 1 kappa is |G|+2+kf1e

%%good blocks is a cell containing the indices of each good block

if nargin<8
    [beta_f, sf0s, sf1s]=POFAC_sub(x0,x1,pfm,pf);
    
    
    good_blocks=SPM_sub(sf0s,sf1s,x0,x1,beta_f,t1,t2,T4);
    
else
    [beta_f, sf0s, sf1s]=POFAC_sub(x0,x1,pfm,pf,prior);
    
    
    good_blocks=SPM_sub(sf0s,sf1s,x0,x1,beta_f,t1,t2,T4,prior);
    
end
