function sel_feats=REMAIN_gaussian(x0,x1,T1,T2,pfm,prior)


%x0 and x1 re the observed sample in classes 0 and 1 respectively. each row is a sample point and each column is a feature (gene), i.e., each observation is a row.
%%T1 and T2 are the threhsolds of REMAIN, we suggest fixing T2 to sample size. T1 is in (0,1). T1 could be small, but it should be larger than 1/(expected number of bad features). 
%%pfm is a matrix that the value in row i and column j is the prior of set{i,j}. Unless good prior information is availbale use an all ones matrix. 
%%if no prior is inputted the code assumes Jeffreys non-informative prior
%%if a prior is being used it has be of the following form:
%%prior is a structure with the following fields
%%prior.sf0 is the s^f_0, the scale paramter of the inverse wishart distribution i.e. it is s^f_0 in the model.
%%prior.sf1 is the s^f_1 and prior.sft is s^f, the scale parameter of bad features
%%prior.kf0, prior.kf1, and prior.kft are the degrees of freedom (kappa)
%%prior.mf0, prior.mf0, and prior.mft are the average mean of features, i.e., m^f_0, m^f_1, and m^f, respectively.
%%vf0, vf1 and vft are nu^f_0, ny^f_1, and nu^f in the model


%%%%%%
%%start with 2MNC robust! we need the un-normlaized posterior matrix! :)

if nargin<6   %%%%if no prior structure is given, useJeffreys non-informative prior
    
    app_post=APP_TMNC_sub(x0,x1,pfm);
    
else   %%%%if we do have a prior, use it!
    
    app_post=APP_TMNC_sub(x0,x1,pfm,prior);
    
end


%%%%%%%%%%%%%%%%%%%%%%
%%%now that we have the matrix call the main remain function

log_app_post=log(app_post);

sel_feats=REMAIN_sub(log_app_post,T1,T2);
