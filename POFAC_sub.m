function [beta_f, sf0s, sf1s]=POFAC_sub(x0,x1,pfm,pf,prior)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%initilize
TF=size(x0,2);

n0=size(x0,1);
n1=size(x1,1);

if nargin<5   %%%%if no prior structure is given, useJeffreys non-informative prior
    
    [app_post_mat, sf0sc, sf1sc]=APP_TMNC_sub2(x0,x1,pfm);
    un_post_f=igib_sub(x0,x1,pf);
    
    
    sf0s=sf0sc/(n0-1);
    sf1s=sf1sc/(n1-1);
    
    
else   %%%%if we do have a prior, use it!
    
    [app_post_mat, sf0sc, sf1sc]=APP_TMNC_sub2(x0,x1,pfm,prior);
    priorp=prior;
    priorp.kf0=priorp.kf0-1;
    priorp.kf1=priorp.kf1-1;
    priorp.kft=priorp.kft-1;
    un_post_f=igib_sub(x0,x1,pf,priorp);
    
    sf0s=sf0sc/(n0-1);
    sf1s=sf1sc/(n1-1);
    
end

beta_mat=app_post_mat./repmat(un_post_f(:),1,TF);

beta_f=sum(beta_mat)/(TF-1);

