function good_blocks=SPM_sub(sf0,sf1,x0,x1,beta_f,t1,t2,T4,prior)

if nargin<9
    JPi=1;
else
    JPi=0;
end

n0=size(x0,1);
n1=size(x1,1);
nt=n0+n1;

TF=size(x0,2);

rfi=ones(1,TF);  %%%%indicator of remaining features, .i.e, features not selected to belong to a good block yet
check=1;


gbp=cell(1,1000);  %%%%%where good blocks will be stored for now

mf=[];

gbcnt=0;  %%%%%%%%%%%good block counter

while (check>0)
    
    r_beta_f=beta_f.*rfi;  %%%%only consider the beta of features which are remining
    
    [max_r_beta_f, u0]=max(r_beta_f); %%%%the maximum beta and the feature corresponding to it being the seed of good seed grower (GSG)
    rfi(u0)=0;
    
    if max_r_beta_f>T4
        
        %%%%%%
        %%find the good block
        if JPi>0
            my_det_block=GSGjp(sf0,sf1,rfi,u0,TF,n0,n1,t1,t2);  %%%%Good seed grower for jeffrey prior
        else
            my_det_block=GSGpp(sf0,sf1,rfi,u0,TF,n0,n1,t1,t2);  %%%%Good seed grower for proper inputted prior
        end
        
        gbcnt=gbcnt+1;
        gbp{gbcnt}=my_det_block;
        
        mf=[mf my_det_block];
        
        rfi(mf)=0;
        
        if sum(rfi)<1
            check=0;
        end
        
    else
        check=0;
    end
end

good_blocks=gbp(1,1:gbcnt);