function mf=REMAIN_sub(log_app_post,T1,T2)


check=1;

check1=0;

mf=[];

posg=exp(log_app_post-max(log_app_post(:)+600));

posg=2*posg./sum(posg(:));

while check
    
    mposg=sum(posg);
    
    found=find(mposg>T1);
    
    if length(found)>0
        
        check1=1;
        
        mf=[mf found];
        
        posg(:,found)=0;
        
        posg=2*posg./sum(posg(:));
        
    else
        
        if check1
            
            log_app_post(mf,:)=-inf;
            log_app_post(:,mf)=-inf;
            
            posg=exp(log_app_post-max(log_app_post(:)+600));
            
            posg=2*posg./sum(posg(:));
            
            check1=0;
            
            [mval dummy]=max(log_app_post);
            if mval<T2
                check=0;
            end
            
        else
            
            check=0;
            
        end
        
    end
    
    
end

