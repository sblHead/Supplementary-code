%%
%---------------Extrapolation of n and k data---- S. Head----------------
% clear all
% clc
%------------------------------------------------------------------------ 
%%

dat = importdata(''); % columns of [wavelength,n,wavelength,k] for each material
lam = (380:1:780);   
N = 2; % number of materials
full_lam = zeros(401,N*2);   % Resolution of 1nm

col1=1;
col2=2;
for n=1:4:N*4
    count = 1;
    count2 =1;
    
    for i=380:1:780
        if dat(count+1,n) > i
            y=((dat(count,n+1)*(dat(count+1,n)-i))+(dat(count+1,n+1)*(i-dat(count,n))))...
                /(dat(count+1,n)-dat(count,n));
            full_lam(i-379,col1)=y;
        else
            count = count+1;
            y=((dat(count,n+1)*(dat(count+1,n)-i))+(dat(count+1,n+1)*(i-dat(count,n))))...
                /(dat(count+1,n)-dat(count,n));
            full_lam(i-379,col1)=y;
        end
        
        if dat(count2+1,n+2) > i
            y2=((dat(count2,n+3)*(dat(count2+1,n+2)-i))+(dat(count2+1,n+3)*(i-dat(count2,n+2))))...
                /(dat(count2+1,n+2)-dat(count2,n+2));
            if isnan(y2)
                y2=0;
            end
            full_lam(i-379,col2)=y2;
            
        else
            count2 = count2+1;
        
            y2=((dat(count2,n+3)*(dat(count2+1,n+2)-i))+(dat(count2+1,n+3)*(i-dat(count2,n+2))))...
                /(dat(count2+1,n+2)-dat(count2,n+2));
            if isnan(y2)
                y2=0;
            end
            full_lam(i-379,col2)=y2;
       
        end
    end  
    col1=col1+2;
    col2=col2+2;
end

condensed = zeros(81,N*2);  

for k =(1:4)
    cond = 1;
    for j=(1:5:401)   % every 5th sample was chosen in this case
       condensed(cond,k)=full_lam(j,k); 
       cond=cond+1;
    end
end

writematrix(condensed,'');  % save output in desired format
