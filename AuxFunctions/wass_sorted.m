function res=wass_sorted(v1,v2) 
%%% this function takes 2 vectors of intrinsic distances from the start of
%%% each trajectory and computes a Wasserstein distance between them (assuming discrete uniform measures)
%%% Note: this function can be used to compute 2-Wasserstein distance
%%% between two uniform measures with supports v1 and v2
v1=sort(v1,'descend');
v2=sort(v2,'descend');
if length(v1)==length(v2)
res=sum((v1-v2).^2)/length(v1); %2-Wasserstein
res=sqrt(res);
else
    N=length(v1);
    M=length(v2);
    
    lam=0;
    my_sum=0;
    for i=1:N
        for j=1:M
            if min(i*M,j*N)>max((i-1)*M,(j-1)*N)
                lam=min(i/N,j/M) - max((i-1)/N,(j-1)/M);
                my_sum=my_sum+lam*(v1(i)-v2(j))^2;   
            end             
 res=sqrt(my_sum);   
        end
    end

end