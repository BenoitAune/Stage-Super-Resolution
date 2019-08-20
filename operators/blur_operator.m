function H=blur_operator(G,tau)

%input: G is the gaussian matrix convolution
%       (tau^2)*1 is the size of the vector in column that we want to blur
%output: H is the blur operator (sparse). Its size is (tau^2)*(tau^2)

sg=size(G,1);  %odd
sg2=floor(size(G,1)/2); 
K=zeros(1,2*tau*sg2+(tau-sg2*2)*2*sg2+sg*sg*((tau-sg2*2)^2));
res=ones(1,2*tau*sg2+(tau-sg2*2)*2*sg2+sg*sg*((tau-sg2*2)^2));
l=zeros(1,2*tau*sg2+(tau-sg2*2)*2*sg2+sg*sg*((tau-sg2*2)^2));
k=0;
for j=1:sg2
    for i=1:tau
       k=k+1;
       K(1,k)=i+(j-1)*tau;
       res(1,k)=1;
       l(1,k)=i+(j-1)*tau;
    end
end   
k=k+1;
for j=sg2+1:tau-sg2  
    for i=1:sg2         
       K(1,k)=i+(j-1)*tau;
       res(1,k)=1;
       l(1,k)=i+(j-1)*tau;
       k=k+1;
    end
    k=k+sg*sg*(tau-sg2*2)+sg2;
end   
k=tau*sg2+sg2+(tau-sg2*2)*sg*sg+1;
for j=sg2+1:tau-sg2  
    for i=tau-sg2+1:tau           
       K(1,k)=i+(j-1)*tau;
       res(1,k)=1;
       l(1,k)=i+(j-1)*tau;
       k=k+1;
    end
    k=k+sg*sg*(tau-sg2*2)+sg2;
end
k=tau*sg2+sg2+1;
for j = sg2+1:tau-sg2
    for i = sg2+1:tau-sg2
        for v=1:sg
            K(1,k:k+sg-1)=i+(j-1)*tau; 
            res(1,k:k+sg-1)=G(v,:);
            l(1,k:k+sg-1)=i+(j-1-sg2+v-1)*tau-sg2:i+(j-1-sg2+v-1)*tau+sg2;
            k=k+sg;
        end

    end
    k=k+sg2*2;
end
k=tau*sg2+(tau-sg2*2)*2*sg2+sg*sg*((tau-sg2*2)^2);
for j=tau-sg2+1:tau
    for i=1:tau
       k=k+1;
       K(1,k)=i+(j-1)*tau;
       res(1,k)=1;
       l(1,k)=i+(j-1)*tau;
    end
end
H=sparse(K,l,res,tau*tau,tau*tau);
end