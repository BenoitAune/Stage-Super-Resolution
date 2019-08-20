function [S,St]=sub_operator(q,tau)

%input: (tau^2*q^2)*1 is the size of the vector in column that we want to sample
%       q is the magnification factor
%output: S is the sub-sampling operator (sparse). Its size is (tau^2)*(tau^2*q^2)
%        St is the over-sampling operator (sparse). Its size is (tau^2*q^2)*(tau^2)


K=1:tau*tau*q*q;
res=ones(1,tau*tau*q*q);
m=zeros(1,tau*tau*q*q);
k=0;
for h=1:tau
for l=1:q
for i=1:tau
    for j=1:q
        k=k+1;
        m(1,k)=(h-1)*tau+i;        
    end
end
end
end
St=sparse(K,m,res,tau*tau*q*q,tau*tau);

K=1:tau*tau;
res=ones(1,tau*tau);
m=zeros(1,tau*tau);
k=0;
for j=1:q:tau*q
    for i=1:q:tau*q
        k=k+1;
        m(1,k)=i+tau*q*(j-1);
    end
end
S=sparse(K,m,res,tau*tau,tau*tau*q*q);
end