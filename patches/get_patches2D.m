function [u] = get_patches2D(U,tau,mask)

% input : (tau * tau) is the size of the patches of U
%         mask gives the position of each center of each patch in the image
% output: u contains all the patches of U
%         the size of u is m * n with
%          m=tau*tau the number of pixels in each patch
%          n the number of patches

[M,N]=size(U);
tau2=floor(tau/2);

[x,y]=meshgrid(1:M,1:N);
x=x'.*mask;
y=y'.*mask;
x=x';
y=y';
x=x(x~=0);
y=y(y~=0);
u = zeros(tau*tau,length(x)); 
for k=1:length(x)
        for s=1:tau
            for d=1:tau
                    u(d+tau*(s-1),k)= U(x(k)-tau2+s-1,y(k)-tau2+d-1); 
            end
        end
end
end
