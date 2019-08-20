function  u=get_image2D(xi,tau,Mx,Nx,k_hat,mu,sigma,mask,gamma)

%input: xi are the patches 
%       (tau * tau) is the size of the patches
%       mu and sigma are the GMM parameters 
%       k_hat gives the position of the gaussian paramater for each patch
%       mask gives the positions of each patch
%       gamma is a weights coeficients for the image reconstruction
%output:u is a new image with a size of (Mx * Nx)  

[~,Nxi]=size(xi);
wi=zeros(1,Nxi);
tau2=floor(tau/2);
[x,y]=meshgrid(1:Mx,1:Nx);
x=x'.*mask;
y=y'.*mask;
x=x';
y=y';
x=x(x~=0);
y=y(y~=0);

u=zeros(Mx,Nx);
Gi=zeros(Mx,Nx);

for l = 1:length(x)
   wi(1,l) = exp((-gamma/2)*((xi(:,l)-mu(:,k_hat(1,l)))'/(sigma(:,:,k_hat(1,l)) ) * (xi(:,l)-mu(:,k_hat(1,l))))); % weights for the patches
   zi=reshape(xi(:,l)',[tau tau]); 
   zi=zi.*wi(1,l);
   if mod(tau,2)==1
         u(x(l)-tau2:x(l)+tau2,y(l)-tau2:y(l)+tau2)=u(x(l)-tau2:x(l)+tau2,y(l)-tau2:y(l)+tau2)+zi;
         Gi(x(l)-tau2:x(l)+tau2,y(l)-tau2:y(l)+tau2)= Gi(x(l)-tau2:x(l)+tau2,y(l)-tau2:y(l)+tau2)+wi(1,l);
   else
        u(x(l)-tau2:x(l)+tau2-1,y(l)-tau2:y(l)+tau2-1)=u(x(l)-tau2:x(l)+tau2-1,y(l)-tau2:y(l)+tau2-1)+zi;
        Gi(x(l)-tau2:x(l)+tau2-1,y(l)-tau2:y(l)+tau2-1)= Gi(x(l)-tau2:x(l)+tau2-1,y(l)-tau2:y(l)+tau2-1)+wi(1,l);
   end
end
u=u./Gi;
u=u';
end