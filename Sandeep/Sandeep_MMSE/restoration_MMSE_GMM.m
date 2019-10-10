function [xi_hat,k_hat]=restoration_MMSE_GMM(y,w,mu,sigma,var,Sub,Subt,H,q,tau)


xi_hat=zeros(tau*tau*q*q,size(y,2));
K=length(w);
pk=zeros(1,size(y,2),K);
Rk=zeros(tau*tau,tau*tau,K);
phimu=zeros(tau*tau,K);
phi=Sub*H;
Dk=zeros(tau*tau*q*q,tau*tau,K);
Mvi = size(mu,1);
Mxi=tau*tau*q*q;
muHk= zeros(Mxi,K);
sigmaHk=zeros(Mxi,Mxi,K);
for k=1:K
muHk(:,k)=mu(1:Mxi,k);
sigmaHk(:,:,k)=sigma(1:Mxi,1:Mxi,k)+eye(tau*tau*q*q)*0;
Dk(:,:,k)=(Sub*(H*(sigmaHk(:,:,k))))';
Rk(:,:,k)=(var)+(Sub*(H*Dk(:,:,k)));
phimu(:,k)=phi*muHk(:,k);
end

model.Sigma=Rk;
model.mu=phimu;
model.w=w; 

[pk(1,:,:),~]=expectation(y(:,:),model);

%with the likelihood, we choose the better parameters for each patch in yi_obs
if K==1
   k_hat=ones(1,size(y,2)); 
else
gamma = zeros(size(y,2),K);
    for k=1:K
          gamma(:,k) = w(k)*mvnpdf(y',mu(Mxi+1:Mvi,k)',sigma(Mxi+1:Mvi,Mxi+1:Mvi,k));   
    end

 [~,k_hat]=max(gamma');
end

for k=1:K       
    g = muHk(:,k)+Dk(:,:,k)*pinv(Rk(:,:,k))*(y(:,:)-phimu(:,k));
    g=g.*pk(1,:,k);
    xi_hat= xi_hat +g;
end
% 
% x=pk.*nk;
% xi_hat(:,:)=sum(x,3);

end
