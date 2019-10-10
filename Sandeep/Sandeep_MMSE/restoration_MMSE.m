function [xi_hat,k_hat] = restoration_MMSE(yi_obs,wVk,muVk,sigmaVk,Mxi) % Using MMSE

%input: wVk, muVk and sigmaVk the K GMM parameters
%       yi_obs the patches we use to choose the parameters
%       Mxi the number of pixels in the patches that we build
%output: xi_hat the new estimation of the patches
%        h_hat the position of the gaussian paramater that we choose for each patch

K= length(wVk);
[Myi_obs,Nyi_obs]=size(yi_obs);
Mvi = Myi_obs + Mxi;

%with the likelihood, we choose the better parameters for each patch in yi_obs
if K==1
   k_hat=ones(1,size(yi_obs,2)); 
else
gamma = zeros(Nyi_obs,K);
    for k=1:K
          gamma(:,k) = wVk(k)*mvnpdf(yi_obs(:,:)',muVk(Mxi+1:Mvi,k)',sigmaVk(Mxi+1:Mvi,Mxi+1:Mvi,k));   
    end

 [~,k_hat]=max(gamma');
end
 % we use MMSE (minimum mean square error) to build the patches of X_hat(estimation of X)
 xi_hat= zeros([Mxi,Nyi_obs]);
for i=1:Nyi_obs
    xi_hat(:,i) = muVk(1:Mxi,k_hat(i)) + sigmaVk(1:Mxi,Mxi+1:Mvi,k_hat(i))/(sigmaVk(Mxi+1:Mvi,Mxi+1:Mvi,k_hat(i)))*(yi_obs(:,i)-muVk(Mxi+1:Mvi,k_hat(i)));
end

end 
