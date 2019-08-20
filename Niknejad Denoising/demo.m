%This code is provided for non-comercial use, corrosponding to the following paper
%"Image Restoration Using Gaussian Mixture Models With Spatially Constrained Patch Clustering" 
%Milad Niknejad, Hossein Rabbani, Massoud Babaie-Zadeh, IEEE Transactions on Image Processing, Vol 24 ,  No. 11, PP. 3624 - 3636 Nov. 2015

%The values of some constant parameters has changed in this version of our code to have an overall better performance.


%A ETE MODIFIE POUR LE FLOU (ne marche pas)


clc;
clear;
tau=9;
j=1;
r=6; %jump size for reference patches
N=30; %size of clustering window
%im=imread('images/house.png');
x_orig=imread('cameraman.tif');
sigma=1; %noise standard deviation
sigma2=2;
x_orig=im2double(x_orig);
x_orig=x_orig*255;
%x_orig = x_orig(32:32+64-1-4,100:163-4);
%y_obs=imgaussfilt(x_orig,sigma2)+sigma*randn(size(x_orig));
l=5;%impair
G=fspecial('gaussian',l,sigma2);
H = eye(tau*tau);

y_obs=imfilter(x_orig,G,'conv')+sigma*randn(size(x_orig));
%y_obs=x_orig+sigma*randn(size(x_orig));
G2=zeros(tau,tau);
G2(tau/2-(l/2)+1:tau/2+(l/2),tau/2-(l/2)+1:tau/2+(l/2))=G;
G2=reshape(G2,tau*tau,1);
H(ceil(tau*tau/2),:)=G2;
for c=1:(floor(l/2))
   H(ceil(tau*tau/2)-c,:)=circshift(G2,-1*c);
    H(ceil(tau*tau/2)+c,:)=circshift(G2,1*c);
end
for c=1:(floor(l/2))
   H(ceil(tau*tau/2)-floor(l/2)-(tau-l)/2  -c*tau :ceil(tau*tau/2)+floor(l/2)+(tau-l)/2  - c*tau,:)= circshift(H(ceil(tau*tau/2)-floor(l/2)-(tau-l)/2:ceil(tau*tau/2)+floor(l/2)+(tau-l)/2,:),-(c*tau),2);
   H(ceil(tau*tau/2)-floor(l/2)-(tau-l)/2  +c*tau :ceil(tau*tau/2)+floor(l/2)+(tau-l)/2  + c*tau,:)= circshift(H(ceil(tau*tau/2)-floor(l/2)-(tau-l)/2:ceil(tau*tau/2)+floor(l/2)+(tau-l)/2,:),+(c*tau),2);
end



colormap gray

m=40;%number of cluster dimension
nmiter=5; %number of iterations of EM-like
%agrregation weights coeficients
gamma=0.1;
imn1=x_orig;

sigmam=125;
sigmaju=20;
sigmam=sigmam+sigmaju;
toiter=floor((sigmam-10)/sigmaju);


tstart = tic;
for iter=1:nmiter
    
 if (iter<toiter || iter==toiter)
        sigma=sigmam-iter*sigmaju;
    else
        sigma=10/(3^(iter-toiter));
    end

if mod(iter,3)==0
    r=5;
elseif mod(iter,3)==1
    r=6;
else
    r=5;
end
[X_hat, wc, Yd]=main(j,r,N,imn1,tau,m,iter,sigma,y_obs,gamma,H);

imn1=X_hat;
end
t = toc(tstart)

colormap gray
figure(2);
colormap gray
psnr2=10*log10(255*255/mean(mean((x_orig-X_hat).^2)))
ssim2=ssim_index(x_orig,X_hat)
subplot(1,3,1);imagesc(x_orig);title('image originale')
subplot(1,3,2);imagesc(y_obs);title('image bruitée')
subplot(1,3,3);imagesc(X_hat);title(sprintf('image reconstruite: PSNR %.1f SSIM %.3f',psnr2,ssim2 ));