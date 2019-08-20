% Niknejad method for the Super resolution.

clc;
clear;

cd(fileparts(mfilename('fullpath')));
addpathrec('.')
deterministic('on');

q=3; % the magnification factor
tau=7; % tau*tau the size of the patches
j=1;  % jump size for all the patches
r=6;  %jump size for reference patches
N=30;  %size of clustering window
img='lena';
% x_orig=ones(512,512)*1;
% x_orig(50:150,50:150)=0;
% x_orig(350:450,350:450)=0;
% for s=0:0.01:2*pi
%    x_orig(-round(150/2*cos(s))+100:round(150/2*cos(s))+100,-round(150/2*sin(s))+400:round(150/2*sin(s))+400)=0;
% end
% for s=50:50:150
%    x_orig(300:500,s-12:s+12)=0; 
%    x_orig(s-12+300:s+12+300,20:200)=0;
% end 
x_orig=imread(['images/',img,'.tif']);
%x_orig=imread(['images/',img,'.jpg']);
%x_orig=imread(['images/',img,'.png']);
sig=2; %noise standard deviation
% x_orig=im2double(x_orig);
x_orig=rgb2gray(im2double(x_orig));
x_orig=x_orig*255; % is the original image

[Mx_orig,Nx_orig]=size(x_orig);
x_orig = x_orig(1:floor(Mx_orig/q)*q,1:floor(Nx_orig/q)*q); % the size of x_orig must be divide by q
[Mx_orig,Nx_orig]=size(x_orig);
B= reshape(x_orig,Mx_orig*Nx_orig,1);

G=fspecial('gaussian',5,0.5); % the gaussian convolution to have the blur operator. Its size is odd
[Subx_orig,~]=sub_operator(q,Mx_orig/q); % the sub-sampling operator for the image
Hx_orig=blur_operator(G,Mx_orig); % the blur operator for the image

y_obs=Subx_orig*(Hx_orig*B);
y_obs=reshape(y_obs,Mx_orig/q,Nx_orig/q)+sig*randn(size(x_orig)/q); % is the observed image
colormap gray
c=40;  %number of cluster dimension
d=3;  %number of iterations
gamma=0.001; %agrregation weights coeficients

xr=imresize(y_obs, size(x_orig), 'bicubic'); %initialisation
imn1=xr;

m=tau*tau*q*q;
[Sub,Subt]=sub_operator(q,sqrt(m)/q);  % the sub-sampling operator for the patches
H=blur_operator(G,sqrt(m));           % the blur operator for the patches

tstart = tic;
for iter=1:d  

if mod(iter,3)==0
    r=5;
elseif mod(iter,3)==1
    r=6;
else
    r=5;
end
[X_hat]=main_niknejad(j,r,N,imn1,tau,c,iter,d,sig,y_obs,gamma,q,H,Subt,Sub);
                                                                                                                                                                                                                                                                                      
imn1=X_hat;
end
t = toc(tstart)

%% Display
colormap gray
figure(2);
colormap gray
psnr2=10*log10(255*255/mean(mean((x_orig-X_hat).^2)))
ssim2=ssim_index(x_orig,X_hat)
subplot(2,2,1);imagesc(x_orig,[0,255]);title('Image originale X')
subplot(2,2,2);imagesc(y_obs,[0,255]);title('Image observée Y')
subplot(2,2,4);imagesc(X_hat,[0,255]);title(sprintf('Image reconstruite X_h (LINC): PSNR %.1f SSIM %.3f',psnr2,ssim2 ));
subplot(2,2,3);imagesc(xr,[0,255]);title(sprintf('Image reconstruite Z (bicubic): PSNR %.1f SSIM %.3f',10*log10(255*255/mean(mean((x_orig-xr).^2))),ssim_index(x_orig,xr)));
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 20]);
saveas(gcf,['fig/',img,'_niknejad_q',num2str(q)],'fig');
saveas(gcf,['png/',img,'_niknejad_q',num2str(q)],'png');

tab_niknejad=[psnr2,ssim2,t];
csvwrite(['tab/tab_',img,'_nikenejad_q',num2str(q),'.csv'],tab_niknejad ) ;