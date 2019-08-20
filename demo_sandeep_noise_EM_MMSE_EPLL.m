% Sandeep method for the denoising. We use EPLL, MMSE with one component
clear all
close all
colormap gray;

%Y=SHX+p
%Y the LR image , X the HR image , p a white gaussian noise ,
%S a subsampling operator and H a blur operator 

cd(fileparts(mfilename('fullpath')));
addpathrec('.')
deterministic('on');

%% Parameters 
 img='cameraman';
sig=40;
h=37;
tau = 8;   % tau * tau the size of a patch in Y
q=1;       % magnification factor
K=300;       % the number of gaussian in the GMM
tau2=floor(tau/2);
tau2q=floor(tau*q/2);
epsilon = 0.000001; 
%% Load and generate Images

x_orig = double(imread('images/lena.png'))/255;
[Mx_orig, Nx_orig] = size(x_orig); 
x_orig = x_orig(1:floor(Mx_orig/q)*q,1:floor(Nx_orig/q)*q);% original image
x_orig(26:27,26:27)=1;
[Mx_orig, Nx_orig] = size(x_orig); 

% y_obs is the subsamping of the HR image x_orig (LR observed image)
sig    = sig/255;
y_obs      = x_orig + sig * randn([Mx_orig, Nx_orig]);

m = min(y_obs(:));
M = max(y_obs(:));
a = 5;
if m < 0 - a*sig || M > 1 + a*sig
    if autonorm
        warning(['Image was not in the range [0 1], ' ...
                 'it has been automatically normalized. ' ...
                 'Note that this may alter the performances. ' ...
                 'For optimal results, please normalized the input ' ...
                 'image correctly in the range [0 1].']);
        shift = m+a*sig;
        scale = (M-m - 2*a*sig);
        sig = sig / scale;
        y_obs     = (y_obs - shift) / scale;
    else
        warning(['Image is not in the range [0 1]. ' ...
                 'Note that this may alter the performances. ' ...
                 'For optimal results, please normalized the input ' ...
                 'image correctly in the range [0 1], or ' ...
                 'rerun GGMM_EPLL with option "autonorm: true"']);
    end
end

[My_obs, Ny_obs] = size(y_obs);
wM=My_obs/2;    
wN=Ny_obs/2; 
randM=1;
randN=1;
init=2;

y = y_obs(randM:randM+wM-1, randN:wN+randN-1);
[My, Ny] = size(y);
wM2 = wM*q;
wN2 = wN*q; % the size of part of the image X that we know
%x is the 'random' selection of a part of the HR image x_orig (HR observed image)
randM2=randM*q-(q-1);  
randN2=randN*q-(q-1);
x = x_orig(randM2:randM2+wM2-1, randN2:wN2+randN2-1);
[Mx, Nx] = size(x);



%% Patches collections

% [Mxi,Nxi]=size(xi);   %Mxi the number of pixels in a patch & Nxi the number of patches
% [Myi,~]=size(yi);     %we want to have the same number of patches

masky   = getmask(My, Ny, tau, 1, false);
yi = getpatches(y, tau, masky);
[Myi,~]=size(yi);
maskx   = getmask(Mx, Nx, tau*q, q, false);
xi = getpatches(x, tau*q, maskx);
[Mxi,Nxi]=size(xi);

vi= [xi;yi]; %the concatenate vector with the same number of patches
Mvi=Mxi+Myi; 

%% GMM parameters initialisation
% 
% muVk2=rand(Mxi,K);          %muVK = [muHk';muLk'] for HR and LR parts respectively
% sigmaVk2=zeros(Mxi,Mxi,K);   %sigmaVk = [sigmaHk sigmaHLk ; sigmaHLk' sigmaLk]

wVk=ones(1,K)*(K/Nxi);
hpatches=randi(Nxi,1,K);
vhnn = knnsearch(vi',vi(:,hpatches)','K',h,'Distance','euclidean');
muVk=zeros(Mvi,K);
sigmaVk=zeros(Mvi,Mvi,K);

for k=1:K
muVk(:,k)=mean(vi(:,vhnn(k,:)),2);    
sigmaVk(:,:,k)=(cov(vi(:,vhnn(k,:))'))+eye(Mvi)*epsilon;
end
model.mu=muVk;
model.Sigma=sigmaVk;
model.w=wVk;


%% EM with MMSE (Sandeep)
tstartEM = tic;
[~,model2,~]=mixGaussEm(vi,model); %EM
tEM = toc(tstartEM);
muVk=model2.mu;
sigmaVk= model2.Sigma;
wVk=model2.w;

tstart = tic;
masky_obs   = getmask(My_obs, Ny_obs, tau,1, false);
yi_obs = getpatches(y_obs, tau, masky_obs);
[Myi_obs,Nyi_obs]= size(yi_obs);
[xi_hat,k_hat] = patches_restoration(yi_obs,K,wVk,muVk,sigmaVk,Mxi); % MMSE

maskx_orig=zeros(My_obs*q,Ny_obs*q);
if mod(tau*q,2)==0
    maskx_orig(tau2q+1:q:My_obs*q-tau2q+1,tau2q+1:q:Ny_obs*q-tau2q+1)=1;
else
    maskx_orig(tau2q+1:q:My_obs*q-tau2q,tau2q+1:q:Ny_obs*q-tau2q)=1;
end
X_hat = get_image2D(xi_hat,tau*q,My_obs*q,Ny_obs*q,k_hat,muVk(1:Mvi-Myi_obs,:),sigmaVk(1:Mvi-Myi_obs,1:Mvi-Myi_obs,:),maskx_orig);
% 
% maskx_orig  = getmask(Mx_orig, Nx_orig, tau*q, q, randommask);
% [Mxi_orig,Nxi_orig]=size(xi);
% X_hat = projpatches(xi_hat,Mx_orig,Nx_orig,maskx_orig); %get image

tMMSE=tEM + toc(tstart);

%% EM with EPLL (Zoran & Weiss)

tstart=tic;
muHk= zeros(Mxi,K);
sigmaHk=zeros(Mxi,Mxi,K);
for k=1:K
muHk(:,k)=muVk(1:Mxi,k);
sigmaHk(:,:,k)=sigmaVk(1:Mxi,1:Mxi,k);
end
prior_model.name='gmm';
nu=ones(tau*tau*q*q,1)*2;
for k=1:K
   [U , S] = eig(sigmaHk(:,:,k));
   prior_model.GS.U{k}=U;
   prior_model.GS.S{k}=S(S~=0);
   prior_model.GS.wts(k)=wVk(k);
   prior_model.GS.mu{k}=muHk(:,k);
   prior_model.GS.nu{k}=nu;
end
prior_model.GS.nmodels=K;
prior_model.GS.dim=tau*tau*q*q;

X_hat2 = ggmm_epll(y_obs,maskx_orig,q,k_hat, sig, prior_model); % EPLL & get image

tEPLL=tEM + toc(tstart);

% X_hat3 = ;
% tFEPLL=toc(tstart);
%% Display & Tab
psnr_X_hat1=psnr(X_hat,x_orig)
ssim_X_hat1=ssim(X_hat,x_orig)
psnr_X_hat2=psnr(X_hat2,x_orig)
ssim_X_hat2=ssim(X_hat2,x_orig)
tstart=tic;
xb = imresize(y_obs, [Mx_orig, Nx_orig], 'bicubic');
txb=toc(tstart);
psnrxb=psnr(xb, x_orig)
ssimxb=ssim(xb, x_orig)

img=[img,'_K', num2str(K),'q',num2str(q),'t',num2str(tau)];
subplot(1,2,1);imagesc(x_orig,[0 1]);title('Image originale X');
subplot(1,2,2);imagesc(x,[0 1]);title('Partie de X (connue)');
colormap gray

figure(2)
subplot(2,2,1);imagesc(y_obs,[0 1]);title('Image observée Y (connue)');
subplot(2,2,2);imagesc(xb,[0 1]);
title(sprintf('Z: PSNR %.1f SSIM %.3f', psnr(xb, x_orig), ssim(xb, x_orig)));
subplot(2,2,3);imagesc(X_hat,[0 1]);title(sprintf('X_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat, x_orig), ssim(X_hat, x_orig)));
subplot(2,2,4);imagesc(X_hat2,[0 1]);title(sprintf('X_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2, x_orig), ssim(X_hat2, x_orig)));
colormap gray
tab_entier=[psnr_X_hat1, ssim_X_hat1, tMMSE; psnr_X_hat2, ssim_X_hat2, tEPLL; psnrxb, ssimxb,txb];

% set(2, 'PaperUnits', 'centimeters');
% set(2, 'PaperPosition', [0 0 25 20]);
% saveas(2,['fig/',img,'_large_',num2str(init)],'fig');
% saveas(2,['png/',img,'_large_',num2str(init)],'png');

figure()
subplot(2,2,1);imagesc(y,[0 1]);title('Partie de Y (connue)');
subplot(2,2,3);imagesc(X_hat(randM2:randM2+wM2-1, randN2:wN2+randN2-1),[0 1]);
title(sprintf('x_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x), ssim(X_hat(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x)));
subplot(2,2,2);imagesc(xb(randM2:randM2+wM2-1, randN2:wN2+randN2-1),[0 1]);
title(sprintf('z: PSNR %.1f SSIM %.3f',psnr(xb(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x), ssim(xb(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x)));
subplot(2,2,4);imagesc(X_hat2(randM2:randM2+wM2-1, randN2:wN2+randN2-1),[0 1]);
title(sprintf('x_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x), ssim(X_hat2(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x)));
colormap gray
tab_partconnue=[psnr(X_hat(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x), ssim(X_hat(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x); ...
psnr(X_hat2(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x), ssim(X_hat2(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x);...
psnr(xb(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x), ssim(xb(randM2:randM2+wM2-1, randN2:wN2+randN2-1), x)];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 20]);
% saveas(gcf,['fig/',img,'_partconnue_',num2str(init)],'fig');
% saveas(gcf,['png/',img,'_partconnue_',num2str(init)],'png');

figure()
posM=1;
posN=256;
%posM2=posM/q;  
posN2=posN/q;
% posM=1;
% posN=1;
posM2=1;
% posN2=1;
subplot(2,2,1);imagesc(y_obs(posM2:posM2+wM-1, posN2:wN+posN2-1),[0 1]);
title('Partie de Y (connue)');
subplot(2,2,3);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,2,2);imagesc(xb(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('z: PSNR %.1f SSIM %.3f',psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,2,4);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
colormap gray
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 20]);
% saveas(gcf,['fig/',img,'_part2_',num2str(init)],'fig');
% saveas(gcf,['png/',img,'_part2_',num2str(init)],'png');

tab_part2=[psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))];

figure()
posM=Mx_orig/2;
posN=Nx_orig/2;
posM2=posM/q;  
posN2=posN/q;
subplot(2,2,1);imagesc(y_obs(posM2:posM2+wM-1, posN2:wN+posN2-1),[0 1]);
title('Partie de Y (connue)');
subplot(2,2,3);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,2,2);imagesc(xb(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('z: PSNR %.1f SSIM %.3f',psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,2,4);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
colormap gray
tab_part3=[psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1));...
psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 20]);
% saveas(gcf,['fig/',img,'_part3_',num2str(init)],'fig');
% saveas(gcf,['png/',img,'_part3_',num2str(init)],'png');

figure()
posM=256;
posN=1;
posM2=posM/q;  
%posN2=posN/q;
% posM=1;
% posN=1;
% posM2=1;
 posN2=1;
subplot(2,2,1);imagesc(y_obs(posM2:posM2+wM-1, posN2:wN+posN2-1),[0 1]);
title('Partie de Y (connue)');
subplot(2,2,3);imagesc(X_hat(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (MMSE): PSNR %.1f SSIM %.3f',psnr(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,2,2);imagesc(xb(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('z: PSNR %.1f SSIM %.3f',psnr(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(xb(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
subplot(2,2,4);imagesc(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1),[0 1]);
title(sprintf('x_h (EPLL): PSNR %.1f SSIM %.3f',psnr(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1)), ssim(X_hat2(posM:posM+wM2-1, posN:wN2+posN-1), x_orig(posM:posM+wM2-1, posN:wN2+posN-1))));
colormap gray
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 25 20]);
% saveas(gcf,['fig/',img,'_part4_',num2str(init)],'fig');
% saveas(gcf,['png/',img,'_part4_',num2str(init)],'png');
