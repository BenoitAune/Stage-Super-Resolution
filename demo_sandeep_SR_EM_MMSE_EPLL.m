% Sandeep method for the Super resolution. We use EPLL, MMSE with one
% component and MMSE with all components
% I compile with differents K numbers of components and differents q
% magnification factors.
% I compile 10 times to have stable results and for 5 images with size
% 512*512



cd(fileparts(mfilename('fullpath')));
addpathrec('.')
deterministic('on');

K=[100, 200, 300];
q=[2,3];
tau=[4,3];
init = 'inithg_MMSEallcomp';

moyenne=10;
for im=1:5
    
tab_moy_entier=zeros(3*size(q,2),3*size(K,2));
tab_moy_parthg=zeros(3*size(q,2),2*size(K,2));
tab_moy_parthd=zeros(3*size(q,2),2*size(K,2));
tab_moy_partbg=zeros(3*size(q,2),2*size(K,2));
tab_moy_partbd=zeros(3*size(q,2),2*size(K,2));
for i=1:moyenne
tab_final_entier=[];
tab_final_parthg=[];
tab_final_parthd=[];
tab_final_partbg=[];
tab_final_partbd=[];
for t=1:2
if im==1
        
    img='camera';
    ext='.jpg';
    x_orig = rgb2gray(double(imread(['images/',img,ext]))/255);
    else
        if im==2
        img='lena';
        ext='.tif';
        x_orig = rgb2gray(double(imread(['images/',img,ext]))/255);
        else
            if im==3
            img='barbara';
            ext='.png'; 
            x_orig =double(imread(['images/',img,ext]))/255;
            else
                if im==4
                img='hill';
                ext='.png';
                x_orig =double(imread(['images/',img,ext]))/255;
                else
                img='shape';
                x_orig=ones(512,512)*1;
                x_orig(50:150,50:150)=0;
                x_orig(350:450,350:450)=0;
                for s=0:0.01:2*pi
                   x_orig(-round(150/2*cos(s))+100:round(150/2*cos(s))+100,-round(150/2*sin(s))+400:round(150/2*sin(s))+400)=0;
                end
                for s=50:50:150
                   x_orig(300:500,s-12:s+12)=0; 
                   x_orig(s-12+300:s+12+300,20:200)=0;
                end   
                end
            end
        end
end
[Mx_orig, Nx_orig] = size(x_orig); 
x_orig = x_orig(1:floor(Mx_orig/q(t))*q(t),1:floor(Nx_orig/q(t))*q(t));% original image
tab1=[];
tab2=[];
tab3=[];
tab4=[];
tab5=[];
for k=K
[tab_entier,tab_parthg,tab_parthd,tab_partbg,tab_partbd]=main_sandeep_SR_EM_MMSE_EPLL(k,q(t),tau(t),x_orig,img,init); %%%%%%
tab1=cat(2,tab1 , tab_entier);
tab2=cat(2,tab2 , tab_parthg);
tab3=cat(2,tab3 , tab_parthd);
tab4=cat(2,tab4 , tab_partbg);
tab5=cat(2,tab5 , tab_partbd);
end
tab_final_entier=cat(1,tab_final_entier, tab1);
tab_final_parthg=cat(1,tab_final_parthg,tab2);
tab_final_parthd=cat(1,tab_final_parthd,tab3);
tab_final_partbg=cat(1,tab_final_partbg,tab4);
tab_final_partbd=cat(1,tab_final_partbd,tab5);
end
tab_moy_entier=tab_moy_entier+tab_final_entier;
tab_moy_parthg=tab_moy_parthg+tab_final_parthg;
tab_moy_parthd=tab_moy_parthd+tab_final_parthd;
tab_moy_partbg=tab_moy_partbg+tab_final_partbg;
tab_moy_partbd=tab_moy_partbd+tab_final_partbd;
end
tab_moy_entier=tab_moy_entier/moyenne;
tab_moy_parthg=tab_moy_parthg/moyenne;
tab_moy_parthd=tab_moy_parthd/moyenne;
tab_moy_partbg=tab_moy_partbg/moyenne;
tab_moy_partbd=tab_moy_partbd/moyenne;

csvwrite(['tab/tab_',img,'_entier_',init,'.csv'],tab_moy_entier ) ;
csvwrite(['tab/tab_',img,'_parthg_',init,'.csv'],tab_moy_parthg) ;
csvwrite(['tab/tab_',img,'_parthd_',init,'.csv'],tab_moy_parthd) ;
csvwrite(['tab/tab_',img,'_partbg_',init,'.csv'],tab_moy_partbg) ;
csvwrite(['tab/tab_',img,'_partbd_',init,'.csv'],tab_moy_partbd) ;
 end
