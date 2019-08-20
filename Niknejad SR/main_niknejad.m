function [X_hat2]= main(j,r,N,X_hat,tau,m,iter,d,sig,y_obs,gamma,q,H,Subt,Sub)

%input: j is the jump size for the patches 
%       r is the jump size for the references patches 
%       N*N is the size of each region
%       X_hat is the reconstructed image at the (iter-1)-th iteration
%       tau*tau is the size of the patches
%       m is the number of cluster dimension
%       d is the number of iterations
%       sig is the variance of noise
%       y_obs is the observed image
%       gamma is a weights coeficients for the image reconstruction
%       q is the magnification factor
%       H is the blur operator
%       Sub & Subt are the sub sampling operators. (Sub to shrink & Subt to extend)
%output: X_hat2 is the image at the iter-th iteration 

dim=size(y_obs,1);
dims2=size(y_obs,2);
X_hat2=zeros(size(X_hat)); % the futur new image 
w=zeros(size(X_hat));      % the weight for the image reconstruction

for rowp=1:r:size(y_obs,1)
    iter/d
    if (size(y_obs,1)-rowp)<tau 
        rowp=size(y_obs,1)-tau+1;end 
    for colp=1:r:size(y_obs,2)
        if (size(y_obs,2)-colp)<tau     
            colp=size(y_obs,2)-tau+1;end  % (rowp,colp) is the position of the reference patch
        refp=X_hat(rowp*q-(q-1):rowp*q-(q-1)+tau*q-1,colp*q-(q-1):colp*q-(q-1)+tau*q-1); %the reference patches in the reconstructed imageX_hat
        refpa=y_obs(rowp:rowp+tau-1,colp:colp+tau-1);                                    %the reference patches in the observed image Y 
        refpv=reshape(refp,[tau*q*q*tau,1]);                                             % in column
        refpva=reshape(refpa,[tau*tau,1]);
        [sr,sc,er,ec]=winsize(rowp,colp,N,dim,dims2); % to have the position at the beginning and at the end of the region
        dimc=0;  
        cs=0;
        for opr=sr:j:er-tau+1
            for opc=sc:j:ec-tau+1
                cs=cs+1; % the number of patches in the region
            end
        end
        Y1=zeros(tau*tau*q*q,cs);       % all the patches of the reconstructed image in the region
        Y1a=zeros(tau*tau,cs);          % all the patches of the observed image in the region
        pos1=zeros(2,cs);               % the positions of the patches in the image
        for opr=sr:j:er-tau+1
            for opc=sc:j:ec-tau+1
                dimc=dimc+1;
                pa=X_hat(opr*q-(q-1):opr*q-(q-1)+tau*q-1,opc*q-(q-1):opc*q-(q-1)+tau*q-1);  
                paa=y_obs(opr:opr+tau-1,opc:opc+tau-1);             
                y=reshape(pa,[tau*tau*q*q,1]);                      % in column
                ya=reshape(paa,[tau*tau,1]);
                
                Y1(:,dimc)=y;
                Y1a(:,dimc)=ya;                     
                pos1(:,dimc)=[opr;opc];
            end
        end
        Y=zeros(tau*tau*q*q,m+1);     %the m (knn) patches of the reconstructed image +  ref patch
        Ya=zeros(tau*tau,m+1);        %the m (knn) patches of the observed image + ref patch
        [Y,Ya,pos]=KNN(refpv,Y1,tau,pos1,m,Y1a);
        dimc2=size(Y,2);
        Y(:,dimc2+1)=refpv;
        Ya(:,dimc2+1)=refpva;
        pos(:,dimc2+1)=[rowp;colp];
        [Yd, wc]=MAP_niknejad(Y,Ya,sig,gamma,H,Subt,Sub); % parameters estimation & patches reconstruction
        clear Y;
        for i=1:size(Yd,2)
            %image reconstruction
            cor1=pos(1,i);
            cor2=pos(2,i);
            blk=reshape(Yd(:,i),[tau*q,tau*q]);
            X_hat2(cor1*q-(q-1):cor1*q-(q-1)+tau*q-1,cor2*q-(q-1):cor2*q-(q-1)+tau*q-1)=X_hat2(cor1*q-(q-1):cor1*q-(q-1)+tau*q-1,cor2*q-(q-1):cor2*q-(q-1)+tau*q-1)+blk.*wc(i,1);
            w(cor1*q-(q-1):cor1*q-(q-1)+tau*q-1,cor2*q-(q-1):cor2*q-(q-1)+tau*q-1)=w(cor1*q-(q-1):cor1*q-(q-1)+tau*q-1,cor2*q-(q-1):cor2*q-(q-1)+tau*q-1)+wc(i,1);

        end
        
             
          
    end
end

X_hat2=X_hat2./w;