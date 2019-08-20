function [X_hat2, wc,Yd]= main(j,r,N,X_hat,s,m,iter,sig,y_obs,gamma,G)
dim=size(X_hat,1);
dims2=size(X_hat,2);
X_hat2=zeros(size(X_hat));
w=zeros(size(X_hat));

for rowp=1:r:size(X_hat,1)
    clc;
    iter
    progress=rowp/dim
    if (size(X_hat,1)-rowp)<s
        rowp=size(X_hat,1)-s+1;end
    for colp=1:r:size(X_hat,2)
        if (size(X_hat,2)-colp)<s
            colp=size(X_hat,2)-s+1;end
        refp=X_hat(rowp:rowp+s-1,colp:colp+s-1);
        refpa=y_obs(rowp:rowp+s-1,colp:colp+s-1);
        refpv=reshape(refp,[s*s,1]);
        refpva=reshape(refpa,[s*s,1]);
        [sr,sc,er,ec]=winsizecal(rowp,colp,N,dim,dims2);
        dimc=0;  
        cs=0;
        for opr=sr:j:er-s+1
            for opc=sc:j:ec-s+1
                cs=cs+1;
            end
        end
        Y1=zeros(s*s,cs);
        Y1a=zeros(s*s,cs);
        pos1=zeros(2,cs);
        for opr=sr:j:er-s+1
            for opc=sc:j:ec-s+1
                dimc=dimc+1;
                pa=X_hat(opr:opr+s-1,opc:opc+s-1);
                paa=y_obs(opr:opr+s-1,opc:opc+s-1);
                y=reshape(pa,[s*s,1]);
                ya=reshape(paa,[s*s,1]);
                Y1(:,dimc)=y;
                Y1a(:,dimc)=ya;
                pos1(:,dimc)=[opr;opc];
            end
        end
        Y=zeros(s*s,m+1);
        Ya=zeros(s*s,m+1);
        [Y,Ya,pos]=KNN(refpv,Y1,s,pos1,m,Y1a);
        dimc2=size(Y,2);
        Y(:,dimc2+1)=refpv;
        Ya(:,dimc2+1)=refpva;
        pos(:,dimc2+1)=[rowp;colp];
        %debruitage
        [Yd wc]=lowrcovsimp2(Y,Ya,sig,gamma,G);
        clear Y;
        for i=1:size(Yd,2)
            cor1=pos(1,i);
            cor2=pos(2,i);
            blk=reshape(Yd(:,i),[s,s]);
            X_hat2(cor1:cor1+s-1,cor2:cor2+s-1)=X_hat2(cor1:cor1+s-1,cor2:cor2+s-1)+blk*wc(i,1);
            w(cor1:cor1+s-1,cor2:cor2+s-1)=w(cor1:cor1+s-1,cor2:cor2+s-1)+wc(i,1);
        end
        
                   
    end
end

X_hat2=X_hat2./w;