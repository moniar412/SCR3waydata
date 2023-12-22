function [U,Mmu,Sig,dif,like,bic,it]=mixhom(X,U,eps,dis)
% 
% homoscedastic Gaussian mixture (ML penalizeded IW)
%
% 05/03/2008
%
[n,J]=size(X);
G=size(U,2);
su=sum(U);
onesn=ones(n,1);
oJ=ones(J,1);
lfig=zeros(n,G);
p=(su/n)';
dif=1;
likeold=-Inf;
it=0;
%
while dif > eps
   it=it+1;
   %
   % update Mmu
   Mmu = diag(1./su)*U'*X;
   %
   % update MOme and sig
   Sig=zeros(J);
   for g=1:G
       Xd=X-onesn*Mmu(g,:);
       Sig=Sig+Xd'*diag(U(:,g))*Xd;
   end
   Sig=Sig/n;
   [Q,L,Q]=svd(Sig);
   Lam1=diag(1./sqrt(diag(L)));
   %
   lam=diag(L);
   MOmeQ=Q;
   %lam=min(oJ*10000,max(oJ*0.0001,lam));
   % update U
   for g=1:G
       lfig(:,g)=-0.5*sum(log(lam))-0.5*sum(((X-onesn*Mmu(g,:))*Q*diag(1./sqrt(lam))).^2,2);
   end
   ind=(lfig<-7.0e2);
   U=exp(-7.0e2*ind+lfig.*(1-ind))*diag(p);
   U=repmat(1./sum(U,2),1,G).*U;
   su=sum(U);
   %
   % update p
   p=(su/n)';
   %
   % stopping rule
   warning off
   UlU=U.*log(U);
   ulu=UlU(:);
   ulu(find(isnan(ulu)))=[];
   sulu=sum(ulu);
   warning on
   like=n*sum(p.*log(p))+sum(sum(U.*lfig))-sulu;
   dif=(like-likeold);
   likeold=like;
   %
end
np=G-1+G*J+(J*J+J)*0.5;
%awe=2*like+2*sulu-2*(3/2+log(n))*np;
awe=1;
bic=2*like-log(n)*np;
aic=2*like-2*np;
if dis==1
    disp(sprintf('mixhom(%g): dif=%g, iter=%g, like=%g, AWE=%g, BIC=%g, AIC=%g',G,dif,it,like,awe,bic,aic))
end
if dif<-0.000001
    disp('------------------------------------------- error -------------------------------------')
    dif
    MOmel
    10000*min(MOmel)
    sum(U)
end 