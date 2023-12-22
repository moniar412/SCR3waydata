function [U,TB,SV,Y,like,bic]=t2mixt(X,U,TB,eps,dis)
% 
% joint Tucker3 model & Gaussian mixture (ML solution)
%
% X = [X_1 ...X_K];   
%
% X  = U*G*kron(C,B)'
% Xy = B*Gy*kron(C,U)'
% Xz = C*Gz*kron(B,U)'
%
% 21/09/2016
%
d1 = size(X,1);
n1 = size(U,2);
[d2,n2] = size(TB);
SV = eye(d2);
su = sum(U);
p = su'/d1;
likeold = -Inf;
dif = 1;
it = 0;
Xbar = diag(1./su)*U'*X;
Y = Xbar*TB;
M = Y*TB';
Xm = zeros(d2,1);
Mm = zeros(d2,1);
lfig = zeros(d1,n1);
onesd1 = ones(d1,1);
od2=ones(d2,1);
%
%
while dif > eps,
   it=it+1;
   %
   % update SV
   SV=0;
   for g=1:n1
       Xd=X-onesd1*M(g,:);
       SV=SV+Xd'*diag(U(:,g))*Xd;
   end
   SV=SV/d1;
   %
   % update TB
   WB=Xbar'*diag(su)*Xbar;
   [P,L,Q]=svd(SV);
   %l=min(od2*10000,max(od2*0.0001,diag(L)));
   l=diag(L);
   SVR=P*diag(sqrt(l))*P';
   iSVR=P*diag(1./sqrt(l))*P';
   [UB,DB,VB]=svd(iSVR*WB*iSVR,0);
   TB=SVR*UB(:,1:n2);
   %
   % update Y
   STm1 = inv(SV);
   Y = Xbar*STm1*TB;
   M = Y*TB';
   %
   % update U
   for g=1:n1
       lfig(:,g)=-0.5*sum( ( (X-onesd1*M(g,:))*P*diag(1./sqrt(l)) ).^2,2);
   end
   ind=(lfig<-7e2);
   U=exp(-7e2*ind+lfig.*(1-ind))*diag(p);
   U=diag(1./sum(U'))*U;
   su=sum(U);
   Xbar=diag(1./su)*U'*X;
   %
   % update p
   p=su'/d1;
   %
   % stopping rule
   warning off
   UlU=U.*log(U);
   ulu=UlU(:);
   ulu(find(isnan(ulu)))=[];
   sulu=sum(ulu);
   warning on
   like = -0.5*d1*log(det(SV)) + d1*sum(p.*log(p)) + sum(sum(U.*lfig)) - sulu;
   dif=(like-likeold);
   likeold=like;
   %
end
if dif<-0.000001
    disp('-------------------------------------- error t2mixt --------------------------------------------')
    dif
    svd(SV)'
end
%
np = n1-1 + d2 + (n1-1)*n2 + (d2-n2)*n2 + (d2*d2+d2)/2 -1;
%
bic = 2*like - log(d1)*np;
if dis==1
    disp(sprintf('T2mix: dif=%g, iter=%g, like=%g, np=%g, BIC=%g',dif,it,like,np,bic))
end