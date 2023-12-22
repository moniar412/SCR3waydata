function [U,TB,TC,SO,SV,Y,like,bic]=t3mixs(X,U,TB,TC,SV,SO,eps,dis)
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
[d3,n3] = size(TC);
su = sum(U);
p = su'/d1;
likeold = -Inf;
dif = 1;
it = 0;
Xbar = diag(1./su)*U'*X;
M = Xbar;
% Y = Xbar*kron(TC,TB);
% [TB,TC,Y]=tuck2a(Xbar,TB,TC,1,0); % rational starting point
% M = Y*kron(TC',TB');
Xm = zeros(d2,d3);
Mm = zeros(d2,d3);
lfig = zeros(d1,n1);
onesd1 = ones(d1,1);
%
%
while dif > eps,
   it=it+1;
   %
   % update SO
   A=0;
   SVm1=inv(SV);
   for i=1:d1
       Xm(:)=X(i,:);
       for g=1:n1
           Mm(:)=M(g,:);
           A=A+U(i,g)*(Xm-Mm)'*SVm1*(Xm-Mm);
       end
   end
   SO=A/(d1*d2);
   %
   % update SV
   A=0;
   SOm1=inv(SO);
   for i=1:d1
       Xm(:)=X(i,:);
       for g=1:n1
           Mm(:)=M(g,:);
           A=A+U(i,g)*(Xm-Mm)*SOm1*(Xm-Mm)';
       end
   end
   SV=A/(d1*d3);
   %
   % update TB
   WB=0;
   SOm1=inv(SO);
   TCTCp=SOm1*TC*TC'*SOm1;
   for g=1:n1,
       Xm(:)=Xbar(g,:);
       WB=WB+su(g)*Xm*TCTCp*Xm';
   end
   [P,L,Q]=svd(SV);
   SVR=P*diag(sqrt(diag(L)))*P';
   iSVR=P*diag(1./sqrt(diag(L)))*P';
   [UB,DB,VB]=svd(iSVR*WB*iSVR,0);
   TB=SVR*UB(:,1:n2);
   %
   %update TC
   WC=0;
   SVm1=inv(SV);
   TBTBp = SVm1*TB*TB'*SVm1;
   for g=1:n1,
       Xm(:)=Xbar(g,:);
       WC=WC+su(g)*Xm'*TBTBp*Xm;
   end
   SOR=sqrtm(SO);
   iSOR=inv(SOR);
   [UC,DC,VC]=svd(iSOR*WC*iSOR,0);
   TC=SOR*UC(:,1:n3);
   %
   % update Y (eta)
   iSO=inv(SO);
   iSV=inv(SV);
   Y = Xbar*kron(iSO*TC,iSV*TB);
   M = Y*kron(TC',TB');
   %
   % update U
   iSORiSVR=kron(iSOR,iSVR);
   for g=1:n1
       lfig(:,g)=-0.5*sum( ( (X-onesd1*M(g,:))*iSORiSVR ).^2,2);
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
   like = -0.5*d1*d2*log(det(SO)) -0.5*d1*d3*log(det(SV))+ d1*sum(p.*log(p)) + sum(sum(U.*lfig)) - sulu;
   dif=(like-likeold);
   likeold=like;
   %
end
if dif<-0.000001
    disp('------------------------------------------- error -------------------------------------------------')
    dif
end
%
np = n1-1 + d2*d3 + (n1-1)*n2*n3 + (d2-n2)*n2+(d3-n3)*n3 + (d2*d2+d2)/2 + (d3*d3+d3)/2 -1;
%
bic = 2*like - log(d1)*np;
if dis==1
    disp(sprintf('T3mix: dif=%g, iter=%g, like=%g, np=%g, BIC=%g',dif,it,like,np,bic))
end