function [X,U,z]=genmixhet(n,p,M,Sig)
%
% generate a random sample from a mixture of Gaussians with
% prior probabilities p, matrix means M and variance matrices Sig.
%
%
[G,J]=size(M);
U=zeros(n,G);
onesn=ones(n,1);
X=zeros(n,J);
%
% generate the number of observations for each class
cp=cumsum(p);
x=rand(n,1);
z=1+sum((x*ones(1,G))>(ones(n,1)*cp'),2);
%
% generate the data
for g=1:G
    ind=find(z==g);
    ng=size(ind,1);
    X(ind,:)=ones(ng,1)*M(g,:)+randn(ng,J)*sqrtm(Sig((1:J)+J*(g-1),:));
end
for g=1:G
    [P,L,Q]=svd(Sig((1:J)+J*(g-1),:));
    lam=diag(L);
    U(:,g)=-0.5*sum(log(lam))-0.5*sum(((X-onesn*M(g,:))*Q*diag(1./sqrt(lam))).^2,2);
end
ind=(U<-7.0e2);
U=exp(-7.0e2*ind+U.*(1-ind))*diag(p);
U=diag(1./sum(U'))*U;