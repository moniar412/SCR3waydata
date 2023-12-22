function [ari,Xt,Utruet]=simula1(N,G,nrep,dgp,ns)
%
% G=#groups, N=#observations, nrep=#starting points, dgp=data generating process
%
%ns=100;
eps=10^-6;
% variables
J=5;
% variables components
Q=2;
% occasions
K=4;
% occasions components
R=2;
%
ari=zeros(ns,3);
Xt=zeros(N,J*K,ns);
Utruet=zeros(N,G,ns);
parfor sa=1:ns
    rng(sa+13);
    pg=rand(G,1); pg=pg/sum(pg);
    Sv=randn(J); Sv=Sv'*Sv; Sv(1:Q,Q+1:J)=0; Sv(Q+1:J,1:Q)=0;
    So=randn(K); So=So'*So; So(1:R,R+1:K)=0; So(R+1:K,1:R)=0;
    S=kron(So,Sv); % model true
    Sf=0.6*randn(J*K); Sf=Sf'*Sf; Sf=(S~=0).*Sf; % model false
    B=eye(J); B=B(:,1:Q);
    C=eye(K); C=C(:,1:R);
    Eta=randn(Q*R,G);
    Mu=(20*kron(C,B)*Eta)'; % model true
    Muf=20*(Mu~=0).*randn(G,J*K); % model false
    switch dgp
        case 1 % model true
        case 2 % model true in S
            Mu=Muf;
        case 3 % model true in Mu
            S=Sf;
        case 4 % model false
            Mu=Muf; S=Sf;
    end
    [X,Utrue,z]=genmixhet(N,pg,Mu,repmat(S,G,1));
    X=preproa(X,1);
    for model=1:3
        blike=-Inf;
        for rep=1:nrep
            rng(10*sa+rep);
            U=rand(N,G); U=diag(1./sum(U'))*U;
            if rep==1
%                 [U,M,~,~,like,~,~]=mixhom(X,U,eps,0);
                %[U,TB,M,ssr]=rkm(X,U,Q*R,0);
                 [U,~,~]=kmeans(X,U,0,0);
%                  U=Utrue;
            end
            switch model
                case 1  % T3
                    Sv=randn(J); Sv=Sv'*Sv; So=randn(K); So=So'*So;
                    TB=rand(J,Q);
                    [A,D,B]=svd(TB'*inv(Sv)*TB); TB=TB*B*diag(1./sqrt(diag(D)));
                    TC=    rand(K,R);
                    [A,D,B]=svd(TC'*inv(So)*TC); TC=TC*B*diag(1./sqrt(diag(D)));
                    %
                    
                    [Ur,TB,TC,~,~,Y,like,~]=t3mixs(X,U,TB,TC,Sv,So,eps,0);
                    %Ur=U; like=1;
                case 2  % T2
                    TB=orth(rand(J*K,Q*R));
                    [Ur,TB,~,Y,like,~]=t2mixt(X,U,TB,eps,0);
                    %Ur=U; like=1;
                case 3  % mixhom
                   
                    [Ur,M,~,~,like,~,~]=mixhom(X,U,eps,0);
                    %Ur=U; like=1;
            end
            if blike<like
                blike=like; bU=Ur;
            end
        end
        ari(sa,model)=mrand(ftoh(Utrue)'*ftoh(bU));
        Xt(:,:,sa)=X;
        Utruet(:,:,sa)=Utrue;
        % los(sa,model)=lossu(Utrue,bU);
    end
    %disp(sa)
    %monia5=[ari(sa,:)  ]
    %save('ARI.txt','X','-ascii','-append')
end
% display results
a=ari;
disp('---------------------------------')
disp('    T3mix     T2mix     Hom')
disp(mean(a,1)), disp(std(a)), disp(prctile(a,25)),  disp(prctile(a,50)), disp(prctile(a,75))
