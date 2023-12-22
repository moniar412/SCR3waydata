function simulaG5ari2
%
%
%
ns=250;
ARI=zeros(ns,3,16);
idx=[1:4;5:8;9:12;13:16];
for G=[5]
    for N=[500,1000]
        for nrep=[1,3]
            if N==500 && nrep==1
                r=1;
            elseif N==500 && nrep==3
                r=2;
              elseif N==1000 && nrep==1
                r=3;
            elseif N==1000 && nrep==3
                r=4;
                       end
            for dgp=1:4
                disp(sprintf('N=%g, G=%g, nrep=%g, dgp=%g',N,G,nrep,dgp))
                %Scenario 2
                ari=simula2(N,G,nrep,dgp,ns);
                %Scenario 1
               % ari=simula1(N,G,nrep,dgp,ns);
                idd=idx(r,dgp);
                ARI(:,:,idd)=ari;
                
                
               


            end
        end
    end
    %
   
    
    %save('xprova.mat','x')
 %save('G5ari2.txt','ARI','-ascii','-append')
 save('G5ari2.mat','ARI')
end

