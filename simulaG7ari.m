function simulaG7ari
%
%
%
ns=250;
ARI=zeros(ns,3,16);
idx=[1:4;5:8;9:12;13:16];
for G=[7]
    for N=[300,500]
        for nrep=[1,3]
           if N==300 && nrep==1
                r=1;
            elseif N==300 && nrep==3
                r=2;
              elseif N==500 && nrep==1
                r=3;
            elseif N==500 && nrep==3
                r=4;
                       end
            for dgp=1:4
                disp(sprintf('N=%g, G=%g, nrep=%g, dgp=%g',N,G,nrep,dgp))
                %Scenario 2
               % ari=simula2(N,G,nrep,dgp,ns);
                %Scenario 1
                ari=simula1(N,G,nrep,dgp,ns);
                idd=idx(r,dgp);
                ARI(:,:,idd)=ari;
                
               


            end
        end
    end
    %
   
    
    %
 save('G7ari.mat','ARI')
end

