function rate_constants = Generate_Symm_Sys(Nspecies,rate_const, W)
%This function generates system parameters for systems with symmetric rate constants for all
%connections

connections=find(abs(W)>0);

%birth rates
Knum_E=zeros(size(W));
Kdenum_E=zeros(size(W));
Knum_I=zeros(size(W));
Kdenum_I=zeros(size(W));
for j=1:size(connections,1)
    i=connections(j);
    if W(i)==1%excitatory
        Knum_E(i) = rate_const;
        Kdenum_E(i) = rate_const;
    else %inhibitory
        Knum_I(i) = rate_const;
        Kdenum_I(i) = rate_const;
    end
end
rate_constants.Knum_E=Knum_E;
rate_constants.Kdenum_E=Kdenum_E;
rate_constants.Knum_I=Knum_I;
rate_constants.Kdenum_I=Kdenum_I;

%% death rates - determined by the rate_constants 
beta=zeros(1,Nspecies);
for i=1:Nspecies
    beta(i)=sum([rate_constants.Knum_E(:,i);rate_constants.Knum_I(:,i)]./[rate_constants.Kdenum_E(:,i)+1;rate_constants.Kdenum_I(:,i)+1]);
end

rate_constants.betas=beta;