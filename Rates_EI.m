function rates = Rates_EI(x,xbar,rate_constants)
%% This function calculates the reaction rates based on the current state x, the equilibrium xbar, and the birth rates (defined by inhibitory, exitatory connections) and death rates rate_constants
Nspecies=length(xbar);

Knum_E=rate_constants.Knum_E;
Kdenum_E=rate_constants.Kdenum_E;
Knum_I=rate_constants.Knum_I;
Kdenum_I=rate_constants.Kdenum_I;
beta=rate_constants.betas;

Exc_rates=sum(diag((x+1e-5)./xbar)*Knum_E*diag(xbar)./(Kdenum_E+x./xbar.*(Kdenum_E>0)+1*(Kdenum_E==0)));
Inh_rates=sum(Knum_I*diag(xbar)./(Kdenum_I+x./xbar.*(Kdenum_I>0)+1*(Kdenum_I==0)));
Death_rates=beta.*x';

rates=zeros(2*Nspecies,1);
for i=1:Nspecies
    rates((i-1)*2+1:i*2)=[Exc_rates(i)+Inh_rates(i);Death_rates(i)];
end
end