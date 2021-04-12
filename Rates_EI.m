function rates = Rates_EI(x,xbar,rate_constants,eps_local)
%% This function calculates the reaction rates based on the current state x, the equilibrium xbar, and the birth rates (defined by inhibitory, exitatory connections) and death rates rate_constants
if ~exist('eps_local','var')
    eps_local=1e-4;
end
Nspecies=length(xbar);
Kdenum_E=rate_constants.Kdenum_E;
Knum_I=rate_constants.Knum_I;
Kdenum_I=rate_constants.Kdenum_I;
Knum_E=Kdenum_E>eps_local;

%% inhibitory and excitatory birth rates - epsilon is added to the excitatory birth rates to avoid complete death at state (0,0)
Exc_rates_mtx_num=((x./xbar+eps_local)).*Knum_E*diag(xbar);
Exc_rates_mtx_denum=Kdenum_E+x./xbar.*(Kdenum_E>0)+1*(Kdenum_E==0);
Exc_rates=sum(Exc_rates_mtx_num./Exc_rates_mtx_denum);
Inh_rates=sum(Knum_I*diag(xbar)./(Kdenum_I+x./xbar.*(Kdenum_I>0)+1*(Kdenum_I==0)));

%% Death rates: we add epsilon to the death rates so xbar is the equilibrium
beta=zeros(size(x))';
for i=1:Nspecies
    temp_vec=Kdenum_E(:,i)>eps_local;
    beta(i)=sum([(1+eps_local)*temp_vec;Knum_I(:,i)]./[Kdenum_E(:,i)+1;Kdenum_I(:,i)+1]);
end
Death_rates=beta.*x';

%% All birth-death rates
rates=zeros(2*Nspecies,1);
for i=1:Nspecies
    rates((i-1)*2+1:i*2)=[Exc_rates(i)+Inh_rates(i);Death_rates(i)];
end
end