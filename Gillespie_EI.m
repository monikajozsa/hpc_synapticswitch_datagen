function [X,T] = Gillespie_EI(xbar,rate_constants,N_steps,N_realisations)

Nspecies=size(xbar,1);

%% Reactions
r=zeros(Nspecies,2*Nspecies);
for i=1:Nspecies
    r(i,(i-1)*2+1)=1;
    r(i,i*2)=-1;
end

%% Generalized I-E Toggle Switch with Gillespie's algorithm
t_0=0;
x_0=xbar;
T=zeros(1,N_steps,N_realisations);
X=zeros(size(r,1),N_steps,N_realisations);

parfor j=1:N_realisations
    rng(j+1)
    [X(:,:,j), T(:,:,j)] = Gillespie_alg_EI(r,xbar,rate_constants,t_0,x_0,N_steps);
end
end

function [X_realisation,T_realisation] = Gillespie_alg_EI(r,xbar,rate_constants,t_0,x_0,N_steps)

    T_realisation=zeros(1,N_steps);
    X_realisation=zeros(size(r,1),N_steps);
        
    t=t_0;
    x=x_0;
    T_realisation(:,1)=t;
    X_realisation(:,1)=x;

    for i=2:N_steps
        rates = Rates_EI(x,xbar,rate_constants);
        tau=-1/sum(rates)*log(rand);
        t=t+tau;
        reac=sum(cumsum(rates/sum(rates))<rand)+1;
        x=x+r(:,reac);
        T_realisation(:,i)=t;
        X_realisation(:,i)=x;
        if any(x>10000) && i<N_steps
            disp('x exceeded allowed threshold')
            break
        end
    end
end