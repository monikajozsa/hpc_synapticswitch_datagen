function [xbar,W,rate_constants] = GenRandChemReac_EI(Nspecies,xbar_interval,rate_const_interval,IE_ratio, W_force)

if ~exist('Nspecies','var')
    Nspecies=2;
end
if ~exist('rate_const_interval','var')
    rate_const_interval=[0.1 0.5];
end
rate_const_intervals.Knum_I=rate_const_interval;
rate_const_intervals.Kdenum_I=rate_const_interval;
rate_const_intervals.Kdenum_E=rate_const_interval;

if ~exist('xbar_interval','var') || isempty(xbar_interval)
    xbar_interval=[20 40];
end

if ~exist('IE_ratio','var')
    IE_ratio=1;
end

%% Connectivity matrix
if exist('W_force', 'var')
    W=W_force;
else
    W=randi([0 4],Nspecies)-2;%randi([0 2],Nspecies)-1;%connectivity matrix: if W(i,j)==-1 then i inhibits j; if W(i,j)==1 then i excites j
    W=sign(W);
    for i=1:Nspecies
        W(i,i)=0;
    end
    %% to avoid complete death of a species (and thus the break of the equilibrium)
    if any(sum(abs(W))==0)
        die_out_vec=find(sum(abs(W))==0);
        for i=die_out_vec
            forced_neighbour=randi([1 Nspecies],1);
            while forced_neighbour==i
                forced_neighbour=randi([1 Nspecies],1);
            end
            W(forced_neighbour,i)=2*(randi([0 1],1)-0.5);
        end
    end

    if exist('IE_ratio','var')
        W=abs(W);
        Connections=find(W==1);
        nC=size(Connections,1);
        nI=round(IE_ratio*nC);
        I_ind=randperm(nC,nI)';
        W(Connections(I_ind))=-ones(size(I_ind));
    end   
end
    
%% generating the rate constants for I and E connections

connections=find(abs(W)>0);
% self_connections=find(eye(Nspecies)==1);
% connections=setdiff(connections,self_connections);
Kdenum_E=zeros(size(W));
Knum_I=zeros(size(W));
Kdenum_I=zeros(size(W));
for j=1:size(connections,1)
    i=connections(j);
    if W(i)==1%excitatory connections
        Kdenum_interval=rate_const_intervals.Kdenum_E;
        Kdenum_E(i) = Kdenum_interval(1)+rand(1)*diff(Kdenum_interval);
    else %inhibitory connections
        Knum_interval=rate_const_intervals.Knum_I;
        Kdenum_interval=rate_const_intervals.Kdenum_I;
        Knum_I(i) = Knum_interval(1)+rand(1)*diff(Knum_interval);
        Kdenum_I(i) = Kdenum_interval(1)+rand(1)*diff(Kdenum_interval);
    end
end

rate_constants.Kdenum_E=Kdenum_E;
rate_constants.Knum_I=Knum_I;
rate_constants.Kdenum_I=Kdenum_I;

%% this symmetry in the inhibitory constants is added as a simplification
rate_constants.Kdenum_I=rate_constants.Knum_I;

%% xbar
for i=1:Nspecies
    xbar=randi(xbar_interval,Nspecies,1);
end

