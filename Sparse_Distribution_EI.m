function [H_sparse,W_sparse,H_sparse_1D,max_XX] = Sparse_Distribution_EI(X,T)
%% Sparse distribution (only stores points with non-zero weight) from time series

%% Discard the beginning of the time series
Xsize=size(X);
Nspecies=Xsize(1);
N_steps=Xsize(2);
N_discarded=max(1000,floor(0.01*N_steps));
if length(Xsize)>2
    N_realisations=Xsize(3);
else
    N_realisations=1;
end
time_int=diff(T);

try %% we use "try" because if there is only one realization then we do not need to reshape X and time_int
    XX=reshape(X(:,N_discarded:(N_steps-1),:),Nspecies,(N_steps-N_discarded)*N_realisations)';
    WW=reshape(time_int(1,N_discarded:end,:),[],1);
catch
    XX=X(:,N_discarded:(N_steps-1))';
    WW=time_int(1,N_discarded:end)';
end

%% Discard points with NaN values - it occurs when time series dies off or reaches upper limit
nan_ind=isnan(XX);
XX=reshape(XX(~nan_ind),[],Nspecies);
WW=WW(~nan_ind(:,1));

XX=XX+1; %% it caused more confusion to start from 0 so we start from 1

%% One-to-one convertion of XX to XX_1D, a one dimensional vector - speeds up "unique" function
% example: if XX(i,:)=[2,3,4] and the maximum coordinates are two digits in each dimensions of X then XX_1D(i)=20304
max_XX=max(XX);
Ndigits=ceil(log10(max(10,abs(max_XX+1))));%% we need max_XX+1 because we will deal with the neighbouring points later
Npowers=flip(cumsum(flip(Ndigits)));
invNpowers=[-Npowers(2:end) 0];
Npowers=[Npowers(2:end) 0];
XX_powers=XX.*10.^Npowers;
XX_1D=sum(XX_powers,2);

%% 1-dimensional sparsely stored support of the distribution
[H_sparse_1D,~,IC]=unique(XX_1D); % store only the visited grid points (in 1D)
%% Nspecies-dimensional sparsely stored support of the distribution
H_sparse=zeros(length(H_sparse_1D),Nspecies);% convert 1D sparse distribution back to Nspecies-dimension
H_sparse_1D_temp=H_sparse_1D;
for i=1:Nspecies
    coord_i=floor(H_sparse_1D_temp*10^invNpowers(i));
    H_sparse(:,i)=coord_i;
    H_sparse_1D_temp=H_sparse_1D_temp-10^Npowers(i)*coord_i;
end

%% Probability weight = time spent in a state
W_sparse = accumarray(IC,WW);
W_sparse = W_sparse./sum(W_sparse);
end
