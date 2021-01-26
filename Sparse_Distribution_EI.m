function [H_sparse,W_sparse,H_sparse_1D,max_XX] = Sparse_Distribution_EI(X,T)
Xsize=size(X);
Nspecies=Xsize(1);
N_steps=Xsize(2);
N_realisations=Xsize(3);

N_discarded=max(1000,floor(0.05*N_steps));%has to be at least 1 because the corresponding time is shifted by one
time_int=diff(T);
XX=reshape(X(:,N_discarded:(N_steps-1),:),Nspecies,(N_steps-N_discarded)*N_realisations)';
WW=reshape(time_int(1,N_discarded:end,:),[],1);

XX=XX+1; %% It caused more confusion to start from 0 so we start from 1

%%STATIONARY DISTRIBUTION
% speeded up version for [H_sparse,~,IC]=unique(XX,'rows'); and sparser
% representation from (Nspecies-D) accumarray

% one-to-one convertion of XX to 1D vector using digits
max_XX=max(XX);
Ndigits=ceil(log10(max(10,abs(max_XX+1))));%% we need max_XX+1 because we will deal with the neaighbouring points
Npowers=flip(cumsum(flip(Ndigits)));
invNpowers=[-Npowers(2:end) 0];
Npowers=[Npowers(2:end) 0];
XX_powers=XX.*10.^Npowers;
XX_1D=sum(XX_powers,2);
% sparse distribution - storing only the visited grid points (in 1D)
[H_sparse_1D,~,IC]=unique(XX_1D);
% Converting 1D sparse distribution back to Nspecies-D
H_sparse=zeros(length(H_sparse_1D),Nspecies);
H_sparse_1D_temp=H_sparse_1D;
for i=1:Nspecies
    coord_i=floor(H_sparse_1D_temp*10^invNpowers(i));
    H_sparse(:,i)=coord_i;
    H_sparse_1D_temp=H_sparse_1D_temp-10^Npowers(i)*coord_i;
end

%time spent in each state
W_sparse = accumarray(IC,WW);
W_sparse = W_sparse./sum(W_sparse);
end
