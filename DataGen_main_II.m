%% This code generates data for 3D birth-death processes 
%% data is generated from 94 connectivity structures from 7 architectures (architecture is used for graphs when connection type is ignored)
%% the simulation is based on exact Gillespie's algorithm

clear all
clc

warning ('off','all')
Nspecies=3;
N_realisations=200; % number of realizations
N_steps=500000; % number of steps in each simulations
%% runtime: it was about 40 second for 1 connectivity with 10 realizations and 10000 steps, so for the whole thing: 40*100(connectivities)*10(100 realizations)*100(1000000 steps)/(3600*24) = it's about 40 days with 6 parallel processors

%% Equilibrium point to be tested
% xbar_min=2;
% xbar_max=20;
% N_xbar=10;
% xbar_vec=round(logspace(log10(xbar_min),log10(xbar_max),N_xbar+1));
xbar_vec=[1 3 5 7 10 13 16 20 24 28 33 38];
N_xbar=length(xbar_vec);


%% rate constants
rate_const_interval_vec=[0.1 0.07 0.05 0.03 0.01 0.007 0.005 0.003];
nrate_const=length(rate_const_interval_vec);

for nA=6:7
    rng(1)
    runtime=tic;  
    W_nA = W_list_v5(Nspecies,nA); % connectivity structure
    for nA_i=1:1%length(W_nA)
        W=W_nA{nA_i};
        %% Variable allocation for sparse distribution
        H_sparse=cell(N_xbar,nrate_const);
        H_sparse_1D=cell(N_xbar,nrate_const);
        W_sparse=cell(N_xbar,nrate_const);
        Weighted_Cov=cell(N_xbar,nrate_const);
        Weighted_Mean=cell(N_xbar,nrate_const);
        for k=1:nrate_const     
            %% System parameters and simulation
            rate_constants = Generate_Symm_Sys(Nspecies,rate_const_interval_vec(k),W);
            for xbar_ind=1:N_xbar
                xbar=xbar_vec(xbar_ind);
                [X,T] = Gillespie_EI(xbar*ones(Nspecies,1),rate_constants,N_steps,N_realisations);
                %% Calculating Sparse Distribution from time series
                [H_sparse{xbar_ind,k},WW_sparse,H_sparse_1D{xbar_ind,k},max_XX]= Sparse_Distribution_EI(X,T);
                W_sparse{xbar_ind,k} = WW_sparse./sum(WW_sparse); %normalization of the weights
                [Weighted_Cov{xbar_ind,k}, Weighted_Mean{xbar_ind,k}] = Sparse_Distribution_weighted_cov(H_sparse{xbar_ind,k},W_sparse{xbar_ind,k});
                clear('X','T','WW_sparse')
            end
            disp(strcat('nA=: ',num2str(nA),' out of 7. W: ',num2str(nA_i),' out of ',num2str(length(W_nA)),'. k=',num2str(k)))
            toc(runtime)
            FileName=strcat('Data_A',num2str(nA),'_W',num2str(nA_i),'.mat');
            save(FileName,'H_sparse','W_sparse','W','H_sparse_1D','Weighted_Cov','Weighted_Mean','rate_const_interval_vec','xbar_vec','N_steps','N_realisations')
        end
        %% Save variables
    end
end
