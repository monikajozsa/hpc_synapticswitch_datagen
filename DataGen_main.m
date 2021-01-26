%% The goal of this code is to find suitable equilibrium points and parameters for each relevant architecture from W_list = W_list(Nspecies,nIE)
% clear all
% clc

warning ('off','all')
Nspecies=3;
N_steps=500000;
N_realisations=100;

%% equilibrium point to be tested
xbar_min=1;
xbar_max=5;
xbar_interval=[xbar_min xbar_max];
N_xbar=10;

%% rate constants to be tested
rate_const_max=0.1;
rate_const_min=0.001;
rate_const_interval=[rate_const_min rate_const_max];

plot_on=0;
gen_max=100;
rep_max=100;

nE=0;
W_all = W_list_v2(Nspecies,nE,plot_on);
for network_i=1:1%length(W_all)
    FileName=strcat('/Users/mj555l/Dropbox (Cambridge University)/Monika-Tim/Synaptic switch paper/Random chemical systems/Data_rnd_sys/E',num2str(nE),'_A',num2str(network_i));
    FileName=strcat(FileName,'.mat');
    H_sparse=cell(N_xbar,rep_max);
    H_sparse_1D=cell(N_xbar,rep_max);
    W_sparse=cell(N_xbar,rep_max);
    rate_constants_all=cell(N_xbar,rep_max);
    xbar_all=cell(N_xbar,rep_max);
    for rep=1:rep_max
        rng(rep)
        tic;       
        %% System parameters and architecture
        W=W_all{network_i};
        for gen_i=1:gen_max
            [xbar,~,rate_constants] = GenRandChemReac_EI(Nspecies,xbar_interval,rate_const_interval,nE/length(W),W);
            [Neq, J] = Num_of_Equilibria(rate_constants,xbar);
            if all(Neq==1) && all(double(real(eig(J)))<0)
                break
            end
        end
        if exist('rate_constants','var')
            for x_k=1:N_xbar
                %% Simulation - change in length as xbar grows (it is multiplied by x_bar_ind)
                N_real=N_realisations+(x_k-1)*5;
                [X,T] = Gillespie_EI(xbar*x_k,rate_constants,N_steps,N_real);
                %% Sparse Distribution
                [H_sparse{x_k,rep},WW_sparse,H_sparse_1D{x_k,rep},max_XX]= Sparse_Distribution_EI(X,T);
                W_sparse{x_k,rep} = WW_sparse./sum(WW_sparse);
                rate_constants_all{x_k,rep}=rate_constants;
                xbar_all{x_k,rep}=xbar*x_k;
                clear('X','T','WW_sparse')
            end        
        else
            disp([nE network_i])
        end
        clear('rate_constants')
        toc
    end
    rep_start=1;
    save(FileName,'H_sparse','W_sparse','W','H_sparse_1D','rate_constants_all','xbar_all','N_steps','N_real','nE','network_i')
end