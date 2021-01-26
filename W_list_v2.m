function [W, break_ind] = W_list_v2(Nspecies,nE,plot_on)
%% This function gives back the relevant network architectures for Nspecies number of species and nIE ratio between inhibitory and exitatory connections
%W(i,j)=1 means that i excites j
if ~exist('plot_on','var')
    plot_on=0;
end

%% 2D
if Nspecies==2 
    W2D=cell(1);
    if nE==0
        W2D{1}=[0, -1; -1, 0];
        disp('The network consists of two I connections.')
    elseif nE==1
        W2D{1}=[0, 1; -1, 0];
        disp('The network consists of one I and one E connection.')
    elseif nE==2
        W2D{1}=[0, 1; 1, 0];
        disp('The network consists of two E connections.')
    end
    W=W2D;
    break_ind=1:3;
end

%% 3D
if Nspecies==3 && nE==0
    W3D=cell(6,1);
    W3D{1}=[0 -1 -1;-1 0 -1; -1 -1 0]; %full inhibitory network
    W3D{2}=[0 0 -1;-1 0 -1; -1 -1 0]; %one missing inhibition
    W3D{3}=[0 -1 0; -1 0 -1; -1 0 0]; %two missing inhibitions, circular inhibition
    W3D{4}=[0 -1 -1; -1 0 -1; 0 0 0]; %two missing inhibitions, no circular inhibition
    W3D{5}=[0 -1 0;-1 0 -1; 0 0 0]; %three missing inhibitions, no circular inhibition
    W3D{6}=[0 -1 0; 0 0 -1; -1 0 0]; %three missing inhibitions, circular inhibition
    W=W3D;
    break_ind=1:6;
    %%
    if plot_on
        figure
        for i=1:6
            subplot(3,2,i)
            G = digraph(W3D{i});
            title(strcat('Architectue:',num2str(i)))
            plot(G)
        end
    end
end

if Nspecies==3 && nE==1
    W3D=cell(16,1);
    %full inhibitory network - one I-E switch
    W3D{1}=[0 1 -1;-1 0 -1; -1 -1 0]; 
    %one missing inhibition - one I-E switch
    W3D{2}=[0 0 1;-1 0 -1; -1 -1 0]; 
    W3D{3}=[0 0 -1;1 0 -1; -1 -1 0]; 
    W3D{4}=[0 0 -1;-1 0 1; -1 -1 0]; 
    W3D{5}=[0 0 -1;-1 0 -1; 1 -1 0]; 
    W3D{6}=[0 0 -1;-1 0 -1; -1 1 0]; 
    %two missing inhibitions, circular inhibition - one I-E switch
    W3D{7}=[0 1 0; -1 0 -1; -1 0 0]; 
    W3D{8}=[0 -1 0; 1 0 -1; -1 0 0]; 
    W3D{9}=[0 -1 0; -1 0 1; -1 0 0]; 
    W3D{10}=[0 -1 0; -1 0 -1; 1 0 0]; 
    %two missing inhibitions, no circular inhibition - one I-E switch
    W3D{11}=[0 1 -1; -1 0 -1; 0 0 0]; 
    W3D{12}=[0 -1 1; -1 0 -1; 0 0 0]; 
    %three missing inhibitions, no circular inhibition
    W3D{13}=[0 1 0;-1 0 -1; 0 0 0]; 
    W3D{14}=[0 -1 0;1 0 -1; 0 0 0]; 
    W3D{15}=[0 -1 0;-1 0 1; 0 0 0]; 
    %three missing inhibitions, circular inhibition
    W3D{16}=[0 1 0; 0 0 -1; -1 0 0];
    W=W3D;
    break_ind=[1,2,7,11,13,16];
    if plot_on
        figure
        for i=1:16
            subplot(4,4,i)
            G = digraph(W3D{i});
            plot(G,'EdgeCData',1+0.2*G.Edges.Weight)
            if ismember(i,break_ind)
                title(strcat('Architectue:',num2str(find(break_ind==i))))
            end
        end
    end
end

if Nspecies==3 && nE==2
    W3D=cell(28,1);
    %full inhibitory network - one I-E switch
    W3D{1}=[0 1 1;-1 0 -1; -1 -1 0];
    W3D{2}=[0 1 -1;1 0 -1; -1 -1 0];
    W3D{3}=[0 -1 -1;1 0 -1; 1 -1 0];
    W3D{4}=[0 -1 1;1 0 -1; -1 -1 0];
    %one missing inhibition - one I-E switch
    W3D{5}=[0 0 1;1 0 -1; -1 -1 0]; 
    W3D{6}=[0 0 1;-1 0 1; -1 -1 0]; 
    W3D{7}=[0 0 1;-1 0 -1; 1 -1 0]; 
    W3D{8}=[0 0 1;-1 0 -1; -1 1 0]; 
    W3D{9}=[0 0 -1;1 0 1; -1 -1 0]; 
    W3D{10}=[0 0 -1;1 0 -1; 1 -1 0]; 
    W3D{11}=[0 0 -1;1 0 -1; -1 1 0]; 
    W3D{12}=[0 0 -1;-1 0 1; 1 -1 0]; 
    W3D{13}=[0 0 -1;-1 0 1; -1 1 0]; 
    W3D{14}=[0 0 -1;-1 0 -1; 1 1 0]; 
    %two missing inhibitions, circular inhibition - one I-E switch
    W3D{15}=[0 1 0; 1 0 -1; -1 0 0];
    W3D{16}=[0 1 0; -1 0 1; -1 0 0];
    W3D{17}=[0 1 0; -1 0 -1; 1 0 0];
    W3D{18}=[0 -1 0; 1 0 1; -1 0 0];
    W3D{19}=[0 -1 0; 1 0 -1; 1 0 0];
    W3D{20}=[0 -1 0; -1 0 1; 1 0 0];
    %two missing inhibitions, no circular inhibition - one I-E switch
    W3D{21}=[0 -1 1; 1 0 -1; 0 0 0]; 
    W3D{22}=[0 -1 1; -1 0 1; 0 0 0];
    W3D{23}=[0 1 -1; 1 0 -1; 0 0 0]; 
    W3D{24}=[0 -1 -1; 1 0 1; 0 0 0];
    %three missing inhibitions, no circular inhibition
    W3D{25}=[0 1 0;1 0 -1; 0 0 0]; 
    W3D{26}=[0 1 0;-1 0 1; 0 0 0]; 
    W3D{27}=[0 -1 0;1 0 1; 0 0 0]; 
    %three missing inhibitions, circular inhibition
    W3D{28}=[0 1 0; 0 0 1; -1 0 0];
    W=W3D;
    break_ind=[1,5,15,21,25,28];
    if plot_on
        figure
        for i=1:28
            subplot(7,4,i)
            G = digraph(W3D{i});
            plot(G,'EdgeCData',1+0.2*G.Edges.Weight)
            if ismember(i,break_ind)
                title(strcat('Architectue:',num2str(find(break_ind==i))))
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nspecies==3 && nE==3
    W3D=cell(23,1);
    %full inhibitory network - one I-E switch
    W3D{1}=[0 1 1;1 0 -1; -1 -1 0];
    W3D{2}=[0 1 -1;1 0 -1; 1 -1 0];
    W3D{3}=[0 1 -1;-1 0 1; 1 -1 0];
    %one missing inhibition - one I-E switch
    W3D{4}=[0 0 1;1 0 1; -1 -1 0]; 
    W3D{5}=[0 0 1;1 0 -1; 1 -1 0]; 
    W3D{6}=[0 0 1;1 0 -1; -1 1 0]; 
    W3D{7}=[0 0 1;-1 0 1; 1 -1 0]; 
    W3D{8}=[0 0 1;-1 0 1; -1 1 0]; 
    W3D{9}=[0 0 1;-1 0 -1; 1 1 0]; 
    W3D{10}=[0 0 -1;1 0 1; 1 -1 0]; 
    W3D{11}=[0 0 -1;1 0 1; -1 1 0]; 
    W3D{12}=[0 0 -1;1 0 -1; 1 1 0]; 
    W3D{13}=[0 0 -1;-1 0 1; 1 1 0]; 
    %two missing inhibitions, circular inhibition - one I-E switch
    W3D{14}=[0 1 0; 1 0 1; -1 0 0];
    W3D{15}=[0 1 0; 1 0 -1; 1 0 0];
    W3D{16}=[0 1 0; -1 0 1; 1 0 0];
    W3D{17}=[0 -1 0; 1 0 1; 1 0 0];
    %two missing inhibitions, no circular inhibition - one I-E switch
    W3D{18}=[0 1 1; 1 0 -1; 0 0 0]; 
    W3D{19}=[0 1 1; -1 0 1; 0 0 0]; 
    W3D{20}=[0 1 -1; 1 0 1; 0 0 0]; 
    W3D{21}=[0 -1 1; 1 0 1; 0 0 0]; 
    %three missing inhibitions, no circular inhibition
    W3D{22}=[0 1 0;1 0 1; 0 0 0]; 
    %three missing inhibitions, circular inhibition
    W3D{23}=[0 1 0; 0 0 1; 1 0 0];
    W=W3D;
    break_ind=[1,4,14,18,22,23];
    if plot_on
        figure
        for i=1:23
            subplot(6,4,i)
            G = digraph(W3D{i});
            plot(G,'EdgeCData',1+0.2*G.Edges.Weight)
            if ismember(i,break_ind)
                title(strcat('Architectue:',num2str(find(break_ind==i))))
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nspecies==3 && nE==4
    W3D=cell(11,1);
    %full inhibitory network - one I-E switch
    W3D{1}=[0 -1 -1;1 0 1; 1 1 0];
    W3D{2}=[0 -1 1;-1 0 1; 1 1 0];
    W3D{3}=[0 1 1;-1 0 1; -1 1 0];
    W3D{4}=[0 1 -1;-1 0 1; 1 1 0];
    %one missing inhibition - one I-E switch
    W3D{5}=[0 0 1;1 0 1; 1 -1 0]; 
    W3D{6}=[0 0 1;1 0 1; -1 1 0]; 
    W3D{7}=[0 0 1;1 0 -1; 1 1 0]; 
    W3D{8}=[0 0 1;-1 0 1; 1 1 0]; 
    W3D{9}=[0 0 -1;1 0 1; 1 1 0]; 
    %two missing inhibitions, circular inhibition - one I-E switch
    W3D{10}=[0 1 0; 1 0 1; 1 0 0];
    %two missing inhibitions, no circular inhibition - one I-E switch
    W3D{11}=[0 1 1; 1 0 1; 0 0 0]; 
    W=W3D;
    break_ind=[1,5,10,11];
    if plot_on
        figure
        for i=1:11
            subplot(3,4,i)
            G = digraph(W3D{i});
            plot(G,'EdgeCData',1+0.2*G.Edges.Weight)
            if ismember(i,break_ind)
                title(strcat('Architectue:',num2str(find(break_ind==i))))
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nspecies==3 && nE==5
    W3D=cell(2,1);
    %full inhibitory network - one I-E switch
    W3D{1}=[0 1 1;1 0 1; 1 -1 0];
    %one missing inhibition - one I-E switch
    W3D{2}=[0 0 1;1 0 1; 1 1 0]; 
    W=W3D;
    break_ind=[1,2];
    if plot_on
        figure
        for i=1:2
            subplot(1,2,i)
            G = digraph(W3D{i});
            plot(G,'EdgeCData',1+0.2*G.Edges.Weight)
            if ismember(i,break_ind)
                title(strcat('Architectue:',num2str(find(break_ind==i))))
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nspecies==3 && nE==6
    W3D=cell(1,1);
    %full inhibitory network - one I-E switch
    W3D{1}=[0 1 1;1 0 1; 1 1 0];
    W=W3D;
    break_ind=1;
    if plot_on
        figure
        G = digraph(W3D{i});
        plot(G,'EdgeCData',1+0.2*G.Edges.Weight)
        title(strcat('Architectue:',num2str(i)))
    end
end

end