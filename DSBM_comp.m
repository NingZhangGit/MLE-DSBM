%% Comparing with existing algorithm on DSBM synthetic dataset


%%
close all
clear all
clc
addpath('../cluster_algs')
addpath('../functions')
addpath('../metrics')
addpath('..')

%% (1) Generate DSBM
% set DSBM parameters
N = 2000;
k = 2;
n1 = 1000;
n2 = 1000;
n = [n1,n2];
m = 1; %number of repeatations

p = 0.02;
q = 0.01;
eps = linspace(0,0.5,6);

ARI_H = [];
ARI_Hrw = [];
ARI_disim = [];
ARI_Bsym = [];
ARI_Asym = [];
ARI_MLEsc = [];
ARI_MLEbm = [];
%%
for iter = 1:m
    for ii=1:length(eps)
        % Sample Ad from DSBM
        [Ad,y] = f_gen_DSBM(N,k,n,p,q,eps(ii)); %Hermitian adj & ture label
        A = Ad+Ad';
        % find largest comp
        [bfs_comp_vertex , comp_number, length_comp] = BFS_connected_components(A);
        id = find(length_comp==max(length_comp));
        cl = find(bfs_comp_vertex==id);
        A = A(cl,cl);
        Ad = Ad(cl,cl); %directed adj
        y = y(cl);
        N_lg = length(cl);

        %  Clustering Algorithms
        % (2.1) Hermitian clustering
        H = 1i*(Ad-Ad');
        y_hat_H = f_Herm(H,k,1);
        % (2.2) Hermitian random walk clustering
        y_hat_Hrw = f_HermRW(Ad,k);

        % (2.3) DI-SIM
        y_hat_disim = f_DI_SIM(Ad,k,1);

        % (2.4) naive symmetrization
        Asym = Ad'+Ad;
        y_hat_Asym = f_SpeCluster(Asym,2);

        % (2.5) B symmetrization
        B_sym = Ad*Ad'+Ad'*Ad;
        y_hat_Bsym = f_SpeCluster(B_sym,2);

        % (2.6) iterative MLE-SC
        y_hat_MLEsc = f_IT_MLE_sc(Ad,k,1,1e-2);
        % (2.7) iterative MLE-BM
        y_hat_MLEbm = f_IT_MLE_bm(Ad,k,1,1e-2);

        % comparision
        ARI_H(ii,iter) = f_ARI(y,y_hat_H);
        ARI_Hrw(ii,iter) = f_ARI(y,y_hat_Hrw);
        ARI_disim(ii,iter) = f_ARI(y,y_hat_disim);
        ARI_Bsym(ii,iter) = f_ARI(y,y_hat_Bsym);
        ARI_Asym(ii,iter) = f_ARI(y,y_hat_Asym);
        ARI_MLEsc(ii,iter) = f_ARI(y,y_hat_MLEsc);
        ARI_MLEbm(ii,iter) = f_ARI(y,y_hat_MLEbm);
    end
    fprintf('Finish %d iterations. \n',iter)
end
%% Visualization
f = figure(1);
hold on
p_ari(1) = plot(eps,mean(ARI_H,2));
p_ari(2) = plot(eps,mean(ARI_Hrw,2),'-.','Color','g');
p_ari(3) = plot(eps,mean(ARI_disim,2),'-^');
p_ari(4) = plot(eps,mean(ARI_Bsym,2),':');
p_ari(5) = plot(eps,mean(ARI_Asym,2),':square');
p_ari(6) = plot(eps,mean(ARI_MLEsc,2),'-o','Color','b');
p_ari(7) = plot(eps,mean(ARI_MLEbm,2),':*','Color','r');

leg = legend(p_ari,'Herm','HermRW','DI-SIM','B-Sym', 'A+A^T', 'MLE-SC', 'MLE-SDP');
set(leg,'Box','off')
hold off
xlabel('$\eta$','Interpreter','latex')
ylabel('ARI','Interpreter','latex')
title('$p=1.0\%, q= 0.5\%$','Interpreter','latex')
xlim([0,0.5])
axis square
box on
set(gcf,'Position',[0 0 220 220])


