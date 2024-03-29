%% DSBM visualize adjacency relation before & after clustering 

close all
clear all
clc
addpath('../cluster_algs')
addpath('../metrics')
addpath('../functions')

%% (1) Generate DSBM
% DSBM parameters
N = 200;
k = 2;
n1 = 100;
n2 = 100;
n = [n1,n2];

%% p=q=1%
ARI_all_H_1 = [];
ARI_all_Hrw_1 = [];
ARI_all_disim_1 = [];
ARI_all_Bsym_1 = [];
ARI_all_MLEsc_1 = [];
ARI_all_MLEbm_1 = [];

p = 0.05;
q = 0.05;
eps = 0.1;


% Sample Ad from DSBM
[Ad,y] = f_gen_DSBM(N,k,n,p,q,eps); %Hermitian adj & ture label
A = Ad+Ad';
% find largest comp
[bfs_comp_vertex , comp_number, length_comp] = BFS_connected_components(A);
id = find(length_comp==max(length_comp));
cl = find(bfs_comp_vertex==id);
A = A(cl,cl);
Ad = Ad(cl,cl); %directed adj
y = y(cl);
N_lg = length(cl);

p = randperm(length(y));
Ad = Ad(p,p);
y = y(p);
%% (2) Spectral Clustering
% (2.1) Hermitian clustering
H = 1i*(Ad-Ad');
y_hat_H = f_Herm(H,k,1);
y_hat_Hrw = f_HermRW(Ad,k);

% (2.2) DI-SIM
y_hat_disim = f_DI_SIM(Ad,k,1);

% (2.3) symmetrization
B_sym = Ad*Ad'+Ad'*Ad;
y_hat_Bsym = f_SpeCluster(B_sym,2);
Asym = Ad'+Ad;
y_hat_Asym = f_SpeCluster(Asym,2);

% (2.5) iterative MLE
y_hat_MLEsc = f_IT_MLE_sc(Ad,k,2,1e-2);
y_hat_MLEbm = f_IT_MLE_bm(Ad,k,2,1e-2);

%% comparision
ARI_H = f_ARI(y,y_hat_H);
ARI_Hrw = f_ARI(y,y_hat_Hrw);
ARI_disim = f_ARI(y,y_hat_disim);
ARI_Bsym = f_ARI(y,y_hat_Bsym);
ARI_Asym = f_ARI(y,y_hat_Asym);
ARI_MLEsc = f_ARI(y,y_hat_MLEsc);
ARI_MLEbm = f_ARI(y,y_hat_MLEbm);

%% order index according to clustering outcomes
c_sc = [find(y_hat_MLEsc==1),find(y_hat_MLEsc==2)];
A_sc = Ad(c_sc,c_sc);

c_bm = [find(y_hat_MLEbm==1),find(y_hat_MLEbm==2)];
A_bm = Ad(c_bm,c_bm);

c_H = [find(y_hat_H==1), find(y_hat_H ==2)];
A_H = Ad(c_H, c_H);

c_Hrw = [find(y_hat_Hrw==1), find(y_hat_Hrw ==2)];
A_Hrw = Ad(c_Hrw, c_Hrw);

c_disim = [find(y_hat_disim==1), find(y_hat_disim==2)];
A_disim = Ad(c_disim, c_disim);

c_Asym = [find(y_hat_Asym==1), find(y_hat_Asym==2)];
A_Asym = Ad(c_Asym,c_Asym);

c_Bsym =[find(y_hat_Bsym==1), find(y_hat_Bsym==2)];
A_Bsym = Ad(c_Bsym,c_Bsym);

%% visualize adjacency relation A-A^T
f = figure(1);
subplot(2,4,1)
imagesc(Ad-Ad')
title('input graph')
axis square

subplot(2,4,2)
imagesc(A_H-A_H')
title(sprintf('Herm, ARI = %.2f',ARI_H))
axis square

subplot(2,4,3)
imagesc(A_Hrw-A_Hrw')
title(sprintf('HermRW, ARI = %.2f',ARI_Hrw))
axis square

subplot(2,4,4)
imagesc(A_disim-A_disim')
title(sprintf('DI-SIM, ARI = %.2f',ARI_disim))
axis square

subplot(2,4,5)
imagesc(A_Bsym-A_Bsym')
title(sprintf('B-Sym, ARI = %.2f',ARI_Bsym))
axis square

subplot(2,4,6)
imagesc(A_Asym-A_Asym')
title(sprintf('A+A^T, ARI = %.2f',ARI_Asym))
axis square

subplot(2,4,7)
imagesc(A_sc-A_sc')
title(sprintf('MLE-SC, ARI = %.2f',ARI_MLEsc))
axis square

subplot(2,4,8)
imagesc(A_bm-A_bm')
title(sprintf('MLE-SDP, ARI = %.2f',ARI_MLEbm))
axis square

mymap = [1 0 0;
    1 1 1;
    0 0 1];
colormap(mymap)

set(gcf,'Position',[0 0 900 400])


