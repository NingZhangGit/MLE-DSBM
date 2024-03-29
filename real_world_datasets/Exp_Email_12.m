%%  Experiment on the Eu-email-core dataset (1st and 2nd)

clear all
close all
clc
addpath('../cluster_algs')
addpath('../metrics')
addpath('../functions')
%% load and preprocess data
data = load("G_email12.mat");
Ad = data.G;
y = data.y;

% randomly permute data
p = randperm(length(y));
Ad = Ad(p,p);
y = y(p);

H = 1i*(Ad-Ad');
assert(ishermitian(H))
A = Ad+Ad'; %symmetric adj (undirected adj)
assert(issymmetric(A))
N = length(Ad);
k = 2;
m = 10; %number of repetitions

% find largest comp
[bfs_comp_vertex , comp_number, length_comp] = BFS_connected_components(A);
id = find(length_comp==max(length_comp));
cl = find(bfs_comp_vertex==id);
A = A(cl,cl);
Ad = Ad(cl,cl);
y = y(cl);
N_lg = length(cl);

%% cluster with different algorithms
for iter = 1:m
    % Herm
    H = 1i*(Ad-Ad');
    y_H = f_Herm(H,k,1);
    % Herm RW
    y_Hrw = f_HermRW(Ad,k);
    % DI-SIM
    davg = sum(A,'all')/N;
    y_disim = f_DI_SIM(Ad,k,davg);
    % A+A'
    y_Asym = f_SpeCluster(A,k);
    %B-Sym
    B_sym = Ad*Ad'+Ad'*Ad;
    y_Bsym = f_SpeCluster(B_sym,2);
    % MLE-SC
    t = 1e-3;
    [y_EMsc,p_sc,q_sc,eps_sc] = f_EM_MLE_sc(Ad,k,3,t);
    % MLE-BM
    [y_EMbm,p_bm,q_bm,eps_bm] = f_EM_MLE_bm(Ad,k,3,t);
   
    % compare with ground truth
    ARI_H(iter) = f_ARI(y,y_H);
    ARI_Hrw(iter) = f_ARI(y,y_Hrw);
    ARI_disim(iter) = f_ARI(y,y_disim);
    ARI_Asym(iter) = f_ARI(y,y_Asym);
    ARI_Bsym(iter) = f_ARI(y,y_Bsym);
    ARI_EMsc(iter) = f_ARI(y,y_EMsc);
    ARI_EMbm(iter) = f_ARI(y,y_EMbm);
    fprintf(' %.3f, %.3f, %.3f, %.3f, %.3f ,%.3f , %.3f\n', ARI_H(iter),...
        ARI_Hrw(iter), ARI_disim(iter), ...
        ARI_Bsym(iter),ARI_Asym(iter),ARI_EMsc(iter),ARI_EMbm(iter))
    fprintf('Iteration: %d finished!\n',iter)

end
%%
fprintf('ARIs: H: %.3f, Hrw: %.3f, DI-SIM:  %.3f ,A-sym: %.3f, B-sym: %.3f, MLE-SC:%.3f , MLE-SDP %.3f\n',...
    mean(ARI_H), mean(ARI_Hrw), mean(ARI_disim), ...
    mean(ARI_Asym),mean(ARI_Bsym),mean(ARI_EMsc),mean(ARI_EMbm))
%% order Adj by clustering results
c_sc = [find(y_EMsc==1),find(y_EMsc==2)];
A_sc = Ad(c_sc,c_sc);

c_bm = [find(y_EMbm==1),find(y_EMbm==2)];
A_bm = Ad(c_bm,c_bm);

c_H = [find(y_H==1), find(y_H ==2)];
A_H = Ad(c_H, c_H);

c_Hrw = [find(y_H==1), find(y_H ==2)];
A_Hrw = Ad(c_Hrw, c_Hrw);

c_disim = [find(y_disim==1), find(y_disim==2)];
A_disim = Ad(c_disim, c_disim);

c_Asym = [find(y_Asym==1), find(y_Asym==2)];
A_Asym = Ad(c_Asym,c_Asym);

c_Bsym =[find(y_Bsym==1), find(y_Bsym==2)];
A_Bsym = Ad(c_Bsym,c_Bsym);

y_true = [find(y==14); find(y==4)];
A_true = Ad(y_true,y_true);

%% visualizing adj after clustering
f = figure(1);
subplot(2,4,1)
imagesc(A_true-A_true')
title('ground truth')
axis square

subplot(2,4,2)
imagesc(A_H-A_H')
title(sprintf('Herm, ARI = %.2f',ARI_H(m)))
axis square

subplot(2,4,3)
imagesc(A_Hrw-A_Hrw')
title(sprintf('HermRW, ARI = %.2f',ARI_Hrw(m)))
axis square

subplot(2,4,4)
imagesc(A_disim-A_disim')
title(sprintf('DI-SIM, ARI = %.2f',ARI_disim(m)))
axis square

subplot(2,4,5)
imagesc(A_Bsym-A_Bsym')
title(sprintf('B-Sym, ARI = %.2f',ARI_Bsym(m)))
axis square

subplot(2,4,6)
imagesc(A_Asym-A_Asym')
title(sprintf('A+A^T, ARI = %.2f',ARI_Asym(m)))
axis square

subplot(2,4,7)
imagesc(A_sc-A_sc')
title(sprintf('MLE-SC, ARI = %.2f',ARI_EMsc(m)))
sc_par = sprintf('p=%.1f %%, q =%.1f %%, eta= %.1f %% ', p_sc*100,q_sc*100,eps_sc*100);
xlabel(sc_par)
axis square

subplot(2,4,8)
imagesc(A_bm-A_bm')
title(sprintf('MLE-SDP, ARI = %.2f',ARI_EMbm(m)))
bm_par = sprintf('p=%.1f %%, q =%.1f %%, eta= %.1f %% ', p_bm*100,q_bm*100,eps_bm*100);
xlabel(bm_par)
axis square

mymap = [1 0 0;
    1 1 1;
    0 0 1];
colormap(mymap)

set(gcf,'Position',[0 0 1000 800])
