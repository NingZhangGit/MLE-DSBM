%% This script is for MLE-driven embedding space visualization

%%
close all
clear all
clc
addpath('functions')
addpath('cluster_algs')
addpath('metrics')
%% Set DSBM parameters
N = 1000;
k = 2;
n1 = 500;
n2 = 500;
n = [n1,n2];
eta = 0.05;
p = 0.05;
q = 0.05;
% sample Ad from DSBM
[Ad,y] = f_gen_DSBM(N,k,n,p,q, eta); 

%% digraph adjacency matrix visualization
figure(1)
imagesc(Ad-Ad')
mymap = [1 0 0;
    1 1 1;
    0 0 1];
colormap(mymap)
colorbar('Ticks',[-1,0,1],...
         'TickLabels',{'-1','0','1'})
axis equal

xlim([1 1000])
ylim([1 1000])


%% Spectral embedding: top eigenvector of H_mle
% compute MLE derived Hermitian matrix
H = 1i*log((1-eta)/eta) *(Ad-Ad') + log(1/(4*eta*(1-eta)))*(Ad+Ad');
[v,val] = eigs(H,1);
figure(2)
scatter(real(v(1:n)), imag(v(1:n)),10,'red')
hold on
scatter(real(v(1+n:N)), imag(v(1+n:N)),'*','blue')
xlim([-0.05 0.05])
ylim([-0.05 0.05])
mymap = [1 0 0;
    0 0 1];
colormap(mymap)
legend('c_1','c_2')
axis equal
%% SDP embedding: top eigenvector of SDP solution 
dim = ceil(sqrt(N));
quiet = 0;
manifold = obliquecomplexfactory(dim,N);
problemBM.M = manifold;
problemBM.cost  = @(x) trace(-x*H*x');
problemBM.egrad = @(x) -2*x*H; 
problemBM.ehess = @(x,u) -2*u*H;
opts.verbosity=quiet; % Set to 0 for no output, 2 for normal output
x = trustregions(problemBM,[],opts);
X = x'*x;
[x_BM, val] = eigs(X,1);

figure(3)
scatter(real(x_BM(1:n1)), imag(x_BM(1:n1)),10,'red')
hold on
scatter(real(x_BM(1+n1:N)), imag(x_BM(1+n1:N)),'*','blue')
xlim([-0.05 0.05])
ylim([-0.05 0.05])
mymap = [1 0 0;
    0 0 1];
colormap(mymap)
legend('c_1','c_2')
axis equal


