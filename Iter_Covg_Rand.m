% The script visualize how the iterative algorithm (Algorithm 4 in our
% paper) updates the DSBM parameters p, q and eta.

%%
close all
clear all
clc

addpath('../cluster_algs')
addpath('../metrics')

%% (1) Generate DSBM
% DSBM parameters
N = 2000;
k = 2;
n1 = 1000;  %|c1|
n2 = 1000;  %|c2|
n = [n1,n2];
p = 0.02;
q = 0.01;
eps = 0.1;

% sample directed graph adj Ad and ture label vector y
[Ad,y] = f_gen_DSBM(N,k,n,p,q,eps); 

% find largest comp
A = Ad+Ad';
[bfs_comp_vertex , comp_number, length_comp] = BFS_connected_components(A);
id = find(length_comp==max(length_comp));
cl = find(bfs_comp_vertex==id);
A = A(cl,cl);
Ad = Ad(cl,cl);
y = y(cl);
N_lg = length(cl);


%% (2) MLE hermitian
if eps==0
    H_mle = -1i*(Ad-Ad') + (Ad+Ad');
elseif eps ==0.5
    H_mle = 2*log(p/q)*(Ad +Ad') +2*log((1-p)/(1-q))*(ones(N_lg,N_lg)-eye(N_lg)-Ad-Ad');
else
    H_mle = 1i*log((1-eps)/eps)*(Ad - Ad') ...
        -log(4*eps*(1-eps))*(Ad +Ad')...
        +2*log(p/q)*(Ad +Ad') ...
        +2*log((1-p)/(1-q))*(ones(N_lg,N_lg)-eye(N_lg)-Ad-Ad');
end
assert(ishermitian(H_mle))

% Spectral Clustering with true p, q, eps
y_MLEsc = f_herm_sc_k2(H_mle,2);
% BM method for SDP clustering with true p, q, eps
dim = ceil(sqrt(N_lg));
y_MLEbm = f_BM_k2(H_mle,dim,0,2);
%% (3) iterative clustering with no side info
count = 1;
diff_sc =[];
diff_bm = [];
max_count = 100;

% randomly initialize p,q and eta for MLE-SC
p0 = rand()*0.3;
q0 = rand()*0.3;
eta0 = rand()*0.3;

eps_val(count) = eta0;
p_val(count) = p0;
q_val(count) = q0;

% randomly initialize p,q and eta for MLE-SDP
eps_val_bm(count) = eta0;
p_val_bm(count) = p0;
q_val_bm(count) = q0;

w1 = -log((1-eps_val(count))/eps_val(count));
w2 = -log(4*eps_val(count)*(1-eps_val(count))) + 2*log(p_val(count)/q_val(count));
w3 =  2*log((1-p_val(count))/(1-q_val(count)));

w1_bm = w1;
w2_bm = w2;
w3_bm = w3;

error = 1;
Niter = 5;

while count < Niter
    count = count+1;
    % compute Herm-MLE matrix
    H_sc = w1*1i*(Ad-Ad') + w2*(Ad+Ad') + w3*(ones(N_lg,N_lg)-eye(N_lg)-Ad-Ad');
    H_bm = w1_bm*1i*(Ad-Ad') + w2_bm*(Ad+Ad') + w3_bm*(ones(N_lg,N_lg)-eye(N_lg)-Ad-Ad');
    % clustering
    y_Hsc = f_herm_sc_k2(H_sc,k);
    y_Hbm = f_BM_k2(H_bm,dim,0,2);
    % comparing
    diff_sc(count) = f_ARI(y_MLEsc, y_Hsc);
    diff_bm(count) = f_ARI(y_MLEbm, y_Hbm);

    %updating p q eps for spectral alg
    y_1 = zeros(N_lg,1);
    y_2 = zeros(N_lg,1);
    y_1(y_Hsc==1)=1;
    y_2(y_Hsc==2)=1;
    len1 = sum(y_1);
    len2 = sum(y_2);
    %count edges
    size1 = 0.5*y_1'*A*y_1; %directed edges within c1
    size2 = 0.5*y_2'*A*y_2; %directed edges within c2
    size12 = y_1'*Ad*y_2;   %edges from c1 to c2
    size21 = y_2'*Ad*y_1;   %edges from c2 to c1
    sizeG = ones(N_lg,1)'*Ad*ones(N_lg,1); % #total edges
    assert(sizeG==size1+size2+size12+size21)
    assert((N_lg-1)*N_lg == len1*(len1-1) + len2*(len2-1) + 2*len1*len2);

    % p = (|E|-TF)/({|c1|/choose2} + {|c2|/choose2})
    p_val(count) = 2*(size1+size2)/(len1*(len1-1) + len2*(len2-1));
    % q = TF/(|c1|*|c2|)
    q_val(count) = (size12+size21)/(len1*len2);
    % eta = min {|c1->c2|/TF, |c2->c1|/TF}
    eps_val(count) = min(size12/(size12+size21),size21/(size12+size21));

    % update weighting parameters 
    w1_new = log((1-eps_val(count))/eps_val(count));
    w2_new = -log(4*eps_val(count)*(1-eps_val(count))) + 2*log(p_val(count)/q_val(count));
    w3_new = 2*log((1-p_val(count))/(1-q_val(count)));
    error = (abs(w1-w1_new) + abs(w2-w2_new) + abs(w3-w3_new));
    fprintf('Iteration %d: update = %.1f%% \n',count,error*100)
    w1 = w1_new;
    w2 = w2_new;
    w3 = w3_new;

    % update parameters for SDP alg
    y_1_bm = zeros(N_lg,1);
    y_2_bm = zeros(N_lg,1);
    y_1_bm(y_Hbm==1)=1;
    y_2_bm(y_Hbm==2)=1;
    len1_bm = sum(y_1_bm);
    len2_bm = sum(y_2_bm);
    %count edges
    size1_bm = 0.5*y_1_bm'*A*y_1_bm; %directed edges within c1
    size2_bm = 0.5*y_2_bm'*A*y_2_bm; %directed edges within c2
    size12_bm = y_1_bm'*Ad*y_2_bm;   %edges from c1 to c2
    size21_bm = y_2_bm'*Ad*y_1_bm;   %edges from c2 to c1
    sizeG_bm = ones(N_lg,1)'*Ad*ones(N_lg,1); % #total edges
    assert(sizeG_bm==size1_bm+size2_bm+size12_bm+size21_bm)
    assert((N_lg-1)*N_lg == len1_bm*(len1_bm-1) + len2_bm*(len2_bm-1) + 2*len1_bm*len2_bm);

    % p = (|E|-TF)/({|c1|/choose2} + {|c2|/choose2})
    p_val_bm(count) = 2*(size1+size2)/(len1*(len1-1) + len2*(len2-1));
    % q = TF/(|c1|*|c2|)
    q_val_bm(count) = (size12+size21)/(len1*len2);
    % eta = min {|c1->c2|/TF, |c2->c1|/TF}
    eps_val_bm(count) = min(size12/(size12+size21),size21/(size12+size21));
    
    % update weighting parameters 
    w1_new_bm = log((1-eps_val_bm(count))/eps_val_bm(count));
    w2_new_bm = -log(4*eps_val_bm(count)*(1-eps_val_bm(count))) + 2*log(p_val_bm(count)/q_val_bm(count));
    w3_new_bm = 2*log((1-p_val_bm(count))/(1-q_val_bm(count)));
    error_bm = (abs(w1_bm-w1_new_bm) + abs(w2_bm-w2_new_bm) + abs(w3_bm-w3_new_bm));
    fprintf('Iteration %d: update = %.1f%% \n',count,error_bm*100)
    w1_bm = w1_new_bm;
    w2_bm = w2_new_bm;
    w3_bm = w3_new_bm;
end

%% visualize the results
f = figure;
subplot(1,3,1)
yline(p*100,'-.','ture $p$','Interpreter','latex')
hold on
p_it(1)= plot(linspace(1,count,count),p_val*100,'-o');
p_it(2)= plot(linspace(1,count,count),p_val_bm*100,'-*');
xlim([1,count])
xlabel('iteration','Interpreter','latex')
ytickformat('%.1f \%')
pmin = min([min(p_val*100), min(p_val_bm*100),p*100])*0.95;
pmax = max([max(p_val*100), max(p_val_bm*100),p*100])*1.05;
ylim([pmin,pmax])
title('$p$','Interpreter','latex','position', [1,pmax])
hold off
axis square

subplot(1,3,2)
hold on
yline(q*100,'-.','ture $q$','Interpreter','latex')
plot(linspace(1,count,count),q_val*100,'-o')
plot(linspace(1,count,count),q_val_bm*100,'-*')
xlabel('iteration','Interpreter','latex')
ytickformat('%.1f \%')
qmin = min([min(q_val*100), min(q_val_bm*100),q*100])*0.95;
qmax = max([max(q_val*100), max(q_val_bm*100),q*100])*1.05;
ylim([qmin,qmax])
xlim([1,count])
title('$q$','Interpreter','latex','Position',[1,qmax])
axis square
hold off


subplot(1,3,3)
hold on
plt(1) = yline(eps,'-.','ture $\eta$','Interpreter','latex');
plt(2) = plot(linspace(1,count,count),eps_val,'-o');
plt(3) = plot(linspace(1,count,count),eps_val_bm,'-*');
xlabel('iteration','Interpreter','latex')
ytickformat('%.2f')
epsmin = min([min(eps_val), min(eps_val_bm),q])*0.95;
epsmax = max([max(eps_val), max(eps_val_bm),q])*1.05;
ylim([epsmin,epsmax])
xlim([1,count])
title('$\eta$','Position',[1,epsmax*1],'Interpreter','latex')
legend(plt,'ture parameters','SC', 'SDP')
axis square
hold off

set(gcf,'Position',[0 0 600 300])

%% Compare ARIs between:
% 1. cluster results from MLE using ture DSBM parameters 
% 2. cluster results from iterative learning

figure()
hold on
cplt(1) = plot(linspace(1,count,count),diff_sc,'-o');
cplt(2) = plot(linspace(1,count,count),diff_bm,'-*');
hold off
legend(cplt,'SC','SDP')
title('ARI')
xlabel('# iterations')

