%% DI_SIM Cluster algorithm (optional:reg)
% Input:
    % Ad: graph adjacency matrix
    % k: number of clusters, in this study k=2
    % s: regularization parameter (0=no reg)
% Output: 
    % y_hat:cluster label vector 

%%
function y_hat = f_DI_SIM(Ad,k,reg)
O = sum(Ad,2);  %sum over rows
P = sum(Ad,1);  %sum over columns
if reg==0
    O_norm = O;
    P_norm = P;
    O_norm(~O_norm) = 1;
    P_norm(~P_norm) = 1;
    A_norm = diag(sqrt(1./O_norm))*Ad*diag(sqrt(1./P_norm));
else
    t = sum(O)/length(Ad);
    O_norm = O+t;
    P_norm = P+t;
    O_norm(~O_norm) = 1;
    P_norm(~P_norm) = 1;
    A_norm = diag(sqrt(1./O_norm))*Ad*diag(sqrt(1./P_norm));
end

% compute leading eigenvector of AA' and A'A
[v_disim_L, val] = eigs(A_norm*A_norm', ceil(k/2));
assert(isreal(val))
[v_disim_R, val] = eigs(A_norm'*A_norm, ceil(k/2));
assert(isreal(val))

% row nomalization
V_disim = [v_disim_L,v_disim_R];
V_disim_rownorm = sqrt(sum(V_disim.^2,2));
V_disim_norm = V_disim./V_disim_rownorm;
% k-means cluster
y_hat = (kmeans(V_disim_norm, k))';
end