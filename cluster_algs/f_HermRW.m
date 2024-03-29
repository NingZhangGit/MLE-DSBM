%% Hermitian cluster on D^{-1}*i(A-Ad)
% For numerical stability: we first compute eigenvectors of Hsym and then 
% normailze and get eigenvectors of Hrw

% Input: 
    % Ad: graph adjacency matrix
    % k: number of clusters, in this study k=2
% Outputs:
    % y_hat: cluster label vector
%%
function y_hat = f_HermRW(Ad,k)

    A = Ad+Ad';
    d = sum(A,1); % degree of the undirected adj
    d(~d) = 1;
    H = 1i*(Ad-Ad');
    H_sym = diag(1./sqrt(d))*H*diag(1./sqrt(d));     % Symmetrically normilized Hermitian adj

    % compute leading eigenvector Hsym
    [V_sym,val] = eigs(H_sym,1);
    assert(isreal(val))
    % compute leading eigenvector Hrw
    v1_rw = diag(1./sqrt(d))*real(V_sym);
    v2_rw = diag(1./sqrt(d))*imag(V_sym);
    % k-means cluster
    y_hat = (kmeans([v1_rw,v2_rw],k))';
end