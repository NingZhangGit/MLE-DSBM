% Spectral cluster on given Hermitian matrix
% Input: 
    % Ad: graph adjacency matrix
    % k: number of clusters, in this study k=2
% Output:
%       y_hat: cluster label vector 
%%
function y_hat = f_herm_sc_k2(H,k)

    % Compute leading eigenvector
    [V_H, val] = eigs(H, 1);
    assert(isreal(val))
    % k-means cluster
    v1 = real(V_H);
    v2 = imag(V_H);
    y_hat = (kmeans([v1,v2],k))';
    
end