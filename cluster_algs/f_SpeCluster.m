%% Spectral clustering algorithm
% Input:
%   M: inpur matrix 
%   k: number of clusters
% Outputs:
    % y_hat: cluster label vector

%%
function y_hat=f_SpeCluster(M,k)
    % compute top-k eigenvectors
    [U,val] = eigs(M,2);
    assert(isreal(val));
    % apply k-means alg
    y_hat = (kmeans(U(:,2),k))';
end