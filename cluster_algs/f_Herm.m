% Hermitian clustering algorithm
% Input: 
%       Hermitian matrix 
%       k: number of clusters
%       rank: rank of the embedding space
% Output:
%       y_hat {1,...,k}

function y_hat = f_Herm(H,k,rank)
    
    [V_H, val] = eigs(H, rank);
    v1 = real(V_H);
    v2 = imag(V_H);
    y_hat = (kmeans([v1,v2],k))';
    assert(isreal(val))

end