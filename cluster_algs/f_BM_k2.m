%% BM method for clustering
% Inputs:
    % H: hermitian matrix i(A-A')
    % dim: dimintion 
    % quiet: 0 for quiet and 2 for display trust regions
    % k: number of clusters, k=2
% Output:
    %y_hat: {1,2}
function y_hat = f_BM_k2(H,dim,quiet,k)
N = length(H);
% define problem structure
manifold = obliquecomplexfactory(dim,N);
problemBM.M = manifold;
problemBM.cost  = @(x) trace(-x*H*x');
problemBM.egrad = @(x) -2*x*H; 
problemBM.ehess = @(x,u) -2*u*H;

% checkgradient(problemBM);
% checkhessian(problemBM);
% solve with ManOpt
opts.verbosity=quiet; % Set to 0 for no output, 2 for normal output
x = trustregions(problemBM,[],opts);
X = x'*x;
[x_BM, val] = eigs(X,1);

assert(isreal(val));
% k-means on x
y_hat = (kmeans([real(x_BM),imag(x_BM)],k))';
