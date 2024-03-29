% % Generate adjancency matrix (sparse) from DSBM 

% % Inputs:
% N : total number of vertices
% k: number of clusters
% n: vector whose entries are community size
% p: probability of edge within a cluster
% q: probability of edge between clusters
% eta: probability on edge from C2 to C1

% % Outputs:
% Ad: adjacency matrix (sparse)
% y: ture cluster label (1,2,3,...,k)

function [Ad,y] = f_gen_DSBM(N, k, n, p, q, eta)
%% Assign cumminity labels
y = [];
for i=1:k
    y = [y,i*ones(1,n(i))];
end

%% Generate undirected adj
% Sample N*N matrix from unif[0,1] 
U = rand(N);
% Sample binary edges
% within edge with prob p
% between cluster edge with prob q
B = sparse(N,N);
B(1:n(1),1:n(1)) = int8(U(1:n(1),1:n(1))<p);
B(n(1):n(1)+n(2),n(1):n(1)+n(2)) = int8(U(n(1):n(1)+n(2),n(1):n(1)+n(2))<p);
B(1:n(1),n(1):n(1)+n(2)) = int8(U(1:n(1),n(1):n(1)+n(2))<q);

G = triu(B,1); %graph with undirected edges
G = G+G';
assert(issymmetric(G))
%% Create direction matrix mask
% Sample N*N matrix from unif[0,1] 
U2 = rand(N);
% Sample edge direction matrix 
% +1 r->c; -1 c->r 
% within edge: prob(r->c)= 0.5
% between cluster edge: prob(r->c)= 1-eta
D = zeros(N,N);
D(1:n(1),1:n(1)) = int8(U2(1:n(1),1:n(1))<0.5)*2-1;
D(n(1):n(1)+n(2),n(1):n(1)+n(2)) = int8(U2(n(1):n(1)+n(2),n(1):n(1)+n(2))<0.5)*2-1;
D(1:n(1),n(1):n(1)+n(2)) = int8(U2(1:n(1),n(1):n(1)+n(2))<1-eta)*2-1;

%% Create directed graph
M = G.*D;
M = triu(M,1);
M = M - M';
assert(sum(M==-M','all')==N*N);

Ad = M>0;
assert(sum(Ad+Ad'==G,'all') == N*N)

end