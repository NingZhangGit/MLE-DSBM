function val = f_ARI(x, y, varargin)

assert(numel(x) == numel(y));
N = numel(x);
x = reshape(x,1,N);
y = reshape(y,1,N);
%% Preliminary computations and cleansing of the partitions
N = length(x);
[~, ~, x] = unique(x);
N1 = max(x);
[~, ~, y] = unique(y);
N2 = max(y);
%% Create the matching matrix
for i=1:1:N1
    for j=1:1:N2
        G1 = find(x==i);
        G2 = find(y==j);
        n(i,j) = length(intersect(G1,G2));
    end
end
%% Otherwise, calculate the adjusted rand index

ssm = 0;
sm1 = 0;
sm2 = 0;
for i=1:1:N1
    for j=1:1:N2
        ssm = ssm + nchoosek2(n(i,j),2);
    end
end
temp = sum(n,2);
for i=1:1:N1
    sm1 = sm1 + nchoosek2(temp(i),2);
end
temp = sum(n,1);
for i=1:1:N2
    sm2 = sm2 + nchoosek2(temp(i),2);
end
NN = ssm - sm1*sm2/nchoosek2(N,2);
DD = (sm1 + sm2)/2 - sm1*sm2/nchoosek2(N,2);
% Special case to handle perfect partitions
if (NN == 0 && DD==0)
    val = 1;
else
    val = NN/DD;
end

%% Special definition of n choose k
    function c = nchoosek2(a,b)
        if a>1
            c = nchoosek(a,b);
        else
            c = 0;
        end
    end
end