
%% Iterative algorithm (Algorithm 4) + spectral clustering (Algorithm1) 
% Inputs:
    % Ad: graph adjacency matrix
    % k: number of clusters, in this study k=2
    % opt: choose initializations for iterative algorithm
        % option 1:  i(A-A^T) + (A+A^T)
        % option 2: Hermitian clustering i(A-A^T)
    % t: threshold
% Outputs:
    % y_bm: cluster label vector
    % p_: estimated DSBM parameter p 
    % q_: estimated DSBM parameter q
    % eta_: estimated DSBM parameter eta 
%%
function [y_sc,p_,q_,eta_] = f_IT_MLE_sc(Ad,k,opt,t)

N_lg = length(Ad);
A = Ad+Ad';
count = 0;
max_count = 50; % maximum num of iterations
% Initialize Hermitian matrix
if opt==1
    w1 = 1;
    w2 = 1;
    w3 = 0;
else
    w1 = 1;
    w2 = 0;
    w3 = 0;
end
error = 100;
p_val =[];
q_val = [];
eta_val =[];
while error>t
    count = count+1;
    % Clustering
    HH = w1*1i*(Ad-Ad') + w2*(Ad+Ad') + w3*(ones(N_lg,N_lg)-eye(N_lg)-Ad-Ad');
    y_sc = f_Herm(HH,k,1);

    % Update p q eta
    y_1 = zeros(N_lg,1);
    y_2 = zeros(N_lg,1);
    y_1(y_sc==1)=1;
    y_2(y_sc==2)=1;
    len1 = sum(y_1);
    len2 = sum(y_2);

    
    if len1==0 || len2==0 || len1==1 || len2==1
       count = count-1;
       continue
    end
    % Count edges
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
    eta_val(count) = min(size12/(size12+size21),size21/(size12+size21));

    if p_val(count)==0 
        p_val(count) = 1e-10;
    end
    if q_val(count)==0 
        q_val(count) = 1e-5;
    end
    if eta_val(count)==0 
        eta_val(count)= 1e-5;
    end

    % Update Hermitian matrix using MLE derived formula
    w1_new = -log((1-eta_val(count))/eta_val(count));
    w2_new = -log(4*eta_val(count)*(1-eta_val(count))) + 2*log(p_val(count)/q_val(count));
    w3_new = 2*log((1-p_val(count))/(1-q_val(count)));
    error = (abs(w1-w1_new) + abs(w2-w2_new) + abs(w3-w3_new))./(abs(w1) + abs(w2) + abs(w3));
    % fprintf('Iteration %d: updates = %.2f \n',count,error)
    w1 = w1_new;
    w2 = w2_new;
    w3 = w3_new;
    
    if count>=max_count
        fprintf('Reach maximum number of iterations! \n')
        break;
    end
end
% Outpute converged p,q,eta values
p_ = p_val(count);
q_ = q_val(count);
eta_ = eta_val(count);
end