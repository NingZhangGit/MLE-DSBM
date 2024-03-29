%% BFS search for finding largest component in a graph
% Inputs:
% G = Adjacency matrix of an undirected graph
% Neighbors{i} = the neighbors of node i. If not given, then it's computed

% Outputs:
% comp_number = number of components  %% Integer
% bfs_comp_vertex(i) = the component to which vertex i belongs. LENGTH=n
% length of each component

%%
function [bfs_comp_vertex , comp_number, length_comp] = BFS_connected_components(G, Neighbors)

n= length(G);
nargin;

if nargin ==1
    %tic;
    Neighbors = cell(n,1);
    for i=1:n
        t = find(G(:,i) == 1);
        Neighbors{i} = t;
        % disp([int2str(i) ' : ' int2str(Neighbors{i}) ', with larger :' int2str(Neighbors_larger{i})]);
    end
    %toc;
    %tmp = toc;
    % disp(['Found the neighborhood lists in ' num2str(tmp) ' seconds.']);
end
    
    
    
%tic;
    
bfs_comp_vertex = zeros(n,1); % = k if this vertex belongs to component k
added_to_queue  = zeros(n,1);
comp_number = 0;
length_comp = zeros(10,1);

go_on = 1;
while go_on == 1
    first_zero = min(find(bfs_comp_vertex==0));
    if isempty(first_zero)
        go_on=0;
    else
        comp_number = comp_number + 1;
%        disp('==========================================================');
%        disp(['Searching for component ' int2str(comp_number)]);
        pointer = first_zero; % the current vertex whose neighbors will be marked
        added_to_queue(pointer) = 1;
        bfs_comp_vertex(pointer) =comp_number;

        clear QUEUE
        QUEUE(1)=pointer;
        head = 0; tail = 1;

        while head<tail
%           disp('*******************************');
            head = head+1;
            tail;
            QUEUE(head : tail);
            pointer = QUEUE(head);
            %%%%%% nbrs_pointer = find( G(pointer,:)==1);
            nbrs_pointer = Neighbors{pointer};
            bfs_comp_vertex(nbrs_pointer)= comp_number;  
            added_or_no_prev = added_to_queue(nbrs_pointer);
            index_vert_to_add = find(added_or_no_prev==0);

            % Add the neighbors to the end of the queue,  if they have never
            % been added previously
            nbrs_to_add = nbrs_pointer(index_vert_to_add);
            if  length(nbrs_to_add) >=1
                QUEUE(tail+1 : tail+length(nbrs_to_add)) = nbrs_to_add;
                tail = tail + length(nbrs_to_add);
                added_to_queue(nbrs_to_add) = 1;
            else
%                disp('nothing to add');
            end
            %input('go on');
        end
        length_comp(comp_number) = head; %% = tail
    end
end
length_comp = length_comp(1:comp_number);
%tmp = toc;
% disp(['Found the components in ' num2str(tmp) ' seconds.']);

end
