function Nodes = Sort_nodes(Nodes,Mu)
%Sort_nodes Sort the nodes according to their cluster, defined by which mu
%they are close to. Then giving the Nodes new indexes. Index ordering within
%a cluster will be consistant with the index from before. 
%  

n_cluster=length(Mu(:,1));  % Number of clusters
n_nodes=length(Nodes);      % Number of nodes

Index_matrix=zeros(n_nodes,n_cluster);  % To store the indeces of teh nodes in aech cluster
Iteration=ones(n_cluster,1);            % To keep track of the current iteration in each cluster
% Iterate through the nodes ang get each nosed index and cluster
for i=1:n_nodes
    cluster=Nodes(i).cluster;
    Index_matrix(Iteration(cluster),cluster)=Nodes(i).index;
    Iteration(cluster)=Iteration(cluster)+1;
end

Nodes2=Node;

iter=1;

% Iterate through the nodes and assign new nodes wit the old ones, by teh
% order from Index_matrix, takin one cluster at the time
for i=1:n_cluster
    j=1;
    while Index_matrix(j,i)>0
        Nodes2(iter)=Nodes(Index_matrix(j,i));
        Nodes2(iter).index=iter;
        iter=iter+1;
        j=j+1;
    end
end

Nodes=Nodes2; % replace the old nodes with the new ones. (same info but resorted) 
end

