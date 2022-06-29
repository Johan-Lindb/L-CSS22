function [Nodes_new,M_new,K_new, P_new] = Restructure_Clusters(Nodes, M, K, P, Mu)
%Restructure_Clusters Restructureing the clusters so that the renewable
%sources shows up in the beginneing of the matrices this will help with the
%plotting to be more comprehensible. 

n_cluster=length(Mu(:,1));  % Number of clusters
n_nodes=length(Nodes);      % Number of nodes
% info will be 1: index, 2: start node, 3: end node, 4: boolean is renewable
Cluster_info=zeros(n_cluster,4);

j=0;

% Iterate through the nodes and fill in the start node and end node of each
% cluser and if it is renewable or not
for i=1:n_cluster
    Cluster_info(i,1)=i;
    prev_cluster_end=j;
    j=j+1;
    if j>n_nodes
        break
    end
    while Nodes(j).cluster==i
        j=j+1;
        if j>n_nodes
            j=n_nodes+1;
            break
        end
    end
    j=j-1;
    Cluster_info(i,2)=prev_cluster_end+1;
    Cluster_info(i,3)=j;
    if contains(Nodes(j).type,'Wind')
        Cluster_info(i,4)=1;
    end
end


K_new=zeros(n_nodes);       % New Laplacian, inertia matrix and power matrix
M_new=zeros(n_nodes);
P_new=zeros(1,n_nodes);
K_temp=K;
Nodes_new=Node;
index_change=zeros(1,n_nodes);
cluster_change=zeros(1,n_nodes);

endpos=0;
cluster=1;
% Start with the clusters with renewables and add them to the new matrices.
% To avoid mistakes the values are then removed from K_temp
for i=1:n_cluster 
    if Cluster_info(i,4)
        start_node=Cluster_info(i,2);
        end_node=Cluster_info(i,3);
        if end_node>=start_node
        M_new(endpos+1:endpos+end_node-start_node+1,endpos+1:endpos+end_node-start_node+1)=M(start_node:end_node,start_node:end_node);
        K_new(endpos+1:endpos+end_node-start_node+1,endpos+1:endpos+end_node-start_node+1)=K(start_node:end_node,start_node:end_node);
        K_temp(start_node:end_node,start_node:end_node)=0;
        P_new(endpos+1:endpos+end_node-start_node+1)=P(start_node:end_node);
        index_change(endpos+1:endpos+end_node-start_node+1)=start_node:1:end_node;
        cluster_change(endpos+1:endpos+end_node-start_node+1)=cluster;
        endpos=endpos+end_node-start_node+1;
        cluster=cluster+1;
        end
    end
end

% Then add the remaining clusters after the reneables (in the order they
% came in before)
for i=1:n_cluster
    if Cluster_info(i,4)<1
        start_node=Cluster_info(i,2);
        end_node=Cluster_info(i,3);
        if start_node>0
        if end_node>=start_node
        M_new(endpos+1:endpos+end_node-start_node+1,endpos+1:endpos+end_node-start_node+1)=M(start_node:end_node,start_node:end_node);
        K_new(endpos+1:endpos+end_node-start_node+1,endpos+1:endpos+end_node-start_node+1)=K(start_node:end_node,start_node:end_node);
        K_temp(start_node:end_node,start_node:end_node)=0;
        P_new(endpos+1:endpos+end_node-start_node+1)=P(start_node:end_node);
        index_change(endpos+1:endpos+end_node-start_node+1)=start_node:1:end_node;
        cluster_change(endpos+1:endpos+end_node-start_node+1)=cluster;
        endpos=endpos+end_node-start_node+1;
        cluster=cluster+1;
        end
        end
    end
end

%fix transmission nodes by iteraing through the remaining K_temp and place
%the values in the new laplacian
for i=1:n_nodes
    for j=i:n_nodes
        if K_temp(i,j)<0
            ind1=find(index_change==i);
            ind2=find(index_change==j);
            K_new(ind1,ind2)=K_temp(i,j);
            K_new(ind2,ind1)=K_temp(i,j);
        end
    end  
end

% Update the index, cluster and the connections
for i=1:n_nodes
    Nodes_new(i)=Nodes(index_change(i));
    Nodes_new(i).index=i;
    Nodes_new(i).cluster=cluster_change(i);
    con=[];
    for j=1:n_nodes
        if K_new(i,j)<0
            con=[con, j];
        end
    end
    Nodes_new(i).connection=con;
end


end

