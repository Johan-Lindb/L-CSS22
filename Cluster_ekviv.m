function [M_red,K_red] = Cluster_ekviv(Nodes,K,Mu)
%Cluster_ekviv Given the Nodes and the Laplacian K for the detailed system
%the reduced mass matrix and the reduced laplacian is generated
%   The detailed system is reduced to one equivalnt node per cluster, with
%   all the power iand inertia of all the nodes in the cluster. 
n_nodes=length(Nodes);      % Number of nodes

n_cluster=length(Mu(:,1));  % Number of clusters

M_red=zeros(n_cluster);     % Empty matrices to store values in

K_red=zeros(n_cluster);

% for each node add the mass to the mass matrix in the position of the
% cluster
for i=1:n_nodes
    cluster=Nodes(i).cluster;
    M_red(cluster,cluster)=M_red(cluster,cluster)+Nodes(i).m;
end

% Construct the reduced Laplacian by adding entries if the clusters are
% different for two nodes. 
for i=1:n_nodes
    for j=i:n_nodes
        if K(i,j)<0
            cluster1=Nodes(i).cluster;
            cluster2=Nodes(j).cluster;
            if cluster1 ~= cluster2
                K_red(cluster1,cluster2)=K_red(cluster1,cluster2)+K(i,j);
                K_red(cluster2,cluster1)=K_red(cluster2,cluster1)+K(i,j);
                K_red(cluster1,cluster1)=K_red(cluster1,cluster1)-K(i,j);
                K_red(cluster2,cluster2)=K_red(cluster2,cluster2)-K(i,j);
            end
        end
    end
end

% Check if a cluster did not have any nodes. Then remove this cluster to
% avoid singularities later. This happend if the number of clusters is  
% close to the number of nodes. 
i=1;
while i<n_cluster+1
    if K_red(i,i)==0
        K_red(i,:)=[];
        K_red(:,i)=[];
        n_cluster=n_cluster-1;
        M_red(i,:)=[];
        M_red(:,i)=[];
    else
        i=i+1;
    end
end

end

