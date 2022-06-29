function [Q, Nodes, resorted] = Connect_lines(Nodes,Mu,k_dist, k_trans, print)
%Connect_lines Connect the lines between the nodes
%   Nodes assosiated with one cluster are first connected through a
%   minimum spanning three. Then the most central nodes of each cluster is
%   picked to be connected to the central nodes of neighboring clusters.
%   This is done to ensure n-1 security between clusters to represent a
%   transmission grid. 
%   The output Q is the weithted Laplacian. The Nodes are returned with the
%   new informationa about ho they are connected. If the nodes were
%   resorted that is returned as a boolean. 
%   The inputs are the Nodes to get information and to add the information
%   about the connections. Mu is used to determine the number of clusters.
%   k_dist and k_trans are the line parameters for the distribution and the
%   transmossion respectively. The name distributon was the working name, but it
%   corresponds more to a sub-transmission grid. If print=1 a map will be
%   drawn with the connections, if print=0, no map will be drawn. 

n_cluster=length(Mu(:,1));      % Number of clusters
n_nodes=length(Nodes);          % Number of nodes
resorted=0;                     % Boolean number to return to tell if the nodes were resorted or not
resort=0;                       % Intermediate boolean variable in the sorting check

% Iterate through the nodes. If one of the nodes have not been assigned a
% cluster assign it the one there it is closest to the midpoint of that
% cluster. Then resort the nodes to get all within a cluster with indexes
% after each other
for i=1:n_nodes
    if Nodes(i).cluster==0
        pos=Nodes(i).pos;
        dist=NaN(1,n_cluster);
        for j=1:n_cluster
            dist(j)=sqrt((pos(1)-Mu(j,2))^2+(pos(2)-Mu(j,1))^2);
        end
        [~, clu]=min(dist);
        Nodes(i).cluster=clu;
        resort=1;
    end
end
if resort
    Nodes=Sort_nodes(Nodes,Mu);
    disp('Nodes Resorted')
    resorted=1;
end


Index_matrix=zeros(n_nodes,n_cluster); % Matrix to store the indeces of each cluster. The indeces of cluster i is strores in column i. 
Iteration=ones(n_cluster,1);           % Vector ro keep track of the current position in each column of Index_matrix
for i=1:n_nodes
    Nodes(i).connection=[];
    cluster=Nodes(i).cluster;
    Index_matrix(Iteration(cluster),cluster)=Nodes(i).index;
    Iteration(cluster)=Iteration(cluster)+1;
end

% The laplacian to start to fill out
Q=[];

% Vector to store the most eigenvalue central nodes of each cluster
Central_nodes=zeros(n_cluster,1);

for i=1:n_cluster % Creating a minimum spanning tree in each cluster
    length_cluster=nnz(Index_matrix(:,i));
    Length_info=[];
    % create a list of all distances between all nodes of a cluster with
    % the indeces of the nodes in the second and third column. Then sort 
    % the list
    for j=1:length_cluster
        pos1=Nodes(Index_matrix(j,i)).pos;
        for k=j+1:length_cluster
            pos2=Nodes(Index_matrix(k,i)).pos;
            dist=sqrt((pos1(1)-pos2(1))^2+(pos1(2)-pos2(2))^2);
            Length_info=[Length_info; [dist Index_matrix(j,i) Index_matrix(k,i)]];
        end
    end
    
    Length_info=sortrows(Length_info);
    
    % q is the local laplacian for the cluster
    q=zeros(length_cluster);
    iter=1;
    j=1;
    if length_cluster==1 % If a cluster has one node that node is the most central the 
        q=0;
        Central_nodes(i)=Index_matrix(1,i);
       
    elseif length_cluster==0 % If a cluster empty the local laplacian is empty
        q=[];
        
    else 
        % If a cluster has more than 1 node a minimum spanning tree is
        % generated and the most cental node is determined
        start_node=min(Length_info(:,2)); %Start val for nodes
        
        % The algoritm for generating the minimum spanning tree is to take
        % the shortest connection and connect that by removing k_dist ro the 
        % positions (j,k) and (k,j) and add k_dist from nodes (j,j) and
        % (k,k). If the rank of the laplacian has increased by 1 then the
        % change is kept, if not take the next connection in the list and
        % try again. After n_cluster-1 iterations a minimum spanning tree
        % will have been genreated. 
        while j<length_cluster
            
            q_temp=q;
            rel_index1=Length_info(iter,2)-start_node+1;
            rel_index2=Length_info(iter,3)-start_node+1;
            
            q_temp(rel_index1,rel_index1)=q_temp(rel_index1,rel_index1)+k_dist;
            q_temp(rel_index2,rel_index2)=q_temp(rel_index2,rel_index2)+k_dist;
            
            q_temp(rel_index1,rel_index2)=q_temp(rel_index1,rel_index2)-k_dist;
            q_temp(rel_index2,rel_index1)=q_temp(rel_index2,rel_index1)-k_dist;
            
            if rank(q_temp)==rank(q)+1
                q=q_temp;
                j=j+1;
            end
            iter=iter+1;
        end
        
        % Algoritm to find the node with highest eigen value centrality
        % Remove the diagonal and find the leading eigenvector, i.e the one
        % with all entries of the same sign. The most central node is the
        % one correspponging to the largers magnitude in the leading
        % eigenvector. 
        q_temp=q-diag(diag(q));
        [v,~]=eig(q_temp);
        for k=1:length(v(:,1))
            if all(v(:,k)>0)
                [~,Cental_node]=max(v(:,end));
                break
            elseif all(v(:,k)<0)
                [~,Cental_node]=min(v(:,end));
                break
            end
        end               
        
        %[~,Cental_node]=max(diag(q)); % alternative way to find the most
        %central node, but more primitive. 
        
        Central_nodes(i)=Cental_node+start_node-1; % Store index of the most central node
    end
    
    % Add the local laplacian to the global
    Q=[Q, zeros(length(Q),length(q)); zeros(length(q),length(Q)) q];
end

%Conncect transmission

% Create a list of all distances between all central nodes with
% the indeces of the nodes in the second and third column. Then sort 
% the list
Cluster_dist_info=[];
for i=1:n_cluster
    if Central_nodes(i)>0
        pos1=Nodes(Central_nodes(i)).pos;
        for j=i+1:n_cluster
            if Central_nodes(j)>0
                pos2=Nodes(Central_nodes(j)).pos;
                dist=sqrt((pos1(1)-pos2(1))^2+(pos1(2)-pos2(2))^2);
                Cluster_dist_info=[Cluster_dist_info; [dist i j]];
            end
        end
    end
    
end

Cluster_dist_info=sortrows(Cluster_dist_info); % Sort the distances
n_cluster_old=n_cluster;

% To take care of situaton when some clusters don't have nodes
n_cluster=length(unique([Cluster_dist_info(:,2); Cluster_dist_info(:,3)])); 

if n_cluster<max(unique([Cluster_dist_info(:,2); Cluster_dist_info(:,3)]))
    Active_clusters=ismember(1:1:n_cluster_old,unique([Cluster_dist_info(:,2); Cluster_dist_info(:,3)]));
    iter=1;
    Cluster_dist_info_index=Cluster_dist_info(:,2:3);
    for i=1:n_cluster_old
        Central_nodes(iter)=Central_nodes(i);
        Cluster_dist_info_index(Cluster_dist_info_index==i)=iter;
        if Active_clusters(i)
            iter=iter+1;
        end
    end
    Cluster_dist_info(:,2:3)=Cluster_dist_info_index;
end

if n_cluster<2 % No transmission network is needed. 
    
else
    Q_cluster=-ones(n_cluster)+(n_cluster)*eye(n_cluster);
    % Try to remove the longest link, according to Cluster_dist_info, Then try
    % to remove every other link (one at the time) and see if rank drops, if it do don't remove
    % that link. Also try to remove every linek conncected to one specifivc 
    % node, If rank drops more than 1 dontremove the link.  Iterate until 
    % you have gone through the whole netwok 
    max_lines_trans=length(Cluster_dist_info(:,1));
    for i=1:max_lines_trans
        removal=1;
        Q_temp=Q_cluster;
        Node1=Cluster_dist_info(max_lines_trans-i+1,2);
        Node2=Cluster_dist_info(max_lines_trans-i+1,3);
        
        Q_temp(Node1,Node1)=Q_temp(Node1,Node1)-1;
        Q_temp(Node2,Node2)=Q_temp(Node2,Node2)-1;
        
        Q_temp(Node1,Node2)=Q_temp(Node1,Node2)+1;
        Q_temp(Node2,Node1)=Q_temp(Node2,Node1)+1;
        
        if rank(Q_temp)<n_cluster-1
            removal=0;
            break            
        end
        
        for j=1:n_cluster
            Q_temp_node_break=Q_temp;
            for k=1:n_cluster 
                if Q_temp(j,k)<0
                    Q_temp_node_break(k,k)=Q_temp_node_break(k,k)-1;
                end
            end
            Q_temp_node_break(j,:)=zeros(n_cluster,1);
            Q_temp_node_break(:,j)=zeros(1,n_cluster);
            
            Rank=rank(Q_temp_node_break);
            if Rank<n_cluster-2
                removal=0;
            end
        end
        
        if removal
            Q_cluster=Q_temp;
        end
        
    end
    
end

% Add the cluster connections to the Global laplacian
for i=1:length(Q_cluster)
    for j=1:length(Q_cluster)
        Q(Central_nodes(i),Central_nodes(j))=Q(Central_nodes(i),Central_nodes(j))+Q_cluster(i,j)*k_trans;       
    end
end

% Add the conncetions to be stored in the Node class for each node
for i=1:length(Q)
    Con=[];
    Nodes(i).connection=[];
    for j=1:length(Q)
        if Q(i,j)<0
            Con=[Con j];
        end
    end
    Nodes(i).connection=Con;
end

% Draw a map with the connctions 
if print
    % Plot the nodes as stars
    for i=1:n_cluster
        length_cluster=nnz(Index_matrix(:,i));
        x=zeros(length_cluster,1);
        y=zeros(length_cluster,1);
        for j=1:length_cluster
            pos=Nodes(Index_matrix(j,i)).pos;
            x(j)=pos(2);
            y(j)=pos(1);
        end
        plot(x,y,'*')
        hold on
    end
    
    % Plot the connections as black lines
    for i=1:length(Q)
        for j=i+1:length(Q)
            if Q(i,j)<0
                pos1=Nodes(i).pos;
                pos2=Nodes(j).pos;
                x=[pos1(2),pos2(2)];
                y=[pos1(1),pos2(1)];
                plot(x,y,'k-')
            end
        end
    end
    xlabel('width')
    ylabel('height')
    title('Map with buses and lines')

    hold off
end

end


