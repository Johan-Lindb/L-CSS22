function [Nodes, Location] = Pick_nodes(M,W,H,n,Mu,Sigma)
%Pick_nodes places nodes at random on the map of power dencity M.
%   n is the number of nodes to be placed. If print=1 a map will be
%   printed. Mu specifies the position and Sigma the varriance of the 
%   distributions in the 2 directions.
%   The output nodes is an array of node element with possition 
%   and cluster assigned. The output Location is a boolean matrix with 1 
%   at the locations of nodes.

h=length(H); % Height and width
w=length(W);

% Normalize the map to get a probability map
M_prob=abs(M)/sum(sum(abs(M)));

M_prob_vec=zeros(h*w,1);

% Stack the columns to a columnn vector
for i=1:w
    M_prob_vec((i-1)*h+1:i*h)=M_prob(:,i);
end

Location=zeros(h,w);

Index=1:1:length(M_prob_vec); % Vector with the index of each poition. Needed when picking a position 
n_picked=sum(sum(Location));  % Is zero here 

while sum(sum(Location))<n % Sometimes the same poition is picked twice. This makes sure n positions are picked
    
    Pick=randsample(Index,n-n_picked,true,M_prob_vec); % Pick n-n_picked nodes with propbability proportional to M_prob_vec
    
    for i=1:n-n_picked % Mark the positions that were picked with a 1 on the Location matrix
        col=floor(Pick(i)/h)+1;
        row=mod(Pick(i),h);
        if row==0
            row=h;
        end
        if col>w % Happens once every 100 try with n_nodes=100 and grid =100x100
            col=w;
        end
        Location(row,col)=1;
    end
    n_picked=sum(sum(Location));
end


Nodes=Node; % Introduce the class Node (a class defined i this project to store data)

row=1;
col=1;

[W2,H2]=meshgrid(W,H);

Map=[W2(:) H2(:)];

% Iterate and assign the properties of the nodes
Pos_map=zeros(h,w);
for i=1:n
    Nodes(i)=Node;
    Nodes(i).index=i;
    while 1 % if statement at the end gives a do while loop
        % Iterate until an active position is found in Location
        row=row+1;
        if row==h+1 % If a the end of a column start at the next one
            col=col+1;
            row=1;
        end
        if col>w % Happens once every 100 try with n_nodes=100 and grid =100x100
            col=w;            
            break
        end
        if Location(row,col)
            break
        end
    end
    
    Nodes(i).pos=[row,col]; % Assign position to the node
    
    cluster=[];
    
    % Assign the couster a node belongs to
    for j=1:length(Mu(:,1))
        M_add=mvnpdf(Map,Mu(j,:),Sigma(2*j-1:2*j,:));
        M_add=reshape(M_add',h,w);
        if M_add(row,col)>0.00001
            cluster=[cluster; j];
        end
    end
    
    % If a node can belong to several clusters, assign it to the one it
    % is the closest to
    if length(cluster)>1
        dist=zeros(length(cluster),1);
        for j=1:length(cluster)
            dist(j)=sqrt((row-Mu(cluster(j),1))^2+(col-Mu(cluster(j),2))^2);
        end
        [~,clu]=min(dist);
        cluster=cluster(clu);
    end
    if length(cluster)<1 % if the position did not belong to a cluster assign it 0
        cluster=0;
    end
    
    Nodes(i).cluster=cluster;
    
end

end

