% Run 100 times with fixed number of nodes in each cluster. 

n_cluster=10;
n_nodes=100;

N_cluster=[7 20 4 8 11, 15 10 5 12 8]; % First 5 clusters are generating, last 5 are loads
w=120;          % Width of map
h=100;          % Height of map

H_big=6;        % Nuclear power 
H_medium=3;     % Hydro power
H_small=0.001;  % Wind/solar Power
f_0=50;         % Hz

Gam_matrix_sum=zeros(50); % Empty matrix to stree the values in

for k=1:100     % Iterate 100 times

Nodes=Node;     % Define Nodes as one instance of ny own clas Node
i=1;            % Start iteration for the nodes
H=zeros(1,n_cluster);
% Define empty matrices to strore data
M=zeros(n_nodes);

for c=1:n_cluster
    location=zeros(h,w);
    % Spread nodes at random all over the map for each clustrer. Since the
    % position of nodes don't effect the Laplacian K and the 
    % Laplacian K don't affect the metric gamma there is no need to have
    % nodes in a "cluster" close to each other. This will not generate as
    % neat maps at the one presented in the paper, but the metric is still
    % correct. The reason for doing this is because it was simpler and
    % could be done with the most randomness posible.
    while sum(sum(location))<N_cluster(c) % Make sure no node ended up in the same position as another one, in thet cluster
        posx=randi(h);
        posy=randi(w);
        location(posx,posy)=1;
    end
    
    row=0;
    col=1;
    for j=1:N_cluster(c)
        Nodes(i)=Node;
        Nodes(i).index=i;
        Nodes(i).cluster=c;
        % Iterate through location until the position of a node is found
        while 1 % if statement at the end gives a do while loop
            row=row+1;
            if row==h+1
                col=col+1;
                row=1;
            end
            if col>w % Happens once every 100 try with n_nodes=100 and grid =100x100
                col=w;
                break
            end
            if location(row,col)
                break
            end
        end
        
        Nodes(i).pos=[row,col];
        
        % To be able to compare the all the 100 iterations the first 2
        % clusters were assigned to be wind/solar, the third conventonal
        % and the two final hudro power. This was to mimik the result of
        % the one example in the paper. Since Conventional power often have
        % higher Power output than hydro power or solar/wind power the
        % P_avr was set to capture this. Powers of the individual nodes was
        % then picked fro a normal distribution. Looking at the P_avr the
        % nodes of wind/solar shpould be understood as a collection of
        % turbines, since these are often connected to a common
        % substation/transformer at a wind or solar park. 
        if c<3
            Nodes(i).type='Wind/Solar';
            H(c)=H_small;
            P_avr=1;
            Nodes(i).P=abs(normrnd(P_avr, abs(P_avr/4)));
        elseif c==3
            Nodes(i).type='Converntional';
            H(c)=H_big;
            P_avr=10;
            Nodes(i).P=abs(normrnd(P_avr, abs(P_avr/4)));
        elseif c<6
            Nodes(i).type='Hydro';
            H(c)=H_medium;
            P_avr=2;
            Nodes(i).P=abs(normrnd(P_avr, abs(P_avr/4)));
        else
            Nodes(i).type='Load';
            P_avr=2.1; % this gives on average a system with the same sum of production and load
            Nodes(i).P=-abs(normrnd(P_avr, abs(P_avr/4)));
        end
        i=i+1; 
    end
end

% Assign the "mass" to the mass matrix M it is proportional to the the ineria
% constant H and the rated power. 
for i=1:n_nodes
    cluster=Nodes(i).cluster;
    Mass=abs(H(cluster)*Nodes(i).P/(pi*f_0)); %M=2*H*P_rated/omega_system
    M(i,i)=Mass;
    Nodes(i).m=Mass;
end

% Get the line parametes, connect the lines and calculate the gamma-matrix
[k_dist, k_trans] = Line_paramters(Nodes,ones(1,10)); 

[K, Nodes, resorted]=Connect_lines(Nodes,zeros(n_cluster,1),1,k_dist, k_trans,0); 

[C, Cl, gam, Gam_matrix]=Metric_check(M,K,2,0); 

Gam_matrix_sum=Gam_matrix_sum+Gam_matrix;
k % print iterations. On my laptop each  iteratin took a few seconds. To make sure nothing had gotten stuck 
end

Gam_matrix=Gam_matrix_sum/100; % Take the average

% Plotting in Matlab
heatmap(log(Gam_matrix))


% Remove the axis values, since it become to cluttered. 
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));


