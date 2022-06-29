function [gam, Gam_matrix, gam_red,Gam_matrix_red,part_renewable,Nodes] = Run_simulation(w,h,n_cluster,n_nodes, Power_amp, range_var, prob_wind,metric_type,print)
%Run_simulation run one instance of the simulation. Create a new model and
%at the end produce a gamma metric and a gamma matrix for either H2 norm
%or H_inf norm
% w is the width of the map, h the height. n_cluster and n_nodes are the
% number of clusters and nodes in the model, Power_amp an amplitude in
% power, can be seen as some sort of per unit base. range_var gives the
% range of varaince for the locations of nodes in a cluster. prob_wind
% gives the probably of a power producing area of being wind or solar 
% power, rather than hydro power. metric type is either 2 or inf, to chose
% between 2-norm or infinity norm. If print=1, plots are drawn, if print=0
% they are not. 

close all

% The function runs in the following way: 
% First generate a dencity map. This
% is representing a country or region, where a positive value corresponds
% to the ability to generate power eg a river that has the potential to
% generate a lot of hydropower or an area with a lot of wind, like an
% offshore area for wind turbines. There will also be generated
% conventrional power producing areas, theet have a very small varaince in
% location in space. The inspiration for the systems is the Nordic
% Synchronous area where the power sourses are (in decreasing order) Hydo
% power, Nuclear power, wind power, other heat power (gas, coal and oil),
% and solar power. 
% Then pick Nodes that will be either loads or generators
% Set the powers of each node and calculate the line parameters needed for
% a network with said nodes. 
% Connect the lines cor construct a transmission and sub-transmission
% networks and get a directed Laplacian for calculations
% Set the masses the get an Inertia matrix
% Sort the clusters to get solar and wind power in the first clusters
% Genrate the gamma metric and the gamma matrix to se how disturbences
% propagate 
% Generate a lumped system where each cluster has been replaced by an
% equivalnt node and get the gamma metric and matrix again. 
% If of interest, get the part of the total amound of power productions 
% that comes from solar/wind power 

% Get the "power volume", positions of the centre and the covariance matrix
% for generation and load areas. 
[Power, Mu, Sigma] = Power_dencity(n_cluster, Power_amp, range_var, w,h);

if print
    figure(1)
end
% Generate a map from the previous results
[Map,W,H]=Create_map(w,h,Power,Mu,Sigma,print);

% Pick Nodes. Higher magnitude in Map gives higher probability of a
% position being picked
[Nodes, Location] = Pick_nodes(Map,W,H,n_nodes,Mu,Sigma);

% Sort the nodes to get nodes in the same cluster to have indeces after to
% each other
Nodes=Sort_nodes(Nodes,Mu);

if print
    figure(2)
end

% Set the power production/consumption of each node
[P, Nodes] = Set_powers(Nodes,Map,print);

% Calculate the Line parameters
[k_dist, k_trans]=Line_paramters(Nodes,Power);

if print
    figure(3)
end

% Connect the lines and get the Laplaceian. This sometimes don't the nodes
% in some clusters, but they are present in the rest of the steps. 
[K, Nodes, resorted]=Connect_lines(Nodes,Mu,k_dist, k_trans, print);

% Set masses/inertia depending on the type of power gerenation
[M,Nodes]=Set_masses(Nodes,Power,Sigma,prob_wind);

% Restructure to get Solar/wind in first clusters
[Nodes,M,K, P] = Restructure_Clusters(Nodes, M, K, P, Mu);

if print
    figure(4)
end

disp("Start optimal control") % To show that the function has rached here and not gotten stuck in the previous steps

% Get the gamma metric and the gamma matrix for 2 or inf norm
[C, Cl, gam, Gam_matrix]=Metric_check(M,K,metric_type,print);

% Create lumed models
[M_red,K_red] = Cluster_ekviv(Nodes,K,Mu);

if print
    figure(5)
end

% Get the gamma metric and the gamma matrix for the lumped system
[C_red, Cl_red, gam_red, Gam_matrix_red]=Metric_check(M_red,K_red,metric_type,print);

% Calcuate the part thet is solar/wind of the total power priduction
part_renewable=Part_wind(Nodes);


end

