% Run script. Does the same as Run_simulation, but as a script

% One instance of the simulation. Create a new model and
% at the end produce a gamma metric and a gamma matrix for either H2 norm
% or H_inf norm with both the full system, the lumed model and with
% increased mass of the renewables

w=120;          % Width of the map
h=100;          % Height of the map
n_cluster=10;   % Number of clusters
n_nodes=100;    % Number of nodes
Power_amp=100;  % Power amplitude, can be seen as some sort of per unit base.
range_var=20;   % Range of variance of positions of nodes in a cluster
prob_wind=0.5;  % Probability of a Cluster bring wind/solar, rather than hydorpower
metric_type=2;  % Type of norm used. 2-norm or inf norm
print=1;        % Print figures or not 1=printing

% The script runs in the following way: 
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

%%

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
%%
if print
    figure(3)
end

% Calculate the Line parameters
[k_dist, k_trans] = Line_paramters(Nodes,Power);

% Connect the lines and get the Laplaceian. This sometimes don't the nodes
% in some clusters, but they are present in the rest of the steps. 
[K, Nodes, resorted]=Connect_lines(Nodes,Mu,k_dist, k_trans, print);

%%

% Set masses/inertia depending on the type of power gerenation
[M,Nodes]=Set_masses(Nodes,Power,Sigma,prob_wind);

% Restructure to get Solar/wind in first clusters
[Nodes,M,K, P] = Restructure_Clusters(Nodes, M, K, P, Mu);

% Calcuate the part thet is solar/wind of the total power priduction
part_wind=Part_wind(Nodes);
%%

% Increase the inertia of wind and solar by a factor
factor=100;
[M_large,Nodes_Large] = Increase_mass_renewable(Nodes, factor);

%%
if print
    figure(4)
end

% Get the gamma metric and the gamma matrix for 2 or inf norm
[C, Cl, gam, Gam_matrix]=Metric_check(M,K,metric_type,print);

% Create lumed models
[M_red,K_red] = Cluster_ekviv(Nodes,K,Mu);
%%
if print
    figure(5)
end

%%

% Get the gamma metric and the gamma matrix for the lumped system
[C_red, Cl_red, gam_red, Gam_matrix_red]=Metric_check(M_red,K_red,metric_type,print);

ratio=gam/gam_red; % The ratio of the gamma parameter between unlumped and lumped
%%
if print
    figure(6)
end

% Get the gamma metric and the gamma matrix for 2 or inf norm fro the
% system with increased (or decreased) inertia
[C_large, Cl_large, gam_large, Gam_matrix_large]=Metric_check(M_large,K,metric_type,print);



