function [k_dist, k_trans] = Line_paramters(Nodes,Power)
%Line_parameter gives the line paramers given how much power is produced or
%consumed in the nodes and clusters
%   The line parametersars are of 2 tyes to represent that that is often
%   the case with stnadardized lines for diffent voltages. the line
%   parameters are calcularted using the formula
%   P_tranmission = V_recieving*V_sending/X_line*sin(phase diffence) =  
%   P_max*sin(phase diffenrce). After linearization this becomes:
%   P_transmission=P_max*cos(Phase_diffence_stationary)*deviation from
%   stationary phase diffence

n_cluster=length(Power);            % Number of clusters
Cluster_power=zeros(1,n_cluster);   % Power in the clusters
n_nodes=length(Nodes);              % Number of nodes

% Calculate the total sum of power in all cluster
for i=1:n_nodes
    Cluster_power(Nodes(i).cluster)=Cluster_power(Nodes(i).cluster)+Nodes(i).P;
end

% Assume operating at 15 deg. In the cluster with the most power 
% production/consumption the line from the cental node must be able to 
% transport the half out/in. For a cluster with one big and all other much
% smaller this will in reality give angle corresponding to much more than
% 15 degrees, but never more than 45, since sin(15 deg)=0.26 and sin(45
% deg)=0.71
k_dist=max(abs(Cluster_power))/2*cos(pi/12); 

% Assume operating at 15 deg. N-1 security gives that all power in/out of 
% cluster must be able to go from one line. For an unblanaced system the
% angle will in reality be larger, is highly unlikely to be over 45 deg at
% normal operation and unlikely at N-1 contingency
k_trans=max(abs(Cluster_power))*cos(pi/12);
end

