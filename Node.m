% Defines a class Node that contains the relevant information for for each
% node. This class is used for storing data about the nodes. 
classdef Node 
    properties 
        index=0;        % Index of the node
        pos=[];         % Position of the node
        m=0;            % Mass of the node in the mass-spring equivalnce. This comes from the ineria in the node
        P=0;            % Power production of node. Consumptions gives neagrive P
        cluster=0;      % The cluster a node bleongs to
        connection=[];  % The nodes this node is connected to
        type=string;    % Type or node. It can be "Conventional" "Hydro", "Wind/solar" or "Load"
    end
end