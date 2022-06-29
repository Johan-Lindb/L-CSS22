function [M,Nodes] = Set_masses(Nodes,Power,Sigma,prob_wind)
%Set_masses Calculate the mass in the mass spring equivalent 
%   Detailed explanation goes here

% Ineria constant unit second
H_big=6;        % Nuclear power 
H_medium=3;     % Hydro power
H_small=0.001;  % Wind/solar Power
f_0=50;         % nominal frequency Hz

n_nodes=length(Nodes);      % Number of nodes
M=zeros(n_nodes,n_nodes);   % Matrix to store the masses

H=zeros(1,length(Power));   % Vector to store the inertia constants
generator_type=strings(1,length(Power)); % Vector to stroe the type of generator/load in each cluster

for i=1:length(H)
    
    if Power(i)>0
        % If very small varance, set it to be a covnetional power plant. 
        if max(max(Sigma(2*i-1:2*i,:)))<1.5
            H(i)=H_big;
            generator_type(i)='Conventional Power';
        else
            r=rand(1);
            % Genrate a random number, if higher than prob_wind set it to
            % hydro power, if lower set it to wind/solar
            if r>prob_wind
                H(i)=H_medium;
                generator_type(i)='Hydro Power';
            else
                H(i)=H_small;
                generator_type(i)='Wind/solar Power';
            end
            
        end
        
    else
        % If power<0 it is a load, with no inertia
        H(i)=0;
        generator_type(i)='Load';
    end
end

% Iterate through the nodes and set the type and the mass of each according 
% to the general formula see for example E. Ørum, M. Kuivaniemi, M. Laasonen, A. I. Bruseth, E. A. Jansson, A. Danell, K. Elkington, and N. Modig, “Future system inertia,” Entsoe, Tech. Rep., 2018.
for i=1:n_nodes
    cluster=Nodes(i).cluster;
    Mass=abs(H(cluster)*Nodes(i).P/(pi*f_0)); %M=2*H*P_rated/omega_system
    M(i,i)=Mass;
    Nodes(i).m=Mass;
    Nodes(i).type=generator_type(cluster);
end
end

