function [part_renewable] = Part_wind(Nodes)
%Part_wind Calcualtes that part of the total power that comes from
%solar or wind. 

n_nodes=length(Nodes); 
Total_power=0;
Renewable_power=0;

for i=1:n_nodes
    p=Nodes(i).P;
    if p>0 % only add the power from gerating nodes. 
        if contains(Nodes(i).type, 'Wind')
            Renewable_power=Renewable_power+p;
        end
        Total_power=Total_power+p;
    end
end

part_renewable=Renewable_power/Total_power;

end

