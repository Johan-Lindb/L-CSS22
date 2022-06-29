function [M,Nodes] = Increase_mass_renewable(Nodes, factor)
%Increase the mass of reneabels with factor. factor<1 gives reduced mass
%   Returns an updated mass matrix M and Nodes with new masses. 
n_nodes=length(Nodes);
M=zeros(n_nodes);
for i=1:n_nodes
    if contains(Nodes(i).type, 'Wind')
        Nodes(i).m=Nodes(i).m*factor;
    end
    M(i,i)=Nodes(i).m;
end

end

