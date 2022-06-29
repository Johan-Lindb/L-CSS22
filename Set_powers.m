function [Powers, Nodes] = Set_powers(Nodes,Map,print)
%Set_powers Sets the powers positive or negative for the nodes, based on
%the power density map M. 

n_nodes=length(Nodes);      % Number of nodes

Powers=zeros(1,n_nodes);    % Vector to stroe the powers of all nodes

tol=0.1;        % Toleance
rate=1;         % 
% To get a balanced network from the randomly gernerated the powers will be
% multiplied with at plosfactor of minus factor (depending on if to is a 
% genrator or load). This is to get sum(production)=sum(consumption)
plus_factor=1;
minus_factor=1;
while rate>tol
    for i=1:n_nodes
        pos=Nodes(i).pos;
        P_avr=Map(pos(1),pos(2));
        if P_avr>0
        Powers(i)=plus_factor*abs(normrnd(P_avr, abs(P_avr/4)));
        else
            Powers(i)=-minus_factor*abs(normrnd(P_avr, abs(P_avr/4)));
        end
    end
    P_sum=sum(Powers);
    P_abs_sum=sum(abs(Powers));
    rate=abs(P_sum)/P_abs_sum;
    if P_sum>0
        minus_factor=minus_factor+0.1;
    else
        plus_factor=plus_factor+0.1;
    end
end

Powers=Powers-sum(Powers)/n_nodes; % to make it completely balanced

% Save the power in the Node
for i=1:n_nodes
    Nodes(i).P=Powers(i);
end

% This gives a rather stranfge looking map, but can be useful to see the
% amplityde of power generation/consumption
if print
    minus_factor % Printed out of curiosity
    plus_factor
    
    h=length(Map(:,1));
    w=length(Map(1,:));
    P_plot=NaN(h+1,w+1);
    
    for i=1:n_nodes
        pos=Nodes(i).pos;
        p=Nodes(i).P;
        P_plot(pos(1),pos(2))=p;
        P_plot(pos(1)+1,pos(2))=p;
        P_plot(pos(1),pos(2)+1)=p;
        P_plot(pos(1)+1,pos(2)+1)=p;
    end
    H=0:1:h;
    W=0:1:w;
    surf(W,H,P_plot)
    xlabel('width')
    ylabel('height')
    zlabel('Power')
    title('Node power map')
    
end
end

