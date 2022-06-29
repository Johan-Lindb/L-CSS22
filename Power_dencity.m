function [Power, Mu, Sigma] = Power_dencity(n_cluster, Power_amp, range_var, w,h)
%Power_dencity Given the number of clusters,a general power amplitude, the
%variance and the size of the map, return volumes, positions and covarnance
%matrices of 2D standard deviations


Power=zeros(1,n_cluster);   % The volumes of the "bells"
Mu=zeros(n_cluster,2);      % The positions of the center
Sigma=zeros(2*n_cluster,2); % The covariance matrices, stacked


for i=1:n_cluster
    Power(i)=sign(n_cluster/2+0.1-i)*Power_amp*rand(1);     % Make half of the clusters generating and half of them loads. If the number of clusters is an odd nomber, make one more generating cluster than loads
    Mu(i,1)=round(w*rand(1));                               % Place the center at random within the map
    Mu(i,2)=round(h*rand(1));
    Sigma(2*i-1:2*i,:)=Covariance(range_var*(i/n_cluster)); % Get a Covaraince. To make sure that some generating clusters are small in size, ie corresponding to conventional power plants, the varainve is shifting between the clusters. Cities are ususlly large in size (compared to power plants and wind farms) and thus have the largest varaince. 
    
end

end

