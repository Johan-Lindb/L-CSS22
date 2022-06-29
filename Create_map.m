function [M,W,H] = Create_map(w,h,Power,Mu,Sigma,print)
%Create_map Creates a map of a power dencity production and consumption
%   w and h specifies the width and height of the map.
%   The vector Power specifies the power production/consumption of an area,
%   Mu specifies the position and Sigma the varriance of the distributions
%   in the 2 directions. If print=1 the map will be plotted

[n,~]=size(Mu); % n is the number of distributions pos and neg.
H=0:1:h-1;      % Create discrete positions in heigth and width and make a grid from it
W=0:1:w-1;

[W2,H2]=meshgrid(W,H);

Map=[W2(:) H2(:)];
M=zeros(h,w);

% Add for each cluster add a normal distrubution to the the map. Mu gives
% the positon and Sigma gives the Covariance matrix. Since it is stacked 
% the correct rows must be picked 
for i=1:n
    M_add=Power(i)*mvnpdf(Map,Mu(i,:),Sigma(2*i-1:2*i,:));
    M=M +reshape(M_add',h,w);
end

% If it should be printed draw a map. For large maps the figure can be
% clutered and thefore the grid is removed in the surf if height or width
% is larger then 200. 
if print
    surf(W,H,M)
    xlabel('width')
    ylabel('height')
    zlabel('Power')
    title('Power dencity map')
    if h >200 || w>200
    shading('flat')
    end
end
