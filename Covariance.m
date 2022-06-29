function [Sigma] = Covariance(range_var)
%Covariance gives a 2x2 positive definite covaraince matrix for a given varianvce 

a=range_var*rand(1); % Picking the variances as a positive random number between 0 and range_var
b=range_var*rand(1);

c=range_var*randn(1); % Picking a random c, this can be both positive and negative.

while a*b-c^2<0 % making sure that the Covariance matrix will be positive definite
    c=randn(1);
end

Sigma=[a c; c b]; % Constructing the covaranve matrix
end

