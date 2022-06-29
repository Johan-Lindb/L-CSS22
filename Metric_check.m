function [C,Cl,gam,Gam_matrix] = Metric_check(M,K,metric_type,print)
%Metric_check given the Inertia matrix M and line matrix K (the laplacian)
%and a metric type of either "2" or "inf" calculates the gamma measure of
%either 2-norm or infinity-norm
%   Returns C as the controller and Cl as the closed loop transfer function
%   gam is the metric for the whole system and Gam_matrix gives the matrix
%   showing how disturbences propagate between nodes. 

n_red=nnz(M);   % Number of generating nodes

M_red=M(1:n_red,1:n_red); % reduced mass matrix onlyy looking at the generating nodes

K_red=K(1:n_red,1:n_red)-(K(1:n_red,n_red+1:end)/K(n_red+1:end,n_red+1:end))*K(n_red+1:end,1:n_red); % Create the Krohn reduced laplacing (Schur comlpiment)

[l,D]=ldl(K_red);   % Decompse K_red to K_red=lDl'

L=l*sqrt(D);        % Get L so that K_red=LL'
L=L(:,1:end-1);     % remove the last line to to remove the eigenvalue corresponding to a uniform rotation of the phases of all nodes. 

% Get the system matrices to fit in the framework of a transfer function
A=[zeros(n_red) -inv(M_red)*L;L' zeros(n_red-1)];

B=[inv(M_red); zeros(n_red-1,n_red)];
C=[eye(n_red) zeros(n_red,n_red-1)];

P=ss(A,[B zeros(2*n_red-1,n_red) B],[C; zeros(n_red, 2*n_red-1);C],[zeros(n_red), zeros(n_red), zeros(n_red); zeros(n_red,2*n_red), eye(n_red) ;zeros(n_red), eye(n_red), zeros(n_red)]);

if metric_type==2 % if 2-norm is of interest
    [C,Cl,gam]=h2syn(P,n_red,n_red); % use h2syn to get the optimal controller and closed loop system
    %gam1= sqrt(2*trace(B*B')) % Can also be used to calculate the gamma metric
    Gam_matrix=zeros(n_red);
    
    % Given the optimal closed loop transfer function, calculate the 2-norm
    % gain from all nodes to all other
    for i=1:n_red
        for j=1:n_red
            Gam_matrix(i,j)=norm(Cl([i;i+n_red],[j;j+n_red]),2);
        end
    end
    
elseif isinf(metric_type) % if inf-norm is of interest
    %[C,Cl,gam]=hinfsyn(P,n_red,n_red);
    Cl=lft(P,-sqrt(2)*eye(n_red)); % using hinfsyn took too long time and the analytical optimal controller were used directly
    Gam_matrix=zeros(n_red);
    
    % Given the optimal closed loop transfer function, calculate the 2-norm
    % gain from all nodes to all other
    for i=1:n_red
        for j=1:n_red
            Gam_matrix(i,j)=norm(Cl([i;i+n_red],[j;j+n_red]),Inf);
        end
    end
    gam=max(max(Gam_matrix));
else  % If no correct norm has ben intered return this
    C=0;
    Cl=0;
    gam=inf;
    Gam_matrix=0;
end

% Print a heatmap showing the disturebence gains. 
if print
    if metric_type==2
        log_gam=log(Gam_matrix);
        log_gam=round(log_gam,1);
        heatmap(log_gam)
        title('Gamma matrix: ln(\gamma_{H_2, ik})')
    elseif isinf(metric_type)
        if length(Gam_matrix)<10
            heatmap(round(Gam_matrix,1))
        else
            heatmap(Gam_matrix)
        end
        title('Gamma matrix: \gamma_{H_\infty, ik}')
    end
    % If the heatmap becomes larger than 10x10 then don't print the 
    % the axis values
    if length(Gam_matrix)>10 
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    end    
end


end

