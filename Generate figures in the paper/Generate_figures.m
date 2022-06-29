
% Script to geneate the figures in the paper. The figures in the paper were
% drawn in Tikz. The figures produced her are the ones in the paper but
% plotted by Matlab. To get the figures in Tikz, save the Gamma_matrix
% files as text-tiles and let Tikz read from them. The White squares and
% labels in figure 2 were added in Tikz, and are not produced here. 

% To get an interpretation of what eveytng is look in the functions, and
% read the comments in Run_script.m or Run_simulation.m

% The full genration of the figures takes a minute or two.

load('Workspace_figures') % Everything here is not used and some things are recalculated when using the functions

figure(1)

% Generate figure 2
[K, Nodes, resorted]=Connect_lines(Nodes,Mu,k_dist, k_trans, print);

colormap parula % the colormap used in the paper. Works good in black and white and can be seen by colorblind
figure(2)

% Generate figure 3a
[C, Cl, gam, Gam_matrix]=Metric_check(M,K,2,1);
colormap parula % the colormap used in the paper. Works good in black and white and can be seen by colorblind

%%
figure(3)
% Generate figure 3b
[C, Cl, gam, Gam_matrix]=Metric_check(M,K,inf,1);
colormap parula % the colormap used in the paper. Works good in black and white and can be seen by colorblind

figure(4)
[M_red,K_red] = Cluster_ekviv(Nodes,K,Mu); % Get the lumped model
%%
% Generate figure 5a
[C_red, Cl_red, gam_red, Gam_matrix_red]=Metric_check(M_red,K_red,2,1);
colormap parula % the colormap used in the paper. Works good in black and white and can be seen by colorblind


figure(5)
% Generate figure 5b
[C_red, Cl_red, gam_red, Gam_matrix_red]=Metric_check(M_red,K_red,inf,1);
colormap parula % the colormap used in the paper. Works good in black and white and can be seen by colorblind


figure(6)
load('Gam_matrix_avr.mat') % Load the average produced from a run of Average_run.m

% Generate figure 4
heatmap(log(Gam_matrix_avr))
colormap parula % the colormap used in the paper. Works good in black and white and can be seen by colorblind

% Remove the axis values, since it becomes to cluttered. 
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
title('Gamma matrix: ln(average(\gamma_{H_2, ik}))')

