function [A] = adjacency_mat(x,y,Gamma,alpha,type)
% Adjacency matrix calculator (2D fluid flow)
% Date Created  : 05/20/2016
% Date Modified : 08/24/2017, alpha, vectorized, type
% By MGM
% Reference:    A. G. Nair & K. Taira
%               "Network-Theoretic Approach to Sparsified Discrete Vortex
%               Dynamics", Journal of Fluid Mechanics, 768, 549-571, 2015
%               (doi.org/10.1017/jfm.2015.97)
%%
% [A] = adjacency_mat(x,y,gamma,alpha,type)
% Function to calculate adjacency matrix given x, y and gamma. 
% A         : adjacency matrix using algebraic mean (default)
% x         : x coordinate vector
% y         : y coordinate vector
% gamma     : circulation vector
% alpha     : network direction parameter (0.5=undirected)
% type      : algebraic or geometric mean ('alg' [default] or 'goem')

if isempty(type)
    type = 'alg';
end
pts = length(x);
x = reshape(x,pts,1); y = reshape(y,pts,1);
r(1:pts,1:2) = [x y];               % position vector matrix, r = (xi,yj)
Gamma = Gamma/(2*pi);               % scalling for Gamma in 2D flow

fprintf('Calculating A matrix of %i nodes...\n',length(x));
A = NaN(pts,pts);
for i = 1:pts   % can use parfor if needed
    dist_ij = r(1:pts,:) - ones(pts,1)*r(i,:);  % distance vector, r_{i->j}
    dist_mag = sqrt(sum(dist_ij.^2,2));         % ||r_{i->j}||
    dist_ij = 1./dist_mag;                      % distance ratio
    dist_ij(isinf(dist_ij)) = 0;                % update r_{i->i} = 0
    u_ij = abs(Gamma(i)) * dist_ij;             % u_{i->j}
    u_ji = abs(Gamma) .* dist_ij;               % u_{j->i}
    % A_{i->j}
    if strcmp(type,'alg') == 1
        A(i,:) = alpha*u_ij + (1-alpha)*u_ji;
    elseif strcmp(type,'geom') == 1
        A(i,:) = alpha*u_ij .* (1-alpha)*u_ji;
    else
        error('Eneter valid network definition "type"')
    end
end
