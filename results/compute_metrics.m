function [epsilon,MSE] = compute_metrics(c,Uavg,Ul,Uu,M,mesh_config)
% COMPUTE_METRICS computes a quantitative metric for assessing the fit
% between the stochastic and deterministic frameworks. For each time step,
% a value of one is assigned to a lattice site if the particle density
% obtained from the determinsitic framework lies outside the 95% confidence
% bounds obtained from the stochastic simulations. These lattice sites are
% then counted and divided by the total number of lattice sites to give a
% quantitative metric.
% Inputs:
%   c - deterministic particle density (matrix)
%   Uavg - stochastic particle density (matrix)
%   Ul - lower 95% CI bound (matrix)
%   Uu - upper 95% CI bound (matrix)
%   mesh_config - mesh configuration (string)
%       Options: 'rectangular', 'flat-top' or 'pointy-top'
% Outputs:
%   epsilon - confidence interval error (scalar)
%   MSE - mean squared error (scalar)

epsilon = zeros(M,1); % initialise metric storage (excluding t = 0)
MSE = zeros(M,1); % mean-squared error

if isequal(mesh_config,'rectangular')
    [Nx,Ny,~] = size(c);
    for n = 2:M+1
        squared_error = zeros(Nx,Ny);
        conf_error = zeros(Nx,Ny);
        p = 0;
        for i = 1:Nx
            for j = 1:Ny
                if Ul(i,j,n) > 0
                    p = p + 1;
                    if c(i,j,n) < Ul(i,j,n) || c(i,j,n) > Uu(i,j,n)
                        conf_error(i,j) = 1;
                    end
                end
                squared_error(i,j) = (c(i,j,n) - Uavg(i,j,n))^2;
            end
        end
        epsilon(n-1) = sum(conf_error,'all') / p;
        MSE(n-1) = mean(squared_error,'all');
    end
else % hexagonal configurations
    Nl = size(c,1); % no. of lattice sites
    for n = 2:M+1
        squared_error = zeros(Nl,1);
        conf_error = zeros(Nl,1);
        p = 0;
        for i = 1:Nl
            if Ul(i,n) > 0
                p = p + 1;
                if c(i,n) < Ul(i,n) || c(i,n) > Uu(i,n)
                    conf_error(i) = 1;
                end
            end
            squared_error(i) = (c(i,n) - Uavg(i,n))^2;
        end
        epsilon(n-1) = sum(conf_error) / p;
        MSE(n-1) = mean(squared_error);
    end
end

% Average over time points
epsilon = mean(epsilon); MSE = mean(MSE);

end