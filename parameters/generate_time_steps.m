function M = generate_time_steps(dx,dy,T,D,mesh_config)
% GENERATE_TIME_STEPS computes the number of time steps required to
% satisfy the monotonicity requirement of the deterministic models and 
% achieve non-negative transition probabilities for the stochastic
% frameworks.
% Inputs:
%   dx - spatial step in x-direction (scalar)
%   dy - spatial step in y-direction (scalar)
%   T - final time (scalar)
%   D - diffusion tensor (matrix)
%   mesh_config - mesh configuration (string)
%       Options: 'rectangular', 'flat-top' or 'pointy-top'
% Output:
%   M - no. of time steps (scalar)

% Extract tensor components
Dxx = D(1,1); Dxy = D(1,2); Dyy = D(2,2);

% Compute upper bound on time step duration
if isequal(mesh_config,'rectangular')
    max_tau = 2*dx^2*dy^2 / (3*dy^2*Dxx + 4*dx*dy*abs(Dxy) + 3*dx^2*Dyy);
elseif isequal(mesh_config,'flat-top')
    dy = 2*dy; % set dy to be full vertical step
    max_tau = min([dx^2*dy^2 / (2*(dy^2*Dxx + 3*dx^2*Dyy)), ...
        dx^2*dy^2/(3*(dy^2*Dxx + 2*dx*dy*abs(Dxy) + dx^2*Dyy)), ...
        2*dx^2*dy^2/(3*(dy^2*Dxx + 2*dx*dy*abs(Dxy) + 5*dx^2*Dyy))]);
else % pointy-top
    dx = 2*dx; % set dx to be full horizontal step
    max_tau = min([dx^2*dy^2 / (2*(dx^2*Dyy + 3*dy^2*Dxx)), ...
        dx^2*dy^2/(3*(dx^2*Dyy + 2*dx*dy*abs(Dxy) + dy^2*Dxx)), ...
        2*dx^2*dy^2/(3*(dx^2*Dyy + 2*dx*dy*abs(Dxy) + 5*dy^2*Dxx))]);
end

% Calculate minimum number of time steps
M = ceil(T / max_tau);

fprintf("No. of time steps: M = %g \n",M)

end