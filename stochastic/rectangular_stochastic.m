function [Uavg,Ul,Uu,Npn] = rectangular_stochastic(xn,yn,dx,dy,V,D,T,M,f,Np,Ns)
% RECTANGULAR_STOCHASTIC conducts a random walk simulation of Np particles 
% on a rectangular lattice in the rectangular domain [0,Lx] x [0,Ly], 
% with no-flux boundary conditions, up to the final time T. Particle 
% movement is described by a symmetric diffusion tensor 
% D = [Dxx Dxy; Dxy Dyy]. After each duration tau, each particle takes a 
% step (in terms of dx and dy) to another fixed lattice site in the domain. 
% A total of Ns simulations are conducted.
% Inputs:
%   x - x-components of lattice sites (vector)
%   y - y-components of lattice sites (vector)
%   dx - horizontal spatial step (scalar)
%   dy - vertical spatial step (scalar)
%   V - control volume areas (vector)
%   D - diffusivity (4x4 matrix)
%   tau - time step duration (scalar)
%   T - final time (scalar)
%   M - no. of time steps (scalar)
%   f - initial condition (vector)
%   Np - no. of particles (scalar)
%   Ns - no. of simulations (scalar)
%   sx - x-index for one-dimensional slice (scalar)
%   sy - y-index for one-dimensional slice (scalar)
% Outputs:
%   U - particle density (matrix)
%   x_CI - confidence intervals for one-dimensional slice at sx (scalar)
%   y_CI - confidence interval for one-dimensional slice at sy (scalar)

% Setup
Nl = length(f); % no. of lattice sites
Nx = length(unique(xn)); % no. of nodes in x-direction
Ny = length(unique(yn)); % no. of nodes in y-direction
Uavg = zeros(Nl,M+1);
Ustd = zeros(Nl,M+1);
tau = T / M; % time step duration
Z = 2.576;

% Initial distribution of particles
Sp = f' * V / Np; % scaling constant
N0 = round(f .* V / Sp); % initial number of particles
Npn = sum(N0); Sp = f' * V / Npn;
P = define_cumulative_probabilities(xn,yn,dx,dy,D,Nl,Nx,Ny,tau);
fprintf('Number of particles scaled to Np = %g \n\n',Npn)

% Run simulations
N = repmat(N0,1,Ns); % initialise no. of nodes
U = N * Sp ./ V; % initial particle density
Uavg(:,1) = mean(U,2); Ustd(:,1) = std(U,[],2);

fprintf('Stochastic progress: ')
k = 0; % counter
for n = 1:M
    if mod(n,floor(M/10)) == 0
        k = k + 1;
        fprintf('%g%% ',k*10)
    end

    for s = 1:Ns
        Nt = zeros(size(N0)); % temporary storage of particles
        for i = 1:Nl
            for p = 1:N(i,s)

                prob = rand;

                if prob <= P(i,1) % i-1,j-1
                    Nt(i - (Nx+1)) = Nt(i - (Nx+1)) + 1;
                elseif prob <= P(i,2) % i,j-1
                    Nt(i - Nx) = Nt(i - Nx) + 1;
                elseif prob <= P(i,3) % i+1,j-1
                    Nt(i - (Nx-1)) = Nt(i - (Nx-1)) + 1;
                elseif prob <= P(i,4) % i-1,j
                    Nt(i - 1) = Nt(i - 1) + 1;
                elseif prob <= P(i,5) % i+1,j
                    Nt(i + 1) = Nt(i + 1) + 1;
                elseif prob <= P(i,6) % i-1,j+1
                    Nt(i + (Nx-1)) = Nt(i + (Nx-1)) + 1;
                elseif prob <= P(i,7) % i,j+1
                    Nt(i + Nx) = Nt(i + Nx) + 1;
                elseif prob <= P(i,8) % i+1,j+1
                    Nt(i + (Nx+1)) = Nt(i + (Nx+1)) + 1;
                else % otherwise, particle stays put
                    Nt(i) = Nt(i) + 1;
                end

            end
        end
        N(:,s) = Nt; % no. of particles at each lattice site at time tn
        U(:,s) = N(:,s) .* Sp ./ V; % particle density at time tn
    end

    Uavg(:,n+1) = mean(U,2); Ustd(:,n+1) = std(U,[],2);

end
fprintf('\n\n')


Ul = Uavg - Z*Ustd / sqrt(Ns); Uu = Uavg + Z*Ustd / sqrt(Ns);

% Compute average density and CI bounds
Uavg = reshape(Uavg,Nx,Ny,M+1); % average density over Ns simulations
Ul = reshape(Ul,Nx,Ny,M+1); % reshaped lower CI bound
Uu = reshape(Uu,Nx,Ny,M+1); % reshaped upper CI bound


end

%% Subfunctions
function P = define_cumulative_probabilities(xn,yn,dx,dy,D,Nl,Nx,Ny,tau)
% RECOVER_INDICES recovers the position indices for each particle recorded
% in the initial condition N0.
% Inputs:
%   xn - x-component of lattice sites (vector)
%   yn - y-component of lattice sites (vector)
%   dx - spatial step in the x-direction (scalar)
%   dy - spatial step in the y-direction (scalar)

Dxx = D(1,1); Dxy = D(1,2); Dyy = D(2,2); % tensor components
P = zeros(Nl,8); % cumulative probabilities for each particle

% General expressions for transition probabilities
Prh = tau*(3*dy^2*Dxx - dx^2*Dyy) / (dx^2*dy^2);
Prv = tau*(3*dx^2*Dyy - dy^2*Dxx) / (dx^2*dy^2);
Prdp = tau*(dy^2*Dxx + 4*dx*dy*Dxy + dx^2*Dyy)/(dx^2*dy^2);
Prdm = tau*(dy^2*Dxx - 4*dx*dy*Dxy + dx^2*Dyy)/(dx^2*dy^2);

for s = 1:Nl
    i = round(xn(s) / dx) + 1; % x-index
    j = round(yn(s) / dy) + 1; % y-index

    if i == 1 % western boundary
        if j == 1 % southwest corner
            P(s,5) = Prh/2; P(s,7) = Prv/2; P(s,8) = Prdp/2;
        elseif j == Ny % northwest corner
            P(s,2) = Prv/2; P(s,3) = Prdm/2; P(s,5) = Prh/2;
        else % western edge
            P(s,2) = Prv/4; P(s,3) = Prdm/4; P(s,5) = Prh/2;
            P(s,7) = Prv/4; P(s,8) = Prdp/4;
        end
    elseif i == Nx % eastern boundary
        if j == 1 % southeast corner
            P(s,4) = Prh/2; P(s,6) = Prdm/2; P(s,7) = Prv/2;
        elseif j == Ny % northeast convex corner
            P(s,1) = Prdp/2; P(s,2) = Prv/2; P(s,4) = Prh/2;
        else % eastern edge
            P(s,1) = Prdp/4; P(s,2) = Prv/4; P(s,4) = Prh/2;
            P(s,6) = Prdm/4; P(s,7) = Prv/4;
        end
    elseif j == 1 % southern boundary
        P(s,4) = Prh/4; P(s,5) = Prh/4; P(s,6) = Prdm/4;
        P(s,7) = Prv/2; P(s,8) = Prdp/4;
    elseif j == Ny % northern boundary
        P(s,1) = Prdp/4; P(s,2) = Prv/2; P(s,3) = Prdm/4; 
        P(s,4) = Prh/4; P(s,5) = Prh/4;
    else % interior
        P(s,1) = Prdp/8; P(s,2) = Prv/4; P(s,3) = Prdm/8; P(s,4) = Prh/4;
        P(s,5) = Prh/4; P(s,6) = Prdm/8; P(s,7) = Prv/4; P(s,8) = Prdp/8;
    end

end

% Define cumulative probabilities
for i = 2:8
    P(:,i) = P(:,i) + P(:,i-1);
end

end