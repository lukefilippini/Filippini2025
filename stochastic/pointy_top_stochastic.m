function [Uavg,Ul,Uu,Npn] = pointy_top_stochastic(xn,yn,dx,dy,V,D,T,M,f,Np,Ns)
% POINTY_TOP_STOCHASTIC conducts a random walk simulation of Np particles 
% on a pointy-top hexagonal lattice in the rectangular domain [0,Lx] x 
% [0,Ly], with no-flux boundary conditions, up to the final time T. 
% Particle movement is described by a symmetric diffusion tensor 
% D = [Dxx Dxy; Dxy Dyy]. After each duration tau, each particle takes a 
% step (in terms of dx and dy) to another fixed lattice site in the domain. 
% A total of Ns simulations are conducted.
% Inputs:
%   Nx - no. of nodes in the x-direction (scalar)
%   Ny - no. of nodes in the y-direction (scalar)
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
% Outputs:
%   U - particle density (matrix: no. of lattice sites by time steps)

% Setup
Nl = length(f); % no. of lattice sites
Nx = length(unique(xn)); % no. of nodes in x-direction
Ny = length(unique(yn)); % no. of nodes in y-direction
Uavg = zeros(Nl,M+1);
Ul = zeros(Nl,M+1);
Uu = zeros(Nl,M+1);
tau = T / M;

% Initial distribution of particles
Sp = f' * V / Np; % scaling constant
N0 = round(f .* V / Sp); % initial number of particles
Npn = sum(N0); Sp = f' * V / Npn;
P = define_cumulative_probabilities(xn,yn,dx,dy,D,Nl,Nx,Ny,tau);
fprintf('No. of scaled particles: Np = %g \n',Npn)

% Run simulations
N = repmat(N0,1,Ns); % initialise no. of nodes
U = N * Sp ./ V; % initial particle density
Uavg(:,1) = mean(U,2);

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
                xi = round(xn(i) / dx) + 1; % x-index

                if prob <= P(i,1) % i-1,j-1
                    if mod(xi,2) == 0
                        Nt(i - Ny/2) = Nt(i - Ny/2) + 1;
                    else
                        Nt(i - (Ny/2+1)) = Nt(i - (Ny/2+1)) + 1;
                    end
                elseif prob <= P(i,2) % i-1,j+1
                    if mod(xi,2) == 0
                        Nt(i - (Ny/2-1)) = Nt(i - (Ny/2-1)) + 1;
                    else
                        Nt(i - Ny/2) = Nt(i - Ny/2) + 1;
                    end
                elseif prob <= P(i,3) % i,j-2
                    Nt(i-1) = Nt(i-1) + 1;
                elseif prob <= P(i,4) % i,j+2
                    Nt(i+1) = Nt(i+1) + 1;
                elseif prob <= P(i,5) % i+1,j-1
                    if mod(xi,2) == 0
                        Nt(i + Ny/2) = Nt(i + Ny/2) + 1;
                    else
                        Nt(i + (Ny/2-1)) = Nt(i + (Ny/2-1)) + 1;
                    end
                elseif prob <= P(i,6) % i+1,j+1
                    if mod(xi,2) == 0
                        Nt(i + (Ny/2+1)) = Nt(i + (Ny/2+1)) + 1;
                    else
                        Nt(i + Ny/2) = Nt(i + Ny/2) + 1;
                    end
                else % otherwise, particle stays put
                    Nt(i) = Nt(i) + 1;
                end

            end
        end

        N(:,s) = Nt; % no. of particles at each lattice site at time tn
        U(:,s) = N(:,s) .* Sp ./ V; % particle density at time tn

    end

    Uavg(:,n+1) = mean(U,2);
    Ul(:,n+1) = quantile(U,0.025,2); Uu(:,n+1) = quantile(U,0.975,2);

end
fprintf('\n')

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

dx = 2*dx;
Dxx = D(1,1); Dxy = D(1,2); Dyy = D(2,2); % tensor components
P = zeros(Nl,6); % cumulative probabilities for each particle

% General expressions for transition probabilities
Ppv = 2*tau*(dx^2*Dyy - dy^2*Dxx) / (dx^2*dy^2);
Ppdp = 2*tau*(dy*Dxx + dx*Dxy)/(dx^2*dy);
Ppdm = 2*tau*(dy*Dxx - dx*Dxy)/(dx^2*dy);

for s = 1:Nl
    i = round(xn(s) / (dx/2)) + 1; % x-index
    j = round(yn(s) / (dy/2)) + 1; % y-index

    if j == 1 % southern concave boundary
        if i == 1 % southwest concave corner
            P(s,6) = 3*Ppdp; P(s,4) = 3*Ppv/2;
        elseif i == Nx % southeast concave corner
            P(s,2) = 3*Ppdm; P(s,4) = 3*Ppv/2;
        else % concave boundary
            P(s,2) = 3*Ppdm/2; P(s,6) = 3*Ppdp/2; P(s,4) = 3*Ppv/2;
        end
    elseif j == Ny % northern concave boundary
        P(s,3) = 3*Ppv/2; P(s,1) = 3*Ppdp/2; P(s,5) = 3*Ppdm/2;
    elseif j == 2 % southern convex boundary
        P(s,1) = 3*Ppdp/4; P(s,5) = 3*Ppdm/4; P(s,2) = 3*Ppdm/2;
        P(s,6) = 3*Ppdp/2; P(s,4) = 3*Ppv/4;
    elseif j == Ny-1 % northern convex boundary
        if i == 1 % northwest convex corner
            P(s,3) = 3*Ppv/4; P(s,5) = 3*Ppdm; P(s,6) = 3*Ppdp/2;
        elseif i == Nx % northeast convex corner
            P(s,3) = 3*Ppv/4; P(s,1) = 3*Ppdp; P(s,2) = 3*Ppdm/2;
        else % convex boundary
            P(s,3) = 3*Ppv/4; P(s,1) = 3*Ppdp/2; P(s,5) = 3*Ppdm/2;
            P(s,2) = 3*Ppdm/4; P(s,6) = 3*Ppdp/4;
        end
    elseif i == 1 % western boundary
        P(s,3) = Ppv/2; P(s,5) = 2*Ppdm; P(s,6) = 2*Ppdp; P(s,4) = Ppv/2;
    elseif i == Nx % eastern boundary
        P(s,3) = Ppv/2; P(s,1) = 2*Ppdp; P(s,2) = 2*Ppdm; P(s,4) = Ppv/2;
    else % interior
        P(s,3) = Ppv/2; P(s,1) = Ppdp; P(s,5) = Ppdm; P(s,2) = Ppdm;
        P(s,6) = Ppdp; P(s,4) = Ppv/2;
    end

end

% Define cumulative probabilities
for i = 2:6
    P(:,i) = P(:,i) + P(:,i-1);
end

end