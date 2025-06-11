function c = compute_numerical_sol(mesh_props,D,T,M,f,mesh_config)
% COMPUTE_NUMERICAL_SOL computes a numerical solution to the diffusion
% equation in two dimensions for a single diffusion tensor (D) and no-flux 
% boundary conditions. An element-based vertex-centred finite volume method 
% is used for the spatial discretisation with forward Euler time-stepping 
% for the temporal discretisation.
% Inputs:
%   mesh_props - mesh properties (struct - see 
%       get_rec_properties.m or get_tri_properties.m)
%   D - diffusion tensor (matrix)
%   T - final time (scalar)
%   M - no. of time steps (scalar)
%   f - initial solution (vector)
%   mesh_config - mesh configuration (string)
%       Options: 'rectangular', 'flat-top' or 'pointy-top'
% Outputs:
%   c - numerical solution (matrix or vector)

% Setup
no_nodes = mesh_props.no_nodes; % number of nodes
no_elements = mesh_props.no_elements; % number of elements
nodes = mesh_props.nodes; % mesh nodes
elements = mesh_props.elements; % vertices of elements (anti-clockwise from bottom-left)
mx = mesh_props.mx; % x-component of control volume edge midpoints
my = mesh_props.my; % y-component of control volume edge midpoints
s = mesh_props.s; % shape function coefficients
V = mesh_props.V; % control volume areas
L = mesh_props.L; % control volume edge lengths

if ~isequal(mesh_config,'rectangular')
    nx = mesh_props.nx; % x component of unit normals
    ny = mesh_props.ny; % y component of unit normals
end

Nx = length(unique(nodes(:,1)));
Ny = length(unique(nodes(:,2)));
tau = T / M; % time step duration

% Spatial discretisation %
A = zeros(no_nodes,no_nodes); % coefficient matrix
for k = 1:no_elements

    v = elements(k,:); % vertices of element

    if isequal(mesh_config,'rectangular')
        
        % Determine nodes associated with flux term contribution
        for j = 1:4 % for each control volume face midpoint (anti-clockwise direction)
            if j == 1
                jnb = 2; n = 1;
                D1 = D(1,1); D2 = D(1,2);
            elseif j == 2
                jnb = 3; n = 1;
                D1 = D(1,2); D2 = D(2,2);
            elseif j == 3
                jnb = 4; n = -1;
                D1 = D(1,1); D2 = D(1,2);
            else
                jnb = 1; n = -1;
                D1 = D(1,2); D2 = D(2,2);
            end
    
            % Coefficients of unknowns in flux terms
            F = zeros(4,1);
            F(1) = n*(D1*(s(k,1) + my(k,j)*s(k,9)) + ... 
                D2*(s(k,5) + mx(k,j)*s(k,9)))*L(k,j);
            F(2) = n*(D1*(s(k,2) + my(k,j)*s(k,10)) + ... 
                D2*(s(k,6) + mx(k,j)*s(k,10)))*L(k,j);
            F(3) = n*(D1*(s(k,3) + my(k,j)*s(k,11)) + ... 
                D2*(s(k,7) + mx(k,j)*s(k,11)))*L(k,j);
            F(4) = n*(D1*(s(k,4) + my(k,j)*s(k,12)) + ... 
                D2*(s(k,8) + mx(k,j)*s(k,12)))*L(k,j);
        
            % Distribute into coefficient matrix
            for i = 1:4
                A(v(j),v(i)) = A(v(j),v(i)) + F(i) / V(v(j));
                A(v(jnb),v(i)) = A(v(jnb),v(i)) - F(i) / V(v(jnb));
            end
        end
    else % hexagonal configurations
        for j = 1:3 % for each control volume face midpoint (anti-clockwise direction)
            if j < 3
                jnb = j + 1;
            else
                jnb = 1;
            end

            Dxx = D(1,1); Dxy = D(1,2); Dyy = D(2,2);
            
            % Coefficients of unknowns in flux terms
            F = zeros(3,1); % store flux terms
            F(1) = (Dxx*s(k,1) + Dxy*s(k,4))*nx(k,j)*L(k,j) + ...
                (Dxy*s(k,1) + Dyy*s(k,4))*ny(k,j)*L(k,j);
            F(2) = (Dxx*s(k,2) + Dxy*s(k,5))*nx(k,j)*L(k,j) + ...
                (Dxy*s(k,2) + Dyy*s(k,5))*ny(k,j)*L(k,j);
            F(3) = (Dxx*s(k,3) + Dxy*s(k,6))*nx(k,j)*L(k,j) + ...
                (Dxy*s(k,3) + Dyy*s(k,6))*ny(k,j)*L(k,j);

            % Distribute into coefficient matrix
            for i = 1:3
                A(v(j),v(i)) = A(v(j),v(i)) + F(i) / V(v(j));
                A(v(jnb),v(i)) = A(v(jnb),v(i)) - F(i) / V(v(jnb));
            end
        end
    end
end

% Convert matrix to sparse format %
A = sparse(A);

% Time stepping %
c = zeros(no_nodes,M+1); % initialise solution
c(:,1) = f; % initial solution
I = speye(no_nodes); % identity matrix
At = I + tau*A; % RHS of forward Euler method

fprintf('Deterministic progress: ')
k = 0; % counter
for n = 1:M
    if mod(n,floor(M/10)) == 0
        k = k + 1;
        fprintf('%g%% ',k*10)
    end
    c(:,n+1) = At * c(:,n);
end
fprintf('\n\n')

if isequal(mesh_config,'rectangular')
    c = reshape(c,Nx,Ny,M+1);
end

end


