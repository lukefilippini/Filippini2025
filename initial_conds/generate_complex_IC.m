function f = generate_complex_IC(xn,yn,Nf,dens,mesh_config)
% GENERATE_COMPLEX_IC generates a complex initial condition of constant 
% density for the deterministic and stochastic frameworks. 
% Inputs:
%   x - x-component of nodes (vector)
%   y - y-component of nodes (vector)
%   Nf - no. of aggregration iterations (scalar)
%   dens - approx. proportion of non-zero density (scalar)
%   mesh_config - mesh configuration (string)
%       Options: 'rectangular', 'flat-top' or 'pointy-top'
% Output:
%   f - initial condition (vector)

x = unique(xn); Nx = length(x); % no. of nodes in x direction
y = unique(yn); Ny = length(y); % no. of nodes in y direction
f = rand(size(xn)); % initial densities

map = @(i,j) map_func(i,j,Nx,Ny,mesh_config); % mapping function

for k = 1:Nf
    f = aggregrate_density(f,Nx,Ny,map,mesh_config);
end

th = (1-dens)*max(f) + dens*min(f); % density threshold
f(f < th) = 0; f(f >= th) = 1; % keep or remove initial density

end

%% Subfunctions
function k = map_func(i,j,Nx,Ny,mesh_config)
% MAP_FUNC defines the mapping function relating node indices i,j to
% position k in the initial condition f. Periodic boundaries are used for
% attempts to index beyond domain boundaries.
% Inputs:
%   i - x-index / column (scalar)
%   j - y-indx / row (scalar)
%   Nx - no. of nodes in the x-direction (scalar)
%   Ny - no. of nodes in the y-direction (scalar)
%   mesh_config - mesh configuration (string)
%       Options: 'rectangular', 'flat-top' or 'pointy-top' 

if isequal(mesh_config,'rectangular')

    if i == 0 % western to eastern
        i = Nx;
    elseif i == Nx+1 % eastern to western
        i = 1;
    end

    if j == 0 % southern to northern
        j = Ny;
    elseif j == Ny+1 % northern to southern
        j = 1;
    end

    k = (j-1)*Nx + i; % mapping function

elseif isequal(mesh_config,'flat-top')

    if mod(j,2) == 0 % convex western and concave eastern
        if i == 0 % convex western (wrap left node)
            i = Nx;
        elseif i == Nx+1 % concave eastern (wrap diagonal nodes)
            i = 1;
        elseif i == Nx+2 % concave eastern (wrap right node)
            i = 2;
        end
    else % concave western and convex eastern
        if i == -1 % concave western (wrap diagonal nodes)
            i = Nx-1;
        elseif i == 0 % concave western (wrap left node)
            i = Nx;
        elseif i == Nx+1 % convex eastern (wrap right node)
            i = 1;
        end
    end

    if j == 0 % southern to northern
        j = Ny;
    elseif j == Ny+1 % northern to southern
        j = 1;
    end

    k = (j-1)*Nx/2 + ceil(i/2);

else % pointy-top

    if mod(i,2) == 0 % convex southern and concave northern
        if j == 0 % convex southern (wrap bottom node)
            j = Ny;
        elseif j == Ny+1 % concave northern (wrap diagonal nodes)
            j = 1;
        elseif j == Ny+2 % concave northern (wrap top node)
            j = 2;
        end
    else % concave southern and convex northern
        if j == -1 % concave southern (wrap bottom node)
            j = Ny-1;
        elseif j == 0 % concave southern (wrap diagonal nodes)
            j = Ny;
        elseif j == Ny+1 % convex northern (wrap top node)
            j = 1;
        end
    end

    if i == 0 % western to eastern
        i = Nx;
    elseif i == Nx+1 % eastern to western
        i = 1;
    end

    k = (i-1)*Ny/2 + ceil(j/2);

end
    
end

function f = aggregrate_density(f,Nx,Ny,map,mesh_config)
% AGGREGATE_DENSITY aggregates density in the initial condition f. The
% density at f(map(i,j)) is modified according to a weighting of densities
% at the current and surrounding nodes. These weightings are inspired by
% those from the lattice Boltzmann method.
% Inputs:
%   Nx - no. of nodes in the x-direction (scalar)
%   Ny - no. of nodes in the y-direction (scalar)
%   map - mapping function (function handle)
%   mesh_config - mesh configuration (string)
%       Options: 'rec', 'flat-top' or 'pointy-top'

if isequal(mesh_config,'rectangular')
    for i = 1:Nx
        for j = 1:Ny
            f(map(i,j)) = 4/9*f(map(i,j)) + 1/9*(f(map(i-1,j)) + ...
                f(map(i+1,j))) + 1/9*(f(map(i,j-1)) + f(map(i,j+1))) + ...
                1/36*(f(map(i-1,j-1)) + f(map(i-1,j+1)) + ...
                f(map(i+1,j-1)) + f(map(i+1,j+1)));
        end
    end
elseif isequal(mesh_config,'flat-top')
    for i = 1:Nx
        for j = 1:Ny
            if mod(i+j,2) == 0 % sum of indices must be even
                f(map(i,j)) = 1/2*f(map(i,j)) + 1/12 * (f(map(i-1,j-1)) + ...
                    f(map(i+1,j-1)) + f(map(i-2,j)) + f(map(i+2,j)) + ...
                    f(map(i-1,j+1)) + f(map(i+1,j+1)));
            end
        end
    end
else % pointy-top
    for i = 1:Nx
        for j = 1:Ny
            if mod(i+j,2) == 0 % sum of indices must be even
                f(map(i,j)) = 1/2*f(map(i,j)) + 1/12 * (f(map(i-1,j-1)) + ...
                    f(map(i+1,j-1)) + f(map(i,j-2)) + f(map(i,j+2)) + ...
                    f(map(i-1,j+1)) + f(map(i+1,j+1)));
            end
        end
    end
end

end