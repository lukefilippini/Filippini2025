function f = generate_simple_IC(xn,yn,Lx,Ly,px,py)
% GENERATE_SIMPLE_IC generates an initial condition consisting of a 
% rectangle of constant density for the deterministic and 
% stochastic frameworks.
% Inputs:
%   x - x-component of nodes (vector)
%   y - y-component of nodes (vector)
%   Lx - length of domain in x-direction (scalar)
%   Ly - length of domain in y-direction (scalar)
%   px - percentage length of rectangle (scalar)
%   py - percentage width of rectangle (scalar)
% Output:
%   f - initial condition (vector)

n_nodes = length(xn);
f = zeros(n_nodes,1);

for n = 1:n_nodes
    % Get current node
    x = xn(n); y = yn(n);

    % Define rectangle
    xb = Lx/2 + [-1,1]*px*Lx/2;
    yb = Ly/2 + [-1,1]*py*Ly/2;
    
    % If node is within rectangle, save in IC
    if (x >= xb(1)) && (x <= xb(2))
        if (y >= yb(1)) && (y <= yb(2))
            f(n) = 1;
        end
    end
end
