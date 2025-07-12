function f = generate_simple_IC(xn,yn,Lx,Ly)
% GENERATE_SIMPLE_IC generates a simple initial condition consisting of a 
% rectangle of constant density for the deterministic and stochastic 
% frameworks. The initial condition is defined as f(x) = 1 for 3Lx/10 <= x
% <= 7Ly/10 and 3Ly/10 <= y <= 7Ly/10 and f(x) = 0 otherwise.
% Inputs:
%   x - x-component of nodes (vector)
%   y - y-component of nodes (vector)
%   Lx - length of domain in x-direction (scalar)
%   Ly - length of domain in y-direction (scalar)
% Output:
%   f - initial condition (vector)

n_nodes = length(xn);
f = zeros(n_nodes,1);

for n = 1:n_nodes
    % Get current node
    x = xn(n); y = yn(n);

    % Define rectangle
    xb = Lx/2 + [-1,1]*Lx/5;
    yb = Ly/2 + [-1,1]*Ly/5;
    
    % If node is within rectangle, save in IC
    if (x >= xb(1)) && (x <= xb(2))
        if (y >= yb(1)) && (y <= yb(2))
            f(n) = 1;
        end
    end
end
