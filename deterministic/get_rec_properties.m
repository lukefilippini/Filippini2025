function mesh_props = get_rec_properties(nodes,elements)
% GET_REC_PROPERTIES computes mesh properties associated with
% each element of a two-dimensional finite volume scheme, using a
% rectangular mesh configuration and no-flux boundary conditions.
% Inputs:
%   nodes - location of each node (matrix : x - first col, y - second col)
%   elements - nodes associated with each element (anti-clockwise from
%       bottom-left; matrix: each row contains an index for each node)
% Ouputs:
%   mesh_props - struct containing mesh properties
%       mesh_props.nodes - location of each node (matrix)
%       mesh_props.elements - nodes associated with each element (matrix)
%       mesh_props.no_nodes - number of nodes (scalar)
%       mesh_props.no_elements - number of elements (scalar)
%       mesh_props.mx - x-component of control volume edge midpoints 
%       (matrix: four columns; one for each midpoint)
%       mesh_props.my - y-component of control volume edge midpoints 
%       (matrix: four columns; one for each midpoint)
%       mesh_props.s - shape function coefficients (matrix: 16 columns)
%       mesh_props.V - control volume areas (vector: no_nodes x 1)
%       mesh_props.L - control volume edge lengths (matrix: four columns; 
%       one for each control volume edge in each element)
%       mesh_props.c - element centroids (vector: no_elements x 1)

% Setup
no_nodes = size(nodes,1); % no. of nodes
no_elements = size(elements,1); % no. of elements

% Initialise
mx = zeros(no_elements,4); % x-component of CV edge midpoints
my = zeros(no_elements,4); % y-component of CV edge midpoints
s = zeros(no_elements,16); % shape function coefficients
V = zeros(no_nodes,1); % control volume areas
L = zeros(no_elements,4); % control volume edge lengths (in each element)
c = zeros(no_elements,2); % centroid of each element

% Generate mesh properties
for k = 1:no_elements

    % Setup
    v = elements(k,:); % element vertices (node numbers)
    xv = nodes(v,1); % x-component of vertices
    yv = nodes(v,2); % y-component of vertices
    
    % Element centroid
    c(k,:) = [mean(xv) mean(yv)];
    
    % Control volume edge lengths
    L(k,1) = c(k,2) - yv(1); 
    L(k,2) = xv(2) - c(k,1);
    L(k,3) = yv(3) - c(k,2);
    L(k,4) = c(k,1) - xv(4);

    % Midpoints of control volume edges
    mx(k,1) = c(k,1); my(k,1) = (yv(1) + c(k,2))/2;
    mx(k,2) = (xv(2) + c(k,1))/2; my(k,2) = c(k,2);
    mx(k,3) = c(k,1); my(k,3) = (yv(3) + c(k,2))/2;
    mx(k,4) = (c(k,1) + xv(4))/2; my(k,4) = c(k,2);

    % Control volume areas
    S = (xv(2) - xv(1))*(yv(3) - yv(2))/4; % sub-control volume area
    V(v) = V(v) + S; % add sub-area to total for each node
    
    % Shape function coefficients
    A = [xv yv xv.*yv ones(4,1)]; detA = det(A); I = eye(4); % determinant and identity
    m = 0; % counter
    for i = 1:4 % for each main coefficient (alpha, beta, etc.)
        for j = 1:4 % for each shape coefficient
            m = m + 1;
            A1 = A; A1(:,i) = I(:,j); % replace columns of A
            s(k,m) = det(A1)/detA; % coefficients
        end
    end
end

% Store mesh properties
mesh_props.nodes = nodes;
mesh_props.elements = elements;
mesh_props.no_nodes = no_nodes;
mesh_props.no_elements = no_elements;
mesh_props.mx = mx;
mesh_props.my = my;
mesh_props.s = s;
mesh_props.V = V;
mesh_props.L = L;
mesh_props.c = c;

end
