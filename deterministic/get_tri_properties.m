function mesh_props = get_tri_properties(nodes,elements)
% GET_TRI_PROPERTIES computes mesh properties associated with
% each element of a two-dimensional finite volume scheme, using an
% flat-top or pointy-top hexagonal mesh (triangular elements) and no-flux 
% boundary conditions.
% Inputs:
%   nodes - location of each node (matrix : x - first col, y - second col)
%   elements - nodes associated with each element (anti-clockwise from
%       bottom-left; matrix: each row contains four node numbers)
% Ouputs:
%   mesh_props - struct containing mesh properties
%       mesh_props.nodes - location of each node (matrix)
%       mesh_props.elements - nodes associated with each element (matrix)
%       mesh_props.no_nodes - number of nodes (scalar)
%       mesh_props.no_elements - number of elements (scalar)
%       mesh_props.mx - x component of control volume face midpoints 
%       (matrix: three columns; one for each midpoint)
%       mesh_props.my - y component of control volume face midpoints 
%       (matrix: three columns; one for each midpoint)
%       mesh_props.nx - x component of unit normal to control volume faces
%       (matrix: three columns; one for each control volume face)
%       mesh_props.ny - y component of unit normal to control volume faces
%       (matrix: three columns; one for each control volume face)
%       mesh_props.s - shape function coefficients (matrix: 16 columns)
%       mesh_props.V - control volume areas (vector: no_nodes x 1)
%       mesh_props.E - element areas (vector: no_elements x 1)
%       mesh_props.L - control volume edge lengths (matrix: three columns; 
%       one for each control volume face in each element)
%       mesh_props.c - element centroids (vector: no_elements x 1)
%       mesh_props.domain_area - unit domain area (scalar)

% Setup
no_nodes = size(nodes,1); % no. of nodes
no_elements = size(elements,1); % no. of elements

% Initialise
mx = zeros(no_elements,3); % x component of CV face midpoints
my = zeros(no_elements,3); % y component of CV face midpoints
nx = zeros(no_elements,3); % x component of unit normals
ny = zeros(no_elements,3); % y component of unit normals
s = zeros(no_elements,9); % shape function coefficients
V = zeros(no_nodes,1); % control volume areas
E = zeros(no_elements,1); % element areas
L = zeros(no_elements,3); % control volume edge lengths (in each element)
c = zeros(no_elements,2); % centroid of each element
domain_area = 0; % unit domain area

% Generate mesh properties
for k = 1:no_elements

    % Setup
    v = elements(k,:); % element vertices (node numbers)
    xv = nodes(v,1); % x component of vertices
    yv = nodes(v,2); % y component of vertices
    
    % Element centroid
    c(k,:) = [mean(xv) mean(yv)];

    % Element edge midpoints
    e1 = [(xv(1) + xv(2))/2 (yv(1) + yv(2))/2];
    e2 = [(xv(2) + xv(3))/2 (yv(2) + yv(3))/2];
    e3 = [(xv(3) + xv(1))/2 (yv(3) + yv(1))/2];

    % Control volume edge vectors
    u1 = c(k,:) - e1;
    u2 = c(k,:) - e2;
    u3 = c(k,:) - e3;

    % Control volume edge lengths
    L(k,1) = norm(u1,2);
    L(k,2) = norm(u2,2);
    L(k,3) = norm(u3,2);

    % Midpoints of control volume edges
    mx(k,1) = (c(k,1) + e1(1))/2; my(k,1) = (c(k,2) + e1(2))/2;
    mx(k,2) = (c(k,1) + e2(1))/2; my(k,2) = (c(k,2) + e2(2))/2;
    mx(k,3) = (c(k,1) + e3(1))/2; my(k,3) = (c(k,2) + e3(2))/2;

    % Unit normals
    nx(k,1) = u1(2)/L(k,1); ny(k,1) = -u1(1)/L(k,1);
    nx(k,2) = u2(2)/L(k,2); ny(k,2) = -u2(1)/L(k,2);
    nx(k,3) = u3(2)/L(k,3); ny(k,3) = -u3(1)/L(k,3);

    % Vectors connecting nodes to centroid
    p1 = c(k,:) - [xv(1) yv(1)];
    p2 = c(k,:) - [xv(2) yv(2)];
    p3 = c(k,:) - [xv(3) yv(3)];

    % Vectors connecting element edge midpoints
    q1 = e1 - e3; q2 = e2 - e1; q3 = e3 - e2;
    
    % Sub-control volume areas
    S1 = (q1(1)*p1(2) - q1(2)*p1(1))/2;
    S2 = (q2(1)*p2(2) - q2(2)*p2(1))/2;
    S3 = (q3(1)*p3(2) - q3(2)*p3(1))/2;
    
    % Control volume areas
    V(v(1)) = V(v(1)) + S1;
    V(v(2)) = V(v(2)) + S2;
    V(v(3)) = V(v(3)) + S3;

    % Element area
    E(k) = S1 + S2 + S3;
    
    % Unit domain area
    domain_area = domain_area + E(k);

    % Shape function coefficients
    A = [xv yv ones(3,1)]; detA = det(A); I = eye(3); % determinant and identity
    m = 0; % counter
    for i = 1:3 % for each main coefficient (alpha, beta, etc.)
        for j = 1:3 % for each shape coefficient
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
mesh_props.nx = nx;
mesh_props.ny = ny;
mesh_props.s = s;
mesh_props.V = V;
mesh_props.E = E;
mesh_props.L = L;
mesh_props.c = c;
mesh_props.domain_area = domain_area;

end
