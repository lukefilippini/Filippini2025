function mesh_props = generate_mesh(Lx,Ly,Ex,Ey,mesh_config)
% GENERATE_MESH generates mesh properties for a structured mesh 
% partitioning the domain (0,Lx) x (0,Ly). The mesh configurations are
% rectangular or hexagonal (flat-top or pointy_top). 
% Inputs:
%   Lx - domain length (scalar)
%   Ly - domain width (scalar)
%   Ex - no. of elements in x direction (scalar)
%   Ey - no. of elements in y direction (scalar)
%   mesh_config - mesh configuration (string)
%       Options: 'rectangular', 'flat-top' or 'pointy-top'
% Outputs:
%   mesh_props - mesh properties (struct - see
%       get_rec_properties.m or get_tri_properties.m)

x = linspace(0,Lx,Ex+1); y = linspace(0,Ly,Ey+1); % node components
nodes = zeros((Ex+1)*(Ey+1),2); % initialise node matrix

if isequal(mesh_config,'rectangular')
    for j = 1:Ey+1
        for i = 1:Ex+1
            indx = (Ex+1)*(j-1) + i; % map node to matrix position
            nodes(indx,1) = x(i);
            nodes(indx,2) = y(j);
        end
    end
    
    elements = zeros(Ex*Ey,4); % initialise element matrix
    vert = [0 1 (Ex+1)+1 (Ex+1)]; % vertices of element (anti-clockwise from bottom-left)
    k = 0; % row of element matrix
    for j = 1:Ey
        for i = 1:Ex
            k = k + 1; % increase row index
            vert = vert + 1; % move horizontally to next element
            elements(k,:) = vert; % store vertices
        end
        vert = vert + 1; % move to next mesh row
    end

    % Generate mesh properties
    mesh_props = get_rec_properties(nodes,elements);

elseif isequal(mesh_config,'flat-top')
    
    delta_x = x(2) - x(1); % spatial step in x direction
           
    for j = 1:Ey+1
        if mod(j,2) == 0 
            xt = x + delta_x/2; % shift every second row
        elseif mod(j-1,2) == 0 
            xt = x; % odd-numbered rows stay the same
        end
    
        for i = 1:Ex+1
            indx = (Ex+1)*(j-1) + i; % map node to matrix position
            nodes(indx,1) = xt(i);
            nodes(indx,2) = y(j);
        end
    end
    
    elements = zeros(2*Ex*Ey,3); % initialise element matrix
    vert_up = [0 1 Ex+1]; % element vertices (upward triangles)
    vert_down = [1 Ex+2 Ex+1]; % element vertices (downward triangles)
    k = 0; % row of element matrix
    for j = 1:Ey
        for i = 1:Ex % upward triangular elements
            k = k + 1;
            vert_up = vert_up + 1;
            elements(k,:) = vert_up;
        end

        for i = 1:Ex % downward triangular elements
            k = k + 1;
            vert_down = vert_down + 1;
            elements(k,:) = vert_down;
        end
        
        % Increase index as normal
        vert_up(1:2) = vert_up(1:2) + 1;
        vert_down(2:3) = vert_down(2:3) + 1;
    
        % Increase index depending on mesh row
        if mod(j,2) == 0
            vert_down(1) = vert_down(1) + 2;        
        else
            vert_up(3) = vert_up(3) + 2;
        end

    end

    % Generate mesh properties
    mesh_props = get_tri_properties(nodes,elements);

else % pointy-top

    delta_y = y(2) - y(1); % spatial step in y direction
           
    for i = 1:Ex+1
        if mod(i,2) == 0 
            yt = y + delta_y/2; % shift every second row
        elseif mod(i-1,2) == 0 
            yt = y; % odd-numbered rows stay the same
        end
    
        for j = 1:Ey+1
            indx = (Ey+1)*(i-1) + j; % map node to matrix position
            nodes(indx,1) = x(i);
            nodes(indx,2) = yt(j);
        end
    end
    
    elements = zeros(2*Ex*Ey,3); % initialise element matrix
    vert_right = [0 1 Ey+1]; % element vertices (triangles pointing to right)
    vert_left = [1 Ey+2 Ey+1]; % element vertices (triangles pointing to left)
    k = 0; % row of element matrix
    for i = 1:Ex
        for j = 1:Ey % right-facing triangular elements
            k = k + 1;
            vert_right = vert_right + 1;
            elements(k,:) = vert_right;
        end

        for j = 1:Ey % downward triangular elements
            k = k + 1;
            vert_left = vert_left + 1;
            elements(k,:) = vert_left;
        end
        
        % Increase index as normal
        vert_right(1:2) = vert_right(1:2) + 1;
        vert_left(2:3) = vert_left(2:3) + 1;
    
        % Increase index depending on mesh row
        if mod(i,2) == 0
            vert_left(1) = vert_left(1) + 2;        
        else
            vert_right(3) = vert_right(3) + 2;
        end

    end
    
    % Generate mesh properties
    mesh_props = get_tri_properties(nodes,elements);

end