function [Ex,Ey] = generate_no_of_elements(Lx,Ly,E,D,mesh_config,res)
% GENERATE_NO_OF_ELEMENTS generates the number of elements in the x (Ex)
% and y (Ey) directions for the mesh configuration 'mesh_config'. Ex and Ey
% are defined to ensure all transition probabilities are non-negative.
% Inputs:
%   Lx - domain length in the x direction (scalar)
%   Ly - domain length in the y direction (scalar)
%   E - minimum number of elements in x or y direction (scalar)
%   D - diffusion tensor (matrix)
%   mesh_config - mesh configuration (string)
%       Options: 'rec', 'flat-top' or 'pointy-top'
%   res - mesh resolution (string)
%       Options: 'low', 'med' or 'high'
% Outputs:
%   Ex - no. of elements in the x direction (scalar)
%   Ey - no. of elements in the y direction (scalar)

Dxx = D(1,1); Dyy = D(2,2); Dxy = D(1,2);

% Generate no. of elements
if isequal(mesh_config,'rectangular')
    % Define no. of elements in x-direction
    Ex = max([ceil(Lx / (Ly * (min(sqrt(3*Dyy/Dxx),Dxx/(2*abs(Dxy))) - ...
        max(sqrt(Dyy/(3*Dxx)),2*abs(Dxy)/Dyy)))),E]);
    
    % Define no. of elements in y-direction (three values)
    min_Ey = max(2*Ex*Ly*abs(Dxy) / (Lx*Dyy),Ex*Ly/Lx * sqrt(Dxx/(3*Dyy)));
    max_Ey = min(Ex*Ly*Dxx / (2*Lx*abs(Dxy)),Ex*Ly/Lx * sqrt(3*Dxx/Dyy));

    if isequal(res,'low')
        Ey = ceil(min_Ey);
    elseif isequal(res,'med')
        Ey = round((min_Ey + max_Ey)/2); 
    elseif isequal(res,'high')
        Ey = floor(max_Ey);
    else
        error("Invalid resolution. 'res' must be 'low', 'med' or 'high'.")
    end
        
    
elseif isequal(mesh_config,'flat-top')
    % Define no. of elements in x-direction
    Ex = max(ceil( Lx*Dyy / (2*Ly*(sqrt(Dxx*Dyy) - abs(Dxy))) ),E);
    
    % Define no. of elements in y-direction (three values)
    min_Ey = 2*Ex*Ly*abs(Dxy) / (Lx*Dyy);
    max_Ey = 2*Ex*Ly/Lx * sqrt(Dxx/Dyy);
    
    if isequal(res,'low')
        Ey = ceil(min_Ey);
        if mod(Ey,2) ~= 0 % must be even
            Ey = Ey + 1;
        end
    elseif isequal(res,'med')
        Ey = round((min_Ey + max_Ey)/2); 
        if mod(Ey,2) ~= 0 % must be even
            Ey = Ey + 1;
        end
    elseif isequal(res,'high')
        Ey = floor(max_Ey);
        if mod(Ey,2) ~= 0 % must be even
            Ey = Ey - 1;
        end
    else
        error("Invalid resolution. 'res' must be 'low', 'med' or 'high'.")
    end

else % pointy-top
    % Define no. of elements in y-direction
    Ey = max(ceil( Ly*Dxx / (2*Lx*(sqrt(Dxx*Dyy) - abs(Dxy))) ),E);
    
    % Define no. of elements in x-direction (three values)
    min_Ex = 2*Ey*Lx*abs(Dxy) / (Ly*Dxx);
    max_Ex = 2*Ey*Lx/Ly * sqrt(Dyy/Dxx);
    
    if isequal(res,'low')
        Ex = ceil(min_Ex);
        if mod(Ex,2) ~= 0 % must be even
            Ex = Ex + 1;
        end
    elseif isequal(res,'med')
        Ex = round((min_Ex + max_Ex)/2); 
        if mod(Ex,2) ~= 0 % must be even
            Ex = Ex + 1;
        end
    elseif isequal(res,'high')
        Ex = floor(max_Ex);
        if mod(Ex,2) ~= 0 % must be even
            Ex = Ex - 1;
        end
    else
        error("Invalid resolution. 'res' must be 'low', 'med' or 'high'.")
    end
end

% fprintf('Dxx = %g, Dxy = %g, Dyy = %g \n\n',Dxx,Dxy,Dyy);
% fprintf("No. of elements: Ex = %g, Ey = %g \n\n",Ex,Ey)

end