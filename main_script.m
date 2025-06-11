%% Filippini et al. 2025 %%
clc, clear, close all
% delete("figures") % remove existing figures

% Include subfolders in directory
folder = fileparts(which("main_script.m"));
addpath(genpath(folder));

% Domain parameters
Lx = 1; Ly = 0.8; % length in x and y directions
E = 61; % minimum no. of elements in x or y direction
% res = 'low'; % mesh resolution
% res = 'med';
res = 'high';

% mesh_config = 'rectangular'; % mesh configuration
mesh_config = 'flat-top';
% mesh_config = 'pointy-top';

% Time parameters
T = 1; % final time
t = [T/32 T/8 T/4 T/2]; % plotting times
Nt = length(t); % no. of plotting times

% Stochastic parameters
Ns = 5; % no. of simulations

% Diffusion tensors 
lambda_x = 0.01; lambda_y = 0.09; % tensor eigenvalues
angles = [pi/24 pi/12 pi/6]; % rotation angles
Nr = length(angles); % no. of rotation angles

% Initial condition
px = 0.4; py = 0.4; % proportion length and width of rectangle (simple)
dens = 0.5; % approx. proportion of domain with non-zero density (complex)
Nf = 20; % no. of aggregrations (complex)
IC_type = 'simple'; % type of initial condition
% IC_type = 'complex'; 

% Visualisation of results
visualise = true;

if visualise
    config_dir = [mesh_config,'_',IC_type,'_',res];
    mkdir("figures",config_dir)
end

colours = [128,0,0; 245,130,48; 0 0 128; 0 130 200]/255;

% Table of quantitative metrics and parameters
output_metrics = true;
output_table = true;

if output_table
    metrics_table = cell(Nr*7,10);
    table_headings = {'epsilon','MSE','theta','Dxx','Dxy','Dyy','Ex', ...
        'Ey','tau(M)','Np'};
    for i = 1:10
        metrics_table{1,i} = table_headings{i};
    end
    m = 1; % counter
end


% Compute for each rotation angle / diffusion tensor
fprintf('Eigenvalues: lambda_x = %g, lambda_y = %g\n',lambda_x,lambda_y)
fprintf('Initial condition: %s\n',IC_type)
for i = 1:Nr
    if visualise
        angle_dir = ['angle_',num2str(i)];
        mkdir(strcat("./figures/",config_dir),angle_dir)
    end

    fprintf('Rotation angle: %g\n',angles(i))

    % Generate diffusion tensor
    theta = angles(i); % rotation angle
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % rotation matrix
    D = R*[lambda_x 0; 0 lambda_y]*R'; % diffusion tensor
    fprintf('Tensor components: Dxx = %g, Dxy = %g, Dyy = %g\n',D(1,1),D(1,2),D(2,2))
    
    % Check rectangular configuration viability
    sb = 2 / abs(lambda_x - lambda_y) * sqrt(lambda_x*lambda_y/3);
    if sb < 1 && isequal(mesh_config,'rectangular')
        theta_crit = 1/2*asin(sb);
        if theta > theta_crit
            fprintf("Rectangular configuration unsuitable for rotation angle %g. \n",theta)
            break
        end
    end

    % Generate no. of elements
    [Ex,Ey] = generate_no_of_elements(Lx,Ly,E,D,mesh_config,res);
    dx = Lx / Ex; dy = Ly / Ey; % spatial steps

    % Generate mesh properties and no. of time steps
    fprintf('---------------------------------------\n')
    fprintf('No. of lattice sites: Ex = %g, Ey = %g \n',Ex,Ey)
    mesh_props = generate_mesh(Lx,Ly,Ex,Ey,mesh_config);
    M = generate_time_steps(dx,dy,T,D,mesh_config);
    
    % Extract properties
    nodes = mesh_props.nodes; xn = nodes(:,1); yn = nodes(:,2); % nodes
    x = unique(xn,'stable'); y = unique(yn,'stable'); % unique nodes
    Nx = length(x); Ny = length(y); % no. of nodes in x and y directions
    elements = mesh_props.elements; % elements
    V = mesh_props.V; % control volumes areas
    tau = T / M; % time step duration

    % Define initial condition
    if isequal(IC_type,'simple')
        f = generate_simple_IC(xn,yn,Lx,Ly,px,py);
    else
        f = generate_complex_IC(xn,yn,Nf,dens,mesh_config);
    end

    % Compute numerical solution
    c = compute_numerical_sol(mesh_props,D,T,M,f,mesh_config);
    
    % Define number of particles
    Np = [25 100 250]*(Ex+1)*(Ey+1); % no. of desired particles
    Nps = length(Np); % no. of particle groups

    % if t(1) == 0
    %     Np_dir = ['Np_0'];
    %     mkdir(strcat('./figures/',config_dir,'/',angle_dir),Np_dir)
    %     path_name = strcat('./figures/',config_dir,'/',angle_dir,'/',Np_dir,'/');
    %     % Deterministic framework
    %     if isequal(mesh_config,'rectangular')
    %         surf(X,Y,c(:,:,1)')
    %         axis([0, Lx, 0, Ly, 0, max(f)])
    %     else % hexagonal configurations
    %         trimesh(elements,xn,yn,c(:,1),'FaceColor','interp','EdgeColor','black')
    % 
    %         if isequal(mesh_config,'flat-top')
    %             axis([0, Lx+dx/2, 0, Ly, 0, max(f)])
    %         else
    %             axis([0, Lx, 0, Ly+dy/2, 0, max(f)])
    %         end
    %     end
    % 
    %     xlab = xlabel('$x$','Interpreter','latex','FontSize',30);
    %     set(xlab, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);
    % 
    %     ylab = ylabel('$y$','Interpreter','latex','FontSize',30);
    %     set(ylab, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    % 
    %     set(gca,'XTick',[0,Lx],'XTickLabel',{'$0$','$L_{x}$'}, ...
    %         'YTick',[0,Ly],'YTickLabel',{'$0$','$L_{y}$'}, ...
    %         'TickLabelInterpreter','latex','FontSize',30)
    % 
    %     shading interp
    %     axis equal
    %     view(0,90)
    %     colormap(cork)
    %     clim([min(f),max(f)])
    % 
    %     print(gcf,[path_name,['2D_',config_dir,'_',angle_dir, ...
    %         '_',Np_dir,'_t0']],'-depsc2')
    % end

    % Run stochastic model for each desired no. of particles
    for j = 1:Nps
        
        % Run stochastic simulations
        if isequal(mesh_config,'rectangular')
            [X,Y] = meshgrid(x,y);
            [Uavg,Ul,Uu,Npn] = rectangular_stochastic(xn,yn,dx,dy,V,D,T,M,f,Np(j),Ns);
        elseif isequal(mesh_config,'flat-top')
            [Uavg,Ul,Uu,Npn] = flat_top_stochastic(xn,yn,dx,dy,V,D,T,M,f,Np(j),Ns);
        else
            [Uavg,Ul,Uu,Npn] = pointy_top_stochastic(xn,yn,dx,dy,V,D,T,M,f,Np(j),Ns);
        end
        
        % Visualisation of results
        if visualise

            Np_dir = ['Np_',num2str(j)];
            mkdir(strcat('./figures/',config_dir,'/',angle_dir),Np_dir)
            path_name = strcat('./figures/',config_dir,'/',angle_dir,'/',Np_dir,'/');
            
            % Two-dimensional comparisons 
            indx = floor(t / tau + 1);
            for n = 1:Nt
                figure
                tiledlayout(1,2,'Padding','compact','TileSpacing','compact')

                % Deterministic framework
                nexttile()
                if isequal(mesh_config,'rectangular')
                    surf(X,Y,c(:,:,indx(n))')
                    axis([0, Lx, 0, Ly, 0, max(Uavg,[],'all')])
                else % hexagonal configurations
                    trimesh(elements,xn,yn,c(:,indx(n)),'FaceColor','interp','EdgeColor','black')
        
                    if isequal(mesh_config,'flat-top')
                        axis([0, Lx+dx/2, 0, Ly, 0, max(Uavg,[],'all')])
                    else
                        axis([0, Lx, 0, Ly+dy/2, 0, max(Uavg,[],'all')])
                    end
                end

                xlab = xlabel('$x$','Interpreter','latex','FontSize',15);
                set(xlab, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);

                ylab = ylabel('$y$','Interpreter','latex','FontSize',15);
                set(ylab, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);

                set(gca,'XTick',[0,Lx],'XTickLabel',{'$0$','$L_{x}$'}, ...
                    'YTick',[0,Ly],'YTickLabel',{'$0$','$L_{y}$'}, ...
                    'TickLabelInterpreter','latex','FontSize',15)

                shading interp
                axis equal
                view(0,90)
                colormap(cork)
                clim([min(f),max(f)])

                % Stochastic framework
                nexttile()
                if isequal(mesh_config,'rectangular')
                    surf(X,Y,Uavg(:,:,indx(n))')
                    axis([0, Lx, 0, Ly, 0, max(Uavg,[],'all')])
                else % hexagonal configurations
                    trimesh(elements,xn,yn,Uavg(:,indx(n)),'FaceColor','interp','EdgeColor','black')
        
                    if isequal(mesh_config,'flat-top')
                        axis([0, Lx+dx/2, 0, Ly, 0, max(Uavg,[],'all')])
                    else
                        axis([0, Lx, 0, Ly+dy/2, 0, max(Uavg,[],'all')])
                    end
                end

                xlab = xlabel('$x$','Interpreter','latex','FontSize',15);
                set(xlab, 'Units', 'Normalized', 'Position', [0.5, -0.04, 0]);

                set(gca,'XTick',[0,Lx],'XTickLabel',{'$0$','$L_{x}$'}, ...
                    'YTick',[],'TickLabelInterpreter','latex', ...
                    'FontSize',15)
                shading interp
                axis equal
                view(0,90)
                colormap(cork)
                clim([min(f),max(f)])

                if n == Nt-1
                    cbh = colorbar;
                    cbh.Layout.Tile = 'south';
                    cbh.FontSize = 10;
                    cbh.TickLabelInterpreter = 'latex';
                end
                
                print(gcf,[path_name,['2D_',config_dir,'_',angle_dir, ...
                    '_',Np_dir,'_t',num2str(n)]],'-depsc2')
            end
    
            % One-dimensional visualisations
            figure
            tiledlayout(2,2,'Padding','loose','TileSpacing','loose')
            for n = 1:Nt
                nexttile()
                hold on
                if isequal(mesh_config,'rectangular')

                    sy = ceil(Ny/2); % slice at approx. Lx/2

                    plot(x,c(:,sy,indx(n))','Color',colours(n,:),'LineWidth',2) % plot deterministic slice
                    % plot(x,Uavg(:,sy,indx(n))','.','Color',colours(n,:),'MarkerSize',10) % plot stochastic average
                    plot(x,Ul(:,sy,indx(n))','Color',colours(n,:),'LineWidth',1) % lower CI bound
                    plot(x,Uu(:,sy,indx(n))','Color',colours(n,:),'LineWidth',1) % upper CI bound

                    for k = 1:Nx % connect lower and upper CI bounds
                        plot([x(k) x(k)],[Ul(k,sy,indx(n)) Uu(k,sy,indx(n))], ...
                            'Color',colours(n,:),'LineWidth',1)
                    end
                    
                    maxUu = max(Uu(:,sy,indx(1)));

                elseif isequal(mesh_config,'flat-top')
                    
                    sy = ceil(Ny/2); % slice at approx. Ly/2

                    ky = (sy-1)*Nx/2 + (1:(Nx/2)); % indexing for slice sx

                    plot(xn(ky),c(ky,indx(n)),'Color',colours(n,:),'LineWidth',2.5) % plot deterministic slice
                    % plot(xn(ky),Uavg(ky,indx(n)),'.','Color',colours(n,:),'MarkerSize',10) % plot stochastic average
                    plot(xn(ky),Ul(ky,indx(n)),'Color',colours(n,:),'LineWidth',1) % lower CI bound
                    plot(xn(ky),Uu(ky,indx(n)),'Color',colours(n,:),'LineWidth',1) % upper CI bound

                    for k = 1:(Nx/2) % connect lower and upper CI bounds
                        plot([xn(ky(k)) xn(ky(k))],[Ul(ky(k),indx(n)) ... 
                            Uu(ky(k),indx(n))],'Color',colours(n,:),'LineWidth',1)
                    end

                    maxUu = max(Uu(ky,indx(1)));

                else

                    sx = ceil(Nx/2); % slice at approx. Lx/2

                    kx = (sx-1)*Ny/2 + (1:(Ny/2)); % indexing for slice sx
                    plot(yn(kx),c(kx,indx(n)),'Color',colours(n,:),'LineWidth',2) % plot deterministic slice
                    plot(yn(kx),Ul(kx,indx(n)),'Color',colours(n,:),'LineWidth',2) % lower CI bound
                    plot(yn(kx),Uu(kx,indx(n)),'Color',colours(n,:),'LineWidth',2) % upper CI bound
                    for k = 1:(Ny/2) % connect lower and upper CI bounds
                        plot([yn(kx(k)) yn(kx(k))],[Ul(kx(k),indx(n)) ... 
                            Uu(kx(k),indx(n))],'Color',colours(n,:),'LineWidth',1);
                    end

                end
                subtitle(strcat('$t=',num2str(round(t(n),2)),'$'), ...
                    'FontSize',25,'Interpreter','latex');
                axis([0 Lx 0 maxUu])
                set(gca,'XTick',[0 Lx/2 Lx],'YTick',[0 max(f)/2 max(f)], ...
                    'FontSize',25)
                

                if mod(n-1,2) == 0
                    ylabel('Particle Density','FontSize',30,'Interpreter','latex')
                end

                if n > Nt-2
                    xlabel('Position ($x$)','Interpreter','latex','FontSize',30)
                end
            end
            
            set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1])
            print(gcf,[path_name,['1D_',config_dir,'_',angle_dir, ...
                '_',Np_dir]],'-depsc2')

            close all

        end

        % Compute quantitative metrics
        if output_metrics
            [epsilon,MSE] = compute_metrics(c,Uavg,Ul,Uu,M,mesh_config);
            fprintf('Confidence interval error: %g\n',epsilon)
            fprintf('Mean squared error: %g\n\n',MSE)
        end
        
        % Input results in a table
        if output_table
            m = m + 1; % next row of table
            current_stats = {epsilon,MSE,theta,D(1,1),D(1,2),D(2,2), ...
                Ex,Ey,[num2str(tau),'(',num2str(M),')'],Npn};
            for k = 1:10
                metrics_table{m,k} = current_stats{k};
            end
        end

    end
    fprintf('---------------------------------------\n\n')

end

if output_table
    save(strcat("figures/",config_dir,'/', ...
        config_dir,"_metrics.mat"),"metrics_table")
end