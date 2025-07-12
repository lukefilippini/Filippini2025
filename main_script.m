%% Filippini et al. 2025 %%
clc, clear, close all

% Include subfolders in directory
folder = fileparts(which("main_script.m"));
addpath(genpath(folder));
cmap = curl;

% Specify output of results
visualise_results = true; % create and save comparison figures
output_results = true; % output and save quantitative results

% Parameters
Lx = 1; Ly = 0.8; % length in x and y directions
E = 100; % no. of elements in x or y direction (depends on configuration)
res = 'high'; % mesh resolution ('low','med','high')
mesh_config = 'flat-top'; % mesh configuration ('rectangular','flat-top','pointy-top')
T = 1; % final time
Ns = 5; % no. of stochastic simulations
Npl = [10 25 50 100]; % no. of particles per lattice site 

% Diffusion tensors 
lambda_x = 0.08; lambda_y = 0.02; % tensor eigenvalues
angles = [pi/24 5*pi/12 3*pi/4]; % rotational angles
Nr = length(angles); % no. of rotational angles

% Initial condition
dens = 0.5; % approx. proportion of domain with non-zero density (complex)
Nf = 25; % no. of aggregrations (complex)
IC_type = 'complex'; % type of initial condition ('simple or 'complex')

% Define directory and plotting times (visualise_results == true)
if visualise_results

    % Plotting times (must be an even number)
    t = [0 T/32 T/8 T/4 T/2];
    Nt = length(t);

    % Plotting colours
    colours = [128,0,0; 245,130,48; 0 0 128; 0 130 200]/255;
    
    % Directory for figure storage
    config_dir = [mesh_config,'_',IC_type,'_',res];
    mkdir("./figures",config_dir)
    
end

% Create file for output of results (output_results == true)
if output_results

    % Create text file
    file_path = ['results/',config_dir,'_results.txt'];
    sim_results = fopen(file_path,'w');
    
    % Record setup for test cases
    fprintf(sim_results,'Filippini_et_al_2025: Output of results\n');
    fprintf(sim_results,'---------------------------------------\n');
    fprintf(sim_results,'Domain: Lx = %g, Ly= %g\n',Lx,Ly);
    fprintf(sim_results,'Mesh configuration: %s\n',mesh_config);
    fprintf(sim_results,'Mesh resolution: %s\n',res);
    fprintf(sim_results,'Eigenvalues: lambda_x = %g, lambda_y = %g\n',lambda_x,lambda_y);
    fprintf(sim_results,'Final time: T = %g\n',T);
    fprintf(sim_results,'No. of stochastic simulations: Ns = %g\n',Ns);
    fprintf(sim_results,'Initial condition: %s\n',IC_type);
    fprintf(sim_results,'---------------------------------------\n');

end

% Compute determinsitic and stochastic models
for i = 1:Nr % for each rotational angle / diffusion tensor
    
    % Generate diffusion tensor
    theta = angles(i); % rotation angle
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % rotation matrix
    D = R*[lambda_x 0; 0 lambda_y]*R'; % diffusion tensor
    
    % Check rectangular configuration viability
    if isequal(mesh_config,'rectangular')
        sb = 2 / abs(lambda_x - lambda_y) * sqrt(lambda_x*lambda_y/3); % bound on sin(2*theta)
        if sb < 1
            theta_crit = 1/2*asin(sb); % critical angle
            if ((theta > theta_crit) && (theta < pi/2 - theta_crit)) || ((theta > pi/2 + theta_crit) && (theta < pi - theta_crit))
                fprintf("Rectangular configuration unsuitable for rotation angle %g. \n",theta)
                break
            end
        end
    end

    % Generate no. of elements
    [Ex,Ey] = generate_no_of_elements(Lx,Ly,E,D,mesh_config,res);
    dx = Lx / Ex; dy = Ly / Ey; % spatial steps

    % Generate mesh properties and no. of time steps
    mesh_props = generate_mesh(Lx,Ly,Ex,Ey,mesh_config);
    M = generate_time_steps(dx,dy,T,D,mesh_config);
    
    % Extract properties
    nodes = mesh_props.nodes; xn = nodes(:,1); yn = nodes(:,2); % nodes
    x = unique(xn,'stable'); y = unique(yn,'stable'); % unique nodes
    Nx = length(x); Ny = length(y); % no. of nodes in x and y directions
    elements = mesh_props.elements; % elements
    V = mesh_props.V; % control volumes areas
    tau = T / M; % time step duration

    if isequal(mesh_config,'rectangular')
        [X,Y] = meshgrid(x,y);
    end
    
    % Define initial condition
    if isequal(IC_type,'simple')
        f = generate_simple_IC(xn,yn,Lx,Ly);
    else
        f = generate_complex_IC(xn,yn,Nf,dens,mesh_config);
    end

    % Compute numerical solution
    c = compute_numerical_sol(mesh_props,D,T,M,f,mesh_config);
    
    % Define number of particles
    Np = Npl*(Ex+1)*(Ey+1); % no. of desired particles
    Nps = length(Np); % no. of test cases for each tensor

    % Create a new folder for current angle (visualise_results == true)
    if visualise_results
        angle_dir = ['angle_',num2str(i)];
        mkdir(['./figures/',config_dir],angle_dir)
    end

    % Record current angle (output_results == true)
    if output_results
        fprintf(sim_results,'\nAngle %g: %g\n',i,angles(i));
        fprintf(sim_results,'---------------------------------------\n');
        fprintf(sim_results,'Tensor components: Dxx = %g, Dxy = %g, Dyy = %g\n',D(1,1),D(1,2),D(2,2));
        fprintf(sim_results,'No. of elements: Ex = %g, Ey = %g \n',Ex,Ey);
        fprintf(sim_results,'No. of lattice sites: N = %g \n',(Ex+1)*(Ey+1));
        fprintf(sim_results,'No. of time steps: M = %g \n',M);
        fprintf(sim_results,'Time step duration: tau = %g \n\n',tau);
    end

    % Run stochastic model for each desired no. of particles
    for j = 1:Nps
        
        % Run stochastic simulations
        if isequal(mesh_config,'rectangular')
            [Uavg,Ul,Uu,Npn] = rectangular_stochastic(xn,yn,dx,dy,V,D,T,M,f,Np(j),Ns);
        elseif isequal(mesh_config,'flat-top')
            [Uavg,Ul,Uu,Npn] = flat_top_stochastic(xn,yn,dx,dy,V,D,T,M,f,Np(j),Ns);
        else
            [Uavg,Ul,Uu,Npn] = pointy_top_stochastic(xn,yn,dx,dy,V,D,T,M,f,Np(j),Ns);
        end
        
        % Visualisation of results (visualise_results == true)
        if visualise_results
            
            % Define path for current test case
            Np_dir = ['Np_',num2str(j)];
            mkdir(['./figures/',config_dir,'/',angle_dir],Np_dir)
            path_name = ['./figures/',config_dir,'/',angle_dir,'/',Np_dir,'/'];
            
            % Plot two-dimensional comparisons
            p = 0; % figure counter
            indx = floor(t / tau + 1);
            for n = indx

                p = p + 1; % increase figure count

                figure
                tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
                max_stoch = max(Uavg,[],'all'); % maximum stochastic particle density

                % Deterministic framework
                nexttile()
                if isequal(mesh_config,'rectangular')

                    surf(X,Y,c(:,:,n)')
                    axis([0, Lx, 0, Ly, 0, max_stoch])

                else % hexagonal configurations
                   
                    trimesh(elements,xn,yn,c(:,n),'FaceColor','interp','EdgeColor','black')

                    if isequal(mesh_config,'flat-top')
                        axis([0, Lx+dx/2, 0, Ly, 0, max_stoch])
                    else
                        axis([0, Lx, 0, Ly+dy/2, 0, max_stoch])
                    end

                end
               
                shading interp
                axis equal
                grid off
                view(0,90)
                colormap(cmap)
                clim([min(f),max(f)])
                
                hd = gca; 
                % set(hd,'PositionConstraint','InnerPosition');
                set(hd,'XColor','none','YColor','none')
                % set(hd,'YTick',[0,Ly],'YTickLabel',{'$0$',['$',num2str(Ly),'$']})
                % text(-0.05,Ly/2,'$y$', ...
                %     'Color','k','FontSize',40,'HorizontalAlignment','center', ...
                %     'VerticalAlignment','middle','Rotation',90,'Interpreter','latex')
                % text(zeros(size(hd.YTick))-0.035*Lx,hd.YTick,hd.YTickLabel, ...
                %     'Color','k','FontSize',40,'HorizontalAlignment','center', ...
                %     'VerticalAlignment','middle','Interpreter','latex')
                
                if n == indx(2)
                    title('Deterministic','FontSize',55,'Interpreter','latex')
                end

                if n == indx(end-1)
                    % set(hd,'XTick',[0,Lx],'XTickLabel',{'$0$',['$',num2str(Lx),'$']})
                    % text(Lx/2,-0.05*Ly,'$x$', ...
                    %     'Color','k','FontSize',40,'HorizontalAlignment','center', ...
                    %     'VerticalAlignment','middle','Interpreter','latex')
                    % text(hd.XTick,zeros(size(hd.XTick))-0.1*Ly,hd.XTickLabel, ...
                    %     'Color','k','FontSize',40,'HorizontalAlignment','center', ...
                    %     'VerticalAlignment','bottom','Interpreter','latex')
                end

                % Stochastic framework
                nexttile()
                if isequal(mesh_config,'rectangular')

                    surf(X,Y,Uavg(:,:,n)')
                    axis([0, Lx, 0, Ly, 0, max_stoch])

                else % hexagonal configurations

                    trimesh(elements,xn,yn,Uavg(:,n),'FaceColor','interp','EdgeColor','black')

                    if isequal(mesh_config,'flat-top')
                        axis([0, Lx+dx/2, 0, Ly, 0, max_stoch])
                    else
                        axis([0, Lx, 0, Ly+dy/2, 0, max_stoch])
                    end

                end
                
                shading interp
                axis equal
                grid off
                view(0,90)
                colormap(cmap)
                clim([min(f),max(f)])
                
                hs = gca; 
                set(hs,'XColor','none','YColor','none')
                % set(h,'Units','normalized','Position',[0 0 0.75 0.75])

                if n == indx(end-1)
                    xlabel('$x$','FontSize',40,'Interpreter','latex')      
    
                    % set(hs,'XTick',[0,Lx],'XTickLabel',{'$0$',['$',num2str(Lx),'$']})
                    % text(Lx/2,-0.05*Ly,'$x$', ...
                    %     'Color','k','FontSize',40,'HorizontalAlignment','center', ...
                    %     'VerticalAlignment','middle','Interpreter','latex')
                    % text(hs.XTick,zeros(size(hs.XTick))-0.1*Ly,hs.XTickLabel, ...
                    %     'Color','k','FontSize',40,'HorizontalAlignment','center', ...
                    %     'VerticalAlignment','bottom','Interpreter','latex')
                end
                
                if n == indx(end)
                    cbh = colorbar;
                    cbh.Layout.Tile = 'south';
                    cbh.FontSize = 40;
                    cbh.TickLabelInterpreter = 'latex';
                    cbh.Position(4) = 5*cbh.Position(4);
                end

                if n == indx(2)
                    title('Stochastic','FontSize',55,'Interpreter','latex')
                end
             
               
                % Save figure in specified directory
                set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1],'Renderer','opengl')
                movegui(gcf,'center')
                print(gcf,[path_name,['2D_',config_dir,'_',angle_dir, ...
                    '_',Np_dir,'_t',num2str(p)]],'-dsvg')
            end
    
            % Plot one-dimensional comparisons
            figure
            tiledlayout((Nt-1)/2,2,'Padding','loose','TileSpacing','loose')
            p = 0; % figure counter
            for n = indx(2:end)

                p = p + 1; % increase figure count

                nexttile()
                hold on
                if isequal(mesh_config,'rectangular')

                    sy = ceil(Ny/2); % row of nodes at Ly/2

                    plot(x,c(:,sy,n)','Color',colours(p,:),'LineWidth',2.5) % plot deterministic density
                    plot(x,Ul(:,sy,n)','Color',colours(p,:),'LineWidth',1) % plot 0.025 quantile
                    plot(x,Uu(:,sy,n)','Color',colours(p,:),'LineWidth',1) % plot 0.975 quantile

                    for k = 1:Nx % connect quantile bounds
                        plot([x(k) x(k)],[Ul(k,sy,n) Uu(k,sy,n)], ...
                            'Color',colours(p,:),'LineWidth',1)
                    end
                    
                    if n == indx(2) % store maximum density for first time point
                        maxUu = max(Uu(:,sy,n));
                    end

                elseif isequal(mesh_config,'flat-top')
                    
                    sy = ceil(Ny/2); % row of nodes at Ly/2
                    ry = (sy-1)*Nx/2 + (1:(Nx/2)); % row indexing

                    plot(xn(ry),c(ry,n),'Color',colours(p,:),'LineWidth',2.5) % plot deterministic density
                    plot(xn(ry),Ul(ry,n),'Color',colours(p,:),'LineWidth',1) % plot 0.025 quantile
                    plot(xn(ry),Uu(ry,n),'Color',colours(p,:),'LineWidth',1) % plot 0.975 quantile

                    for k = 1:(Nx/2) % connect quantile bounds
                        plot([xn(ry(k)) xn(ry(k))],[Ul(ry(k),n) ... 
                            Uu(ry(k),n)],'Color',colours(p,:),'LineWidth',1)
                    end
                    
                    if n == indx(2) % store maximum density for first time point
                        maxUu = max(Uu(ry,n));
                    end

                else % pointy-top configuration

                    sx = ceil(Nx/2); % column of nodes at Lx/2
                    cx = (sx-1)*Ny/2 + (1:(Ny/2)); % column indexing

                    plot(yn(cx),c(cx,n),'Color',colours(p,:),'LineWidth',2.5) % plot deterministic density
                    plot(yn(cx),Ul(cx,n),'Color',colours(p,:),'LineWidth',1) % plot 0.025 quantile
                    plot(yn(cx),Uu(cx,n),'Color',colours(p,:),'LineWidth',1) % plot 0.975 quantile

                    for k = 1:(Ny/2) % connect quantile bounds
                        plot([yn(cx(k)) yn(cx(k))],[Ul(cx(k),n) ... 
                            Uu(cx(k),n)],'Color',colours(p,:),'LineWidth',1);
                    end
                    
                    if n == indx(2) % store maximum density for first time point
                        maxUu = max(Uu(cx,n));
                    end

                end

                subtitle(strcat('$t=',num2str(round(t(p+1),2)),'$'), ...
                    'FontSize',25,'Interpreter','latex');
                axis([0 Lx 0 max(1.2,maxUu)])
                set(gca,'XTick',[0 Lx/2 Lx],'YTick',[0 max(f)/2 max(f)], ...
                    'FontSize',25,'TickLabelInterpreter','latex')
                
                % Plot y label for left column of figures
                if mod(p-1,2) == 0
                    ylabel('Particle Density','FontSize',30,'Interpreter','latex')
                end
                
                % Plot x label for bottom row of figures
                if p >= Nt-1
                    xlabel('Position ($x$)','Interpreter','latex','FontSize',30)
                end
            end
            
            set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1])
            print(gcf,[path_name,['1D_',config_dir,'_',angle_dir, ...
                '_',Np_dir]],'-depsc2')

            close all

        end
               
        % Compute quantitative metrics and output results
        if output_results
            MSE = compute_metrics(c,Uavg,Ul,Uu,M,mesh_config);
            
            fprintf(sim_results,'No. of desired particles: Np = %d (N x %g)\n',Np(j),Npl(j));
            fprintf(sim_results,'No. of scaled particles: Np = %d\n',Npn);
            %fprintf(sim_results,'Epsilon: %g\n',epsilon);
            fprintf(sim_results,'MSE: %g\n\n',MSE);

        end
       
    end

    fprintf(sim_results,'---------------------------------------\n\n');

end

% Close output of results (output_results == true)
if output_results
    sim_results = fclose('all');
end
