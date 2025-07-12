figure

% Deterministic framework
if isequal(mesh_config,'rectangular')

    surf(X,Y,c(:,:,1)')
    axis([0, Lx, 0, Ly, 0, max_stoch])

else % hexagonal configurations
   
    trimesh(elements,xn,yn,c(:,1),'FaceColor','interp','EdgeColor','black')

    if isequal(mesh_config,'flat-top')
        axis([0, Lx+dx/2, 0, Ly, 0, max(f)])
    else
        axis([0, Lx, 0, Ly+dy/2, 0, max(f)])
    end

end

shading interp
axis off; axis equal
grid off
view(0,90)
colormap(cmap)
clim([min(f),max(f)])

%  hd = gca; 
% set(hd,'PositionConstraint','InnerPosition');
% set(hd,'XColor','none','YColor','none')
% set(hd,'YTick',[0,Ly],'YTickLabel',{'$0$',['$',num2str(Ly),'$']})
% text(-0.065,Ly/2,'$y$', ...
%     'Color','k','FontSize',40,'HorizontalAlignment','center', ...
%     'VerticalAlignment','middle','Rotation',90,'Interpreter','latex')
% text(zeros(size(hd.YTick))-0.05*Lx,hd.YTick,hd.YTickLabel, ...
%     'Color','k','FontSize',40,'HorizontalAlignment','center', ...
%     'VerticalAlignment','middle','Interpreter','latex')