function PlotCoGs(clusters,coord,COGs,similold,tri,folderpath)
    %PlotFileCoGs Plots persistence clusters with their CoGs
    %   PlotFileCoGs(FOLDERPATH) Plots persistence clusters with their CoGs
    %   for part belonging to to Output\FOLDERPATH\CoGClusteredPersistance-%.fig

    % Define x, y and z coordinates
    x=COGs(:,1,:);
    y=COGs(:,2,:);
    z=COGs(:,3,:);

    % Plot and save CoG figures
    figure()
    i = 2;
    trimesh(tri,coord(:,1),coord(:,2),coord(:,3),clusters(2*i,:),'edgecolor',...
        [0 0 0],'marker','o','markerfacecolor','flat','markersize',1);
    axis equal
    alpha 0;
    hold on
    plot3(x(:,i),y(:,i),z(:,i),'o','MarkerFaceColor',[0 0 0],'MarkerSize',18);
    hold off
    colorbar  
end