function createfigure(clscls,clusters,coord,flclsind,tri)
cluster=clusters(2*flclsind-1,:);
pts=[];
clus=find(clscls==0);
for i=1:length(clus)
    pts=[pts find(cluster==clus(i))];
end

coord2=coord;
coord2(pts,:)=NaN;
trisurf(tri,coord2(:,1),coord2(:,2),coord2(:,3),clusters(2*flclsind,:));
axis equal
title('Tips plot')
colorbar
end