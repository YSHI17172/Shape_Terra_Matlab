function [abc,cg]=FindProtrusions(folderpath)
%FindProtrusions Finds the protrusion/pocket axis and cross section of features.
%   [ABC,CG] = FindProtrusions(FOLDERPATH) Finds the protrusion/pocket axis 
%   of the features of FOLDERPATH part and their cross section, plots figures
%   and saves these as FeaturesAxes.fig in Output\FOLDERPATH

tic
% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load clusters and other variables for feature identification
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'features','coord','clusters','tri');

% Identify protrusion/pocket direction and cross sections
ptft=features;
n=max(ptft);
t=-6:0.1:6;
% Display screen comments
DispStr = 'Identifying protrusion/pocket direction and cross-section';
ScreenComment(DispStr,[DispStr ' for:'])

coord2=coord;
coord2(features==0)=NaN;
abc=[];
cg=[];
% Generate plot with features + axes

h=figure;
trisurf(tri,coord2(:,1),coord2(:,2),coord2(:,3),clusters(size(clusters,1),:));
hold on
zm=max(coord(:,3));
axis equal

for i=1:n
    ScreenComment('',['Identifying ' partname ' feature #' num2str(i) ...
        ' direction and cross-section']);
    
    pts=find(ptft==i);
    ptscoord2=coord(pts,:);
    cluster=clusters(size(clusters,1)-1,:);
    clusk=unique(cluster(pts));

    [Vi,cgi]= FindProtrusionAxis(ptscoord2);
    abci=Vi(:,3);
    ptchk=coord(find(cluster==min(clusk),1,'first'),:)-cgi;
    dir=ptchk*abci;
    if dir<0
        Vi=-Vi;
        abci=-abci;
    end
    abc(i,:)=abci';
    cg(i,:)=cgi;
    
    xyz=abci*t+cgi'*ones(1,length(t));
    xyz=xyz';
    %hold on
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3));
    [npts ori]=FindProtrusionCrossection(ptscoord2,Vi,zm);
    scatter3(npts(:,1),npts(:,2),npts(:,3));
    clear clus pts ptscoord ptscoord2 x y F X Y Z
end
title(strcat(partname,' features+axes'))
colorbar
% Save figure and close
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'FeatureAxes.fig']);
close(h);

% Save abc, cg and ori variables to partrecord.mat
save(file,'abc','cg','ori','-append');

t_elapsed = toc;
% Display screen comment level two
ScreenComment('',['Elapsed time finding protrusions: '...
    num2str(t_elapsed) '[s]']);
end
    
function [W,cg]=FindProtrusionAxis(ptscoord)
cg=mean(ptscoord);
ptscoord(:,1)=ptscoord(:,1)-cg(1);
ptscoord(:,2)=ptscoord(:,2)-cg(2);
ptscoord(:,3)=ptscoord(:,3)-cg(3);
[U,S,V]=svd(ptscoord);
W=[V(:,3),V(:,2),V(:,1)];
end

function [npts ori]=FindProtrusionCrossection(pts,V,zm)
abc=V(:,3);
ori=1;
% v2=[0;0;0];
% v3=v2;
% v2(1)=sqrt(abc(3)^2/(abc(1)^2+abc(3)^2));
% v2(3)=-abc(1)*v2(1)/abc(3);
% v3(2)=sqrt(abc(3)^2/(abc(2)^2+abc(3)^2));
% v3(3)=-abc(2)*v2(2)/abc(3);
% W=[v2 v3 abc];
W=V;
cg=mean(pts);
d=-sum(abc'.*cg);

pts2(:,1)=pts(:,1)-cg(1);
pts2(:,2)=pts(:,2)-cg(2);
pts2(:,3)=pts(:,3)-cg(3);

dotps=pts2*abc;
npts2=pts(find(dotps>=0),:);

h=(abs(npts2*abc+d)/sum(abc.^2))';
nptscoord=npts2-(abc*h)';
crssn=nptscoord;
nptscoord(:,1)=nptscoord(:,1);
nptscoord(:,2)=nptscoord(:,2);
nptscoord(:,3)=nptscoord(:,3)-cg(3);
if dot(abc,[0;0;1])<0;
   nptscoord(:,3)=nptscoord(:,3)+zm+0.5;
   ori=0;
end
    

npts=nptscoord;
% npts(:,1)=npts(:,1)+cg(1);
% npts(:,2)=npts(:,2)+cg(2);

end