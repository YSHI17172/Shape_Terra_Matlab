function []=FindPlanarSurf(folderpath,r)
%FindPlanarSurf Finds planar surfaces by disctintion of number of R*100 points
%   FindPlanarSurf(FOLDERPATH,R) finds planar surfaces of part, with data in
%   FOLDERPATH/partrecord.mat, that contain more than R*100 % of the total
%   points in the part, plots them and saves figures as FilteredSurfaces.fig
%   in Output\FOLDERPATH

tic
% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load clusters and mesh coordinates for this part
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'coord','tri','clusters');

plnsrf=zeros(max(max(tri)),1);
normv=zeros(size(coord));
tritri=[tri;tri(:,2),tri(:,3),tri(:,1);tri(:,3),tri(:,1),tri(:,2)]; 

ScreenComment('Filtering large planar surfaces',...
    'Filtering large planar surfaces');
for i=1:size(coord,1)
    trind=find(tritri(:,1)==i);
    
    adjpt1=tritri(trind,2);
    adjpt2=tritri(trind,3);
    
    coord1=coord(adjpt1,:);
    coord2=coord(adjpt2,:);
    
    v1=coord1-ones(size(coord1,1),1)*coord(i,:);
    v2=coord2-ones(size(coord2,1),1)*coord(i,:);
    
    normvis=cross(v1,v2);
    
    normvi=sum(normvis,1);
    
    normv(i,:)=normvi./sqrt(sum(normvi.^2));
    clear trind adjpt1 adjpt2 coord1 coord2 v1 v2 normvis normvi
end

surfn=zeros(size(coord,1),1);
pln=[];
count=1;

while any(surfn==0)
    i=find(surfn==0,1,'first');
    abc=normv(i,:);
    d=-sum(abc.*coord(i,:));
    
    dist=abs(sum((ones(size(coord,1),1)*abc).*coord,2)+d*ones(size(coord,1),1));
    npts=find(dist==0);
    surfn(npts)=count;
    pln(count)=length(npts);
    count=count+1;
    clear i abc d dist npts
end

pln=pln./size(coord,1);
plns=find(pln<r);
pts=[];
for i=1:length(plns)
    pts=[pts find(surfn==plns(i))'];
end
plnsrf(pts)=1;

% Save planar surface results to pathrecord.mat in Output\folderpath folder
save(file,'plnsrf','-append');

% Create figure and save the results in Output\folderpath
h=figure;
coord(plnsrf==0,:)=NaN;
trisurf(tri,coord(:,1),coord(:,2),coord(:,3),clusters(size(clusters,1)));
axis equal
title([partname ' filtered surfaces'])
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'FilteredSurfaces.fig']);
close(h);

t_elapsed = toc;
% Display screen comment level two
ScreenComment('',['Elapsed time filtering large planar surfaces: '...
    num2str(t_elapsed) '[s]'])
end