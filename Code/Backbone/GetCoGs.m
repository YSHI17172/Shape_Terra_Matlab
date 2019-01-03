% TODO
% Update matrix size COG and COGPs in function description
function [COGs,COGPs]=GetCoGs(folderpath)
%GetCoGs Calculates and saves CoGs of persistence clusters
%   [COGs,COGPs]=GetCoGs(FILENAME) calculates and saves CoGs of persistence
%   clusters for part FILENAME to Output\filename\partrecord.mat 
%   COGs and COGPs are matrices of size [?]x[?]x[#clusters]

tic
% Load global file separator fs
global fs

% Load files from output folder for folderpath
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'coord','clusters','similold','tri');

n=size(clusters,1);
l=n/2;
inds=1:2:n-1;

cluster=clusters(inds,:);

m=max(max(cluster));

COGs=NaN(m,3,l);
COGPs=NaN(m,3,l);

ScreenComment('Computing CoGs','Computing CoGs of persistence clusters')
for i=1:l
    [COGsi,COGPsi]=CoG(coord,cluster(i,:));
    d=size(COGsi,1);
    COGs(1:d,:,i)=COGsi(:,:);
    COGPs(1:d,:,i)=COGPsi(:,:);
    clear COGsi COGPsi d
end
% Save results to partrecord.mat
save(file,'COGs','COGPs','-append');

% Display level 2 ScreenComment elapsed time
t_elapsed = toc;
ScreenComment('',['Elapsed time computing CoGs: ' num2str(t_elapsed) '[s]']);

% Plot and save COGs
PlotCoGs(clusters,coord,COGs,similold,tri,folderpath);
end


function [COG,COGP]=CoG(coord,clusters)

n=max(clusters);

COG=zeros(n,3);
COGP=zeros(n,3);

for i=1:n
    cluster=find(clusters==i);
    m=length(cluster);
    x=sum(coord(cluster,1))/m;
    y=sum(coord(cluster,2))/m;
    z=sum(coord(cluster,3))/m;
    COG(i,:)=[x y z];
    
    coords=coord(cluster,:);
    cogx=x*ones(m,1);
    cogy=y*ones(m,1);
    cogz=z*ones(m,1);
    cogcoords=[cogx cogy cogz];
    dists=sum((coords-cogcoords).^2,2);
    dis=min(dists);
    clospt=find(dists==dis,1);
    COGP(i,:)=coords(clospt,:);
    
    
    clear cluster m x y z coords cogx cogy cogz cogcoords dists dis clospt  
end
end


function PlotCoGs(clusters,coord,COGs,similold,tri,folderpath)
%PlotFileCoGs Plots persistence clusters with their CoGs
%   PlotFileCoGs(FOLDERPATH) Plots persistence clusters with their CoGs
%   for part belonging to to Output\FOLDERPATH\CoGClusteredPersistance-%.fig

% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Define x, y and z coordinates
x=COGs(:,1,:);
y=COGs(:,2,:);
z=COGs(:,3,:);

% Plot and save CoG figures
ScreenComment('Saving CoGs plots','')
h=figure;
for i=1:length(similold)
    trimesh(tri,coord(:,1),coord(:,2),coord(:,3),clusters(2*i,:),'edgecolor',...
        [0 0 0],'marker','o','markerfacecolor','flat','markersize',7);
    axis equal
    alpha 0;
    hold on
    plot3(x(:,i),y(:,i),z(:,i),'o','MarkerFaceColor',[0 0 0],'MarkerSize',18);
    hold off
    title(strcat(partname,' ',num2str(similold(i)*100),'%'))
    colorbar 
    ScreenComment('',['Saving CoG ' num2str(similold(i)*100) ...
        '%-peristence similarity cluster plot'])
    saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'CoG-ClusteredPersistence-'...
        num2str(similold(i)*100) '.fig']);  
end

close(h);
end