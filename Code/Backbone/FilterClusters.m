function []=FilterClusters(folderpath)
%FilterClusters Filter clusters for different persistence similarities
%   FilterClusters(FOLDERPATH) filter clusters of part belonging to folderpath
%   for different persistence similarities, plots the filtered part at various
%   similarity percentages and saves figures as FilteredPart-%.fig 
%   in Output\FOLDERPATH

tic
% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load clusters, similarity precentages and mesh coordinates for this part
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'clusters','tri','similold','coord');

% Determine number of filtering runs, one for each similarity level
n=length(similold);
fltpts=zeros(max(max(tri)),n);

% Create figure handle for new plots
h=figure;

% Filter clusters, create plots and save results
ScreenComment('Filtering Clusters','');
val = zeros(size(coord,1), n);
label = cell(n,1);
alpha = zeros(size(coord,1), n);
for i=1:n
    coord2=coord;
    ScreenComment('',['Filtering clusters for ' num2str(similold(i)*100)...
        '% similarity']);
    rempts=FilterClustersi(clusters,tri,i);    
    fltpts(rempts,i)=1;
    
    alpha(rempts, i) = 1;
    val(:, i) = clusters(2*i,:)';
    label{i} = ['Filtered clusters for ', num2str(similold(i)*100), '% similarity'];
    
%     coord2(:,:)=NaN;
%     coord2(rempts,:)=coord(rempts,:);
%     trisurf(tri,coord2(:,1),coord2(:,2),coord2(:,3),clusters(2*i,:));
%     axis equal
%     title([partname ' filtered for ' num2str(similold(i)*100) '% similarity']);
%     colorbar
%     saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'FilteredPart-'...
%         num2str(similold(i)*100) '.fig']);    
    clear rempts
end

h = SliderFigure(coord, tri, val, alpha, label);
title('Filtered Clusters');
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'FilteredPart.fig']); 

% Save filtered points to partrecord.mat
save(file,'fltpts','-append');
close(h);
end


% Filtering functions
function rempt=FilterClustersi(clusters,tri,k)
i=1;
n=2*i-1;
inicl=clusters(n,:);
cl=FindBiggestCluster(inicl);

n=2*k-1;
cluster=clusters(n,:);

flt=(inicl==cl)*1;
ncluster=cluster.*flt;
lst=FindBiggestCluster(ncluster);

clmat=GetClusterAdjMatrix(cluster,tri);
clmat=clmat+eye(size(clmat));
nclmat=Filter(clmat,lst);

del=find(sum(nclmat,2)<1);
nclmat(del,:)=0;
nclmat(:,del)=0;
[clsi,~]=find(nclmat);
cls=unique(clsi);
rempt=[];
for i=1:length(cls)
    rempt=[rempt find(cluster==cls(i))];
end

end

function adjclm=GetClusterAdjMatrix(cluster,tri)
conn=[tri(:,1) tri(:,2);tri(:,2) tri(:,3);tri(:,3) tri(:,1)];
clusadj=zeros(3*size(tri,1),2);
clusadj(:,1) = cluster(conn(:,1));
clusadj(:,2) = cluster(conn(:,2));

adjclm=zeros(max(cluster));
for l=1:size(clusadj,1)   
    a=clusadj(l,1);
    b=clusadj(l,2);
    if (adjclm(a,b)==0) && (adjclm(b,a)==0) && (a~=b)
        adjclm(a,b)=1;
        adjclm(b,a)=1;
    end
end
end

function cl=FindBiggestCluster(cluster)
i=1;
j=max(cluster);
sumcl=0;
cl=0;
for k=i:j
    nsumcl=sum(cluster==k);
    if nsumcl>sumcl
        sumcl=nsumcl;
        cl=k;
    end
end
end

function clmat=Filter(clmat,lst)

for i=1:length(lst)
    l=find(clmat(:,lst(i)));
    clmat(:,lst(i))=0;
    nlst=find(clmat(lst(i),:));    
    if length(nlst)>1
        clmat(lst(i),:)=0;
        % ambiguition
        %clmat(:,nlst)=0;
        %
        clmat=Filter(clmat,nlst);
    else
        clmat(l,lst(i))=1;
%         clmat(:,lst(i))=zeros(size(clmat,1),1);
    end
        
end
end