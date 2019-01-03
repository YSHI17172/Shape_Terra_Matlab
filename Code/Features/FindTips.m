function clscls=FindTips(folderpath)
%FindTips Finds tips which are remaining clusters
%   CLSCLS = FindPlanarSurf(FOLDERPATH) Finds the clusters remaining after 
%   the filtering, plots the tips and saves figures as 
%   Tips.fig in Output\FOLDERPTAH                         

t_start = tic;
% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load clusters, mesh coordinates, filtered points and surfaces and
% similarity percentages for this part from partrecord.mat
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'clusters','tri','fltpts','coord','plnsrf','similold');
% flclsind=size(clusters,1)/2;
%h=figure;

DispStr = ['Finding ' partname ' tips'];
ScreenComment(DispStr,DispStr);
features2=prod(fltpts,2).*plnsrf;
% clsind=find(features2==0);
maxi=0;
clscls=[];
for flclsind=1:length(similold)
    ScreenComment('',['Finding tips for ' ...
        num2str(100*similold(flclsind)) '% similarity']);
    delind=unique(clusters(2*flclsind-1,features2==0));
    
    adjclm=GetClusterAdjMatrix(clusters(2*flclsind-1,:),tri);

    adjclm(delind,:)=0;
    adjclm(:,delind)=0;

    clscls2=ClusterCls(adjclm);
    if max(clscls2)>maxi
        maxi=max(clscls2);
        clscls=clscls2;
        ind=flclsind;
    end
end
    
flclsind=ind; % error here for some parts
% Save tip calculation results to partrecord.mat
save(file,'clscls','flclsind','-append');
ScreenComment('',['Number of tips found: ' num2str(max(clscls))])
cluster=clusters(2*flclsind-1,:);
pts=[];
clus=find(clscls==0);
for i=1:length(clus)
    pts=[pts find(cluster==clus(i))];
end

% Create tips plot and save result
h=figure;
coord2=coord;
coord2(pts,:)=NaN;
trisurf(tri,coord2(:,1),coord2(:,2),coord2(:,3),clusters(2*flclsind,:));
axis equal
title(strcat(partname,' tips'))
colorbar
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'Tips.fig']);
close(h);

t_elapsed = toc(t_start);
ScreenComment('',['Finding tips elapsed time is ' num2str(t_elapsed) '[s]'])
end

function clscls=ClusterCls(nclmat)
nclmat2=nclmat+eye(size(nclmat,1));
count=1;
i=1;
vnew=nclmat2(i,:);
vold=zeros(1,length(vnew));
clscls=zeros(size(nclmat2,1),1);

while any(any(nclmat2))
    if sum(vnew)>1
        while ~isequal(vnew,vold)
            vold=vnew;
            vnew=vold*nclmat2;
            vnew=(vnew>0)*1;
        end
        del=find(vnew);
        clscls(del)=count;
        count=count+1;
    else
        del=i;
    end
    nclmat2(del,:)=0;
    nclmat2(:,del)=0;
    summ=sum(nclmat2,2);
    i=find(summ,1,'first');
    vnew=nclmat2(i,:);
    vold=zeros(1,length(vnew));
    clear del
end
end
function adjclm=GetClusterAdjMatrix(cluster,tri)
conn=[tri(:,1) tri(:,2);tri(:,2) tri(:,3);tri(:,3) tri(:,1)];
clusadj=zeros(3*size(tri,1),2);

for k=1:length(cluster)
    inds1=find(conn(:,1)==k);
    inds2=find(conn(:,2)==k);
    clusno=cluster(k);
    clusadj(inds1,1)=clusno;
    clusadj(inds2,2)=clusno;
end

adjclm=zeros(max(cluster));
for l=1:size(clusadj,1)   
    a=clusadj(l,1);
    b=clusadj(l,2);
    if ~(a==b)
        adjclm(a,b)=1;
    end
end

end