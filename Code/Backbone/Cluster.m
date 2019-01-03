function [] = Cluster(folderpath,simil)
%Cluster Clusters points and saves generated clusters
%   v = Cluster(FOLDERPATH,SIMIL) cluster points for part specified 
%   similarity levels SIMIL and saves the generated clusters as 'clusters' 
%   Clusters is a matrix of size [2*size(simil)]x[#points] and save results
%   in folderpath\partrecord.mat

t_start = tic;
% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load variables needed for clustering from partrecord.mat
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'coord','tri','v','hthreshold','stepold','iterold');

simil=unique(simil); % Remove duplicate entries in simil
% Check if (same) clusters exist in partrecord.mat
[CLTR_match,CLTR_flag] = CheckCluster(folderpath,iterold,stepold,hthreshold);

% If not or if different: generate new clusters
if isempty(find(CLTR_flag,1))
    clusters=zeros(2*length(simil),size(coord,1));
    ScreenComment('No clusters exist, generating new clusters',['No '...
        partname ' clusters exist']);
    for l=1:length(simil)
        tic
        ScreenComment('',['Clustering for ',num2str(simil(l)*100),'% similarity']);
        newcluster=Cluster5(coord,tri,simil(l),v);
        clusters(2*l-1:2*l,:)=newcluster;
        t_elapsed = toc;
        ScreenComment('',[num2str(simil(l)*100),'% similarity cluster calculation in '...
            num2str(t_elapsed) '[s]']);
    end
elseif isempty(CLTR_match) % No matching clusters found
    ScreenComment(['Clusters persistence level threshold has changed.' ...
        10 ' Clustering for new threshold'],'');
    clear similold clusters
    clusters=[];
    clusters=zeros(2*length(simil),size(coord,1));
    for l=1:length(simil)
        tic
        ScreenComment('',['Clustering for ',num2str(simil(l)*100),'% similarity']);
        newcluster=Cluster5(coord,tri,simil(l),v);
        clusters(2*l-1:2*l,:)=newcluster;
        t_elapsed = toc;
        ScreenComment('',[num2str(simil(l)*100),'% similarity cluster calculation in '...
            num2str(t_elapsed) '[s]']);
    end   
else
    % Clusters with same threshold level exist, check if there are new 
    % entries in simil which do not exist in partrecord.mat:
    partdata=['ShapeTerra' fs 'Output' fs tradename fs partname fs ...
        meshname fs CLTR_match{end} fs 'partrecord.mat'];
    load(partdata,'clusters','similold');
        
    % if so run cluster calculation for new level
    clusters2=zeros(2*length(simil),size(coord,1));
    Cmnt1_flag = 0; % Flag to display level 1 comment only once 
    for l=1:length(simil)
        % Check simil level in partrecord.mat
        chk=find(similold==simil(l),1);
        % If there is not match for this simil level run again
        if isempty(chk)
            % Display level1 screencomment only once
            if Cmnt1_flag == 0
                ScreenComment('Create new clusters for new similarity values','')
            end
            Cmnt1_flag = 1;
            tic
            newcluster=Cluster5(coord,tri,simil(l),v);
            clusters2(2*l-1:2*l,:)=newcluster;
            t_elapsed = toc;
            ScreenComment('',['Clustering for ' num2str(simil(l)*100)...
                '% similarity done in ' num2str(t_elapsed) '[s]']);
        else % Clusters with same simil and threshold exist
            % Display level1 screencomment only once
            if Cmnt1_flag == 0
                ScreenComment('No clusters with new similarity values','')
            end
            Cmnt1_flag = 1;
            ScreenComment('',['Cluster for ' num2str(simil(l)*100)...
                '% similarity already exists'])
            newcluster=clusters(2*chk-1:2*chk,:);
            clusters2(2*l-1:2*l,:)=newcluster;
        end
        clear chk
    end
    clusters=[];
    clusters=clusters2;
end
%hthreshclus=hthreshold;
clear similold
similold=simil;
% Save variables to Output\folderpath\partrecord.mat
save(file,'clusters','similold','-append');  

% Display level 2 ScreenComment elapsed time
t_elapsed = toc(t_start);
ScreenComment('',['Elapsed time clustering: ' num2str(t_elapsed) '[s]']);

% Plot and save clusters
PlotCluster(clusters,coord,similold,tri,folderpath);
end


% Check if same clusters exists in folderpath partrecord
function [CLTR_match,CLTR_flag] = CheckCluster(folderpath,iters,steps,thresh)
% Load global file separator fs
global fs

% Retreive partname and meshfolder from folderpath
indsep = strfind(folderpath,fs);
tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Define folder where you want so search timestamp subfolders in
FoInfo = dir(['ShapeTerra' fs 'Output' fs tradename fs partname fs meshname]);
FoInfo(~[FoInfo.isdir]) = []; % Remove non-folder entries in folder info
RemIdx(1) = find(strcmp({FoInfo.name},{'.'})==1);
RemIdx(2) = find(strcmp({FoInfo.name},{'..'})==1);
FoInfo(sort(RemIdx)) = []; % Remove . and .. entries from folder info

% Search for CLTR_match in Output\partname\meshsize timestamp subfolders
CLTR_match = {''}; % Create empty cell array to store names of 
                   % subsubfolders with similar clusters in it
CLTR_flag = zeros(1,length(FoInfo));
for FoIdx = 1:length(FoInfo)
    % Clear HKS from previous loop iteration if it exists
    if exist('clusters')
        clear clusters
    end
    filetest=['ShapeTerra' fs 'Output' fs tradename fs partname fs meshname...
        fs FoInfo(FoIdx).name fs 'partrecord.mat'];
    %load(filetest,'HKS');
    % Check if clusters can be found in partrecord.mat
    try
        CLTRtest = whos('-file',filetest,'-regexp','clusters');
    catch err
        ScreenComment(['No partrecord.mat found at ' filetest],...
            ['No partrecord.mat found at ' filetest])
        CLTRtest = [];
    end

    if ~isempty(CLTRtest)
       % Clear hthreshclus from prev loop iter if they exist
        if exist('hthreshold')
            clear hthreshold
        end

        % clusters  and hthreshclus exists 
        CLTR_flag(FoIdx) = 1;
        load(filetest,'clusters','hthreshold','iterold','stepold');
        
        % Now check persistence level threshold has changed
        if (steps~=stepold) || (iters~=iterold) || (hthreshold~=thresh)
            % different cluster persistence levels found and or different
            % iter + step size for peristence calculations
        else
            % clusters found with same persistence level threshold
            CLTR_match(FoIdx) = {FoInfo(FoIdx).name};
        end
    else % clusters do not exist
        CLTR_flag(FoIdx) = 0;
    end
end

% Remove empty entries from CLTR_match and sort
CLTR_match(cellfun(@isempty,CLTR_match)) = [];
CLTR_match = sort(CLTR_match);
end

function clusters = Cluster5(coord,tri,simil,v)
notri = size(tri,1);
nopts = size(coord,1);
neighbor_similarity = simil;

p1=v(:,2);

I = [tri(:,1);tri(:,2);tri(:,3)];
J = [tri(:,2);tri(:,3);tri(:,1)];
adj = sparse(I,J,ones(3*notri,1));
[x,i] = sort(p1,'descend');
clusters = zeros(2,size(coord,1));
count=1;
visited = zeros(1,nopts);
while (isempty(i) == 0)
    list = i(1);
    rootval = x(1);
    while (isempty(list) == 0)
        if visited(list(1)) == 0
        adjpts = find(adj(list(1),:)==1);
        adjpts = adjpts(find(visited(adjpts)==0));
        else
            adjpts = [];
        end
%         adjpts

%         adj(:,list(1)) = 0;
%         adj(list(1),:) = 0;
        visited(list(1)) = 1;
        curpos = find(i==list(1));
        banana = p1(adjpts);
        citrus = neighbor_similarity*rootval;
        addpts = find(p1(adjpts)>(neighbor_similarity*rootval));
        list = [list;adjpts(addpts)'];
        clusters(1,list(1)) = count;
        clusters(2,list(1)) = rootval;
        list(1) = [];
        i(curpos) = []; x(curpos) = [];
         if (isempty(list) == 1)
              count=count+1;
         end
    end
end
end


function PlotCluster(clusters,coord,similold,tri,folderpath)
%PlotCluster Plots and save clusters with same similarity persistence levels 
%   PlotCluster(FOLDERPATH) plots clustered persistence similarity levels for 
%   specific part and saves the generated similarity figures as
%   folderpath\ClusteredPersistence-%.fig

% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Plot and save figures
ScreenComment('Saving Cluster plots','')
h=figure();

n=length(similold);
val = zeros(size(coord,1), n);
label = cell(n,1);
for i=1:n
    val(:, i) = clusters(2*i,:)';
    label{i} = ['Clusters for ', num2str(similold(i)*100), '% similarity'];
end
h = SliderFigure(coord, tri, val, label);
title('Clusters');
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'ClusteredPersistence.fig']); 
% for i=1:length(similold)
%     trisurf(tri,coord(:,1),coord(:,2),coord(:,3),clusters(2*i,:));
%     axis equal
%     title(strcat(partname,' ',num2str(similold(i)*100),'%'))
%     colorbar
%     ScreenComment('',['Saving ' num2str(similold(i)*100) '%' ...
%         ' peristence similarity cluster plot'])
%     saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'ClusteredPersistence-'...
%         num2str(similold(i)*100) '.fig']);    
% end

close(h);
end