function [v] = Persistence(folderpath,per)
%Persistence Generates & saves persistence level & value for defined threshold
%   v = Persistence(FOLDERPATH,PER) generates the level of persistence for
%   strucutral part under investigation. v is a matrix of size
%   [#points x #2] and consists of peristence value and level per point
%
%   FOLDERPATH is the partrecord path where calculation results are stored
%   PER is the persistence threshold defined by the user, if not given it
%   is set to the maximum HKS value at the last iteration

tic
% Load global file separator fs
global fs

% Load partrecord file of current part run
indsep = strfind(folderpath,fs);
tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
meshname = folderpath(indsep(2)+1:indsep(3)-1);
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'HKS','stepold','iterold');

% Set/determine persistence level
if nargin==1
    thresh=max(HKS(:,size(HKS,2)-1));
else
    thresh=per;
end

% Check if persistence calculation results exist
[per_match,per_flag] = CheckPersistence(folderpath,iterold,stepold,thresh);

% If not or if changes calculate peristence value and level
if isempty(find(per_flag,1))
    ScreenComment(['Persistence levels not found. Generating '...
        'Persistence levels'],'Persistence levels not found.');
    tic

    n=size(HKS,1);
    v=zeros(n,2);
    for i=1:1:n
        k=find(HKS(i,:)<=thresh,1,'first')-1;
        if isempty(k)
            k=size(HKS,2);
        end
        s=sum(HKS(i,1:k));
        v(i,1)=k;
        v(i,2)=s;
    end
    hthreshold=thresh;
    %save(file,'v','hthreshold','-append');
    t_elapsed = toc;
    ScreenComment('',...
        ['Persistence levels generated in ' num2str(t_elapsed) '[s]']);
else % v variable found in previous calculation results
    if isempty(per_match)
        DispStr = 'Persistence levels found with different threshold level';
        ScreenComment([DispStr 10 'Generating Persistence levels with '...
            'new parameters'],DispStr);
        tic
        n=size(HKS,1);
        v=zeros(n,2);
        for i=1:1:n
            k=find(HKS(i,:)<=thresh,1,'first')-1;
            if isempty(k)
                k=size(HKS,2);
            end
            s=sum(HKS(i,1:k));
            v(i,1)=k;
            v(i,2)=s;
        end
        hthreshold=thresh;
        t_elapsed = toc;
        ScreenComment('',['Persistence levels with new parameters '...
            'generated in ' num2str(t_elapsed) '[s]']);
    else % ~isempty(Per_match)
        DispStr = 'Persistence levels found with matching parameters';
        ScreenComment(DispStr,DispStr);
        
        % Load found v and htreshold data from most recent same partrecord.mat 
        partdata=['ShapeTerra' fs 'Output' fs tradename fs partname fs ...
            meshname fs per_match{end} fs 'partrecord.mat'];
        load(partdata,'v','hthreshold');
    end
end
% Save v and hthreshold results to current partrecord.mat in folderpath
save(file,'v','hthreshold','-append');
end

% Check if same persistence levels exist in folderpath with same steps+iters
function [per_match,per_flag] = CheckPersistence(folderpath,iters,steps,thresh)
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

% Search for Per_match in Output\partname\meshsize timestamp subfolders
per_match = {''}; % Create empty cell array to store names of 
                   % subsubfolders with similar v matrices in it
per_flag = zeros(1,length(FoInfo));
for FoIdx = 1:length(FoInfo)
    % Clear v from previous loop iteration if it exists
    if exist('v')
        clear v
    end
    filetest=['ShapeTerra' fs 'Output' fs tradename fs partname fs meshname...
        fs FoInfo(FoIdx).name fs 'partrecord.mat'];
    %load(filetest,'HKS');
    % Check if HKS can be found in partrecord.mat
    try
        per_test = whos('-file',filetest,'-regexp','v');
    catch err
        per_test = [];
    end

    if ~isempty(per_test) % So if v exists
       % Clear iterold, stepold and hthreshold from previous loop iteration
        if exist('stepold') || exist('iterold') || exist('hthreshold')
            clear stepold iterold hthreshold
        end

        % v exists
        per_flag(FoIdx) = 1;
        load(filetest,'v','stepold','iterold','hthreshold');
        
        % Check if number of iters, steps and persistence levels are the same
        if (steps~=stepold) || (iters~=iterold) || (thresh~=hthreshold)
            % different v/htreshold found
        else
            % v found with same iter and step and htreshold
            per_match(FoIdx) = {FoInfo(FoIdx).name};
        end
    else % v does not exist
        per_flag(FoIdx) = 0;
    end
end

% Remove empty entries from per_match and sort
per_match(cellfun(@isempty,per_match)) = [];
per_match = sort(per_match);
end

    