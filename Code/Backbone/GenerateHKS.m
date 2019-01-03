function HKS = GenerateHKS(folderpath,stepin,iterin)
%GenerateHKS Generates and saves Heat Kernel Signature. Usually
%   HKS = GenerateHKS(FOLDERPATH,STEPSIZE,#ITERATIONS) generates the HKS 
%   of the  structural part under investigation. HKS is a matrix of size
%   [#points x #iterations] and consists of point heat values per point timestep
%
%   FOLDERPATH points to .mat file where calculations are stored
%   STEPSIZE is the timestep value in seconds
%   #ITERATIONS is the number of timesteps that are calculated
%   If step STEPSIZE and #ITERATIONS not specified use default values
%   STEPSIZE = 0.001,  #ITERATIONS = 1000

tic
% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
meshname = folderpath(indsep(2)+1:indsep(3)-1);
varfile=(['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat']);

% Determine if number of iterations is given, if not use default time step
% and number of iteraions
if nargin==3
    steps=stepin;
    iters=iterin;
else
    steps=0.001;
    iters=1000;
end

% Call CheckHKS to see if HKS results already exist
[HKS_match,HKS_flag] = CheckHKS(folderpath,steps,iters);

% Check if HKS calculation has been done before, if not run HKS calculation
if isempty(find(HKS_flag,1))
    ScreenComment('No HKS found. Generating HKS',...
        ['No HKS found: ShapeTerra' fs 'Output' fs tradename fs partname ...
        fs 'meshfolder' fs ' does not contain HKS data' 10 ...
        'Generating new HKS matrix']);
    
    % Generate new HKS matrix with calcHKS function
    [HKS,stepold,iterold] = CalcHKS(varfile,steps,iters);

else % HKS file(s) for PARTNAME found in one of the subsubfolders of output
    if isempty(HKS_match)
        ScreenComment('Different HKS found, generating new HKS',...
        ['HKS found with different step and iter parameters. Generating HKS'...
            ' for new parameters']);
        
        % Calculate HKS matrix for new parameters
        [HKS,stepold,iterold] = CalcHKS(varfile,steps,iters);
    
    else
        ScreenComment('HKS found with matching parameters',...
        'HKS found with matching parameters');
        
        % Load found HKS iter and step data from most recent same partrecord.mat 
        partdata=['ShapeTerra' fs 'Output' fs tradename fs partname fs ...
            meshname fs HKS_match{end} fs 'partrecord.mat'];
        load(partdata,'HKS','stepold','iterold');
        
    end 
end

%%%%%%%%%%%%%%%%%%%
%--> FIGURE
load(varfile,'coord','tri');
HKS_labels = cell(size(HKS,2),1);
for i=1:size(HKS,2)
    HKS_labels{i} = ['time = ' num2str(1000*i*steps, '%.1f') 'ms'];
end
%h = SliderFigure(coord, tri, HKS, HKS_labels);
%saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'HKS.fig']); 
%<--
%%%%%%%%%%%%%%%%%%%

% Save results to partrecord.mat
save(varfile,'HKS','stepold','iterold','-append');
end

% Calculate new HKS matrix
function [HKS,stepold,iterold] = CalcHKS(varfile,steps,iters)
    % Load variables coord and tri from partrecord
    load(varfile,'coord','tri');
    nopts=size(coord,1);
    notri=size(tri,1);
    
    % Compute Laplacian
    ScreenComment('Computing Laplacian Eigenvalues..','Computing Laplacian Eigenvalues..');
    tic
    
    [L,M] = Laplacian(coord,tri);
    eigno = 300;
    options=struct('disp',0);
    [V,D,flag] = eigs(M,L,eigno);
    if flag ~= 0
       errrrrrrror=1
    end
   
    i1 = tri(:,1); i2 = tri(:,2); i3 = tri(:,3);
    v1 = coord(i3,:) - coord(i2,:);  v2 = coord(i1,:) - coord(i3,:); v3 = coord(i2,:) - coord(i1,:);
    n  = cross(v1,v2,2); 
    dblA = (sqrt(sum((n').^2)))';
    totalarea = sum(dblA)/2;

    V = real(V);
    D = real(D);
    D = diag(D);
    D = diag((totalarea/4/pi)*ones(size(D))./D);
    
    clear flag eigno options
    
    t_elapsed = toc;
    ScreenComment('',['Laplacian computed in ' num2str(t_elapsed) '[s]']);
    
    % Start iterations for HKS generation
    t =  steps;                                                                
    ScreenComment('','Iterating for HKS matrix..')
    t_start = tic;
    
    HKS = zeros(size(V,1), iters);
    diagD = diag(D);
    
    for i = 1:iters
        % Record at what time first iteration begins
        if i==1
            t_iter = tic;
        end
        % Display a level 2 ScreenComment every 25 iters
        if rem(i,25)==0
            t_elapsed = toc(t_iter);
            ScreenComment('',['HKS iteration ' num2str(i) ' elapsed time '...
                num2str(t_elapsed) '[s]']);
            % Record at what time new batch of 25 iterations begins
            t_iter = tic;
        end
        
        % Calculate HKS values
        H = (V.^2)*exp(diagD*t);
        HKS(:, i) = H;
        
        t = t + steps;
        
        if i==iters
            t_elapsed = toc(t_iter);
            ScreenComment('','Last HKS iteration completed');
        end
    end 
    stepold=steps;
    iterold=iters;
    
    % Display level 2 ScreenComment elapsed time
    t_elapsed = toc(t_start);
    ScreenComment('HKS computation finished',...
        ['Total HKS computation time: ' num2str(t_elapsed) '[s]']);
end

% Check if same HKS exists in folderpath
function [HKS_match,HKS_flag] = CheckHKS(folderpath,steps,iters)
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

% Search for HKS_match in Output\partname\meshsize timestamp subfolders
HKS_match = {''}; % Create empty cell array to store names of 
                   % subsubfolders with similar HKS matrices in it
HKS_flag = zeros(1,length(FoInfo));
for FoIdx = 1:length(FoInfo)
    % Clear HKS from previous loop iteration if it exists
    if exist('HKS')
        clear HKS
    end
    filetest=['ShapeTerra' fs 'Output' fs tradename fs partname fs meshname...
        fs FoInfo(FoIdx).name fs 'partrecord.mat'];
    %load(filetest,'HKS');
    % Check if HKS can be found in partrecord.mat
    try
        HKStest = whos('-file',filetest,'-regexp','HKS');
    catch err
        HKStest = [];
    end

    if ~isempty(HKStest)
       % Clear iterold and stepold from previous loop iteration if they exist
        if exist('stepold') || exist('iterold')
            clear stepold iterold
        end

        % HKS exists
        HKS_flag(FoIdx) = 1;
        load(filetest,'HKS','stepold','iterold');
        
        % Now check if number of iterations and stepsize HKS are the same
        if (steps~=stepold) || (iters~=iterold)
            % different HKS found
        else
            % HKS found with same iter and step
            HKS_match(FoIdx) = {FoInfo(FoIdx).name};
        end
    else % HKS does not exist
        HKS_flag(FoIdx) = 0;
    end
end

% Remove empty entries from HKS_match and sort
HKS_match(cellfun(@isempty,HKS_match)) = [];
HKS_match = sort(HKS_match);
end