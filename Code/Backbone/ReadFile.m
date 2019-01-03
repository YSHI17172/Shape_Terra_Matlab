function [coord,tri]=ReadFile(folderpath,filename,Nmax)
%ReadFile Extract mesh coordinates and tri data from mesh file.
%   [COORD,TRI] = ReadFile(FOLDERPATH,FILENAME,NMAX) reads mesh COORD and TRI 
%   data from FILENAME.dat in the input folder and returns COORD and TRI 
%   mesh variables. NMAX is maximum number of mesh points allowed to continue
%
%   ReadFile checks if the COORD and TRI already exist in all output subfolders
%   of this part. If exact COORD and tri data are found the existing COORD
%   and TRI data are used. Else, the COORD and TRI data is (re)generated for
%   the new calculation run. Coord and Tri saved in new partrecord.mat

tic
% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load input files and set path to new partrecord file
file=strcat(['ShapeTerra' fs 'Input' fs 'Database' fs filename]);
varfile=(['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat']);

% Search output folder for subfolders with filename and same mesh,
% then check for partrecord.mat in that subsubfolder
[mesh_match,part_flag] = CheckMesh(folderpath);
if isempty(find(part_flag,1))
    % No partrecord found for this meshsize
    DispStr = [partname ' partrecord.mat not found'];
    ScreenComment([DispStr ', generating Coord and Tri'],...
        [DispStr ' in Output' fs tradename fs partname fs meshname fs '..'...
        10 'Generating mesh variables Coord and Tri']);
    s1=strcat(file,'.dat');
    s2=strcat(file,'.off');
    disp(file)
    if exist(s1)>0
        [coord,tri]=ReadDAT(s1);
        save(varfile,'coord','tri');
        ScreenComment('',['Read ' partname ' mesh data from ' s1]);
    elseif exist(s2)>0
        [coord,tri]=ReadOFF(s2); 
        save(varfile,'coord','tri');
        ScreenComment('',['Read ' partname ' mesh data from ' s2]);
    end
    t_elapsed = toc;
    ScreenComment('',['Tri and Coord generated in ' num2str(t_elapsed) '[s]']);
else % partrecord for FILENAME found in one of the subsubfolders of output
    if isempty(mesh_match)
        % Partrecord found but no similar mesh
        ScreenComment('No similar mesh found, generating Coord and Tri',...
            [partname ' partrecord.mat found in Output' fs meshname 10 ...
            'No similar mesh found, generating Coord and Tri']);
        % Read new Coord and Tri data from either .dat or .off file
        s1=strcat(file,'.dat');
        s2=strcat(file,'.off');
        if exist(s1)>0       
            [coord,tri]=ReadDAT(s1);
        elseif exist(s2)>0
            [coord,tri]=ReadOFF(s2); 
        end
        t_elapsed = toc;
        ScreenComment('',['Tri and Coord generated in ' num2str(t_elapsed) '[s]']);
    else
        % ScreenComment partrecord and similar mesh found
        ScreenComment('Similar mesh found, Coord and Tri found',...
            [partname ' partrecord.mat found in ' mesh_match{end} 10 ...
            'Similar mesh found, Coord and Tri found']);
        % Load Coord and Tri data from last timestamp subfolder mesh_match{end}
        partdata=['ShapeTerra' fs 'Output' fs tradename fs partname fs meshname...
            fs mesh_match{end} fs 'partrecord.mat'];
        load(partdata,'coord','tri');
    end
    
    % Display warning and possibly cut of here if #points > 30k
    Npoints = size(coord,1);
    if Npoints > Nmax
        ScreenComment(['Mesh contains more than maximum ' num2str(Nmax) ...
            ' points'],['Mesh number of points ' num2str(Npoints) ...
            ' exceeds maximum of ' num2str(Nmax)])
        return
    else
        ScreenComment('','Mesh number of points lower than maximum # points')
    end
    
    % Save new or existing mesh to partrecord in new Output subfolder
    % folderpath
    save(varfile,'coord','tri');
    % The following is for testing purposes
    %dlmwrite([partname '_coord.csv'],coord,'precision',15);
    %dlmwrite([partname '_tri.csv'],tri,'precision',15)
end

% Generate and save mesh plots
PlotMesh(coord,tri,folderpath);
end


% Function to read Mesh data from .off file
% Not yet tested!!
function [coord,tri]=ReadOFF(path)
fid=fopen(path);
fgetl(fid);
prop=str2num(fgetl(fid));
nopts=prop(1);
notri=prop(2);

coord=zeros(nopts,3);
tri=zeros(notri,3);

for i=1:nopts
    val=str2num(fgetl(fid));
    coord(i,:)=val;
    clear val;
end

for i=1:notri
    val=str2num(fgetl(fid));
    tri(i,:)=val(2:4);
    clear val;
end
tri=tri+ones(notri,3);
fclose(fid);
end


function [coord,tri,cfile,names] =ReadDAT(path)
[names]=textscan(fopen(path),'%s');
% Determine if .dat file is in CATIA R19 or R20 format
% 21st place in R19 file is '*', in R20 this is at the 22nd place
% If both tests fail assume R19 file format
if strcmp(names{1}(21),'*')
    CATIAv = 'R19';
elseif strcmp(names{1}(22),'*')
    CATIAv = 'R20';
else
    CATIAv = 'R?-unknown';
end

% Display for cmnt_level what .dat CATIA mesh version is detected
ScreenComment('',['Mesh .dat file-version detected CATIA V5' CATIAv]);

% Create vector with converted double values names
cfile = str2double(names{:});

% Start reading in file where XYZ coordinates are
i=18;
j=1;

% Read XYZ coordinates from file
while ~isnan(cfile(i))

    if strcmp(CATIAv,'R20')
        coord(j,1)=cfile(i+1);
        coord(j,2)=cfile(i+2);
        coord(j,3)=cfile(i+6);
        i = i + 8; 

    else % 'R19'
        coord(j,1)=cfile(i+1);
        coord(j,2)=cfile(i+2);
        coord(j,3)=cfile(i+4);
        i=i+6;
    end
   j = j + 1;
end

% Determine i counter to jump to start of tri numbers in the .dat file
% Step will be deleted in the future when own mesh generation is implemented
dollarind = find(strcmp(names{:},'$..')==1);

% Jump to second last $.. index + 2, here the tri part starts
i = dollarind(end-1)+2;
l=1;

% Read trimetric mesh data from file
while ~isnan(cfile(i))

    tri(l,1)=cfile(i+2);
    tri(l,2)=cfile(i+3);
    tri(l,3)=cfile(i+4);
    i=i+6;
    l=l+1;
end
end


% Function to check if there already is a same mesh in folderpath
function [mesh_match,part_flag] = CheckMesh(folderpath)
% Load global file separator fs
global fs

% Retrieve tradename, partname and meshname from folderpath
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

% Search for mesh_match in Output\partname\meshsize timestamp subfolders
mesh_match = {''}; % Create empty cell array to store names of 
                   % subsubfolders with similar meshes in it
part_flag = zeros(1,length(FoInfo));
for FoIdx = 1:length(FoInfo)
    %if strfind(FoInfo(FoIdx).name,filename)
    filetest=['ShapeTerra' fs 'Output' fs tradename fs partname fs meshname...
        fs FoInfo(FoIdx).name fs 'partrecord.mat'];
    if exist(filetest)
        % Partrecord exists
        part_flag(FoIdx) = 1;
        load(filetest,'coord','tri');
        if exist('coord') && exist('tri')
            % Partrecord found, coord and tri found
            mesh_match(FoIdx) = {FoInfo(FoIdx).name};
        end
    else
        part_flag(FoIdx) = 0;
    end
end

% Remove empty entries from mesh_match and sort
mesh_match(cellfun(@isempty,mesh_match)) = [];
mesh_match = sort(mesh_match);
end


function PlotMesh(coord,tri,folderpath)
%PlotFileMesh Plot the trimetric CATIA part mesh
%   PlotFileMesh(COORD,TRI,FOLDERPATH) will plot trimetric CATIA part mesh
%   data COORD and TRI, plot is saved in folderpath\mesh.fig

% Load global file separator fs
global fs

% Determine partname folder
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Create plot
h = figure();
trisurf(tri,coord(:,1),coord(:,2),coord(:,3),0.5*ones(size(coord,1),1));
title([partname ' trimetric mesh'])
colorbar
axis equal
% Save plot in Output folder and close figure
saveas(h,(['ShapeTerra' fs 'Output' fs folderpath fs 'mesh.fig']));
close(h);

% Display extensive comments only
ScreenComment('',['Plot and save ' partname ' mesh']);
end