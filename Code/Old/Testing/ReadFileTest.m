function [coord,tri,cfile,names]=ReadFileTest
%ReadFile Extract mesh coordinates and tri data from mesh file.
%   [COORD,TRI] = ReadFile(FOLDERPATH,FILENAME) reads mesh COORD and TRI 
%   data from FILENAME.dat in the input folder and returns COORD and TRI 
%   mesh variables.
%
%   ReadFile checks if the COORD and TRI already exist in all output subfolders
%   of this part. If exact COORD and tri data are found the existing COORD
%   and TRI data are used. Else, the COORD and TRI data is (re)generated for
%   the new calculation run. Coord and Tri saved in new partrecord.mat

tic
% Depending on operating system set file separator
OS = computer;
if strcmp(OS, 'PCWIN') || strcmp(OS, 'PCWIN64')
    % MATLAB operating on a windows PC
    filesep = '\';
else % GLNX86 || GLNXA64 || MACI64 (wasn't tested on Linux)    
    % MATLAB operating on Mac oor other UNIX system
    filesep = '/';
end
global fs
fs = filesep;

addpath(['ShapeTerra' fs 'Code' fs 'Backbone' fs]);
filename = 'TPockets02MultiHeights02';
testpath = ['TPockets' fs '02MultiHeights' fs 'Mesh-2mm' fs '13-Oct-2014-16hr39'];
%mkdir(['ShapeTerra' fs 'Output' fs testpath]);
% Retreive names from folderpath
indsep = strfind(testpath,fs);
tradename = testpath(1:indsep(1)-1);
partname = testpath(indsep(1)+1:indsep(2)-1);
meshname = testpath(indsep(2)+1:indsep(3)-1);

% Load input files and set path to new partrecord file
file=strcat(['ShapeTerra' fs 'Database' fs 'Input' fs filename]);
varfile=(['ShapeTerra' fs 'Output' fs testpath fs 'partrecord.mat']);

global cmnt_level
commentlevel = 2;
cmnt_level = commentlevel;

% Search output folder for subfolders with filename and same mesh,
% then check for partrecord.mat in that subsubfolder
[mesh_match,part_flag] = CheckMesh(testpath);
if isempty(find(part_flag,1))
    % No partrecord found for this meshsize
    DispStr = [partname ' partrecord.mat not found'];
    ScreenComment([DispStr ', generating Coord and Tri'],...
        [DispStr ' in Output' fs tradename fs partname fs meshname ...
        ' subfolders, generating Coord and Tri']);
    s1=strcat(file,'.dat')
    s2=strcat(file,'.off');
    if exist(s1)>0
        %try % Catch error if necessary
            [coord,tri]=ReadDAT(s1);
            save(varfile,'coord','tri');
        %catch err
            % Display custom error message
         %   error([pwd fs 'ReadFile:ReadDAT error' 10 'Check if .dat file is'...
          %      ' in correct Catia V5-R19 or -R20 format'])
        %end  % end try/catch
        ScreenComment('',['Read ' partname ' mesh data from ' s1]);
    elseif exist(s2)>0
        try % Catch error if necessary
            [coord,tri]=ReadOFF(s2); 
            save(varfile,'coord','tri');
        catch err
            % Display custom error message
            error([pwd fs 'ReadFile:ReadOFF error' 10 'Check if .dat file is'...
                ' in correct Catia V5-R19 or -R20 format'])
        end  % end try/catch
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
        s1=strcat(file,'.dat')
        s2=strcat(file,'.off');
        if exist(s1)>0
            [coord,tri,cfile,names]=ReadDAT(s1)
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
    
    % Save new or existing mesh to partrecord in new Output subfolder
    % folderpath
    save(varfile,'coord','tri');
end
end

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
    CATIAv = 'R19';
end

% Display for cmnt_level what .dat CATIA mesh version is detected
ScreenComment('',['Mesh .dat file-version detected CATIA V5' CATIAv]);

% cfile=zeros(length(names{1}),1);
% for i=1:length(names{1})
%     cfile(i)=str2double(names{1}(i));
% end
cfile = str2double(names{:});

i=18;
j=1;

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
%i=i+14;
i = dollarind(end-1)+2;
l=1;

while ~isnan(cfile(i))

    tri(l,1)=cfile(i+2);
    tri(l,2)=cfile(i+3);
    tri(l,3)=cfile(i+4);
    i=i+6;
    l=l+1;
end


end