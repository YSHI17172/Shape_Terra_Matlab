function [coord,tri]=ReadMesh(filename)
%ReadMesh Extract mesh coordinates and tri data from mesh file.
%   [COORD,TRI] = ReadFile(FILENAME) reads mesh COORD and TRI 
%   data from FILENAME.dat in the input folder and returns COORD and TRI 
%   mesh variables.

tic
% Load global file separator fs
global fs

% Load input files and set path to new partrecord file
file=strcat(['..' fs '..' fs 'Input' fs 'Database' fs filename]);

% Generate new mesh to check if it can be done
s1=strcat(file,'.dat');
s2=strcat(file,'.off');
ScreenComment('Checking list file','Checking list file');
if exist(s1)>0
    try % Catch error if necessary
        [coord,tri,~,~]=ReadDAT(s1);
    catch err
        % Display custom error message
        error([pwd fs 'Code' fs 'Checking' fs 'ReadFile:ReadDAT format error'...
            10 'Mesh .dat file is not in correct Catia V5-R19 or -R20 format'])
    end  % end try/catch
    ScreenComment('',['Read ' filename ' mesh data from ' s1]);
elseif exist(s2)>0
    try % Catch error if necessary
        [coord,tri,~,~]=ReadOFF(s2); 
    catch err
        % Display custom error message
        error([pwd fs 'Code' fs 'Checking' fs 'ReadFile:ReadOFF format error'...
            10 'Mesh .dat file is not in correct Catia V5-R19 or -R20 format'])
    end  % end try/catch
    ScreenComment('',['Read ' filename ' mesh data from ' s2]);
else
    coord = [];
    tri = [];
    DispStr = ['Cannot read ' filename ' mesh data from either ' s1 ' or ' s2];
    ScreenComment(DispStr,DispStr);
end
t_elapsed = toc;
ScreenComment('',['Tri and Coord generated in ' num2str(t_elapsed) '[s]']);
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