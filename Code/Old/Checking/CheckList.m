function CheckList(listfile,commentlevel)
clc, close all, clearvars -except listfile commentlevel
Nmax = 13000; % Maximum number of points allowed per part

% Read listfile name
% Check if listfile is given, if not set listfile to 'standard'
if ~exist('listfile','var') || isempty(listfile)
  listfile = 'standard';
end

% Set comment level:    0 show no progress comments
%                       1 show short progress comments
%                       2 show extensive progress comments
% Check if comment level cmt_level is given, if not set level to one
if ~exist('commentlevel','var') || isempty(commentlevel)
  commentlevel = 1;
end
global cmnt_level
cmnt_level = commentlevel;

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

% Make checking directory in output
mkdir(['..' fs '..' fs 'Output' fs 'Checking' fs]);
% Determine path to list file
path=['..' fs '..' fs 'Input' fs 'Lists' fs listfile '.txt'];
if exist(path)==0
    disp(['List file ' listfile ' does not exist!'])
else
    fid=fopen(path);
    i=1;
    C={};
    while ~feof(fid)
        clear filename
        filename=fgetl(fid);
        C=[C;filename];
        i=i+1;
    end
    i=i-1;
    disp(['List file ' listfile ' found'])
    disp(['Files to be checked: ' num2str(i)])
    points_flag = zeros(1,i);
    for j=1:i
        clear filename coord tri
        filename=cell2mat(C(j));
        disp(['File: ' filename])
        %Reads from catia, generates and saves coord and tri  
        [coord,tri] = ReadMesh(filename);    
        Npoints = size(coord,1);
        if Npoints > Nmax
            ScreenComment('Npoints > Nmax',['Npoints > Nmax: '...
                num2str(Npoints) ' > ' num2str(Nmax)]);
            points_flag(i) = 1;
        else
            points_flag(i) = 0;
        end
        % Plots the initial mesh
        PlotMesh(filename,coord,tri);     
        %
        disp('Done')
        disp(' ')
        %SendEmail(strcat(listfile,' File Completion'),strcat(filename,' is done processing'));
    end
    disp('All files processed succesfully')
    %SendEmail(strcat(listfile,' all Files Completion'),strcat('All files of',' ',listfile,'have been processed succesfully'));
end
end