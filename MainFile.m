function MainFile(listfile,commentlevel)
%MainFile Run ShapeTerra analysis on specified parts in a list.
%   MainFile(LISTFILE,COMMENTLEVEL) wil run analysis on filenames in the
%   LISTFILE list. Depending on value of COMMENTLEVEL display show:
%   0: No on-screen comments
%   1: Concise on-screen comments
%   2: Extensive on-screen comments
%
%   Output of MainFile will be stored in the Output folder with a different
%   folder for each part named "[FILENAME-TIMESTAMP]"
%   result consisting of MATLAB figures and calculation .mat file
% start using repo
clc, close all, clearvars -except listfile commentlevel

% Depending on operating system set file separator
OS = computer;
if strcmp(OS, 'PCWIN') || strcmp(OS, 'PCWIN64')
    % MATLAB operating on a windows PC
    filesep = '\';
else % GLNX86 || GLNXA64 || MACI64 (wasn't tested on Linux)    
    % MATLAB operating on Mac or other UNIX system
    filesep = '/';
end
global fs
fs = filesep;

% Add code folders to path, to call all functions in the Code subfolders
addpath(['ShapeTerra' fs 'Code' fs 'Backbone' fs]);
addpath(['ShapeTerra' fs 'Code' fs 'Thinness' fs]);
addpath(['ShapeTerra' fs 'Code' fs 'Features' fs]);
addpath(['ShapeTerra' fs 'Code' fs 'Additive' fs]);
addpath(['ShapeTerra' fs 'Code' fs 'Subtractive' fs]);
addpath(['ShapeTerra' fs 'Code' fs 'Plotting' fs]); % For HKS plotting

% Get day and tome for output folder name with timestamp
day = date;
time= fix(clock);

% Read listfile name
% Check if listfile is given, if not set listfile to 'standard'
if ~exist('listfile','var') || isempty(listfile)
  listfile = 'standard';
end

% Set global variable cmnt_level once by user for screen comment display
% Set comment level:    0 show no progress comments
%                       1 show short progress comments
%                       2 show extensive progress comments
% Check if COMMENTLEVEL is given, if not set global cmnt_level to one
global cmnt_level
if ~exist('commentlevel','var') || isempty(commentlevel)
    cmnt_level = 1;
else
    cmnt_level = commentlevel;
end

% Determine path to list file
path=(['ShapeTerra' fs 'Input' fs 'Lists' fs listfile '.txt']);
ind=6;
if exist(path)==0
    ScreenComment('List File does not exist!',...
        ['List File in path ' path ' does not exist!']);
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
    ScreenComment('List File found, starting computations..',...
        ['List File found, Files to be processed: ' num2str(i) ...
        ', starting computations..']);
    
	% Loop over different partnames in listfile
    for j=1:i
        clear filename timestamp
        filename=cell2mat(C(j));
        % Create unique folder for this calculation run of the part in Output
        % with subfolders for meshsize containing subfolders for timestamp
        timestamp = [day '-' num2str(time(4)) 'hr' num2str(time(5))];
        
        % Detect meshsize info in filename and convert to mm value
        indnum = regexp(filename,'\d'); % Find positions numbers in filename
        diffindnum = diff(indnum); % Determine position differences in indnum
        numdiv = find(diffindnum>1); % Determine separation of numbers related 
                                     % to part number and meshsize 
        indmeshnum = indnum(numdiv+1:end);
        indpartnum = indnum(1:numdiv);
        %partnum = str2double(filename(indpartnum))
        meshsize = str2double(filename(indmeshnum));
        
        % Create folderpath with tradename, partname and meshsize value
        tradename = filename(1:indpartnum(1)-1);
        % Remove mesh size info from filename
        partname = filename(indpartnum(1):indmeshnum(1)-1); 
        folderpath = [tradename fs partname  fs 'Mesh-' ...
            num2str(meshsize,indnum) 'mm' fs timestamp];
        mkdir(['ShapeTerra' fs 'Output' fs folderpath]);
        ScreenComment(['Computations for part: ' filename],...
            ['Computations for part: ' filename]);
        
        % Try to run calculations from MainCalc for this part, else skip
      %  try
            MainCalc(folderpath,filename);
        %end
%         catch err
%             disp(['Cannot perform calculations for ' filename])
%             continue
%         end
        
        disp(['Done' 10])
        %SendEmail(['ShapeTerra ' listfile ' File Processing Completed'],...
        %    [filename ' ShapeTerra processing is done']);
    end
    ScreenComment('All files processed',...
        ['All files of ' listfile ' have been processed']);
    %SendEmail(['ShapeTerra' listfile ' Files All Processed'],...
    %    ['All files of ' listfile ' have been processed']);
end

% Remove Code subfolder from path to clean up
rmpath(['ShapeTerra' fs 'Code' fs 'Backbone' fs]);
rmpath(['ShapeTerra' fs 'Code' fs 'Thinness' fs]);
rmpath(['ShapeTerra' fs 'Code' fs 'Features' fs]);
rmpath(['ShapeTerra' fs 'Code' fs 'Additive' fs]);
rmpath(['ShapeTerra' fs 'Code' fs 'Subtractive' fs]);
rmpath(['ShapeTerra' fs 'Code' fs 'Plotting' fs]);
end

% Function to send an email once the ShapeTerra calculations are done
function[]=SendEmail(subj,msg)
myaddress = 'persistent.heat@gmail.com';
mypassword = ''; % Not set now for privacy reasons, ask ramy for pwd

setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(myaddress, subj, msg);
clear
end
