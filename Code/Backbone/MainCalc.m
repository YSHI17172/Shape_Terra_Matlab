function MainCalc(folderpath,filename)
%MainCalc calls ShapeTerra calculation functions and plot functions
%   MainCalc(FOLDERPATH,FILENAME) is the main ShapeTerra calculation file
%   It calls all other ShapeTerra calculation functions.
%   FOLDERPATH is name of the subfolder in the output folder 
%   of the part FILENAME for which calculations are performed
%   MainCalc conists of the following blocks:
%   - Backbone
%       - ReadFile
%       - PlotFileMesh
%       - GenerateHKS        
%       - Persistence
%       - Cluster
%       - PlotFile
%       - GetCoGs
%       - PlotFileCoGs
%       - FilterClusters
%   - Features
%       - FindPlanarSurf
%       - FindTips
%       - FindFeatures
%       - FindProtrusions
%       - PlotFileLO
%   - Thinness
%       - FindIfThin
%       - PlotFileSDDr
%       - PlotSDDrrate2

%% General parameters
% Set calculation parameters - Move to MainFile - Barend
simil=[0.7 0.75 0.8 0.85 0.9];  %Persistence level similarities
r=0.1; %Percentage of points in a single planar surface to be filtered
per=0.0005; %heat value threshold to compute persistence level and value
Nmax = 40000; % Maximum number of points allowed per mesh file

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call ShapeTerra backbone calculations functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Try to read mesh coordinates Catia .dat/.off file 
%   If possible generate and save coord and tri, plot mesh 
%   else return to MainFile loop
try
    ReadFile(folderpath,filename,Nmax);      
catch err
   disp(['Cannot generate ' filename ' mesh'])
   return
end   

% 2. Generate and save the Heat Kernel Signature (HKS)
%   Usually GenerateHKS(filename,step_size,number_of_iteraions)
%   If step size and number of iterations not specified, use default values:
%   step_size=0.001, number_of_iterations=1000;
GenerateHKS(folderpath);
% The following is for testing purposes
% PlotHKS(folderpath);
% return
% 3. Generates and saves persistence level and value for defined threshold
Persistence(folderpath,per);

% 4. Cluster points and plots for specified similarity level, save the 
%   generated clusters in folderpath/partecord.mat and save figures
Cluster(folderpath,simil);

% 5. Find and save the Center of Gravity of the different clusters
GetCoGs(folderpath);

% 6. Filter clusters for different persistence similarities, save results
%   in folderpath and plot the filtered part at different similarity
FilterClusters(folderpath);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call ShapeTerra Features calculation functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Find planar surfaces that contain more than r*100 % of the total points
%  in the part, save them and plot the filtered part
FindPlanarSurf(folderpath,r); 

% 2. Find the clusters remaining after the filtering, plot the tips
%   and save results to folderpath
% try
    FindTips(folderpath);  
% catch err
%     error(err)
%     return
% end   

% 3. Find the features from the tips, plot them and save results
FindFeatures(folderpath);

% 4. Find the protrusion/pocket axis of the features and their cross section
%   plot them and save plots in folderpath, also plot leftovers of part
FindProtrusions(folderpath);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call ShapeTerra Thinness calculation functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Find thin features
FindIfThin(folderpath);    

% 2. Plot SDDr figures and save them
PlotFileSDDr(folderpath);         

% 3. Plot SDDrate2 figures and save
PlotSDDrrate2(folderpath);        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show final ShapeTerra results                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FeatureFigure(folderpath);
end