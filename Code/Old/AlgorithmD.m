% TODO
% Add more elaborate comments last three functions called by AlgorithmD
function AlgorithmD(folderpath,filename)
%AlgorithmD calls ShapeTerra calculation functions and plot functions
%   AlgorithmD(FOLDERPATH,FILENAME) calls ShapeTerra functions to perform 
%   calculations. FOLDERPATH is name of the subfolder in the output folder 
%   of the part FILENAME for which calculations are performed
%   AlgorithmD calls the following functions:
%       - ReadFile
%       - PlotFileMesh
%       - GenerateHKS        
%       - Persistence
%       - Cluster
%       - PlotFile
%       - GetCoGs
%       - PlotFileCoGs
%       - FilterClusters
%       - FindPlanarSurf
%       - FindTips
%       - FindFeatures
%       - FindProtrusions
%       - PlotFileLO
%       - FindIfThin
%       - PlotFileSDDr
%       - PlotSDDrrate2

% Set calculation parameters
simil=[0.7 0.75 0.8 0.85 0.9];  %Persistence level similarities
r=0.1; %Percentage of points in a single planar surface to be filtered
per=0.0005; %heat value threshold to compute persistence level and value

% Call ShapeTerra calculations functions
ReadFile(folderpath,filename);      %Reads from Catia, generates and saves coord and tri  
%
%
PlotFileMesh(folderpath);         %Plots the initial mesh and saves it in "Mesh" Folder
%
%
GenerateHKS(folderpath);          %Generates and saves HKS. Usually
%                               GenerateHKS(filename,step_size,number_of_iteraions)
%                               If step size and number of iterations not specified, 
%                               Use default values:
%                               step_size=0.001
%                               number_of_iterations=1000;
%
Persistence(folderpath,per);      %Generates and saves persistence level and
%                               value for defined threshold
%
%
Cluster(folderpath,simil);        %Cluster points for specified similarity level
%                               %and save the generated clusters
% 
%
PlotFile(folderpath);             %Plots the mesh clustered at different persistence
%                               leves and saves it in the "Plots" folder
%
GetCoGs(folderpath);              %Finds and saves the CoGs of the clusters
%
% 
PlotFileCoGs(folderpath);         %Plots the mesh clustered at different persistence
%                               leves along with their COGs and saves it in the 
%                               "Plots + COGs" folder
%
% 
FilterClusters(folderpath);       %Filters clusters for different persistence
%                               similarities, saves the results and plots the
%                               filtered part at different similarity and saves
%                               in the "Filtered Parts" folder
%
% 
FindPlanarSurf(folderpath,r);     %Finds planar surfaces that contain more
%                               than r*100 % of the total points in the part,
%                               saves them and plots the filtered part and 
%                               saves it in the "Filtered Surfaces" folder 
% 
%
FindTips(folderpath);             %Finds the clusters remaining after the filtering,
%                               saves them and plots the tips in the "Tips" folder
%
% 
FindFeatures(folderpath);         %Finds the features from the tips, saves them,
%                               plots them and saves the plots in the
%                               "Features" Folder
%
% 
FindProtrusions(folderpath);      %Finds the protrusion/pocket axis of the
%                               features and their cross section, saves
%                               them, plots them and saves the plots in the
%                               "Features + axes" folder
%
% 
PlotFileLO(folderpath);           %Plots the leftovers of the part and saves them in
%                               the "Left overs" folder
%
%
FindIfThin(folderpath);           % Find thin features
%
%
PlotFileSDDr(folderpath);         % Plot SDDr figures and save them
%
%
PlotSDDrrate2(folderpath);        % Plot SDDrate2 figures and save
end