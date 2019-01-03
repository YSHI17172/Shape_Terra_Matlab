function []=PlotSDDrrate2(folderpath)
%PlotSDDrrate2 Plot and save SDDrate2 figure.
%   PlotSDDrrate2(FOLDERPATH)) plots the SDDrrate2 figures of part beloning
%   to FOLDERPATH\partrecord.mat and saves plots as SDDrrate2.fig in 
%   Output\FOLDERPATH

% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load variables for SDDrrate2 plot
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'colorv','features','coord','tri','abc','cg');

% Generate level 2 screen comment
ScreenComment('',['Generating ' partname ' SDDrrate2 figure']);

% Loop over the number of features
for i=1:max(features)
    ind=find(features==i);
    % New coordinate system for each feature with z longest direction
    D=coord(ind,:)*abc(i,:)';
%     t=-((abc(i,:)*cg(i,:)')*ones(size(D))+D)./norm(abc(i,:))^2;
%     maxp=ind(find(t==max(t)),1);
%     minp=ind(find(t==min(t)),1);
%     
%     h(i)=sqrt(sum((coord(maxp,:)-coord(minp,:)).^2));
    
    h(i) = max(D(:)) - min(D(:));
end

h=h./min(h);
colorf=zeros(size(features));

for i=1:max(features)
    ind=find(features==i);
    colorf(ind)=h(i)*mean(colorv(ind));
    feature_thinness(i) = h(i)*mean(colorv(ind));
end
save(file, 'feature_thinness', '-append');

ind=find(features==0);
coord(ind,:)=NaN;

% Plot and save SDDRrrate2 figure
h=figure;
trisurf(tri,coord(:,1),coord(:,2),coord(:,3),colorf);
axis equal
title([partname ' SDDRrate2'])
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'SDDRrate2.fig']);
close(h);
end
    