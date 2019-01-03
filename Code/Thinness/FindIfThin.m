function [SDDr]=FindIfThin(folderpath)
%FindIfThin Find thin features.
%   [SDDR] = FindIfThin(FOLDERPATH) finds thin features of part with 
%   associated data in FOLDERPATH\partrecord.mat

t_start = tic;
% Load global file separator fs
global fs

% Load coord, tri, features and other for idenfication thin features
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
%load(file,'coord','tri','features','conn','SDDr');
load(file,'coord','tri','features');

% Check if conn can be found in partrecord.mat
conntest = whos('-file',file,'-regexp','conn');

% ScreenComment on what is happening
DispStr = 'Finding thin features';
ScreenComment(DispStr,[DispStr ' ' partname '..'])
% if ~exist('SDDr')
    colorv=zeros(size(features));
    %tic
    %if ~exist('conn')
    if isempty(conntest)
        conn=ConnectionMatrix(tri);
        save(file,'conn','-append');
    end   
    %toc
    n=max(features);
    SDDr=zeros(2,n);
    MinD=zeros(2,n);
    
    for i=1:n
        clear t_elapsed
        t_iter = tic;
        feati=find(features==i);
        
        coord2=coord(feati,:);
        conn2=conn(feati,feati);
        
        dist=FindDistance(coord2);
        dist2=dist+eye(size(dist));
        
        SCD=FindSCD(conn2);
      
        SDDri=SCD./dist2-eye(size(dist));
        colorvi=max(SDDri,[],2);
        rmax=max(max(SDDri));
        [a b]=find(SDDri==rmax,1,'first');
        SDDr(1,i)=feati(a);   
        SDDr(2,i)=feati(b);
        dist2=dist+max(max(dist))*eye(size(dist));
        
        mindi=min(min(dist2));
        [c d]=find(dist2==mindi,1,'first');
        MinD(1,i)=feati(c);   
        MinD(2,i)=feati(d);
        
        colorv(feati)=colorvi;
        t_elapsed = toc(t_iter);
        ScreenComment('',['Elapsed time finding thin feature #' num2str(i)...
            ': ' num2str(t_elapsed) '[s]']);
    end
    % Save results to partrecord.mat
    save(file,'SDDr','MinD','colorv','-append');

    t_total = toc(t_start);
    % Display screen comment level two
    ScreenComment('',['Elapsed total time finding thin features: '...
        num2str(t_total) '[s]']);
end

function [dist]=FindDistance(coord)
dist=[];
H=ones(size(coord,1),1);

for i=1:size(coord,1)
    x=coord(i,1);
    y=coord(i,2);
    z=coord(i,3);
    d=sum((coord-[x.*H y.*H z.*H]).^2,2).^0.5;
    dist=[dist d];
end
end


function [D]=FindSCD(conn2)
B=conn2;
C=B;
E=C;

D=ones(size(conn2));
i=1;
while any(any(E==0))
    if i==size(conn2,1)
        break
    end
    E=C;
    D=D+(C==0);
    C=C*B;
    i=i+1;
end
end

% function [conn]=ConnectionMatrix(tri)
% n=max(max(tri));
% conn=zeros(n);
% for i=1:size(tri,1)
%     conn(tri(i,1),tri(i,2))=1;
%     conn(tri(i,1),tri(i,3))=1;
%     conn(tri(i,2),tri(i,3))=1;
% end
% conn=conn+conn'+eye(n);
% end

function [conn]=ConnectionMatrix(tri)
sptri=[tri(:,[1 2]);tri(:,[2 3]);tri(:,[3 1])];
v=[1:max(max(tri))]';
sptri=[sptri; [v v]];
conn=sparse(sptri(:,1),sptri(:,2),ones(size(sptri,1),1));
end