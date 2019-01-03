% TODO
% toc display work needed
function HKS = GenerateHKS_old(folderpath,stepin,iterin)
%GenerateHKS_V1 Generates and saves Heat Kernel Signature. Usually
%   HKS = GenerateHKS(FILENAME,STEPSIZE,#ITERATIONS) generates the HKS 
%   of the  structural part under investigation. HKS is a matrix of size
%   [#points x #iterations] and consists of point heat values per point timestep
%
%   FILENAME is the filemesh path for which calculations are performed
%   STEPSIZE is the timestep value in seconds
%   #ITERATIONS is the number of timesteps that are calculated
%   If step STEPSIZE and #ITERATIONS not specified use dfault values
%   STEPSIZE = 0.001,  #ITERATIONS = 1000
%
%   Further output of GenerateHKS is stored in FILENAME.mat under 
%   variable names: STEPOLD ITEROLD

tic
% Retreive filename from folderpath
indsep = strfind(folderpath,'\');
partname = folderpath(1:indsep(1)-1);
meshfolder = folderpath(indsep(1)+1:indsep(2)-1);
varfile=strcat('ShapeTerra\Output\',folderpath,'\partrecord.mat');

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
        ['No HKS found: ShapeTerra\Output\' partname '\meshfolder\ does not contain HKS data' 10 'Generating new HKS matrix']);
    load(varfile,'coord','tri');
    nopts=size(coord,1);
    notri=size(tri,1);
    ScreenComment('','Computing Laplacian');
    tic
    [V,D,L] = Laplacian(coord,tri,nopts,notri);
    toc
    t =  steps;                                                                
    ScreenComment('Iterating for HKS matrix..','Iterating for HKS matrix..')
    tic
    for i = 1:iters
        if i==1
            tic
        end
        if rem(i,25)==0
            ScreenComment('',['Iteration ' num2str(i)]);
            toc
            tic
        end
        H = V*expm(D*t)*V';                                                   
        H = real(diag(H));                                                    
        if (i == 1)                                                           
            Hrec = H;                                                         
        else                                                                  
            Hrec = [Hrec,H];                                                  
        end                                                                   

        t = t + steps;
        
        if i==iters
            toc % Put this in a ScreenComment one day
        end
    end 
    HKS=Hrec;
    stepold=steps;
    iterold=iters;
    save(varfile,'HKS','stepold','iterold','-append');
    t_elapsed = toc;
    ScreenComment(['Total HKS computation time: ' num2str(t_elapsed) '[s]'],...
        ['Total HKS computation time: ' num2str(t_elapsed) '[s]']);
else % HKS file(s) for PARTNAME found in one of the subsubfolders of output
    if isempty(HKS_match)
        ScreenComment('Different HKS found, generating new HKS',...
        ['HKS found with different step and iter parameters. Generating HKS'...
        ' for new parameters']);
        
        load(varfile,'coord','tri');
        nopts=size(coord,1);
        notri=size(tri,1);
        ScreenComment('','Computing Laplacian');
        HKS=[];
        tic
        [V,D,L] = Laplacian(coord,tri,nopts,notri); 
        toc
        t =  steps;                                                                
        ScreenComment('Iterating for HKS matrix..','Iterating for HKS matrix..')
        tic
        for i = 1:iters
            if i==1
                ScreenComment('',['Iteration ' num2str(i)]);
                tic
            end 

            if rem(i,25)==0     
                ScreenComment('',['Iteration ' num2str(i)]);
                toc % Put toc output in ScreenComment or store as variable that is not displayed
                tic     
            end

            H = V*expm(D*t)*V';                                                   
            H = real(diag(H));                                                    
            if (i == 1)                                                           
                Hrec = H;                                                         
            else                                                                  
                Hrec = [Hrec,H];                                                  
            end                                                                   

            t = t + steps;

            if i==iters
                toc
            end       
        end
        stepold=steps;
        iterold=iters;
        HKS=Hrec;
        t_elapsed = toc;
        ScreenComment(['Total HKS computation time: ' num2str(t_elapsed) '[s]']);
    else
        ScreenComment('HKS found with matching parameters',...
        'HKS found with matching parameters');
        % Load found HKS iter and step data from most recent same partrecord.mat 
        partdata=['ShapeTerra\Output\' partname '\' meshfolder '\' ...
            HKS_match{end} '\partrecord.mat'];
        load(partdata,'HKS','stepold','iterold');
    end
    % Save (new) HKS results to partrecord.mat in folderpath
    save(varfile,'HKS','stepold','iterold','-append');
end
end

% Check if same HKS exists in folderpath
function [HKS_match,HKS_flag] = CheckHKS(folderpath,steps,iters)
% Retreive partname and meshfolder from folderpath
indsep = strfind(folderpath,'\');
partname = folderpath(1:indsep(1)-1);
meshfolder = folderpath(indsep(1)+1:indsep(2)-1);

% Define folder where you want so search timestamp subfolders in
FoInfo = dir(['ShapeTerra\Output\' partname '\' meshfolder]);
FoInfo(~[FoInfo.isdir]) = []; % Remove non-folder entries in folder info
RemIdx(1) = find(strcmp({FoInfo.name},{'.'})==1);
RemIdx(2) = find(strcmp({FoInfo.name},{'..'})==1);
FoInfo(sort(RemIdx)) = []; % Remove . and .. entries from folder info

% Search for mesh_match in Output\partname\meshsize timestamp subfolders
HKS_match = {''}; % Create empty cell array to store names of 
                   % subsubfolders with similar meshes in it
HKS_flag = zeros(1,length(FoInfo));
for FoIdx = 1:length(FoInfo)
    % Clear HKS from previous loop iteration if it exists
    if exist('HKS')
        clear HKS
    end
    filetest=['ShapeTerra\Output\' partname '\' meshfolder '\' ...
        FoInfo(FoIdx).name '\partrecord.mat'];
    %load(filetest,'HKS');
    % Check if HKS can be found in partrecord.mat
    HKStest = whos('-file',filetest,'-regexp','HKS');
    clear stepold iterold
    if ~isempty(HKStest)
       % Clear iterold and stepold from previous loop iteration if they exist
        if exist('stepold') || exist('iterold')
            clear stepold iterold
        end

        % HKS exists
        HKS_flag(FoIdx) = 1;
        load(filetest,'HKS','stepold','iterold');
        
        % Now check if number of iterations and stepsize HKS are the same
        if ~((steps==stepold) || (iters==iterold))
            % different HKS found
        else
            % HKS found with same iter and step
            HKS_match(FoIdx) = {FoInfo(FoIdx).name};
        end
    else % HKS does not exist
        HKS_flag(FoIdx) = 0;
    end
end

% Remove empty entries from mesh_match and sort
HKS_match(cellfun(@isempty,HKS_match)) = [];
HKS_match = sort(HKS_match);
end

function [V,Dee,L]=Laplacian(coord,tri,nopts,notri)

[V,Dee,L] = CotanLaplacian5(coord,tri,nopts,notri);

end

function [V,D,L] = CotanLaplacian5(coord,tri,nopts,notri)
eigno = 300;  
    v1 = coord(tri(:,1),:)-coord(tri(:,2),:);
    v2 = coord(tri(:,1),:)-coord(tri(:,3),:);
    v3 = coord(tri(:,2),:)-coord(tri(:,3),:);   
    d1=sqrt(sum((v1).^2,2));
    d2=sqrt(sum((v2).^2,2));
    d3=sqrt(sum((v3).^2,2));
    d = [d2;d3;d1];
    
    adj = zeros(nopts,nopts);
    dt = zeros(nopts,1);
    for i = 1:nopts
        temp = find(tri(:,1)==i|tri(:,2)==i|tri(:,3)==i);
        dtot = 0;
        for j = 1:length(temp)
            index = find(tri(temp(j,1),:)==i);
            if (index == 1), dtot = dtot + d1(temp(j,1),1) + d2(temp(j,1),1);, end;
            if (index == 2), dtot = dtot + d1(temp(j,1),1) + d3(temp(j,1),1);, end;
            if (index == 3), dtot = dtot + d2(temp(j,1),1) + d3(temp(j,1),1);, end;
        end
        dt(i,1) = dtot/(2*length(temp));
    end
    dtotal = [dt(tri(:,1),1);dt(tri(:,3),1);dt(tri(:,2),1)];

    angle1=acos(dot(v1,v2,2)./(d1.*d2));
    angle2=acos(dot(v1,-v3,2)./(d1.*d3));
    angle3=pi-angle1-angle2;
    angles=[angle2;angle1;angle3];

   
    TR1=tri(:,[1 3 2]); 
    TR2=tri(:,[3 2 1]); 
    L=sparse(double(TR1(:)),double(TR2(:)),cot(angles).*(1./d).*dtotal,nopts,nopts);
    L=(L+L')./2;
    D=sum(L,2);
    D(D==0)=1;
    L=L-spdiags(D,0,nopts,nopts);
    
    R=d1./(2*sin(angle3));

    RA1=(d1./8).*sqrt(4*R.^2-d1.^2);
    RA2=(d3./8).*sqrt(4*R.^2-d3.^2);
    RA3=(d2./8).*sqrt(4*R.^2-d2.^2);
    Areas=1/2.*sqrt( (d1.^2).*(d3.^2)-((d1.^2+d3.^2-d2.^2)./2).^2);
    totalarea=sum(Areas);

    RA1(angle1>pi/2)=1/4.*Areas(angle1>pi/2);    
    RA2(angle1>pi/2)=0;
    RA3(angle1>pi/2)=RA1(angle1>pi/2);
    
    RA2(angle2>pi/2)=1/4.*Areas(angle2>pi/2);
    RA3(angle2>pi/2)=0;
    RA1(angle2>pi/2)=RA2(angle2>pi/2);
    
    RA3(angle3>pi/2)=1/4.*Areas(angle3>pi/2);    
    RA1(angle3>pi/2)=0;
    RA2(angle3>pi/2)=RA3(angle3>pi/2);
    
    RA=[RA1;RA2;RA3];
    RA=[RA;RA];
    TR1=tri(:,[1 2 3 2 3 1]); 
    TR2=tri(:,[2 3 1 3 1 2]);

    AM=sparse(double(TR1(:)),double(TR2(:)),RA,nopts,nopts);    
    AM=full(sum(AM,2));
    AM(AM<1e-6)=mean(AM(AM>0));
       
M=spdiags(AM,0,nopts,nopts);

options=struct('disp',0);
[V,D,flag] = eigs(L,4*pi*M/totalarea,eigno,0.1,options);
if flag ~= 0
   errrrrrrror=1
end
end