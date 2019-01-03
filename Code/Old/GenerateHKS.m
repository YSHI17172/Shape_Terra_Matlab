% TODO
% toc display work needed
function HKS = GenerateHKS(filename,stepin,iterin)
%GenerateHKS Generates and saves Heat Kernel Signature. Usually
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

% Determine if number of iterations is given, if not use default time step
% and number of iteraions
if nargin==3
    step=stepin;
    iter=iterin;
else
    step=0.001;
    iter=1000;
end

% Load calculation results for this part
file=strcat('ShapeTerra\Output\',filename,'\partrecord.mat');
load(file,'HKS','stepold','iterold');

% Check if HKS calculation has been done before, if not run HKS calculation
if ~exist('HKS')
    ScreenComment('No HKS found. Generating HKS',...
        ['No HKS found: ' file '.HKS does not exist' 10 'Generating new HKS matrix']);
    load(file,'coord','tri');
    nopts=size(coord,1);
    notri=size(tri,1);
    ScreenComment('','Computing Laplacian');
    tic
    [V,D,L] = Laplacian(coord,tri,nopts,notri);
    toc
    t =  step;                                                                
    ScreenComment('Iterating for HKS matrix..','Iterating for HKS matrix..')
    tic
    for i = 1:iter
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

        t = t + step;
        
        if i==iter
            toc % Put this in a ScreenComment one day
        end
    end 
    HKS=Hrec;
    stepold=step;
    iterold=iter;
    save(file,'HKS','stepold','iterold','-append');
    t_elapsed = toc;
    ScreenComment(['Total HKS computation time: ' num2str(t_elapsed) '[s]'],...
        ['Total HKS computation time: ' num2str(t_elapsed) '[s]']);
elseif ~((step==stepold)||(iter==iterold))
    ScreenComment('Different HKS found, generating new HKS',...
        'HKS found with different step and iter parameters. Generating HKS for new parameters');
    nopts=size(coord,1);
    notri=size(tri,1);
    ScreenComment('','Computing Laplacian');
    HKS=[];
    tic
    [V,D,L] = Laplacian(coord,tri,nopts,notri); 
    toc
    t =  step;                                                                
    ScreenComment('Iterating for HKS matrix..','Iterating for HKS matrix..')
    tic
    for i = 1:iter
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

        t = t + step;
        
        if i==iter
            toc
        end       
    end
    stepold=step;
    iterold=iter;
    HKS=Hrec;
    save(file,'HKS','stepold','iterold','-append');
    t_elapsed = toc;
    ScreenComment(['Total HKS computation time: ' num2str(t_elapsed) '[s]']);
     
else
    ScreenComment('HKS found with matching parameters',...
        'HKS found with matching parameters');
end

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