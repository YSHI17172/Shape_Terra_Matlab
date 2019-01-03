function HKS = GenerateHKS_SB(filename,stepin,iterin)

if nargin==3
    step=stepin;
    iter=iterin;
else
    step=0.001;
    iter=1000;
end

%--> Check OS (10/6/2014 by Stephen Baek)
OS = computer;
if strcmp(OS, 'PCWIN') || strcmp(OS, 'PCWIN64')
    file=strcat('FileRecord\',filename,'.mat')
else % GLNX86 || GLNXA64 || MACI64 (wasn't tested on Linux)    
    file=strcat('FileRecord/',filename,'.mat')
end

load(file,'HKS','stepold','iterold');
if ~exist('HKS')
    disp('No HKS found. Generating HKS')
    load(file,'coord','tri');
    nopts=size(coord,1);
    notri=size(tri,1);
    disp('Computing Laplacian')
    tic
    [V,D,L] = Laplacian(coord,tri,nopts,notri);
    toc
    
    t =  step;                                                                
    disp('Iterating...')
    tic
    HKS = zeros(size(V,1), iter);
    
    diagD = diag(D);
    
    for i = 1:iter
        if i==1
            tic
        end
        if rem(i,25)==0
            disp(strcat('Iteration',' ',num2str(i)))
            toc
            tic
        end
        H = (V.^2)*exp(diagD*t);
        
        HKS(:, i) = H;                                                        

        t = t + step;
        
        if i==iter
            toc
        end
    end 
    stepold=step;
    iterold=iter;
    save(file,'HKS','stepold','iterold','-append');
    disp('Total HKS computation time: ')
    toc
elseif ~((step==stepold)&&(iter==iterold)) 
    disp('HKS found with different parameters. Generating HKS for new parameters') 
    load(file,'coord','tri'); 
    nopts=size(coord,1);
    notri=size(tri,1);
    disp('Computing Laplacian')
    tic
    [V,D,L] = Laplacian(coord,tri,nopts,notri); 
    toc
    t =  step;   
    HKS = zeros(size(V,1), iter);  
    diagD = diag(D);
                                                               
    disp('Iterating...')
    tic
    for i = 1:iter
        if i==1
            disp(strcat('Iteration',' ',num2str(i)))
            tic
        end 
        
        if rem(i,25)==0     
            disp(strcat('Iteration',' ',num2str(i)))
            toc
            tic     
        end
        
        H = (V.^2)*exp(diagD*t);
        HKS(:, i) = H;                                                                

        t = t + step;
        
        if i==iter
            toc
        end       
    end
    stepold=step;
    iterold=iter;
    save(file,'HKS','stepold','iterold','-append');
    disp('Total HKS computation time: ')
    toc
     
else
    disp('HKS found with matching parameters')
end

end

function [V,Dee,L]=Laplacian(coord,tri,nopts,notri)

%[V,Dee,L] = CotanLaplacian5(coord,tri,nopts,notri);
[V,Dee,L] = CotanLaplacian_SB(coord,tri);

end

%--> New way of calculating cotan Laplacian (10/6/2014 by Stephen Baek)
function [V,D,L] = CotanLaplacian_SB(coord,tri)
eigno = 300;
    i1 = tri(:,1); i2 = tri(:,2); i3 = tri(:,3);
    v1 = coord(i3,:) - coord(i2,:);  v2 = coord(i1,:) - coord(i3,:); v3 = coord(i2,:) - coord(i1,:);

    n  = cross(v1,v2,2); 
    dblA = (sqrt(sum((n').^2)))';
    totalarea = sum(dblA/2);

    cot12 = -dot(v1,v2,2)./dblA/2; cot23 = -dot(v2,v3,2)./dblA/2; cot31 = -dot(v3,v1,2)./dblA/2;
    diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;

    i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
    j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
    v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];

    L = sparse(i,j,v,size(coord,1),size(coord,1));

%     i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
%     j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
%     offd_v = dblA/24.;
%     diag_v = dblA/12.;
%     v = [offd_v,offd_v, offd_v,offd_v, offd_v,offd_v, diag_v,diag_v,diag_v];  
    i = [i1 i2 i3];
    j = [i1 i2 i3];
    diag_v = dblA/6.;
    v = [diag_v,diag_v,diag_v];
    M = sparse(i,j,v,size(coord,1),size(coord,1));
    
    options=struct('disp',0);
    [V,D,flag] = eigs(M,L,eigno);
    if flag ~= 0
       errrrrrrror=1
    end
    
    V = real(V);
    D = real(D);
    D = diag(D);
    D = diag((totalarea/4/pi)*ones(size(D))./D);
end

%<--


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