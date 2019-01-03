function [coord tri]=ReadFile2(filename)
filetest=strcat('FileRecord\',filename,'.mat');
file=strcat('Input\',filename);
varfile=strcat('FileRecord\',filename,'.mat');
if ~exist(filetest)
    disp('No file found, generating Coord and Tri')
    tic
    s1=strcat(file,'.dat');
    s2=strcat(file,'.off');
    if exist(s1)>0
        [coord tri]=ReadDAT(s1);
    elseif exist(s2)>0
        [coord tri]=ReadOFF(s3);
    end
    save(varfile,'coord','tri');
    toc
else
    load(filetest,'coord','tri');
    if ~(exist('coord') && exist('tri'))
        disp('File found, no Coord nor Tri found, generating Coord and Tri')
        tic
        s1=strcat(file,'.dat');
        s2=strcat(file,'.off');
        if exist(s1)>0
            [coord tri]=ReadDAT(s1);
        elseif exist(s2)>0
            [coord tri]=ReadOFF(s3);
        end
        save(varfile,'coord','tri','-append');
        tic
    else
        disp('File found, Coord and Tri found')
    end
end

end

function [coord tri]=ReadOFF(path)
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

function [coord tri] =ReadDAT(path)
[names]=textread(path,'%s');
cfile=zeros(length(names),1);
for i=1:length(names)
    cfile(i)=str2double(names(i));
end
    
i=18;
j=1;

while ~isnan(cfile(i))

    coord(j,1)=cfile(i+1);
    coord(j,2)=cfile(i+2);
    coord(j,3)=cfile(i+4);
    i=i+6;
    j=j+1;    
end

i=i+14;
l=1;

while ~isnan(cfile(i))

    tri(l,1)=cfile(i+2);
    tri(l,2)=cfile(i+3);
    tri(l,3)=cfile(i+4);
    tri(l,4)=cfile(i+5);
    i=i+7;
    l=l+1;
end
end