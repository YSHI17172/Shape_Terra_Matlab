function [conn]=ConnMat(filename)
path=strcat('FileRecord\',filename,'.mat');
load(path,'tri');

n=max(max(tri));

% conn=zeros(n);
v=[];
for i=1:size(tri,1)
    v(1,1)=1;
    v(1,1)=1;
    v(1,1)=1;
%     conn(tri(i,1),tri(i,2))=1;
%     conn(tri(i,1),tri(i,3))=1;
%     conn(tri(i,2),tri(i,3))=1;
end

% conn=conn+conn'+eye(n);
conn=v;
end