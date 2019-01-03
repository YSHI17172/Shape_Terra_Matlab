function [conn]=ConnMatTest
path=strcat('ShapeTerra\Output\thinpart\partrecord.mat');
load(path,'tri');

n=max(max(tri));

conn=zeros(n);

for i=1:size(tri,1)

    conn(tri(i,1),tri(i,2))=1;
    conn(tri(i,1),tri(i,3))=1;
    conn(tri(i,2),tri(i,3))=1;
end

conn=conn+conn'+eye(n);
conn=sparse(conn);
end