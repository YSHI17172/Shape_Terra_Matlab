function [mesh_dist,avg_dist]=MeshSize(n)
path=strcat('..\ShapeTerra\Output\Old\FileRecord\MultiheightpocketB.mat');

% First check if connection matrix exists
load(path,'conn','coord');
% if ~exist(conn)
%     conn = ConnMatrix(path);
%     save(path,'conn','-append');
% end

for i = 1:n
    clear conn_pts
    conn_pts = find(conn(i,:));
    for j = 1:length(conn_pts)
        if i~=j
            z = conn_pts(j);
            clear X Y Z
            %disp(['i = ' num2str(i)])
            %disp(['z = ' num2str(z)])
            X = coord(i,1) - coord(z,1);
            %disp(['i = ' num2str(i)])
            %disp(['z = ' num2str(z)])
            Y = coord(i,2) - coord(z,2);
            %disp(['i = ' num2str(i)])
            %disp(['z = ' num2str(z)])
            Z = coord(i,3) - coord(z,3);
            mesh_dist(i,j) = sqrt(X^2+Y^2+Z^2);    
        end
    end
    avg_dist(i) = sum(mesh_dist(i,:))/(length(conn_pts)-1);
end
end

function [conn]=ConnMatrix(path)
load(path,'tri');

n=max(max(tri));

conn=zeros(n);

for i=1:size(tri,1)

    conn(tri(i,1),tri(i,2))=1;
    conn(tri(i,1),tri(i,3))=1;
    conn(tri(i,2),tri(i,3))=1;
end

conn=conn+conn'+eye(n);
end