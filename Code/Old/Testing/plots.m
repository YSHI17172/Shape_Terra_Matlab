function plots
path=strcat('ShapeTerra\Output\Old\FileRecord\thinpart.mat');
% First check if connection matrix exists
load(path,'coord','conn');

figure()
hold on
grid on
plot3(coord(:,1),coord(:,2),coord(:,3),'.k');
redpoints1 = find(conn(1,:));
plot3(coord(redpoints1,1),coord(redpoints1,2),coord(redpoints1,3),'.r')
title('Mesh point plot')
hold off
end

