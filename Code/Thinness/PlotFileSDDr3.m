function []=PlotFileSDDr3(filename)

file=strcat('FileRecord\',filename,'.mat');
load(file,'tri','coord','features','SDDr','colorv','clusters');
mkdir('Plots + SDDr');
h=figure;
n=9;
colorv2=zeros(size(colorv));

cluster=clusters(7,:);

for i=1:max(cluster)
    ind=find(cluster==i);
    val=mean(colorv(ind));
    colorv2(ind)=val;
end


    trisurf(tri,coord(:,1),coord(:,2),coord(:,3),colorv2,'EdgeAlpha',0)
    colormap gray
    alpha('color')
    %trisurf(tri,coord(:,1),coord(:,2),coord(:,3),'FaceAlpha','flat','AlphaDataMapping','scaled','AlphaData',colorv2,'FaceColor','black');
    axis equal
    title(filename)
    saveas(h,strcat('Plots + SDDr\',filename,' ','SDDr2','.fig'));    

close(h);
end