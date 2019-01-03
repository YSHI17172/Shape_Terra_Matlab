function []=PlotSDDrrate(filename,flag)
path=strcat('FileRecord\',filename,'.mat');
load(path,'clusters','colorv','features','coord','abc','cg','COGs');
cluster=clusters(size(clusters,1)-1,:);
val=[];

for i=1:max(cluster)
    ind=find(cluster==i);
    if nargin==1
        val(i)=mean(colorv(ind));
    else
        val(i)=max(colorv(ind));
    end
end

mkdir('SDDr Analysis',filename);
for i=1:max(features)
    ind=find(features==i);
    clscls=unique(cluster(ind));
    
    D=coord(ind,:)*abc(i,:)';
    t=-((abc(i,:)*cg(i,:)')*ones(size(D))+D)./norm(abc(i,:))^2;

    maxp=ind(find(t==max(t)),1);
    minp=ind(find(t==min(t)),1);

    cl1=cluster(maxp);
    cl2=cluster(minp);
    val1=clusters(size(clusters,1),maxp);
    val2=clusters(size(clusters,1),minp);
    
    if val1>=val2
        cl=cl1;
    else
        cl=cl2;
    end
    

    D=[];
    t=[];
    xyz=[];
    
    D=COGs(clscls,:,5)*abc(i,:)';
    t=-((abc(i,:)*cg(i,:)')*ones(size(D))+D)./norm(abc(i,:))^2;    
    xyz(:,1)=abc(i,1).*t+cg(i,1);
    xyz(:,2)=abc(i,2).*t+cg(i,2);
    xyz(:,3)=abc(i,3).*t+cg(i,3);
    
    xyz0=xyz(find(clscls==cl),:);
    
    Xs=sqrt(sum((xyz-ones(size(xyz,1),1)*xyz0).^2,2));
    vals=val(clscls);
    A=[length(Xs) sum(Xs);sum(Xs) sum(Xs.^2)];
    B=[sum(vals);sum(vals.*Xs')];
    ab=inv(A)*B;
    a=ab(2);
    b=ab(1);
    Xsi=linspace(min(Xs),max(Xs),1000);
    valsi=a*Xsi+b;
    mkdir(strcat('SDDr Analysis\',filename),strcat('Feature ',num2str(i)));
    h=figure;
    scatter(Xs,vals)
    hold
    plot(Xsi,valsi)
    saveas(h,strcat('SDDr Analysis\',filename,'\Feature ',num2str(i),'\graph.fig'));
    h=figure;
    scatter3(coord(ind,1),coord(ind,2),coord(ind,3))
    axis equal
    saveas(h,strcat('SDDr Analysis\',filename,'\Feature ',num2str(i),'\feature.fig'));
    close(h);
end
end
    
    
    
    
    
    

