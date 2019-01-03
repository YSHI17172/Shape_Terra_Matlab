function ptft=FindFeatures(folderpath) % Barend
%FindFeatures Finds the features from the tips.
%   FindFeatures(FOLDERPATH) finds features of part with data in 
%   FOLDERPATH/partrecord.mat from the tips, plots the different features 
%   and saves as Features.fig in Output\FOLDERPATH

t_start = tic;
% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load clusters and other variables for feature identification
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'clusters','tri','v','clscls','coord','plnsrf','flclsind');

% flclsind=size(clusters,1)/2;
cluster=clusters(2*flclsind-1,:);

val=v(:,2)';
n=size(clusters,2)/2;
I = [tri(:,1);tri(:,2);tri(:,3)];
J = [tri(:,2);tri(:,3);tri(:,1)];
conn = sparse(I,J,ones(3*size(tri,1),1));

ptft=zeros(size(conn,1),1);
conn3=conn;

conn3(plnsrf==0,:)=0;
conn3(:,plnsrf==0)=0;

% Identify features and display various levels of screen comments
ScreenComment('Identifying features',['Identifying features ' partname]);
for i=1:max(clscls)
    ScreenComment('',['Identifying ' partname ' feature #' num2str(i)]);
    conn2=conn3;
    del1=find(clscls~=i);
    del2=find(clscls);
    del=intersect(del1,del2);
    delpts=[];
    for j=1:length(del)
        delpts=[delpts find(cluster==del(j))];
    end
        
    conn2(delpts,:)=0;
    conn2(:,delpts)=0;
    
    inclus=find(clscls==i);
    
    pts=[];
    for j=1:length(inclus)
        pts=[pts find(cluster==inclus(j))];
    end
    
    ptft(pts)=i;
    list = pts;
    minval = val(list(1));
    rep=[];
    while ~isempty(list)
        if ~ismember(list(1),rep)
            adjpts = find(conn2(list(1),:));
            conn2(:,list(1)) = 0;
            conn2(list(1),:) = 0;
            addpts1 = find(val<=minval);
            addpts=intersect(addpts1,adjpts);
            ptft(addpts)=i;
            list = [list,addpts];
            rep=[rep list(1)];
        end
        list(1) = [];
        if ~isempty(list)
            minval=val(list(1));
        end
        
    end
        
    ind=find(ptft==i);
    conn3(ind,:)=0;
    conn3(:,ind)=0;
    clear ind
end

features=ptft;

% Start to refine features, screencomment on what is going on
ScreenComment('Refining features',['Refining features ' partname])

% begin code SB
t_refine = tic;

% Detect high curvature points and divide the whole part into several
% segments by using those high curvature points as a cutting boundary.
label = segmentByCurvature(coord, tri);
nl = max(label);

% Build halfedge structure for an efficient computation
M = HalfEdge(coord, tri);

% For each segments, see if the vertices therin are a part of features.
% If more than 30 percent of vertices in a segment are belonging to a
% certain feature, then merge every other vertex in the same segment to that
% feature.
newfeatures = zeros(size(features));
for i=1:nl
    vtx_in_label = M.facet(find(label==i),:);
    vtx_in_label = unique(vtx_in_label(:));
    
    [feature_number, ia, ic] = unique(sort(features(vtx_in_label)));
    feature_vote = [ia(2:length(ia)); length(vtx_in_label)] - ia;
    [maxval max_id] = max(feature_vote);
    fn = feature_number(max_id);
    if fn == 0 && length(feature_vote)~=1
        feature_vote(1) = 0;
        [maxval max_id] = max(feature_vote);
        if maxval/length(vtx_in_label) > 0.3
            fn = feature_number(max_id);
        end
    end
    newfeatures(vtx_in_label) = fn;
end
% For some cases, small sized features are discarded (or merged into the
% features) during the refinement. For those cases, remap the feature
% indices in a right order.
%   ex) Feature idx before refinement = {0, 1, 2, 3, 4}
%       After discarding feature #2 = {0, 1, 3, 4}
%       After remapping the indices = {0, 1, 2, 3} (3-->2, 4-->3)
nzidx = find(newfeatures~=0);
un = unique(newfeatures(nzidx));
label_map(un) = 1:length(un);
newfeatures(nzidx) = label_map(newfeatures(nzidx));

% Update features with the new ones.
features = newfeatures;

t_elapsed = toc(t_refine);
ScreenComment('',['Features of ' partname ' refined in ' ...
    num2str(t_elapsed) '[s]']);
% end code SB

% Save features to partrecord.mat in Output folder
save(file,'features','-append');

% ScreenComment on plotting tasks
ScreenComment('','Generating features- and left-over-plots');

% Plot and save features figure
h=figure;
coord2=coord;
coord2(ptft==0,:)=NaN;
trisurf(tri,coord2(:,1),coord2(:,2),coord2(:,3),clusters(size(clusters,1),:));
axis equal 
title(strcat(partname,' features'))
colorbar
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'Features.fig']);
close(h);

% Create and save left-over plot
h=figure;
i=size(clusters,1)/2;
coord(features>0,:)=NaN;
trisurf(tri,coord(:,1),coord(:,2),coord(:,3),clusters(2*i,:));
axis equal
title([partname ' left-overs from features'])
colorbar
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'LeftOvers.fig']);    
close(h);

t_elapsed = toc(t_start);
% Display screen comment level two
% Level 2 screen comment when done
ScreenComment('',['All features of ' partname ' identified and plotted in ' ...
    num2str(t_elapsed) '[s]']);
end


% Function to detect planar surfaces and label them
% It will separate curved surfaces from planar surfaces
function label = segmentByCurvature(coord, tri)
    M = HalfEdge(coord, tri);
    
    nv = size(M.vertex,1);
    ne = size(M.edge,1);
    nf = size(M.facet,1);
    
    L = Laplacian(coord, tri);
    MeanCurv = L*coord;
    MeanCurv = sqrt(sum(MeanCurv.^2,2));

    % label = kmeans(MeanCurv, [max(MeanCurv);min(MeanCurv)]);
    % label = label-1;
%     label = MeanCurv < 0.01;
%     label = (1:n)'.*label;
    bndry_vtx = MeanCurv > 0.01;
    bndry_edge = bndry_vtx(M.edge(:,1)) & bndry_vtx(M.edge(:,2));
    
    nei = cell(nf,1);
    for i=1:nf
        fh = M.facet_child_halfedge(i);
        start_fh = fh;
        
        pf = M.halfedge_parent_facet(M.halfedge_opposite_halfedge(fh));
        if pf == 0 || bndry_edge(M.halfedge_parent_edge(fh))==1
            pf = [];
        end
        nei{i} = [nei{i}; pf];
        fh = M.halfedge_next_halfedge(fh);
        while fh ~= start_fh
            pf = M.halfedge_parent_facet(M.halfedge_opposite_halfedge(fh));
            if pf == 0 || bndry_edge(M.halfedge_parent_edge(fh))==1
                pf = [];
            end
            nei{i} = [nei{i}; pf];
            fh = M.halfedge_next_halfedge(fh);
        end
    end
    
    label = (1:nf)';
    old_label = zeros(size(label));
    while sum(label~=old_label)
        old_label = label;
        label_map = (1:nf)';
        for i = 1:nf
            nei_label = label([i; nei{i}]);
            min_nei_label = min(label_map(nei_label));
            label_map(nei_label) = min_nei_label;
        end
        label = label_map(label);
    end
    un = unique(label(find(label~=0)));
    label_map(un) = 1:length(un);
    label = label_map(label);

    [sorted_label ii] = sort(label);
    [unique_label ia ic] = unique(sorted_label);
    label_cnt = [ia(2:length(ia)); length(label)+1] - ia;
    
    one_el_label = find(label_cnt == 1);
    
    label_map = ones(size(unique_label));
    label_map(one_el_label) = 0;
    bw_label = label_map(label);
    
    i = find(bw_label == 0);
    while ~isempty(i)
        mask = FloodFill(M, i(1), bw_label);
        if length(mask) > 2
            label(mask) = max(label)+1;
        end
        bw_label(mask) = 1;
        i = find(bw_label == 0);
    end
    un = unique(label(find(label~=0)));
    label_map(un) = 1:length(un);
    label = label_map(label);
    
    
    [sorted_label ii] = sort(label);
    [unique_label ia ic] = unique(sorted_label);
    label_cnt = [ia(2:length(ia)); length(label)+1] - ia;
    
    one_el_label = find(label_cnt == 1);
    
    for i=1:length(one_el_label)
        fid = find(label==one_el_label(i));
        v1 = M.vertex(M.facet(fid,2),:) - M.vertex(M.facet(fid,1),:);
        v2 = M.vertex(M.facet(fid,3),:) - M.vertex(M.facet(fid,1),:);
        my_n = cross(v1, v2);
        my_n = my_n / sqrt(sum(my_n.^2));
        
        fh = M.facet_child_halfedge(fid);
        start_fh = fh;
        
        pf = M.halfedge_parent_facet(M.halfedge_opposite_halfedge(fh));
        v1 = M.vertex(M.facet(pf,2),:) - M.vertex(M.facet(pf,1),:);
        v2 = M.vertex(M.facet(pf,3),:) - M.vertex(M.facet(pf,1),:);
        n = cross(v1, v2);
        n = n / sqrt(sum(n.^2));
        nei_n = n;
        nei_id = pf;
        fh = M.halfedge_next_halfedge(fh);
        
        while fh ~= start_fh
            pf = M.halfedge_parent_facet(M.halfedge_opposite_halfedge(fh));
            v1 = M.vertex(M.facet(pf,2),:) - M.vertex(M.facet(pf,1),:);
            v2 = M.vertex(M.facet(pf,3),:) - M.vertex(M.facet(pf,1),:);
            n = cross(v1, v2);
            n = n / sqrt(sum(n.^2));
            nei_n = [nei_n; n];
            nei_id = [nei_id; pf];
            fh = M.halfedge_next_halfedge(fh);
        end
        
        [maxval max_id] = max(nei_n*my_n');
        label(fid) = label(nei_id(max_id));
    end
    un = unique(label(find(label~=0)));
    label_map(un) = 1:length(un);
    label = label_map(label);
    
end


% Function to propagate a seed on a surface until boundaries are reached
% Used for surface detection by segmentByCurvature
function mask = FloodFill(he_mesh, seed, label)
    n = size(he_mesh.vertex, 1);
    nei = cell(1,n);
    M = he_mesh;
%     for i=1:n
%         oh = M.vertex_outgoing_halfedge(i);
%         start_oh = oh;
%         nei{i} = [nei{i}; M.halfedge(oh,2)];
%         oh = M.halfedge_next_halfedge(M.halfedge_opposite_halfedge(oh));
%         while oh ~= start_oh
%             nei{i} = [nei{i}, M.halfedge(oh,2)];
%             oh = M.halfedge_next_halfedge(M.halfedge_opposite_halfedge(oh));
%         end
%     end
    
    list = seed;
    mask = [];
    visited = zeros(size(M.facet,1),1);
    while ~isempty(list)
        curr = list(1);
        if visited(curr) == 1
            list(1) = [];
            continue;
        end
        mask = [mask; curr];
        visited(curr) = 1;
        list(1) = [];
        
        fh = 0;
        start_fh = M.facet_child_halfedge(curr);
        while start_fh ~= fh
            if fh == 0
                fh = start_fh;
            end
            oh = M.halfedge_opposite_halfedge(fh);
            nf = M.halfedge_parent_facet(oh);
            if nf~=0 && label(nf)==0 && visited(nf)==0
                list = [list; nf];
            end
            fh = M.halfedge_next_halfedge(fh);
        end
    end
end