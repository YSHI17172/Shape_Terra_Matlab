function h = FeatureFigure(folderpath)

% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load clusters, mesh coordinates, filtered points and surfaces and
% similarity percentages for this part from partrecord.mat
path=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];

    load(path,'coord','tri','features','feature_thinness');
    label = features;
    
    nv = size(coord,1);
    nf = length(unique(sort(label)))-1;
    adata = ones(size(coord,1), nf+1);
    adata(find(label==0),1) = 0;
    feature_idx = cell(nf, 1);
    popup_legend = cell(nf+1,1);
    popup_legend{1} = 'All Features';
    for i=1:nf
        popup_legend{i+1} = ['Feature #', num2str(i)];
        adata(find(label~=i), i+1) = 0;
        feature_idx{i} = find(label==i);
    end
    
    
    i1 = tri(:,1); i2 = tri(:,2); i3 = tri(:,3);
    v1 = coord(i3,:) - coord(i2,:);  v2 = coord(i1,:) - coord(i3,:); v3 = coord(i2,:) - coord(i1,:);

    face_normal  = cross(v1,v2,2); 
    face_normal = face_normal./repmat(sqrt(sum(face_normal.^2,2)),1,3);
    
    vtx_normal = full([sparse(i1, 1, face_normal(:,1), nv, 1) sparse(i1, 1, face_normal(:,2), nv, 1) sparse(i1, 1, face_normal(:,3), nv, 1)])+...
                 full([sparse(i2, 1, face_normal(:,1), nv, 1) sparse(i2, 1, face_normal(:,2), nv, 1) sparse(i2, 1, face_normal(:,3), nv, 1)])+...
                 full([sparse(i3, 1, face_normal(:,1), nv, 1) sparse(i3, 1, face_normal(:,2), nv, 1) sparse(i3, 1, face_normal(:,3), nv, 1)]);
    vtx_normal = vtx_normal./repmat(sqrt(sum(vtx_normal.^2,2)),1,3);

                    
    for i=1:nf
        [U, S, V] = svd(vtx_normal(feature_idx{i},:));
        feature_pts = coord(feature_idx{i},:);
        feature_trans = feature_pts*V;
        feature_cg(i,:) = sum(feature_pts)./length(feature_idx{i});
        
        feature_offset = feature_pts - repmat(feature_cg(i,:), size(feature_pts,1), 1);
        
        feature_perp = feature_offset - repmat(sum(feature_offset.*repmat(V(3,:), size(feature_pts,1),1), 2),1,3).*repmat(V(3,:), size(feature_pts,1),1);
        
        pocket = sum(feature_perp.*vtx_normal(feature_idx{i},:),2) < 0;
        feature_pocketness(i) = sum(pocket) > 0.5*length(pocket);
        feature_dir(i,:) = V(3,:);
        
        
        
        feature_dim(i,:) = max(feature_trans) - min(feature_trans);
        chull = convhull(feature_trans(:,1:2));
        feature_cross_section{i} = [feature_trans(chull,1:2), zeros(length(chull),1)]*V';
        
        feature_cross_section_area{i} = 0;
        for j=3:length(chull)
            v1 = feature_trans(chull(j-1),1:2) - feature_trans(chull(1),1:2);
            v2 = feature_trans(chull(j  ),1:2) - feature_trans(chull(1),1:2);
            ca = v1(1)*v2(2)-v1(2)*v2(1);
            
            feature_cross_section_area{i} = feature_cross_section_area{i} + abs(ca)*0.5;
        end
    
        label_string{i+1} = {'FEATURE INFORMATION',...
                             '====================',...
                             '   Type:','','',...
                             '   Section Area:','','',...
                             '   Axis Direction:','','',...
                             '   Feature Dimension:','','',...
                             '   Thinness: ',''};
        if feature_pocketness(i)
            label_string{i+1}{4} = '            Pocket';
        else
            label_string{i+1}{4} = '            Protrusion';
        end
        label_string{i+1}{7} = ['           ', num2str(feature_cross_section_area{i}, '%.6f')];
        label_string{i+1}{10} = ['          ', num2str(feature_dir(i,1), '%.3f'), ', ', num2str(feature_dir(i,2), '%.3f'), ', ', num2str(feature_dir(i,3), '%.3f')];
        label_string{i+1}{13} = ['          ', num2str(feature_dim(i,1), '%.2f'), ' X ', num2str(feature_dim(i,2), '%.2f'), ' X ', num2str(feature_dim(i,3), '%.2f')];
        label_string{i+1}{16} = ['          ', num2str(feature_thinness(i), '%.6f')];
    end
    
    label_string{1} = {'FEATURE SUMMARY',...
                       '====================',...
                       '   Num. of Features:',['          ', num2str(nf)],'',...
                       '   Num. of Protrusions:','','',...
                       '   Num. of Pockets:','','',...
                       '   Num. of Thin Features:','','',...
                       '   Information Here: ',''};
    
    h = figure;
    t = trisurf(tri,coord(:,1),coord(:,2),coord(:,3));
    set(t, 'FaceAlpha', 'interp', 'FaceVertexAlphaData',adata(:,1));
    hPopup = uicontrol('Style', 'popup', 'String', popup_legend);
    hLabel = uicontrol('Style','text','String',label_string{1},'HorizontalAlignment','left');
    
    myhandles = guihandles(h); 
    myhandles.hPopup = hPopup;
    myhandles.hLabel = hLabel;
    myhandles.t = t;
    myhandles.hObj = h;
    myhandles.adata = adata;
    myhandles.label_string = label_string;
    guidata(h, myhandles);
    % set(hPopup
    %            'Position', [20 340 100 50],...
    
    hold on
        
    for i=1:nf
        arrow = [feature_cg(i,:) - 0.7*feature_dir(i,:)*feature_dim(i,3);...
                 feature_cg(i,:) + 0.7*feature_dir(i,:)*feature_dim(i,3)];
        if feature_pocketness(i)
            plot3(arrow(:,1), arrow(:,2), arrow(:,3), '-b','LineWidth',3);
        else
            plot3(arrow(:,1), arrow(:,2), arrow(:,3), '-g','LineWidth',3);
        end
        plot3(feature_cross_section{i}(:,1), feature_cross_section{i}(:,2), feature_cross_section{i}(:,3), '-k', 'LineWidth',3);
    end

    axis equal 
    

    popup_callback = ['data = guidata(gcbo);',...
        't = data.t;',...
        'hPopup = data.hPopup;',...
        'hLabel = data.hLabel;',...
        'hObj = data.hObj;',...
        'adata = data.adata;',...
        'label_string = data.label_string;',...
        'val = get(hPopup,''Value'');',...
        'set(hLabel,''String'', label_string{val});',...
        'set(t, ''FaceAlpha'', ''interp'', ''FaceVertexAlphaData'', adata(:,val));'];

    set(hPopup,'Callback',popup_callback);
    
    
    
    resize_callback = ['data = guidata(gcbo);',...
    'hLabel = data.hLabel;',...
    'hPopup = data.hPopup;',...
    'hObj = data.hObj;',...
    'pos = get(hObj, ''Position'');',...
    'margin = 10/pos(3);',...
    'if(pos(3) > pos(4)) ',...
    'x = 0.5+pos(4)*0.6/pos(3);',...
    'else ',...
    'x = 1-16*margin;',...
    'end;',...
    'if(x > 1-16*margin)',...
    'x = 1-16*margin;',...
    'end;',...
    'y = 0.91-225/pos(4);',...
    'set(hLabel,''Units'',''Normalized'',''Position'',[x y 15*margin 225/pos(4)]);',...
    'set(hPopup,''Units'',''Normalized'', ''Position'',[x-1/pos(3) y+230/pos(4) 15*margin+2/pos(3) 20/pos(4)]);'];
    set(h,'ResizeFcn',resize_callback);
    
    title('Final Result');
    saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'FinalResult.fig']); 
    
end
