function varargout = SliderFigure(varargin)
%% Draw a Figure with a Slider Bar
% h = SliderFigure(coord, tri, cdata)
% h = SliderFigure(coord, tri, cdata, adata)
% h = SliderFigure(coord, tri, cdata, labels)
% h = SliderFigure(coord, tri, cdata, adata, labels)
%
% coord: (NV X 3) matrix containing coordinate values of vertices
% tri:   (NF X 3) matrix containing the triangulation
% cdata: (NV X m) matrix containing color values (or level of physical properties, e.g., HKS, Persistence, ...)
%         cf) m is a number of slider ticks.
% adata: (NV X m) matrix containing the alpha values (=transparancy values: alpha==0 means fully transparent, alpha==1 means not transparent)
% label: (m X 1) cell with strings (description of slider figures)

    narginchk(3,5);         
    coord = varargin{1};
    tri = varargin{2};
    cdata = varargin{3};

    hasAlpha = false;
    hasLabel = false;
    if nargin == 5
        hasAlpha = true;
        hasLabel = true;
        adata = varargin{4};
        label_string = varargin{5};
    elseif nargin == 4
        if isfloat(varargin{4})
            hasAlpha = true;
            adata = varargin{4};
        else
            hasLabel = true;
            label_string = varargin{4};
        end
    end


    % CHECK DATA VALIDITY
    if size(cdata, 1) < size(cdata, 2)
        cdata = cdata';
    end
    if hasAlpha
        if size(adata, 1) < size(adata, 2)
            adata = adata';
        end
    end
    nv = size(coord,1);
    nf = size(tri, 1);
    if size(cdata,1) ~= nv
        error('val must be a matrix of size (N_VERTEX)x(N_DATA) or (N_DATA)x(N_VERTEX)')
    end
    nd = size(cdata,2);
   
%     if hasTitle
%         if length(title_string) ~= nd
%             error('
%         end
%     end
    
    
    h = figure;
    t = trisurf(tri,coord(:,1),coord(:,2),coord(:,3));
    set(t, 'CData', cdata(:,1));
    if hasAlpha
        set(t, 'FaceAlpha', 'interp', 'FaceVertexAlphaData', adata(:,1));
    end
    
    axis equal
    
    colorbar
    
    if hasLabel
        hLabel = uicontrol('Style','text','String',label_string{1});
    else
        hLabel = uicontrol('Style','text','String','Level');
    end
    set(hLabel,'BackgroundColor',get(h,'Color'));
    set(hLabel,'HorizontalAlignment','left');
    
    hSlider = uicontrol('Style', 'slider',...
        'Min',1,'Max',nd,'Value',1,'SliderStep', [1.0/nd,1.0/nd]);
    
    
    
    
    myhandles = guihandles(h); 
    myhandles.hSlider = hSlider;
    myhandles.hLabel = hLabel;
    myhandles.hObj = h;
    myhandles.t = t;
    myhandles.cdata = cdata;
    if hasAlpha
        myhandles.adata = adata;
    end
    if hasLabel
        myhandles.label_string = label_string;
    end
    guidata(h, myhandles);
    
    
    
    
    
    
    if hasAlpha & hasLabel
        slider_callback = ['data = guidata(gcbo);',...
        'hSlider = data.hSlider;',...
        't = data.t;',...
        'cdata = data.cdata;',...
        'adata = data.adata;',...
        'label_string = data.label_string;',...
        'sliderVal = round(get(hSlider,''Value''));',...
        'set(hSlider,''Value'',sliderVal);',...
        'set(t, ''CData'', cdata(:,sliderVal));',...
        'set(t, ''FaceAlpha'', ''interp'', ''FaceVertexAlphaData'', adata(:,sliderVal));',...
        'set(hLabel,''String'', label_string{sliderVal});',...
        'drawnow'];      
        set(hSlider,'Callback',slider_callback);  
%         if all(version('-release') >= '2014a')
%             addlistener(hSlider,'ContinuousValueChange',@(hObject, event) slidercallback_all(hObject,event,t,cdata,adata,hLabel,label_string));
%         else
%             addlistener(hSlider,'ActionEvent',@(hObject, event) slidercallback_all(hObject,event,t,cdata,adata,hLabel,label_string));
%         end
    elseif hasLabel
        slider_callback = ['data = guidata(gcbo);',...
        'hSlider = data.hSlider;',...
        't = data.t;',...
        'cdata = data.cdata;',...
        'label_string = data.label_string;',...
        'sliderVal = round(get(hSlider,''Value''));',...
        'set(hSlider,''Value'',sliderVal);',...
        'set(t, ''CData'', cdata(:,sliderVal));',...
        'set(hLabel,''String'', label_string{sliderVal});',...
        'drawnow'];  
        set(hSlider,'Callback',slider_callback);  
%         if all(version('-release') >= '2014a')
%             addlistener(hSlider,'ContinuousValueChange',@(hObject, event) slidercallback_label(hObject,event,t,cdata,hLabel,label_string));
%         else
%             addlistener(hSlider,'ActionEvent',@(hObject, event) slidercallback_label(hObject,event,t,cdata,hLabel,label_string));
%         end
    elseif hasAlpha
        slider_callback = ['data = guidata(gcbo);',...
        'hSlider = data.hSlider;',...
        't = data.t;',...
        'cdata = data.cdata;',...
        'adata = data.adata;',...
        'sliderVal = round(get(hSlider,''Value''));',...
        'set(hSlider,''Value'',sliderVal);',...
        'set(t, ''CData'', cdata(:,sliderVal));',...
        'set(t, ''FaceAlpha'', ''interp'', ''FaceVertexAlphaData'', adata(:,sliderVal));',...
        'drawnow'];  
        set(hSlider,'Callback',slider_callback); 
%         if all(version('-release') >= '2014a')
%             addlistener(hSlider,'ContinuousValueChange',@(hObject, event) slidercallback_alpha(hObject,event,t,cdata,adata));
%         else
%             addlistener(hSlider,'ActionEvent',@(hObject, event) slidercallback_alpha(hObject,event,t,cdata,adata));
%         end
    else
        slider_callback = ['data = guidata(gcbo);',...
        'hSlider = data.hSlider;',...
        't = data.t;',...
        'cdata = data.cdata;',...
        'sliderVal = round(get(hSlider,''Value''));',...
        'set(hSlider,''Value'',sliderVal);',...
        'set(t, ''CData'', cdata(:,sliderVal));',...
        'drawnow'];
        
%     hhSlider = handle(hSlider);
%     hProp = findprop(hhSlider,'Value');  % a schema.prop object
%     hListener = handle.listener(hhSlider,hProp,'PropertyPreSet',slider_callback);
%     setappdata(hSlider,'sliderListener',hListener);
%         set(hSlider,'Callback',@(src,data) disp(src.UserData * src.Value));
% hSlider.addlistener('Value','PostSet',@(src,data) data.AffectedObject.Callback(data.AffectedObject,struct('Source',data.AffectedObject,'ActionEvent','Action')));


%         hListener = handle.listener(hSlider,'ActionEvent',slider_callback);
%         setappdata(hSlider,'sliderListener',hListener);  % this is important - read below


%         if all(version('-release') >= '2014a')
            set(hSlider,'Callback',slider_callback);
%             addlistener(hSlider,'ContinuousValueChange',@(hObject, event) slidercallback(hObject,event,t,cdata));
%             addlistener(hSlider,'ContinuousValueChange',@(hObject, event) guidata(gcbo));
%         else
%             addlistener(hSlider,'ActionEvent',@(hObject, event) slidercallback(hObject,event,t,cdata));
%         end
    end
%     set(h, 'ResizeFcn', {@figureResize});
    
    resize_callback = ['data = guidata(gcbo);',...
    'hSlider = data.hSlider;',...
    'hLabel = data.hLabel;',...
    'hObj = data.hObj;',...
    'pos = get(hObj, ''Position'');',...
    'margin = 10/pos(3);',...
    'set(hSlider,''Units'',''Normalized'',''Position'',[margin margin 1.0-2*margin 20.0/pos(4)]);',...
    'set(hLabel,''Units'',''Normalized'',''Position'',[margin 20.0/pos(4)+margin 1 20.0/pos(4)]);'];
    set(h,'ResizeFcn',resize_callback);
    
    varargout{1} = h;
    
    
%     function slidercallback(hObj,event,t,cdata)
%         sliderVal = round(get(hObj,'Value'));
%         set(hObj,'Value',sliderVal);
%         set(t, 'CData', cdata(:,sliderVal));
%     end
end

% function slidercallback_alpha(hObj,event,t,cdata,adata)
%     sliderVal = round(get(hObj,'Value'));
%     set(hObj,'Value',sliderVal);
%     set(t, 'CData', cdata(:,sliderVal));
%     set(t, 'FaceAlpha', 'interp', 'FaceVertexAlphaData', adata(:,sliderVal));
% end
% function slidercallback_label(hObj,event,t,cdata,hLabel,label_string)
%     sliderVal = round(get(hObj,'Value'));
%     set(hObj,'Value',sliderVal);
%     set(t, 'CData', cdata(:,sliderVal));
%     set(hLabel,'String',label_string{sliderVal});
% end
% function slidercallback_all(hObj,event,t,cdata,adata,hLabel,label_string)
%     sliderVal = round(get(hObj,'Value'));
%     set(hObj,'Value',sliderVal);
%     set(t, 'CData', cdata(:,sliderVal));
%     set(t, 'FaceAlpha', 'interp', 'FaceVertexAlphaData', adata(:,sliderVal));
%     set(hLabel,'String',label_string{sliderVal});
% end