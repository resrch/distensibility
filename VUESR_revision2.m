% changed: shearpoints from 15 to total in wall shear function
classdef VUESR_revision2 < handle

    %This is a dicom image slice viewer with built in scroll, contrast, zoom,
    %,ROI, and image measurement tools.
    %
    %   Use this class to place a self-contained image viewing panel within
    %   a GUI (or any figure). Similar to imtool but with slice scrolling.
    %   Only designed to view grayscale (intensity) images. Use the mouse
    %   to control how the image is displayed. A left click allows window
    %   ad leveling, a right click is for panning, and a middle click is
    %   for zooming. Also the scroll wheel can be used to scroll through
    %   slices.
    %----------------------------------------------------------------------
    %Inputs:
    %
    %   I           An m x n x k image array of grayscale values. Default
    %               is a 100x100x3 random noise image.
    %   position    The position of the panel containing the image and all
    %               the tools. Format is [xmin ymin width height]. Default
    %               position is [0 0 1 1] (units = normalized). See the
    %               setPostion and setUnits methods to change the postion
    %               or units.
    %   h           Handle of the parent figure. If no handles is provided,
    %               a new figure will be created.
    %   range       The display range of the image. Format is [min max].
    %               The range can be adjusted with the contrast tool or
    %               with the setRange method. Default is [min(I) max(I)].
    %----------------------------------------------------------------------
    %Output:
    %
    %   tool        The VUESR object. Use this object as input to the
    %               class methods described below.
    %----------------------------------------------------------------------
    %Constructor Syntax
    %
    %tool = VUESR() creates an VUESR panel in the current figure with
    %a random noise image. Returns the VUESR object.
    %
    %tool = VUESR(I) sets the image of the VUESR panel.
    %
    %tool = VUESR(I,position) sets the position of the VUESR panel
    %within the current figure. The default units are normalized.
    %
    %tool = VUESR(I,position,h) puts the VUESR panel in the figure
    %specified by the handle h.
    %
    %tool = VUESR(I,position,h,range) sets the display range of the
    %image according to range=[min max].
    %
    %Note that you can pass an empty matrix for any input variable to have
    %the constructor use default values. ex. tool=VUESR([],[],h,[]).
    %----------------------------------------------------------------------
    %Methods:
    %
    %   setImage(tool, I) displays a new image.
    %
    %   I = getimage(tool) returns the image being shown by the tool
    %
    %   setPostion(tool,position) sets the position of tool.
    %
    %   position = getPosition(tool) returns the position of the tool
    %   relative to its parent figure.
    %
    %   setUnits(tool,Units) sets the units of the position of tool. See
    %   uipanel properties for possible unit strings.
    %
    %   units = getUnits(tool) returns the units of used for the position
    %   of the tool.
    %
    %   handles = getHandles(tool) returns a structured variable, handles,
    %   which contains all the handles to the various objects used by
    %   VUESR.
    %
    %   setDisplayRange(tool,range) sets the display range of the image.
    %   see the 'Clim' property of an Axes object for details.
    %
    %   range=getDisplayRange(tool) returns the current display range of
    %   the image.
    %
    %   setWindowLevel(tool,W,L) sets the display range of the image in
    %   terms of its window (diff(range)) and level (mean(range)).
    %
    %   [W,L] = getWindowLevel(tool) returns the display range of the image
    %   in terms of its window (W) and level (L)
    %
    %   ROI = getcurrentROI(tool) returns info about the currently selected
    %   region of interest (ROI). If no ROI is currently selected, the
    %   method returns an empty matrix. ROI is a structured variable with
    %   the following fields:
    %       -ROI.mask is a binary mask that defines the pixels within the
    %       ROI.
    %       -ROI.stats is a structured variable containing stats about the
    %       ROI. Included stats are, Area, Perimeter, MaxIntensity,
    %       MinIntensity, MeanIntensity, and STD.
    %
    %   setCurrentSlice(tool,slice) sets the current displayed slice.
    %
    %   slice = getCurrentSlice(tool) returns the currently displayed
    %   slice.
    %
    %----------------------------------------------------------------------
    %Notes:
    %
    %   Based on imtool3D v2.1 by Justin Solomon
    %   Modified by Barry Belmont and Luc Hildebrand to different ends
    %
    %   Requires the image processing toolbox
    
    properties (SetAccess = private, GetAccess = private)
        I               % Image data (MxNxK) matrix of image data (double)
        handles         % Structured variable with all the handles
        handlesROI      % list of ROI handles
        currentROI      % Currently selected ROI
        centers         % list of bin centers for histogram
        calibration     % Converts pixels to mm
        pointLog        % Log of all the points being tracked
        pointLog2       % Log points for the 2 point tracker
        pointLogScores  % Validity scores for first point log
        pixelDensity    % Amount of pixels being tracked
        accFrames       % Frames used to find accumulated strain
        map             % Mean Arterial blood pressure
        fRate           % Frame rate in Hz
        hRate           % Heart rate in BPM
        fName           % File opened
        zlocation       % Measurement location on the arm
        zdiameter       % Average measured diameter
        zdistensibility % Average measured distensibility
        zelasticity     % Average measured elasticity
        zpatientID
        pointDistCm
        studydate
        
    end
    
    methods
        
        function tool = VUESR_revision2(varargin)  %Constructor
            %%
            %Check the inputs and set things appropriately
            tool.fName = [];
            switch nargin
                case 0  %tool = VUESR()
                    [fileName, filePath] = uigetfile('*.DCM;*.dcm;*.mat;*', ...
                        'Choose DICOM images to import', pwd, ...
                        'MultiSelect', 'off');
                    if filePath(1) == 0
                        disp('No files chosen, exiting function')
                        return;
                    else
                        fileName = fullfile(filePath, fileName);       %WDR fixed
                        disp(['User selected: ', fullfile(fileName)]);
                        [~, ~, ext] = fileparts(fileName);
                        if strcmp(ext,'.DCM') || strcmp(ext,'.dcm') || strcmp(ext,'')
                            tool.I = permute(dicomread(fileName),[1, 2, 4, 3]);
                            myslice = tool.I(:,:,10,1);
                        else
                            load( fileName);
                            tool.I = permute(image_change,[1 2 4 3]);
                        end
                        position=[0, 0, 1, 1];
                        heightHistogram=  figure('Position', [400 200 650 600],'Name','VUESR Imaging Toolbox - alpha','NumberTitle','off');
                        set(heightHistogram,'Toolbar','none','Menubar','none');                
                        pixelValueRange = [min(tool.I(:)), max(tool.I(:))];
                        tool.fName = fileName;
                    end
                case 1  %tool = VUESR(I)
                    tool.I = varargin{1}; position=[0 0 1 1];
                    heightHistogram= figure('Position', [400 200 600 600]);
                    set(heightHistogram,'Toolbar','none','Menubar','none')
                    pixelValueRange = [min(tool.I(:)), max(tool.I(:))];
                case 2  %tool = VUESR(I,position)
                    tool.I=varargin{1}; position=varargin{2};
                    heightHistogram=figure;
                    set(heightHistogram,'Toolbar','none','Menubar','none')
                    pixelValueRange=[min(tool.I(:)), max(tool.I(:))];
                case 3  %tool = VUESR(I,position,h)
                    tool.I=varargin{1};
                    position=varargin{2};
                    heightHistogram=varargin{3};
                    pixelValueRange=[min(tool.I(:)), max(tool.I(:))];
                case 4  %tool = VUESR(I,position,h,range)
                    tool.I=varargin{1};
                    position=varargin{2};
                    heightHistogram=varargin{3};
                    pixelValueRange=varargin{4};
            end
            
            if isempty(tool.I)
                tool.I=random('unif',-50,50,[100 100 3]);
            end
            %oldI = tool.I; %Save backup of original image data so that it can be recalled later.
            if isempty(position)
                position=[0 0 1 1];
            end
            
            if isempty(heightHistogram)
                heightHistogram=figure;
            end
            
            if isempty(pixelValueRange)
                pixelValueRange=[min(tool.I(:)), max(tool.I(:))];
            end
               
            %Make the aspect ratio of the figure match that of the image
            if nargin<3
                set(heightHistogram,'Units','Pixels');
                pos=get(heightHistogram,'Position');
                Af=pos(3)/pos(4);   %Aspect Ratio of the figure
                AI=size(tool.I,2)/size(tool.I,1); %Aspect Ratio of the image
                if Af>AI    %Figure is too wide, make it taller to match
                    pos(4)=pos(3)/AI;
                elseif Af<AI    %Figure is too long, make it wider to match
                    pos(3)=AI*pos(4);
                end
                set(heightHistogram,'Position',pos)
                set(heightHistogram,'Units','normalized');
            end
            tool.handles.fig = heightHistogram;
            
            initializeVUESR(tool);
            
            %%
            % Create the panels and slider
            widthSidePanel = 30; %Pixel width of the side panels
            heightHistogram = 110; %Pixel height of the histogram panel
            widthButtons = 20; %Pixel size of the buttons
            
            createSliderandPanels(tool,position,pixelValueRange,widthSidePanel,heightHistogram);            
            
            %%
            % Set up mouse button controls
            fun=@(hObject,eventdata) imageButtonDownFunction(tool,hObject,eventdata);
            set(tool.handles.I,'ButtonDownFcn',fun)
            
            % Create the tool buttons
            wp = widthSidePanel;
            widthSidePanel = widthButtons;
            buff = (wp - widthSidePanel)/2;
            
            %%
            %Create the histogram plot
            tool.handles.HistAxes = ...
                axes('Position',[0.025, 0.15, 0.95, 0.55],...
                'Parent',tool.handles.Panels.Hist);
            im = tool.I(:,:,1);
            
            tool.centers = linspace(min(double(tool.I(:))),max(double(tool.I(:))),256);
            nElements = hist(im(:),tool.centers);
            nElements = nElements./max(nElements);
            tool.handles.HistLine=plot(tool.centers,nElements,'-w','LineWidth',1);
            set(tool.handles.HistAxes,'Color','none',...
                'XColor','w','YColor','w','FontSize',9,'YTick',[])
            axis on
            hold on
            axis fill
            xlim(get(gca,'Xlim'))
            
            tool.handles.Histrange(1) = plot([pixelValueRange(1),...
                pixelValueRange(1) pixelValueRange(1)],[0, 0.5, 1],'.-r');
            tool.handles.Histrange(2) = plot([pixelValueRange(2),...
                pixelValueRange(2) pixelValueRange(2)],[0, 0.5, 1],'.-r');
            tool.handles.Histrange(3) = plot([mean(pixelValueRange),...
                mean(pixelValueRange), mean(pixelValueRange)],[0, 0.5, 1],'.--r');
            
            tool.handles.HistImageAxes = axes('Position',[0.025 0.75 0.95 0.2],...
                'Parent',tool.handles.Panels.Hist);
            set(tool.handles.HistImageAxes,'Units','Pixels');
            pos = get(tool.handles.HistImageAxes,'Position');
            set(tool.handles.HistImageAxes,'Units','Normalized');
            
            tool.handles.HistImage = ...
                imshow(repmat(tool.centers,[round(pos(4)) 1]),pixelValueRange);
            set(tool.handles.HistImageAxes,'XColor','w','YColor','w',...
                'XTick',[],'YTick',[])
            axis on;
            axis normal;
            box on;
            
            tool.centers = tool.centers;
            fun = @(hObject,evnt)histogramButtonDownFunction(tool,hObject,evnt,1);
            set(tool.handles.Histrange(1),'ButtonDownFcn',fun);
            fun = @(hObject,evnt)histogramButtonDownFunction(tool,hObject,evnt,2);
            set(tool.handles.Histrange(2),'ButtonDownFcn',fun);
            fun = @(hObject,evnt)histogramButtonDownFunction(tool,hObject,evnt,3);
            set(tool.handles.Histrange(3),'ButtonDownFcn',fun);
            
            % Create histogram checkbox
            tool.handles.Tools.Hist = uicontrol(tool.handles.Panels.Tools,...
                'Style','Checkbox','String','Hist?',...
                'Position',[buff, buff, 2.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Show Histogram',...
                'BackgroundColor','k','ForegroundColor','w');
            fun=@(hObject,evnt) ShowHistogram(tool,hObject,evnt,wp,heightHistogram);
            set(tool.handles.Tools.Hist,'Callback',fun)
            lp = buff+2.5*widthSidePanel;
            
            %%
            % Set up the resize function
            fun=@(x,y) panelResizeFunction(tool,x,y,wp,heightHistogram,widthButtons);
            set(tool.handles.Panels.Large,'ResizeFcn',fun)
            
            
            % Create window and level boxes
            tool.handles.Tools.TW = ...
                uicontrol(tool.handles.Panels.Tools,...
                'Style','text','String','W','Position',...
                [lp+buff buff widthSidePanel widthSidePanel],...
                'BackgroundColor','k','ForegroundColor','w',...
                'TooltipString','Window Width');
            tool.handles.Tools.W = ...
                uicontrol(tool.handles.Panels.Tools,...
                'Style','Edit','String',...
                num2str(pixelValueRange(2)-pixelValueRange(1)),...
                'Position',[lp+buff+widthSidePanel buff 2*widthSidePanel widthSidePanel],...
                'TooltipString','Window Width');
            tool.handles.Tools.TL = ...
                uicontrol(tool.handles.Panels.Tools,...
                'Style','text','String','L','Position',...
                [lp+2*buff+3*widthSidePanel, buff, widthSidePanel, widthSidePanel],...
                'BackgroundColor','k','ForegroundColor','w',...
                'TooltipString','Window Level');
            tool.handles.Tools.L = ...
                uicontrol(tool.handles.Panels.Tools,'Style','Edit','String',...
                num2str(mean(pixelValueRange)),'Position',...
                [lp+2*buff+4*widthSidePanel, buff, 2*widthSidePanel, widthSidePanel],...
                'TooltipString','Window Level');
            lp = lp+buff+7*widthSidePanel;
            
            % Create window and level callbacks
            fun = @(hobject,evnt) WindowLevel_callback(tool,hobject,evnt);
            set(tool.handles.Tools.W,'Callback',fun);
            set(tool.handles.Tools.L,'Callback',fun);
            
            % Create view restore button
            tool.handles.Tools.ViewRestore = ...
                uicontrol(tool.handles.Panels.Tools,'Style','pushbutton',...
                'String','','Position',[lp, buff, widthSidePanel, widthSidePanel],...
                'TooltipString','Reset Pan and Zoom');
            [iptdir, MATLABdir] = ipticondir; %#ok<ASGLU>
            icon_save = makeToolbarIconFromPNG([iptdir '/overview_zoom_in.png']);
            set(tool.handles.Tools.ViewRestore,'CData',icon_save);
            fun = @(hobject,evnt) resetViewCallback(tool,hobject,evnt);
            set(tool.handles.Tools.ViewRestore,'Callback',fun)
            lp = lp + widthSidePanel + 2*buff;
            
            % Create grid checkbox and grid lines
            axes(tool.handles.Axes)
            tool.handles.Tools.Grid = ...
                uicontrol(tool.handles.Panels.Tools,'Style','checkbox',...
                'String','Grid?','Position',...
                [lp, buff, 2.5*widthSidePanel, widthSidePanel],...
                'BackgroundColor','k','ForegroundColor','w');
            nGrid = 7;
            nMinor = 4;
            x = linspace(1,size(tool.I,2),nGrid);
            y = linspace(1,size(tool.I,1),nGrid);
            hold on;
            tool.handles.grid=[];
            gColor=[255, 38, 38]./256;
            mColor=[255, 102, 102]./256;
            for i = 1:nGrid
                tool.handles.grid(end+1) = ...
                    plot([0.5 size(tool.I,2)-0.5],[y(i) y(i)],'-',...
                    'LineWidth',1.2,'HitTest','off','Color',gColor);
                tool.handles.grid(end+1) = ...
                    plot([x(i) x(i)],[0.5 size(tool.I,1)-0.5],'-',...
                    'LineWidth',1.2,'HitTest','off','Color',gColor);
                
                if i < nGrid
                    xm = linspace(x(i),x(i+1),nMinor+2);
                    xm = xm(2:end-1);
                    ym = linspace(y(i),y(i+1),nMinor+2);
                    ym = ym(2:end-1);
                    
                    for j = 1:nMinor
                        tool.handles.grid(end+1) = ...
                            plot([.5 size(tool.I,2)-.5],[ym(j) ym(j)],'-r',...
                            'LineWidth',.9,'HitTest','off','Color',mColor);
                        tool.handles.grid(end+1) = ...
                            plot([xm(j) xm(j)],[.5 size(tool.I,1)-.5],'-r',...
                            'LineWidth',.9,'HitTest','off','Color',mColor);
                    end
                end
            end
            tool.handles.grid(end+1) = ...
                scatter(0.5+size(tool.I,2)/2,0.5+size(tool.I,1)/2,'r','filled');
            set(tool.handles.grid,'Visible','off')
            fun = @(hObject,evnt) toggleGrid(tool,hObject,evnt);
            set(tool.handles.Tools.Grid,'Callback',fun)
            set(tool.handles.Tools.Grid,'TooltipString','Toggle Gridlines')
            lp = lp+3*widthSidePanel;
            
            % Create colormap pulldown menu
            mapNames = ...
                {'Gray','Hot','Jet','HSV','Cool',...
                'Spring','Summer','Autumn','Winter',...
                'Bone','Copper','Pink','Lines',...
                'colorcube','flag','prism','white'};
            tool.handles.Tools.Color = ...
                uicontrol(tool.handles.Panels.Tools,...
                'Style','popupmenu',...
                'String',mapNames,...
                'Position',[lp, buff, 4*widthSidePanel, widthSidePanel]);
            fun = @(hObject,evnt) changeColormap(tool,hObject,evnt);
            set(tool.handles.Tools.Color,'Callback',fun)
            set(tool.handles.Tools.Color,'TooltipString','Select a colormap')
            lp = lp + 3*widthSidePanel;
            
            createButtons(tool,lp,buff,widthSidePanel);
            
        end

        function setPosition(tool,~)
            set(tool.handles.Panels.Large,'Position',newPosition)
        end
        
        function position = getPosition(tool)
            position = get(tool.handles.Panels.Large,'Position');
        end
        
        function setUnits(tool,units)
            set(tool.handles.Panels.Large,'Units',units)
        end
        
        function units = getUnits(tool)
            units = get(tool.handles.Panels.Large,'Units');
        end
        
        function setImage(varargin)
            switch nargin
                case 1
                    tool=varargin{1}; tool.I=random('unif',-50,50,[100 100 3]);
                    range=[-50 50];
                case 2
                    tool=varargin{1}; tool.I=varargin{2};
                    range=[min(tool.I(:)) max(tool.I(:))];
                case 3
                    tool=varargin{1}; tool.I=varargin{2};
                    range=varargin{3};
            end
            
            if isempty(tool.I)
                tool.I=random('unif',-50,50,[100 100 3]);
            end
            if isempty(range)
                range=[min(tool.I(:)) max(tool.I(:))];
            end
            
            
            
            % Update the histogram
            im=tool.I(:,:,1);
            tool.centers=linspace(min(tool.I(:)),max(tool.I(:)),256);
            nelements=hist(im(:),tool.centers);
            nelements=nelements./max(nelements);
            set(tool.handles.HistLine,'XData',tool.centers,'YData',nelements);
            axes(tool.handles.HistAxes);
            xlim([tool.centers(1) tool.centers(end)])
            axis fill
            
            %Update the window and level
            setWL(tool,diff(range),mean(range))
            
            %Update the image
            set(tool.handles.I,'CData',im)
            axes(tool.handles.Axes);
            xlim([0 size(tool.I,2)])
            ylim([0 size(tool.I,1)])
            
            %Update the gridlines
            axes(tool.handles.Axes);
            delete(tool.handles.grid)
            nGrid=7;
            nMinor=4;
            x=linspace(1,size(tool.I,2),nGrid);
            y=linspace(1,size(tool.I,1),nGrid);
            hold on;
            tool.handles.grid=[];
            gColor=[255 38 38]./256;
            mColor=[255 102 102]./256;
            for i=1:nGrid
                tool.handles.grid(end+1)=plot([.5 size(tool.I,2)-.5],[y(i) y(i)],'-','LineWidth',1.2,'HitTest','off','Color',gColor);
                tool.handles.grid(end+1)=plot([x(i) x(i)],[.5 size(tool.I,1)-.5],'-','LineWidth',1.2,'HitTest','off','Color',gColor);
                if i<nGrid
                    xm=linspace(x(i),x(i+1),nMinor+2); xm=xm(2:end-1);
                    ym=linspace(y(i),y(i+1),nMinor+2); ym=ym(2:end-1);
                    for j=1:nMinor
                        tool.handles.grid(end+1)=plot([.5 size(tool.I,2)-.5],[ym(j) ym(j)],'-r','LineWidth',.9,'HitTest','off','Color',mColor);
                        tool.handles.grid(end+1)=plot([xm(j) xm(j)],[.5 size(tool.I,1)-.5],'-r','LineWidth',.9,'HitTest','off','Color',mColor);
                    end
                end
            end
            tool.handles.grid(end+1)=scatter(.5+size(tool.I,2)/2,.5+size(tool.I,1)/2,'r','filled');
            toggleGrid(tool,tool.handles.Tools.Grid,[])
            
            %Update the slider
            setupSlider(tool)
            
            %Show the first slice
            showSlice(tool)
            
            
        end
        
        function I = getImage(tool)
            I=tool.I;
        end
        
        function handles=getHandles(tool)
            handles=tool.handles;
        end
        
        function setDisplayRange(tool,range)
            W=diff(range);
            L=mean(range);
            setWL(tool,W,L);
        end
        
        function range=getDisplayRange(tool)
            range=get(tool.handles.Axes,'Clim');
        end
        
        function setWindowLevel(tool,W,L)
            setWL(tool,W,L);
        end
        
        function [W,L] = getWindowLevel(tool)
            range=get(tool.handles.Axes,'Clim');
            W=diff(range);
            L=mean(range);
        end
        
        function ROI = getcurrentROI(tool)
            if ~isempty(tool.currentROI)
                if isvalid(tool.currentROI)
                    mask = createMask(tool.currentROI);
                    im=get(tool.handles.I,'CData');
                    stats= regionprops(mask,im,'Area','Perimeter','MaxIntensity','MinIntensity','MeanIntensity');
                    stats.STD=std(im(mask));
                    ROI.mask=mask;
                    ROI.stats=stats;
                end
            else
                ROI=[];
            end
        end
        
        function setCurrentSlice(tool,slice)
            showSlice(tool,slice)
        end
        
        function slice = getCurrentSlice(tool)
            slice=round(get(tool.handles.Slider,'value'));
        end
        
    end
    
    methods (Access = private)
        
        function addhandlesROI(tool,h)
            tool.handlesROI{end+1}=h;
        end
        
        function scrollWheel(tool,~,evnt)
            %Check to see if the mouse is hovering over the axis
            units=get(tool.handles.fig,'Units');
            set(tool.handles.fig,'Units','Pixels')
            point=get(tool.handles.fig, 'CurrentPoint');
            set(tool.handles.fig,'Units',units)
            
            units=get(tool.handles.Panels.Large,'Units');
            set(tool.handles.Panels.Large,'Units','Pixels')
            pos_p=get(tool.handles.Panels.Large,'Position');
            set(tool.handles.Panels.Large,'Units',units)
            
            units=get(tool.handles.Panels.Image,'Units');
            set(tool.handles.Panels.Image,'Units','Pixels')
            pos_a=get(tool.handles.Panels.Image,'Position');
            set(tool.handles.Panels.Image,'Units',units)
            
            xmin=pos_p(1)+pos_a(1); xmax=xmin+pos_a(3);
            ymin=pos_p(2)+pos_a(2); ymax=ymin+pos_a(4);
            
            if point(1)>=xmin && point(1)<=xmax && point(2)>=ymin && point(2)<=ymax
                newSlice=get(tool.handles.Slider,'value')-evnt.VerticalScrollCount;
                if newSlice>=1 && newSlice <=size(tool.I,3)
                    set(tool.handles.Slider,'value',newSlice);
                    showSlice(tool)
                end
            end
            
        end
        
        function showSlice(varargin)
            switch nargin
                case 1
                    tool=varargin{1};
                    n=round(get(tool.handles.Slider,'value'));
                case 2
                    tool=varargin{1};
                    n=varargin{2};
                otherwise
                    tool=varargin{1};
                    n=round(get(tool.handles.Slider,'value'));
            end
            
            if n < 1
                n=1;
            end
            
            if n > size(tool.I,3)
                n=size(tool.I,3);
            end
            
            set(tool.handles.I,'CData',tool.I(:,:,n))
            set(tool.handles.SliceText,'String',[num2str(n) '/' num2str(size(tool.I,3))])
            if get(tool.handles.Tools.Hist,'value')
                im=tool.I(:,:,n);
                nelements=hist(im(:),tool.centers); nelements=nelements./max(nelements);
                set(tool.handles.HistLine,'YData',nelements);
            end
            
        end
        
        function setupSlider(tool)
            n=size(tool.I,3);
            if n==1
                set(tool.handles.Slider,'visible','off');
            else
                set(tool.handles.Slider,'visible','on');
                set(tool.handles.Slider,'min',1,'max',size(tool.I,3),'value',1)
                set(tool.handles.Slider,'SliderStep',[1/(size(tool.I,3)-1) 1/(size(tool.I,3)-1)])
                fun=@(hobject,eventdata)showSlice(tool,[],hobject,eventdata);
                set(tool.handles.Slider,'Callback',fun);
            end
            
        end
        
        function setWL(tool,W,L)
            set(tool.handles.Axes,'Clim',[L-W/2 L+W/2])
            set(tool.handles.Tools.W,'String',num2str(W));
            set(tool.handles.Tools.L,'String',num2str(L));
            set(tool.handles.HistImageAxes,'Clim',[L-W/2 L+W/2])
            set(tool.handles.Histrange(1),'XData',[L-W/2 L-W/2 L-W/2])
            set(tool.handles.Histrange(2),'XData',[L+W/2 L+W/2 L+W/2])
            set(tool.handles.Histrange(3),'XData',[L L L])
        end
        
        function WindowLevel_callback(tool,~,~)
            range=get(tool.handles.Axes,'Clim');
            Wold=range(2)-range(1); Lold=mean(range);
            W=str2num(get(tool.handles.Tools.W,'String')); %#ok<ST2NM>
            if isempty(W) || W<=0
                W=Wold;
                set(tool.handles.Tools.W,'String',num2str(W))
            end
            L=str2num(get(tool.handles.Tools.L,'String')); %#ok<ST2NM>
            if isempty(L)
                L=Lold;
                set(tool.handles.Tools.L,'String',num2str(L))
            end
            setWL(tool,W,L)
        end
        
        function imageButtonDownFunction(tool,~,~)
            cp=get(tool.handles.Axes,'CurrentPoint');
            cp=[cp(1,1) cp(1,2)];
            switch get(tool.handles.fig,'SelectionType')
                case 'normal'   %Adjust window and level
                    CLIM=get(tool.handles.Axes,'Clim');
                    W=CLIM(2)-CLIM(1);
                    L=mean(CLIM);
                    fun=@(src,evnt) adjustContrastMouse(tool,src,evnt,cp,tool.handles.Axes,W,L);
                    fun2=@(src,evnt) buttonUpFunction(tool,src,evnt);
                    set(tool.handles.fig,'WindowButtonMotionFcn',fun,'WindowButtonUpFcn',fun2)
                case 'extend'  %Zoom
                    fun=@(src,evnt) adjustZoomMouse(tool,src,evnt,cp,tool.handles.Axes);
                    fun2=@(src,evnt) buttonUpFunction(tool,src,evnt);
                    set(tool.handles.fig,'WindowButtonMotionFcn',fun,'WindowButtonUpFcn',fun2)
                case 'alt'
                    xlims=get(tool.handles.Axes,'Xlim');
                    ylims=get(tool.handles.Axes,'Ylim');
                    fun=@(src,evnt) adjustPanMouse(tool,src,evnt,cp,tool.handles.Axes,xlims,ylims);
                    fun2=@(src,evnt) buttonUpFunction(tool,src,evnt);
                    set(tool.handles.fig,'WindowButtonMotionFcn',fun,'WindowButtonUpFcn',fun2)
            end
        end
        
        function histogramButtonDownFunction(tool,~,~,line)
            
            switch line
                case 1 %Lower limit of range
                    fun=@(src,evnt) newLowerRangePosition(tool,src,evnt,tool.handles.HistAxes);
                    fun2=@(src,evnt) buttonUpFunction(tool,src,evnt);
                    set(tool.handles.fig,'WindowButtonMotionFcn',fun,'WindowButtonUpFcn',fun2)
                case 2 %Upper limt of range
                    fun=@(src,evnt) newUpperRangePosition(tool,src,evnt,tool.handles.HistAxes);
                    fun2=@(src,evnt) buttonUpFunction(tool,src,evnt);
                    set(tool.handles.fig,'WindowButtonMotionFcn',fun,'WindowButtonUpFcn',fun2)
                case 3 %Middle line
                    fun=@(src,evnt) newLevelRangePosition(tool,src,evnt,tool.handles.HistAxes);
                    fun2=@(src,evnt) buttonUpFunction(tool,src,evnt);
                    set(tool.handles.fig,'WindowButtonMotionFcn',fun,'WindowButtonUpFcn',fun2)
            end
        end
        
        function toggleGrid(tool, hObject,~)
            if get(hObject,'Value')
                set(tool.handles.grid,'Visible','on')
            else
                set(tool.handles.grid,'Visible','off')
            end
        end
        
        function changeColormap(tool,hObject,~)
            n=get(hObject,'Value');
            maps=get(hObject,'String');
            colormap(tool.handles.Axes, maps{n});
            colormap(tool.handles.HistAxes, maps{n});
            colormap(tool.handles.HistImageAxes, maps{n});
        end
        
        function exportROI(tool,~,~)
            if ~isempty(tool.currentROI)
                if isvalid(tool.currentROI)
                    mask = createMask(tool.currentROI);
                    im=get(tool.handles.I,'CData');
                    stats= regionprops(mask,im,'Area','Perimeter','MaxIntensity','MinIntensity','MeanIntensity');
                    stats.STD=std(im(mask));
                    ROI.mask=mask;
                    ROI.stats=stats;
                    name = inputdlg('Enter variable name');
                    name=name{1};
                    assignin('base', name, ROI)
                end
            end
        end
        
        function measureImageCallback(tool,~,~,type)
            
            switch type
                case 'ellipse'
                    fcn = makeConstrainToRectFcn('imellipse',[1 size(tool.I,2)],[1 size(tool.I,1)]);
                    h = imellipse(tool.handles.Axes,'PositionConstraintFcn',fcn);
                    addhandlesROI(tool,h)
                    fcn=@(pos) newROIposition(tool,pos,h);
                    addNewPositionCallback(h,fcn);
                    setPosition(h,getPosition(h));
                    tool.pointLog = [];
                            
                case 'rectangle'
                    fcn = makeConstrainToRectFcn('imrect',[1 size(tool.I,2)],[1 size(tool.I,1)]);
                    h = imrect(tool.handles.Axes,'PositionConstraintFcn',fcn);
                    addhandlesROI(tool,h)
                    fcn=@(pos) newROIposition(tool,pos,h);
                    addNewPositionCallback(h,fcn);
                    setPosition(h,getPosition(h));
                    tool.pointLog = [];
                    
                case 'polygon'
                    fcn = makeConstrainToRectFcn('impoly',[1 size(tool.I,2)],[1 size(tool.I,1)]);
                    h = impoly(tool.handles.Axes,'PositionConstraintFcn',fcn);
                    addhandlesROI(tool,h)
                    fcn=@(pos) newROIposition(tool,pos,h);
                    addNewPositionCallback(h,fcn);
                    setPosition(h,getPosition(h));
                    tool.pointLog = [];
                     
                case 'ruler'
                    h = imdistline(tool.handles.Axes);
                    fcn = makeConstrainToRectFcn('imline',[1 size(tool.I,2)],[1 size(tool.I,1)]);
                    setPositionConstraintFcn(h,fcn);
           
                case 'profile'
                    axes(tool.handles.Axes);
                    improfile(); grid on;
                otherwise
            end
            
            
        end
        
        function deletecurrentROI(tool,~,~)
            %if ~ISEMPTY(tool.currentROI)
                if isvalid(tool.currentROI)
                    delete(tool.currentROI)
                    set(tool.handles.ROIinfo,'String','STD:                    Mean:                    ');
                    tool.pointLog = [];
                end
            %end
        end
        
        function displayHelp(~,~,~)
            
            message={'Welcome to VUESR', ...
                '',...
                'Left Mouse Button: Window and Level', ...
                'Right Mouse Button: Pan', ...
                'Middle Mouse Button: Zoom', ...
                'Scroll Wheel: Change Slice',...
                '',...
                };
            
            msgbox(message)
        end
        
        function CropImageCallback(tool,~,~)
            [~, rect] = imcrop(tool.handles.Axes);
            rect=round(rect);
            setImage(tool, tool.I(rect(2):rect(2)+rect(4)-1,rect(1):rect(1)+rect(3)-1,:))
            tool.pointLog = [];
        end
        
        function AutoCropCallback(tool,~,~)
           %Crop the image to only show ultrasound measurements
           x = dicominfo(tool.fName);
           if  isfield(x, 'SequenceOfUltrasoundRegions') && ...             %WDR fixed
                isfield(x.SequenceOfUltrasoundRegions, 'Item_1') && ...
                isfield(x.SequenceOfUltrasoundRegions.Item_1, 'RegionLocationMinX0')
                xmin = x.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0;
                xmax = x.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1;
                ymin = x.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0;
                ymax = x.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1;

                setImage(tool, tool.I(ymin:ymax,xmin:xmax,:));
           else
               warndlg('Autocrop only works for Ultrasound');
           end
        end
        
        function CropFramesCallback(tool,~, ~)
            g = figure('Visible','off','Position',[200,200,150,120],'NumberTitle', 'off','ToolBar','none','MenuBar','none','Name','Edit/Export Data');
                                  
              frameRange = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
                           'String', 'New frame range:', 'Position', [25 80 120 20], 'BackgroundColor', [.8,.8,.8]);
              currentRange = uicontrol('Style', 'text','FontSize',6,'HorizontalAlignment','left',...
                           'String', ['(Current Frame Range: ' '1 - ' num2str(size(tool.I,3)) ')'] , 'Position', [12 30 150 20], 'BackgroundColor', [.8,.8,.8]);
              minFrame = uicontrol('Style', 'edit','FontSize', 10,...
                         'String', num2str(tool.accFrames(1)), 'Position', [30 60 30 20]);
               maxFrame = uicontrol('Style', 'edit','FontSize', 10,...
                         'String', num2str(tool.accFrames(2)), 'Position', [ 80 60 30 20]);
               to = uicontrol('Style', 'text','FontSize', 9,'HorizontalAlignment','left',...
                           'String', 'to', 'Position', [64 60 12 20], 'BackgroundColor', [.8,.8,.8]);
                Crop = uicontrol('Style', 'pushbutton', 'String', 'Crop',...
                    'FontSize', 10,'Position', [10 10 60 20],'Callback', @(hObject,evnt) Crop_callback(tool,hObject,evnt));
                Cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel',...
                    'FontSize', 10,'Position', [80 10 60 20], 'Callback', @(hObject,evnt) Cancel_callback(tool,hObject,evnt));

                set(g,'Name','Crop Frames')
                set([g,frameRange, currentRange, minFrame, to, maxFrame, Crop, Cancel], 'Units','normalized') 
                % Move the GUI to the center of the screen. 
                movegui(g,'center')
                % Make the GUI visible. 
                set(g,'Visible','on','Resize','off');

                function Cancel_callback(~,~,~)
                     %Closes window without saving changes
                      close(g);
                end
                function Crop_callback(tool,~,~)
%%%%%%%%%%%%%%%%%%%%%%%%% 
filename = 'SSP.xlsx';
BBB=[ str2double(get(minFrame,'String')) str2double(get(maxFrame,'String')) ]; 
sheet = 2;
xlRange = 'A1';
xlswrite(filename,BBB,sheet,xlRange);
%%%%%%%%%%%%%%%%%%%%%%%%%
                    setImage(tool, tool.I(:,:,round(str2double(get(minFrame,'String'))):round(str2double(get(maxFrame,'String')))))
                    tool.pointLog = [];
                    tool.accFrames = [1 size(tool.I,3)-1];
                    close(g);
                    msgbox('New frame range set')
                end
        end
        
        function openImageCallback(tool, ~,~)
                [fileName, filePath] = uigetfile('*.DCM;*.dcm;*.mat;*', ...
                    'Choose DICOM images to import', pwd, ...
                    'MultiSelect', 'off');
                if filePath(1) == 0
                    disp('No files chosen, exiting function')
                    return;
                else
                    fileName = fullfile(filePath, fileName);               %WDR fixed
                    disp(['User selected: ', fullfile(fileName)]);
                    [~, ~, ext] = fileparts(fileName);
                    if strcmp(ext,'.DCM') || strcmp(ext,'.dcm') || strcmp(ext,'')
                        tool.I = permute(dicomread(fileName),[1, 2, 4, 3]);
                    else
                        load( fileName);
                        tool.I = permute(image_change,[1 2 4 3]);
                    end       
                    tool.fName = fileName;
                end
                initializeVUESR(tool);
                setImage(tool,tool.I);
                initializeVUESR(tool);
                
        end
        
        function resetViewCallback(tool,~,~)
            set(tool.handles.Axes,'Xlim',get(tool.handles.I,'XData'))
            set(tool.handles.Axes,'Ylim',get(tool.handles.I,'YData'))
        end
        
        function createSliderandPanels(tool,position,pixelValueRange,widthSidePanel,heightHistogram)
            tool.handles.Panels.Large = ...
                uipanel(tool.handles.fig,'Position',position,'Title','');
            set(tool.handles.Panels.Large,'Units','Pixels');
            pos = get(tool.handles.Panels.Large,'Position');
            set(tool.handles.Panels.Large,'Units','normalized');
            
            tool.handles.Panels.Hist = ...
                uipanel(tool.handles.Panels.Large,'Units','Pixels', ...
                'Position',[widthSidePanel, pos(4)-widthSidePanel-heightHistogram, ...
                pos(3)-2*widthSidePanel, heightHistogram],'Title','');
            
            tool.handles.Panels.Image = ...
                uipanel(tool.handles.Panels.Large,'Units','Pixels',...
                'Position',[widthSidePanel, widthSidePanel, ...
                pos(3)-3.8*widthSidePanel, pos(4)-2*widthSidePanel],'Title','');
            
            tool.handles.Panels.Tools = ....
                uipanel(tool.handles.Panels.Large,'Units','Pixels', ....
                'Position',[0, pos(4)-widthSidePanel, ...
                pos(3), widthSidePanel],'Title','');
            
            tool.handles.Panels.ROItools = ...
                uipanel(tool.handles.Panels.Large,'Units','Pixels',...
                'Position',[pos(3)-2.8*widthSidePanel,  widthSidePanel, ...
               2.8*widthSidePanel, pos(4)-2*widthSidePanel],'Title','');
            
            tool.handles.Panels.Slider = ...
                uipanel(tool.handles.Panels.Large,'Units','Pixels',...
                'Position',[0, widthSidePanel, ...
                widthSidePanel, pos(4)-2*widthSidePanel],'Title','');
            
            tool.handles.Panels.Info = ...
                uipanel(tool.handles.Panels.Large,'Units','Pixels',...
                'Position',[0, 0, pos(3), widthSidePanel],'Title','');
            
            pan_cell = struct2cell(tool.handles.Panels);                    %WDR fix
            set( vertcat(pan_cell{:}), ...
                'BackgroundColor','k','ForegroundColor','w','HighlightColor','k')
            
            %%
            % Create slider to scroll through image stack
            tool.handles.Slider = ...
                uicontrol(tool.handles.Panels.Slider,'Style','Slider',...
                'Units','Normalized','Position',[0, 0, 1, 1],...
                'TooltipString','Change Slice (can use scroll wheel also)');
            setupSlider(tool)
            fun = @(scr,evnt) scrollWheel(tool,scr,evnt);
            set(tool.handles.fig,'WindowScrollWheelFcn',fun);
            
            % Create image axis
            tool.handles.Axes = ...
                axes('Position', [0, 0, 1, 1], ...
                'Parent',tool.handles.Panels.Image,'Color','none');
            tool.handles.I = imshow(tool.I(:,:,1),pixelValueRange);
            set(tool.handles.Axes,'Position',[0, 0, 1, 1],...
                'Color','none','XColor','r','YColor','r',...
                'GridLineStyle','--','LineWidth',1.5,...
                'XTickLabel','','YTickLabel','');
            axis off
            grid off
            axis fill
            
            % Set up image info display
            tool.handles.Info = uicontrol(tool.handles.Panels.Info,...
                'Style','text','String','(x,y) val','Units','Normalized',...
                'Position',[0, 0.1, 0.5, 0.8],'BackgroundColor','k',...
                'ForegroundColor','w','FontSize',12,...
                'HorizontalAlignment','Left');
            tool.handles.ROIinfo = uicontrol(tool.handles.Panels.Info,...
                'Style','text','String',...
                'STD:                    Mean:                    ',...
                'Units','Normalized','Position',[0.5, 0.1, 0.5, 0.8],...
                'BackgroundColor','k','ForegroundColor','w',....
                'FontSize',12,'HorizontalAlignment','Right');
            fun = @(src,evnt)getImageInfo(tool,src,evnt);
            set(tool.handles.fig,'WindowButtonMotionFcn',fun);
            tool.handles.SliceText=uicontrol(tool.handles.Panels.Tools,...
                'Style','text','String',['1/' num2str(size(tool.I,3))], ...
                'Units','Normalized','Position',[0.5, 0.1, 0.48, 0.8], ...
                'BackgroundColor','k','ForegroundColor','w', ...
                'FontSize',12,'HorizontalAlignment','Right');
        end
        
        % ***********************NEW***********************************************************
        function pixelCalibrateCallback(tool,~,~)
            %Select Points to Track
            uiwait(msgbox('Select two points with a know distance, then hit "Enter"'));
   
            figHandle = gcf;
            [poiX, poiY] = getpts(figHandle);

            poiX = round(poiX);     poiY = round(poiY);
            %Calculate the distance in pixels
            pixels = sqrt ((poiX(1) - poiX(2))^2+(poiY(1) - poiY(2))^2);
            
            %Have user input distance in cm
            prompt = {'Input distance between 2 points (cm):'};
                    dlg_title = 'Input';
                    num_lines = 1;
                    default = {'1'};
                    options.Resize='on';
                    options.WindowStyle='normal';
                    answer = inputdlg(prompt,dlg_title,num_lines,default,options);
                    cm = str2double(answer{1,1});
            
           tool.calibration = cm/pixels;
           disp(tool.calibration)
           msgbox('Calibration Complete')
                
        end
        
        function pixelSettingsCallback(tool,~,~)
            
            prompt = {'cm/pixel calibration factor:', '% of pixels analyzed (1-100%):','Frame range for accumulated strain (separated by a comma):','Mean Arterial Blood Pressure (mmHg)','Frame Rate (Hz):'};
                    dlg_title = 'Settings';
                    num_lines = [1, 60; 1, 60; 1, 60; 1, 60; 1,60];
                    cal = num2str(tool.calibration);
                    pixels = num2str(tool.pixelDensity);
                    oldpixels = tool.pixelDensity;
                    accRange = strcat([num2str((tool.accFrames(1))),',',num2str((tool.accFrames(2)))]);
                    oldRange = tool.accFrames;
                    bloodp = num2str(tool.map);
                    frameRate = num2str(tool.fRate);
                    %accRange = '1,10';
                    %disp(accRange); disp(size(accRange)); disp(cal); disp(size(cal)); disp(pixels);
                    default = {cal,pixels,accRange,bloodp,frameRate};
                    options.Resize='on';
                    options.WindowStyle='normal';
                    answer = inputdlg(prompt,dlg_title,num_lines,default,options);
                    if ~isempty(answer)
                        tool.calibration = str2double(answer{1,1});
                        tool.pixelDensity = str2double(answer{2,1});
                        tool.accFrames = str2num(answer{3,1}); %#ok<ST2NM>
                        tool.map = str2num(answer{4,1}); %#ok<ST2NM>
                        tool.fRate = str2num(answer{5,1}); %#ok<ST2NM>
                    end
                    if oldpixels ~= tool.pixelDensity || ~isequal(oldRange, tool.accFrames)
                        tool.pointLog = [];
                    end
        end
        
        function pixelTrack2Callback(tool,~,~)
              
                    dicomFrames = size(tool.I,3);
                    newI = uint8(tool.I); 
                    J = uint8(tool.I);
                    % Get region of interest
                    framenum = 1;
                    objectFrame = J(:,:,framenum);
                   
                    %Select Points to Track
                    uiwait(msgbox('Select two points to track, then hit "Enter"'));
                    
                    figHandle = gcf;
                    [poiX, poiY] = getpts(figHandle);

                    poiX = round(poiX);     poiY = round(poiY);
%%%%%%%%%%%%%%%%%%%%%%%%% RECORD IN EXCEL FILE THE VALUES OF THE POINTS TO
%%%%%%%%%%%%%%%%%%%%%%%%% BE TRACKED....
filename = 'SSP.xlsx';
AAA=poiX; AAAA=poiY;
sheet = 1;
xlRange = 'A1';
xlswrite(filename,AAA,sheet,xlRange);
xlRange = 'B1';
xlswrite(filename,AAAA,sheet,xlRange);
%%%%%%%%%%%%%%%%%%%%%%%%%
                    nPoints = size(poiX,1);
                   
                    if  nPoints  ~= 2
                        warndlg('Please select only 2 points. Exiting function', 'Warning');
                        return;
                    end
                    tool.pointLog2 = zeros(nPoints, 2, dicomFrames);
                    points = [poiX, poiY];
                    pointImage = insertMarker(objectFrame, points, '+', 'Color', 'white');
                    newI(:,:,1) = pointImage(:,:,1);
                    pointDist = zeros(3,dicomFrames);
                    
                    % Create object tracker
                    tracker = vision.PointTracker('MaxBidirectionalError', 3);

                    % Initialize object tracker
                    initialize(tracker, points(:,:,1), objectFrame);
                    % Show the points getting tracked
                    while framenum <= dicomFrames
                         %Track the points     
                          frame =J(:,:,framenum);
                          [points, validity,scores] = step(tracker, frame);
                          tool.pointLog2(:,:,framenum) = points;
                          out = insertMarker(frame, points(validity, :), '+', 'Color', 'white');
                          newI(:,:,framenum) = out(:,:,1);
                          scoreshold(:,framenum)= scores;
                          if any(validity==0) || any(scores<0.5)
                              disp('The current points being tracked have been lost')
                              return
                          end
                          %Compute the distance between the 2 points and
                          %store in a matrix 1: total distance 2: x distance
                          %3: y distance
                          pointDist(1,framenum) = sqrt((tool.pointLog2(1,1,framenum) - tool.pointLog2(2,1,framenum)).^2+(tool.pointLog2(1,2,framenum) - tool.pointLog2(2,2,framenum)).^2);
                          pointDist(2,framenum) = (tool.pointLog2(2,1,framenum) - tool.pointLog2(1,1,framenum));
                          pointDist(3,framenum) = (tool.pointLog2(1,2,framenum) - tool.pointLog2(2,2,framenum));
                          framenum = framenum + 1;
                    end
                  %  disp(mean(mean(scoreshold)))
                    %Convert pixels to cm and percent
                    cm = tool.calibration;
                    tool.pointDistCm = pointDist.*cm;
                    pointDistPercent = pointDist.*100./max(max(pointDist));
                    
                    %Distensibility = dVolume / (Original V * dPressure)
                    %Compliance = dVolume / dPressure
                    pressure = tool.map*133.32; %mmHg to Pa
                    %volume = pi.*(tool.pointDistCm./200).^2; %Volume of circle in m2
                    %deltaV = max(volume) - min(volume);
                    %compliance = deltaV/(min(volume)*pressure);
                    
                    %hoopStress = pressure*mean(tool.pointDistCm)/(2*thickness);
                    %radialStress = -pressure/2;
                    
                    pointDistTp = tool.pointDistCm - min(tool.pointDistCm);
                    pointDistTpPercent = pointDistPercent-min(pointDistPercent);
                    
                    tool.zdiameter = mean(tool.pointDistCm);
                    tool.zdistensibility = mean(pointDistTp);
                    %tool.zcompliance = mean(compliance);
                    
                    sampleFreq = tool.fRate;
                    time = (1:dicomFrames)/sampleFreq;
                    
                    %Find heart rate by looking at mins in cardiac cycle
 %                   [~,indd]=lmin(tool.pointDistCm,1);
  %                  indd1 = indd/sampleFreq;
                    
   %                 avgbeat = mean(diff(indd1));
    %                tool.hRate = 60/avgbeat; %Heart rate in bpm
                    
%                     ImageViewer(newI); "Imageviewer function isn't often
%                     working in many computers. Instead of using this
%                     function parameters are saved for post processing."
                        save newI.mat newI;
                    %Average cardio cycle
%                    mymatrix = [];
 %                   mymatrix2 = [];
  %                  for p = 1:length(indd)-1
                        %not done
                        %use a matrix here taking the mode of the lengths
                        %between the mins and use that mode length to cut
                        %the other cardio cycles to get an average of a set
                        %length becaus otherwise the lengths will vary and
                        %cause dimensional mismatches and an unclear frame.
   %                     modefinder(p,1) = indd(p+1) - indd(p);
    %                end
     %               [mymode,freq] = mode(modefinder);
      %              if freq > 1
       %                 for p = 1:length(indd)-1
        %                   mymatrix(p,1:mymode) = tool.pointDistCm(indd(p):indd(p+1));
         %               end
          %                 avgcycle = mean(mymatrix);
           %         elseif freq == 1
            %            mymooder something;
             %       end
              %      figure
               %     plot(time(1:length(indd(p+1):indd(p))),avgcycle)
                    % Plot the tracked pixel movement in a switchable GUI
                    % Create and then hide the GUI as it is being constructed. 
                    f = figure('Visible','off','Position',[360,500,525,350]); %Left bottom width height

                    % Construct the components. 

                    hpopup = uicontrol('Style','popupmenu',... 
                        'String',{'Distance (cm)','Distance (%)','Distensibility (cm)','Distensibility (%)'},... 
                        'Position',[25,320,100,25],...
                        'Callback',{@(hObject,evnt) popup_menu_Callback(tool,hObject,evnt)});

                    hdrawrange = uicontrol('Style','pushbutton',... 
                        'String','Select Strain Range','Position',[150,320,125,25],...
                        'Callback',{@(hObject,evnt) drawrange_button_Callback(tool,hObject,evnt)});
                    
                    hcropframes =uicontrol('Style','pushbutton',... 
                        'String','Crop # frames','Position',[300,320,100,25],...
                        'Callback',{@(hObject,evnt) cropframes_button_Callback(tool,hObject,evnt)});
                    
                    
                    ha = axes('Units','pixels','Position',[25,25,425,270]); 

                    % Change units to normalized so components resize automatically. 
                    set([f,ha,hpopup,hdrawrange, hcropframes],'Units','normalized');

                    % Create a plot in the axes. 
%%%%%%%%%%%ADDED EL BARTO 02112019%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%ADDED EL BARTO 02112019%%%%%%%%%%%%%%%%%%%%%%%%%%
                    DDD=[time', tool.pointDistCm(1,:)'];
%%%%%%%%%%%ADDED EL BARTO 02112019%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%ADDED EL BARTO 02112019%%%%%%%%%%%%%%%%%%%%%%%%%%
                    plot(time, tool.pointDistCm(1,:))
                    xlabel('Time (s)'); ylabel('Distance (cm)')
                    title('Distance between 2 points')

                    % Assign the GUI a name to appear in the window title. 
                    set(f,'Name','Distance Between 2 Points')
                    % Move the GUI to the center of the screen. 
                    movegui(f,'center')
                     
                    %Comment out the above plot command and these can be
                    %used to make 4 subplots of the x and y distances
                    %between two points and the associated x and ystrains 
                    %for the material between them
                    
%                    subplot(2,2,1)
 %                   plot(time, tool.pointDistCm(2,:)),
  %                  title('X Distance')
   %                xlabel('Time (s)'); ylabel('Distance (cm)');
    %                subplot(2,2,2)
     %               for j = 2:length(tool.pointDistCm(2,:))
      %                  xstrain2pts(j-1) = (tool.pointDistCm(2,j) - tool.pointDistCm(2,1))./tool.pointDistCm(2,1);
       %             end
        %            plot(time(2:length(xstrain2pts)+1), xstrain2pts)
         %           title('X Strain'); ylabel('Strain'); xlabel('Time (s)');
          %          subplot(2,2,3)
           %         plot(time, tool.pointDistCm(3,:))
            %        title('Y Distance'); xlabel('Time (s)'); ylabel('Distance (cm)');
             %       subplot(2,2,4)
              %      for j = 2:length(tool.pointDistCm(3,:))
               %         ystrain2pts(j-1) = (tool.pointDistCm(3,j) - tool.pointDistCm(3,1))./tool.pointDistCm(3,1);
                %    end
                 %   plot(time(2:length(ystrain2pts)+1), ystrain2pts)
                  %  title('Y Strain'); ylabel('Strain'); xlabel('Time (s)');
                  % Make the GUI visible.
                    set(f,'Visible','on');
                    
                    % Pop-up menu callback. Read the pop-up menu and display property
%%%%%%%%%%%%%%% MODIFICATION BARTO 02112019%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODIFICATION BARTO 02112019%%%%%%%%%%%%%%%%%%%%%  
% FLAG ed=1; %El barto esto no hace nada...
uiwait(msgbox('Save the picture in Tiff format first, please hit Enter to continue'));                    
uiwait(msgbox('Select the first 5 Max and after the 5 min,please hit Enter to continue'));                    
figHandle = gcf;
[posicionX, posicionY] = getpts(figHandle);
filename = 'SSP.xlsx';
CCC=[posicionX, posicionY];
sheet = 3;
xlRange = 'A1';
xlswrite(filename,CCC,sheet,xlRange);
%DDD=[time', tool.pointDistCm(1,:)'];
 sheet = 4;
 xlRange = 'A1';
 xlswrite(filename,DDD,sheet,xlRange);
uiwait(msgbox('The data has been saved in the file: SSP.xlsx'));                    
%%%%%%%%%%%%%%% MODIFICATION BARTO 02112019%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%% MODIFICATION BARTO 02112019%%%%%%%%%%%%%%%%%%%%%  
                        function popup_menu_Callback(tool,source,~) 
                            % Determine the selected data set. 
                            str = get(source, 'String'); 
                            val = get(source,'Value'); 

                            % Set current data to the selected data set.
                            switch str{val} 
                                case 'Distance [cm]' % User selects cm
                                    plot(time, tool.pointDistCm)
                                    xlabel('Time [s]'); ylabel('Distance [cm]')
                                    title('Distance between 2 points')
                                case 'Distensibility [cm]' % User selects cm
                                    plot(time, pointDistTp)
                                    xlabel('Time [s]'); ylabel('Distance [cm]')
                                    title('Distensibility between 2 points')    
                                case 'Distance [%]'   % User selects percent
                                    plot(time, pointDistPercent)
                                    xlabel('Time [s]'); ylabel('Distance [%]')
                                    title('Distance between 2 points')
                                case 'Distensibility [%]' % User selects cm
                                    plot(time, pointDistTpPercent)
                                    xlabel('Time [s]'); ylabel('Distance [%]')
                                    title('Distensibility between 2 points')   
                            end
                        end

                         function drawrange_button_Callback(tool,~,~)
                                uiwait(msgbox({'Draw line for the lower limit of the accumulated strain range' 'Double click to save positions'}));
                                hmin = imline;
                                pos1 = wait(hmin);
                                uiwait(msgbox({'Draw line for the upper limit of the accumulated strain range.'  'Double click to save positions'}));
                                hmax = imline;
                                pos2 = wait(hmax);
                                tool.accFrames = [round((pos1(1,1)+pos1(2,1))*sampleFreq/2) , round((pos2(1,1)+pos2(2,1))*sampleFreq/2)];
                                msgbox(['Frame range set as: ' num2str(tool.accFrames(1)) '-' num2str(tool.accFrames(2))]);    
                         end
                         
                         function cropframes_button_Callback(tool,~,~)
                                  g = figure('Visible','off','Position',[200,200,150,120],'NumberTitle', 'off','ToolBar','none','MenuBar','none','Name','Edit/Export Data');
                                  
                                  frameRange = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
                                               'String', 'New frame range:', 'Position', [25 80 120 20], 'BackgroundColor', [.8,.8,.8]);
                                  currentRange = uicontrol('Style', 'text','FontSize',6,'HorizontalAlignment','left',...
                                               'String', ['(Current Frame Range: ' '1 - ' num2str(size(tool.I,3)) ')'] , 'Position', [12 30 150 20], 'BackgroundColor', [.8,.8,.8]);
                                  minFrame = uicontrol('Style', 'edit','FontSize', 10,...
                                             'String', num2str(tool.accFrames(1)), 'Position', [30 60 30 20]);
                                   maxFrame = uicontrol('Style', 'edit','FontSize', 10,...
                                             'String', num2str(tool.accFrames(2)), 'Position', [ 80 60 30 20]);
                                   to = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
                                               'String', 'to', 'Position', [64 60 12 20], 'BackgroundColor', [.8,.8,.8]);
                                    Crop = uicontrol('Style', 'pushbutton', 'String', 'Crop',...
                                        'FontSize', 10,'Position', [10 10 60 20],'Callback', @(hObject,evnt) Crop_callback(tool,hObject,evnt));
                                    Cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel',...
                                        'FontSize', 10,'Position', [80 10 60 20], 'Callback', @(hObject,evnt) Cancel_callback(tool,hObject,evnt));
                         
                                    set(g,'Name','Crop Frames')
                                    set([g,frameRange, currentRange, minFrame, to, maxFrame, Crop, Cancel], 'Units','normalized') 
                                    % Move the GUI to the center of the screen. 
                                    movegui(g,'center')
                                    % Make the GUI visible. 
                                    set(g,'Visible','on','Resize','off');
                                    
                                    function Cancel_callback(~,~,~)
                                         %Closes window without saving changes
                                          close(g);
                                    end
                                    function Crop_callback(tool,~,~)
                                        setImage(tool, tool.I(:,:,round(str2double(get(minFrame,'String'))):round(str2double(get(maxFrame,'String')))))
                                        tool.pointLog = [];
                                        tool.accFrames = [1 size(tool.I,3)-1];
                                        close(g);
                                        msgbox('New frame range set')
                                    end
                                    
                         end
                         
                        function [lmval,indd]=lmin(xx,filt)
                            %Find the local minima in a data set, excludes the first
                            %and last points.
                            % Created by Serge Koptenko, Guigne International Ltd.
                            x=xx;
                            len_x = length(x);
                                fltr=[1 1 1]/3;
                              if nargin <2, filt=0; 
                                else
                            x1=x(1); x2=x(len_x); 

                                for jj=1:filt
                                c=conv(fltr,x);
                                x=c(2:len_x+1);
                                x(1)=x1;  
                                    x(len_x)=x2; 
                                end
                              end

                            lmval=[];
                            indd=[];
                            i=2;		% start at second data point in time series

                                while i < len_x-1
                                if x(i) < x(i-1)
                                   if x(i) < x(i+1)	% definite min
                            lmval =[lmval x(i)];
                            indd = [ indd i];

                                   elseif x(i)==x(i+1)&&x(i)==x(i+2)	% 'long' flat spot
                            %lmval =[lmval x(i)];	%1   comment these two lines for strict case 
                            %indd = [ indd i];	%2 when only  definite min included
                            i = i + 2;  		% skip 2 points

                                   elseif x(i)==x(i+1)	% 'short' flat spot
                            %lmval =[lmval x(i)];	%1   comment these two lines for strict case
                            %indd = [ indd i];	%2 when only  definite min included
                            i = i + 1;		% skip one point
                                   end
                                end
                                i = i + 1;
                                end
                           if filt>0 && ~isempty(indd)
                                if (indd(1)<= 3)||(indd(length(indd))+2>length(xx)) 
                                   rng=1;	%check if index too close to the edge
                                else
                                    rng=2;
                                end

                                   for ii=1:length(indd) 
                                       a=0;b=0;
                                    [val(ii), iind(ii)] = min(xx(indd(ii) -rng:indd(ii) +rng));
                                    iind(ii)=indd(ii) + iind(ii)  -rng-1;
                                   end
                              indd=iind; lmval=val;
                            else
                            end
                        end     
                        
                        function [lmval,indd]=lmax(xx,filt)
                            %Find the local minima in a data set, excludes the first
                            %and last points.
                            % Created by Serge Koptenko, Guigne International Ltd.
                            x=xx;
                            len_x = length(x);
                                fltr=[1 1 1]/3;
                              if nargin <2, filt=0; 
                                else
                            x1=x(1); x2=x(len_x); 

                                for jj=1:filt
                                c=conv(fltr,x);
                                x=c(2:len_x+1);
                                x(1)=x1;  
                                    x(len_x)=x2; 
                                end
                              end

                            lmval=[];
                            indd=[];
                            i=2;		% start at second data point in time series

                                while i < len_x-1
                                if x(i) > x(i-1)
                                   if x(i) > x(i+1)	% definite min
                            lmval =[lmval x(i)];
                            indd = [ indd i];

                                   elseif x(i)==x(i+1)&&x(i)==x(i+2)	% 'long' flat spot
                            %lmval =[lmval x(i)];	%1   comment these two lines for strict case 
                            %indd = [ indd i];	%2 when only  definite min included
                            i = i + 2;  		% skip 2 points

                                   elseif x(i)==x(i+1)	% 'short' flat spot
                            %lmval =[lmval x(i)];	%1   comment these two lines for strict case
                            %indd = [ indd i];	%2 when only  definite min included
                            i = i + 1;		% skip one point
                                   end
                                end
                                i = i + 1;
                                end
                           if filt>0 && ~isempty(indd)
                                if (indd(1)<= 3)||(indd(length(indd))+2>length(xx)) 
                                   rng=1;	%check if index too close to the edge
                                else
                                    rng=2;
                                end

                                   for ii=1:length(indd) 
                                    [val(ii), iind(ii)] = max(xx(indd(ii) -rng:indd(ii) +rng));
                                    iind(ii)=indd(ii) + iind(ii)  -rng-1;
                                   end
                              indd=iind; lmval=val;
                            else
                            end
                        end     
                        
                         
        end
            
        function pixelStrainCallback(tool,~,~)
                 tool.pointLog = tool.pointLog;
                 if 1%(isempty(tool.pointLog))
                     tool.pointLog=[];
                     msgbox('Tracking all points'); pause(1); close;
                     TrackAll(tool);
                 end
                   dicomFrames = size(tool.I,3);
                   choice = questdlg('Which measurement would you like to display?', ...
                        'Select strain type', 'Displacement','Strain','Strain');
                   
                    if ~isempty(tool.currentROI)                 
                                  if isvalid(tool.currentROI)
                                        pos = round(getPosition(tool.currentROI));
                                        pixelsX = pos(3); pixelsY = pos(4);
                                  else
                                      pixelsX =size(tool.I,2); pixelsY = size(tool.I,1);    
                                  end
                     else
                         pixelsX =size(tool.I,2); pixelsY = size(tool.I,1);    
                     end
                      
                    pixelsXtracked = round(pixelsX*tool.pixelDensity/100);
                    pixelsYtracked = round(pixelsY*tool.pixelDensity/100);
                    trackedPixels = pixelsXtracked*pixelsYtracked;
                    
                    startFrame = 1;
                    endFrame = size(tool.pointLog,3);
                    
%%   DaeWoo
                    switch choice
                       case 'Displacement'
                           if endFrame <= 1
                               warndlg('Need multiple frames to measure displacement');
                               return;
                           end
                           count = 1;
                           counter = 1;
                           counter1 = 1;
                           pixelsYtracked = size(tool.pointLog(tool.pointLog(:,1,1)==tool.pointLog(1,1,1)));
                           rectpoints= [];
                           %Here the points are separated out into a
                           %rectangular matrix structured more like the
                           %image so they are more intuitive. X and Y 
                           %displacements are then stored 
                           for indFrame = startFrame:endFrame
                               for ind = 1:pixelsYtracked:size(tool.pointLog)
                                   rectpoints(:,counter,indFrame) = tool.pointLog(ind:(ind+pixelsYtracked-1),1,indFrame);
                                   rectpoints(:,counter+1,indFrame) = tool.pointLog(ind:(ind+pixelsYtracked-1),2,indFrame);
                                   %the scores can be used to filter the
                                   %tracked points for quality
                                   myboxscores(:,counter1,indFrame) = tool.pointLogScores(ind:(ind+pixelsYtracked-1),1,indFrame);
                                   counter = counter+2;
                                   counter1 = counter1+1;
                               end
                               counter1=1;
                               counter = 1;
                               if indFrame > 1
                                  xdisplacement(:,:,indFrame) = rectpoints(:,1:2:end,indFrame) - rectpoints(:,1:2:end,1);
                                  ydisplacement(:,:,indFrame) = rectpoints(:,2:2:end,indFrame) - rectpoints(:,2:2:end,1);
                               end
                           end
           %                cm = tool.calibration;
            %               xdisplacement = xdisplacement.*cm;
             %              ydisplacement = ydisplacement.*cm;
                           totaldisplacement = sqrt(xdisplacement.^2 + ydisplacement.^2);
                           counter = size(rectpoints,3);%this plots only the final displacement without having it as a fucntion of time
                           figure
                       %    subplot(1,3,1)
                        %   imagesc([rectpoints(1,1,counter),rectpoints(end,end-1,counter)],[rectpoints(1,2,counter),rectpoints(end,end,counter)], totaldisplacement(:,:,counter))
                         %  axis off
                        %   title('Total Displacement')
                          % colorbar
                       %    subplot(1,3,2)
                           imagesc([rectpoints(1,1,counter),rectpoints(end,end-1,counter)],[rectpoints(1,2,counter),rectpoints(end,end,counter)], xdisplacement(:,:,counter))
                           axis off
                           title('X Displacement')
                           colorbar
                  %         subplot(1,3,3)
                   %        imagesc([rectpoints(1,1,counter),rectpoints(end,end-1,counter)],[rectpoints(1,2,counter),rectpoints(end,end,counter)], ydisplacement(:,:,counter))
                    %       axis off
                     %      title('Y Displacement')
                      %     colorbar
                          
                          % save Displacement.mat;
%                        case 'OldFrame to Frame'
%                             % Simple difference
%                             for indFrames = startFrame:endFrame
%                                 pointLogDiff(:,:,indFrames) = tool.pointLog(:,:,indFrames+1) ...
%                                     - tool.pointLog(:,:,indFrames);       
%                             end
                        case 'Strain'
                           if endFrame <= 1
                               warndlg('Need multiple frames to measure strain');
                               return;
                           end

                           %Here the points are separated out into a
                           %rectangular matrix structured more like the
                           %image so they are more intuitive. 
                           counter = 1;
                           myboxes=[];
                           pixelsYtracked = size(tool.pointLog(tool.pointLog(:,1,1)==tool.pointLog(1,1,1)));
                           counter1=1;
                           for indFrame = startFrame:endFrame
                               for ind = 1:pixelsYtracked:size(tool.pointLog,1)
                                   myboxes(:,counter,indFrame) = tool.pointLog(ind:(ind+pixelsYtracked-1),1,indFrame);
                                   myboxes(:,counter+1,indFrame) = tool.pointLog(ind:(ind+pixelsYtracked-1),2,indFrame);
                                   counter = counter+2;
                                   myboxscores(:,counter1,indFrame) = tool.pointLogScores(ind:(ind+pixelsYtracked-1),1,indFrame);
                                   counter1= counter1+1;
                               end
                               counter = 1;
                               counter1=1;
                           end

                           %Strain Looper 
                           %For each frame, column, and row of boxes, the
                           %midpoints of the sizes are determined to get x
                           %and y dimensions of the boxes as they deform.
                           %These are used to calculate the strain. Vectors
                           %are calculated from the bottome left corner of
                           %each box (P2o) to get angles for shear strain. The
                           %lengths are recorded of the box diagonal (xls etc.) to
                           %later be aggregated into circumference. The
                           %mean of the corner points is also calculated to
                           %act as the center of the box which will be used
                           %to later plot representative box strain values
                           finalLs= size(myboxes,3);
                           for framecount = 1:size(myboxes,3)-1
                            for colcount = 1:2:size(myboxes,2)-2
                               for rowcount = 1:size(myboxes,1)-1
                                   topl1 = [myboxes(rowcount,colcount,1),myboxes(rowcount,colcount+1,1)];
                                   btml1 = [myboxes(rowcount+1,colcount,1),myboxes(rowcount+1,colcount+1,1)];
                                   btmr1 = [myboxes(rowcount+1,colcount+2,1),myboxes(rowcount+1,colcount+3,1)];
                                   topl2 = [myboxes(rowcount,colcount,framecount+1),myboxes(rowcount,colcount+1,framecount+1)];
                                   btml2 = [myboxes(rowcount+1,colcount,framecount+1),myboxes(rowcount+1,colcount+1,framecount+1)];
                                   btmr2 = [myboxes(rowcount+1,colcount+2,framecount+1),myboxes(rowcount+1,colcount+3,framecount+1)];
                                   xStrain(rowcount,(colcount+1)./2,framecount+1) = ((((myboxes(rowcount,colcount+2,framecount+1) + myboxes(rowcount+1,colcount+2,framecount+1))./2) - ((myboxes(rowcount,colcount,framecount+1) + myboxes(rowcount+1,colcount,framecount+1))./2))  ...
                                     - (((myboxes(rowcount,colcount+2,1)+myboxes(rowcount+1,colcount+2,1))./2) - ((myboxes(rowcount,colcount,1) + myboxes(rowcount+1,colcount,1))./2)))./(((myboxes(rowcount,colcount+2,1)+myboxes(rowcount+1,colcount+2,1))./2) - ((myboxes(rowcount,colcount,1) + myboxes(rowcount+1,colcount,1))./2));        
                                   yStrain(rowcount,(colcount+1)./2,framecount+1) = ((((myboxes(rowcount+1,colcount+1,framecount+1)+myboxes(rowcount+1,colcount+3,framecount+1))./2) - ((myboxes(rowcount,colcount+3,framecount+1) +myboxes(rowcount,colcount+1,framecount+1))./2))...
                                     - (((myboxes(rowcount+1,colcount+1,1)+myboxes(rowcount+1,colcount+3,1))./2) - ((myboxes(rowcount,colcount+3,1) + myboxes(rowcount,colcount+1,1))./2)))./(((myboxes(rowcount+1,colcount+1,1)+myboxes(rowcount+1,colcount+3,1))./2) - ((myboxes(rowcount,colcount+3,1) + myboxes(rowcount,colcount+1,1))./2));
                                   ShearStrain(rowcount,(colcount+1)./2,framecount+1) = acos(dot((topl1-btml1),(btmr1-btml1))./norm(topl1-btml1)./norm(btmr1-btml1))-acos(dot((topl2-btml2),(btmr2-btml2))./norm(topl2-btml2)./norm(btmr2-btml2));
                                   xcorners= [myboxes(rowcount,colcount,framecount+1),myboxes(rowcount+1,colcount,framecount+1),myboxes(rowcount,colcount+2,framecount+1),myboxes(rowcount+1,colcount+2,framecount+1)];
                                   ycorners= [myboxes(rowcount,colcount+1,framecount+1),myboxes(rowcount+1,colcount+1,framecount+1),myboxes(rowcount,colcount+3,framecount),myboxes(rowcount+1,colcount+3,framecount+1)];
                                   center = means(xcorners,ycorners);
                                   %modification starts
                                   xStrain(1,:,:) = NaN;
                                   xStrain(size(myboxes,1)-1,:,:) = NaN;
                                   
                                   
                                   yStrain(1,:,:) = NaN;
                                   yStrain(size(myboxes,1)-1,:,:) = NaN;
                                   
                                   
                                   ShearStrain(1,:,:) = NaN;
                                   ShearStrain(size(myboxes,1)-1,:,:) = NaN;
                                   
                                   
                                   %modification ends
                                   xdim(rowcount,(colcount+1)./2,framecount)=(((myboxes(rowcount,colcount+2,framecount)+myboxes(rowcount+1,colcount+2,framecount))./2) - ((myboxes(rowcount,colcount,framecount) + myboxes(rowcount+1,colcount,framecount))./2));
                                   ydim(rowcount,(colcount+1)./2,framecount)=(((myboxes(rowcount+1,colcount+1,framecount)+myboxes(rowcount+1,colcount+3,framecount))./2) - ((myboxes(rowcount,colcount+3,framecount) + myboxes(rowcount,colcount+1,framecount))./2));
                                   xmidpoints(rowcount,(colcount+1)./2,framecount+1)= center(1);
                                   ymidpoints(rowcount,(colcount+1)./2,framecount+1)= center(2);
                                   xdim(rowcount,(colcount+1)./2,finalLs)=(((myboxes(rowcount,colcount+2,finalLs)+myboxes(rowcount+1,colcount+2,finalLs))./2) - ((myboxes(rowcount,colcount,finalLs) + myboxes(rowcount+1,colcount,finalLs))./2));
                                   ydim(rowcount,(colcount+1)./2,finalLs)=(((myboxes(rowcount+1,colcount+1,finalLs)+myboxes(rowcount+1,colcount+3,finalLs))./2) - ((myboxes(rowcount,colcount+3,finalLs) + myboxes(rowcount,colcount+1,finalLs))./2));
                               end
                               %modification starts
                               xStrain(:,1,:) = NaN;
                               %xStrain(:,size(myboxes,2)-2,:) = NaN;
                               yStrain(:,1,:) = NaN;
                               %yStrain(:,size(myboxes,2)-2,:) = NaN;
                               ShearStrain(:,1,:) = NaN;
                               %ShearStrain(:,size(myboxes,2)-2,:) = NaN;
                               %modificiation ends
                            end
                            xStrain(:,(colcount+1)./2,:) = NaN;
                            yStrain(:,(colcount+1)./2,:) = NaN;
                            ShearStrain(:,(colcount+1)./2,:) = NaN;
                           end
                           
                           StrainMagnitude = (xStrain.^2 +yStrain.^2).^0.5;
                           diaglength= (xdim.^2+ydim.^2).^0.5;
                           %The user is prompted to create a new ROI to be able to
                           %filter out strains within the vessel to be
                           %zero.
                          imgholder = tool.I;
                           choice = questdlg('Would you like to filter out undesired elements?', ...
                        'Filter Vessel', 'Yes','No','No');
                    switch choice
                        case 'Yes'
                            uiwait(msgbox('Change the ROI to select undesired element. Then, click "OK"'));
                            if ~isempty(tool.currentROI)                 
                                  if isvalid(tool.currentROI)
                                        newROIpos = round(getPosition(tool.currentROI));
                                  else
                                      newROIpos = [0,0,0,0];    
                                  end
                            else
                                 newROIpos = [0,0,0,0];
                            end
                           imgholder(newROIpos(2)-5:newROIpos(2)+newROIpos(4)+5,newROIpos(1)-5:newROIpos(1)+newROIpos(3)+5,:)= 0;
                        case 'No'
                    end
                    prompt = {'Input filter tolerance:'};
                    dlg_title = 'Input';
                    num_lines = 1;
                    default = {'30'};
                    options.Resize='on';
                    options.WindowStyle='normal';
                    answer = inputdlg(prompt,dlg_title,num_lines,default,options);
                    filt = str2double(answer{1,1});
            %FAILED LOOP FOR CIRCUMFERENTIAL STRAIN. FAILURE TO REFERERENCE
            %ORIGINAL FRAME. AQUIRES NEW POINTS AS IT GOES.
                           %instead of prompting the user, we set the
                           %strain inside the vessel to zero by defining a
                           %tolerance for pixel values to determine which
                           %pixels are inside the vessel. Then, strains and
                           %lengths around the edge of that intensity
                           %defined vessel wall will be aggregated to get a
                           %full perimeter value.
%                           vesselpoints={[],[]};circum0=0; circumend=0; strainlines={[],[]};selected={[],[]};
 %                          xcircstrain= zeros(1,size(xmidpoints,3));ycircstrain= zeros(1,size(xmidpoints,3));magcircstrain=zeros(1,size(xmidpoints,3));
  %                         for framecounter = 2:size(xmidpoints,3)
   %                         i=1;
    %                        j=1;
                            %the segment below uses the size of a matrix to
                            %get offsets that can be coupled with an index
                            %to provide the indexed value and the points
                           %directly surrounding it. a multiplier was used
                           %to expand this region.
%                            xmid= xmidpoints(:,:,framecounter);
%                            ymid= ymidpoints(:,:,framecounter);
%                            s=size(xmidpoints(:,:,1));
%                            N=length(s);
%                            [c1{1:N}]=ndgrid(1:3);
%                            c2(1:N)={2};
%                            offsets=(sub2ind(s,c1{:}) - sub2ind(s,c2{:}));
%                            image= imgholder(:,:,framecounter);
%                            s1=size(image);
%                           N1=length(s);
%                             [c11{1:N1}]=ndgrid(1:3);% MAYBE ALTER THE 3 TO GET BIGGER NEIGHBORS
%                             c21(1:N1)={2};
%                             offsets1=(sub2ind(s1,c11{:}) - sub2ind(s1,c21{:})).*10;
%                             c=1;
%                             gunpoint=[];
%                            for colcounter = 2:size(xmidpoints,2)-1
%                                 for rowcounter = 2:size(xmidpoints,1)-1
%                                     x= uint16(round(xmidpoints(rowcounter,colcounter,framecounter)));
%                                     y= uint16(round(ymidpoints(rowcounter,colcounter,framecounter)));
%                                     neighbors = image(sub2ind(s1,y,x)+offsets1);
%                                     neighborsx= xmid(sub2ind(s,rowcounter,colcounter)+offsets);
%                                     neighborsy= ymid(sub2ind(s,rowcounter,colcounter)+offsets);
%                                     if imgholder(y,x,framecounter)<filt && sum(sum(neighbors<filt)) >=2
%                                         xStrain(rowcounter,colcounter,framecounter) = 0;
%                                         yStrain(rowcounter,colcounter,framecounter) = 0;
%                                         ShearStrain(rowcounter,colcounter,framecounter) = 0;
%                                         StrainMagnitude(rowcounter,colcounter,framecounter) = 0;
%                                    end
%                                    if imgholder(y,x,framecounter)~=0 && sum(sum(neighbors<filt))>= 1 && sum(sum(neighbors<filt)) <= 5
%                                       if i>2
%                                        hits=0;
%                                        tol = eps(.5);
%                                        for p= 1:size(neighborsx,1)*size(neighborsx,2)
%                                           [in,~]= find(abs(vesselpoints{framecounter}(:,1)-neighborsx(p))<=tol);
%                                           if ~isempty(in)
%                                               h=0;
%                                           end
%                                           if abs(vesselpoints{framecounter}(in,2)-neighborsy(p))<=tol
%                                           hits= hits+1;
%                                           end
%                                        end
%                                         if hits>=3 
%                                            gunpoint(c,1:2)= [xmidpoints(rowcounter,colcounter,framecounter),ymidpoints(rowcounter,colcounter,framecounter)];
%                                            c=c+1;
%                                           continue 
%                                         end
%                                        end
%                                        if framecounter == 2
%                                             circum0 = circum0+diaglength(rowcounter,colcounter,1);
%                                             selected{2}(rowcounter,colcounter)=diaglength(rowcounter,colcounter,1);
%                                             selected{1}(i,:)= diaglength(rowcounter,colcounter,1);
%                                             strainlines{framecounter}(i,:)= [xmidpoints(rowcounter,colcounter,framecounter)+xdim(rowcounter,colcounter,1)/2,ymidpoints(rowcounter,colcounter,framecounter)-ydim(rowcounter,colcounter,1)/2,...
%                                                 xmidpoints(rowcounter,colcounter,framecounter)-xdim(rowcounter,colcounter,1)/2,ymidpoints(rowcounter,colcounter,framecounter)+ydim(rowcounter,colcounter,1)];
%                                        end
%                                        if framecounter == size(xmidpoints,3)
%                                             circumend = circumend+diaglength(rowcounter,colcounter,end);
%                                             selected{framecounter}(rowcounter,colcounter)= diaglength(rowcounter,colcounter,end);
%                                             selected{framecounter-1}(i,:)= diaglength(rowcounter,colcounter,end);
%                                             strainlines{framecounter}(i,:)= [xmidpoints(rowcounter,colcounter,framecounter)+xdim(rowcounter,colcounter,end)/2,ymidpoints(rowcounter,colcounter,framecounter)-ydim(rowcounter,colcounter,end)/2,...
%                                                 xmidpoints(rowcounter,colcounter,framecounter)-xdim(rowcounter,colcounter,end)/2,ymidpoints(rowcounter,colcounter,framecounter)+ydim(rowcounter,colcounter,end)];
%                                        end
%                                        xcircstrain(framecounter) = xcircstrain(framecounter) + xStrain(rowcounter,colcounter,framecounter);
%                                        ycircstrain(framecounter) = ycircstrain(framecounter) + yStrain(rowcounter,colcounter,framecounter);
%                                        magcircstrain(framecounter) = magcircstrain(framecounter) + StrainMagnitude(rowcounter,colcounter,framecounter);
%                                        vesselpoints{framecounter}(i,:) = [xmidpoints(rowcounter,colcounter,framecounter),ymidpoints(rowcounter,colcounter,framecounter)];
% 
%                                        i=i+1;
%                                        j=j+1;
%                                    end
%                                 end
%                            end
%                            end
%                            ystrainav=0;
%                            xstrainav=0;
                %SUCCESSFUL LOOP FOR CIRCUMFERENTIAL STRAINS
                          vesselpoints={[],[]};circum0=0; circumend=0; strainlines={[],[]};selected={[],[]};
                          selectedrows= [];
                          selectedcols= [];
                          framecounter=2;
                          i=1;
                           %the segment below uses the size of a matrix to
                           %get offsets that can be coupled with an index
                           %to provide the indexed value and the points
                           %directly surrounding it. a multiplier was used
                           %to expand this region.
                          xmid= xmidpoints(:,:,framecounter);
                          ymid= ymidpoints(:,:,framecounter);

                          s=size(xmidpoints(:,:,1));
                          N=length(s);
                          [c1{1:N}]=ndgrid(1:3);
                          c2(1:N)={2};
                          offsets=(sub2ind(s,c1{:}) - sub2ind(s,c2{:}));
                          image= imgholder(:,:,framecounter);
                          s1=size(image);
                          N1=length(s);
                          [c11{1:N1}]=ndgrid(1:19);
                          c21(1:N1)={10};
                          offsets1=(sub2ind(s1,c11{:}) - sub2ind(s1,c21{:}));
                          c=1;
                          o=1;
                          gunpoint=[];
                          for colcounter = 2:size(xmidpoints,2)-1
                               for rowcounter = 2:size(xmidpoints,1)-1
                                   x= uint16(round(xmidpoints(rowcounter,colcounter,framecounter)));
                                   y= uint16(round(ymidpoints(rowcounter,colcounter,framecounter)));
                                   neighbors = image(sub2ind(s1,y,x)+offsets1);
                                   neighborsx= xmid(sub2ind(s,rowcounter,colcounter)+offsets);
                                   neighborsy= ymid(sub2ind(s,rowcounter,colcounter)+offsets);
                                   if imgholder(y,x,framecounter)<filt && sum(sum(neighbors<filt)) >=2
                                       xStrain(rowcounter,colcounter,end) = 0;
                                       yStrain(rowcounter,colcounter,end) = 0;
                                       ShearStrain(rowcounter,colcounter,end) = 0;
                                       StrainMagnitude(rowcounter,colcounter,end) = 0;
                                    end
                                   if imgholder(y,x,framecounter)~=0 && sum(sum(neighbors<filt))>= 21 && sum(sum(neighbors<filt)) <= 200
                   %                    boxscore= [myboxscores(rowcounter,colcounter,2:end),myboxscores(rowcounter+1,colcounter,2:end),...
                    %                       myboxscores(rowcounter,colcounter+1,2:end), myboxscores(rowcounter+1,colcount+1,2:end)];
                     %                  if any(boxscore==0)
                      %                     cutpoint(o,1:2)=[xmidpoints(rowcounter,colcounter,framecounter),ymidpoints(rowcounter,colcounter,framecounter)];
                       %                    o=o+1;
                        %               end
                                       if i>2
                                       hits=0;
                                       tol = eps(.5);
                                       for p= 1:size(neighborsx,1)*size(neighborsx,2)
                                          [in,~]= find(abs(vesselpoints{framecounter}(:,1)-neighborsx(p))<=tol);
                                          if ~isempty(in)
                                              h=0;
                                          end
                                          if abs(vesselpoints{framecounter}(in,2)-neighborsy(p))<=tol
                                          hits= hits+1;
                                          end
                                       end
                                       if hits>=3 
                                          gunpoint(c,1:2)= [xmidpoints(rowcounter,colcounter,framecounter),ymidpoints(rowcounter,colcounter,framecounter)];
                                          c=c+1;
                                         continue 
                                       end
                                      end
                                      circum0 = circum0+diaglength(rowcounter,colcounter,1);
                                      selected{2}(rowcounter,colcounter)=diaglength(rowcounter,colcounter,1);
                                      selected{1}(i,:)= diaglength(rowcounter,colcounter,1);
                                      selectedrows(i)= rowcounter;
                                      selectedcols(i)= colcounter;
                                      strainlines{framecounter}(i,:)= [xmidpoints(rowcounter,colcounter,framecounter)+xdim(rowcounter,colcounter,1)/2,ymidpoints(rowcounter,colcounter,framecounter)-ydim(rowcounter,colcounter,1)/2,...
                                         xmidpoints(rowcounter,colcounter,framecounter)-xdim(rowcounter,colcounter,1)/2,ymidpoints(rowcounter,colcounter,framecounter)+ydim(rowcounter,colcounter,1)];
                                      vesselpoints{framecounter}(i,:) = [xmidpoints(rowcounter,colcounter,framecounter),ymidpoints(rowcounter,colcounter,framecounter)];
                                      i=i+1;
                                  end
                               end
                          end
                          i=1;
                          selectedrows=selectedrows.';
                          selectedcols=selectedcols.';
                          selpts=[selectedrows,selectedcols];
                          xstrainav=0;ystrainav=0;
                          framecounter = size(xmidpoints,3);
                          
                          for cnt= 1:size(selpts,1)
                                      circumend = circumend+diaglength(selpts(cnt,1),selpts(cnt,2),end);
                                      selected{framecounter}(selpts(cnt,1),selpts(cnt,2))= diaglength(selpts(cnt,1),selpts(cnt,2),end);
                                      selected{framecounter-1}(i,:)= diaglength(selpts(cnt,1),selpts(cnt,2),end);
                                      strainlines{framecounter}(i,:)= [xmidpoints(selpts(cnt,1),selpts(cnt,2),framecounter)+xdim(selpts(cnt,1),selpts(cnt,2),end)/2,ymidpoints(selpts(cnt,1),selpts(cnt,2),framecounter)-ydim(selpts(cnt,1),selpts(cnt,2),end)/2,...
                                        xmidpoints(selpts(cnt,1),selpts(cnt,2),framecounter)-xdim(selpts(cnt,1),selpts(cnt,2),end)/2,ymidpoints(selpts(cnt,1),selpts(cnt,2),framecounter)+ydim(selpts(cnt,1),selpts(cnt,2),end)];
                                      vesselpoints{framecounter}(i,:) = [xmidpoints(selpts(cnt,1),selpts(cnt,2),framecounter),ymidpoints(selpts(cnt,1),selpts(cnt,2),framecounter)];
                                      xstrainav=xstrainav+xStrain(selpts(cnt,1),selpts(cnt,2),end);
                                      ystrainav=ystrainav+yStrain(selpts(cnt,1),selpts(cnt,2),end);
                                      i=i+1;
                          end
                          strainsmatrix= zeros(size(StrainMagnitude,1),size(StrainMagnitude,2),size(StrainMagnitude,3),4);
                          strainsmatrix(:,:,:,1) = StrainMagnitude;
                          strainsmatrix(:,:,:,2) = xStrain;
                          strainsmatrix(:,:,:,3) = yStrain;
                          strainsmatrix(:,:,:,4) = ShearStrain;
                          ystrainav=ystrainav/(i-1);
                          xstrainav=xstrainav/(i-1);
                    %The following code is used to create images showing
                    %where the algorithm has chosen points to aggregate
                    %circumfrential strains and lengths around the vessel
                    dicomFrames = size(tool.I,3);
                    newI = uint8(tool.I); 
                    J = uint8(tool.I);
                    v={};
                    for t= 1:size(vesselpoints{2},1)
                        v{t}= 'cyan';
                    end
                    if ~isempty(gunpoint)
                    for u= 1:size(gunpoint,1)
                       v{t+u}= 'red';
                    end
                    vesselpoints{2}=[vesselpoints{2};gunpoint];
                    end
        %            if ~isempty(cutpoint)
         %          for r= 1:size(cutpoint,1)
          %             v{t+u+r}= 'magenta';
           %        end
            %       vesselpoints{2} = [vesselpoints{2};cutpoint];
             %       end
               %     pointImage= insertShape(J(:,:,1),'Line',strainlines{2} ,'Color','blue');
                    pointImage = insertMarker(J(:,:,1), vesselpoints{2}, '+', 'Color', v);
                    pointImage2 = insertMarker(J(:,:,end), vesselpoints{end}, '+','Color','white');
                %    pointImage2= insertShape(J(:,:,end),'Line',strainlines{end} ,'Color','blue');
                    figure
                    imshow(pointImage(:,:,:))
                    title('startframe')
                    figure
                    imshow(pointImage2(:,:,:))
                    title('endframe')
%                     g=1;
%                     for j = 1:size(xmidpoints,2)
%                         midplot(g:g+size(xmidpoints,1)-1,1)= xmidpoints(:,j,2);
%                         g=g+size(xmidpoints,1);
%                     end
%                     g=1;
%                     for j = 1:size(ymidpoints,2)
%                         midplot(g:g+size(ymidpoints,1)-1,2) = ymidpoints(:,j,2);
%                         g= g+size(ymidpoints,1);
%                     end
%                     pointImage3 = insertMarker(J(:,:,2), midplot, '+','Color','white');
%                     figure
%                     imshow(pointImage3(:,:,1))
%                     title('xmid')
                            DeltL = [circum0,circumend,(circumend-circum0)/circum0.*100];
                            strainav=[xstrainav,ystrainav];
                            %The user will select which data they want then
                            %imagesc will plot a color coded image of the
                            %representative strain values for each box over the
                            %selected region. These will be the resultant total
                            %strain from the first to last frame of the
                            %image. 
                            save('Strain','strainsmatrix','DeltL','strainav')
                            choice2 = questdlg('Which measurement would you like to display?', ...
                        'Select strain type', 'Magnitude','Specific Strain','All','All');
                            switch choice2
                                case 'Magnitude'
                                    chosenmatrix = StrainMagnitude;
                                    mytitle = 'Strain Magnitude';
                                    counter = size(StrainMagnitude,3);
                                    figure
                                    imagesc([xmidpoints(1,1,counter),xmidpoints(end,end,counter)],[ymidpoints(1,1,counter),ymidpoints(end,end,counter)], chosenmatrix(:,:,counter))
                                    axis off
                                    title(mytitle)
                                    colorbar
                                case 'Specific Strain'
                                    choice3 = questdlg('Which strain would you like to display?',...
                                        'Select Strain Type','X','Y','Shear','Shear');
                                    switch choice3
                                        case 'X'
                                                chosenmatrix= xStrain;
                                                mytitle = 'X Strain';
                                        case 'Y'
                                                 chosenmatrix= yStrain;
                                                 mytitle = 'Y Strain';
                                       case 'Shear'
                                                chosenmatrix = ShearStrain;
                                                mytitle = 'Shear Strain';
                                    end
                                    counter = size(StrainMagnitude,3);
                                    figure
                                    imagesc([xmidpoints(1,1,counter),xmidpoints(end,end,counter)],[ymidpoints(1,1,counter),ymidpoints(end,end,counter)], chosenmatrix(:,:,counter))
                                    axis off
                                    title(mytitle)
                                    colorbar
                                case 'All'
                                    myfig = figure; set(myfig,'Visible','off');
                                    subplot(2,2,1)
                                    subnum=1;
                                        for index = 1:4
                                        subplot(2,2,subnum)
                                        imagesc([xmidpoints(2,2,end),xmidpoints(end-1,end-1,end)],[ymidpoints(2,2,end),ymidpoints(end-1,end-1,end)], strainsmatrix(:,:,end,index)) %,[min(min(min(min(strainmatrix(:,:,:,2:4))))),max(max(max(max(strainmatrix(:,:,:,2:4)))))])
                                        axis off
                                        colorbar
                                        subnum = subnum + 1;
                                        if isequal(index,1)
                                            title('Magnitude')
                                        elseif isequal(index,2)
                                            title('X Strain')
                                        elseif isequal(index,3)
                                            title ('Y Strain')
                                        elseif isequal(index,4)
                                            title('Shear Strain')
                                        end
                                        end
                                        %this code and the part above after imagesc above
                                        %will create one colorbar for
                                        %the entire figure however it does
                                        %not differentiate the small
                                        %strains as well
                                        %hp4 = get(subplot(2,2,4),'Position');
                                        %colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.1  hp4(2)+hp4(3)*2.1])
                                        %colorbar('Location','eastoutside')
                                        %caxis([min(min(min(min(strainmatrix(:,:,:,2:4))))),max(max(max(max(strainmatrix(:,:,:,2:4)))))])
                                        prompt = {'Patient Number:', 'Vessel Type:'};
                                        dlg_title = 'Plot Data';
                                        num_lines = [1, 60; 1, 60;];
                                        default = {'#','Blood'};
                                        options.Resize='on';
                                        options.WindowStyle='normal';
                                        answer = inputdlg(prompt,dlg_title,num_lines,default,options);
                                        if ~isempty(answer)
                                            set(myfig,'Visible','on')
                                            suptitle(['Patient ',answer{1,1},' ',answer{2,1}])
                                        end
                                        set(myfig,'Visible','on')
                                        
                            end

                            %The following commented out can be used
                            %instead of the above plotting and plotting
                            %loops. It will plot a figure for the strain
                            %magnitude on each box for all of the frames up to
                            %the current frame as it cycles through all of
                            %the frames thus displaying strain over time.
                            %However a reasonable manner to display those
                            %figures was not found for this language.
                            
                            %StrainFrames = cell(size(overallStrain),3);
                            %mywaitbar = waitbar(0,'Processing Strain');
                            %for counter = 2:size(overallStrain,3)
                            %myfig = figure; set(gcf,'Visible','off');
                            %color = overallStrain(:,:,counter)./max(max(overallStrain)).*255;
                            %imagesc([xmidpoints(1,1,counter),xmidpoints(end,end,counter)],[ymidpoints(1,1,counter),ymidpoints(end,end,counter)], chosenmatrix(:,:,counter))
                            %axis([pos(1),pos(1)+pos(3),pos(2),pos(2)+pos(4)])
                            
                            %THIS IS WHERE THE imagesc PLOTTER GOES
                            
                            %L = getframe;
                            %[StrainFrames(:,:,:,counter),~] = frame2im(L);
                            %close
                            %waitbar(counter/size(overallStrain,3))
                            %end
                            %close(mywaitbar)


                    end

                
                         function TrackAll(tool)

                            
                            firstFrame = tool.accFrames(1); 
                            lastFrame = tool.accFrames(2);
                            dicomFrames = lastFrame - firstFrame + 1;
                            newI = uint8(tool.I); 
                             J = uint8(tool.I);
                             H = zeros(size(J,1),size(J,2),3,lastFrame);
                            if ~isempty(tool.currentROI)                 
                                  if isvalid(tool.currentROI)
                                        pos = round(getPosition(tool.currentROI));
                                  end
                            end
                            %Create grid of points on the image                    
                            if tool.pixelDensity >100
                                tool.pixelDensity = 100;
                            elseif tool.pixelDensity <=0
                                tool.pixelDensity = 1;
                            end

                            if ~isempty(tool.currentROI)                 
                                  if isvalid(tool.currentROI)
                                      %If there is an roi get dimensions
                                      pixelsX = pos(3); pixelsY = pos(4);
                                      offsetX = round(pos(1)); offsetY = round(pos(2));
                                  else
                                      %If no roi get dimensions for entire image
                                      pixelsX =size(tool.I,2); pixelsY = size(tool.I,1);
                                      offsetX = 0.001; offsetY = 0.001;                              
                                  end
                            else
                                %If no roi get dimensions for entire image
                                pixelsX =round(size(tool.I,2)); pixelsY = round(size(tool.I,1));
                                offsetX = .0001; offsetY = .0001;   
                            end
                            % Find pixel spacing using pixel density factor (tool.pixelDensity)
                            pixelsBetweenX = 5;%(pixelsX-1)/round((pixelsX-1)*tool.pixelDensity/100); 
                            pixelsBetweenY = 5;%(pixelsY-1)/round((pixelsY-1)*tool.pixelDensity/100); 
                            count = 1;
                            countX = 1+offsetX;
                            % We get a grid that is %PixelDensity^2*(pixelsX*pixelsY)
                            points = zeros(1,2);
                            while countX <= pixelsX+offsetX
                                countY=1+offsetY;
                                while countY <= pixelsY+offsetY
                                    points(count,:) = [round(countX) round(countY)];
                                    countY = countY + pixelsBetweenY;
                                    count = count+1;
                                end
                                countX = countX + pixelsBetweenX;
                            end
                            nPoints = count - 1;
                            tool.pointLog = zeros(nPoints, 2, dicomFrames);
                            tool.pointLog(:,:,1) = points;
                            framenum = firstFrame+1;
                            objectFrame = newI(:,:,firstFrame);
                            pointImage = insertMarker(objectFrame, points, '+', 'Color', 'white');
                            newI(:,:,firstFrame) = pointImage(:,:,1);
                            quality = ones(1,dicomFrames);
                            tool.pointLogScores = zeros(nPoints,1,dicomFrames);
                            % Create object tracker
                            tracker = vision.PointTracker('MaxBidirectionalError', 3);
                            ii = 2;
                            % Initialize object tracker
                            initialize(tracker, points(:,:,1), objectFrame);
                            h = waitbar(0,'Running pixel tracker...');
                            % Show the points getting tracked
                            filtered=zeros(1,lastFrame);
                            losses=zeros(1,lastFrame);
                            po = 1;
                            workingDir = tempname;
                            mkdir(workingDir)
                            mkdir(workingDir,'images')
                            while framenum <= lastFrame 
                                 %Track the points and set the ones that are invalid to never move  
                                  frame =J(:,:,framenum);
                                  [points, validity,scores] = step(tracker, frame);
                                  tool.pointLog(:,:,ii) = points;
                                  tool.pointLogScores(:,1,ii) = scores;
                                  if framenum~=1
                                  [a,~]= find(validity==0);
                                  tool.pointLog(a,1,ii) = tool.pointLog(a,1,ii-1);
                                  tool.pointLog(a,2,ii) = tool.pointLog(a,2,ii-1);
                                  filtered(framenum)= filtered(framenum)+ sum(validity==0);
                                  [b,~]= find(scores<0.5);
                                  tool.pointLog(b,1,ii) = tool.pointLog(a,1,ii-1);
                                  tool.pointLog(b,2,ii)= tool.pointLog(b,2,ii-1);
                                  losses(framenum)= losses(framenum)+sum(scores<0.5);
                                  end
                                  out = insertMarker(frame, points(validity, :), '+', 'Color', 'white');
                                  newI(:,:,framenum) = out(:,:,1);
                                  framenum = framenum + 1;  ii = ii +1;
                                  quality(ii) = sum(validity)/length(validity);
                                  waitbar(ii/dicomFrames)
                                  filename = [sprintf('%03d',po) '.jpg'];
                                  fullname = fullfile(workingDir,'images',filename);
                                  imwrite(insertMarker(J(:,:,framenum), points,'+','Color','red'),fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
                                  po = po+1;
                            end
                            close(h)
                            choice = questdlg('Would you like to generate a video of the points tracked?', ...
                                'Point Video', 'Yes','No','No');
                            %figure
                            %scatter(tool.pointLog(:,1,1),tool.pointLog(:,2,1))
%                             ImageViewer(newI);
                            frames = (1:dicomFrames);
                            quality = quality*100;
                             %figure;
                             %plot(frames, quality)
                             %xlabel('Frames'); ylabel('% of Points Tracked')
                             %title('Tracking Quality')
                switch choice
                    case 'Yes'
                             imageNames = dir(fullfile(workingDir,'images','*.jpg'));
                             imageNames = {imageNames.name}';
                             outputVideo = VideoWriter(fullfile('C:\Users\Leo\Documents\Summer\Summer Research\UltrasoundMeasurement\VUESR_ Kai','points.avi'));
                             outputVideo.FrameRate = 30;
                             open(outputVideo)
                             for ii = 1:length(imageNames)
                                 img = imread(fullfile(workingDir,'images',imageNames{ii}));
                                 writeVideo(outputVideo,img)
                             end
                             close(outputVideo)
                             pointsAvi = VideoReader(fullfile('C:\Users\Leo\Documents\Summer\Summer Research\UltrasoundMeasurement\VUESR_ Kai','points.avi'));
                             ii = 1;
                             while hasFrame(pointsAvi)
                             mov(ii) = im2frame(readFrame(pointsAvi));
                             ii = ii+1;
                             end
                             figure
                             imshow(mov(1).cdata, 'Border', 'tight')
                             movie(mov,1)
                    case 'No'
                        return
                end
                         end
                         
        end
        
        function pixelShearCallback(tool,~,~)
%             if ~isempty(tool.currentROI)                 
%                   if isvalid(tool.currentROI)
%                        if (isempty(tool.pointLog))
%                            msgbox('Please run pixel tracking first to get strain')
%                        else
                          
                            J = uint8(tool.I);
                            dicomFrames = size(J,3);
                            if dicomFrames <= 1
                                warndlg('Need multiple frames to calculate sheer', 'Warning');
                                return;
                            end
                            %Select Points to Track
                            uiwait(msgbox('Select 2 points, 1 on  the vessel edge, and 1 near the vessel center, then hit "Enter"'));
                            figHandle = gcf;
                            [poiX, poiY] = getpts(figHandle);
                            poiX = round(poiX);     poiY = round(poiY);
                            point = [poiX(1), poiY(1)];                       
                            nPoints = size(poiX,1);
                            if  nPoints  ~= 2
                                warndlg('Please select only 2 points. Exiting function', 'Warning');
                                return;
                            end
                            % Create object tracker
                            tracker = vision.PointTracker('MaxBidirectionalError', 3);

                            % Initialize object tracker
                            framenum=1;
                            objectFrame = J(:,:,framenum);
                            initialize(tracker, point(:,:,1), objectFrame);

                            % Show the points getting tracked
                            while framenum <= dicomFrames
                                 %Track the points     
                                  frame =J(:,:,framenum);
                                  [point, validity] = step(tracker, frame);
                                  points(:,:,framenum) = point;
                                  framenum = framenum + 1;
                            end
                            
                            slope = (poiY(2)-poiY(1))/(poiX(2)-poiX(1));
                            %shearPoints = 15;
                            shearPoints = sqrt((poiX(2)-poiX(1))^2 + (poiY(2)-poiY(1))^2);
                            if abs(slope) >= 1.5
                                voffset = 1;
                                if poiY(1) > poiY(2)
                                    voffset = voffset * -1;
                                end
                                hoffset = voffset/slope;
                            else
                                hoffset = 1;
                                if poiX(1) > poiX(2)
                                    hoffset = hoffset * -1;
                                end
                                voffset = hoffset*slope;
                            end
                            for ind = 1:dicomFrames
                                for count = 2:shearPoints
                                    points(count,1,ind) = points(count-1,1,ind)+hoffset;
                                    points(count,2,ind) = points(count-1,2,ind)+voffset;
                                end
                            end
                            
                            pointsTracked = size(points,1);
                            shear = J;
                            % Create new image showing shear rate magnitudes
                            h = waitbar(0,'Calculating wall shear rate...');
                            for indFrame = 1:dicomFrames-1
                                waitbar(indFrame/dicomFrames)
                                for ind = 1:pointsTracked
%                                     IX = J(:,:,indFrame);                  %Frame 1
%                                     IY = J(:,:,indFrame+1);              %Frame 2
                                    %FILT = ones(5);                           %Filter matrix
                                    KRNL_LMT = [3 3];                   %Group of pixels you're trying to find in next image
                                    SRCH_LMT = [2 2];                   %Region
                                    POS = round(points(ind,:,indFrame));  %Origin of krnl and srch
                                    %[RHO]=corr2D(IX,IY,FILT,KRNL_LMT,SRCH_LMT,POS);
                                    [RHO]=newcorr2D(KRNL_LMT,SRCH_LMT,POS,indFrame);
                                    rho(ind,indFrame) = mean(mean(RHO));
                                    rho2(ind,indFrame) = max(max(RHO));
                                    %shear(POS(2)-2:POS(2)+2,POS(1)-2:POS(1)+2,indFrame) = max(max(RHO));
                                end
                            end
                            close(h);
                            %Normalize rho values to display as image
                            %intensities
                            for indFrame = 1:dicomFrames-1
                                for ind = 1:pointsTracked
                                    POS = round(points(ind,:,indFrame)); 
                                    maxrho = max(max(rho));
                                    newrho = 200.*rho./maxrho;
                                    shear(POS(2),POS(1),indFrame) = newrho(ind,indFrame);
                                end
                            end
%                             ImageViewer(shear);
                                   
                            for ind = 1:shearPoints
                                wallshear(ind) = mean(rho(ind,:));
                                wallshear2(ind) = mean(rho2(ind,:));
                            end
                            
                           
                            time = 1:dicomFrames-1;
                            z1 = rho(1,:); z2=rho(2,:); z3=rho(3,:); z4 = rho(4,:); z5 = rho(5,:);
                            figure;
                            d1(1:dicomFrames-1) = 1; d2(1:dicomFrames-1) = 2; d3(1:dicomFrames-1) = 3;
                            d4(1:dicomFrames-1) = 4; d5(1:dicomFrames-1) = 5;
                            plot3(d1,time,z1,d2,time,z2,d3,time,z3,d4,time,z4,d5,time,z5);
                            %surf(rho);
                            xlabel('Distance from Wall'); ylabel('Time'); zlabel('Wall Shear (corr coeff)');
                            figure;
                            dist = 1:shearPoints;
                            plot(dist, wallshear2)
                            xlabel('Distance from Wall'); ylabel('Wall Shear (max corr coeff)'); 
                            title('Wall Shear')
                            figure;
                            imagesc(rho);
                            hold on;
                            colorbar;
                            xlabel('Time');ylabel('Distance from wall');title('Average corr coeff mapping')
                            save shear.mat
                            function correlation = newcorr2D(krnl,srch,pos,frame)
%
                                correlation = zeros(2*srch(1)+1,2*srch(2)+1);
                                for i = 1:(2*srch(1)+1)
                                    for j = 1:(2*srch(2)+1)
                                        x = -srch(1)+i-1;
                                        y = -srch(2)+j-1;
                                        kernel = J(pos(1)-krnl(1):pos(1)+krnl(1),pos(2)-krnl(2):pos(2)+krnl(2),frame);
                                        search = J(pos(1)-krnl(1)+x:pos(1)+krnl(1)+x,pos(2)-krnl(2)+y:pos(2)+krnl(2)+y,frame+1);
                                        correlation(i,j) = corr2(kernel,search);
                                    end
                                end
                            end
%                        end   
%                   else
%                        msgbox('Please select a region of interest');
%                        return;
%                   end
%             else
%                   msgbox('Please select a region of interest');
%                   return;  
%             end

        end
        
        function pixelWallCallback(tool,~,~)
            %Wall Strain
            J = uint8(tool.I);
           dicomFrames = size(tool.I,3);

            %Select Points to Track
            uiwait(msgbox('Select 2 points, 1 on  the vessel edge, and 1 near the vessel center, then hit "Enter"'));
            figHandle = gcf;
            [poiX, poiY] = getpts(figHandle);
            poiX = round(poiX);     poiY = round(poiY);
            point = [poiX(1), poiY(1)];
            
            slope = (poiY(2)-poiY(1))/(poiX(2)-poiX(1));
            if abs(slope) >= 1.5
                voffset = 4;
                if poiY(1) > poiY(2)
                    voffset = voffset * -1;
                end
                hoffset = voffset/slope;
            else
                hoffset = 4;
                if poiX(1) > poiX(2)
                    hoffset = hoffset * -1;
                end
                voffset = hoffset*slope;
            end
            point(2,1) = point(1,1)+hoffset;
            point(2,2) = point(1,2)+voffset;
            point(3,1) = point(1,1)-hoffset;
            point(3,2) = point(1,2)-voffset;
            
            points(:,:,1) = point;
            distX(1) = points(2,1,1)-points(3,1,1);
            distY(1) = points(2,2,1)-points(3,2,1);
            dist(1) = sqrt(distX(1)^2+distY(1)^2);

            % Create object tracker
            tracker = vision.PointTracker('MaxBidirectionalError', 1);
            
            % Initialize object tracker
            framenum=1;
            objectFrame = J(:,:,framenum);
            initialize(tracker, point(:,:,1), objectFrame);
            pointImage = insertMarker(objectFrame, points, '+', 'Color', 'white');
            newI(:,:,1) = pointImage(:,:,1);

            % Show the points getting tracked
            while framenum < dicomFrames
                 %Track the points     
                  frame =J(:,:,framenum);
                  [point, validity] = step(tracker, frame);
                  framenum = framenum + 1;
                  points(:,:,framenum) = point;
                  out = insertMarker(frame, point(validity, :), '+', 'Color', 'white');
                  newI(:,:,framenum) = out(:,:,1);
                  distX(framenum) = points(2,1,framenum)-points(3,1,framenum);
                  distY(framenum) = points(2,2,framenum)-points(3,2,framenum);
                  dist(framenum) = sqrt(distX(framenum)^2+distY(framenum)^2);
            end

            avgdist = mean(dist);
            for ind = 1:dicomFrames
                strain(ind) = (dist(ind)-avgdist)/avgdist;
            end
%             ImageViewer(newI);
            frame = 1:dicomFrames;
            figure;
            plot(frame, strain)
            xlabel('Distance from wall'); ylabel('Wall Strain')
            title('Wall Strain')            
        end
        
        function pixelEdgeCallback(tool,~,~)
            if ~isempty(tool.currentROI)                 
                  if isvalid(tool.currentROI)
                       %imageROI = tool.currentROI;
                       disp(tool.currentROI);
                      % GRAYTHRESH EDGE DETECT
                        indFrame = 1;
                        imageROI = uint8(tool.I);
                        nFrames = size(tool.I,3);
                        while indFrame <= nFrames
                            imageROI_adjusted(:,:,indFrame) = imadjust(imageROI(:,:,indFrame));
                            imageROI_level(indFrame) = graythresh(imageROI_adjusted(:,:,indFrame));
                            imageROI_BW(:,:,indFrame) = im2bw(imageROI_adjusted(:,:,indFrame),...
                                imageROI_level(indFrame)*.3);
                            indFrame = indFrame + 1;
                        end
                        imageROI_BW = uint8(imageROI_BW);
%                        ImageViewer(imageROI_BW)
                        % An interesting result
                        time = (1:nFrames)./16;
                        figure;
                        plot(time,imageROI_level)
                        xlabel('Time [s]')
                        ylabel('Graythreshold')
                        title('A heartbeat measure by image contrast')
                        else
                            msgbox('Please select a region of interest');
                            return;
                   end
             else
                  msgbox('Please select a region of interest');
                  return;
             end

        end
        
        function filterImageCallback(tool,~,~)
            for i=1:size(tool.I,3)
                    a(:,:,i) = medfilt2(tool.I(:,:,i));
            end
            tool.I = a;
        end
        
    end
 end

 %Leos code


function newLowerRangePosition(tool,~,~,hObject)
cp = get(hObject,'CurrentPoint'); cp=[cp(1,1) cp(1,2)];
range=get(tool.handles.Axes,'Clim');
Xlims=get(hObject,'Xlim');
range(1)=cp(1);
W=diff(range);
L=mean(range);
if W>0 && range(1)>=Xlims(1)
    setWL(tool,W,L)
end
end

function newUpperRangePosition(tool,~,~,hObject)
cp = get(hObject,'CurrentPoint'); cp=[cp(1,1) cp(1,2)];
range=get(tool.handles.Axes,'Clim');
Xlims=get(hObject,'Xlim');
range(2)=cp(1);
W=diff(range);
L=mean(range);
if W>0 && range(2)<=Xlims(2)
    setWL(tool,W,L)
end
end

function newLevelRangePosition(tool,~,~,hObject)
cp = get(hObject,'CurrentPoint'); cp=[cp(1,1) cp(1,2)];
range=get(tool.handles.Axes,'Clim');
Xlims=get(hObject,'Xlim');
L=cp(1);
W=diff(range);
if L>=Xlims(1) && L<=Xlims(2)
    setWL(tool,W,L)
end
end

function newROIposition(tool,~,hObject)
handlesROI=tool.handlesROI;
for i=1:length(handlesROI)
    if isvalid(handlesROI{i})
        setColor(handlesROI{i},'b');
    end
end
setColor(hObject,'r');
mask = createMask(hObject);
im=get(tool.handles.I,'CData');
m=mean(im(mask));
noise=std(im(mask));
set(tool.handles.ROIinfo,'String',['STD:' num2str(noise,'%+.4f') '   Mean:' num2str(m,'%+.4f')])
tool.currentROI=hObject;
end

function adjustContrastMouse(tool,~,~,dp,hObject,W,L)
cp = get(hObject,'CurrentPoint'); cp=[cp(1,1) cp(1,2)];
d=round(cp-dp);
W2=W+d(1); L=L-d(2);
if W2>=1
    W=W2;
end
setWL(tool,W,L)
end

function adjustZoomMouse(tool,~,~,dp,hObject)
cp = get(hObject,'CurrentPoint'); cp=[cp(1,1) cp(1,2)];
d=cp(2)-dp(2);
zFactor=.025;
if d>0
    zoom(1+zFactor)
elseif d<0
    zoom(1-zFactor)
end
fun=@(Newsrc,Newevnt) adjustZoomMouse(tool,Newsrc,Newevnt,cp,tool.handles.Axes);
set(tool.handles.fig,'WindowButtonMotionFcn',fun)
axis fill

end

function adjustPanMouse(~,~,~,map,hObject,xlims,ylims)
cp = get(hObject,'CurrentPoint'); cp=[cp(1,1) cp(1,2)];
d=(map-cp)/1.25;
set(hObject,'Xlim',xlims+d(1),'Ylim',ylims+d(2))
end

function buttonUpFunction(tool,src,~)

fun=@(src,evnt)getImageInfo(tool,src,evnt);
set(src,'WindowButtonMotionFcn',fun);

end

function getImageInfo(tool,~,~)
pos=round(get(tool.handles.Axes,'CurrentPoint'));
pos=pos(1,1:2);
Xlim=get(tool.handles.Axes,'Xlim');
Ylim=get(tool.handles.Axes,'Ylim');
n=round(get(tool.handles.Slider,'value'));
if n == 0
    n = 1;
end
if pos(1)>0 && pos(1)<=size(tool.I,2) && pos(1)>=Xlim(1) && pos(1) <=Xlim(2) && pos(2)>0 && pos(2)<=size(tool.I,1) && pos(2)>=Ylim(1) && pos(2) <=Ylim(2)
    set(tool.handles.Info,'String',['(' num2str(pos(1)) ',' num2str(pos(2)) ') ' num2str(tool.I(pos(2),pos(1),n))])
else
    set(tool.handles.Info,'String','(x,y) val')
end



end

function panelResizeFunction(tool,~,~,w,h,~)
units=get(tool.handles.Panels.Large,'Units');
set(tool.handles.Panels.Large,'Units','Pixels')
pos=get(tool.handles.Panels.Large,'Position');
set(tool.handles.Panels.Large,'Units',units)
if get(tool.handles.Tools.Hist,'value')
    set(tool.handles.Panels.Image,'Position',[w w pos(3)-3.8*w pos(4)-2*w-h])
else
    set(tool.handles.Panels.Image,'Position',[w w pos(3)-3.8*w pos(4)-2*w])
end
%set(tool.handles.Panels.Image,'Position',[w w pos(3)-2*w pos(4)-2*w])
set(tool.handles.Panels.Hist,'Position',[w pos(4)-w-h pos(3)-2*w h])
set(tool.handles.Panels.Tools,'Position',[0 pos(4)-w pos(3) w])
set(tool.handles.Panels.ROItools,'Position',[pos(3)-2.8*w  w 2.8*w pos(4)-2*w])
set(tool.handles.Panels.Slider,'Position',[0 w w pos(4)-2*w])
set(tool.handles.Panels.Info,'Position',[0 0 pos(3) w])
axis(tool.handles.Axes,'fill');
%buff=(w-wbutt)/2;
%pos=get(tool.handles.Panels.ROItools,'Position');


end

function icon = makeToolbarIconFromPNG(filename)
% makeToolbarIconFromPNG  Creates an icon with transparent
%   background from a PNG image.

%   Copyright 2004 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2004/08/10 01:50:31 $

% Read image and alpha channel if there is one.
[icon,map,alpha] = imread(filename);

% If there's an alpha channel, the transparent values are 0.  For an RGB
% image the transparent pixels are [0, 0, 0].  Otherwise the background is
% cyan for indexed images.
if (ndims(icon) == 3) % RGB
    
    idx = 0;
    if ~isempty(alpha)
        mask = alpha == idx;
    else
        mask = icon==idx;
    end
    
else % indexed
    
    % Look through the colormap for the background color.
    for i=1:size(map,1)
        if all(map(i,:) == [0 1 1])
            idx = i;
            break;
        end
    end
    
    mask = icon==(idx-1); % Zero based.
    icon = ind2rgb(icon,map);
    
end

% Apply the mask.
icon = im2double(icon);

for p = 1:3
    
    tmp = icon(:,:,p);
    if ndims(mask)==3
        tmp(mask(:,:,p))=NaN;
    else
        tmp(mask) = NaN;
    end
    icon(:,:,p) = tmp;
    
end

end

function saveImage(tool,~,~)
cmap = colormap(tool.handles.Axes);
switch get(tool.handles.Tools.SaveOptions,'value')
    case 1 %Save just the current slice
        I=get(tool.handles.I,'CData'); lims=get(tool.handles.Axes,'CLim');
        I=gray2ind(mat2gray(I,lims),256);
        [FileName,PathName] = uiputfile({'*.png';'*.tif';'*.jpg';'*.bmp';'*.gif';'*.hdf'; ...
            '*.jp2';'*.pbm';'*.pcx';'*.pgm'; ...
            '*.pnm';'*.ppm';'*.ras';'*.xwd'},'Save Image');
        
        if FileName == 0
        else
            imwrite(I,cmap,[PathName FileName])
        end
    case 2 %Save entire dicom stack
        lims=get(tool.handles.Axes,'CLim');
        [FileName,PathName] = uiputfile({'*.dcm';'*.tif'},'Save Image Stack');
        [~, ~, ext] = fileparts(FileName);

        if FileName == 0
        else
            if strcmp(ext,'.tif') || strcmp(ext,'.TIF')
                for i=1:size(tool.I,3)
                    imwrite(gray2ind(mat2gray(tool.I(:,:,i),lims),256),cmap, [PathName FileName], 'WriteMode', 'append',  'Compression','none');
                end
            elseif strcmp(ext,'.DCM') || strcmp(ext,'.dcm')
                   dicomwrite(permute(uint8(tool.I), [1 2 4 3]), FileName);
            else
            end
        end
end
end

function exportCallback(tool,~,~)
     
    f = figure('Visible','off','Position',[360,500,320,280],'NumberTitle', 'off','ToolBar','none','MenuBar','none','Name','Edit/Export Data');
    
    % Create headers and text boxes
    patientIDText = uicontrol('Style', 'edit','FontSize', 10,...
           'String', num2str(tool.zpatientID), 'Position', [140 240 120 20]);
   dicomFileText = uicontrol('Style', 'text','FontSize', 10,'BackgroundColor', [.8,.8,.8],...
           'String', tool.fName, 'Position', [140 220 120 20]);
    locationText = uicontrol('Style', 'edit','FontSize', 7,...
           'String', cellstr(tool.zlocation), 'Position', [140 200 120 20]);
    diameterText = uicontrol('Style', 'edit','FontSize', 10,...
           'String', num2str(tool.zdiameter), 'Position', [140 180 120 20]);
    distensibilityText = uicontrol('Style', 'edit','FontSize', 10,...
           'String', num2str(tool.zdistensibility), 'Position', [140 160 120 20]);
    elasticityText = uicontrol('Style', 'edit','FontSize', 10,...
           'String', num2str(tool.zelasticity), 'Position', [140 120 120 20]);
    mapText = uicontrol('Style', 'edit','FontSize', 10,...
           'String', num2str(tool.map), 'Position', [140  80 120 20]);
    hrText = uicontrol('Style', 'edit','FontSize', 10,...
           'String', num2str(tool.hRate), 'Position', [140 60 120 20]);

     patientIDHeader = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'Patient ID:', 'Position', [10 240 128 20], 'BackgroundColor', [.8,.8,.8]);
     dicomFileHeader = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'DICOM File:', 'Position', [10 220 128 20], 'BackgroundColor', [.8,.8,.8]);
     locationHeader = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'Blood Vessel Type:', 'Position', [10 200 128 20], 'BackgroundColor', [.8,.8,.8]);
      diameterHeader = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'Blood Vessel Diam:', 'Position', [10 180 128 20], 'BackgroundColor', [.8,.8,.8]);
      distensibilityHeader = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'Distensibility:', 'Position', [10 160 128 20], 'BackgroundColor', [.8,.8,.8]);
      elasticityHeader = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'Elasticity:', 'Position', [10 120 128 20], 'BackgroundColor', [.8,.8,.8]);
      mapHeader = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'Mean Art. Pressure:', 'Position', [10 80 128 20], 'BackgroundColor', [.8,.8,.8]);
      hrHeader = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'Heart Rate:', 'Position', [10 60 128 20], 'BackgroundColor', [.8,.8,.8]);
       
      diameterUnits = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'cm', 'Position', [262 180 30 20], 'BackgroundColor', [.8,.8,.8]);
      distensibilityUnits = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'cm', 'Position', [262 160 30 20], 'BackgroundColor', [.8,.8,.8]);
      elasticityUnits = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', '', 'Position', [262 120 30 20], 'BackgroundColor', [.8,.8,.8]);
      mapUnits = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'mmHg', 'Position', [262 80 40 20], 'BackgroundColor', [.8,.8,.8]);
      hrUnits = uicontrol('Style', 'text','FontSize', 10,'HorizontalAlignment','left',...
           'String', 'bpm', 'Position', [262 60 30 20], 'BackgroundColor', [.8,.8,.8]);
       
       
     % Create push buttons
    Export = uicontrol('Style', 'pushbutton', 'String', 'Export',...
        'FontSize', 11,'Position', [10 10 65 25],'Callback', @(hObject,evnt) Export_callback(tool,hObject,evnt));
     Ok = uicontrol('Style', 'pushbutton', 'String', 'Ok',...
        'FontSize', 11,'Position', [160 10 70 25],'Callback', @(hObject,evnt) OK_callback(tool,hObject,evnt));
    Cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel',...
        'FontSize', 11,'Position', [240 10 80 25], 'Callback', @(hObject,evnt) Cancel_callback(tool,hObject,evnt));
    
    %Create check boxes
    pulsgraph = uicontrol('Style', 'checkbox', 'FontSize', 9,'HorizontalAlignment','left',...
        'String', 'Include pulsatility graph', 'Position', [100 141 200 18], 'BackgroundColor', [.8,.8,.8]);
    elastgraph = uicontrol('Style', 'checkbox', 'FontSize', 9,'HorizontalAlignment','left',...
        'String', 'Include elasticity graph', 'Position', [100 101 200 18], 'BackgroundColor', [.8,.8,.8]);
    
    % Assign the GUI a name to appear in the window title. 
    set(f,'Name','Edit or Export Data')
     set([f, patientIDText, dicomFileText, locationText, diameterText, distensibilityText, elasticityText, mapText, hrText,...
              patientIDHeader, dicomFileHeader, locationHeader, diameterHeader, distensibilityHeader, elasticityHeader, mapHeader, hrHeader,...
              diameterUnits, distensibilityUnits, elasticityUnits, mapUnits, hrUnits, Export, Ok, Cancel, pulsgraph,elastgraph],'Units','normalized');
    
    % Move the GUI to the center of the screen. 
    movegui(f,'center')
    % Make the GUI visible. 
    set(f,'Visible','on','Resize','off');
    
    function Cancel_callback(~,~,~)
         %Closes window without saving changes
          close;
    end
    function OK_callback(tool,~,~)
          
          choice = questdlg('Are you sure you want to overwrite data?', ...
                'Warning', ...
                'Yes','No','No');
            % Handle response
            switch choice
                case 'Yes'
                      %Saves changes, then closes window
                      tool.zpatientID = str2double(get(patientIDText,'String'));
                      tool.fName = get(dicomFileText,'String');
                      tool.zlocation =get(locationText, 'String');
                      tool.zdiameter =str2num(get(diameterText, 'String')); %#ok<ST2NM>
                      tool.zdistensibility =str2num(get(distensibilityText, 'String')); %#ok<ST2NM>
                      tool.map = str2num(get(mapText, 'String')); %#ok<ST2NM>
                      tool.zelasticity =str2num(get(elasticityText,'String')); %#ok<ST2NM>
                      tool.hRate = str2num(get(hrText, 'String')); %#ok<ST2NM>
                      close;
                case 'No'
                    %Don't save changes
            end
    end

    function Export_callback(tool,~,~)
          disp('Export')
          %Exports data to excel
          patientID = cellstr(get(patientIDText, 'String'));
          dicomfile = cellstr(get(dicomFileText, 'String'));
          location =  cellstr(get(locationText, 'String'));
          diameter =  cellstr(get(diameterText, 'String'));
          distensibility =  cellstr(get(distensibilityText, 'String'));
          MAP =  cellstr(get(mapText, 'String'));
          elasticity =  cellstr(get(elasticityText, 'String'));
          heartrate =  cellstr(get(hrText, 'String'));
          time = transpose((1:size(tool.I,3))/tool.fRate);
          pulsatility = transpose(tool.pointDistCm);
          disp(size(time)); disp(size(pulsatility));
          T1 = table(patientID, dicomfile, location, diameter, distensibility, MAP, elasticity, heartrate);     

          %Save data to .csv file
          fileName = strcat('patient', patientID, '_', location, '.xlsx');
          [filename,pathname] = uiputfile(fileName,'Save data file');                      
          if isequal(filename,0) || isequal(pathname,0)
             disp('User selected Cancel')
          else
             disp(['User selected ',fullfile(pathname,filename)])
             writetable(T1,fullfile(pathname,filename),'Sheet',1);
             if length(pulsatility) == length(time) 
                 if get(pulsgraph,'Value') == 1 && get(elastgraph,'Value') == 1
                       T2 = table(time, pulsatility);
                       writetable(T2,fullfile(pathname,filename),'Sheet',2);
                 elseif get(pulsgraph,'Value') == 1
                      T2 = table(time, pulsatility);
                      writetable(T2,fullfile(pathname,filename), 'Sheet',2);
                 elseif get(elastgraph,'Value') == 1
                 end
             else
                 msgbox('Variables are not the same size')
             end
             msgbox('Data export successful')      
          end
          
    end
        
end

function ShowHistogram(tool,~,~,w,h)
set(tool.handles.Panels.Large,'Units','Pixels')
pos=get(tool.handles.Panels.Large,'Position');
set(tool.handles.Panels.Large,'Units','normalized')

if get(tool.handles.Tools.Hist,'value')
    set(tool.handles.Panels.Image,'Position',[w w pos(3)-2*w pos(4)-2*w-h])
else
    set(tool.handles.Panels.Image,'Position',[w w pos(3)-2*w pos(4)-2*w])
end
axis(tool.handles.Axes,'fill');
showSlice(tool)

end

function initializeVUESR(tool)
          
           % Initialize Variables
            tool.handlesROI = [];
            tool.currentROI = [];
            %tool.calibration = 1;
            tool.pixelDensity = 10;
            tool.accFrames = [1 size(tool.I,3)-1];
            tool.map = 50;
            tool.hRate = [];
            %tool.fRate = [];
            tool.I = double(tool.I);
            %tool.zlocation = []; 
            tool.zdiameter = [];
            tool.zdistensibility = []; 
            tool.zelasticity =[]; 
            tool.pointDistCm = [];
            
            %Do optical character recognition to find measurement location
            bottom = size(tool.I,1);
            BW = im2bw(uint8(tool.I(floor(bottom/1.5):bottom,:,1)), 0.5);  %WDR fixed
            
            str = '';
            has_ocr = ismember(exist('ocr'), [2 3 4 5 6 8]);               %WDR fixed
            if has_ocr
                try
                  ocrResults = ocr(BW);
                  str = ocrResults.Text;
                  str = regexprep(str,'[^a-zA-Z]','');
                catch
                    has_ocr = false;
                end
            end
            tool.zlocation =str;
            
            %Do optical character recognition to find the frame rate in hz
            right = size(tool.I,2);
            BW = im2bw(uint8(tool.I(:,round(right/1.5):right,1)), 0.5);
            if has_ocr                                                  %WDR fixed
              ocrResults = ocr(BW);
              results = char(ocrResults.Text);
               k = strfind(results,'Hz');
               if length(k) > 1
                  k = k(2);
               end
            else
                k = [];
            end
            if ~isempty(k)
                 hz = results(k-3:k-2);
                 tool.fRate = str2double(hz);
            else
                disp('No frame rate detected, asuming 1 Hz')
                tool.fRate = 1;
            end
            
            %Use dicominfo to read important metadata from file
            x = dicominfo(tool.fName);
            
            if  isfield(x, 'SequenceOfUltrasoundRegions') && ...             %WDR fixed
                isfield(x.SequenceOfUltrasoundRegions, 'Item_1') && ...
                isfield(x.SequenceOfUltrasoundRegions.Item_1, 'PhysicalDeltaX')
              tool.calibration = x.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX;
            elseif isfield(x, 'PixelSpacing')
                tool.calibration = x.PixelSpacing(1);
            else
                disp('No pixel spacing detected, assuming 0.3')
                tool.calibration = 0.3;
            end
            tool.studydate = x.StudyDate;
           
           %x.MechanicalIndex
            %tool.fRate = ocrResults.Text;
end

function createButtons(tool,lp,buff,widthSidePanel)
% BUTTONS
           [iptdir, MATLABdir] = ipticondir;
            % Create save button
            tool.handles.Tools.Save = ...
                uicontrol(tool.handles.Panels.Tools,...
                'Style','pushbutton',...
                'String','',...
                'Position', [lp+1.1*widthSidePanel, buff, widthSidePanel, widthSidePanel]);
            icon_save = makeToolbarIconFromPNG([MATLABdir '/file_save.png']);
            set(tool.handles.Tools.Save,'CData',icon_save);
            lp = lp+2*widthSidePanel;
            
            tool.handles.Tools.SaveOptions = ...
                uicontrol(tool.handles.Panels.Tools,...
                'Style','popupmenu',...
                'String',{'as single slice','as entire stack'},...
                'Position',[lp+buff, buff, 5.5*widthSidePanel, widthSidePanel]);
            fun = @(hObject,evnt) saveImage(tool,hObject,evnt);
            set(tool.handles.Tools.Save,'Callback',fun)
            set(tool.handles.Tools.Save,'TooltipString','Save image as slice or entire stack')
            lp = lp+5.7*widthSidePanel;
            
             % Export button
            tool.handles.Tools.Export = ...
                uicontrol(tool.handles.Panels.Tools,...
                'Style','pushbutton',...
                'String','Edit/Export Data',...
                'Position', [lp+buff, buff, 5*widthSidePanel, widthSidePanel],...
                'TooltipString','Export Data To Excel');
            fun = @(hObject,evnt) exportCallback(tool,hObject,evnt);
            set(tool.handles.Tools.Export,'Callback',fun)
            lp = lp+5.2*widthSidePanel;
            
            % Open New File
            tool.handles.Tools.OpenImage = ...
                uicontrol(tool.handles.Panels.Tools,...
                'Style','pushbutton',...
                'String','New Image',...
                'Position',[lp+buff, buff, 4*widthSidePanel, widthSidePanel],...
                'TooltipString','Open New Image to Analyze');
            fun=@(hObject,evnt) openImageCallback(tool,hObject,evnt);
            set(tool.handles.Tools.OpenImage ,'Callback',fun)

            % Create Circle ROI button
            tool.handles.Tools.CircleROI = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','',...
                'Position',[buff+2.25*widthSidePanel, buff+2*widthSidePanel, widthSidePanel, widthSidePanel],...
                'TooltipString','Create Elliptical Region of Interest');
            icon_ellipse = makeToolbarIconFromPNG([MATLABdir '/tool_shape_ellipse.png']);
            set(tool.handles.Tools.CircleROI,'Cdata',icon_ellipse)
            fun = @(hObject,evnt) measureImageCallback(tool,hObject,evnt,'ellipse');
            set(tool.handles.Tools.CircleROI,'Callback',fun)
            
            % Create Square ROI button
            tool.handles.Tools.SquareROI = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','',...
                'Position',[buff+.25*widthSidePanel, buff+2*widthSidePanel, widthSidePanel, widthSidePanel],...
                'TooltipString','Create Rectangular Region of Interest');
            icon_rect = makeToolbarIconFromPNG([MATLABdir '/tool_shape_rectangle.png']);
            set(tool.handles.Tools.SquareROI,'Cdata',icon_rect)
            fun = @(hObject,evnt) measureImageCallback(tool,hObject,evnt,'rectangle');
            set(tool.handles.Tools.SquareROI,'Callback',fun)
            
            % Create Polygon ROI button
            tool.handles.Tools.PolyROI = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','\_/',...
                'Position',[buff+1.25*widthSidePanel, buff+2*widthSidePanel, widthSidePanel, widthSidePanel],...
                'TooltipString','Create Polygon Region of Interest');
            fun = @(hObject,evnt) measureImageCallback(tool,hObject,evnt,'polygon');
            set(tool.handles.Tools.PolyROI,'Callback',fun)
            
            % Create Delete Button
            tool.handles.Tools.DeleteROI = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Delete ROI',...
                'Position',[buff, buff+widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Delete Region of Interest');
            fun = @(hObject,evnt) deletecurrentROI(tool,hObject,evnt);
            set(tool.handles.Tools.DeleteROI,'Callback',fun)
            
            % Create Export ROI Button
            tool.handles.Tools.ExportROI = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Export ROI',...
                'Position',[buff, buff, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Export Region of Interest to Workspace',...
                'ForegroundColor','k');
            fun = @(hObject,evnt) exportROI(tool,hObject,evnt);
            set(tool.handles.Tools.ExportROI,'Callback',fun)
            
            % Create Ruler button
            tool.handles.Tools.Ruler = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','',...
                'Position',[buff+2.25*widthSidePanel, buff+6.25*widthSidePanel, widthSidePanel, widthSidePanel],...
                'TooltipString','Measure Distance');
            icon_distance = makeToolbarIconFromPNG([MATLABdir '/tool_line.png']);
            set(tool.handles.Tools.Ruler,'CData',icon_distance);
            fun = @(hObject,evnt) measureImageCallback(tool,hObject,evnt,'ruler');
            set(tool.handles.Tools.Ruler,'Callback',fun)
            
            % Create Line Profile button
            tool.handles.Tools.Profile = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','',...
                'Position',[buff+1.25*widthSidePanel, buff+6.25*widthSidePanel, widthSidePanel, widthSidePanel],...
                'TooltipString','Get Line Profile');
            icon_profile = makeToolbarIconFromPNG([iptdir '/profile.png']);
            set(tool.handles.Tools.Profile,'Cdata',icon_profile)
            fun = @(hObject,evnt) measureImageCallback(tool,hObject,evnt,'profile');
            set(tool.handles.Tools.Profile,'Callback',fun)
            
            % Create Crop tool button
            tool.handles.Tools.Crop = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','',...
                'Position',[buff+.25*widthSidePanel, buff+6.25*widthSidePanel, widthSidePanel, widthSidePanel],...
                'TooltipString','Crop Image');
            icon_profile = makeToolbarIconFromPNG([iptdir '/crop_tool.png']);
            set(tool.handles.Tools.Crop ,'Cdata',icon_profile)
            fun=@(hObject,evnt) CropImageCallback(tool,hObject,evnt);
            set(tool.handles.Tools.Crop ,'Callback',fun)
            
            % Create Crop # frames button
            tool.handles.Tools.Cropframes = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','# frames',...
                'Position',[buff+.25*widthSidePanel, buff + 4.25*widthSidePanel, 3*widthSidePanel, widthSidePanel],...
                'TooltipString','Crop # of frames in the cine loop',...
                'ForegroundColor','k');
            fun=@(hObject,evnt) CropFramesCallback(tool,hObject,evnt);
            set(tool.handles.Tools.Cropframes ,'Callback',fun)
            
            % Create Auto Crop button
            tool.handles.Tools.AutoCrop = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Auto Crop',...
                'Position',[buff+.25*widthSidePanel, buff + 5.25*widthSidePanel, 3*widthSidePanel, widthSidePanel],...
                'TooltipString','Automatically crop image',...
                'ForegroundColor','k');
            fun=@(hObject,evnt) AutoCropCallback(tool,hObject,evnt);
            set(tool.handles.Tools.AutoCrop ,'Callback',fun)
            
            
            % ********************************NEW STUFF**************************

            %Filter image
            tool.handles.Tools.FilterImage = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Filter',...
                'Position',[buff, buff+20*widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Use median filter to smooth image');
            fun=@(hObject,evnt) filterImageCallback(tool,hObject,evnt);
            set(tool.handles.Tools.FilterImage ,'Callback',fun)
            
            % Create 2 pixel tracking button
            tool.handles.Tools.Track2 = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Distensibility',...   %'String','Pulsatility',...
                'Position',[buff, buff+12*widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Track 2 Points in the Image ');
            fun=@(hObject,evnt) pixelTrack2Callback(tool,hObject,evnt);
            set(tool.handles.Tools.Track2 ,'Callback',fun)
            
            %Create Strain button
            tool.handles.Tools.Strain = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Strain',...
                'Position',[buff, buff+11*widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Find Strain for Entire Region');
            fun=@(hObject,evnt) pixelStrainCallback(tool,hObject,evnt);
            set(tool.handles.Tools.Strain ,'Callback',fun)
            
            %Create Wall Shear button
            tool.handles.Tools.Shear = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Shear',...
                'Position',[buff, buff+10*widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Find Wall Shear Rate in Region of Interest ');
            fun=@(hObject,evnt) pixelShearCallback(tool,hObject,evnt);
            set(tool.handles.Tools.Shear ,'Callback',fun)

            %Create Wall Strain button
            tool.handles.Tools.Wall = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Wall Strain',...
                'Position',[buff, buff+9*widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Find strain near the vessel wall');
            fun=@(hObject,evnt) pixelWallCallback(tool,hObject,evnt);
            set(tool.handles.Tools.Wall ,'Callback',fun)
             
            %Create Calibration Button
            tool.handles.Tools.Calibrate = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Calibrate',...
                'Position',[buff, buff+16*widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Calibrate image pixels to cm ');
            fun = @(hObject,evnt) pixelCalibrateCallback(tool,hObject,evnt);
            set(tool.handles.Tools.Calibrate,'Callback',fun)
            
            %Create Settings Button
            tool.handles.Tools.Settings = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Settings',...
                'Position',[buff, buff+17*widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','File Settings');
            fun = @(hObject,evnt) pixelSettingsCallback(tool,hObject,evnt);
            set(tool.handles.Tools.Settings,'Callback',fun)
            
            % Create Help Button
            tool.handles.Tools.Help = ...
                uicontrol(tool.handles.Panels.ROItools,...
                'Style','pushbutton',...
                'String','Help',...
                'Position',[buff, buff+15*widthSidePanel, 3.5*widthSidePanel, widthSidePanel],...
                'TooltipString','Help with VUESR');
            fun = @(hObject,evnt) displayHelp(tool,hObject,evnt);
            set(tool.handles.Tools.Help,'Callback',fun)
                      
             %Create text boxes for user guidance
             tool.handles.Tools.SetUp = uicontrol(tool.handles.Panels.ROItools,'Style','text',...
                'String','Set Up','BackgroundColor','k','ForegroundColor','w','FontWeight','bold',...
                'Position',[buff, buff+18*widthSidePanel, 3.5*widthSidePanel, widthSidePanel]);
            tool.handles.Tools.ROI = uicontrol(tool.handles.Panels.ROItools,'Style','text',...
                'String','ROI Tools','BackgroundColor','k','ForegroundColor','w','FontWeight','bold',...
                'Position',[buff, buff+3*widthSidePanel, 3.5*widthSidePanel, 1*widthSidePanel]);
            tool.handles.Tools.Image = uicontrol(tool.handles.Panels.ROItools,'Style','text',...
                'String','Image Tools','BackgroundColor','k','ForegroundColor','w','FontWeight','bold',...
                'Position',[buff, buff+7.25*widthSidePanel, 3.5*widthSidePanel, widthSidePanel]);        
            tool.handles.Tools.Analyze = uicontrol(tool.handles.Panels.ROItools,'Style','text',...
                'String','Analysis','BackgroundColor','k','ForegroundColor','w','FontWeight','bold',...
                'Position',[buff, buff+13*widthSidePanel, 3.5*widthSidePanel, widthSidePanel]);        
%             tool.handles.Tools.ImageTracking = uicontrol(tool.handles.Panels.ROItools,'Style','text',...
%                 'String', 'Image Tracking','BackgroundColor','k','ForegroundColor','w','FontWeight','bold',...
%                 'Position',[buff, buff+15*widthSidePanel, 3.5*widthSidePanel, 2*widthSidePanel]);
            
            % Set font size of all the tool objects
            
            tools_cell = struct2cell(tool.handles.Tools);                    %WDR fix
            set(vertcat(tools_cell{:}), ...
                'FontSize',9,...
                'Units','Pixels')
            
end
function geom = means( x, y ) 
%
% H.J. Sommer III, Ph.D., Professor of Mechanical Engineering, 337 Leonhard Bldg
% The Pennsylvania State University, University Park, PA  16802
% (814)863-8997  FAX (814)865-9693  hjs1-at-psu.edu  www.mne.psu.edu/sommer/
 
% begin function POLYGEOM
 
% check if inputs are same size
if ~isequal( size(x), size(y) )
  error( 'X and Y must be the same size');
end

%means of corner points are taken
xm = mean(x);
ym = mean(y);
geom= [xm,ym];
end