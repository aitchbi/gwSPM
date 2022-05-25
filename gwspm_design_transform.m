function varargout = gwspm_design_transform(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gwspm_design_transform_OpeningFcn, ...
    'gui_OutputFcn',  @gwspm_design_transform_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gwspm_design_transform is made visible.
function gwspm_design_transform_OpeningFcn(hObject, ~, handles, varargin)

if  ~isempty(varargin) && isstruct(varargin{1})
    handles.inputOptions = varargin{1};
    handles.directCall = 0;
else
    handles.directCall = 1;
end

dummy = load(strcat(fileparts(mfilename('fullpath')),filesep,'private',filesep,'flowLut.mat'));
handles.flowLut = dummy.flowLut;

handles.nans = 0; 

if isfield(handles,'inputOptions')
    handles.resolution = handles.inputOptions.resolution;
    switch handles.resolution
        case 2
            set(handles.rb_20mm, 'Value',1);
            set(handles.rb_25mm, 'Value',0);
            set(handles.rb_30mm, 'Value',0);
        case 2.5
            set(handles.rb_20mm, 'Value',0);
            set(handles.rb_25mm, 'Value',1);
            set(handles.rb_30mm, 'Value',0);
        case 3
            set(handles.rb_20mm, 'Value',0);
            set(handles.rb_25mm, 'Value',0);
            set(handles.rb_30mm, 'Value',1);
    end
else
    handles.resolution = 3;
    set(handles.rb_20mm, 'Value',0);
    set(handles.rb_25mm, 'Value',0);
    set(handles.rb_30mm, 'Value',1);
end

handles.firstResChange = 1;

handles.currentOrientation = 2; 
set(handles.rb_axial, 'Value',0);
set(handles.rb_sagital, 'Value',0);
set(handles.rb_coronal, 'Value',1);

handles.firstLoad = 1;

handles = updateData(handles);

set(handles.s_sliceNumber, 'Min',1);

handles.changeView = 0;

handles = updateCurrentSliceNumbers(handles);
handles = updateCurrentSlices(handles);
updatePlots(handles);


handles.frame.choice = 'cbr';
set(handles.rb_choice_cbr, 'Value',1);
set(handles.rb_choice_cbl, 'Value',0);

handles.frame.wav_scales = 2;
set(handles.rb_2scales, 'Value',1);
set(handles.rb_3scales, 'Value',0);
set(handles.rb_4scales, 'Value',0);

maxShift =30;
set(handles.s_shift,'Max',maxShift);
set(handles.s_shift, 'Min',1);
set(handles.s_shift,'SliderStep',[1/maxShift, 0.1000]);
handles.frame.shift = 20;
set(handles.s_shift, 'Value',handles.frame.shift);

maxOrder =200;
set(handles.s_chebyOrder,'Max',maxOrder);
set(handles.s_chebyOrder, 'Min',10);
set(handles.s_chebyOrder,'SliderStep',[1/maxOrder, 0.1000]);
handles.frame.chebyOrder = 50;
set(handles.s_chebyOrder, 'Value',handles.frame.chebyOrder);
set(handles.t_chebyOrder,'String',handles.frame.chebyOrder);

nChunks =10;
set(handles.s_zoom,'Max',nChunks);
set(handles.s_zoom, 'Min',1);
set(handles.s_zoom,'SliderStep',[1/nChunks, 0.1000]);
handles.frame.zoom = 1;
set(handles.s_zoom, 'Value',handles.frame.zoom);

handles = updateFrame(handles,'initialize');

handles.currentNode = [];

handles.parOption = 0; % paralle pool check
handles.poolOn=0;
handles.firstEst = 1;

handles.firstNodeChoice = 1;
handles.firstChebyChange = 1;

handles.output = hObject;

guidata(hObject, handles);

uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gwspm_design_transform_OutputFcn(~, ~, handles)

varargout{1} = handles.frame;
varargout{2} = handles.graph;
varargout{3} = handles.resolution;
varargout{4} = handles.parOption;
varargout{5} = get(handles.cb_proceed,'value');

delete(handles.figure1);


% --- Executes on slider movement.
function s_sliceNumber_Callback(hObject, ~, handles)

set(hObject,'Max',size(handles.GM,handles.currentOrientation));
switch handles.currentOrientation
    case 1
        set(hObject,'SliderStep',[1/handles.sz1, 0.1000]);
    case 2
        set(hObject,'SliderStep',[1/handles.sz2, 0.1000]);
    case 3
        set(hObject,'SliderStep',[1/handles.sz3, 0.1000]);
end
handles.currentSliceNumber = round(get(hObject,'Value'));

handles = updateCurrentSlices(handles);
updatePlots(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function s_sliceNumber_CreateFcn(hObject, ~, ~)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in rb_axial.
function rb_axial_Callback(hObject, ~, handles)

set(handles.rb_sagital, 'Value',0);
set(handles.rb_coronal, 'Value',0);
handles.currentOrientation = 3; % 1:sagital, 2:coronal, 3:axial,

handles.changeView = 1;
handles = updateCurrentSliceNumbers(handles);
handles.changeView = 0;

handles = updateCurrentSlices(handles);
updatePlots(handles);
guidata(hObject, handles);


% --- Executes on button press in rb_coronal.
function rb_coronal_Callback(hObject, ~, handles)

set(handles.rb_axial, 'Value',0);
set(handles.rb_sagital, 'Value',0);
handles.currentOrientation = 2; % 1:sagital, 2:coronal, 3:axial,

handles.changeView = 1;
handles = updateCurrentSliceNumbers(handles);
handles.changeView = 0;

handles = updateCurrentSlices(handles);
updatePlots(handles);
guidata(hObject, handles);


% --- Executes on button press in rb_sagital.
function rb_sagital_Callback(hObject, ~, handles)

set(handles.rb_axial, 'Value',0);
set(handles.rb_coronal, 'Value',0);
handles.currentOrientation = 1; % 1:sagital, 2:coronal, 3:axial,

handles.changeView = 1;
handles = updateCurrentSliceNumbers(handles);
handles.changeView = 0;

handles = updateCurrentSlices(handles);
updatePlots(handles);
guidata(hObject, handles);


% --- Executes on button press in rb_20mm.
function rb_20mm_Callback(hObject, ~, handles)

if handles.firstResChange && isfield(handles,'inputOptions')
    if handles.inputOptions.resolution == 2.5
        dummy = ' 2.5';
    else
        dummy = num2str(handles.inputOptions.resolution);
    end
    
    message = strcat('Note that the contrast maps using which you',...
        sprintf(' have performed second level analysis are in %s',dummy),...
        ' mm cubic resolution.',...
        sprintf('\n\nHowever, here in this GUI, you may explore a different resolution.'),...
        sprintf('\n\n''BUT'' you may not proceed to the next step'),...
        sprintf(' of the algorithm with a resolution other than %s',dummy),' mm cubic');
    uiwait(msgbox(message));
    
    handles.firstResChange = 0;
    set(handles.cb_proceed,'Value',0)

elseif get(handles.cb_proceed,'Value') &&...
        isfield(handles,'inputOptions') &&...
        handles.inputOptions.resolution~=2
    
    if handles.inputOptions.resolution == 2.5
        dummy = ' 2.5';
    else
        dummy = num2str(handles.inputOptions.resolution);
    end
    
    message = strcat('Note that the contrast maps using which you',...
        sprintf(' have performed second level analysis are in %s',dummy),...
        ' mm cubic resolution.',...
        sprintf('\n\nHowever, here in this GUI, you may explore a different resolution.'),...
        sprintf('\n\n''BUT'' you may not proceed to the next step'),...
        sprintf(' of the algorithm with a resolution other than %s',dummy),' mm cubic',...
        sprintf('\n\nThe ''check'' for finalizing will now be removed.\n\n'));
    uiwait(msgbox(message));
    
    set(handles.cb_proceed,'Value',0)
end

set(handles.rb_20mm, 'Value',1);
set(handles.rb_25mm, 'Value',0);
set(handles.rb_30mm, 'Value',0);
handles.resolution = 2;

handles = updateData(handles);
handles = updateCurrentSliceNumbers(handles);
handles = updateCurrentSlices(handles);
updatePlots(handles);

handles = updateFrame(handles,'resolutionChange');
guidata(hObject, handles);


% --- Executes on button press in rb_25mm.
function rb_25mm_Callback(hObject, ~, handles)

if handles.firstResChange && isfield(handles,'inputOptions')
    if handles.inputOptions.resolution == 2.5
        dummy = ' 2.5';
    else
        dummy = num2str(handles.inputOptions.resolution);
    end
    
    message = strcat('Note that the contrast maps using which you',...
        sprintf(' have performed second level analysis are in %s',dummy),...
        ' mm cubic resolution.',...
        sprintf('\n\nHowever, here in this GUI, you may explore a different resolution.'),...
        sprintf('\n\n''BUT'' you may not proceed to the next step'),...
        sprintf(' of the algorithm with a resolution other than %s',dummy),' mm cubic');
    uiwait(msgbox(message));
    
    handles.firstResChange = 0;
    set(handles.cb_proceed,'Value',0)

elseif get(handles.cb_proceed,'Value') &&...
        isfield(handles,'inputOptions') &&...
        handles.inputOptions.resolution~=2.5
    
    if handles.inputOptions.resolution == 2.5
        dummy = ' 2.5';
    else
        dummy = num2str(handles.inputOptions.resolution);
    end
    
    message = strcat('Note that the contrast maps using which you',...
        sprintf(' have performed second level analysis are in %s',dummy),...
        ' mm cubic resolution.',...
        sprintf('\n\nHowever, here in this GUI, you may explore a different resolution.'),...
        sprintf('\n\n''BUT'' you may not proceed to the next step'),...
        sprintf(' of the algorithm with a resolution other than %s',dummy),' mm cubic',...
        sprintf('\n\nThe ''check'' for finalizing will now be removed.\n\n'));    uiwait(msgbox(message));
    set(handles.cb_proceed,'Value',0)
end

set(handles.rb_20mm, 'Value',0);
set(handles.rb_25mm, 'Value',1);
set(handles.rb_30mm, 'Value',0);
handles.resolution = 2.5;

handles = updateData(handles);
handles = updateCurrentSliceNumbers(handles);
handles = updateCurrentSlices(handles);
updatePlots(handles);

handles = updateFrame(handles,'resolutionChange');
guidata(hObject, handles);


% --- Executes on button press in rb_30mm.
function rb_30mm_Callback(hObject, ~, handles)

if handles.firstResChange && isfield(handles,'inputOptions')
    if handles.inputOptions.resolution == 2.5
        dummy = ' 2.5';
    else
        dummy = num2str(handles.inputOptions.resolution);
    end
    
    message = strcat('Note that the contrast maps using which you',...
        sprintf(' have performed second level analysis are in %s',dummy),...
        ' mm cubic resolution.',...
        sprintf('\n\nHowever, here in this GUI, you may explore a different resolution.'),...
        sprintf('\n\n''BUT'' you may not proceed to the next step'),...
        sprintf(' of the algorithm with a resolution other than %s',dummy),' mm cubic');
    
    uiwait(msgbox(message));
    
    handles.firstResChange = 0;
    set(handles.cb_proceed,'Value',0)

elseif get(handles.cb_proceed,'Value') &&...
        isfield(handles,'inputOptions') &&...
        handles.inputOptions.resolution~=3
    
    if handles.inputOptions.resolution == 2.5
        dummy = ' 2.5';
    else
        dummy = num2str(handles.inputOptions.resolution);
    end
    message = strcat('Note that the contrast maps using which you',...
        sprintf(' have performed second level analysis are in %s',dummy),...
        ' mm cubic resolution.',...
        sprintf('\n\nHowever, here in this GUI, you may explore a different resolution.'),...
        sprintf('\n\n''BUT'' you may not proceed to the next step'),...
        sprintf(' of the algorithm with a resolution other than %s',dummy),' mm cubic',...
        sprintf('\n\nThe ''check'' for finalizing will now be removed.\n\n'));

    uiwait(msgbox(message));
    set(handles.cb_proceed,'Value',0)
end

set(handles.rb_20mm, 'Value',0);
set(handles.rb_25mm, 'Value',0);
set(handles.rb_30mm, 'Value',1);
handles.resolution = 3;

handles = updateData(handles);
handles = updateCurrentSliceNumbers(handles);
handles = updateCurrentSlices(handles);
updatePlots(handles);

handles = updateFrame(handles,'resolutionChange');
guidata(hObject, handles);


% --- Executes on button press in rb_choice_cbr.
function rb_choice_cbr_Callback(hObject, ~, handles)

handles.frame.choice = 'cbr';
set(handles.rb_choice_cbl,'Value',0);
handles = updateFrame(handles,'choiceChange');
guidata(hObject, handles);


% --- Executes on button press in rb_choice_cbl.
function rb_choice_cbl_Callback(hObject, ~, handles)

handles.frame.choice = 'cbl';
set(handles.rb_choice_cbr,'Value',0);
handles = updateFrame(handles,'choiceChange');
guidata(hObject, handles);


% --- Executes on button press in rb_2scales.
function rb_2scales_Callback(hObject, ~, handles)

set(handles.rb_2scales, 'Value',1);
set(handles.rb_3scales, 'Value',0);
set(handles.rb_4scales, 'Value',0);
handles.frame.wav_scales = 2;

handles = updateFrame(handles,'scalesChange');
handles = updateAtoms(handles,'scalesChange');

guidata(hObject, handles);


% --- Executes on button press in rb_3scales.
function rb_3scales_Callback(hObject, ~, handles)

set(handles.rb_2scales, 'Value',0);
set(handles.rb_3scales, 'Value',1);
set(handles.rb_4scales, 'Value',0);
handles.frame.wav_scales = 3;

handles = updateFrame(handles,'scalesChange');
handles = updateAtoms(handles,'scalesChange');

guidata(hObject, handles);


% --- Executes on button press in rb_4scales.
function rb_4scales_Callback(hObject, ~, handles)

set(handles.rb_2scales, 'Value',0);
set(handles.rb_3scales, 'Value',0);
set(handles.rb_4scales, 'Value',1);
handles.frame.wav_scales = 4;

handles = updateFrame(handles,'scalesChange');
handles = updateAtoms(handles,'scalesChange');

guidata(hObject, handles);


% --- Executes on slider movement.
function s_shift_Callback(hObject, ~, handles)

handles.frame.shift = round(get(hObject,'Value'));

handles = updateFrame(handles,'shiftChange');
handles = updateAtoms(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function s_shift_CreateFcn(hObject, ~, ~)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function s_chebyOrder_Callback(hObject, ~, handles)

if handles.firstChebyChange
    message = strcat('Select an order, which makes the ''dashed'' line approximately straight.',...
        sprintf('\n\n(i.e. the dashed line that is seen in the plot of the Chebyshev polynomial approximaton of kernels)'));
    uiwait(msgbox(message));
    handles.firstChebyChange = 0;
end

handles.frame.chebyOrder = round(get(hObject,'Value'));
set(handles.t_chebyOrder,'String',round(get(hObject,'Value')));

handles = updateFrame(handles,'chebyChange');

handles = updateAtoms(handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function s_chebyOrder_CreateFcn(hObject, ~, ~)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function s_zoom_Callback(hObject, ~, handles)

handles.frame.zoom = round(get(hObject,'Value'));
handles = updateFrame(handles,'zoomChange');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function s_zoom_CreateFcn(hObject, ~, ~)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pb_selectNode.
function pb_selectNode_Callback(hObject, ~, handles)

set(handles.txt_clickHere,'String','Here!')

if handles.firstNodeChoice
    message = sprintf('''double click'' on a voxel within the gray matter.\n\nThe point should be chosen in the figure tagged with ''Here!''.\n\n');
    uiwait(msgbox(message));
    handles.firstNodeChoice = 0;
end

[x,y] = getpts(handles.axes7);

set(handles.txt_clickHere,'String','')

xx = round(y); % due to MATLAB's conventions
y = round(x);
x =xx;

if numel(x)>1 || numel(y)>1
    message = sprintf('More than one point was selectted.\n\nRetry selecting a node.\n\nOnly do a single ''double click'' somewhere on the gray matter.');
    uiwait(msgbox(message));
    return
end

dummy2 = size(squeeze(handles.GM_mask_currentSlice),2);

% since we transpose and flipud when plotting. we do the reverse here:
x = dummy2-x+1; % flipud
xx=y; % & then transpose
yy=x;

% Check whether selected ppoints are too close to border ----------
tooCloseToImageEdges = 0;
switch handles.currentOrientation
    case 1
        tempX =handles.sz2;
        tempY =handles.sz3;
    case 2
        tempX =handles.sz1;
        tempY =handles.sz3;
    case 3
        tempX =handles.sz1;
        tempY =handles.sz2;
end

check1 = xx>0;
check2 = yy>0;

check3 = xx<=tempX;
check4 = yy<=tempY;

check =  ~check1+~check2+~check3+~check4;
if check
    tooCloseToImageEdges=1;
end

if ~tooCloseToImageEdges
    
    switch handles.currentOrientation
        case 1
            node = sub2ind([handles.sz1,handles.sz2,handles.sz3],...
                handles.currentSliceNumber,xx,yy);
        case 2
            node = sub2ind([handles.sz1,handles.sz2,handles.sz3],...
                xx,handles.currentSliceNumber,yy);
        case 3
            node = sub2ind([handles.sz1,handles.sz2,handles.sz3],...
                xx,yy,handles.currentSliceNumber);
    end
    
    if ismember(node,handles.graph.indGM)
        [~,indCenter] = ismember(node,handles.graph.indCbr);
        if indCenter
            ind = handles.graph.indCbr;
            atoms = sgwt_cheby_op(sgwt_delta(numel(ind),indCenter),...
                handles.graph.L_cbr,handles.frame.c_cbr,...
                handles.frame.arange_cbr);
            indElse = handles.graph.indCbl;
            
        else
            [~,indCenter] = ismember(node,handles.graph.indCbl);
            ind = handles.graph.indCbl;
            atoms = sgwt_cheby_op(sgwt_delta(numel(ind),indCenter),...
                handles.graph.L_cbl,handles.frame.c_cbl,...
                handles.frame.arange_cbl);
            indElse = handles.graph.indCbr;
        end
        
        handles.wave0(ind) = atoms{1};
        handles.wave0(indElse) = 0;
        handles.wave1(ind) = atoms{2};
        handles.wave1(indElse) = 0;
        handles.wave2(ind) = atoms{3};
        handles.wave2(indElse) = 0;
        
        if handles.frame.wav_scales == 3
            handles.wave3(ind) = atoms{4};
            handles.wave3(indElse) = 0;
        elseif handles.frame.wav_scales == 4
            handles.wave3(ind) = atoms{4};
            handles.wave3(indElse) = 0;
            handles.wave4(ind) = atoms{5};
            handles.wave4(indElse) = 0;
        else
            
        end
        
        handles.currentNode = node; 
        
    end    
else
    message = sprintf('You have either click outside the specified figure, or you have chosen a point too close to the border.\n\nPlease select another point.\n\n');
    uiwait(msgbox(message));
end

handles = updateCurrentSlices(handles);
updatePlots(handles);
guidata(hObject, handles);


% --- Executes on button press in pb_estimateTime.
function pb_estimateTime_Callback(hObject, ~, handles)

set(handles.t_timeInfo1,'String',[])
set(handles.t_timeInfo2,'String',[])

if handles.firstEst
    message = strcat(sprintf('The estimation takes a few minutes.'),...
        sprintf('\n\nMeanwhile, please do not change the settings.'),...
        sprintf('\n\nAlso note that this estimate is based on the current status of your computer.'),...
        sprintf('\n\nAs such, if you run the algorithm (step. 5) over a night, or a weekened,'),...
        ' the actual time can be close to the minimal estimate.');
    uiwait(msgbox(message));
    handles.firstEst = 0;
end

set(handles.t_waitEstInProg,'String','Wait! estimation in progres..')
pause(1)

trans.wav_scales = handles.frame.wav_scales;
trans.wav_dim = size(handles.GM);
trans.cbr.arange = handles.frame.arange_cbr;
trans.cbr.c = handles.frame.c_cbr;
trans.cbr.L = handles.graph.L_cbr;
trans.cbr.indices = handles.graph.indCbr;

trans.cbl.arange = handles.frame.arange_cbl;
trans.cbl.c = handles.frame.c_cbl;
trans.cbl.L = handles.graph.L_cbl;
trans.cbl.indices = handles.graph.indCbl;

switch handles.parOption
    case 0
        option = [];
    case 1
        option = 'parallel';
        
       if ~handles.poolOn
             parpool();
             handles.poolOn=1;
        end
end

[tMx,tMn,option] = gwspm_estimate_computeTime(trans,24,48,option);

if strcmp(option,'parallel')
    handles.parOption = 1;
else
    handles.parOption = 0;
end

set(handles.t_waitEstInProg,'String',[])

if handles.resolution~=2.5
    dummy1 = sprintf('For the current settings, i.e. GM template resolution: %d mm cubic, number of wavelet scales: %d, Chebyshev polynomial oder: %d',...
        handles.resolution,handles.frame.wav_scales,handles.frame.chebyOrder);
else
    dummy1 = sprintf('For the current settings, i.e. GM template resolution: 2.5 mm cubic, number of wavelet scales: %d, Chebyshev polynomial oder: %d',...
        handles.frame.wav_scales,handles.frame.chebyOrder);
end
dummy = strcat(sprintf('%d to %d',round(tMn),round(tMx)),sprintf(' hours'));
set(handles.t_timeInfo1,'String',strcat(dummy1))
set(handles.t_timeInfo2,'String',dummy)

guidata(hObject, handles);


% --- Executes on button press in cb_proceed.
function cb_proceed_Callback(hObject, ~, handles)

if isfield(handles,'inputOptions') && handles.inputOptions.resolution~=handles.resolution
    
    if handles.inputOptions.resolution == 2.5
        dummy = ' 2.5';
    else
        dummy = num2str(handles.inputOptions.resolution);
    end
    message = strcat('Note that the contrast maps using which you',...
    sprintf(' have performed second level analysis are in %s',dummy),...
    ' mm cubic resolution.',...
        sprintf('\n\nHowever, here in this GUI, you may explore a different resolution.'),...
        sprintf('\n\n''BUT'' you may not proceed to the next step'),...
        sprintf(' of the algorithm with a resolution other than %s',dummy),' mm cubic');
    uiwait(msgbox(message));
    
    set(hObject,'Value',0)
    
else
    
    if get(hObject,'Value')
        handles.graph.L_cbr = sgwt_laplacian(...
            gwspm_compute_adjacency(handles.cbr_mask,26,'Weight','no'),...
            'opt','normalized');
        
        handles.graph.L_cbl = sgwt_laplacian(...
            gwspm_compute_adjacency(handles.cbl_mask,26,'Weight','no'),...
            'opt','normalized');
    else
        handles.graph.L_cbr = [];
        handles.graph.L_cbl = [];
    end
    
end

guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, handles)

if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
      delete(hObject);
end


% --- Executes on button press in pb_saveAtoms.
function pb_saveAtoms_Callback(hObject, ~, handles)

try
    dummy1 = sprintf('sgwt_sc%d_sh%d_or%d',handles.wav_scales,handles.shift,handles.chebyOrder);
    
    currTransDir  = strcat(handles.inputOptions.atomsDir,filesep,dummy1);
    
    if ~exist(currTransDir,'dir')
        mkdir(currTransDir)
    end
    
    currNodeDir = strcat(currTransDir,filesep,'node_index_',num2str(handles.currentNode));
    if ~exist(currNodeDir,'dir')
        mkdir(currNodeDir)
    end
    
    V = handles.volInfo;
    V.fname= strcat(currNodeDir,filesep,'atom_scaling.nii');
    spm_write_vol(V,handles.wave0);
    
    V.fname= strcat(currNodeDir,filesep,'atom_wav_scale_1.nii');
    spm_write_vol(V,handles.wave1);
    
    V.fname= strcat(currNodeDir,filesep,'atom_wav_scale_2.nii');
    spm_write_vol(V,handles.wave2);
    
    if handles.wav_scales>2
        V.fname= strcat(currNodeDir,filesep,'atom_wav_scale_3.nii');
        spm_write_vol(V,handles.wave3);
    end
    if handles.wav_scales>3
        V.fname= strcat(currNodeDir,filesep,'atom_wav_scale_4.nii');
        spm_write_vol(V,handles.wave4);
    end
catch
    try
        fprintf('Problem in saving atoms. No atoms saved.\n')
    catch
        
    end
end


function handles = updateAtoms(handles,varargin) 

if nargin>1
    option = varargin{1};
else
    option = ' ';
end

if strcmp(option,'scalesChange')
    switch handles.frame.wav_scales
        case 2
            handles.wave3 = zeros(size(handles.GM));
            handles.wave4 = zeros(size(handles.GM));
        case 3
            handles.wave4 = zeros(size(handles.GM));
    end
end

if ~isempty(handles.currentNode)
    if ismember(handles.currentNode,handles.graph.indGM)
        [~,indCenter] = ismember(handles.currentNode,handles.graph.indCbr);
        if indCenter
            ind = handles.graph.indCbr;
            atoms = sgwt_cheby_op(sgwt_delta(numel(ind),indCenter),...
                handles.graph.L_cbr,handles.frame.c_cbr,...
                handles.frame.arange_cbr);
        else
            [~,indCenter] = ismember(handles.currentNode,handles.graph.indCbl);
            ind = handles.graph.indCbl;
            atoms = sgwt_cheby_op(sgwt_delta(numel(ind),indCenter),...
                handles.graph.L_cbl,handles.frame.c_cbl,...
                handles.frame.arange_cbl);
        end
        
        handles.wave0(ind) = atoms{1};
        handles.wave1(ind) = atoms{2};
        handles.wave2(ind) = atoms{3};
        if handles.frame.wav_scales == 3
            handles.wave3(ind) = atoms{4};
        elseif handles.frame.wav_scales == 4
            handles.wave3(ind) = atoms{4};
            handles.wave4(ind) = atoms{5};
        end
    else
        error('Something is fishy. Check why ended up here..')
    end
end
handles = updateCurrentSlices(handles);
updatePlots(handles);


function handles = updateData(handles)

switch handles.resolution
    case 2
        dummy = 'w2_';
    case 2.5
        dummy = 'w2point5_';
    case 3
        dummy = 'w3_';
end

if isfield(handles,'inputOptions')
    
    Im_GM_cbr = strcat(handles.inputOptions.templateDirRoot,...
        filesep,handles.inputOptions.templateDirName,...
        filesep,dummy, 'GM_cerebrum_adjusted.nii');
    Im_GM_cbl = strcat(handles.inputOptions.templateDirRoot,...
        filesep,handles.inputOptions.templateDirName,...
        filesep,dummy,'GM_cerebellum_adjusted.nii');
    
    Im_GM_cbr_mask = strcat(handles.inputOptions.templateDirRoot,...
        filesep,handles.inputOptions.templateDirName,...
        filesep,dummy,'GM_cerebrum_mask_adjusted.nii');
    Im_GM_cbl_mask = strcat(handles.inputOptions.templateDirRoot,...
        filesep,handles.inputOptions.templateDirName,...
        filesep,dummy,'GM_cerebellum_mask_adjusted.nii');
    
else
    
    try
        dummyPath = ' '; % '.../gwspm_templates/';
    catch
        fprintf('please define the path where the GM_cerebrum_adjusted.nii is located, in the blank space above. ')
    end
        
    Im_GM_cbr = strcat(dummyPath,dummy,'GM_cerebrum_adjusted.nii');
    
    Im_GM_cbl = strcat(dummyPath,dummy,'GM_cerebellum_adjusted.nii');
    
    Im_GM_cbr_mask = strcat(dummyPath,dummy,'GM_cerebrum_mask_adjusted.nii');
    
    Im_GM_cbl_mask = strcat(dummyPath,dummy,'GM_cerebellum_mask_adjusted.nii');
end

if strcmp(Im_GM_cbr(end),'1')
    handles.volInfo = spm_vol(Im_GM_cbr);
    
    handles.GM = spm_read_vols(spm_vol(Im_GM_cbr))...
        +spm_read_vols(spm_vol(Im_GM_cbl));
    
    handles.cbr_mask = spm_read_vols(spm_vol(Im_GM_cbr_mask));
    handles.cbl_mask = spm_read_vols(spm_vol(Im_GM_cbl_mask));
else
    handles.volInfo = spm_vol(strcat(Im_GM_cbr,',1'));
    
    handles.GM = spm_read_vols(spm_vol(strcat(Im_GM_cbr,',1')))...
        +spm_read_vols(spm_vol(strcat(Im_GM_cbl,',1')));
    
    handles.cbr_mask = spm_read_vols(spm_vol(strcat(Im_GM_cbr_mask,',1')));
    handles.cbl_mask = spm_read_vols(spm_vol(strcat(Im_GM_cbl_mask,',1')));
end

handles.graph.indCbr = find(handles.cbr_mask);
handles.graph.indCbl = find(handles.cbl_mask);

temp = numel(handles.graph.indCbl);
handles.graph.indCbl = handles.graph.indCbl(...
    ~ismember(handles.graph.indCbl,...
    handles.graph.indCbr));
handles.cbl_mask = zeros(size(handles.cbl_mask));
handles.cbl_mask(handles.graph.indCbl) = 1;

if temp ~= numel(handles.graph.indCbl)
    sprintf('There were %d voxels overlapping between the cerebrum and cerebellum masks.\n',...
        temp - numel(handles.graph.indCbl))
    sprintf('These voxels were removed from the cerebellum mask.\n')
end

handles.GM_mask = handles.cbr_mask+handles.cbl_mask;
handles.graph.indGM = find(handles.GM_mask);

handles.wave0 = zeros(size(handles.GM));
handles.wave1 = zeros(size(handles.GM));
handles.wave2 = zeros(size(handles.GM));
handles.wave3 = zeros(size(handles.GM));
handles.wave4 = zeros(size(handles.GM));

% Laplacian matrices
handles.graph.L_cbr = sgwt_laplacian(gwspm_compute_adjacency(...
    handles.cbr_mask,26,'Weight','no'),'opt','normalized');
handles.graph.L_cbl = sgwt_laplacian(gwspm_compute_adjacency(...
    handles.cbl_mask,26,'Weight','no'),'opt','normalized');

if ~handles.firstLoad
    switch handles.currentOrientation
        case 1
            handles.previousSz = handles.sz1;
        case 2
            handles.previousSz = handles.sz2;
        case 3
            handles.previousSz = handles.sz3;
    end
else
    handles.firstLoad = 0;
end

[sz1,sz2,sz3] = size(handles.GM);
handles.sz1 = sz1;
handles.sz2 = sz2;
handles.sz3 = sz3;


function handles = updateCurrentSliceNumbers(handles)

switch handles.currentOrientation
    case 1
        if isfield(handles,'previousSz') && ~handles.changeView
            handles.currentSliceNumber = floor(handles.sz1*(...
                handles.currentSliceNumber/handles.previousSz));
        else
            handles.currentSliceNumber = floor(handles.sz1/3);
        end
        set(handles.s_sliceNumber, 'Max',handles.sz1);
    case 2
        if isfield(handles,'previousSz') && ~handles.changeView
            handles.currentSliceNumber = floor(handles.sz2*(...
                handles.currentSliceNumber/handles.previousSz));
        else
            handles.currentSliceNumber = floor(handles.sz2/3);
        end
        set(handles.s_sliceNumber, 'Max',handles.sz2);
    case 3
        if isfield(handles,'previousSz') && ~handles.changeView
            handles.currentSliceNumber = floor(handles.sz3*(...
                handles.currentSliceNumber/handles.previousSz));
        else
            handles.currentSliceNumber = floor(handles.sz3/3);
        end
        set(handles.s_sliceNumber, 'Max',handles.sz3);
end
set(handles.s_sliceNumber, 'Value',handles.currentSliceNumber);


function handles = updateCurrentSlices(handles)

switch handles.currentOrientation
    case 1
        handles.GM_currentSlice = handles.GM(handles.currentSliceNumber,:,:);
        handles.GM_mask_currentSlice = handles.GM_mask(handles.currentSliceNumber,:,:);
        handles.wave0_currentSlice = handles.wave0(handles.currentSliceNumber,:,:);
        handles.wave1_currentSlice = handles.wave1(handles.currentSliceNumber,:,:);
        handles.wave2_currentSlice = handles.wave2(handles.currentSliceNumber,:,:);
        handles.wave3_currentSlice = handles.wave3(handles.currentSliceNumber,:,:);
        handles.wave4_currentSlice = handles.wave4(handles.currentSliceNumber,:,:);
    case 2
        handles.GM_currentSlice = handles.GM(:,handles.currentSliceNumber,:);
        handles.GM_mask_currentSlice = handles.GM_mask(:,handles.currentSliceNumber,:);
        handles.wave0_currentSlice = handles.wave0(:,handles.currentSliceNumber,:);
        handles.wave1_currentSlice = handles.wave1(:,handles.currentSliceNumber,:);
        handles.wave2_currentSlice = handles.wave2(:,handles.currentSliceNumber,:);
        handles.wave3_currentSlice = handles.wave3(:,handles.currentSliceNumber,:);
        handles.wave4_currentSlice = handles.wave4(:,handles.currentSliceNumber,:);
    case 3
        handles.GM_currentSlice = handles.GM(:,:,handles.currentSliceNumber);
        handles.GM_mask_currentSlice = handles.GM_mask(:,:,handles.currentSliceNumber);
        handles.wave0_currentSlice = handles.wave0(:,:,handles.currentSliceNumber);
        handles.wave1_currentSlice = handles.wave1(:,:,handles.currentSliceNumber);
        handles.wave2_currentSlice = handles.wave2(:,:,handles.currentSliceNumber);
        handles.wave3_currentSlice = handles.wave3(:,:,handles.currentSliceNumber);
        handles.wave4_currentSlice = handles.wave4(:,:,handles.currentSliceNumber);
end


function updatePlots(handles)

dummy = squeeze(handles.GM_mask_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes7),
colormap(handles.axes7,gray)

dummy = squeeze(handles.GM_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes8),
colormap(handles.axes8,gray)


dummy = squeeze(handles.wave0_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes1),
colormap(handles.axes1,handles.flowLut)
colorbar(handles.axes1)
dummy = max(abs(handles.wave0(:)));
if ~isnan(dummy) && dummy~=0
    set(handles.axes1,'CLim',[-dummy,dummy])
end
dummy = squeeze(handles.wave1_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes2),
colormap(handles.axes2,handles.flowLut)
colorbar(handles.axes2)
dummy = max(abs(handles.wave1(:)));
if ~isnan(dummy) && dummy~=0
    set(handles.axes2,'CLim',[-dummy,dummy])
end

dummy = squeeze(handles.wave2_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes3),
colormap(handles.axes3,handles.flowLut)
colorbar(handles.axes3)
dummy = max(abs(handles.wave2(:)));
if ~isnan(dummy) && dummy~=0
    set(handles.axes3,'CLim',[-dummy,dummy])
end


dummy = squeeze(handles.wave3_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes4),
colormap(handles.axes4,handles.flowLut)
colorbar(handles.axes4)
dummy = max(abs(handles.wave3(:)));
if ~isnan(dummy) && dummy~=0
    set(handles.axes4,'CLim',[-dummy,dummy])
end

dummy = squeeze(handles.wave4_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes5),
colormap(handles.axes5,handles.flowLut)
colorbar(handles.axes5)
dummy = max(abs(handles.wave4(:)));
if ~isnan(dummy) && dummy~=0
    set(handles.axes5,'CLim',[-dummy,dummy])
end

for i=[1:5,7,8]
    eval(strcat('set(handles.axes',num2str(i),...
        ',''XtickLabel'',[],''YtickLabel'',[],''Xtick'',[],''Ytick'',[])'))
end


function handles = updateFrame(handles,option)

if ~strcmp(option,'zoomChange')
    
    % arange --------------------------------------------------------------
    if strcmp(option, 'initialize') || strcmp(option, 'resolutionChange')
       
        handles.frame.arange_cbr(1) = 0;
        handles.frame.arange_cbr(2) =sgwt_rough_lmax(...
            sgwt_laplacian(gwspm_compute_adjacency(...
            handles.cbr_mask,26,'Weight','no'),'opt','normalized'));
        
        handles.frame.arange_cbl(1) = 0;
        handles.frame.arange_cbl(2) =sgwt_rough_lmax(...
            sgwt_laplacian(gwspm_compute_adjacency(...
            handles.cbl_mask,26,'Weight','no'),'opt','normalized'));
    end
    
    % Spectral kernels ----------------------------------------------------
    if strcmp(option, 'initialize') || strcmp(option,'resolutionChange')...
            || strcmp(option,'scalesChange') || strcmp(option,'shiftChange')
        
        handles.frame.g_cbr=gwspm_filter_design(handles.frame.arange_cbr(2),...
            handles.frame.wav_scales,handles.frame.shift,'designtype','mey',...
            'wav_type','meyer','lpfactor',0);
        handles.frame.g_cbl=gwspm_filter_design(handles.frame.arange_cbl(2),...
            handles.frame.wav_scales,handles.frame.shift,'designtype','mey',...
            'wav_type','meyer','lpfactor',0);
    end
    
    % Chebyshev poly. coeffs. ---------------------------------------------
    if strcmp(option, 'initialize') || strcmp(option,'resolutionChange')...
            || strcmp(option,'scalesChange') || strcmp(option,'shiftChange')...
            || strcmp(option,'chebyChange')
        
        for k=1:numel(handles.frame.g_cbr)
            handles.frame.c_cbr{k}=sgwt_cheby_coeff(handles.frame.g_cbr{k},...
                handles.frame.chebyOrder,handles.frame.chebyOrder+1,...
                handles.frame.arange_cbr);
            handles.frame.c_cbl{k}=sgwt_cheby_coeff(handles.frame.g_cbl{k},...
                handles.frame.chebyOrder,handles.frame.chebyOrder+1,...
                handles.frame.arange_cbl);
        end
    end
    
    if strcmp(handles.frame.choice,'cbr')
        dummy1 = handles.frame.g_cbr;
        dummy2 = handles.frame.arange_cbr;
    else
        dummy1 = handles.frame.g_cbl;
        dummy2 = handles.frame.arange_cbl;
    end
    
    % Frame plots ---------------------------------------------------------
    if strcmp(option, 'initialize') || strcmp(option,'resolutionChange')...
            || strcmp(option,'scalesChange') || strcmp(option,'shiftChange')...
            || strcmp(option,'choiceChange')
        
       gwspm_view_design(dummy1,dummy2,...
            'showLegend','no','plotLineWidth',2,'guiHandle',handles.axes10);
        
        gwspm_view_design(dummy1,dummy2,...
            'showLegend','no','plotLineWidth',2,'guiHandle',handles.axes11);
        
        % keep viewing max range fixed.. for better visual comparison
        dummy3 = [0, max(handles.frame.arange_cbr(2),handles.frame.arange_cbl(2))];
        
        set(handles.axes10,'XLim',dummy3,'YLim',[0,1.2],...
            'YTick',[0 1],'Box','off');
        
        set(handles.axes11,'XLim',[...
            dummy3(2)/10*(handles.frame.zoom-1),...
            dummy3(2)/10*handles.frame.zoom],...
            'YLim',[0,1.2],'YTick',[],'Box','off');
    end
    
    % Cheby. Poly. Approx. plots ------------------------------------------
    if strcmp(option, 'initialize') || strcmp(option,'resolutionChange')...
            || strcmp(option,'scalesChange') || strcmp(option,'shiftChange')...
            || strcmp(option,'chebyChange') || strcmp(option,'choiceChange')
        
        gwspm_view_design(dummy1,dummy2,...
            'showLegend','no','plotLineWidth',2,'guiHandle',handles.axes12,...
            'chebyOrder',handles.frame.chebyOrder,'onlyChebyApprox','yes');
        
        gwspm_view_design(dummy1,dummy2,...
            'showLegend','no','plotLineWidth',2,'guiHandle',handles.axes13,...
            'chebyOrder',handles.frame.chebyOrder,'onlyChebyApprox','yes');
        
        dummy3 = [0, max(handles.frame.arange_cbr(2),handles.frame.arange_cbl(2))];
        
        set(handles.axes12,'XLim',dummy3,'YLim',[0,1.2],...
            'YTick',[0 1],'Box','off');
        
        set(handles.axes13,'XLim',[...
            dummy3(2)/10*(handles.frame.zoom-1),...
            dummy3(2)/10*handles.frame.zoom],...
            'YLim',[0,1.2],'YTick',[],'Box','off');
        
    end
else
    % Adjust zoom region --------------------------------------------------
    dummy3 = [0, max(handles.frame.arange_cbr(2),handles.frame.arange_cbl(2))];
    set(handles.axes11,'XLim',[...
        dummy3(2)/10*(handles.frame.zoom-1),...
        dummy3(2)/10*handles.frame.zoom],...
        'YLim',[0,1.2],'YTick',[],'Box','off');
    
    set(handles.axes13,'XLim',[...
        dummy3(2)/10*(handles.frame.zoom-1),...
        dummy3(2)/10*handles.frame.zoom],...
        'YLim',[0,1.2],'YTick',[],'Box','off');
end
