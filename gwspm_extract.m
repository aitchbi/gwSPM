% This function is part of the toolbox:
%       gwSPM: Graph-based, Wavelet-based Statistical Parametric Mapping
%       (v1.00)
%
% 	Author: Hamid Behjat
% 
%   Biomedical Signal Processing Group, 
%   Dept. of Biomedical Engineering,
%   Lund University, Sweden
% 
%   June 2015
%
function varargout = gwspm_extract(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gwspm_extract_OpeningFcn, ...
    'gui_OutputFcn',  @gwspm_extract_OutputFcn, ...
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


% --- Executes just before gwspm_extract is made visible.
function gwspm_extract_OpeningFcn(hObject, ~, handles, varargin)

handles.initialTH =0.1;

if  ~isempty(varargin) && isstruct(varargin{1})
    handles.inputOptions = varargin{1};
    Im = strcat(handles.inputOptions.templateDirRoot,filesep,...
        handles.inputOptions.templateDirName,filesep,...
        handles.inputOptions.templateName);
    handles.directCall = 0;
else
    try
        crapDir = pwd;
        cd('/Volumes/ADATA UFD/OHBM_2016/MATLAB/gwspm_templates/');
        Im = '/Volumes/ADATA UFD/OHBM_2016/MATLAB/gwspm_templates/w1Template_6.nii,1';
        cd(crapDir)
    catch
        [Im,~] = spm_select(1,'^w1Template_6','Load the MNI normalized Dartel Template (w1Template_6.nii)'); %'^wTemplate_6\.nii$'
    end
    handles.directCall = 1;
end

if strcmp(Im(end),'1')
    handles.volInfo = spm_vol(Im);
else
    handles.volInfo = spm_vol(strcat(Im,',1'));
end

dummy = spm_read_vols(handles.volInfo);
dummy(dummy<handles.initialTH)=0;
handles.GMvolume = dummy;


I{1} = strcat(handles.volInfo.fname(1:end-4),'_cerebrum.nii');
I{2} = strcat(handles.volInfo.fname(1:end-4),'_cerebellum.nii');
clear dummy
dummy = cell(2,1);
for count=1:2
    try
        dummy{count} = spm_read_vols(spm_vol(I{count}));
        volExists(count) = 1;
    catch
        volExists(count) = 0;
    end
end

if sum(volExists)==2
    
    message = sprintf('Previously editted cerebrum and cerebellum masks were available for the selected template.\n\nThese are also uploaded, together with the selected whole brain template.\n\nIf desired, these editions can be removed by pressing the ''Reset'' button.');
    uiwait(msgbox(message));
    
    handles.cbrGMvolume = dummy{1};
    handles.cbrMask = zeros(size(handles.GMvolume));
    temp = dummy{1};
    temp(temp<=0.5)=0;
    handles.cbrMask(find(temp))=1;
    
    handles.cblGMvolume = dummy{2};
    handles.cblMask = zeros(size(handles.GMvolume));
    temp = dummy{2};
    temp(temp<=0.5)=0;
    handles.cblMask(find(temp))=1;
    
    handles.GMmask = zeros(size(handles.GMvolume));
    temp = handles.GMvolume;
    temp(temp<=0.5) = 0;
    handles.GMmask(find(temp))=1;
    
else
    
    handles.cbrGMvolume = handles.GMvolume;
    handles.cbrMask = zeros(size(handles.GMvolume));
    temp = handles.GMvolume;
    temp(temp<=0.5)=0;
    handles.cbrMask(find(temp))=1;
    
    handles.cblMask = zeros(size(handles.GMvolume));
    handles.cblGMvolume = zeros(size(handles.GMvolume));
    
    handles.GMmask = handles.cbrMask;
end

[sz1,sz2,sz3] = size(handles.GMvolume);
handles.sz1 = sz1;
handles.sz2 = sz2;
handles.sz3 = sz3;

handles.currentOrientation = 2; 
set(handles.cb_axial, 'Value',0);
set(handles.cb_sagital, 'Value',0);
set(handles.cb_coronal, 'Value',1);

set(handles.s_sliceNumber, 'Min',1);
handles = updateCurrentSliceNumbers(handles);
handles = updateCurrentSlices(handles);

handles.markerSize = 3;
set(handles.cb_3by3, 'Value',1);
set(handles.cb_5by5, 'Value',0);
set(handles.cb_15by15, 'Value',0);
set(handles.cb_1by1, 'Value',0);

handles.markerDepth = 'onlyCurrentSlice';
set(handles.cb_currentSlice, 'Value',1);
set(handles.cb_pmOneSlice, 'Value',0);
set(handles.cb_pmThreeSlices, 'Value',0);
set(handles.cb_pmFiveSlices, 'Value',0);

handles.editCount = 0;

handles.backupsRemoved = 0;

handles.backup = [];

handles.resetCount = 0;

updatePlots(handles)

handles.output = hObject;

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = gwspm_extract_OutputFcn(hObject, ~, handles)

varargout{1} = handles.output;  


% --- Executes on slider movement.
function s_sliceNumber_Callback(hObject, ~, handles)

set(hObject,'Max',size(handles.GMvolume,handles.currentOrientation));
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


% --- Executes on button press in pb_toggleOff.
function pb_toggleOff_Callback(hObject, ~, handles)

if handles.editCount < 3
    if handles.editCount == 0
        message = sprintf('Mark out the cerebelum by ''left clicking'' on positions of interest.\n\nNote that all pixels within a range specified by the marker size will be selected.\n\nThe selected pixels will be toggled ''OFF''.\n\n''Double Click'' to apply changes and update the figures.');
        uiwait(msgbox(message));
    end
    handles.backup = handles;
    handles.editCount = handles.editCount+1;
elseif handles.editCount == 3
    handles.backup.backup.backup = [];
    handles.backup = handles;
    handles.backupsRemoved = 1;
end

[x,y] = getpts(handles.axes1);
xx = round(y); 
y = round(x);
x =xx;
x=x';
y=y';

switch handles.markerSize
    case 1
        x=x';
        y=y';
    case 3
        x = [x,x,x,x+1,x+1,x+1,x-1,x-1,x-1]';
        y = [y,y+1,y-1,y,y+1,y-1,y,y+1,y-1]';
    case 5
        x = [repmat(x,1,5),...
            repmat(x+1,1,5),repmat(x+2,1,5),...
            repmat(x-1,1,5),repmat(x-2,1,5),]';
        y = repmat([y,y+1,y+2,y-1,y-2],1,5)';
    case 15
        x = [repmat(x,1,15),...
            repmat(x+1,1,15),repmat(x+2,1,15),repmat(x+3,1,15),repmat(x+4,1,15),...
            repmat(x+5,1,15),repmat(x+6,1,15),repmat(x+7,1,15),...
            repmat(x-1,1,15),repmat(x-2,1,15),repmat(x-3,1,15),repmat(x-4,1,15),...
            repmat(x-5,1,15),repmat(x-6,1,15),repmat(x-7,1,15)]';
        y = repmat([y,y+1,y+2,y+3,y+4,y+5,y+6,y+7,...
            y-1,y-2,y-3,y-4,y-5,y-6,y-7],1,15)';
    
end

dummy1 = squeeze(handles.GMvolume_currentSlice);
dummy2 = size(dummy1,2);

x = dummy2-x+1; 
xx=y; 
yy=x;
N_pixels = numel(xx);

% Check whether selected ppoints are too close to border ------------------
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

check1 = numel(find(xx>0))==N_pixels;
check2 = numel(find(yy>0))==N_pixels;

check3 = numel(find(xx<=tempX))==N_pixels;
check4 = numel(find(yy<=tempY))==N_pixels;

check =  ~check1+~check2+~check3+~check4;
if check
    tooCloseToImageEdges=1;
end
% -------------------------------------------------------------------------

if ~tooCloseToImageEdges
    
    ind = sub2ind(size(dummy1),xx,yy);
    
    switch handles.markerDepth
        case 'onlyCurrentSlice'
            slicesToUpdate = 0;
        case '+-1'
            slicesToUpdate = [0,1,-1];
        case '+-3'
            slicesToUpdate = [0,1,2,3,-1,-2,-3];
        case '+-5'
            slicesToUpdate = [0,1,2,3,4,5,-1,-2,-3,-4,-5];
    end
    
    currSliceNum = handles.currentSliceNumber;
    
    % Check whether selected slice is too close to volume borders ---------
    tooCloseToVolumeBorder = 0;
    switch handles.currentOrientation
        case 1
            temp0 =handles.sz1;
        case 2
            temp0 =handles.sz2;
        case 3
            temp0 =handles.sz3;
    end
    check =  numel(find(currSliceNum+slicesToUpdate >temp0 ))...
        + numel(find(currSliceNum+slicesToUpdate < 1));
    
    if check
        tooCloseToVolumeBorder = 1;
    end
    % ---------------------------------------------------------------------
    
    if ~tooCloseToVolumeBorder
        
        
        for iDepth=slicesToUpdate
            handles.currentSliceNumber = currSliceNum+iDepth;
            handles = updateCurrentSlices(handles);
            dummy3 = squeeze(handles.cbrGMvolume_currentSlice);
            dummy4 = squeeze(handles.cbrMask_currentSlice);
            dummy5 = squeeze(handles.cblMask_currentSlice);
            dummy6 = squeeze(handles.cblGMvolume_currentSlice);
            dummy7 = squeeze(handles.GMvolume_currentSlice);
            dummy8 = squeeze(handles.GMmask_currentSlice);
            
            switch handles.currentOrientation
                case 1
                    dummy5(ind) =dummy8(ind);
                    dummy6(ind) =dummy7(ind);
                    
                    handles.cblMask(handles.currentSliceNumber,:,:) = dummy5;
                    handles.cblGMvolume(handles.currentSliceNumber,:,:) = dummy6;
                    
                    dummy3(ind) =0;
                    dummy4(ind) =0;
                    handles.cbrGMvolume(handles.currentSliceNumber,:,:) = dummy3;
                    handles.cbrMask(handles.currentSliceNumber,:,:) = dummy4;
                    
                case 2
                    dummy5(ind) =dummy8(ind);
                    dummy6(ind) =dummy7(ind);
                    
                    handles.cblMask(:,handles.currentSliceNumber,:) = dummy5;
                    handles.cblGMvolume(:,handles.currentSliceNumber,:) = dummy6;
                    
                    dummy3(ind) =0;
                    dummy4(ind) =0;
                    
                    handles.cbrGMvolume(:,handles.currentSliceNumber,:) = dummy3;
                    handles.cbrMask(:,handles.currentSliceNumber,:) = dummy4;
                    
                case 3
                    dummy5(ind) =dummy8(ind);
                    dummy6(ind) =dummy7(ind);
                    
                    handles.cblMask(:,:,handles.currentSliceNumber) = dummy5;
                    handles.cblGMvolume(:,:,handles.currentSliceNumber) = dummy6;
                    
                    dummy3(ind) =0;
                    dummy4(ind) =0;
                    handles.cbrGMvolume(:,:,handles.currentSliceNumber) = dummy3;
                    handles.cbrMask(:,:,handles.currentSliceNumber) = dummy4;
            end
        end
        handles.currentSliceNumber = currSliceNum; 
        handles = updateCurrentSlices(handles); 
        guidata(hObject, handles);
        
        updatePlots(handles);
        
    else
        fprintf('-------------------------------------------------------------\n')
        fprintf('Current slice is too close to the border of the volume for the selected marker depth. \n')
        fprintf('No update done.\n')
        fprintf('Either change the slice to a slice further towards the center of the volume.. \n')
        fprintf('or choose a smaller marker depth. \n')
        fprintf('-------------------------------------------------------------\n')
    end
else
    fprintf('-------------------------------------------------------------\n')
    fprintf('Selected points were to close to image borders for this ''marker size''. \n')
    fprintf('No update done.\n')
    fprintf('Please select points further away from the image border.\n')
    fprintf('-------------------------------------------------------------\n')
end


% --- Executes on button press in pb_toggleOn.
function pb_toggleOn_Callback(hObject, ~, handles)

if handles.editCount < 3
    if handles.editCount == 0
        message = sprintf('Toggle ''ON'' voxels that have been toggled off by ''left clicking'' on positions of interest.\n\nNote that all pixels within a range specified by the marker size will be selected.\n\nThose voxels that are already ON will be unaffected.\n\n''Right click'' to apply changes and update the figures.');
        uiwait(msgbox(message));
    end
    handles.backup = handles;
    handles.editCount = handles.editCount+1;
elseif handles.editCount == 3
    handles.backup.backup.backup = [];
    handles.backup = handles;
    handles.backupsRemoved = 1;
end

[x,y] = getpts(handles.axes1);
xx = round(y); 
y = round(x);
x =xx;

x=x';
y=y';

switch handles.markerSize
    
    case 1
        x=x';
        y=y';
    case 3
        x = [x,x,x,x+1,x+1,x+1,x-1,x-1,x-1]';
        y = [y,y+1,y-1,y,y+1,y-1,y,y+1,y-1]';
    case 5
        x = [repmat(x,1,5),...
            repmat(x+1,1,5),repmat(x+2,1,5),...
            repmat(x-1,1,5),repmat(x-2,1,5),]';
        y = repmat([y,y+1,y+2,y-1,y-2],1,5)';
    case 15
        x = [repmat(x,1,15),...
            repmat(x+1,1,15),repmat(x+2,1,15),repmat(x+3,1,15),repmat(x+4,1,15),...
            repmat(x+5,1,15),repmat(x+6,1,15),repmat(x+7,1,15),...
            repmat(x-1,1,15),repmat(x-2,1,15),repmat(x-3,1,15),repmat(x-4,1,15),...
            repmat(x-5,1,15),repmat(x-6,1,15),repmat(x-7,1,15)]';
        y = repmat([y,y+1,y+2,y+3,y+4,y+5,y+6,y+7,...
            y-1,y-2,y-3,y-4,y-5,y-6,y-7],1,15)';
end

dummy1 = squeeze(handles.GMvolume_currentSlice);
dummy2 = size(dummy1,2);

x = dummy2-x+1; 
xx=y;
yy=x;
N_pixels = numel(xx);

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

check1 = numel(find(xx>0))==N_pixels;
check2 = numel(find(yy>0))==N_pixels;

check3 = numel(find(xx<=tempX))==N_pixels;
check4 = numel(find(yy<=tempY))==N_pixels;

check =  ~check1+~check2+~check3+~check4;
if check
    tooCloseToImageEdges=1;
end
% ---------------------------------------------------------

if ~tooCloseToImageEdges
    
    ind = sub2ind(size(dummy1),xx,yy);
    
    switch handles.markerDepth
        case 'onlyCurrentSlice'
            slicesToUpdate = 0;
        case '+-1'
            slicesToUpdate = [0,1,-1];
        case '+-3'
            slicesToUpdate = [0,1,2,3,-1,-2,-3];
        case '+-5'
            slicesToUpdate = [0,1,2,3,4,5,-1,-2,-3,-4,-5];
    end
    
    currSliceNum = handles.currentSliceNumber;
    
    % Check whether selected slice is too close to volume borders -----
    tooCloseToVolumeBorder = 0;
    switch handles.currentOrientation
        case 1
            temp0 =handles.sz1;
        case 2
            temp0 =handles.sz2;
        case 3
            temp0 =handles.sz3;
    end
    check =  numel(find(currSliceNum+slicesToUpdate >temp0 ))...
        + numel(find(currSliceNum+slicesToUpdate < 1));
    
    if check
        tooCloseToVolumeBorder = 1;
    end
    % ---------------------------------------------------------
    
    if ~tooCloseToVolumeBorder
        for iDepth=slicesToUpdate
            handles.currentSliceNumber = currSliceNum+iDepth;
            handles = updateCurrentSlices(handles);
            dummy3 = squeeze(handles.cbrGMvolume_currentSlice);
            dummy4 = squeeze(handles.cbrMask_currentSlice);
            dummy5 = squeeze(handles.cblMask_currentSlice);
            dummy6 = squeeze(handles.cblGMvolume_currentSlice);
            dummy7 = squeeze(handles.GMvolume_currentSlice);
            dummy8 = squeeze(handles.GMmask_currentSlice);
            
            switch handles.currentOrientation
                case 1
                    dummy5(ind) =0;
                    dummy6(ind) =0;
                    
                    handles.cblMask(handles.currentSliceNumber,:,:) = dummy5;
                    handles.cblGMvolume(handles.currentSliceNumber,:,:) = dummy6;
                    
                    dummy3(ind) =dummy7(ind);
                    dummy4(ind) =dummy8(ind);
                    handles.cbrGMvolume(handles.currentSliceNumber,:,:) = dummy3;
                    handles.cbrMask(handles.currentSliceNumber,:,:) = dummy4;
                    
                case 2
                    dummy5(ind) =0;
                    dummy6(ind) =0;
                    
                    handles.cblMask(:,handles.currentSliceNumber,:) = dummy5;
                    handles.cblGMvolume(:,handles.currentSliceNumber,:) = dummy6;
                    
                    dummy3(ind) =dummy7(ind);
                    dummy4(ind) =dummy8(ind);
                    
                    handles.cbrGMvolume(:,handles.currentSliceNumber,:) = dummy3;
                    handles.cbrMask(:,handles.currentSliceNumber,:) = dummy4;
                    
                case 3
                    dummy5(ind) =0;
                    dummy6(ind) =0;
                    
                    handles.cblMask(:,:,handles.currentSliceNumber) = dummy5;
                    handles.cblGMvolume(:,:,handles.currentSliceNumber) = dummy6;
                    
                    dummy3(ind) =dummy7(ind);
                    dummy4(ind) =dummy8(ind);
                    handles.cbrGMvolume(:,:,handles.currentSliceNumber) = dummy3;
                    handles.cbrMask(:,:,handles.currentSliceNumber) = dummy4;
            end
        end
        handles.currentSliceNumber = currSliceNum; 
        handles = updateCurrentSlices(handles); 
        guidata(hObject, handles);
        
        updatePlots(handles);
    else
        fprintf('-------------------------------------------------------------\n')
        fprintf('Current slice is too close to the border of the volume for the selected marker depth. \n')
        fprintf('No update done.\n')
        fprintf('Either change the slice to a slice further towards the center of the volume.. \n')
        fprintf('or choose a smaller marker depth. \n')
        fprintf('-------------------------------------------------------------\n')
    end
else
    fprintf('-------------------------------------------------------------\n')
    fprintf('Selected points were to close to image borders for this ''marker size''. \n')
    fprintf('No update done.\n')
    fprintf('Please select points further away from the image border.\n')
    fprintf('-------------------------------------------------------------\n')
end


% --- Executes on button press in pb_draw.
function pb_draw_Callback(hObject, ~, handles)

if handles.editCount < 3
    if handles.editCount == 0
        message = sprintf('''left click'' and hold on to mark an enclosed region.\n\nThe selected pixels within the contour will be toggled ''OFF''.\n\nRelease the mouse key to apply changes and update the figures.');
        uiwait(msgbox(message));
    end
    handles.backup = handles;
    handles.editCount = handles.editCount+1;
elseif handles.editCount == 3
    handles.backup.backup.backup = [];
    handles.backup = handles;
    handles.backupsRemoved = 1;
end

dummy1 = squeeze(handles.GMvolume_currentSlice);

hFH = imfreehand(handles.axes1);
bIm = hFH.createMask();
bIm = flipud(bIm);
ind = find(bIm');

switch handles.markerDepth
    case 'onlyCurrentSlice'
        slicesToUpdate = 0;
    case '+-1'
        slicesToUpdate = [0,1,-1];
    case '+-3'
        slicesToUpdate = [0,1,2,3,-1,-2,-3];
    case '+-5'
        slicesToUpdate = [0,1,2,3,4,5,-1,-2,-3,-4,-5];
end

currSliceNum = handles.currentSliceNumber;

for iDepth=slicesToUpdate
    
    handles.currentSliceNumber = currSliceNum+iDepth;
    handles = updateCurrentSlices(handles);
    
    dummy3 = squeeze(handles.cbrGMvolume_currentSlice);
    dummy4 = squeeze(handles.cbrMask_currentSlice);
    dummy5 = squeeze(handles.cblMask_currentSlice);
    dummy6 = squeeze(handles.cblGMvolume_currentSlice);
    dummy7 = squeeze(handles.GMvolume_currentSlice);
    dummy8 = squeeze(handles.GMmask_currentSlice);
    
    switch handles.currentOrientation
        case 1
            dummy5(ind) =dummy8(ind);
            dummy6(ind) =dummy7(ind);
            
            handles.cblMask(handles.currentSliceNumber,:,:) = dummy5;
            handles.cblGMvolume(handles.currentSliceNumber,:,:) = dummy6;
            
            dummy3(ind) =0;
            dummy4(ind) =0;
            handles.cbrGMvolume(handles.currentSliceNumber,:,:) = dummy3;
            handles.cbrMask(handles.currentSliceNumber,:,:) = dummy4;
            
        case 2
            dummy5(ind) =dummy8(ind);
            dummy6(ind) =dummy7(ind);
            
            handles.cblMask(:,handles.currentSliceNumber,:) = dummy5;
            handles.cblGMvolume(:,handles.currentSliceNumber,:) = dummy6;
            
            dummy3(ind) =0;
            dummy4(ind) =0;
            
            handles.cbrGMvolume(:,handles.currentSliceNumber,:) = dummy3;
            handles.cbrMask(:,handles.currentSliceNumber,:) = dummy4;
            
        case 3
            dummy5(ind) =dummy8(ind);
            dummy6(ind) =dummy7(ind);
            
            handles.cblMask(:,:,handles.currentSliceNumber) = dummy5;
            handles.cblGMvolume(:,:,handles.currentSliceNumber) = dummy6;
            
            dummy3(ind) =0;
            dummy4(ind) =0;
            handles.cbrGMvolume(:,:,handles.currentSliceNumber) = dummy3;
            handles.cbrMask(:,:,handles.currentSliceNumber) = dummy4;
    end
end
handles.currentSliceNumber = currSliceNum; 
handles = updateCurrentSlices(handles); 
guidata(hObject, handles);
updatePlots(handles);


% --- Executes on button press in cb_axial.
function cb_axial_Callback(hObject, ~, handles)

set(handles.cb_sagital, 'Value',0);
set(handles.cb_coronal, 'Value',0);
handles.currentOrientation = 3; 
handles = updateCurrentSliceNumbers(handles);
handles = updateCurrentSlices(handles);
updatePlots(handles);
guidata(hObject, handles);


% --- Executes on button press in cb_sagital.
function cb_sagital_Callback(hObject, ~, handles)

set(handles.cb_axial, 'Value',0);
set(handles.cb_coronal, 'Value',0);
handles.currentOrientation = 1; 
handles = updateCurrentSliceNumbers(handles);
handles = updateCurrentSlices(handles);
updatePlots(handles);
guidata(hObject, handles);


% --- Executes on button press in cb_sagital.
function cb_coronal_Callback(hObject, ~, handles)

set(handles.cb_axial, 'Value',0);
set(handles.cb_sagital, 'Value',0);
handles.currentOrientation = 2; 
handles = updateCurrentSliceNumbers(handles);
handles = updateCurrentSlices(handles);
updatePlots(handles);
guidata(hObject, handles);



% --- Executes on button press in cb_3by3.
function cb_3by3_Callback(hObject, ~, handles)

set(handles.cb_5by5, 'Value',0);
set(handles.cb_15by15, 'Value',0);
set(handles.cb_1by1, 'Value',0);
handles.markerSize=3;

updatePlots(handles)

guidata(hObject, handles);


% --- Executes on button press in cb_5by5.
function cb_5by5_Callback(hObject, ~, handles)

set(handles.cb_3by3, 'Value',0);
set(handles.cb_15by15, 'Value',0);
set(handles.cb_1by1, 'Value',0);
handles.markerSize=5;

updatePlots(handles)

guidata(hObject, handles);

% --- Executes on button press in cb_15by15.
function cb_15by15_Callback(hObject, ~, handles)

set(handles.cb_3by3, 'Value',0);
set(handles.cb_5by5, 'Value',0);
set(handles.cb_1by1, 'Value',0);
handles.markerSize=15;

updatePlots(handles)

guidata(hObject, handles);


% --- Executes on button press in cb_1by1.
function cb_1by1_Callback(hObject, ~, handles)

set(handles.cb_3by3, 'Value',0);
set(handles.cb_5by5, 'Value',0);
set(handles.cb_15by15, 'Value',0);
handles.markerSize=1;

updatePlots(handles)

guidata(hObject, handles);




% --- Executes on button press in cb_currentSlice.
function cb_currentSlice_Callback(hObject, ~, handles)

set(handles.cb_pmOneSlice, 'Value',0);
set(handles.cb_pmThreeSlices, 'Value',0);
set(handles.cb_pmFiveSlices, 'Value',0);
handles.markerDepth = 'onlyCurrentSlice';
guidata(hObject, handles);


% --- Executes on button press in cb_pmOneSlice.
function cb_pmOneSlice_Callback(hObject, ~, handles)

set(handles.cb_currentSlice, 'Value',0);
set(handles.cb_pmThreeSlices, 'Value',0);
set(handles.cb_pmFiveSlices, 'Value',0);
handles.markerDepth = '+-1';
guidata(hObject, handles);


% --- Executes on button press in cb_pmThreeSlices.
function cb_pmThreeSlices_Callback(hObject, ~, handles)

set(handles.cb_currentSlice, 'Value',0);
set(handles.cb_pmOneSlice, 'Value',0);
set(handles.cb_pmFiveSlices, 'Value',0);
handles.markerDepth = '+-3';
guidata(hObject, handles);


% --- Executes on button press in cb_pmFiveSlices.
function cb_pmFiveSlices_Callback(hObject, ~, handles)

set(handles.cb_currentSlice, 'Value',0);
set(handles.cb_pmOneSlice, 'Value',0);
set(handles.cb_pmThreeSlices, 'Value',0);
handles.markerDepth = '+-5';
guidata(hObject, handles);


% --- Executes on button press in pb_undo1step.
function pb_undo1step_Callback(hObject, ~, handles)

if isempty(handles.backup)
    if handles.backupsRemoved
        message = sprintf('No further undos can be done.\n\nA maximum of 3 undos are saved in memory.\n\n');
        uiwait(msgbox(message));
    else
        message = sprintf('There is no undo to be done.\n\nIt seems you have not done any further edits to be undone.\n\n');
        uiwait(msgbox(message));
    end
else
    currSl = handles.currentSliceNumber;
    currOr = handles.currentOrientation;
    handles = handles.backup;
    handles.currentSliceNumber = currSl;
    handles.currentOrientation = currOr;
    
    handles = updateCurrentSlices(handles);
    updatePlots(handles);
end
guidata(hObject, handles);


% --- Executes on button press in pb_save.
function pb_save_Callback(~, ~, handles)

V=handles.volInfo;
skipp = 0;
[pathstr,name,~] = fileparts(handles.volInfo.fname);
for i=1:2
    if i==2
        if ~handles.directCall
            name = 'GM';
        else
            skipp=1;
        end
    end
    if ~skipp
        V.fname =strcat(pathstr,filesep,name,'_cerebrum.nii');
        spm_write_vol(V,handles.cbrGMvolume);
        V.fname =strcat(pathstr,filesep,name,'_cerebrum_mask.nii');
        spm_write_vol(V,handles.cbrMask);
        V.fname = strcat(pathstr,filesep,name,'_cerebellum.nii');
        spm_write_vol(V,handles.cblGMvolume);
        V.fname = strcat(pathstr,filesep,name,'_cerebellum_mask.nii');
        spm_write_vol(V,handles.cblMask);
    end
end


% --- Executes on button press in pb_resetEditting.
function pb_resetEditting_Callback(hObject, ~, handles)

if handles.resetCount == 0
    message = sprintf('By pressing ''Reset'' all edits will be removed.\n\nIf this is what you want, press ''Reset'' again.');
    uiwait(msgbox(message));
    handles.resetCount = 1;
elseif handles.resetCount == 1
    handles.resetCount = 0;
    
    handles.cbrGMvolume = handles.GMvolume;
    handles.cbrMask = zeros(size(handles.GMvolume));
    handles.cbrMask(find(handles.GMvolume))=1;
    
    handles.cblMask = zeros(size(handles.GMvolume));
    handles.cblGMvolume = zeros(size(handles.GMvolume));
    handles = updateCurrentSlices(handles); 
    updatePlots(handles);
end
guidata(hObject, handles);


% --- Executes on button press in pb_proceed.
function pb_proceed_Callback(~, ~, handles)

close(handles.figure1)


function handles = updateCurrentSlices(handles)

switch handles.currentOrientation
    case 1
        handles.GMvolume_currentSlice = handles.GMvolume(handles.currentSliceNumber,:,:);
        handles.cbrGMvolume_currentSlice = handles.cbrGMvolume(handles.currentSliceNumber,:,:);
        handles.cbrMask_currentSlice = handles.cbrMask(handles.currentSliceNumber,:,:);
        handles.cblGMvolume_currentSlice = handles.cblGMvolume(handles.currentSliceNumber,:,:);
        handles.cblMask_currentSlice = handles.cblMask(handles.currentSliceNumber,:,:);
        handles.GMmask_currentSlice = handles.GMmask(handles.currentSliceNumber,:,:);
        
    case 2
        handles.GMvolume_currentSlice = handles.GMvolume(:,handles.currentSliceNumber,:);
        handles.cbrGMvolume_currentSlice = handles.cbrGMvolume(:,handles.currentSliceNumber,:);
        handles.cbrMask_currentSlice = handles.cbrMask(:,handles.currentSliceNumber,:);
        handles.cblGMvolume_currentSlice = handles.cblGMvolume(:,handles.currentSliceNumber,:);
        handles.cblMask_currentSlice = handles.cblMask(:,handles.currentSliceNumber,:);
        handles.GMmask_currentSlice = handles.GMmask(:,handles.currentSliceNumber,:);
        
    case 3
        handles.GMvolume_currentSlice = handles.GMvolume(:,:,handles.currentSliceNumber);
        handles.cbrGMvolume_currentSlice = handles.cbrGMvolume(:,:,handles.currentSliceNumber);
        handles.cbrMask_currentSlice = handles.cbrMask(:,:,handles.currentSliceNumber);
        handles.cblGMvolume_currentSlice = handles.cblGMvolume(:,:,handles.currentSliceNumber);
        handles.cblMask_currentSlice = handles.cblMask(:,:,handles.currentSliceNumber);
        handles.GMmask_currentSlice = handles.GMmask(:,:,handles.currentSliceNumber);
end

function handles = updateCurrentSliceNumbers(handles)

switch handles.currentOrientation
    case 1
        handles.currentSliceNumber = floor(handles.sz1/3);
        set(handles.s_sliceNumber, 'Max',handles.sz1);
    case 2
        handles.currentSliceNumber = floor(handles.sz2/3);
        set(handles.s_sliceNumber, 'Max',handles.sz2);
    case 3
        handles.currentSliceNumber = floor(handles.sz3/3);
        set(handles.s_sliceNumber, 'Max',handles.sz3);
end
set(handles.s_sliceNumber, 'Value',handles.currentSliceNumber);


function updatePlots(handles)

distToEdges = 1; 
if isempty(handles.markerSize)
    dummy = squeeze(handles.cbrGMvolume_currentSlice);
    imagesc(flipud(dummy'),'Parent', handles.axes1), colormap(gray)
else
    dummy = squeeze(handles.cbrGMvolume_currentSlice);
    dummy = flipud(dummy');
    dummy1 = zeros(size(dummy));
    dummy1(distToEdges+1:distToEdges+handles.markerSize,end-(distToEdges+handles.markerSize-1):end-distToEdges)=1;
    dummy = dummy+dummy1;
    imagesc(dummy,'Parent', handles.axes1), colormap(gray)
end

dummy = squeeze(handles.cblMask_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes5),colormap(gray)
dummy = squeeze(handles.GMvolume_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes3), colormap(gray)
dummy = squeeze(handles.cblGMvolume_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes2), colormap(gray)
dummy = squeeze(handles.cbrMask_currentSlice);
imagesc(flipud(dummy'),'Parent', handles.axes4), colormap(gray)

set(handles.axes1,'XtickLabel',[],'YtickLabel',[],'Xtick',[],'Ytick',[]);
set(handles.axes2,'XtickLabel',[],'YtickLabel',[],'Xtick',[],'Ytick',[]);
set(handles.axes3,'XtickLabel',[],'YtickLabel',[],'Xtick',[],'Ytick',[]);
set(handles.axes4,'XtickLabel',[],'YtickLabel',[],'Xtick',[],'Ytick',[]);
set(handles.axes5,'XtickLabel',[],'YtickLabel',[],'Xtick',[],'Ytick',[]);
