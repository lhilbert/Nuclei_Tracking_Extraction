function varargout = TrackReview_GUI(varargin)
% TRACKREVIEW_GUI MATLAB code for TrackReview_GUI.fig
%      TRACKREVIEW_GUI, by itself, creates a new TRACKREVIEW_GUI or raises the existing
%      singleton*.
%
%      H = TRACKREVIEW_GUI returns the handle to a new TRACKREVIEW_GUI or the handle to
%      the existing singleton*.
%
%      TRACKREVIEW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKREVIEW_GUI.M with the given input arguments.
%
%      TRACKREVIEW_GUI('Property','Value',...) creates a new TRACKREVIEW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrackReview_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrackReview_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrackReview_GUI

% Last Modified by GUIDE v2.5 21-Sep-2016 23:09:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrackReview_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @TrackReview_GUI_OutputFcn, ...
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


% --- Executes just before TrackReview_GUI is made visible.
function TrackReview_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrackReview_GUI (see VARARGIN)

% Choose default command line output for TrackReview_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TrackReview_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TrackReview_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in showReject_checkbox.
function showReject_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showReject_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showReject_checkbox

redraw(hObject, eventdata, handles);



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[load_sourceFile,load_sourceDir] = uigetfile('*.*');
thisPath = fullfile(load_sourceDir,load_sourceFile);

data = load(thisPath);

data.loadPath = thisPath;

set(handles.frameSlider,'Max',data.numFrames,...
    'Min',1,'SliderStep',[1./(data.numFrames-1),10./data.numFrames],'Value',1);

data.activeFrame = get(handles.frameSlider,'Value');

if ~isfield(data,'rejectFlags')
    data.rejectFlags = false(1,numel(data.tracks));
end

if ~isfield(data,'minRatio')
    data.minRatio = str2num(get(handles.minRatioEdit,'String'));
else
    set(handles.minRatioEdit,'String',num2str(data.minRatio));
end
if ~isfield(data,'Vmin')
    data.Vmin = str2num(get(handles.minVolEdit,'String'));
else
    set(handles.minVolEdit,'String',num2str(data.Vmin));
end
if ~isfield(data,'Vmax')
    data.Vmax = str2num(get(handles.maxVolEdit,'String'));
else
    set(handles.maxVolEdit,'String',num2str(data.Vmax));
end
if ~isfield(data,'minLength')
    data.minLength = str2num(get(handles.minLengthEdit,'String'));
else
    set(handles.minLengthEdit,'String',num2str(data.minLength));
end
if ~isfield(data,'minX')
    handles.minX = str2num(get(handles.minXedit,'String'));
else
    set(handles.minXedit,'String',num2str(data.minX));
end
if ~isfield(data,'maxX')
    handles.maxX = str2num(get(handles.maxXedit,'String'));
else
    set(handles.maxXedit,'String',num2str(data.maxX));
end
if ~isfield(data,'minY')
    handles.minY = str2num(get(handles.minYedit,'String'));
else
    set(handles.minYedit,'String',num2str(data.minY));
end
if ~isfield(data,'maxY')
    handles.maxY = str2num(get(handles.maxYedit,'String'));
else
    set(handles.maxYedit,'String',num2str(data.maxY));
end
if ~isfield(data,'minZ')
    minZ = str2num(get(handles.minZedit,'String'));
else
    set(handles.minZedit,'String',num2str(data.minZ));
end
if ~isfield(data,'maxZ')
    maxZ = str2num(get(handles.maxZedit,'String'));
else
    set(handles.maxZedit,'String',num2str(data.maxZ));
end

handles.data = data;

set(handles.filename_textBox,'String',sprintf('Track set file:\n%s',...
    data.loadPath))

guidata(hObject,handles);

updateLimiters(hObject, eventdata, handles);





% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);

inLimitsFlags = handles.data.inLimitsFlags;
rejectFlags = handles.data.rejectFlags;

minRatio = str2num(get(handles.minRatioEdit,'String'));
Vmin = str2num(get(handles.minVolEdit,'String'));
Vmax = str2num(get(handles.maxVolEdit,'String'));
minLength = str2num(get(handles.minLengthEdit,'String'));

minX = str2num(get(handles.minXedit,'String'));
maxX = str2num(get(handles.maxXedit,'String'));
minY = str2num(get(handles.minYedit,'String'));
maxY = str2num(get(handles.maxYedit,'String'));
minZ = str2num(get(handles.minZedit,'String'));
maxZ = str2num(get(handles.maxZedit,'String'));

save(handles.data.loadPath,'inLimitsFlags','rejectFlags',...
    'minRatio','Vmin','Vmax','minLength',...
    'minX','maxX','minY','maxY','minZ','maxZ','-append');

% --- Executes on slider movement.
function frameSlider_Callback(hObject, eventdata, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderVal = get(hObject,'Value');
sliderVal = round(sliderVal);
set(hObject,'Value',sliderVal);

handles.data.activeFrame = get(hObject,'Value');

guidata(hObject,handles);

redraw(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function frameSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





function redraw(hObject, eventdata, handles)

handles = guidata(hObject);

showRejected = get(handles.showReject_checkbox,'Value');

if showRejected
    mapRejectFlags = ...
        ~handles.data.rejectFlags(handles.data.inLimitsInds);
else
    mapRejectFlags = ...
        handles.data.rejectFlags(handles.data.inLimitsInds);
end

ff = handles.data.activeFrame;

drawInds = find(cellfun(@(track)ismember(ff,track.timeInd),...
    handles.data.inLimitsTracks) & ~mapRejectFlags);
drawTracks = handles.data.inLimitsTracks(drawInds);
timePullInds = cellfun(@(track)find(track.timeInd==ff),drawTracks,...
    'UniformOutput',false);

yyCoords = cellfun(@(track,pullInd)track.centroid{pullInd}(1),...
    drawTracks,timePullInds);
xxCoords = cellfun(@(track,pullInd)track.centroid{pullInd}(2),...
    drawTracks,timePullInds);
zzCoords = cellfun(@(track,pullInd)track.centroid{pullInd}(3),...
    drawTracks,timePullInds);

handles.data.drawInds = drawInds;
handles.data.drawTracks = drawTracks;
handles.data.timePullInds = timePullInds;

handles.data.xxCoords = xxCoords;
handles.data.yyCoords = yyCoords;
handles.data.zzCoords = zzCoords;

guidata(hObject,handles);

% Saturate the brightest pixels of the image
intLim = prctile(handles.data.zmaxProj_cell{ff}(:),[1,99]);

axes(handles.axes_xy)
imagesc(handles.axes_xy,...
    [0,floor(handles.data.rawStackSizeY./handles.data.binning).*handles.data.voxelSizeY],...
    [0,floor(handles.data.rawStackSizeX./handles.data.binning).*handles.data.voxelSizeX],...
    handles.data.zmaxProj_cell{ff}',intLim);
hold on
if ~showRejected
    plot(handles.axes_xy,xxCoords,yyCoords,'ro')
else
    plot(handles.axes_xy,xxCoords,yyCoords,'bo')
end
set(get(handles.axes_xy,'Children'),'ButtonDownFcn',@axes_xy_ButtonDownFcn)
xlabel('x [\mum]')
ylabel('y [\mum]')
colormap(gray)
axis equal
axis tight
hold off

axes(handles.axes_xz)
imagesc(handles.axes_xz,[0,floor(handles.data.rawStackSizeY./handles.data.binning).*handles.data.voxelSizeY],...
    [0,floor(handles.data.rawStackSizeZ./handles.data.binning).*handles.data.voxelSizeZ],...
    handles.data.xmaxProj_cell{ff}.',intLim)
hold on
if ~showRejected
    plot(handles.axes_xz,xxCoords,zzCoords,'ro')
else
    plot(handles.axes_xz,xxCoords,zzCoords,'bo')
end
set(get(handles.axes_xz,'Children'),'ButtonDownFcn',@axes_xz_ButtonDownFcn)
xlabel('x [\mum]')
ylabel('z [\mum]')
colormap(gray)
axis equal
axis tight
set(gca,'YDir','normal')
hold off


axes(handles.axes_yz)
imagesc(handles.axes_yz,[0,floor(handles.data.rawStackSizeX./handles.data.binning).*handles.data.voxelSizeX],...
    [0,floor(handles.data.rawStackSizeZ./handles.data.binning).*handles.data.voxelSizeZ],...
    handles.data.ymaxProj_cell{ff}.',intLim)
hold on
if ~showRejected
    plot(handles.axes_yz,yyCoords,zzCoords,'ro')
else
    plot(handles.axes_yz,yyCoords,zzCoords,'bo')
end
set(get(handles.axes_yz,'Children'),'ButtonDownFcn',@axes_yz_ButtonDownFcn)
xlabel('y [\mum]')
ylabel('z [\mum]')
colormap(gray)
axis equal
axis tight
set(gca,'YDir','normal')
hold off



function updateLimiters(hObject, eventdata, handles)

handles.data.minRatio = str2num(get(handles.minRatioEdit,'String'));
handles.data.Vmin = str2num(get(handles.minVolEdit,'String'));
handles.data.Vmax = str2num(get(handles.maxVolEdit,'String'));
handles.data.minLength = str2num(get(handles.minLengthEdit,'String'));

handles.data.minX = str2num(get(handles.minXedit,'String'));
handles.data.maxX = str2num(get(handles.maxXedit,'String'));
handles.data.minY = str2num(get(handles.minYedit,'String'));
handles.data.maxY = str2num(get(handles.maxYedit,'String'));
handles.data.minZ = str2num(get(handles.minZedit,'String'));
handles.data.maxZ = str2num(get(handles.maxZedit,'String'));

intChannel = handles.data.intChannel;

% Still add x, y, and z limits to this

inLimFlags = ...
    cellfun(@(elmt)numel(elmt.time),handles.data.tracks)>=handles.data.minLength ...
    & cellfun(@(elmt)...
    max(elmt.nucIntensity{intChannel} ...
    ./elmt.cytoIntensity{intChannel}),handles.data.tracks) ...
    >=handles.data.minRatio ...
    & cellfun(@(elmt)max(elmt.volume),handles.data.tracks)<=handles.data.Vmax ...
    & cellfun(@(elmt)min(elmt.volume),handles.data.tracks)>=handles.data.Vmin ...
    & cellfun(@(elmt)min(cellfun(@(centroid)centroid(1),elmt.centroid)),...
    handles.data.tracks)>=handles.data.minY...
    & cellfun(@(elmt)max(cellfun(@(centroid)centroid(1),elmt.centroid)),...
    handles.data.tracks)<=handles.data.maxY...
    & cellfun(@(elmt)min(cellfun(@(centroid)centroid(2),elmt.centroid)),...
    handles.data.tracks)>=handles.data.minX...
    & cellfun(@(elmt)max(cellfun(@(centroid)centroid(2),elmt.centroid)),...
    handles.data.tracks)<=handles.data.maxX...
    & cellfun(@(elmt)min(cellfun(@(centroid)centroid(3),elmt.centroid)),...
    handles.data.tracks)>=handles.data.minZ...
    & cellfun(@(elmt)max(cellfun(@(centroid)centroid(3),elmt.centroid)),...
    handles.data.tracks)<=handles.data.maxZ;


handles.data.inLimitsFlags = inLimFlags;

handles.data.inLimitsInds = find(inLimFlags);

handles.data.inLimitsTracks = ...
    handles.data.tracks(handles.data.inLimitsInds);

guidata(hObject,handles);

redraw(hObject, eventdata, handles);



function minRatioEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minRatioEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minRatioEdit as text
%        str2double(get(hObject,'String')) returns contents of minRatioEdit as a double


% --- Executes during object creation, after setting all properties.
function minRatioEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minRatioEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minVolEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minVolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minVolEdit as text
%        str2double(get(hObject,'String')) returns contents of minVolEdit as a double


% --- Executes during object creation, after setting all properties.
function minVolEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minVolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxVolEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxVolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxVolEdit as text
%        str2double(get(hObject,'String')) returns contents of maxVolEdit as a double


% --- Executes during object creation, after setting all properties.
function maxVolEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxVolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minLengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minLengthEdit as text
%        str2double(get(hObject,'String')) returns contents of minLengthEdit as a double



% --- Executes during object creation, after setting all properties.
function minLengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes_xy_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_xy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
coordinates = get(handles.axes_xy,'CurrentPoint');

rejectTrack(hObject,coordinates,1);
redraw(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes_xz_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_xz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
coordinates = get(handles.axes_xz,'CurrentPoint');

rejectTrack(hObject,coordinates,2);
redraw(hObject, eventdata, handles);

% --- Executes on mouse press over axes background.
function axes_yz_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_yz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
coordinates = get(handles.axes_yz,'CurrentPoint');

rejectTrack(hObject,coordinates,3);
redraw(hObject, eventdata, handles);


function rejectTrack(hObject,coordinates,axisNum)

handles = guidata(hObject);

xxCoords = handles.data.xxCoords;
yyCoords = handles.data.yyCoords;
zzCoords = handles.data.zzCoords;

if axisNum == 1
    dists = ...
        sqrt((xxCoords - coordinates(1,1)).^2 ...
        + (yyCoords - coordinates(1,2)).^2);
elseif axisNum == 2
    dists = ...
        sqrt((xxCoords - coordinates(1,1)).^2 ...
        + (zzCoords - coordinates(1,2)).^2);
elseif axisNum == 3
    dists = ...
        sqrt((yyCoords - coordinates(1,1)).^2 ...
        + (zzCoords - coordinates(1,2)).^2);
end

% find closest coordinate on given axis
[~,minInd] = min(dists);

% add track to rejection set
rejectInd = handles.data.inLimitsInds(handles.data.drawInds(minInd));

handles.data.rejectFlags(rejectInd) = ...
    ~get(handles.showReject_checkbox,'Value');

guidata(hObject,handles);



function minXedit_Callback(hObject, eventdata, handles)
% hObject    handle to minXedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minXedit as text
%        str2double(get(hObject,'String')) returns contents of minXedit as a double



% --- Executes during object creation, after setting all properties.
function minXedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minXedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxXedit_Callback(hObject, eventdata, handles)
% hObject    handle to maxXedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxXedit as text
%        str2double(get(hObject,'String')) returns contents of maxXedit as a double


% --- Executes during object creation, after setting all properties.
function maxXedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxXedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minYedit_Callback(hObject, eventdata, handles)
% hObject    handle to minYedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minYedit as text
%        str2double(get(hObject,'String')) returns contents of minYedit as a double


% --- Executes during object creation, after setting all properties.
function minYedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minYedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxYedit_Callback(hObject, eventdata, handles)
% hObject    handle to maxYedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxYedit as text
%        str2double(get(hObject,'String')) returns contents of maxYedit as a double


% --- Executes during object creation, after setting all properties.
function maxYedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxYedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minZedit_Callback(hObject, eventdata, handles)
% hObject    handle to minZedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minZedit as text
%        str2double(get(hObject,'String')) returns contents of minZedit as a double


% --- Executes during object creation, after setting all properties.
function minZedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minZedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxZedit_Callback(hObject, eventdata, handles)
% hObject    handle to maxZedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxZedit as text
%        str2double(get(hObject,'String')) returns contents of maxZedit as a double


% --- Executes during object creation, after setting all properties.
function maxZedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxZedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.rejectFlags(:) = false;

redraw(hObject, eventdata, handles);

guidata(hObject,handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updateLimiters(hObject, eventdata, handles);
redraw(hObject, eventdata, handles);
