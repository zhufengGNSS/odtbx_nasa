function varargout = meas_sched_mock_kenny(varargin)
% MEAS_SCHED_MOCK_KENNY MATLAB code for meas_sched_mock_kenny.fig
%      MEAS_SCHED_MOCK_KENNY, by itself, creates a new MEAS_SCHED_MOCK_KENNY or raises the existing
%      singleton*.
%
%      H = MEAS_SCHED_MOCK_KENNY returns the handle to a new MEAS_SCHED_MOCK_KENNY or the handle to
%      the existing singleton*.
%
%      MEAS_SCHED_MOCK_KENNY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEAS_SCHED_MOCK_KENNY.M with the given input arguments.
%
%      MEAS_SCHED_MOCK_KENNY('Property','Value',...) creates a new MEAS_SCHED_MOCK_KENNY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before meas_sched_mock_kenny_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to meas_sched_mock_kenny_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help meas_sched_mock_kenny

% Last Modified by GUIDE v2.5 28-Sep-2012 09:22:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @meas_sched_mock_kenny_OpeningFcn, ...
                   'gui_OutputFcn',  @meas_sched_mock_kenny_OutputFcn, ...
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


% --- Executes just before meas_sched_mock_kenny is made visible.
function meas_sched_mock_kenny_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to meas_sched_mock_kenny (see VARARGIN)
global boxes;
clear global boxes; % Get rid of data from previous runs

global meas_add_remove;
clear global meas_add_remove;
meas_add_remove = 0;
set(handles.meas_schedule_mode,'SelectedObject',[handles.Add]);
% get(handles.meas_schedule_mode, 'SelectedObject')
% meas_add_remove

% Choose default command line output for meas_sched_mock_kenny
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Scroll bar code
handles.slide_handles = get(handles.ground_stations, 'Children');
% Work around for uicontrols not clipping correctly
handles.slide_labels = [handles.gs_label1, handles.gs_label2, ...
    handles.gs_label3, handles.gs_label4, handles.gs_label5, ...
    handles.gs_label6, handles.gs_label7, handles.gs_label8, ...
    handles.gs_label9, handles.gs_label10, handles.gs_label11, ...
    handles.gs_label12, handles.gs_label13, handles.gs_label14, ...
    handles.gs_label15, handles.gs_label16, handles.gs_label17, ...
    handles.gs_label18, handles.gs_label19, handles.gs_label20];

% Get original positions of objects
handles.slide_pos = get(handles.slide_handles, 'position');

% Update handle structure
guidata(hObject, handles);

% Initially, we only want the first six label buttons visible (work around
% for uicontrols not clipping correctly)
set(handles.slide_labels, 'Visible', 'off');
initial_slide_labels = [handles.gs_label1, handles.gs_label2, ...
    handles.gs_label3, handles.gs_label4, handles.gs_label5, ...
    handles.gs_label6];
set(initial_slide_labels, 'Visible', 'on');

% Make the axes uniform
% Link all the axes
linkaxes([handles.axes1, handles.axes2, handles.axes3, handles.axes4, handles.axes5...
    handles.axes6, handles.axes7, handles.axes8, handles.axes9, handles.axes10, ...
    handles.axes11, handles.axes12, handles.axes13, handles.axes14, handles.axes15, ...
    handles.axes16, handles.axes17, handles.axes18, handles.axes19, handles.axes20, ...
    handles.meas_total], 'xy');

% Set the size of the axes (will replace with dates)
set(handles.axes1,'XLim',[0 10]);
set(handles.axes1,'YLim',[0 1]);

% Set default add/remove measurements state to "add"
set(handles.meas_schedule_mode, 'SelectedObject', handles.Add);

% UIWAIT makes meas_sched_mock_kenny wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = meas_sched_mock_kenny_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Axes1
% % Visibility
% visa1 = rectangle('Position',[1,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% visa2 = rectangle('Position',[5,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% set(visa1, 'parent', handles.axes1);
% set(visa2, 'parent', handles.axes1);
% 

% % Axes2
% % Visibility
% visb1 = rectangle('Position',[3,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% set(visb1, 'parent', handles.axes2);
% 

% % Axes3
% % Visibility
% visc1 = rectangle('Position',[.5,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% visc2 = rectangle('Position',[3.5,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% visc3 = rectangle('Position',[7,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% set(visc1, 'parent', handles.axes3);
% set(visc2, 'parent', handles.axes3);
% set(visc3, 'parent', handles.axes3);
 
% % Axes4
% % Visibility
% visd1 = rectangle('Position',[2,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% set(visd1, 'parent', handles.axes4);


% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function edit_Callback(hObject, eventdata, handles)
% hObject    handle to edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function satellite_Callback(hObject, eventdata, handles)
% hObject    handle to satellite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% slide_max = get(hObject, 'Max')
% slide_min = get(hObject, 'Min')
slide_val = get(hObject, 'Value');
pos = cellfun(@(y) {[0 45 0 0]+y-[0 get(handles.slider1,'value') 0 0]}, handles.slide_pos);
set(handles.slide_handles, {'position'}, pos);

% Workaround to make uicontrols disappear
panel_pos = get(handles.ground_stations, 'position')

for i = 1:length(handles.slide_labels)
    button_pos = get(handles.slide_labels(i), 'position')
    if (((button_pos(2) + button_pos(4)) < panel_pos(4)) && (button_pos(2) > 0))
        set(handles.slide_labels(i), 'Visible', 'on');
    else
        set(handles.slide_labels(i), 'Visible', 'off');
    end
end
    

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Max',45);
set(hObject, 'Min', -40);
set(hObject, 'Value', get(hObject, 'Max'));



% Functions for creating boxes when the axes are clicked on

% --- Executes when selected object is changed in meas_schedule_mode.
function meas_schedule_mode_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in meas_schedule_mode 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global meas_add_remove;

mode = get(eventdata.NewValue, 'String');
if (strcmp(mode, 'Add'))
    meas_add_remove = 0; % Global variable to determine what mode the gui is in
elseif (strcmp(mode, 'Remove'))
    meas_add_remove = 1; % Global variable to determine what mode the gui is in
else
    meas_add_remove = -1; % Global variable to determine what mode the gui is in
end
% get(handles.meas_schedule_mode, 'SelectedObject')
meas_add_remove;


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes7_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes8_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes9_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes10_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes11_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes12_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes13_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes14_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes15_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes16_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes17_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes18_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes19_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


% --- Executes on mouse press over axes background.
function axes20_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
schedule_measurements(hObject, eventdata, handles);


function schedule_measurements(hObject, eventdata, handles)

persistent add_coords;
global meas_add_remove;
global boxes;

if (isempty(meas_add_remove))
    meas_add_remove = 0;
end
% meas_add_remove

if (isempty(boxes))
    boxes = struct('ground_station', [], ...
        'type', [], ...
        'x', [0, 0], ...
        'handle', [], ...
        'total_handle', []);
end

mousepos = get(hObject, 'currentpoint');
% screenpos = get(handles.axes1, 'position')

if (meas_add_remove == 0) % Add boxes to axes
    if isempty(add_coords) % First mouse click
        add_coords(1,1:2) = mousepos(1,1:2);
    else % Second mouse click
        add_coords(2,1:2) = mousepos(1,1:2);
        boxes = create_a_box(add_coords, [hObject, handles.meas_total], boxes); 
        clear add_coords;
    end
elseif (meas_add_remove == 1) % Remove boxes from axes
    clear add_cords;
    remove_coords(1,1) = mousepos(1,1);
    boxes = remove_a_box(remove_coords, [hObject, handles.meas_total], boxes);
else
    clear add_cords;
end

% Debugging
assignin('base', 'boxes', boxes(:));


function [boxes] = create_a_box(coords, axes_handles, boxes)

if (coords(2,1) < coords(1,1)) % Rectangles can only have positive deltas
    coords_temp = coords(1,1);
    coords(1,1) = coords(2,1);
    coords(2,1) = coords_temp;
end
x = coords(1,1);
% y = coords(1,2);
dx = coords(2,1) - coords(1,1);
% dy = coords(2,2) - coords(1,2);

% Plot a box on parent axes
meas = rectangle('Position',[x,0,dx,1], 'EdgeColor','g','LineWidth', 5);%,'FaceColor',[175/255 1 175/255]);
set(meas, 'parent', axes_handles(1));

% Plot a box on total axes
meas_on_total = rectangle('Position',[x,0,dx,1], 'EdgeColor','g','LineWidth', 5);%,'FaceColor','g');
set(meas_on_total, 'parent', axes_handles(2));

% Save a box in memory
% for i = 1:boxes(end+1)

boxes(end+1) = struct('ground_station', get(axes_handles(1), 'Tag'), ...
    'type', 'measurement', ...
    'x', [coords(1,1) coords(2,1)], ...
    'handle', meas, ...
    'total_handle', meas_on_total);


function [boxes] = remove_a_box(coords, axes_handles, boxes)

i = 1;
while (i <= length(boxes)) % While loop, *not* for loop (we need length recalculated every iteration)
    if (isequal(axes_handles(1), get(boxes(i).handle, 'parent')))
        if ((boxes(i).x(1) <= coords(1,1)) && (coords(1,1) <= boxes(i).x(2)))

            % Delete boxes
            delete(boxes(i).handle);
            delete(boxes(i).total_handle);
            
            % Delete the record
            boxes(i) = [];
        end
    end
    i = i + 1;
end