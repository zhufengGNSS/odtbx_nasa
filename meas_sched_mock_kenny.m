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

% Last Modified by GUIDE v2.5 17-Sep-2012 16:04:13

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

% Choose default command line output for meas_sched_mock_kenny
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Link all the axes
linkaxes([handles.axes1, handles.axes2, handles.axes3, handles.axes4], 'y')
linkaxes([handles.axes1, handles.axes2, handles.axes3, handles.axes4, handles.axes5], 'x')

% Set the size of the axes (will replace with dates)
set(handles.axes1,'XLim',[0 10])
set(handles.axes1,'YLim',[0 1])

% UIWAIT makes meas_sched_mock_kenny wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = meas_sched_mock_kenny_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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
% % Measurements
% measa1 = rectangle('Position',[1.5,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% measa2 = rectangle('Position',[2.5,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% measa3 = rectangle('Position',[5.5,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% measa4 = rectangle('Position',[6.5,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% set(measa1, 'parent', handles.axes1);
% set(measa2, 'parent', handles.axes1);
% set(measa3, 'parent', handles.axes1);
% set(measa4, 'parent', handles.axes1);
% 
% % Axes2
% % Visibility
% visb1 = rectangle('Position',[3,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% set(visb1, 'parent', handles.axes2);
% 
% % Measurements
% measb1 = rectangle('Position',[3.25,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% measb2 = rectangle('Position',[4.5,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% set(measb1, 'parent', handles.axes2);
% set(measb2, 'parent', handles.axes2);
% 
% % Axes3
% % Visibility
% visc1 = rectangle('Position',[.5,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% visc2 = rectangle('Position',[3.5,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% visc3 = rectangle('Position',[7,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% set(visc1, 'parent', handles.axes3);
% set(visc2, 'parent', handles.axes3);
% set(visc3, 'parent', handles.axes3);
% 
% % Measurements
% measc1 = rectangle('Position',[1,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% measc2 = rectangle('Position',[2,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% measc3 = rectangle('Position',[3.75,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% measc4 = rectangle('Position',[5.25,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% measc5 = rectangle('Position',[8.25,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% set(measc1, 'parent', handles.axes3);
% set(measc2, 'parent', handles.axes3);
% set(measc3, 'parent', handles.axes3);
% set(measc4, 'parent', handles.axes3);
% set(measc5, 'parent', handles.axes3);
% 
% 
% % Axes4
% % Visibility
% visd1 = rectangle('Position',[2,0,2.5,1], 'EdgeColor','r','LineWidth', 5,'FaceColor','w')
% set(visd1, 'parent', handles.axes4);
% 
% % Measurements
% measd1 = rectangle('Position',[3,0,.5,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w')
% set(measd1, 'parent', handles.axes4);


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
function ground_station_Callback(hObject, eventdata, handles)
% hObject    handle to ground_station (see GCBO)
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


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
persistent coords;

mousepos = get(handles.axes1, 'currentpoint')
% screenpos = get(handles.axes1, 'position')
if isempty(coords)
    coords(1,1:2) = mousepos(1,1:2);
else
    coords(2,1:2) = mousepos(1,1:2);
    create_a_box(coords, handles.axes1);
    clear coords
end

function create_a_box(coords, axes_handle)

if (coords(2,1) < coords(1,1)) % Rectangles can only have positive deltas
    coords_temp = coords(1,1);
    coords(1,1) = coords(2,1);
    coords(2,1) = coords_temp;
end
x = coords(1,1);
% y = coords(1,2);
dx = coords(2,1) - coords(1,1);
% dy = coords(2,2) - coords(1,2);

meas = rectangle('Position',[x,0,dx,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','w');
set(meas, 'parent', axes_handle);
