% PATTERN_ADD This is a function of MEAS_SCHED used to create repeating
% measurements.
%
% (This file is part of ODTBX, The Orbit Determination Toolbox, and is
%  distributed under the NASA Open Source Agreement.  See file source for
%  more details.)

% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2011 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.


function varargout = pattern_add(varargin)
% PATTERN_ADD MATLAB code for pattern_add.fig
%      PATTERN_ADD, by itself, creates a new PATTERN_ADD or raises the existing
%      singleton*.
%
%      H = PATTERN_ADD returns the handle to a new PATTERN_ADD or the handle to
%      the existing singleton*.
%
%      PATTERN_ADD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PATTERN_ADD.M with the given input arguments.
%
%      PATTERN_ADD('Property','Value',...) creates a new PATTERN_ADD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pattern_add_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pattern_add_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pattern_add

% Last Modified by GUIDE v2.5 23-Oct-2012 11:37:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pattern_add_OpeningFcn, ...
                   'gui_OutputFcn',  @pattern_add_OutputFcn, ...
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


% --- Executes just before pattern_add is made visible.
function pattern_add_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pattern_add (see VARARGIN)

% Choose default command line output for pattern_add
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

dontOpen = false;
mainGuiInput = find(strcmp(varargin, 'meas_sched'));
if (isempty(mainGuiInput)) ...
    || (length(varargin) <= mainGuiInput) ...
    || (~ishandle(varargin{mainGuiInput+1}))
    dontOpen = true;
else
    % Remember the handle, and adjust our position
    handles.meas_sched_Main = varargin{mainGuiInput+1};
    handles.time_in = varargin{mainGuiInput+3};
    
    % Obtain handles using GUIDATA with the caller's handle 
    mainHandles = guidata(handles.meas_sched_Main);
    
    set(handles.date_begin_string, 'String', datestr(handles.time_in(1), 'mm/dd/yyyy HH:MM:SS'));
    set(handles.date_final_string, 'String', datestr(handles.time_in(2), 'mm/dd/yyyy HH:MM:SS'));
    
    visible_time = get(mainHandles.meas_total, 'XLim');
    set(handles.repeat_until_string, 'String', datestr(visible_time(2), 'mm/dd/yyyy HH:MM:SS'));

    % Determine the position of the dialog - centered on the callback figure
    % if available, else, centered on the screen
    FigPos=get(0,'DefaultFigurePosition');
    OldUnits = get(hObject, 'Units');
    set(hObject, 'Units', 'pixels');
    OldPos = get(hObject,'Position');
    FigWidth = OldPos(3);
    FigHeight = OldPos(4);
    if isempty(gcbf)
        ScreenUnits=get(0,'Units');
        set(0,'Units','pixels');
        ScreenSize=get(0,'ScreenSize');
        set(0,'Units',ScreenUnits);

        FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
        FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
    else
        GCBFOldUnits = get(gcbf,'Units');
        set(gcbf,'Units','pixels');
        GCBFPos = get(gcbf,'Position');
        set(gcbf,'Units',GCBFOldUnits);
        FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                       (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
    end
    FigPos(3:4)=[FigWidth FigHeight];
    set(hObject, 'Position', FigPos);
    set(hObject, 'Units', OldUnits);

    % Make the GUI modal
    set(handles.figure1,'WindowStyle','modal')
end

% Update handles structure
guidata(hObject, handles);


if dontOpen
   disp('-----------------------------------------------------');
   disp('Improper input arguments. Pass a property value pair') 
   disp('whose name is "changeme_main" and value is the handle')
   disp('to the changeme_main figure, e.g:');
   disp('   x = changeme_main()');
   disp('   changeme_dialog(''changeme_main'', x)');
   disp('-----------------------------------------------------');
else    
% UIWAIT makes pattern_add wait for user response (see UIRESUME)
    uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = pattern_add_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.date_begin_string;
varargout{3} = handles.date_final_string;
varargout{4} = handles.repeat_freq;
varargout{5} = handles.repeat_until;

% The figure can be deleted now
delete(handles.figure1);



function date_begin_string_Callback(hObject, eventdata, handles)
% hObject    handle to date_begin_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of date_begin_string as text
%        str2double(get(hObject,'String')) returns contents of date_begin_string as a double


% --- Executes during object creation, after setting all properties.
function date_begin_string_CreateFcn(hObject, eventdata, handles)
% hObject    handle to date_begin_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function date_final_string_Callback(hObject, eventdata, handles)
% hObject    handle to date_final_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of date_final_string as text
%        str2double(get(hObject,'String')) returns contents of date_final_string as a double


% --- Executes during object creation, after setting all properties.
function date_final_string_CreateFcn(hObject, eventdata, handles)
% hObject    handle to date_final_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(hObject,'String');
handles.date_begin_string = 0;
handles.date_final_string = 0;
handles.repeat_freq = 0;
handles.repeat_until = 0;

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes on button press in ok_button.
function ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Bring in the changed values
handles.date_begin_string = datenum(get(handles.date_begin_string, 'String'));
handles.date_final_string = datenum(get(handles.date_final_string, 'String'));

repeat_every = str2num(get(handles.repeat_every_string, 'String'));

time_unit = get(handles.units, 'Value');
switch time_unit
    case 1 % Minutes
        time_multiplier = 1/24*1/60;
    case 2 % Hours
        time_multiplier = 1/24;
    case 3 % Days
        time_multiplier = 1;
    case 4 % Weeks
        time_multiplier = 7;
end

handles.repeat_freq = repeat_every * time_multiplier;
handles.repeat_until = datenum(get(handles.repeat_until_string, 'String'));

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = 'Cancel';
handles.date_begin_string = 0;
handles.date_final_string = 0;
handles.repeat_freq = 0;
handles.repeat_until = 0;

guidata(hObject, handles)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end



function repeat_every_string_Callback(hObject, eventdata, handles)
% hObject    handle to repeat_every_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of repeat_every_string as text
%        str2double(get(hObject,'String')) returns contents of repeat_every_string as a double


% --- Executes during object creation, after setting all properties.
function repeat_every_string_CreateFcn(hObject, eventdata, handles)
% hObject    handle to repeat_every_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in units.
function units_Callback(hObject, eventdata, handles)
% hObject    handle to units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns units contents as cell array
%        contents{get(hObject,'Value')} returns selected item from units


% --- Executes during object creation, after setting all properties.
function units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function repeat_until_string_Callback(hObject, eventdata, handles)
% hObject    handle to repeat_until_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of repeat_until_string as text
%        str2double(get(hObject,'String')) returns contents of repeat_until_string as a double


% --- Executes during object creation, after setting all properties.
function repeat_until_string_CreateFcn(hObject, eventdata, handles)
% hObject    handle to repeat_until_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
