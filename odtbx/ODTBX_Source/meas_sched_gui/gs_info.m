% GS_INFO This is a function of MEAS_SCHED used to select ground stations.
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


function varargout = gs_info(varargin)
% GS_INFO MATLAB code for gs_info.fig
%      GS_INFO, by itself, creates a new GS_INFO or raises the existing
%      singleton*.
%
%      H = GS_INFO returns the handle to a new GS_INFO or the handle to
%      the existing singleton*.
%
%      GS_INFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GS_INFO.M with the given input arguments.
%
%      GS_INFO('Property','Value',...) creates a new GS_INFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gs_info_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gs_info_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gs_info

% Last Modified by GUIDE v2.5 01-Nov-2012 17:07:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gs_info_OpeningFcn, ...
                   'gui_OutputFcn',  @gs_info_OutputFcn, ...
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
end


% --- Executes just before gs_info is made visible.
function gs_info_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gs_info (see VARARGIN)

% Choose default command line output for gs_info
handles.output = hObject;

% So that we don't have to rebuild this over and over again
global measOptions;
handles.gsList = createGroundStationList();
measOptions = setOdtbxOptions(measOptions, 'gsList', handles.gsList);

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
    
    % Obtain handles using GUIDATA with the caller's handle 
    mainHandles = guidata(handles.meas_sched_Main);
    
    % Import the data from main function
    set(handles.gs_name, 'String', ...
        get(mainHandles.gs_label_current, 'String'));
    gsECEF_local = get(mainHandles.gs_label_current, 'UserData');
    if (~isempty(gsECEF_local))
        set(handles.gs_pos_x, 'String', num2str(gsECEF_local(1)));
        set(handles.gs_pos_y, 'String', num2str(gsECEF_local(2)));
        set(handles.gs_pos_z, 'String', num2str(gsECEF_local(3)));
    end
    
    set(mainHandles.gs_label_current, 'UserData', ...
        [str2num(get(handles.gs_pos_x, 'String')); ... 
        str2num(get(handles.gs_pos_y, 'String')); ...
        str2num(get(handles.gs_pos_z, 'String'))]);

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
   disp('Improper input arguments: gs_info.m ') 
   disp('-----------------------------------------------------');
else    
% UIWAIT makes gs_info wait for user response (see UIRESUME)
    uiwait(handles.figure1);
end

end


% --- Outputs from this function are returned to the command line.
function varargout = gs_info_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);
end



function gs_name_Callback(hObject, eventdata, handles)
% hObject    handle to gs_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gs_name as text
%        str2double(get(hObject,'String')) returns contents of gs_name as a double
end


% --- Executes during object creation, after setting all properties.
function gs_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gs_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);
end



% --- Executes on button press in ok_button.
function ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set the ground station text value on the main GUI
main = handles.meas_sched_Main;
% Obtain handles using GUIDATA with the caller's handle 
if(ishandle(main))
    mainHandles = guidata(main);
    % Change the visible label
    set(mainHandles.gs_label_current, 'String', ...
        get(handles.gs_name, 'String'));
    set(mainHandles.gs_label_current, 'UserData', ...
        [str2num(get(handles.gs_pos_x, 'String')); ... 
        str2num(get(handles.gs_pos_y, 'String')); ...
        str2num(get(handles.gs_pos_z, 'String'))]);
%     set(mainHandles.gs_axes_current, 'UserData', ...
%         get(handles.gs_name, 'String'));
end

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

end


function gs_pos_x_Callback(hObject, eventdata, handles)
% hObject    handle to gs_pos_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gs_pos_x as text
%        str2double(get(hObject,'String')) returns contents of gs_pos_x as a double

end


% --- Executes during object creation, after setting all properties.
function gs_pos_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gs_pos_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end



function gs_pos_y_Callback(hObject, eventdata, handles)
% hObject    handle to gs_pos_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gs_pos_y as text
%        str2double(get(hObject,'String')) returns contents of gs_pos_y as a double

end


% --- Executes during object creation, after setting all properties.
function gs_pos_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gs_pos_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end



function gs_pos_z_Callback(hObject, eventdata, handles)
% hObject    handle to gs_pos_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gs_pos_z as text
%        str2double(get(hObject,'String')) returns contents of gs_pos_z as a double

end


% --- Executes during object creation, after setting all properties.
function gs_pos_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gs_pos_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on selection change in gs_name_list.
function gs_name_list_Callback(hObject, eventdata, handles)
% hObject    handle to gs_name_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gs_name_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gs_name_list

end


% --- Executes during object creation, after setting all properties.
function gs_name_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gs_name_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in find_gs_button.
function find_gs_button_Callback(hObject, eventdata, handles)
% hObject    handle to find_gs_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Test epoch - DELETE THIS VALUE
epoch = datenum('Jan 1 2006');

% Take the value and find the corresponding ECEF coordinates
% gsList = createGroundStationList();
gsID = { char(get(handles.gs_name_list, 'String')) };
gsECEF = zeros(3,1);

try
    gsECEF(:,1) = getGroundStationInfo(handles.gsList,gsID{1},'ecefPosition',epoch);
catch exception
    % Warn people if their ground station wasn't in the list
    errordlg(exception.message, 'Not Found!');
end

% Display the name in the name box
set(handles.gs_name, 'String', gsID{1});

% Display the coordinates in the ECEF boxes
set(handles.gs_pos_x, 'String', gsECEF(1,1));
set(handles.gs_pos_y, 'String', gsECEF(2,1));
set(handles.gs_pos_z, 'String', gsECEF(3,1));

end


% --- Executes during object creation, after setting all properties.
function find_gs_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to find_gs_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end
