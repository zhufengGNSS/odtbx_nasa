function varargout = satellite_edit(varargin)
    % SATELLITE_EDIT MATLAB code for satellite_edit.fig
    %      SATELLITE_EDIT, by itself, creates a new SATELLITE_EDIT or raises the existing
    %      singleton*.
    %
    %      H = SATELLITE_EDIT returns the handle to a new SATELLITE_EDIT or the handle to
    %      the existing singleton*.
    %
    %      SATELLITE_EDIT('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in SATELLITE_EDIT.M with the given input arguments.
    %
    %      SATELLITE_EDIT('Property','Value',...) creates a new SATELLITE_EDIT or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before satellite_edit_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to satellite_edit_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help satellite_edit

    % Last Modified by GUIDE v2.5 01-Nov-2012 12:33:32

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @satellite_edit_OpeningFcn, ...
                       'gui_OutputFcn',  @satellite_edit_OutputFcn, ...
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

% --- Executes just before satellite_edit is made visible.
function satellite_edit_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to satellite_edit (see VARARGIN)
    global time_sim;
    global sat_state_sim;
    global propagator;

    set(handles.x_pos, 'String', num2str(sat_state_sim.pos_x));
    set(handles.y_pos, 'String', num2str(sat_state_sim.pos_y));
    set(handles.z_pos, 'String', num2str(sat_state_sim.pos_z));
    set(handles.x_vel, 'String', num2str(sat_state_sim.vel_x));
    set(handles.y_vel, 'String', num2str(sat_state_sim.vel_y));
    set(handles.z_vel, 'String', num2str(sat_state_sim.vel_z));

    set(handles.time_begin, 'String', num2str(time_sim.begin));
    set(handles.time_inc, 'String', num2str(time_sim.increment));
    set(handles.time_end, 'String', num2str(time_sim.end));

    set(handles.prop_handle, 'String', propagator);


    % Choose default command line output for satellite_edit
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

        % Obtain handles using GUIDATA with the caller's handle 
        mainHandles = guidata(handles.meas_sched_Main);

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
       disp('Improper input arguments: satellite_edit.m ') 
       disp('-----------------------------------------------------');
    else    
    % UIWAIT makes satellite_edit wait for user response (see UIRESUME)
        uiwait(handles.figure1);
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = satellite_edit_OutputFcn(hObject, eventdata, handles) 
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

    global time_sim;
    global sat_state_sim;
    global propagator;

    state_values = struct('X_Position', str2num(get(handles.x_pos, 'String')), ...
                          'Y_Position', str2num(get(handles.y_pos, 'String')), ...
                          'Z_Position', str2num(get(handles.z_pos, 'String')), ...
                          'X_Velocity', str2num(get(handles.x_vel, 'String')), ...
                          'Y_Velocity', str2num(get(handles.y_vel, 'String')), ...
                          'Z_Velocity', str2num(get(handles.z_vel, 'String')));
    
    time_values = struct('Begin_Time', str2num(get(handles.time_begin, 'String')), ...
                         'Time_Increment', str2num(get(handles.time_inc, 'String')), ...
                         'End_Time', str2num(get(handles.time_end, 'String')));
    
    prop_values = struct('Propagator', get(handles.prop_handle, 'String'));
        
    error = check_entry_validity(state_values, time_values, prop_values);
    
    if (~error)
        % Assign the values to be passed out
        sat_state_sim.pos_x = state_values.X_Position;
        sat_state_sim.pos_y = state_values.Y_Position;
        sat_state_sim.pos_z = state_values.Z_Position;
        sat_state_sim.vel_x = state_values.X_Velocity;
        sat_state_sim.vel_y = state_values.Y_Velocity;
        sat_state_sim.vel_z = state_values.Z_Velocity;

        time_sim.begin = time_values.Begin_Time;
        time_sim.increment = time_values.Time_Increment;
        time_sim.end = time_values.End_Time;

        propagator = prop_values.Propagator;


        handles.output = get(hObject,'String');

        % Update handles structure
        guidata(hObject, handles);

        % Use UIRESUME instead of delete because the OutputFcn needs
        % to get the updated handles structure.
        uiresume(handles.figure1);
    end
end


function [error] = check_entry_validity(state_values, time_values, prop_values)
    % This function keeps the GUI from exiting with bad value    

    error = 0;
    
    % Check the state for non-numbers
    state_field_names = fieldnames(state_values);
    for i = 1:length(state_field_names)
        if (isempty(state_values.(state_field_names{i})))
            errordlg('Must be a number!', state_field_names{i});
            error = 1;
            return;
        end
    end
    
    % Check the time interval for non-numbers
    time_field_names = fieldnames(time_values);
    for i = 1:length(time_field_names)
        if (isempty(time_values.(time_field_names{i})))
            errordlg('Must be a number!', time_field_names{i});
            error = 1;
            return;
        end
    end
    
    % Make sure propagator exists
    % The user has to actually enter something
    prop_field_names = fieldnames(prop_values);
    try
        % Warnings appear informing us of the next error check (so we don't
        % really want to see them. In the future this will throw an error,
        % not a warning.
        warning off
        function_handle = str2func(prop_values.(prop_field_names{1}));
        warning on
    catch exception
        errordlg('Is not a valid function!', prop_field_names{1});
        error = 1;
        return;
    end
    % And whatever they enter should be a defined function
    info = functions(function_handle);
    if (isempty(info.file))
        errordlg('Function does not exist!', prop_field_names{1});
        error = 1;
        return
    end
    
end


% --- Executes on button press in check_range.
function check_range_Callback(hObject, eventdata, handles)
% hObject    handle to check_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_range
end


% --- Executes on button press in check_ea.
function check_ea_Callback(hObject, eventdata, handles)
% hObject    handle to check_ea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_ea
end


% --- Executes on button press in check_visibility.
function check_visibility_Callback(hObject, eventdata, handles)
% hObject    handle to check_visibility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_visibility
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


function x_pos_Callback(hObject, eventdata, handles)
    % hObject    handle to x_pos (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of x_pos as text
    %        str2double(get(hObject,'String')) returns contents of x_pos as a double
end


% --- Executes during object creation, after setting all properties.
function x_pos_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to x_pos (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function y_pos_Callback(hObject, eventdata, handles)
    % hObject    handle to y_pos (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of y_pos as text
    %        str2double(get(hObject,'String')) returns contents of y_pos as a double
end


% --- Executes during object creation, after setting all properties.
function y_pos_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to y_pos (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function z_pos_Callback(hObject, eventdata, handles)
    % hObject    handle to z_pos (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of z_pos as text
    %        str2double(get(hObject,'String')) returns contents of z_pos as a double
end

% --- Executes during object creation, after setting all properties.
function z_pos_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to z_pos (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function x_vel_Callback(hObject, eventdata, handles)
    % hObject    handle to x_vel (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of x_vel as text
    %        str2double(get(hObject,'String')) returns contents of x_vel as a double
end


% --- Executes during object creation, after setting all properties.
function x_vel_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to x_vel (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function y_vel_Callback(hObject, eventdata, handles)
    % hObject    handle to y_vel (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of y_vel as text
    %        str2double(get(hObject,'String')) returns contents of y_vel as a double

end


% --- Executes during object creation, after setting all properties.
function y_vel_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to y_vel (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function z_vel_Callback(hObject, eventdata, handles)
    % hObject    handle to z_vel (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of z_vel as text
    %        str2double(get(hObject,'String')) returns contents of z_vel as a double
end


% --- Executes during object creation, after setting all properties.
function z_vel_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to z_vel (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function time_begin_Callback(hObject, eventdata, handles)
    % hObject    handle to time_begin (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of time_begin as text
    %        str2double(get(hObject,'String')) returns contents of time_begin as a double
end


% --- Executes during object creation, after setting all properties.
function time_begin_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to time_begin (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function time_inc_Callback(hObject, eventdata, handles)
    % hObject    handle to time_inc (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of time_inc as text
    %        str2double(get(hObject,'String')) returns contents of time_inc as a double
end


% --- Executes during object creation, after setting all properties.
function time_inc_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to time_inc (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function time_end_Callback(hObject, eventdata, handles)
    % hObject    handle to time_end (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of time_end as text
    %        str2double(get(hObject,'String')) returns contents of time_end as a double
end


% --- Executes during object creation, after setting all properties.
function time_end_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to time_end (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end



function edit12_Callback(hObject, eventdata, handles)
    % hObject    handle to edit12 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit12 as text
    %        str2double(get(hObject,'String')) returns contents of edit12 as a double
end


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit12 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function edit13_Callback(hObject, eventdata, handles)
    % hObject    handle to edit13 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit13 as text
    %        str2double(get(hObject,'String')) returns contents of edit13 as a double
end


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit13 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function edit14_Callback(hObject, eventdata, handles)
    % hObject    handle to edit14 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of edit14 as text
    %        str2double(get(hObject,'String')) returns contents of edit14 as a double
end


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to edit14 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function prop_handle_Callback(hObject, eventdata, handles)
    % hObject    handle to prop_handle (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of prop_handle as text
    %        str2double(get(hObject,'String')) returns contents of prop_handle as a double
end


% --- Executes during object creation, after setting all properties.
function prop_handle_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to prop_handle (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end
