% MEAS_SCHED This is a graphical tool used for scheduling ground station
% measurements
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

function varargout = meas_sched(varargin)
    % MEAS_SCHED MATLAB code for meas_sched.fig
    %      MEAS_SCHED, by itself, creates a new MEAS_SCHED or raises the existing
    %      singleton*.
    %
    %      H = MEAS_SCHED returns the handle to a new MEAS_SCHED or the handle to
    %      the existing singleton*.
    %
    %      MEAS_SCHED('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in MEAS_SCHED.M with the given input arguments.
    %
    %      MEAS_SCHED('Property','Value',...) creates a new MEAS_SCHED or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before meas_sched_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to meas_sched_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to calculate (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % visualization the above text to modify the response to help meas_sched

    % Last Modified by GUIDE v2.5 20-Nov-2012 22:31:51

    % Begin initialization code - DO NOT VISUALIZATION
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @meas_sched_OpeningFcn, ...
                       'gui_OutputFcn',  @meas_sched_OutputFcn, ...
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
    % End initialization code - DO NOT VISUALIZATION
end



% --- Executes just before meas_sched is made visible.
function meas_sched_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to meas_sched (see VARARGIN)

    global meas_add_remove;
    meas_add_remove = 1;
    set(handles.meas_schedule_mode,'SelectedObject',[handles.Add]);

    global measOptions;
    measOptions = odtbxOptions('measurement');
    
    global boxes;
    if (~isempty(boxes)) % Get rid of data from previous runs
    boxes = struct('ground_station', [], ...
        'type', [], ...
        'x', [0, 0], ...
        'handle', [], ...
        'total_handle', []);
    end
    
%     global time_prop;
%     clear global time_sim;
    
%     global sat_state_prop;
%     clear global sat_state_sim;
    
%     global propagator;
%     clear global propagator;
    
    % Choose default command line output for meas_sched
    handles.output = hObject;

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

    % Handle to the current gs_label button
    handles.gs_label_current = handles.gs_label1;

    % Initially, we only want the first six label buttons visible (work around
    % for uicontrols not clipping correctly)
    set(handles.slide_labels, 'Visible', 'off');
    initial_slide_labels = [handles.gs_label1, handles.gs_label2, ...
        handles.gs_label3, handles.gs_label4, handles.gs_label5, ...
        handles.gs_label6];
    set(initial_slide_labels, 'Visible', 'on');

    % Make the axes uniform
    % Link all the axes
    handles.axes_handles = [handles.axes1, handles.axes2, handles.axes3, handles.axes4, handles.axes5...
        handles.axes6, handles.axes7, handles.axes8, handles.axes9, handles.axes10, ...
        handles.axes11, handles.axes12, handles.axes13, handles.axes14, handles.axes15, ...
        handles.axes16, handles.axes17, handles.axes18, handles.axes19, handles.axes20, ...
        handles.meas_total];
    linkaxes(handles.axes_handles, 'xy');
    
    % Give the axes a default name
    set(handles.axes_handles, 'UserData', 0);

    % Set the size of the axes (will replace with dates)
    set(handles.axes_handles,'YLim',[0 1]);
    % Select a starting date:
    startDate = datenum('07-03-2012');
    % Select an ending date:
    endDate = datenum('07-27-2012');
    % Create xdata to correspond to the number of 
    % units between the start and end dates:
    set(handles.axes_handles, 'XLim', [startDate endDate])
    xData = linspace(startDate,endDate,7);
    % Set the number of XTicks to the number of points in xData:
    set(handles.axes_handles,'XTick',xData)
    datetick('x', 'mm/dd/yy', 'keepticks');

    % Update handle structure
    guidata(hObject, handles);

    % Set default add/remove measurements state to "add"
    set(handles.meas_schedule_mode, 'SelectedObject', handles.Add);

    % UIWAIT makes meas_sched wait for user response (see UIRESUME)
%     uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = meas_sched_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    % varargout{1} = handles.output;
end


% --- Executes on button press in export_button.
function export_button_Callback(hObject, eventdata, handles)
    % hObject    handle to export_button (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    export_options_to_workspace();
end


% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
    % hObject    handle to file (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function simulation_Callback(hObject, eventdata, handles)
% hObject    handle to simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end



% --------------------------------------------------------------------
function visualization_Callback(~, eventdata, handles)
    % hObject    handle to visualization (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
    % hObject    handle to help (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
    % hObject    handle to about (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    about('meas_sched', handles.figure1);
end


% --------------------------------------------------------------------
function satellite_Callback(hObject, eventdata, handles)
    % hObject    handle to satellite (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    change_satellite(hObject, eventdata, handles);
end


% --------------------------------------------------------------------
function options_Callback(hObject, eventdata, handles)
% hObject    handle to options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    gs_options('meas_sched', handles.figure1);
end


% --------------------------------------------------------------------
function export_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to export_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    export_options_to_workspace()
end


% --------------------------------------------------------------------
function export_csv_Callback(hObject, eventdata, handles)
    % hObject    handle to export_csv (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    export_schedule_to_csv();
end


% --------------------------------------------------------------------
function import_options_Callback(hObject, eventdata, handles)
% hObject    handle to import_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    import_options_from_workspace(hObject, eventdata, handles);
end


% --------------------------------------------------------------------
function quit_Callback(hObject, eventdata, handles)
    % hObject    handle to quit (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    clear_data();
    % The figure can be deleted now
    delete(handles.figure1);
end


% --------------------------------------------------------------------
function time_Callback(hObject, eventdata, handles)
    % hObject    handle to time (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    change_time(hObject, eventdata, handles);
end


% --------------------------------------------------------------------
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function generate_Callback(hObject, eventdata, handles)
% hObject    handle to generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global measOptions;
    global T;
    
    if (~isempty(T))
        make_meas();

        % Adjust the time range we can see on the plots to be the same as the
        % time span over which the orbit was propagated

        % Pull in the current display data
        xData = get(handles.meas_total, 'Xtick');
        increments = length(xData);

        % Calculate the new time range
        new_begin = getOdtbxOptions(measOptions, 'epoch');
        new_end = new_begin + T*1/60*1/60*1/24;
        xData = linspace(new_begin,new_end(end),increments);

        % Set the new limits on all of the axes
        set(handles.axes_handles,'XLim',[xData(1) xData(end)]);

        % Set the number of XTicks to the number of points in xData:
        set(handles.axes_handles,'XTick',xData);

        datetick('x', 'mm/dd/yy', 'keepticks');

        % Redraw all the boxes on the potentially new axes scale
        plot_meas(hObject, eventdata, handles);
        redraw_boxes(hObject, eventdata, handles);
    end
end


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
    panel_pos = get(handles.ground_stations, 'position');

    for i = 1:length(handles.slide_labels)
        button_pos = get(handles.slide_labels(i), 'position');
        if (((button_pos(2) + button_pos(4)) < panel_pos(4)) && (button_pos(2) > 0))
            set(handles.slide_labels(i), 'Visible', 'on');
        else
            set(handles.slide_labels(i), 'Visible', 'off');
        end
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
end


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
        meas_add_remove = 1; % Global variable to determine what mode the gui is in
    elseif (strcmp(mode, 'Remove'))
        meas_add_remove = -1; % Global variable to determine what mode the gui is in
    elseif (strcmp(mode, 'Pattern Add'))
        meas_add_remove = 2; % Global variable to determine what mode the gui is in
    elseif (strcmp(mode, 'Edit'))
        meas_add_remove = 0; % Global variable to determine what mode the gui is in
    end
end



% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes2 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes3 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes4 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes5_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes5 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes6_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes6 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes7_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes7 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes8_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes8 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes9_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes9 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes10_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes10 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes11_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes11 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes12_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes12 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes13_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes13 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes14_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes14 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes15_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes15 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes16_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes16 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes17_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes17 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes18_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes18 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes19_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes19 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


% --- Executes on mouse press over axes background.
function axes20_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to axes20 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    schedule_measurements(hObject, eventdata, handles);
end


function schedule_measurements(hObject, eventdata, handles)

    persistent add_coords;
    global meas_add_remove;
    global boxes;

    if (isempty(meas_add_remove))
        meas_add_remove = 1;
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

    switch meas_add_remove
        case 1 % Add boxes to axes
            if isempty(add_coords) % First mouse click
                add_coords(1) = mousepos(1,1);

            else % Second mouse click
                add_coords(2) = mousepos(1,1);
                create_a_box(add_coords, [hObject, handles.meas_total]); 
                [boxes(end).handle, boxes(end).total_handle] = ...
                    draw_a_box(add_coords, [hObject, handles.meas_total]);
                clear add_coords;

            end

        case -1 % Remove boxes from axes
            clear add_cords;
            remove_coords(1) = mousepos(1,1);
            remove_a_box(remove_coords, [hObject, handles.meas_total]);

        case 2 % Add boxes in a pattern
            if isempty(add_coords) % First mouse click
                add_coords(1) = mousepos(1,1);

            else % Second mouse click
                add_coords(2) = mousepos(1,1);
                create_many_much_boxen(add_coords, [hObject, handles.meas_total], handles); 
                clear add_coords;

            end

        case 0 % Edit box information
            clear add_coords;
            edit_coords(1) = mousepos(1,1);
            edit_a_box(edit_coords, [hObject, handles.meas_total], handles);

        otherwise 
            clear add_coords;

    end

    % Debugging
    assignin('base', 'boxes', boxes(:));
end


function create_a_box(coords, axes_handles)
    global boxes;
    
    if (isempty(boxes))
        boxes = struct('ground_station', [], ...
            'type', [], ...
            'x', [0, 0], ...
            'handle', [], ...
            'total_handle', []);
    end
    
    % Rectangles can only have positive deltas
    if (coords(2) < coords(1)) 
        coords_temp = coords(1);
        coords(1) = coords(2);
        coords(2) = coords_temp;
    end

%     [meas, meas_on_total] = draw_a_box(coords, axes_handles);

    % Save a box in memory
    boxes(end+1) = struct('ground_station', get(axes_handles(1), 'UserData'), ...
        'type', 'measurement', ...
        'x', [coords(1), coords(2)], ...
        'handle', [], ...
        'total_handle', []);
end


function remove_a_box(coords, axes_handles)
    global boxes;

    i = 1;
    while (i <= length(boxes)) % While loop, *not* for loop (we need length recalculated every iteration)
        if (isequal(axes_handles(1), get(boxes(i).handle, 'parent')))
            if ((boxes(i).x(1) <= coords(1)) && (coords(1) <= boxes(i).x(2)))

                % Delete boxes
                delete(boxes(i).handle);
                delete(boxes(i).total_handle);

                % Delete the record
                boxes(i) = [];
            end
        end
        i = i + 1;
    end
end


function edit_a_box(coords, axes_handles, handles)
    global boxes;

    i = 1;
    % While loop, *not* for loop (we need length recalculated every iteration)
    while (i <= length(boxes)) 
        if (isequal(axes_handles(1), get(boxes(i).handle, 'parent')))
            if ((boxes(i).x(1) <= coords(1)) && (coords(1) <= boxes(i).x(2)))

                % Open visualization window
                [returned_handle, edited_times(1), edited_times(2)] = ...
                    meas_edit('meas_sched', handles.figure1, ...
                    'meas_times', [boxes(i).x(1), boxes(i).x(2)]);

                if (~strcmp(returned_handle, 'Cancel'))
                    % Assign new values into boxes
                    boxes(i).x(1) = edited_times(1);
                    boxes(i).x(2) = edited_times(2);

                    % Redraw all the boxes
                    redraw_boxes(hObject, eventdata, handles);

                end
            end
        end
        i = i + 1;
    end
end


function create_many_much_boxen(coords, axes_handles, handles)
    global boxes;

    % Bring up GUI to visualization original information and set repeat
    [returned_handle, edited_times(1), edited_times(2), repeat_freq, repeat_until] = ...
                    pattern_add('meas_sched', handles.figure1, ...
                    'meas_times', [coords(1), coords(2)]);
     
    if (~strcmp(returned_handle, 'Cancel'))
        duration = edited_times(2) - edited_times(1);

        % Add in the repeated boxes
        for new_start = edited_times(1):abs(repeat_freq):repeat_until
            series_coords(1,1) = new_start;
            series_coords(2,1) = new_start + duration;
            
            create_a_box(series_coords, axes_handles);
            [boxes(end).handle, boxes(end).total_handle] = ...
                draw_a_box(series_coords, [axes_handles(1), axes_handles(2)]);
        end
        
    else
        % If they cancel out, delete the original box they made
        remove_a_box(coords, axes_handles)
    end
end


function [meas, meas_on_total] = draw_a_box(coords, axes_handles)
    % Rectangles can only have positive deltas
    if (coords(2) < coords(1)) 
        coords_temp = coords(1);
        coords(1) = coords(2);
        coords(2) = coords_temp;
    end

    % Prevent boxes from showing up out of chart bounds
    visible_time = get(axes_handles(1), 'XLim');
    if (visible_time(2) < coords(2))
        x = coords(1);
        dx = visible_time(2) - coords(1);
    elseif (visible_time(1) > coords(1))
        x = visible_time(1);
        datestr(visible_time(1))
        dx = coords(2) - visible_time(1);
    else
        x = coords(1);
        dx = coords(2) - coords(1);
    end

    % Plot a box on parent axes
    meas = rectangle('Position',[x,0,dx,1], 'EdgeColor','g',...
        'LineWidth', 3);%,'FaceColor',[175/255 1 175/255]);
    set(meas, 'parent', axes_handles(1));

    % Plot a box on total axes
    meas_on_total = rectangle('Position',[x,0,dx,1], 'EdgeColor','g',...
        'LineWidth', 5,'FaceColor','g');
    set(meas_on_total, 'parent', axes_handles(2));
end


function redraw_boxes(hObject, eventdata, handles)
    % This function redraws *all* the boxes
    % Box images have already been deleted (probably by cla)
    
    global boxes;

    for axes_num = 1:length(handles.axes_handles)-1
        i = 2; % The first box is a decoy structure box
        while (i <= length(boxes)) 
            if (get(handles.axes_handles(axes_num), 'UserData') == ...
                    boxes(i).ground_station)
                
                % Draw the boxes and return the handles to the drawings
                [meas, meas_on_total] = draw_a_box(boxes(i).x, ...
                    [handles.axes_handles(axes_num), handles.meas_total]);
                
                % Save the new rectangles to the data structure
                boxes(i).handle = meas;
                boxes(i).total_handle = meas_on_total;
            end

            i = i + 1;
        end
    end
end


function harvest_gs(hObject, eventdata, handles)
    % This function will go through all the defined ground stations and
    % save their data to the options structure
    
    global measOptions;
%     measOptions = odtbxOptions('measurement');
    default_options = odtbxOptions('measurement');
    
    % Clear out the entries in the old structure (this function will
    % completely rebuild the entries)
    
    measOptions = setOdtbxOptions(measOptions, 'gsID', ...
        getOdtbxOptions(default_options, 'gsID'));
    measOptions = setOdtbxOptions(measOptions, 'gsECEF', ...
        getOdtbxOptions(default_options, 'gsECEF'));
    
    for i = 1:length(handles.slide_labels)
        if (~strcmp(get(handles.slide_labels(i), 'String'), '[ ]'))
            % Pull in the full variables from the options structure
            local_gsID = getOdtbxOptions(measOptions, 'gsID');
            local_gsECEF = getOdtbxOptions(measOptions, 'gsECEF');
            
            % Get the new values to add to the structure
            new_gsID = get(handles.slide_labels(i), 'String');
            new_gsECEF = get(handles.slide_labels(i), 'UserData');
            
            % Append the new data to the end of the current data
            local_gsID{end+1} = new_gsID;
            local_gsECEF(1:3, end+1) = new_gsECEF;
            
            % Save it back to the options structure
            measOptions = setOdtbxOptions(measOptions, 'gsID', local_gsID);
            measOptions = setOdtbxOptions(measOptions, 'gsECEF', local_gsECEF);
            
            % Save the plot number to the corresponding axes_handle
            set(handles.axes_handles(i), 'UserData', i);
        
        else
%             set(handles.axes_handles(i), 'UserData', []);
        end
    end
    
%     % Write out variable
    assignin('base', 'optOut', measOptions);
end



function change_gs(axes_handle)
    % This function changes the ground station of a box
    global boxes;
    
    % Change axes_handle ground station
%     set(axes_handle, 'UserData', new_gs_name);
    
    i = 2; % The first box is a decoy structure box
    while (i <= length(boxes)) % While loop, *not* for loop (we need length recalculated every iteration)
        % If it's on the correct axes
        if (axes_handle == get(boxes(i).handle, 'Parent'))
            % Change boxes(i).ground_station
            boxes(i).ground_station = get(axes_handle, 'UserData');
        end
        i = i + 1;
    end
end


% --- Executes on mouse press over axes background.
function meas_total_ButtonDownFcn(hObject, eventdata, handles)
    % hObject    handle to meas_total (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    change_time(hObject, eventdata, handles);
end


function change_time(hObject, eventdata, handles)
    % hObject    handle to time (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Calls GUI to get new information for time
    % This function also adjusts the axes scaling (but not the boxes
    % themselves)
    time_info('meas_sched', handles.figure1);
    datetick('x', 'mm/dd/yy', 'keepticks');
    % Redraw all the boxes on the potentially new axes scale
    plot_meas(hObject, eventdata, handles);
    redraw_boxes(hObject, eventdata, handles);
    
end


% --- Executes on button press in gs_label1.
function gs_label1_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label1;
    handles.gs_axes_current = handles.axes1;
    
    % Update handle structure
    guidata(hObject, handles);
   
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label2.
function gs_label2_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label2 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label2;
    handles.gs_axes_current = handles.axes2;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label3.
function gs_label3_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label3 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label3;
    handles.gs_axes_current = handles.axes3;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label4.
function gs_label4_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label4 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label4;
    handles.gs_axes_current = handles.axes4;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label5.
function gs_label5_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label5 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label5;
    handles.gs_axes_current = handles.axes5;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label6.
function gs_label6_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label6 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label6;
    handles.gs_axes_current = handles.axes6;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label7.
function gs_label7_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label7 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label7;
    handles.gs_axes_current = handles.axes7;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label8.
function gs_label8_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label8 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label8;
    handles.gs_axes_current = handles.axes8;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label9.
function gs_label9_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label9 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
   
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label9;
    handles.gs_axes_current = handles.axes9;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label10.
function gs_label10_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label10 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label10;
    handles.gs_axes_current = handles.axes10;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label11.
function gs_label11_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label11 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label11;
    handles.gs_axes_current = handles.axes11;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label12.
function gs_label12_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label12 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label12;
    handles.gs_axes_current = handles.axes12;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label13.
function gs_label13_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label13 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label13;
    handles.gs_axes_current = handles.axes13;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label14.
function gs_label14_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label14 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label14;
    handles.gs_axes_current = handles.axes14;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label15.
function gs_label15_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label15 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label15;
    handles.gs_axes_current = handles.axes15;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label16.
function gs_label16_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label16 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label16;
    handles.gs_axes_current = handles.axes16;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label17.
function gs_label17_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label17 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label17;
    handles.gs_axes_current = handles.axes17;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label18.
function gs_label18_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label18 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label18;
    handles.gs_axes_current = handles.axes18;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label19.
function gs_label19_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label19 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label19;
    handles.gs_axes_current = handles.axes19;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


% --- Executes on button press in gs_label20.
function gs_label20_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label20 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Change the handle to reference the button that has been pressed (this
    % will be used in gs_info) 
    handles.gs_label_current = handles.gs_label20;
    handles.gs_axes_current = handles.axes20;
    
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Harvest the new data (take it from the button and collect it to the options structure)
    harvest_gs(hObject, eventdata, handles);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.gs_axes_current);
end


function export_schedule_to_csv()
    % Write all of the measurements to a file
    global boxes;

    [name, path] = uiputfile('*.csv','Export To','measurement_schedule.csv');
    if (name ~= 0)
        filename = strcat(path, name);
        fid = fopen(filename, 'w');

        if (fid ~= -1) % Make sure the file opens correctly
            % Print a header
            fprintf(fid, 'Measurement Schedule\n\n');
            fprintf(fid, 'Ground Station, Event Type, Start Date/Time, Finish Date/Time,\n');

            i = 2; % The first box is a decoy structure box
                while (i <= length(boxes)) % While loop, *not* for loop (we need length recalculated every iteration)
                    fprintf(fid, '%s, %s, %s, %s\n', ...
                                boxes(i).ground_station, boxes(i).type, ...
                                datestr(boxes(i).x(1), 'mm/dd/yyyy HH:MM:SS'), ...
                                datestr(boxes(i).x(2), 'mm/dd/yyyy HH:MM:SS'));
                    i = i + 1;
                end
%             close(fid);
        end
    end
end


function export_options_to_workspace()
    global measOptions;
    global boxes;
    
    % Prompt user for new name
    prompt = {'Enter target variable name:'};
    title = 'Export to Workspace';
    lines = 1;
    def = {'measOptions'};
    answer = inputdlg(prompt, title, lines, def);
    
    if (isempty(answer))
        % If the dialog box is cancelled, we don't do anything.
        return;
    end
    
    % Change schedule to the appropriate format
    % gsmeas() expects the [ID, ti, tf] where the start and stop times are
    % in seconds from epoch.
    epoch = getOdtbxOptions(measOptions, 'epoch');
    schedule = [];
    for i = 2:length(boxes)
        ti = (boxes(i).x(1) - epoch) * 24 * 60 * 60; % Was in days
        tf = (boxes(i).x(2) - epoch) * 24 * 60 * 60; % Now in seconds
        new_entry = [boxes(i).ground_station, ti, tf];
        if isempty(schedule)
            schedule = new_entry;
        else
            schedule = vertcat(schedule, new_entry);
        end
    end
    
    % Assign schedule to measOptions
    exportedOptions = setOdtbxOptions(measOptions, 'Schedule', schedule);

    % Write out variable
    assignin('base', answer{1}, exportedOptions);
%     save(answer{1}, 'measOptions');

end


function import_options_from_workspace(hObject, eventdata, handles)
    global measOptions;
    global boxes;
    global T;
    global X;
    
    % Clear any current values from measOptions
    measOptions = odtbxOptions('measurement');
    
    % Prompt user for new name
    prompt = {'Enter workspace variable name:'};
    title = 'Import from Workspace';
    lines = 1;
    def = {'measOptions'};
    answer = inputdlg(prompt, title, lines, def);
    
    if (isempty(answer))
        % If the dialog box is cancelled, we don't do anything.
        return;
    end
    
    has_options = 1;
    % Read in options structure
    try
        measOptions = evalin('base', answer{1});
    catch exception
        errordlg(exception.message, 'Import error!');
        has_options = 0;
    end
    
    if (has_options)
        % Create the imported ground stations
        % Delete all the old ground stations (if there were any)
        for i = 1:length(handles.slide_labels)
            if (~strcmp(get(handles.slide_labels(i), 'String'), '[ ]'))
                set(handles.slide_labels(i), 'String', '[ ]');
                set(handles.slide_labels(i), 'UserData', '');
            end
        end

        % Write all the new GS data (seed)
        gsID_local = getOdtbxOptions(measOptions, 'gsID');
        gsECEF_local = getOdtbxOptions(measOptions, 'gsECEF');
        for j = 1:length(gsID_local)
            set(handles.slide_labels(j), 'String', gsID_local{j});
            set(handles.slide_labels(j), 'UserData', gsECEF_local(:,j));
            set(handles.axes_handles(j), 'UserData', j);
        end

        % Update the box structure
        for k = 1:length(handles.axes_handles)-1
            change_gs(handles.axes_handles(k));
        end

        % We don't know if there will be valid data to propagate a state or
        % not. If there is, then great. We'll trust the user knows whether or
        % not the schedule and what they've propagated match up. If not, just
        % get rid of the schedule data and start fresh.
        if (~isempty(T) && ~isempty(X))
            % Clear box data that already exists in the workspace
            boxes = [];
            
            % Import schedule
            schedule = getOdtbxOptions(measOptions, 'Schedule');
            epoch = getOdtbxOptions(measOptions, 'epoch');

            % Import schedule information
            if (~isempty(epoch))
                for mysched = 1:length(schedule)

                    % Convert times
                    add_coords(1) = (schedule(mysched,2) / 60 / 60 / 24 ) + epoch; % Was in seconds
                    add_coords(2) = (schedule(mysched,3) / 60 / 60 / 24 ) + epoch; % Now in days

                    % Find correct axes handles
                    gs = schedule(mysched, 1);
                    for myaxes = 1:length(handles.axes_handles)-1
                        axes_gs = get(handles.axes_handles(myaxes), 'UserData');
                        if (gs == axes_gs)
                            display_axes = handles.axes_handles(myaxes);
                            break;
                        end
                    end

                    % Create and display box
                    create_a_box(add_coords, [display_axes, handles.meas_total]); 
    %                 [boxes(end).handle, boxes(end).total_handle] = ...
    %                                 draw_a_box(add_coords, [display_axes, handles.meas_total]);
                end
            end

            % Pop up a window telling the user about the assumptions we made with
            % regards to the propagation time and interval
            errordlg({'Functions propagated with current satellite data.', ...
                'May not be data originally used for scheduling.'}, ...
                'Warning!');
        else
            % Inform the user as to why no data was shown
            errordlg({'No propagated satellite data available.', ...
                'Please enter this data in Satellite menu.'}, ...
                'Information');
        end

        % Delete schedule in imported measOptions (we don't want them limiting
        % what gsmeas generates)
        measOptions = setOdtbxOptions(measOptions, 'Schedule', []);
        generate_Callback(hObject, eventdata, handles)
    end
end


function change_satellite(hObject, eventdata, handles)

    % Set up the necessary structures (if they haven't been created
    % already)
    global time_prop;
    if (isempty(time_prop))
        time_prop = struct('begin', [], ...
                          'increment', [], ...
                          'end', []);
    end
    
    global sat_state_prop;
    if (isempty(sat_state_prop))
        sat_state_prop = struct('pos_x', [], ...
                               'pos_y', [], ...
                               'pos_z', [], ...
                               'vel_x', [], ...
                               'vel_y', [], ...
                               'vel_z', []);
    end
    
    % GUI will gather the data and save it to the structures
    output = satellite_edit('meas_sched', handles.figure1);
    if (~strcmp(output, 'Cancel'))
        % Do the actual propagation of the cartesian state using the specified
        % function
        global propagator;
        time = time_prop.begin:time_prop.increment:time_prop.end;
        coords = [sat_state_prop.pos_x;
                  sat_state_prop.pos_y;
                  sat_state_prop.pos_z;
                  sat_state_prop.vel_x;
                  sat_state_prop.vel_y;
                  sat_state_prop.vel_z];

        % Set numerical integration tolerances
        opts = odeset('reltol',1e-9,'abstol',1e-9);      
        mu = 3.986e5;           % Pancake gravitational parameter
        
        global T;
        global X;
        try
            [T,X] = integ(propagator, time, coords, opts, mu);
        catch exceptions
            errordlg(exceptions.message, 'Propagation Error!');
        end
    end
end


function make_meas()
    global T;
    global X;
    global measOptions;
    
    try
        [y, H, R] = gsmeas(T, X, measOptions);
        separate_measurements(y);
    catch exceptions
        errordlg(exceptions.message, 'Measurement Error!');
    end
end


function separate_measurements(meas_in)
    global measurements;
    global measOptions;

    % Create the structure template
    measurements = struct('plot', [], ...
        'ground_station', [], ...
        'type', [], ...
        'data', []);

    % Figure out how many ground stations we have measurements for
    ground_stations = getOdtbxOptions(measOptions, 'gsID');
    num_gs = length(ground_stations);
    
    % Figure out what measurements we have
    options = {};
    option_size = [];
    if (getOdtbxOptions(measOptions, 'useRange'))
        options{end+1} = 'useRange';
        option_size(end+1) = 1;
    end
    if (getOdtbxOptions(measOptions, 'useRangeRate'))
        options{end+1} = 'useRangeRate';
        option_size(end+1) = 1;
    end
    if (getOdtbxOptions(measOptions, 'useDoppler'))
        options{end+1} = 'useDoppler';
        option_size(end+1) = 1;
    end
    if (getOdtbxOptions(measOptions, 'useUnit'))
        options{end+1} = 'useUnit';
        option_size(end+1) = 3;
    end
    if (getOdtbxOptions(measOptions, 'useAngles'))
        options{end+1} = 'useAngles';
        option_size(end+1) = 2;
    end
    num_options = length(options);
    
    % Save the data into a more useable form
    row_current = 1;
    for gs = 1:num_gs
        for opt = 1:num_options
            measurements(end+1) = struct('plot', gs, ...
                'ground_station', ground_stations(gs), ...
                'type', options(opt), ...
                'data', meas_in(row_current:row_current+option_size(opt)-1,:));
            row_current = row_current + option_size(opt);
        end
    end
end


function plot_meas(hObject, eventdata, handles)
    global measurements;
    global T;
    global measOptions;
    
    if (~isempty(T))
        % Convert relative time to absolute time
        time_abs = getOdtbxOptions(measOptions, 'epoch') + T*1/60*1/60*1/24;

        % In the future, the controls that dictate which plots will be shown
        % will go in here. For now, we show all the plots.
        for axes_loop = 1:length(handles.axes_handles)-1
            user_data = get(handles.axes_handles(axes_loop), 'UserData');
            if (~isempty(user_data) && user_data > 0)
                cla(handles.axes_handles(axes_loop));
                for meas_loop = 2:length(measurements)
                    plot_num = measurements(meas_loop).plot;
                    if (user_data == plot_num)       
                        % Set up color scheme and normalizations for different data here
                        switch measurements(meas_loop).type
                            case 'useRange'
                                color = 'r';
                                % Scale the data to fit the axes, take the max
                                % to be 1
                                max_val = max(max(abs(measurements(meas_loop).data)));
                                scaled_data = measurements(meas_loop).data / (2 * max_val) +.5
                            case 'useRangeRate'
                                color = 'b';
                                % Scale the data to fit the axes, take the max
                                % to be 1
                                max_val = max(max(abs(measurements(meas_loop).data)));
                                scaled_data = measurements(meas_loop).data / (2 * max_val) +.5
                            case 'useDoppler'
                                color = 'k';
                                % Scale the data to fit the axes, take the max
                                % to be one
                                max_val = max(max(abs(measurements(meas_loop).data)));
                                scaled_data = measurements(meas_loop).data / (2 * max_val) +.5
                            case 'useUnit'
                                color = 'c';
                                % These are unit vectors, they should already
                                % be scaled
                                scaled_data = measurements(meas_loop).data / 2 + .5
                            case 'useAngles'
                                color = 'm';
                                % Scale the data to fit the axes
                                measurements(meas_loop).data
                                scaled_data(1,:) = measurements(meas_loop).data(1,:) / (2*pi)
                                scaled_data(2,:) = measurements(meas_loop).data(2,:) / (pi/2)
                            otherwise
                                color = 'y';
                                % Scale the data to fit the axes
                                max_val = max(max(measurements(meas_loop).data));
                                scaled_data = measurements(meas_loop).data / (2 * max_val) +.5
                        end

                        line(time_abs, scaled_data, 'Color', color, 'LineWidth', 2, ...
                            'Parent', handles.axes_handles(axes_loop));
                        hold on;
                    end
                end
            end
            hold off;
        end
    end
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear_data();
% Hint: delete(hObject) closes the figure
delete(hObject);

end


function clear_data()
clear global T;
clear global X;
clear global measOptions;
clear global measurements;
clear global boxes;
clear global meas_add_remove;
clear global sat_state_prop;
clear global time_prop;
clear global propagator;
end
