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
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % visualization the above text to modify the response to help meas_sched

    % Last Modified by GUIDE v2.5 07-Nov-2012 13:37:15

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
    global boxes;
    clear global boxes; % Get rid of data from previous runs

    global meas_add_remove;
%     clear global meas_add_remove;
    meas_add_remove = 1;
    set(handles.meas_schedule_mode,'SelectedObject',[handles.Add]);

    global measOptions;
%     clear global measOptions;
    measOptions = odtbxOptions('measurement');
    
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
    set(handles.axes_handles, 'UserData', ['Unspecified']);

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
    global measOptions;
%     measOptions = setOdtbxOptions(measOptions, 'useRange', true);
%     measOptions = setOdtbxOptions(measOptions, 'rangeType', '2way');
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
    import_options_from_workspace();
end


% --------------------------------------------------------------------
function quit_Callback(hObject, eventdata, handles)
    % hObject    handle to quit (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
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

    if (coords(2) < coords(1)) % Rectangles can only have positive deltas
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
    while (i <= length(boxes)) % While loop, *not* for loop (we need length recalculated every iteration)
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
                    redraw_boxes();

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
    % Prevent boxes from showing up out of chart bounds
    visible_time = get(axes_handles(1), 'XLim');
    if (visible_time(2) < coords(2))
        x = coords(1);
        dx = visible_time(2) - coords(1);
    elseif (visible_time(1) > coords(1))
        x = visible_time(1);
        dx = coords(2) - visible_time(1);
    else
        x = coords(1);
        dx = coords(2) - coords(1);
    end
    
    % Plot a box on parent axes
    meas = rectangle('Position',[x,0,dx,1], 'EdgeColor','g','LineWidth', 3);%,'FaceColor',[175/255 1 175/255]);
    set(meas, 'parent', axes_handles(1));

    % Plot a box on total axes
    meas_on_total = rectangle('Position',[x,0,dx,1], 'EdgeColor','g','LineWidth', 5,'FaceColor','g');
    set(meas_on_total, 'parent', axes_handles(2));
end


function redraw_boxes()
    % This function redraws *all* the boxes

    global boxes;

    i = 2; % The first box is a decoy structure box
    while (i <= length(boxes)) % While loop, *not* for loop (we need length recalculated every iteration)

        % Find the parent handles (so we know where to redraw the boxes)
        axes_handle = get(boxes(i).handle, 'parent');
        axes_total_handle = get(boxes(i).total_handle, 'parent');

        % Erase the old boxes by their handles
        delete(boxes(i).handle);
        delete(boxes(i).total_handle);
            
        % Draw the boxes and return the handles to the drawings
        [meas, meas_on_total] = draw_a_box(boxes(i).x, [axes_handle, axes_total_handle]);

        % Save the new rectangles to the data structure
        boxes(i).handle = meas;
        boxes(i).total_handle = meas_on_total;

        i = i + 1;
    end
end


function change_gs(axes_handle, new_gs_name)
    % This function changes the ground station of a box
    global boxes;
    
    % Change axes_handle ground station
    set(axes_handle, 'UserData', [new_gs_name]);
    
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
    redraw_boxes();
end


% --- Executes on button press in gs_label1.
function gs_label1_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label1;
    % Update handle structure
    guidata(hObject, handles);
   
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes1, get(hObject, 'String'));
end


% --- Executes on button press in gs_label2.
function gs_label2_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label2 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label2;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes2, get(hObject, 'String'));
end


% --- Executes on button press in gs_label3.
function gs_label3_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label3 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label3;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes3, get(hObject, 'String'));
end


% --- Executes on button press in gs_label4.
function gs_label4_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label4 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label4;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes4, get(hObject, 'String'));
end


% --- Executes on button press in gs_label5.
function gs_label5_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label5 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label5;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes5, get(hObject, 'String'));
end


% --- Executes on button press in gs_label6.
function gs_label6_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label6 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label6;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes6, get(hObject, 'String'));
end


% --- Executes on button press in gs_label7.
function gs_label7_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label7 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label7;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes7, get(hObject, 'String'));
end


% --- Executes on button press in gs_label8.
function gs_label8_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label8 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label8;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes8, get(hObject, 'String'));
end


% --- Executes on button press in gs_label9.
function gs_label9_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label9 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label9;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes9, get(hObject, 'String'));
end


% --- Executes on button press in gs_label10.
function gs_label10_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label10 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label10;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes10, get(hObject, 'String'));
end


% --- Executes on button press in gs_label11.
function gs_label11_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label11 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label11;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes11, get(hObject, 'String'));
end


% --- Executes on button press in gs_label12.
function gs_label12_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label12 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label12;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes12, get(hObject, 'String'));
end


% --- Executes on button press in gs_label13.
function gs_label13_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label13 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label13;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes13, get(hObject, 'String'));
end


% --- Executes on button press in gs_label14.
function gs_label14_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label14 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label14;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes14, get(hObject, 'String'));
end


% --- Executes on button press in gs_label15.
function gs_label15_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label15 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label15;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes15, get(hObject, 'String'));
end


% --- Executes on button press in gs_label16.
function gs_label16_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label16 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label16;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes16, get(hObject, 'String'));
end


% --- Executes on button press in gs_label17.
function gs_label17_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label17 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label17;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes17, get(hObject, 'String'));
end


% --- Executes on button press in gs_label18.
function gs_label18_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label18 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label18;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes18, get(hObject, 'String'));
end


% --- Executes on button press in gs_label19.
function gs_label19_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label19 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label19;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes19, get(hObject, 'String'));
end


% --- Executes on button press in gs_label20.
function gs_label20_Callback(hObject, eventdata, handles)
    % hObject    handle to gs_label20 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    handles.gs_label_current = handles.gs_label20;
    % Update handle structure
    guidata(hObject, handles);
       
    % Call the GUI to get updated info
    gs_info('meas_sched', handles.figure1);
    
    % Change the boxes structures to reflect new changes
    change_gs(handles.axes20, get(hObject, 'String'));
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
    
    % Prompt user for new name
    prompt = {'Enter target variable name:'};
    title = 'Export to Workspace';
    lines = 1;
    def = {'measOptions'};
    answer = inputdlg(prompt, title, lines, def);

    % Write out variable
    assignin('base', answer{1}, measOptions);

end


function import_options_from_workspace()
    global measOptions;
    
    % Prompt user for new name
    prompt = {'Enter workspace variable name:'};
    title = 'Import from Workspace';
    lines = 1;
    def = {'measOptions'};
    answer = inputdlg(prompt, title, lines, def);

    % Read in variable
    try
        measOptions = evalin('base', answer{1});
    catch exception
        errordlg(exception.message, 'Enter a function name!');
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

        try
            [T,X] = integ(propagator, time, coords, opts, mu);
        catch exceptions
            errordlg(exceptions.message, 'Propagation Error!');
        end
    end
end
