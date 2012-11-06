function varargout = gs_options(varargin)
% GS_OPTIONS MATLAB code for gs_options.fig
%      GS_OPTIONS, by itself, creates a new GS_OPTIONS or raises the existing
%      singleton*.
%
%      H = GS_OPTIONS returns the handle to a new GS_OPTIONS or the handle to
%      the existing singleton*.
%
%      GS_OPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GS_OPTIONS.M with the given input arguments.
%
%      GS_OPTIONS('Property','Value',...) creates a new GS_OPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gs_options_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gs_options_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gs_options

% Last Modified by GUIDE v2.5 02-Nov-2012 13:38:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gs_options_OpeningFcn, ...
                   'gui_OutputFcn',  @gs_options_OutputFcn, ...
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


% --- Executes just before gs_options is made visible.
function gs_options_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gs_options (see VARARGIN)

% Choose default command line output for gs_options
handles.output = hObject;

% Set up the display with the correct values
global measOptions;
if getOdtbxOptions(measOptions, 'useRange')
    set(handles.range, 'Value', getOdtbxOptions(measOptions, 'useRange'));
end
if getOdtbxOptions(measOptions, 'rangeType')
    set(handles.rangeType, 'String', getOdtbxOptions(measOptions, 'rangeType'));
end
if getOdtbxOptions(measOptions, 'useRangeRate')
    set(handles.range_rate, 'Value', getOdtbxOptions(measOptions, 'useRangeRate'));
end
if getOdtbxOptions(measOptions, 'useDoppler')
    set(handles.doppler, 'Value', getOdtbxOptions(measOptions, 'useDoppler'));
end
if getOdtbxOptions(measOptions, 'useUnit')
    set(handles.unit_vec, 'Value', getOdtbxOptions(measOptions, 'useUnit'));
end
if getOdtbxOptions(measOptions, 'useAngles')
    set(handles.angles, 'Value', getOdtbxOptions(measOptions, 'useAngles'));
end
if getOdtbxOptions(measOptions, 'epoch')
    set(handles.epoch, 'String', getOdtbxOptions(measOptions, 'epoch'));
end
if getOdtbxOptions(measOptions, 'gsElevationConstraint')
    set(handles.elevationconstraint, 'String', getOdtbxOptions(measOptions, 'gsElevationConstraint'));
end
if getOdtbxOptions(measOptions, 'frequencyTransmit')
    set(handles.trans_freq, 'String', getOdtbxOptions(measOptions, 'frequencyTransmit'));
end
if getOdtbxOptions(measOptions, 'rSigma')
    set(handles.meas_cov, 'String', getOdtbxOptions(measOptions, 'rSigma'));
end
if getOdtbxOptions(measOptions, 'useLightTime')
    set(handles.light_time, 'Value', getOdtbxOptions(measOptions, 'useLightTime'));
end
if getOdtbxOptions(measOptions, 'useGPSIonosphere')
    set(handles.gps_ionosphere, 'Value', getOdtbxOptions(measOptions, 'useGPSIonosphere'));
end
if getOdtbxOptions(measOptions, 'useIonosphere')
    set(handles.ionosphere, 'Value', getOdtbxOptions(measOptions, 'useIonosphere'));
end
if getOdtbxOptions(measOptions, 'useTroposphere')
    set(handles.troposphere, 'Value', getOdtbxOptions(measOptions, 'useTroposphere'));
end
if getOdtbxOptions(measOptions, 'EarthAtmMaskRadius')
    set(handles.earth_atm_mask_radius, 'String', getOdtbxOptions(measOptions, 'EarthAtmMaskRadius'));
end
if getOdtbxOptions(measOptions, 'PrecnNutnExpire')
    set(handles.Precn_Nutn_Expire, 'String', getOdtbxOptions(measOptions, 'PrecnNutnExpire'));
end

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
   disp('Improper input arguments: gs_options.m');
   disp('-----------------------------------------------------');
else    
% UIWAIT makes gs_options wait for user response (see UIRESUME)
    uiwait(handles.figure1);
end

end

% --- Outputs from this function are returned to the command line.
function varargout = gs_options_OutputFcn(hObject, eventdata, handles) 
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
global measOptions;

% Save all the changes from the form to the structurecd
measOptions = setOdtbxOptions(measOptions, 'useRange', get(handles.range, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'rangeType', get(handles.rangeType, 'String'));
measOptions = setOdtbxOptions(measOptions, 'useRangeRate', get(handles.range_rate, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'useDoppler', get(handles.doppler, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'useUnit', get(handles.unit_vec, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'useAngles', get(handles.angles, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'epoch', get(handles.epoch, 'String'));
measOptions = setOdtbxOptions(measOptions, 'gsElevationConstraint', get(handles.elevationconstraint, 'String'));
measOptions = setOdtbxOptions(measOptions, 'frequencyTransmit', get(handles.trans_freq, 'String'));
measOptions = setOdtbxOptions(measOptions, 'rSigma', get(handles.meas_cov, 'String'));
measOptions = setOdtbxOptions(measOptions, 'useLightTime', get(handles.light_time, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'useGPSIonosphere', get(handles.gps_ionosphere, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'useIonosphere', get(handles.ionosphere, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'useTroposphere', get(handles.troposphere, 'Value'));
measOptions = setOdtbxOptions(measOptions, 'EarthAtmMaskRadius', get(handles.earth_atm_mask_radius, 'String'));
measOptions = setOdtbxOptions(measOptions, 'PrecnNutnExpire', get(handles.Precn_Nutn_Expire, 'String'));

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);
end


% --- Executes on button press in range.
function range_Callback(hObject, eventdata, handles)
% hObject    handle to range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of range
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


function elevationconstraint_Callback(hObject, eventdata, handles)
% hObject    handle to elevationconstraint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elevationconstraint as text
%        str2double(get(hObject,'String')) returns contents of elevationconstraint as a double
end



function epoch_Callback(hObject, eventdata, handles)
% hObject    handle to epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epoch as text
%        str2double(get(hObject,'String')) returns contents of epoch as a double
end


% --- Executes during object creation, after setting all properties.
function epoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function rangeType_Callback(hObject, eventdata, handles)
% hObject    handle to rangeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rangeType as text
%        str2double(get(hObject,'String')) returns contents of rangeType as a double
end


% --- Executes during object creation, after setting all properties.
function rangeType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in range_rate.
function range_rate_Callback(hObject, eventdata, handles)
% hObject    handle to range_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of range_rate
end


% --- Executes on button press in doppler.
function doppler_Callback(hObject, eventdata, handles)
% hObject    handle to doppler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doppler
end


% --- Executes on button press in unit_vec.
function unit_vec_Callback(hObject, eventdata, handles)
% hObject    handle to unit_vec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of unit_vec
end


function trans_freq_Callback(hObject, eventdata, handles)
% hObject    handle to trans_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans_freq as text
%        str2double(get(hObject,'String')) returns contents of trans_freq as a double
end


% --- Executes during object creation, after setting all properties.
function trans_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function meas_cov_Callback(hObject, eventdata, handles)
% hObject    handle to meas_cov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meas_cov as text
%        str2double(get(hObject,'String')) returns contents of meas_cov as a double
end


% --- Executes during object creation, after setting all properties.
function meas_cov_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meas_cov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in light_time.
function light_time_Callback(hObject, eventdata, handles)
% hObject    handle to light_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of light_time
end


% --- Executes on button press in gps_ionosphere.
function gps_ionosphere_Callback(hObject, eventdata, handles)
% hObject    handle to gps_ionosphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gps_ionosphere
end


% --- Executes on button press in ionosphere.
function ionosphere_Callback(hObject, eventdata, handles)
% hObject    handle to ionosphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ionosphere
end


% --- Executes on button press in troposphere.
function troposphere_Callback(hObject, eventdata, handles)
% hObject    handle to troposphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of troposphere
end


% --- Executes on button press in angles.
function angles_Callback(hObject, eventdata, handles)
% hObject    handle to angles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of angles
end


% --- Executes on button press in charged_particle.
function charged_particle_Callback(hObject, eventdata, handles)
% hObject    handle to charged_particle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of charged_particle
end


function earth_atm_mask_radius_Callback(hObject, eventdata, handles)
% hObject    handle to earth_atm_mask_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of earth_atm_mask_radius as text
%        str2double(get(hObject,'String')) returns contents of earth_atm_mask_radius as a double
end


% --- Executes during object creation, after setting all properties.
function earth_atm_mask_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to earth_atm_mask_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function Precn_Nutn_Expire_Callback(hObject, eventdata, handles)
% hObject    handle to Precn_Nutn_Expire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Precn_Nutn_Expire as text
%        str2double(get(hObject,'String')) returns contents of Precn_Nutn_Expire as a double
end


% --- Executes during object creation, after setting all properties.
function Precn_Nutn_Expire_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Precn_Nutn_Expire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
