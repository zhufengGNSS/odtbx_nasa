function GetMousePosition(varargin)
%This function works with subplots as well as individual axes in one figure
%or multiple figure, 
%it requires handle to figure object that the functionality need to be
%associated with.
%multiple figure handles can be passed to it. it then finds the axes
%automatically if there are any. it returns with no error otherwise if it
%can not find at least one axes object in a vector of figure/s supplied to
%it.
%it can be set off, if varargin{2}='' or [] or {} or 'off'.
%if nargin = 1 then it asigns the functionality with all the figure/s
%supplied to it.
%Examples:
%
%GetMousePosition([f1, f2, f3, f4, f5, f6], 'on')
%set it off on some.
%GetMousePosition([f1, f3, f5], '')
%or
%GetMousePosition([f1, f3, f5], [])
%GetMousePosition([f1, f3, f5], {})
%
%Auther: Rasam Aliazizi (email:rasama@hotmail.co.uk)-COPY RIGHT 2011
%last tested: 03/10/2011-OK.
%
clc,
%check input, return if input is not, 
%             1) handle to a figure or a vector of handles to a figure, if
%             a vector is supplied then all its values must be handles to
%             figure objects.
%             2) if in the figure object/s atleast one axes object is not
%             found.
narginchk(1, 2),    
if ~all(arrayfun(@ishandle, varargin{1}))
    error(['Function expects its argument to be handle to a valid figure object.' ...
        'If multiple an array of handles are supplied then all must valid handle of figure objects.']);
else
   InputType = get(varargin{1}, 'type');
end
if ~all(strcmp(InputType, 'figure'))
    error(['Function expects its argument to be handle to a valid figure object.' ...
        'If multiple an array of handles are supplied then all must valid handle of figure objects.']);
else
    S.fig = varargin{1};
end
if nargin == 2
    switch varargin{2}
        case 'on'
            set(S.fig, 'windowbuttonmotionfcn', @motion_callback);
        case 'off' 
            set(S.fig, 'windowbuttonmotionfcn', @do_nothing);
            return,
        case ''
            set(S.fig, 'windowbuttonmotionfcn', @do_nothing);
            return,
        case []
            set(S.fig, 'windowbuttonmotionfcn', @do_nothing);
            return,
        case {}
            set(S.fig, 'windowbuttonmotionfcn', @do_nothing);
            return,
        otherwise
            set(S.fig, 'windowbuttonmotionfcn', @motion_callback);
    end
else
    set(S.fig, 'windowbuttonmotionfcn', @motion_callback);
end
S.ax = findobj(S.fig, 'type', 'axes');
if isempty(S.ax)
    return,
else
    set(S.ax,'unit','pix');
    S.AXP = get(S.ax,'pos');
    for i = 1:length(S.fig)
        textBx(i) = uicontrol(S.fig(i), 'style', 'text', 'backgroundcolor', [0 1 1],...
        'foregroundcolor', [1 1 1], 'unit', 'pix', 'position', [10 10 10 10], 'visible', 'off');
    end
end

    function motion_callback(varargin)      
        %make sure the text box fits within the figure so that user can
        %read it.
        screenpos = get(gcf, 'position');
        mouseposition = get(gcf, 'currentpoint');mouseposition(3)=150;mouseposition(4)=30;
        if (mouseposition(3)+mouseposition(1))-screenpos(3)>0
            offset2 = (mouseposition(3)+mouseposition(1))-screenpos(3);
        else
            offset2 = 0;
        end
        if (mouseposition(4)+mouseposition(2))-screenpos(4)>0
            offset1 = (mouseposition(4)+mouseposition(2))-screenpos(4);
        else
            offset1 = 0;
        end
        S.Info = get(gca, {'position', 'xlim', 'ylim'});
        %check if mouse is within an axes
        S.tf = mouseposition(1) >=  S.Info{1}(1) && mouseposition(1) <= S.Info{1}(1) + S.Info{1}(3) && ...
            mouseposition(2) >= S.Info{1}(2) && mouseposition(2) <= S.Info{1}(2) + S.Info{1}(4);
        if S.tf
            %activate a particular text box.
            set(textBx(find(gcf==S.fig)), 'visible', 'on', 'position',...
                [mouseposition(1)-offset2 mouseposition(2)-offset1 mouseposition(3) mouseposition(4)],...
             'string',...
             ['(X,Y) Position:(' num2str(((mouseposition(1)-S.Info{1}(1))*(S.Info{2}(2)-S.Info{2}(1)))/S.Info{1}(3) + S.Info{2}(1))...
             ', ' num2str(((mouseposition(2)-S.Info{1}(2))*(S.Info{3}(2)-S.Info{3}(1)))/S.Info{1}(4) + S.Info{3}(1)) ')']),
        else
            if strcmp(get(textBx(find(gcf==S.fig)), 'visible'), 'on')
                %clean screen.
                set(textBx(find(gcf==S.fig)), 'visible', 'off')
            end
        end
    end
    function do_nothing(varargin)
        %do nothing
    end

end