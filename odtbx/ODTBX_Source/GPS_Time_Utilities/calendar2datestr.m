% CALENDAR2DATESTR - reformat 1x6 time vector into formatted string
%
% datestr = calendar2datestr(calendardate,output_format)
% 
% Input:    calendardate : [2002,6,1,12,0,0.00000000000]
%           output_format: 0-   6/1/2002 12:00:00  (with zero padding on times)
%                          1-   1-JUN-2002 12:00:00  (with zero padding on times)
%                          2-   1 Jun 2002 12:00:00.000000000  (with zero padding on times)
%                          3-   YYYYMMDD.HHMMSS  (with zero padding)
% Output:   datestr:       formatted string
%
% See also: DATESTR2CALENDAR, DISPLAYTIME
% Created: M.Moreau  March 2005

function datestr = calendar2datestr(calendardate,output_format)

if nargin <2
    output_format = 0;
elseif (output_format ~= 0) & (output_format ~= 1) & (output_format ~= 2) & (output_format ~= 3)
    error('Invalid output format.')
end


switch output_format
case 0
	timevector = calendardate(4:6);
	
	if timevector(1)<10 
       hh = ['0',num2str(timevector(1))];
	else
       hh = num2str(timevector(1));
	end
	
	if timevector(2)<10 
       mm = ['0',num2str(timevector(2))];
	else
       mm = num2str(timevector(2));
	end
	
	if timevector(3)<10 
       ss = ['0',sprintf('%1d',timevector(3))];
	else
       ss = sprintf('%2d',timevector(3));
	end
	
	datestr = sprintf('%d/%d/%d %s:%s:%s',calendardate(2:3),calendardate(1),hh,mm,ss);
case 1
	monthstr = ['JAN'; 'FEB'; 'MAR'; 'APR'; 'MAY'; 'JUN'; 'JUL'; 'AUG'; 'SEP'; 'OCT'; 'NOV'; 'DEC'];
	datestr = sprintf('%2d-%s-%4d',calendardate(3),monthstr(calendardate(2),:),calendardate(1));
	
	timevector = calendardate(4:6);
	
	if timevector(1)<10 
       hh = ['0',num2str(timevector(1))];
	else
       hh = num2str(timevector(1));
	end
	
	if timevector(2)<10 
       mm = ['0',num2str(timevector(2))];
	else
       mm = num2str(timevector(2));
	end
	
	if timevector(3)<10 
       ss = ['0',sprintf('%1d',timevector(3))];
	else
       ss = sprintf('%2d',timevector(3));
	end
	
	datestr = [datestr,' ',hh,':',mm,':',ss];

case 2
    monthstr = ['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'];
	datestr = sprintf('%2d %s %4d',calendardate(3),monthstr(calendardate(2),:),calendardate(1));
	
	timevector = calendardate(4:6);
	
	if timevector(1)<10 
       hh = ['0',num2str(timevector(1))];
	else
       hh = num2str(timevector(1));
	end
	
	if timevector(2)<10 
       mm = ['0',num2str(timevector(2))];
	else
       mm = num2str(timevector(2));
	end
	
	if timevector(3)<10 
       ss = ['0',sprintf('%0.9f',timevector(3))];
	else
       ss = sprintf('%0.9f',timevector(3));
	end
	
	datestr = [datestr,' ',hh,':',mm,':',ss];
case 3
    datevector = calendardate(1:3);
    
    yyyy = num2str(datevector(1));
    
	if datevector(2)<10 
       month = ['0',num2str(datevector(2))];
	else
       month = num2str(datevector(2));
	end

  	if datevector(3)<10 
       day = ['0',sprintf('%1d',datevector(3))];
	else
       day = sprintf('%2d',datevector(3));
	end

    timevector = calendardate(4:6);
	
	if timevector(1)<10 
       hh = ['0',num2str(timevector(1))];
	else
       hh = num2str(timevector(1));
	end
	
	if timevector(2)<10 
       mm = ['0',num2str(timevector(2))];
	else
       mm = num2str(timevector(2));
	end
	
	if timevector(3)<10 
       ss = ['0',sprintf('%1d',timevector(3))];
	else
       ss = sprintf('%2d',timevector(3));
	end
	
	datestr = [yyyy,month,day,'.',hh,mm,ss];
    
end
