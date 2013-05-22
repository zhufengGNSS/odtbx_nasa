% create_html    Creates html files from Contents.m files in ODTBX_Source and ODTBX_Examples
%
% Note: the help2html function is unsupported by Mathworks and subject to
% change. In addition, it does not seem to properly respect the "pagetitle"
% argument (2nd argument), which is why that argument is not formatted
% here.
% 
% Create by
%     Keith Speckman    03/15/2008      Original
%     Ravi Mathur       11/07/2012      New HTML folder location
%     Ravi Mathur       05/22/2013      Add meas_sched_gui

f = filesep;

% Generate ODTBX function listing from ODTBX_Source folder
fid_f = fopen(['..',f,'HTML',f,'odtbx_functionlist.html'],'w');
fwrite(fid_f,help2html('ODTBX_Source','ODTBX_Source'),'char');
fclose(fid_f);

% Generate ODTBX tutorials/demos listing from ODTBX_Examples folder
fid_e = fopen(['..',f,'HTML',f,'odtbx_examplelist.html'],'w');
fwrite(fid_e,help2html('ODTBX_Examples','ODTBX_Examples'),'char');
fclose(fid_e);

% Generate JAT Adapters listing from JAT_Adapters folder
fid_j = fopen(['..',f,'HTML',f,'jat_adapters.html'],'w');
fwrite(fid_f,help2html('JAT_Adapters','JAT_Adapters'),'char');
fclose(fid_j);

% Generate Measurement Scheduling GUI Adapters listing from meas_sched_gui folder
fid_j = fopen(['..',f,'HTML',f,'meassched_functionlist.html'],'w');
fwrite(fid_f,help2html('meas_sched_gui','meas_sched_gui'),'char');
fclose(fid_j);