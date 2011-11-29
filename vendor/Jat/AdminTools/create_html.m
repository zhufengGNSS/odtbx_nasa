% create_html    Creates html files from Contents.m files in ODTBX_Source and ODTBX_Examples
%
% Note: the help2html function is unsupported by Mathworks and subject to change
% Create by
%     Keith Speckman    03/15/2008      Original

f = filesep;
% ODTBX_Source
fid_f = fopen(['..',f,'ODTBX_Source',f,'HTML',f,'odtbx_functionlist.html'],'w');
fwrite(fid_f,help2html('ODTBX_Source','ODTBX_Source'),'char');
fclose(fid_f);

% ODTBX_Examples
fid_e = fopen(['..',f,'ODTBX_Source',f,'HTML',f,'odtbx_examplelist.html'],'w');
fwrite(fid_e,help2html('ODTBX_Examples','ODTBX_Examples'),'char');
fclose(fid_e);

% JAT_Adapters
fid_j = fopen(['..',f,'ODTBX_Source',f,'HTML',f,'jat_adapters.html'],'w');
fwrite(fid_f,help2html('JAT_Adapters','JAT_Adapters'),'char');
fclose(fid_j);
