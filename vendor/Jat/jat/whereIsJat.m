function JATPath = whereIsJat()

fp = mfilename('fullpath');
mn = mfilename;
pathlength=length(fp)-length(mn);
JATPath=fp(1:pathlength); %Gets path to m-file's directory