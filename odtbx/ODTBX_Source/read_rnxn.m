function [eph,KlobucharCoefs] = read_rnxn(rnxname)
% READ_RNXN Read rinex formated GPS broadcast ephemeris data
%
% [eph] = read_rnxn(fn)
%
% Reads rinex navigation data file into nx34 matrixs 'eph'
% Ignores, the data contained in the header!
% Returns: 
%     eph [n,38]   - ephemeris records
%
% Written by Mike Moreau: 12/14/2000
% Modified 2/2004

% eph is a nx34 matix made up of the following columns:
%   1 short int prn;                                      
%   2 logical vflg;                                   0=No valid data
%   3 long TofXmission;  Time of subframe 1 transmission, sec of week
%   4 short int s1hlth;                        Subframe 1 health code
%   5 short int codeL2;                               Code on L2 flag
%   6 short int wkn;         GPS week at time of subframe 1 reception
%   7 short int L2Pdata;                               L2 P data flag
%   8 short int ura;                             Satellite's URA code
%   9 short int iodc;                            Issue of data, clock
%  10  double tgd;                               Group delay parameter
%  11  short int tocwk;                  GPS week corresponding to toc
%  12  double toc;          Reference time of clock data parameter set
%  13  double af0;             Clock correction polynomial coefficient
%  14  double af1;             Clock correction polynomial coefficient
%  15  double af2;             Clock correction polynomial coefficient
%  16  short int iode;                        Issue of data, ephemeris
%  17  double crs;          Sine harmonic correction to orbital radius
%  18  double deltan;            Mean motion delta from computed value
%  19  double m0;                                  Mean anomaly at TOE
%  20  double cuc;        Cosine harmonic correction to orbital radius
%  21  double ecc;                                        Eccentricity
%  22  double cus;          Sine harmonic corr to argument of latitude
%  23  double sqrta;                     Square root of semimajor axis
%  24  short int toewk;                  GPS week corresponding to toe
%  25  double toe;                Reference time of ephemeris data set
%  26  short int fti;                                     Fit interval
%  27  double cic;                 Cosine harmonic corr to inclination
%  28  double om0;                              Right ascension at TOE
%  29  double cis;                   Sine harmonic corr to inclination
%  30  double in0;                                  Inclination at TOE
%  31  double crc;        Cosine harmonic correction to orbital radius
%  32  double olc;                          Argument of perigee at TOE
%  33  double omd;                             Rate of right ascension
%  34  double idot;                                Rate of inclination
%
% Mapping of GPS broadcast ephemeris parameters to rinex record:
%   prn toc(yy,mm,dd,hh,mm,ss)       af0                   af1                   af2
%      iode                          crs                deltan                    m0 
%      cuc                           ecc                   cus                 sqrta
%      toesec                        cic                   om0                   cis
%      in0                           crc                   olc                   omd
%      idot                       codeL2               toeweek               L2Pdata
%      ura                      s1heatlh                   tgd                  iodc
%      TofXmission                   fti                spare1                spare2
%
% Relation of rinex formatted data and columns of 'ephem'
%   1prn,2yy,3mm,4dd,5hh,6mm,7ss,8svbias,             9svdrift,        10svdriftrate
%   11IODE,                        12Crs,            13Delta n,                 14M0 
%   15Cuc,                           16e,                17Cus,              18roota
%   19Toe,                         20Cic,              21OMEGA,                22Cis
%   23i,                           24Crc,              25omega,           26OMEGADOT
%   27IDOT,                     28codeL2,               29week,            30L2Pdata
%   31SVacc,                  32SVhealth,                33Tgd,               34IODC
%   35transtime,                   36fit,             37spare1,             38spare2
%
% Note: MATLAB does not recognize the exponential format used in RINEX: -0.5222959D-11
%   This must be replaced with -0.5222959E-11 in the text file for this script to work.

gps_constants
format long e

if nargin == 0
    defname = ('R2_n.txt');
    rnxname = input(sprintf('\nEnter name of RINEX observation file (%s default): ',defname),'s');
    if isempty(rnxname)
        rnxname = defname;
    end
end

fid = fopen(rnxname,'r');
clear defname

EOH = 0;
% Read header (but currently don't do anything with header data)
while (feof(fid) == 0) & (EOH == 0)
    line = fgetl(fid);
    if strcmp(line(61:65),'ION A')
        Alphas = line(1:60);
        Alphas = strrep(Alphas,'D','E');
        KlobucharCoefs.A = sscanf(Alphas,'%f');
    elseif strcmp(line(61:65),'ION B')
        Betas = line(1:60);
        Betas = strrep(Betas,'D','E');
        KlobucharCoefs.B = sscanf(Betas,'%f');
    end
    labela = line(61:length(line));
    label = '                    ';
    label(1:length(labela)) = labela;
    EOH = sum(label(1:13)=='END OF HEADER') >= 13;
end  % at end of loop, line pointer at first line of file body
clear EOH label labela line 

step = 1;

ephem = zeros(1000,38);

% Start loop through time
while feof(fid) == 0
   % Read first record from file (8 lines)
   line = fgetl(fid);  % 1prn, 2yy, 3mm, 4dd, 5hh, 6mm, 7ss, 8svbias, 9svdrift, 10svdriftrate
   line = strrep(line,'D','E');
   [line0,ct] = sscanf(line,'%f');
   line = strrep(line,'D','E');
   line = fgetl(fid);          % 11IODE, 12Crs, 13Delta n, 14M0 
   line = strrep(line,'D','E');
   [line1,ct] = sscanf(line,'%f');
   line = fgetl(fid);          % 15Cuc, 16e, 17Cus, 18roota
   line = strrep(line,'D','E');
   [line2,ct] = sscanf(line,'%f');
   line = fgetl(fid);          % 19Toe, 20Cic, 21OMEGA, 22Cis
   line = strrep(line,'D','E');
   [line3,ct] = sscanf(line,'%f');
   line = fgetl(fid);          % 23i, 24Crc, 25omega, 26OMEGADOT
   line = strrep(line,'D','E');
   [line4,ct] = sscanf(line,'%f');
   line = fgetl(fid);          % 27IDOT, 28codeL2, 29week, 30L2Pdata
   line = strrep(line,'D','E');
   [line5,ct] = sscanf(line,'%f');
   line = fgetl(fid);          % 31SVacc, 32SVhealth, 33Tgd, 34IODC
   line = strrep(line,'D','E');
   [line6,ct] = sscanf(line,'%f');
   line = fgetl(fid);          % 35transtime, 36fit, 37spare1, 38spare2
   line = strrep(line,'D','E');
   [line7,ct] = sscanf(line,'%f');

   % Checks if proper number of records were read in for first three lines
   % (known problem with pivot rinex files will cause an error here)
   if (length(line0) < 10) | (length(line1) < 4) | (length(line2) < 4)
       msg = sprintf('\nMissing data in first three lines of record %d,\n possibly caused by out-of-place header fields (from PiVoT)\n',step);
       fclose(fid);
       error(msg); 
   elseif length([line0',line1',line2',line3',line4',line5',line6',line7(1)']) < 35
       % Check to make sure first 35 records were read in correctly
       % **** PiVoT currently only outputs 35 records ****
       msg = sprintf('\nPossible missing data in record %d,\n\n',step);
       warning(msg);      
   end
   
   % Write record to matrix
   record_length = length([line0',line1',line2',line3',line4',line5',line6',line7']);
   ephem(step,1:record_length) = [line0',line1',line2',line3',line4',line5',line6',line7'];
   
   step = step+1;
   clear line line0 line1 line2 line3 line4 line5 line6 line7 ct
end

ephem = ephem(1:(step-1),:);

fclose(fid);
clear step msg fid rnxname

toc_gpstime = utc2gps(ephem(:,2:7),0);
utctime = sec2utc(utc2sec(ephem(:,2:7)));
time = gpst2sec(utc2gps(utctime));
prn = ephem(:,1);

% Map original ephem matrix (38,n) to new eph matrix (n,34)
eph(:,1)  = ephem(:,1);  % prn
eph(:,2)  = ones(size(ephem(:,1)));           % vflg
eph(:,3)  = ephem(:,35);  % TofXmission
eph(:,4)  = ephem(:,32);  % s1hlth
eph(:,5)  = ephem(:,28);  % codeL2
eph(:,6)  = ephem(:,29);  % wkn  **** This is redundant, set equal to toewk
eph(:,7)  = ephem(:,30);  % L2Pdata
eph(:,8)  = ephem(:,31);  % ura
eph(:,9)  = ephem(:,34);  % iodc
eph(:,10) = ephem(:,33);  % tgd
eph(:,11) = toc_gpstime(:,1);  % tocwk
eph(:,12) = toc_gpstime(:,2);  % tocsec
eph(:,13) = ephem(:,8);  % af0
eph(:,14) = ephem(:,9);  % af1
eph(:,15) = ephem(:,10);  % af2
eph(:,16) = ephem(:,11);  % iode
eph(:,17) = ephem(:,12);  % crs 
eph(:,18) = ephem(:,13);  % deltan
eph(:,19) = ephem(:,14);  % m0
eph(:,20) = ephem(:,15);  % cuc
eph(:,21) = ephem(:,16);  % ecc
eph(:,22) = ephem(:,17);  % cus
eph(:,23) = ephem(:,18);  % sqrta
eph(:,24) = ephem(:,29);  % toewk
eph(:,25) = ephem(:,19);  % toe
eph(:,26) = ephem(:,36);  % fti
eph(:,27) = ephem(:,20);  % cic
eph(:,28) = ephem(:,21);  % om0
eph(:,29) = ephem(:,22);  % cis
eph(:,30) = ephem(:,23);  % in0
eph(:,31) = ephem(:,24);  % crc
eph(:,32) = ephem(:,25);  % olc
eph(:,33) = ephem(:,26);  % omd
eph(:,34) = ephem(:,27);  % idot

% Check time of transmission week number
resetweek = abs(eph(:,25) - eph(:,3)) >= (SECONDS_IN_WEEK/2);
eph(:,6) = eph(:,6) - resetweek;

% Summary info to display
ind_first = min(find(time==min(time)));
ind_last = max(find(time==max(time)));
for j=1:32
    byprn(j) = sum(prn==ones(length(prn),1)*j);
end
fprintf('\nBroadcast Records Read:     %20d\n',length(prn));
fprintf('Number With Invalid Health:   %18d\n', sum(ephem(:,32)~=0));
fprintf('Earliest Epoch:            %s\n',calendar2datestr(utctime(ind_first,:)))
fprintf('Lastest Epoch:             %s\n',calendar2datestr(utctime(ind_last,:)))
fprintf('Records By PRN:\n');
fprintf('%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n',[1:16]);
fprintf('%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n',byprn(1:16));
fprintf('\n%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n',[17:32]);
fprintf('%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n',byprn(17:32));
