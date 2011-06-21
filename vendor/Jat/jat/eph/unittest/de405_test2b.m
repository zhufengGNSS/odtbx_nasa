% Tests JAT DE405 position output vs SPICE.
% (Assumes the spice DE405 SPK kernel is loaded)
% Creates output files suitable for use by DE405Test.java and creates
% position plots for visual checking.

close all;

cspice_furnsh('naif0009.tls');
cspice_furnsh('de405.bsp');
cspice_furnsh('pck00009.tpc');
de = jat.eph.DE405;

% JAT DE405 args:
jb = {jat.eph.DE405_Body.MERCURY, jat.eph.DE405_Body.VENUS, ...
    jat.eph.DE405_Body.EARTH, jat.eph.DE405_Body.MARS, ...
    jat.eph.DE405_Body.JUPITER, jat.eph.DE405_Body.SATURN, ...
    jat.eph.DE405_Body.URANUS, jat.eph.DE405_Body.NEPTUNE, ...
    jat.eph.DE405_Body.PLUTO, ...
    jat.eph.DE405_Body.GEOCENTRIC_SUN, jat.eph.DE405_Body.GEOCENTRIC_MOON, ...
    jat.eph.DE405_Body.EM_BARY, ...
    jat.eph.DE405_Body.SUN, jat.eph.DE405_Body.MOON, ...
    jat.eph.DE405_Body.SOLAR_SYSTEM_BARY
    };

% SPICE args:
% Note, the JAT designations vary for the SPICE reference and observation 
% locations.  Some are the barycenter while others are the body centers!
planet = {'MERCURY BARYCENTER', 'VENUS BARYCENTER', ...
    'EARTH', 'MARS BARYCENTER', ...
    'JUPITER BARYCENTER', 'SATURN BARYCENTER', ...
    'URANUS BARYCENTER', 'NEPTUNE BARYCENTER', ...
    'PLUTO BARYCENTER', ...
    'SUN', 'MOON', ...
    'EARTH BARYCENTER' ...
    'SUN', 'MOON', ...
    'SOLAR SYSTEM BARYCENTER'};
obs = {'SOLAR SYSTEM BARYCENTER', 'SOLAR SYSTEM BARYCENTER', ...
    'SOLAR SYSTEM BARYCENTER', 'SOLAR SYSTEM BARYCENTER', ...
    'SOLAR SYSTEM BARYCENTER', 'SOLAR SYSTEM BARYCENTER', ...
    'SOLAR SYSTEM BARYCENTER', 'SOLAR SYSTEM BARYCENTER', ...
    'SOLAR SYSTEM BARYCENTER', ...
    'EARTH', 'EARTH', ...
    'SOLAR SYSTEM BARYCENTER', ...
    'SOLAR SYSTEM BARYCENTER', 'SOLAR SYSTEM BARYCENTER', ...
    'SOLAR SYSTEM BARYCENTER'};

% rough comparative orbit radius, just for sanity check:
jat_au = 149597870000.0; % IAY 1976 (m)
orbrad = [0.39, .72, ...
    1.0, 1.52, ...
    5.20, 9.54, ...
    19.22, 30.06, ...
    (29.7+49.3)/2, ...
    1.0, 0.00257, ...
    1.0, ...
    0.0001, 1.0, ...
    0] * jat_au/1000; %km

maxdifftbl = [];

for p = 1:length(planet)

    fd = fopen(['de405_' planet{p} '_data.dat'],'w+t');
    if fd == -1
        disp('Can''t open file');
    end
    
    years = 1:10;
    months = 1:12;
    days = 1:4:28;
    sz = length(years)*length(months)*length(days);

    t_tt = zeros(sz,1);
    r_sun = cell(sz,1);
    r_sun_mag = zeros(sz,1);
    sstate = cell(sz,1);
    r_mag_spice = r_sun_mag;
    i = 1;
    
    % write the file header
    fprintf(fd,'# SPICE DE405 data for %s \n',[planet{p} ' from '  obs{p}]);
    fprintf(fd,'# For the period %d-1-1-12:00:00 to %d-%d-%d-12:00:00\n',...
        2000, 2000+max(years)-1,max(months),max(days));
    fprintf(fd,'# 7 items per row: time (TT, mjd), XYZ position (m), XYZ velocity (m/s)\n');
    fprintf(fd,'# For use with JAT: jat.eph.unittest.DE405Test.readTargetPosVels().\n');
    fprintf(fd,'# NOTE: The majority of the JAT-SPICE DE405 differences in this data\n');
    fprintf(fd,'#       are due to an unexplained -80 microsec (-8.046852e-5 sec) time difference\n');
    fprintf(fd,'#       between JAT and SPICE.  When this difference is artifically compensated for,\n');
    fprintf(fd,'#       the position results millimeters or less (at numerical limits) for all bodies.\n');
    
    t_mjd_tt_j2k = jat.spacetime.TimeUtils.MJD_J2000; % TT, MJD, J2K epoch
    t_mjd_tdb_j2k = jat.spacetime.TimeUtils.TTtoTDB(t_mjd_tt_j2k); % TDB, MJD
    t_jd_tdb_j2k = jat.spacetime.TimeUtils.MJDtoJD(t_mjd_tdb_j2k); % TDB, JD, J2K epoch
    
    for y = years
        for m = months
            for d = days
                t_tt(i) = datenum(2000+y-1, m, d, 12, 0, 0); % TT
                t_mjd_tt = matlabTime2MJD(t_tt(i)); % TT, MJD

                %JAT DE405:
                r_sun{i} = de.get_planet_posvel(jb{p}, t_mjd_tt); % m
                r_sun_mag(i) = sqrt(sum(r_sun{i}.toDouble().*r_sun{i}.toDouble()))/1000; % km

                % SPICE DE405:
                t_mjd_tdb = jat.spacetime.TimeUtils.TTtoTDB(t_mjd_tt); % TDB, MJD
                t_jd_tdb = jat.spacetime.TimeUtils.MJDtoJD(t_mjd_tdb); % TDB, JD
                
                % SPICE ephemeris time is seconds past the J2000 epoch
                % (TDB).
                % NOTE: There is an unexplained 80microsec time difference
                % between the SPICE results and the JAT results.
                et = (t_jd_tdb-t_jd_tdb_j2k)*jat.spacetime.TimeUtils.days2sec;
                sstate{i} = cspice_spkezr( planet{p}, et, 'J2000', 'NONE', obs{p} ); % km
                r_mag_spice(i) = sqrt(sum(sstate{i}(1:3).*sstate{i}(1:3)));
                
                % write SPICE data to file:
                % time  (mjd) r_x r_y r_z (m) v_x v_y v_z (m/s) 
                sd = sstate{i}*1000;
                fprintf(fd,'%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n', ...
                    t_mjd_tt, sd(1), sd(2), sd(3),sd(4), sd(5), sd(6));
                
                i = i + 1;
            end
        end   
    end
    
    fclose(fd);

    jat_au = 149597870000.0; % IAY 1976

    jatx = zeros(length(r_sun),1);
    jaty = jatx;
    jatz = jatx;
    spicex = jatx;
    spicey = jatx;
    spicez = jatx;
    jatvx = jatx;
    jatvy = jatx;
    jatvz = jatx;
    spicevx = jatx;
    spicevy = jatx;
    spicevz = jatx;
    for i = 1:length(r_sun)
        r = (r_sun{i}.toDouble)';
        jatx(i) = r(1)/1000;
        jaty(i) = r(2)/1000;
        jatz(i) = r(3)/1000;
        jatvx(i) = r(4)/1000;
        jatvy(i) = r(5)/1000;
        jatvz(i) = r(6)/1000;
        spicex(i) = sstate{i}(1);
        spicey(i) = sstate{i}(2);
        spicez(i) = sstate{i}(3);
        spicevx(i) = sstate{i}(4);
        spicevy(i) = sstate{i}(5);
        spicevz(i) = sstate{i}(6);
    end
    
    
    figure;
    subplot(2,1,1)
    plot(t_tt-t_tt(1),r_mag_spice,'ro',t_tt-t_tt(1),r_sun_mag,'b.-');
    ylabel('range (km)');
    xlabel('days');
    title([planet{p} ' from '  obs{p}]);
    legend('JAT DE405','SPICE');
    subplot(2,1,2)
    hold on;
    plot(spicex, spicey, 'ro');
    plot(jatx, jaty, 'b.-');
    orx = orbrad(p)*cos(0:.1:(2*pi));
    ory = orbrad(p)*sin(0:.1:(2*pi));
    plot(orx, ory, 'g-');
    axis square;
    xlabel('X (km)');
    ylabel('Y (km)');
    hold off;

    figure;
    subplot(3,1,1);
    plot(t_tt-t_tt(1),spicex,'ro',t_tt-t_tt(1),jatx,'b.-');
    ylabel('J2K X (km)');
    title([planet{p} ' from '  obs{p}]);
    subplot(3,1,2);
    plot(t_tt-t_tt(1),spicey,'ro',t_tt-t_tt(1),jaty,'b.-');
    ylabel('J2K Y (km)');
    subplot(3,1,3);
    plot(t_tt-t_tt(1),spicez,'ro',t_tt-t_tt(1),jatz,'b.-');
    ylabel('J2K Y (km)');

    dx = abs(jatx - spicex);
    dy = abs(jaty - spicey);
    dz = abs(jatz - spicez);
    dvx = abs(jatvx - spicevx);
    dvy = abs(jatvy - spicevy);
    dvz = abs(jatvz - spicevz);
    fprintf(1,'%s max X difference mag: %f (km)\n',[planet{p} ' from '  obs{p}],max(dx));
    fprintf(1,'%s max Y difference mag: %f (km)\n',[planet{p} ' from '  obs{p}],max(dy));
    fprintf(1,'%s max Z difference mag: %f (km)\n',[planet{p} ' from '  obs{p}],max(dz));
    fprintf(1,'%s max VX difference mag: %f (km)\n',[planet{p} ' from '  obs{p}],max(dvx));
    fprintf(1,'%s max VY difference mag: %f (km)\n',[planet{p} ' from '  obs{p}],max(dvy));
    fprintf(1,'%s max VZ difference mag: %f (km)\n',[planet{p} ' from '  obs{p}],max(dvz));
    maxdifftbl(1+end,:) = [max(dx) max(dy) max(dz) max(dvx) max(dvy) max(dvz)];
    
    deltatx = (jatx-spicex)./spicevx;
    deltaty = (jaty-spicey)./spicevy;
    deltatz = (jatz-spicez)./spicevz;
    fprintf(1,'%s mean delta-T X difference: %e (s)\n',[planet{p} ' from '  obs{p}],mean(deltatx));
    fprintf(1,'%s mean delta-T Y difference: %e (s)\n',[planet{p} ' from '  obs{p}],mean(deltaty));
    fprintf(1,'%s mean delta-T Z difference: %e (s)\n',[planet{p} ' from '  obs{p}],mean(deltatz));
    
    figure;
    subplot(4,1,1);
    plot(t_tt-t_tt(1),jatx-spicex,'k.');
    ylabel('Diff X (km)');
    title(['JAT-SPICE ' planet{p} ' from '  obs{p}]);
    subplot(4,1,2);
    plot(t_tt-t_tt(1),jaty-spicey,'k.');
    ylabel('Diff Y (km)');
    subplot(4,1,3);
    plot(t_tt-t_tt(1),jatz-spicez,'k.');
    ylabel('Diff Y (km)');
    subplot(4,1,4);
    plot(t_tt-t_tt(1),deltatx,'k.', t_tt-t_tt(1),deltaty,'b+',t_tt-t_tt(1),deltatz,'ro');
    ylabel('time diff (s)');

end

% max tolernace data to the screen and the file
fprintf(1,'\nMaximum differences for comparison (m, m/s):\n');
fprintf(1,'0 %e %e %e %e %e %e\n',(maxdifftbl*1000)');
fd = fopen('de405_TOL_data.dat','w+t');
if fd == -1
    disp('Can''t open file');
end
fprintf(fd,'# DE405 comparison tolerance data for each planet and in each axis\n');
fprintf(fd,'# time (empty), position tolerance (m), velocity tolerance (m/s)\n');
fprintf(fd,'# For use with JAT: jat.eph.unittest.DE405Test.readTargetPosVels().\n');
fprintf(fd,'# NOTE: The majority of the JAT-SPICE DE405 differences in this data\n');
fprintf(fd,'#       are due to an unexplained -80 microsec (-8.046852e-5 sec) time difference\n');
fprintf(fd,'#       between JAT and SPICE.  When this difference is artifically compensated for,\n');
fprintf(fd,'#       the position results millimeters or less (at numerical limits) for all bodies.\n');
fprintf(fd,'0 %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n',(maxdifftbl*1000)');
fclose(fd);
