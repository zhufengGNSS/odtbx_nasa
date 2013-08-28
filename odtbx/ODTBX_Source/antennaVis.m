function [AntVis,HVIScases_Hvis, HVIScases_Hvis_CNOdyn] = antennaVis(health, Hvis_earth, AntLB, betar_gss, ...
    Hvis_atm, GPS_SIZE, CN0_lim, dyn_range, HVIS, HVISdyn, ANT, HRP, HCN0)

    AntVis = struct('Hvis_gss',[],'Hel',[],'Hvis',[],...
        'Hvis_CN0dyn',[],'Hvisdyn',[]);
    
    % Compute visibility for each antenna
    % Visibility for the GSS simulator, takes the form:
    % HXvis_gss = Health & Hvis_earth & (HXalpha_r <= betar_gss);
%    eval(sprintf('H%dvis_gss = health'' & Hvis_earth & (H%dalpha_r <= betar_gss);',ANT*ones(1,2)));
    AntVis.Hvis_gss = health' & Hvis_earth & (AntLB.Halpha_r <= betar_gss);

    % User elevation angle with respect to antenna boresite
    % Set non-existent/unhealthy SVs to -90 deg
%    eval(sprintf('H%del = (pi/2) - H%dalpha_r;',ANT*ones(1,2)));    % (nn,GPS_SIZE)
    AntVis.Hel = (pi/2) - AntLB.Halpha_r;    % (nn,GPS_SIZE)
%    eval(sprintf('H%del(~health'') = -pi/2;',ANT));
    AntVis.Hel(~health'') = -pi/2;
    
    % Compute sat_pos in lla
    % (commented out because results were unused)
    %[sat_lla(1,:),sat_lla(2,:),sat_lla(3,:)] = ecef2LLA(sat_pos);

    % Compute visibility based on dynamic CN0 threshold (imposes limit on range of simultaneous C/No values)
    % Compute maximum CN0 for VISIBLE SV at each time step

    % Compute normal visibility for each antenna
    %finalstring = sprintf('H%dvis = health''.*Hvis_earth.*H%dvis_beta.*Hvis_atm.*H%dvis_CN0',ANT*ones(1,3));
    %   H1vis = Health.*Hvis_earth.*H1vis_beta.*Hvis_atm.*H1vis_CN0;   % (GPS_SIZE,mm)
    %eval([finalstring,';']);
    AntVis.Hvis = health' .* Hvis_earth .* AntLB.Hvis_beta .* ...
        Hvis_atm .* AntLB.Hvis_CN0;   % (GPS_SIZE,mm)

%    eval(sprintf('maxCN0 = max(H%dvis.*H%dCN0);',ANT*ones(1,2)));
    maxCN0 = max(AntVis.Hvis .* AntLB.HCN0);
    
    % Compute the visibility based on the dynamic CN0 limit at each time step
    CN0_lim_dyn = ones(GPS_SIZE,1)*max(CN0_lim,(maxCN0 - dyn_range));
    %eval(sprintf('H%dvis_CN0dyn = H%dCN0 >= CN0_lim_dyn;',ANT*ones(1,2)));
    AntVis.Hvis_CN0dyn = AntLB.HCN0 >= CN0_lim_dyn;
    %clear maxCN0 CN0_lim_dyn

    % Compute visibility for each antenna using dynamic tracking threshold
    %finalstring = sprintf('H%dvisdyn = health''.*Hvis_earth.*H%dvis_beta.*Hvis_atm.*H%dvis_CN0dyn',ANT*ones(1,3));
    %   H1visdyn = Health.*Hvis_earth.*H1vis_beta.*Hvis_atm.*H1vis_CN0dyn;   % (GPS_SIZE,mm)
    %eval([finalstring,';']);
    AntVis.Hvisdyn = health' .* Hvis_earth .* AntLB.Hvis_beta...
        .* Hvis_atm .* AntVis.Hvis_CN0dyn;   % (GPS_SIZE,mm)

    %eval(sprintf('HVIS = HVIS | H%dvis;',ANT));
    HVIS = HVIS | AntVis.Hvis;

    %eval(sprintf('HVISdyn = HVISdyn | H%dvis_CN0dyn;',ANT));
    HVISdyn = HVISdyn | AntVis.Hvis_CN0dyn;
    
    %HVIScases=eval(sprintf('union(HVIScases, ''H%dvis'');',ANT));
    %HVIScases=eval(sprintf('union(HVIScases, ''AntVis{%d}.Hvis'');',ANT));
    %HVIScases = union(HVIScases, ['AntVis{',num2str(ANT),'}.Hvis']);
    HVIScases_Hvis = AntVis.Hvis;
    
    %HVIScases=eval(sprintf('union(HVIScases, ''H%dvis_CN0dyn'');',ANT));
    %HVIScases=eval(sprintf('union(HVIScases, ''AntVis{%d}.Hvis_CN0dyn'');',ANT));
    %HVIScases = union(HVIScases, ['AntVis{',num2str(ANT),'}.Hvis_CN0dyn']);
    HVIScases_Hvis_CNOdyn = AntVis.Hvis_CN0dyn;

    % Compute composite RP, C/No across all simulated antennas (take the max C/No)
    %eval(sprintf('HRP = max(HRP,H%dRP);',ANT));
    HRP = max(HRP, AntLB.HRP);
    
    % eval(sprintf('HCN0 = max(HCN0,H%dCN0);',ANT));
    HCN0 = max(HCN0, AntLB.HCN0);

end