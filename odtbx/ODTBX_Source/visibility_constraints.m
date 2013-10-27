function [HVIS] = visibility_constraints(AntLB, options, health, Hvis_earth, Hvis_atm)

link_budget = getOdtbxOptions(options, 'linkbudget', []);

% Error check inputs
link_budget = linkbudget_default(link_budget, 'RecAcqThresh', 32 ); % dB-Hz, Receiver acquisition threshold
link_budget = linkbudget_default(link_budget, 'RecTrackThresh', link_budget.RecAcqThresh ); % dB-Hz, Receiver tracking threshold
link_budget = linkbudget_default(link_budget, 'DynamicTrackRange', 15 ); % dB
    % Dynamic tracking range of receiver, or maximum difference  
    % in power levels tracked simultaneously. If the difference
    % in snrs between two satellites is more that link_budget.DynamicTrackRange,
    % the weaker of the two will not be considered visible. 
link_budget = linkbudget_default(link_budget, 'AntennaPattern', {'sensysmeas_ant.txt','sensysmeas_ant.txt'});
    %  Specify antenna pattern for each antenna, existing antennas are:
    %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
    %     omni.txt                  - zero dB gain,  180 degree half beamwidth
    %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
    %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
    %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
num_ant = length(link_budget.AntennaPattern); %hasn't been tested for >4 antennas


TARGET_SIZE = size(AntLB{1}.Halpha_r,1);
                
% Set receiver antenna loop number
loop = max([1,num_ant]);

%% Initialize - Combine results across multiple antennas

d2r = pi/180;

HVIS = zeros(size(health'));
HCN0 = ones(size(health')).*-300;
HRP = ones(size(health')).*-300;
HVISdyn = zeros(size(health'));

% HVIScases is a cell array of matrices that contain
% visibility information for both antennas.  It is sized 2+(2*n) where n is
% the number of antennas on the satellite.
%HVIScases={'HVIS','HVISdyn'};
HVIScases=cell(1, 2*loop+2);

betar_gss = 180*d2r; %Aperature mask (half angle) used with rcv antennas in GSS simulator


% Compute visibility and other derived parameters
% These can be re-computed quickly and therefore are not saved

AntVis = cell(loop,1); % cell array of structs to hold visibility data
                      % for each antenna
for ANT=1:loop

    AntVis{ANT} = struct('Hvis_gss',[],'Hel',[],'Hvis',[],...
        'Hvis_CN0dyn',[],'Hvisdyn',[]);
    
    % Compute visibility for each antenna
    % Visibility for the GSS simulator, takes the form:
    % HXvis_gss = Health & Hvis_earth & (HXalpha_r <= betar_gss);
%    eval(sprintf('H%dvis_gss = health'' & Hvis_earth & (H%dalpha_r <= betar_gss);',ANT*ones(1,2)));
    AntVis{ANT}.Hvis_gss = health' & Hvis_earth & (AntLB{ANT}.Halpha_r <= betar_gss);

    % User elevation angle with respect to antenna boresite
    % Set non-existent/unhealthy SVs to -90 deg
%    eval(sprintf('H%del = (pi/2) - H%dalpha_r;',ANT*ones(1,2)));    % (nn,GPS_SIZE)
    AntVis{ANT}.Hel = (pi/2) - AntLB{ANT}.Halpha_r;    % (nn,GPS_SIZE)
%    eval(sprintf('H%del(~health'') = -pi/2;',ANT));
    AntVis{ANT}.Hel(~health'') = -pi/2;
    
    % Compute sat_pos in lla
    % (commented out because results were unused)
    %[sat_lla(1,:),sat_lla(2,:),sat_lla(3,:)] = ecef2LLA(sat_pos);

    % Compute visibility based on dynamic CN0 threshold (imposes limit on range of simultaneous C/No values)
    % Compute maximum CN0 for VISIBLE SV at each time step

    % Compute normal visibility for each antenna
    %finalstring = sprintf('H%dvis = health''.*Hvis_earth.*H%dvis_beta.*Hvis_atm.*H%dvis_CN0',ANT*ones(1,3));
    %   H1vis = Health.*Hvis_earth.*H1vis_beta.*Hvis_atm.*H1vis_CN0;   % (GPS_SIZE,mm)
    %eval([finalstring,';']);
    AntVis{ANT}.Hvis = health' .* Hvis_earth .* AntLB{ANT}.Hvis_beta .* ...
        Hvis_atm .* AntLB{ANT}.Hvis_CN0;   % (GPS_SIZE,mm)

%    eval(sprintf('maxCN0 = max(H%dvis.*H%dCN0);',ANT*ones(1,2)));
    maxCN0 = max(AntVis{ANT}.Hvis .* AntLB{ANT}.HCN0);
    
    % Compute the visibility based on the dynamic CN0 limit at each time step
    CN0_lim_dyn = ones(TARGET_SIZE,1)*max(link_budget.RecTrackThresh,(maxCN0 - link_budget.DynamicTrackRange));
    %eval(sprintf('H%dvis_CN0dyn = H%dCN0 >= CN0_lim_dyn;',ANT*ones(1,2)));
    AntVis{ANT}.Hvis_CN0dyn = AntLB{ANT}.HCN0 >= CN0_lim_dyn;
    %clear maxCN0 CN0_lim_dyn

    % Compute visibility for each antenna using dynamic tracking threshold
    %finalstring = sprintf('H%dvisdyn = health''.*Hvis_earth.*H%dvis_beta.*Hvis_atm.*H%dvis_CN0dyn',ANT*ones(1,3));
    %   H1visdyn = Health.*Hvis_earth.*H1vis_beta.*Hvis_atm.*H1vis_CN0dyn;   % (GPS_SIZE,mm)
    %eval([finalstring,';']);
    AntVis{ANT}.Hvisdyn = health' .* Hvis_earth .* AntLB{ANT}.Hvis_beta...
        .* Hvis_atm .* AntVis{ANT}.Hvis_CN0dyn;   % (GPS_SIZE,mm)

    %eval(sprintf('HVIS = HVIS | H%dvis;',ANT));
    HVIS = HVIS | AntVis{ANT}.Hvis;

    %eval(sprintf('HVISdyn = HVISdyn | H%dvis_CN0dyn;',ANT));
    HVISdyn = HVISdyn | AntVis{ANT}.Hvis_CN0dyn;
    
    %HVIScases=eval(sprintf('union(HVIScases, ''H%dvis'');',ANT));
    %HVIScases=eval(sprintf('union(HVIScases, ''AntVis{%d}.Hvis'');',ANT));
    %HVIScases = union(HVIScases, ['AntVis{',num2str(ANT),'}.Hvis']);
    HVIScases{2*ANT-1} = AntVis{ANT}.Hvis;
    
    %HVIScases=eval(sprintf('union(HVIScases, ''H%dvis_CN0dyn'');',ANT));
    %HVIScases=eval(sprintf('union(HVIScases, ''AntVis{%d}.Hvis_CN0dyn'');',ANT));
    %HVIScases = union(HVIScases, ['AntVis{',num2str(ANT),'}.Hvis_CN0dyn']);
    HVIScases{2*ANT} = AntVis{ANT}.Hvis_CN0dyn;

    % Compute composite RP, C/No across all simulated antennas (take the max C/No)
    %eval(sprintf('HRP = max(HRP,H%dRP);',ANT));
    HRP = max(HRP, AntLB{ANT}.HRP);
    
    % eval(sprintf('HCN0 = max(HCN0,H%dCN0);',ANT));
    HCN0 = max(HCN0, AntLB{ANT}.HCN0);

end
HVIScases{2*loop+1} = HVIS;
HVIScases{2*loop+2} = HVISdyn;

%% Apply Aquisition constraint to visibility (HCNO with HVIS and HVISdyn)
for HVc=1:length(HVIScases)    
    %VISCN0=eval(['HCN0.*' HVIScases{HVc}]);
    VISCN0=HCN0.*HVIScases{HVc};
    VISCN0acq=zeros(size(VISCN0));
    for n=1:size(VISCN0,1)
        Vacq=find(VISCN0(n,:)>=link_budget.RecAcqThresh);
        Vzero=union(find(VISCN0(n,:)==0),size(VISCN0,2)+1);
        if isempty(Vacq)
            VStartIndex=[];
        else
            VStartIndex=1;
        end
        while ~isempty(VStartIndex)
            VStart=Vacq(VStartIndex);
            VStop=Vzero(find(Vzero>=VStart,1,'first'))-1;
            VISCN0acq(n,VStart:VStop)=VISCN0(n,VStart:VStop);
            VStartIndex=VStartIndex+find((Vacq(VStartIndex+1:end)-Vacq(VStartIndex:end-1))>1,1,'first');
        end
    end
    %eval([HVIScases{HVc} '=(VISCN0acq~=0);']);
    HVIScases{HVc} = (VISCN0acq~=0);
    clear n VISCNO VISCN0acq Vacq Vzero VStartIndex VStart VStop
end


%HRR = sqrt(sum(sat_pos.^2));   % [1,nn]  magnitude of user vehicle position vector

HVIS = HVIS.*health';  % Make sure all non-existent SVs are not visible
%HCN0(~health') = -300;  % Set C/No for all non-existent SVs to -300
%HRP(~health') = -300;  % Set RP for all non-existent SVs to -300
end