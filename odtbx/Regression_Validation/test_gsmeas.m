classdef test_gsmeas < matlab.unittest.TestCase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        epoch
        t
        x
        gsList
        gsID
        gsECEF
        measOptions_noLinkBudget
        measOptions_LinkBudget
        measOptions_2D_Antenna_LinkBudget
    end
    
    methods (TestMethodSetup)
        function create_epoch(testCase)
            testCase.epoch = datenum('Jan 1 2006');
        end
        
        function create_tandx(testCase)
            testCase.t = [0 60 120 180];

            testCase.x = [
                5158.712338     -1514.921648      3777.724544       2.754779       9.375862      -0.001963;
                5310.829717      -949.016866      3768.064835       2.313986       9.479618      -0.319662;
                5436.235143      -378.334492      3739.451984       1.865419       9.535073      -0.633026;
                5534.648919       194.234583      3692.270615       1.415291       9.542710      -0.937961]';
        end
        
        function create_gsList_ID_ECEF(testCase)
            testCase.gsList = createGroundStationList();
            testCase.gsID   = {  
                        'DS12'
                        'DS16' 
                        'DS46'
                        'DS66' };
            nGS    = length(testCase.gsID);
            testCase.gsECEF = zeros(3,nGS);
            for n=1:nGS
                testCase.gsECEF(:,n) = getGroundStationInfo(...
                    testCase.gsList,testCase.gsID{n},'ecefPosition',testCase.epoch);
            end
        end
        
        function create_measOptions_noLinkBudget(testCase)
            testCase.measOptions_noLinkBudget = odtbxOptions('measurement');
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'epoch',testCase.epoch);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useRange', true);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useRangeRate', true);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useDoppler', false);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'frequencyTransmit', JATConstant('L1Frequency'));
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'rangeType','2way');
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'gsElevationConstraint',0);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'gsECEF',testCase.gsECEF);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'gsID', testCase.gsID);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useLightTime',false);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useIonosphere',false);
            testCase.measOptions_noLinkBudget = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useTroposphere',false);
        end
        
        function create_measOptions_LinkBudget(testCase)
            d2r = pi/180;
            %% Link budget information
            % Parameters specific to link budget
            link_budget.AntennaPattern    = {'omni.txt'};           % Specify receive antenna pattern for each antenna
                %  Specify antenna pattern for each antenna, existing antennas are:
                %     sensysmeas_ant.txt        - hemi antenna, 4 dB peak gain, 157 degree half beamwidth
                %     omni.txt                  - zero dB gain,  180 degree half beamwidth
                %     trimblepatch_ant.txt      - hemi antenna, 4.5 dB gain, 90 deg half beamwidth
                %     ballhybrid_10db_60deg.txt - high gain, 10 db peak gain, 60 degree half-beamwidth
                %     ao40_hga_measured_10db.txt- another 10 dB HGA with 90 deg beamwidth
            link_budget.RXAntennaMask     = 180*d2r;     % Cut off angle for the receive antenna
            link_budget.AtmosphereMask    = 50;                      % Mask altitude, km
            link_budget.NoiseTemp         = 300;                    % System noise temp of receive antenna (K)
                % System noise temp [K], space pointing antenna = 290
                % System noise temp [K], earth pointing antenna = 300
            link_budget.AtmAttenuation    = 0.0;                    % Attenuation due to atmosphere, should be negative (dB)
            link_budget.TransPowerLevel   = 2;                      % Transmitter power level (1=min, 2=typical, 3=max)
            link_budget.TransPowerOffset  = 0;                      % Global transmitter power offset (dB)
            link_budget.TXAntennaMask     = 90 * d2r;               % Cut off angle for the transmit antenna (rad)
                %  The actual mask used is the lesser of this mask and the limit of the defined pattern
                %  Note:  mask = 70 deg includes entire defined pattern
                %         mask = 42 deg includes only main and first side lobes
                %         mask = 26 deg includes only main lobe
            link_budget.ReceiverNoise     = -3;                     % Noise figure of receiver/LNA (dB)
            link_budget.RecConversionLoss = 0;                      % Receiver implementation, A/D conversion losses (dB)
            link_budget.SystemLoss        = 0;                      % System losses, in front of LNA (dB)
            link_budget.LNAGain           = 40;                     % Gain provided by the LNA (dB)
            link_budget.CableLoss         = -2;                     % Cable losses after LNA (dB)
            link_budget.RecAcqThresh      = 34;                     % Receiver acquisition threshold (dB-Hz)
            link_budget.RecTrackThresh    = 34;                     % Receiver tracking threshold (dB-Hz)
            link_budget.DynamicTrackRange = 15;                     % Receiver dynamic range (dB)
            %link_budget.RXpattern = 'ao40_hga_measured_10db.txt';
            link_budget.TXpattern = 'ao40_hga_measured_10db.txt';
            link_budget.TX_AntennaPointing= 1; % 1 for zenith pointing, -1 for nadir pointing
            testCase.measOptions_LinkBudget = testCase.measOptions_noLinkBudget;
            testCase.measOptions_LinkBudget = setOdtbxOptions(testCase.measOptions_LinkBudget, 'linkbudget', link_budget);
        end
        
        function create_measOptions_2D_Antenna_LinkBudget(testCase)
            testCase.measOptions_2D_Antenna_LinkBudget = testCase.measOptions_LinkBudget;
            link_budget = getOdtbxOptions(testCase.measOptions_2D_Antenna_LinkBudget,'linkbudget');
            link_budget.AntennaPattern    = {'sensysmeas_ant_2D.txt'};
            link_budget.TXpattern = 'sensysmeas_ant_2D.txt';
            testCase.measOptions_2D_Antenna_LinkBudget = setOdtbxOptions(testCase.measOptions_2D_Antenna_LinkBudget, 'linkbudget', link_budget);
        end
            
    end
    
    methods (Test)
        function test_noLinkBudget(testCase)
            % run gsmeas with no LinkBudget and range and rangerate
            % measurements
            
            % Run gsmeas
            [y,H,R] = gsmeas(testCase.t, testCase.x, testCase.measOptions_noLinkBudget);
            % Check y
            y_expected = 1.0e3*[
               0.199032532782455   0.591898264330166   1.131011079112153   1.678243829361832
               0.000001650017813   0.008741270829023   0.009100998874876   0.009120816998485
               0.199213890231063   0.597831760859305   1.137200002369821   1.684481624233693
               0.000293647613281   0.008751923697204   0.009102497539121   0.009121215711152
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
                0.783398512203442        -0.229870841447615         0.577447978028843                         0                         0                         0
                0.0132960278492823        0.0452769539688087     -1.42643643439925e-05         0.783398512203442        -0.229870841447615         0.577447978028843
                0.804421968278982        -0.203607181003063         0.558076529513788                         0                         0                         0
                0.0121027687146992        0.0455355407620577     -0.000832090074627532         0.804421968278982        -0.203607181003063         0.558076529513788
                0.653140372345824         0.317407825572411         0.687502673651057                         0                         0                         0
               -5.56343720182588e-05       0.00075071336040374     -0.000293737389887119         0.653140372345824         0.317407825572411         0.687502673651057
                0.666193352210258        -0.744728490074234       -0.0395713727663203                         0                         0                         0
                0.00075052704805523       0.00067261244066317     -2.31863897485416e-05         0.666193352210258        -0.744728490074234       -0.0395713727663203];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);

        end
        
        function test_noLinkBudget_no_gsECEF(testCase)
            % run gsmeas with no LinkBudget and range and rangerate
            % measurements
            measOptions = setOdtbxOptions(testCase.measOptions_noLinkBudget,'gsECEF',[]);
            % Run gsmeas
            [y,H,R] = gsmeas(testCase.t, testCase.x, measOptions);
            % Check y
            y_expected = 1.0e3*[
               0.199032532782455   0.591898264330166   1.131011079112153   1.678243829361832
               0.000001650017813   0.008741270829023   0.009100998874876   0.009120816998485
               0.199213890231063   0.597831760859305   1.137200002369821   1.684481624233693
               0.000293647613281   0.008751923697204   0.009102497539121   0.009121215711152
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
                0.783398512203442        -0.229870841447615         0.577447978028843                         0                         0                         0
                0.0132960278492823        0.0452769539688087     -1.42643643439925e-05         0.783398512203442        -0.229870841447615         0.577447978028843
                0.804421968278982        -0.203607181003063         0.558076529513788                         0                         0                         0
                0.0121027687146992        0.0455355407620577     -0.000832090074627532         0.804421968278982        -0.203607181003063         0.558076529513788
                0.653140372345824         0.317407825572411         0.687502673651057                         0                         0                         0
               -5.56343720182588e-05       0.00075071336040374     -0.000293737389887119         0.653140372345824         0.317407825572411         0.687502673651057
                0.666193352210258        -0.744728490074234       -0.0395713727663203                         0                         0                         0
                0.00075052704805523       0.00067261244066317     -2.31863897485416e-05         0.666193352210258        -0.744728490074234       -0.0395713727663203];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);

        end
        
        function test_noLinkBudget_withIonosphere(testCase)
            % run gsmeas with no LinkBudget and range and rangerate
            % measurements with Ionosphere
            measOptions = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useIonosphere',true);
            % Run gsmeas
            [y,H,R] = gsmeas(testCase.t, testCase.x, measOptions);
            % Check y
            y_expected = 1.0e3*[
               0.199033819191285   0.591901226220171   1.131014867419668   1.678247840232908
               0.000004389458457   0.008742330388856   0.009101226998172   0.009120816898797
               0.199215177572410   0.597834740268890   1.137203796624833   1.684485636694389
               0.000296387199687   0.008752966980883   0.009102721129734   0.009121215614482
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
                0.783398512203442        -0.229870841447615         0.577447978028843                         0                         0                         0
                0.0132960278492823        0.0452769539688087     -1.42643643439925e-05         0.783398512203442        -0.229870841447615         0.577447978028843
                0.804421968278982        -0.203607181003063         0.558076529513788                         0                         0                         0
                0.0121027687146992        0.0455355407620577     -0.000832090074627532         0.804421968278982        -0.203607181003063         0.558076529513788
                0.653140372345824         0.317407825572411         0.687502673651057                         0                         0                         0
               -5.56343720182588e-05       0.00075071336040374     -0.000293737389887119         0.653140372345824         0.317407825572411         0.687502673651057
                0.666193352210258        -0.744728490074234       -0.0395713727663203                         0                         0                         0
                0.00075052704805523       0.00067261244066317     -2.31863897485416e-05         0.666193352210258        -0.744728490074234       -0.0395713727663203];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);

        end
        
        function test_noLinkBudget_withTroposphere(testCase)
            % run gsmeas with no LinkBudget and range and rangerate
            % measurements with Troposphere
            measOptions = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useTroposphere',true);
            % Run gsmeas
            [y,H,R] = gsmeas(testCase.t, testCase.x, measOptions);
            % Check y
            y_expected = 1.0e3*[
               0.199034774665485   0.591905388372341   1.131028024139526   1.678282955094066
               0.000001881114885   0.008742662315549   0.009104226309806   0.009131699208765
               0.199216133918880   0.597838974945844   1.137217146379245   1.684521380162365
               0.000293908839859   0.008753327094060   0.009105774051251   0.009132365545611
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
                0.783398512203442        -0.229870841447615         0.577447978028843                         0                         0                         0
                0.0132960278492823        0.0452769539688087     -1.42643643439925e-05         0.783398512203442        -0.229870841447615         0.577447978028843
                0.804421968278982        -0.203607181003063         0.558076529513788                         0                         0                         0
                0.0121027687146992        0.0455355407620577     -0.000832090074627532         0.804421968278982        -0.203607181003063         0.558076529513788
                0.653140372345824         0.317407825572411         0.687502673651057                         0                         0                         0
               -5.56343720182588e-05       0.00075071336040374     -0.000293737389887119         0.653140372345824         0.317407825572411         0.687502673651057
                0.666193352210258        -0.744728490074234       -0.0395713727663203                         0                         0                         0
                0.00075052704805523       0.00067261244066317     -2.31863897485416e-05         0.666193352210258        -0.744728490074234       -0.0395713727663203];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);

        end
        
        function test_noLinkBudget_withLightTime(testCase)
            % run gsmeas with no LinkBudget and range and rangerate
            % measurements using LightTime
            measOptions = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useLightTime',true);
            % Run gsmeas
            [y,H,R] = gsmeas(testCase.t, testCase.x, measOptions);
            % Check y
            y_expected = 1.0e3*[
               0.199032532782620   0.591898264309697   1.131011079192571   1.678243829364543
               0.000001650017817   0.008741270828982   0.009100998874894   0.009120816998485
               0.199213890231232   0.597831760838404   1.137200002451057   1.684481624236411
               0.000293647613285   0.008751923697164   0.009102497539139   0.009121215711152
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
               0.783398869395764  -0.229869625116898   0.577447977637861                   0                   0                   0
               0.013295980949963   0.045276967708025  -0.000014298864235   0.783398869395764  -0.229869625116898   0.577447977637861
               0.804422293995511  -0.203605957446998   0.558076506417436                   0                   0                   0
               0.012102719707128   0.045535549256131  -0.000832123267837   0.804422293995511  -0.203605957446998   0.558076506417436
               0.653141152169440   0.317407173316794   0.687502233938057                   0                   0                   0
              -0.000055634404980   0.000750713359407  -0.000293736770650   0.653141152169440   0.317407173316794   0.687502233938057
               0.666192635083809  -0.744729130164298  -0.039571399328553                   0                   0                   0
               0.000750527756830   0.000672611778259  -0.000023186456647   0.666192635083809  -0.744729130164298  -0.039571399328553];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);

        end
        
        function test_noLinkBudget_withLightTime_andIonosphere(testCase)
            % run gsmeas with no LinkBudget and range and rangerate
            % measurements using LightTime and Ionosphere
            measOptions = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useLightTime',true);
            measOptions = setOdtbxOptions(measOptions,'useIonosphere',true);
            % Run gsmeas
            [y,H,R] = gsmeas(testCase.t, testCase.x, measOptions);
            % Check y
            y_expected = 1.0e3*[
               0.199033819191449   0.591901226199701   1.131014867500087   1.678247840235619
               0.000004389458462   0.008742330388815   0.009101226998190   0.009120816898797
               0.199215177572579   0.597834740247988   1.137203796706070   1.684485636697107
               0.000296387199692   0.008752966980842   0.009102721129752   0.009121215614481
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
               0.783398869398071  -0.229869625109037   0.577447977637861                   0                   0                   0
               0.013295980949660   0.045276967708114  -0.000014298864458   0.783398869398071  -0.229869625109037   0.577447977637861
               0.804422293997615  -0.203605957439092   0.558076506417287                   0                   0                   0
               0.012102719706811   0.045535549256186  -0.000832123268051   0.804422293997615  -0.203605957439092   0.558076506417287
               0.653141152170607   0.317407173315818   0.687502233937399                   0                   0                   0
              -0.000055634404980   0.000750713359407  -0.000293736770649   0.653141152170607   0.317407173315818   0.687502233937399
               0.666192635083651  -0.744729130164438  -0.039571399328559                   0                   0                   0
               0.000750527756831   0.000672611778259  -0.000023186456647   0.666192635083651  -0.744729130164438  -0.039571399328559];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);

        end
        
        function test_noLinkBudget_withLightTime_andTroposphere(testCase)
            % run gsmeas with no LinkBudget and range and rangerate
            % measurements using LightTime and Troposphere
            measOptions = setOdtbxOptions(testCase.measOptions_noLinkBudget,'useLightTime',true);
            measOptions = setOdtbxOptions(measOptions,'useTroposphere',true);
            % Run gsmeas
            [y,H,R] = gsmeas(testCase.t, testCase.x, measOptions);
            % Check y
            y_expected = 1.0e3*[
               0.199034774665650   0.591905388351870   1.131028024219949   1.678282955096779
               0.000001881114890   0.008742662315508   0.009104226309825   0.009131699208765
               0.199216133919049   0.597838974924941   1.137217146460486   1.684521380165084
               0.000293908839863   0.008753327094019   0.009105774051269   0.009132365545611
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
               0.783398869399785  -0.229869625103198   0.577447977637859                   0                   0                   0
               0.013295980949435   0.045276967708180  -0.000014298864624   0.783398869399785  -0.229869625103198   0.577447977637859
               0.804422293999178  -0.203605957433218   0.558076506417178                   0                   0                   0
               0.012102719706576   0.045535549256227  -0.000832123268211   0.804422293999178  -0.203605957433218   0.558076506417178
               0.653142296215773   0.317406216372359   0.687501588872470                   0                   0                   0
              -0.000055634453305   0.000750713357960  -0.000293735862169   0.653142296215773   0.317406216372359   0.687501588872470
               0.666191310868355  -0.744730312123536  -0.039571448392541                   0                   0                   0
               0.000750529065631   0.000672610555102  -0.000023186580184   0.666191310868355  -0.744730312123536  -0.039571448392541];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);

        end
        
        function test_LinkBudget(testCase)
            % run gsmeas with a LinkBudget and range and rangerate
            % measurements
            
            % Run gsmeas
            [y,H,R,AntLB] = testCase.verifyWarning(@()gsmeas(testCase.t, testCase.x, testCase.measOptions_LinkBudget),...
                'ODTBX:GSMEAS:noBodyQuat');
            % Check y
            y_expected = 1.0e3*[
               0.199032532782455   0.591898264330166   1.131011079112153                 NaN
               0.000001650017813   0.008741270829023   0.009100998874876                 NaN
               0.199213890231063   0.597831760859305   1.137200002369821                 NaN
               0.000293647613281   0.008751923697204   0.009102497539121                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
                0.783398512203442        -0.229870841447615         0.577447978028843                         0                         0                         0
                0.0132960278492823        0.0452769539688087     -1.42643643439925e-05         0.783398512203442        -0.229870841447615         0.577447978028843
                0.804421968278982        -0.203607181003063         0.558076529513788                         0                         0                         0
                0.0121027687146992        0.0455355407620577     -0.000832090074627532         0.804421968278982        -0.203607181003063         0.558076529513788
                0.653140372345824         0.317407825572411         0.687502673651057                         0                         0                         0
               -5.56343720182588e-05       0.00075071336040374     -0.000293737389887119         0.653140372345824         0.317407825572411         0.687502673651057
                0.666193352210258        -0.744728490074234       -0.0395713727663203                         0                         0                         0
                0.00075052704805523       0.00067261244066317     -2.31863897485416e-05         0.666193352210258        -0.744728490074234       -0.0395713727663203];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.Halpha_r
            Halpha_r_expected = [
           0.003115704800299   1.168260941956579   1.275128535884713   1.274785638792920
           0.037196767999882   1.171438738118571   1.275774240924258   1.274913097277424
           0.582968622078663   0.545602755137455   0.508188440538280   0.471279381977418
           0.833955946910231   0.862165276398346   0.887662419802437   0.909884804974218];
            ABSTOL = 10*eps(max(max(abs(Halpha_r_expected))));
            testCase.verifyEqual(AntLB{1}.Halpha_r,Halpha_r_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HCN0
            HCN0_expected = 1.0e2*[
           0.684544605279052   0.427418255761894   0.366928073181581   0.332256035399327
           0.683590650852834   0.426404347525161   0.366437253323645   0.331934315096683
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HCN0_expected))));
            testCase.verifyEqual(AntLB{1}.HCN0,HCN0_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAd
            HAd_expected = 1.0e2*[
          -1.423741917056861  -1.518406516420129  -1.574650474966648  -1.608928114934059
          -1.423821026412351  -1.519272899922783  -1.575124473513245  -1.609250358666021
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAd_expected))));
            testCase.verifyEqual(AntLB{1}.HAd,HAd_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAr
            HAr_expected = 1.0e2*[
             0     0     0     0
             0     0     0     0
             0     0     0     0
             0     0     0     0];
            ABSTOL = 10*eps(max(max(abs(HAr_expected))));
            testCase.verifyEqual(AntLB{1}.HAr,HAr_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAP
            HAP_expected = 1.0e2*[
          -1.323743269248981  -1.580869618766140  -1.641359801346453  -1.676031839128706
          -1.324697223675199  -1.581883527002873  -1.641850621204388  -1.676353559431351
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAP_expected))));
            testCase.verifyEqual(AntLB{1}.HAP,HAP_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HRP
            HRP_expected = 1.0e2*[
          -1.323743269248981  -1.580869618766140  -1.641359801346453  -1.676031839128706
          -1.324697223675199  -1.581883527002873  -1.641850621204388  -1.676353559431351
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HRP_expected))));
            testCase.verifyEqual(AntLB{1}.HRP,HRP_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAt
            HAt_expected = 1.0e2*[
           0.099998647807879  -0.062463102346010  -0.066709326379805  -0.067103724194647
           0.099123802737152  -0.062610627080090  -0.066726147691143  -0.067103200765329
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAt_expected))));
            testCase.verifyEqual(AntLB{1}.HAt,HAt_expected,'AbsTol',ABSTOL);
        end
        
        function test_LinkBudget_wLightTime(testCase)
            % run gsmeas with a LinkBudget and range and rangerate
            % measurements
            measOptions = setOdtbxOptions(testCase.measOptions_LinkBudget,'useLightTime',true);
            % Run gsmeas
            [y,H,R,AntLB] = testCase.verifyWarning(@()gsmeas(testCase.t, testCase.x, measOptions),...
                'ODTBX:GSMEAS:noBodyQuat');
            % Check y
            y_expected = 1.0e3*[
               0.199032532782620   0.591898264309697   1.131011079192571                 NaN
               0.000001650017817   0.008741270828982   0.009100998874894                 NaN
               0.199213890231232   0.597831760838404   1.137200002451057                 NaN
               0.000293647613285   0.008751923697164   0.009102497539139                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);
            
            % Check H
            % only check first part of H to save space and time
            H_expected = [
           0.783398869395764  -0.229869625116898   0.577447977637861                   0                   0                   0
           0.013295980949963   0.045276967708025  -0.000014298864235   0.783398869395764  -0.229869625116898   0.577447977637861
           0.804422293995511  -0.203605957446998   0.558076506417436                   0                   0                   0
           0.012102719707128   0.045535549256131  -0.000832123267837   0.804422293995511  -0.203605957446998   0.558076506417436
           0.653141152169440   0.317407173316794   0.687502233938057                   0                   0                   0
          -0.000055634404980   0.000750713359407  -0.000293736770650   0.653141152169440   0.317407173316794   0.687502233938057
           0.666192635083809  -0.744729130164298  -0.039571399328553                   0                   0                   0
           0.000750527756830   0.000672611778259  -0.000023186456647   0.666192635083809  -0.744729130164298  -0.039571399328553];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.Halpha_r
            Halpha_r_expected = [
           0.003115704800299   1.168260941956579   1.275128535884713   1.274785638792920
           0.037196767999882   1.171438738118571   1.275774240924258   1.274913097277424
           0.582968622078663   0.545602755137455   0.508188440538280   0.471279381977418
           0.833955946910231   0.862165276398346   0.887662419802437   0.909884804974218];
            ABSTOL = 10*eps(max(max(abs(Halpha_r_expected))));
            testCase.verifyEqual(AntLB{1}.Halpha_r,Halpha_r_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HCN0
            HCN0_expected = 1.0e2*[
           0.684544595969717   0.427418237195692   0.366928069198436   0.332256035557153
           0.683590614486464   0.426404330079141   0.366437249445125   0.331934315301062
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HCN0_expected))));
            testCase.verifyEqual(AntLB{1}.HCN0,HCN0_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAd
            HAd_expected = 1.0e2*[
          -1.423741917056861  -1.518406516420129  -1.574650474966648  -1.608928114934059
          -1.423821026412351  -1.519272899922783  -1.575124473513245  -1.609250358666021
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAd_expected))));
            testCase.verifyEqual(AntLB{1}.HAd,HAd_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAr
            HAr_expected = 1.0e2*[
             0     0     0     0
             0     0     0     0
             0     0     0     0
             0     0     0     0];
            ABSTOL = 10*eps(max(max(abs(HAr_expected))));
            testCase.verifyEqual(AntLB{1}.HAr,HAr_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAP
            HAP_expected = 1.0e2*[
          -1.323743278558316  -1.580869637332342  -1.641359805329598  -1.676031838970881
          -1.324697260041570  -1.581883544448893  -1.641850625082909  -1.676353559226972
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAP_expected))));
            testCase.verifyEqual(AntLB{1}.HAP,HAP_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HRP
            HRP_expected = 1.0e2*[
          -1.323743278558316  -1.580869637332342  -1.641359805329598  -1.676031838970881
          -1.324697260041570  -1.581883544448893  -1.641850625082909  -1.676353559226972
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HRP_expected))));
            testCase.verifyEqual(AntLB{1}.HRP,HRP_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAt
            HAt_expected = 1.0e2*[
           0.099998638498544  -0.062463120912213  -0.066709330362950  -0.067103724036822
           0.099123766370781  -0.062610644526110  -0.066726151569663  -0.067103200560950
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAt_expected))));
            testCase.verifyEqual(AntLB{1}.HAt,HAt_expected,'AbsTol',ABSTOL);
        end
        
        function test_LinkBudget_2D_Antenna_with_Quat(testCase)
            % run gsmeas with a LinkBudget and range and rangerate
            % measurements
            qatt = repmat([0;0;0;1],1,length(testCase.t));
            % Run gsmeas
            [y,H,R,AntLB] = gsmeas(testCase.t, testCase.x, testCase.measOptions_2D_Antenna_LinkBudget,qatt);
            % Check y
            y_expected = 1.0e3*[
               0.199032532782455   0.591898264330166   1.131011079112153                 NaN
               0.000001650017813   0.008741270829023   0.009100998874876                 NaN
               0.199213890231063   0.597831760859305   1.137200002369821                 NaN
               0.000293647613281   0.008751923697204   0.009102497539121                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN
                             NaN                 NaN                 NaN                 NaN];
            ABSTOL = 10*eps(max(max(y_expected)));
            testCase.verifyEqual(y,y_expected,'AbsTol',ABSTOL);

            % Check H
            % only check first part of H to save space and time
            H_expected = [
                0.783398512203442        -0.229870841447615         0.577447978028843                         0                         0                         0
                0.0132960278492823        0.0452769539688087     -1.42643643439925e-05         0.783398512203442        -0.229870841447615         0.577447978028843
                0.804421968278982        -0.203607181003063         0.558076529513788                         0                         0                         0
                0.0121027687146992        0.0455355407620577     -0.000832090074627532         0.804421968278982        -0.203607181003063         0.558076529513788
                0.653140372345824         0.317407825572411         0.687502673651057                         0                         0                         0
               -5.56343720182588e-05       0.00075071336040374     -0.000293737389887119         0.653140372345824         0.317407825572411         0.687502673651057
                0.666193352210258        -0.744728490074234       -0.0395713727663203                         0                         0                         0
                0.00075052704805523       0.00067261244066317     -2.31863897485416e-05         0.666193352210258        -0.744728490074234       -0.0395713727663203];
            ABSTOL = 10*eps(max(max(H_expected)));
            testCase.verifyEqual(H(:,:,1),H_expected,'AbsTol',ABSTOL);
            % Check R
            R_expected = repmat(1e-6*eye(8),[1 1 4]);
            ABSTOL = 10*eps(max(max(max(abs(R_expected)))));
            testCase.verifyEqual(R,R_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.Halpha_r
            Halpha_r_expected = [
           2.186395708928785   1.749609010554037   1.638635019157889   1.588369514848663
           2.162862292752860   1.741439416839951   1.634956732029307   1.586075376313282
           2.328840774007467   2.304386536620801   2.279665098056304   2.254813938013354
           1.531214619318362   1.528636157135541   1.523545687653922   1.515756382082139];
            ABSTOL = 10*eps(max(max(abs(Halpha_r_expected))));
            testCase.verifyEqual(AntLB{1}.Halpha_r,Halpha_r_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HCN0
            HCN0_expected = 1.0e2*[
           0.545189938464646   0.422546574695219   0.346528304808443   0.299461781540165
           0.548593904464370   0.421843121665536   0.345875889088881   0.299074790716980
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HCN0_expected))));
            testCase.verifyEqual(AntLB{1}.HCN0,HCN0_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAd
            HAd_expected = 1.0e2*[
          -1.423741917056861  -1.518406516420129  -1.574650474966648  -1.608928114934059
          -1.423821026412351  -1.519272899922783  -1.575124473513245  -1.609250358666021
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAd_expected))));
            testCase.verifyEqual(AntLB{1}.HAd,HAd_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAr
            HAr_expected = 1.0e2*[
          -0.079355937953287  -0.058289587857354  -0.050873919192214  -0.049775266289088
          -0.075871974077950  -0.057790752337911  -0.050704102713239  -0.049788870177261
          -1.000000000000000  -1.000000000000000  -1.000000000000000  -1.000000000000000
          -1.000000000000000  -1.000000000000000  -1.000000000000000  -1.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAr_expected))));
            testCase.verifyEqual(AntLB{1}.HAr,HAr_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAP
            HAP_expected = 1.0e2*[
          -1.383741998110100  -1.527451711975461  -1.610885650527376  -1.659050826698781
          -1.383821995985714  -1.528654000524587  -1.611707882725913  -1.659424213633793
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAP_expected))));
            testCase.verifyEqual(AntLB{1}.HAP,HAP_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HRP
            HRP_expected = 1.0e2*[
          -1.463097936063387  -1.585741299832815  -1.661759569719591  -1.708826092987869
          -1.459693970063664  -1.586444752862497  -1.662411985439153  -1.709213083811054
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HRP_expected))));
            testCase.verifyEqual(AntLB{1}.HRP,HRP_expected,'AbsTol',ABSTOL);
            % Check AntLB{1}.HAt
            HAt_expected = 1.0e2*[
           0.039999918946760  -0.009045195555332  -0.036235175560728  -0.050122711764722
           0.039999030426637  -0.009381100601804  -0.036583409212668  -0.050173854967772
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000
          -3.000000000000000  -3.000000000000000  -3.000000000000000  -3.000000000000000];
            ABSTOL = 10*eps(max(max(abs(HAt_expected))));
            testCase.verifyEqual(AntLB{1}.HAt,HAt_expected,'AbsTol',ABSTOL);
        end
        
        function test_LinkBudget_2D_Antenna_with_BadQuat(testCase)
            % run gsmeas with a LinkBudget and range and rangerate
            % measurements
            qatt = repmat([0;0;1],1,length(testCase.t));
            % Run gsmeas
            testCase.verifyError(@()gsmeas(testCase.t, testCase.x, testCase.measOptions_2D_Antenna_LinkBudget,qatt),...
                'ODTBX:GSMEAS:badQuat');
            
        end
    end
    
end

