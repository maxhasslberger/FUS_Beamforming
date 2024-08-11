classdef simulationApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        SlicexLabel                     matlab.ui.control.Label
        UpdateSliceButton               matlab.ui.control.Button
        TabGroup                        matlab.ui.container.TabGroup
        InitTab                         matlab.ui.container.Tab
        UpdateButtonInit                matlab.ui.control.Button
        MediumPanel                     matlab.ui.container.Panel
        hetPanel                        matlab.ui.container.Panel
        homPanel                        matlab.ui.container.Panel
        HomzEditField                   matlab.ui.control.NumericEditField
        zEditField_2Label               matlab.ui.control.Label
        HomyEditField                   matlab.ui.control.NumericEditField
        yEditField_2Label               matlab.ui.control.Label
        xHomEditField                   matlab.ui.control.NumericEditField
        xEditField_2Label               matlab.ui.control.Label
        GridSizemmLabel                 matlab.ui.control.Label
        GreensFunctionbasedCheckBox     matlab.ui.control.CheckBox
        ScanSpatialResolutionmmEditField  matlab.ui.control.NumericEditField
        ScanSpatialResolutionmmEditFieldLabel  matlab.ui.control.Label
        CTfilenameDropDown              matlab.ui.control.DropDown
        CTfilenameDropDownLabel         matlab.ui.control.Label
        T1wfilenameDropDown             matlab.ui.control.DropDown
        T1wfilenameDropDownLabel        matlab.ui.control.Label
        scanyEditField                  matlab.ui.control.NumericEditField
        yEditFieldLabel                 matlab.ui.control.Label
        scanzEditField                  matlab.ui.control.NumericEditField
        zEditFieldLabel                 matlab.ui.control.Label
        scanxEditField                  matlab.ui.control.NumericEditField
        xEditFieldLabel                 matlab.ui.control.Label
        MinScanIndexineachDimensionLabel  matlab.ui.control.Label
        SwitchLabel                     matlab.ui.control.Label
        MediumSwitch                    matlab.ui.control.Switch
        DisplayAdvancedSettingsCheckBox  matlab.ui.control.CheckBox
        GeneralPanel                    matlab.ui.container.Panel
        SliceIndexEditField             matlab.ui.control.NumericEditField
        DSliceIndexEditFieldLabel       matlab.ui.control.Label
        SliceDirectionDropDown          matlab.ui.control.DropDown
        SliceDirectionDropDownLabel     matlab.ui.control.Label
        SimSpatialResolutionmmEditField  matlab.ui.control.NumericEditField
        SimSpatialResolutionmmEditFieldLabel  matlab.ui.control.Label
        CenterFreqkHzEditField          matlab.ui.control.NumericEditField
        CenterFreqkHzEditFieldLabel     matlab.ui.control.Label
        DimSwitch                       matlab.ui.control.Switch
        Label                           matlab.ui.control.Label
        SaveSimulationResultsCheckBox   matlab.ui.control.CheckBox
        RealdxLabel                     matlab.ui.control.Label
        AdvancedGridMediumSettingsPanel  matlab.ui.container.Panel
        PMLinsideCheckBox               matlab.ui.control.CheckBox
        PMLsizeEditField                matlab.ui.control.NumericEditField
        PMLsizeEditFieldLabel           matlab.ui.control.Label
        hounsfieldmaxEditField          matlab.ui.control.NumericEditField
        hounsfieldmaxEditFieldLabel     matlab.ui.control.Label
        hounsfieldminEditField          matlab.ui.control.NumericEditField
        hounsfieldminEditFieldLabel     matlab.ui.control.Label
        cflEditField                    matlab.ui.control.NumericEditField
        cflEditFieldLabel               matlab.ui.control.Label
        ppwEditField                    matlab.ui.control.NumericEditField
        ppwEditFieldLabel               matlab.ui.control.Label
        alphapowerEditField             matlab.ui.control.NumericEditField
        alphapowerEditFieldLabel        matlab.ui.control.Label
        alphamaxdBMHzycmEditField       matlab.ui.control.NumericEditField
        alphamaxdBMHzycmLabel           matlab.ui.control.Label
        alphamindBMHzycmEditField       matlab.ui.control.NumericEditField
        alphamindBMHzycmLabel           matlab.ui.control.Label
        alphawaterdBMHzycmEditField     matlab.ui.control.NumericEditField
        alphawaterdBMHzycmLabel         matlab.ui.control.Label
        rhomaxkgm3EditField             matlab.ui.control.NumericEditField
        rhomaxkgm3EditFieldLabel        matlab.ui.control.Label
        rho0kgm3EditField               matlab.ui.control.NumericEditField
        rho0kgm3EditFieldLabel          matlab.ui.control.Label
        cmaxmsEditField                 matlab.ui.control.NumericEditField
        cmaxmsEditFieldLabel            matlab.ui.control.Label
        c0msEditField                   matlab.ui.control.NumericEditField
        c0msEditFieldLabel              matlab.ui.control.Label
        TransducersTab                  matlab.ui.container.Tab
        GeneralPanel_2                  matlab.ui.container.Panel
        WidthmmEditField                matlab.ui.control.NumericEditField
        WidthmmEditFieldLabel           matlab.ui.control.Label
        LengthmmEditField               matlab.ui.control.NumericEditField
        LengthmmEditFieldLabel          matlab.ui.control.Label
        ElementGeometrySwitch           matlab.ui.control.Switch
        ElementGeometrySwitchLabel      matlab.ui.control.Label
        SparsityfilenameDropDown        matlab.ui.control.DropDown
        SparsityfilenameDropDownLabel   matlab.ui.control.Label
        ArrayElementsPositionsfilenameDropDown  matlab.ui.control.DropDown
        ArrayElementsPositionsfilenameDropDownLabel  matlab.ui.control.Label
        PropagationMatrixAfilenameDropDown  matlab.ui.control.DropDown
        PropagationMatrixAfilenameDropDownLabel  matlab.ui.control.Label
        UpdateButtonTransducer          matlab.ui.control.Button
        RemoveTransducerButton          matlab.ui.control.Button
        AddTransducerButton             matlab.ui.control.Button
        Transducer1Panel                matlab.ui.container.Panel
        gammaEditField                  matlab.ui.control.NumericEditField
        gammaEditFieldLabel             matlab.ui.control.Label
        betaEditField                   matlab.ui.control.NumericEditField
        betaEditFieldLabel              matlab.ui.control.Label
        alphaEditField                  matlab.ui.control.NumericEditField
        alphaEditFieldLabel             matlab.ui.control.Label
        RotationdegxyzintrinsicLabel    matlab.ui.control.Label
        trPoszEditField                 matlab.ui.control.NumericEditField
        zEditField_3Label               matlab.ui.control.Label
        trPosyEditField                 matlab.ui.control.NumericEditField
        yEditField_3Label               matlab.ui.control.Label
        trPosxEditField                 matlab.ui.control.NumericEditField
        xEditField_3Label               matlab.ui.control.Label
        FacePositionCenterLabel         matlab.ui.control.Label
        TransducerDropDown              matlab.ui.control.DropDown
        TransducerDropDownLabel         matlab.ui.control.Label
        TargetingTab                    matlab.ui.container.Tab
        MaxPressurekPaEditField_2       matlab.ui.control.NumericEditField
        MaxPressurekPaEditField_2Label  matlab.ui.control.Label
        LimitSkullPressureCheckBox      matlab.ui.control.CheckBox
        LimitIntracranialOffTargetPressureCheckBox  matlab.ui.control.CheckBox
        MaxPressurekPaEditField         matlab.ui.control.NumericEditField
        MaxPressurekPaEditFieldLabel    matlab.ui.control.Label
        UpdateButtonTargeting           matlab.ui.control.Button
        RegionTargetDropDown            matlab.ui.control.DropDown
        RegionTargetDropDownLabel       matlab.ui.control.Label
        RemoveRegionTargetButton        matlab.ui.control.Button
        AddRegionTargetButton           matlab.ui.control.Button
        Transducer1Panel_3              matlab.ui.container.Panel
        PressureAmplitudekPaEditField_2  matlab.ui.control.NumericEditField
        PressureAmplitudekPaEditField_2Label  matlab.ui.control.Label
        todeclaredPressureSwitch        matlab.ui.control.Switch
        todeclaredPressureSwitchLabel   matlab.ui.control.Label
        MinPointDistancemmEditField_2   matlab.ui.control.NumericEditField
        MinPointDistancemmEditField_2Label  matlab.ui.control.Label
        BrainRegionDropDown             matlab.ui.control.DropDown
        BrainRegionDropDownLabel        matlab.ui.control.Label
        ManualTargetDropDown            matlab.ui.control.DropDown
        ManualTargetDropDownLabel       matlab.ui.control.Label
        RemoveManualTransducerButton    matlab.ui.control.Button
        AddManualTargetButton           matlab.ui.control.Button
        Transducer1Panel_2              matlab.ui.container.Panel
        todeclaredPressureSwitch_2      matlab.ui.control.Switch
        todeclaredPressureSwitch_2Label  matlab.ui.control.Label
        MinPointDistancemmEditField     matlab.ui.control.NumericEditField
        MinPointDistancemmEditFieldLabel  matlab.ui.control.Label
        PressureAmplitudekPaEditField   matlab.ui.control.NumericEditField
        PressureAmplitudekPaEditFieldLabel  matlab.ui.control.Label
        FocusRadiusmmEditField          matlab.ui.control.NumericEditField
        FocusRadiusmmEditFieldLabel     matlab.ui.control.Label
        focuszEditField                 matlab.ui.control.NumericEditField
        zEditField_4Label               matlab.ui.control.Label
        focusyEditField                 matlab.ui.control.NumericEditField
        yEditField_4Label               matlab.ui.control.Label
        focusxEditField                 matlab.ui.control.NumericEditField
        xEditField_4Label               matlab.ui.control.Label
        FocusPositionLabel              matlab.ui.control.Label
        OptimizeTab                     matlab.ui.container.Tab
        GroundTruthResolutionFactorEditField  matlab.ui.control.NumericEditField
        GroundTruthResolutionFactorEditFieldLabel  matlab.ui.control.Label
        PlotSkullCheckBox               matlab.ui.control.CheckBox
        PlotEntireDomainCheckBox        matlab.ui.control.CheckBox
        OptimizationModeButtonGroup     matlab.ui.container.ButtonGroup
        OptimizeforPhasesand1AmpperTransducerButton  matlab.ui.control.RadioButton
        OptimizeforTransducerPhasesandAmplitudesButton  matlab.ui.control.RadioButton
        OptimizeButton                  matlab.ui.control.Button
        ComputeInitialSolutionButton    matlab.ui.control.Button
        DoGroundTruthSimulationButton   matlab.ui.control.Button
    end

    
    properties (Access = public)
        n_dim
        kgrid
        medium
        sensor
        sensor_mask
        dx_factor
        grid_size
        t1w_filename
        dx_scan
        plot_offset
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: UpdateButtonInit
        function UpdateButtonInitPushed(app, event)
            f0 = app.CenterFreqkHzEditField.Value * 1e3; % Hz
            app.n_dim = 2 + (app.DimSwitch.Value == "3D");
            app.plot_offset = [-app.scanxEditField.Value, -app.scanyEditField.Value, -app.scanzEditField.Value] + 1;

            %% Obtain Constants
            const.c0 = app.c0msEditField.Value; % m/s - water
            const.c_max = app.cmaxmsEditField.Value; % m/s - skull
            
            const.rho0 = app.rho0kgm3EditField.Value; % kg/m^3
            const.rho_max = app.rhomaxkgm3EditField.Value; % kg/m^3
            
            const.alpha_coeff_water = app.alphawaterdBMHzycmEditField.Value; % dB/(MHz^y cm)
            const.alpha_coeff_min = app.alphamindBMHzycmEditField.Value; % dB/(MHz^y cm)
            const.alpha_coeff_max = app.alphamaxdBMHzycmEditField.Value; % dB/(MHz^y cm)
            
            const.alpha_power = app.alphapowerEditField.Value;
            
            const.hu_min = app.hounsfieldminEditField.Value; % Hounsfield units
            const.hu_max = app.hounsfieldmaxEditField.Value;
            
            const.ppw = app.ppwEditField.Value; % >= 2
            const.cfl = app.cflEditField.Value;

            %% Define Grid and Medium
            SimSpatialResolutionmmEditFieldValueChanged(app);

            if strcmp(app.MediumSwitch, 'Homogeneous')
                app.grid_size = [app.xHomEditField, app.HomyEditField, app.HomzEditField] * 1e-3;
                app.t1w_filename = [];
                ct_filename = [];
                app.dx_scan = 1e-3; % m
            else
                app.grid_size = [];
                app.t1w_filename = fullfile('..', 'Scans', strcat(app.T1wfilenameDropDown.Value, '.nii'));
                ct_filename = fullfile('..', 'Scans', strcat(app.CTfilenameDropDown.Value, '.nii'));
                app.dx_scan = app.ScanSpatialResolutionmmEditField.Value;
            end

            [app.kgrid, app.medium, ppp] = init_grid_medium(f0, app.grid_size, 'n_dim', app.n_dim, 'dx_factor', app.dx_factor, ...
                'ct_scan', ct_filename, 'slice_idx', round(app.plot_offset(2) + app.SliceIndexEditField.Value), ...
                'dx_scan', app.dx_scan, 'constants', const);
            [app.sensor, app.sensor_mask] = init_sensor(app.kgrid, ppp);

            app.dx_factor = app.dx_scan / app.kgrid.dx;
        end

        % Value changed function: SimSpatialResolutionmmEditField
        function SimSpatialResolutionmmEditFieldValueChanged(app, event)
            value = app.SimSpatialResolutionmmEditField.Value * 1e-3;
            
            dx_std = 1500 / (app.CenterFreqkHzEditField.Value * 1e3) / 3; % = (c0 / f0 / ppw)
            if value <= 0
                app.dx_factor = 1;
                dx = dx_std;
            else
                app.dx_factor = dx_std / value;
                dx = value;
            end

            app.RealdxLabel.Text = strcat("-> Real: ", string(dx * 1e3), " mm");
        end

        % Value changed function: DisplayAdvancedSettingsCheckBox
        function DisplayAdvancedSettingsCheckBoxValueChanged(app, event)
            value = app.DisplayAdvancedSettingsCheckBox.Value;
            
            if value
                app.AdvancedGridMediumSettingsPanel.Visible = true;
            else
                app.AdvancedGridMediumSettingsPanel.Visible = false;
            end
        end

        % Value changed function: MediumSwitch
        function MediumSwitchValueChanged(app, event)
            value = app.MediumSwitch.Value;
            
            if strcmp(value, "Homogeneous")
                app.homPanel.Visible = true;
            else
                app.homPanel.Visible = false;

                % Create t1w and ct dropdown menu
                items = dir(fullfile('..', 'Scans'));
                t1_pattern = {'T1', 't1'};
                ct_pattern = {'ct', 'CT'};
                scan_pattern = '.nii';

                % Check if each item matches any of the patterns
                t1w_matches = false(size(items))';
                ct_matches = false(size(items))';
                for i = 1:length(t1_pattern)
                    t1w_matches = t1w_matches | (contains({items.name}, t1_pattern{i}) & contains({items.name}, scan_pattern));
                    ct_matches = ct_matches | (contains({items.name}, ct_pattern{i}) & contains({items.name}, scan_pattern));
                end

                % Assign to t1 and ct dropdown menus
                app.T1wfilenameDropDown.Items = {items(t1w_matches).name};
                app.CTfilenameDropDown.Items = {items(ct_matches).name};
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1064 550];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 727 550];

            % Create InitTab
            app.InitTab = uitab(app.TabGroup);
            app.InitTab.Title = 'Init';

            % Create AdvancedGridMediumSettingsPanel
            app.AdvancedGridMediumSettingsPanel = uipanel(app.InitTab);
            app.AdvancedGridMediumSettingsPanel.Title = 'Advanced Grid/Medium Settings';
            app.AdvancedGridMediumSettingsPanel.Visible = 'off';
            app.AdvancedGridMediumSettingsPanel.Position = [33 13 678 158];

            % Create c0msEditFieldLabel
            app.c0msEditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.c0msEditFieldLabel.HorizontalAlignment = 'right';
            app.c0msEditFieldLabel.Position = [54 104 48 22];
            app.c0msEditFieldLabel.Text = 'c0 (m/s)';

            % Create c0msEditField
            app.c0msEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.c0msEditField.Position = [113 104 45 22];
            app.c0msEditField.Value = 1500;

            % Create cmaxmsEditFieldLabel
            app.cmaxmsEditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.cmaxmsEditFieldLabel.HorizontalAlignment = 'right';
            app.cmaxmsEditFieldLabel.Position = [35 70 68 22];
            app.cmaxmsEditFieldLabel.Text = 'c max (m/s)';

            % Create cmaxmsEditField
            app.cmaxmsEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.cmaxmsEditField.Position = [114 70 45 22];
            app.cmaxmsEditField.Value = 3100;

            % Create rho0kgm3EditFieldLabel
            app.rho0kgm3EditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.rho0kgm3EditFieldLabel.HorizontalAlignment = 'right';
            app.rho0kgm3EditFieldLabel.Position = [23 42 79 22];
            app.rho0kgm3EditFieldLabel.Text = 'rho0 (kg/m^3)';

            % Create rho0kgm3EditField
            app.rho0kgm3EditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.rho0kgm3EditField.Position = [112 42 48 22];
            app.rho0kgm3EditField.Value = 1000;

            % Create rhomaxkgm3EditFieldLabel
            app.rhomaxkgm3EditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.rhomaxkgm3EditFieldLabel.HorizontalAlignment = 'right';
            app.rhomaxkgm3EditFieldLabel.Position = [5 11 98 22];
            app.rhomaxkgm3EditFieldLabel.Text = 'rho max (kg/m^3)';

            % Create rhomaxkgm3EditField
            app.rhomaxkgm3EditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.rhomaxkgm3EditField.Position = [113 11 48 22];
            app.rhomaxkgm3EditField.Value = 1900;

            % Create alphawaterdBMHzycmLabel
            app.alphawaterdBMHzycmLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.alphawaterdBMHzycmLabel.HorizontalAlignment = 'right';
            app.alphawaterdBMHzycmLabel.Position = [173 104 152 22];
            app.alphawaterdBMHzycmLabel.Text = 'alpha water (dB/MHz^y cm)';

            % Create alphawaterdBMHzycmEditField
            app.alphawaterdBMHzycmEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.alphawaterdBMHzycmEditField.Position = [337 104 39 22];
            app.alphawaterdBMHzycmEditField.Value = 0.75;

            % Create alphamindBMHzycmLabel
            app.alphamindBMHzycmLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.alphamindBMHzycmLabel.HorizontalAlignment = 'right';
            app.alphamindBMHzycmLabel.Position = [183 71 142 22];
            app.alphamindBMHzycmLabel.Text = 'alpha min (dB/MHz^y cm)';

            % Create alphamindBMHzycmEditField
            app.alphamindBMHzycmEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.alphamindBMHzycmEditField.Position = [337 71 39 22];
            app.alphamindBMHzycmEditField.Value = 4;

            % Create alphamaxdBMHzycmLabel
            app.alphamaxdBMHzycmLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.alphamaxdBMHzycmLabel.HorizontalAlignment = 'right';
            app.alphamaxdBMHzycmLabel.Position = [183 42 145 22];
            app.alphamaxdBMHzycmLabel.Text = 'alpha max (dB/MHz^y cm)';

            % Create alphamaxdBMHzycmEditField
            app.alphamaxdBMHzycmEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.alphamaxdBMHzycmEditField.Position = [340 42 39 22];
            app.alphamaxdBMHzycmEditField.Value = 8.7;

            % Create alphapowerEditFieldLabel
            app.alphapowerEditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.alphapowerEditFieldLabel.HorizontalAlignment = 'right';
            app.alphapowerEditFieldLabel.Position = [257 11 70 22];
            app.alphapowerEditFieldLabel.Text = 'alpha power';

            % Create alphapowerEditField
            app.alphapowerEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.alphapowerEditField.Position = [337 11 42 22];
            app.alphapowerEditField.Value = 1.43;

            % Create ppwEditFieldLabel
            app.ppwEditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.ppwEditFieldLabel.HorizontalAlignment = 'right';
            app.ppwEditFieldLabel.Position = [455 11 27 22];
            app.ppwEditFieldLabel.Text = 'ppw';

            % Create ppwEditField
            app.ppwEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.ppwEditField.Position = [490 11 38 22];
            app.ppwEditField.Value = 3;

            % Create cflEditFieldLabel
            app.cflEditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.cflEditFieldLabel.HorizontalAlignment = 'right';
            app.cflEditFieldLabel.Position = [454 42 25 22];
            app.cflEditFieldLabel.Text = 'cfl';

            % Create cflEditField
            app.cflEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.cflEditField.Position = [487 42 40 22];
            app.cflEditField.Value = 0.3;

            % Create hounsfieldminEditFieldLabel
            app.hounsfieldminEditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.hounsfieldminEditFieldLabel.HorizontalAlignment = 'right';
            app.hounsfieldminEditFieldLabel.Position = [390 105 82 22];
            app.hounsfieldminEditFieldLabel.Text = 'hounsfield min';

            % Create hounsfieldminEditField
            app.hounsfieldminEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.hounsfieldminEditField.Position = [482 105 48 22];
            app.hounsfieldminEditField.Value = 300;

            % Create hounsfieldmaxEditFieldLabel
            app.hounsfieldmaxEditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.hounsfieldmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.hounsfieldmaxEditFieldLabel.Position = [386 71 86 22];
            app.hounsfieldmaxEditFieldLabel.Text = 'hounsfield max';

            % Create hounsfieldmaxEditField
            app.hounsfieldmaxEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.hounsfieldmaxEditField.Position = [482 71 48 22];
            app.hounsfieldmaxEditField.Value = 2000;

            % Create PMLsizeEditFieldLabel
            app.PMLsizeEditFieldLabel = uilabel(app.AdvancedGridMediumSettingsPanel);
            app.PMLsizeEditFieldLabel.HorizontalAlignment = 'right';
            app.PMLsizeEditFieldLabel.Position = [549 104 54 22];
            app.PMLsizeEditFieldLabel.Text = 'PML size';

            % Create PMLsizeEditField
            app.PMLsizeEditField = uieditfield(app.AdvancedGridMediumSettingsPanel, 'numeric');
            app.PMLsizeEditField.Position = [608 104 35 22];
            app.PMLsizeEditField.Value = 10;

            % Create PMLinsideCheckBox
            app.PMLinsideCheckBox = uicheckbox(app.AdvancedGridMediumSettingsPanel);
            app.PMLinsideCheckBox.Text = 'PML inside';
            app.PMLinsideCheckBox.Position = [562 71 81 22];
            app.PMLinsideCheckBox.Value = true;

            % Create GeneralPanel
            app.GeneralPanel = uipanel(app.InitTab);
            app.GeneralPanel.Title = 'General';
            app.GeneralPanel.Position = [31 220 298 292];

            % Create RealdxLabel
            app.RealdxLabel = uilabel(app.GeneralPanel);
            app.RealdxLabel.Position = [47 58 81 22];
            app.RealdxLabel.Text = '-> Real: 1 mm';

            % Create SaveSimulationResultsCheckBox
            app.SaveSimulationResultsCheckBox = uicheckbox(app.GeneralPanel);
            app.SaveSimulationResultsCheckBox.Text = 'Save Simulation Results';
            app.SaveSimulationResultsCheckBox.Position = [74 9 152 30];

            % Create Label
            app.Label = uilabel(app.GeneralPanel);
            app.Label.HorizontalAlignment = 'center';
            app.Label.Position = [40 240 26 22];
            app.Label.Text = 'Dim';

            % Create DimSwitch
            app.DimSwitch = uiswitch(app.GeneralPanel, 'slider');
            app.DimSwitch.Items = {'2D', '3D'};
            app.DimSwitch.Position = [30 210 45 20];
            app.DimSwitch.Value = '2D';

            % Create CenterFreqkHzEditFieldLabel
            app.CenterFreqkHzEditFieldLabel = uilabel(app.GeneralPanel);
            app.CenterFreqkHzEditFieldLabel.HorizontalAlignment = 'right';
            app.CenterFreqkHzEditFieldLabel.Position = [126 210 101 22];
            app.CenterFreqkHzEditFieldLabel.Text = 'Center Freq (kHz)';

            % Create CenterFreqkHzEditField
            app.CenterFreqkHzEditField = uieditfield(app.GeneralPanel, 'numeric');
            app.CenterFreqkHzEditField.Position = [239 210 46 22];
            app.CenterFreqkHzEditField.Value = 500;

            % Create SimSpatialResolutionmmEditFieldLabel
            app.SimSpatialResolutionmmEditFieldLabel = uilabel(app.GeneralPanel);
            app.SimSpatialResolutionmmEditFieldLabel.HorizontalAlignment = 'right';
            app.SimSpatialResolutionmmEditFieldLabel.Position = [41 79 157 22];
            app.SimSpatialResolutionmmEditFieldLabel.Text = 'Sim Spatial Resolution (mm)';

            % Create SimSpatialResolutionmmEditField
            app.SimSpatialResolutionmmEditField = uieditfield(app.GeneralPanel, 'numeric');
            app.SimSpatialResolutionmmEditField.ValueChangedFcn = createCallbackFcn(app, @SimSpatialResolutionmmEditFieldValueChanged, true);
            app.SimSpatialResolutionmmEditField.Position = [213 79 45 22];

            % Create SliceDirectionDropDownLabel
            app.SliceDirectionDropDownLabel = uilabel(app.GeneralPanel);
            app.SliceDirectionDropDownLabel.HorizontalAlignment = 'right';
            app.SliceDirectionDropDownLabel.Position = [41 160 90 22];
            app.SliceDirectionDropDownLabel.Text = 'Slice (Direction)';

            % Create SliceDirectionDropDown
            app.SliceDirectionDropDown = uidropdown(app.GeneralPanel);
            app.SliceDirectionDropDown.Items = {'X', 'Y', 'Z'};
            app.SliceDirectionDropDown.Position = [146 160 45 22];
            app.SliceDirectionDropDown.Value = 'Y';

            % Create DSliceIndexEditFieldLabel
            app.DSliceIndexEditFieldLabel = uilabel(app.GeneralPanel);
            app.DSliceIndexEditFieldLabel.HorizontalAlignment = 'right';
            app.DSliceIndexEditFieldLabel.Position = [46 126 82 22];
            app.DSliceIndexEditFieldLabel.Text = '2D Slice Index';

            % Create SliceIndexEditField
            app.SliceIndexEditField = uieditfield(app.GeneralPanel, 'numeric');
            app.SliceIndexEditField.Position = [144 126 48 22];
            app.SliceIndexEditField.Value = 30;

            % Create DisplayAdvancedSettingsCheckBox
            app.DisplayAdvancedSettingsCheckBox = uicheckbox(app.InitTab);
            app.DisplayAdvancedSettingsCheckBox.ValueChangedFcn = createCallbackFcn(app, @DisplayAdvancedSettingsCheckBoxValueChanged, true);
            app.DisplayAdvancedSettingsCheckBox.Text = 'Display Advanced Settings';
            app.DisplayAdvancedSettingsCheckBox.Position = [32 182 164 22];

            % Create MediumPanel
            app.MediumPanel = uipanel(app.InitTab);
            app.MediumPanel.Title = 'Medium';
            app.MediumPanel.Position = [379 220 323 292];

            % Create MediumSwitch
            app.MediumSwitch = uiswitch(app.MediumPanel, 'slider');
            app.MediumSwitch.Items = {'Homogeneous', 'Heterogeneous'};
            app.MediumSwitch.ValueChangedFcn = createCallbackFcn(app, @MediumSwitchValueChanged, true);
            app.MediumSwitch.Position = [139 242 45 20];
            app.MediumSwitch.Value = 'Homogeneous';

            % Create SwitchLabel
            app.SwitchLabel = uilabel(app.MediumPanel);
            app.SwitchLabel.HorizontalAlignment = 'center';
            app.SwitchLabel.Position = [148 216 25 22];
            app.SwitchLabel.Text = ' ';

            % Create hetPanel
            app.hetPanel = uipanel(app.MediumPanel);
            app.hetPanel.Position = [1 0 322 233];

            % Create MinScanIndexineachDimensionLabel
            app.MinScanIndexineachDimensionLabel = uilabel(app.hetPanel);
            app.MinScanIndexineachDimensionLabel.HorizontalAlignment = 'center';
            app.MinScanIndexineachDimensionLabel.Position = [111 104 104 30];
            app.MinScanIndexineachDimensionLabel.Text = {'Min Scan Index in '; 'each Dimension'};

            % Create xEditFieldLabel
            app.xEditFieldLabel = uilabel(app.hetPanel);
            app.xEditFieldLabel.HorizontalAlignment = 'right';
            app.xEditFieldLabel.Position = [22 71 25 22];
            app.xEditFieldLabel.Text = 'x';

            % Create scanxEditField
            app.scanxEditField = uieditfield(app.hetPanel, 'numeric');
            app.scanxEditField.Position = [62 71 42 22];
            app.scanxEditField.Value = -96;

            % Create zEditFieldLabel
            app.zEditFieldLabel = uilabel(app.hetPanel);
            app.zEditFieldLabel.HorizontalAlignment = 'right';
            app.zEditFieldLabel.Position = [225 71 25 22];
            app.zEditFieldLabel.Text = 'z';

            % Create scanzEditField
            app.scanzEditField = uieditfield(app.hetPanel, 'numeric');
            app.scanzEditField.Position = [265 71 42 22];
            app.scanzEditField.Value = -126;

            % Create yEditFieldLabel
            app.yEditFieldLabel = uilabel(app.hetPanel);
            app.yEditFieldLabel.HorizontalAlignment = 'right';
            app.yEditFieldLabel.Position = [119 71 25 22];
            app.yEditFieldLabel.Text = 'y';

            % Create scanyEditField
            app.scanyEditField = uieditfield(app.hetPanel, 'numeric');
            app.scanyEditField.Position = [159 71 42 22];
            app.scanyEditField.Value = -127;

            % Create T1wfilenameDropDownLabel
            app.T1wfilenameDropDownLabel = uilabel(app.hetPanel);
            app.T1wfilenameDropDownLabel.HorizontalAlignment = 'right';
            app.T1wfilenameDropDownLabel.Position = [27 195 76 22];
            app.T1wfilenameDropDownLabel.Text = 'T1w filename';

            % Create T1wfilenameDropDown
            app.T1wfilenameDropDown = uidropdown(app.hetPanel);
            app.T1wfilenameDropDown.Items = {};
            app.T1wfilenameDropDown.Position = [118 195 189 22];
            app.T1wfilenameDropDown.Value = {};

            % Create CTfilenameDropDownLabel
            app.CTfilenameDropDownLabel = uilabel(app.hetPanel);
            app.CTfilenameDropDownLabel.HorizontalAlignment = 'right';
            app.CTfilenameDropDownLabel.Position = [34 160 69 22];
            app.CTfilenameDropDownLabel.Text = 'CT filename';

            % Create CTfilenameDropDown
            app.CTfilenameDropDown = uidropdown(app.hetPanel);
            app.CTfilenameDropDown.Items = {};
            app.CTfilenameDropDown.Position = [118 160 189 22];
            app.CTfilenameDropDown.Value = {};

            % Create ScanSpatialResolutionmmEditFieldLabel
            app.ScanSpatialResolutionmmEditFieldLabel = uilabel(app.hetPanel);
            app.ScanSpatialResolutionmmEditFieldLabel.HorizontalAlignment = 'right';
            app.ScanSpatialResolutionmmEditFieldLabel.Position = [47 18 164 22];
            app.ScanSpatialResolutionmmEditFieldLabel.Text = 'Scan Spatial Resolution (mm)';

            % Create ScanSpatialResolutionmmEditField
            app.ScanSpatialResolutionmmEditField = uieditfield(app.hetPanel, 'numeric');
            app.ScanSpatialResolutionmmEditField.Position = [215 18 30 22];
            app.ScanSpatialResolutionmmEditField.Value = 1;

            % Create homPanel
            app.homPanel = uipanel(app.hetPanel);
            app.homPanel.Position = [-1 14 308 219];

            % Create GreensFunctionbasedCheckBox
            app.GreensFunctionbasedCheckBox = uicheckbox(app.homPanel);
            app.GreensFunctionbasedCheckBox.Text = 'Green''s Function based';
            app.GreensFunctionbasedCheckBox.Position = [89 160 149 22];
            app.GreensFunctionbasedCheckBox.Value = true;

            % Create GridSizemmLabel
            app.GridSizemmLabel = uilabel(app.homPanel);
            app.GridSizemmLabel.Position = [126 92 86 22];
            app.GridSizemmLabel.Text = 'Grid Size (mm)';

            % Create xEditField_2Label
            app.xEditField_2Label = uilabel(app.homPanel);
            app.xEditField_2Label.HorizontalAlignment = 'right';
            app.xEditField_2Label.Position = [34 62 25 22];
            app.xEditField_2Label.Text = 'x';

            % Create xHomEditField
            app.xHomEditField = uieditfield(app.homPanel, 'numeric');
            app.xHomEditField.Limits = [0 Inf];
            app.xHomEditField.Position = [69 62 37 22];
            app.xHomEditField.Value = 192;

            % Create yEditField_2Label
            app.yEditField_2Label = uilabel(app.homPanel);
            app.yEditField_2Label.HorizontalAlignment = 'right';
            app.yEditField_2Label.Position = [114 62 25 22];
            app.yEditField_2Label.Text = 'y';

            % Create HomyEditField
            app.HomyEditField = uieditfield(app.homPanel, 'numeric');
            app.HomyEditField.Limits = [0 Inf];
            app.HomyEditField.Position = [149 62 37 22];
            app.HomyEditField.Value = 256;

            % Create zEditField_2Label
            app.zEditField_2Label = uilabel(app.homPanel);
            app.zEditField_2Label.HorizontalAlignment = 'right';
            app.zEditField_2Label.Position = [198 62 25 22];
            app.zEditField_2Label.Text = 'z';

            % Create HomzEditField
            app.HomzEditField = uieditfield(app.homPanel, 'numeric');
            app.HomzEditField.Limits = [0 Inf];
            app.HomzEditField.Position = [233 62 37 22];
            app.HomzEditField.Value = 256;

            % Create UpdateButtonInit
            app.UpdateButtonInit = uibutton(app.InitTab, 'push');
            app.UpdateButtonInit.ButtonPushedFcn = createCallbackFcn(app, @UpdateButtonInitPushed, true);
            app.UpdateButtonInit.Position = [602 182 100 23];
            app.UpdateButtonInit.Text = 'Update';

            % Create TransducersTab
            app.TransducersTab = uitab(app.TabGroup);
            app.TransducersTab.Title = 'Transducer(s)';

            % Create TransducerDropDownLabel
            app.TransducerDropDownLabel = uilabel(app.TransducersTab);
            app.TransducerDropDownLabel.HorizontalAlignment = 'right';
            app.TransducerDropDownLabel.Position = [48 311 75 22];
            app.TransducerDropDownLabel.Text = 'Transducer #';

            % Create TransducerDropDown
            app.TransducerDropDown = uidropdown(app.TransducersTab);
            app.TransducerDropDown.Items = {'1'};
            app.TransducerDropDown.Position = [138 311 100 22];
            app.TransducerDropDown.Value = '1';

            % Create Transducer1Panel
            app.Transducer1Panel = uipanel(app.TransducersTab);
            app.Transducer1Panel.Title = 'Transducer 1';
            app.Transducer1Panel.Position = [294 159 383 229];

            % Create FacePositionCenterLabel
            app.FacePositionCenterLabel = uilabel(app.Transducer1Panel);
            app.FacePositionCenterLabel.Position = [127 179 125 22];
            app.FacePositionCenterLabel.Text = 'Face Position (Center)';

            % Create xEditField_3Label
            app.xEditField_3Label = uilabel(app.Transducer1Panel);
            app.xEditField_3Label.HorizontalAlignment = 'right';
            app.xEditField_3Label.Position = [55 149 25 22];
            app.xEditField_3Label.Text = 'x';

            % Create trPosxEditField
            app.trPosxEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.trPosxEditField.Position = [90 149 37 22];
            app.trPosxEditField.Value = -59;

            % Create yEditField_3Label
            app.yEditField_3Label = uilabel(app.Transducer1Panel);
            app.yEditField_3Label.HorizontalAlignment = 'right';
            app.yEditField_3Label.Position = [135 149 25 22];
            app.yEditField_3Label.Text = 'y';

            % Create trPosyEditField
            app.trPosyEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.trPosyEditField.Position = [170 149 37 22];
            app.trPosyEditField.Value = 30;

            % Create zEditField_3Label
            app.zEditField_3Label = uilabel(app.Transducer1Panel);
            app.zEditField_3Label.HorizontalAlignment = 'right';
            app.zEditField_3Label.Position = [219 149 25 22];
            app.zEditField_3Label.Text = 'z';

            % Create trPoszEditField
            app.trPoszEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.trPoszEditField.Position = [254 149 37 22];
            app.trPoszEditField.Value = 68;

            % Create RotationdegxyzintrinsicLabel
            app.RotationdegxyzintrinsicLabel = uilabel(app.Transducer1Panel);
            app.RotationdegxyzintrinsicLabel.Position = [100 77 183 22];
            app.RotationdegxyzintrinsicLabel.Text = 'Rotation (deg) - x-y''-z'''' -> intrinsic';

            % Create alphaEditFieldLabel
            app.alphaEditFieldLabel = uilabel(app.Transducer1Panel);
            app.alphaEditFieldLabel.HorizontalAlignment = 'right';
            app.alphaEditFieldLabel.Position = [42 47 34 22];
            app.alphaEditFieldLabel.Text = 'alpha';

            % Create alphaEditField
            app.alphaEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.alphaEditField.Position = [86 47 37 22];

            % Create betaEditFieldLabel
            app.betaEditFieldLabel = uilabel(app.Transducer1Panel);
            app.betaEditFieldLabel.HorizontalAlignment = 'right';
            app.betaEditFieldLabel.Position = [137 47 28 22];
            app.betaEditFieldLabel.Text = 'beta';

            % Create betaEditField
            app.betaEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.betaEditField.Position = [175 47 37 22];
            app.betaEditField.Value = 45;

            % Create gammaEditFieldLabel
            app.gammaEditFieldLabel = uilabel(app.Transducer1Panel);
            app.gammaEditFieldLabel.HorizontalAlignment = 'right';
            app.gammaEditFieldLabel.Position = [233 47 45 22];
            app.gammaEditFieldLabel.Text = 'gamma';

            % Create gammaEditField
            app.gammaEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.gammaEditField.Position = [288 47 37 22];
            app.gammaEditField.Value = 180;

            % Create AddTransducerButton
            app.AddTransducerButton = uibutton(app.TransducersTab, 'push');
            app.AddTransducerButton.Position = [129 251 122 23];
            app.AddTransducerButton.Text = 'Add Transducer';

            % Create RemoveTransducerButton
            app.RemoveTransducerButton = uibutton(app.TransducersTab, 'push');
            app.RemoveTransducerButton.Position = [128 221 123 23];
            app.RemoveTransducerButton.Text = 'Remove Transducer';

            % Create UpdateButtonTransducer
            app.UpdateButtonTransducer = uibutton(app.TransducersTab, 'push');
            app.UpdateButtonTransducer.Position = [574 116 100 23];
            app.UpdateButtonTransducer.Text = 'Update';

            % Create PropagationMatrixAfilenameDropDownLabel
            app.PropagationMatrixAfilenameDropDownLabel = uilabel(app.TransducersTab);
            app.PropagationMatrixAfilenameDropDownLabel.HorizontalAlignment = 'right';
            app.PropagationMatrixAfilenameDropDownLabel.Position = [39 119 174 22];
            app.PropagationMatrixAfilenameDropDownLabel.Text = 'Propagation Matrix (A) filename';

            % Create PropagationMatrixAfilenameDropDown
            app.PropagationMatrixAfilenameDropDown = uidropdown(app.TransducersTab);
            app.PropagationMatrixAfilenameDropDown.Items = {''};
            app.PropagationMatrixAfilenameDropDown.Position = [228 119 100 22];
            app.PropagationMatrixAfilenameDropDown.Value = '';

            % Create GeneralPanel_2
            app.GeneralPanel_2 = uipanel(app.TransducersTab);
            app.GeneralPanel_2.Title = 'General';
            app.GeneralPanel_2.Position = [28 414 649 104];

            % Create ArrayElementsPositionsfilenameDropDownLabel
            app.ArrayElementsPositionsfilenameDropDownLabel = uilabel(app.GeneralPanel_2);
            app.ArrayElementsPositionsfilenameDropDownLabel.HorizontalAlignment = 'right';
            app.ArrayElementsPositionsfilenameDropDownLabel.Position = [7 46 188 22];
            app.ArrayElementsPositionsfilenameDropDownLabel.Text = 'Array Elements Positions filename';

            % Create ArrayElementsPositionsfilenameDropDown
            app.ArrayElementsPositionsfilenameDropDown = uidropdown(app.GeneralPanel_2);
            app.ArrayElementsPositionsfilenameDropDown.Items = {};
            app.ArrayElementsPositionsfilenameDropDown.Position = [210 46 100 22];
            app.ArrayElementsPositionsfilenameDropDown.Value = {};

            % Create SparsityfilenameDropDownLabel
            app.SparsityfilenameDropDownLabel = uilabel(app.GeneralPanel_2);
            app.SparsityfilenameDropDownLabel.HorizontalAlignment = 'right';
            app.SparsityfilenameDropDownLabel.Position = [98 17 97 22];
            app.SparsityfilenameDropDownLabel.Text = 'Sparsity filename';

            % Create SparsityfilenameDropDown
            app.SparsityfilenameDropDown = uidropdown(app.GeneralPanel_2);
            app.SparsityfilenameDropDown.Items = {};
            app.SparsityfilenameDropDown.Position = [210 17 100 22];
            app.SparsityfilenameDropDown.Value = {};

            % Create ElementGeometrySwitchLabel
            app.ElementGeometrySwitchLabel = uilabel(app.GeneralPanel_2);
            app.ElementGeometrySwitchLabel.HorizontalAlignment = 'center';
            app.ElementGeometrySwitchLabel.Position = [399 49 105 22];
            app.ElementGeometrySwitchLabel.Text = 'Element Geometry';

            % Create ElementGeometrySwitch
            app.ElementGeometrySwitch = uiswitch(app.GeneralPanel_2, 'slider');
            app.ElementGeometrySwitch.Items = {'Disc', 'Rect'};
            app.ElementGeometrySwitch.Position = [429 22 45 20];
            app.ElementGeometrySwitch.Value = 'Rect';

            % Create LengthmmEditFieldLabel
            app.LengthmmEditFieldLabel = uilabel(app.GeneralPanel_2);
            app.LengthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.LengthmmEditFieldLabel.Position = [529 47 73 22];
            app.LengthmmEditFieldLabel.Text = 'Length (mm)';

            % Create LengthmmEditField
            app.LengthmmEditField = uieditfield(app.GeneralPanel_2, 'numeric');
            app.LengthmmEditField.Position = [607 47 33 22];
            app.LengthmmEditField.Value = 3;

            % Create WidthmmEditFieldLabel
            app.WidthmmEditFieldLabel = uilabel(app.GeneralPanel_2);
            app.WidthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.WidthmmEditFieldLabel.Position = [536 17 67 22];
            app.WidthmmEditFieldLabel.Text = 'Width (mm)';

            % Create WidthmmEditField
            app.WidthmmEditField = uieditfield(app.GeneralPanel_2, 'numeric');
            app.WidthmmEditField.Position = [608 17 33 22];
            app.WidthmmEditField.Value = 3;

            % Create TargetingTab
            app.TargetingTab = uitab(app.TabGroup);
            app.TargetingTab.Title = 'Targeting';

            % Create Transducer1Panel_2
            app.Transducer1Panel_2 = uipanel(app.TargetingTab);
            app.Transducer1Panel_2.Title = 'Transducer 1';
            app.Transducer1Panel_2.Position = [291 278 383 229];

            % Create FocusPositionLabel
            app.FocusPositionLabel = uilabel(app.Transducer1Panel_2);
            app.FocusPositionLabel.Position = [149 179 84 22];
            app.FocusPositionLabel.Text = 'Focus Position';

            % Create xEditField_4Label
            app.xEditField_4Label = uilabel(app.Transducer1Panel_2);
            app.xEditField_4Label.HorizontalAlignment = 'right';
            app.xEditField_4Label.Position = [55 149 25 22];
            app.xEditField_4Label.Text = 'x';

            % Create focusxEditField
            app.focusxEditField = uieditfield(app.Transducer1Panel_2, 'numeric');
            app.focusxEditField.Position = [90 149 37 22];
            app.focusxEditField.Value = -18;

            % Create yEditField_4Label
            app.yEditField_4Label = uilabel(app.Transducer1Panel_2);
            app.yEditField_4Label.HorizontalAlignment = 'right';
            app.yEditField_4Label.Position = [135 149 25 22];
            app.yEditField_4Label.Text = 'y';

            % Create focusyEditField
            app.focusyEditField = uieditfield(app.Transducer1Panel_2, 'numeric');
            app.focusyEditField.Position = [170 149 37 22];
            app.focusyEditField.Value = 30;

            % Create zEditField_4Label
            app.zEditField_4Label = uilabel(app.Transducer1Panel_2);
            app.zEditField_4Label.HorizontalAlignment = 'right';
            app.zEditField_4Label.Position = [219 149 25 22];
            app.zEditField_4Label.Text = 'z';

            % Create focuszEditField
            app.focuszEditField = uieditfield(app.Transducer1Panel_2, 'numeric');
            app.focuszEditField.Position = [254 149 37 22];
            app.focuszEditField.Value = -27;

            % Create FocusRadiusmmEditFieldLabel
            app.FocusRadiusmmEditFieldLabel = uilabel(app.Transducer1Panel_2);
            app.FocusRadiusmmEditFieldLabel.HorizontalAlignment = 'right';
            app.FocusRadiusmmEditFieldLabel.Position = [43 80 110 22];
            app.FocusRadiusmmEditFieldLabel.Text = 'Focus Radius (mm)';

            % Create FocusRadiusmmEditField
            app.FocusRadiusmmEditField = uieditfield(app.Transducer1Panel_2, 'numeric');
            app.FocusRadiusmmEditField.Position = [164 80 40 22];
            app.FocusRadiusmmEditField.Value = 7;

            % Create PressureAmplitudekPaEditFieldLabel
            app.PressureAmplitudekPaEditFieldLabel = uilabel(app.Transducer1Panel_2);
            app.PressureAmplitudekPaEditFieldLabel.HorizontalAlignment = 'right';
            app.PressureAmplitudekPaEditFieldLabel.Position = [12 47 141 22];
            app.PressureAmplitudekPaEditFieldLabel.Text = 'Pressure Amplitude (kPa)';

            % Create PressureAmplitudekPaEditField
            app.PressureAmplitudekPaEditField = uieditfield(app.Transducer1Panel_2, 'numeric');
            app.PressureAmplitudekPaEditField.Position = [164 47 40 22];
            app.PressureAmplitudekPaEditField.Value = 300;

            % Create MinPointDistancemmEditFieldLabel
            app.MinPointDistancemmEditFieldLabel = uilabel(app.Transducer1Panel_2);
            app.MinPointDistancemmEditFieldLabel.HorizontalAlignment = 'right';
            app.MinPointDistancemmEditFieldLabel.Position = [7 17 140 22];
            app.MinPointDistancemmEditFieldLabel.Text = 'Min. Point Distance (mm)';

            % Create MinPointDistancemmEditField
            app.MinPointDistancemmEditField = uieditfield(app.Transducer1Panel_2, 'numeric');
            app.MinPointDistancemmEditField.Position = [163 17 43 22];
            app.MinPointDistancemmEditField.Value = 5;

            % Create todeclaredPressureSwitch_2Label
            app.todeclaredPressureSwitch_2Label = uilabel(app.Transducer1Panel_2);
            app.todeclaredPressureSwitch_2Label.HorizontalAlignment = 'center';
            app.todeclaredPressureSwitch_2Label.Position = [251 39 116 22];
            app.todeclaredPressureSwitch_2Label.Text = 'to declared Pressure';

            % Create todeclaredPressureSwitch_2
            app.todeclaredPressureSwitch_2 = uiswitch(app.Transducer1Panel_2, 'slider');
            app.todeclaredPressureSwitch_2.Items = {'Force', 'Limit'};
            app.todeclaredPressureSwitch_2.Position = [285 67 45 20];
            app.todeclaredPressureSwitch_2.Value = 'Force';

            % Create AddManualTargetButton
            app.AddManualTargetButton = uibutton(app.TargetingTab, 'push');
            app.AddManualTargetButton.Position = [112 401 122 23];
            app.AddManualTargetButton.Text = 'Add Manual Target';

            % Create RemoveManualTransducerButton
            app.RemoveManualTransducerButton = uibutton(app.TargetingTab, 'push');
            app.RemoveManualTransducerButton.Position = [90 371 166 23];
            app.RemoveManualTransducerButton.Text = 'Remove Manual Transducer';

            % Create ManualTargetDropDownLabel
            app.ManualTargetDropDownLabel = uilabel(app.TargetingTab);
            app.ManualTargetDropDownLabel.HorizontalAlignment = 'right';
            app.ManualTargetDropDownLabel.Position = [15 461 91 22];
            app.ManualTargetDropDownLabel.Text = 'Manual Target #';

            % Create ManualTargetDropDown
            app.ManualTargetDropDown = uidropdown(app.TargetingTab);
            app.ManualTargetDropDown.Items = {'1'};
            app.ManualTargetDropDown.Position = [121 461 100 22];
            app.ManualTargetDropDown.Value = '1';

            % Create Transducer1Panel_3
            app.Transducer1Panel_3 = uipanel(app.TargetingTab);
            app.Transducer1Panel_3.Title = 'Transducer 1';
            app.Transducer1Panel_3.Visible = 'off';
            app.Transducer1Panel_3.Position = [292 84 383 175];

            % Create BrainRegionDropDownLabel
            app.BrainRegionDropDownLabel = uilabel(app.Transducer1Panel_3);
            app.BrainRegionDropDownLabel.HorizontalAlignment = 'right';
            app.BrainRegionDropDownLabel.Position = [79 115 74 22];
            app.BrainRegionDropDownLabel.Text = 'Brain Region';

            % Create BrainRegionDropDown
            app.BrainRegionDropDown = uidropdown(app.Transducer1Panel_3);
            app.BrainRegionDropDown.Items = {};
            app.BrainRegionDropDown.Position = [168 115 100 22];
            app.BrainRegionDropDown.Value = {};

            % Create MinPointDistancemmEditField_2Label
            app.MinPointDistancemmEditField_2Label = uilabel(app.Transducer1Panel_3);
            app.MinPointDistancemmEditField_2Label.HorizontalAlignment = 'right';
            app.MinPointDistancemmEditField_2Label.Position = [15 34 140 22];
            app.MinPointDistancemmEditField_2Label.Text = 'Min. Point Distance (mm)';

            % Create MinPointDistancemmEditField_2
            app.MinPointDistancemmEditField_2 = uieditfield(app.Transducer1Panel_3, 'numeric');
            app.MinPointDistancemmEditField_2.Position = [171 34 43 22];
            app.MinPointDistancemmEditField_2.Value = 5;

            % Create todeclaredPressureSwitchLabel
            app.todeclaredPressureSwitchLabel = uilabel(app.Transducer1Panel_3);
            app.todeclaredPressureSwitchLabel.HorizontalAlignment = 'center';
            app.todeclaredPressureSwitchLabel.Position = [253 38 116 22];
            app.todeclaredPressureSwitchLabel.Text = 'to declared Pressure';

            % Create todeclaredPressureSwitch
            app.todeclaredPressureSwitch = uiswitch(app.Transducer1Panel_3, 'slider');
            app.todeclaredPressureSwitch.Items = {'Force', 'Limit'};
            app.todeclaredPressureSwitch.Position = [287 67 45 20];
            app.todeclaredPressureSwitch.Value = 'Force';

            % Create PressureAmplitudekPaEditField_2Label
            app.PressureAmplitudekPaEditField_2Label = uilabel(app.Transducer1Panel_3);
            app.PressureAmplitudekPaEditField_2Label.HorizontalAlignment = 'right';
            app.PressureAmplitudekPaEditField_2Label.Position = [20 65 141 22];
            app.PressureAmplitudekPaEditField_2Label.Text = 'Pressure Amplitude (kPa)';

            % Create PressureAmplitudekPaEditField_2
            app.PressureAmplitudekPaEditField_2 = uieditfield(app.Transducer1Panel_3, 'numeric');
            app.PressureAmplitudekPaEditField_2.Position = [172 65 40 22];
            app.PressureAmplitudekPaEditField_2.Value = 300;

            % Create AddRegionTargetButton
            app.AddRegionTargetButton = uibutton(app.TargetingTab, 'push');
            app.AddRegionTargetButton.Position = [113 153 122 23];
            app.AddRegionTargetButton.Text = 'Add Region Target';

            % Create RemoveRegionTargetButton
            app.RemoveRegionTargetButton = uibutton(app.TargetingTab, 'push');
            app.RemoveRegionTargetButton.Position = [105 123 138 23];
            app.RemoveRegionTargetButton.Text = 'Remove Region Target';

            % Create RegionTargetDropDownLabel
            app.RegionTargetDropDownLabel = uilabel(app.TargetingTab);
            app.RegionTargetDropDownLabel.HorizontalAlignment = 'right';
            app.RegionTargetDropDownLabel.Position = [17 213 90 22];
            app.RegionTargetDropDownLabel.Text = 'Region Target #';

            % Create RegionTargetDropDown
            app.RegionTargetDropDown = uidropdown(app.TargetingTab);
            app.RegionTargetDropDown.Items = {};
            app.RegionTargetDropDown.Position = [122 213 100 22];
            app.RegionTargetDropDown.Value = {};

            % Create UpdateButtonTargeting
            app.UpdateButtonTargeting = uibutton(app.TargetingTab, 'push');
            app.UpdateButtonTargeting.Position = [574 24 100 23];
            app.UpdateButtonTargeting.Text = 'Update';

            % Create MaxPressurekPaEditFieldLabel
            app.MaxPressurekPaEditFieldLabel = uilabel(app.TargetingTab);
            app.MaxPressurekPaEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxPressurekPaEditFieldLabel.Position = [66 24 114 22];
            app.MaxPressurekPaEditFieldLabel.Text = 'Max. Pressure (kPa)';

            % Create MaxPressurekPaEditField
            app.MaxPressurekPaEditField = uieditfield(app.TargetingTab, 'numeric');
            app.MaxPressurekPaEditField.Position = [191 24 44 22];
            app.MaxPressurekPaEditField.Value = 200;

            % Create LimitIntracranialOffTargetPressureCheckBox
            app.LimitIntracranialOffTargetPressureCheckBox = uicheckbox(app.TargetingTab);
            app.LimitIntracranialOffTargetPressureCheckBox.Text = 'Limit Intracranial Off-Target Pressure';
            app.LimitIntracranialOffTargetPressureCheckBox.Position = [50 50 219 22];
            app.LimitIntracranialOffTargetPressureCheckBox.Value = true;

            % Create LimitSkullPressureCheckBox
            app.LimitSkullPressureCheckBox = uicheckbox(app.TargetingTab);
            app.LimitSkullPressureCheckBox.Text = 'Limit Skull Pressure';
            app.LimitSkullPressureCheckBox.Position = [333 50 128 22];

            % Create MaxPressurekPaEditField_2Label
            app.MaxPressurekPaEditField_2Label = uilabel(app.TargetingTab);
            app.MaxPressurekPaEditField_2Label.HorizontalAlignment = 'right';
            app.MaxPressurekPaEditField_2Label.Visible = 'off';
            app.MaxPressurekPaEditField_2Label.Position = [333 24 114 22];
            app.MaxPressurekPaEditField_2Label.Text = 'Max. Pressure (kPa)';

            % Create MaxPressurekPaEditField_2
            app.MaxPressurekPaEditField_2 = uieditfield(app.TargetingTab, 'numeric');
            app.MaxPressurekPaEditField_2.Visible = 'off';
            app.MaxPressurekPaEditField_2.Position = [458 24 44 22];
            app.MaxPressurekPaEditField_2.Value = 1000;

            % Create OptimizeTab
            app.OptimizeTab = uitab(app.TabGroup);
            app.OptimizeTab.Title = 'Optimize';

            % Create DoGroundTruthSimulationButton
            app.DoGroundTruthSimulationButton = uibutton(app.OptimizeTab, 'push');
            app.DoGroundTruthSimulationButton.Position = [545 478 164 23];
            app.DoGroundTruthSimulationButton.Text = 'Do Ground Truth Simulation';

            % Create ComputeInitialSolutionButton
            app.ComputeInitialSolutionButton = uibutton(app.OptimizeTab, 'push');
            app.ComputeInitialSolutionButton.Position = [329 478 142 23];
            app.ComputeInitialSolutionButton.Text = 'Compute Initial Solution';

            % Create OptimizeButton
            app.OptimizeButton = uibutton(app.OptimizeTab, 'push');
            app.OptimizeButton.Position = [329 423 142 45];
            app.OptimizeButton.Text = 'Optimize';

            % Create OptimizationModeButtonGroup
            app.OptimizationModeButtonGroup = uibuttongroup(app.OptimizeTab);
            app.OptimizationModeButtonGroup.Title = 'Optimization Mode';
            app.OptimizationModeButtonGroup.Position = [12 412 224 106];

            % Create OptimizeforTransducerPhasesandAmplitudesButton
            app.OptimizeforTransducerPhasesandAmplitudesButton = uiradiobutton(app.OptimizationModeButtonGroup);
            app.OptimizeforTransducerPhasesandAmplitudesButton.Text = {'Optimize for Transducer '; 'Phases and Amplitudes'};
            app.OptimizeforTransducerPhasesandAmplitudesButton.Position = [11 52 153 30];
            app.OptimizeforTransducerPhasesandAmplitudesButton.Value = true;

            % Create OptimizeforPhasesand1AmpperTransducerButton
            app.OptimizeforPhasesand1AmpperTransducerButton = uiradiobutton(app.OptimizationModeButtonGroup);
            app.OptimizeforPhasesand1AmpperTransducerButton.Text = {'Optimize for Phases '; 'and 1 Amp per Transducer'};
            app.OptimizeforPhasesand1AmpperTransducerButton.Position = [11 15 163 30];

            % Create PlotEntireDomainCheckBox
            app.PlotEntireDomainCheckBox = uicheckbox(app.OptimizeTab);
            app.PlotEntireDomainCheckBox.Text = 'Plot Entire Domain';
            app.PlotEntireDomainCheckBox.Position = [15 377 122 22];

            % Create PlotSkullCheckBox
            app.PlotSkullCheckBox = uicheckbox(app.OptimizeTab);
            app.PlotSkullCheckBox.Text = 'Plot Skull';
            app.PlotSkullCheckBox.Position = [163 377 72 22];

            % Create GroundTruthResolutionFactorEditFieldLabel
            app.GroundTruthResolutionFactorEditFieldLabel = uilabel(app.OptimizeTab);
            app.GroundTruthResolutionFactorEditFieldLabel.HorizontalAlignment = 'right';
            app.GroundTruthResolutionFactorEditFieldLabel.Position = [562 435 99 30];
            app.GroundTruthResolutionFactorEditFieldLabel.Text = {'Ground Truth '; 'Resolution Factor'};

            % Create GroundTruthResolutionFactorEditField
            app.GroundTruthResolutionFactorEditField = uieditfield(app.OptimizeTab, 'numeric');
            app.GroundTruthResolutionFactorEditField.Position = [673 443 34 22];
            app.GroundTruthResolutionFactorEditField.Value = 1;

            % Create UpdateSliceButton
            app.UpdateSliceButton = uibutton(app.UIFigure, 'push');
            app.UpdateSliceButton.Position = [868 25 100 23];
            app.UpdateSliceButton.Text = 'Update Slice';

            % Create SlicexLabel
            app.SlicexLabel = uilabel(app.UIFigure);
            app.SlicexLabel.Position = [898 71 40 22];
            app.SlicexLabel.Text = 'Slice x';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = simulationApp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end