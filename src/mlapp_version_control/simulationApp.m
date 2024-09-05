classdef simulationApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        SliceDirectionDropDown_2        matlab.ui.control.DropDown
        SliceDirectionLabel             matlab.ui.control.Label
        UpdateSliceButton               matlab.ui.control.Button
        CloseallPlotsButton             matlab.ui.control.Button
        Slice30Label                    matlab.ui.control.Label
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
        LoadResultsfromFileDropDown     matlab.ui.control.DropDown
        LoadResultsfromFileDropDownLabel  matlab.ui.control.Label
        SaveResultsCheckBox             matlab.ui.control.CheckBox
        SliceIndexEditField             matlab.ui.control.NumericEditField
        SliceIndexLabel                 matlab.ui.control.Label
        SliceDimDropDown                matlab.ui.control.DropDown
        SliceDimDropDownLabel           matlab.ui.control.Label
        SimSpatialResolutionmmEditField  matlab.ui.control.NumericEditField
        SimSpatialResolutionmmEditFieldLabel  matlab.ui.control.Label
        CenterFreqkHzEditField          matlab.ui.control.NumericEditField
        CenterFreqkHzEditFieldLabel     matlab.ui.control.Label
        DimSwitch                       matlab.ui.control.Switch
        Label                           matlab.ui.control.Label
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
        TransducerPreviewButton         matlab.ui.control.Button
        RemoveTransducerButton          matlab.ui.control.Button
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
        AddTransducerButton             matlab.ui.control.Button
        Transducer1Panel                matlab.ui.container.Panel
        ConfirmButton                   matlab.ui.control.Button
        TransducerLengthmmEditField     matlab.ui.control.NumericEditField
        TransducerLengthmmEditFieldLabel  matlab.ui.control.Label
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
        MaxPressurekPaSkullEditField    matlab.ui.control.NumericEditField
        MaxPressurekPaSkullEditFieldLabel  matlab.ui.control.Label
        LimitSkullPressureCheckBox      matlab.ui.control.CheckBox
        LimitIntracranialOffTargetPressureCheckBox  matlab.ui.control.CheckBox
        MaxPressurekPaEditField         matlab.ui.control.NumericEditField
        MaxPressurekPaEditFieldLabel    matlab.ui.control.Label
        UpdateTargetingButton           matlab.ui.control.Button
        RegionTargetDropDown            matlab.ui.control.DropDown
        RegionTargetDropDownLabel       matlab.ui.control.Label
        RemoveRegionTargetButton        matlab.ui.control.Button
        AddRegionTargetButton           matlab.ui.control.Button
        TargetRegPanel                  matlab.ui.container.Panel
        ConfirmRegTargetButton          matlab.ui.control.Button
        PressureAmplitudekPaEditFieldReg  matlab.ui.control.NumericEditField
        PressureAmplitudekPaEditField_2Label  matlab.ui.control.Label
        todeclaredPressureSwitchReg     matlab.ui.control.Switch
        todeclaredPressureSwitchLabel   matlab.ui.control.Label
        MinPointDistancemmEditFieldReg  matlab.ui.control.NumericEditField
        MinPointDistancemmEditField_2Label  matlab.ui.control.Label
        BrainRegionDropDown             matlab.ui.control.DropDown
        BrainRegionDropDownLabel        matlab.ui.control.Label
        ManualTargetDropDown            matlab.ui.control.DropDown
        ManualTargetDropDownLabel       matlab.ui.control.Label
        RemoveManualTargetButton        matlab.ui.control.Button
        AddManualTargetButton           matlab.ui.control.Button
        TargetManPanel                  matlab.ui.container.Panel
        ConfirmManTargetButton          matlab.ui.control.Button
        todeclaredPressureSwitch        matlab.ui.control.Switch
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
        PlotResultsButton               matlab.ui.control.Button
        opt_mode                        matlab.ui.control.CheckBox
        PlotPanel                       matlab.ui.container.Panel
        MaskPressurePlotkPaEditField    matlab.ui.control.NumericEditField
        MaskPressurePlotkPaEditFieldLabel  matlab.ui.control.Label
        LimittoSkullCheckBox            matlab.ui.control.CheckBox
        LimittoIntracranialFieldCheckBox  matlab.ui.control.CheckBox
        AdvancedOptimizationOptionsPanel  matlab.ui.container.Panel
        LxNormEditField                 matlab.ui.control.NumericEditField
        LxNormEditFieldLabel            matlab.ui.control.Label
        MinimizeArrrayAmplitudesCheckBox  matlab.ui.control.CheckBox
        MaxNoofIterationsEditField      matlab.ui.control.NumericEditField
        MaxNoofIterationsEditFieldLabel  matlab.ui.control.Label
        ConstraintToleranceEditField    matlab.ui.control.NumericEditField
        ConstraintToleranceLabel        matlab.ui.control.Label
        RelFunctionToleranceEditField   matlab.ui.control.NumericEditField
        RelFunctionToleranceEditFieldLabel  matlab.ui.control.Label
        DisplayAdvancedOptimizationOptionsCheckBox  matlab.ui.control.CheckBox
        GroundTruthResolutionFactorEditField  matlab.ui.control.NumericEditField
        GroundTruthResolutionFactorEditFieldLabel  matlab.ui.control.Label
        OptimizationModeButtonGroup     matlab.ui.container.ButtonGroup
        OptimizeforPhasesand1AmpButton  matlab.ui.control.RadioButton
        OptimizeforPhasesand1AmpperTransducerButton  matlab.ui.control.RadioButton
        OptimizeforPhasesandAmplitudesButton  matlab.ui.control.RadioButton
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
        dx_factor_init
        dx_factor
        grid_size
        t1w_filename
        ct_filename
        dx_scan
        plot_offset
        tr_offset_karr 
        segment_ids
        segment_labels
        slice_grid_2D
        logical_dom_ids
        input_args
        use_greens_fctn
        current_datetime
        t_pos
        t_rot
        tr_len
        el_per_t
        t_mask_ps
        mask2el
        karray_t
        preplot_arg
        active_ids
        point_pos_m
        point_pos
        focus_radius
        des_pressures
        min_dist
        force_pressures
        tar_reg_labels
        des_pressures_reg
        min_dist_reg
        force_pressures_reg
        b_mask
        full_bmask
        b_des
        b_ip_des
        ip
        A
        skull_ids
        vol_ids
        init_ids
        p_curr
        plot_title_curr
        slice_dim
        dims_2D
        sv_obj
    end
    
    methods (Access = private)
        
        function entries = getDropdownEntries(app, path, patterns, file_end)
            items = dir(path);

            % Check if each item matches any of the patterns
            matches = false(size(items))';
            for i = 1:length(patterns)
                matches = matches | (contains({items.name}, patterns{i}) & contains({items.name}, file_end));
            end

            % List entries
            entries = {items(matches).name};
        end
        
        function evaluate_results(app, p, plot_title)
            %% Obtain pressure distribution
            b = app.A * p;
            b = reshape(b, size(app.kgrid.k));

            %% Plot results
            b_lim = b;
            if app.LimittoIntracranialFieldCheckBox.Value || app.LimittoSkullCheckBox.Value
                excl_ids = ~(app.LimittoIntracranialFieldCheckBox.Value & app.logical_dom_ids) & ...
                ~(app.LimittoSkullCheckBox.Value & app.skull_ids);

                b_lim(excl_ids) = 0.0;
            end

            plot_thr = app.MaskPressurePlotkPaEditField.Value * 1e3;
            save_results = app.SaveResultsCheckBox.Value;

            app.sv_obj = plot_results(app.kgrid, p, b_lim, plot_title, app.mask2el, app.t1w_filename, app.plot_offset, ...
                app.grid_size, app.dx_factor, save_results, app.current_datetime, 'slice', app.SliceIndexEditField.Value, ...
                'slice_dim', app.SliceDimDropDown.Value, 'axes', []);%app.UIAxesParam);

            % Plot mask with pressure above off-target limit
            masked_b = abs(b_lim);
            masked_b(masked_b <= plot_thr) = 0.0;
            plot_results(app.kgrid, [], masked_b, strcat(plot_title, ' Mask'), app.mask2el, app.t1w_filename, ...
                app.plot_offset, app.grid_size, app.dx_factor, save_results, app.current_datetime, 'slice', ...
                app.SliceIndexEditField.Value, 'fig_pos', {[], [1400 450 475 525], [1660 55 250 300]}, ...
                'slice_dim', app.SliceDimDropDown.Value);

            %% Evaluate pressure distribution
            app.p_curr = p;
            app.plot_title_curr = plot_title;

            evaluate_metrics(app, plot_title, b, app.vol_ids, app.logical_dom_ids, app.skull_ids, app.init_ids);

            % Save results in mat file
            app.ip.b = b;
            save_results_mat(app);
        end
        
        function get_initial_solution(app)
            app.ip.p_init = pinv(app.A(app.init_ids, :)) * app.b_ip_des(app.init_ids);
        end
        
        function [kgrid_out, medium_out, sensor_out, sensor_mask_out, dx_factor_out, seg_ids_out, log_dom_ids_out, ...
                tr_offset_karr_out, slice_grid_2D_out] = init_general(app, dx_factor_in)
            f0 = app.CenterFreqkHzEditField.Value * 1e3; % Hz
            app.n_dim = 2 + (app.DimSwitch.Value == "3D");

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
            if isempty(dx_factor_in)
                SimSpatialResolutionmmEditFieldValueChanged(app);
                dx_factor_in = app.dx_factor_init;
            end

            if strcmp(app.MediumSwitch.Value, 'Homogeneous')
                app.grid_size = [app.xHomEditField.Value, app.HomyEditField.Value, app.HomzEditField.Value] * 1e-3;
                app.t1w_filename = [];
                app.ct_filename = [];
                app.dx_scan = 1e-3; % m

                app.plot_offset = app.grid_size / app.dx_scan / 2 + 1; % Offset to center
            else
                app.grid_size = [];
                app.t1w_filename = fullfile('..', 'Scans', app.T1wfilenameDropDown.Value);
                app.ct_filename = fullfile('..', 'Scans', app.CTfilenameDropDown.Value);
                app.dx_scan = app.ScanSpatialResolutionmmEditField.Value * 1e-3;

                app.plot_offset = [-app.scanxEditField.Value, -app.scanyEditField.Value, -app.scanzEditField.Value] + 1;
            end

            app.slice_dim = dim2num(app.SliceDimDropDown.Value);
            scan_slice = round(app.plot_offset(app.slice_dim) + app.SliceIndexEditField.Value);
            app.dims_2D = exclude_dim(app.slice_dim);

            [kgrid_out, medium_out, app.grid_size, ppp] = init_grid_medium(f0, app.grid_size, app.SliceDimDropDown.Value, ...
                'n_dim', app.n_dim, 'dx_factor', dx_factor_in, 'ct_scan', app.ct_filename, ...
                'slice_idx', scan_slice, 'dx_scan', app.dx_scan, 'constants', const);
            [sensor_out, sensor_mask_out] = init_sensor(kgrid_out, ppp);

            dx_factor_out = app.dx_scan / kgrid_out.dx;

            %% Segment the brain
            if ~isempty(app.t1w_filename)
                [seg_ids_out] = segment_space(app.t1w_filename, app.dx_scan);
                app.segment_labels = unique(seg_ids_out(:));

                if abs(dx_factor_out) ~= 1.0
                    % Interpolate to adapt to grid size
                    grid_sz_vox = size(kgrid_out.k);
                    seg_sz = size(seg_ids_out);
                    [uniqueStrings, ~, seg_nums] = unique(seg_ids_out);
                    seg_nums = reshape(seg_nums, size(seg_ids_out)); % Ensure it has the same shape as the original 3D array
                
                    if app.n_dim == 2
                        grid_sz_vox = grid_sz_vox(app.dims_2D);
                        seg_sz = seg_sz(app.dims_2D);
                        seg_nums = index2Dto3D(seg_nums, app.SliceDimDropDown.Value, scan_slice);

                        [X, Z] = meshgrid(1:seg_sz(1), 1:seg_sz(2));
                        [Xq, Zq] = meshgrid(linspace(1, seg_sz(1), grid_sz_vox(1)), linspace(1, seg_sz(2), grid_sz_vox(2)));
                        seg_nums = interp2(X, Z, double(seg_nums)', Xq, Zq, "nearest")';
                    else
                        [X, Y, Z] = meshgrid(1:seg_sz(1), 1:seg_sz(2), 1:seg_sz(3));
                        [Xq, Yq, Zq] = meshgrid(linspace(1, seg_sz(1), grid_sz_vox(1)), linspace(1, seg_sz(2), grid_sz_vox(2)), ...
                            linspace(1, seg_sz(3), grid_sz_vox(3)));
                        seg_nums = permute(interp3(X, Y, Z, permute(double(seg_nums), [2 1 3]), Xq, Yq, Zq, "nearest"), [2 1 3]);
                    end
                
                    % Map back to strings
                    seg_nums = round(seg_nums); % Ensure indices are integers
                    seg_ids_out = uniqueStrings(seg_nums);
                end
                domain_ids = seg_ids_out ~= "background"; % Mask entire brain
            else
                seg_ids_out = [];
                domain_ids = ones(size(kgrid_out.k));
            end

            log_dom_ids_out = false(numel(medium_out.sound_speed), 1);
            if app.n_dim == 3
                tr_offset_karr_out = ((app.plot_offset - 1) * app.dx_scan - app.grid_size / 2)'; % karray offset in m

                slice_grid_2D_out = [];
                log_dom_ids_out(domain_ids) = true;
            else
                tr_offset_karr_out = [];

                slice_grid_2D_out = round((app.plot_offset(app.slice_dim) + app.SliceIndexEditField.Value) * dx_factor_out);
                log_dom_ids_out(index2Dto3D(domain_ids, app.SliceDimDropDown.Value, slice_grid_2D_out)) = true;
            end

            %% Set global options
            app.input_args = {'PMLSize', app.PMLsizeEditField.Value, 'PMLInside', app.PMLinsideCheckBox.Value, ...
                'PlotPML', true, 'DisplayMask', 'off', 'RecordMovie', false};
            app.use_greens_fctn = app.GreensFunctionbasedCheckBox.Value & max(medium_out.sound_speed(:)) ...
                == min(medium_out.sound_speed(:)); % Update green's fctn flag

            disp("Init successful")
        end
        
        function [t_mask_ps_out, karray_t_out, el_per_t_out, active_ids_out] = transducer_geometry_init(app, kgrid_in, ...
                dx_factor_in, tr_offset_karr_in)

            n_trs = length(app.TransducerDropDown.Items);

            if app.n_dim == 2
                spacing = 1;
            
                t_mask_ps_out = false(kgrid_in.Nx, kgrid_in.Ny);
                el_per_t_out = zeros(1, n_trs);
                t_ids = [];
                tr_len_grid = app.tr_len * 1e-3 / kgrid_in.dx;

                for i = 1:n_trs
                    x_offset = round((app.plot_offset(app.dims_2D(1)) + app.t_pos(app.dims_2D(1), i)) * dx_factor_in); % grid points
                    y_offset = round((app.plot_offset(app.dims_2D(2)) + app.t_pos(app.dims_2D(2), i)) * dx_factor_in); 
                    % tangential shift in grid points
                
                    new_arr = create_linear_array(kgrid_in, tr_len_grid(i), x_offset, y_offset, spacing, app.t_rot(app.slice_dim, i));
            
                    el_per_t_out(i) = sum(new_arr(:));
                    t_ids = [t_ids; find(new_arr)];
                    t_mask_ps_out = t_mask_ps_out | logical(new_arr);
                end
                
                [~, el2mask_ids] = sort(t_ids);
                [~, app.mask2el] = sort(el2mask_ids);
            
                karray_t_out = [];
                active_ids_out = [];
            else
                % Planar Array
                t_name = app.ArrayElementsPositionsfilenameDropDown.Value(1:end-4);
                sparsity_name = app.SparsityfilenameDropDown.Value(1:end-4);

                t_pos_3D = app.t_pos * app.dx_scan;
                active_tr_ids = 1:n_trs;

                if strcmp(app.ElementGeometrySwitch.Value, 'Rect')
                    el_sz = [app.LengthmmEditField.Value, app.WidthmmEditField.Value];
                else
                    el_sz = app.LengthmmEditField.Value;
                end
            
                [karray_t_out, t_mask_ps_out, active_ids_out, num_elements, app.mask2el] = create_transducer(kgrid_in, ...
                    app.plot_offset, tr_offset_karr_in, t_name, sparsity_name, t_pos_3D, app.t_rot, active_tr_ids, el_sz * 1e-3);
            
                el_per_t_out = num_elements * ones(1, length(active_tr_ids));
            end
        end
        
        function save_results_mat(app)
            res_filename = "results";
            if app.SaveResultsCheckBox.Value

                % Std param
                dim_switch = app.DimSwitch.Value;
                f0 = app.CenterFreqkHzEditField.Value;
                slice_direction = app.SliceDimDropDown.Value;
                slice_index = app.SliceIndexEditField.Value;
                dx = app.SimSpatialResolutionmmEditField.Value;

                medium_type = app.MediumSwitch.Value; MediumSwitchValueChanged(app);
                green_function = app.GreensFunctionbasedCheckBox.Value;
                hom_coord = [app.xHomEditField.Value, app.HomyEditField.Value, app.HomzEditField.Value];
                t1w_file = app.T1wfilenameDropDown.Value;
                ct_file = app.CTfilenameDropDown.Value;
                het_coord = [app.scanxEditField.Value, app.scanyEditField.Value, app.scanzEditField.Value];
                dx_t1w = app.ScanSpatialResolutionmmEditField.Value;
                seg_labels = app.segment_labels;
                
                % Transducer param
                TransducersTabButtonDown(app);
                dropdowns = [];
                if ~isempty(app.ArrayElementsPositionsfilenameDropDown.Items)
                    dropdowns.array_pos = app.ArrayElementsPositionsfilenameDropDown.Value;
                end
                if ~isempty(app.SparsityfilenameDropDown.Items)
                    dropdowns.sparsity = app.SparsityfilenameDropDown.Value;
                end

                array.geom = app.ElementGeometrySwitch.Value;
                array.len = app.LengthmmEditField.Value;
                array.wid = app.WidthmmEditField.Value;
                tr.pos = app.t_pos;
                tr.rot = app.t_rot;
                tr.length = app.tr_len;
                matrix_name = app.PropagationMatrixAfilenameDropDown.Value;

                % Target param
                man.points = app.point_pos_m;
                man.radius = app.focus_radius;
                man.pressure = app.des_pressures;
                man.min_dist = app.min_dist;
                man.force = app.force_pressures;
                reg.labels = app.tar_reg_labels;
                reg.pressure = app.des_pressures_reg;
                reg.min_dist = app.min_dist_reg;
                reg.force = app.force_pressures_reg;

                ineq_flag = app.LimitIntracranialOffTargetPressureCheckBox.Value;
                ineq_val = app.MaxPressurekPaEditField.Value;
                skull_flag = app.LimitSkullPressureCheckBox.Value;
                skull_val = app.MaxPressurekPaSkullEditField.Value;

                % Optimization param
                opt.phaseAmp = app.OptimizeforPhasesandAmplitudesButton.Value;
                opt.phaseOnly = app.OptimizeforPhasesand1AmpperTransducerButton.Value;
                opt.intracranial = app.LimittoIntracranialFieldCheckBox.Value;
                opt.skull = app.LimittoSkullCheckBox.Value;
                gt_dx_factor = app.GroundTruthResolutionFactorEditField.Value;

                term.advanced_opt = app.DisplayAdvancedOptimizationOptionsCheckBox.Value;
                term.min_arr_amp = app.MinimizeArrrayAmplitudesCheckBox.Value;
                term.fun_tol = app.RelFunctionToleranceEditField.Value;
                term.constr_tol = app.ConstraintToleranceEditField.Value;
                term.max_iter = app.MaxNoofIterationsEditField.Value;
                term.algorithm = app.AlgorithmDropDown.Value;

                ip_sol = app.ip;

                save(fullfile("..", "Results", app.current_datetime + "_" + res_filename + ".mat"), "dim_switch", "f0", "slice_direction", "slice_index", "dx", ...
                    "medium_type", "green_function", "hom_coord", "t1w_file", "ct_file", "het_coord", "dx_t1w", "dropdowns", "array", "tr", "matrix_name", ...
                    "man", "reg", "ineq_flag", "ineq_val", "skull_flag", "skull_val", "opt", "gt_dx_factor", "term", "ip_sol", "seg_labels");
                
            end
        end
        
        function evaluate_metrics(app, plot_title, b, vol_ids, logical_dom_ids, skull_ids, init_ids)
            real_ip = abs(reshape(b, [], 1));
            
            offTar_real_ip = real_ip;
            offTar_real_ip(vol_ids | ~logical_dom_ids) = [];
            
            skull_real_ip = real_ip(skull_ids);
            
            init_real_ip = real_ip(init_ids);


            fprintf(strcat("\n", plot_title, " Excitation L2-norm (kPa):\n"))
            disp(norm(app.p_curr) * 1e-3)

            fprintf(strcat("\n", plot_title, " max. Pressure (kPa):\n"))
            disp(max(real_ip) * 1e-3)

            fprintf(strcat("\n", plot_title, " max. Skull Pressure (kPa):\n"))
            disp(max(skull_real_ip) * 1e-3)
            
            fprintf(strcat("\n", plot_title, " Forced Points (kPa):\n"))
            disp(init_real_ip' * 1e-3)
            
            fprintf(strcat("\n", plot_title, " max. Off-Target Pressure (kPa):\n"))
            disp(max(offTar_real_ip) * 1e-3)
        end
        
        function labels2D_visible(app, visible_3D, dim)
            
            if dim == 1
                app.xHomEditField.Visible = visible_3D;
                app.xEditField_2Label.Visible = visible_3D;
        
                app.trPosxEditField.Visible = visible_3D;
                app.xEditField_3Label.Visible = visible_3D;
                app.betaEditField.Visible = visible_3D;
                app.betaEditFieldLabel.Visible = visible_3D;
                app.gammaEditField.Visible = visible_3D;
                app.gammaEditFieldLabel.Visible = visible_3D;
        
                app.focusxEditField.Visible = visible_3D;
                app.xEditField_4Label.Visible = visible_3D;
        
            elseif dim == 2
                app.HomyEditField.Visible = visible_3D;
                app.yEditField_2Label.Visible = visible_3D;
        
                app.trPosyEditField.Visible = visible_3D;
                app.yEditField_3Label.Visible = visible_3D;
                app.alphaEditField.Visible = visible_3D;
                app.alphaEditFieldLabel.Visible = visible_3D;
                app.gammaEditField.Visible = visible_3D;
                app.gammaEditFieldLabel.Visible = visible_3D;
        
                app.focusyEditField.Visible = visible_3D;
                app.yEditField_4Label.Visible = visible_3D;
        
            else
                app.HomzEditField.Visible = visible_3D;
                app.zEditField_2Label.Visible = visible_3D;
        
                app.trPoszEditField.Visible = visible_3D;
                app.zEditField_3Label.Visible = visible_3D;
                app.alphaEditField.Visible = visible_3D;
                app.alphaEditFieldLabel.Visible = visible_3D;
                app.betaEditField.Visible = visible_3D;
                app.betaEditFieldLabel.Visible = visible_3D;
        
                app.focuszEditField.Visible = visible_3D;
                app.zEditField_4Label.Visible = visible_3D;
            end
        
            % General properties 
            app.UpdateSliceButton.Visible = visible_3D;
            app.SliceDirectionDropDown_2.Visible = visible_3D;
            app.SliceDirectionLabel.Visible = visible_3D;
            app.SliceDimDropDown.Visible = ~visible_3D;
            app.SliceDimDropDownLabel.Visible = ~visible_3D;

            app.TransducerLengthmmEditField.Visible = ~visible_3D;
            app.TransducerLengthmmEditFieldLabel.Visible = ~visible_3D;
        
            app.GeneralPanel_2.Visible = visible_3D;
            app.GroundTruthResolutionFactorEditField.Visible = visible_3D;
            app.GroundTruthResolutionFactorEditFieldLabel.Visible = visible_3D;
        end

    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: UpdateButtonInit
        function UpdateButtonInitPushed(app, event)
            app.current_datetime = string(datestr(now, 'yyyymmddHHMMSS'));
            [app.kgrid, app.medium, app.sensor, app.sensor_mask, app.dx_factor, app.segment_ids, app.logical_dom_ids, ...
                app.tr_offset_karr, app.slice_grid_2D] = init_general(app, []);

            %% Show Skull and Scan in Preview
            if ~strcmp(app.MediumSwitch.Value, 'Homogeneous')
                skull_arg = app.medium.sound_speed - min(app.medium.sound_speed(:));
                skull_arg = skull_arg / max(skull_arg(:));

                app.sv_obj = plot_results(app.kgrid, [], skull_arg, 'Scan/Skull Preview', [], app.t1w_filename, ...
                    app.plot_offset, app.grid_size, app.dx_factor, false, [], 'slice', app.SliceIndexEditField.Value, ...
                    'colorbar', false, 'cmap', hot(), 'slice_dim', app.SliceDimDropDown.Value);
            end
        end

        % Value changed function: SimSpatialResolutionmmEditField
        function SimSpatialResolutionmmEditFieldValueChanged(app, event)
            value = app.SimSpatialResolutionmmEditField.Value * 1e-3;
            
            dx_std = app.c0msEditField.Value / (app.CenterFreqkHzEditField.Value * 1e3) / app.ppwEditField.Value;
            if value <= 0 || dx_std < value
                app.dx_factor_init = 1;
                dx = dx_std;
                if dx_std < value
                    app.RealdxLabel.Text = strcat("Spatial Aliasing Warning! -> Real: ", num2str(dx * 1e3, 3), " mm");
                    return;
                end
            else
                app.dx_factor_init = dx_std / value;
                dx = value;
            end

            app.RealdxLabel.Text = strcat("-> Real: ", num2str(dx * 1e3, 3), " mm");
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

                % Exclude Region targets
                app.RegionTargetDropDown.Visible = false;
                app.RegionTargetDropDownLabel.Visible = false;
                app.AddRegionTargetButton.Visible = false;
                app.RemoveRegionTargetButton.Visible = false;
            else
                app.homPanel.Visible = false;

                % Create t1w and ct dropdown menu
                scan_path = fullfile('..', 'Scans');
                t1_pattern = {'T1', 't1'};
                ct_pattern = {'ct', 'CT'};
                scan_pattern = '.nii';

                % Assign to t1 and ct dropdown menus
                app.T1wfilenameDropDown.Items = getDropdownEntries(app, scan_path, t1_pattern, scan_pattern);
                app.CTfilenameDropDown.Items = getDropdownEntries(app, scan_path, ct_pattern, scan_pattern);

                % Include Region targets
                app.RegionTargetDropDown.Visible = true;
                app.RegionTargetDropDownLabel.Visible = true;
                app.AddRegionTargetButton.Visible = true;
                app.RemoveRegionTargetButton.Visible = true;
            end
        end

        % Value changed function: CenterFreqkHzEditField
        function CenterFreqkHzEditFieldValueChanged(app, event)
            SimSpatialResolutionmmEditFieldValueChanged(app);
        end

        % Button down function: TransducersTab
        function TransducersTabButtonDown(app, event)
            % Get Transducer files for dropdown menu
            trFile_path = fullfile('..', 'Array_Positions');
            tr_pattern = {'std'};
            sparsity_pattern = {'spars'};
            file_ending = '.mat';

            app.ArrayElementsPositionsfilenameDropDown.Items = getDropdownEntries(app, trFile_path, tr_pattern, file_ending);

            sparse_entries = getDropdownEntries(app, trFile_path, sparsity_pattern, file_ending);
            app.SparsityfilenameDropDown.Items = [{''}, sparse_entries(:)'];

            % Get Propagation matrix files for dropdown menu
            A_path = fullfile('..', 'Lin_Prop_Matrices');
            A_pattern = {'A'};
            file_ending = '.mat';

            A_entries = getDropdownEntries(app, A_path, A_pattern, file_ending);
            app.PropagationMatrixAfilenameDropDown.Items = [{''}, A_entries(:)'];
        end

        % Value changed function: ElementGeometrySwitch
        function ElementGeometrySwitchValueChanged(app, event)
            value = app.ElementGeometrySwitch.Value;
            if strcmp(value, 'Rect')
                app.WidthmmEditField.Visible = true;
                app.WidthmmEditFieldLabel.Visible = true;
            else
                app.WidthmmEditField.Visible = false;
                app.WidthmmEditFieldLabel.Visible = false;
            end
        end

        % Button pushed function: AddTransducerButton
        function AddTransducerButtonPushed(app, event)
            app.ConfirmButtonPushed();
            
            new_item = num2str(length(app.TransducerDropDown.Items) + 1);
            app.TransducerDropDown.Items = [app.TransducerDropDown.Items(:)', {new_item}];

            app.TransducerDropDown.Value = {new_item};

            app.Transducer1Panel.Title = strcat("Transducer ", num2str(str2num(app.TransducerDropDown.Value)));
            app.Transducer1Panel.Visible = true;
        end

        % Button pushed function: UpdateButtonTransducer
        function UpdateButtonTransducerPushed(app, event)

            TransducerPreviewButtonPushed(app);

            %% Obtain or load propagation matrix
            if strcmp(app.PropagationMatrixAfilenameDropDown.Value, "")
                get_current_A = false;
            else
                get_current_A = app.PropagationMatrixAfilenameDropDown.Value(1:end-4);
            end
            app.A = obtain_linear_propagator(app.kgrid, app.medium, app.sensor, app.sensor_mask, app.input_args, ...
                app.t_mask_ps, app.karray_t, app.CenterFreqkHzEditField.Value * 1e3, get_current_A, app.use_greens_fctn, ...
                'active_ids', app.active_ids);

            disp('Transducer init successful')
        end

        % Value changed function: TransducerDropDown
        function TransducerDropDownValueChanged(app, event)
            value = str2num(app.TransducerDropDown.Value);

            app.trPosxEditField.Value = app.t_pos(1, value);
            app.trPosyEditField.Value = app.t_pos(2, value);
            app.trPoszEditField.Value = app.t_pos(3, value);
            
            app.alphaEditField.Value = app.t_rot(1, value);
            app.betaEditField.Value = app.t_rot(2, value);
            app.gammaEditField.Value = app.t_rot(3, value);

            app.TransducerLengthmmEditField.Value = app.tr_len(value);

            app.Transducer1Panel.Title = strcat("Transducer ", num2str(value));
            app.Transducer1Panel.Visible = true;
        end

        % Button pushed function: ConfirmButton
        function ConfirmButtonPushed(app, event)
            value = app.TransducerDropDown.Value;

            if ~isempty(value)
                value = str2num(value);

                % Add new column for new transducer
                if size(app.t_pos, 2) < value
                    app.t_pos = [app.t_pos, inf(3, 1)];
                    app.t_rot = [app.t_rot, inf(3, 1)];
                    app.tr_len = [app.tr_len, inf];
                end
    
                % Assign values to global transducer vars
                app.t_pos(:, value) = [app.trPosxEditField.Value, app.trPosyEditField.Value, app.trPoszEditField.Value];
                app.t_rot(:, value) = [app.alphaEditField.Value, app.betaEditField.Value, app.gammaEditField.Value];
                app.tr_len(value) = app.TransducerLengthmmEditField.Value;
            end
        end

        % Button pushed function: AddManualTargetButton
        function AddManualTargetButtonPushed(app, event)
            app.ConfirmManTargetButtonPushed();
            
            new_item = num2str(length(app.ManualTargetDropDown.Items) + 1);
            app.ManualTargetDropDown.Items = [app.ManualTargetDropDown.Items(:)', {new_item}];

            app.ManualTargetDropDown.Value = {new_item};
            app.TargetManPanel.Title = strcat("Target ", num2str(str2num(app.ManualTargetDropDown.Value)));

            app.TargetManPanel.Visible = true;
        end

        % Button pushed function: ConfirmManTargetButton
        function ConfirmManTargetButtonPushed(app, event)
            value = app.ManualTargetDropDown.Value;

            if ~isempty(value)
                value = str2num(value);
    
                % Add new column for new target
                if length(app.des_pressures) < value
                    app.point_pos_m = [app.point_pos_m, inf(3, 1)];
    
                    app.focus_radius = [app.focus_radius, inf];
                    app.des_pressures = [app.des_pressures, inf];
                    app.min_dist = [app.min_dist, inf];
    
                    app.force_pressures = [app.force_pressures, false];
                end
    
                % Assign values to global target vars
                app.point_pos_m(:, value) = [app.focusxEditField.Value, app.focusyEditField.Value, app.focuszEditField.Value];
    
                app.focus_radius(value) = app.FocusRadiusmmEditField.Value;
                app.des_pressures(value) = app.PressureAmplitudekPaEditField.Value;
                app.min_dist(value) = app.MinPointDistancemmEditField.Value;
    
                app.force_pressures(value) = app.todeclaredPressureSwitch.Value == "Force";
            end
        end

        % Button pushed function: AddRegionTargetButton
        function AddRegionTargetButtonPushed(app, event)
            app.ConfirmRegTargetButtonPushed();

            new_item = num2str(length(app.RegionTargetDropDown.Items) + 1);
            app.RegionTargetDropDown.Items = [app.RegionTargetDropDown.Items(:)', {new_item}];

            app.RegionTargetDropDown.Value = {new_item};
            app.TargetRegPanel.Title = strcat("Target ", num2str(str2num(app.RegionTargetDropDown.Value)));

            app.TargetRegPanel.Visible = true;
        end

        % Button pushed function: ConfirmRegTargetButton
        function ConfirmRegTargetButtonPushed(app, event)
            value = app.RegionTargetDropDown.Value;
            if ~isempty(value)
                value = str2num(value);

                % Add new column for new target
                if length(app.des_pressures_reg) < value
                    app.tar_reg_labels = [app.tar_reg_labels, ""];
    
                    app.des_pressures_reg = [app.des_pressures_reg, inf];
                    app.min_dist_reg = [app.min_dist_reg, inf];
    
                    app.force_pressures_reg = [app.force_pressures_reg, false];
                end
    
                % Assign values to global target vars
                app.tar_reg_labels(value) = app.BrainRegionDropDown.Value;
    
                app.des_pressures_reg(value) = app.PressureAmplitudekPaEditFieldReg.Value;
                app.min_dist_reg(value) = app.MinPointDistancemmEditFieldReg.Value;
    
                app.force_pressures_reg(value) = app.todeclaredPressureSwitchReg.Value == "Force";
            end
        end

        % Value changed function: ManualTargetDropDown
        function ManualTargetDropDownValueChanged(app, event)
            value = str2num(app.ManualTargetDropDown.Value);
            
            app.focusxEditField.Value = app.point_pos_m(1, value);
            app.focusyEditField.Value = app.point_pos_m(2, value);
            app.focuszEditField.Value = app.point_pos_m(3, value);

            app.FocusRadiusmmEditField.Value = app.focus_radius(value);
            app.PressureAmplitudekPaEditField.Value = app.des_pressures(value);
            app.MinPointDistancemmEditField.Value = app.min_dist(value);

            if app.force_pressures(value)
                app.todeclaredPressureSwitch.Value = "Force";
            else
                app.todeclaredPressureSwitch.Value = "Limit";
            end
                
            app.TargetManPanel.Title = strcat("Target ", num2str(value));
        end

        % Value changed function: RegionTargetDropDown
        function RegionTargetDropDownValueChanged(app, event)
            value = str2num(app.RegionTargetDropDown.Value);

            app.BrainRegionDropDown.Value = app.tar_reg_labels(value);

            app.PressureAmplitudekPaEditFieldReg.Value = app.des_pressures_reg(value);
            app.MinPointDistancemmEditFieldReg.Value = app.min_dist_reg(value);

            if app.force_pressures_reg(value)
                app.todeclaredPressureSwitchReg.Value = "Force";
            else
                app.todeclaredPressureSwitchReg.Value = "Limit";
            end
            
            app.TargetRegPanel.Title = strcat("Target ", num2str(value));
        end

        % Button down function: TargetingTab
        function TargetingTabButtonDown(app, event)
            % Get region labels for dropdown menu
            if ~isempty(app.segment_labels)
                app.BrainRegionDropDown.Items = app.segment_labels;
            end
        end

        % Button pushed function: UpdateTargetingButton
        function UpdateTargetingButtonPushed(app, event)

            ConfirmManTargetButtonPushed(app);
            ConfirmRegTargetButtonPushed(app);

            if app.n_dim == 2

                % Define targets
                if ~isempty(app.ManualTargetDropDown.Items)
                    plot_offset_rep = repmat(app.plot_offset, length(app.ManualTargetDropDown.Items), 1);
                    app.point_pos = round((plot_offset_rep' + app.point_pos_m) * app.dx_factor);

                    amp_in = app.des_pressures' * 1e3; % Pa
                
                    % Assign amplitude acc. to closest position
                    idx = sub2ind([app.kgrid.Nx, app.kgrid.Ny], app.point_pos(app.dims_2D(1), :), app.point_pos(app.dims_2D(2), :));
                    [~, order] = sort(idx);
                    amp_in = amp_in(order);
                else
                    app.point_pos = [];
                end
            
                app.b_mask = zeros(app.kgrid.Nx, app.kgrid.Ny, length(app.des_pressures) + length(app.tar_reg_labels));
                
                amp_vol = -1 * ones(numel(app.kgrid.k), 1);
            
                % Stimulate Disc pattern
                for i = 1:length(app.des_pressures)
                    disc = makeDisc(app.kgrid.Nx, app.kgrid.Ny, app.point_pos(app.dims_2D(1), i), app.point_pos(app.dims_2D(2), i), ...
                        round(app.focus_radius(i) * 1e-3 / app.kgrid.dx), false);
                    amp_vol(logical(disc)) = amp_in(i) * ones(sum(disc(:)), 1);
                    app.b_mask(:, :, i) = disc;
                end
            
                if ~isempty(app.t1w_filename) && ~isempty(app.RegionTargetDropDown.Items)
                    % Stimulate brain region
                    amp_in_reg = app.des_pressures_reg' * 1e3; % Pa
                    stim_regions = index2Dto3D(app.segment_ids, app.SliceDimDropDown.Value, app.slice_grid_2D);

                    for j = 1:length(app.tar_reg_labels)
                        reg_mask = stim_regions == app.tar_reg_labels(j);
                        amp_vol(logical(reg_mask)) = amp_in_reg(j) * ones(sum(reg_mask(:)), 1);
                        app.b_mask(:, :, length(app.des_pressures) + j) = reg_mask;
                    end
                end
            
                b_cross = amp_vol / max(amp_vol);
                b_cross(b_cross < 0.0) = 0.0;
                amp_in = amp_vol(amp_vol >= 0);
            
            else
                
                % Define targets
                if ~isempty(app.ManualTargetDropDown.Items)
                    plot_offset_rep = repmat(app.plot_offset, length(app.ManualTargetDropDown.Items), 1);
                    app.point_pos = round((plot_offset_rep' + app.point_pos_m) * app.dx_factor);

                    amp_in = app.des_pressures' * 1e3; % Pa
                
                    % Assign amplitude acc. to closest position
                    idx = sub2ind([app.kgrid.Nx, app.kgrid.Ny, app.kgrid.Nz], ...
                        app.point_pos(1, :), app.point_pos(2, :), app.point_pos(3, :));
                    [~, order] = sort(idx);
                    amp_in = amp_in(order);
                else
                    app.point_pos = [];
                end
            
                app.b_mask = zeros(app.kgrid.Nx, app.kgrid.Ny, app.kgrid.Nz, length(app.des_pressures) + length(app.tar_reg_labels));
            
                amp_vol = -1 * ones(numel(app.kgrid.k), 1);
            
                % Stimulate Disc pattern
                for i = 1:length(app.des_pressures)
                    ball = makeBall(app.kgrid.Nx, app.kgrid.Ny, app.kgrid.Nz, app.point_pos(1, i), app.point_pos(2, i), ...
                        app.point_pos(3, i), round(app.focus_radius(i) * 1e-3 / app.kgrid.dx), false);
                    amp_vol(logical(ball)) = amp_in(i) * ones(sum(ball(:)), 1);
                    app.b_mask(:, :, :, i) = ball;
                end
            
                if ~isempty(app.t1w_filename)
                    
                    % Stimulate brain region
                    amp_in_reg = app.des_pressures_reg' * 1e3; % Pa
                    stim_regions = app.segment_ids;
                    for j = 1:length(app.tar_reg_labels)
                        reg_mask = stim_regions == app.tar_reg_labels(j);
                        amp_vol(logical(reg_mask)) = amp_in_reg(j) * ones(sum(reg_mask(:)), 1);
                        app.b_mask(:, :, :, length(app.des_pressures) + j) = reg_mask;
                    end
                end
            
                b_cross = amp_vol / max(amp_vol);
                b_cross(b_cross < 0.0) = 0.0;
                amp_in = amp_vol(amp_vol >= 0);
            end
            app.b_mask = logical(app.b_mask);
            
            %% Create signal vector
            phase = zeros(length(amp_in), 1); % Zero phase for entire observation plane
            
            app.b_des = amp_in .* exp(1j*phase); % only observed elements
            app.full_bmask = sum(app.b_mask, app.n_dim + 1);
            app.full_bmask = logical(app.full_bmask);
            
%             b_max = max(abs(app.b_des));
            app.b_ip_des = app.MaxPressurekPaEditField.Value * 1e3 * ones(numel(app.kgrid.k), 1); % Entire plane max amp
            app.b_ip_des(app.full_bmask) = app.b_des; % Target amp

            skullMask = app.medium.sound_speed > min(app.medium.sound_speed(:));
            app.skull_ids = logical(reshape(skullMask, [], 1));
            app.b_ip_des(app.skull_ids & ~app.logical_dom_ids) ...
                = app.MaxPressurekPaSkullEditField.Value * 1e3; % Skull max amp

            %% Introduce constraints and prepare optimization
            app.vol_ids = reshape(logical(app.full_bmask), numel(app.full_bmask), 1); % Indices that correspond to the target volume(s)
            
            [app.init_ids, ~, b_mask_plot] = get_init_ids(app.kgrid, ...
                [app.min_dist, app.min_dist_reg] * 1e-3, app.b_mask, ...
                find([app.force_pressures, app.force_pressures_reg])); % Indices where pressure values forced

            % Update displayed slice based on first init_id
            if app.n_dim == 3
                [dispX, dispY, dispZ] = ind2sub(size(app.kgrid.k), find(app.init_ids));
                disp_id = [dispX(1), dispY(1), dispZ(1)];
                app.SliceIndexEditField.Value = round(disp_id(app.slice_dim) / app.dx_factor - app.plot_offset(app.slice_dim));
                SliceIndexEditFieldValueChanged(app);
            end
            
            % Create preview plot
            b_mask_plot = b_mask_plot + app.full_bmask;
            preplot_arg2 = app.preplot_arg + b_mask_plot / max(b_mask_plot(:));

            app.sv_obj = plot_results(app.kgrid, [], preplot_arg2, 'Target Preview', app.mask2el, app.t1w_filename, ...
                app.plot_offset, app.grid_size, app.dx_factor, false, [], 'slice', app.SliceIndexEditField.Value, ...
                'colorbar', false, 'cmap', hot(), 'slice_dim', app.SliceDimDropDown.Value);

            disp('Target init successful')
        end

        % Value changed function: LimitSkullPressureCheckBox
        function LimitSkullPressureCheckBoxValueChanged(app, event)
            value = app.LimitSkullPressureCheckBox.Value;
            
            if value
                app.MaxPressurekPaSkullEditField.Visible = true;
                app.MaxPressurekPaSkullEditFieldLabel.Visible = true;
            else
                app.MaxPressurekPaSkullEditField.Visible = false;
                app.MaxPressurekPaSkullEditFieldLabel.Visible = false;
            end
        end

        % Value changed function: 
        % LimitIntracranialOffTargetPressureCheckBox
        function LimitIntracranialOffTargetPressureCheckBoxValueChanged(app, event)
            value = app.LimitIntracranialOffTargetPressureCheckBox.Value;
            
            if value
                app.MaxPressurekPaEditField.Visible = true;
                app.MaxPressurekPaEditFieldLabel.Visible = true;
            else
                app.MaxPressurekPaEditField.Visible = false;
                app.MaxPressurekPaEditFieldLabel.Visible = false;
            end
        end

        % Button pushed function: ComputeInitialSolutionButton
        function ComputeInitialSolutionButtonPushed(app, event)
            get_initial_solution(app);
            evaluate_results(app, app.ip.p_init, "Initial Solution");
        end

        % Button pushed function: OptimizeButton
        function OptimizeButtonPushed(app, event)
            get_initial_solution(app);

            %% Get ids considered in inequality constraints
            dom_active = app.LimitIntracranialOffTargetPressureCheckBox.Value;
            skull_active = app.LimitSkullPressureCheckBox.Value;

            ineq_active = skull_active || dom_active;
            cons_ids = (skull_active & app.skull_ids) | (dom_active & app.logical_dom_ids);

            %% Get Optimization Options
            term.fun_tol = app.RelFunctionToleranceEditField.Value;
            term.constr_tol = app.ConstraintToleranceEditField.Value;
            term.iter_tol = 10;
            term.iter_lim = max([term.iter_tol, app.MaxNoofIterationsEditField.Value - term.iter_tol]);
            term.norm_val = app.LxNormEditField.Value;

            app.ip.beta = double(app.MinimizeArrrayAmplitudesCheckBox.Value);
            
            %% Optimize
            tic
            if app.OptimizeforPhasesandAmplitudesButton.Value
                app.ip.p = solvePhasesAmp(app.opt_mode.Value, app.A, app.b_ip_des, cons_ids, app.vol_ids, app.ip.p_init, ...
                    app.init_ids, app.ip.beta, ineq_active, term);
            else
                app.ip.p = solvePhasesOnly(app.opt_mode.Value, app.A, app.b_ip_des, cons_ids, app.vol_ids, app.ip.p_init, ...
                    app.init_ids, app.ip.beta, ineq_active, term, ...
                    app.mask2el, app.el_per_t, true, app.OptimizeforPhasesand1AmpButton.Value);
            end

            app.ip.t_solve = toc;
            disp("Time until solver converged: " + string(app.ip.t_solve / 60) + " min")

            app.ip.b = app.A * app.ip.p;
            app.ip.b = reshape(app.ip.b, size(app.kgrid.k));
            evaluate_results(app, app.ip.p, "Inverse Problem");
        end

        % Value changed function: MaxPressurekPaEditField
        function MaxPressurekPaEditFieldValueChanged(app, event)
            value = app.MaxPressurekPaEditField.Value;
            
        end

        % Button pushed function: DoGroundTruthSimulationButton
        function DoGroundTruthSimulationButtonPushed(app, event)
            f0 = app.CenterFreqkHzEditField.Value * 1e3;

            %% Compute GT solution
            [kgridP, mediumP, sensorP, sensor_maskP, dx_factorP, app.segment_ids, ~, tr_offset_karrP, ~] ...
                = init_general(app, app.dx_factor_init * app.GroundTruthResolutionFactorEditField.Value);

            [t_mask_psP, karray_tP, ~, ~] = transducer_geometry_init(app, kgridP, dx_factorP, tr_offset_karrP);

            if app.use_greens_fctn
                [amp_in, phase_in] = get_amp_phase_mask(app.kgrid, f0, app.ip.p, t_mask_psP, karray_tP);
                b_gt = acousticFieldPropagator(amp_in, phase_in, app.kgrid.dx, f0, app.medium.sound_speed);
            else
                b_gt = sim_exe(kgridP, mediumP, sensorP, f0, app.ip.p, t_mask_psP, sensor_maskP, true, app.input_args, ...
                    'karray_t', karray_tP);
                b_gt = reshape(b_gt, size(app.kgrid.k));
            end

            %% Plot
            save_results = app.SaveResultsCheckBox.Value;
            plot_thr = app.MaskPressurePlotkPaEditField.Value * 1e3;

            plot_title = 'Ground Truth';

            b_gt_lim = b_gt;
            if app.GroundTruthResolutionFactorEditField.Value == 1 % && ...
                   % app.LimittoIntracranialFieldCheckBox.Value || app.LimittoSkullCheckBox.Value

                excl_ids = (app.LimittoIntracranialFieldCheckBox.Value || app.LimittoSkullCheckBox.Value) & ...
                    ~(app.LimittoIntracranialFieldCheckBox.Value & app.logical_dom_ids) & ...
                    ~(app.LimittoSkullCheckBox.Value & app.skull_ids);

                b_gt_lim(excl_ids) = 0.0;

                % Evaluate Pressure Distribution
                evaluate_metrics(app, plot_title, b_gt_lim, app.vol_ids, app.logical_dom_ids, app.skull_ids, app.init_ids);
            end

            app.sv_obj = plot_results(kgridP, [], b_gt_lim, plot_title, app.mask2el, app.t1w_filename, app.plot_offset, app.grid_size, ...
                dx_factorP, save_results, app.current_datetime, 'slice', app.SliceIndexEditField.Value, ...
                'slice_dim', app.SliceDimDropDown.Value);

            masked_b = abs(b_gt_lim);
            masked_b(masked_b <= plot_thr) = 0.0;
            plot_results(kgridP, [], masked_b, strcat(plot_title, ' Mask'), app.mask2el, app.t1w_filename, ...
                app.plot_offset, app.grid_size, dx_factorP, save_results, app.current_datetime, 'slice', ...
                app.SliceIndexEditField.Value, 'fig_pos', {[], [1400 450 475 525], [1660 55 250 300]}, ...
                'slice_dim', app.SliceDimDropDown.Value);

            % Save results in mat file
            app.ip.b_gt = b_gt;
            save_results_mat(app);
        end

        % Button down function: OptimizeTab
        function OptimizeTabButtonDown(app, event)
            app.MaskPressurePlotkPaEditField.Value = app.MaxPressurekPaEditField.Value + 1;
        end

        % Callback function
        function UpdateSliceButtonPushed(app, event)
            evaluate_results(app, app.p_curr, app.plot_title_curr);
        end

        % Value changed function: DimSwitch
        function DimSwitchValueChanged(app, event)
            app.n_dim = 2 + (app.DimSwitch.Value == "3D");
            dim = dim2num(app.SliceDimDropDown.Value);

            visible_3D = app.n_dim == 3;
            labels2D_visible(app, visible_3D, dim);
        end

        % Callback function
        function RemoveTransducerButtonPushed(app, event)

        end

        % Button pushed function: RemoveTransducerButton
        function RemoveTransducerButtonPushed2(app, event)
            n_trs = length(app.TransducerDropDown.Items) - 1;

            if n_trs >= 0
                curr_item = str2num(app.TransducerDropDown.Value);
    
                % Delete item in global variables
                app.t_pos(:, curr_item) = [];
                app.t_rot(:, curr_item) = [];
                app.tr_len(curr_item) = [];
    
                if n_trs > 0
                    % Delete item in Dropdown and assign new indices from 1 to N
                    rem_items = num2str( (1:n_trs)' );
                    app.TransducerDropDown.Items = cellstr(rem_items)';
        
                    app.TransducerDropDown.Value = '1';
                    TransducerDropDownValueChanged(app);
                else
                    app.TransducerDropDown.Items = {};
                    app.Transducer1Panel.Visible = false;
                end
            end
        end

        % Button pushed function: RemoveManualTargetButton
        function RemoveManualTargetButtonPushed(app, event)
            n_man_tars = length(app.ManualTargetDropDown.Items) - 1;

            if n_man_tars >= 0
                curr_item = str2num(app.ManualTargetDropDown.Value);
    
                % Delete item in global variables
                app.point_pos_m(:, curr_item) = [];
    
                app.focus_radius(curr_item) = [];
                app.des_pressures(curr_item) = [];
                app.min_dist(curr_item) = [];
    
                app.force_pressures(curr_item) = [];
    
                if n_man_tars > 0
                    % Delete item in Dropdown and assign new indices from 1 to N
                    rem_items = num2str( (1:n_man_tars)' );
                    app.ManualTargetDropDown.Items = cellstr(rem_items)';
        
                    app.ManualTargetDropDown.Value = '1';
                    ManualTargetDropDownValueChanged(app);
                else
                    app.ManualTargetDropDown.Items = {};
                    app.TargetManPanel.Visible = false;
                end
            end
        end

        % Button pushed function: RemoveRegionTargetButton
        function RemoveRegionTargetButtonPushed(app, event)
            n_reg_tars = length(app.RegionTargetDropDown.Items) - 1;

            if n_reg_tars >= 0
                curr_item = str2num(app.RegionTargetDropDown.Value);
    
                % Delete item in global variables
                app.tar_reg_labels(curr_item) = [];
                app.des_pressures_reg(curr_item) = [];
                app.min_dist_reg(curr_item) = [];
    
                app.force_pressures_reg(curr_item) = [];
    
                if n_reg_tars > 0
                    % Delete item in Dropdown and assign new indices from 1 to N
                    rem_items = num2str( (1:n_reg_tars)' );
                    app.RegionTargetDropDown.Items = cellstr(rem_items)';
        
                    app.RegionTargetDropDown.Value = '1';
                    RegionTargetDropDownValueChanged(app);
                else
                    app.RegionTargetDropDown.Items = {};
                    app.TargetRegPanel.Visible = false;
                end
            end
        end

        % Clicked callback: LoadResultsfromFileDropDown
        function LoadResultsfromFileDropDownClicked(app, event)
            item = event.InteractionInformation.Item;
            
            % Get results file names for dropdown menu
            res_path = fullfile('..', 'Results');
            res_pattern = {'results'};
            file_ending = '.mat';

            result_files = getDropdownEntries(app, res_path, res_pattern, file_ending);
            app.LoadResultsfromFileDropDown.Items = [{''}, result_files(:)'];
        end

        % Value changed function: LoadResultsfromFileDropDown
        function LoadResultsfromFileDropDownValueChanged(app, event)
            value = app.LoadResultsfromFileDropDown.Value;
            if ~strcmp(value, '')
                filename = fullfile("..", "Results", value);
                data = load(filename);

                % Std param
                app.DimSwitch.Value = data.dim_switch; DimSwitchValueChanged(app);
                app.CenterFreqkHzEditField.Value = data.f0;
                app.SimSpatialResolutionmmEditField.Value = data.dx; CenterFreqkHzEditFieldValueChanged(app);

                app.SliceDimDropDown.Value = data.slice_direction;
                app.SliceIndexEditField.Value = data.slice_index; SliceIndexEditFieldValueChanged(app);

                app.MediumSwitch.Value = data.medium_type; MediumSwitchValueChanged(app);
                app.GreensFunctionbasedCheckBox.Value = data.green_function;
                app.xHomEditField.Value = data.hom_coord(1);
                app.HomyEditField.Value = data.hom_coord(2);
                app.HomzEditField.Value = data.hom_coord(3);
                app.T1wfilenameDropDown.Value = data.t1w_file;
                app.CTfilenameDropDown.Value = data.ct_file;
                app.scanxEditField.Value = data.het_coord(1);
                app.scanyEditField.Value = data.het_coord(2);
                app.scanzEditField.Value = data.het_coord(3);
                app.ScanSpatialResolutionmmEditField.Value = data.dx_t1w;
                app.segment_labels = data.seg_labels;

                % Transducer param
                TransducersTabButtonDown(app);
                app.ElementGeometrySwitch.Value = data.array.geom; ElementGeometrySwitchValueChanged(app);
                app.LengthmmEditField.Value = data.array.len;
                app.WidthmmEditField.Value = data.array.wid;
                app.t_pos = data.tr.pos;
                app.t_rot = data.tr.rot;
                app.tr_len = data.tr.length; 
                
                app.Transducer1Panel.Visible = ~isempty(app.t_pos);
                if ~isempty(app.t_pos)
                    app.TransducerDropDown.Items = cellstr( num2str((1:size(app.t_pos, 2))') );
                    app.TransducerDropDown.Value = '1';

                    app.trPosxEditField.Value = app.t_pos(1, 1);
                    app.trPosyEditField.Value = app.t_pos(2, 1);
                    app.trPoszEditField.Value = app.t_pos(3, 1);

                    app.alphaEditField.Value = app.t_rot(1, 1);
                    app.betaEditField.Value = app.t_rot(2, 1);
                    app.gammaEditField.Value = app.t_rot(3, 1);

                    app.TransducerLengthmmEditField.Value = app.tr_len(1);
                end
                
                if isfield(data, 'dropdowns') && isfield(data.dropdowns, 'array_pos')
                    app.ArrayElementsPositionsfilenameDropDown.Value = data.dropdowns.array_pos;
                end
                if isfield(data, 'dropdowns') && isfield(data.dropdowns, 'sparsity')
                    app.SparsityfilenameDropDown.Value = data.dropdowns.sparsity;
                end
                
                app.PropagationMatrixAfilenameDropDown.Value = data.matrix_name;

                % Target param
                TargetingTabButtonDown(app)
                app.point_pos_m = data.man.points;
                app.focus_radius = data.man.radius;
                app.des_pressures = data.man.pressure;
                app.min_dist = data.man.min_dist;
                app.force_pressures = data.man.force;
                
                app.TargetManPanel.Visible = ~isempty(app.des_pressures);
                if ~isempty(app.des_pressures)
                    app.ManualTargetDropDown.Items = cellstr( num2str((1:length(app.des_pressures))') );
                    app.ManualTargetDropDown.Value = '1';

                    app.focusxEditField.Value = app.point_pos_m(1, 1);
                    app.focusyEditField.Value = app.point_pos_m(2, 1);
                    app.focuszEditField.Value = app.point_pos_m(3, 1);

                    app.FocusRadiusmmEditField.Value = app.focus_radius(1);
                    app.PressureAmplitudekPaEditField.Value = app.des_pressures(1);
                    app.MinPointDistancemmEditField.Value = app.min_dist(1);

                    if app.force_pressures(1)
                        switch_val_man = 'Force';
                    else
                        switch_val_man = 'Limit';
                    end
                    app.todeclaredPressureSwitch.Value = switch_val_man;
                end
                
                app.tar_reg_labels = data.reg.labels;
                app.des_pressures_reg = data.reg.pressure;
                app.min_dist_reg = data.reg.min_dist;
                app.force_pressures_reg = data.reg.force;

                app.TargetRegPanel.Visible = ~isempty(app.des_pressures_reg);
                if ~isempty(app.des_pressures_reg)
                    app.RegionTargetDropDown.Items = cellstr( num2str((1:length(app.des_pressures_reg))') );
                    app.RegionTargetDropDown.Value = '1';

                    app.BrainRegionDropDown.Value = convertStringsToChars(app.tar_reg_labels(1));

                    app.PressureAmplitudekPaEditFieldReg.Value = app.des_pressures_reg(1);
                    app.MinPointDistancemmEditFieldReg.Value = app.min_dist_reg(1);

                    if app.force_pressures_reg(1)
                        switch_val_reg = 'Force';
                    else
                        switch_val_reg = 'Limit';
                    end
                    app.todeclaredPressureSwitchReg.Value = switch_val_reg;
                end
                
                app.LimitIntracranialOffTargetPressureCheckBox.Value = data.ineq_flag;
                app.MaxPressurekPaEditField.Value = data.ineq_val;
                app.LimitSkullPressureCheckBox.Value = data.skull_flag;
                app.MaxPressurekPaSkullEditField.Value = data.skull_val;

                % Optimization param
                OptimizeTabButtonDown(app);
                app.OptimizeforPhasesandAmplitudesButton.Value = data.opt.phaseAmp;
                app.OptimizeforPhasesand1AmpperTransducerButton.Value = data.opt.phaseOnly;
                app.LimittoIntracranialFieldCheckBox.Value = data.opt.intracranial;
                app.LimittoSkullCheckBox.Value = data.opt.skull;
                app.GroundTruthResolutionFactorEditField.Value = data.gt_dx_factor;

                app.DisplayAdvancedOptimizationOptionsCheckBox.Value = data.term.advanced_opt;
                app.MinimizeArrrayAmplitudesCheckBox.Value = data.term.min_arr_amp;
                app.RelFunctionToleranceEditField.Value = data.term.fun_tol;
                app.ConstraintToleranceEditField.Value = data.term.constr_tol;
                app.MaxNoofIterationsEditField.Value = data.term.max_iter;
                app.AlgorithmDropDown.Value = data.term.algorithm;

                app.ip = data.ip_sol;
            end
        end

        % Value changed function: SliceIndexEditField
        function SliceIndexEditFieldValueChanged(app, event)
            value = app.SliceIndexEditField.Value;
            app.Slice30Label.Text = strcat("Slice ", num2str(value));
        end

        % Value changed function: 
        % DisplayAdvancedOptimizationOptionsCheckBox
        function DisplayAdvancedOptimizationOptionsCheckBoxValueChanged(app, event)
            value = app.DisplayAdvancedOptimizationOptionsCheckBox.Value;

            if value
                app.AdvancedOptimizationOptionsPanel.Visible = true;
            else
                app.AdvancedOptimizationOptionsPanel.Visible = false;
            end
        end

        % Button pushed function: CloseallPlotsButton
        function CloseallPlotsButtonPushed(app, event)
            close all;
        end

        % Value changed function: MinimizeArrrayAmplitudesCheckBox
        function MinimizeArrrayAmplitudesCheckBoxValueChanged(app, event)
            value = app.MinimizeArrrayAmplitudesCheckBox.Value;
            app.ip.beta = double(value);

            % Update default optimization options
            if value
                app.MaxNoofIterationsEditField.Value = 1000;
            else
                app.MaxNoofIterationsEditField.Value = 200;
            end
        end

        % Value changed function: SliceDimDropDown
        function SliceDimDropDownValueChanged(app, event)
            value = app.SliceDimDropDown.Value;
            
            persistent prev_dim
            if isempty(prev_dim)
                prev_dim = 2;
            end
            
            labels2D_visible(app, true, prev_dim); % Make dimensions of previous dim visible
            DimSwitchValueChanged(app);

            prev_dim = dim2num(value);

            app.SliceDirectionDropDown_2.Value = value;
            UpdateSliceButtonPushed2(app);
        end

        % Button pushed function: UpdateSliceButton
        function UpdateSliceButtonPushed2(app, event)
            if ~isempty(app.sv_obj) && app.n_dim == 3
                if isvalid(app.sv_obj)
                    currentSliceIndex = app.sv_obj.SliceNumber;
                    app.SliceIndexEditField.Value = round(currentSliceIndex / app.dx_factor - app.plot_offset(app.slice_dim));
                    SliceIndexEditFieldValueChanged(app);
                end
            end
        end

        % Value changed function: SliceDirectionDropDown_2
        function SliceDirectionDropDown_2ValueChanged(app, event)
            value = app.SliceDirectionDropDown_2.Value;
            app.SliceDimDropDown.Value = value;
            SliceDimDropDownValueChanged(app);

            app.slice_dim = dim2num(app.SliceDimDropDown.Value);
            app.dims_2D = exclude_dim(app.slice_dim);
            UpdateSliceButtonPushed2(app);
        end

        % Button pushed function: TransducerPreviewButton
        function TransducerPreviewButtonPushed(app, event)
            ConfirmButtonPushed(app);
            [app.t_mask_ps, app.karray_t, app.el_per_t, app.active_ids] = transducer_geometry_init(app, app.kgrid, ...
                app.dx_factor, app.tr_offset_karr);

            %% Create preview plot
            app.preplot_arg = zeros(size(app.kgrid.k));
            app.preplot_arg(logical(app.t_mask_ps)) = 1.0;
            
            if ~isscalar(app.medium.sound_speed)
                skull_arg = app.medium.sound_speed / max(app.medium.sound_speed(:));
                skull_arg = skull_arg - min(skull_arg(:));
                app.preplot_arg = app.preplot_arg + skull_arg;
            end

            app.sv_obj = plot_results(app.kgrid, [], app.preplot_arg, 'Transducer Preview', app.mask2el, app.t1w_filename, ...
                app.plot_offset, app.grid_size, app.dx_factor, false, [], 'slice', app.SliceIndexEditField.Value, ...
                'colorbar', false, 'cmap', hot(), 'slice_dim', app.SliceDimDropDown.Value);
        end

        % Drop down opening function: RegionTargetDropDown
        function RegionTargetDropDownOpening(app, event)
            ConfirmRegTargetButtonPushed(app);
        end

        % Drop down opening function: ManualTargetDropDown
        function ManualTargetDropDownOpening(app, event)
            ConfirmManTargetButtonPushed(app);
        end

        % Drop down opening function: TransducerDropDown
        function TransducerDropDownOpening(app, event)
            ConfirmButtonPushed(app);
        end

        % Clicked callback: RegionTargetDropDown
        function RegionTargetDropDownClicked(app, event)
%             item = event.InteractionInformation.Item;
        end

        % Button pushed function: PlotResultsButton
        function PlotResultsButtonPushed(app, event)
            if isfield(app.ip, 'p')
                evaluate_results(app, app.ip.p, "Inverse Problem");
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [25 445 727 581];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 32 727 550];

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
            app.PMLsizeEditField.RoundFractionalValues = 'on';
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
            app.GeneralPanel.Position = [31 221 298 291];

            % Create RealdxLabel
            app.RealdxLabel = uilabel(app.GeneralPanel);
            app.RealdxLabel.Position = [47 57 238 22];
            app.RealdxLabel.Text = '-> Real: 1 mm';

            % Create Label
            app.Label = uilabel(app.GeneralPanel);
            app.Label.HorizontalAlignment = 'center';
            app.Label.Position = [40 239 26 22];
            app.Label.Text = 'Dim';

            % Create DimSwitch
            app.DimSwitch = uiswitch(app.GeneralPanel, 'slider');
            app.DimSwitch.Items = {'2D', '3D'};
            app.DimSwitch.ValueChangedFcn = createCallbackFcn(app, @DimSwitchValueChanged, true);
            app.DimSwitch.Position = [30 209 45 20];
            app.DimSwitch.Value = '2D';

            % Create CenterFreqkHzEditFieldLabel
            app.CenterFreqkHzEditFieldLabel = uilabel(app.GeneralPanel);
            app.CenterFreqkHzEditFieldLabel.HorizontalAlignment = 'right';
            app.CenterFreqkHzEditFieldLabel.Position = [126 209 101 22];
            app.CenterFreqkHzEditFieldLabel.Text = 'Center Freq (kHz)';

            % Create CenterFreqkHzEditField
            app.CenterFreqkHzEditField = uieditfield(app.GeneralPanel, 'numeric');
            app.CenterFreqkHzEditField.ValueChangedFcn = createCallbackFcn(app, @CenterFreqkHzEditFieldValueChanged, true);
            app.CenterFreqkHzEditField.Position = [239 209 46 22];
            app.CenterFreqkHzEditField.Value = 500;

            % Create SimSpatialResolutionmmEditFieldLabel
            app.SimSpatialResolutionmmEditFieldLabel = uilabel(app.GeneralPanel);
            app.SimSpatialResolutionmmEditFieldLabel.HorizontalAlignment = 'right';
            app.SimSpatialResolutionmmEditFieldLabel.Position = [41 78 157 22];
            app.SimSpatialResolutionmmEditFieldLabel.Text = 'Sim Spatial Resolution (mm)';

            % Create SimSpatialResolutionmmEditField
            app.SimSpatialResolutionmmEditField = uieditfield(app.GeneralPanel, 'numeric');
            app.SimSpatialResolutionmmEditField.ValueChangedFcn = createCallbackFcn(app, @SimSpatialResolutionmmEditFieldValueChanged, true);
            app.SimSpatialResolutionmmEditField.Position = [213 78 45 22];

            % Create SliceDimDropDownLabel
            app.SliceDimDropDownLabel = uilabel(app.GeneralPanel);
            app.SliceDimDropDownLabel.HorizontalAlignment = 'right';
            app.SliceDimDropDownLabel.Position = [75 159 56 22];
            app.SliceDimDropDownLabel.Text = 'Slice Dim';

            % Create SliceDimDropDown
            app.SliceDimDropDown = uidropdown(app.GeneralPanel);
            app.SliceDimDropDown.Items = {'X', 'Y', 'Z'};
            app.SliceDimDropDown.ValueChangedFcn = createCallbackFcn(app, @SliceDimDropDownValueChanged, true);
            app.SliceDimDropDown.Position = [146 159 45 22];
            app.SliceDimDropDown.Value = 'Y';

            % Create SliceIndexLabel
            app.SliceIndexLabel = uilabel(app.GeneralPanel);
            app.SliceIndexLabel.HorizontalAlignment = 'right';
            app.SliceIndexLabel.Position = [64 125 64 22];
            app.SliceIndexLabel.Text = 'Slice Index';

            % Create SliceIndexEditField
            app.SliceIndexEditField = uieditfield(app.GeneralPanel, 'numeric');
            app.SliceIndexEditField.RoundFractionalValues = 'on';
            app.SliceIndexEditField.ValueChangedFcn = createCallbackFcn(app, @SliceIndexEditFieldValueChanged, true);
            app.SliceIndexEditField.Position = [144 125 48 22];
            app.SliceIndexEditField.Value = 30;

            % Create SaveResultsCheckBox
            app.SaveResultsCheckBox = uicheckbox(app.GeneralPanel);
            app.SaveResultsCheckBox.Text = 'Save Results';
            app.SaveResultsCheckBox.Position = [14 34 93 22];

            % Create LoadResultsfromFileDropDownLabel
            app.LoadResultsfromFileDropDownLabel = uilabel(app.GeneralPanel);
            app.LoadResultsfromFileDropDownLabel.HorizontalAlignment = 'right';
            app.LoadResultsfromFileDropDownLabel.Position = [9 8 125 22];
            app.LoadResultsfromFileDropDownLabel.Text = 'Load Results from File';

            % Create LoadResultsfromFileDropDown
            app.LoadResultsfromFileDropDown = uidropdown(app.GeneralPanel);
            app.LoadResultsfromFileDropDown.Items = {''};
            app.LoadResultsfromFileDropDown.ValueChangedFcn = createCallbackFcn(app, @LoadResultsfromFileDropDownValueChanged, true);
            app.LoadResultsfromFileDropDown.ClickedFcn = createCallbackFcn(app, @LoadResultsfromFileDropDownClicked, true);
            app.LoadResultsfromFileDropDown.Position = [149 8 143 22];
            app.LoadResultsfromFileDropDown.Value = '';

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
            app.scanxEditField.RoundFractionalValues = 'on';
            app.scanxEditField.Position = [62 71 42 22];
            app.scanxEditField.Value = -96;

            % Create zEditFieldLabel
            app.zEditFieldLabel = uilabel(app.hetPanel);
            app.zEditFieldLabel.HorizontalAlignment = 'right';
            app.zEditFieldLabel.Position = [225 71 25 22];
            app.zEditFieldLabel.Text = 'z';

            % Create scanzEditField
            app.scanzEditField = uieditfield(app.hetPanel, 'numeric');
            app.scanzEditField.RoundFractionalValues = 'on';
            app.scanzEditField.Position = [265 71 42 22];
            app.scanzEditField.Value = -126;

            % Create yEditFieldLabel
            app.yEditFieldLabel = uilabel(app.hetPanel);
            app.yEditFieldLabel.HorizontalAlignment = 'right';
            app.yEditFieldLabel.Position = [119 71 25 22];
            app.yEditFieldLabel.Text = 'y';

            % Create scanyEditField
            app.scanyEditField = uieditfield(app.hetPanel, 'numeric');
            app.scanyEditField.RoundFractionalValues = 'on';
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
            app.homPanel.Position = [1 13 308 219];

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
            app.xHomEditField.RoundFractionalValues = 'on';
            app.xHomEditField.Position = [69 62 37 22];
            app.xHomEditField.Value = 192;

            % Create yEditField_2Label
            app.yEditField_2Label = uilabel(app.homPanel);
            app.yEditField_2Label.HorizontalAlignment = 'right';
            app.yEditField_2Label.Visible = 'off';
            app.yEditField_2Label.Position = [114 62 25 22];
            app.yEditField_2Label.Text = 'y';

            % Create HomyEditField
            app.HomyEditField = uieditfield(app.homPanel, 'numeric');
            app.HomyEditField.Limits = [0 Inf];
            app.HomyEditField.RoundFractionalValues = 'on';
            app.HomyEditField.Visible = 'off';
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
            app.HomzEditField.RoundFractionalValues = 'on';
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
            app.TransducersTab.ButtonDownFcn = createCallbackFcn(app, @TransducersTabButtonDown, true);

            % Create TransducerDropDownLabel
            app.TransducerDropDownLabel = uilabel(app.TransducersTab);
            app.TransducerDropDownLabel.HorizontalAlignment = 'right';
            app.TransducerDropDownLabel.Position = [48 311 75 22];
            app.TransducerDropDownLabel.Text = 'Transducer #';

            % Create TransducerDropDown
            app.TransducerDropDown = uidropdown(app.TransducersTab);
            app.TransducerDropDown.Items = {};
            app.TransducerDropDown.DropDownOpeningFcn = createCallbackFcn(app, @TransducerDropDownOpening, true);
            app.TransducerDropDown.ValueChangedFcn = createCallbackFcn(app, @TransducerDropDownValueChanged, true);
            app.TransducerDropDown.Position = [138 311 100 22];
            app.TransducerDropDown.Value = {};

            % Create Transducer1Panel
            app.Transducer1Panel = uipanel(app.TransducersTab);
            app.Transducer1Panel.Title = 'Transducer 1';
            app.Transducer1Panel.Visible = 'off';
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
            app.trPosxEditField.RoundFractionalValues = 'on';
            app.trPosxEditField.Position = [90 149 37 22];
            app.trPosxEditField.Value = -59;

            % Create yEditField_3Label
            app.yEditField_3Label = uilabel(app.Transducer1Panel);
            app.yEditField_3Label.HorizontalAlignment = 'right';
            app.yEditField_3Label.Visible = 'off';
            app.yEditField_3Label.Position = [135 149 25 22];
            app.yEditField_3Label.Text = 'y';

            % Create trPosyEditField
            app.trPosyEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.trPosyEditField.RoundFractionalValues = 'on';
            app.trPosyEditField.Visible = 'off';
            app.trPosyEditField.Position = [170 149 37 22];
            app.trPosyEditField.Value = 30;

            % Create zEditField_3Label
            app.zEditField_3Label = uilabel(app.Transducer1Panel);
            app.zEditField_3Label.HorizontalAlignment = 'right';
            app.zEditField_3Label.Position = [219 149 25 22];
            app.zEditField_3Label.Text = 'z';

            % Create trPoszEditField
            app.trPoszEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.trPoszEditField.RoundFractionalValues = 'on';
            app.trPoszEditField.Position = [254 149 37 22];
            app.trPoszEditField.Value = 60;

            % Create RotationdegxyzintrinsicLabel
            app.RotationdegxyzintrinsicLabel = uilabel(app.Transducer1Panel);
            app.RotationdegxyzintrinsicLabel.Position = [100 105 183 22];
            app.RotationdegxyzintrinsicLabel.Text = 'Rotation (deg) - x-y''-z'''' -> intrinsic';

            % Create alphaEditFieldLabel
            app.alphaEditFieldLabel = uilabel(app.Transducer1Panel);
            app.alphaEditFieldLabel.HorizontalAlignment = 'right';
            app.alphaEditFieldLabel.Visible = 'off';
            app.alphaEditFieldLabel.Position = [42 75 34 22];
            app.alphaEditFieldLabel.Text = 'alpha';

            % Create alphaEditField
            app.alphaEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.alphaEditField.Visible = 'off';
            app.alphaEditField.Position = [86 75 37 22];

            % Create betaEditFieldLabel
            app.betaEditFieldLabel = uilabel(app.Transducer1Panel);
            app.betaEditFieldLabel.HorizontalAlignment = 'right';
            app.betaEditFieldLabel.Position = [137 75 28 22];
            app.betaEditFieldLabel.Text = 'beta';

            % Create betaEditField
            app.betaEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.betaEditField.Position = [175 75 37 22];
            app.betaEditField.Value = 45;

            % Create gammaEditFieldLabel
            app.gammaEditFieldLabel = uilabel(app.Transducer1Panel);
            app.gammaEditFieldLabel.HorizontalAlignment = 'right';
            app.gammaEditFieldLabel.Visible = 'off';
            app.gammaEditFieldLabel.Position = [233 75 45 22];
            app.gammaEditFieldLabel.Text = 'gamma';

            % Create gammaEditField
            app.gammaEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.gammaEditField.Visible = 'off';
            app.gammaEditField.Position = [288 75 37 22];
            app.gammaEditField.Value = 180;

            % Create TransducerLengthmmEditFieldLabel
            app.TransducerLengthmmEditFieldLabel = uilabel(app.Transducer1Panel);
            app.TransducerLengthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.TransducerLengthmmEditFieldLabel.Position = [36 23 137 22];
            app.TransducerLengthmmEditFieldLabel.Text = 'Transducer Length (mm)';

            % Create TransducerLengthmmEditField
            app.TransducerLengthmmEditField = uieditfield(app.Transducer1Panel, 'numeric');
            app.TransducerLengthmmEditField.Position = [183 23 55 22];
            app.TransducerLengthmmEditField.Value = 70;

            % Create ConfirmButton
            app.ConfirmButton = uibutton(app.Transducer1Panel, 'push');
            app.ConfirmButton.ButtonPushedFcn = createCallbackFcn(app, @ConfirmButtonPushed, true);
            app.ConfirmButton.Visible = 'off';
            app.ConfirmButton.Position = [275 11 100 23];
            app.ConfirmButton.Text = 'Confirm';

            % Create AddTransducerButton
            app.AddTransducerButton = uibutton(app.TransducersTab, 'push');
            app.AddTransducerButton.ButtonPushedFcn = createCallbackFcn(app, @AddTransducerButtonPushed, true);
            app.AddTransducerButton.Position = [129 251 122 23];
            app.AddTransducerButton.Text = 'Add Transducer';

            % Create UpdateButtonTransducer
            app.UpdateButtonTransducer = uibutton(app.TransducersTab, 'push');
            app.UpdateButtonTransducer.ButtonPushedFcn = createCallbackFcn(app, @UpdateButtonTransducerPushed, true);
            app.UpdateButtonTransducer.Position = [576 63 100 23];
            app.UpdateButtonTransducer.Text = 'Update';

            % Create PropagationMatrixAfilenameDropDownLabel
            app.PropagationMatrixAfilenameDropDownLabel = uilabel(app.TransducersTab);
            app.PropagationMatrixAfilenameDropDownLabel.HorizontalAlignment = 'right';
            app.PropagationMatrixAfilenameDropDownLabel.Position = [46 62 174 22];
            app.PropagationMatrixAfilenameDropDownLabel.Text = 'Propagation Matrix (A) filename';

            % Create PropagationMatrixAfilenameDropDown
            app.PropagationMatrixAfilenameDropDown = uidropdown(app.TransducersTab);
            app.PropagationMatrixAfilenameDropDown.Items = {''};
            app.PropagationMatrixAfilenameDropDown.Position = [235 62 174 22];
            app.PropagationMatrixAfilenameDropDown.Value = '';

            % Create GeneralPanel_2
            app.GeneralPanel_2 = uipanel(app.TransducersTab);
            app.GeneralPanel_2.Title = 'General';
            app.GeneralPanel_2.Visible = 'off';
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
            app.ElementGeometrySwitch.ValueChangedFcn = createCallbackFcn(app, @ElementGeometrySwitchValueChanged, true);
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

            % Create RemoveTransducerButton
            app.RemoveTransducerButton = uibutton(app.TransducersTab, 'push');
            app.RemoveTransducerButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveTransducerButtonPushed2, true);
            app.RemoveTransducerButton.Position = [129 213 122 23];
            app.RemoveTransducerButton.Text = 'Remove Transducer';

            % Create TransducerPreviewButton
            app.TransducerPreviewButton = uibutton(app.TransducersTab, 'push');
            app.TransducerPreviewButton.ButtonPushedFcn = createCallbackFcn(app, @TransducerPreviewButtonPushed, true);
            app.TransducerPreviewButton.Position = [555 123 121 23];
            app.TransducerPreviewButton.Text = 'Transducer Preview';

            % Create TargetingTab
            app.TargetingTab = uitab(app.TabGroup);
            app.TargetingTab.Title = 'Targeting';
            app.TargetingTab.ButtonDownFcn = createCallbackFcn(app, @TargetingTabButtonDown, true);

            % Create TargetManPanel
            app.TargetManPanel = uipanel(app.TargetingTab);
            app.TargetManPanel.Title = 'Target 1';
            app.TargetManPanel.Visible = 'off';
            app.TargetManPanel.Position = [291 278 383 229];

            % Create FocusPositionLabel
            app.FocusPositionLabel = uilabel(app.TargetManPanel);
            app.FocusPositionLabel.Position = [149 179 84 22];
            app.FocusPositionLabel.Text = 'Focus Position';

            % Create xEditField_4Label
            app.xEditField_4Label = uilabel(app.TargetManPanel);
            app.xEditField_4Label.HorizontalAlignment = 'right';
            app.xEditField_4Label.Position = [55 149 25 22];
            app.xEditField_4Label.Text = 'x';

            % Create focusxEditField
            app.focusxEditField = uieditfield(app.TargetManPanel, 'numeric');
            app.focusxEditField.RoundFractionalValues = 'on';
            app.focusxEditField.Position = [90 149 37 22];
            app.focusxEditField.Value = -18;

            % Create yEditField_4Label
            app.yEditField_4Label = uilabel(app.TargetManPanel);
            app.yEditField_4Label.HorizontalAlignment = 'right';
            app.yEditField_4Label.Visible = 'off';
            app.yEditField_4Label.Position = [135 149 25 22];
            app.yEditField_4Label.Text = 'y';

            % Create focusyEditField
            app.focusyEditField = uieditfield(app.TargetManPanel, 'numeric');
            app.focusyEditField.RoundFractionalValues = 'on';
            app.focusyEditField.Visible = 'off';
            app.focusyEditField.Position = [170 149 37 22];
            app.focusyEditField.Value = 30;

            % Create zEditField_4Label
            app.zEditField_4Label = uilabel(app.TargetManPanel);
            app.zEditField_4Label.HorizontalAlignment = 'right';
            app.zEditField_4Label.Position = [219 149 25 22];
            app.zEditField_4Label.Text = 'z';

            % Create focuszEditField
            app.focuszEditField = uieditfield(app.TargetManPanel, 'numeric');
            app.focuszEditField.RoundFractionalValues = 'on';
            app.focuszEditField.Position = [254 149 37 22];
            app.focuszEditField.Value = -27;

            % Create FocusRadiusmmEditFieldLabel
            app.FocusRadiusmmEditFieldLabel = uilabel(app.TargetManPanel);
            app.FocusRadiusmmEditFieldLabel.HorizontalAlignment = 'right';
            app.FocusRadiusmmEditFieldLabel.Position = [43 80 110 22];
            app.FocusRadiusmmEditFieldLabel.Text = 'Focus Radius (mm)';

            % Create FocusRadiusmmEditField
            app.FocusRadiusmmEditField = uieditfield(app.TargetManPanel, 'numeric');
            app.FocusRadiusmmEditField.Position = [164 80 40 22];
            app.FocusRadiusmmEditField.Value = 7;

            % Create PressureAmplitudekPaEditFieldLabel
            app.PressureAmplitudekPaEditFieldLabel = uilabel(app.TargetManPanel);
            app.PressureAmplitudekPaEditFieldLabel.HorizontalAlignment = 'right';
            app.PressureAmplitudekPaEditFieldLabel.Position = [12 47 141 22];
            app.PressureAmplitudekPaEditFieldLabel.Text = 'Pressure Amplitude (kPa)';

            % Create PressureAmplitudekPaEditField
            app.PressureAmplitudekPaEditField = uieditfield(app.TargetManPanel, 'numeric');
            app.PressureAmplitudekPaEditField.Position = [164 47 40 22];
            app.PressureAmplitudekPaEditField.Value = 300;

            % Create MinPointDistancemmEditFieldLabel
            app.MinPointDistancemmEditFieldLabel = uilabel(app.TargetManPanel);
            app.MinPointDistancemmEditFieldLabel.HorizontalAlignment = 'right';
            app.MinPointDistancemmEditFieldLabel.Position = [7 17 140 22];
            app.MinPointDistancemmEditFieldLabel.Text = 'Min. Point Distance (mm)';

            % Create MinPointDistancemmEditField
            app.MinPointDistancemmEditField = uieditfield(app.TargetManPanel, 'numeric');
            app.MinPointDistancemmEditField.Position = [163 17 43 22];
            app.MinPointDistancemmEditField.Value = 5;

            % Create todeclaredPressureSwitch_2Label
            app.todeclaredPressureSwitch_2Label = uilabel(app.TargetManPanel);
            app.todeclaredPressureSwitch_2Label.HorizontalAlignment = 'center';
            app.todeclaredPressureSwitch_2Label.Position = [249 69 116 22];
            app.todeclaredPressureSwitch_2Label.Text = 'to declared Pressure';

            % Create todeclaredPressureSwitch
            app.todeclaredPressureSwitch = uiswitch(app.TargetManPanel, 'slider');
            app.todeclaredPressureSwitch.Items = {'Force', 'Limit'};
            app.todeclaredPressureSwitch.Position = [283 97 45 20];
            app.todeclaredPressureSwitch.Value = 'Force';

            % Create ConfirmManTargetButton
            app.ConfirmManTargetButton = uibutton(app.TargetManPanel, 'push');
            app.ConfirmManTargetButton.ButtonPushedFcn = createCallbackFcn(app, @ConfirmManTargetButtonPushed, true);
            app.ConfirmManTargetButton.Visible = 'off';
            app.ConfirmManTargetButton.Position = [260 17 100 23];
            app.ConfirmManTargetButton.Text = 'Confirm';

            % Create AddManualTargetButton
            app.AddManualTargetButton = uibutton(app.TargetingTab, 'push');
            app.AddManualTargetButton.ButtonPushedFcn = createCallbackFcn(app, @AddManualTargetButtonPushed, true);
            app.AddManualTargetButton.Position = [112 401 122 23];
            app.AddManualTargetButton.Text = 'Add Manual Target';

            % Create RemoveManualTargetButton
            app.RemoveManualTargetButton = uibutton(app.TargetingTab, 'push');
            app.RemoveManualTargetButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveManualTargetButtonPushed, true);
            app.RemoveManualTargetButton.Position = [104 371 139 23];
            app.RemoveManualTargetButton.Text = 'Remove Manual Target';

            % Create ManualTargetDropDownLabel
            app.ManualTargetDropDownLabel = uilabel(app.TargetingTab);
            app.ManualTargetDropDownLabel.HorizontalAlignment = 'right';
            app.ManualTargetDropDownLabel.Position = [15 461 91 22];
            app.ManualTargetDropDownLabel.Text = 'Manual Target #';

            % Create ManualTargetDropDown
            app.ManualTargetDropDown = uidropdown(app.TargetingTab);
            app.ManualTargetDropDown.Items = {};
            app.ManualTargetDropDown.DropDownOpeningFcn = createCallbackFcn(app, @ManualTargetDropDownOpening, true);
            app.ManualTargetDropDown.ValueChangedFcn = createCallbackFcn(app, @ManualTargetDropDownValueChanged, true);
            app.ManualTargetDropDown.Position = [121 461 100 22];
            app.ManualTargetDropDown.Value = {};

            % Create TargetRegPanel
            app.TargetRegPanel = uipanel(app.TargetingTab);
            app.TargetRegPanel.Title = 'Target 1';
            app.TargetRegPanel.Visible = 'off';
            app.TargetRegPanel.Position = [292 84 383 175];

            % Create BrainRegionDropDownLabel
            app.BrainRegionDropDownLabel = uilabel(app.TargetRegPanel);
            app.BrainRegionDropDownLabel.HorizontalAlignment = 'right';
            app.BrainRegionDropDownLabel.Position = [30 115 74 22];
            app.BrainRegionDropDownLabel.Text = 'Brain Region';

            % Create BrainRegionDropDown
            app.BrainRegionDropDown = uidropdown(app.TargetRegPanel);
            app.BrainRegionDropDown.Items = {};
            app.BrainRegionDropDown.Position = [119 115 237 22];
            app.BrainRegionDropDown.Value = {};

            % Create MinPointDistancemmEditField_2Label
            app.MinPointDistancemmEditField_2Label = uilabel(app.TargetRegPanel);
            app.MinPointDistancemmEditField_2Label.HorizontalAlignment = 'right';
            app.MinPointDistancemmEditField_2Label.Position = [15 34 140 22];
            app.MinPointDistancemmEditField_2Label.Text = 'Min. Point Distance (mm)';

            % Create MinPointDistancemmEditFieldReg
            app.MinPointDistancemmEditFieldReg = uieditfield(app.TargetRegPanel, 'numeric');
            app.MinPointDistancemmEditFieldReg.Position = [171 34 43 22];
            app.MinPointDistancemmEditFieldReg.Value = 5;

            % Create todeclaredPressureSwitchLabel
            app.todeclaredPressureSwitchLabel = uilabel(app.TargetRegPanel);
            app.todeclaredPressureSwitchLabel.HorizontalAlignment = 'center';
            app.todeclaredPressureSwitchLabel.Position = [252 50 116 22];
            app.todeclaredPressureSwitchLabel.Text = 'to declared Pressure';

            % Create todeclaredPressureSwitchReg
            app.todeclaredPressureSwitchReg = uiswitch(app.TargetRegPanel, 'slider');
            app.todeclaredPressureSwitchReg.Items = {'Force', 'Limit'};
            app.todeclaredPressureSwitchReg.Position = [286 79 45 20];
            app.todeclaredPressureSwitchReg.Value = 'Force';

            % Create PressureAmplitudekPaEditField_2Label
            app.PressureAmplitudekPaEditField_2Label = uilabel(app.TargetRegPanel);
            app.PressureAmplitudekPaEditField_2Label.HorizontalAlignment = 'right';
            app.PressureAmplitudekPaEditField_2Label.Position = [20 65 141 22];
            app.PressureAmplitudekPaEditField_2Label.Text = 'Pressure Amplitude (kPa)';

            % Create PressureAmplitudekPaEditFieldReg
            app.PressureAmplitudekPaEditFieldReg = uieditfield(app.TargetRegPanel, 'numeric');
            app.PressureAmplitudekPaEditFieldReg.Position = [172 65 40 22];
            app.PressureAmplitudekPaEditFieldReg.Value = 300;

            % Create ConfirmRegTargetButton
            app.ConfirmRegTargetButton = uibutton(app.TargetRegPanel, 'push');
            app.ConfirmRegTargetButton.ButtonPushedFcn = createCallbackFcn(app, @ConfirmRegTargetButtonPushed, true);
            app.ConfirmRegTargetButton.Visible = 'off';
            app.ConfirmRegTargetButton.Position = [259 10 100 23];
            app.ConfirmRegTargetButton.Text = 'Confirm';

            % Create AddRegionTargetButton
            app.AddRegionTargetButton = uibutton(app.TargetingTab, 'push');
            app.AddRegionTargetButton.ButtonPushedFcn = createCallbackFcn(app, @AddRegionTargetButtonPushed, true);
            app.AddRegionTargetButton.Visible = 'off';
            app.AddRegionTargetButton.Position = [113 153 122 23];
            app.AddRegionTargetButton.Text = 'Add Region Target';

            % Create RemoveRegionTargetButton
            app.RemoveRegionTargetButton = uibutton(app.TargetingTab, 'push');
            app.RemoveRegionTargetButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveRegionTargetButtonPushed, true);
            app.RemoveRegionTargetButton.Visible = 'off';
            app.RemoveRegionTargetButton.Position = [105 123 138 23];
            app.RemoveRegionTargetButton.Text = 'Remove Region Target';

            % Create RegionTargetDropDownLabel
            app.RegionTargetDropDownLabel = uilabel(app.TargetingTab);
            app.RegionTargetDropDownLabel.HorizontalAlignment = 'right';
            app.RegionTargetDropDownLabel.Visible = 'off';
            app.RegionTargetDropDownLabel.Position = [17 213 90 22];
            app.RegionTargetDropDownLabel.Text = 'Region Target #';

            % Create RegionTargetDropDown
            app.RegionTargetDropDown = uidropdown(app.TargetingTab);
            app.RegionTargetDropDown.Items = {};
            app.RegionTargetDropDown.DropDownOpeningFcn = createCallbackFcn(app, @RegionTargetDropDownOpening, true);
            app.RegionTargetDropDown.ValueChangedFcn = createCallbackFcn(app, @RegionTargetDropDownValueChanged, true);
            app.RegionTargetDropDown.Visible = 'off';
            app.RegionTargetDropDown.ClickedFcn = createCallbackFcn(app, @RegionTargetDropDownClicked, true);
            app.RegionTargetDropDown.Position = [122 213 100 22];
            app.RegionTargetDropDown.Value = {};

            % Create UpdateTargetingButton
            app.UpdateTargetingButton = uibutton(app.TargetingTab, 'push');
            app.UpdateTargetingButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateTargetingButtonPushed, true);
            app.UpdateTargetingButton.Position = [574 24 100 23];
            app.UpdateTargetingButton.Text = 'Update';

            % Create MaxPressurekPaEditFieldLabel
            app.MaxPressurekPaEditFieldLabel = uilabel(app.TargetingTab);
            app.MaxPressurekPaEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxPressurekPaEditFieldLabel.Position = [66 24 114 22];
            app.MaxPressurekPaEditFieldLabel.Text = 'Max. Pressure (kPa)';

            % Create MaxPressurekPaEditField
            app.MaxPressurekPaEditField = uieditfield(app.TargetingTab, 'numeric');
            app.MaxPressurekPaEditField.ValueChangedFcn = createCallbackFcn(app, @MaxPressurekPaEditFieldValueChanged, true);
            app.MaxPressurekPaEditField.Position = [191 24 44 22];
            app.MaxPressurekPaEditField.Value = 200;

            % Create LimitIntracranialOffTargetPressureCheckBox
            app.LimitIntracranialOffTargetPressureCheckBox = uicheckbox(app.TargetingTab);
            app.LimitIntracranialOffTargetPressureCheckBox.ValueChangedFcn = createCallbackFcn(app, @LimitIntracranialOffTargetPressureCheckBoxValueChanged, true);
            app.LimitIntracranialOffTargetPressureCheckBox.Text = 'Limit Intracranial Off-Target Pressure';
            app.LimitIntracranialOffTargetPressureCheckBox.Position = [50 50 219 22];
            app.LimitIntracranialOffTargetPressureCheckBox.Value = true;

            % Create LimitSkullPressureCheckBox
            app.LimitSkullPressureCheckBox = uicheckbox(app.TargetingTab);
            app.LimitSkullPressureCheckBox.ValueChangedFcn = createCallbackFcn(app, @LimitSkullPressureCheckBoxValueChanged, true);
            app.LimitSkullPressureCheckBox.Text = 'Limit Skull Pressure';
            app.LimitSkullPressureCheckBox.Position = [317 50 128 22];

            % Create MaxPressurekPaSkullEditFieldLabel
            app.MaxPressurekPaSkullEditFieldLabel = uilabel(app.TargetingTab);
            app.MaxPressurekPaSkullEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxPressurekPaSkullEditFieldLabel.Visible = 'off';
            app.MaxPressurekPaSkullEditFieldLabel.Position = [333 24 114 22];
            app.MaxPressurekPaSkullEditFieldLabel.Text = 'Max. Pressure (kPa)';

            % Create MaxPressurekPaSkullEditField
            app.MaxPressurekPaSkullEditField = uieditfield(app.TargetingTab, 'numeric');
            app.MaxPressurekPaSkullEditField.Visible = 'off';
            app.MaxPressurekPaSkullEditField.Position = [458 24 44 22];
            app.MaxPressurekPaSkullEditField.Value = 1000;

            % Create OptimizeTab
            app.OptimizeTab = uitab(app.TabGroup);
            app.OptimizeTab.Title = 'Optimize';
            app.OptimizeTab.ButtonDownFcn = createCallbackFcn(app, @OptimizeTabButtonDown, true);

            % Create DoGroundTruthSimulationButton
            app.DoGroundTruthSimulationButton = uibutton(app.OptimizeTab, 'push');
            app.DoGroundTruthSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @DoGroundTruthSimulationButtonPushed, true);
            app.DoGroundTruthSimulationButton.Position = [534 496 164 23];
            app.DoGroundTruthSimulationButton.Text = 'Do Ground Truth Simulation';

            % Create ComputeInitialSolutionButton
            app.ComputeInitialSolutionButton = uibutton(app.OptimizeTab, 'push');
            app.ComputeInitialSolutionButton.ButtonPushedFcn = createCallbackFcn(app, @ComputeInitialSolutionButtonPushed, true);
            app.ComputeInitialSolutionButton.Position = [311 495 142 23];
            app.ComputeInitialSolutionButton.Text = 'Compute Initial Solution';

            % Create OptimizeButton
            app.OptimizeButton = uibutton(app.OptimizeTab, 'push');
            app.OptimizeButton.ButtonPushedFcn = createCallbackFcn(app, @OptimizeButtonPushed, true);
            app.OptimizeButton.Position = [311 441 142 45];
            app.OptimizeButton.Text = 'Optimize';

            % Create OptimizationModeButtonGroup
            app.OptimizationModeButtonGroup = uibuttongroup(app.OptimizeTab);
            app.OptimizationModeButtonGroup.Title = 'Optimization Mode';
            app.OptimizationModeButtonGroup.Position = [14 387 224 132];

            % Create OptimizeforPhasesandAmplitudesButton
            app.OptimizeforPhasesandAmplitudesButton = uiradiobutton(app.OptimizationModeButtonGroup);
            app.OptimizeforPhasesandAmplitudesButton.Text = {'Optimize for Phases '; 'and Amplitudes'};
            app.OptimizeforPhasesandAmplitudesButton.Position = [11 78 133 30];
            app.OptimizeforPhasesandAmplitudesButton.Value = true;

            % Create OptimizeforPhasesand1AmpperTransducerButton
            app.OptimizeforPhasesand1AmpperTransducerButton = uiradiobutton(app.OptimizationModeButtonGroup);
            app.OptimizeforPhasesand1AmpperTransducerButton.Text = {'Optimize for Phases '; 'and 1 Amp per Transducer'};
            app.OptimizeforPhasesand1AmpperTransducerButton.Position = [11 41 163 30];

            % Create OptimizeforPhasesand1AmpButton
            app.OptimizeforPhasesand1AmpButton = uiradiobutton(app.OptimizationModeButtonGroup);
            app.OptimizeforPhasesand1AmpButton.Text = {'Optimize for Phases '; 'and 1 Amp'};
            app.OptimizeforPhasesand1AmpButton.Position = [11 7 133 30];

            % Create GroundTruthResolutionFactorEditFieldLabel
            app.GroundTruthResolutionFactorEditFieldLabel = uilabel(app.OptimizeTab);
            app.GroundTruthResolutionFactorEditFieldLabel.HorizontalAlignment = 'right';
            app.GroundTruthResolutionFactorEditFieldLabel.Position = [551 455 99 30];
            app.GroundTruthResolutionFactorEditFieldLabel.Text = {'Ground Truth '; 'Resolution Factor'};

            % Create GroundTruthResolutionFactorEditField
            app.GroundTruthResolutionFactorEditField = uieditfield(app.OptimizeTab, 'numeric');
            app.GroundTruthResolutionFactorEditField.Position = [662 463 34 22];
            app.GroundTruthResolutionFactorEditField.Value = 1;

            % Create DisplayAdvancedOptimizationOptionsCheckBox
            app.DisplayAdvancedOptimizationOptionsCheckBox = uicheckbox(app.OptimizeTab);
            app.DisplayAdvancedOptimizationOptionsCheckBox.ValueChangedFcn = createCallbackFcn(app, @DisplayAdvancedOptimizationOptionsCheckBoxValueChanged, true);
            app.DisplayAdvancedOptimizationOptionsCheckBox.Text = 'Display Advanced Optimization Options';
            app.DisplayAdvancedOptimizationOptionsCheckBox.Position = [22 217 234 22];

            % Create AdvancedOptimizationOptionsPanel
            app.AdvancedOptimizationOptionsPanel = uipanel(app.OptimizeTab);
            app.AdvancedOptimizationOptionsPanel.Title = 'Advanced Optimization Options';
            app.AdvancedOptimizationOptionsPanel.Visible = 'off';
            app.AdvancedOptimizationOptionsPanel.Position = [25 17 236 187];

            % Create RelFunctionToleranceEditFieldLabel
            app.RelFunctionToleranceEditFieldLabel = uilabel(app.AdvancedOptimizationOptionsPanel);
            app.RelFunctionToleranceEditFieldLabel.HorizontalAlignment = 'right';
            app.RelFunctionToleranceEditFieldLabel.Position = [6 106 131 22];
            app.RelFunctionToleranceEditFieldLabel.Text = 'Rel. Function Tolerance';

            % Create RelFunctionToleranceEditField
            app.RelFunctionToleranceEditField = uieditfield(app.AdvancedOptimizationOptionsPanel, 'numeric');
            app.RelFunctionToleranceEditField.Position = [142 106 52 22];
            app.RelFunctionToleranceEditField.Value = 1e-08;

            % Create ConstraintToleranceLabel
            app.ConstraintToleranceLabel = uilabel(app.AdvancedOptimizationOptionsPanel);
            app.ConstraintToleranceLabel.HorizontalAlignment = 'right';
            app.ConstraintToleranceLabel.Position = [22 76 115 22];
            app.ConstraintToleranceLabel.Text = 'Constraint Tolerance';

            % Create ConstraintToleranceEditField
            app.ConstraintToleranceEditField = uieditfield(app.AdvancedOptimizationOptionsPanel, 'numeric');
            app.ConstraintToleranceEditField.Position = [142 76 50 22];
            app.ConstraintToleranceEditField.Value = 0.001;

            % Create MaxNoofIterationsEditFieldLabel
            app.MaxNoofIterationsEditFieldLabel = uilabel(app.AdvancedOptimizationOptionsPanel);
            app.MaxNoofIterationsEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxNoofIterationsEditFieldLabel.Position = [22 46 116 22];
            app.MaxNoofIterationsEditFieldLabel.Text = 'Max. No of Iterations';

            % Create MaxNoofIterationsEditField
            app.MaxNoofIterationsEditField = uieditfield(app.AdvancedOptimizationOptionsPanel, 'numeric');
            app.MaxNoofIterationsEditField.RoundFractionalValues = 'on';
            app.MaxNoofIterationsEditField.Position = [143 46 50 22];
            app.MaxNoofIterationsEditField.Value = 100;

            % Create MinimizeArrrayAmplitudesCheckBox
            app.MinimizeArrrayAmplitudesCheckBox = uicheckbox(app.AdvancedOptimizationOptionsPanel);
            app.MinimizeArrrayAmplitudesCheckBox.ValueChangedFcn = createCallbackFcn(app, @MinimizeArrrayAmplitudesCheckBoxValueChanged, true);
            app.MinimizeArrrayAmplitudesCheckBox.Text = 'Minimize Arrray Amplitudes';
            app.MinimizeArrrayAmplitudesCheckBox.Position = [13 136 167 22];

            % Create LxNormEditFieldLabel
            app.LxNormEditFieldLabel = uilabel(app.AdvancedOptimizationOptionsPanel);
            app.LxNormEditFieldLabel.HorizontalAlignment = 'right';
            app.LxNormEditFieldLabel.Position = [86 12 50 22];
            app.LxNormEditFieldLabel.Text = 'Lx Norm';

            % Create LxNormEditField
            app.LxNormEditField = uieditfield(app.AdvancedOptimizationOptionsPanel, 'numeric');
            app.LxNormEditField.RoundFractionalValues = 'on';
            app.LxNormEditField.Position = [143 12 50 22];
            app.LxNormEditField.Value = 10;

            % Create PlotPanel
            app.PlotPanel = uipanel(app.OptimizeTab);
            app.PlotPanel.Title = 'Plot';
            app.PlotPanel.Position = [14 253 309 119];

            % Create LimittoIntracranialFieldCheckBox
            app.LimittoIntracranialFieldCheckBox = uicheckbox(app.PlotPanel);
            app.LimittoIntracranialFieldCheckBox.Text = 'Limit to Intracranial Field';
            app.LimittoIntracranialFieldCheckBox.Position = [24 68 153 22];

            % Create LimittoSkullCheckBox
            app.LimittoSkullCheckBox = uicheckbox(app.PlotPanel);
            app.LimittoSkullCheckBox.Text = 'Limit to Skull';
            app.LimittoSkullCheckBox.Position = [199 69 90 22];

            % Create MaskPressurePlotkPaEditFieldLabel
            app.MaskPressurePlotkPaEditFieldLabel = uilabel(app.PlotPanel);
            app.MaskPressurePlotkPaEditFieldLabel.HorizontalAlignment = 'right';
            app.MaskPressurePlotkPaEditFieldLabel.Position = [19 21 141 22];
            app.MaskPressurePlotkPaEditFieldLabel.Text = 'Mask Pressure Plot (kPa)';

            % Create MaskPressurePlotkPaEditField
            app.MaskPressurePlotkPaEditField = uieditfield(app.PlotPanel, 'numeric');
            app.MaskPressurePlotkPaEditField.Position = [169 21 52 22];

            % Create opt_mode
            app.opt_mode = uicheckbox(app.OptimizeTab);
            app.opt_mode.Text = 'Use double-iterative mode';
            app.opt_mode.Position = [312 410 162 22];
            app.opt_mode.Value = true;

            % Create PlotResultsButton
            app.PlotResultsButton = uibutton(app.OptimizeTab, 'push');
            app.PlotResultsButton.ButtonPushedFcn = createCallbackFcn(app, @PlotResultsButtonPushed, true);
            app.PlotResultsButton.Position = [335 296 118 34];
            app.PlotResultsButton.Text = 'Plot Results';

            % Create Slice30Label
            app.Slice30Label = uilabel(app.UIFigure);
            app.Slice30Label.Position = [653 6 67 22];
            app.Slice30Label.Text = 'Slice 30';

            % Create CloseallPlotsButton
            app.CloseallPlotsButton = uibutton(app.UIFigure, 'push');
            app.CloseallPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @CloseallPlotsButtonPushed, true);
            app.CloseallPlotsButton.Position = [8 6 100 23];
            app.CloseallPlotsButton.Text = 'Close all Plots';

            % Create UpdateSliceButton
            app.UpdateSliceButton = uibutton(app.UIFigure, 'push');
            app.UpdateSliceButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateSliceButtonPushed2, true);
            app.UpdateSliceButton.Visible = 'off';
            app.UpdateSliceButton.Position = [543 5 100 23];
            app.UpdateSliceButton.Text = 'Update Slice';

            % Create SliceDirectionLabel
            app.SliceDirectionLabel = uilabel(app.UIFigure);
            app.SliceDirectionLabel.HorizontalAlignment = 'right';
            app.SliceDirectionLabel.Visible = 'off';
            app.SliceDirectionLabel.Position = [390 6 82 22];
            app.SliceDirectionLabel.Text = 'Slice Direction';

            % Create SliceDirectionDropDown_2
            app.SliceDirectionDropDown_2 = uidropdown(app.UIFigure);
            app.SliceDirectionDropDown_2.Items = {'X', 'Y', 'Z'};
            app.SliceDirectionDropDown_2.ValueChangedFcn = createCallbackFcn(app, @SliceDirectionDropDown_2ValueChanged, true);
            app.SliceDirectionDropDown_2.Visible = 'off';
            app.SliceDirectionDropDown_2.Position = [487 6 45 22];
            app.SliceDirectionDropDown_2.Value = 'Y';

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