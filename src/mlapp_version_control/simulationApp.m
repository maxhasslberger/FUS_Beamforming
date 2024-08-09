classdef simulationApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        TabGroup                      matlab.ui.container.TabGroup
        InitTab                       matlab.ui.container.Tab
        SaveResultsCheckBox           matlab.ui.control.CheckBox
        GreensFunctionbasedCheckBox   matlab.ui.control.CheckBox
        MediumSwitchLabel             matlab.ui.control.Label
        MediumSwitch                  matlab.ui.control.Switch
        RealxmmLabel                  matlab.ui.control.Label
        SpatialResolutionmmEditField  matlab.ui.control.NumericEditField
        SpatialResolutionmmEditFieldLabel  matlab.ui.control.Label
        CenterFreqkHzEditField        matlab.ui.control.NumericEditField
        CenterFreqkHzEditFieldLabel   matlab.ui.control.Label
        DimSwitch                     matlab.ui.control.Switch
        Label                         matlab.ui.control.Label
        TransducerTab                 matlab.ui.container.Tab
        TargetingTab                  matlab.ui.container.Tab
        InspectResultsTab             matlab.ui.container.Tab
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1064 572];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [2 1 715 572];

            % Create InitTab
            app.InitTab = uitab(app.TabGroup);
            app.InitTab.Title = 'Init';

            % Create Label
            app.Label = uilabel(app.InitTab);
            app.Label.HorizontalAlignment = 'center';
            app.Label.Position = [54 488 26 22];
            app.Label.Text = 'Dim';

            % Create DimSwitch
            app.DimSwitch = uiswitch(app.InitTab, 'slider');
            app.DimSwitch.Items = {'2D', '3D'};
            app.DimSwitch.Position = [44 458 45 20];
            app.DimSwitch.Value = '2D';

            % Create CenterFreqkHzEditFieldLabel
            app.CenterFreqkHzEditFieldLabel = uilabel(app.InitTab);
            app.CenterFreqkHzEditFieldLabel.HorizontalAlignment = 'right';
            app.CenterFreqkHzEditFieldLabel.Position = [184 456 101 22];
            app.CenterFreqkHzEditFieldLabel.Text = 'Center Freq (kHz)';

            % Create CenterFreqkHzEditField
            app.CenterFreqkHzEditField = uieditfield(app.InitTab, 'numeric');
            app.CenterFreqkHzEditField.Position = [297 456 46 22];
            app.CenterFreqkHzEditField.Value = 500;

            % Create SpatialResolutionmmEditFieldLabel
            app.SpatialResolutionmmEditFieldLabel = uilabel(app.InitTab);
            app.SpatialResolutionmmEditFieldLabel.HorizontalAlignment = 'right';
            app.SpatialResolutionmmEditFieldLabel.Position = [23 276 133 22];
            app.SpatialResolutionmmEditFieldLabel.Text = 'Spatial Resolution (mm)';

            % Create SpatialResolutionmmEditField
            app.SpatialResolutionmmEditField = uieditfield(app.InitTab, 'numeric');
            app.SpatialResolutionmmEditField.Position = [171 276 45 22];

            % Create RealxmmLabel
            app.RealxmmLabel = uilabel(app.InitTab);
            app.RealxmmLabel.Position = [27 255 80 22];
            app.RealxmmLabel.Text = '-> Real: x mm';

            % Create MediumSwitch
            app.MediumSwitch = uiswitch(app.InitTab, 'slider');
            app.MediumSwitch.Items = {'Homogeneous', 'Heterogeneous'};
            app.MediumSwitch.Position = [508 457 45 20];
            app.MediumSwitch.Value = 'Homogeneous';

            % Create MediumSwitchLabel
            app.MediumSwitchLabel = uilabel(app.InitTab);
            app.MediumSwitchLabel.HorizontalAlignment = 'center';
            app.MediumSwitchLabel.Position = [507 488 48 22];
            app.MediumSwitchLabel.Text = 'Medium';

            % Create GreensFunctionbasedCheckBox
            app.GreensFunctionbasedCheckBox = uicheckbox(app.InitTab);
            app.GreensFunctionbasedCheckBox.Text = {'Green''s Function'; 'based'};
            app.GreensFunctionbasedCheckBox.Position = [424 396 113 30];
            app.GreensFunctionbasedCheckBox.Value = true;

            % Create SaveResultsCheckBox
            app.SaveResultsCheckBox = uicheckbox(app.InitTab);
            app.SaveResultsCheckBox.Text = 'Save Results';
            app.SaveResultsCheckBox.Position = [27 174 93 22];

            % Create TransducerTab
            app.TransducerTab = uitab(app.TabGroup);
            app.TransducerTab.Title = 'Transducer';

            % Create TargetingTab
            app.TargetingTab = uitab(app.TabGroup);
            app.TargetingTab.Title = 'Targeting';

            % Create InspectResultsTab
            app.InspectResultsTab = uitab(app.TabGroup);
            app.InspectResultsTab.Title = 'Inspect Results';

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