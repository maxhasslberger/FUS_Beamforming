classdef simulationApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        TabGroup           matlab.ui.container.TabGroup
        InitTab            matlab.ui.container.Tab
        Switch             matlab.ui.control.Switch
        Label              matlab.ui.control.Label
        TransducerTab      matlab.ui.container.Tab
        TargetingTab       matlab.ui.container.Tab
        InspectResultsTab  matlab.ui.container.Tab
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
            app.Label.Position = [54 469 25 22];
            app.Label.Text = ' ';

            % Create Switch
            app.Switch = uiswitch(app.InitTab, 'slider');
            app.Switch.Items = {'2D', '3D'};
            app.Switch.Position = [43 506 45 20];
            app.Switch.Value = '2D';

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