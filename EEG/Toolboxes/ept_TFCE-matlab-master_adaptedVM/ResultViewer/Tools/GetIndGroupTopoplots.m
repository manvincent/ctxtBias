% Load the results file...

sample = 92;

e_loc=Info.Electrodes.e_loc;

a = cellfun(@(x) x(:,:,sample), Data, 'UniformOutput', false);
b = cellfun(@mean, a, 'UniformOutput', false);

% New dataset for that particular sample
nData = (cell2mat(b(:,1))+ cell2mat(b(:,2)))/2;


H.Figure = figure;
set(H.Figure,...
    'Name',             'ERP Topoplots'    ,...
    'NumberTitle',      'off'               ,...
    'Units',            'pixels'            ,...
    'Position',         [100 100 1000 500]   ,...
    'Color',            'w'                 );

% H.Axes = axes;
% set(H.Axes,...
%     'FontName',         'Helvetica'         ,...
%     'ColorOrder',       flipud(hsv(nL))     ,...     % set the new ColorOrder to the first three scaled colours of HSV colormap
%     'NextPlot',         'replacechildren'   );

% Plot the actual amplitude topoplots for each group  
topoplot(nData(3,:), e_loc,...
    'electrodes',           'on',           ...         % display markers ("labels" shows electrode names
    'whitebk',              'on',           ...
    'colormap',             hot,            ...
    'shading',              'interp',       ...         % (flat or interp) useless with style is "fill"
    'style',                'fill',         ...
    'numcontour',           10,             ...         % Increase for more contour details
    'maplimits',            [-10,5],         ...
    'plotrad',              max(abs(cell2mat({e_loc.radius})))      ); % 'map' 'contour' 'both' 'fill' 'blank'
    
% Plot the TFCE_Obs Values
topoplot(Results.TFCE_Obs.A(:,sample)', e_loc,...
        'electrodes',           'on',           ...         % display markers ("labels" shows electrode names
        'whitebk',              'on',           ...
        'colormap',             hot,            ...
        'shading',              'interp',       ...         % (flat or interp) useless with style is "fill"
        'style',                'fill',         ...
        'numcontour',           6,             ...         % Increase for more contour details
        'plotrad',              max(abs(cell2mat({e_loc.radius})))      ); % 'map' 'contour' 'both' 'fill' 'blank'