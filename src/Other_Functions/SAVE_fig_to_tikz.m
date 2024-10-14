function [ ] = SAVE_fig_to_tikz( filename, cmap )
%UNTITLED Summary of this function goes here
% Detailed explanation goes here
set(gcf,'colormap', cmap);

for i = 1:2
matlab2tikz([filename,'.tex'],'width','\fw','height','\fh','extraaxisoptions',['title style={font=\small},'...
'xlabel style={font=\small},'...
'ylabel style={font=\small},',...
'legend style={font=\small},',...
'ticklabel style={font=\small}']);
fprintf('\n\nFINISHED!\n');

set(gcf,'colormap', gray());
filename = [filename,'_gray'];
end
end
