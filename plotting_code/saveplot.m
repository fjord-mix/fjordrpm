function saveplot(fw, fh, res, filename)

% SAVEPLOT Save a plot.
%   SAVEPLOT(FW, FH, RES, FILENAME) saves a plot with paper size FW, FH,
%   resolution RES and filename FILENAME.

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'inverthardcopy', 'off');
set(gcf, 'PaperPosition', [0 0 fw fh]);
set(gcf, 'PaperSize', [fw fh]);
set(gcf, 'units', 'centimeters', 'position', [0 0 fw fh]);
set(gcf, 'color', 'w');
res_string = ['-r', num2str(res)];
print('-dpng', res_string, filename);

end