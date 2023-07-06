function animate(name, nframes)

% ANIMATE Animate the results of a box model simulation.
%   ANIMATE(NAME, NFRAMES) creates a video mp4 file of a box model
%   simulation showing the temperature and salinity in each layer, the
%   potential and the subglacial discharge.

%% File setup
% Load the results of the required box model simulation.
load(['./output_', name, '/run_params.mat'], 'p');
load(['./output_', name, '/out_boxmodel.mat'], 's', 'f');

% Set location for animation files.
addpath(['output_', name])
outputfile = ['output_', name, '/output_animation'];

%% Plotting features
% Change default plot line width to 1.
set(0, 'DefaultLineLineWidth', 1);

% Change all default interpreters to latex.
list_factory = fieldnames(get(groot, 'factory'));
index_interpreter = find(contains(list_factory, 'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)}, 'factory', 'default');
    set(groot, default_name, 'latex');
end

% Axis limits, legend, subglacial discharge.
Slims = [min([min(s.S(:)), min(f.Ss(:))]), max([max(s.S(:)), max(f.Ss(:))])];
Tlims = [min([min(s.T(:)), min(f.Ts(:))]), max([max(s.T(:)), max(f.Ts(:))])];
zlims = [-p.H, 0];
x0 = 1500;
sf = x0*0.2/max([max(abs(s.QVs(:))), max(abs(s.QVg(:))), max(abs(s.QVb(:)))]);
legstr = cell(1, size(s.H, 1)); 
for i = 1:size(s.H, 1)
    legstr{i} = [legstr, num2str(i)];
end
f_qsg = interp1(linspace(0, max(s.t), length(f.Qsg)), f.Qsg, s.t, 'linear');

% Plot positions and parameters
fs = 10;
fs2 = 6;
lspace = 0.06;
hspace1 = 0.03;
hspace2 = 0.08;
rspace = 0.02;
pw1 = 0.14;
pw2 = pw1/2;
pw3 = 1-lspace-rspace-2*pw1-2*pw2-3*hspace1-hspace2;
bspace = 0.12;
tspace = 0.06;
vspace = 0.05;
ph = 1-bspace-tspace;
ph3 = (1-bspace-tspace-2*vspace)/3;
cby = 0.075;
cbh = 0.02;
cbw = 0.2;

%% Create animation files
for i = 1:round((length(s.t)-1)/nframes):length(s.t)-1
    figure(); 
    set(gcf, 'Visible', 'off');

    ints = [0; cumsum(s.H(:, i))];
    y = -0.5*(ints(1:end-1)+ints(2:end));
    y2 = cumsum(s.H(1:end-1, i));

    % Box model plots
    a1 = axes('position', [lspace, bspace, pw1, ph]); 
    hold on;
    pcolor([0, x0], -[0;cumsum(s.H(:, i))], [[s.S(:, i);0], [s.S(:, i);0]]);
    set(gca, 'xtick', [], 'xticklabel', {});
    q1 = quiver(0*y+x0, y, sf*s.QVs(:, i), 0*s.QVs(:, i), 'autoscale', 'off');
    set(q1, 'color', 'k', 'linewidth', 2);
    q2 = quiver(0*y, y, sf*s.QVg(:, i), 0*s.QVg(:, i), 'autoscale', 'off');
    set(q2, 'color', 'k', 'linewidth', 2);
    q3 = quiver(0*y2+x0/2, -y2, 0*y2, sf*cumsum(s.QVk(1:end-1, i)), 'autoscale', 'off');
    set(q3, 'color', 'k', 'linewidth', 2);
    q4 = quiver(0*y2+x0/3, -y2, 0*y2, sf*cumsum(s.QVb(1:end-1, i)), 'autoscale', 'off');
    set(q4, 'color', 'r', 'linewidth', 2);
    set(gca, 'clipping', 'off', 'box', 'on', 'fontsize', fs);
    clim(a1, Slims); xlim([0, x0]); ylim(zlims);
    text(0.02*x0, -800, ['t = ', num2str(0.01*round(100*s.t(i))), ' days'], ...
        'fontsize', fs, 'VerticalAlignment', 'bottom');
    title('Salinity: fjord', 'fontsize', fs);
    ylabel('depth (m)', 'fontsize', fs);

    a2 = axes('position', [lspace+pw1+hspace1, bspace, pw2, ph]);
    pcolor([0, x0], f.zs, [f.Ss(:, i), f.Ss(:, i)]); 
    shading flat;
    xlim([0, x0]); ylim(zlims);
    set(gca, 'box', 'on', 'fontsize', fs, 'xtick', [], 'ytick', [], 'yticklabel', {}, 'clipping', 'off');
    title('shelf', 'fontsize', fs);
    colorbar(a2, 'southoutside', 'position', [lspace+0.5*(pw1+hspace1+pw2)-cbw/2, cby, cbw, cbh], 'fontsize', fs);
    clim(a2, Slims);

    a3 = axes('position', [lspace+pw1+2*hspace1+pw2, bspace, pw1, ph]); 
    hold on;
    pcolor([0, x0], -[0;cumsum(s.H(:, i))], [[s.T(:, i);0], [s.T(:, i);0]]);
    set(gca, 'xtick', [], 'xticklabel', {});
    q1 = quiver(0*y+x0, y, sf*s.QVs(:, i), 0*s.QVs(:, i), 'autoscale', 'off');
    set(q1, 'color', 'k', 'linewidth', 2);
    q2 = quiver(0*y, y, sf*s.QVg(:, i), 0*s.QVg(:, i), 'autoscale', 'off');
    set(q2, 'color', 'k', 'linewidth', 2);
    q3 = quiver(0*y2+x0/2, -y2, 0*y2, sf*cumsum(s.QVk(1:end-1, i)), 'autoscale', 'off');
    set(q3, 'color', 'k', 'linewidth', 2);
    q4 = quiver(0*y2+x0/3, -y2, 0*y2, sf*cumsum(s.QVb(1:end-1, i)), 'autoscale', 'off');
    set(q4, 'color', 'r', 'linewidth', 2);
    set(gca, 'clipping', 'off', 'box', 'on', 'fontsize', fs, 'xtick', [], 'ytick', []);
    clim(gca, Tlims); xlim([0, x0]); ylim(zlims);
    title('Temperature: fjord', 'fontsize', fs);
    colormap(a3, jet);

    a4 = axes('position', [lspace+2*pw1+pw2+3*hspace1, bspace, pw2, ph]);
    pcolor([0, x0], f.zs, [f.Ts(:, i), f.Ts(:, i)]); 
    shading flat;
    xlim([0, x0]); ylim(zlims);
    set(gca, 'box', 'on', 'fontsize', fs, 'xtick', [], 'ytick', [], 'yticklabel', {}, 'clipping', 'off');
    title('shelf', 'fontsize', fs);
    colorbar(a4, 'southoutside', 'position', [lspace+pw1+2*hspace1+pw2+0.5*(pw1+hspace1+pw2)-cbw/2, cby, cbw, cbh], 'fontsize', fs);
    colormap(a4, jet); 
    clim(a4, Tlims);

    % Illustration of sill
    if p.sill==1
        annotation('line', (lspace+pw1+hspace1/2)*[1, 1], bspace+[0, ph*(1-p.silldepth/p.H)], ...
            'color', 'k', 'linewidth', 2);
        annotation('line', (lspace+2*pw1+pw2+2.5*hspace1)*[1, 1], bspace+[0, ph*(1-p.silldepth/p.H)], ...
            'color', 'k', 'linewidth', 2);
    end

    % Time series
    axes('position', [lspace+2*pw1+2*pw2+3*hspace1+hspace2, bspace+2*(vspace+ph3), pw3, ph3]); hold on;
    plot(s.t, f_qsg, 'linewidth', 1);
    plot(s.t(i), f_qsg(i), 'k.');
    ylabel('subglacial discharge (m$^3$/s)', 'fontsize', fs2);
    set(gca, 'box', 'on', 'fontsize', fs2);

    axes('position', [lspace+2*pw1+2*pw2+3*hspace1+hspace2, bspace+1*(vspace+ph3), pw3, ph3]); hold on;
    plot(s.t, s.H, 'linewidth', 1);
    plot(s.t(i), s.H(:, i), 'k.');
    ylabel('$H$ (m)', 'fontsize', fs2);
    set(gca, 'box', 'on', 'fontsize', fs2);
    legend(legstr, 'location', 'north', 'orientation', 'horizontal', 'fontsize', 4, 'NumColumns', 2);

    axes('position', [lspace+2*pw1+2*pw2+3*hspace1+hspace2, bspace+0*(vspace+ph3), pw3, ph3]); hold on;
    plot(s.t, s.phi, 'linewidth', 1);
    plot(s.t(i), s.phi(:, i), 'k.');
    xlabel('$t$ (days)', 'fontsize', fs2);
    ylabel('$\phi$ (shelf-fjord pressure difference)', 'fontsize', fs2);
    set(gca, 'box', 'on', 'fontsize', fs2);

    % Save files
    if i < 10
        savenum = ['00', num2str(i)];
    elseif i<100
        savenum = ['0', num2str(i)];
    else
        savenum = num2str(i);
    end
    saveplot(25, 10, 300, [outputfile, '_', savenum, '.png']);
    pause
    close all;

end

%% Write video
video = VideoWriter(['output_', name, '/animation_', name, '.mp4'], 'MPEG-4');
video.FrameRate = 10;
open(video);
for i = 1:round((length(s.t)-1)/nframes):length(s.t)-1

    if i < 10
        savenum = ['00', num2str(i)];
    elseif i < 100
        savenum = ['0', num2str(i)];
    else
        savenum = num2str(i);
    end

    I = imread([outputfile, '_', savenum, '.png']);
    writeVideo(video, I);

end
close(video);
delete([outputfile, '*.png'])

end
