function plotrpm(p, s)

% PLOTRPM makes summary plots after the model has run.
%   PLOTRPM(p, s) makes example plots of temperature, salinity,
%   fluxes and melt rates to show the sort of thing that can be done.
%   Inputs are the parameters structure p and solution structure s.

% time indices to plot
ip = [1:round(length(s.t)/50):length(s.t)];

% colourmaps
cmapt = parula(length(ip));
cmapz = parula(length(s.z));

% line width
lw = 2;

%% temperature
figure();

tcl = tiledlayout(2,2);
title(tcl,'TEMPERATURES');

% fjord profiles
nexttile; hold on;
for i=1:length(ip)
    plot(s.T(:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('temperature (C)'); ylabel('depth (m)');
h = colorbar; colormap(gca,cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)'); ylim([-p.H,0]);
set(gca,'box','on'); grid on; title('FJORD');

% fjord time series
nexttile; hold on;
for i=1:length(s.z)
    plot(s.t,s.T(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('temperature (C)');
h = colorbar; colormap(gca,cmapz); caxis([-p.H,0]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on; title('FJORD');

% shelf profiles
nexttile; hold on;
for i=1:length(ip)
    plot(s.Ts(:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('temperature (C)'); ylabel('depth (m)');
h = colorbar; colormap(gca,cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)'); ylim([-p.H,0]);
set(gca,'box','on'); grid on; title('SHELF');

% shelf time series
nexttile; hold on;
for i=1:length(s.z)
    plot(s.t,s.Ts(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('temperature (C)');
h = colorbar; colormap(gca,cmapz); caxis([-p.H,0]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on; title('SHELF');

%% salinity
figure();

tcl = tiledlayout(2,2);
title(tcl,'SALINITIES');

% fjord profiles
nexttile; hold on;
for i=1:length(ip)
    plot(s.S(:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('salinity'); ylabel('depth (m)');
h = colorbar; colormap(gca,cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)'); ylim([-p.H,0]);
set(gca,'box','on'); grid on; title('FJORD');

% fjord time series
nexttile; hold on;
for i=1:length(s.z)
    plot(s.t,s.S(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('salinity');
h = colorbar; colormap(gca,cmapz); caxis([-p.H,0]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on; title('FJORD');

% shelf profiles
nexttile; hold on;
for i=1:length(ip)
    plot(s.Ss(:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('salinity'); ylabel('depth (m)');
h = colorbar; colormap(gca,cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)'); ylim([-p.H,0]);
set(gca,'box','on'); grid on; title('SHELF');

% shelf time series
nexttile; hold on;
for i=1:length(s.z)
    plot(s.t,s.Ss(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('salinity');
h = colorbar; colormap(gca,cmapz); caxis([-p.H,0]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on; title('SHELF');

%% volume fluxes
figure();

tcl = tiledlayout(2,2);
title(tcl,'VOLUME FLUXES');

% plume
nexttile; hold on;
if size(s.QVp,1)==1
    for i=1:length(ip)
        plot(squeeze(s.QVp(1,:,ip(i))),s.z,'color',cmapt(i,:),'linewidth',lw);
    end
else
    for i=1:length(ip)
        plot(squeeze(sum(s.QVp(:,:,ip(i)))),s.z,'color',cmapt(i,:),'linewidth',lw);
    end
end
xlabel('volume flux (m$^3$/s)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)'); ylim([-p.H,0]);
set(gca,'box','on'); grid on; title('PLUME');

% shelf
nexttile; hold on;
for i=1:length(ip)
    plot(s.QVs(:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('volume flux (m$^3$/s)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)'); ylim([-p.H,0]);
set(gca,'box','on'); grid on; title('SHELF');

% icebergs
% note for these vertical fluxes we plot the layer interface fluxes
% Q_{j+1,j} rather than the net fluxes Q_j
ints = -[0;cumsum(s.H(1:end-1))];
nexttile; hold on;
for i=1:length(ip)
    Qint = -cumsum(s.QVi(:,ip(i)),'reverse');
    plot(Qint,ints,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('volume flux (m$^3$/s)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)'); ylim([-p.H,0]);
set(gca,'box','on'); grid on; title('ICEBERGS');

% vertical advection
% note for these vertical fluxes we plot the layer interface fluxes
% Q_{j+1,j} rather than the net fluxes Q_j
ints = -[0;cumsum(s.H(1:end-1))];
nexttile; hold on;
for i=1:length(ip)
    Qint = -cumsum(s.QVv(:,ip(i)),'reverse');
    plot(Qint,ints,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('volume flux (m$^3$/s)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)'); ylim([-p.H,0]);
set(gca,'box','on'); grid on; title('VERTICAL ADVECTIVE');

%% freshwater
figure();

tcl = tiledlayout(2,3);
title(tcl,'FRESHWATER MELT');

% iceberg melt rate
nexttile; hold on;
for i=1:length(ip)
    plot(s.mi(:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('iceberg melt rate (m/d)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

% iceberg melt flux
nexttile; hold on;
for i=1:length(ip)
    plot(s.QMi(:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('iceberg melt flux (m$^3$/s)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

% total iceberg melt flux
nexttile; hold on;
plot(s.t,sum(s.QMi),'k','linewidth',lw);
xlabel('time (days)'); ylabel('total iceberg melt flux (m$^3$/s)');
set(gca,'box','on'); grid on;

% plume melt rate
nexttile; hold on;
for i=1:length(ip)
    if size(s.QVp,1)==1
        plot(s.mp(1,:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
    else
        plot(mean(s.mp(:,:,ip(i))),s.z,'color',cmapt(i,:),'linewidth',lw);
    end
end
xlabel('plume melt rate (m/d)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

% plume melt flux
nexttile; hold on;
for i=1:length(ip)
    if size(s.QVp,1)==1
        plot(s.QMp(1,:,ip(i)),s.z,'color',cmapt(i,:),'linewidth',lw);
    else
        plot(sum(s.QMp(:,:,ip(i))),s.z,'color',cmapt(i,:),'linewidth',lw);
    end
end
xlabel('plume melt flux (m$^3$/s)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

% total plume melt flux
nexttile; hold on;
if size(s.QVp,1)==1
    plot(s.t,squeeze(sum(s.QMp)),'k','linewidth',lw);
else
    plot(s.t,squeeze(sum(sum(s.QMp))),'k','linewidth',lw);
end
xlabel('time (days)'); ylabel('total plume melt flux (m$^3$/s)');
set(gca,'box','on'); grid on;

end