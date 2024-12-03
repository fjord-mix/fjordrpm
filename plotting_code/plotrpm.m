function plotrpm(outputfile)

% a script to demonstrate some possible plots from model output

% load data
load([outputfile,'.mat']);

% colourmaps
cmapt = parula(length(s.t));
cmapz = parula(length(s.z));
lw = 1;

%% temperature
figure();
subplot(1,2,1); hold on;
for i=1:length(s.t)
    plot(s.T(:,i),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('temperature (C)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

subplot(1,2,2); hold on;
for i=1:length(s.z)
    plot(s.t,s.T(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('temperature (C)');
h = colorbar; colormap(cmapz); caxis([min(s.z),max(s.z)]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on;

%% salinity
figure();
subplot(1,2,1); hold on;
for i=1:length(s.t)
    plot(s.S(:,i),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('salinity'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

subplot(1,2,2); hold on;
for i=1:length(s.z)
    plot(s.t,s.S(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('salinity');
h = colorbar; colormap(cmapz); caxis([min(s.z),max(s.z)]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on;

%% iceberg melt flux
figure();
subplot(1,2,1); hold on;
for i=1:length(s.t)
    plot(s.QMi(:,i),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('iceberg melt flux (m$^3$/s)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

subplot(1,2,2); hold on;
for i=1:length(s.z)
    plot(s.t,s.QMi(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('iceberg melt flux (m$^3$/s)');
h = colorbar; colormap(cmapz); caxis([min(s.z),max(s.z)]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on;

%% iceberg melt rate
figure();
subplot(1,2,1); hold on;
for i=1:length(s.t)
    plot(s.icebergmeltrate(:,i),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('iceberg melt rate (m/d)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

subplot(1,2,2); hold on;
for i=1:length(s.z)
    plot(s.t,s.icebergmeltrate(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('iceberg melt rate (m/d)');
h = colorbar; colormap(cmapz); caxis([min(s.z),max(s.z)]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on;

%% shelf exchange
figure();
subplot(1,2,1); hold on;
for i=1:length(s.t)
    plot(s.QVs(:,i),s.z,'color',cmapt(i,:),'linewidth',lw);
end
xlabel('fjord-shelf exchange (m$^3$/s)'); ylabel('depth (m)');
h = colorbar; colormap(cmapt); caxis([min(s.t),max(s.t)]);
ylabel(h,'time (days)');
set(gca,'box','on'); grid on;

subplot(1,2,2); hold on;
for i=1:length(s.z)
    plot(s.t,s.QVs(i,:),'color',cmapz(length(s.z)-i+1,:),'linewidth',lw);
end
xlabel('time (days)'); ylabel('fjord-shelf exchange (m$^3$/s)');
h = colorbar; colormap(cmapz); caxis([min(s.z),max(s.z)]);
ylabel(h,'depth (m)');
set(gca,'box','on'); grid on;

end