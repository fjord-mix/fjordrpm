function animate(outputfile,nframes)

% load data
load([outputfile,'.mat']);

% delete existing video if exists
warning off; delete([outputfile,'.mp4']); warning on;

% colour maps from https://uk.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
cmap_vel = cmocean('balance');
cmap_sal = cmocean('haline');
cmap_temp = cmocean('thermal');
cmap_anom = cmocean('delta');

% colour limits
u0 = [s.QVs./(p.W*s.H)];
ulims = max(abs(u0(:)))*[-1,1];
w0 = cumsum(s.QVp+s.QVs)/(p.W*p.L);
wlims = max(abs(w0(:)))*[-1,1];
Slims = [min(s.S(:)),max(s.S(:))];
Tlims = [min(s.T(:)),max(s.T(:))];
Sanomlims = max(abs(s.S(:)-s.Ss(:)))*[-1,1];
Tanomlims = max(abs(s.T(:)-s.Ts(:)))*[-1,1];

% grids
nx = 15;
x = linspace(0,p.L,nx)/1e3;
xg = linspace(0,p.L,nx+1)/1e3;
nz = length(s.z);
zc = s.z;
zg = -[0;cumsum(s.H)];

% sill patch
if p.sill==1
    xsill = [p.L,p.L,0.98*p.L,0.98*p.L]/1e3;
    zsill = -[p.H,p.Hsill,p.Hsill,p.H];
else
    xsill = NaN;
    zsill = NaN;
end
scolor = 0.5*[1,1,1];
swidth = 1;

% time steps to plot
inx = round(linspace(1,length(s.t),nframes));

% loop over time steps to plot
for k=1:length(inx)

    i = inx(k);
    figure(); set(gcf,'Visible','off');
    
    % plotting data for this timestep
    u = interp1([0,p.L]/1e3,[s.QVp(:,i),-s.QVs(:,i)]',x)'./(p.W*s.H);
    w = -[0;cumsum(s.QVp(:,i)+s.QVs(:,i))/(p.W*p.L)];
    w = 0.5*(w(1:end-1)+w(2:end));
    w = interp1([0,p.L]/1e3,[w,w]',x)';
    S = interp1([0,p.L]/1e3,[s.S(:,i),s.S(:,i)]',x)';
    Ss = interp1([0,p.L]/1e3,[s.Ss(:,i),s.Ss(:,i)]',x)';
    T = interp1([0,p.L]/1e3,[s.T(:,i),s.T(:,i)]',x)';
    Ts = interp1([0,p.L]/1e3,[s.Ts(:,i),s.Ts(:,i)]',x)';
    
    % pad data to deal with pcolor
    uplot = [[u,NaN(nz,1)];NaN(1,nx+1)];
    wplot = [[w,NaN(nz,1)];NaN(1,nx+1)];
    Splot = [[S,NaN(nz,1)];NaN(1,nx+1)];
    Sanomplot = [[S-Ss,NaN(nz,1)];NaN(1,nx+1)];
    Tplot = [[T,NaN(nz,1)];NaN(1,nx+1)];
    Tanomplot = [[T-Ts,NaN(nz,1)];NaN(1,nx+1)];
    
    % plot along-fjord velocity
    subplot(2,3,1); hold on;
    pcolor(xg,zg,uplot); shading flat;
    patch(xsill,zsill,scolor,'edgecolor','k','linewidth',swidth);
    colorbar; colormap(gca,cmap_vel); caxis(ulims);
    xlim([0,p.L]/1e3); ylim([-p.H,0]);
    title('along-fjord velocity (m/s)');
    xlabel('x (km)'); ylabel('depth (m)');
    set(gca,'box','on','layer','top');
    
    % plot vertical velocity
    subplot(2,3,4); hold on;
    pcolor(xg,zg,wplot); shading flat;
    patch(xsill,zsill,scolor,'edgecolor','k','linewidth',swidth);
    colorbar; colormap(gca,cmap_vel); caxis(wlims);
    xlim([0,p.L]/1e3); ylim([-p.H,0]);
    title('vertical velocity (m/s)');
    xlabel('x (km)'); ylabel('depth (m)');
    set(gca,'box','on','layer','top');
    
    % plot salinity
    subplot(2,3,2); hold on;
    pcolor(xg,zg,Splot); shading flat;
    patch(xsill,zsill,scolor,'edgecolor','k','linewidth',swidth);
    colorbar; colormap(gca,cmap_sal); caxis(Slims);
    xlim([0,p.L]/1e3); ylim([-p.H,0]);
    title('fjord salinity');
    xlabel('x (km)'); ylabel('depth (m)');
    set(gca,'box','on','layer','top');
    
    % plot salinity anomaly
    subplot(2,3,5); hold on;
    pcolor(xg,zg,Sanomplot); shading flat;
    patch(xsill,zsill,scolor,'edgecolor','k','linewidth',swidth);
    colorbar; colormap(gca,cmap_anom); caxis(Sanomlims);
    xlim([0,p.L]/1e3); ylim([-p.H,0]);
    title('fjord-shelf salinity anomaly');
    xlabel('x (km)'); ylabel('depth (m)');
    set(gca,'box','on','layer','top');
    
    % plot temperature
    subplot(2,3,3); hold on;
    pcolor(xg,zg,Tplot); shading flat;
    patch(xsill,zsill,scolor,'edgecolor','k','linewidth',swidth);
    colorbar; colormap(gca,cmap_temp); caxis(Tlims);
    xlim([0,p.L]/1e3); ylim([-p.H,0]);
    title('fjord temperature (C)');
    xlabel('x (km)'); ylabel('depth (m)');
    set(gca,'box','on','layer','top');
    
    % plot temperature anomaly
    subplot(2,3,6); hold on;
    pcolor(xg,zg,Tanomplot); shading flat;
    patch(xsill,zsill,scolor,'edgecolor','k','linewidth',swidth);
    colorbar; colormap(gca,cmap_anom); caxis(Tanomlims);
    xlim([0,p.L]/1e3); ylim([-p.H,0]);
    title('fjord-shelf temperature anomaly (C)');
    xlabel('x (km)'); ylabel('depth (m)');
    set(gca,'box','on','layer','top');
    
    % time step
    sgtitle(['t = ',num2str(0.01*round(100*s.t(i))),' days']);
    
    % save plot
    if k<10
        savenum = ['00', num2str(k)];
    elseif k<100
        savenum = ['0', num2str(k)];
    else
        savenum = num2str(k);
    end
    saveplot(40,15,300,[outputfile,'_',savenum,'.png']);
    close all;

end

% write video
video = VideoWriter([outputfile,'.mp4'],'MPEG-4');
video.FrameRate = 5;
open(video);
pngs = dir(fullfile([outputfile,'*.png']));
for k = 1:length(pngs)
    I = imread(pngs(k).name);
    writeVideo(video, I);
end
close(video);
delete([outputfile,'*.png']);

end
