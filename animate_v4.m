% function to animate box model results
function [] = animate_v4(outputfile,fname,nframes)

clearvars -except outputfile fname nframes; close all;
load(outputfile,'fjord_model');

names = fieldnames(fjord_model);
for i=1:length(names)
    eval([names{i} '=fjord_model.' names{i} ]);
end

orig_dir = pwd;
folder = dir(outputfile);
eval(['cd ',folder.folder,'/plots/']);

Slims = [min([min(s.S(:)),min(f.Ss(:))]),...
         max([max(s.S(:)),max(f.Ss(:))])];
zlims = [-p.H,0];
x0 = 1500;
sf = x0*0.2/max([max(abs(s.QVs(:))),max(abs(s.QVg(:)))]);

fs = 12;
lspace = 0.08;
hspace1 = 0.05;
hspace2 = 0.1;
rspace = 0.02;
cby = 0.13;
cbh = 0.02;
cbw = 0.4;
pw1 = 0.15;
pw3 = 0.25;
pw2 = 1-lspace-rspace-hspace1-hspace2-pw1-pw3;
bspace = 0.19;
tspace = 0.06;
vspace = 0.05;
ph = 1-bspace-tspace;
ph2 = (1-bspace-tspace-vspace)/2;
fs2 = 6;

for i=1:round((length(s.t)-1)/nframes):length(s.t)-1,

    figure(); set(gcf,'Visible','off');

    ints = [0;cumsum(s.H(:,i))];
    y = -0.5*(ints(1:end-1)+ints(2:end));
    y2 = cumsum(s.H(1:3,i));

    a3 = axes('position',[lspace+pw2+hspace1+hspace2+pw1,bspace+ph2+vspace,pw3,ph2]); hold on;
    plot(s.t,s.H,'linewidth',1);
    plot(s.t(i),s.H(:,i),'k.');
    ylabel('H (m)','fontsize',fs2);
    set(gca,'box','on','fontsize',fs2);
    legend({'1','2','3','4'},'location','north','orientation','horizontal','fontsize',4);

    a4 = axes('position',[lspace+pw2+hspace1+hspace2+pw1,bspace,pw3,ph2]); hold on;
    plot(s.t,s.phi,'linewidth',1);
    plot(s.t(i),s.phi(:,i),'k.');
    xlabel('t (days)','fontsize',fs2);
    ylabel('\phi (shelf-fjord pressure difference)','fontsize',fs2);
    set(gca,'box','on','fontsize',fs2);

    a2 = axes('position',[lspace+pw2+hspace1,bspace,pw1,ph]);
    pcolor([0,x0],f.zs,[f.Ss(:,i),f.Ss(:,i)]); shading flat;
    caxis(Slims);
    xlim([0,x0]); ylim(zlims);
    set(gca,'box','on','fontsize',fs);
    set(gca,'xtick',[],'ytick',zlims,'yticklabel',{});
    title('shelf','fontsize',fs);

    a1 = axes('position',[lspace,bspace,pw2,ph]); hold on;
    pcolor([0,x0],-[0;cumsum(s.H(:,i))],[[s.S(:,i);0],[s.S(:,i);0]]);
    q1 = quiver(0*y+x0,y,sf*s.QVs(:,i),0*s.QVs(:,i),'autoscale','off');
    set(q1,'color','k','linewidth',2);
    q2 = quiver(0*y,y,sf*s.QVg(:,i),0*s.QVg(:,i),'autoscale','off');
    set(q2,'color','k','linewidth',2);
    q3 = quiver(0*y2+x0/2,-y2,0*y2,sf*[s.QVk(1,i),s.QVk(1,i)+s.QVk(2,i),0]','autoscale','off');
    set(q3,'color','k','linewidth',2);
    q4 = quiver(0*y2+x0/3,-y2,0*y2,sf*[s.QVb(1,i),s.QVb(1,i)+s.QVb(2,i),s.QVb(1,i)+s.QVb(2,i)+s.QVb(3,i)]','autoscale','off');
    set(q4,'color','r','linewidth',2);

    h = colorbar('southoutside');
    set(h,'position',[0.5-cbw/2,cby,cbw,cbh],'fontsize',fs);
    xlabel(h,'salinity','fontsize',fs);
    caxis(Slims);
    title('fjord','fontsize',fs);
    xlim([0,x0]); ylim(zlims);
    set(gca,'clipping','off','box','on','fontsize',fs);
    set(gca,'xtick',[],'ytick',zlims);
    text(0.01*x0,-800,['t = ',num2str(0.01*round(100*s.t(i))),' days'],...
        'fontsize',fs,'VerticalAlignment','bottom');

    if i<10, savenum = ['00',num2str(i)];
    elseif i<100, savenum = ['0',num2str(i)];
    else savenum = num2str(i);
    end

    saveplot(25,10,300,[fname,'_',savenum,'.png']);
    close all;

end

% write video
video = VideoWriter([fname,'.mp4'],'MPEG-4');
video.FrameRate = 10;
open(video);
for i=1:round((length(s.t)-1)/nframes):length(s.t)-1,

    if i<10, savenum = ['00',num2str(i)];
    elseif i<100, savenum = ['0',num2str(i)];
    else savenum = num2str(i);
    end

    I = imread([fname,'_',savenum,'.png']);
    writeVideo(video,I);

end
close(video);
delete([fname,'*.png'])
eval(['cd ',orig_dir]); 

end