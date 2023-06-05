% function to animate box model results - updated to show more variables as well
function [] = animate_v4p1(fjord_model,outputdir,outputfile,nframes)

% clearvars -except fjord_model outputfile nframes; close all;
% load(outputfile,'fjord_model');

names = fieldnames(fjord_model);
for i=1:length(names)
    eval([names{i} '=fjord_model.' names{i},';']);
end

orig_dir = pwd;
eval(['cd ',outputdir,'/plots/']);

Slims = [min([min(s.S(:)),min(f.Ss(:))]), max([max(s.S(:)),max(f.Ss(:))])];
Tlims = [min([min(s.T(:)),min(f.Ts(:))]), max([max(s.T(:)),max(f.Ts(:))])];
zlims = [-p.H,0];
x0 = 1500;
sf = x0*0.2/max([max(abs(s.QVs(:))),max(abs(s.QVg(:)))]);

fs = 10;
lspace = 0.08;
hspace1 = 0.03;
hspace2 = 0.08;
rspace = 0.02;
cby = 0.075;
cbh = 0.02;
cbw = 0.2;
pw1 = 0.05;
pw3 = 0.25;
pw2 = 1-lspace-rspace-hspace1-hspace2-8*pw1-pw3; %0.15;
pw2b = 1-lspace-rspace-hspace1-hspace2-8*pw1-pw3;
bspace = 0.12;
tspace = 0.06;
vspace = 0.05;
ph = 1-bspace-tspace;
ph3 = (1-bspace-tspace-vspace)/3;
fs2 = 6;

f_qsg = interp1(linspace(0,max(s.t),length(f.Qsg)),f.Qsg,s.t,'linear');

for i=1:round((length(s.t)-1)/nframes):length(s.t)-1

    figure(); set(gcf,'Visible','off');

    ints = [0;cumsum(s.H(:,i))];
    y = -0.5*(ints(1:end-1)+ints(2:end));
    y2 = cumsum(s.H(1:3,i));

    %% Time series
    a3 = axes('position',[lspace+2*pw2+2*hspace1+hspace2+2*pw1,bspace,pw3,ph3]); hold on;
    plot(s.t,s.phi,'linewidth',1);
    plot(s.t(i),s.phi(:,i),'k.');
    xlabel('t (days)','fontsize',fs2);
    ylabel('\phi (shelf-fjord pressure difference)','fontsize',fs2);
    set(gca,'box','on','fontsize',fs2);
    
    a4 = axes('position',[lspace+2*pw2+2*hspace1+hspace2+2*pw1,bspace+ph3+vspace,pw3,ph3]); hold on;
    plot(s.t,s.H,'linewidth',1);
    plot(s.t(i),s.H(:,i),'k.');
    ylabel('H (m)','fontsize',fs2);
    set(gca,'box','on','fontsize',fs2);
    legend({'1','2','3','4'},'location','north','orientation','horizontal','fontsize',4,'NumColumns',2);
    
    a5 = axes('position',[lspace+2*pw2+2*hspace1+hspace2+2*pw1,bspace+2.2*ph3+vspace,pw3,ph3]); hold on;
    plot(s.t,f_qsg,'linewidth',1);
    plot(s.t(i),f_qsg(i),'k.');    
    ylabel('Subgl. discharge (m^3s^{-1})','fontsize',fs2);
    set(gca,'box','on','fontsize',fs2);


    %% box model plots

    a1 = axes('position',[lspace,bspace,pw2,ph]); hold on;
    text(0.01*x0,-800,['t = ',num2str(0.01*round(100*s.t(i))),' days'],...
        'fontsize',fs,'VerticalAlignment','bottom');
    
    title('Salinity: fjord','fontsize',fs);
    pcolor([0,x0],-[0;cumsum(s.H(:,i))],[[s.S(:,i);0],[s.S(:,i);0]]);
    set(gca,'xtick',[],'xticklabel',{});
    q1 = quiver(0*y+x0,y,sf*s.QVs(:,i),0*s.QVs(:,i),'autoscale','off');
    set(q1,'color','k','linewidth',2);
    q2 = quiver(0*y,y,sf*s.QVg(:,i),0*s.QVg(:,i),'autoscale','off');
    set(q2,'color','k','linewidth',2);
    q3 = quiver(0*y2+x0/2,-y2,0*y2,sf*[s.QVk(1,i),s.QVk(1,i)+s.QVk(2,i),0]','autoscale','off');
    set(q3,'color','k','linewidth',2);
    q4 = quiver(0*y2+x0/3,-y2,0*y2,sf*[s.QVb(1,i),s.QVb(1,i)+s.QVb(2,i),s.QVb(1,i)+s.QVb(2,i)+s.QVb(3,i)]','autoscale','off');
    set(q4,'color','r','linewidth',2);
    set(a1,'clipping','off','box','on','fontsize',fs);
    clim(a1,Slims);
    xlim([0,x0]); ylim(zlims);

    a2 = axes('position',[lspace+pw2+hspace1,bspace,pw1,ph]);    
    pcolor([0,x0],f.zs,[f.Ss(:,i),f.Ss(:,i)]); shading flat;    
    xlim([0,x0]); ylim(zlims);
    set(gca,'box','on','fontsize',fs);
    set(gca,'xtick',[],'ytick',zlims,'yticklabel',{});    
    title('shelf','fontsize',fs);
    clim(a2,Slims);
    
    h = colorbar('southoutside');
    set(h,'position',[0.2-cbw/2,cby,cbw,cbh],'fontsize',fs);
    % xlabel(h,'salinity','fontsize',fs);
    clim(Slims);    
    xlim([0,x0]); ylim(zlims);    
    set(gca,'xtick',[],'ytick',zlims);    


    a1b = axes('position',[lspace+pw1+pw2+2*hspace1,bspace,pw2b,ph]); hold on;
    title('Temp.: fjord','fontsize',fs);
    pcolor([0,x0],-[0;cumsum(s.H(:,i))],[[s.T(:,i);0],[s.T(:,i);0]]);
    set(gca,'xtick',[],'ytick',zlims,'yticklabel',{});
    q1 = quiver(0*y+x0,y,sf*s.QVs(:,i),0*s.QVs(:,i),'autoscale','off');
    set(q1,'color','k','linewidth',2);
    q2 = quiver(0*y,y,sf*s.QVg(:,i),0*s.QVg(:,i),'autoscale','off');
    set(q2,'color','k','linewidth',2);
    q3 = quiver(0*y2+x0/2,-y2,0*y2,sf*[s.QVk(1,i),s.QVk(1,i)+s.QVk(2,i),0]','autoscale','off');
    set(q3,'color','k','linewidth',2);
    q4 = quiver(0*y2+x0/3,-y2,0*y2,sf*[s.QVb(1,i),s.QVb(1,i)+s.QVb(2,i),s.QVb(1,i)+s.QVb(2,i)+s.QVb(3,i)]','autoscale','off');
    set(q4,'color','r','linewidth',2);
    set(a1b,'clipping','off','box','on','fontsize',fs);
    clim(a1b,Tlims);
    xlim([0,x0]); ylim(zlims);

    a2b = axes('position',[lspace+pw1+pw2+2*hspace1+pw2+hspace1,bspace,pw1,ph]);    
    pcolor([0,x0],f.zs,[f.Ts(:,i),f.Ts(:,i)]); shading flat;    
    xlim([0,x0]); ylim(zlims);
    set(gca,'box','on','fontsize',fs);
    set(gca,'xtick',[],'ytick',zlims,'yticklabel',{});        
    title('shelf','fontsize',fs);
    clim(a2b,Tlims);

    h2 = colorbar('southoutside');
    set(h2,'position',[0.46-cbw/2,cby,cbw,cbh],'fontsize',fs);
    % xlabel(h2,'Temperature (^oC)','fontsize',fs);
    clim(Tlims);    
    xlim([0,x0]); ylim(zlims);    
    set(gca,'xtick',[],'ytick',zlims);            

    if i<10, savenum = ['00',num2str(i)];
    elseif i<100, savenum = ['0',num2str(i)];
    else 
        savenum = num2str(i);
    end

    saveplot(25,10,300,[outputfile,'_',savenum,'.png']);
    close all;

end

% write video
video = VideoWriter([outputfile,'.mp4'],'MPEG-4');
video.FrameRate = 10;
open(video);
for i=1:round((length(s.t)-1)/nframes):length(s.t)-1

    if i<10, savenum = ['00',num2str(i)];
    elseif i<100, savenum = ['0',num2str(i)];
    else 
        savenum = num2str(i);
    end

    I = imread([outputfile,'_',savenum,'.png']);
    writeVideo(video,I);

end
close(video);
delete([outputfile,'*.png'])
eval(['cd ',orig_dir]); 

end
