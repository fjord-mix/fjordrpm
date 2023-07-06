function hf = show_box_model(hf,i,t,H,T,S,QVs,QVg,QVk,QVb,f)

first_time=0;
if size(H,2) > 1
    ints = [0; cumsum(H(:,i))];
    y2 = cumsum(H(1:3,i));
else
    ints = [0;cumsum(H,1)];
    y2 = cumsum(H(1:3));
end
y = -0.5*(ints(1:end-1)+ints(2:end));


Slims = [min([min(S(:)),min(f.Ss(:))]),...
         max([max(S(:)),max(f.Ss(:))])];
zlims = [-ints(end),0];
x0 = 1500;
if ~(isempty(QVs) || isempty(QVg) || isempty(QVk) || isempty(QVb))
    sf = x0*0.2/max([max(abs(QVs(:))),max(abs(QVg(:)))]);
end

fs = 12;
lspace = 0.12; % 0.08
hspace1 = 0.05;
hspace2 = 0.1;
rspace = 0.02;
cby = 0.13;
cbh = 0.02;
cbw = 0.4;
pw1 = 0.25;
pw3 = 0.;
pw2 = 1-lspace-rspace-hspace1-hspace2-pw1-pw3;
bspace = 0.19;
tspace = 0.06;
vspace = 0.05;
ph = 1-bspace-tspace;

if isempty(hf)   
    hf = figure(); %set(hf,'Visible','off');
    first_time=1;
else
    figure(hf);    
end


a2 = axes('position',[lspace+pw2+hspace1,bspace,pw1,ph]);
pcolor(a2,[0,x0],f.zs,[f.Ss(:,i),f.Ss(:,i)]); shading flat;
clim(a2,Slims);
xlim(a2,[0,x0]); ylim(a2,zlims);
set(a2,'xtick',[],'ytick',zlims,'yticklabel',{});
set(a2,'box','on','fontsize',fs);
title(a2,'shelf','fontsize',fs);

a1 = axes('position',[lspace,bspace,pw2,ph]); hold on;
pcolor(a1,[0,x0],-ints,[[S(:,i);0],[S(:,i);0]]); 
if ~(isempty(QVs) || isempty(QVg) || isempty(QVk) || isempty(QVb))
    q1 = quiver(a1,0*y+x0,y,sf*QVs(:,i),0*QVs(:,i),'autoscale','off');
    set(q1,'color','k','linewidth',2);    
    q2 = quiver(a1,0*y,y,sf*QVg(:,i),0*QVg(:,i),'autoscale','off');
    set(q2,'color','k','linewidth',2);
    q3 = quiver(a1,0*y2+x0/2,-y2,0*y2,sf*[QVk(1,i),QVk(1,i)+QVk(2,i),0]','autoscale','off');
    set(q3,'color','k','linewidth',2);
    q4 = quiver(a1,0*y2+x0/3,-y2,0*y2,sf*[QVb(1,i),QVb(1,i)+QVb(2,i),QVb(1,i)+QVb(2,i)+QVb(3,i)]','autoscale','off');
    set(q4,'color','r','linewidth',2);
end

if first_time
    h = colorbar('southoutside');
    set(h,'position',[0.5-cbw/2,cby,cbw,cbh],'fontsize',fs);
    xlabel(h,'salinity','fontsize',fs);
    clim(Slims);
end
title(a1,'fjord','fontsize',fs);
xlim(a1,[0,x0]); ylim(a1,zlims);
set(a1,'clipping','off','box','on','fontsize',fs);
set(a1,'xtick',[],'ytick',zlims);
text(a1,0.01*x0,-800,['t = ',num2str(0.01*round(100*t(i))),' days'],...
    'fontsize',fs,'VerticalAlignment','bottom');
hold off


end