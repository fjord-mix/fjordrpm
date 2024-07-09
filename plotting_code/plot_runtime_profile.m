function plot_runtime_profile(i, p, t, f, s)

% Preallocate variables
H = s.H; S = s.S; Se = s.Se; QVs = s.QVs; z = f.zs;
Sfjord = zeros(size(H,1));
Sext = zeros(size(H,1));
Uext = zeros(size(H,1));

% Find the box and shelf salinity and fjord-shelf velocity in each layer.
ints = -[0;cumsum(H(:,i))];
for k=1:size(H,1)
    inds=find(z<=ints(k)&z>=ints(k+1));
    Sfjord(inds)=S(k,i);
    Sext(inds)=Se(k,i);
    Uext(inds)=QVs(k,i)/(p.W*H(k,i));
end

Sfjord(1) = S(end,i);
Sext(1) = Se(end,i);
Uext(1) = QVs(end,i)/(p.W*H(end,i));

%  Plot the shelf salinity and the shelf/box salinity in each layer.
subplot(1,2,1);
plot(f.Ss(:,i),z,Sext,z,Sfjord,z);
set(gca,'box','on'); grid on;
xlabel('salinity');
ylabel('depth (m)');
legend('true shelf','shelf in box','fjord','location','southwest');
drawnow;

% Plot the fjord-shelf velocity in each layer.
subplot(1,2,2);
plot(Uext,z);
set(gca,'box','on'); grid on;
xlabel('exchange velocity (m/s)');
ylabel('depth (m)');
text(gca,0.05,0.05,['t = ',num2str(0.01*round(100*t(i+1))),' days'],'fontsize',12,'VerticalAlignment','bottom','units','normalized');
drawnow;

end