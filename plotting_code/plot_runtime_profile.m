function plot_runtime_profile(i, p, t, s)

% PLOT_RUNTIME_PROFILE makes basic plots as the model is solving.
%   PLOT_RUNTIME_PROFILE(i, p, t, s) makes plots of temperature, salinity
%   and fjord-shelf exchange as the model is running. Inputs are the time
%   step i, parameters structure p, time vector t and solution structure s.
%   Called from run_model if p.plot_runtime = 1.

% layer interfaces
ints = -[0;cumsum(s.H)];

% plot temperatures
subplot(1,3,1); cla; hold on;
stairs([s.Ts(1,i);s.Ts(:,i)],ints,'linewidth',2);
stairs([s.T(1,i);s.T(:,i)],ints,'linewidth',2);
set(gca,'box','on'); grid on; ylim([-p.H,0]);
xlabel('temperature (C)');
ylabel('depth (m)');
legend('shelf','fjord','location','southwest');

% plot salinities
subplot(1,3,2); cla; hold on;
stairs([s.Ss(1,i);s.Ss(:,i)],ints,'linewidth',2);
stairs([s.S(1,i);s.S(:,i)],ints,'linewidth',2);
set(gca,'box','on'); grid on; ylim([-p.H,0]);
xlabel('salinity');
ylabel('depth (m)');
title(['t = ',num2str(0.01*round(100*t(i))),' days']);

% plot fjord-shelf exchange velocity
subplot(1,3,3);
stairs(-[s.QVs(1,i)/(p.W*s.H(1));s.QVs(:,i)./(p.W*s.H)],ints,'k','linewidth',2);
set(gca,'box','on'); grid on; ylim([-p.H,0]);
xlabel('fjord-shelf exchange velocity (m/s)');
ylabel('depth (m)');
title('positive values directed out of fjord');
drawnow;

end