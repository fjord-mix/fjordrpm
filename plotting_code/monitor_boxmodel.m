function hf_track = monitor_boxmodel(hf_track,i,H,T,S,f)

if ~isempty(hf_track)
    figure(hf_track);    
    subplot(1,3,1); hold on
    title(['time step ',num2str(i+1)])
    t_plot = (i+1).*ones(size(H(:,1))); H_plot = -cumsum(H(:,i+1));
    scatter(t_plot,H_plot,30,'black','filled'); box on
    ylim([-1e3 0]); 
    
    subplot(1,3,2); hold on
    plot(T(:,i+1),H_plot); 
    plot(f.Ts(:,i+1),f.zs); 
    xlim([-2 5]); ylim([-1e3 0]); 
    
    subplot(1,3,3); hold on
    plot(S(:,i+1),H_plot);
    plot(f.Ss(:,i+1),f.zs);
    xlim([0 35]); ylim([-1e3 0]); 
    
    legend('model','shelf','location','southwest')
else
    hf_track = figure();
    subplot(1,3,1);
    scatter(1,-cumsum(H),30,'black','filled');
    xlabel('time step'); ylabel('depth (m)')

    subplot(1,3,2); hold on
    plot(f.Ts(:,1),f.zs)
    plot(T,-cumsum(H))
    xlabel('Temperature (^oC)'); 
    
    subplot(1,3,3); hold on
    plot(S,-cumsum(H))
    plot(f.Ss(:,1),f.zs)
    xlabel('Salinity')
    legend('model','shelf','location','southwest')
end


end