function plot_debug_profile(i,t,f,p,H,S,s_bnds)
    % this is to monitor the how the fjord profile fares compared with the shelf
    % mainly for debugging purposes, which is why it's not a refined plot
    
    ints = [0;-cumsum(H(:,i+1))]; z_levels=0.5*(ints(2:end)+ints(1:end-1));
    figure(99); 
    plot(f.Ss(:,i+1),f.zs); hold on; ylim([-p.H 0]); xlim(s_bnds);% to monitor interpolation results
    scatter(S(:,i+1),z_levels); yline(-cumsum(H(:,i+1)),':k','linewidth',0.5); % to monitor interpolation results  
    text(gca,0.05,0.05,['t = ',num2str(0.01*round(100*t(i+1))),' days'],'fontsize',12,'VerticalAlignment','bottom','units','normalized');
    hold off; 
end