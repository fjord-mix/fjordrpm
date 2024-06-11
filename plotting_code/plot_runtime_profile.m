function plot_runtime_profile(i,t,f,p,H,S,T,Se,QVs)


s_bnds = [min(f.Ss(:)) max(f.Ss(:))+0.1];

   ints = -[0;cumsum(H(:,i))];
        z = f.zs;
        Sshelf = f.Ss(:,i);
        for k=1:size(H,1)
            inds=find(z<=ints(k)&z>=ints(k+1));
            Sfjord(inds)=S(k,i);
            Sext(inds)=Se(k,i);
            Uext(inds)=QVs(k,i)/(p.W*H(k,i));
        end
        Sfjord(1) = S(end,i);
        Sext(1) = Se(end,i);
        Uext(1) = QVs(end,i)/(p.W*H(end,i));
        % hf_track = monitor_boxmodel(hf_track,i,H,T,S,f);
        % hf_track = show_boxmodel([],i,t,H,T,S,QVs,QVg,QVk,QVb,f);
%         plot_debug_profile(i,t,f,p,H,S,[]);
        subplot(1,2,1);
        plot(f.Ss(:,i),z,Sext,z,Sfjord,z); 
          set(gca,'box','on'); grid on;
        xlabel('salinity');
        ylabel('depth (m)');
        legend('true shelf','shelf in box','fjord','location','southwest');       
        drawnow;
        subplot(1,2,2);
        plot(Uext,z); 
        set(gca,'box','on'); grid on;
        xlabel('exchange velocity (m/s)');
        ylabel('depth (m)');
        text(gca,0.05,0.05,['t = ',num2str(0.01*round(100*t(i+1))),' days'],'fontsize',12,'VerticalAlignment','bottom','units','normalized');
        drawnow;
end