function plot_outputs(fjord_model1,plot_forcing)

if ~isfield(fjord_model1,'m'), fjord_model1.m.name='Unamed'; end
if ~isfield(fjord_model1.m,'name'), fjord_model1.m.name='Unamed'; end
model_runtime1 = fjord_model1.s.t(1:size(fjord_model1.s.H,2));
if isfield(fjord_model1.m,'time_axis')
    runtime_axis = fjord_model1.m.time_axis;
    t0 = convertTo(runtime_axis(1),'datenum');
    taxis1 = NaT([size(fjord_model1.s.H,2),1]);
    isdays='';
    for i_time=1:length(taxis1)
        taxis1(i_time) = datetime(t0+model_runtime1(i_time),'ConvertFrom','datenum');
    end
else
    taxis1=fjord_model1.s.t;
    isdays = ' (days)';
end

layer_lbls =cell(size(fjord_model1.a.H0));
for i=1:fjord_model1.p.N+fjord_model1.p.sill
    layer_lbls{i}=num2str(i);
end

Hlims = [0, max(fjord_model1.s.H(:))+100];
% Tlims = [min(fjord_model1.s.T(:))-0.5, max(fjord_model1.s.T(:))+0.5];
% Slims = [min(fjord_model1.s.S(:))-0.5, max(fjord_model1.s.S(:))+0.5];
Tlims = [-0.5 6.5];
Slims = [27 36];

m=3; n=1;
figure('Position',[20 20 800 600],'Name',fjord_model1.m.name); 
subplot(m,n,1); plot(taxis1,fjord_model1.s.H,'linewidth',1.5); hold on; ylim(Hlims)
ylabel('Thickness (m)'); 
subplot(m,n,2); plot(taxis1,fjord_model1.s.T,'linewidth',1.5); hold on; ylim(Tlims)
ylabel('Temperature ($^o$C)','interpreter','latex'); 
subplot(m,n,3); plot(taxis1,fjord_model1.s.S,'linewidth',1.5); hold on; ylim(Slims)
ylabel('Salinity'); xlabel(['Model time',isdays])
hl = legend(layer_lbls,'Location','southeast'); title(hl,'Layer'); hl.NumColumns=2;
if nargin > 1 && plot_forcing
    fjord_temp=NaN(size(fjord_model1.s.T));
    fjord_salt=NaN(size(fjord_model1.s.S));
    for i_time=1:length(taxis1)
        [fjord_temp(:,i_time), fjord_salt(:,i_time)] = bin_ocean_profiles(fjord_model1.f.Ts(:,i_time),fjord_model1.f.Ss(:,i_time),-fjord_model1.f.zs,fjord_model1.s.H(:,i_time)',fjord_model1.p);
    end
    subplot(m,n,2); plot(taxis1,fjord_temp,'linewidth',1.5,'linestyle',':')
    subplot(m,n,3); plot(taxis1,fjord_salt,'linewidth',1.5,'linestyle',':')
end

end