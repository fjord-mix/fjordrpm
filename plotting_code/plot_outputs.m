function plot_outputs(fjord_model1,fjord_model2)

if ~isfield(fjord_model1.m,'name'), fjord_model1.m.name='Unamed'; end
model_runtime1 = fjord_model1.s.t(1:size(fjord_model1.s.H,2));
runtime_axis = fjord_model1.m.time_axis;
t0 = convertTo(runtime_axis(1),'datenum');
taxis1 = NaT([size(fjord_model1.s.H,2),1]);
for i_time=1:length(taxis1)
    taxis1(i_time) = datetime(t0+model_runtime1(i_time),'ConvertFrom','datenum');
end

Hlims = [0, max(fjord_model1.s.H(:))+100];
% Tlims = [min(fjord_model1.s.T(:))-0.5, max(fjord_model1.s.T(:))+0.5];
% Slims = [min(fjord_model1.s.S(:))-0.5, max(fjord_model1.s.S(:))+0.5];
Tlims = [-0.5 6.5];
Slims = [27 36];

m=3; n=1;
figure('Position',[20 20 800 600],'Name',fjord_model1.m.name); 
subplot(m,n,1); plot(taxis1,fjord_model1.s.H,'linewidth',1.5); hold on; ylim(Hlims)
ylabel('Thickness (m)'); legend('1','2','3','4','Location','east');
subplot(m,n,2); plot(taxis1,fjord_model1.s.T,'linewidth',1.5); hold on; ylim(Tlims)
ylabel('Temperature (^oC)'); 
subplot(m,n,3); plot(taxis1,fjord_model1.s.S,'linewidth',1.5); hold on; ylim(Slims)
ylabel('Salinity'); xlabel('Model time (days)')

if nargin > 1
    model_runtime2 = fjord_model2.s.t(1:size(fjord_model1.s.H,2));
    runtime_axis = fjord_model2.m.time_axis;
    t0 = convertTo(runtime_axis(1),'datenum');
    taxis2 = NaT([size(fjord_model2.s.H,2),1]);
    for i_time=1:length(taxis2)
        taxis2(i_time) = datetime(t0+model_runtime2(i_time),'ConvertFrom','datenum');
    end

    subplot(m,n,1); plot(taxis2,fjord_model2.s.H,'linewidth',1.5,'LineStyle','--'); 
    subplot(m,n,2); plot(taxis2,fjord_model2.s.T,'linewidth',1.5,'LineStyle','--');
    subplot(m,n,3); plot(taxis2,fjord_model2.s.S,'linewidth',1.5,'LineStyle','--');
end

end