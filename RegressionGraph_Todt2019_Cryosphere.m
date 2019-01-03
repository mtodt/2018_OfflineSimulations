tic
clear
close all
%{
Multi-linear regression for all single-PFT Toy Model sites.
%}

%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
load('MLRdata_Alptal.mat','PAI_Alptal','LWR_Atmosphere_Alptal_2005',...
    'LWR_Subcanopy_Observations_Alptal_2005','LWR_Subcanopy_CLM_Alptal_2005',...
    'SWR_Incoming_Alptal_2005','Sky_Emissivity_Alptal_2005');
load('MLRdata_Cherskiy.mat','PAI_Cherskiy','LWR_Atmosphere_Cherskiy',...
    'LWR_Subcanopy_Observations_Cherskiy','LWR_Subcanopy_CLM_Cherskiy',...
    'SWR_Incoming_Cherskiy','Sky_Emissivity_Cherskiy');
load('MLRdata_Seehornwald.mat','PAI_Seehornwald',...
    'LWR_Atmosphere_Seehornwald_2009','LWR_Subcanopy_Observations_Seehornwald_2009',...
    'LWR_Subcanopy_CLM_Seehornwald_2009','SWR_Incoming_Seehornwald_2009',...
    'Sky_Emissivity_Seehornwald_2009');
load('MLRdata_Sodankyla.mat','PAI_Sodankyla','LWR_Atmosphere_Sodankyla',...
    'LWR_Subcanopy_Observations_Sodankyla','LWR_Subcanopy_CLM_Sodankyla',...
    'SWR_Incoming_Sodankyla','Sky_Emissivity_Sodankyla');


%-------------------------------------------------------------------------%
%-------------------------  calculation of bias  -------------------------%
%-------------------------------------------------------------------------%
% calculate LWR from vegetation using vegetation fraction from CLM

% Alptal
f_veg = 1 - exp(-PAI_Alptal);
LWR_Vegetation_CLM_Alptal_2005 = (LWR_Subcanopy_CLM_Alptal_2005 ...
    - (1-f_veg)*LWR_Atmosphere_Alptal_2005)/f_veg;
LWR_Vegetation_Obs_Alptal_2005 = (LWR_Subcanopy_Observations_Alptal_2005 ...
    - (1-f_veg)*LWR_Atmosphere_Alptal_2005)/f_veg;
LWR_Vegetation_Ratio_Alptal_2005 = LWR_Vegetation_CLM_Alptal_2005 ...
    ./ LWR_Vegetation_Obs_Alptal_2005;
% Cherskiy
f_veg = 1 - exp(-PAI_Cherskiy);
LWR_Vegetation_CLM_Cherskiy = (LWR_Subcanopy_CLM_Cherskiy ...
    - (1-f_veg)*LWR_Atmosphere_Cherskiy)/f_veg;
LWR_Vegetation_Obs_Cherskiy = (LWR_Subcanopy_Observations_Cherskiy ...
    - (1-f_veg)*LWR_Atmosphere_Cherskiy)/f_veg;
LWR_Vegetation_Ratio_Cherskiy ...
    = LWR_Vegetation_CLM_Cherskiy./LWR_Vegetation_Obs_Cherskiy;

% Seehornwald
f_veg = 1 - exp(-PAI_Seehornwald);
LWR_Vegetation_CLM_Seehornwald_2009 = (LWR_Subcanopy_CLM_Seehornwald_2009 ...
    - (1-f_veg)*LWR_Atmosphere_Seehornwald_2009)/f_veg;
LWR_Vegetation_Obs_Seehornwald_2009 = (LWR_Subcanopy_Observations_Seehornwald_2009 ...
    - (1-f_veg)*LWR_Atmosphere_Seehornwald_2009)/f_veg;
LWR_Vegetation_Ratio_Seehornwald_2009 = LWR_Vegetation_CLM_Seehornwald_2009...
    ./ LWR_Vegetation_Obs_Seehornwald_2009;

% Sodankyla
LWR_Subcanopy_CLM_Sodankyla_avg = nan(size(LWR_Atmosphere_Sodankyla));
LWR_Subcanopy_Observations_Sodankyla_avg = nan(size(LWR_Atmosphere_Sodankyla));
for l=1:length(LWR_Atmosphere_Sodankyla)
    LWR_Subcanopy_CLM_Sodankyla_avg(l) = mean(LWR_Subcanopy_CLM_Sodankyla(l,:));
    LWR_Subcanopy_Observations_Sodankyla_avg(l) ...
        = mean(LWR_Subcanopy_Observations_Sodankyla(l,:));
end
PAI_Sodankyla_avg = mean(PAI_Sodankyla);
f_veg = 1 - exp(-PAI_Sodankyla_avg);
LWR_Vegetation_CLM_Sodankyla_avg = (LWR_Subcanopy_CLM_Sodankyla_avg ...
    - (1-f_veg)*LWR_Atmosphere_Sodankyla)/f_veg;
LWR_Vegetation_Obs_Sodankyla_avg = (LWR_Subcanopy_Observations_Sodankyla_avg ...
    - (1-f_veg)*LWR_Atmosphere_Sodankyla)/f_veg;
LWR_Vegetation_Ratio_Sodankyla_avg = LWR_Vegetation_CLM_Sodankyla_avg...
    ./ LWR_Vegetation_Obs_Sodankyla_avg;


%-------------------------------------------------------------------------%
%-----------------------  test regression estimate  -----------------------
%-------------------------------------------------------------------------%

em_sky = 0.5:0.05:1.1;
SWin = 0:200:800;
RegrTest = nan(length(em_sky),length(SWin));
for i=1:length(em_sky)
    RegrTest(i,1) = 0.7582 + 0.2342*em_sky(i) + 0*SWin(1) ...
        + 0*SWin(1)*em_sky(i);
    for j=2:length(SWin)
        RegrTest(i,j) = 0.8685 + 0.1223*em_sky(i) ...
            + 5.2627*10^(-4)*SWin(j) - 3.6065*10^(-4)*SWin(j)*em_sky(i);
    end
end
             
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
sp = subplot(1,1,1);
hold on
scatter(Sky_Emissivity_Alptal_2005,LWR_Vegetation_Ratio_Alptal_2005,13,...
    SWR_Incoming_Alptal_2005)
scatter(Sky_Emissivity_Seehornwald_2009,LWR_Vegetation_Ratio_Seehornwald_2009,13,...
    SWR_Incoming_Seehornwald_2009)
scatter(Sky_Emissivity_Cherskiy,LWR_Vegetation_Ratio_Cherskiy,13,...
    SWR_Incoming_Cherskiy)
scatter(Sky_Emissivity_Sodankyla,LWR_Vegetation_Ratio_Sodankyla_avg,13,...
    SWR_Incoming_Sodankyla)
cm = colormap(hot(100));
caxis([0 1100])
plot(em_sky,RegrTest(:,1),'Color',cm(1,:))
plot(em_sky,RegrTest(:,2),'Color',cm(round(SWin(2)/(1101/100)),:))
plot(em_sky,RegrTest(:,3),'Color',cm(round(SWin(3)/(1101/100)),:))
plot(em_sky,RegrTest(:,4),'Color',cm(round(SWin(4)/(1101/100)),:))
plot(em_sky,RegrTest(:,5),'Color',cm(round(SWin(5)/(1101/100)),:))
hold off
cb=colorbar('SouthOutside');
% ylim([-25 25])
xlim([0.4 1.2])
set(gca,'XTick',0.4:0.2:1.2)
xlabel('Effective Emissivity of the Sky','FontSize',15,'FontWeight','bold')
ylabel('Ratio of LWR from Vegetation','FontSize',15,'FontWeight','bold')
ylabel(cb,'Insolation [W m^{-2}]','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
box on
% size & position
set(sp,'Position',[0.11,0.32,0.85,0.66])
set(cb,'Position',[0.11,0.13,0.85,0.04])
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 5];
fig.PaperSize = [6 5];
print(fig,'-dpng','-r600','MLR_Explanation.png')
set(gcf, 'Renderer', 'opengl')
print(fig,'-dpdf','-r600','MLR_Explanation.pdf')

toc