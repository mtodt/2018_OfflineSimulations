tic

clear
close all
% Create maps for analysis of global impact by correction factors.

%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
%-------------------------------------------------------------------------%

if exist('fig_num','var') == 0
    fig_num = 1;
end

% create colormaps for land without snow
bla = parula(10000);
beige   = [0.94 0.87 0.8];
    amarula = vertcat(beige,bla(2:10000,:));
iceblue = [0.1 0.9 1]; steelblue = [0.55 0.65 0.85]; blue = [0 0.447 0.741];
imperial= [0.1 0.1 0.45]; violet = [0.5 0.1 0.5];
coolio = nan(10000,3);
for c=1:3
    coolio(1:2500,c) = linspace(iceblue(c),steelblue(c),2500);
    coolio(2500:5000,c) = linspace(steelblue(c),blue(c),2501);
    coolio(5000:7500,c) = linspace(blue(c),imperial(c),2501);
    coolio(7500:9999,c) = linspace(imperial(c),violet(c),2500);
end
coolio(10000,:) = beige;

%--------------------------------------------------------------------------
%--------------------  basic variables and parameters  --------------------
%--------------------------------------------------------------------------
load('Lightbulb_Global_dimensions.mat')
    lon_str = {'30°E' '60°E' '90°E' '120°E' '150°E' '180°' '150°W' ...
        '120°W' '90°W' '60°W' '30°W'};
    lat_str = {'20°N' '40°N' '60°N' '80°N'};

% snow mask
load('Lightbulb_Global_SnowOff_control.mat','mask_perennial')

% defining snow years - dates not really important overall
Jan1 = nan(length(time)/365-1,1);      % without 1981
SYstart = nan(length(time)/365-1,1);   % without 2016
SYend = nan(length(time)/365-1,1);     % without 1981
for t=1:length(time)
    for year=1982:2016
        if time(t) == datenum(year,1,1)
            Jan1(year-1981) = t;
        elseif time(t) == datenum(year-1,10,1)
            SYstart(year-1981) = t;
        elseif time(t) == datenum(year,9,30)
            SYend(year-1981) = t;
        end
    end
end

% surface and PFT data
surfacefile = 'surfdata_0.9x1.25_simyr2000_c130418.nc';
    sinfo = ncinfo(surfacefile);

elev = ncread(surfacefile,'TOPO');  % elevation on land

NEBTs = 3-1;
NETTs = 2-1;
vegetated = 1;

PFTperc = ncread(surfacefile,'PCT_PFT');
PFTperc_name = ncreadatt(surfacefile,'PCT_PFT','long_name');
PFTperc_NETs = squeeze(PFTperc(:,:,NEBTs+1)+PFTperc(:,:,NETTs+1));

% PAIs
bla = ncread(surfacefile,'MONTHLY_LAI');
    LAI_JFM = squeeze((bla(:,:,:,1) + bla(:,:,:,2) + bla(:,:,:,3)))/3;
bla = ncread(surfacefile,'MONTHLY_SAI');
    SAI_JFM = squeeze((bla(:,:,:,1) + bla(:,:,:,2) + bla(:,:,:,3)))/3;
PAI_JFM = LAI_JFM + SAI_JFM;
clear LAI_JFM SAI_JFM

% weighted PAI
PAI_JFM_NETs = zeros(size(PFTperc_NETs));
for x=1:length(lon)
    for y=1:length(lat)
        if PFTperc_NETs(x,y) > 0
            PAI_JFM_NETs(x,y) = (PAI_JFM(x,y,2)*PFTperc(x,y,2) ...
                + PAI_JFM(x,y,3)*PFTperc(x,y,3))/PFTperc_NETs(x,y);
        end
        if isnan(mask_perennial(x,y))
            PAI_JFM_NETs(x,y) = nan;
            PFTperc_NETs(x,y) = nan;
            elev(x,y) = nan;
        end
    end
end

% rearrange maps for Europe as centre
lon_fig = nan(size(lon));
lon_fig(1:length(lon_fig)/2) = lon(length(lon_fig)/2+1:end)-360;
lon_fig(length(lon_fig)/2+1:end) = lon(1:length(lon_fig)/2);
lon_fig = vertcat(lon_fig,180);
lon_str_old = lon_str;
lon_str = {'150°W' '120°W' '90°W' '60°W' '30°W' '0°' '30°E' '60°E'...
    '90°E' '120°E' '150°E'};

elev_fig = nan(size(elev));
elev_fig(1:length(lon)/2,:) = elev(length(lon)/2+1:end,:);
elev_fig(length(lon)/2+1:end,:) = elev(1:length(lon)/2,:);

PFTperc_NETs_fig = nan(size(PFTperc_NETs));
PFTperc_NETs_fig(1:length(lon)/2,:) = PFTperc_NETs(length(lon)/2+1:end,:);
PFTperc_NETs_fig(length(lon)/2+1:end,:) = PFTperc_NETs(1:length(lon)/2,:);

PAI_JFM_NETs_fig = nan(size(PAI_JFM_NETs));
PAI_JFM_NETs_fig(1:length(lon)/2,:) = PAI_JFM_NETs(length(lon)/2+1:end,:);
PAI_JFM_NETs_fig(length(lon)/2+1:end,:) = PAI_JFM_NETs(1:length(lon)/2,:);


%--------------------------------------------------------------------------
%-----------------------  averages and rearranging  -----------------------
%--------------------------------------------------------------------------

%---------------------------------  LWE  ---------------------------------%
load Lightbulb_Global_LWEavg_control.mat
load Lightbulb_Global_LWEavg_lightbulb.mat

LWE_ctrl_NETs_DJFMAM = zeros(size(LWE_ctrl_NEBTs_DJFMAM));
LWE_lbulb_NETs_DJFMAM = zeros(size(LWE_lbulb_NEBTs_DJFMAM));
LWE_ctrl_NETs_DJFMAM_avg = zeros(size(mask_perennial)).*mask_perennial;
LWE_diff_NETs_DJFMAM_avg = zeros(size(mask_perennial)).*mask_perennial;
for x=1:length(lon)
    for y=1:length(lat)
        LWE_ctrl_NETs_DJFMAM(x,y,:) = LWE_ctrl_NETs_DJFMAM(x,y,:)*mask_perennial(x,y);
        LWE_lbulb_NETs_DJFMAM(x,y,:) = LWE_lbulb_NETs_DJFMAM(x,y,:)*mask_perennial(x,y);
        if PFTperc_NETs(x,y) > 0
            if PFTperc(x,y,NEBTs+1) > 0
                LWE_ctrl_NETs_DJFMAM(x,y,:) = LWE_ctrl_NETs_DJFMAM(x,y,:) ...
                    + LWE_ctrl_NEBTs_DJFMAM(x,y,:)*PFTperc(x,y,NEBTs+1)/PFTperc_NETs(x,y);
                LWE_lbulb_NETs_DJFMAM(x,y,:) = LWE_lbulb_NETs_DJFMAM(x,y,:) ...
                    + LWE_lbulb_NEBTs_DJFMAM(x,y,:)*PFTperc(x,y,NEBTs+1)/PFTperc_NETs(x,y);
            end
            if PFTperc(x,y,NETTs+1) > 0
                LWE_ctrl_NETs_DJFMAM(x,y,:) = LWE_ctrl_NETs_DJFMAM(x,y,:) ...
                    + LWE_ctrl_NETTs_DJFMAM(x,y,:)*PFTperc(x,y,NETTs+1)/PFTperc_NETs(x,y);
                LWE_lbulb_NETs_DJFMAM(x,y,:) = LWE_lbulb_NETs_DJFMAM(x,y,:) ...
                    + LWE_lbulb_NETTs_DJFMAM(x,y,:)*PFTperc(x,y,NETTs+1)/PFTperc_NETs(x,y);
            end       
            LWE_ctrl_NETs_DJFMAM_avg(x,y) = mean(squeeze(LWE_ctrl_NETs_DJFMAM(x,y,:)));
            LWE_diff_NETs_DJFMAM_avg(x,y) ...
                = mean(squeeze(LWE_lbulb_NETs_DJFMAM(x,y,:)-LWE_ctrl_NETs_DJFMAM(x,y,:)));
        end
    end
end

% rearrange maps for Europe as centre
LWE_ctrl_NETs_DJFMAM_fig = nan(size(LWE_ctrl_NETs_DJFMAM_avg));
    LWE_ctrl_NETs_DJFMAM_fig(1:length(lon)/2,:)...
        = LWE_ctrl_NETs_DJFMAM_avg(length(lon)/2+1:end,:);
    LWE_ctrl_NETs_DJFMAM_fig(length(lon)/2+1:end,:)...
        = LWE_ctrl_NETs_DJFMAM_avg(1:length(lon)/2,:);
LWE_diff_NETs_DJFMAM_fig = nan(size(LWE_diff_NETs_DJFMAM_avg));
    LWE_diff_NETs_DJFMAM_fig(1:length(lon)/2,:)...
        = LWE_diff_NETs_DJFMAM_avg(length(lon)/2+1:end,:);
    LWE_diff_NETs_DJFMAM_fig(length(lon)/2+1:end,:)...
        = LWE_diff_NETs_DJFMAM_avg(1:length(lon)/2,:);

LWE_ctrl_NETs_DJFMAM_fig ...
    = vertcat(LWE_ctrl_NETs_DJFMAM_fig,LWE_ctrl_NETs_DJFMAM_fig(1,:));
LWE_diff_NETs_DJFMAM_fig ...
    = vertcat(LWE_diff_NETs_DJFMAM_fig,LWE_diff_NETs_DJFMAM_fig(1,:));


%-----------------------  snow surface temperature  -----------------------
load('Lightbulb_Global_TsnotAvg_control.mat')
load('Lightbulb_Global_TsnotAvg_lightbulb.mat')

Tsnot_ctrl_avg = nan(size(mask_perennial));
Tsnot_diff_avg = nan(size(mask_perennial));
for x=1:length(lon)
    for y=1:length(lat)
        if mask_perennial(x,y) == 1
            for t=1:length(Jan1)
                if Tsnot_ctrl_SYavg(x,y,t) > 275
                    Tsnot_ctrl_SYavg(x,y,t) = nan;
                end
                if Tsnot_lbulb_SYavg(x,y,t) > 275
                    Tsnot_lbulb_SYavg(x,y,t) = nan;
                end
            end
            Tsnot_ctrl_avg(x,y) ...
                = nanmean(squeeze(Tsnot_ctrl_SYavg(x,y,:)));
            Tsnot_diff_avg(x,y) ...
                = nanmean(squeeze(Tsnot_lbulb_SYavg(x,y,:)-Tsnot_ctrl_SYavg(x,y,:)));
        elseif mask_perennial(x,y) == 0
            Tsnot_ctrl_avg(x,y) = 273.15+1;
            Tsnot_diff_avg(x,y) = 0;
        end
    end
end

% rearrange maps for Europe as centre
Tsnot_ctrl_avg_fig = nan(size(Tsnot_ctrl_avg));
    Tsnot_ctrl_avg_fig(1:length(lon)/2,:)...
        = Tsnot_ctrl_avg(length(lon)/2+1:end,:);
    Tsnot_ctrl_avg_fig(length(lon)/2+1:end,:)...
        = Tsnot_ctrl_avg(1:length(lon)/2,:);
Tsnot_diff_avg_fig = nan(size(Tsnot_diff_avg));
    Tsnot_diff_avg_fig(1:length(lon)/2,:)...
        = Tsnot_diff_avg(length(lon)/2+1:end,:);
    Tsnot_diff_avg_fig(length(lon)/2+1:end,:)...
        = Tsnot_diff_avg(1:length(lon)/2,:);

Tsnot_ctrl_avg_fig = vertcat(Tsnot_ctrl_avg_fig,Tsnot_ctrl_avg_fig(1,:));
Tsnot_diff_avg_fig = vertcat(Tsnot_diff_avg_fig,Tsnot_diff_avg_fig(1,:));


%-----------------------------  cold content  -----------------------------
load('Lightbulb_Global_CCsnowAvg_control.mat')
load('Lightbulb_Global_CCsnowAvg_lightbulb.mat')

CCsnow_ctrl_avg = nan(size(mask_perennial));
CCsnow_diff_avg = nan(size(mask_perennial));
CCsnow_ratio_avg = nan(size(mask_perennial));
for x=1:length(lon)
    for y=1:length(lat)
        if mask_perennial(x,y) == 1
            CCsnow_ctrl_avg(x,y) = nanmean(squeeze(CCsnow_ctrl_SYavg(x,y,:)));
            CCsnow_diff_avg(x,y) ...
                = nanmean(squeeze(CCsnow_lbulb_SYavg(x,y,:)-CCsnow_ctrl_SYavg(x,y,:)));
            CCsnow_ratio_avg(x,y) ...
                = nanmean(squeeze(CCsnow_lbulb_SYavg(x,y,:)-CCsnow_ctrl_SYavg(x,y,:)))...
                /CCsnow_ctrl_avg(x,y);
            
        elseif mask_perennial(x,y) == 0
            CCsnow_ctrl_avg(x,y) = -1;  % values set for maps
            CCsnow_diff_avg(x,y) = 0;  % values set for maps
            CCsnow_ratio_avg(x,y) = 0;  % values set for maps
        end
    end
end

% rearrange maps for Europe as centre
CCsnow_ctrl_avg_fig = nan(size(CCsnow_ctrl_avg));
    CCsnow_ctrl_avg_fig(1:length(lon)/2,:)...
        = CCsnow_ctrl_avg(length(lon)/2+1:end,:);
    CCsnow_ctrl_avg_fig(length(lon)/2+1:end,:)...
        = CCsnow_ctrl_avg(1:length(lon)/2,:);
CCsnow_ratio_avg_fig = nan(size(CCsnow_ratio_avg));
    CCsnow_ratio_avg_fig(1:length(lon)/2,:)...
        = CCsnow_ratio_avg(length(lon)/2+1:end,:);
    CCsnow_ratio_avg_fig(length(lon)/2+1:end,:)...
        = CCsnow_ratio_avg(1:length(lon)/2,:);

CCsnow_ctrl_avg_fig = vertcat(CCsnow_ctrl_avg_fig,CCsnow_ctrl_avg_fig(1,:));
CCsnow_ratio_avg_fig = vertcat(CCsnow_ratio_avg_fig,CCsnow_ratio_avg_fig(1,:));


%----------------------------  snow-off date  ----------------------------%
load('Lightbulb_Global_SnowOff_control.mat')
load('Lightbulb_Global_SnowOff_lightbulb.mat')

SnowOff_ctrl_avg = nan(size(mask_perennial));
SnowOff_diff_avg = nan(size(mask_perennial));
for x=1:length(lon)
    for y=1:length(lat)
        if mask_perennial(x,y) == 0
            SnowOff_ctrl(x,y,:) = -1;       % values set for maps
            SnowOff_lbulb(x,y,:) = -1;      % values set for maps
        end
        SnowOff_ctrl_avg(x,y) = mean(squeeze(SnowOff_ctrl(x,y,:)));
        SnowOff_diff_avg(x,y) ...
            = mean(squeeze(SnowOff_lbulb(x,y,:)-SnowOff_ctrl(x,y,:)));
    end
end

% rearrange maps for Europe as centre
SnowOff_ctrl_avg_fig = nan(size(SnowOff_ctrl_avg));
    SnowOff_ctrl_avg_fig(1:length(lon)/2,:)...
        = SnowOff_ctrl_avg(length(lon)/2+1:end,:);
    SnowOff_ctrl_avg_fig(length(lon)/2+1:end,:)...
        = SnowOff_ctrl_avg(1:length(lon)/2,:);
SnowOff_diff_avg_fig = nan(size(SnowOff_diff_avg));
    SnowOff_diff_avg_fig(1:length(lon)/2,:)...
        = SnowOff_diff_avg(length(lon)/2+1:end,:);
    SnowOff_diff_avg_fig(length(lon)/2+1:end,:)...
        = SnowOff_diff_avg(1:length(lon)/2,:);

SnowOff_ctrl_avg_fig ...
    = vertcat(SnowOff_ctrl_avg_fig,SnowOff_ctrl_avg_fig(1,:));
SnowOff_diff_avg_fig ...
    = vertcat(SnowOff_diff_avg_fig,SnowOff_diff_avg_fig(1,:));



%-------------------------------------------------------------------------%
%----------------------------  paper figures  ----------------------------%
%-------------------------------------------------------------------------%

elev_fig = vertcat(elev_fig,elev_fig(1,:));
PFTperc_NETs_fig = vertcat(PFTperc_NETs_fig,PFTperc_NETs_fig(1,:));
PAI_JFM_NETs_fig = vertcat(PAI_JFM_NETs_fig,PAI_JFM_NETs_fig(1,:));
for x=1:length(lon_fig)
    for y=1:length(lat)
        if PFTperc_NETs_fig(x,y) == 0
            PAI_JFM_NETs_fig(x,y) = -1;     % values set for maps
        end
    end
end

fig=figure(fig_num);fig_num = fig_num+1;
% PFT percentage for NEBTs + NETTs
sp1=subplot(3,1,1);
hold on
tx1 = text(-1,1,'a','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,PFTperc_NETs_fig')
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-');
colormap(sp1,x2b2y('single',0,0,100,'vegetation'))
cb1=colorbar('EastOutside');
ylabel(cb1,'Plant Functional Type coverage [%]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
% PAI weighted for NEBTs and NETTs
sp2=subplot(3,1,2);
hold on
tx2 = text(-1,1,'b','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,PAI_JFM_NETs_fig')
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-');
colormap(sp2,amarula)
cb2=colorbar('EastOutside');
caxis([-0.001 6])
ylabel(cb2,'Plant Area Index [m^2 m^{-2}]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
% elevation
sp3=subplot(3,1,3);
hold on
tx3 = text(-1,1,'c','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,elev_fig')
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-');
colormap(sp3,flipud(x2b2y('single',0,1500,3000,'vegetation')))
cb3=colorbar('EastOutside');
ylabel(cb3,'Elevation [m]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
% position & size
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6.432 15];
fig.PaperSize = [6.432 15];
set(sp1,'Position',[0.03,0.69,0.75,0.29])
set(cb1,'Position',[0.83,0.7,0.03,0.27])
set(sp2,'Position',[0.03,0.355,0.75,0.29])
set(cb2,'Position',[0.83,0.365,0.03,0.27])
set(sp3,'Position',[0.03,0.02,0.75,0.29])
set(cb3,'Position',[0.83,0.03,0.03,0.27])
print(fig,'-dpng','-r600','InfoMaps_Paper_vertical.png')
print(fig,'-dpdf','-r600','InfoMaps_Paper_vertical.pdf')


%------------------------  deficiencies concerns  ------------------------%
load('FaultFrequences.mat')
FaultFrequenceDiff = nan(size(FaultFrequenceDiff_NEBTs));
FaultFrequenceDiff(~isnan(mask_perennial)) = -1;
for x=1:length(lon)
    for y=1:length(lat)
        if PFTperc_NETs(x,y) > 0
            FaultFrequenceDiff(x,y) ...
                = nanmean([FaultFrequenceDiff_NEBTs(x,y) FaultFrequenceDiff_NETTs(x,y)]);
        end
    end
end

FaultFrequenceDiff_fig = nan(size(FaultFrequenceDiff));
FaultFrequenceDiff_fig(1:length(lon)/2,:) = FaultFrequenceDiff(length(lon)/2+1:end,:);
FaultFrequenceDiff_fig(length(lon)/2+1:end,:) = FaultFrequenceDiff(1:length(lon)/2,:);
FaultFrequenceDiff_fig = vertcat(FaultFrequenceDiff_fig,FaultFrequenceDiff_fig(1,:));


fig=figure(fig_num);fig_num = fig_num+1;
sp=subplot(1,1,1);
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,FaultFrequenceDiff_fig')
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',8);
colormap(sp,amarula)
caxis([-0.0001 0.5])
cb=colorbar('SouthOutside');
ylabel(cb,'Frequency','FontSize',10)
set(gca,'FontSize',10)
set(gcf, 'Renderer', 'opengl')
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 8.3 9.4];
fig.PaperSize = [8.3 9.4];
set(sp,'Position',[0.075,0.15,0.85,0.85])
set(cb,'Position',[0.1,0.11,0.8,0.03])
print(fig,'-dpng','-r600','FaultFrequency.png')
print(fig,'-dpdf','-r600','FaultFrequency.pdf')

% contour line for 0.1 threshold
threshold = 0.1;
bla = contour(lon_fig,lat,FaultFrequenceDiff_fig',[threshold threshold],...
    'LineWidth',2,'LineColor','k');
BoundFault = nan(size(bla));
for b=1:length(bla(1,:))
    if round(bla(1,b),2) ~= threshold
        BoundFault(1,b) = bla(1,b);
        BoundFault(2,b) = bla(2,b);
    end
end
blablub = 10^6*ones(length(bla(1,:)),1);
line(BoundFault(1,300:1100),BoundFault(2,300:1100),blablub(300:1100),...
    'Color','k','LineWidth',2)


fig=figure(fig_num);fig_num = fig_num+1;
    % control LWE
sp1=subplot(4,2,1);
hold on
tx1 = text(-1,1,'a','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,LWE_ctrl_NETs_DJFMAM_fig')
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',7);
colormap(sp1,amarula)
caxis([0.9 1.5])
cb1=colorbar('WestOutside');
ylabel(cb1,'Longwave Enhancement','FontSize',12)
t1 = title('CTRL','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
    % difference in LWE
sp2=subplot(4,2,2);
hold on
tx2 = text(-1,1,'b','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,LWE_diff_NETs_DJFMAM_fig')
m_line(BoundFault(1,300:1100),BoundFault(2,300:1100),'Color','k','LineWidth',0.5)
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',7);
colormap(sp2,x2b2y('single',-0.06,0,0.06,'hotncold')) 
cb2=colorbar('EastOutside');
set(cb2,'YTick',-0.06:0.03:0.06)
ylabel(cb2,'Longwave Enhancement','FontSize',12)
t2 = title('CORR - CTRL','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
    % control snow surface temperature
sp3=subplot(4,2,3);
hold on
tx3 = text(-1,1,'c','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,Tsnot_ctrl_avg_fig'-273.15)
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',7);
colormap(sp3,coolio)
caxis([-40 0.01])
cb3=colorbar('WestOutside');
ylabel(cb3,'Snow Surface Temperature [°C]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
    % difference in snow surface temperature
sp4=subplot(4,2,4);
hold on
tx4 = text(-1,1,'d','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,Tsnot_diff_avg_fig')
m_line(BoundFault(1,300:1100),BoundFault(2,300:1100),'Color','k','LineWidth',0.5)
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',7);
colormap(sp4,x2b2y('single',-3,0,3,'hotncold'))
cb4=colorbar('EastOutside');
set(cb4,'YTick',-3:1:3)
ylabel(cb4,'Snow Surface Temperature [°C]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
    % control cold content
sp5=subplot(4,2,5);
hold on
tx5 = text(-1,1,'e','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,CCsnow_ctrl_avg_fig')
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',7);
colormap(sp5,flipud(coolio))
caxis([-0.001 5])
cb5=colorbar('WestOutside');
ylabel(cb5,'Cold Content [MJ m^{-2}]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
    % difference in cold content
sp6=subplot(4,2,6);
hold on
tx6 = text(-1,1,'f','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,CCsnow_ratio_avg_fig'*100)
m_line(BoundFault(1,300:1100),BoundFault(2,300:1100),'Color','k','LineWidth',0.5)
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',7);
colormap(sp6,flipud(x2b2y('single',-50,0,50,'hotncold')))
cb6=colorbar('EastOutside');
set(cb6,'YTick',-50:25:50)
ylabel(cb6,'Cold Content [%]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
    % control snow off date
sp7=subplot(4,2,7);
hold on
tx7 = text(-1,1,'g','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,SnowOff_ctrl_avg_fig')
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',7);
colormap(sp7,amarula)
caxis([-0.1 200])
cb7=colorbar('WestOutside');
ylabel(cb7,'Snow Off Date [DoY]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
    % difference in snow off date
sp8=subplot(4,2,8);
hold on
tx8 = text(-1,1,'h','FontSize',12,'FontWeight','bold');
m_proj('stereographic','lat',90,'long',0,'radius',50);
m_pcolor(lon_fig,lat,SnowOff_diff_avg_fig')
m_line(BoundFault(1,300:1100),BoundFault(2,300:1100),'Color','k','LineWidth',0.5)
m_grid('xtick',12,'tickdir','out','ytick',[50 60 70 80],'linest','-','fontsize',7);
colormap(sp8,flipud(x2b2y('single',-10,0,10,'hotncold')))
cb8=colorbar('EastOutside');
ylabel(cb8,'Snow Off Date [days]','FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Renderer', 'opengl')
    % position & size
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 15 18];
fig.PaperSize = [15 18];
set(t1,'Position',[0,1.18,-5e+15])
set(t2,'Position',[0,1.18,-5e+15])
set(sp1,'Position',[0.19,0.726,0.25,0.25])
set(cb1,'Position',[0.075,0.751,0.03,0.2])
set(sp2,'Position',[0.54,0.726,0.25,0.25])
set(cb2,'Position',[0.88,0.751,0.03,0.2])
set(sp3,'Position',[0.19,0.484,0.25,0.25])
set(cb3,'Position',[0.075,0.509,0.03,0.2])
set(sp4,'Position',[0.54,0.484,0.25,0.25])
set(cb4,'Position',[0.88,0.509,0.03,0.2])
set(sp5,'Position',[0.19,0.242,0.25,0.25])
set(cb5,'Position',[0.075,0.267,0.03,0.2])
set(sp6,'Position',[0.54,0.242,0.25,0.25])
set(cb6,'Position',[0.88,0.267,0.03,0.2])
set(sp7,'Position',[0.19,0,0.25,0.25])
set(cb7,'Position',[0.075,0.025,0.03,0.2])
set(sp8,'Position',[0.54,0,0.25,0.25])
set(cb8,'Position',[0.88,0.025,0.03,0.2])
print(fig,'-dpng','-r600','ResultsMaps_Paper_contour.png')
print(fig,'-dpdf','-r600','ResultsMaps_Paper_contour.pdf')



%----------------------  explanation scatter graph  ----------------------%
fig=figure(fig_num);fig_num = fig_num+1;
% snow off date
sp1=subplot(1,3,1);
hold on
text(-30+0.05*40,-8+0.95*16,'a','FontSize',15,'FontWeight','bold')
for x=33:113
scatter(CCsnow_ratio_avg(x,144:171)*100,SnowOff_diff_avg(x,144:171),17,...
    SnowOff_ctrl_avg(x,144:171))
caxis([0 200])
end
xlim([-30 10])
ylim([-8 8]) %ylim([-6 4])
xlabel('Change in Cold Content [%]','FontSize',15)
ylabel('Change in Snow Off Date [days]','FontSize',15)
cb1=colorbar('SouthOutside');
ycb1 = ylabel(cb1,'Snow Off Date [DoY]','FontSize',15);
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
set(gcf, 'Renderer', 'opengl')
box on
% cold content
sp2=subplot(1,3,2);
hold on
text(-30+0.05*40,-8+0.95*16,'b','FontSize',15,'FontWeight','bold')
for x=33:113
scatter(CCsnow_ratio_avg(x,144:171)*100,SnowOff_diff_avg(x,144:171),17,...
    CCsnow_ctrl_avg(x,144:171))
caxis([-0.001 3])
end
colormap(sp2,flipud(coolio))
xlim([-30 10])
ylim([-8 8]) %ylim([-6 4])
xlabel('Change in Cold Content [%]','FontSize',15)
cb2=colorbar('SouthOutside');
ycb2 = ylabel(cb2,'Cold Content [MJ m^{-2}]','FontSize',15);
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
set(gcf, 'Renderer', 'opengl')
box on
get(ycb2,'position')%1.4995 -1.1333 0 - but after printing: 1.4995 -2.3538 0
set(ycb2,'position',[1.4995 -1.53 0])
% elevation
sp3=subplot(1,3,3);
hold on
text(-30+0.05*40,-8+0.95*16,'c','FontSize',15,'FontWeight','bold')
for x=33:113
scatter(CCsnow_ratio_avg(x,144:171)*100,SnowOff_diff_avg(x,144:171),17,...
    elev(x,144:171))
caxis([0 1000])
end
colormap(sp3,flipud(x2b2y('single',0,500,1000,'vegetation')))
xlim([-30 10])
ylim([-8 8]) %ylim([-6 4])
xlabel('Change in Cold Content [%]','FontSize',15)
cb3=colorbar('SouthOutside');
ycb3 = ylabel(cb3,'Elevation [m]','FontSize',15);
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
set(gcf, 'Renderer', 'opengl')
box on
% position & size
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 15 6];
fig.PaperSize = [15 6];
set(sp1,'Position',[0.04,0.265,0.29,0.72])
set(cb1,'Position',[0.05,0.1075,0.27,0.03])
set(sp2,'Position',[0.37,0.265,0.29,0.72])
set(cb2,'Position',[0.38,0.1075,0.27,0.03])
set(sp3,'Position',[0.7,0.265,0.29,0.72])
set(cb3,'Position',[0.71,0.1075,0.27,0.03])
print(fig,'-dpng','-r600','ExplainerScatter_Siberia_Paper.png')
print(fig,'-dpdf','-r600','ExplainerScatter_Siberia_Paper.pdf')

toc