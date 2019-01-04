tic
%{
Create ultimate script to create paper(-quality) graphs for lightbulb
analysis.
%}

%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
%-------------------------------------------------------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end

% site colours
GlacialGrey = [0.4 0.7 1];
MagicMaroon = [0.65 0.32 0.35];
EcstaticEmerald = [0.17 0.52 0.5];
GorgeousGold = [1 0.85 0];
BastilleBleu = [0.12 0.34 0.56];
RussianRed = [1 0.032 0.064];
ViolentViolet = [0.6 0 0.7];

% met colours
GargantuanGreen = [0.1 0.5 0.2];
matlabblue = [0 0.447 0.741];
matlabred = [0.85 0.325 0.098];
matlabyellow = [0.929 0.694 0.125];


%--------------------------------------------------------------------------
%----------------------------  site locations  ----------------------------
%--------------------------------------------------------------------------
Abisko_loc = [16 169];      % 18.8206E, 68.3556N
Alptal_loc = [8 146];       % 8.8E, 47.05N
Borden_loc = [225 144];     % 79.9333W == 280.0667E, 44.3167N
Cherskiy_loc = [130 169];   % 161.45482E, 68.75574N
Seehornwald_loc = [9 146];  % 9.8558E, 46.8153N
Sodankyla_loc = [22 168];   % 26.6333E, 67.3667N
Yakutsk_loc = [105 163];    % 129.6189E, 62.255N
% sparse site close-ish to Sodankyla
Sparse_loc = [25 172];      % 30E, 71.15N


%--------------------------------------------------------------------------
%------------------  import of global simulation output  ------------------
%--------------------------------------------------------------------------
load('LightbulbAnalysis_Alptal_Meteorology.mat')
load('LightbulbAnalysis_Alptal_Control.mat')
load('LightbulbAnalysis_Alptal_Lightbulb_LWsub.mat')


%--------------------------------------------------------------------------
%--------------------  import forest stand-scale data  --------------------
%--------------------------------------------------------------------------
load('MLRdata_Alptal.mat','Hourly_Time_Alptal_2004','LWR_Atmosphere_Alptal_2004',...
    'LWR_Subcanopy_Observations_Alptal_2004','Sky_Emissivity_Alptal_2004',...
    'Hourly_Time_Alptal_2005','LWR_Atmosphere_Alptal_2005',...
    'LWR_Subcanopy_Observations_Alptal_2005','Sky_Emissivity_Alptal_2005',...
    'Hourly_Time_Alptal_2006','LWR_Atmosphere_Alptal_2006',...
    'LWR_Subcanopy_Observations_Alptal_2006','Sky_Emissivity_Alptal_2006',...
    'Hourly_Time_Alptal_2007','LWR_Atmosphere_Alptal_2007',...
    'LWR_Subcanopy_Observations_Alptal_2007','Sky_Emissivity_Alptal_2007')

bla = horzcat(LWR_Atmosphere_Alptal_2004,LWR_Subcanopy_Observations_Alptal_2004,...
    Sky_Emissivity_Alptal_2004);
[Time_obs_Alp_04,Obs_Alp_04] = NaNgapper(24,Hourly_Time_Alptal_2004,bla);
clear bla

bla = horzcat(LWR_Atmosphere_Alptal_2005,LWR_Subcanopy_Observations_Alptal_2005,...
    Sky_Emissivity_Alptal_2005);
[Time_obs_Alp_05,Obs_Alp_05] = NaNgapper(24,Hourly_Time_Alptal_2005,bla);
clear bla

bla = horzcat(LWR_Atmosphere_Alptal_2006,LWR_Subcanopy_Observations_Alptal_2006,...
    Sky_Emissivity_Alptal_2006);
[Time_obs_Alp_06,Obs_Alp_06] = NaNgapper(24,Hourly_Time_Alptal_2006,bla);
clear bla

bla = horzcat(LWR_Atmosphere_Alptal_2007,LWR_Subcanopy_Observations_Alptal_2007,...
    Sky_Emissivity_Alptal_2007);
[Time_obs_Alp_07,Obs_Alp_07] = NaNgapper(24,Hourly_Time_Alptal_2007,bla);


%-------------------------------------------------------------------------%
%--------------------------  analysis & graphs  --------------------------%
%-------------------------------------------------------------------------%
% sub-canopy LWR
%{
Fortunately, Alptal is in CET (UTC+1 in winter), so that CLM4.5 and
observations both start at 1:00 CET or 0:00 UTC.
%}
LWsub_ctrl_Alp_04_dc = nan(24,1); LWsub_lbulb_Alp_04_dc = nan(24,1);
LWsub_obs_Alp_04_dc = nan(24,1);
LWsub_ctrl_Alp_05_dc = nan(24,1); LWsub_lbulb_Alp_05_dc = nan(24,1);
LWsub_obs_Alp_05_dc = nan(24,1);
LWsub_ctrl_Alp_06_dc = nan(24,1); LWsub_lbulb_Alp_06_dc = nan(24,1);
LWsub_obs_Alp_06_dc = nan(24,1);
LWsub_ctrl_Alp_07_dc = nan(24,1); LWsub_lbulb_Alp_07_dc = nan(24,1);
LWsub_obs_Alp_07_dc = nan(24,1);
for h=1:24
    LWsub_ctrl_Alp_04_dc(h) = mean(LWsub_ctrl_Alp_04(h:24:end-24+h));
    LWsub_lbulb_Alp_04_dc(h) = mean(LWsub_lbulb_Alp_04(h:24:end-24+h));
    LWsub_obs_Alp_04_dc(h) ...
        = mean(LWR_Subcanopy_Observations_Alptal_2004(h:24:end-24+h));
    LWsub_ctrl_Alp_05_dc(h) = mean(LWsub_ctrl_Alp_05(h:24:end-24+h));
    LWsub_lbulb_Alp_05_dc(h) = mean(LWsub_lbulb_Alp_05(h:24:end-24+h));
    LWsub_obs_Alp_05_dc(h) ...
        = mean(LWR_Subcanopy_Observations_Alptal_2005(h:24:end-24+h));
    LWsub_ctrl_Alp_06_dc(h) = mean(LWsub_ctrl_Alp_06(h:24:end-24+h));
    LWsub_lbulb_Alp_06_dc(h) = mean(LWsub_lbulb_Alp_06(h:24:end-24+h));
    LWsub_obs_Alp_06_dc(h) ...
        = mean(LWR_Subcanopy_Observations_Alptal_2006(h:24:end-24+h));
    LWsub_ctrl_Alp_07_dc(h) = mean(LWsub_ctrl_Alp_07(h:24:end-24+h));
    LWsub_lbulb_Alp_07_dc(h) = mean(LWsub_lbulb_Alp_07(h:24:end-24+h));
    LWsub_obs_Alp_07_dc(h) ...
        = mean(LWR_Subcanopy_Observations_Alptal_2007(h:24:end-24+h));
end


% longwave enhancement
LWE_ctrl_Alp_04 = LWsub_ctrl_Alp_04./LWatm_Alp_04;
LWE_lbulb_Alp_04 = LWsub_lbulb_Alp_04./LWatm_Alp_04;
LWE_obs_Alp_04 = squeeze(Obs_Alp_04(:,2))./squeeze(Obs_Alp_04(:,1));
LWE_ctrl_Alp_05 = LWsub_ctrl_Alp_05./LWatm_Alp_05;
LWE_lbulb_Alp_05 = LWsub_lbulb_Alp_05./LWatm_Alp_05;
LWE_obs_Alp_05 = squeeze(Obs_Alp_05(:,2))./squeeze(Obs_Alp_05(:,1));
LWE_ctrl_Alp_06 = LWsub_ctrl_Alp_06./LWatm_Alp_06;
LWE_lbulb_Alp_06 = LWsub_lbulb_Alp_06./LWatm_Alp_06;
LWE_obs_Alp_06 = squeeze(Obs_Alp_06(:,2))./squeeze(Obs_Alp_06(:,1));
LWE_ctrl_Alp_07 = LWsub_ctrl_Alp_07./LWatm_Alp_07;
LWE_lbulb_Alp_07 = LWsub_lbulb_Alp_07./LWatm_Alp_07;
LWE_obs_Alp_07 = squeeze(Obs_Alp_07(:,2))./squeeze(Obs_Alp_07(:,1));
%{
Fortunately, Alptal is in CET (UTC+1 in winter), so that CLM4.5 and
observations both start at 1:00 CET or 0:00 UTC.
%}
LWE_ctrl_Alp_04_dc = nan(24,1); LWE_lbulb_Alp_04_dc = nan(24,1);
LWE_obs_Alp_04_dc = nan(24,1);
LWE_ctrl_Alp_05_dc = nan(24,1); LWE_lbulb_Alp_05_dc = nan(24,1);
LWE_obs_Alp_05_dc = nan(24,1);
LWE_ctrl_Alp_06_dc = nan(24,1); LWE_lbulb_Alp_06_dc = nan(24,1);
LWE_obs_Alp_06_dc = nan(24,1);
LWE_ctrl_Alp_07_dc = nan(24,1); LWE_lbulb_Alp_07_dc = nan(24,1);
LWE_obs_Alp_07_dc = nan(24,1);
for h=1:24
    LWE_ctrl_Alp_04_dc(h) = mean(LWE_ctrl_Alp_04(h:24:end-24+h));
    LWE_lbulb_Alp_04_dc(h) = mean(LWE_lbulb_Alp_04(h:24:end-24+h));
    bla = LWR_Subcanopy_Observations_Alptal_2004./LWR_Atmosphere_Alptal_2004;
    LWE_obs_Alp_04_dc(h) = mean(bla(h:24:end-24+h));
    LWE_ctrl_Alp_05_dc(h) = mean(LWE_ctrl_Alp_05(h:24:end-24+h));
    LWE_lbulb_Alp_05_dc(h) = mean(LWE_lbulb_Alp_05(h:24:end-24+h));
    bla = LWR_Subcanopy_Observations_Alptal_2005./LWR_Atmosphere_Alptal_2005;
    LWE_obs_Alp_05_dc(h) = mean(bla(h:24:end-24+h));
    LWE_ctrl_Alp_06_dc(h) = mean(LWE_ctrl_Alp_06(h:24:end-24+h));
    LWE_lbulb_Alp_06_dc(h) = mean(LWE_lbulb_Alp_06(h:24:end-24+h));
    bla = LWR_Subcanopy_Observations_Alptal_2006./LWR_Atmosphere_Alptal_2006;
    LWE_obs_Alp_06_dc(h) = mean(bla(h:24:end-24+h));
    LWE_ctrl_Alp_07_dc(h) = mean(LWE_ctrl_Alp_07(h:24:end-24+h));
    LWE_lbulb_Alp_07_dc(h) = mean(LWE_lbulb_Alp_07(h:24:end-24+h));
    bla = LWR_Subcanopy_Observations_Alptal_2007./LWR_Atmosphere_Alptal_2007;
    LWE_obs_Alp_07_dc(h) = mean(bla(h:24:end-24+h));
end


time_ticks = [0 15 31 44 59 73 89];
time_str = {'1 Jan','16 Jan','1 Feb','14 Feb','1 Mar','15 Mar','31 Mar'};
fig=figure(fig_num);fig_num = fig_num+1;
    % time series of sub-canopy LWR at Alptal 2006
ts1=subplot(2,2,1);
hold on
text(0+90*0.05/(0.72/0.12),150+0.95*250,'a','FontSize',15,'FontWeight','bold')
plot(time_06-datenum(2006,01,01,00,00,00),LWsub_ctrl_Alp_06,'k')
plot(time_06-datenum(2006,01,01,00,00,00),LWsub_lbulb_Alp_06,'r')
plot(Time_obs_Alp_06-datenum(2006,01,01,00,00,00),Obs_Alp_06(:,2),...
    'Color',EcstaticEmerald)
xlim([0 90])
ylim([150 400])
set(gca,'XTick',time_ticks,'XTickLabel',time_str)
ylabel('sub-canopy LWR [W m^{-2}]','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
ctrl=plot(-1,-1,'MarkerEdgeColor','k','MarkerFaceColor','k',...
    'Marker','o','MarkerSize',10,'LineStyle','none');
lbulb=plot(-1,-1,'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'Marker','o','MarkerSize',10,'LineStyle','none');
obs=plot(-1,-1,'MarkerEdgeColor',EcstaticEmerald,'MarkerFaceColor',...
    EcstaticEmerald,'Marker','o','MarkerSize',10,'LineStyle','none');
legend([ctrl,lbulb,obs],'CTRL','CORR','OBS',...
    'location','southwest','orientation','horizontal',...
    'FontSize',15,'FontWeight','bold')
legend('boxoff')
box on
    % diurnal cycle of sub-canopy LWR at Alptal 2006
dc1=subplot(2,2,2);
hold on
text(1+23*0.05,225+0.95*125,'b','FontSize',15,'FontWeight','bold')
plot([1 24],[mean(LWsub_ctrl_Alp_06) mean(LWsub_ctrl_Alp_06)],':k')
plot([1 24],[mean(LWsub_lbulb_Alp_06) mean(LWsub_lbulb_Alp_06)],':r')
    bla = mean(LWR_Subcanopy_Observations_Alptal_2006);
plot([1 24],[bla bla],'Color',EcstaticEmerald,'LineStyle',':')
plot(LWsub_ctrl_Alp_06_dc,'k','LineWidth',2)
plot(LWsub_lbulb_Alp_06_dc,'r','LineWidth',2)
plot(LWsub_obs_Alp_06_dc,'Color',EcstaticEmerald,'LineWidth',2)
xlim([1 24])
set(gca,'XTick',6:6:18)
ylim([225 350])
set(gca,'YTick',225:25:350)
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
box on
    % time series of LW enhancement at Alptal 2006
ts2=subplot(2,2,3);
hold on
text(0+90*0.05/(0.72/0.12),0.8+0.95*1,'c','FontSize',15,'FontWeight','bold')
plot(time_06-datenum(2006,01,01,00,00,00),LWE_ctrl_Alp_06,'k')
plot(time_06-datenum(2006,01,01,00,00,00),LWE_lbulb_Alp_06,'r')
plot(Time_obs_Alp_06-datenum(2006,01,01,00,00,00),...
    LWE_obs_Alp_06,'Color',EcstaticEmerald)
xlim([0 90])
ylim([0.8 1.8])
set(gca,'XTick',time_ticks,'XTickLabel',time_str)
xlabel('2006','FontSize',15,'FontWeight','bold')
ylabel('Longwave Enhancement','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
box on
    % diurnal cycle of LW enhancement at Alptal 2006
dc2=subplot(2,2,4);
hold on
text(1+23*0.05,1.1+0.95*0.4,'d','FontSize',15,'FontWeight','bold')
plot([1 24],[mean(LWE_ctrl_Alp_06) mean(LWE_ctrl_Alp_06)],':k')
plot([1 24],[mean(LWE_lbulb_Alp_06) mean(LWE_lbulb_Alp_06)],':r')
    bla = mean(LWR_Subcanopy_Observations_Alptal_2006./LWR_Atmosphere_Alptal_2006);
plot([1 24],[bla bla],'Color',EcstaticEmerald,'LineStyle',':')
plot(LWE_ctrl_Alp_06_dc,'k','LineWidth',2)
plot(LWE_lbulb_Alp_06_dc,'r','LineWidth',2)
plot(LWE_obs_Alp_06_dc,'Color',EcstaticEmerald,'LineWidth',2)
xlim([1 24])
set(gca,'XTick',6:6:18)
ylim([1.1 1.5])
xlabel('Hour of Day','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
box on
% position & size
set(ts1,'Position',[0.07,0.56,0.7,0.41])
set(ts2,'Position',[0.07,0.075,0.7,0.41])
set(dc1,'Position',[0.83,0.56,0.16,0.41])
set(dc2,'Position',[0.83,0.075,0.16,0.41])
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 11.7 8.3]; % A4
fig.PaperSize = [11.7 8.3]; % A4
print(fig,'-dpng','-r600','Lightbulb_DC_Alptal_Paper.png')
print(fig,'-dpdf','-r600','Lightbulb_DC_Alptal_Paper.pdf')

fig=figure(fig_num);fig_num = fig_num+1;
sp=subplot(1,1,1);
hold on
plot(Sky_Emissivity_Alptal_2006,LWR_Subcanopy_Observations_Alptal_2006./LWR_Atmosphere_Alptal_2006,'.','Color',EcstaticEmerald)
plot(Esky_Alp_06,LWE_ctrl_Alp_06,'.k')
plot(Esky_Alp_06,LWE_lbulb_Alp_06,'.r')
xlim([0.5 1.2])
ylim([0.8 1.8])
ylabel('Longwave Enhancement','FontSize',15,'FontWeight','bold')
xlabel('Effective Emissivity of the Sky','FontSize',15,'FontWeight','bold')
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2)
ctrl=plot(-1,-1,'MarkerEdgeColor','k','MarkerFaceColor','k',...
    'Marker','o','MarkerSize',10,'LineStyle','none');
lbulb=plot(-1,-1,'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'Marker','o','MarkerSize',10,'LineStyle','none');
obs=plot(-1,-1,'MarkerEdgeColor',EcstaticEmerald,'MarkerFaceColor',...
    EcstaticEmerald,'Marker','o','MarkerSize',10,'LineStyle','none');
legend([ctrl,lbulb,obs],'CTRL','CORR','OBS',...
    'location','northeast','orientation','vertical',...
    'FontSize',15,'FontWeight','bold')
legend('boxoff')
box on
set(sp,'Position',[0.11,0.123,0.85,0.86])
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 5];
fig.PaperSize = [6 5];
print(fig,'-dpng','-r600','Scatter_Esky_LWE_2006.png')
print(fig,'-dpdf','-r600','Scatter_Esky_LWE_2006.pdf')

toc