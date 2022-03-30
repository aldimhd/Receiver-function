% Stack Align

close all
clear all

load Z35S2pick
load Z20S2pick
load Z05S2pick
load NGNSpick

load PPNGpick
load DSAOpick
load SLBU2pick
load MLKRpick
load BSATpick

% w=[ Z35S_LQTtp_ts(500:2500) Z20S_LQTtp_ts(500:2500) Z05S_LQTtp_ts(500:2500) NGNS_LQTtp_ts(500:2500)] %siberut iterative
% w=[ Z35S_LQTfp_ts(1000:5000) Z20S_LQTfp_ts(1000:5000) Z05S_LQTfp_ts(30:4030) NGNS_LQTfp_ts(1000:5000)] %siberut water level 100 hz

% w=[ PPNG_LQTtp_ts(500:2500) DSAO_LQTtp_ts(500:2500) SLBU_LQTtp_ts(500:2500) MLKR_LQTtp_ts(500:2500) BSAT_LQTtp_ts(500:2500) ]
w=[ PPNG_LQTfp_ts(1560:3560) DSAO_LQTfp_ts(500:2500) SLBU_LQTfp_ts(500:2500) MLKR_LQTtp_ts(500:2500) BSAT_LQTfp_ts(500:2500)]
maxi= 2001
fs= 50
% maxi= 4001
% fs= 100
taxis= (0:(maxi-1))/fs;

w=normalization3(w);
% x=[0 24.97 54.94 81.33] %Siberut
x=[0 46.7  96.66  124.66 140.49 142.39]; %Sipora dan pagai

figure(1)
wigbC4(w,0.6437,x,taxis);
xlim([-20 100]) %siberut
xlim([-20 170])  %sipora dan pagai
ylabel('Time (s)')
xlabel('Distance (km)')
set(gca,'yTicklabel',[-10 -5 0 5 10 15 20 25 30])

% Pulau Siberut
% annotation('textbox',[0.22 0.89 .1 .1],'String','Z35S','FitBoxToText','on');
% annotation('textbox',[0.38 0.89 .1 .1],'String','Z20S','FitBoxToText','on');
% annotation('textbox',[0.57 0.89 .1 .1],'String','Z05S','FitBoxToText','on');
% annotation('textbox',[0.74 0.89 .1 .1],'String','NGNS','FitBoxToText','on');

% Pulau Sipora&pagai
annotation('textbox',[0.17 0.89 .1 .1],'String','PPNG','FitBoxToText','on');
annotation('textbox',[0.35 0.825 .1 .1],'String','DSAO','FitBoxToText','on');
annotation('textbox',[0.55 0.89 .1 .1],'String','SLBU','FitBoxToText','on');
annotation('textbox',[0.67 0.825 .1 .1],'String','MLKR','FitBoxToText','on');
annotation('textbox',[0.74 0.89 .1 .1],'String','BSAT','FitBoxToText','on');


% sta = 'Z40S Z35S Z20S Z05S'
% str= sprintf('Station %s',sta)
%  help wigbC4

% for u=1:f2
%     
%     
%  [ BSATpick(u,:), zaxis ] =depthconverter3(dum1(1015:6015,u),0.04,100,500);
% end
% 
% %stack data in depth
% d5=ME43pickdepth';
% normd5=normalization3(d5);
% 
% d6=sum(normd5,2);
% PPNG_LQTfp_ts=normalization3(d6);
%     %root stack
%     nthrootstackME43=rootstack(d5',2)'
%     nthrootstacknorME43=normalization3(nthrootstackME43)
% figure(5)
% 
% plot(PPNG_LQTfp_ts)
% figure (6)
% plot(nthrootstackME43)
