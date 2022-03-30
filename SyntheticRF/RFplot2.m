clear all; clc; close all;
data=load('lowfreq.mat');
data2=load('hifreq.mat');
cd ('E:\Kuliah\Semester 8\Tugas Akhir\Pengolahan data\Data RF TA Kepulauan Mentawai')
load Z20S2pick

maxi=length(Z20S_LQTfp_ts(:,1));
fs= 100
taxis=(-2000:(maxi-2000))/fs;

nr = Z20S_LQTfp_ts(2101:2701);
nfw = data.cc4norm(9:63);
% nfw2 = data.cc4norm(17:113);

nrp = normalization3(nr);
nfwp = normalization3((nfw)');
% nfwp2 = normalization3((nfw)');

nfwp3= nfwp-0.4;
plot(taxis(2101:2701),nrp(:),'black','linewidth',3);hold on
plot(data.tt(9:63),nfwp3(:),'linewidth',3); hold on

% title('Synthetic Receiver Function', 'fontsize',20);
ylim([-1 0.6]);
annotation('textbox',[0.15 0.79 .1 .1],'String','Stasiun Z20S','FitBoxToText','on');
legend('-Hz','1 Hz','4 Hz');
xlabel('Time (s)');
ylabel('Amplitude');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load BSATpick
% 
% maxi=length(BSAT_LQTfp_ts(:,1));
% fs= 50
% taxis=(-1000:(maxi-1000))/fs;
% 
% nr = BSAT_LQTfp_ts(1055:1401);
% nfw = data.cc4norm(9:63);
% 
% nrp = normalization3(nr);
% nfwp = normalization3((nfw)');
% 
% nfwp3= nfwp-0.7;
% plot(taxis(1055:1401),nrp(:),'black','linewidth',3);hold on
% plot(data.tt(9:63),nfwp3(:),'linewidth',3); hold on
% 
% ylim([-1 0.4]);
% annotation('textbox',[0.15 0.79 .1 .1],'String','Stasiun BSAT','FitBoxToText','on');
% legend('-Hz','1 Hz','4 Hz');
% xlabel('Time (s)');
% ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%  DIMULAI 1 SEKON  %%%%%%%%%%%%%%%%%
% nr = LT34_LQTfp_ts(2100:3000);
% nfw = data.cc4norm(9:81);
% nrp = normalization3(nr);
% nfwp = normalization3((nfw)');
% % plot(nfwp)
% plot(taxis(2101:3001),nrp(:),'black','linewidth',3);hold on
% plot(data.tt(9:81),nfwp(:),'linewidth',3); hold on
% % plot(data2.tt(11:71),data2.cc4norm(11:71),'linewidth',3); hold off
% ylim([-1 1]);
% legend('real data','1 Hz', 'fontsize',20);
% title('Synthetic Receiver Function', 'fontsize',20);
% xlabel('Time (s)');
% ylabel('Amplitude');
% % set(gca,'fontsize',20);