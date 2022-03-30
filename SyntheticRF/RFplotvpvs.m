clear all; clc; close all;

data=load('lowfreqvpvs1.mat');
data2=load('lowfreqvpvs2.mat');
data3=load('lowfreqvpvs3.mat');
data4=load('lowfreqsimplemod.mat');
plot(data4.tt(1:181),data4.cc4norm(1:181),'-b','linewidth',3); hold on
plot(data.tt(1:181),data.cc4norm(1:181),'--r','linewidth',3); hold on
plot(data2.tt(1:181),data2.cc4norm(1:181),'-k','linewidth',3); hold on
plot(data3.tt(1:181),data3.cc4norm(1:181),'--m','linewidth',3); hold off
ylim([-0.5 1]);
legend('simple model','model 1 (vp meningkat)','model 2 (vs menurun)','model 3 (vp vs menurun)');
set(gca,'fontsize',20);
% title('Synthetic Receiver Function');
xlabel('Time (s)');
ylabel('Amplitude');