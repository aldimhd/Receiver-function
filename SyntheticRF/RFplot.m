clear all; clc; close all;

data=load('lowfreqLT34lio.mat');
% data2=load('hifreq.mat');
plot(data.tt(1:181)*100,data.cc4norm(1:181),'linewidth',3); hold on
% plot(data2.tt(1:481)*100,data2.cc4norm(1:481),'linewidth',3); hold off
ylim([-0.5 1]);
legend('1.5 Hz','4 Hz');
% title('Synthetic Receiver Function');
xlabel('Time (s)');
ylabel('Amplitude');
set(gca,'fontsize',20);