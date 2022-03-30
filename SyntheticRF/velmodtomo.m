clear all; clc; close all;

data=dlmread('velmod_ME03tomo.txt', ' ', 2, 0);
data(:,3)=[];
data(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data(1:2:end-1,:)=[];       %delete odd rows
depth=data(:,1);
vp=data(:,2);
vs=data(:,3);

data2=dlmread('velmod_ME10tomo.txt', ' ', 2, 0);
data2(:,3)=[];
data2(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data2(1:2:end-1,:)=[];       %delete odd rows
depth2=data2(:,1);
vp2=data2(:,2);
vs2=data2(:,3);

data3=dlmread('velmod_ME16tomo.txt', ' ', 2, 0);
data3(:,3)=[];
data3(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data3(1:2:end-1,:)=[];       %delete odd rows
depth3=data3(:,1);
vp3=data3(:,2);
vs3=data3(:,3);

data4=dlmread('velmod_ME17tomo.txt', ' ', 2, 0);
data4(:,3)=[];
data4(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data4(1:2:end-1,:)=[];       %delete odd rows
depth4=data4(:,1);
vp4=data4(:,2);
vs4=data4(:,3);

data5=dlmread('velmod_ME20tomo.txt', ' ', 2, 0);
data5(:,3)=[];
data5(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data5(1:2:end-1,:)=[];       %delete odd rows
depth5=data5(:,1);
vp5=data5(:,2);
vs5=data5(:,3);

data6=dlmread('velmod_ME24tomo.txt', ' ', 2, 0);
data6(:,3)=[];
data6(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data6(1:2:end-1,:)=[];       %delete odd rows
depth6=data6(:,1);
vp6=data6(:,2);
vs6=data6(:,3);

data7=dlmread('velmod_ME25tomo.txt', ' ', 2, 0);
data7(:,3)=[];
data7(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data7(1:2:end-1,:)=[];       %delete odd rows
depth7=data7(:,1);
vp7=data7(:,2);
vs7=data7(:,3);

data8=dlmread('velmod_ME28tomo.txt', ' ', 2, 0);
data8(:,3)=[];
data8(:,3)=[];
data8(:,3)=[];
% data(:,4)=[];
data8(1:2:end-1,:)=[];       %delete odd rows
depth8=data8(:,1);
vp8=data8(:,2);
vs8=data8(:,3);

data9=dlmread('velmod_ME32tomo.txt', ' ', 2, 0);
data9(:,3)=[];
data9(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data9(1:2:end-1,:)=[];       %delete odd rows
depth9=data9(:,1);
vp9=data9(:,2);
vs9=data9(:,3);

data10=dlmread('velmod_ME36tomo.txt', ' ', 2, 0);
data10(:,3)=[];
data10(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data10(1:2:end-1,:)=[];       %delete odd rows
depth2=data10(:,1);
vp10=data10(:,2);
vs10=data10(:,3);

data11=dlmread('velmod_ME39tomo.txt', ' ', 2, 0);
data11(:,3)=[];
data11(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data11(1:2:end-1,:)=[];       %delete odd rows
depth11=data11(:,1);
vp11=data11(:,2);
vs11=data11(:,3);

data12=dlmread('velmod_ME42tomo.txt', ' ', 2, 0);
data12(:,3)=[];
data12(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data12(1:2:end-1,:)=[];       %delete odd rows
depth12=data12(:,1);
vp12=data12(:,2);
vs12=data12(:,3);

data13=dlmread('velmod_ME53tomo.txt', ' ', 2, 0);
data13(:,3)=[];
data13(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data13(1:2:end-1,:)=[];       %delete odd rows
depth13=data13(:,1);
vp13=data13(:,2);
vs13=data13(:,3);

plot(vp,depth,'linewidth',3,'b'); hold on
plot(vp2,depth,'linewidth',3,'k'); hold on
plot(vp3,depth,'linewidth',3,'r'); hold on
plot(vp4,depth,'linewidth',3,'g'); hold on
plot(vp5,depth,'linewidth',3,'y'); hold on
plot(vp6,depth,'linewidth',3,'c'); hold on
plot(vp7,depth,'linewidth',3,'m'); hold on
plot(vp8,depth,'linewidth',3); hold on
plot(vp9,depth,'linewidth',3); hold on
plot(vp10,depth,'linewidth',3); hold on
plot(vp11,depth,'linewidth',3); hold on
plot(vp12,depth,'linewidth',3); hold on
plot(vp13,depth,'linewidth',3); hold on
plot(vs,depth,'linewidth',3,'--b'); hold on
plot(vs2,depth,'linewidth',3,'--k'); hold on
plot(vs3,depth,'linewidth',3,'--r'); hold on
plot(vs4,depth,'linewidth',3,'--g'); hold on
plot(vs5,depth,'linewidth',3,'--y'); hold on
plot(vs6,depth,'linewidth',3,'--c'); hold on
plot(vs7,depth,'linewidth',3,'--m'); hold on
plot(vs8,depth,'linewidth',3); hold on
plot(vs9,depth,'linewidth',3); hold on
plot(vs10,depth,'linewidth',3); hold on
plot(vs11,depth,'linewidth',3); hold on
plot(vs12,depth,'linewidth',3); hold on
plot(vs13,depth,'linewidth',3); hold off

set(gca,'ydir','reverse','fontsize',20);
ylim([0 50000]);
hAx=gca;                    % handle to current axes
hAx.YAxis.Exponent=3;       % don't use exponent
% hAx.YAxis.Format='%d';    % use integer format, no decimal points
hAx.XAxis.Exponent=3;
% title('Velocity Model');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
% legend('Vp','Vs');