clear all; clc; close all;

data=dlmread('Aldi2_No_LVZ.txt', ' ', 2, 0);
data(:,3)=[];
data(:,3)=[];
% data(:,3)=[];
% data(:,4)=[];
data(1:2:end-1,:)=[];       %delete odd rows

depth=data(:,1);
vp=data(:,2);
vs=data(:,3);

plot(vp,depth,'linewidth',3); hold on
plot(vs,depth,'linewidth',3); hold off
set(gca,'ydir','reverse','fontsize',20);
% ylim([0 50000]);
hAx=gca;                    % handle to current axes
hAx.YAxis.Exponent=3;       % don't use exponent
% hAx.YAxis.Format='%d';    % use integer format, no decimal points
hAx.XAxis.Exponent=3;
% title('Velocity Model');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
legend('Vp','Vs');