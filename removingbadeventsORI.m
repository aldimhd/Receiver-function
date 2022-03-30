close all
clear all
clc
load Z05S

%MAke plot dan remove the bad one

[x y]=size(Z05S_header)
flow=0.1;
fhigh=1;
fs=50;
[B1,A1] = butter(2,[flow/(fs/2) fhigh/(fs/2)]);
for i=1:y
        
    figure(1)
        
        dat1 = filtfilt(B1,A1,detrend(Z05S_vert_ori(:,i)));
        dat2 = normalization3(dat1)
        d1=max(Z05S_vert_ori(:,i));
        d2=min(Z05S_vert_ori(:,i));
        d3=max(dat1);
        d4=min(dat1);
 subplot(2,1,1),plot(Z05S_vert_ori(:,i));line([1000 1000],[d1 d2],'color','red');grid on
%  subplot(2,1,1),plot(detrend(dat1));line([1000 1000],[d3 d4],'color','red');grid on
 ylabel('Amplitudo (m)')
 xlabel('npts')
%  set(gca,'xTicklabel',[-10 -5 0 5 10 15 20 25 30])
 subplot(2,1,2),plot(detrend(dat1));line([1000 1000],[d3 d4],'color','red');grid on
 ylabel('Amplitudo (m)')
 xlabel('npts')
%  set(gca,'xTicklabel',[-10 -5 0 5 10 15 20 25 30])
%  subplot(4,1,3),plot(Z05S_LQTf(1500:6500,i));line([500 500],[1 -1],'color','red');grid on
%  subplot(4,1,4),plot(Z05S_LQTt(:,i));line([500 500],[1 -1],'color','red');grid on

%       title(sprintf('file number %d, of %d',i,y))
             [xx yy]=ginput(1);
             if xx<1000
                 Z05S_vert(:,i)=0;
                 Z05S_vert_ori(:,i)=0;
                 Z05S_ew(:,i)=0;
                 Z05S_ns(:,i)=0;
                 Z05S_LQTf(:,i)=0;
                 Z05S_LQTt(:,i)=0;
                Z05S_header(:,i)=0;
             
             else
                 Z05S_vert(:,i)=Z05S_vert(:,i);
                 Z05S_ew(:,i)=Z05S_ew(:,i);
                 Z05S_ns(:,i)=Z05S_ns(:,i);
                 Z05S_vert_ori(:,i)=Z05S_vert_ori(:,i);
                 Z05S_LQTf(:,i)=Z05S_LQTf(:,i);
                 Z05S_LQTt(:,i)=Z05S_LQTt(:,i);
                 Z05S_header(:,i)=Z05S_header(:,i);
                             
             end
     end


%Z05Spick=zeroout(Z05Sall,2);
  Z05S_vertp=zeroout(Z05S_vert,2);
  Z05S_ewp=zeroout(Z05S_ew,2);
  Z05S_nsp=zeroout(Z05S_ns,2);
  Z05S_vert_orip=zeroout(Z05S_vert_ori,2);
  Z05S_LQTfp=zeroout(Z05S_LQTf,2);
  Z05S_LQTtp=zeroout(Z05S_LQTt,2);
  Z05S_headerp=zeroout(Z05S_header,2);

[f1 f2]=size(Z05S_LQTfp);
% dum1=Z05Spick;g
%zaxis=1:1:200

%calculate stack pick data in time

%stack all data in time
%repick data to move in 
d6=sum(Z05S_LQTfp,2);
Z05S_LQTfp_ts=normalization3(d6);

d7=sum(Z05S_LQTtp,2);
Z05S_LQTtp_ts=normalization3(d7);



figure(77);
plot(Z05S_LQTfp_ts(1:5000),'red');hold on;plot(Z05S_LQTtp_ts(1:5000),'black')


% save Z05Spick Z05S_vertp Z05S_vert_orip Z05S_nsp Z05S_ewp Z05S_LQTfp Z05S_LQTtp Z05S_headerp Z05S_LQTfp_ts Z05S_LQTtp_ts
% %stack all data in depth
% for u=1:f2
%     
%     
%  [ Z05Spickdepth(u,:), zaxis ] =depthconverter(dum1(1015:6015,u),0.04,100,500);
% end
% 
% 
% %stack data in depth
% d5=Z05Spickdepth';
% normd5=normalization3(d5);
% 
% d6=sum(normd5,2);
% Z05S_LQTfp_ts=normalization3(d6);
%     %root stack
%     nthrootstackZ05S=rootstack(d5',2)'
%     nthrootstacknorZ05S=normalization3(nthrootstackZ05S)
% figure(5)
% 
% plot(Z05S_LQTfp_ts)
% figure (6)
% plot(nthrootstackZ05S)
