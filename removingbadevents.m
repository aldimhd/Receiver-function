close all;
clear all;
clc
load Z05STNG
%Make plot dan remove the bad one

[x,y]=size(Z05S_header); %mengetahui ukuran data
%Butterworth(bandpass filter)
flow=0.1;
fhigh=2;

for i=1:y
    A = Z05S_A(:,i);
    fs= A.Vert.stats.sampling_rate;
    taxis = 1:1:8500;
    maxi = 8501;
    taxis1 = (0:(maxi-1))/fs;
    [B1,A1] = butter(2,[flow/(fs/2) fhigh/(fs/2)]); %hitung koefisien,B1= orde, A1= cut off
    figure(1)
        %Zero-phase digital filtering
        dat1 = filtfilt(B1,A1,detrend(Z05S_vert_ori(:,i)));
        d1=max(Z05S_vert_ori(:,i));
        d2=min(Z05S_vert_ori(:,i));
        d3=max(dat1);
        d4=min(dat1);
%  subplot(4,1,1),plot(Z05S_vert_ori(:,i));line([1000 1000],[d1 d2],'color','red');grid on %Vertika Ori
%  subplot(4,1,2),plot(detrend(dat1));line([1000 1000],[d3 d4],'color','red');grid on %vertikal ori tanpa noise
 if fs > 50
 subplot(2,1,1),plot(taxis1(1:4001),Z05S_LQTf(1000:5000,i));line([10 10],[1 -1],'color','red');grid on %Rf frekuensi 100 hZ
 title('Dekonvolusi Water Level')
 ylabel('Amplitudo (m)')
 xlabel('Time (s)')
 set(gca,'xTicklabel',[-10 -5 0 5 10 15 20 25 30])
 subplot(2,1,2),plot(taxis1(1:4001),Z05S_LQTt(:,i));line([10 10],[1 -1],'color','red');grid on %Rf Time
 title('Dekonvolusi Iterative')
 ylabel('Amplitudo (m)')
 xlabel('Time (s)')
 set(gca,'xTicklabel',[-10 -5 0 5 10 15 20 25 30])
 else
 subplot(2,1,1),plot(taxis1(1:2001), Z05S_LQTf(500:2500,i));line([10 10],[1 -1],'color','red');grid on %Rf frekuensi 50 hZ
 title('Dekonvolusi Water Level')
 ylabel('Amplitudo (m)')
 xlabel('Time (s)')
 set(gca,'xTicklabel',[-10 -5 0 5 10 15 20 25 30])
 subplot(2,1,2),plot(taxis1(1:2001),Z05S_LQTt(500:2500,i));line([10 10],[1 -1],'color','red');grid on %Rf Time 
 ylabel('Amplitudo (m)')
 xlabel('Time (s)')
 set(gca,'xTicklabel',[-10 -5 0 5 10 15 20 25 30])
 end
 
      title({sprintf('File Number %d, of %d , Fs= %d Hz',i,y,fs);'Dekonvolusi Iterative'})
%       title('Dekonvolusi Iterative')
      %sprintf menampilkan output data numeric dalam format string
             [xx,yy]=ginput(1);
             if xx<10
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

%membuat tempat
Z05S_vertp=zeroout(Z05S_vert,2);
Z05S_ewp=zeroout(Z05S_ew,2);
Z05S_nsp=zeroout(Z05S_ns,2);
Z05S_vert_orip=zeroout(Z05S_vert_ori,2);
Z05S_LQTfp=zeroout(Z05S_LQTf,2);
Z05S_LQTtp=zeroout(Z05S_LQTt,2);
Z05S_headerp=zeroout(Z05S_header,2);

[f1,f2]=size(Z05S_LQTfp);
% dum1=ME43pick;g
%zaxis=1:1:200

%calculate stack pick data in time

%stack all data in time
%repick data to move in 
d6=sum(Z05S_LQTfp,2);
Z05S_LQTfp_ts=normalization3(d6);

d7=sum(Z05S_LQTtp,2);
Z05S_LQTtp_ts=normalization3(d7);

figure(77);
if fs > 50
    plot(taxis1(1:4001),Z05S_LQTfp_ts(1000:5000),'red');hold on;plot(taxis1(1:4001),Z05S_LQTtp_ts,'black') 
    ylabel('Amplitudo')
    xlabel('Time (s)')
    set(gca,'xTicklabel',[-10 -5 0 5 10 15 20 25 30])
    legend('Dekonvolusi Water Level','Dekonvolusi Iterative')
else
    plot(taxis1(1:2001),Z05S_LQTfp_ts(500:2500),'red');hold on;plot(taxis1(1:2001),Z05S_LQTtp_ts(500:2500),'black')
    ylabel('Amplitudo')
    xlabel('Time (s)')
    set(gca,'xTicklabel',[-10 -5 0 5 10 15 20 25 30])
    legend('Dekonvolusi Water Level','Dekonvolusi Iterative')
end

% save Z05S3pick Z05S_vertp Z05S_vert_orip Z05S_nsp Z05S_ewp Z05S_LQTfp Z05S_LQTtp Z05S_headerp Z05S_LQTfp_ts Z05S_LQTtp_ts Z05S_A

%stack all data in depth
% for u=1:f2
%     
%     
%  [ ME43pickdepth(u,:), zaxis ] =depthconverter(dum1(1015:6015,u),0.04,100,500);
% end
% 
% 
% %stack data in depth
% d5=ME43pickdepth';
% normd5=normalization3(d5);
% 
% d6=sum(normd5,2);
% Z05S_LQTfp_ts=normalization3(d6);
%     %root stack
%     nthrootstackME43=rootstack(d5',2)'
%     nthrootstacknorME43=normalization3(nthrootstackME43)
% figure(5)
% 
% plot(Z05S_LQTfp_ts)
% figure (6)
% plot(nthrootstackME43)
