close all;
clear all;
clc
load Z05S2pick
%Make plot dan remove the bad one

[x,y]=size(Z05S_headerp); %mengetahui ukuran data
%Butterworth(bandpass filter)
flow=0.1;
fhigh=1;

for i=1:y
    A = Z05S_A(:,i);
    fs= A.Vert.stats.sampling_rate;
    [B1,A1] = butter(2,[flow/(fs/2) fhigh/(fs/2)]); %hitung koefisien,B1= orde, A1= cut off
    figure(1)
        %Zero-phase digital filtering
        dat1 = filtfilt(B1,A1,detrend(Z05S_vert_orip(:,i)));
        d1=max(Z05S_vert_orip(:,i));
        d2=min(Z05S_vert_orip(:,i));
        d3=max(dat1);
        d4=min(dat1);
 subplot(2,1,1),plot(Z05S_vert_ori(:,i));line([1000 1000],[d1 d2],'color','red');grid on %Vertika Ori
 subplot(2,1,2),plot(detrend(dat1));line([1000 1000],[d3 d4],'color','red');grid on %vertikal ori tanpa noise
%  subplot(2,1,1),plot(Z05S_LQTfp(1000:5000,i));line([1000 1000],[1 -1],'color','red');grid on %Rf frekuensi
%  ylabel('Amplitudo')
%  xlabel('ntps')
%  subplot(2,1,2),plot(Z05S_LQTtp(:,i));line([1000 1000],[1 -1],'color','red');grid on %Rf Time
%  ylabel('Amplitudo')
%  xlabel('ntps')
      title(sprintf('file number %d, of %d , fs= %d Hz',i,y,fs))
      %sprintf menampilkan output data numeric dalam format string
             [xx,yy]=ginput(1);
             if xx<1000
                 Z05S_vertp(:,i)=0;
                 Z05S_vert_orip(:,i)=0;
                 Z05S_ewp(:,i)=0;
                 Z05S_nsp(:,i)=0;
                 Z05S_LQTfp(:,i)=0;
                 Z05S_LQTtp(:,i)=0;
                 Z05S_headerp(:,i)=0;
                  
             else
                 Z05S_vertp(:,i)=Z05S_vertp(:,i);
                 Z05S_ewp(:,i)=Z05S_ewp(:,i);
                 Z05S_nsp(:,i)=Z05S_nsp(:,i);
                 Z05S_vert_orip(:,i)=Z05S_vert_orip(:,i);
                 Z05S_LQTfp(:,i)=Z05S_LQTfp(:,i);
                 Z05S_LQTtp(:,i)=Z05S_LQTtp(:,i);
                 Z05S_headerp(:,i)=Z05S_headerp(:,i);
                                             
             end
     
end


%membuat tempat
Z05S_vertp=zeroout(Z05S_vertp,2);
Z05S_ewp=zeroout(Z05S_ewp,2);
Z05S_nsp=zeroout(Z05S_nsp,2);
Z05S_vert_orip=zeroout(Z05S_vert_orip,2);
Z05S_LQTfp=zeroout(Z05S_LQTfp,2);
Z05S_LQTtp=zeroout(Z05S_LQTtp,2);
Z05S_headerp=zeroout(Z05S_headerp,2);

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
    plot(Z05S_LQTfp_ts(1000:5000),'red');hold on;plot(Z05S_LQTtp_ts,'black') 
else
    plot(Z05S_LQTfp_ts(1:5000),'red');hold on;plot(Z05S_LQTtp_ts,'black')
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
