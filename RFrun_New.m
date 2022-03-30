clear all   %menghapus data di memori Matlab
close all   %menghapus semua gambar yang tampil sebelumnya

file_out=('E:\Kuliah\Semester 8\Tugas Akhir\Pengolahan data\Data RF TA Kepulauan Mentawai\Sipura Island\SOBS\*Signal*SOBS*.mat');
file_out2=dir(file_out); %untuk mengetahui file apa saja yang ada di current directory
%Butterworth(bandpass filter)
flow=0.1;
fhigh=2;
fs=50;
[B1,A1] = butter(2,[flow/(fs/2) fhigh/(fs/2)]); %hitung koefisien,B1= orde, A1= cut off
% freqz(B1,A1)
     alph=4; %lebar gausian filter pada water level
     wl=0.01; %fraksi water level
     DT=0.02; %sample interval (s)
     MINLAG=0; %Min rf lag time  (usually 0 for p->s, -ve for s->p) 
     MAXLAG=60;%Max rf lag time  (usually +ve for p->s, 0 for s->p)
     TSHIFT=20; %Time until beginning of receiver function (s)
     F0=2.5; %width of gaussian filter
     ITMAX=1000; %membatasi jumlah iterasi
     MINDERR=0.0005;%Min change in error required for stopping iterations
     ISVB=1;%flag for verbose output

for n=1:length(file_out2)
    clc % menghapus layar di command window
    filename1=file_out2(n).name;
     A=importdata(filename1);
        %sort data
        A.vert=A.Vert.data;
        A.vert_ori=A.vert;
        A.East=A.East.data;
        A.North=A.North.data;
        %Zero-phase digital filtering
        A.vert=filtfilt(B1,A1,detrend(A.vert));
        A.East=filtfilt(B1,A1,detrend(A.East));
        A.North=filtfilt(B1,A1,detrend(A.North));
%         subplot(3,1,1),plot(A.vert);title('vertical data');grid on
%         subplot(3,1,2),plot(A.East);title('timur data');grid on
%         subplot(3,1,3),plot(A.North);title('utara data');grid on
    %menghitung azimut 
    A.Header(10)
    baz = azimuth(A.Header(2),A.Header(1),A.Header(4),A.Header(5)); %back azimut
    iang = A.Header(10); %sudut datang yang diukur dari arah vertikal
    
    %Rotasi
    datL=cosd(iang).*A.vert - sind(iang)*sind(baz).*A.East - sind(iang)*cosd(baz).*A.North; %Selaras dengan arah perambatan gelombang P
    datQ=sind(iang).*A.vert + cosd(iang)*sind(baz).*A.East + cosd(iang)*cosd(baz).*A.North; %Selaras dengan arah gerakan fase SV
    datT=sind(baz).*A.North - cosd(baz)*A.East; %Disejajarkan dengan arah gerakan fase SH
    datQ=-1*datQ; %
  
%     figure(11);plot(datL(1500:6500));title('L comp self calculated')
%     figure(12);plot(datQ(1500:6500));title('Q comp self calculated')
%     figure(13);plot(datT(1500:7500));

%RF frequency domain, water level method
    % [ recfun, taxis ] = rf_min(z,r,fs)
    % z,r are vertical and radial-horizontal component data(or radial and tangential data)
    % fs is sampling frequency
    [RFtmp3, RFtaxis3] = rf_minwl(datL,datQ,wl,alph,fs,20); 
%         figure(5)
%         plot(RFtmp3(1:8000));
        
        UIN=datQ; %numerator (radial for PdS)
        WIN=datL; %denominator (vertical component for PdS)
    [RFd, RMS]=makeRFitdecon(UIN,WIN,DT,...
               MINLAG,MAXLAG,TSHIFT,F0,...
               ITMAX,MINDERR,ISVB)
        
     %Normalisasi
     nRFf=normalization3(RFtmp3);
     nRFd=normalization3(RFd');
%      figure(6)
%      plot(nRFf(1:8000));
        
     %Staking
     SOBS_vert(:,n)=A.vert;
     SOBS_ns(:,n)=A.North;
     SOBS_ew(:,n)=A.East;
     SOBS_vert_ori(:,n)=A.vert_ori';
     SOBS_LQTf(:,n)=nRFf;
     SOBS_LQTt(:,n)=nRFd;
     SOBS_header(:,n)=A.Header;
     
   %  [[stlon], [stlat], [stel], [lat[j]], [lon[j]], [depth[j]], [mag[j]], [arc_distance], [ray_param], [angle]]
end
    
save SOBS SOBS_vert SOBS_vert_ori SOBS_ns SOBS_ew SOBS_LQTf SOBS_LQTt SOBS_header
