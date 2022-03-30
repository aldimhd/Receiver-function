clear all
close all
% cd ("E:\Kuliah\Semester 8\Tugas Akhir\Pengolahan data\Data TA")
file_out=('E:\Kuliah\Semester 8\Tugas Akhir\Pengolahan data\Data RF TA Kepulauan Mentawai\Siberut Island\Z05S\*Signal*Z05S*.mat');
file_out2=dir(file_out);
flow=0.1;
fhigh=1;
for n=1:length(file_out2)
    filename= file_out2.name
    A = importdata(filename)
    Z05S_A(:,n)= A;
    A.Vert.stats.sampling_rate
    
    fs= A.Vert.stats.sampling_rate;
end
[B1,A1] = butter(2,[flow/(fs/2) fhigh/(fs/2)]); %hitung koefisien
     alph=2.5;
     wl= 0.001;
          
            DT=0.02
            MINLAG=0
            MAXLAG=60
            TSHIFT=20 %waktu tunda 0s akan berada pada waktu tshift sesudahnya
            F0=2.5
            ITMAX=1000 %membatasi jumlah iterasi / lonjakan untuk ditambahkan
            MINDERR=0.0005 %hentikan pengulangan saat perubahan kesalahan dari penambahan lonjakan lain turun di bawah ambang batas ini
            ISVB=1
for n=1:length(file_out2)
    
        clc % menghapus layar di command window
       filename1=file_out2(n).name;
        A=importdata(filename1);
%         a1=rdsac('C:/rf/01.SAC') ; %HHT = SH wave
%         b1=rdsac('C:/rf/02.SAC') ; %HHQ = SV wave
%         c1=rdsac('C:/rf/03.SAC') ; %HHL = P wave
%          figure(2);plot(c1.d);title('dari code stream')
%          figure(3);plot(b1.d);title('Q dari code stream')
%         figure(4);plot(a1.d)
          
            A.vert=A.Vert.data;
            A.vert_ori=A.vert;
            A.East=A.East.data;
            A.North=A.North.data;
            A.vert=filtfilt(B1,A1,detrend(A.vert));
            A.East=filtfilt(B1,A1,detrend(A.East));
            A.North=filtfilt(B1,A1,detrend(A.North));
            
%             figure(1)
%             plot(A.vert(1500:6500));title('original vertical data')
%         dat=A.vert;
%         dat1 = filtfilt(B1,A1,dat);
%         dat1=detrend(dat1);
%         figure(1);plot(A.North);
%         figure(2);plot(dat1);
        
        baz=azimuth(A.Header(2),A.Header(1),A.Header(4),A.Header(5))
%         baz2=azimuth(A.Header(4),A.Header(5),A.Header(2),A.Header(1))
         iang=A.Header(10)
        
         datL=cosd(iang).*A.vert - sind(iang)*sind(baz).*A.East - sind(iang)*cosd(baz).*A.North;
         datQ=sind(iang).*A.vert + cosd(iang)*sind(baz).*A.East + cosd(iang)*cosd(baz).*A.North;
         datT=sind(baz).*A.North - cosd(baz)*A.East;
         datQ=-1*datQ;
         
         
%         datR2 = -A.East.*sind(baz) - A.North.*cosd(baz);
%         datT2 =  - A.East.*sind((baz-90)) - A.North.*cosd((baz-90));
%         datZ3 =  A.vert;
%         figure(11);plot(datL(1500:6500));title('L comp self calculated')
%          figure(12);plot(datQ(1500:6500));title('Q comp self calculated')
% %         figure(13);plot(datT(1500:7500));

%RF frequency domain, water level method
        [RFtmp3, RFtaxis3] = rf_minwl(datL,datT,wl,alph,fs,20);
%         figure(5)
%         plot(RFtmp3(1:8000));

            UIN=datT;
            WIN=datL;
        [RFd RMS]=makeRFitdecon(UIN,WIN,DT,...
				   MINLAG,MAXLAG,TSHIFT,F0,...
				   ITMAX,MINDERR,ISVB);

   
%      datZrot =  datZ3.*cosd(iang) + datR2.*sind(iang); %L = P wave 
%      datRrot =  - datZ3.*sind(iang) + datR2.*cosd(iang); % Q = SV
%             figure(7)
%            plot(datZrot(1500:6500),'black');title('dave P wave comp')
%             figure(8)
%             plot(datRrot(1500:6500),'blue');title('dave Rrot wave comp');hold on;plot(datQ(1500:6500),'*red')
%             figure(9)
%             plot(datR2(1500:6500));title('dave Radial')

%                [RFtmp4, RFtaxis4] = rf_minwl(datZrot,datRrot,wl,alph,fs,20);
%                [RFd2, RMS2]=makeRFitdecon(datRrot,datZrot,DT,...
% 				   MINLAG,MAXLAG,TSHIFT,F0,...
% 				   ITMAX,MINDERR,ISVB);
%          figure(66);plot(RFtmp3,'red');hold on;plot(RFtmp4,'black')
%           figure(67);plot(RFtaxis4,RFtmp4,'red');hold on; plot(RFd2)
%           figure(68);plot(RFd,'black')

          
          
%           nc=normalization3(b1.d);figure(71);plot(nc)
          nRFf=normalization3(RFtmp3);
%           nRFf2=normalization3(RFtmp4);
          nRFd=normalization3(RFd');
%           nRFd2=normalization3(RFd2');
%           figure(81);plot(nc,'red');hold on;plot(nRFf(1500:6500),'black');hold on;plot(nRFd(1500:6500),'blue')
%           figure(82);plot(nRFf(1500:6500),'black');hold on;plot(nRFd1500:6500),'red')
%          figure(6);plot(RFd,'red');hold on;plot(RFd2,'black')
        %calculate azimuth az=azimuth(stla,stlo,evla,avlo)
        %perform seismic rotation using BZ
        %pefform RF
        %perform rotation with incidence angle
        %perform RFkem
        Z05S_vert(:,n)=A.vert;
        Z05S_ns(:,n)=A.North;
        Z05S_ew(:,n)=A.East;
        Z05S_vert_ori(:,n)=A.vert_ori';
        Z05S_LQTf(:,n)=nRFf;
        Z05S_LQTt(:,n)=nRFd;
        Z05S_header(:,n)=A.Header;
      
       %  [[stlon], [stlat], [stel], [lat[j]], [lon[j]], [depth[j]], [mag[j]], [arc_distance], [ray_param], [angle]]
end

save Z05STNG Z05S_vert Z05S_vert_ori Z05S_ns Z05S_ew Z05S_LQTf Z05S_LQTt Z05S_header Z05S_A


