%%Zhu kanamori%%
%zhu and kanamori data stack
clear all
close all
clc
% cd c:/RF9

%Header = [[stlon], [stlat], [stel], [lat[j]], [lon[j]], [depth[j]], [mag[j]], [arc_distance], [ray_param], [angle]]
load PRKBpick
delay=20;

% average crustal p-wave velocity to use
vp=6.5;

% crustal thickness estimate min and max
Hthick1=22; Hthick2=40; 

% Vp/Vs min and max to search over
vpvs1=1.70; vpvs2=2.20; 

% number of grid points to use in grid search
numpts=55;

flow = 1/20;
fhigh = 2;
flow2 = .025;
fhigh2 = .5;
fs=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% no modifications should be needed below this line.

tstart=delay*fs;
vpvs = (vpvs1:(vpvs2-vpvs1)/(numpts-1):vpvs2);
Hthick = Hthick1:(Hthick2-Hthick1)/(numpts-1):Hthick2;
solgrid = zeros(numpts,numpts);
stdgrid = solgrid;
                 [B1,A1] = butter(2,[flow/(fs/2) fhigh/(fs/2)]);
                 [B2,A2] = butter(2,[flow2/(fs/2) fhigh2/(fs/2)]);
                 
[r1, r2]=size(PRKB_nsp);      

  for n = 1:r2,
      
            dat = detrend(PRKB_LQTtp(:,n)); %ambil Rf time domain
            adist = PRKB_headerp(8,n);  %arc distance

                 dat1 = filtfilt(B1,A1,dat);
dat1=detrend(dat1);
                 dat2 = filtfilt(B2,A2,dat);
dat2=detrend(dat2);
         rayp(n)=sind(PRKB_headerp(10,n))/vp; %menghitung raypat dengan input angle
           adists(n) = adist;
           thissta(:,n) = dat1(tstart:min(tstart+2500,length(dat1)));
           thisstalf(:,n) = dat2(tstart:min(tstart+2500,length(dat1)));
thissta(:,n) = thissta(:,n)./sqrt(sum(thissta(:,n).^2));
thisstalf(:,n) = thisstalf(:,n)./sqrt(sum(thisstalf(:,n).^2));

%end of loop through files
  end    
  
  maxi=length(thissta(:,1));
	   taxis=(0:(maxi-1))/fs;
           thissta(maxi,:) = thissta(maxi,:).*0;
           thisstalf(maxi,:) = thisstalf(maxi,:).*0;

            [Y,I] = sort(rayp);
[M,N] = size(thissta);
thissta2f= thissta(:,I);
thisstalf= thisstalf(:,I);

for m = 1:numpts,
 for n = 1:numpts,
            VS = vp / vpvs(m);
            tstimes = Hthick(n) * ( VS^-2 - Y.^2).^.5;
            tptimes = Hthick(n) * ( vp^-2 - Y.^2).^.5;
            tps = tstimes - tptimes;
            tpsi = min(ceil((tps)*fs),maxi); %ceil pembulatan ke plus tak hingga
            t2p1s = tstimes+tptimes;
            t2p1si = min(ceil((t2p1s)*fs),maxi);
            t2s1p = 2.*tstimes;
            t2s1pi = min(ceil((t2s1p)*fs),maxi);
            for nn = 1:length(Y),
             solvector(nn) = .5*thissta2f(tpsi(nn),nn) + .3*thisstalf(t2p1si(nn),nn) - .2*thisstalf(t2s1pi(nn),nn);
            end
            solgrid(m,n) = sum(solvector);
            stdgrid(m,n) = std(solvector);

 end
end

figure(1); clf
[cc,hh] =contour(Hthick,vpvs,solgrid./max(max(solgrid)),[ 0 .25 .5 .75 .95 .99]);
set(hh,'linewidth',2)
clabel(cc,hh)
grid on
xlabel('thickness (km)')
ylabel('Vp/Vs ratio')

%ambil nilai max dari grid
[Yn, Ih] = max(max(solgrid)); %Yn= nilai maksimum, Ih= urutan kolom
[Ym, Iv] = max(max(solgrid')); %%Ym= nilai maksimum, Iv= urutan kolom setelah di transpose

% this may need to be reversed
[Fx,Fy] = gradient(solgrid,(Hthick2-Hthick1)/(numpts-1),(vpvs2-vpvs1)/(numpts-1));
Fx2 = gradient(Fx,(Hthick2-Hthick1)/(numpts-1));
[Fxx, Fy2] = gradient(Fy,(vpvs2-vpvs1)/(numpts-1));

stdH = abs(2*stdgrid(Iv,Ih)^0.5 / Fx2(Iv,Ih)).^0.5;
stdV = abs(2*stdgrid(Iv,Ih)^0.5 / Fy2(Iv,Ih)).^0.5;

%disp([ vp Hthick(Ih) stdH vpvs(Iv) stdV ])
 h=Hthick(Ih);
 VpVs=vpvs(Iv);
 herr=stdH;
 vpvserr=stdV;

            tstimes =  Hthick(Ih)* ( (vp/vpvs(Iv))^-2 - Y.^2).^.5;
            tptimes =  Hthick(Ih)* ( vp^-2 - Y.^2).^.5;
            tps = tstimes - tptimes;
            t2p1s = tstimes + tptimes;
            hhh = 2*tstimes;
            kkk = mean(hhh);

     figure(2); clf
     wigbC4(thissta2f,1.5,adists(I),taxis) %tanda (I)= mengurutkan dari yang terbesar hingga terkecil
axis([min(adists(I))+10 max(adists(I))+5 -2 25])
figure(2)
plot(adists(I),tps,'k-','linewidth',2)
plot(adists(I),t2p1s,'b-','linewidth',2)
plot(adists(I),2*tstimes,'r-','linewidth',2)
     ylabel('time (s)')
     xlabel('epicentral distance (deg)')
     title ('Stasiun PRKB')
   saveas(gcf,'PRKB_RF_witharcdistance.png') 
sta = 'PRKB'
figure(1);
xlabel('crustal thickness (km)')
     title(sprintf('Station %s: Vp=%2.2f, VpVs=%1.2f+-%1.2f, H=%2.2f+-%1.2f',sta,vp,VpVs,vpvserr,h,herr))

     figure(2);
     title(sprintf('Station %s',sta))
     saveas(gcf,'PRKBstackcontourzhukanaa.png')
     







