%%Zhu kanamori%%
%zhu and kanamori data stack
clear all
close all
clc

%Header = [[stlon], [stlat], [stel], [lat[j]], [lon[j]], [depth[j]], [mag[j]], [arc_distance], [ray_param], [angle]]
load SMYApick
delay=20;

% average crustal p-wave velocity to use
vp=5.50; 

% crustal thickness estimate min and max
% Hthick1=13.67; Hthick2=35.8; %Water level
Hthick1=13; Hthick2=30; %Iterative

% Vp/Vs min and max to search over
vpvs1=1.55; vpvs2=2.2; %water level
% vpvs1=1.56; vpvs2=2.2; %Iterative

% number of grid points to use in grid search
numpts=55;

flow = 1/20;
fhigh = 1;
flow2 = .025;
fhigh2 = .5;
[r1, r2]=size(SMYA_nsp);%NortSouth data
% for i=1:r2
%     A = SMYA_A(:,i);
%     fs= A.Vert.stats.sampling_rate;
% end
fs = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% no modifications should be needed below this line.

tstart=delay*fs;
vpvs = (vpvs1:(vpvs2-vpvs1)/(numpts-1):vpvs2);
Hthick = Hthick1:(Hthick2-Hthick1)/(numpts-1):Hthick2;
solgrid = zeros(numpts,numpts);
stdgrid = solgrid;
                 [B1,A1] = butter(2,[flow/(fs/2) fhigh/(fs/2)]);
                 [B2,A2] = butter(2,[flow2/(fs/2) fhigh2/(fs/2)]);

  for n = 1:r2,

%             dat = detrend(SMYA_LQTfp(1000:8500,n)); % khusus fs 50 hz Rf Water level 
            dat = detrend(SMYA_LQTfp(:,n)); %Rf time domain/Water level
            adist = SMYA_headerp(8,n); %arc disatance

                 dat1 = filtfilt(B1,A1,dat);
dat1=detrend(dat1);
                 dat2 = filtfilt(B2,A2,dat);
dat2=detrend(dat2);
         %   iang = asin(raypS(1)*rvel);
         rayp(n)=sind(SMYA_headerp(10,n))/vp;
%           rayp(n) = (SMYA_headerp(9,n))/111.2;
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
            tpsi = min(ceil((tps)*fs),maxi);
            t2p1s = tstimes+tptimes;
            t2p1si = min(ceil((t2p1s)*fs),maxi);
            t2s1p = 2.*tstimes;
            t2s1pi = min(ceil((t2s1p)*fs),maxi);
            for nn = 1:length(Y),
             solvector(nn) = .7*thissta2f(tpsi(nn),nn) + .2*thisstalf(t2p1si(nn),nn) - .1*thisstalf(t2s1pi(nn),nn);
            end
            solgrid(m,n) = sum(solvector);
            stdgrid(m,n) = std(solvector);
 end
end

figure(21); clf
[cc,hh] =contour(Hthick,vpvs,solgrid./max(max(solgrid)),[ 0 .25 .5 .75 .95 .99]);
set(hh,'linewidth',2)
clabel(cc,hh)
grid on
xlabel('thickness (km)')
ylabel('Vp/Vs ratio')
%title('station SPA')
[Yn, Ih] = max(max(solgrid));
[Ym, Iv] = max(max(solgrid'));
% this may need to be reversed
[Fx,Fy] = gradient(solgrid,(Hthick2-Hthick1)/(numpts-1),(vpvs2-vpvs1)/(numpts-1));
Fx2 = gradient(Fx,(Hthick2-Hthick1)/(numpts-1));
[Fxx, Fy2] = gradient(Fy,(vpvs2-vpvs1)/(numpts-1));

stdH = abs(2*stdgrid(Iv,Ih)^.5 / Fx2(Iv,Ih)).^.5;
stdV = abs(2*stdgrid(Iv,Ih)^.5 / Fy2(Iv,Ih)).^.5;

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
     wigbC4(thissta2f,1.5,adists(I),taxis)
axis([min(adists(I))+-5 max(adists(I))+5 0 25])
figure(2)
ax = plot(adists(I),tps,'k-','linewidth',2,'Displayname','tPs');
legend(ax,'Location','southeast')
plot(adists(I),t2p1s,'b-','linewidth',2,'Displayname','tPpPs');
plot(adists(I),2*tstimes,'r-','linewidth',2,'Displayname','tPpSs+PsPs')
     ylabel('time (s)')
     xlabel('epicentral distance (deg)')
     title ('Stasiun SMYA')
%    saveas(gcf,'SMYA_RF_witharcdistance.png') 
sta = 'SMYA'
figure(21);
xlabel('H(km)')
     vp = 5.80; 
     title(sprintf('Station %s: Vp=%2.2f, VpVs=%1.2f+-%1.2f, H=%2.2f+-%1.2f',sta,vp,VpVs,vpvserr,h,herr))
     figure(2);
     title(sprintf('Station %s',sta))
%      saveas(gcf,'SMYAstackcontourzhukanaa.png')

%  figure(91);clf
AA=[rayp;thissta2f];
[w1, w2]=size(AA);
dzummy=zeros(w1-1,w2);

%calculate n
for g=1:8
   for f=1:w2
       if(rayp(f)>(g*0.005+0.035) && rayp(f)<(g*0.005+0.04)),
       dumx(f)=1;
       else
       dumx(f)=0;
       end
   end
    nn(g)=sum(dumx);
    clear dumx
end

for ii=1:8
    for jj=1:w2
        if ((AA(1,jj)>(ii*0.005+0.035)) && (AA(1,jj)<(ii*0.005+0.04))),
        dzummy(:,jj)=AA(2:w1,jj);
        else
        dzummy(:,jj)=0;
        end
    end
     
        newthis(:,ii)=detrend(sum(dzummy,2));
        %newthis(:,ii)=normalization(newthis(:,ii)');
end

for t=1:8
    newthis(:,t)=(1/nn(t).*newthis(:,t));
end

testa=newthis';
%testa=normalization3(testa);
testa=testa';
testa=detrend(testa)
%newthis=newthis+0.000001;
%newthis=normalization(newthis);
%newthis(:,3)=newthis(:,3);
angdis2=(0.0425:0.005:0.0775)
% wigbC4(testa,0.9,angdis2,taxis);
% figure(91);
% plot(rayp(I),tps,'g-','linewidth',2)
% plot(rayp(I),t2p1s,'r-','linewidth',2)
% plot(rayp(I),2*tstimes,'b-','linewidth',2)
%      ylabel('time (s)')
%      xlabel('ray parameter (deg)')
% title('SMYA stack bin')
% saveas(gcf,'SMYAstackdepthwithrayparameter.png')

% ave = mean2(tps);
% perhitungan stack data by azimuth termasuk juga zhu and kanamori method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% normd5=normalization3(SMYAtPSV);
% d6=sum(normd5,2);
% SMYAtPSVazimfilt_time_stack=normalization3(d6);
% figure(791);
% plot(SMYAtPSVazimfilt_time_stack);title('SMYA azimuth time stack')
% 
% 
% 
% [g1 g2]=size(SMYAtPSV);
% dum1=SMYAtPSV;
% 
% for u=1:g2
%     
%     
%  [ SMYAtPSV_azifilt_depth(:,u), zaxis ] =depthconverter3(dum1(1000:4000,u),SMYAheader(4,u)/111.12,100,elevation);
%  
% end
% 
% 
% 
% normd5=normalization3(SMYAtPSV_azifilt_depth);
% d6=sum(normd5,2);
% SMYAtPSV_azifilt_depth_stack=normalization3(d6);
%     root stack
% nthrootstackSMYAPSV2_azifilt_depth=rootstack(normd5',2)'
%   
% 
% figure(221)