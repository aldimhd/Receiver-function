close all
clear all
%
%
global model model2 stiff
% ======== lines 48-69
%
%  we reduce the condition numbers of matrices by
%  normalizing physical quantities to make them dimensionless
%  Ill use the normal-mode normalizations
%  which are a little peculiar for the crust, but what the hey!
rbar=5.515e3;
ren=1.075190645e-3;
radi=6.371e6;
vbar=ren*radi;
con=rbar*vbar^2;
z1=complex(1,0);
z0=complex(0,0);
%     read the frequency you want to work with,
%     and decide how many frequency points you need.
nfrq=512;
gainfac = 16;
%  gainfac above will compensate timeseries values for
% the effects of spectral smoothing.
%  it includes a factor of 2 to remove the effect of the cos**2 taper
% and a halved ratio of freq. series length to the number of frequencies
% for nfrq=512 and npad = 8192 gainfac= 8
% ======== lines 70-93
reply=[];
% while isempty(reply)
%     reply = input('enter maximum frequency:','s');
% end
% frqmax = str2double(reply);
frqmax=1.5;

if frqmax > 1
    nfrq=1024;
    gainfac=8;
end
%	read backazimuth (this only makes sence for anisotropic models)
reply=[];
% while isempty(reply)
%     reply = input('backazimuth:','s');
% end
% baz = str2double(reply);
baz=165;
% input file
reply=[];
% while isempty(reply)
%     reply = input('Enter model file name: ','s');
% end
% filein = reply;
filein='tomovelME32.txt'

fh = fopen(filein,'r');
if fh == -1
    disp('Can not open the input file')
    return;
end
title = fgetl(fh);
% title = fscanf(fh,'%s');
disp(title);
nl = fscanf(fh,'%d',1);
%  read in theta,phi in degrees - polar coords of fast axis
nlp = nl+1;
nlm = nl-1;
w = zeros(3,nlp);
z = zeros(nlp,1);
vp = z;
vp2 = vp;
vp4 = vp;
vs = z;
vs2 = vs;
rho = z;

%
%  94-127
for i=1:nlp
    % angle conventions for anisotropic tensor orientations
    % numbers supplied in the input file are :
    % phig - azimuth of symmetry axis, theta - tilt from vertical
    % internal computations use phi measured
    % counterclockwise from X (=radial=baz-180)
    data = fscanf(fh,'%f %f',[2,1]);
    theta = data(1);
    phig = data(2);
    % fscanf(fh,'%f %f',theta,phig);
    %-----correct for backazimuth
    %        phi = baz-180-phig
    % change to:		12/13/96 BUG DISCOVERED THE NIGHT BEFORE AGU...
    phi = phig-baz;
    %  end Change
    w(1,i)=sind(theta)*cosd(phi);
    w(2,i)=sind(theta)*sind(phi);
    w(3,i)=cosd(theta);
    fprintf('%f  %f  %f\n',w(1,i), w(2,i), w(3,i));
    %  read depth to ith interface, vp (m/sec), pk-to-pk cos(2th) relative P pert
    %  pk-to-pk cos(4th) relative P pert, v_s,  pk-to-pk cos(2th) relative S pert
    %  density (kg/m**3)
    data = fscanf(fh,'%f %f %f %f %f %f %f',[7,1]);
    z(i) = data(1);
    vp(i) = data(2);
    vp2(i) = data(3);
    vp4(i) = data(4);
    vs(i) = data(5);
    vs2(i) = data(6);
    rho(i) = data(7);
    % fscanf(fh,'%f %f %f %f %f %f %f',z(i),vp(i),vp2(i),vp4(i),vs(i),vs2(i),rho(i));
    %  recall that we interpret fractional values of b,c,e
    %  as peak-to-peak relative velocity perts.
    %  therefore, e=0.02 is 2% pert to mu from slowest to fastest
end
fclose(fh);
xmu=rho.*vs.^2/con;
xmu2=vs2.*xmu;
xla=rho.*vp.^2/con;
xla2=vp2.*xla;
xla4=vp4.*xla;
vs=vs./vbar;
vp=vp./vbar;
rho=rho/rbar;
z=z/radi;

% line 130
dz = diff(z(1:end-1));
dz = [z(1) ; dz];
model.dz = dz;
model.z = z;

% line 134-166
%  print the organ-pipe mode count for 1Hz
%  the lowest layer (nl+1) is taken as evanescent region.
sn = sum(dz./vs(1:nl))/ren;
pn = sum(dz./vp(1:nl))/ren;

disp('organ-pipe mode count at 1 Hz in stack of layers: S & P')
fprintf('%10.1f  %10.1f\n', sn, pn);
%  search for cmin, cmax is halfspace velocity
cmin = max(vs);
csmin=vs(end)*vbar/1000.;
cpmin=vp(end)*vbar/1000.;
fprintf('minimum phase velocity for S wave (km/sec) %f\n', csmin)
fprintf('minimum phase velocity for P wave (km/sec) %f\n', cpmin)
reply=[];
% while isempty(reply)
%     reply = input('enter phase velocity of incident wave (km/sec): ','s');
% end
% % line 157
% cc = str2double(reply);

cc=20;
if cc <= 0
    return;
end
cc=cc*1000/vbar;
%  calc the eigenvector decomps for the layers
%  need to identify upgoing P,SV,SH in the halfspace
% TODO matget
model.rho = rho;
model2.xla = xla;


model2.xla2 = xla2;
model2.xla4 = xla4;
model2.xmu = xmu;
model2.xmu2 = xmu2;

stiff.w = w;

% TODO
matget(nl,cc)

%  loop over frequency, calc reflection/transmission matrices
%  calc 3-comp transfer fct response at surface
df=frqmax/nfrq;

frqq = df*(1:nfrq);
resp = zeros(3,3,nfrq);
for jf=1:nfrq
    om=2*pi*jf*df/ren;
    temp = respget(nl,om,cc);
    resp(:,:,jf) =  temp;
end

npad=8192;
dt=1/(npad*df);

% line 177
zz4(1:npad) = 0;

%  zero the DC and Nyquist
%  ick switches the sign of y and z components to transverse & vertical
cc4(1:2)=0;
cc5(1:2)=0;

%      open(10,file = "RRF.aspec")
%      open(11,file = "fac")

for jf=1:nfrq
    %    figure out the value of the cosine-squared scaler
    fac=cos(0.5*pi*(jf-1)/nfrq)^2;
    
    zz= -resp(1,1,jf) / resp(3,1,jf);
    cc4(2*jf+1)=real(zz)*fac;
    cc4(2*jf+2)=imag(zz)*fac;
    %
    %   ampsp =  cc4(2*jf+1)**2 + cc4(2*jf+2)**2
    %   write (10,332) jf*df, ampsp
    %   write (11,332) jf*df, fac
    %   cc4 has spectra of P-SV reciever function
    zz= resp(2,1,jf) / resp(3,1,jf);
    cc5(2*jf+1)=real(zz)*fac;
    cc5(2*jf+2)=imag(zz)*fac;
    
    %   ampsp = cc5(2*jf+1)**2 + cc5(2*jf+2)**2
    %   write (11,332) jf*df, ampsp
    
    %   cc5 has a spectra of P-SH reciever function
end


% line 228
cc4(2*nfrq+3 : npad)=0;
cc5(2*nfrq+3 : npad)=0;

cc4r = cc4(1:2:end);
cc4i = cc4(2:2:end);
cc4 = refft(cc4r,cc4i);
% cc4p= complex(cc4r,cc4i);

cc5r = cc5(1:2:end);
cc5i = cc5(2:2:end);
cc5 = refft(cc5r,cc5i);
% cc5p= complex(cc5r,cc5i);
% TODO refft
%%% call refft(cc4,npad,-1,-1)
%%% call refft(cc5,npad,-1,-1)
%   cc4 and cc5 have time domain reciever functions
%        now mult by 2 since cosine**2 taper is applied to spectral RF
%        also mult by gainfac to compensate for losing high freqs
% cc4 = fft(cc4);
% cc5 = fft(cc5);
cc4(1: npad)=cc4(1 : npad)*gainfac;
cc5(1: npad)=cc5(1 : npad)*gainfac;

% line 245
amx = max(cc4);
amx1 = max(cc5);
%   define new trace length - don't need a long tail of zeroes and junk
%   output only first 100 seconds
npad = 1600;
%   output RFs in ascii form, appropriate for use, e.g., with
%   GMT command pswiggle
% TODO write output to files RRF.table & TRF.table
fhr = fopen('RRF.table', 'w');
fht = fopen('TRF.table', 'w');
for i = 1:npad
    timestep = (i-1)*dt;    
    fprintf(fhr, '%10.4f  %8.3f %11.8f \n', timestep, baz, cc4(i));
    fprintf(fht, '%10.4f  %8.3f %11.8f \n', timestep, baz, cc5(i));
end
fclose all;
t=(0:npad-1)*dt;
subplot(2,1,1)
plot(t,cc4(1:npad));
subplot(2,1,2)
plot(t,cc5(1:npad));

B=[t(1:50)' cc4(1:50)]

%if time axis is not right
%normalisasi data
vv=max(abs(cc4));
vvv=find(cc4==vv);

v1=length(cc4)
for i=1:v1
    cc4norm(i)=cc4(i)/vv;
end


tt=t*2

cc=find(tt>10)
figure(2)
subplot(2,1,1)
h=500
plot(tt(1:cc(1)),cc4(1:cc(1)));
subplot(2,1,2)
plot(tt(1:cc(1)),cc5(1:cc(1)));

figure(99999)
plot(tt(1:cc(1)),cc4norm(1:cc(1)));
% open(10,file = "RRF.table")
% open(11,file = "TRF.table")
% for i=1:npad

%     timestep = (i-1)*dt;
%     write (10,333) timestep, baz, cc4(i)
%     write (11,333) timestep, baz, cc5(i)
% end
% f10.4,1x,f8.3,1x,f11.8)

