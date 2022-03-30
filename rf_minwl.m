function [ recfun, taxis, rfspec, faxis ] = rf_min(z,r,wl,a,fs,delay)
% rf_min
% compute receiver function by building tradeoff curve between
% spectral division data fit and solution roughness defined as
% the norm of the first derivative.
%
% usage:
%      [ recfun, taxis ] = rf_min(z,r,fs)
%
%  z,r are vertical and radial-horizontal component data
% (or radial and tangential data), fs is sampling frequency

[s1, s2] = size(z);
if s1>s2, z=z'; end
[s1, s2] = size(r);
if s1>s2, r=r'; end

datz  = [ z zeros(1,length(z))];
datr = [ r zeros(1,length(r))];

n = length(datz);
df = fs/n;
w = [(0:df:fs/2) -((fs/2-df):-df:df)].*(2*pi);
G =  exp(-1*w.^2 ./ (4*a^2));
fdatz = fft(datz);
fdatr = fft(datr);

magdatz = fdatz .* conj(fdatz);
fampmax = max(magdatz);
phi = sqrt(max(magdatz,(wl * fampmax)));
angles = angle(fdatz);
     FdatZwl = phi.*cos(angles) + i*phi.*sin(angles);

% spectral division
clear Frf Frfmat
T1=now;

p=length(fdatz);
clear AA ii jj kk
    ii=[];jj=[];kk=[];
    realvec=real(FdatZwl);
    imagvec=imag(FdatZwl);
    ii = [ ((1:p)*2-1) (1:p)*2];
    jj = [ ((1:p)*2-1) (1:p)*2];
    kk = [ realvec realvec ];
    ii = [ ii (1:p)*2 ((1:p)*2-1) ];
    jj = [ jj ((1:p)*2-1) (1:p)*2 ];
    kk = [ kk imagvec -imagvec];
    AA = sparse(ii,jj,kk);
    dd = zeros(p*2,1);
    dd(1:2:2*p-1) = real(fdatr.');
    dd(2:2:2*p) = imag(fdatr.');

    ii=[];jj=[];kk=[];
    ii = [  ((1:p)*2-1) (1:p)*2 ];
    jj = [  ((1:p)*2-1) (1:p)*2  ];
    kk = [ w w];
    BB = sparse(ii,jj,kk);
    dd0 = zeros(2*p,1);
maxval=max([realvec imagvec]);
alphacount=1;
clear model  m
for expon=-4:.25:1,
  alpha2 = maxval.*10^expon;
  BB2=BB .* alpha2;

  m=cgls3([ AA; BB2],[ dd; dd0],100) ;
  model(:,alphacount) = m;
  misfit(alphacount)=norm((AA*m-dd));
  roughness(alphacount)=norm(((BB*m)));
  alphaval(alphacount)=alpha2;
  alphacount = alphacount+1;

end
mexact = (fdatr./fdatz).';
 m(1:2:2*p-1,1)=real(mexact);
m(2:2:2*p,1)=imag(mexact);
 maxmis = norm((dd));
 maxrough= norm(((BB*m)));

[YY2,II2]=min((misfit/max(maxmis)).^2+(roughness/max(maxrough)).^2);


Frf=model(:,II2);
Frf2 = Frf(1:2:2*p-1,1).' + sqrt(-1).*Frf(2:2:2*p,1).';

recfun=real(ifft(Frf2.*G.*exp(sqrt(-1).*w.*-delay)));
taxis=(1:length(recfun))/fs -delay;

recfun=recfun(1:length(z)).';
taxis=taxis(1:length(z));

return;

