function [ data2, zaxis ] =depthconverter3(thissta,rayp,fs1,stel)
%
%function to convert time to depth
% input mexico depth profile
invel='merapi_velfile2.dat';
zaxis=0:1:300;
velin=load(invel);
     a=6371;
     d2km=(2*pi*a/360);
     dz=min(diff(zaxis))/3;
     zv=0:dz:(max(zaxis)+10);
     vp = interp1(velin(:,1),velin(:,2),zv);
     vs = interp1(velin(:,1),velin(:,3),zv);
     zf = -a .* log((a-zv)./a);
     Vpf = (a./(a-zv)).*vp;
     Vsf = (a./(a-zv)).*vs;
    lz = length(zf);
    lz2=length(zaxis);
    h = diff(zf);
    
    AngS = asin(min(rayp.*Vsf(1:lz-1),1));
    realang=[ 1, 1-floor(AngS/(pi/2))];
    xS = real(tan(AngS).*h);
    Vsf2 = max(Vsf(1:lz-1).^-2,rayp^2);
    Vpf2 = max(Vpf(1:lz-1).^-2,rayp^2);
    ts = h.*(Vsf2 - rayp^2).^.5;
    tp = h.*(Vpf2 - rayp^2).^.5;
    Tsp = [ 0,  cumsum(ts)-cumsum(tp) ];
    XS = [ 0,  cumsum(xS)];
    data = thissta(min(ceil(Tsp*fs1 +1),length(thissta)));
    XS2=interp1([-10, zv-stel/1000],[0, XS],zaxis);
    data2=interp1([-10, zv-stel/1000],[0; data],zaxis);
    data2=detrend(data2);
    