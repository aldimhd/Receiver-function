function resp = respget(nl,om,cc)
    global model propag stiff defect defect1 pstff qstff
%  returns surface response for a stack of anisotropic layers
%  incident p,sv,sh waves with freq om and phase velocity cc
%  iev=1 if the waves are evanescent in the top layer, iev=0 otherwise
% implicit real*8 (a-h,o-z)
% implicit integer*4 (i-n)
% complex*16 pp,u0,ee,z1,z0,xnu,Eye,e1,e2,zla,rtm
% complex*16 rt,tt,rt0,trc,xl,pfac,u,co,resp,ur
% common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
% x                     			 r(3,3),x(3),y(3)
% common/model/z(100),dz(100),rho(101),vp(101),vs(101),vp2(101),
% x               vp4(101),vs2(101),vss(101)
% common/model2/xmu(101),xla(101),xmu2(101),xla2(101),xla4(101)
% common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
% common/defect/idfct(4,101),adf(2,101)
% common/rrt/rtm(6,6,100)
% common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
% common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
% common/pstf2/co(6,101),ur(3)
% common/qstff/qi(6,6),xr(6),xi(6),yr(6),yi(6),ips(3)
% common/rstff/rt(3,3,101),tt(3,3,101),rt0(3,3),trc(3,3)
% dimension resp(3,3)
% data pi/3.14159265358979d0/,eps/1.d-6/,tol/1.d-7/

%  set iev=1
%  toggle to iev=0 if there is a purely propagating wave in the top layer n=1
iev=1;
z1=complex(1,0);
z0=complex(0,0);
Eye=complex(0,1);
rbar=5.515e3;
ren=1.075190645e-3;
radi=6.371e6;
vbar=radi*ren;
con=rbar*radi*radi*ren*ren;

resp = zeros(3,3);
zla = zeros(6,1);

nlp=nl+1;
nlm=nl-1;
e1 = pstff.e1;
xnu = propag.xnu;
dz = model.dz;
z = model.z;
ee = pstff.ee;
idfct = defect.idfct;
ips = qstff.ips;
%  first calculate vertical wavenumbers and propagating waves for each layer
%    an eigenvector problem was solved in prior subroutine
%   results stored in array ee(6,6,101)
%  in general, the evanescent vertical wavenumbers have nonzero real parts
%   complex exponential fct is used to avoid endless branching
%
%  calculate modified R/T coefficients
%  first calc the propagation factors
%  note that for dipping fast axes the upgoing and downgoing wavenumbers are
%  independent, so we must calc all to be safe
xl = zeros(6,nl);
rtm = zeros(6,6,nl);
rt0 = zeros(3,3);
trc = zeros(3,3);
for n=1:nl
    for k=1:3
        xl(k,n)=exp(om*xnu(k,n)*dz(n));		% downgoing
        xl(k+3,n)=exp(-om*xnu(k+3,n)*dz(n))	;% upgoing
    end
end
%      for i=1:6
%        print 1002,xnu(i,3),xl(i,3)
%      end
%  1002 format('i*k_z:',2g15.6,',  propfac is',2g15.6)
%  calculate modified R/T coefficients at each interface
e1 = zeros(6,6);
e2 = zeros(6,6);
for n=1:nl
    %  rearrange to e1: waves approaching and e2: waves leaving an interface
    for k=1:3
        for i=1:6
            e1(i,k)=ee(i,k,n+1);
            e2(i,k)=ee(i,k,n);
            e1(i,k+3)=-ee(i,k+3,n);
            e2(i,k+3)=-ee(i,k+3,n+1);
        end
        zla(k)=xl(k,n);
        if(n < nl)
            zla(k+3)=xl(k+3,n+1);
        else
            % c  reference the upcoming wave amplitude to the top of halfspace
            % c  therefore no propagation factor, not valid for evanescent waves in halfspace
            % c  in surface wave code this is zero, so that upgoing evanescent waves vanish
            zla(k+3)=1;
        end
    end
    % mult the columns of e2
    for k=1:6
        for i=1:6
            e2(i,k)=e2(i,k)*zla(k);
        end
    end
    %  the possibility of defective matrices must be contemplated here
    %  k=1,2,3 columns are downgoing in nth layer
    %  k=4,5,6 columns are upgoing in (n+1)th layer
    %  the vector e2(.,k1) has already been multiplied by exponential factor zla
    if(idfct(1,n) ~= 0)
        k1=idfct(1,n);
        k2=idfct(2,n);
        for i=1:6
            e2(i,k2)=e2(i,k2)+adf(1,n)*dz(n)*e2(i,k1);
        end
    end
    %  the sign change on dz is for upgoing waves
    if(idfct(3,n+1) ~= 0)
        k1=idfct(3,n+1);
        k2=idfct(4,n+1);
        for i=1:6
            e2(i,k2)=e2(i,k2)-adf(2,n+1)*dz(n+1)*e2(i,k1);
        end
    end
    %  in order to use csolve to invert e1, must separate into real/imag parts
    %  its clumsy, but im lazy
    %  we calc e1^{-1}\cdot e2\cdot \Gamma one column at a time
    % e1\e2 e1\Gamma
    %for k=1:6
        rtm(:,:,n) = e1\e2;
        % for i=1:6
        %     qq(i,k)=real(e1(i,k));
        %     qi(i,k)=imag(e1(i,k));
        % end
    %end
    % nn=6;
    % for k=1:6
    %     for i=1:6
    %         yr(i)=real(e2(i,k));
    %         yi(i)=imag(e2(i,k));
    %     end
    %     call csolve(nn,qq,qi,xr,xi,yr,yi)
    %     nn=-6
    %     for i=1:6
    %         rtm(i,k,n)=complex(xr(i),xi(i))
    %     end
    % end
end
s = zeros(3,3);
t = zeros(3,3);
%  calc R_ud at the free surface
%  note that first two factors in Chen (20) dont collapse
%  mult by inv-matrix one column at a time
for k=1:3
    for i=1:3
        rt0(i,k)=ee(i+3,k+3,1)*xl(k+3,1);
        s(i,k)=real(ee(i+3,k,1));
        t(i,k)=imag(ee(i+3,k,1));
    end
end
st = complex(s,t);

%  the possibility of defective matrices must be contemplated here
%  these waves are upgoing in 1st layer
%  the sign change on dz is for upgoing waves, and xl(k1,1)=xl(k2,1)
if(idfct(3,1) ~= 0)
    k1=idfct(3,1);
    k2=idfct(4,1)-3;
    for i=1:3
        rt0(i,k2)=rt0(i,k2)-adf(2,1)*dz(1)*ee(i+3,k1,1)*xl(k1,1);
    end
end
rt0 = -st\rt0;

% nn=3;
% for k=1:3
%     for i=1:3
%         yr(i)=real(rt0(i,k));
%         yi(i)=imag(rt0(i,k));
%     end
%     call csolve(nn,s,t,xr,xi,yr,yi)
%     nn=-3
%     for i=1:3
%         rt0(i,k)=-complex(xr(i),xi(i));
%     end
% end
%  recursive calc of generalized R/T coefs:
%  in contrast to the surface-wave code, we start from the top layer and
%  iterate down to the halfspace
%  we also uses submatrices of generalized R/T matrix in different order
rt = zeros(3,3,nl);
tt = zeros(3,3,nl);

for n=1:nl
    %  first the generalized upward-transmission coef:
    for k=1:3
        for i=1:3
            trc(i,k)=z0;
            if(n > 1)
                for j=1:3
                    trc(i,k)=trc(i,k)-rtm(i+3,j,n)*rt(j,k,n-1);
                end
            else
                %c  use free-surface reflection matrix in top layer (interface "zero")
                for j=1:3
                    trc(i,k)=trc(i,k)-rtm(i+3,j,n)*rt0(j,k);
                end
            end
        end
        trc(k,k)=trc(k,k)+z1;
    end
    for k=1:3
        for i=1:3
            s(i,k)=real(trc(i,k));
            t(i,k)=imag(trc(i,k));
        end
    end

%     for k=1:3
%         for i=1:3
%             yr(i)=real(rtm(i+3,k+3,n));
%             yi(i)=imag(rtm(i+3,k+3,n));
%         end
%         
%         tt(:,k,n) = trc\complex(yr,yi);
%         % call csolve(nn,s,t,xr,xi,yr,yi)
%         % nn=-3
%         % for i=1:3
%         %     tt(i,k,n)=complex(xr(i),xi(i));
%         % end
%     end
    tt(:,:,n) = trc\rtm(4:6,4:6,n);
    %tt(:,:,n)
    %tt(:,:,n)
    %c  next the generalized reflection coef:
    for k=1:3
        for i=1:3
            trc(i,k)=z0;
            if(n > 1)
                for j=1:3
                    trc(i,k)=trc(i,k)+rt(i,j,n-1)*tt(j,k,n);
                end
            else
                % c  use free-surface reflection matrix in top layer (interface "zero")
                for j=1:3
                    trc(i,k)=trc(i,k)+rt0(i,j)*tt(j,k,n);
                end
            end
        end
    end
    for k=1:3
        for i=1:3
            rt(i,k,n)=rtm(i,k+3,n);
            for j=1:3
                rt(i,k,n)=rt(i,k,n)+rtm(i,j,n)*trc(j,k);
            end
        end
    end
end
% 1001 format(6f14.6)
% c      print *,'free-surface reflection'
% c      print 1001,((rt0(i,j),j=1,3),i=1,3)
% c      do n=1,nl
% c        print *,'interface',n
% c        print 1001,((rt(i,j,n),j=1,3),i=1,3)
% c        print 1001,((tt(i,j,n),j=1,3),i=1,3)
% c      end
% c  using the p,sv,sh identification, we propagate upward to the surface,
% c  calculate the wave coefs in the top layer, then the particle displacement
co = zeros(6,nlp);
%rt
for iup=1:3
    for i=1:3
        co(i+3,nlp)=z0;
    end
    co(iup+3,nlp)=z1;
    % c  from upgoing coefs in the n+1 layer, calculate
    % c  upgoing coefs in the nth layer, downgoing coefs in the n+1 layer
    for n=nl:-1:1
        for i=1:3
            co(i+3,n)=z0;
            co(i,n+1)=z0;
            for j=1:3
                co(i+3,n)=co(i+3,n)+tt(i,j,n)*co(j+3,n+1);
                co(i,n+1)=co(i,n+1)+rt(i,j,n)*co(j+3,n+1);
            end
        end
    end
    % c  then downgoing coefs in the top layer:
    for i=1:3
        co(i,1)=z0;
        for j=1:3
            co(i,1)=co(i,1)+rt0(i,j)*co(j+3,1);
        end
    end
    % c        print *,'upgoing coefs'
    % c        print 1001,((co(j+3,n),j=1,3),n=1,nlp)
    % c        print *,'downgoing coefs'
    % c        print 1001,((co(j,n),j=1,3),n=1,nlp)
    % c  calc the surface displacement
    h1=0;
    h2=z(1);
    ur = zeros(3,1);
    
    for i=1:3
        % ur(i)=z0;
        for k=1:3
            ur(i)=ur(i)+co(k,1)*ee(i,k,1) ...
                +co(k+3,1)*ee(i,k+3,1)*(exp(om*xnu(k+3,1)*(-h2)));
        end
        % c  check for the xtra terms associated with defective matrices
        if(idfct(1,1) ~= 0)
            ii=idfct(1,1);
            jj=idfct(2,1);
            ur(i)=ur(i)+co(jj,1)*adf(1,1)*ee(i,ii,1)*(-h1);
        end
        if(idfct(3,1) ~= 0)
            ii=idfct(3,1);
            jj=idfct(4,1);
            ur(i)=ur(i) +co(jj,1)*adf(2,1)*ee(i,ii,1)*(-h2)*(exp(om*xnu(ii,1)*(-h2)));
        end
    end
    
    %  copy the surface displacement into the response matrix
    for i=1:3
        resp(i,ips(iup))=ur(i);
    end
end
% disp('response');
% disp(resp)
pstff.e1 = e1;
pstff.ee = ee;
propag.xnu =xnu;
end