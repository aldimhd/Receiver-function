function  adf = defective(i,j,n,a,b,c,d,e,px)
global propag stiff defect defect1 pstff
% call defective(i,j,n,adf(1,n),a,b,c,d,e,px)
% call defective(i,j,n,adf(2,n),a,b,c,d,e,px)
%subroutine defective(i,j,n,adf,a,b,c,d,e,px)
%  kluge for dealing with nearly defective propagator matrices
%  in which the eigenvectors,
%  which represent the particle motion of upgoing and downgoing waves
%  become nearly parallel.
%  in this case the solution for system of ODEs is
%  a_1 \bf_1 e^xnu*(z-z0) + a_2*(\bf_2 + adf*(z-z0)*\bf_1)e^xnu*(z-z0)
%
%   implicit real*8 (a-h,o-z)
%   implicit integer*4 (i-n)
%   complex*16 pp,u0,ee,z1,z0,znu,xnu,e1,e2,zla,xl,u,pfac,Eye
%   complex*16 zq1,zq2,u1,u2,zq3,xee
%   common/stfff/w(3,101),t(3,3),ttl(3,3),s(3,3),stl(3,3),
%  x                     			 r(3,3),x(3),y(3)
%   common/propag/xnu(6,101),xl(6,100),pfac(6,3),u(3,6)
%   common/defect1/zq1(3,3),zq2(3,3),u1(3),u2(3),zq3(2,2),xee(3)
%   common/defect2/edr(6),edi(6),qdr(5,5),qdi(5,5),ydr(6),ydi(6)
%   common/defect3/q1r(3,3),q1i(3,3),q2r(3,3),q2i(3,3),fv2(3),fv3(3)
%   common/mstff/qq(6,6),wr(6),wi(6),zr(6,6),zi(6,6),iv(6),fv(6)
%   common/pstff/pp(3),u0(3),ee(6,6,101),e1(6,6),e2(6,6),zla(6)
%   common/qstff/qi(6,6),xr(6),xi(6),yr(6),yi(6),ips(3)
z1=complex(1,0);
z0=complex(0,0);
Eye=complex(0,1);
zq1 = zeros(3,3);
zq2 = zeros(3,3);

q1r = zeros(3,3);
q1i = zeros(3,3);

q1= zeros(3,3);
q2 = zeros(3,3);


q2r = zeros(3,3);
q2i = zeros(3,3);

xnu = propag.xnu;
w = stiff.w;
ee = zeros(6,6,101);
% TODO xnu , global


%  for the extravector, need to solve system of equations
%  based on original 6x6 Q matrix
%  the plane-wave solutions generalize to the form
%  u0*e^{i*nu*(z-z0)}  and  u1*e^{i*nu*(z-z0)} + adf* u0*(z-z0)*e^{i*nu*(z-z0)}
%  u1 is the solution to
%  (\bTtil + nu*\bStil + nu^2*\bI).u1=i*adf*(\bStil + 2*nu*\bI).u0
%  in practice, we absorb the adf factor into u1, then normalize
%  (\bTtil + nu*\bStil + nu^2*\bI).(u1/adf)=i*(\bStil + 2*nu*\bI).u0
%  since nu is the known eigenvalue of u0, the solution is easier
%  form the matrices on either side
znu=-Eye*xnu(i,n);
for ii=1:3
    for jj=1:3
        zq1(jj,ii)=complex(ttl(jj,ii),0)+znu*stl(jj,ii);
        zq2(jj,ii)=complex(stl(jj,ii),0);
    end
    zq1(ii,ii)=zq1(ii,ii)+znu*znu;
    zq2(ii,ii)=zq2(ii,ii)+2*znu;
end
%  we wish to find the eigenvector of the near-defective matrix
%  in the region where its eigenvectors are numerically unstable
%   we explicitly calculate the eigenvector with smallest right-eigenvalue of
%  (\bTtil + nu*\bStil + nu^2*\bI)=zq1
%  copy into real, imag matrices
for ii=1:3
    for jj=1:3
        q1r(jj,ii)=real(zq1(jj,ii));
        q1i(jj,ii)=imag(zq1(jj,ii));
    end
end


[wr, wi, q2r, q2i, ierr] = cbal(q1r, q1i);
q1r
% %  into eispack
% call cbal(3,3,q1r,q1i,low,igh,fv)
% call corth(3,3,low,igh,q1r,q1i,fv2,fv3)
% call comqr2(3,3,low,igh,fv2,fv3,q1r,q1i,wr,wi,q2r,q2i,ierr)
% if ierr ~= 0
%     disp('Error');
%     return;
% end
% call cbabk2(3,3,low,igh,fv,3,q2r,q2i)

[~ , ij] = min(wr.^2+wi.^2);
% amn=wr(1)^2+wi(1)^2;
% ij=1;
% for ii=2:3
%     amm=wr(ii)^2+wi(ii)^2;
%     if amm < amn
%         ij=ii;
%         amn=amm;
%     end
% end
% q2r = real(q2);
% q2i = imag(q2);
u0 = complex(q2r(:,ij),q2i(:,ij));
u0 = u0/norm(u0);
% sum=0
% for ii=1:3
%     u0(ii)=complex(q2r(ii,ij),q2i(ii,ij));
%     sum=sum+abs(u0(ii))^2;
% end
% sum=sqrt(sum);
% for ii=1:3
%     u0(ii)=u0(ii)/sum;
% end
%  assemble the ith stress-displacement vector
%  calculate the traction components, with i removed
pp(1)=complex(px,0);
pp(2)=z0;
pp(3)=znu;
abcde=a-b+c-2*d+2*e;
bce=b-4*c-4*e;
de=d-e;
pu = sum(pp.*u0);
pw = sum(pp.*w(:,n));
uw = sum(u0.*w(:,n));
% for ii=1:3
%     pu=pu+pp(ii)*u0(ii);
%     pw=pw+pp(ii)*w(ii,n);
%     uw=uw+u0(ii)*w(ii,n);
% end
for ii=1:3
    ee(ii,i,n)= u0(ii);
    ee(ii+3,i,n)=w(ii,n)*(pu*w(3,n)*bce+8*pw*uw*w(3,n)*c ...
        +2*(pw*u0(3)+uw*pp(3))*e);
    ee(ii+3,i,n)=ee(ii+3,i,n)+pp(ii)*(u0(3)*de+2*uw*w(3,n)*e);
    ee(ii+3,i,n)=ee(ii+3,i,n)+u0(ii)*(pp(3)*de+2*pw*w(3,n)*e);
end
ee(6,i,n)=ee(6,i,n)+pu*abcde+pw*uw*bce;
%  almost lastly, mult traction by i
for ii=1:3
    ee(ii+3,i,n)=Eye*ee(ii+3,i,n);
end
%  extract u0 from ee(*,i,n) use it to calculate the additional traction terms
%  and store in ee(*,j,n)
%  additional traction terms involve gradient of (z-z0)
%  so can be calculated from standard formulas with \bk=zhat
%  we dont multiply by i
pp(1)=z0;
pp(2)=z0;
pp(3)=z1;
abcde=a-b+c-2*d+2*e;
bce=b-4*c-4*e;
de=d-e;

% for ii=1:3
%     u0(ii)=ee(ii,i,n);
%     pu=pu+pp(ii)*u0(ii);
%     pw=pw+pp(ii)*w(ii,n);
%     uw=uw+u0(ii)*w(ii,n);
% end

pu = sum(pp(:).*ee(:,i,n));
pw = sum(pp(:).*w(:,n));
uw = sum(ee(:,i,n).*w(:,n));

xee = zeros(3,1);
for ii=1:3
    xee(ii)=w(ii,n)*(pu*w(3,n)*bce+8*pw*uw*w(3,n)*c ...
        +2*(pw*u0(3)+uw*pp(3))*e);
    xee(ii)=xee(ii)+pp(ii)*(u0(3)*de+2*uw*w(3,n)*e);
    xee(ii)=xee(ii)+u0(ii)*(pp(3)*de+2*pw*w(3,n)*e);
end
xee(3)=xee(3)+pu*abcde+pw*uw*bce;
%  extract u0 from ee(*,i,n), mult by i*(\bStil + 2*nu*\bI), replace in u0
for ii=1:3
    u0(ii) = Eye *sum(zq2(ii,:).*ee(:,i,n));
    % for jj=1:3
    %     u0(ii)=u0(ii)+zq2(ii,jj)*ee(jj,i,n)
    % end
    % u0(ii)=Eye*u0(ii);
end
% 1002 format(3(2g14.6,3x))
%  for znu NOT an eigenvalue,
%  but rather the average of closely-space eigenvalues
%  in this case, zq1 is nonsingular, and we just solve for u1
% for ii=1:3
%     yr(ii)=real(u0(ii));
%     yi(ii)=imag(u0(ii));
%     for jj=1:3
%         q1r(jj,ii)=real(zq1(jj,ii));
%         q1i(jj,ii)=imag(zq1(jj,ii));
%     end
% end
u1 = zq1\u0;
% call csolve(3,q1r,q1i,xr,xi,yr,yi)
% for ii=1:3
%     u1(ii)=complex(xr(ii),xi(ii));
% end
%   End, different tactic

%  normalize
u1 = u1/norm(u1);
%  adf is the normalization constant
adf=1/norm(u1);
%  calculate the traction
%  and place the new stress-displacement vector in column j
%  pp is the wavenumber vector, and first two components are already in place
pp(1)=complex(px,0);
pp(2)=z0;
pp(3)=znu;
abcde=a-b+c-2*d+2*e;
bce=b-4*c-4*e;
de=d-e;
pu = sum(pp(:).*u1(:)); 
pw = sum(pp(:).*w(:,n));
uw = sum(u1(:).*w(:,n));
% for ii=1:3
%     pu=pu+pp(ii)*u1(ii);
%     pw=pw+pp(ii)*w(ii,n);
%     uw=uw+u1(ii)*w(ii,n);
% end
for ii=1:3
    ee(ii,j,n)=u1(ii);
    ee(ii+3,j,n)=w(ii,n)*(pu*w(3,n)*bce+8*pw*uw*w(3,n)*c ...
        +2*(pw*u1(3)+uw*pp(3))*e);
    ee(ii+3,j,n)=ee(ii+3,j,n)+pp(ii)*(u1(3)*de+2*uw*w(3,n)*e);
    ee(ii+3,j,n)=ee(ii+3,j,n)+u1(ii)*(pp(3)*de+2*pw*w(3,n)*e);
end
ee(6,j,n)=ee(6,j,n)+pu*abcde+pw*uw*bce;
%  almost lastly, mult traction by i
%  and add extra traction from (z-z0) term (not mult by i)
%  TEST - mult xee by zero, see if it is important --- it IS important
for ii=1:3
    ee(ii+3,j,n)=Eye*ee(ii+3,j,n)+adf*xee(ii);
end
stiff.w = w;

defect1.xee = xee;
pstff.ee = ee;
end
