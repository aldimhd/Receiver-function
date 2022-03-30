function matget(nl,cc)
global model model2 stiff propag pstff qstff defect

tol = 1e-7;


xla = model2.xla;
xla2 = model2.xla2;
xla4 = model2.xla4;
xmu = model2.xmu;
xmu2 = model2.xmu2;

w = stiff.w;
rho = model.rho;

rbar=5.515e3;
ren=1.075190645e-3;
radi=6.371e6;
vbar=ren*radi;
con = rbar*radi*radi*ren*ren;

nlp=nl+1;
nlm=nl-1;

ee = zeros(6,6,nlp);
xnu = zeros(6,nlp);
idfct = zeros(4,nlp);
adf = zeros(2,nlp);

z1 = complex(1,0);
z0 = complex(0,0);
Eye = complex(0,1);

r=zeros(3,3);
s=zeros(3,3);
ttl=zeros(3,3);
stl=zeros(3,3);
t=zeros(3,3);
x=zeros(3,1);
y=zeros(3,1);

q=zeros(6,6);

ee = zeros(6,6,nlp);
pstff.ee = ee;

% line 312
%    first calculate vertical wavenumbers and propagating waves for each layer
%     requires an eigenvector problem be solved
%    in general, the evanescent vertical wavenumbers have nonzero real parts
%     complex exponential fct is used to avoid endless branching
%    horizontal slowness p_x
px=1/cc;
for n=1 : nlp
    a=xla(n);
    b=xla2(n);
    c=xla4(n);
    d=xmu(n);
    e=xmu2(n);
    fact=8*w(1,n)*w(1,n)*c+2*e;
    facs=16*w(1,n)*w(3,n)*c;
    facr=8*w(3,n)*w(3,n)*c+2*e;
    for i=1:3
        %    first the what-0-what tensor
        for j=1:3
            t(j,i)=fact*w(j,n)*w(i,n);
            s(j,i)=facs*w(j,n)*w(i,n);
            r(j,i)=facr*w(j,n)*w(i,n);
        end
        %    next the identity tensor - correct an error on 7/6/95
        t(i,i)=t(i,i)+d+e*(2*w(1,n)*w(1,n)-1);
        s(i,i)=s(i,i)+4*e*w(1,n)*w(3,n);
        r(i,i)=r(i,i)+d+e*(2*w(3,n)*w(3,n)-1);
    end
    
    fac=b-4*c-2*e;
    %    next the what-0-xhat and what-0-zhat tensors
    for i=1:3
        t(1,i)=t(1,i)+fac*w(1,n)*w(i,n);
        t(i,1)=t(i,1)+fac*w(1,n)*w(i,n);
        s(1,i)=s(1,i)+fac*w(3,n)*w(i,n);
        s(i,1)=s(i,1)+fac*w(3,n)*w(i,n);
        s(3,i)=s(3,i)+fac*w(1,n)*w(i,n);
        s(i,3)=s(i,3)+fac*w(1,n)*w(i,n);
        r(3,i)=r(3,i)+fac*w(3,n)*w(i,n);
        r(i,3)=r(i,3)+fac*w(3,n)*w(i,n);
    end
    fac=a-b+c-d+e;
    %    finally the xhat-0-xhat, zhat-0-zhat, xhat-0-zhat, zhat-0-xhat tensors
    t(1,1)=t(1,1)+fac;
    s(3,1)=s(3,1)+fac;
    s(1,3)=s(1,3)+fac;
    r(3,3)=r(3,3)+fac;
    %    mult by horizontal slowness and calc the modified T-matrix
    for i=1:3
        for j=1:3
            t(j,i)=t(j,i)*px*px;
            s(j,i)=s(j,i)*px;
        end
        t(i,i)=t(i,i)-rho(n);
    end
    %    calculate R ^(-1).S, R ^(-1).T, using routine solve
    stl = r\s;
    ttl = r\t;
    % nn=3;
    % for i=1:3
    %     for j=1:3
    %         y(j)=s(j,i);
    %     end
    %     %% TODO call solve
    %     call solve(nn,r,x,y)
    %     for j=1:3
    %         stl(j,i)=x(j);
    %     end
    %     nn=-3;
    % end
    % for i=1:3
    %     for j=1:3
    %         y(j)=t(j,i);
    %     end
    %     %% TODO cal solve
    %     call solve(nn,r,x,y)
    %     for j=1:3
    %         ttl(j,i)=x(j);
    %     end
    % end
    %    fill the 6x6 Q-matrix
    qq = zeros(6,6);
    for i=1:3
        for j=1:3
            qq(j,i)=-stl(j,i);
            qq(j,i+3)=-ttl(j,i);
            qq(j+3,i)=0;
            qq(j+3,i+3)=0;
        end
        qq(i+3,i)=1;
    end
    
    % line 402
    %
    %  solve eigenvalue problem for polarization vectors and vertical slownesses
    %  matrix system is nonsymmetric real valued
    %  solution from the eispack guide
    
    
    % call balanc(6,6,qq,is1,is2,fv)
    zi = zeros(6,6);
    [iv, wr, wi, zr, ierr] = balanc(qq);
    if ierr ~= 0
        disp('error');
        return;
    end
    
    % call elmhes(6,6,is1,is2,qq,iv)
    % call eltran(6,6,is1,is2,qq,iv,zr)
    % call hqr2(6,6,is1,is2,qq,wr,wi,zr,ierr)
    % if ierr ~= 0
    %     print *, ierr,'   error!'
    %     return
    % end
    % call balbak(6,6,is1,is2,fv,6,zr)
    %        print *,'for layer',n
    %        print *, 'for phase velocity',cc,'  the vertical slownesses are'
    %        print 101,(wr(i),wi(i),i=1,6)
    %        pause
    % 101 format(6g12.4)
    %  eigenvector unpacking, see EISPACK guide, page 88
    %  bad eigenvector order is flagged by wi(i)>0. for odd i
    
    % ww = eig(qq);
    % wr = real(ww);
    % wi = imag(ww);
    % [zz , ~] = eig(qq);
    % disp(zz)
    % for i=1:6
    %     smq = norm(zz(4:6,i));
    %     zz(:,i)=zz(:,i)/smq;
    % end
    % zr = real(zz);
    % zi = imag(zz);
    iflag =0;
    for i=1:6
        if wi(i)==0
            if n==1
                iev=0;
            end
            zi(:,i)=0;
        elseif wi(i) > 0
            %  bad eigenvector order is flagged by wi(i)>0 for even i
            if mod(i,2)==0
                iflag=iflag+1;
                iv(iflag)=i;
            end
            zi(:,i)=zr(:,i+1);
            % line 436
            % for j=1:6
            %     zi(j,i)=zr(j,i+1)
            % end
        else
            % line 450
            zi(:,i)=-zi(:,i-1);
            zr(:,i)=zr(:,i-1) ;
        end
        %  normalize by the last three indices
        smq = sqrt(sum(zr(4:6,i).^2+zi(4:6,i).^2));
             
        % for j=4:6
        %     sum=sum+zr(j,i) ^2+zi(j,i) ^2
        % end;
        zr(:,i)=zr(:,i)/smq;
        zi(:,i)=zi(:,i)/smq;
        % line 451
        % for j=1:6
        %     zr(j,i)=zr(j,i)/sum
        %     zi(j,i)=zi(j,i)/sum
        % end
    end
    %  assemble the stress-displacement vectors
    %  calculate the traction components, with i removed
    
    u0 = zeros(3,1);
    e1 = zeros(6,6);
    pp(1)=complex(px,0);
    pp(2)=z0;
    for k=1:6
        pp(3)=complex(wr(k),wi(k));
        for i=1:3
            u0(i)=complex(zr(i+3,k),zi(i+3,k));
        end
        pu=z0;
        pw=z0;
        uw=z0;
        abcde=a-b+c-2*d+2*e;
        bce=b-4*c-4*e;
        de=d-e;
        for i=1:3
            pu=pu+pp(i)*u0(i);
            pw=pw+pp(i)*w(i,n);
            uw=uw+u0(i)*w(i,n);
        end
        for i=1:3
            e1(i,k)=u0(i);
            e1(i+3,k)=w(i,n)*(pu*w(3,n)*bce+8*pw*uw*w(3,n)*c ...
                +2*(pw*u0(3)+uw*pp(3))*e);
            e1(i+3,k)=e1(i+3,k)+pp(i)*(u0(3)*de+2*uw*w(3,n)*e);
            e1(i+3,k)=e1(i+3,k)+u0(i)*(pp(3)*de+2*pw*w(3,n)*e);
        end
        e1(6,k)=e1(6,k)+pu*abcde+pw*uw*bce;
        %  almost lastly, mult traction by i
        for i=1:3
            e1(i+3,k)=Eye*e1(i+3,k);
        end
    end
    %  reorder into upgoing and downgoing waves
    %  we use the exp(-i*omega*t) convention with z increasing downward
    %  so downgoing oscillatory waves have p_z>0, k_z real
    %  downgoing evanescent waves have Im(p_z)>0
    %  if the axis of symmetry is tilted, there are cases where a pair of
    %  near-horizontal plane waves will be both upgoing or both downgoing
    %  since Chen's algorithm depends on a 3,3 split, we must adopt a kluge
    %  similarly, there are cases where the EISPACK routines dont return
    %  the vertical wavenumbers in ordered pairs, but mix them up a bit
    %  this seems to cause problems, so a fix is necessary
    %
    %  first, test for bad eigenvector order, switch k-1->k+1, k->k-1, k+1->k
    %   worst case is iflag=2, real,imag1+,imag1-,imag2+,imag2-,real
    if iflag > 0
        for i=1:iflag
            k=iv(i);
            wrr=wr(k-1);
            wii=wi(k-1);
            wr(k-1)=wr(k);
            wi(k-1)=wi(k);
            wr(k)=wr(k+1);
            wi(k)=wi(k+1);
            wr(k+1)=wrr;
            wi(k+1)=wii;
            for j=1:6
                pu=e1(j,k-1);
                e1(j,k-1)=e1(j,k);
                e1(j,k)=e1(j,k+1);
                e1(j,k+1)=pu;
            end
        end
    end
    %  second, divide into upgoing and downgoing waves
    iv = zeros(6,1);
    isum=0;
    for k=1:6
        iv(k)=0;
        if wi(k)==0 && wr(k) > 0
            iv(k)=1;
        end
        if wi(k) > 0
            iv(k)=1;
        end
        isum=isum+iv(k);
    end
    %  if up and downgoing cohorts are not equal, switch the sense of the
    %  pure-oscillatory wave with smallest wavenumber
    % line 530
    % 140   continue
    while isum ~= 3
        wr0=max(abs(wr));
        % for k=1:6
        %     wr0=dmax1(wr0,dabs(wr(k)))
        % end
        for k = 1:6
            if wi(k)==0
                if abs(wr(k)) < wr0
                    wr0=abs(wr(k));
                    kk=k;
                end
            end
        end
        if iv(kk)==0
            iv(kk)=1;
        else
            iv(kk)=0;
        end
        %  check that we have equal up/down cohorts
        isum=sum(iv);
        % line 552
        % for k=1:6
        %     isum=isum+iv(k)
        % end
        % go to 140
    end
    % line 556
    jdown=1;
    jup=4;
    %        print *,'for layer',n,'  the vert wavenums are (0=up,1=dn)'
    % 1001 format(i2,2g15.6)
    for k=1:6
        if iv(k)==1
            ki=jdown;
            jdown=jdown+1;
        else
            ki=jup;
            jup=jup+1;
        end
        for i=1:6
            ee(i,ki,n)=e1(i,k);
        end
        %  incorporate the factor of i into the stored vertical slowness
        xnu(ki,n)=complex(-wi(k),wr(k));
    end
    %  OK, here's where we check whether two downgoing stress-disp vectors
    %  are nearly parallel - we check the dotproducts of displacement components
    % line 576
    idfct = zeros(4,nlp);
    adf = zeros(2,nlp);
    for i=1:2
        for j=i+1:3
            r1=0;
            r2=0;
            zz=z0;
            for k=1:3
                r1=r1+abs(ee(k,i,n))^2;
                r2=r2+abs(ee(k,j,n))^2;
                zz=zz+ee(k,j,n)*conj(ee(k,i,n));
            end
            qqq=1-abs(zz)/sqrt(r1*r2);
            propag.xnu = xnu;
            % line 591
            if qqq < tol
                %ccc=cc*vbar;
                idfct(1,n)=i;
                idfct(2,n)=j;
                %              print 1008,'vert slownesses',xnu(i,n),' and',xnu(j,n)
                %  we average eigenvalues (vert slownesses)
                %  and solve for eigenvector in subroutine defective
                xnu(i,n)=(xnu(i,n)+xnu(j,n))/2;
                xnu(j,n)=xnu(i,n);
                %  calculate the extravector for defective repeated eigenvalue
                propag.xnu = xnu;
                pstff.ee = ee;
                adf(1,n) =  defective(i,j,n,a,b,c,d,e,px);
                %              print *,i,j,n,ccc,qqq,adf(1,n)
                xnu = propag.xnu;
                ee = pstff.ee;
                
            end
        end
    end
    % 1008 format(a,2g15.6,a,2g15.6)
    %  OK, here's where we check whether two upgoing stress-disp vectors
    %  are nearly parallel - we check the dotproducts of displacement components
    % line 609
    for i=4:5
        for j=i+1:6
            r1=0;
            r2=0;
            zz=z0;
            for k=1:3
                r1=r1+abs(ee(k,i,n))^2;
                r2=r2+abs(ee(k,j,n)) ^2;
                zz=zz+ee(k,j,n)*conj(ee(k,i,n));
            end
            qqq=1-abs(zz)/sqrt(r1*r2);
            if qqq < tol
                %ccc=cc*vbar
                idfct(3,n)=i;
                idfct(4,n)=j;
                %              print 1008,'vert slownesses',xnu(i,n),' and',xnu(j,n)
                %  we average the eigenvalues
                xnu(i,n)=(xnu(i,n)+xnu(j,n))/2;
                xnu(j,n)=xnu(i,n);
                propag.xnu = xnu;
                pstff.ee = ee;
                %  calculate the extravector for defective repeated eigenvalue
                %  as well as coefficient (adf) of linear (z-z0) term
                adf(2,n) =  defective(i,j,n,a,b,c,d,e,px);
                %              print *,i,j,n,ccc,qqq,adf(2,n)
                ee = pstff.ee;
                xnu = propag.xnu;
            end
        end
    end
end
%  now, must identify which upgoing waves in the halfspace are P,SV,SH
%  crud, this goes back to array ee
%  3: SH is y-motion
%  2: SV is (-sqrt((1/vs) ^2-p_x ^2),0,-p_x)  ! recall that z points down
%  1: P is (p_x,0,-sqrt((1/vp) ^2-p_x ^2)
%  so we branch on size of u_y, and relative sign of u_x and u_z

disp('in the halfspace:');
for i=4:6
    fprintf('for i*k_z= (%f , %f)  the the disp-stress vector is \n', real(xnu(i,nlp)), imag(xnu(i,nlp)));
    % print *,'for i*k_z=',xnu(i,nlp),', the disp-stress vector is'
    xr = real(ee(:,i,nlp));
    xi = imag(ee(:,i,nlp));
    % for j=1:6
    %     xi(j)=imag(ee(j,i,nlp));
    %     xr(j)=real(ee(j,i,nlp));
    % end
    for jj = 1:6
        fprintf('%f  %f\n', xr(jj), xi(jj));
    end
    % print 101,(xr(j),j=1,6),(xi(j),j=1,6)
end
ips = zeros(3,1);
for i=4:6
    ips(i-3)=3;
    if abs(ee(2,i,nlp)) < sqrt(tol) % not SH
        test=real(ee(1,i,nlp))/real(ee(3,i,nlp));
        if test > 0
            ips(i-3)=2;
        else
            ips(i-3)=1;
        end
    end
end
qstff.ips = ips;

propag.xnu = xnu;
defect.idfct = idfct;
pstff.ee = ee;
pstff.e1 = e1;
model.rho = rho;
fprintf('wave types: %d  %d %d\n', ips(1), ips(2), ips(3));
% print *,'wave types:',(ips(i),i=1,3)
