function fx= refft(xr,xi)
n2= length(xr);
n = n2*2;
sn = 0.0;
cn = 1.0;
arg = pi / n;
aa = sin (arg);
cd = 2.0 * aa * aa;
sd = -sin (arg+arg);

aa = xr(1);
bb = xi(1);
xr(1) = (aa + bb)/2;
xi(1) = (bb - aa)/2;

i = 2;
j = n2;
while i<=j
    aa = cd * cn + sd * sn;
    sn = sn + sd * cn - cd * sn;
    cn = cn-aa;
    aa = (xr(i) + xr(j)) * .5;
    ab = (xr(i) - xr(j)) * .5;
    ba = (xi(i) + xi(j)) * .5;
    bb = (xi(i) - xi(j)) * .5;
    p_real = cn * ba + sn * ab;
    p_imag = sn * ba - cn * ab;
    xi(j) = p_imag - bb;
    xi(i) = p_imag + bb;
    xr(j) = aa - p_real;
    xr(i) = aa + p_real;
    i=i+1;
    j=j-1;
end
i=0;
j=0;
step=0;
sintab = sin(pi./2.^(1:20));

i=0;
% xp = x(1), x(2), ...s
for xp=1:n2
    yp = i+1;
    if xp<yp
        temp = xr(yp);
        xr(yp) = xr(xp);
        xr(xp) = temp;
        temp = xi(yp);
        xi(yp) = xi(xp);
        xi(xp) = temp;
    end
    j=n2/2;
    while j>=1 && i>=j
        i = i-j;
        j = j/2;
    end
    i = i+j;
end
xr=xr/n2;
xi=xi/n2;
for xp=1:2:n2
    yp=xp+1;
    temp = xr(yp);
    xr(yp) = xr(xp)-temp;
    xr(xp)=xr(xp)+temp;
    temp = xi(yp);
    xi(yp) = xi(xp)-temp;
    xi(xp)=xi(xp)+temp;
end
[~,p]=find(2.^(2:20)==n);
si=1;
for step=2.^(2:p)
    i=step/2;
    sd=-sintab(si);
    temp=sintab(si+1);
    si=si+1;
    cd=2*temp*temp;
    cn=1;
    sn=0;
    for j=0:i-1
        for xp=j:step:n2-1
            yp=xp+i;
            rl = cn*xr(yp+1)-sn*xi(yp+1);
            im = sn*xr(yp+1)+cn*xi(yp+1);
            xr(yp+1)=xr(xp+1)-rl;
            xi(yp+1)=xi(xp+1)-im;
            xr(xp+1)=xr(xp+1)+rl;
            xi(xp+1)=xi(xp+1)+im;
        end
        temp = cd * cn + sd * sn;
        sn = sn+ sd * cn - cd * sn;
        cn =cn-temp;
    end
end
xi(1:n2)=-xi(1:n2);
fx=[xr xi]';
fx=fx(:);
end

