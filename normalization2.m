function [normz] = normalization2(P)
%this function will normalize data by divide to its maximum/minimum value.
[x y]=size(P)
%format long 
for r=1:y
    vm1=max(P(:,r));
    vm2=min(P(:,r));
    if vm1>(-1*vm2),
    P(:,r)=(1/vm1).*P(:,r);
    else
    P(:,r)=(1/(-1*vm2)).*P(:,r);    
    end
end
normz=P;