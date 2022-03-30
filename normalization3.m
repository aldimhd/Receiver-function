function [normz] = normalization3(A)
%this function will normalize data by divide to its maximum/minimum value.
[x y]=size(A)
%format long 
for r=1:y
    vm1=max(A(:,r));
    vm2=min(A(:,r));
    if vm1>(-1*vm2),
    A(:,r)=(1/vm1).*A(:,r);
    else
    A(:,r)=(1/(-1*vm2)).*A(:,r);    
    end
end
normz=A;