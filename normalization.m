function [normz] = normalization(RFtmp3f)
%this function will normalize data by divide to its maximum/minimum value.
[x y]=size(RFtmp3f)
%format long 
for r=1:y
    vm1=max(RFtmp3f(:,r));
    vm2=min(RFtmp3f(:,r));
    if vm1>(-1*vm2),
    RFtmp3f(:,r)=(1/vm1).*RFtmp3f(:,r);
    else
    RFtmp3f(:,r)=(1/(-1*vm2)).*RFtmp3f(:,r);    
    end
end
normz=RFtmp3f;