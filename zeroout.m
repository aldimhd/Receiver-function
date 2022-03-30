function [B]=zeroout(A,C)
if C~=1,
A=A';
end
B=A;
B(~any(A,2),:)=[];
if C~=1
    B=B';
end
