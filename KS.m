function [B] = KS(A)
n=size(A,1);
for i=1:n
    r=floor(rand*(n+1-i))+1;
    a=A(n+1-i,:);
    b=A(r,:);
    A(r,:)=a;
    A(n+1-i,:)=b;
end
B=A;
end

