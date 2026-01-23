function V=n3pdfe(x, bins)
%
% histogram estimations by using Parzen triangle window
%   
%   V=n3pdfe(x, bins)
%
%   x: value vector
%   bins: bin vector

N=length(x);
n=length(bins);
h=bins(2)-bins(1);
V=zeros(size(bins));
for i=1:n
    s=(bins(i)-x)/h;
    phi=zeros(size(s));
    ind=find(abs(s)<1);
    phi(ind)=1-abs(s(ind));
    V(i)=sum(phi);
end

