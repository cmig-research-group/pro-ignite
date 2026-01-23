function X=gaussian(fwhm, tsize)
% DblMat X(size,1,0.0);
%   
%   double factor = 4.0*log(2.0)/(fwhm*fwhm);
%   double scale = 2.0*sqrt(log(2.0)/M_PI) / fwhm;
%   X(0,0) = scale;
%   
%   for(int i = 1; i <= (size-1)/2; i++)
%     {
%       X(i,0) = X(size-i,0) = scale*exp(-i*i*factor);
%     }
%   
%   if((size-1)/2 != size/2) // if size is even
%     // fill in middle value
%     X(size/2,0) = scale*exp(-size*size*factor/4.0);
% 
%   return(X);

X=zeros(tsize,1);
factor = 4.0*log(2.0)/(fwhm*fwhm);
scale = 2.0*sqrt(log(2.0)/pi) / fwhm;
X(1)=scale;

for i=1:(tsize-1)/2
    X(i+1)=scale*exp(-i*i*factor);
    X(tsize-i+1)=X(i+1);
end
if (floor((tsize-1)/2) ~= floor(tsize/2))
    X(tsize/2+1) = scale*exp(-tsize*tsize*factor/4.0);
end