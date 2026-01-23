function vols=sharpen_volume(vol, volm, Nbins, fwhm, noise)
% 
% Sharpen Volume by deconvolve with fixed guassian kernel
%
%   vols=sharpen_volume(vol, volm, Nbins, fwhm, noise)
%
%       vol: volume strcture
%       volm: mask volume
%       Nbins: number of bins for histogram
%       fwhm: fwhm for gaussian kernel
%       noise: noise for wiener filtering

% find mask indices
ind=find(volm.imgs>0);


% building histogram
u_h=vol.imgs(ind);
des=[min(u_h) max(u_h) Nbins];
dbin=(des(2)-des(1))/(Nbins-1);
bins=des(1):dbin:des(2);
U=n3pdfe(u_h,bins);

% Sharpen Volume
padded_size= floor(2^(ceil(log(Nbins)/log(2.0))+1)+0.5);
offset=(padded_size - Nbins)/2;

% Create blur kernel
blur=fft(gaussian(fwhm/dbin, padded_size));

% create wiener filter kenerl
wfilter=conj(blur)./(conj(blur).*blur+noise);

U_pad=zeros(padded_size,1);
U_pad(offset+1:offset+Nbins)=U;

% estimate U
f=real(ifft(fft(U_pad).*wfilter));

indf=find(f<0);
f(indf)=0;


% sharpen U
moment=(bins(1) + ([0:padded_size-1]'-offset)*dbin).*f;
Y_pad=real(ifft(fft(moment).*blur))./real(ifft(fft(f).*blur));
Y=Y_pad(offset+1:offset+Nbins);
indy=find(Y==nan);
Y(indy)=0;
u_h_e=interp1(bins, Y, u_h,'linear');
 
% create sharpen volume
vols=vol;
vols.imgs(ind)=u_h_e;
[vols.maxI vols.minI]=maxmin(vols.imgs);