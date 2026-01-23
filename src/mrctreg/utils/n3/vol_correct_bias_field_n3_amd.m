function [vol_bc, f_vol]=vol_correct_bias_field_n3(volo,estop,maxiter,Nbins,fwhm,noise,distance,shrink,th)
%
% N3 Bias field correction
% 
%  REF:  J. G. Sled, A. P. Zijdenbos, and A. C. Evans, 
%  'A non-parametric method for automatic correction of intensity 
%   non-uniformity in MRI data,'
%   IEEE Transactions on Medical Imaging, vol. 17, pp. 87-97, February 1998.
%
%  [vol_bc, f_vol]=vol_correct_bias_field_n3(vol, [estop], [maxiter], [Nbins]
%                                       [fwhm], [noise], [distance], [shrink]
%
%  INPUT:
%
%   vol: Volume strcuture
%   estop: stop criterion (default=0.001)
%   maxiter: maximum iteration (default=50)
%   Nbins: number of bins for histogram (default=200)
%   fwhm: Full-Width Half-Maximum for gaussian kernel (default = 0.15)
%   noise: noise for wiener filerting (default = 0.01)
%   distance: spline knot distance unit:mm (default = 200)
%   shrink: shrinking factor for working volume (default = 4)
%
%  OUPUT:
%   
%   vol_bc: bias field corrected volume
%   f_vol: estimated bias field volume

if ~exist('estop','var') || isempty(estop)
  estop=0.001;
end

if ~exist('maxiter','var') || isempty(maxiter)
  maxiter=50;
end

if ~exist('Nbins','var') || isempty(Nbins)
  Nbins=200;
end

if ~exist('fwhm','var') || isempty(fwhm)
  fwhm=0.15;
end

if ~exist('noise','var') || isempty(noise)
  noise=0.01;
end

if ~exist('distance','var') || isempty(distance)
  distance=200;
end

if ~exist('shrink','var') || isempty(shrink)
  shrink=4;
end

% Shrinking volume
vol=shrink_volume(volo, shrink, 0);

%clamping 
vol=volvxlclamp(vol, 1, 1.7e308);

%log
vol_log=vol;
vol_log.imgs=log(vol.imgs);
[vol_log.maxI vol_log.minI]=maxmin(vol_log.imgs);

%Getting threshhold and creating mask;
if ~exist('th','var') || isempty(th)
  th=vol_getBiModalThreshold(vol);
end

vol_mask=vol;
vol_mask.imgs=zeros(size(vol.imgs));
ind=find(vol.imgs>th);
vol_mask.imgs(ind)=1;
[vol_mask.maxI vol_mask.minI]=maxmin(vol_mask.imgs);


%apply mask to log volume
vol_log.imgs=vol_log.imgs.*vol_mask.imgs;
[vol_log.maxI vol_log.minI]=maxmin(vol_log.imgs);



% Initialization
vol_u_h=vol_log;
vol_f_h=vol_log;
vol_f_h.imgs=zeros(size(vol_f_h.imgs));
vol_f_h_s=vol_f_h;
lamda=1/shrink^3;

for i=1:maxiter
    % sharpen volume
    vol_u_h.imgs=vol_log.imgs-vol_f_h_s.imgs;
    vol_u_h_e=sharpen_volume(vol_u_h, vol_mask, Nbins, fwhm, noise);

    % Calculate f_h
    vol_f_h.imgs(ind)=vol_log.imgs(ind)-vol_u_h_e.imgs(ind);
    [vol_f_h.maxI vol_f_h.minI]=maxmin(vol_f_h.imgs);

    % Smoothing f_h Volume with tensor cubic spline
    vol_f_h_s_old=vol_f_h_s;
    vol_f_h_s=spline_smooth(vol_f_h, vol_mask, lamda, distance);
    e=std(vol_f_h_s_old.imgs(ind)-vol_f_h_s.imgs(ind));
    %[i e]
    if e <estop
        break
    end
end

% Get Bias field
vol_f=vol_f_h_s;
vol_f.imgs=exp(vol_f_h_s.imgs);
[vol_f.maxI vol_f.minI]=maxmin(vol_f.imgs);
[vol_f_s, spt]=spline_smooth(vol_f, vol_mask, 1, distance);
f_vol=sptvol(volo, spt.domain, spt.distance, spt.coef);

% Correct Volume
vol_bc=volo;
ind = find(f_vol.imgs>1e-5); % AMD: should make this a parameter?
vol_bc.imgs(ind)=volo.imgs(ind)./f_vol.imgs(ind);
%vol_bc.imgs=volo.imgs./f_vol.imgs;

%[vol_bc.maxI vol_bc.minI]=maxmin(vol_bc.imgs);
