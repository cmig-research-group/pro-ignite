function volc=volvxlclamp(vol, Imin, Imax)
% Volume Vxl Clamping
% volc=volvxlclamp(vol, min, max)

volc=vol;
clear vol;
ind=find(volc.imgs<Imin);
volc.imgs(ind)=Imin;
ind=find(volc.imgs>Imax);
volc.imgs(ind)=Imax;
[volc.maxI volc.minI]=maxmin(volc.imgs);
