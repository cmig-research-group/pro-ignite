function vol=shrink_volume(volo, shrink, varargin)
%
%   Shrinking volume  
%     vols=shrink_volume(vol, shrink, [interpm], [padding])
%
%       vol: volume struct
%       shrink: shrink factor
%       interpm : 0:Nearest 1:Linear(default) 2:cubic 
%                 3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%       padding : half width of number of points used in the interpolation,
%                  only for 4 and 5. For example, 3 means it would use
%                  6*6*6 points for interpolation.


interpm = 1;
if nargin >= 3
  interpm = varargin{1};
end


padding=3;
if nargin >= 4
  padding = varargin{2};
end

if (interpm==0) % Nearest Neighbor padding=0
    padding=0;
elseif (interpm==1) % Lineae  padding=1
    padding=1;
elseif (interpm==2) % cubic and key's padding =2
    padding=2;
elseif (interpm==3)
    padding=2;
end

% Kepp all geometric info and shrink the dim
vol.dimc=ceil((volo.dimc-1)/shrink)+1;
vol.dimr=ceil((volo.dimr-1)/shrink)+1;
vol.dimd=ceil((volo.dimd-1)/shrink)+1;
vol.vx=volo.vx*shrink;
vol.vy=volo.vy*shrink;
vol.vz=volo.vz*shrink;
vol.DirRow=volo.DirRow;
vol.DirCol=volo.DirCol;
vol.DirDep=volo.DirDep;
vol.lphcent=volo.lphcent;
Mvxl2lph=eye(4,4);
% Row direction vector
Mvxl2lph(:,2)=vol.vx*[vol.DirRow;0];
%Column Direction vector
Mvxl2lph(:,1)=vol.vy*[vol.DirCol;0];
%Slice direction vector
Mvxl2lph(:,3)=(vol.vz)*[vol.DirDep;0];
%For matlab Index sake: start from 1
T=Mvxl2lph*[-(vol.dimr+1)/2; -(vol.dimc+1)/2; -(vol.dimd+1)/2;1];
Mvxl2lph(:,4)=T+[vol.lphcent;0];
vol.Mvxl2lph=Mvxl2lph;
vol.imgs= zeros(vol.dimr,vol.dimc,vol.dimd);
vol=vol_resample(volo,vol,eye(4,4), interpm, padding);
