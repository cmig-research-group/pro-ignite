function [vols spt]=spline_smooth(vol, volm, lamda, distance)


vols=vol;
vols.imgs=zeros(size(vols.imgs));

%determine domain from size of volume in world coordinates
separations=zeros(3,1);
separations(1)=vol.vx;
separations(2)=vol.vy;
separations(3)=vol.vz;

sizes=size(vol.imgs);

domain=zeros(3,2);
for i=1:3
    domain(i,1)=-0.5*separations(i);
    domain(i,2)=(sizes(i)-0.5)*separations(i);
end

four=4^3;
%domain

% Initial splines


% sd=size(domain);
% nDim = sd(1);
% start=zeros(nDim,1);
% scale = 1.0/(distance*distance*distance);
% 
% 
% % determine number of basis functions in each dimension
% n=zeros(nDim,1);
% for i=1:nDim
%     n(i)= ceil((domain(i,2) - domain(i,1))/distance) + 3;
% end    
% %n
% 
% knots=zeros(max(n)+4,nDim);
% for i=1:nDim
%     starts=0.5*(domain(i,1) + domain(i,2) - distance*(n(i)+3));
%     for j=1:n(i)+4
%         knots(j,i) = starts + distance*(j-1);
%     end
% end
% %knots
% 
% zerok=knots(4,:);
% 
% nProduct=1;four=1;smallN=zeros(4,1);
% for i=1:nDim
%     nProduct=nProduct*n(i);
%     four=four*4;
%     smallN(i)=4;
% end
% %nProduct
% %four
% %smallN
% 
% AtA=zeros(nProduct, nProduct);
% AtF=zeros(nProduct,1);
% coef=zeros(nProduct,1);
% dloci=zeros(four,1);
% dlocj=zeros(four,1);
% 
% % Create Lookup Table
% spline1d=zeros(nDim, max(sizes), 4);
% offset=zeros(nDim, max(sizes),4);
% 
% for i=1:nDim
%     offsetStep=1;
%     for j=1:nDim-i
%         offsetStep=offsetStep*n(nDim-j);
%     end
%     [i offsetStep]
%     
%     if (domain(i,1)>start(i))
%         first=ceil((domain(i,1) - start(i))/separations(i))+1
%     else
%         first=1;
%     end
%     %first
%     if (domain(i,2) <(start(i)+(sizes(i)-1)*separations(i)))
%         last=floor((domain(i,2)-start(i))/separations(i))+1;
%     else
%         last=sizes(i);
%     end
%     %last
%     [first last]
%     x = start(i)+(first-1)*separations(i);
%     
%     for j=first:last
%     
%         blockIndex= ceil((x-zerok(i))/distance)-1;
%         if (blockIndex <0)
%             blockIndex=0;
%         elseif (blockIndex>(n(i)-4))
%             blockIndex=n(i)-4;
%         end
%         
%         currentOffset = blockIndex*offsetStep;
%          offset(i,j,1)=currentOffset;
%         for k=2:4
%             currentOffset=currentOffset+offsetStep;
%             offset(i,j,k)=currentOffset;
%         end
%         
%         temp=scale*(knots(blockIndex+5,i)-x)^3;
%         spline1d(i,j,1)=temp;
%         spline1d(i,j,2)=scale*(knots(blockIndex+6,i)-x)^3-4*temp;
%         temp=scale*(x-knots(blockIndex+4,i))^3;
%         spline1d(i,j,3)=scale*(x-knots(blockIndex+3,i))^3-4*temp;
%         spline1d(i,j,4)=temp;
%         x=x+separations(i);
%     end
%     
% end

[spline1d offset nProduct n]=getsplinesMEX(domain, separations, distance, sizes);

%% Add data point
%ind=find(volm.imgs>0);
%nL=length(ind);
% lower=zeros(3,1);
% upper=zeros(3,1);
% for i=1:nDim
%     lower(i)=ceil(domain(i,1)/separations(i))+1;
%     upper(i)=floor(domain(i,2)/separations(i))+1;
% end
% 
% for x=lower(1):upper(1)
%    for y=lower(2):upper(2)
%        for z=lower(3):upper(3)
%            if (volm.imgs(y,x,z)>0)
%                val=vol.imgs(y, x, z);
%                value=zeros(four,1);
%                location=zeros(four,1);
%                xp=squeeze(spline1d(1,x,:));
%                yp= squeeze(spline1d(2,y,:));
%                zp= squeeze(spline1d(3,z,:));
%                xop= squeeze(offset(1,x,:));
%                yop= squeeze(offset(2,y,:));
%                zop= squeeze(offset(3,z,:));
%                
%                for i=1:4
%                    xv=xp(i);
%                    xo=xop(i);
%                    for j=1:4
%                        xyv=xv*yp(j);
%                        xyo=xo+yop(j);
%                        for k=1:4
%                            value(16*(i-1)+4*(j-1)+k)=xyv*zp(k);
%                            location(16*(i-1)+4*(j-1)+k)=xyo+zop(k);
%                        end
%                    end
%                end
%                dloci(1:four-1)=diff(location);
%                dlocj(1:four-1)=nProduct*dloci(1:four-1);
%                
%                for k=four:-1:1
%                    AtF(location(k)+1)=AtF(location(k)+1)+val*value(k);
%                    AtA(location(k)+1, location(k)+1)=AtA(location(k)+1, location(k)+1)+value(k)^2;
%                    for l=k-1:-1:1
%                        AtA(location(k)+1, location(l)+1)=AtA(location(k)+1, location(l)+1)+value(k)*value(l);
%                    end
%                   AtA(1:location(k), location(k)+1)=AtA(1:location(k), location(k)+1)+value(k)*value(1:k-1);
%                end
%                AtF
%                AtA;
%            end
%        end
%    end
% end

[AtA, AtF]=adddataMEX(domain, separations, four, sizes, vol.imgs, volm.imgs, size(spline1d), spline1d, offset, nProduct);

%Create Bending Energy
%J=bendingEnergyTensor(n);

[J, BM]=bendingEnergyTensorMEX(n, max(n), nProduct);
% Calculate coeff.
%A=csy(AtA)+lamda*J;
A=AtA+lamda*J;
coef=A\AtF;

% calcualte volumeoutputCompactField
% for x=1:sizes(2)
%     for y=1:sizes(1)
%         for z=1:sizes(3)
%             if (volm.imgs(y,x,z)>0)
%                xp=squeeze(spline1d(1,x,:));
%                yp= squeeze(spline1d(2,y,:));
%                zp= squeeze(spline1d(3,z,:));
%                xop= squeeze(offset(1,x,:));
%                yop= squeeze(offset(2,y,:));
%                zop= squeeze(offset(3,z,:));
%                val=0;
%                for i=1:4
%                    xv=xp(i);
%                    xo=xop(i);
%                    for j=1:4
%                        xyv=xv*yp(j);
%                        xyo=xo+yop(j);
%                        for k=1:4
%                            val=val+xyv*zp(k)*coef(xyo+zop(k)+1);
%                        end
%                    end
%                end
%                vols.imgs(y,x,z)=val;
%             end
%         end
%     end
% end
% [vols.maxI vols.minI]=maxmin(vols.imgs);

bmask=1;
vols.imgs=getbiasfieldvolMEX(sizes, coef, size(spline1d), spline1d, offset, bmask, volm.imgs);
spt.coef=coef;
spt.domain=domain;
spt.distance=distance;

