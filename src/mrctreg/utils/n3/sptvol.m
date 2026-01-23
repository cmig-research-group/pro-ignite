function vols=sptvol(vol, domain, distance, coef)

sizes=size(vol.imgs);
vols=vol;
vols.imgs=zeros(size(vol.imgs));

separations=zeros(3,1);
separations(1)=vol.vx;
separations(2)=vol.vy;
separations(3)=vol.vz;

[spline1d offset nProduct n]=getsplinesMEX(domain, separations, distance, sizes);
bmask=0;
vols.imgs=getbiasfieldvolMEX(sizes, coef, size(spline1d), spline1d, offset, bmask, []);
[vols.maxI vols.minI]=maxmin(vols.imgs);
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
% % Create Lookup Table
% spline1d=zeros(nDim, max(sizes), 4);
% offset=zeros(nDim, max(sizes),4);
% 
% for i=1:nDim
%     offsetStep=1;
%     for j=1:nDim-i
%         offsetStep=offsetStep*n(nDim-j);
%     end
%     %offsetStep
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
% 
% 
% for x=1:sizes(2)
%     for y=1:sizes(1)
%         for z=1:sizes(3)
%             
%             xp=squeeze(spline1d(1,x,:));
%             yp= squeeze(spline1d(2,y,:));
%             zp= squeeze(spline1d(3,z,:));
%             xop= squeeze(offset(1,x,:));
%             yop= squeeze(offset(2,y,:));
%             zop= squeeze(offset(3,z,:));
%             val=0;
%             for i=1:4
%                 xv=xp(i);
%                 xo=xop(i);
%                 for j=1:4
%                     xyv=xv*yp(j);
%                     xyo=xo+yop(j);
%                     for k=1:4
%                         val=val+xyv*zp(k)*coef(xyo+zop(k)+1);
%                     end
%                 end
%             end
%             vols.imgs(y,x,z)=val;
%         end
%     end
% end
%[vols.maxI vols.minI]=maxmin(vols.imgs);