function vol_resampled = resfun_cropped(ctx_src, ctx_target, x)

ctx_resampled = vol_resample(ctx_src, ctx_target, Mtrans(x(1),x(2),x(3)));
vol_resampled = ctx_resampled.imgs;

[rows, cols, ~] = size(vol_resampled);
numel_half_slice = (rows*cols)/2;
mask_zeros = (vol_resampled == 0);
CC = bwconncomp(mask_zeros, 6);
numPixels = cellfun(@numel, CC.PixelIdxList);
[biggest, idx] = max(numPixels);
if ~isempty(biggest) && (biggest > numel_half_slice)
  mask_bg = false(size(vol_resampled));
  mask_bg(CC.PixelIdxList{idx}) = true;
  vol_resampled(mask_bg) = -1000;
end

end
