% data = smoothdatavecs(data,fwhm,mask,voxsize)
function data = smoothdatavecs(data,fwhm,mask,voxsize)

if fwhm==0
    % quick and easy
    return
end
parfor dat = 1:size(data,1)
    % make matrix with NaNs outside mask
    datmat = datavec2mat(data(dat,:),mask,'NaN');
    % use Kendrick's smoothing fun that handles NaNs properly
    datmat = smoothvolumes(datmat,voxsize,[fwhm fwhm fwhm]);
    % return to data
    data(dat,:) = datmat(mask);
end
