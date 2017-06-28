% extract age and gender from dicom in input dcmdir (or perhaps from one of the
% series sub-dirs inside dcmdir, if subdirpattern is not empty).
%
% [age,gender] = dcm2demographic(dcmdir,[subdirpattern='Series_*'])
function [age,gender] = dcm2demographic(dcmdir,subdirpattern)

if ~exist('subdirpattern','var')
    subdirpattern = 'Series_*';
end

if ~isempty(subdirpattern)
    subdirhit = dir(fullfile(dcmdir,subdirpattern));
    % pick one
    dcmdir = fullfile(dcmdir,subdirhit(1).name);
end

dcms = dir(fullfile(dcmdir,'*.dcm'));
assert(~isempty(dcms),'no dcm files found in %s',dcmdir);
% pick one
h = spm_dicom_headers(fullfile(dcmdir,dcms(1).name));
gender = h{1}.PatientSex;
age = h{1}.PatientAge;
