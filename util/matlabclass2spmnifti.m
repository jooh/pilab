function t = matlabclass2spmnifti(data)

% figure out what precision nifti we need
if isnumeric(data)
    switch class(data)
        case 'double'
            t = 'float64';
        case 'single'
            t = 'float32';
        otherwise
            error('cannot write data of class: %s',class(data));
    end
else
    t = 'int32';
end
t = spm_type(t);
