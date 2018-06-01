% Converts a set of volumes (spm V struct) to an animated gif.
%
% V is a cell/char array of paths to SPM-readable files OR spm_vol structs
% outfile is where the gif gets saved
% fps (default 6) is the framerate
% view (default axi or axial) the slice plane (cor / sag / axi)
%
% volumes2gifanimation(V,outfile,[fps],[view])
function volumes2gifanimation(V,outfile,fps,view)

if iscell(V)
    % attempt to concatenate to char
    V = vertcat(V{:});
end

if ischar(V)
    % load volume
    V = spm_vol(V);
end
% ... otherwise we assume you provided a volume

if ~exist('view','var') || isempty(view)
    view = 'axial';
end

if ~exist('fps','var') || isempty(fps)
    fps = 5;
end

% now save memory by only loading one vol at a time
for v = 1:length(V);
    xyz = spm_read_vols(V(v));
    switch view(1:3)
        case 'sag'
            xyz = flipdim(permute(xyz,[3 2 1]),1);
        case 'cor'
            xyz = flipdim(permute(xyz,[3 1 2]),1);
        otherwise
            assert(strcmp(view(1:3),'axi'),'unknown view: %s',view);
    end
    im = makeimagestack(xyz,-2,2,[],2);
    im = uint8(im * 256);
    [imind,cm] = rgb2ind(cat(3,im,im,im),256);
    if v == 1
        % initialise gif
        imwrite(imind,cm,outfile,'gif','loopcount',0,'delaytime',1/fps);
    else
        % append frame to gif
        imwrite(imind,cm,outfile,'gif','writemode','append',...
            'delaytime',1/fps);
    end
end
