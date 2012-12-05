classdef Searchlight < handle
    % sl = Searchlight(maskpath,mapmode,radvox) given a maskpath, a mapmode
    % ('radius' or 'nvox') and a parameter for the mapmode (mm for radius),
    % Searchlight constructs a mapper that can then be queried with the
    % mapcoords(xyz) or mapinds(n) methods to return linear indices to a
    % searchlight for a given location (or optionally a binary mask instead
    % if outputmask=1)
    %
    % the nvox mapmode operates by fminsearch and so runs slower than
    % radius-based searchlights.
    properties
        mapmode = 'radius' % radius or nvox
        radius = []; % searchlight radius (mm)
        radius_vox % searchlight radius (voxels)
        nvox = []; % number of voxels in searchlight
        nwantedvox = []; % only defined if doing mapping by nvox
        nwantedtolerance % accept searchlights that this far from nwanted
        stepsize_radius = .05; % find target in these increments
        searchsphere % searchsphere for current radius in binary form
        spherecoords % coordinates for searchsphere centered on 0
        spherenvox % n voxels in ideal sphere (not masked)
        outputmask = 0 % methods return column indices or binary mask
        vol % Volume instance for mask
        xyz % coordinates of last sphere
        optimoptions = struct('TolX',1e-100,'TolFun',1e-100,'MaxIter',1e3);
    end

    methods
        function sl = Searchlight(maskpath,mapmode,radvox)
        % sl = Searchlight(maskpath,mapmode,radvox)
            sl.vol = Volume([],maskpath);
            % decide on how we are going to map
            sl.mapmode = lower(mapmode);
            switch sl.mapmode
                case 'radius'
                    sl.radius = radvox;
                case 'nvox'
                    % always define a radius because otherwise we can't do
                    % the initial makesphere to initialise the searchsphere
                    sl.radius = 5;
                    sl.nwantedvox = radvox;
                    sl.nwantedtolerance = ceil(radvox * .1);
                otherwise
                    error('unknown mapmode: %s',mapmode)
            end
            % make initial searchsphere. If we are doing radius-based
            % mapping we will never have to call this method again
            sl.makesphere;
        end

        function makesphere(self)
        % makesphere(self)
            self.radius_vox = self.radius ./ self.vol.voxsize;
            % by thresholding a meshgrid we obtain a discontinuous
            % searchsphere
            % (or well, a cube with 1 inside searchsphere)
            [x,y,z] = meshgrid(-self.radius_vox(1):self.radius_vox(1), ...
                -self.radius_vox(2):self.radius_vox(2), ...
                -self.radius_vox(3):self.radius_vox(3));
            self.searchsphere = (x*self.vol.voxsize(1)).^2 + ...
                (y*self.vol.voxsize(2)).^2 + ...
                (z*self.vol.voxsize(3)).^2 <= self.radius^2;
            % get size of searchsphere (cube) (extra dribbling here is
            % needed to prevent Matlab from squeezing last dim for very
            % small (z=1) spheres)
            sphsize_vox = [size(self.searchsphere), ...
                ones(1,3-ndims(self.searchsphere))];
            % get coordinates of the searchsphere inside the cube
            [sphereSUBx,sphereSUBy,sphereSUBz]=ind2sub(sphsize_vox, ...
                find(self.searchsphere)); 
            sphereSUBs=[sphereSUBx,sphereSUBy,sphereSUBz];
            % find indices for centre
            ctrSUB=sphsize_vox/2+[.5 .5 .5]; 
            % shift coordinates so centre is 0
            % (and floor to avoid half voxels)
            self.spherecoords=floor(sphereSUBs - ...
                ones(size(sphereSUBs,1),1)*ctrSUB); 
            self.spherenvox = size(self.spherecoords,1);
        end

        function out = mapcoords(self,xyz)
        % Create a searchlight sphere for a given xyz coordinate (voxel
        % index). Returns column indices into a mapped dataset by default,
        % optionally a 3D mask instead if self.outputmask=1.
        % out = mapcoords(self,xyz)
        %
            assert(self.vol.mask(xyz(1),xyz(2),xyz(3))==1,...
                'coordinates %d are outside mask',xyz);
            self.xyz = xyz;
            % find voxel coordinates for current searchlight (add sphere
            % coordinates to xyz)
            coords = repmat(xyz,[self.spherenvox 1])+self.spherecoords;
            % find searchlight coordinates outside volume limits
            insidev = all(coords>=1,2) & all(coords<=repmat(...
                self.vol.V.dim,[self.spherenvox 1]),2);
            coords_insidev = coords(insidev,:);
            % find coordinates inside mask
            % first back to linear indices
            lininds_insidev = self.vol.coord2linind(coords_insidev);
            % then restrict to mask lininds
            lininds_out = intersect(lininds_insidev,self.vol.lininds);
            % final estimate of how many voxels we ended up with
            self.nvox = length(lininds_out);
            % maybe we have to recurse in here (isselfcall stops infinity)
            if ~isempty(self.nwantedvox) && ...
                    (abs(self.nvox-self.nwantedvox)>0) && ~isselfcall
                opts = optimset(self.optimoptions);
                [r,errmargin,exflag,output] = fminsearch(@self.findrfornvox,...
                    self.radius,opts);
                assert(exflag==1,'optimisation did not converge!')
                % and generate output
                out = self.mapcoords(xyz);
                assert(abs(self.nvox-self.nwantedvox) < ...
                    self.nwantedtolerance,'optimisation failed!');
                %% find the right radius
                %ofun = @(m) abs(-self.nwantedvox);
                %if self.nvox < self.nwantedvox
                    %self.radius = self.radius + self.stepsize_radius;
                %else
                    %self.radius = self.radius - self.stepsize_radius;
                %end
                % re-do the sphere
                %self.makesphere;
                % and recursively call again
                % nb, if you get crashes due to recursion limit here you
                % likely need to have a higher nwantedtolerance
                % A far more efficient solution here would be some kind of
                % gradient descent. The current linear search is quite
                % dumb.
                %out = self.mapcoords(xyz);
            else
                if self.outputmask
                    % make the binary mask and populate with the final
                    % sphere
                    out = self.vol.data2mat(lininds_out);
                else
                    % Find the feature indices 
                    out = self.vol.linind2featind(lininds_out);
                end
            end
        end

        function out = mapinds(self,ind)
        % out = mapinds(self,ind)
            % wrapper around mapcoords
            out = self.mapcoords(self.vol.linind2coord(self.vol.lininds(ind)));
        end

        function ndiff = findrfornvox(self,r)
        % ndiff = findrfornvox(self,r)
        % used internally with optimset
            self.radius = r;
            self.makesphere;
            out = self.mapcoords(self.xyz);
            ndiff = abs(self.nwantedvox-self.nvox);
        end
    end
end
