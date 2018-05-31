classdef Searchlight < handle
    % sl = Searchlight(maskpath,mapmode,radvox) given a maskpath, a mapmode
    % ('radius' or 'nvox') and a parameter for the mapmode (mm for radius),
    % Searchlight constructs a mapper that can then be queried with the
    % mapcoords(xyz) or mapinds(n) methods to return linear indices to a
    % searchlight for a given location.
    % TODO: eventually sub-class for MriSearchlight and other potential
    % lights
    properties
        mapmode = 'radius' % radius or nvox
        radius = []; % searchlight radius (mm)
        nwantedvox % searchlight size (voxels)
        lastspheren % store the last n necessary to achieve masked nvox
        mapcoords % handle to either mapr or mapn methods
        nvox = []; % number of voxels in actual searchlight (after mask)
        vol % MriVolume or sub-class for mask
        distances % mm distances for spheres with <min(vol.header.dim) radius
        xyz % voxel coordinates (centered on 0) for distances
        rmax % largest possible radius (limited by volume size)
        nmax % largest possible nvox
    end

    methods
        function sl = Searchlight(maskpath,mapmode,radvox)
        % sl = Searchlight(maskpath,mapmode,radvox)
            if ischar(maskpath) 
                sl.vol = SPMVolume([],maskpath);
            else
                sl.vol = MriVolume([],maskpath);
            end
            % decide on how we are going to map
            sl.mapmode = lower(mapmode);
            switch sl.mapmode
                case 'radius'
                    sl.radius = radvox;
                    sl.mapcoords = @sl.mapr;
                case 'nvox'
                    sl.nwantedvox = radvox;
                    % initialise with nwanted
                    sl.lastspheren = sl.nwantedvox;
                    sl.mapcoords = @sl.mapn;
                otherwise
                    error('unknown mapmode: %s',mapmode)
            end
            % find coordinate / distance mapping 
            % we assume that you'd never want a radius that exceeds the
            % length of the shortest dimension in your volume
            dims = repmat(min(sl.vol.header.dim)...
              +isodd(min(sl.vol.header.dim)),[1 3]);
            % build xyz matrices for estimating distance from 0
            [x,y,z] = meshgrid(-dims(1):dims(1),-dims(2):dims(2),...
                -dims(3):dims(3));
            % euclidean distance of each dimension (scaled to mm by
            % voxsize) from 0
            distancemat = sqrt((x*sl.vol.voxsize(1)).^2 + ...
                (y*sl.vol.voxsize(2)).^2 + ...
                (z*sl.vol.voxsize(3)).^2);
            % linear indices for distancemat
            dist_lininds = find(distancemat);
            % vector of distances
            sl.distances = distancemat(dist_lininds);
            % % extract voxel coordinates
            [coords_x,coords_y,coords_z] = ind2sub(size(distancemat),...
                dist_lininds);
            % make 3 by n matrix and shift to center on 0
            % (we then add current centre coordinates to get to the right
            % place for a given searchlight sphere)
            sl.xyz = (([coords_x coords_y coords_z]) - repmat(dims+1,...
                [length(sl.distances) 1]))';
            % sort by distance - can now index the first n sl.coords to get
            % a searchlight of nvox size (may have to go into a while loop
            % to find nvox that are also inside mask), or for radius-based
            % mapping, index the sl.xyz where sl.distances <= r
            [sl.distances,inds] = sort(sl.distances);
            sl.xyz = sl.xyz(:,inds);
            % insert central voxel
            % nb, because searchlight maps often end up being sparse we
            % insert a very small value rather than 0 to avoid problems
            % later.
            sl.distances = [realmin('single'); sl.distances];
            sl.xyz = [[0 0 0]' sl.xyz];
            % this approach only supports searchlights up to a certain
            % (ludicrous) size
            sl.rmax = max(sl.distances);
            sl.nmax = min([sum(sl.distances < sl.rmax) sl.vol.nfeatures]);
        end

        function out = mapr(self,xyz)
        % out = mapr(self,xyz)
        % Create a searchlight sphere for a given xyz coordinate (voxel
        % index) based on self.radius. Returns column indices into a mapped
        % dataset 
            assert(self.vol.mask(xyz(1),xyz(2),xyz(3))==1,...
                'coordinates %s are outside mask',mat2str(xyz));
            assert(self.radius < self.rmax,...
                'current radius (%f) exceeds rmax (%f)',self.radius,...
                self.rmax);
            % get coordinates inside radius
            maskxyz = self.xyz(:,self.distances <= self.radius);
            % pass off to general method
            out = self.mapsphere(xyz,maskxyz);
        end

        function out = mapn(self,xyz)
        % out = mapn(self,xyz)
            assert(self.vol.mask(xyz(1),xyz(2),xyz(3))==1,...
                'coordinates %d are outside mask',xyz);
            assert(self.nwantedvox < self.nmax,...
                'current nwanted (%d) exceeds nmax (%d)',...
                self.nwantedvox,self.nmax);
            % need to iteratively scale searchlight to find a size that
            % achieves the desired n
            currentn = self.nwantedvox; %self.lastspheren;
            done = 0;
            step = 100;
            niter = 0;
            while ~done
                % get first n coordinates
                coords = self.xyz(:,1:currentn);
                % this search method is accurate but slow
                out = self.mapsphere(xyz,coords);
                if self.nvox < self.nwantedvox
                    currentn = currentn+step;
                elseif self.nvox > self.nwantedvox
                    % overshot the mark. Step back, reduce step size and
                    % try again
                    currentn = currentn-step;
                    assert(step~=1,'should never happen')
                    step = ceil(step/2);
                else
                    % must be equal
                    done = 1;
                    % store the radius / sphere n necessary to achieve this
                    self.radius = self.distances(currentn);
                    self.lastspheren = currentn;
                end
                niter = niter+1;
                assert(niter<10e3,'iteration limit exceeded');
            end
        end

        function out = mapsphere(self,xyz,coords)
        % out = mapsphere(xyz,coords)
        % Find the volume feature indices centered on xyz corresponding to
        % the entered _searchlight_ coords. Used internally by mapr/mapn.
        % You should probably call these directly instead.
            spherenvox = size(coords,2);
            % shift to current position
            coords = coords + repmat(xyz,[1 spherenvox]);
            % find coordinates outside volume limits
            insidev = all(coords>=1,1) & all(coords<=repmat(...
                self.vol.header.dim',[1 spherenvox]),1);
            coords_insidev = coords(:,insidev);
            % find coordinates inside mask
            % first back to linear indices
            lininds_insidev = self.vol.coord2linind(coords_insidev);
            % then restrict to mask lininds
            lininds_out = intersect(lininds_insidev,self.vol.linind);
            % final estimate of how many voxels we ended up with
            self.nvox = length(lininds_out);
            % convert linear indices to feature indices
            out = self.vol.linind2featind(lininds_out);
        end

        function out = mapinds(self,ind)
        % out = mapinds(self,n)
        % Get the nth searchlight sphere.
            out = self.mapcoords(self.vol.linind2coord(...
                self.vol.linind(ind)));
        end
    end
end
