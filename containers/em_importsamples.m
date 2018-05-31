% Import samples for data analysis. Optionally enter a parameters file to
% produce consistent trigger coding across datasets - also checks that
% other parameters are consistent across datasets.
% [samples,parameters] = em_parsesamplestxt(samplests,[parameters])
function [samples, par] = em_importsamples(samplestxt,expar)

% Parse the sample text file
% [samples, par] = em_parsesamplestxt(samplestxt)
fid = fopen(samplestxt);
% Read header info
nhlines = 21;
theader = textscan(fid,'%s',nhlines,'delimiter','\n','whitespace','');
theader = theader{1};
% Read data
tdata = textscan(fid,'%s','delimiter','\n','whitespace','');
tdata = tdata{1};
fclose(fid);
% Now to make sense of things...

% Extract parameters from this session
datetime = textscan(theader{3}, '%*s %*s %s %s');
par.rec_date = datetime{1}{1};
par.rec_time = datetime{2}{1};
txtpath = textscan(theader{2},'%*s %*s %*s %s');
par.txt_path = txtpath{1}{1};
samplingr = textscan(theader{5},'%*s %*s %*s %f');
par.sampling_rate_hz = samplingr{1};
filename = textscan(theader{7},'%*s %*s %s');
par.filename = filename{1}{1};
screenres = textscan(theader{11},'%*s %*s %*s %f %f');
par.screen_res_px = cell2mat(screenres);
screenwidth = textscan(theader{13},'%*s %*s %*s %*s %f %*d');
par.screenwidth = screenwidth{1};
viewdistance = textscan(theader{14},'%*s %*s %*s %*s %f');
par.viewdistance = viewdistance{1};
% Detect binocular acquisitions
par.binocular = ~isempty(strfind(theader{19},'RIGHT'));

% Marker code start number
mc = 10e3;

% No previous markers to consider, of checks to make
if ~exist('expar','var') || ~isfield(expar,'markers')
	% Initialise empty markers struct
	expar.markers = struct;
end

% Make some basic consistency checks
for fn = {'sampling_rate_hz','screenwidth','viewdistance','binocular'}
	fn = fn{1};
	if isfield(expar,fn)
		if ~par.(fn) == expar.(fn)
			error('Inconsistent parameters for %s %s', ...
				par.filename,expar.filename)
		end
	end
end

% Figure out how many markers we've already used
codes = cell2mat(structfun(@(x){x},expar.markers));
% And ensure we don't duplicate
mc = max([max(codes) mc]);

% Now update par with expar
par = catstruct(par,expar);

% NB These fields are not actually used now - just reserved for later use
% in event extraction
if ~isfield(par.markers,'Saccade');
	mc = mc + 1e3;
	par.markers.Saccade = mc;
	mc = mc + 1e3;
	par.markers.Blink = mc;
	mc = mc + 1e3;
	par.markers.Fixation = mc;
end

% Initialise output matrix
nsamples = length(tdata);
% row: [timestamp eventflag Lx Ly [Rx Ry]]

% textscan setup for messages
msgstr = '%f %*s %*d %*s %*s %s';
if par.binocular
	samples = zeros(nsamples,6);
	% textscan for binocular data
	scanstr = '%f %*s %*d %f %f %f %f';
	linewidth = 6;
else
	samples = zeros(nsamples,4);
	% textscan for monocular data
	scanstr = '%f %*s %*d %f %f';
	linewidth = 4;
end
samples(:) = NaN;

% Now parse samples
for s = 1:length(tdata);
	% Behaviour depends on line type
	switch isempty(strfind(tdata{s},'MSG'));
		case 1 % regular data line (SMP)
			% 3 or 5 entries
			lineout = cell2mat(textscan(tdata{s},scanstr));
			% Column indices to insert data into
			lineIs = [1 3:linewidth];
		case 0 % message line
			linedata = textscan(tdata{s},msgstr);
			% Strip any extensions
			[x,fn,x] = fileparts(linedata{2}{1});
			if ~isfield(par.markers,fn)
				mc = mc + 1e3;
				par.markers.(fn) = mc;
			end
			lineout = [linedata{1} par.markers.(fn)];
			lineIs = [1 2];
	end
	% Add data to outmat
	samples(s,lineIs) = lineout;
end

% Detect sorting errors in sample log
if ~issorted(samples(:,1))
	fprintf('Repaired sample sort error in %s', ...
		par.filename)
	[x,I] = sort(samples(:,1));
	samples = samples(I,:);
end

% Rescore time to ms
samples(:,1) = samples(:,1) * .001;
