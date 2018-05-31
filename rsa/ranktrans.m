% fast, vectorised rank transform of the rows in data input x. 
%
% adapted from ranktrans_ function on file exchange:
% http://www.mathworks.com/matlabcentral/fileexchange/34560-ranktrans--x-dim-
%
% r = ranktrans(x)
function r = ranktrans(x)

siz = size(x);
wasnan = isnan(x);
% rescore nans to a number for speed (each individual nan gets a unique
% rank otherwise)
x(wasnan) = max(x(:)) + 1;

%[Step 1]: sort and get sorting indicies
[x,Ind] = sort(x,1);

%[Step 2]: create matrix [D], which has +1 at the start of
%consecutive runs and -1 at the end, with zeros elsewhere.
D = zeros(siz,'int8');
D(2:end-1,:) = diff(x(1:end-1,:) == x(2:end,:));
D(1,:) = x(1,:) == x(2,:);
D(end,:) = -( x(end,:) == x(end-1,:) );

%[Step 3]: calculate the averaged rank for each consecutive
%run
[a,~] = find(D);
a = reshape(a,2,[]);
h = sum(a,1)/2;

%[Step 4]: insert the troublseome ranks in the relevant
% places
L = zeros(siz);
L(D==1) = h;
L(D==-1) = -h;
L = cumsum(L);
L(D==-1) = h; %cumsum set these ranks to zero, but we
% wanted them to be h

%[Step 5]: insert the simple ranks (i.e. the ones that
% didn't clash)
[L(~L),~] = find(~L);

%[Step 6]: assign the ranks to the relevant position in
% the matrix
r = NaN(siz,class(x));
r(bsxfun(@plus,Ind,(0:siz(2)-1)*siz(1))) = L;
% reintroduce NaNs
r(wasnan) = NaN;
