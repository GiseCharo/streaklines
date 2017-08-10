function df=high_order_diff_along_streak(f,N,dt0)
% Compute the derivative of f along the first dimension
% with a N-point scheme (N=<7 and N odd)
% The result will be of size size(f)
%

P = (N-1)/2;
if N < 3 || N > 7 || floor(P)~=P
    error('incorrect order in the high-order finite difference scheme');
end

d=ndims(f);
Mt0 = size(f,1);

idx='';
for k_dim=2:d
    idx = [idx ',:'];
end
df = zeros(size(f)); % m x d2 x Mt0 x MY

MY= size(f);
MY= MY(2:end);

if (Mt0 < 12 && floor((Mt0+1)/2)==(Mt0+1)/2) || Mt0 < 2
    error(['There are not enough points along the streakline in order ' ...
        'to compute the Laplacien.' ]);
end


switch N
    case 7
        % Filter for interrior
        filter = [ -1/60  3/20 -3/4 0 3/4 -3/20 1/60 ]';
        % Filters for the left boundary
        filterleft_bound = [ -49/20 6 -15/2 20/3 -15/4 6/5 -1/6 ]';
    case 5
        % Filter for interrior
        filter = [ 1/12  -2/3 0 2/3 -1/12 ]';
        % Filters for the left boundary
        filterleft_bound = [ -25/12 4 -3 4/3 -1/4 ]';
    case 3
        % Filter for interrior
        filter = [ -1/2 0 1/2 ]';
        % Filters for the left boundary
        filterleft_bound = [-3/2 2 -1/2 ]';
    otherwise
        error('This number of points is not encoded.')
end
% Filters for the right boundary
filterright_bound = - filterleft_bound(end:-1:1) ;

% Boundaries

% From 1 to P pixels far from the left boundary
for k=1:min(P,Mt0/2)
    % Get neighborhoods
    len_filter = N;
    %     len_filter = 2*(k-1)+1;
    neighb1 = zeros([ MY len_filter]);
    eval(['neighb1(:' idx ')= permute( f(k:(k-1+len_filter)' idx ') ,' ...
        '[ 2:ndims(f) 1]);']);
    filter1=permute(filterleft_bound,[2:ndims(neighb1) 1 ]); % 1 x (1 x)*length(MY) x len_filter
    da1 = bsxfun(@times,neighb1, filter1); % MY x len_filter
    clear neighb1;
    da1= sum(da1,ndims(da1)); % MY
    eval(['df(k' idx ')= 1/dt0 * da1;']);
end

% From 1 to P pixels far from the right boundary
for k=1:min(P,Mt0/2)
    % Get neighborhoods
    len_filter = N;
    %     len_filter = 2*(k-1)+1;
    neighb1 = zeros([ MY len_filter]);
    eval(['neighb1(:' idx ')= permute( f((1-k+end-len_filter+1):(1-k+end)' idx ') ,' ...
        '[ 2:ndims(f) 1]);']);
    filter1=permute(filterright_bound,[2:ndims(neighb1) 1 ]); % 1 x (1 x)*length(MY) x len_filter
    da1 = bsxfun(@times,neighb1, filter1); % 2 x MY x len_filter
    clear neighb1;
    da1= sum(da1,ndims(da1)); % 2 x MY
    eval(['df( Mt0-(k-1) ' idx ')= 1/dt0 * da1;']);
end

if Mt0 >= N
    % Interior
    % Get neighborhoods
    len_filter = N;
    neighb1 = get_neighborhood_mat_1D( f, len_filter,idx); % Mt0 x MY x len_filter
    eval(['neighb1 = neighb1( (P+1):(end-P)' idx ',:);']); % (Mt0-12) x MY x len_filter
    filter1=permute(filter,[2:ndims(neighb1) 1 ]); % 1 x (1 x)*length(MY) x len_filter
    da1 = bsxfun(@times,neighb1, filter1); % (Mt0-12) x MY x len_filter
    clear neighb1;
    da1= sum(da1,ndims(da1));  % (Mt0-12) x MY
    eval(['df((P+1):(end-P)' idx ')= 1/dt0 * da1 ;']);
end

end


function neighb = get_neighborhood_mat_1D(A,len,idx)
% Get the len^2-1 neighborhood of each points
% A must have the size :  Mt0 x MY (2D!!!)
%

s=size(A);
neighb = zeros([s len]);% Mt0 x MY x len

le=ceil((len-1)/2);

for i=-le:0
    eval(['neighb(1:end+i' idx ',le+1+i)=A(1-i:end' idx ');']);
end
for i =1:le
    eval(['neighb((1+i):end' idx ',le+1+i)=A(1:end-i' idx ');']);
end

eval(['neighb=neighb(:' idx ',end:-1:1);']);

end
