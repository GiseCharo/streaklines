function df=high_order_diff_along_streak(f,N,dt0)
% Compute the derivative of f along the first dimension 
% at the order N (<12)
% The result will be of size size(f)
%

P = (N-1)/2;
if N > 13 || floor(P)~=P
    error('incorrect order in the high-order finite difference scheme');
end

type = 'standard';
%type = 'optimized';

d=ndims(f);
% d=ndims(f)-2;
% d2=size(f,2);
Mt0 = size(f,1);

idx='';
for k_dim=2:d
    idx = [idx ',:'];
end
%
% f = permute(f,[(1:ndims(f)-d) ndims(f)-d+1+mod((l-1):(l-2+d),d)]); % Put the k coordinate in first place
% % m x d2 x Mt0 x MY
df = zeros(size(f)); % m x d2 x Mt0 x MY

% if d==2
%     [ m, ~, Mt0, My] = size(f);
%     MY=My;
% else
%     [ m, ~, Mt0, My, Mz] = size(f);
%     MY=[My Mz];
% end
MY= size(f);
MY= MY(2:end);
% dt0=dX(l);

if (Mt0 < 12 && floor((Mt0+1)/2)==(Mt0+1)/2) || Mt0 < 2
    error(['There are not enough points along the streakline in order ' ...
        'to compute the Laplacien.' ]);
end


if strcmp(type,'standard')
    % Filter for interrior
    FD13p = [ 1/5544 -3/1155 1/56 -5/63 15/56 -6/7 0 ...
        6/7 -15/56 5/63 -1/56 3/1155 -1/5544 ]';
    % Filters for the boundaries
    FD11p = [ -1/1260 5/504 -5/84 5/21 -5/6 0 5/6 -5/21 5/84 -5/504 1/1260 ]';
    FD9p = [ 1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280 ]';
    FD7p = [ -1/60  3/20 -3/4 0 3/4 -3/20 1/60 ]';
    FD5p = [ 1/12  -2/3 0 2/3 -1/12 ]';
    FD3p = [ -1/2 0 1/2 ]';
    FD2p = [ -1 1]';
elseif strcmp(type,'optimized')
    % Filter for interrior
    FD13p = [ +0.001456501759 ...
        -0.011169294114 ...
        +0.045246480208 ...
        -0.133442885327 ...
        +0.337048393268 ...
        -0.907646591371 ...
        0              ...
        +0.907646591371 ...
        -0.337048393268 ...
        +0.133442885327 ...
        -0.045246480208 ...
        +0.011169294114 ...
        -0.001456501759 ]';
    
    % Filters for the boundaries
    FD11p = [ -0.002484594688 ...
        +0.020779405824 ...
        -0.090320001280 ...
        +0.286511173973 ...
        -0.872756993962 ...
        0              ...
        +0.872756993962 ...
        -0.286511173973 ...
        +0.090320001280 ...
        -0.020779405824 ...
        +0.002484594688 ]';
    
    FD9p = [  +0.007650904064 ...
        -0.059463584768 ...
        +0.244678631765 ...
        -0.841570125482 ...
        0               ...
        +0.841570125482 ...
        -0.244678631765 ...
        +0.059463584768 ...
        -0.007650904064 ]';
    
    FD7p = [ -1/60  3/20 -3/4 0 3/4 -3/20 1/60 ]';
    FD5p = [ 1/12  -2/3 0 2/3 -1/12 ]';
    FD3p = [ -1/2 0 1/2 ]';
    FD2p = [ -1 1]';
end
filters={FD2p FD3p FD5p FD7p FD9p FD11p FD13p};

% Boundaries

% Boundaries pixels
neighb1 = zeros([ 2 MY 2]);

eval(['neighb1(1,:' idx ')= permute( f([1 2]' idx ') ,' ...
    '[ ndims(f)+1 2:ndims(f) 1]);']);
eval(['neighb1(2,:' idx ')= permute( f([end-1 end]' idx ') ,' ...
    '[ ndims(f)+1 2:ndims(f) 1]);']);
filter1=permute(filters{1},[2:ndims(neighb1) 1 ]); % 1 x (1 x)*length(MY) x len_filter
da1 = bsxfun(@times,neighb1, filter1); % 2 x MY x len_filter
clear neighb1;
da1= sum(da1,ndims(da1)); % 2 x MY
eval(['df([1 end]' idx ')= 1/dt0 * da1;']);

% From 2 to P pixels far from the boundaries
for k=2:min(P,Mt0/2)
    % Get neighborhoods
    len_filter = 2*(k-1)+1;
    neighb1 = zeros([ 2 MY len_filter]);
    eval(['neighb1(1,:' idx ')= permute( f(1:len_filter' idx ') ,' ...
        '[ ndims(f)+1 2:ndims(f) 1]);']);
    eval(['neighb1(2,:' idx ')= permute( f((end-len_filter+1):end' idx ') ,' ...
        '[ ndims(f)+1 2:ndims(f) 1]);']);
    filter1=permute(filters{k},[2:ndims(neighb1) 1 ]); % 1 x (1 x)*length(MY) x len_filter
    da1 = bsxfun(@times,neighb1, filter1); % 2 x MY x len_filter
    clear neighb1;
    da1= sum(da1,ndims(da1)); % 2 x MY
    eval(['df([k (Mt0-(k-1)) ]' idx ')= 1/dt0 * da1;']);
end

if Mt0 >= N
    % Interior
    % Get neighborhoods
    len_filter = N;
    neighb1 = get_neighborhood_mat_1D( f, len_filter,idx); % Mt0 x MY x len_filter
    eval(['neighb1 = neighb1( (P+1):(end-P)' idx ',:);']); % (Mt0-12) x MY x len_filter
    filter1=permute(filters{P+1},[2:ndims(neighb1) 1 ]); % 1 x (1 x)*length(MY) x len_filter
    da1 = bsxfun(@times,neighb1, filter1); % (Mt0-12) x MY x len_filter
    clear neighb1;
    da1= sum(da1,ndims(da1));  % (Mt0-12) x MY
    eval(['df((P+1):(end-P)' idx ')= 1/dt0 * da1 ;']);
end

% df = permute(df,[(1:ndims(df)-d) ndims(df)-d+1+mod((1-l):(d-l),d)]); % Put the first coordinate in first place again

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
