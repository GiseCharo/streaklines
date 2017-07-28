function Hstreak=interpstreak_any(stline,H,X,Y)
%%interpolates the H along the streakline
%

% Reshape
s = size(H);
N = prod(s(3:end));
H = reshape(H,[s(1:2) N]);

% Initialisation
Hstreak=zeros([size(stline,1),N]);
for k=1:N % Loop along the streakline
    Hstreak(:,k) =interp2(X(:,1),Y(1,:),H(:,:,k)',stline(:,1),stline(:,2));
end

% Reshape
if length(s)<3
    s(3)=1;
end
Hstreak = reshape(Hstreak,[size(stline,1) s(3:end)]);