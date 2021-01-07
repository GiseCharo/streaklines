function v=projection_streak(f,s)
%%proyection along s=dtophi/norm(dtophi)
%%size(f)=(lengthstreakline,2)
%%size(s)=(lengthstreakline,2)
%f N 1 d
f=permute(f,[1 3 2]); %N d 1 
v=f.*s;
v=sum(v,2);

end