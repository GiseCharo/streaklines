function norm2v=norm2(v)
%%compute the (norm v)^2
dimv=size(v,1)*size(v,2);
if dimv>2
    norm2v=v(:,:,1).*v(:,:,1)+v(:,:,2).*v(:,:,2);
else
    norm2v=v(1).*v(1)+v(2).*v(2);
end
end