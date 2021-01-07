function fperp=ortogonal(f)
fperp=zeros(length(f),2);
fperp(:,1)=f(:,2);
fperp(:,2)=-f(:,1);
end