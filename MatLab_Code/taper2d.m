function gt=taper2d(g,npts,ny,nx,ydiff,xdiff)
gt=zeros(npts); gt(ydiff+1:ydiff+ny,xdiff+1:xdiff+nx)=g;
gp=g(:,1:3);     [gpx1,~]=gradient(gp);  % sides
gp=g(:,nx-2:nx); [gpx2,~]=gradient(gp);
x1=0; x2=(2*xdiff)+1;
x=[1 1 0 0;x1 x2 1 1; x1^2 x2^2 2*x1 2*x2; x1^3 x2^3 3*x1^2 3*x2^2];
for I=1:ny
    y=[g(I,nx) g(I,1) gpx2(I,3) gpx1(I,1)];
    c=y/x;
    for J=1:xdiff
        gt(I+ydiff,J)=c(1)+(J+xdiff)*c(2)+c(3)*(J+xdiff)^2+c(4)*(J+xdiff)^3;
        gt(I+ydiff,J+nx+xdiff)=c(1)+J*c(2)+c(3)*J^2+c(4)*J^3;
    end
end
gp=g(1:3,:);     [~,gpx1]=gradient(gp);  % top and bottom
gp=g(ny-2:ny,:); [~,gpx2]=gradient(gp);
x1=0; x2=(2*ydiff)+1;
x=[1 1 0 0;x1 x2 1 1; x1^2 x2^2 2*x1 2*x2; x1^3 x2^3 3*x1^2 3*x2^2];
for J=1:nx
    y=[g(ny,J) g(1,J) gpx2(3,J) gpx1(1,J)];
    c=y/x;
    for I=1:ydiff
        gt(I,J+xdiff)=c(1)+(I+ydiff)*c(2)+c(3)*(I+ydiff)^2+c(4)*(I+ydiff)^3;
        gt(I+ydiff+ny,J+xdiff)=c(1)+I*c(2)+c(3)*I^2+c(4)*I^3;
    end
end
for I=ydiff+ny+1:npts
    for J=xdiff+nx+1:npts
        if (I-ny-ydiff)>(J-nx-xdiff)
            gt(I,J)=gt(I,nx+xdiff);
        else
            gt(I,J)=gt(ny+ydiff,J);
        end
    end
end
for I=1:ydiff
    for J=1:xdiff
        if I>J
            gt(I,J)=gt(ydiff+1,J);
        else
            gt(I,J)=gt(I,xdiff+1);
        end
    end
end
for I=1:ydiff
    for J=xdiff+nx+1:npts
        if I>(npts-J)
            gt(I,J)=gt(ydiff+1,J);
        else
            gt(I,J)=gt(I,xdiff+nx);
        end
    end
end
for I=ydiff+ny+1:npts
    for J=1:xdiff
        if (npts-I)>J
            gt(I,J)=gt(ydiff+ny,J);
        else
            gt(I,J)=gt(I,xdiff+1);
        end
    end
end