if ( 0 )
scale = 0.2;
wd = (pi^2) / (8 * scale^2); wa = 1;

x = [ 0 : 0.001 : scale ];
y = [ 0 : pi/length(x) : pi ];
[xx,yy] = meshgrid(x,abs(min(y,repmat(pi,1,size(y,2))-y)));
z = wd * xx.^2 + wa * yy.^2;
surf(x,y,z); 
xlabel('dist'), ylabel('angle');
end

%  ellipse
figure();

range_mult = 3;
ang_thresh = 0.3;
wd = 0.33;

x = [ 0 : (scale * range_mult)/100 : scale*range_mult ];
y = [ 0 : ang_thresh*range_mult/length(x) : ang_thresh*range_mult ];
[xx,yy] = meshgrid(x,y);
z = ((wd * xx).^2) / (scale*scale) + (yy.^2) / (ang_thresh*ang_thresh);
z( z > 1 ) = 2;
mesh(x,y,z); 
xlabel('dist'), ylabel('angle');
%set(gca,'YTickLabel',y)
return;


% hyperbole
figure();

scale = 0.05;
range_mult = 2;
ang_thresh = 0.3;
wd = sqrt(scale);

x = [ 0 : (scale * range_mult)/100 : scale*range_mult ];
y = [ 0 : ang_thresh*range_mult/length(x) : ang_thresh*range_mult ];
[xx,yy] = meshgrid(x,y);
z = ((2*ang_thresh-yy).^2) / (ang_thresh*ang_thresh) - ((xx).^2) / (scale*scale);
z( z <= 1 ) = -2;
mesh(x,y,z); 
xlabel('dist'), ylabel('angle');
%set(gca,'YTickLabel',y*180/pi)
return;