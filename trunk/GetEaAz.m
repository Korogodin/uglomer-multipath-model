function up = GetEaAz(user, st)

  llh = ConvXYZ(user);
  north.X = -cos(llh.Y)*sin(llh.X);
  north.Y = -sin(llh.Y)*sin(llh.X);
  north.Z = cos(llh.X);
  east.X = -sin(llh.Y);
  east.Y = cos(llh.Y);
  east.Z = 0.0;
  up.X = cos(llh.Y)*cos(llh.X);
  up.Y = sin(llh.Y)*cos(llh.X);
  up.Z = sin(llh.X);

  xls = st.X-user.X;
  yls = st.Y-user.Y;
  zls = st.Z-user.Z;
  range1 = sqrt(xls*xls+yls*yls+zls*zls);

	tdot = (up.X*xls+up.Y*yls+up.Z*zls)/range1;
	xls = xls/range1;
	yls = yls/range1;
	zls = zls/range1;

	if tdot >= 1.00 
        satang = pi/2.0;
    elseif tdot <= -1.00 
        satang = -pi/2.0; 
    else
        satang = asin(tdot);
    end

	xaz = east.X*xls+east.Y*yls;
	yaz = north.X*xls+north.Y*yls+north.Z*zls;
	
    if ((xaz ~= 0.0) || (yaz ~= 0.0)) 
        az2 = atan2(xaz,yaz); 
    else
        az2 = 0.0;
    end
    
    if (az2 < 0) 
        az2 = az2 + 2*pi;
    end
	up.X = satang; %// угол места
	up.Y  = az2;   %// азимут
%     figure(1); plot(user.X, user.Y, '*', st.X, st.Y, '*'); grid on;
%     figure(2); polar(up.Y, pi/2 - up.X, '*');
%     figure(3); plot(llh.Y, llh.X, '*'); grid on;
    
function resa = ConvXYZ(coord)

    Ra = 6378137.0;
	f = 1.0/298.257223563;     %//reciprocal flattening
	b = Ra*(1.0-f);            % //semi-minor axis
	e2 = 2.0*f-f*f;            %//first eccentricity squared
	ep2 = f*(2.0-f)/((1.0-f)*(1.0-f));% //second eccentricity squared

    if ((coord.X == 0)&&(coord.Y == 0)&&(coord.Z == 0))
        resa.X = 0; resa.Y = 0; resa.Z = -Ra;
    else
	r2 = coord.X*coord.X + coord.Y*coord.Y;
	r = sqrt(r2);
	E3 = Ra*Ra - b*b;
	F1 = 54.0*b*b*coord.Z*coord.Z;
	G = r2 + (1.0-e2)*coord.Z*coord.Z - e2*E3;
	c1 = (e2*e2*F1*r2)/(G*G*G);
	s = ( 1.0 + c1 + sqrt(c1*c1 + 2.0*c1) )^(1.0/3.0);
	P = F1/(3.0*(s+1.0/s+1.0)*(s+1.0/s+1.0)*G*G);
	Q = sqrt(1.0+2.0*e2*e2*P);
	ro = -(e2*P*r)/(1.0+Q) + sqrt((Ra*Ra/2.0)*(1.0+1.0/Q) - ...
		((1.0-e2)*P*coord.Z*coord.Z)/(Q*(1.0+Q)) - P*r2/2.0);
	tmp = (r - e2*ro)*(r - e2*ro);
	U = sqrt( tmp + coord.Z*coord.Z );
	V = sqrt( tmp + (1.0-e2)*coord.Z*coord.Z );
	zo = (b*b*coord.Z)/(Ra*V);
	h = U*( 1.0 - b*b/(Ra*V));
	phi = atan( (coord.Z + ep2*zo)/r );
	lambda = atan2(coord.Y,coord.X);

	resa.X = phi;    %// широта latitude
	resa.Y = lambda; %// долгота longitude
	resa.Z = h;     %// высота height
    end
