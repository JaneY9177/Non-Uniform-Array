%% Configuration
f0=24e9;
wavelength=1;
beta=2*pi/wavelength;
theta=0:0.05:360;
mainbeam=5;
n=9;
d=[0.25 + (2.5-0.25).*rand(1),0.5 + (5-0.5).*rand(1,7)];

%% 
%A=[2*cos(beta*d0*cosd(theta)) 2*cos(beta*d1*cosd(theta)) 2*cos(beta*d2*cosd(theta)) 2*cos(beta*d3*cosd(theta))];
A=[2*cos(beta*d(1)*cosd(theta))
    2*cos(beta*(d(2)+d(1))*cosd(theta))
    2*cos(beta*(d(3)+d(2)+d(1))*cosd(theta))
    2*cos(beta*(d(4)+d(3)+d(2)+d(1))*cosd(theta))
    2*cos(beta*(d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta))
    2*cos(beta*(d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta))
    2*cos(beta*(d(7)+d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta))
    2*cos(beta*(d(8)+d(7)+d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta))
    ]';

theta1=0:0.05:90-mainbeam/2;
theta2=90+mainbeam/2:0.05:180;

A1=[2*cos(beta*d(1)*cosd(theta1))
    2*cos(beta*(d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(7)+d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(8)+d(7)+d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    ]';

b=[-1*ones(1,length(theta1))
    2*cos(beta*d(1)*cosd(theta1))
    2*cos(beta*(d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(7)+d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    2*cos(beta*(d(8)+d(7)+d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    ]';
% 
c=[-1*ones(1,length(theta1))
    -2*cos(beta*d(1)*cosd(theta1))
    -2*cos(beta*(d(2)+d(1))*cosd(theta1))
    -2*cos(beta*(d(3)+d(2)+d(1))*cosd(theta1))
    -2*cos(beta*(d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    -2*cos(beta*(d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    -2*cos(beta*(d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    -2*cos(beta*(d(7)+d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    -2*cos(beta*(d(8)+d(7)+d(6)+d(5)+d(4)+d(3)+d(2)+d(1))*cosd(theta1))
    ]';

cvx_begin
    variable x(n)
    minimize( abs(x(1)) )
    subject to
        [0 2 2 2 2 2 2 2 2]*x==1;
        b*x<=0;
        c*x<=0;
cvx_end

% n=8;
% lambda = 1;           % wavelength
% theta_tar = 90;       % target direction
% min_sidelobe = -10;   % maximum sidelobe level in dB
% 
% max_half_beam = 50;   % starting half beamwidth (must be feasible)
% halfbeam_bot = 1;
% halfbeam_top = max_half_beam;
% while( halfbeam_top - halfbeam_bot > 1)
%   % try to find a feasible design for given specs
%   halfbeam_cur = ceil( (halfbeam_top + halfbeam_bot)/2 );
% 
%   % create optimization matrices for the stopband
%   ind = find(theta <= (theta_tar-halfbeam_cur) | ...
%              theta >= (theta_tar+halfbeam_cur) );
% 
%   % formulate and solve the feasibility antenna array problem
%   cvx_begin quiet
%     variable w(n)
%     % feasibility problem
%     [2 2 2 2 2 2 2 2]*w == 1;
%     abs(A1*w) <= 10^(min_sidelobe/20);
%   cvx_end
% 
%   % bisection
%   if strfind(cvx_status,'Solved') % feasible
%     fprintf(1,'Problem is feasible for half beam-width = %d degress\n',...
%                halfbeam_cur);
%     halfbeam_top = halfbeam_cur;
%   else % not feasible
%     fprintf(1,'Problem is not feasible for half beam-width = %d degress\n',...
%                halfbeam_cur);
%     halfbeam_bot = halfbeam_cur;
%   end
% end

plot(theta,20*log10(abs(A*x(2:9))));
axis([0,180,-30,0]);
