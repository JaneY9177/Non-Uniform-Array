clear;
%% Configuration
wavelength=1;
beta=2*pi/wavelength;   % wave number
theta=0:0.05:360;
elementNumber=16;
spacingMIN=0.5;
spacingMAX=5;

count=1;

%%
for m=1:10000
    d=[spacingMIN/2 + (spacingMAX/2-spacingMIN/2).*rand(1), spacingMIN + (spacingMAX-spacingMIN).*rand(1,7)];

    mainbeam=5;

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
    
    n=elementNumber/2+1;


    cvx_begin
        variable x(n)
            minimize( abs(x(1)) )
        subject to
            [0 2 2 2 2 2 2 2 2]*x==1;
            b*x<=0;
            c*x<=0;
    cvx_end
    
    if x(1)<10^(-10/20)
        resultd(count,:)=d;
        resultx(count,:)=x';
        count=count+1;
        %plot(theta,20*log10(abs(A*x(2:9))));
        %axis([0,180,-70,0]);
        %hold on;
    end
end
        
% plot(theta,20*log10(abs(A*x(2:9)))-max(20*log10(abs(A*x(2:9)))));
% axis([0,180,-70,0]);
