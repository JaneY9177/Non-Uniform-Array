clear;
%% Configuration
wavelength=1;
k=2*pi/wavelength;   % wave number
angleStep=0.05;
theta=0:angleStep:360;
elementNumber=16;
spacingMIN=0.5;
spacingMAX=3;
mainbeam1=8;
% mainbeam2=20;
thetaM1=90;
% thetaM2=45;

count=1;

%%
for m=1:10000
    d=[spacingMIN/2 + (spacingMAX/2-spacingMIN/2).*rand(1), spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber/2-1)];
    d1=fliplr(d);
    d=[0,d1(1:length(d1)-1),d1(length(d1))+d(1),d(2:length(d))];
    %d=[0.25,0.5*ones(1,elementNumber/2-1)];
    for nn=2:length(d)
        d(nn)=d(nn-1)+d(nn);
    end
    %theta=90;
    
    A=zeros(length(theta),elementNumber);
    for nn=1:length(d)
        A(:,nn)=exp(1i*k*d(nn)*cosd(theta));
    end
    
    A_M=zeros(1,elementNumber);
    for nn=1:length(d)
        A_M(:,nn)=exp(1i*k*d(nn)*cosd(thetaM1));
    end
    %w=ones(16,1);
    
    theta_SL1=[0:angleStep:thetaM1-mainbeam1/2,thetaM1+mainbeam1/2:angleStep:180];
%     theta_SL2=[0:angleStep:thetaM2-mainbeam2/2,thetaM2+mainbeam2/2:angleStep:180];
    %theta_ML=90+mainbeam/2:angleStep:180;
    
    A_SL1=zeros(length(theta_SL1),elementNumber);
    for nn=1:length(d)
        A_SL1(:,nn)=exp(1i*k*d(nn)*cosd(theta_SL1));
    end
    
%     A_SL2=zeros(length(theta_SL2),elementNumber);
%     for nn=1:length(d)
%         A_SL2(:,nn)=exp(1i*k*d(nn)*cosd(theta_SL2));
%     end
    
    cvx_begin
    variable w(elementNumber) complex
    minimize( max(abs(A_SL1*w)) )
    subject to
    A_M*w==1;
    cvx_end
    
    if max(abs(A_SL1*w))<10^(-10/20)
        resultd(count,:)=d;
        resultw(count,:)=w';
        resultSLL(count)=max(abs(A_SL1*w));
        save('resultd.mat', 'resultd');
        save('resultw.mat', 'resultw');
        save('SLL.mat', 'resultSLL');
        count=count+1;
        %plot(theta,20*log10(abs(A*x(2:9))));
        %axis([0,180,-70,0]);
        %hold on;
    end
end

% plot(theta,20*log10(abs(A*w))-max(20*log10(abs(A*w))));
% axis([0,180,-30,0]);
