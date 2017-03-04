clear;

%% Configuration
wavelength=1;
k=2*pi/wavelength;   % wave number
angleStep=0.05;
theta=0:angleStep:360;
elementNumber=16;
spacingMIN=0.5;
spacingMAX=5;
mainbeam=5;

%%

d=[spacingMIN/2 + (spacingMAX/2-spacingMIN/2).*rand(1), spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber/2-1)];
%d=[0.25,0.5*ones(1,elementNumber/2-1)];

for nn=2:length(d)
    d(nn)=d(nn-1)+d(nn);
end

A=zeros(length(theta),elementNumber);

for nn=1:length(d)
    A(:,nn)=exp(-1i*k*d(length(d)-nn+1)*cosd(theta));
end

for nn=1:length(d)
    A(:,nn+length(d))=exp(1i*k*d(nn)*cosd(theta));
end

theta_SL=0:angleStep:90-mainbeam/2;
theta_ML=90+mainbeam/2:angleStep:180;

A_SL=zeros(length(theta_SL),elementNumber);

for nn=1:length(d)
    A_SL(:,nn)=exp(-1i*k*d(length(d)-nn+1)*cosd(theta_SL));
end

for nn=1:length(d)
    A_SL(:,nn+length(d))=exp(1i*k*d(nn)*cosd(theta_SL));
end

n=elementNumber;

%%
cvx_begin
variable w(n)
minimize( max(abs(A_SL*w)) )
subject to
ones(1,elementNumber)*w==1;
cvx_end

%%
plot(theta,20*log10(abs(A*w))-max(20*log10(abs(A*w))));
axis([0,180,-30,0]);
