% linear network simulation function
function [r,int]=linear_network_simulation_neg(w,t,tau_s,init,samplenum)

noiseamp=0.01;


delt=mean(diff(t));
n=size(w,2);
r=zeros(length(t),n);  % zero out the r vector.  will record something at each time step

r(1,:)=init;  % the first r is the one passed in, which is random

for i=2:length(t);
% the r(i-1,:)*W is actually like w_ij*r_j, since the r(i-1,:) is a "column vector" its confusing
    a=r(i-1,:)+delt/tau_s*(-r(i-1,:)+r(i-1,:)*w) + noiseamp*sqrt(delt)*randn(1,n);
% a constraint take out negative firing rates:  not used here
%     a(a<0)=0;
    r(i,:)=a;
end


% use all but the last vector of r, since its doing a difference (and it will has length of one less with the difference)
int=tau_s*diff(r)/delt+r(1:end-1,:);
r=r(1:end-1,:);


end
