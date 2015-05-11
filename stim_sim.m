clear all; clf

duration=6;tau_s=1;delt=.1;

t=0:delt:duration;

 

% generate random connectivity
n=100
w=randn(n,n);w=real(w/max((eig(w))))*.9;
n=size(w,2);

dub=w; % setting the "dub" variable equal to the matrix of random weights just made

 

nums=[1:15];  %number of perturbations

 
figure(1)
clf(1)

 

for M=1:length(nums)

    r=[];int=[];

    for k=1:M %run stimulation sequence M times
	percentageperturb=1;
        init(:,k)=1*randn(n,1).*(rand(n,1)<percentageperturb); % make a "column vector" of length n to pass into the simulation
	 	

	% perform a "linear network simulation" with the weight matrix dub, for times t, with constant tau_s, and with the random initialized column vector
        [R,INT]=linear_network_simulation_neg(dub,t,tau_s,init(:,k)); 
	% add the results of the simulation to the stored variables.  it stores "r" and "int" 
	% both of these are stored as functions of time.  So they will be long!
        r=[r;R];int=[int;INT];
	% singular value decomposition of the matrix r.  Not sure why this is done here, inthe loop.  It looks 
        [u,s,v]=svd(r);

    end

   

   
	% take the pseudoinverse of "r".  r is Mxn, and int is Mxn, so r^-1 times int is an n by n matrix
% so , r is some measure of activity output from the simulation, and so is int... not sure what these are
% r is probably average rates of each one? (rates of each neuron in the nxn network.  maybe int is the output?
    beta=pinv(r)*int;  %connectivity derived from activity

   

    %assess quality of fit

    d=(dub-beta).^2;a=dub.^2;a=mean(a(:));

    d=sqrt(d/a)*100;  % this is a normalization, since a is a constant

% use the median different as a measure of the goodness of fit
    D(M)=median(d(:));  %real vs derived difference

% use SVD to get "eigenvalues" of the activity modes.  
% the magnitude is approximated as the real part of the sqrt of the eigenvalues of r'*r,
    sv(:,M)=flipud(real(sqrt(eig(r'*r))));  %singular values of activity modes

% visualize the differences
    figure(1); imagesc((dub-beta),[-0.1 0.1]);title(num2str(M)); pause(0.1)

end

% so, the above loop seems like its applying a sort of random stimulation to a network.  Then, its using an eigenvalue decomposition to see how close an estimated connectivity resembles the true connectivity (which here was generated as random).  Based on what Emre said, after about 10 of these random stimulations (perturbations), estimated connectivity based on using the pinv(r)*int gets close to the actual connectivity, dub
% plot(D) would be useful, to see how the median measure of error of estimating the connectivity decreases with more random perturbations

 

 
% visualize the vectors of activity after each random perturbation
% x:  neuron number
% y:  activity measure
figure(2); plot(r)

 

%%

figure(3);clf

% sv contains the eigenvalues of the activity modes, calculated as more perturbations are applied
% normalize this 
cn=max(sv)./min(sv);
% take out infinite values and get some type of log of the max value, rounded up to the nearest integer. ....
g=cn;g(g==inf)=0;g=ceil(log10(max(g)));
% make infinite ones large, for displaying
cn(cn==inf)=10^9;

 

% set(gca,'yticklabel',{'2','4','6','8','\infty'})

% this plots the decrease of D, the median of the difference between the "true" and "estimated" connectivity matrices, as a function of the number of perturbations.  It probably decreases very sharply at first and then not so much around k=10, because that is where Emre said the fit starts to level off
figure(3);
plot(D,'k.-','linewidth',2,'markersize',20);
set(gca,'box','off');axis square
title('Matrix Difference Measure')
xlabel('number of stims')

