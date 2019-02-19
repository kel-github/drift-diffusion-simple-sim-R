function [resp,RT,choiceRT,mrt] = simDiff;

 

%Quick and dirty method of generating a bunch of simulated choice-RT

 

%The next 4 lines set the values for the core parameters in the diffusion model: Drift rate, boundary separation, non-decision time, and the start-point of evidence accumulation.

v = 0.2; %drift rate

a = .04; %boundary separation

Ter = 0.2; %non-decision time

z = a/2; %unbiased start point

 

%The next line sets the value of standard deviation of within-trial noise in the evidence accumulation process. This is traditionally fixed to 0.1 (following Ratcliff, 1978), but others set s = 1. Either is fine, as this parameter “scales” the other parameters in the model (i.e., diffusion model parameters are only defined on a ratio scale, and so you can get the same predictions from different combinations of parameter values---if one parameter were to be mulitplied by a constant, the same predictions would obtain if all other parameters were multiplied (or divided) by the same value).

s = 0.1; %diffusion coefficient

 

dt = .001; %Time Step in Seconds

 

%This simulation uses Euler’s method to generate sample paths of the diffusion process. The simulation is constructed as a loop that will generate 10000 simulated sample paths. These paths provide the basis for generating predictions about proportions of correct/error responses and RT distributions.

nTrials = 10000; %number of simulated trials

t0 = 0;

 

%This is the loop that controls the evidence accumulation process. At each time step, the current evidence total (x) is compared against the value of the two absorbing boundaries (a and 0). If the evidence total is sitting between these to values the process has not yet terminated (i.e., neither boundary has been reached) and so more information is sampled. Each sample of information changes the evidence total by an amount that is fixed across all time steps (v * dt; or the mean Drift Rate multiplied by the size of the time step, here set to 1 ms) plus a (scaled) random draw from a normal distribution with mean 0 and standard deviation equal to s. If the process has reached/exceeded one of the absorbing boundaries, the number of time steps required to reach that threshold is recorded along with the absorbing boundary that was reached.

for i = 1:nTrials

    %initialize evidence

    x = z;

    t = t0;

    while x < a && x > 0

        x = x + v*dt + s*randn*sqrt(dt);

        t = t+dt;

    end

    if x > a;

        resp(i) = 1;

    elseif x < 0;

        resp(i) = 0;

    end

    RT(i) = t + Ter;

end

 

%These lines collate the simulated responses and RTs into a single variable, then generate histograms that show the distributions of correct and error RTs.

choiceRT = [resp',RT'];

figure;

hist(choiceRT(choiceRT(:,1)==1,2),100); %histogram of corrects

hold on;

hist(choiceRT(choiceRT(:,1)==0,2),100); %histogram of errors

 

%These lines summarize the mean RTs for the correct and error responses---this can be useful to summarize the model predictions in a more conventional manner.

mrt_c = mean(choiceRT(choiceRT(:,1)==1,2));

mrt_e = mean(choiceRT(choiceRT(:,1)==0,2));

mrt = [mrt_c, mrt_e];