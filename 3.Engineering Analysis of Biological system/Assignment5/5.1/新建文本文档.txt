%% HW 5
%% read in data
format compact
data =readmatrix('MM_data.xlsx');
%set up xobs, yobs into vectors


%% parameter guesses; define beta0

fnameFOR=@MM2; %%  whatever function name you wish to use

%% Scaled sensitivity coefficients before running inverse problem
%to determine a) which parameters can be estimated, and b) which parameters
%wil be most accurate (lowest relative error)
xs=min(t):1:max(t); %xs are the times for SSCs to make a smooth curve.
xs=xs';%make xs a column.
ns=length(xs);%length of xs for plotting
Xp=SSC_V3(beta0,xs,fnameFOR);

%% plot X' for each dependent variable
%plot for V
cmap = ['r' 'g' 'b' 'c' 'y'  'm' 'k' ]';
figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
%plot V vs t to know the total span

ypred=fnameFOR(beta,xs);
h2(1)=plot(xs(1:ns),ypred(1:ns),'-','color',cmap(1,:),'LineWidth',3.5); %plot the predicted C to compare to SSCs
for i=1:p
    h2(i+1) = plot(xs(1:ns),Xp(1:ns,i),'-','color',cmap(i+1,:),'LineWidth',2);
end
legend('X''_V','X''_V_{max}','X''_{K_m}')


%nlinfit returns parameters, residuals, Jacobian (sensitivity coefficient matrix), 
%covariance matrix, and mean square error
fnameINV= %%put the name of your inverse function here
[beta,resids,J,COVB,mse] = nlinfit(,fnameINV,beta0);
Vmax=beta(1)
Km=beta(2)
rmse=sqrt(mse)
rmseunits='mmolarpermin'
condX=cond(J)
detXTX=det(J'*J)
%R is the correlation matrix for the parameters, sigma is the standard deviation vector
[R,sigma]=corrcov(COVB);
R
corr12=R(2,1) %correlation coefficient between the 2 parameters
sigma
relerr=sigma./beta'
relerr1=100*relerr(1) %percent error beta1
relerr2=100*relerr(2) %percent error beta2

ypredp=MM2(beta,xs);
yspan=range(ypredp)
relrmse=rmse/yspan

%% confidence intervals for parameters
%use nlparci



%nonlinear regression confidence intervals-- 'on' means simultaneous
%bounds; 'off' is for nonsimultaneous bounds; must use 'curve' for
%regression line, 'observation' for prediction interval
[ypred, delta] = nlpredci(fnameINV,t,param,resids,J,0.05,'on','curve'); %confidence band for regression line
[ypred, deltaob] =nlpredci(fnameINV,t,param,resids,J,0.05,'on','observation');%prediction band for individual points

%% simultaneous confidence bands for regression line
CBu=Cpred+delta;
CBl=Cpred-delta;

%simultaneous prediction bands for regression line
PBu=Cpred+deltaob;
PBl=Cpred-deltaob;

V1000=fnameFOR(beta,1000); %ypred @ 1000
S1000=find(xobs==1000,1);
CBu1000=CBu(S1000)
PBl1000=PBl(S1000)

%plot yobs, ypred line, confidence band for regression line
ypredp=fnameFOR(beta,xs); %gives smooth line for ypred 

h1(1)=plot(xobs,yobs,'square', 'Markerfacecolor', 'b','markersize',10); %observed
h1(2) = plot(xobs,ypredp,'-k','LineWidth',2); %predicted
h1(3) = plot(xobs,CBu,'--r','LineWidth',2); %upper CB
plot(t,CBl,'--r','LineWidth',2); %lower CB doesn't need a legend

%plot prediction band for regression line
h1(4) = plot(xobs,PBu,'-.g','LineWidth',2); %upper PB
plot(xobs,PBl,'-.g','LineWidth',2);  %lower PB doesn't need a legend

meanres=mean(resids);
resid10=           ;
%% Scaled sensitivity coefficients before running inverse problem
%to determine a) which parameters can be estimated, and b) which parameters
%wil be most accurate (lowest relative error)

Xp=SSC_V3(beta,xs,fnameFOR);

%% plot X' for each dependent variable
%plot for C
cmap = ['r' 'g' 'b' 'c' 'y'  'm' 'k' ]';
figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
%plot C vs t to know the total span
ypred=fnameFOR(beta,xs);
h2(1)=plot(xs(1:ns),ypred(1:ns),'-','color',cmap(1,:),'LineWidth',3.5); %plot the predicted C to compare to SSCs
for i=1:p
    h2(i+1) = plot(xs(1:ns),Xp(1:ns,i),'-','color',cmap(i+1,:),'LineWidth',2);
end
legend('X''_C','X''_V_{max}','X''_{K_m}')
xlabel('time'); ylabel('scaled sensitivity coefficient for C, gmol/L');
t1000=find(xs==1000)
SSC1_1000=Xp(t1000,1)
SSC2_1000=Xp(t1000,2)

%% residual scatter plot
figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
plot(t, resids, 'square', 'Markerfacecolor', 'b')
plot([0,max(t)],[0,0], 'R')
ylabel('Observed y - Predicted y','fontsize',16,'fontweight','bold')
xlabel('S (mmol)','fontsize',16,'fontweight','bold')


%% residuals histogram--same as dfittool, but no curve fit here
 figure
h=histogram(resids)
 xlabel('Y_{observed} - Y_{predicted}','fontsize',16,'fontweight','bold')
 ylabel('Frequency','fontsize',16,'fontweight','bold')
set(gca, 'fontsize',20,'fontweight','bold');

%% standard statistical assumptions

function Vpred = MM2(beta,S)
%Michaelis-Menten model
Vpred=
end

function Xp=SSC_V3(beta,x,yfunc)
%Computes scaled sensitivity coefficients =Xp, nxp matrix
%can have k dependent variables that are stacked in a column vector
%all y1s, then all y2s, ...last are yks
%n is the number of data
%p is the number of parameters
%Xp1 = dY/dbeta1~[y(beta1(1+d), beta2,...betap) - y(beta1,
%beta2,...betap)]/d...
%d is the arbitrary delta, usually 0.001
%beta is the p x 1 parameter vector
%yhat is nx1 vector, the y values when only one parameter has been successively perturbed by d
%ypred is nx1 vector,  the y values when parameters are set to beta
%betain is px1 vector, the parameter values with only one parameter perturbed by d
%x are the independent variables (can be one or more independent variables)
%yfunc is a function (m file or an anonymous) defined by the user outside
%of this file

%% X' = scaled sensitivity coefficients using forward-difference
% This is a forward problem with known approximate parameters
d=0.001;
ypred=yfunc(beta,x);
for i = 1:length(beta)  %scaled sens coeff for forward problem
    betain = beta; %reset beta
    betain(i) = beta(i)*(1+d);
    yhat{i} = yfunc(betain,x);
    SSC{i} = (yhat{i}-ypred)/d;%scaled sens coeff for ith parameter
    Xp(:,i)=SSC{i}; %extract from cell array to 2D array
end
end
