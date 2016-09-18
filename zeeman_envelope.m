clear all
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Painting the Gain Curve               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width = zeros(1,5);
FWHM = zeros(1,5);

for i = 1:5
    file = ['data\ZE_160202_env_',num2str(i),'.csv'];
    A = csvread(file,1,0,[1 0 370 2]);
    voltages = A(:,2);
    voltages = voltages(2:2:end,:);
    empty = 1:2:2*length(voltages);
    
    f_model = fittype('a*exp(-(((x-c1)/l)^2)) + a*exp(-(((x-c2)/l)^2)) + d',...   % Define your fitting function
                'independent',{'x'},...                                             % Define the independent variable
                'coefficients',{'a','c1','c2','l','d'});                  % Define the fit parameters
            
    opt = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0.035 112 345 21 .04],...                   % Lower bounds on parameters (order = a, b, c1, c2, l1, l2, d)
                'Upper',[0.095 135 365 50 .075],...                   % Upper bounds
                'StartPoint',[0.075 125 355 35 0.063]); 
            
    f_fit = fit(empty',voltages,f_model,opt);
    %subplot(2,3,i);
    plot(voltages,'k.')
    hold on
    plot(f_fit(empty),'b')
    ylabel('Voltage (V)')
    xlabel('Data Scale')
    hold off
    
    c = coeffvalues(f_fit);
    c_err = confint(f_fit,0.95);
    uncert = (c_err(2,:) - c_err(1,:))/2;
    
    FSR = c(3) - c(2); % Define FSR
    FSR_unc = sqrt(uncert(3)^2 + uncert(2)^2);
    
    FSR_scaler = 8/FSR;
    
    st_dev = c(4)/sqrt(2); % standard deviation of gaussian
    st_dev_unc = uncert(4)/sqrt(2);
    
    FWHM(i) = 2.3548*st_dev*FSR_scaler; %equation that relates stan dev to full width half max
    FWHM_unc(i) = 2.3548*st_dev_unc*FSR_scaler
    
end

% We now have unc and spacing so we solve using
% Lyons equation

inv_sigma = 0;
mean_sigma = 0;
for i = 1:5
    inv_sigma = inv_sigma + 1/(FWHM_unc(i)^2);
    mean_sigma = mean_sigma + FWHM(i)/(FWHM_unc(i)^2);
end

mean = mean_sigma/inv_sigma
unc = sqrt(1/inv_sigma)
    





    