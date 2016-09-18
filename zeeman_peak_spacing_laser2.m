clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak Spacing for laser 2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1 = zeros(1,5);c2 = zeros(1,5);c3 = zeros(1,5);
c1_unc = zeros(1,5);c2_unc = zeros(1,5);c3_unc = zeros(1,5);

for i = 81:1:85
    file = ['Data\ZE_160204_0',num2str(i),'.csv'];
    A = csvread(file,1,0, [1 0 350 2]);
    voltages = A(:,2);
    empty = 1:length(voltages);
    
    f_model = fittype('a*exp(-(((x-c1)/l)^2)) + a*exp(-(((x-c2)/l)^2)) + a*exp(-(((x-c3)/l)^2)) + d',...   % Define your fitting function
                'independent',{'x'},...                                             % Define the independent variable
                'coefficients',{'a','c1','c2','c3','l','d'});                  % Define the fit parameters
            
    opt = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0.07 -10 165 350 30 .025],...                   % Lower bounds on parameters (order = a, b, c1, c2, l1, l2, d)
                'Upper',[0.11 10 195 390 80 .075],...                  % Upper bounds
                'StartPoint',[0.09 0 180 370 60 0.05]); 
            
    f_fit = fit(empty',voltages,f_model,opt);
    subplot(1,5,(i-80));
    plot(voltages,'.')
    hold on
    plot(f_fit(empty),'r')
    hold off
    
    c = coeffvalues(f_fit)
    c_err = confint(f_fit,0.95);
    uncert = (c_err(2,:) - c_err(1,:))/2;
    
    c1_unc(i-80) = uncert(2);c2_unc(i-80) = uncert(3);c3_unc(i-80) = uncert(4);    
    c1(i-80) = c(2);c2(i-80) = c(3);c3(i-80) = c(4);
end

inv_sigma1 = 0;inv_sigma2 = 0;inv_sigma3 = 0;

for i=1:5
    inv_sigma1 = inv_sigma1 + (1/(c1_unc(i)^2));
    inv_sigma2 = inv_sigma2 + (1/(c2_unc(i)^2));
    inv_sigma3 = inv_sigma3 + (1/(c3_unc(i)^2));
end

mean_sigma1 = 0;mean_sigma2 = 0;mean_sigma3 = 0;

for i=1:5
    mean_sigma1 = mean_sigma1 + c1(i)/(c1_unc(i)^2);
    mean_sigma2 = mean_sigma2 + c2(i)/(c2_unc(i)^2);
    mean_sigma3 = mean_sigma3 + c3(i)/(c3_unc(i)^2);
end

means = [mean_sigma1 / inv_sigma1, mean_sigma2 / inv_sigma2, mean_sigma3 / inv_sigma3] 
uncs = [1/sqrt(inv_sigma1), 1/sqrt(inv_sigma2), 1/sqrt(inv_sigma3)];

FSR = means(3) - means(1)
FSR_unc = sqrt((uncs(1)^2)+(uncs(3)^2))

FSR_scaler = 8/FSR;

MOD1 = means(2)-means(1);
MOD2 = means(3)-means(2);

MOD1_unc = (uncs(1)^2)+(uncs(2)^2);
MOD2_unc = (uncs(2)^2)+(uncs(3)^2);

spacing = (((MOD1/MOD1_unc)+(MOD2/MOD2_unc))/((1/MOD1_unc)+(1/MOD2_unc)))*FSR_scaler
spacing_unc = (1/((1/MOD1_unc)+(1/MOD2_unc)))*FSR_scaler