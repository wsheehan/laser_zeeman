clear all
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak Spacing for laser 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1 = zeros(1,10);c2 = zeros(1,10);c3 = zeros(1,10);c4 = zeros(1,10);
c1_unc = zeros(1,10);c2_unc = zeros(1,10);c3_unc = zeros(1,10);c4_unc = zeros(1,10);

for i = 1:10
    file = ['Data\ZE_160202_0',num2str(i),'.csv'];
    A = csvread(file,1,0, [1 0 350 2]);
    voltages = A(:,2);
    empty = 1:length(voltages);
    
    f_model = fittype('a1*exp(-(((x-c1)/l)^2)) + a2*exp(-(((x-c2)/l)^2)) + a1*exp(-(((x-c3)/l)^2)) + a2*exp(-(((x-c4)/l)^2)) + d',...   % Define your fitting function
                'independent',{'x'},...                                             % Define the independent variable
                'coefficients',{'a1','a2','c1','c2','c3','c4','l','d'});                  % Define the fit parameters
            
    opt = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0.04 0.02 35 50 240 280 5 .025],...                   % Lower bounds on parameters (order = a, b, c1, c2, l1, l2, d)
                'Upper',[0.07 0.05 55 80 270 300 30 .075],...                   % Upper bounds
                'StartPoint',[0.055 0.035 45 60 250 290 15 0.05]); 
            
    f_fit = fit(empty',voltages,f_model,opt);
    %subplot(2,5,i);
    plot(voltages,'k.')
    hold on
    plot(f_fit(empty),'b')
    ylabel('Voltage (V)')
    xlabel('Frame of View')
    hold off
    
    c = coeffvalues(f_fit);
    c_err = confint(f_fit,0.95);
    uncert = (c_err(2,:) - c_err(1,:))/2;
    
    c1_unc(i) = uncert(3);c2_unc(i) = uncert(4);c3_unc(i) = uncert(5);c4_unc(i) = uncert(6);  
    c1(i) = c(3);c2(i) = c(4);c3(i) = c(5);c4(i) = c(6);
end

inv_sigma1 = 0;inv_sigma2 = 0;inv_sigma3 = 0;inv_sigma4 = 0;

for i=1:10
    inv_sigma1 = inv_sigma1 + (1/(c1_unc(i)^2));
    inv_sigma2 = inv_sigma2 + (1/(c2_unc(i)^2));
    inv_sigma3 = inv_sigma3 + (1/(c3_unc(i)^2));
    inv_sigma4 = inv_sigma4 + (1/(c4_unc(i)^2));
end

mean_sigma1 = 0;mean_sigma2 = 0;mean_sigma3 = 0;mean_sigma4 = 0;

for i=1:10
    mean_sigma1 = mean_sigma1 + c1(i)/(c1_unc(i)^2);
    mean_sigma2 = mean_sigma2 + c2(i)/(c2_unc(i)^2);
    mean_sigma3 = mean_sigma3 + c3(i)/(c3_unc(i)^2);
    mean_sigma4 = mean_sigma4 + c4(i)/(c4_unc(i)^2);
end

means = [mean_sigma1 / inv_sigma1, mean_sigma2 / inv_sigma2, mean_sigma3 / inv_sigma3, mean_sigma4 / inv_sigma4] 
uncs = [1/sqrt(inv_sigma1), 1/sqrt(inv_sigma2), 1/sqrt(inv_sigma3), 1/sqrt(inv_sigma4)];

FSR = means(3) - means(1)
FSR_unc = sqrt((uncs(1)^2)+(uncs(3)^2))

FSR_scaler = 8/FSR;

peak1 = means(2)-means(1);
peak2 = means(4)-means(3);

peak1_unc = (uncs(1)^2)+(uncs(2)^2);
peak2_unc = (uncs(3)^2)+(uncs(4)^2);

mean_spacing = ((peak1/peak1_unc)+(peak2/peak2_unc))/((1/peak1_unc)+(1/peak2_unc))*FSR_scaler
spacing_unc = 1/((1/peak1_unc)+(1/peak2_unc))*FSR_scaler



