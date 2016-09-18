clear all
clf

del_freq = zeros(1,14);
del_freq_unc = zeros(1,14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zeeman Splitting Effect  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 2:2:28    
    c1 = zeros(1,5);c2 = zeros(1,5);c3 = zeros(1,5);c4 = zeros(1,5);
    c1_unc = zeros(1,5);c2_unc = zeros(1,5);c3_unc = zeros(1,5);c4_unc = zeros(1,5);
    
    for i = 1:5
        file = ['Data\ZE_160204_',num2str(j),'_',num2str(i),'.csv'];
        A = csvread(file,1,0, [1 0 350 2]);
        voltages = A(:,2);
        empty = 1:length(voltages);
        
        scaler = j*3;
        
        if j > 16
            offset = 20;
        else
            offset = 0;
        end

        f_model = fittype('a1*exp(-(((x-c1)/l)^2)) + a2*exp(-(((x-c2)/l)^2)) + a1*exp(-(((x-c3)/l)^2)) + a2*exp(-(((x-c4)/l)^2)) + d',...   % Define your fitting function
                    'independent',{'x'},...                                             % Define the independent variable
                    'coefficients',{'a1','a2','c1','c2','c3','c4','l','d'});                  % Define the fit parameters

        opt = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',[0.04 0.02 0 10 150 150 5 .025],...                   % Lower bounds on parameters (order = a, b, c1, c2, l1, l2, d)
                    'Upper',[0.07 0.05 50 (70+scaler+offset) (230+offset) (300+scaler+offset) 20 .075],...                   % Upper bounds
                    'StartPoint',[0.055 0.035 20 (40+scaler+offset) (190+offset) (220+scaler+offset) 10 0.05]); 

        f_fit = fit(empty',voltages,f_model,opt);
        
        c = coeffvalues(f_fit);
        c_err = confint(f_fit,0.95);
        uncert = (c_err(2,:) - c_err(1,:))/2;
        
        c1_unc(i) = uncert(3);c2_unc(i) = uncert(4);c3_unc(i) = uncert(5);c4_unc(i) = uncert(6);
        c1(i) = c(3);c2(i) = c(4);c3(i) = c(5);c4(i) = c(6);
    end    
    inv_sigma1 = 0;inv_sigma2 = 0;inv_sigma3 = 0;inv_sigma4 = 0;
    for k=1:5
        inv_sigma1 = inv_sigma1 + (1/(c1_unc(k)^2));
        inv_sigma2 = inv_sigma2 + (1/(c2_unc(k)^2));
        inv_sigma3 = inv_sigma3 + (1/(c3_unc(k)^2));
        inv_sigma4 = inv_sigma4 + (1/(c4_unc(k)^2));
    end
    mean_sigma1 = 0;mean_sigma2 = 0;mean_sigma3 = 0;mean_sigma4 = 0;
    for k=1:5
        mean_sigma1 = mean_sigma1 + c1(k)/(c1_unc(k)^2);
        mean_sigma2 = mean_sigma2 + c2(k)/(c2_unc(k)^2);
        mean_sigma3 = mean_sigma3 + c3(k)/(c3_unc(k)^2);
        mean_sigma4 = mean_sigma4 + c4(k)/(c4_unc(k)^2);
    end
    means = [mean_sigma1 / inv_sigma1, mean_sigma2 / inv_sigma2, mean_sigma3 / inv_sigma3, mean_sigma4 / inv_sigma4]; 
    uncs = [1/sqrt(inv_sigma1), 1/sqrt(inv_sigma2), 1/sqrt(inv_sigma3), 1/sqrt(inv_sigma4)];
    FSR = means(3) - means(1);
    FSR_unc = sqrt((uncs(1)^2)+(uncs(3)^2));
    FSR_scaler = 8/FSR;
    space1 = means(2)-means(1);
    space2 = means(4)-means(3);
    space1_unc = (uncs(1)^2)+(uncs(2)^2);
    space2_unc = (uncs(3)^2)+(uncs(4)^2);
    del_freq(j/2) = (((space1/space1_unc)+(space2/space2_unc))/((1/space1_unc)+(1/space2_unc)))*FSR_scaler;    
    del_freq_unc(j/2) = (1/((1/space1_unc)+(1/space2_unc)))*FSR_scaler;
end

currents = [0,2,4,6,8,10,12,14,16,18,20]; 
B_current1 = [currents; -.005,.172,.369,.559,.695,.799,.896,.984,1.066,1.145,1.216];
B_current2 = [currents; -.004,.187,.373,.550,.691,.801,.895,.983,1.066,1.142,1.214];
B_current = zeros(2,11);
mag_unc = .1 / sqrt(11)

for i = 1:11
    B_current(2,i) = (B_current2(2,i) + B_current1(2,i))/2;
    B_current(1,i) = (i-1)*2;
end

Fit = polyfit(B_current(1,:),B_current(2,:),1); % fit B field versus current

B = zeros(1,14);
for i = 2:2:28
    B(i/2) = Fit(1)*i + Fit(2);
end

plot(B,del_freq,'b.')
hold on
Final_fit = polyfit(B,del_freq,1) % Fit B with delta frequency
unc = mag_unc*Final_fit(1)
y_fit = polyval(Final_fit,B)
plot(B,y_fit)
xlabel('B')
ylabel('\Delta \nu')
hold off

