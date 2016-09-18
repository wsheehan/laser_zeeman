clf
file = ['data\ZE_160202_env_',num2str(1),'.csv'];
A = csvread(file,1,0);
A = A(
voltages=A(:,2)
PZT=A(:,3)
plot(voltages)
hold on
plot(PMT)
xlabel('Frame of View')
ylabel('Voltage (V)')
hold off