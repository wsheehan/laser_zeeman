clf
file = ['data\ZE_160202_env_',num2str(1),'.csv'];
A = csvread(file,1,0,[1 0 370 2]);
plot(voltages,'.r')
xlabel('Data Scale')
ylabel('Voltage (V)')
hold off