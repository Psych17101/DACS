data=readtable("C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\AE\45\AE_45_2.csv");
data2=readtable("C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\AE\90\AE_90_3.csv");
%data.Properties.DimensionNames
%[number, Time, Chan, Load [N], E [eu]]
%plot(data{:,1}, data{:,5}, "*b"); xlabel("number"); ylabel("E");
%plot(data2{:,1}, data2{:,5}, "xr"); xlabel("number"); ylabel("E");
mean1=mean(data.E_eu_);
std1=std(data.E_eu_);

x=mean1-2*std1:1:mean1+2*std1;
y=normpdf(x, mean1, std1);
plot(x, y)

