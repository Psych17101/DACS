%% geometry, from the geometryfile (hardcoded for speed)
geomUD=[  250.1800   14.9867    0.8383
  250.1450   15.1367    0.9917
  250.1400   15.1767    1.0400
  250.1450   15.1100    1.0100
  250.1500   15.0900    0.9867
  250.1600   15.0967    0.9850
  250.1700   15.1400    1.0000
  250.1400   15.1467    1.0017
  250.1200   15.0800    1.0000
  250.1100   15.1567    0.9950
  250.1200   15.0600    0.9967
  250.1650   15.1233    0.9867
  250.1600   15.1433    0.9917
  250.1300   15.1467    1.0033
  250.1100   15.0833    1.0400
  250.1000   15.2333    1.0217
  250.0850   15.2167    0.8983];
geomUD=geomUD.*1e-3; %convert to [m]

geom90=[
  175.1150   25.0700    2.1233
  175.0950   24.9600    2.1450
  175.1150   25.0633    2.1133
  175.1000   25.0467    2.1133
  175.1000   25.0733    2.1250
  175.0950   25.0333    2.1217
  175.0800   25.0367    2.1217
  175.0850   25.0600    2.1400
  175.0650   25.0200    2.1383];
geom90=geom90.*1e-3;

geom45=[  250.0950   25.0367    2.0733
  250.0850   25.0433    2.0700
  250.0700   25.0200    2.0467
  250.0600   25.0300    2.0600
  250.0700   25.0833    2.0483
  250.0700   25.0800    2.0483
  250.0700   25.0500    2.0417
  250.0600   24.9467    2.0700
  250.0650   25.0467    2.0500];
geom45=geom45.*1e-3;

% in the following code, go to either MTS or DIC folder and run it before 
%working only with the sections below it
%% MTS folder
folder_path_UD = 'C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\MTS\UD\'; % Replace with your folder path
file_list = dir([folder_path_UD '*.csv']); % Get a list of all CSV files in the folder

table_listUD = cell(length(file_list), 1); % Initialize a cell array to store the table data

for i = 1:length(file_list)
    filename = [folder_path_UD file_list(i).name]; % Get the full file path
    table_data = readtable(filename); % Read the CSV file using readtable
    table_listUD{i} = table_data; % Store the table data in the cell array
    % Do something with the table data
end
samplesUD=[1,2,3,4,5,6,8,9,11,12,13];
%table_list{i} gives the table if index i, in this case the Samples used
%are 1,2,3,4,5,6,8,9,11,12,13
% its variable names are {'Var1'}{'Time_s_'}{'Load_N_'}{'Displacement_mm_'}

%45 folder
folder_path_45 = 'C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\MTS\45\'; % Replace with your folder path
file_list = dir([folder_path_45 '*.csv']); % Get a list of all CSV files in the folder

table_list45 = cell(length(file_list), 1); % Initialize a cell array to store the table data

for i = 1:length(file_list)
    filename = [folder_path_45 file_list(i).name]; % Get the full file path
    table_data = readtable(filename); % Read the CSV file using readtable
    table_list45{i} = table_data; % Store the table data in the cell array
    % Do something with the table data
end
samples45=[2,3,6,7,8];

%90 folder
folder_path_90 = 'C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\MTS\90\'; % Replace with your folder path
file_list = dir([folder_path_90 '*.csv']); % Get a list of all CSV files in the folder

table_list90 = cell(length(file_list), 1); % Initialize a cell array to store the table data
%length(file_list)
for i = 1:length(file_list)
    filename = [folder_path_90 file_list(i).name]; % Get the full file path
    table_data = readtable(filename); % Read the CSV file using readtable
    table_list90{i} = table_data; % Store the table data in the cell array
    % Do something with the table data
end
samples90=[3,4,5,6,7,8,9];

%Quasi folder
folder_path_Quasi = 'C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\MTS\Quasi\'; % Replace with your folder path
file_list = dir([folder_path_Quasi '*.csv']); % Get a list of all CSV files in the folder

table_listQuasi = cell(length(file_list), 1); % Initialize a cell array to store the table data

for i = 1:length(file_list)
    filename = [folder_path_Quasi file_list(i).name]; % Get the full file path
    table_data = readtable(filename); % Read the CSV file using readtable
    table_listQuasi{i} = table_data; % Store the table data in the cell array
    % Do something with the table data
end
samplesQuasi=[1,2,3,5,6,7,8,9,10];


%% calculating E1 by using the UD

AreasUD=[];
LengthsUD=[];
for i= 1:length(samplesUD)
    AreasUD=[AreasUD, geomUD(i, 2).*geomUD(i, 3)]; %add only the areas of the samples used
    LengthsUD=[LengthsUD, geomUD(:, 1)];
end

%AtimesLUD=AreasUD.*LengthsUD; %This is the denominator when calculating E1, different for each sample

E1table= cell(length(AreasUD), 1); %for each sample, there will be calculated E1
%variable names are {'Var1'}{'Time_s_'}{'Load_N_'}{'Displacement_mm_'}
%we cut off the first 100 and linear behaviour until index=14000
%we remove the values because we're interested in the linear behaviour of
%E1
for i=1:length(AreasUD)
    E1=(table_listUD{i}{100:14000, 3}.*LengthsUD(i))./(table_listUD{i}{100:14000, 4}.*1e-3*AreasUD(i));    
    E1table{i}=E1;
end

%% calculating the gaussian distribution variables for E1
%meanE1=mean(E1table{1}); 
%stdE1=std(E1table{1});
%using all values:
E1tot=[];
for i=1:length(AreasUD)
    E1tot=[E1tot; E1table{i}];
end
meanE1s=[];
for i=1:length(AreasUD)
    meanE1s=[meanE1s, mean(E1table{i})];
end

meanE1=mean(meanE1s); %result:165.22 GPa
stdE1=std(meanE1s); %16.922 GPa


%x=meanE1-2*stdE1:1e8:meanE1+2*stdE1;
%plot(x./1e9, normpdf(x, meanE1, stdE1));xlabel("E [Gpa]");


%% calculating E2 by using the 90

Areas90=[];
Lengths90=[];
for i= 1:length(samples90)
    Areas90=[Areas90, geom90(i, 2).*geom90(i, 3)]; %add only the areas of the samples used
    Lengths90=[Lengths90, geom90(:, 1)];
end

%AtimesLUD=AreasUD.*LengthsUD; %This is the denominator when calculating E1, different for each sample

E2table= cell(length(Areas90), 1); %for each sample, there will be calculated E1
%variable names are {'Var1'}{'Time_s_'}{'Load_N_'}{'Displacement_mm_'}
%we cut off the first 100 and linear behaviour until index=14000
%we remove the values because we're interested in the linear behaviour of
%E1
for i=1:length(Areas90)
    E2=(table_list90{i}{100:7000, 3}.*Lengths90(i))./(table_list90{i}{100:7000, 4}.*1e-3*Areas90(i));
    E2table{i}=E2;
end

%% calculating the gaussian distribution variables for E2
%meanE2=mean(E2table{1}); %mean in Gpa
%stdE2=std(E2table{1});
E2tot=[];
for i=1:length(Areas90)
    E2tot=[E2tot; E2table{i}];
end
meanE2s=[];
for i=1:length(Areas90)
    meanE2s=[meanE2s, mean(E2table{i})];
end
meanE2=mean(meanE2s); %result:8.444 GPa
stdE2=std(meanE2s); %1.077 GPa
%x=meanE2-2*stdE2:1e8:meanE2+2*stdE2;
%plot(x./1e9, normpdf(x, meanE2, stdE2));xlabel("E [Gpa]");

%% calculating Ex to obtain E2=E1 and G12 for the 45

Areas45=[];
Lengths45=[];
for i= 1:length(samples45)
    Areas45=[Areas45, geom45(i, 2).*geom45(i, 3)]; %add only the areas of the samples used
    Lengths45=[Lengths45, geom45(:, 1)];
end

%AtimesLUD=AreasUD.*LengthsUD; %This is the denominator when calculating E1, different for each sample

E21table= cell(length(Areas45), 1); %for each sample, there will be calculated E1
G12table=cell(length(Areas45), 1);
%variable names are {'Var1'}{'Time_s_'}{'Load_N_'}{'Displacement_mm_'}
%we cut off the first 100 and linear behaviour until index=14000
%we remove the values because we're interested in the linear behaviour of
%E1
for i=1:length(Areas45)
    E21=(table_list45{i}{100:14000, 3}.*Lengths45(i))./(table_list45{i}{100:14000, 4}.*1e-3*Areas45(i));
    E21table{i}=E21./(1-meanvxy); %using elastic properties for specified lamina
    G12table{i}=E21./(2.*(1+meanvxy));%the value of E21 is equal to G12 of the UD as the 45 experiences pure shear
end

%using basic transformation matrix with theta=45 and only sigmax non-zero
%gives E1=E2=Ex/(1-vxy)

E21=E21./(1-meanvxy);
%% calculating the gaussian distribution variables for E21
%this is basically souble checking as the E21 can be calculated from E1

%meanE21=mean(E21table{1}); %mean in Gpa
%stdE21=std(E21table{1});

E21tot=[];
for i=1:length(Areas45)
    E21tot=[E21tot; mean(E21table{i})];
end

meanE21=mean(E21tot); %result: 79.317 GPa
stdE21=std(E21tot); %13.708 GPa 

%x=meanE21-2*stdE21:1e8:meanE21+2*stdE21;
%plot(x./1e9, normpdf(x, meanE21, stdE21));xlabel("E [Gpa]");

%% calculating the gaussian distribution variables for G12
G12tot=[];
for i=1:length(Areas45)
    G12tot=[G12tot; mean(G12table{i})];
end

meanG12=mean(G12tot); %result: 6.4167 GPa
stdG12=std(G12tot); %1.1089 GPa

%x=meanE21-2*stdE21:1e8:meanE21+2*stdE21;
%plot(x./1e9, normpdf(x, meanE21, stdE21));xlabel("E [Gpa]");



%% calculating S12 using +-45
S12table= [];
for i =1:length(samples45)
    maxload=max(table_list45{i}{:, 3});%retrieves the maximum load
    S12table=[S12table, maxload./2./Areas45(i)];
end

meanS12=mean(S12table); %result:152.4 MPa
stdS12=std(S12table); %1.7844 MPa





%% 
for i=1:length(v21table)
    mean(v21table{i})
end


%% calculating Xt for UD laminate, using UD
%plot(1:1:length(table_listUD{1}{:, 3}), table_listUD{1}{:, 3})
Xttable= [];
for i =1:length(samplesUD)
    maxload=max(table_listUD{i}{:, 3});%retrieves the maximum load
    Xttable=[Xttable, maxload./AreasUD(i)];
end

meanXt=mean(Xttable); %result:1923.7 MPa
stdXt=std(Xttable); %108.65

%% calculating Yt for UD laminate, using 90
Yttable= [];
for i =1:length(samples90)
    maxload=max(table_list90{i}{:, 3});%retrieves the maximum load
    Yttable=[Yttable, maxload./Areas90(i)];
end
meanYt=mean(Yttable); %result:107.2 MPa
stdYt=std(Yttable); %9.3507 MPa


%% calculating Xt=Yt for 45, using 45
XYttable= [];
for i =1:length(samples45)
    maxload=max(table_list45{i}{:, 3});%retrieves the maximum load
    XYttable=[XYttable, maxload./Areas45(i)];
end
meanXYt=mean(XYttable); %result:305 MPa
stdXYt=std(XYttable);



%% DIC folder
folder_path_UD = 'C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\DIC\UD\'; % Replace with your folder path
file_list = dir([folder_path_UD '*.csv']); % Get a list of all CSV files in the folder

table_listUD = cell(length(file_list), 1); % Initialize a cell array to store the table data

for i = 1:length(file_list)
    filename = [folder_path_UD file_list(i).name]; % Get the full file path
    table_data = readtable(filename); % Read the CSV file using readtable
    table_listUD{i} = table_data; % Store the table data in the cell array
    % Do something with the table data
end
samplesUD=[1,2,3,4,5,6,8,9,11,12,13];
%table_list{i} gives the table if index i, in this case the Samples used
%are 1,2,3,4,5,6,8,9,11,12,13
% its variable names are {'Var1'}{'Time_s_'}{'Load_N_'}{'Displacement_mm_'}

%45 folder
folder_path_45 = 'C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\DIC\45\'; % Replace with your folder path
file_list = dir([folder_path_45 '*.csv']); % Get a list of all CSV files in the folder

table_list45 = cell(length(file_list), 1); % Initialize a cell array to store the table data

for i = 1:length(file_list)
    filename = [folder_path_45 file_list(i).name]; % Get the full file path
    table_data = readtable(filename); % Read the CSV file using readtable
    table_list45{i} = table_data; % Store the table data in the cell array
    % Do something with the table data
end
samples45=[2,3,6,7,8];

%90 folder
folder_path_90 = 'C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\DIC\90\'; % Replace with your folder path
file_list = dir([folder_path_90 '*.csv']); % Get a list of all CSV files in the folder

table_list90 = cell(length(file_list), 1); % Initialize a cell array to store the table data

for i = 1:length(file_list)
    filename = [folder_path_90 file_list(i).name]; % Get the full file path
    table_data = readtable(filename); % Read the CSV file using readtable
    table_list90{i} = table_data; % Store the table data in the cell array
    % Do something with the table data
end
samples90=[3,4,5,6,7,8,9];

%% calculating v12 major poisson ratio using UD

v12table= cell(length(samplesUD), 1);


for i=1:length(samplesUD)
    v12=-1.*str2double(table_listUD{i}{17:end-15, 2})./str2double(table_listUD{i}{17:end-15, 1});
    v12table{i}=v12;
end
%% calc mean
v12tot=[];
for i=1:length(samplesUD)
    v12tot=[v12tot; v12table{i}];
end

meanv12=mean(v12tot); %result:0.35318
stdv12=std(v12tot); %0.1809


%% calculating v21 minor poisson ratio using 90

v21table= cell(length(samples90), 1);


for i=1:length(samples90)
    v21=-1.*str2double(table_list90{i}{17:end-15, 2})./str2double(table_list90{i}{17:end-15, 1});
    v21table{i}=abs(v21);
end
%% calc mean
v21tot=[];
for i=1:length(samples90)
    v21tot=[v21tot; v21table{i}];
end

meanv21=mean(v21tot); %result:0.018489
stdv21=std(v21tot); %0.031465

%% validation

meanv21theory=meanv12*meanE2/meanE1;
relError=(meanv21theory-meanv21)/meanv21;%2percent error, almost too good



%% calculating vxy=vyx for the 45
vxytable= cell(length(samples45), 1);


for i=1:length(samples45)
    vxy=-1*str2double(table_list45{i}{10:end-15, 2})./str2double(table_list45{i}{10:end-15, 1});
    vxytable{i}=vxy;
end
%% calc mean
vxytot=[];
for i=1:length(samples45)
    vxytot=[vxytot; vxytable{i}];
end

meanvxy=mean(vxytot); %result:0.72147
stdvxy=std(vxytot); %result: 0.10113

%% calculating G12 of UD
G12UDtable= cell(length(samplesUD), 1);


for i=1:length(samplesUD)
    A=table_listUD{i}.Variables;
    G12UDtable{i}=str2double(A(10:end-15, 6))./(log(1+abs(str2double(A(10:end-15, 3)))).*AreasUD(i));
    %G12UDtable{i}=abs(2*str2double(A(10:end-15, 6)).*1e3./((str2double(A(10:end-15, 1))-str2double(A(10:end-15, 2))).*geomUD(i,3)));
end
%plot(1:1:length(G12UDtable{1}),G12UDtable{1}) 
%% calculating means and std
meanG12UD=[];
    
for i =1:length(samplesUD)
    meanG12UD=[meanG12UD, mean(G12UDtable{i})];
end
meanG12UDtot=mean(meanG12UD);
stdG12UDtot=std(meanG12UD);
%% same with G12 of 90 to double check

G1290table= cell(length(samples90), 1);


for i=1:length(samples90)
    A=table_list90{i}.Variables;
    G1290table{i}=str2double(A(10:end-15, 6))./(str2double(A(10:end-15, 3)).*geomUD(i,3));
    %G1290table{i}=(2*str2double(A(10:end-15, 6)).*1e3./((str2double(A(10:end-15, 1))-str2double(A(10:end-15, 2))).*geomUD(i,3)));
end

%% calculating means and std
meanG1290=[];
    
for i =1:length(samples90)
    meanG1290=[meanG1290, mean(G1290table{i})];
end
meanG1290tot=mean(meanG1290);
stdG1290tot=std(meanG1290);


%%

Vxy=meanE1/2/(1+0.35);
Vyx=meanE2/2/(1+0.0185);
b=1-(0.35*0.0185);
G12UDtable= cell(length(samplesUD), 1);


for i=1:length(samplesUD)
    A=table_listUD{i}.Variables;
    G12UDtable{i}=str2double(A(10:end-15, 6))./AreasUD(i).*str2double(A(10:end-15, 3)).*b./(Vxy.*str2double(A(10:end-15, 1))+Vyx.*str2double(A(10:end-15, 2)));
    %G12UDtable{i}=abs(2*str2double(A(10:end-15, 6)).*1e3./((str2double(A(10:end-15, 1))-str2double(A(10:end-15, 2))).*geomUD(i,3)));
end
meanG12UD=[];
    
for i =1:length(samplesUD)
    meanG12UD=[meanG12UD, mean(G12UDtable{i})];
end
meanG12UDtot=mean(meanG12UD);
stdG12UDtot=std(meanG12UD);