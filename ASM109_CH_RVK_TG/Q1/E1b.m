geom=[  250.0350   25.0600    2.0050
  250.0450   25.0467    2.0783
  250.0450   25.0433    2.0467
  250.0600   25.0633    2.0467
  250.0700   25.0633    2.0517
  250.0800   25.0667    2.0483
  250.0800   25.0033    2.0433
  250.0800   25.1067    2.0483
  250.1150   25.0567    2.0750
  250.1350   25.0667    2.0117];
geom=geom./1e3;
folder_path = 'C:\Users\rkerk\OneDrive\Dokument\GitHub\DACS\ASM109_2021_data\MTS\Quasi\'; % Replace with your folder path
file_list = dir([folder_path '*.csv']); % Get a list of all CSV files in the folder

table_list = cell(length(file_list), 1); % Initialize a cell array to store the table data

for i = 1:length(file_list)
    filename = [folder_path file_list(i).name]; % Get the full file path
    table_data = readtable(filename); % Read the CSV file using readtable
    table_list{i} = table_data; % Store the table data in the cell array
    % Do something with the table data
end
samples=[1,2,3,5,6,8,9,10];
Areas=[];
for i= 1:length(samples)
    Areas=[Areas, geom(i, 2).*geom(i, 3)]; %add only the areas of the samples used
end

%% calculating Xt 
%plot(1:1:length(table_listUD{1}{:, 3}), table_listUD{1}{:, 3})
Xttable= [];
for i =1:length(samples)
    maxload=max(table_list{i}{:, 3});%retrieves the maximum load
    Xttable=[Xttable, maxload./Areas(i)];
end
meanXt=mean(Xttable)
maxloadtot=mean(maxload)
stdXt=std(Xttable);
stdmaxload=std(maxload)


x=meanXt-3*stdXt:1e5:meanXt+3*stdXt;
plot(x./1e6, normpdf(x, meanXt, stdXt));xlabel("E [MPa]"); grid on;
