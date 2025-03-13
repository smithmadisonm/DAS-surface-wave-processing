% batch write of netCDF files from processed RBR .mat results

clear all, close all

flist = dir('*.mat');

for fi=1:length(flist)

    load(flist(fi).name)
    variable = who('RBR*');
    outputfilename = [flist(fi).name(1:end-3) '.nc']

    RBRresults2NC(eval(variable{1}), outputfilename)

end