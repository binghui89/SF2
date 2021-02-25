function add_pvlib()
% Add path and Site ID mapping
if ispc
    addpath('C:\Users\bxl180002\git\MATLAB_PV_LIB');
    addpath('C:\Users\bxl180002\git\MATLAB_PV_LIB\Example Data');
    addpath('C:\Users\bxl180002\git\MATLAB_PV_LIB\Required Data');
elseif isunix
    addpath('/home/bxl180002/git/MATLAB_PV_LIB');
    addpath('/home/bxl180002/git/MATLAB_PV_LIB/Example Data');
    addpath('/home/bxl180002/git/MATLAB_PV_LIB/Required Data');
end
end