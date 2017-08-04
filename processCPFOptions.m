% Do some pre-processing based on CPFOptions
%
%
function CPFOptions = processCPFOptions(CPFOptions)

define_constants;

% function handle to case file
CPFOptions.caseFunction = str2func(CPFOptions.caseFile);
mpc = CPFOptions.caseFunction();

% count number of contingencies
CPFOptions.nContingencies = 1;
if CPFOptions.tripAllLines
    CPFOptions.nLines = size(mpc.branch,1);
    CPFOptions.tripLines = 1:CPFOptions.nLines;
else
    CPFOptions.nLines = length(CPFOptions.tripLines);
end
CPFOptions.nContingencies = CPFOptions.nContingencies + CPFOptions.nLines;

CPFOptions.nGenerators = length(CPFOptions.tripGenerators);
CPFOptions.nContingencies = CPFOptions.nContingencies + CPFOptions.nGenerators;

CPFOptions.nWindPoints = length(CPFOptions.pWind);

% default is to increase load at all buses
if isempty(CPFOptions.loadIncreaseBuses)
    CPFOptions.loadIncreaseBuses = 1:size(mpc.bus,1);
end
