params = {{'L'}    {'R'}    {'h'}    {'rho'}    {'E'}    {'v'}    {'RW1'}    {'RW2'}    {'Cwk'}};
params = convertCharsToStrings([params{:}]);

T = table('Size', [3, length(params)],'VariableNames', params, 'VariableTypes', repmat("double", 1, length(params)));

for i = 1:length(params)
    [Optimised_Value, Percentage_Error] = ccm.optimiseParamtot(params(i));
    T.(params(i))(1) = ccm.(params(i));
    T.(params(i))(2) = Optimised_Value;
    T.(params(i))(3) = round(Percentage_Error, 2);
end

writetable(T,'Data/CCA_GlobalOpt.csv');