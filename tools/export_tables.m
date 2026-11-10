% EXPORT_TABLES  Write the nine tables of WAMCES as plain structs into raw_tables.mat.
%
% The published .mat files hold MATLAB tables. Python cannot read them. This script
% converts each table into a struct of columns. Scipy can read the file raw_tables.mat.
% The script tools/make_data.py then writes the CSV files of the package from it.
%
% Run it in a folder that holds the files of the Zenodo record 10.5281/zenodo.18721887:
%   matlab -batch "export_tables"
% The environment variable WAMCES_REF_OUT names the output folder. The default is the
% current folder.

out = getenv('WAMCES_REF_OUT');
if isempty(out), out = pwd; end
raw = struct();
for f = {'data_bus','data_lin','data_sm','data_avr','data_pss','data_tur','data_gfm','data_gfl','data_pmu'}
    S = load([f{1} '.mat']);
    name = fieldnames(S);
    T = S.(name{1});
    raw.(name{1}) = table2struct(T, 'ToScalar', true);
    raw.([name{1} '_columns']) = T.Properties.VariableNames;
end
save(fullfile(out, 'raw_tables.mat'), '-struct', 'raw', '-v7');
