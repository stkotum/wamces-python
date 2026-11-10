% MAKE_REFERENCE  Reference results of the published Simulink model for the tests.
%
% Run it in a folder that holds the files of the Zenodo record 10.5281/zenodo.18721887:
%   matlab -batch "make_reference"
% It writes wamces_reference.mat into the folder named by the environment variable
% WAMCES_REF_OUT, or into the current folder. The tests read it from tests/reference/.
% With WAMCES_FULL=1 it also writes every state at every step into states_<case>.mat.
% These files are large and serve a full comparison only.
%
% Two cases are simulated:
%   published  the published script and model, 50 s, traces at every tenth step
%   variant    the wide-area gain Kw = 20, a load step of 1000 MW at bus 1 and the
%              power setpoints of the grid-following converters raised by 0.1 pu, 10 s
% The variant drives the parts of the model that rest in the published case. Bus 1 has
% no source. The step there changes the reduced network in a second way.
%
% The script changes none of the published files. It writes a copy of the published
% script with the edits of each case, runs it and deletes it.

out = getenv('WAMCES_REF_OUT');
if isempty(out), out = pwd; end
full = strcmp(getenv('WAMCES_FULL'), '1');

cases = struct( ...
    'name', {'published', 'variant'}, ...
    'stop', {50, 10}, ...
    'every', {10, 1}, ...
    'snap', {[0 10 100 1000 5000], [0 1 2 3 10 100 1000]}, ...
    'edits', {{}, {'Tgfm.Kw(:) = 0;', 'Tgfm.Kw(:) = 20;'; ...
                   'idis = 2178;', 'idis = 1;'; ...
                   'Tgfl.Q0 = Tgfl.Q;', sprintf('Tgfl.Q0 = Tgfl.Q;\nTgfl.P0 = Tgfl.P0 + 0.1;')}});

REF = struct();
REF.release = version('-release');
for kc = 1:numel(cases)
    c = cases(kc);
    clearvars -except REF out full cases kc c

    % ----- Initialisation as published, with the edits of the case
    txt = fileread('wamces_script.m');
    txt = strrep(txt, 'clearvars;', '');
    txt = strrep(txt, 'clc;', '');
    txt = strrep(txt, 'sim(''wamces_model.mdl'');', '');
    for ke = 1:size(c.edits, 1)
        assert(contains(txt, c.edits{ke, 1}), 'edit not found: %s', c.edits{ke, 1});
        txt = strrep(txt, c.edits{ke, 1}, c.edits{ke, 2});
    end
    fid = fopen('wamces_init_tmp.m', 'w'); fwrite(fid, txt); fclose(fid);
    wamces_init_tmp;
    delete('wamces_init_tmp.m');

    % ----- Simulation with every state at every step
    load_system('wamces_model');
    tic;
    so = sim('wamces_model', 'StopTime', num2str(c.stop), 'SaveState', 'on', ...
        'StateSaveName', 'xout', 'SaveFormat', 'StructureWithTime', 'SaveTime', 'on', ...
        'TimeSaveName', 'tout', 'SaveOutput', 'off', 'LimitDataPoints', 'off', ...
        'Decimation', '1', 'ReturnWorkspaceOutputs', 'on');
    wall = toc;
    x = so.get('xout');
    t = so.get('tout');

    % ----- Snapshots of all states and three columns of each block over time.
    % The traces of the published case keep every tenth step only.
    R = struct();
    R.wall_time = wall;
    R.solver = so.SimulationMetadata.ModelInfo.SolverInfo.Solver;
    R.every = c.every;
    R.t = t(1:c.every:end);
    R.snap_steps = c.snap;
    R.block = {x.signals.blockName};
    for k = 1:numel(x.signals)
        v = x.signals(k).values;
        cols = unique(round(linspace(1, size(v, 2), 3)));
        R.(sprintf('snap%02d', k)) = v(c.snap + 1, :);
        R.(sprintf('cols%02d', k)) = cols;
        R.(sprintf('trace%02d', k)) = v(1:c.every:end, cols);
    end
    REF.(c.name) = R;
    if full
        F = struct('t', t, 'block', {{x.signals.blockName}});
        for k = 1:numel(x.signals)
            F.(sprintf('s%03d', k)) = x.signals(k).values;
        end
        save(fullfile(out, ['states_' c.name '.mat']), '-struct', 'F', '-v7.3');
    end
    close_system('wamces_model', 0);
end
save(fullfile(out, 'wamces_reference.mat'), '-struct', 'REF', '-v7');
