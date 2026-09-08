function run_all()
% run_all  Run the full load-shedding + intermittency suite, in order.
%
% Each stage is self-contained: it loads case14 straight from MATPOWER, depends
% on NO prior .mat file, writes its own <stage>.log via diary, and saves
% its own <stage>_results.mat.  Run in sequence they reproduce the whole
% argument -- base-case load flow, the N-1 critical-line screen, the 15% shed,
% the three contingency states, the ICPSO-vs-PSO study, and the intermittency
% layer that finally gives the optimiser something to do.
%
% Requires MATPOWER on the path (loadcase, runpf) and powerflow_lib.m in this
% folder.  No Optimization or Statistics toolbox is used anywhere.  The full
% protocol runs many thousands of load flows and takes ~30-60 min, dominated by
% Stage C; set QUICK=true at the top of re_intermittency.m for a ~5 min
% smoke test that still exercises every block.
%
% Ajiroye Victor Olusegun (190403034) -- FYP, supervised by Prof. S.O. Adetona

here = fileparts(mfilename('fullpath'));
cd(here);
if exist('loadcase', 'file') ~= 2
    error(['MATPOWER not found on the path (loadcase is missing). ' ...
           'Add MATPOWER with addpath, then rerun RUN_ALL.']);
end

stages = { ...
    'basecase_loadflow'; ...
    'critical_line_screen'; ...
    'load_shedding'; ...
    'contingency_scenarios'; ...
    'icpso_vs_pso'; ...
    're_intermittency'};

fprintf('\n============================================================\n');
fprintf(' LOAD-SHEDDING & ICPSO SIMULATION SUITE -- running %d stages in order\n', numel(stages));
fprintf('============================================================\n');

secs = zeros(numel(stages), 1);
okf  = false(numel(stages), 1);
for i = 1:numel(stages)
    fprintf('\n>>> [%d/%d] %s\n', i, numel(stages), stages{i});
    t = tic;
    okf(i) = runstage(stages{i});
    secs(i) = toc(t);
    if okf(i), tag = 'done'; else, tag = 'FAILED'; end
    fprintf('<<< [%d/%d] %s -- %s in %.1f s\n', i, numel(stages), stages{i}, tag, secs(i));
end

fprintf('\n============================================================\n');
fprintf(' SUMMARY\n');
fprintf('============================================================\n');
for i = 1:numel(stages)
    if okf(i), tag = 'ok    '; else, tag = 'FAILED'; end
    fprintf('  %-26s %s  %8.1f s\n', stages{i}, tag, secs(i));
end
fprintf('  %-26s %s  %8.1f s\n', '(total)', '      ', sum(secs));
if all(okf)
    fprintf('\nAll stages completed.  See each *.log and *_results.mat.\n');
else
    fprintf(2, '\nSome stages FAILED -- see the messages above.\n');
end
end

function ok = runstage(name)
% Run one stage script inside THIS helper's workspace, so the script's
% top-of-file CLEAR cannot reach RUN_ALL's loop variables.  Never aborts the
% batch: a failing stage is reported and the next one still runs.  'ok' is set
% only AFTER run() returns, because the script's CLEAR wipes this workspace.
try
    run(name);
    ok = true;
catch e
    ok = false;
    fprintf(2, '  ERROR in %s: %s\n', name, e.message);
end
end
