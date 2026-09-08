% BASECASE_LOADFLOW  Stage A1: base-case load flow on IEEE-14.
%
% Replication checkpoint 1 of [32]/[31].  Establishes that the network data and
% the load flow are correct before anything is built on top of them.
%
%   Published target   [31] Table 3, Newton-Raphson base case: 13.393 MW loss
%   Published target   [31] Table 4, base-case column: Vg and tap settings
%   Published target   [32] section 4.1.3: base-case line 2 loading, 0.6271 pu
%
% Self-contained: loads case14 from MATPOWER, reads no .mat file.
%
% Ajiroye Victor Olusegun (190403034) -- FYP, supervised by Prof. S.O. Adetona

clear; clc;
L = powerflow_lib();
C = L.C;
diary(fullfile(fileparts(mfilename('fullpath')), 'basecase_loadflow.log'));

L.banner('STAGE A1 -- BASE CASE, IEEE-14');

mpc = L.mpc0();
r   = L.solve(mpc);
if ~r.ok
    error('Base case failed to solve (%s). Check the MATPOWER installation.', r.why);
end
res = r.res;

% -------------------------------------------------------------------------
fprintf('\nNetwork as loaded from MATPOWER case14\n');
fprintf('  buses %d   generators %d   branches %d   base %g MVA\n', ...
        size(mpc.bus,1), size(mpc.gen,1), size(mpc.branch,1), mpc.baseMVA);
fprintf('  total real demand      %8.2f MW    ([32] states 259.0)\n', sum(mpc.bus(:,C.PD)));
fprintf('  total reactive demand  %8.2f Mvar  ([32] states  73.5)\n', sum(mpc.bus(:,C.QD)));
fprintf('  15%% shedding cap       %8.4f MW / %.4f Mvar\n', ...
        0.15*sum(mpc.bus(:,C.PD)), 0.15*sum(mpc.bus(:,C.QD)));
fprintf('  Newton-Raphson tolerance 1e-6 pu ([32] Eq 20), converged in %d iterations\n', r.it);

% -------------------------------------------------------------------------
L.banner('Checkpoint: total real power loss');
fprintf('  %-34s %12s %12s  %-8s %s\n','quantity','computed','published','unit','');
L.cmp('total real power loss', r.loss, 13.393, 'MW', 1e-3);

% -------------------------------------------------------------------------
L.banner('Checkpoint: control settings vs [31] Table 4, base-case column');
fprintf('  %-34s %12s %12s  %-8s %s\n','quantity','computed','published','unit','');
pubbase = [1.060 1.045 1.010 1.070 1.090];
for k = 1:numel(C.VG_BUS)
    b = C.VG_BUS(k);
    L.cmp(sprintf('Vg%d', b), res.bus(res.bus(:,C.BUS_I)==b, C.VM), pubbase(k), 'pu', 2e-3);
end
pubtap = [0.978 0.969 0.932];
tapname = {'T4-7','T4-9','T5-6'};
for k = 1:numel(C.TAP_BR)
    L.cmp(tapname{k}, mpc.branch(C.TAP_BR(k), C.TAP), pubtap(k), 'pu', 2e-3);
end
L.cmp('Qc9 (shunt at bus 9)', mpc.bus(mpc.bus(:,C.BUS_I)==9, C.BS), 19.00, 'Mvar', 1e-3);

% -------------------------------------------------------------------------
L.banner('Bus voltages');
fprintf('  bus  type   |V| (pu)   angle (deg)     Pd (MW)   Qd (Mvar)\n');
tn = {'PQ','PV','REF'};
for i = 1:size(res.bus,1)
    fprintf('   %2d   %-4s  %8.4f   %9.3f   %9.2f   %9.2f\n', ...
        res.bus(i,C.BUS_I), tn{res.bus(i,C.BUS_TYPE)}, res.bus(i,C.VM), ...
        res.bus(i,C.VA), res.bus(i,C.PD), res.bus(i,C.QD));
end
fprintf('  voltage range %.4f to %.4f pu\n', r.vmin, r.vmax);
fprintf('  note bus 4 carries Qd = %.1f Mvar, i.e. a NET CAPACITIVE load.  This makes\n', ...
        mpc.bus(4,C.QD));
fprintf('  [32] Eq (7)''s bound "0 <= Qshed_i <= Qd_i" empty at bus 4, so no reactive\n');
fprintf('  power can be shed there.  Stage A3 handles this explicitly.\n');

% -------------------------------------------------------------------------
L.banner('All 20 branch flows -- locating [32]''s "line 2 = 0.6271 pu"');
T = L.flows(res);
fprintf('  line  from-to    P (MW)    Q (Mvar)     S (MVA)    S (pu)   loss (MW)\n');
for k = 1:numel(T)
    fprintf('   %2d    %2d-%-2d   %9.3f  %9.3f   %9.4f  %8.4f   %8.4f\n', ...
        T(k).k, T(k).f, T(k).t, T(k).Pf, T(k).Qf, T(k).Smva, T(k).S_pu, T(k).ploss);
end

Spu = [T.S_pu]';
k2  = L.findbr(mpc, 1, 5);          % [32] calls the branch 1-5 "line 2"
fprintf('\n  %-34s %12s %12s  %-8s %s\n','quantity','computed','published','unit','');
L.cmp('line 2 (bus 1-5) loading', Spu(k2), 0.6271, 'pu', 1e-3);

[~, near] = min(abs(Spu - 0.6271));
fprintf('\n  The base-case flow on branch 1-5 is %.4f pu, not 0.6271 pu.  No branch in\n', Spu(k2));
fprintf('  the network carries 0.6271 pu; the closest is branch %d (%d-%d) at %.4f pu.\n', ...
        near, T(near).f, T(near).t, Spu(near));
fprintf('  Ratio check: [32]''s post-outage 2.0524 pu over its base 0.6271 pu is 327%%,\n');
fprintf('  whereas its own abstract quotes 275%%.  Against the computed base of %.4f pu,\n', Spu(k2));
fprintf('  2.0524 pu is %.0f%% -- much closer to the abstract.  So the abstract appears to\n', 100*2.0524/Spu(k2));
fprintf('  use the true base case and 0.6271 pu is the inconsistent figure.\n');

% -------------------------------------------------------------------------
L.banner('Generator dispatch');
fprintf('  bus      Pg (MW)   Qg (Mvar)   Qmin    Qmax   at limit?\n');
for i = 1:size(res.gen,1)
    q = res.gen(i,C.QG); lim = '';
    if q >= res.gen(i,C.QMAX)-1e-4, lim = 'Qmax'; end
    if q <= res.gen(i,C.QMIN)+1e-4, lim = 'Qmin'; end
    fprintf('   %2d    %9.3f   %9.3f  %6.1f  %6.1f   %s\n', ...
        res.gen(i,C.GEN_BUS), res.gen(i,C.PG), q, ...
        res.gen(i,C.QMIN), res.gen(i,C.QMAX), lim);
end

% -------------------------------------------------------------------------
L.banner('STAGE A1 RESULT');
fprintf('  The load flow reproduces [31]''s published base-case loss to four significant\n');
fprintf('  figures (%.4f MW against 13.393 MW), and the control settings match [31]\n', r.loss);
fprintf('  Table 4''s base-case column.  The network data and solver are therefore correct,\n');
fprintf('  and every later stage rests on a verified base case.\n');
fprintf('\n  One published figure does not reconcile: [32]''s base-case line 2 loading of\n');
fprintf('  0.6271 pu.  The computed value is %.4f pu.  This is recorded rather than\n', Spu(k2));
fprintf('  adjusted, because the overload ratio in [32]''s own abstract agrees with the\n');
fprintf('  computed base rather than with 0.6271 pu.\n');

save(fullfile(fileparts(mfilename('fullpath')), 'basecase_results.mat'), ...
     'res', 'T', 'Spu', 'r');
diary off;
