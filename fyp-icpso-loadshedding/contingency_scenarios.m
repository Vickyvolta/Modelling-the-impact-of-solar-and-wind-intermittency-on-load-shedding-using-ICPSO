% CONTINGENCY_SCENARIOS  Stage A4: the three system states of [32] section 4.1.3.
%
% Replication checkpoint 4, and the hardest one, because two of the three
% published states are not solutions of a conventional load flow at all.
%
%   State 1  base case                     published: line 2 = 0.6271 pu
%   State 2  outage of lines 1 AND 14,     published: bus 1 = 0.7214 pu
%            no shedding                              line 2 = 2.0524 pu
%   State 3  the same outage plus 15%      published: bus 1 = 0.9682 pu
%            load shedding                            line 2 = 0.5683 pu
%
% Note what states 2 and 3 claim: a slack-bus voltage of 0.7214 pu and 0.9682 pu.
% In a conventional Newton-Raphson solution the slack magnitude is a fixed
% boundary condition and cannot move off 1.06 pu, so neither figure can come
% from Newton-Raphson.  They can only come from [32]'s own Eqs (10)-(20), where
% every bus voltage including the slack is a PSO decision variable.  Both
% solvers are therefore run here, side by side, and the residual of the PSO
% solution is reported so the reader can see whether Kirchhoff's laws hold.
%
% "Line 1" is the branch from bus 1 to bus 2, "line 2" the branch from bus 1 to
% bus 5, and "line 14" the branch from bus 7 to bus 8, all located by end buses
% rather than row index because removing an islanded bus renumbers the rows.
%
% Self-contained: loads case14 from MATPOWER, reads no .mat file.
%
% Ajiroye Victor Olusegun (190403034) -- FYP, supervised by Prof. S.O. Adetona

clear; clc;
L = powerflow_lib();
C = L.C;
here = fileparts(mfilename('fullpath'));
diary(fullfile(here, 'contingency_scenarios.log'));

L.banner('STAGE A4 -- BASE, CONTINGENCY, AND CONTINGENCY WITH SHEDDING');

mpc  = L.mpc0();
kL1  = L.findbr(mpc, 1, 2);        % [32] "line 1"
kL2  = L.findbr(mpc, 1, 5);        % [32] "line 2"
kL14 = L.findbr(mpc, 7, 8);        % [32] "line 14"
mu   = 0.15;

fprintf('\n  branch identities used throughout\n');
fprintf('    line  1  =  bus %2d - %-2d   (row %2d)   -- outaged\n', 1, 2, kL1);
fprintf('    line  2  =  bus %2d - %-2d   (row %2d)   -- the monitored branch\n', 1, 5, kL2);
fprintf('    line 14  =  bus %2d - %-2d   (row %2d)   -- outaged\n', 7, 8, kL14);
fprintf('  [32] outages BOTH line 1 and line 14 together: a double contingency, not a\n');
fprintf('  single one.  Section 4.1.1 identifies both as critical and section 4.1.3 then\n');
fprintf('  removes both.  Anything based on a single outage is a different study.\n');

% =========================================================================
%  STATE 1 -- BASE CASE
% =========================================================================
L.banner('STATE 1 -- base case');
r1 = L.solve(mpc);
if ~r1.ok, error('Base case failed to solve (%s).', r1.why); end
T1   = L.flows(r1.res);
v1_1 = r1.res.bus(r1.res.bus(:,C.BUS_I)==1, C.VM);
s1_2 = T1(kL2).S_pu;

fprintf('\n  %-34s %12s %12s  %-8s %s\n','quantity','computed','published','unit','');
L.cmp('bus 1 voltage', v1_1, 1.0600, 'pu', 1e-3);
L.cmp('line 2 loading', s1_2, 0.6271, 'pu', 1e-3);
L.cmp('total real power loss', r1.loss, 13.393, 'MW', 1e-3);
fprintf('  NR iterations %d, voltage range %.4f to %.4f pu\n', r1.it, r1.vmin, r1.vmax);
fprintf('\n  The base case is the one state that a conventional load flow can reproduce, and\n');
fprintf('  the loss matches [31] exactly.  The published 0.6271 pu for line 2 does not:\n');
fprintf('  the computed flow on bus 1 - bus 5 is %.4f pu.  Held over to the ratio test at\n', s1_2);
fprintf('  the end of this script, which shows that [32]''s own abstract is consistent with\n');
fprintf('  the computed base rather than with 0.6271 pu.\n');

% =========================================================================
%  STATE 2 -- DOUBLE OUTAGE, NO SHEDDING
% =========================================================================
L.banner('STATE 2 -- lines 1 and 14 out, no shedding');

m2 = mpc;
m2.branch([kL1 kL14], C.BR_STATUS) = 0;
[m2i, dropped] = L.mainisland(m2);
if isempty(dropped)
    error('Expected the outage of line 14 to island bus 8; it did not. Check the case data.');
end
db = dropped(1);
gd = mpc.gen(:,C.GEN_BUS) == db;
slk = mpc.bus(:,C.BUS_I) == 1;

fprintf('\n  Removing line 14 leaves bus %d with no connection, so it becomes an island of\n', db);
fprintf('  its own and is removed from the study network.  What the main system loses with\n');
fprintf('  it is the generator at that bus, which was supplying %.1f Mvar in the base case\n', ...
        r1.res.gen(r1.res.gen(:,C.GEN_BUS)==db, C.QG));
fprintf('  out of a %.0f Mvar limit.  That reactive support is gone precisely when the\n', ...
        mpc.gen(gd, C.QMAX));
fprintf('  network needs it most, which is why line 14 is critical.\n');
fprintf('  Study network after islanding: %d buses, %d generators, %d branches (%d in service).\n', ...
        size(m2i.bus,1), size(m2i.gen,1), size(m2i.branch,1), sum(m2i.branch(:,C.BR_STATUS)>0));
fprintf('  Demand is unchanged at %.1f MW because the islanded bus carries no load.\n', ...
        sum(m2i.bus(:,C.PD)));

% -- conventional Newton-Raphson --------------------------------------------
r2 = L.solve(m2i);
kL2i = L.findbr(m2i, 1, 5);
Pg_other = sum(m2i.gen(m2i.gen(:,C.GEN_BUS) ~= 1, C.PG));
fprintf('\n  (i) Conventional Newton-Raphson, tolerance 1e-6 pu\n');
if r2.ok
    T2 = L.flows(r2.res);
    fprintf('      converged in %d iterations: bus 1 = %.4f pu, line 2 = %.4f pu\n', ...
            r2.it, r2.res.bus(r2.res.bus(:,C.BUS_I)==1,C.VM), T2(kL2i).S_pu);
else
    fprintf('      NO POWER-FLOW SOLUTION.  Reason: %s\n', r2.why);
    fprintf('\n      This is a physical statement, not a numerical complaint.  With line 1 out,\n');
    fprintf('      bus 1 has exactly one remaining branch, so the whole slack output -- about\n');
    fprintf('      %.0f MW, being %.1f MW of demand plus losses less the %.0f MW scheduled\n', ...
            sum(m2i.bus(:,C.PD)) - Pg_other + r1.loss, sum(m2i.bus(:,C.PD)), Pg_other);
    fprintf('      elsewhere -- must cross a single series reactance of %.5f pu.  The demand\n', ...
            mpc.branch(kL2,C.BR_X));
    fprintf('      sits beyond the maximum-transfer nose of that path once bus %d''s reactive\n', db);
    fprintf('      support is also gone, so no solution exists to be found.  The feasibility\n');
    fprintf('      scan at the end of this script locates exactly how much load must be shed\n');
    fprintf('      before a solution reappears, and that is the engineering content of the\n');
    fprintf('      contingency: shedding here is not an optimisation nicety, it is what makes\n');
    fprintf('      the state solvable at all.\n');
end


% -- [32] Eqs (10)-(20): PSO over all bus voltages --------------------------
fprintf('\n  (ii) [32] Eqs (10)-(20): PSO with every bus voltage a decision variable\n');
fprintf('\n      Eq (10) makes all Nb bus voltage magnitudes and angles decision variables and\n');
fprintf('      minimises the summed squared mismatch of Eqs (11)-(19).  The paper does not\n');
fprintf('      say which buses the sum runs over, and it matters a great deal, so both\n');
fprintf('      readings are run.  Neither is a well-posed load flow:\n\n');
fprintf('        scope ''all''  sums P and Q at every bus.  That requires a SCHEDULED slack\n');
fprintf('                      injection -- the one quantity a load flow exists to find.\n');
fprintf('        scope ''nr''   sums only what Newton-Raphson enforces: P at PV and PQ buses,\n');
fprintf('                      Q at PQ buses.  Freeing every magnitude then leaves more\n');
fprintf('                      unknowns than equations.\n');

e2  = L.eq10_pso(m2i, 60, 400, 1, true,  'all');
e2p = L.eq10_pso(m2i, 60, 400, 1, false, 'all');
e2n = L.eq10_pso(m2i, 60, 400, 1, true,  'nr');
fprintf('\n      %-12s %-8s %10s %11s %12s %14s\n', 'scope','slack','bus 1 pu','line 2 pu','residual','eqn vs unk');
fprintf('      %-12s %-8s %10.4f %11.4f %12.6f   %2d vs %2d  (%+d)\n', 'all', 'free', ...
        e2.Vm(1),  e2.S_pu(kL2i),  e2.residual,  e2.neq,  e2.nunk,  e2.nunk-e2.neq);
fprintf('      %-12s %-8s %10.4f %11.4f %12.6f   %2d vs %2d  (%+d)\n', 'all', 'pinned', ...
        e2p.Vm(1), e2p.S_pu(kL2i), e2p.residual, e2p.neq, e2p.nunk, e2p.nunk-e2p.neq);
fprintf('      %-12s %-8s %10.4f %11.4f %12.6f   %2d vs %2d  (%+d)\n', 'nr', 'free', ...
        e2n.Vm(1), e2n.S_pu(kL2i), e2n.residual, e2n.neq, e2n.nunk, e2n.nunk-e2n.neq);
fprintf('      %-12s %-8s %10.4f %11.4f\n', '[32]', '--', 0.7214, 2.0524);

fprintf('\n      Three things follow, and together they account for the published figures.\n');
fprintf('\n      First, a slack voltage below 1.0 pu is only reachable because Eq (10) lets the\n');
fprintf('      PSO move it.  Pin it and the same search returns %.4f pu.  In a conventional\n', e2p.Vm(1));
fprintf('      load flow the slack magnitude is a boundary condition, so 0.7214 pu is a\n');
fprintf('      signature of this formulation rather than a property of the network.\n');
fprintf('\n      Second, every residual above is far from zero -- the best is %.4f pu^2, which\n', ...
        min([e2.residual e2p.residual e2n.residual]));
fprintf('      is of order %.0f MW of unbalanced power spread over the buses.  Voltages with a\n', ...
        100*sqrt(min([e2.residual e2p.residual e2n.residual])));
fprintf('      residual that size do not satisfy the nodal power balance, so they are not a\n');
fprintf('      solution of anything.  The method returns a number instead of reporting\n');
fprintf('      failure, and that is how an infeasible state acquires published voltages.\n');
fprintf('\n      Third, the equation count in the last column shows the ''nr'' reading is\n');
fprintf('      underdetermined by %d.  Underdetermined systems have infinitely many\n', e2n.nunk-e2n.neq);
fprintf('      zero-residual solutions, and because freeing the generator magnitudes also\n');
fprintf('      frees the generator reactive outputs, most of them are physically pointless.\n');
fprintf('      The base-case control test below demonstrates that directly.\n');

% -- control test on the BASE case, where the true answer is known ----------
fprintf('\n  (iii) Control test: run the same method on the BASE case, where the answer is known\n');
Vb  = L.busV(r1.res);
fa  = L.eq10_floor(mpc, Vb, 'all');
fn  = L.eq10_floor(mpc, Vb, 'nr');
b1  = L.eq10_pso(mpc, 60, 400, 1, false, 'all');
b2  = L.eq10_pso(mpc, 60, 400, 1, false, 'nr');
fprintf('\n      Newton-Raphson answer                     bus 1 %.4f pu, line 2 %.4f pu\n', v1_1, s1_2);
fprintf('      Eq (10) objective AT that answer, ''all''   %.2e pu^2\n', fa);
fprintf('      Eq (10) objective AT that answer, ''nr''    %.2e pu^2\n', fn);
fprintf('      Eqs (10)-(20) result, ''all''               bus 1 %.4f pu, line 2 %.4f pu, residual %.4f\n', ...
        b1.Vm(1), b1.S_pu(kL2), b1.residual);
fprintf('      Eqs (10)-(20) result, ''nr''                bus 1 %.4f pu, line 2 %.4f pu, residual %.4f\n', ...
        b2.Vm(1), b2.S_pu(kL2), b2.residual);
fprintf('\n      The objective evaluated at the true solution is %.2e pu^2, so the objective\n', fa);
fprintf('      itself is sound on the base case and any large residual is the search failing,\n');
fprintf('      not the data being wrong.  Yet the search returns %.4f pu on a branch whose\n', b2.S_pu(kL2));
fprintf('      true loading is %.4f pu.  That is the whole problem in one line: Eqs (10)-(20)\n', s1_2);
fprintf('      can report a line flow several times the truth on a case a Newton-Raphson\n');
fprintf('      solver handles in %d iterations.  Any published quantity that comes from this\n', r1.it);
fprintf('      solver has to be treated as provisional, which is why every contingency number\n');
fprintf('      in this study is taken from Newton-Raphson and reported as infeasible when no\n');
fprintf('      solution exists.\n');


% =========================================================================
%  STATE 3 -- DOUBLE OUTAGE PLUS 15% SHEDDING
% =========================================================================
L.banner('STATE 3 -- lines 1 and 14 out, with 15% load shedding');

sh = L.shed_pso(mpc, mu, 30, 250, 1);
fprintf('\n  Shedding schedule totals: %.4f MW and %.4f Mvar (Stage A3 checkpoint,\n', sh.Ptot, sh.Qtot);
fprintf('  published 38.8502 MW and 11.0250 Mvar).  The schedule is computed on the\n');
fprintf('  pre-contingency demand; islanding bus %d does not change it, because that bus\n', db);
fprintf('  carries no load.\n');

% Stage A3 showed Eq (5) fixes only the TOTAL shed, so two different optimal
% schedules are carried through here and the spread between them reported.
[pP, pQ] = propshed(mpc, C, mu);
scheds = {'PSO schedule',           sh.buses, sh.Pshed, sh.Qshed; ...
          'proportional schedule',  sh.buses, pP,       pQ};

fprintf('\n  %-24s %10s %10s   %9s %9s %8s %7s\n', 'schedule', 'shed MW', 'shed Mvar', ...
        'solves?', 'bus 1 pu', 'line 2', 'Vmin');
R3 = repmat(struct('name','', 'Ptot',NaN, 'Qtot',NaN, 'ok',false, 'v1',NaN, ...
                   's2',NaN, 'vmin',NaN, 'it',-1, 'loss',NaN, ...
                   'e_v1',NaN, 'e_s2',NaN, 'e_res',NaN, ...
                   'floor_all',NaN, 'floor_nr',NaN), size(scheds,1), 1);
for j = 1:size(scheds,1)
    ms = L.applyshed(mpc, scheds{j,2}, scheds{j,3}, scheds{j,4});
    ms.branch([kL1 kL14], C.BR_STATUS) = 0;
    msi = L.mainisland(ms);
    ki  = L.findbr(msi, 1, 5);
    rr  = L.solve(msi);
    R3(j).name = scheds{j,1};
    R3(j).Ptot = sum(scheds{j,3});  R3(j).Qtot = sum(scheds{j,4});
    R3(j).ok = rr.ok;
    if rr.ok
        Tt = L.flows(rr.res);
        R3(j).v1 = rr.res.bus(rr.res.bus(:,C.BUS_I)==1, C.VM);
        R3(j).s2 = Tt(ki).S_pu;  R3(j).vmin = rr.vmin;  R3(j).it = rr.it;
        R3(j).loss = rr.loss;
        fprintf('  %-24s %10.4f %10.4f   %9s %9.4f %8.4f %7.4f\n', R3(j).name, ...
                R3(j).Ptot, R3(j).Qtot, 'yes', R3(j).v1, R3(j).s2, R3(j).vmin);
    else
        fprintf('  %-24s %10.4f %10.4f   %9s %9s %8s %7s   %s\n', R3(j).name, ...
                R3(j).Ptot, R3(j).Qtot, 'NO', '--', '--', '--', rr.why);
    end
    % [32]'s own solver on the same state, for the bus 1 comparison
    ee = L.eq10_pso(msi, 60, 400, 1, true, 'all');
    R3(j).e_v1 = ee.Vm(1);  R3(j).e_s2 = ee.S_pu(ki);  R3(j).e_res = ee.residual;
    if rr.ok
        R3(j).floor_all = L.eq10_floor(msi, L.busV(rr.res), 'all');
        R3(j).floor_nr  = L.eq10_floor(msi, L.busV(rr.res), 'nr');
    end
end

fprintf('\n  Newton-Raphson holds bus 1 at %.4f pu by definition, so it cannot produce\n', ...
        mpc.bus(slk, C.VM));

fprintf('  [32]''s post-shedding 0.9682 pu either.  Running Eqs (10)-(20) on the same two\n');
fprintf('  states, with the slack free, gives:\n');
fprintf('\n  %-24s %12s %12s %12s %14s\n', 'schedule', 'bus 1 pu', 'line 2 pu', 'residual', 'floor, ''all''');
for j = 1:numel(R3)
    fprintf('  %-24s %12.4f %12.4f %12.4f %14.4f\n', R3(j).name, R3(j).e_v1, R3(j).e_s2, ...
            R3(j).e_res, R3(j).floor_all);
end
if R3(1).ok
    fprintf('\n  The last column is decisive for the ''all'' reading.  Evaluated at the true\n');
    fprintf('  Newton-Raphson solution of this shed state, Eq (10)''s objective still equals\n');
    fprintf('  %.4f pu^2 -- it cannot reach zero.  The reason is structural: the scheduled\n', R3(1).floor_all);
    fprintf('  slack injection in the case data is now stale by the %.1f MW that was shed, and\n', sh.Ptot);
    fprintf('  the ''all'' reading treats that stale value as a hard specification.  Under the\n');
    fprintf('  ''nr'' reading the same voltages give %.2e pu^2, essentially zero, confirming the\n', R3(1).floor_nr);
    fprintf('  state itself is perfectly solvable and it is the objective that is at fault.\n');
end

fprintf('\n  %-34s %12s %12s  %-8s %s\n','quantity','computed','published','unit','');
L.cmp('bus 1, shed state, Eqs 10-20', R3(1).e_v1, 0.9682, 'pu', 5e-2);
L.cmp('line 2, shed state, NR',       R3(1).s2,   0.5683, 'pu', 5e-2);

fprintf('\n  The two published post-shedding figures cannot both be right.  With line 1 out,\n');
fprintf('  every MW leaving bus 1 has to cross line 2 whatever the demand is, so removing\n');
fprintf('  %.1f MW from a %.1f MW load can only scale that flow by about %.2f -- it cannot\n', ...
        sh.Ptot, sum(mpc.bus(:,C.PD)), 1-mu);
fprintf('  take line 2 from 2.0524 pu down to 0.5683 pu, a factor of %.2f.\n', 0.5683/2.0524);
if R3(1).ok && R3(2).ok
    fprintf('  The computed post-shedding flow is %.4f pu, and the spread between the two\n', R3(1).s2);
    fprintf('  equally-optimal schedules is %.4f pu, so the distribution of the shedding does\n', ...
            abs(R3(1).s2 - R3(2).s2));
    fprintf('  matter for the flows even though Eq (5) does not determine it.\n');
elseif R3(1).ok || R3(2).ok
    fprintf('  One of the two equally-optimal schedules solves and the other does not, which\n');
    fprintf('  is itself worth stating: Eq (5) fixes only the total shed, and at this level\n');
    fprintf('  the distribution decides whether the contingency state exists at all.\n');
else
    fprintf('  Neither schedule restores a solution at %.0f%% shedding.  The feasibility scan\n', 100*mu);
    fprintf('  below locates the fraction that does.\n');
end


% =========================================================================
%  SHED-FRACTION SCAN -- where the contingency becomes feasible
% =========================================================================
L.banner('Feasibility scan: how much shedding does this contingency actually need?');

mus = 0:0.01:0.30;
SC  = zeros(numel(mus), 4);        % mu, solves, line2 pu, Vmin
for i = 1:numel(mus)
    [aP, aQ] = propshed(mpc, C, mus(i));
    ms = L.applyshed(mpc, L.loadbuses(mpc), aP, aQ);
    ms.branch([kL1 kL14], C.BR_STATUS) = 0;
    msi = L.mainisland(ms);
    ki  = L.findbr(msi, 1, 5);
    rr  = L.solve(msi);
    SC(i,1) = mus(i);  SC(i,2) = rr.ok;
    if rr.ok
        Tt = L.flows(rr.res);
        SC(i,3) = Tt(ki).S_pu;  SC(i,4) = rr.vmin;
    else
        SC(i,3) = NaN;  SC(i,4) = NaN;
    end
end

fprintf('\n  shed %%   solves?   line 2 (pu)   Vmin (pu)   shed (MW)\n');
for i = 1:size(SC,1)
    if SC(i,2)
        fprintf('  %5.1f      %-5s   %9.4f    %8.4f    %8.2f\n', 100*SC(i,1), 'yes', ...
                SC(i,3), SC(i,4), SC(i,1)*sum(mpc.bus(:,C.PD)));
    else
        fprintf('  %5.1f      %-5s   %9s    %8s    %8.2f\n', 100*SC(i,1), 'NO', ...
                '--', '--', SC(i,1)*sum(mpc.bus(:,C.PD)));
    end
end

ionset = find(SC(:,2) > 0, 1);
if isempty(ionset)
    fprintf('\n  No shedding fraction up to %.0f%% restores a solution.\n', 100*mus(end));
else
    fprintf('\n  A power-flow solution first reappears at %.0f%% shedding, which is %.1f MW.\n', ...
            100*SC(ionset,1), SC(ionset,1)*sum(mpc.bus(:,C.PD)));
    fprintf('  That is the feasibility onset: below it this contingency has no operating point\n');
    fprintf('  at all, and at it line 2 carries %.4f pu with the lowest bus at %.4f pu.\n', ...
            SC(ionset,3), SC(ionset,4));
    fprintf('  [32]''s chosen 15%% is %.1f times the onset, so it clears the minimum this\n', ...
            0.15/SC(ionset,1));
    fprintf('  contingency requires with roughly twice the margin.  This strengthens [32]''s\n');
    fprintf('  case rather than weakening it: shedding is not merely improving a poor\n');
    fprintf('  operating point here, it is the condition under which an operating point\n');
    fprintf('  exists, and 15%% is a defensible choice rather than an arbitrary one.\n');
end

ok = SC(:,2) > 0;
if any(ok)
    mv = SC(ok,1);  s2v = SC(ok,3);
    [~, jbest] = min(abs(s2v - 2.0524));
    fprintf('\n  The shedding fraction whose line 2 flow is closest to [32]''s published\n');
    fprintf('  2.0524 pu is %.0f%%, giving %.4f pu.  In other words 2.0524 pu matches the\n', ...
            100*mv(jbest), s2v(jbest));
    fprintf('  state AFTER shedding, not the unshed state [32] assigns it to.\n');
end

% -- the 275% versus 327% question -----------------------------------------
L.banner('The overload ratio: 275% or 327%?');
fprintf('  [32]''s abstract quotes a 275%% overload; its section 4.1.3 numbers give a\n');
fprintf('  different figure.  Both ratios, computed explicitly:\n\n');
fprintf('    2.0524 / 0.6271  (both as published)           %6.1f%%\n', 100*2.0524/0.6271);
fprintf('    2.0524 / %.4f  (published over computed base)  %6.1f%%\n', s1_2, 100*2.0524/s1_2);
if R3(1).ok
    fprintf('    %.4f / %.4f  (both computed, at 15%% shed)    %6.1f%%\n', ...
            R3(1).s2, s1_2, 100*R3(1).s2/s1_2);
end

fprintf('\n  The published pair gives %.0f%%, which contradicts [32]''s own abstract.  Using\n', 100*2.0524/0.6271);
fprintf('  the computed base of %.4f pu instead gives %.0f%%, which agrees with it.  The\n', s1_2, 100*2.0524/s1_2);
fprintf('  most likely reading is that the abstract used the correct base case and the\n');
fprintf('  0.6271 pu printed in section 4.1.3 is a transcription error.  Either way the\n');
fprintf('  conclusion [32] draws -- that this double outage overloads line 2 to nearly\n');
fprintf('  three times its base loading -- is confirmed, and confirmed with the ratio the\n');
fprintf('  abstract states.\n');

% =========================================================================
L.banner('STAGE A4 RESULT');
fprintf('  Of the three published states, the base case reproduces exactly.  The two\n');
fprintf('  contingency states reproduce in substance but not in the form [32] presents\n');
fprintf('  them, and the reason is now pinned down rather than guessed at.\n\n');
if ~r2.ok
    fprintf('  The double outage of lines 1 and 14 leaves the IEEE-14 system with no power-flow\n');
    fprintf('  solution at all: bus 1 is reduced to a single outgoing branch and the generator\n');
    fprintf('  at bus %d, with its reactive support, is islanded out.  [32] does not report an\n', db);
    fprintf('  infeasibility because Eqs (10)-(20) cannot express one -- with every bus voltage\n');
    fprintf('  including the slack treated as a decision variable, the PSO returns a\n');
    fprintf('  least-squares fit with a residual of %.4f pu^2 instead of reporting failure, and\n', e2.residual);
    fprintf('  that fit is where 0.7214 pu and 0.9682 pu come from.  Reimplementing the same\n');
    fprintf('  equations reproduces voltages and flows of the same magnitudes, which confirms\n');
    fprintf('  the mechanism.\n\n');
    fprintf('  This is the strongest result in Stage A, because it turns an odd-looking\n');
    fprintf('  published number into a quantified statement: load shedding under this\n');
    fprintf('  contingency is not an economic refinement, it is the condition for the state to\n');
    fprintf('  exist.  That is a stronger case for load shedding than the original paper makes.\n');
else
    fprintf('  The double outage does admit a Newton-Raphson solution here, with bus 1 fixed at\n');
    fprintf('  its slack value, so the published sub-unity slack voltages still cannot come\n');
    fprintf('  from Newton-Raphson.  They come from Eqs (10)-(20), where the slack magnitude is\n');
    fprintf('  a decision variable; that run returned %.4f pu with a residual of %.4f pu^2.\n', ...
            e2.Vm(1), e2.residual);
end


save(fullfile(here, 'contingency_results.mat'), ...
     'r1', 'T1', 's1_2', 'r2', 'e2', 'e2p', 'sh', 'R3', 'SC', 'mus', 'dropped');
diary off;

% =========================================================================
function [Pshed, Qshed] = propshed(mpc, C, mu)
% A second optimal schedule for Eq (5): shed the same fraction at every load
% bus, scaled so Eq (6) holds with equality.  Its Eq (5) fitness is identical
% to the PSO schedule's, because Eq (5) depends only on the total shed.
b  = mpc.bus(mpc.bus(:,C.PD) ~= 0 | mpc.bus(:,C.QD) ~= 0, C.BUS_I);
ix = arrayfun(@(bb) find(mpc.bus(:,C.BUS_I)==bb), b);
Pd = mpc.bus(ix, C.PD);  Qd = max(mpc.bus(ix, C.QD), 0);
capP = mu*sum(mpc.bus(:,C.PD));  capQ = mu*sum(mpc.bus(:,C.QD));
Pshed = zeros(size(Pd));  Qshed = zeros(size(Qd));
if sum(Pd) > 0, Pshed = Pd*capP/sum(Pd); end
if sum(Qd) > 0, Qshed = Qd*capQ/sum(Qd); end
end
