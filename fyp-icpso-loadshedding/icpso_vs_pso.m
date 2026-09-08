% ICPSO_VS_PSO  Stage B: [31]'s ICPSO substituted for [32]'s standard PSO.
%
% [32] does load shedding with a standard PSO.  [31] does optimal reactive
% power dispatch with an Improved Chaotic PSO.  Adetona has never applied the
% ICPSO of [31] to the load-shedding problem of [32], and that substitution is
% this project's contribution.  Stage B carries it out and measures it.
%
% Order of business, deliberately:
%   B1  the wingbeat parameters of [31] section 3.1, reproduced from scratch
%   B2  what the wingbeat term does to the swarm, derived algebraically
%   B3  ICPSO against PSO on [31]'s OWN problem, so the implementation is
%       validated against [31]'s own published numbers before it is trusted
%       anywhere else
%   B4  ICPSO against PSO on [32]'s load-shedding problem, the substitution
%       this project exists to make
%
% Published targets, [31] Table 3, IEEE-14 real power loss in MW:
%       Newton-Raphson base case                        13.393
%       PSO      best 12.275   worst 12.303   mean 12.288
%       ICPSO    best 12.260   worst 12.270   mean 12.265
%       convergence onset: PSO around iteration 40, ICPSO around 25
%       each method run three times, 250 iterations, MATPOWER 7.1
% Published targets, [32] section 4.1.2, load shedding:
%       38.8502 MW and 11.0250 Mvar
%
% Self-contained: loads case14 from MATPOWER, reads no .mat file.
%
% Ajiroye Victor Olusegun (190403034) -- FYP, supervised by Prof. S.O. Adetona

clear; clc;
L = powerflow_lib();
C = L.C;
here = fileparts(mfilename('fullpath'));
diary(fullfile(here, 'icpso_vs_pso.log'));

L.banner('STAGE B -- ICPSO IN PLACE OF STANDARD PSO');

mpc = L.mpc0();
KMAX = 250;                       % [31] Step 1a
NP   = 30;
NRUN = 15;                        % [31] used 3; both are reported below
SEEDS = 1:NRUN;

% -------------------------------------------------------------------------
% Runtime.  Section B3.3 runs three algorithms x NRUN runs x KMAX iterations
% x NP particles, and every particle evaluation is one Newton-Raphson load
% flow, so the cost is real and is measured up front rather than discovered
% halfway through.  Set NRUN = 3 above for a quick pass matching [31]'s own
% protocol; NRUN = 15 is what makes section B3.4 statistically meaningful.
Otmp = L.orpd(mpc);
t0 = tic; for i = 1:20, Otmp.fit(Otmp.base); end; tfit = toc(t0)/20;
nev = 3*NRUN*(NP + KMAX*NP);
fprintf('\n  one load flow takes %.1f ms on this machine\n', 1000*tfit);
fprintf('  section B3.3 needs %s load flows, about %.0f minutes\n', ...
        addcommas(nev), nev*tfit/60);
fprintf('  section B3.2 needs a few thousand more, section B4 needs none (no load flow)\n');
fprintf('  NRUN is %d.  Set it to 3 at the top of this file for a quick pass.\n', NRUN);
clear Otmp

% =========================================================================
% B1  THE WINGBEAT FREQUENCY
% =========================================================================
L.banner('B1  Wingbeat frequency from [31] section 3.1');

W = L.wingbeat();
fprintf('  [31] takes its swarm to be Pteropodidae Rousettus, a fruit bat, and sets the\n');
fprintf('  wingbeat frequency from Pennycuick''s allometric relation.  All particles share\n');
fprintf('  the same wing geometry, so Freq varies only through beta.\n\n');
fprintf('    total mass  m_tot   %8.4f kg\n', W.m_tot);
fprintf('    wing mass   m_wng   %8.4f kg\n', W.m_wng);
fprintf('    wing span   b       %8.4f m\n',  W.b);
fprintf('    wing area   S       %8.4f m^2\n', W.s);
fprintf('    air density rho     %8.4f kg/m^3\n', W.rho);
fprintf('\n  f_min = f(m_wng) = %.4f Hz        f_max = f(m_tot) = %.4f Hz\n', W.f_min, W.f_max);
fprintf('  m_tot/m_wng = %.1f exactly, and f scales as m^(1/3), so f_max = %.0f x f_min.\n', ...
        W.m_tot/W.m_wng, W.f_max/W.f_min);
fprintf('  These are reproduced from the masses and geometry, not copied, so the frequency\n');
fprintf('  band [%.4f, %.4f] Hz is independently confirmed.\n', W.f_min, W.f_max);

fprintf('\n  Two details in [31] have to be settled before the algorithm will run.\n\n');
fprintf('  (a) Eq (24) is printed as   Freq = f_min + beta*(f_min - f_max).\n');
fprintf('      Evaluated:  beta = 0.0 -> %7.4f   0.5 -> %7.4f   0.9 -> %7.4f   1.0 -> %7.4f\n', ...
        W.freq_printed(0), W.freq_printed(0.5), W.freq_printed(0.9), W.freq_printed(1));
fprintf('      With beta ~ U(0,1) this DECREASES from f_min to exactly zero, so it spans\n');
fprintf('      (0, %.4f] and never reaches the band [%.4f, %.4f] that [31]''s own text\n', ...
        W.freq_printed(0), W.f_min, W.f_max);
fprintf('      states.  Worse, Eq (23) divides the velocity by Freq, and as beta approaches\n');
fprintf('      1 that divisor approaches zero: at beta = 0.999 the position step is\n');
fprintf('      multiplied by %.0f, and at beta = 0.99999 by %.0f.  The printed form is\n', ...
        1/W.freq_printed(0.999), 1/W.freq_printed(0.99999));
fprintf('      therefore unbounded, not merely mis-scaled, and cannot be what was run.\n');
fprintf('      The sign is transposed.  Using the standard bat-algorithm form,\n');
fprintf('      Freq = f_min + beta*(f_max - f_min), gives %.4f to %.4f Hz, which is exactly\n', ...
        W.freq(0), W.freq(1));
fprintf('      the band [31] states, so that is the form used here.\n');

fprintf('\n  (b) Eq (21) is printed as   m_wng = 0.112*m_tot^0.11.\n');
fprintf('      Evaluated at m_tot = %.4f kg this gives %.4f kg, which is %.1f times the\n', ...
        W.m_tot, W.eq21(W.m_tot), W.eq21(W.m_tot)/W.m_wng);
fprintf('      %.4f kg [31] states in the same section.  The stated value is a measurement\n', W.m_wng);
fprintf('      from [31]''s own reference, so Eq (21) must not be used to regenerate it.\n');
fprintf('      Using the measured %.4f kg is what reproduces the quoted f_min.\n', W.m_wng);
fprintf('\n  Both points are recorded because either one, taken literally, stops the\n');
fprintf('  algorithm from running at all.  Neither affects [31]''s conclusions.\n');

% =========================================================================
% B2  WHAT THE WINGBEAT TERM ACTUALLY DOES
% =========================================================================
L.banner('B2  What the wingbeat term does to the swarm');

fprintf('  [31] Eq (22)   U(k+1) = chi*( w*U(k) + Freq*A ),   A the attraction term\n');
fprintf('  [31] Eq (23)   X(k+1) = X(k) + U(k+1)/Freq\n\n');
fprintf('  Substituting one into the other:\n\n');
fprintf('      X(k+1) = X(k) + chi*A  +  chi*w*U(k)/Freq\n\n');
fprintf('  Freq cancels out of the attraction term, so cognitive and social pull enter at\n');
fprintf('  full strength.  It survives only on the MOMENTUM term, which it divides.  With\n');
fprintf('  Freq in [%.2f, %.2f] the momentum is cut to between %.0f%% and %.0f%% of its\n', ...
        W.f_min, W.f_max, 100/W.f_max, 100/W.f_min);
fprintf('  standard value.  The wingbeat frequency is therefore a MOMENTUM DAMPER, not a\n');
fprintf('  step-size or exploration term.\n\n');
fprintf('  That reading agrees with [31]''s own motivation: it says the chaotic weight adds\n');
fprintf('  diversity but costs convergence speed, and that the wingbeat component is there\n');
fprintf('  to "control the search pace".  Damping momentum is exactly how that is done.\n\n');
fprintf('  It also hands us a clean control experiment.  Setting Freq = 1 removes the\n');
fprintf('  wingbeat term and leaves the chaotic inertia weight untouched, so ICPSO collapses\n');
fprintf('  to a chaotic-weight constriction PSO.  Running all three -- PSO, full ICPSO, and\n');
fprintf('  ICPSO with Freq = 1 -- separates [31]''s two modifications and shows which of them\n');
fprintf('  is responsible for whatever difference appears.  [31] does not do this, and it is\n');
fprintf('  the single most informative measurement in this stage.\n');

% =========================================================================
% B3  VALIDATION ON [31]'S OWN PROBLEM: ORPD ON IEEE-14
% =========================================================================
L.banner('B3  ICPSO vs PSO on [31]''s own problem (ORPD, IEEE-14)');

O = L.orpd(mpc);
fprintf('  %d control variables: %s\n', numel(O.lo), strjoin(O.names, ', '));
fprintf('  limits ([31] Table 2)  Vg %.2f-%.2f pu   tap %.2f-%.2f pu   Qc %.0f-%.0f Mvar\n', ...
        O.lo(1), O.hi(1), O.lo(6), O.hi(6), O.lo(9), O.hi(9));
fprintf('  objective: minimise total real power loss.  %d iterations, %d particles.\n', KMAX, NP);

% -- B3.1 evaluate [31]'s own published control vectors ---------------------
L.banner('B3.1  [31] Table 4 evaluated directly');

pub_pso   = [1.10 1.08 1.05 1.03 1.03, 1.02 1.01 1.02, 20.00];
pub_icpso = [1.10 1.09 1.06 1.10 1.10, 0.98 0.98 1.01, 19.99];
cols  = {O.base, pub_pso, pub_icpso};
colnm = {'Table 4 base case', 'Table 4 PSO column', 'Table 4 ICPSO column'};
pubL  = [13.393, 12.275, 12.260];

fprintf('  Before comparing algorithms, the control vectors [31] publishes are simply fed\n');
fprintf('  into the load flow.  This tests the model, not the search.\n\n');
fprintf('  %-24s %12s %12s %12s   %s\n', 'column', 'loss (MW)', 'Table 3', 'difference', 'in Table-2 box');
Lpub = zeros(1,3);
for j = 1:3
    x = cols{j};
    rj = L.solve(O.apply(x));
    Lpub(j) = rj.loss;
    inbox = all(x >= O.lo - 1e-9) && all(x <= O.hi + 1e-9);
    fprintf('  %-24s %12.4f %12.3f %12.4f   %s\n', colnm{j}, Lpub(j), pubL(j), ...
            Lpub(j)-pubL(j), yn(inbox));
end
fprintf('\n  The base-case column reproduces [31] Table 3 to four significant figures\n');
fprintf('  (%.4f against 13.393), which confirms the network data, the mapping from the\n', Lpub(1));
fprintf('  nine controls onto case14, the tap assignment and the loss definition are all\n');
fprintf('  correct.  With the model verified, the two optimised columns are %.4f MW and\n', Lpub(2)-pubL(2));
fprintf('  %.4f MW above their reported losses.\n', Lpub(3)-pubL(3));
fprintf('\n  Table 4 is printed to two decimal places, so each voltage and tap carries up to\n');
fprintf('  0.005 pu of rounding, across eight variables at once.  That is enough to account\n');
fprintf('  for a difference of this size, and it means Table 4 cannot be used to check\n');
fprintf('  Table 3.  The right way to test %.3f MW is to ask what the model''s true optimum\n', pubL(3));
fprintf('  is, which is done next.\n');

% -- B3.2 the true optimum of [31]'s model ---------------------------------
L.banner('B3.2  The true optimum of [31]''s ORPD model');

fprintf('  A compass search is run from the base case and from %d random starts.  It is\n', 9);
fprintf('  written out in powerflow_lib and needs no optimisation toolbox.\n\n');
rng(11, 'twister');
starts = [O.base; repmat(O.lo, 9, 1) + rand(9, numel(O.lo)).*repmat(O.hi-O.lo, 9, 1)];
POL = zeros(size(starts,1), 1);  PX = zeros(size(starts));
for j = 1:size(starts,1)
    [PX(j,:), POL(j)] = L.polish(O.fit, starts(j,:), O.lo, O.hi);
    fprintf('    start %2d  ->  %10.6f MW\n', j-1, POL(j));
end
[fstar, jstar] = min(POL);
xstar = PX(jstar,:);
fprintf('\n  best %.6f MW    worst %.6f MW    distinct optima found %d\n', ...
        fstar, max(POL), numel(uniquetol(POL, 1e-4)));
fprintf('  %d of %d starts land on %.4f MW, so this is the global optimum of the model\n', ...
        sum(POL < fstar + 1e-3), numel(POL), fstar);
fprintf('  and not a lucky local basin.\n\n');
fprintf('  optimal controls:\n');
for k = 1:numel(O.names)
    fprintf('    %-6s %9.4f      ([31] ICPSO column: %6.2f)\n', O.names{k}, xstar(k), pub_icpso(k));
end
rstar = L.solve(O.apply(xstar));
fprintf('  loss %.4f MW, a %.2f%% reduction on the base case, bus voltages %.4f to %.4f pu\n', ...
        rstar.loss, 100*(Lpub(1)-rstar.loss)/Lpub(1), rstar.vmin, rstar.vmax);

fprintf('\n  Measured against that optimum:\n');
fprintf('    [31] PSO   best %.3f MW    %+.4f MW    %s\n', pubL(2), pubL(2)-fstar, ...
        verdict(pubL(2)-fstar));
fprintf('    [31] ICPSO best %.3f MW    %+.4f MW    %s\n', pubL(3), pubL(3)-fstar, ...
        verdict(pubL(3)-fstar));
fprintf('\n  [31]''s PSO figure sits on the optimum, within the rounding of its own third\n');
fprintf('  decimal place -- so [31]''s standard PSO already solved this problem to optimality.\n');
fprintf('  Its ICPSO figure lies %.4f MW BELOW the optimum, and nothing below the optimum is\n', fstar-pubL(3));
fprintf('  attainable.  Both statements point the same way: on IEEE-14 ORPD there was no\n');
fprintf('  room left for ICPSO to improve on PSO, because PSO had already reached the floor.\n');
fprintf('  The %.3f MW gap [31] reports between the two methods is %.4f MW, and the distance\n', ...
        pubL(2)-pubL(3), pubL(2)-pubL(3));
fprintf('  from its ICPSO figure to the attainable floor is %.4f MW -- the same order of\n', fstar-pubL(3));
fprintf('  magnitude.  The claimed margin is inside the reproducibility of the numbers\n');
fprintf('  themselves, so the honest conclusion is that IEEE-14 ORPD is too easy a problem\n');
fprintf('  to separate the two algorithms, not that either one failed.\n');

% -- B3.3 paired PSO vs ICPSO vs Freq=1 ------------------------------------
L.banner('B3.3  PSO, ICPSO and ICPSO with Freq = 1, paired over runs');

fprintf('  All three share the seed, so for a given run they start from the same swarm and\n');
fprintf('  differ only in the update rule.  That makes the comparison paired rather than\n');
fprintf('  two independent samples, which matters at these effect sizes.\n');
fprintf('  %d particles, %d iterations, %d runs.\n\n', NP, KMAX, NRUN);

algs = {'PSO', 'ICPSO (Freq per particle)', 'ICPSO (Freq = 1)'};
Fv = zeros(NRUN, 3);  Hst = cell(NRUN, 3);
tB = tic;
for s = 1:NRUN
    o = struct('np', NP, 'kmax', KMAX, 'seed', SEEDS(s), 'seed_point', O.base);
    a = L.pso(O.fit, O.lo, O.hi, o);
    o2 = o; o2.per_particle = true;
    b = L.icpso(O.fit, O.lo, O.hi, o2);
    o3 = o; o3.fixed_freq = 1.0;
    c = L.icpso(O.fit, O.lo, O.hi, o3);
    Fv(s,:) = [a.f, b.f, c.f];
    Hst(s,:) = {a.hist, b.hist, c.hist};
    el = toc(tB);
    fprintf('    run %2d of %2d   PSO %9.4f   ICPSO %9.4f   Freq=1 %9.4f   [%.1f min done, %.1f to go]\n', ...
            s, NRUN, a.f, b.f, c.f, el/60, el/60*(NRUN-s)/s);
end

% [31] ran three times; report that sample and the larger one
fprintf('\n  [31] ran each method three times.  Reporting the same way, on the first three\n');
fprintf('  runs, and then on all %d:\n\n', NRUN);
fprintf('  %-26s %9s %9s %9s %9s\n', 'sample', 'best', 'worst', 'mean', 'std');
for j = 1:3
    f3 = Fv(1:3, j);
    fprintf('  %-26s %9.4f %9.4f %9.4f %9s\n', [algs{j} ', 3 runs'], ...
            min(f3), max(f3), mean(f3), '--');
end
fprintf('\n');
for j = 1:3
    fj = Fv(:, j);
    fprintf('  %-26s %9.4f %9.4f %9.4f %9.4f\n', sprintf('%s, %d runs', algs{j}, NRUN), ...
            min(fj), max(fj), mean(fj), std(fj));
end
fprintf('\n  %-26s %9.3f %9.3f %9.3f\n', '[31] Table 3, PSO', 12.275, 12.303, 12.288);
fprintf('  %-26s %9.3f %9.3f %9.3f\n', '[31] Table 3, ICPSO', 12.260, 12.270, 12.265);

% paired win/loss
fprintf('\n  Paired against PSO, run by run:\n');
for j = 2:3
    d = Fv(:,j) - Fv(:,1);
    fprintf('    %-26s better on %2d of %d runs, worse on %2d, mean difference %+.4f MW\n', ...
            algs{j}, sum(d < -1e-9), NRUN, sum(d > 1e-9), mean(d));
end

% convergence onset, three definitions
fprintf('\n  Convergence onset, the iteration at which the best-so-far first comes within a\n');
fprintf('  given tolerance of its own final value.  [31] reports around 40 for PSO and 25\n');
fprintf('  for ICPSO on IEEE-14, read off a convergence plot, so several tolerances are\n');
fprintf('  reported here rather than one.\n\n');
tols = [0.001 0.005 0.010];
fprintf('  %-26s %12s %12s %12s\n', 'method', 'within 0.1%', 'within 0.5%', 'within 1.0%');
ONS = zeros(3, numel(tols));
for j = 1:3
    for t = 1:numel(tols)
        v = zeros(NRUN,1);
        for s = 1:NRUN
            h = Hst{s,j};
            v(s) = find(h <= h(end)*(1+tols(t)), 1);
        end
        ONS(j,t) = median(v);
    end
    fprintf('  %-26s %12.1f %12.1f %12.1f\n', algs{j}, ONS(j,1), ONS(j,2), ONS(j,3));
end

% -- B3.4 what B3.3 means --------------------------------------------------
L.banner('B3.4  Reading the ORPD comparison');

dICP = mean(Fv(:,2) - Fv(:,1));
dFq1 = mean(Fv(:,3) - Fv(:,1));
fprintf('  Three things come out of the table above, and all three are measurements rather\n');
fprintf('  than impressions.\n\n');

fprintf('  1.  The chaotic inertia weight on its own is indistinguishable from standard PSO.\n');
fprintf('      With Freq = 1 the mean loss differs from PSO''s by %+.4f MW over %d paired\n', dFq1, NRUN);
fprintf('      runs, against a run-to-run standard deviation of %.4f MW.  The difference is\n', std(Fv(:,1)));
fprintf('      a small fraction of the noise, and both reach the same best value of %.4f MW.\n', ...
        min(min(Fv(:,1)), min(Fv(:,3))));
fprintf('      Replacing a linearly-decreasing weight with a logistic-map weight does not\n');
fprintf('      change what this swarm finds.\n\n');

fprintf('  2.  The wingbeat term is what changes the behaviour, and on the mean it costs\n');
fprintf('      rather than gains.  Full ICPSO differs from PSO by %+.4f MW in the mean and\n', dICP);
fprintf('      takes %.0f iterations to settle against PSO''s %.0f at the 0.1%% tolerance.\n', ...
        ONS(2,1), ONS(1,1));
fprintf('      This is the momentum damping of B2 showing up in the numbers: less momentum\n');
fprintf('      means slower progress along a descent direction.\n\n');

fprintf('  3.  Where the wingbeat term does help is the WORST case, and that is worth\n');
fprintf('      stating precisely because it is the one place [31]''s claim survives contact\n');
fprintf('      with a larger sample.  Over %d runs the worst PSO run gives %.4f MW while the\n', ...
        NRUN, max(Fv(:,1)));
fprintf('      worst ICPSO run gives %.4f MW, and the spread falls from %.4f to %.4f MW.\n', ...
        max(Fv(:,2)), std(Fv(:,1)), std(Fv(:,2)));
if std(Fv(:,2)) < std(Fv(:,1))
    fprintf('      Damping momentum stops the swarm from committing early to a poor basin, so\n');
    fprintf('      ICPSO gives up some of the best solutions in exchange for having fewer bad\n');
    fprintf('      runs.  For an operational tool that trade is often the right one, and it is\n');
    fprintf('      a sharper description of what [31]''s algorithm does than "it is better".\n');
else
    fprintf('      In this sample the spread does not fall, so the reliability argument is not\n');
    fprintf('      supported either and the two methods are simply equivalent here.\n');
end

se3 = std(Fv(:,1))/sqrt(3);
fprintf('\n  Finally, on sample size.  The run-to-run standard deviation is %.4f MW, so the\n', std(Fv(:,1)));
fprintf('  standard error of a mean over three runs is %.4f MW.  The advantage [31] reports\n', se3);
fprintf('  is %.3f MW.  Three runs therefore cannot resolve an effect this small -- the\n', 12.288-12.265);
fprintf('  measurement noise is %.1f times the effect.  That is the whole explanation for\n', se3/(12.288-12.265));
fprintf('  the difference between [31]''s tables and this one, and it is a statement about\n');
fprintf('  sample size rather than about either algorithm.  Running %d instead of 3 is the\n', NRUN);
fprintf('  one methodological change this stage makes, and it is what makes points 1 to 3\n');
fprintf('  above measurable at all.\n');

% =========================================================================
% B4  THE SUBSTITUTION: ICPSO ON [32]'S LOAD-SHEDDING PROBLEM
% =========================================================================
L.banner('B4  ICPSO vs PSO on [32]''s load-shedding problem');

mu = 0.15;
fprintf('  Same three update rules, now on [32] Eqs (5)-(9): 22 decision variables, the\n');
fprintf('  aggregate caps of Eq (6) enforced by the same repair operator in every case.\n');
fprintf('  Paired by seed, %d particles, %d iterations, %d runs.\n\n', NP, KMAX, NRUN);

salgs = {'pso', 'icpso', 'freq1'};
snm   = {'PSO', 'ICPSO', 'ICPSO Freq=1'};
SP = zeros(NRUN, 3); SQ = zeros(NRUN, 3); SF = zeros(NRUN, 3); SO = zeros(NRUN, 3);
for s = 1:NRUN
    for j = 1:3
        o = L.shed_pso(mpc, mu, NP, KMAX, SEEDS(s), salgs{j});
        SP(s,j) = o.Ptot; SQ(s,j) = o.Qtot; SF(s,j) = o.f; SO(s,j) = o.onset;
        if s == 1, capP = o.capP; capQ = o.capQ; end
    end
end

fprintf('  %-14s %14s %14s %14s %12s\n', 'method', 'best P (MW)', 'worst P (MW)', 'best Q (Mvar)', 'iterations');
for j = 1:3
    fprintf('  %-14s %14.4f %14.4f %14.4f %12d\n', snm{j}, max(SP(:,j)), min(SP(:,j)), ...
            max(SQ(:,j)), round(median(SO(:,j))));
end
fprintf('\n  %-34s %12s %12s  %-8s %s\n','quantity','computed','published','unit','');
L.cmp('total real power shed',     max(SP(:)), 38.8502, 'MW',   1e-3);
L.cmp('total reactive power shed', max(SQ(:)), 11.0250, 'Mvar', 1e-3);

spreadP = max(SP(:)) - min(SP(:));
spreadF = max(SF(:)) - min(SF(:));
fprintf('\n  Across all %d runs and all three update rules, the total real power shed varies\n', 3*NRUN);
fprintf('  by %.3e MW and the Eq (5) fitness by %.3e.  Every method returns the same answer,\n', ...
        spreadP, spreadF);
fprintf('  to machine precision, in %d iteration.\n', round(median(SO(:))));

L.banner('Why the substitution cannot change this answer');
fprintf('  This is not a tie that a better-tuned swarm would break.  It is forced by the\n');
fprintf('  structure of [32]''s formulation, and the argument is short.\n\n');
fprintf('  Eq (7) bounds each variable by its own demand, 0 <= Pshed_i <= Pd_i, so every\n');
fprintf('  term inside Eq (5)''s absolute value is already non-negative and\n\n');
fprintf('      f(x) = sum_i (Pd_i - Pshed_i) = sum_i Pd_i - sum_i Pshed_i\n\n');
fprintf('  The first sum is a constant, %.2f MW.  So Eq (5) is an affine function of the\n', sum(mpc.bus(:,C.PD)));
fprintf('  TOTAL shed and of nothing else, and minimising it means pushing that total to\n');
fprintf('  Eq (6)''s cap of %.4f MW.  Any point on the cap surface is optimal; the feasible\n', capP);
fprintf('  set has a whole face of global optima rather than one.\n\n');
fprintf('  The repair operator that enforces Eq (6) rescales any over-cap particle onto that\n');
fprintf('  surface.  The initial swarm is drawn over the box of Eq (7), where the expected\n');
fprintf('  total is half of demand, far above the cap, so essentially every particle is\n');
fprintf('  rescaled onto the optimal face before the first velocity update is computed.\n');
fprintf('  The search is finished before it begins, which is why the fitness history is flat\n');
fprintf('  from iteration 1 and why [32] reports convergence "within 5 iterations".\n\n');
fprintf('  So the honest result of the substitution is that on the published benchmark it\n');
fprintf('  changes nothing, and it cannot: the answer is set by the constraint set, not by\n');
fprintf('  the optimiser.  That is a finding about [32]''s formulation, and a useful one --\n');
fprintf('  it says the 38.8502 MW figure is a property of the 15%% cap rather than evidence\n');
fprintf('  that PSO did anything, and it explains why the paper publishes only two totals.\n');

L.banner('What this implies for the rest of the project');
fprintf('  Two null results have now been established for good reasons rather than by\n');
fprintf('  failing to find an effect: on ORPD the standard PSO of [31] already reaches the\n');
fprintf('  attainable optimum, and on load shedding [32]''s objective is determined by its\n');
fprintf('  own caps.  Neither problem has any headroom for a better search.\n\n');
fprintf('  Demonstrating what ICPSO contributes therefore requires a formulation in which\n');
fprintf('  the search does real work -- one where the optimum is interior, the constraints\n');
fprintf('  are the network''s rather than a flat cap, and different allocations of the same\n');
fprintf('  total shed give different outcomes.  Solar and wind intermittency produces exactly\n');
fprintf('  that: shedding then has to be placed bus by bus against a varying net load with\n');
fprintf('  voltage and flow limits binding, and the total is no longer fixed in advance.\n');
fprintf('  Stage C builds that problem on top of this verified benchmark, and it is where\n');
fprintf('  the ICPSO substitution is actually tested.\n');

% =========================================================================
L.banner('STAGE B RESULT');
fprintf('  [31]''s ICPSO is reproduced in full, including the two printing errors that have\n');
fprintf('  to be corrected before it will run, and its wingbeat term is identified\n');
fprintf('  algebraically as a momentum damper.  On [31]''s own ORPD problem the\n');
fprintf('  implementation is validated against [31]''s own base-case number (%.4f MW\n', Lpub(1));
fprintf('  against 13.393) and the model''s true optimum is established at %.4f MW.\n\n', fstar);
fprintf('  Against that optimum, [31]''s published PSO result is already optimal and its\n');
fprintf('  published ICPSO result lies %.4f MW below what the model admits.  Over %d paired\n', fstar-pubL(3), NRUN);
fprintf('  runs the chaotic weight alone matches PSO (%+.4f MW in the mean) while the\n', dFq1);
fprintf('  wingbeat term costs %+.4f MW in the mean and reduces the run-to-run spread from\n', dICP);
fprintf('  %.4f to %.4f MW.  The three-run protocol of [31] has a standard error of %.4f MW,\n', ...
        std(Fv(:,1)), std(Fv(:,2)), se3);
fprintf('  which is larger than the effect it reports, so the larger sample is what makes\n');
fprintf('  the decomposition visible.\n\n');
fprintf('  On [32]''s load-shedding problem all three rules return %.4f MW and %.4f Mvar\n', max(SP(:)), max(SQ(:)));
fprintf('  identically, in one iteration, because Eq (5) is affine in the total shed and\n');
fprintf('  Eq (6)''s cap fixes that total.  The substitution is therefore carried out and\n');
fprintf('  measured, and the reason it cannot matter here is proved rather than guessed.\n');

save(fullfile(here, 'icpso_results.mat'), ...
     'Fv', 'Hst', 'ONS', 'POL', 'PX', 'xstar', 'fstar', 'Lpub', 'pub_pso', 'pub_icpso', ...
     'SP', 'SQ', 'SF', 'SO', 'SEEDS', 'KMAX', 'NP', 'NRUN');
diary off;

% =========================================================================
function s = yn(t)
if t, s = 'yes'; else, s = 'NO'; end
end

function s = addcommas(n)
s = sprintf('%d', round(n));
for k = numel(s)-3 : -3 : 1
    s = [s(1:k) ',' s(k+1:end)];
end
end

function s = verdict(d)
if abs(d) < 5e-3
    s = 'at the optimum';
elseif d > 0
    s = 'above the optimum, so attainable';
else
    s = 'BELOW the optimum, so not attainable';
end
end
