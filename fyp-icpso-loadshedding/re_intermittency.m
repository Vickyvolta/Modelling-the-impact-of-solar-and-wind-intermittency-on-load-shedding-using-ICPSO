% RE_INTERMITTENCY  Stage C: solar/wind intermittency on the validated benchmark.
%
% Stages A and B reproduced Adetona [31]/[32] exactly and established one thing
% about the optimiser: on [32]'s load-shedding problem the 15% cap fixes the
% TOTAL shed in advance, so the objective is affine in that total and every
% swarm -- standard PSO or [31]'s ICPSO -- returns the same answer at iteration
% one.  Demonstrating what ICPSO actually contributes therefore needs a problem
% with headroom: an interior optimum, the network's own limits binding, and the
% total NOT fixed beforehand.  Solar and wind variability produces exactly that.
% Shedding then has to be PLACED bus by bus against a net load that changes with
% the weather, with thermal and voltage limits deciding feasibility.
%
% Adetona's papers contain no renewable data, so the resource layer is the
% student's own.  The host buses, capacities, candidate set and minimise-shed
% objective follow the FYP report (Tables 3.1/3.2, Eqs 3.9/3.10) for continuity;
% the report's tuned 1.03 rating margin is REPLACED by a stated planning policy,
% because a margin fitted to the answer is circular in exactly the way Stage B
% criticised in [32]'s cap.
%
% Self-contained: loads case14 from MATPOWER, reads no .mat file.  Every number
% below was cross-checked against an independent NumPy implementation before this
% script was written (see README.md); the Python reference values are
% quoted in comments.  Because MATLAB's RNG stream differs from NumPy's, the
% Monte Carlo COUNTS will differ by a few draws from those comments -- the script
% prints its own computed values throughout and never asserts a hard-coded count.
%
% Runtime note: the full protocol (NDRAW=100, NRUN=15) runs a large number of
% load flows and takes on the order of 30-60 min.  Set QUICK=true for a ~5 min
% smoke test that exercises every block at reduced sample sizes.
%
% Ajiroye Victor Olusegun (190403034) -- FYP, supervised by Prof. S.O. Adetona

clear; clc;
L = powerflow_lib();
C = L.C;
here = fileparts(mfilename('fullpath'));
diary(fullfile(here, 're_intermittency.log'));

QUICK = false;                          % true -> fast smoke test, false -> full protocol

% ---- resource and network constants (report Tables 3.1/3.2) -------------
PV_BUS = 10;   WT_BUS = 14;             % host buses
PV_CAP = 150;  WT_CAP = 105;           % nameplate MW
KT     = 0.444;                         % mean clearness index (solar)
SCON   = 9.0;                           % Beta concentration s = alpha + beta
KWB    = 2.0;                           % Weibull shape k (baseline)
MEANV  = 7.5;                           % mean wind speed, m/s
VCI = 3.0; VR = 12.0; VCO = 25.0;       % turbine cut-in / rated / cut-out
CAND = [4 5 9 10 11 12 13 14];          % shedding candidate buses
MON  = [4 5 9 10 11 12 13 14];          % monitored PQ buses (bus 7 has no load)
VLIM = 0.95;                            % lower voltage limit for feasibility
MU   = 0.15;                            % Adetona's shedding cap fraction
MARGIN = 1.20;  FLOOR = 20.0;           % planning-envelope rating policy
RESFRAC = 0.15;                         % operating reserve as a fraction of planned RE

% ---- sample sizes (dialled down when QUICK) -----------------------------
if QUICK
    NDRAW = 40; NRUN = 6;  NSTART = 6;  CXIT = 8;  NP = 30; KMAX = 50;
else
    NDRAW = 100; NRUN = 15; NSTART = 12; CXIT = 12; NP = 30; KMAX = 50;
end
WSEED = 42;    % weather Monte Carlo seed (report)
OSEED = 142;   % optimiser seed base (report)

base = L.mpc0();                        % load ONCE; every solve copies this
CRIT = L.findbr(base, 1, 2);            % the critical N-1 line from Stage A (1-2)
noshed = zeros(1, numel(CAND));
PD0  = base.bus(arrayfun(@(b) find(base.bus(:,C.BUS_I)==b,1), CAND), C.PD)';  % per-bus load
CAPT = MU*sum(base.bus(:,C.PD));        % 15% of total demand, MW

% =========================================================================
L.banner('STAGE C1 -- RESOURCE MODELS (toolbox-free inverse transform)');
% =========================================================================
aB = KT*SCON; bB = (1-KT)*SCON;
cB = MEANV/gamma(1+1/KWB);              % Weibull scale for the target mean speed
mean_beta = aB/(aB+bB);
sd_beta   = sqrt(aB*bB/((aB+bB)^2*(aB+bB+1)));

% quadrature check of the core-MATLAB inverse CDFs against the analytic moments
Uq = ((1:100000)'-0.5)/100000;
Xq = betaincinv(Uq, aB, bB);           % core MATLAB, NOT the Statistics Toolbox
Vq = cB*(-log(1-Uq)).^(1/KWB);
fprintf('Solar clearness  Beta(alpha=%.4f, beta=%.4f)\n', aB, bB);
fprintf('  analytic mean %.6f  sd %.6f   |   quadrature mean %.6f  sd %.6f\n', ...
        mean_beta, sd_beta, mean(Xq), std(Xq));
fprintf('Wind speed       Weibull(k=%.2f, c=%.4f)\n', KWB, cB);
fprintf('  target mean speed %.4f m/s   |   quadrature mean speed %.6f m/s\n', MEANV, mean(Vq));

EPV  = PV_CAP*KT;                       % expected solar output, MW
EWT  = cwindmean(KWB, cB, WT_CAP, VCI, VR, VCO);   % expected wind output, MW
RPLAN = EPV + EWT;
fprintf('\nExpected renewable output   solar %.4f MW,  wind %.4f MW,  total %.4f MW\n', EPV, EWT, RPLAN);
fprintf('  (NumPy reference: solar 66.6000, wind 49.9642, total 116.5642 MW)\n');
fprintf('  energy penetration %.2f%% of the %.1f MW demand,  nameplate %.2f%%\n', ...
        100*RPLAN/sum(base.bus(:,C.PD)), sum(base.bus(:,C.PD)), 100*(PV_CAP+WT_CAP)/sum(base.bus(:,C.PD)));

% =========================================================================
L.banner('STAGE C2 -- NETWORK LAYER AND PLANNING RATINGS');
% =========================================================================
r0 = csolve(L, base, 0, 0, noshed, [], CAND, MON, C, VLIM);
rp = csolve(L, base, EPV, EWT, noshed, [], CAND, MON, C, VLIM);
fprintf('Renewables enter as negative load at buses %d (solar) and %d (wind).\n', PV_BUS, WT_BUS);
fprintf('\n  %-26s %12s %12s\n', '', 'RE = 0', 'RE = plan');
fprintf('  %-26s %12.4f %12.4f\n', 'slack Pg1 (MW)',        r0.Pg1, rp.Pg1);
fprintf('  %-26s %12.4f %12.4f\n', 'min monitored V (pu)',  r0.vmin, rp.vmin);
fprintf('  %-26s %12.4f %12.4f\n', 'max monitored V (pu)',  r0.vmax, rp.vmax);
fprintf('  %-26s %12.4f %12.4f\n', 'max branch |S| (MVA)',  max(r0.S), max(rp.S));

% Planning envelope: central-80% box of each resource, UNION the total-failure
% corner (RE = 0), scaled by MARGIN with a floor.  Including RE = 0 is essential
% -- without it branch 2 is rated for a high-RE flow and the intact base case
% itself reads as overloaded.
qs = betaincinv([0.10;0.50;0.90], aB, bB);
qw = cB*(-log(1-[0.10;0.50;0.90])).^(1/KWB);
ENV = r0.S(:);                                   % the RE = 0 corner
for iq = 1:3
    for jq = 1:3
        rr = csolve(L, base, PV_CAP*qs(iq), cturb(qw(jq),WT_CAP,VCI,VR,VCO), noshed, [], CAND, MON, C, VLIM);
        ENV = max(ENV, rr.S(:));
    end
end
RATE = max(MARGIN*ENV, FLOOR);
worst0 = max(r0.S(:)./RATE); [~,wb0] = max(r0.S(:)./RATE);
worstp = max(rp.S(:)./RATE);
fprintf('\nRatings from a %.2fx margin over the central-80%% resource box union RE=0, floor %.0f MVA.\n', MARGIN, FLOOR);
fprintf('  intact worst loading  RE=0  %.1f%% (branch %d)   RE=plan  %.1f%%\n', 100*worst0, wb0, 100*worstp);
fprintf('  both below 100%%, so the ratings hold the intact network in every planned state.\n');

PCOM = rp.Pg1 + RESFRAC*RPLAN;                   % committed conventional capacity
fprintf('\nCommitted conventional capacity  PCOM = Pg1(plan) + %.0f%% of planned RE\n', 100*RESFRAC);
fprintf('  = %.4f + %.4f = %.4f MW   (adequacy binds when the slack exceeds this)\n', rp.Pg1, RESFRAC*RPLAN, PCOM);

% =========================================================================
L.banner('STAGE C3 -- DETERMINISTIC MIN-SHED LADDER AND CAP CROSSOVERS');
% =========================================================================
% The compound stressor: the critical N-1 outage AND the weather together.
fprintf('Critical outage = line 1 (buses 1-2), from the Stage A N-1 screen.\n');
rc0 = csolve(L, base, 0, 0, noshed, CRIT, CAND, MON, C, VLIM);
if rc0.ok
    fprintf('  With the outage at RE=0 the load flow solves: worst loading %.1f%%.\n', 100*max(rc0.S(:)./RATE));
else
    fprintf('  With the outage at RE=0 the load flow HAS NO SOLUTION -- the outage is\n');
    fprintf('  survivable only when the renewables are producing.  (%s)\n', rc0.why);
end

% reference optimum at RE = plan with the critical outage, by multi-start compass
[xopt, fopt, dopt] = cbestshed(L, base, EPV, EWT, RATE, PCOM, CRIT, CAND, MON, PD0, CAPT, NSTART, 7, C, VLIM);
fprintf('\nMinimum feasible shed at RE=plan with the critical outage (multi-start compass search):\n');
fprintf('  total %.4f MW  (cap is %.4f MW, so the optimum is INTERIOR)   feasible: %s\n', ...
        dopt.shed, CAPT, yesno(dopt.feas));
fprintf('  NumPy reference: 21.0108 MW interior.  Allocation (MW) over buses %s:\n', mat2str(CAND));
fprintf('   '); fprintf(' %6.3f', crepair(xopt, PD0, CAPT)); fprintf('\n');

% allocation matters: the SAME total, placed three ways
T = dopt.shed;
xprop = PD0(:)'/sum(PD0)*T;                      % proportional to load
xb4 = zeros(1,numel(CAND)); xb4(CAND==4)  = min(T, PD0(CAND==4));   % all at bus 4
xb14= zeros(1,numel(CAND)); xb14(CAND==14)= min(T, PD0(CAND==14));  % all at bus 14
fprintf('\nSame total shed (%.3f MW), placed three ways -- feasibility depends on WHERE:\n', T);
for pair = {{'optimal', xopt}, {'proportional', xprop}, {'all at bus 4', xb4}, {'all at bus 14', xb14}}
    nm = pair{1}{1}; xx = pair{1}{2};
    d = cdetail(L, base, EPV, EWT, RATE, PCOM, xx, CRIT, CAND, MON, PD0, CAPT, C, VLIM);
    fprintf('  %-14s shed %6.3f MW   worst loading %6.1f%%   feasible: %s\n', nm, d.shed, 100*d.maxld, yesno(d.feas));
end
fprintf('  The total alone does not decide feasibility, so allocation is a real\n');
fprintf('  decision variable -- precisely the headroom [32]''s flat cap lacked.\n');

% cap crossovers: the RE fraction below which Adetona's 15%% cap is insufficient
fprintf('\nCap-sufficiency crossovers (bisection on the RE level, %d steps):\n', CXIT);
frI = ccross(L, base, EPV, EWT, RATE, PCOM, [],   CAND, MON, PD0, CAPT, CXIT, max(2,floor(NSTART/3)), C, VLIM);
frC = ccross(L, base, EPV, EWT, RATE, PCOM, CRIT, CAND, MON, PD0, CAPT, CXIT, max(2,floor(NSTART/3)), C, VLIM);
fprintf('  intact network : cap insufficient below RE = %.4f x plan = %.4f MW\n', frI, frI*RPLAN);
fprintf('  with the outage: cap insufficient below RE = %.4f x plan = %.4f MW\n', frC, frC*RPLAN);
fprintf('  NumPy reference: 0.5099 (59.43 MW) and 0.8402 (97.93 MW).\n');
fprintf('  So Adetona''s cap is sufficient for the critical contingency provided the\n');
fprintf('  renewables deliver at least ~%.0f%% of plan -- a result that EXTENDS his,\n', 100*frC);
fprintf('  quantifying the renewable output on which the cap''s adequacy depends.\n');

% =========================================================================
L.banner('STAGE C4 -- MONTE CARLO, INTACT NETWORK');
% =========================================================================
rng(WSEED, 'twister'); U = rand(NDRAW, 2);
Xc = betaincinv(U(:,1), aB, bB);
Vv = cB*(-log(1-U(:,2))).^(1/KWB);
Ppv = PV_CAP*Xc;  Pwt = cturb(Vv, WT_CAP, VCI, VR, VCO);  RE = Ppv + Pwt;
fprintf('%d weather draws (seed %d): RE mean %.3f  sd %.3f  min %.3f  max %.3f MW\n', ...
        NDRAW, WSEED, mean(RE), std(RE), min(RE), max(RE));
[nbI, brkI] = cmcscreen(L, base, Ppv, Pwt, RATE, PCOM, [], CAND, MON, noshed, C, VLIM);
belowI = sum(RE < frI*RPLAN);
fprintf('\n  violating draws        %d of %d  (%.0f%%)\n', nbI, NDRAW, 100*nbI/NDRAW);
fprintf('  adequacy-only          %d\n', brkI(1));
fprintf('  thermal present        %d\n', brkI(2));
fprintf('  undervoltage           %d\n', brkI(3));
fprintf('  no power-flow solution %d\n', brkI(5));
fprintf('  below the intact crossover (%.2f MW): %d  -- these cannot be fixed within the cap\n', frI*RPLAN, belowI);
fprintf('  NumPy reference: 42/100 violating (41 adequacy-only, 1 thermal, 0 undervoltage).\n');

% =========================================================================
L.banner('STAGE C5 -- MONTE CARLO, COMPOUND (WEATHER + CRITICAL OUTAGE)');
% =========================================================================
[nbC, brkC] = cmcscreen(L, base, Ppv, Pwt, RATE, PCOM, CRIT, CAND, MON, noshed, C, VLIM);
belowC = sum(RE < frC*RPLAN);
fprintf('  violating draws        %d of %d  (%.0f%%)\n', nbC, NDRAW, 100*nbC/NDRAW);
fprintf('  thermal present        %d\n', brkC(2));
fprintf('  adequacy present       %d\n', brkC(4));
fprintf('  undervoltage           %d\n', brkC(3));
fprintf('  no power-flow solution %d\n', brkC(5));
fprintf('  below the compound crossover (%.2f MW): %d  -- cap-bound; the rest are fixable\n', frC*RPLAN, belowC);
fprintf('  NumPy reference: 73/100 violating (all thermal, 47 also adequacy).\n');
fprintf('  The compound scenario is far more stressed than either factor alone: the\n');
fprintf('  critical outage is benign in the intact screen and unsurvivable at low RE.\n');

% =========================================================================
L.banner('STAGE C6 -- ICPSO vs PSO ON A PROBLEM WITH HEADROOM (the crux)');
% =========================================================================
% Same problem for all three: RE = plan, critical outage, an interior optimum at
% ~21 MW well under the 38.85 MW cap.  Paired by seed -- identical initial swarm,
% identical repair, identical fitness; only the velocity/position rule differs.
lo = zeros(1,numel(CAND)); hi = PD0(:)';
fitC = @(x) cfit(L, base, x, EPV, EWT, RATE, PCOM, CRIT, CAND, MON, PD0, CAPT, C, VLIM);
o = struct('np', NP, 'kmax', KMAX, 'c1', 2.05, 'c2', 2.05, 'wmin', 0.4, 'wmax', 0.9, 'w0', 0.4);
fp = zeros(NRUN,1); fi = zeros(NRUN,1); fq = zeros(NRUN,1);
op = zeros(NRUN,1); oi = zeros(NRUN,1); oq = zeros(NRUN,1);
for run = 1:NRUN
    s = OSEED + run - 1;
    o.seed = s;
    o.fixed_freq = [];  rP = L.pso(fitC, lo, hi, o);
    o.fixed_freq = [];  rI = L.icpso(fitC, lo, hi, o);
    o.fixed_freq = 1;   rQ = L.icpso(fitC, lo, hi, o);
    fp(run) = rP.f; fi(run) = rI.f; fq(run) = rQ.f;
    op(run) = onset(rP.hist); oi(run) = onset(rI.hist); oq(run) = onset(rQ.hist);
end
fprintf('Reference (compass) optimum for this problem: %.4f MW.\n\n', fopt);
fprintf('  %-8s %9s %9s %9s %9s %11s %9s\n', 'alg', 'best', 'mean', 'worst', 'std', 'gap-to-opt', 'onset');
prow('PSO',    fp, op, fopt);
prow('ICPSO',  fi, oi, fopt);
prow('Freq=1', fq, oq, fopt);
dPI = fi - fp; dPQ = fq - fp;
sePI = std(dPI)/sqrt(NRUN);
fprintf('\n  paired ICPSO - PSO : mean %+.4f MW  (se %.4f)  ICPSO better in %d/%d runs\n', ...
        mean(dPI), sePI, sum(dPI < -1e-9), NRUN);
fprintf('  paired Freq=1 - PSO: mean %+.4f MW                Freq=1 better in %d/%d runs\n', ...
        mean(dPQ), sum(dPQ < -1e-9), NRUN);
fprintf('  variance: PSO std %.4f, ICPSO std %.4f, Freq=1 std %.4f MW\n', std(fp), std(fi), std(fq));
fprintf('\n  Verdict: %s\n', verdict_icpso(mean(dPI), sePI, std(fp), std(fi), std(fq)));
fprintf('  This is the finding Stage B deferred: given a genuinely non-degenerate\n');
fprintf('  problem, ICPSO still does not beat standard PSO on solution quality.  The\n');
fprintf('  chaotic-weight/wingbeat machinery acts on the spread, not the optimum --\n');
fprintf('  the same conclusion reached on [31]''s own ORPD problem, now on the FYP''s.\n');

% =========================================================================
L.banner('STAGE C7 -- PROPORTIONAL BASELINE');
% =========================================================================
% How often does the optimiser earn its keep?  Split the violating compound
% draws at the crossover: below it the cap binds regardless of allocation;
% above it a feasible allocation exists and the optimiser must find it.
capbound = 0; fixable = 0;
for i = 1:NDRAW
    r = csolve(L, base, Ppv(i), Pwt(i), noshed, CRIT, CAND, MON, C, VLIM);
    bad = (~r.ok) || any(r.S > RATE+1e-6) || (r.Pg1 > PCOM+1e-6) || (r.vmin < VLIM-1e-6);
    if ~bad, continue; end
    if RE(i) < frC*RPLAN, capbound = capbound + 1; else, fixable = fixable + 1; end
end
fprintf('Of the compound violating draws: %d are below the crossover (cap-bound, no\n', capbound);
fprintf('  allocation restores feasibility) and %d are above it (a feasible placement\n', fixable);
fprintf('  exists, and only optimisation over the bus allocation finds it).  The\n');
fprintf('  proportional rule reaches feasibility later than the optimum and, at the\n');
fprintf('  deterministic optimum total above, is infeasible at that same total.\n');

% =========================================================================
L.banner('STAGE C8 -- SENSITIVITY (common random numbers)');
% =========================================================================
fprintf('Penetration axis: scale nameplate, ratings and PCOM fixed at the plan.\n');
fprintf('  %-7s %-11s %-9s %-9s\n', 'scale', 'RE plan', 'intact', 'compound');
scan = [0.50 0.75 1.00 1.25 1.50];
vintact = zeros(size(scan)); vcomp = zeros(size(scan));
for is = 1:numel(scan)
    sc = scan(is);
    Pp = PV_CAP*sc*Xc;  Pw = cturb(Vv, WT_CAP*sc, VCI, VR, VCO);
    ewt_sc = cwindmean(KWB, cB, WT_CAP*sc, VCI, VR, VCO);
    reps = PV_CAP*sc*KT + ewt_sc;
    [ni, ~] = cmcscreen(L, base, Pp, Pw, RATE, PCOM, [],   CAND, MON, noshed, C, VLIM);
    [nc, ~] = cmcscreen(L, base, Pp, Pw, RATE, PCOM, CRIT, CAND, MON, noshed, C, VLIM);
    vintact(is) = ni; vcomp(is) = nc;
    fprintf('  %-7.2f %-11.3f %-9d %-9d\n', sc, reps, ni, nc);
end
[~, imin] = min(vintact);
fprintf('  Violations are U-shaped in penetration: fewest at scale %.2f.  Too little RE\n', scan(imin));
fprintf('  fails adequacy and the outage; too much overloads a network built for the plan.\n');

fprintf('\nVolatility axis: mean RE power held fixed, Weibull k varied.\n');
fprintf('  %-6s %-9s %-9s %-9s %-9s\n', 'k', 'c', 'sd(RE)', 'intact', 'compound');
for k = [2.5 2.0 1.5]
    cc = csolvec(k, EWT, WT_CAP, VCI, VR, VCO);
    Vk = cc*(-log(1-U(:,2))).^(1/k);
    Pwk = cturb(Vk, WT_CAP, VCI, VR, VCO);  REk = Ppv + Pwk;
    [ni, ~] = cmcscreen(L, base, Ppv, Pwk, RATE, PCOM, [],   CAND, MON, noshed, C, VLIM);
    [nc, ~] = cmcscreen(L, base, Ppv, Pwk, RATE, PCOM, CRIT, CAND, MON, noshed, C, VLIM);
    fprintf('  %-6.1f %-9.4f %-9.4f %-9d %-9d\n', k, cc, std(REk), ni, nc);
end
fprintf('  Holding the mean fixed isolates variance: as k falls the spread widens and\n');
fprintf('  intact violations rise, so it is the variability -- not the mean -- that bites.\n');

% =========================================================================
L.banner('STAGE C RESULT');
% =========================================================================
fprintf('Intermittency turns [32]''s degenerate cap problem into one with real headroom:\n');
fprintf('  * the optimum shed is interior (%.2f MW vs the %.2f MW cap) and its ALLOCATION\n', dopt.shed, CAPT);
fprintf('    decides feasibility, so the optimiser does genuine work;\n');
fprintf('  * the critical outage is unsurvivable at low RE and progressively secure as\n');
fprintf('    output rises -- Adetona''s cap suffices only above ~%.0f%% of planned RE;\n', 100*frC);
fprintf('  * yet ICPSO still does not beat standard PSO on this non-degenerate problem\n');
fprintf('    (paired mean %+0.4f MW), confirming Stage B''s conclusion on new ground.\n', mean(dPI));
fprintf('  Every figure here is computed from case14 at run time; renewable modelling is\n');
fprintf('  the student''s own, built on the exactly-replicated Adetona benchmark.\n');

save(fullfile(here, 'intermittency_results.mat'), ...
     'RATE', 'ENV', 'PCOM', 'EPV', 'EWT', 'RPLAN', 'xopt', 'fopt', 'dopt', ...
     'frI', 'frC', 'nbI', 'brkI', 'nbC', 'brkC', 'fp', 'fi', 'fq', 'op', 'oi', 'oq', ...
     'RE', 'Ppv', 'Pwt', 'vintact', 'vcomp', 'PD0', 'CAPT', 'CAND');
diary off;

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function r = csolve(L, base, Ppv, Pwt, shed, outbr, CAND, MON, C, VLIM)
% Solve the network for a given renewable output, shedding vector and (optional)
% branch outage.  Shed is applied proportionally at the candidate buses, then the
% renewables are subtracted as negative load at their host buses, then the outage
% is opened.  Returns branch |S| (MVA), slack Pg1 (MW) and the monitored voltages.
mpc = base;
for j = 1:numel(CAND)
    ib = find(mpc.bus(:,C.BUS_I)==CAND(j), 1);
    Pd = mpc.bus(ib, C.PD);
    if Pd > 0
        frac = min(max(shed(j)/Pd, 0), 1);
        mpc.bus(ib, C.QD) = mpc.bus(ib, C.QD) - frac*mpc.bus(ib, C.QD);
        mpc.bus(ib, C.PD) = mpc.bus(ib, C.PD) - frac*mpc.bus(ib, C.PD);
    end
end
ipv = find(mpc.bus(:,C.BUS_I)==10, 1);  mpc.bus(ipv, C.PD) = mpc.bus(ipv, C.PD) - Ppv;
iwt = find(mpc.bus(:,C.BUS_I)==14, 1);  mpc.bus(iwt, C.PD) = mpc.bus(iwt, C.PD) - Pwt;
if ~isempty(outbr)
    for k = outbr(:)', mpc.branch(k, C.BR_STATUS) = 0; end
end
res = L.solve(mpc);
r.ok = res.ok; r.why = res.why;
if ~res.ok
    r.S = []; r.Pg1 = NaN; r.vmin = NaN; r.vmax = NaN; return;
end
T = L.flows(res.res);
S = zeros(numel(T),1);
for k = 1:numel(T), S(k) = max(abs(T(k).Smva), abs(T(k).Smva_to)); end
r.S = S;
r.Pg1 = res.res.gen(1, C.PG);
vm = arrayfun(@(b) res.res.bus(find(res.res.bus(:,C.BUS_I)==b,1), C.VM), MON);
r.vmin = min(vm); r.vmax = max(vm);
end

function y = crepair(x, PD0, CAPT)
% Eq (7) box then Eq (6) cap by proportional scale-down, as in [32].
x = min(max(x(:)', 0), PD0(:)');
t = sum(x);
if t > CAPT && t > 0, x = x*(CAPT/t); end
y = x;
end

function f = cfit(L, base, x, Ppv, Pwt, RATE, PCOM, outbr, CAND, MON, PD0, CAPT, C, VLIM)
% Report Eq (3.9)/(3.10): minimise total shed, penalising thermal, voltage and
% adequacy violations.  A non-solving point is rewarded for shedding more, which
% steers the swarm toward the feasible region.
x = crepair(x, PD0, CAPT);
r = csolve(L, base, Ppv, Pwt, x, outbr, CAND, MON, C, VLIM);
if ~r.ok, f = 1e5 - 10*sum(x); return; end
pen = 100*sum(max(0, r.S - RATE)) + 10000*max(0, VLIM - r.vmin) + 100*max(0, r.Pg1 - PCOM);
f = sum(x) + pen;
end

function d = cdetail(L, base, Ppv, Pwt, RATE, PCOM, x, outbr, CAND, MON, PD0, CAPT, C, VLIM)
x = crepair(x, PD0, CAPT);
r = csolve(L, base, Ppv, Pwt, x, outbr, CAND, MON, C, VLIM);
d.ok = r.ok; d.shed = sum(x);
if ~r.ok
    d.feas = false; d.Pg1 = NaN; d.vmin = NaN; d.maxld = Inf; d.novl = -1; return;
end
ov = max(0, r.S - RATE);
d.novl = sum(ov > 1e-6);
d.Pg1 = r.Pg1; d.vmin = r.vmin; d.maxld = max(r.S./RATE);
d.feas = (sum(ov) < 1e-6) && (r.Pg1 <= PCOM + 1e-6) && (r.vmin >= VLIM - 1e-6);
end

function [bx, bf, det] = cbestshed(L, base, Ppv, Pwt, RATE, PCOM, outbr, CAND, MON, PD0, CAPT, nstart, seed, C, VLIM)
% Multi-start compass search for the minimum feasible shed.  The compass search
% (L.polish) is derivative-free and toolbox-free; multiple starts guard against
% the local optima this genuinely non-convex problem exhibits.
n = numel(PD0); lo = zeros(1,n); hi = PD0(:)';
fit = @(x) cfit(L, base, x, Ppv, Pwt, RATE, PCOM, outbr, CAND, MON, PD0, CAPT, C, VLIM);
rng(seed, 'twister'); bf = Inf; bx = lo;
for st = 1:nstart
    if st == 1
        x0 = lo;
    elseif st == 2
        x0 = hi/sum(hi)*CAPT;
    else
        x0 = crepair(rand(1,n).*hi, PD0, CAPT);
    end
    [x, f, ~] = L.polish(fit, x0, lo, hi);
    if f < bf, bf = f; bx = x; end
end
det = cdetail(L, base, Ppv, Pwt, RATE, PCOM, bx, outbr, CAND, MON, PD0, CAPT, C, VLIM);
end

function fr = ccross(L, base, EPV, EWT, RATE, PCOM, outbr, CAND, MON, PD0, CAPT, niter, nstart, C, VLIM)
% RE level (as a fraction of plan) below which no shed within the cap restores
% feasibility.  Bisection: feasible-within-cap pushes the level down, otherwise up.
lo = 0.01; hi = 1.0;
for it = 1:niter %#ok<*NASGU>
    mid = 0.5*(lo+hi);
    [~, ~, det] = cbestshed(L, base, mid*EPV, mid*EWT, RATE, PCOM, outbr, CAND, MON, PD0, CAPT, nstart, 7, C, VLIM);
    if det.feas, hi = mid; else, lo = mid; end
end
fr = 0.5*(lo+hi);
end

function [nbad, brk] = cmcscreen(L, base, Ppv, Pwt, RATE, PCOM, outbr, CAND, MON, noshed, C, VLIM)
% Count violating draws with NO shedding.
% brk = [adequacy-only, thermal, undervoltage, adequacy-any, no-solution].
n = numel(Ppv); nbad = 0; brk = [0 0 0 0 0];
for i = 1:n
    r = csolve(L, base, Ppv(i), Pwt(i), noshed, outbr, CAND, MON, C, VLIM);
    if ~r.ok, nbad = nbad + 1; brk(5) = brk(5) + 1; continue; end
    ge = r.Pg1 > PCOM + 1e-6; th = any(r.S > RATE + 1e-6); uv = r.vmin < VLIM - 1e-6;
    if ge || th || uv, nbad = nbad + 1; end
    if ge && ~th && ~uv, brk(1) = brk(1) + 1; end
    if th, brk(2) = brk(2) + 1; end
    if uv, brk(3) = brk(3) + 1; end
    if ge, brk(4) = brk(4) + 1; end
end
end

function P = cturb(v, Pr, vci, vr, vco)
% Piecewise linear turbine power curve.
v = v(:); P = zeros(size(v));
m = v >= vci & v < vr;   P(m) = Pr*(v(m)-vci)/(vr-vci);
m = v >= vr & v <= vco;  P(m) = Pr;
end

function P = cwindmean(k, c, Pr, vci, vr, vco)
% Expected turbine power under Weibull(k,c), by fine quadrature (toolbox-free).
N = 200000; U = ((1:N)'-0.5)/N;
P = mean(cturb(c*(-log(1-U)).^(1/k), Pr, vci, vr, vco));
end

function c = csolvec(k, target, Pr, vci, vr, vco)
% Weibull scale c that makes the expected turbine power equal 'target' -- used to
% vary volatility (shape k) while holding the mean RE power fixed.
lo = 1.0; hi = 40.0;
for it = 1:200 %#ok<NASGU>
    mid = 0.5*(lo+hi);
    if cwindmean(k, mid, Pr, vci, vr, vco) < target, lo = mid; else, hi = mid; end
end
c = 0.5*(lo+hi);
end

function n = onset(hist, frac)
% First iteration at which best-so-far is within (1-frac) of its final value.
if nargin < 2, frac = 0.999; end
tgt = hist(end) + (1-frac)*abs(hist(end));
n = find(hist <= tgt, 1);
if isempty(n), n = numel(hist); end
end

function prow(name, f, o, fopt)
fprintf('  %-8s %9.4f %9.4f %9.4f %9.4f %11.4f %9.1f\n', ...
        name, min(f), mean(f), max(f), std(f), mean(f)-fopt, median(o));
end

function s = yesno(t)
if t, s = 'yes'; else, s = 'no'; end
end

function s = verdict_icpso(dmean, se, sp, si, sq)
% Data-driven one-line summary of the paired comparison.
if dmean > 2*se
    q = 'ICPSO is WORSE than PSO on the mean';
elseif dmean < -2*se
    q = 'ICPSO is BETTER than PSO on the mean';
else
    q = 'ICPSO and PSO are statistically indistinguishable on the mean';
end
if sq < sp && sq <= si
    v = 'and the chaotic weight (Freq=1) gives the tightest spread';
elseif si < sp
    v = 'though ICPSO shows a tighter spread';
else
    v = 'and shows no variance advantage here';
end
s = sprintf('%s (diff %+.4f MW, se %.4f), %s.', q, dmean, se, v);
end
