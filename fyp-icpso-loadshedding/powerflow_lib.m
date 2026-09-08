function L = powerflow_lib()
% POWERFLOW_LIB  Shared helpers for replicating Adetona et al. [31] and [32] on IEEE-14.
%
%   L = powerflow_lib();  then call e.g.  mpc = L.mpc0();  r = L.solve(mpc);
%
% Self-contained: loads case14 straight from MATPOWER and depends on no .mat file.
% Every routine here was validated against an independent NumPy implementation of
% Newton-Raphson before being written; see README.md for the checkpoints.
%
% References
%   [31] Adetona et al. (2022), Improved Chaotic PSO with wingbeat frequency,
%        Acta Marisiensis 19(2).  Solves optimal reactive power dispatch.
%   [32] Adetona et al. (2025), PSO-Based Maximum Constrained Load-Shedding with
%        Critical Severity Index, J. Engineering (Univ. Baghdad) 31(9).  THE BENCHMARK.

L = struct();

% ---- data and power flow -------------------------------------------------
L.C          = constants();
L.mpc0       = @mpc0;
L.opt        = @pfopt;
L.solve      = @solve;
L.mainisland = @mainisland;
L.flows      = @flows;
L.findbr     = @findbr;
L.loss       = @loss;
L.busV       = @busV;

% ---- [32] machinery ------------------------------------------------------
L.loadbuses  = @loadbuses;
L.dV         = @dV;
L.eq3_pso    = @eq3_pso;
L.n1screen   = @n1screen;
L.applyshed  = @applyshed;
L.shed_pso   = @shed_pso;
L.eq10_pso   = @eq10_pso;
L.eq10_floor = @eq10_floor;

% ---- [31] machinery ------------------------------------------------------
L.wingbeat   = @wingbeat;
L.orpd       = @orpd_setup;
L.pso        = @pso;
L.icpso      = @icpso;
L.polish     = @polish;

% ---- reporting -----------------------------------------------------------
L.cmp        = @cmp;
L.banner     = @banner;
end

% =========================================================================
%  DATA AND POWER FLOW
% =========================================================================

function C = constants()
% MATPOWER column indices, hard-coded so nothing depends on define_constants.
% Cached, for the same reason as pfopt: this is called on every load flow.
persistent C_cached
if ~isempty(C_cached), C = C_cached; return; end
C.BUS_I=1; C.BUS_TYPE=2; C.PD=3; C.QD=4; C.GS=5; C.BS=6; C.BUS_AREA=7;
C.VM=8; C.VA=9; C.BASE_KV=10; C.ZONE=11; C.VMAX=12; C.VMIN=13;
C.GEN_BUS=1; C.PG=2; C.QG=3; C.QMAX=4; C.QMIN=5; C.VG=6; C.MBASE=7;
C.GEN_STATUS=8; C.PMAX=9; C.PMIN=10;
C.F_BUS=1; C.T_BUS=2; C.BR_R=3; C.BR_X=4; C.BR_B=5; C.RATE_A=6;
C.RATE_B=7; C.RATE_C=8; C.TAP=9; C.SHIFT=10; C.BR_STATUS=11;
C.ANGMIN=12; C.ANGMAX=13; C.PF=14; C.QF=15; C.PT=16; C.QT=17;
C.PQ=1; C.PV=2; C.REF=3;
% [31] Table 4 control set for IEEE-14: 5 gen voltages, 3 taps, 1 shunt.
C.VG_BUS  = [1 2 3 6 8];
C.TAP_BR  = [8 9 10];        % branches 4-7, 4-9, 5-6
C.SHUNT_BUS = 9;
C_cached = C;
end

function mpc = mpc0()
% Unmodified MATPOWER case14.  This is the single source of network data.
mpc = loadcase('case14');
end

function o = pfopt()
% [32] Eq (20): convergence test max(|dP|,|dQ|) <= 1e-6 pu.
% MATPOWER's default is 1e-8, so this MUST be set explicitly to match the paper.
%
% Cached in a persistent, because mpoption builds and validates a large struct
% and Stage B calls the load flow a few hundred thousand times.  The options
% never change, so this is a pure speed-up with no effect on any result.
% pf.enforce_q_lims is deliberately LEFT OFF here and replaced by the
% PV-bus-only enforcement in qlim_pv below.  The reason is specific and
% important: in case14 the slack generator's own reactive limits are
% Qmin = 0, Qmax = 10 Mvar, while the solved base case needs Qg1 = -16.55
% Mvar.  MATPOWER's built-in enforcement therefore judges the SLACK to be in
% violation, switches it off, folds its 232.4 MW into the bus-1 load, converts
% bus 1 to PQ and promotes bus 2 to slack.  That returns 13.2506 MW of loss
% instead of 13.3933 MW, and it silently replaces the fixed slack voltage that
% Adetona's formulation assumes.  Both papers report 13.393, so the reading
% that reproduces them is the one that leaves the slack alone.
persistent o_cached
if isempty(o_cached)
    o_cached = mpoption('pf.tol', 1e-6, 'pf.enforce_q_lims', 0, ...
                        'verbose', 0, 'out.all', 0);
end
o = o_cached;
end

function r = solve(mpc)
% Newton-Raphson load flow with [32]'s tolerance.  Never errors: returns a
% struct whose .ok field says whether a usable solution was found.
%
% "ok = false" carries real physical meaning here, not just numerical trouble:
% under the double outage of lines 1 and 14 the IEEE-14 power flow HAS NO
% SOLUTION, and this is what reports that honestly.
C = constants();
r = struct('ok', false, 'res', [], 'it', -1, 'Vm', [], 'Va', [], ...
           'loss', NaN, 'vmin', NaN, 'vmax', NaN, 'why', '', 'qpinned', []);

% An out-of-service branch that leaves a bus with no connection makes the
% Jacobian singular.  Detect that first so the failure is labelled, not thrown.
[iso, msg] = isolated_buses(mpc);
if ~isempty(iso)
    r.why = sprintf('bus(es) %s isolated -> %s', mat2str(iso), msg);
end

ws = warning('off', 'all');
try
    [res, ok, qp] = qlim_pv(mpc);
catch err
    warning(ws); r.why = ['runpf errored: ' err.message]; return;
end
warning(ws);
r.qpinned = qp;

r.res = res; r.Vm = res.bus(:,C.VM); r.Va = res.bus(:,C.VA);
if isfield(res, 'iterations'), r.it = res.iterations;
elseif isfield(res,'raw') && isfield(res.raw,'iterations'), r.it = res.raw.iterations;
end

% Guard against a "converged" flag on a physically absurd point.
bad = ~ok || any(~isfinite(r.Vm)) || max(r.Vm) > 2.0 || min(r.Vm) < 0.1;
if bad
    if isempty(r.why), r.why = 'Newton-Raphson did not converge to a feasible point'; end
    return;
end
r.ok    = true;
r.vmin  = min(r.Vm);
r.vmax  = max(r.Vm);
r.loss  = loss(res);
end

function [res, ok, pinned] = qlim_pv(mpc)
% Reactive limit enforcement on the PV-bus generators ONLY, never on the slack.
%
% Each pass solves the load flow, then any PV-bus generator whose computed Qg
% has run past a limit is pinned at that limit and its bus reclassified PQ, so
% its voltage becomes free.  Iterated to a fixed point, at most a few passes.
%
% Why not mpoption('pf.enforce_q_lims',1): see the note in pfopt.  MATPOWER's
% version also judges the slack, and on case14 it moves the slack off bus 1,
% which changes the answer and the boundary condition the papers assume.
% Leaving the slack out is also the physically right reading -- the slack's Qg
% is the network's residual, not a schedulable quantity.
C = constants();
work = mpc;  pinned = [];
for pass = 1:12
    [res, ok] = runpf(work, pfopt());
    if ~ok, return; end
    ng   = size(res.gen,1);
    move = false(ng,1);  lim = zeros(ng,1);
    for g = 1:ng
        if res.gen(g,C.GEN_STATUS) <= 0, continue; end
        ib = find(work.bus(:,C.BUS_I) == res.gen(g,C.GEN_BUS), 1);
        if isempty(ib) || work.bus(ib,C.BUS_TYPE) ~= C.PV, continue; end
        if     res.gen(g,C.QG) > res.gen(g,C.QMAX) + 1e-6
            move(g) = true;  lim(g) = res.gen(g,C.QMAX);
        elseif res.gen(g,C.QG) < res.gen(g,C.QMIN) - 1e-6
            move(g) = true;  lim(g) = res.gen(g,C.QMIN);
        end
    end
    if ~any(move), return; end
    for g = find(move)'
        ib = find(work.bus(:,C.BUS_I) == res.gen(g,C.GEN_BUS), 1);
        work.bus(ib,C.BUS_TYPE) = C.PQ;   % voltage now free
        work.gen(g,C.QG)        = lim(g); % reactive output pinned at the limit
        pinned = [pinned; res.gen(g,C.GEN_BUS) lim(g)]; %#ok<AGROW>
    end
end
end

function [iso, msg] = isolated_buses(mpc)
% Buses left with no in-service branch, and buses whose removal splits the graph.
C = constants();
nb  = size(mpc.bus,1);
on  = mpc.branch(:,C.BR_STATUS) > 0;
f   = mpc.branch(on, C.F_BUS); t = mpc.branch(on, C.T_BUS);
deg = accumarray([f;t], 1, [max(mpc.bus(:,C.BUS_I)) 1]);
iso = mpc.bus(deg(mpc.bus(:,C.BUS_I)) == 0, C.BUS_I);
msg = 'no in-service branch';
if isempty(iso)
    % connectivity sweep from the slack
    ref = mpc.bus(mpc.bus(:,C.BUS_TYPE)==C.REF, C.BUS_I);
    seen = false(max(mpc.bus(:,C.BUS_I)),1); stack = ref(1); seen(ref(1)) = true;
    while ~isempty(stack)
        b = stack(end); stack(end) = [];
        nb2 = [t(f==b); f(t==b)];
        for k = 1:numel(nb2)
            if ~seen(nb2(k)), seen(nb2(k)) = true; stack(end+1) = nb2(k); end %#ok<AGROW>
        end
    end
    iso = mpc.bus(~seen(mpc.bus(:,C.BUS_I)), C.BUS_I);
    msg = 'not connected to the slack bus (islanded)';
end
end

function k = findbr(mpc, fb, tb)
% Locate a branch by its end buses rather than by row index.  Necessary because
% MAINISLAND deletes rows, which shifts the numbering: "line 2" must stay the
% branch from bus 1 to bus 5 regardless of what row it now occupies.
C = constants();
k = find((mpc.branch(:,C.F_BUS)==fb & mpc.branch(:,C.T_BUS)==tb) | ...
         (mpc.branch(:,C.F_BUS)==tb & mpc.branch(:,C.T_BUS)==fb), 1);
end

function [m, dropped] = mainisland(mpc)
% Keep only the island that contains the slack bus.
%
% Needed for the outage of line 14 (7-8): bus 8 is a generator bus with no
% load and no other connection, so removing line 14 leaves it as an island of
% its own.  Deleting it is the physically correct treatment -- and the reason
% line 14 is critical is exactly that the main system loses bus 8's reactive
% support (Qmax 24 Mvar, 17.4 Mvar in the base case) when that happens.
C = constants(); m = mpc; dropped = [];
on = m.branch(:,C.BR_STATUS) > 0;
f  = m.branch(on,C.F_BUS); t = m.branch(on,C.T_BUS);
ref = m.bus(m.bus(:,C.BUS_TYPE)==C.REF, C.BUS_I); ref = ref(1);
seen = false(max(m.bus(:,C.BUS_I)),1); seen(ref) = true; stack = ref;
while ~isempty(stack)
    b = stack(end); stack(end) = [];
    nbr = [t(f==b); f(t==b)];
    for k = 1:numel(nbr)
        if ~seen(nbr(k)), seen(nbr(k)) = true; stack(end+1) = nbr(k); end %#ok<AGROW>
    end
end
keep    = seen(m.bus(:,C.BUS_I));
dropped = m.bus(~keep, C.BUS_I);
if isempty(dropped), return; end

gkeep = ismember(m.gen(:,C.GEN_BUS), m.bus(keep,C.BUS_I));
if isfield(m,'gencost') && size(m.gencost,1) == size(m.gen,1)
    m.gencost = m.gencost(gkeep,:);
end
m.gen    = m.gen(gkeep,:);
bkeep    = ismember(m.branch(:,C.F_BUS), m.bus(keep,C.BUS_I)) & ...
           ismember(m.branch(:,C.T_BUS), m.bus(keep,C.BUS_I));
m.branch = m.branch(bkeep,:);
m.bus    = m.bus(keep,:);
end

function T = flows(res)
% Per-branch flows.  S_pu is the apparent power at the FROM end in per-unit on
% the 100 MVA base -- that is the quantity [32] quotes as "line 2 = 0.6271 pu".
C = constants();
n = size(res.branch,1);
T = struct('k', num2cell((1:n)'), 'f', num2cell(res.branch(:,C.F_BUS)), ...
           't', num2cell(res.branch(:,C.T_BUS)));
Sf = hypot(res.branch(:,C.PF), res.branch(:,C.QF));
St = hypot(res.branch(:,C.PT), res.branch(:,C.QT));
for k = 1:n
    T(k).Pf = res.branch(k,C.PF);  T(k).Qf = res.branch(k,C.QF);
    T(k).Smva = Sf(k);  T(k).S_pu = Sf(k)/res.baseMVA;
    T(k).Smva_to = St(k);
    T(k).on = res.branch(k,C.BR_STATUS) > 0;
    T(k).ploss = res.branch(k,C.PF) + res.branch(k,C.PT);
end
end

function p = loss(res)
% Total real power loss in MW.
C = constants();
on = res.gen(:,C.GEN_STATUS) > 0;
p  = sum(res.gen(on,C.PG)) - sum(res.bus(:,C.PD));
end

function V = busV(res)
% Complex bus voltages of a solved case, as a column vector in bus-row order.
C = constants();
V = res.bus(:,C.VM).*exp(1j*res.bus(:,C.VA)*pi/180);
end

% =========================================================================
%  [32] CRITICAL LINE IDENTIFICATION  (Eqs 3-4 and the N-1 screen)
% =========================================================================

function d = dV(res, mode)
% Kernel of [32] Eq (3): |V_from,k - V_to,k| per branch.
%
% The paper does not say whether this is a magnitude difference or a phasor
% difference, and the two give different answers, so both are provided:
%   'mag'    sum over 20 branches = 0.4179  -> E[f] at a random swarm = 0.2090
%   'phasor' sum over 20 branches = 1.1094  -> E[f] at a random swarm = 0.5547
% [32] reports a CSI fitness of "approximately 0.5", which identifies the
% phasor form evaluated before the search has moved.  Default is 'phasor'.
if nargin < 2, mode = 'phasor'; end
C = constants();
f = res.branch(:,C.F_BUS); t = res.branch(:,C.T_BUS);
Vm = res.bus(:,C.VM); Va = res.bus(:,C.VA)*pi/180;
switch lower(mode)
    case 'mag',    d = abs(Vm(f) - Vm(t));
    case 'phasor', V = Vm.*exp(1j*Va); d = abs(V(f) - V(t));
    otherwise, error('dV: mode must be ''mag'' or ''phasor''');
end
end

function out = eq3_pso(d, sense, NP, NIT, seed)
% [32] Eq (3):  f(x) = sum_k x_k * |dV_k|   subject only to  0 <= x_k <= 1.
%
% Implemented exactly as printed, which means it is DEGENERATE: the objective
% is linear in x over a box, so minimising drives every x_k to 0 (f = 0) and
% maximising drives every x_k to 1 (f = sum|dV_k|).  Either way no line is
% singled out.  The N-1 screen below is what actually reproduces [32]'s
% critical set {line 1, line 14}; this routine exists to document that the
% equation as published cannot.
if nargin < 2 || isempty(sense), sense = 'min'; end
if nargin < 3 || isempty(NP),    NP = 30;  end
if nargin < 4 || isempty(NIT),   NIT = 250; end
if nargin < 5 || isempty(seed),  seed = 1;  end
rng(seed, 'twister');
D = numel(d); d = d(:);
sg = 1; if strcmpi(sense,'max'), sg = -1; end
X = rand(NP,D); V = zeros(NP,D);
F = sg*(X*d); pb = X; pbF = F;
[gbF, gi] = min(pbF); gb = pb(gi,:);
hist = zeros(NIT,1);
for k = 1:NIT
    V = 0.729*(V + 2.05*rand(NP,D).*(pb-X) + 2.05*rand(NP,D).*(repmat(gb,NP,1)-X));
    V = max(min(V, 0.2), -0.2);
    X = max(min(X + V, 1), 0);
    F = sg*(X*d);
    up = F < pbF; pbF(up) = F(up); pb(up,:) = X(up,:);
    [m, gi] = min(pbF);
    if m < gbF, gbF = m; gb = pb(gi,:); end
    hist(k) = sg*gbF;
end
out = struct('x', gb, 'f', sg*gbF, 'hist', hist, ...
             'sum_dV', sum(d), 'E_random', sum(d)/2);
end

function S = n1screen(mpc)
% Outage each branch in turn and record whether the network still solves.
% This is the operative critical-line test.  On IEEE-14 it returns exactly
% {line 1 (1-2), line 14 (7-8)}, matching [32] -- 18 of the other 20 outages
% solve cleanly in 3 to 8 Newton-Raphson iterations.
C  = constants();
nl = size(mpc.branch,1);
S  = repmat(struct('k',0,'f',0,'t',0,'ok',false,'radial',false,'it',-1, ...
                   'vmin',NaN,'vmax',NaN,'loss',NaN,'sev',Inf,'why',''), nl, 1);
on  = mpc.branch(:,C.BR_STATUS) > 0;
fa  = mpc.branch(on,C.F_BUS); ta = mpc.branch(on,C.T_BUS);
deg = accumarray([fa;ta], 1, [max(mpc.bus(:,C.BUS_I)) 1]);
for k = 1:nl
    m = mpc; m.branch(k, C.BR_STATUS) = 0;
    fb = mpc.branch(k,C.F_BUS); tb = mpc.branch(k,C.T_BUS);
    S(k).k = k; S(k).f = fb; S(k).t = tb;
    S(k).radial = (deg(fb) == 1) || (deg(tb) == 1);
    r = solve(m);
    S(k).ok = r.ok; S(k).it = r.it; S(k).why = r.why;
    if r.ok
        S(k).vmin = r.vmin; S(k).vmax = r.vmax; S(k).loss = r.loss;
        S(k).sev  = sum(max(0, 0.95 - r.Vm).^2 + max(0, r.Vm - 1.05).^2);
    end
end
end

% =========================================================================
%  [32] TARGETED LOAD SHEDDING  (Eqs 5-9)
% =========================================================================

function b = loadbuses(mpc)
% [32] sheds across ALL load buses.  On IEEE-14 that is the 11 buses
% {2,3,4,5,6,9,10,11,12,13,14}.  Note bus 4 has Qd = -3.9 Mvar, which makes
% Eq (7)'s "0 <= Qshed_i <= Qd_i" empty there, so no Q can be shed at bus 4.
C = constants();
b = mpc.bus(mpc.bus(:,C.PD) ~= 0 | mpc.bus(:,C.QD) ~= 0, C.BUS_I);
end

function m = applyshed(mpc, buses, Pshed, Qshed)
% Eq (8): Pd_new = Pd - Pshed, Qd_new = Qd - Qshed, with Eq (9) Pd_new >= 0.
C = constants(); m = mpc;
for i = 1:numel(buses)
    r = find(m.bus(:,C.BUS_I) == buses(i));
    m.bus(r,C.PD) = max(0, m.bus(r,C.PD) - Pshed(i));
    if m.bus(r,C.QD) >= 0
        m.bus(r,C.QD) = max(0, m.bus(r,C.QD) - Qshed(i));
    end   % bus 4's negative Qd is a capacitive injection: leave it alone
end
end

function out = shed_pso(mpc, mu, NP, NIT, seed, alg)
% [32] Eqs (5)-(9).  Eq (5) minimises sum|Pd_i - Pshed_i|, which drives every
% Pshed_i UP toward Pd_i, so the objective MAXIMISES shedding and Eq (6)'s
% cap is what stops it.  The cap therefore binds by construction, which is
% exactly why [32]'s answer lands on 15% of demand to four significant
% figures.  (Section 3.3 of [32] says "minimize the load shedding", but the
% section title, Eq (5) and the published result all say maximise; the
% maximise reading is the only self-consistent one.)
%
% Benchmark: mu = 0.15 gives 38.8500 MW and 11.0250 Mvar,
%            against [32]'s published 38.8502 MW and 11.0250 Mvar.
%
% alg selects the swarm update, everything else held identical so the two are
% a properly paired comparison -- same seed, same initial swarm, same repair
% operator, same fitness, only the velocity/position rule differs:
%   'pso'    [32]'s standard constriction PSO                        (default)
%   'icpso'  [31]'s chaotic weight + per-particle wingbeat frequency
%   'freq1'  [31]'s chaotic weight with Freq = 1, i.e. the wingbeat term
%            switched off.  This isolates which half of ICPSO does the work.
if nargin < 2 || isempty(mu),   mu = 0.15; end
if nargin < 3 || isempty(NP),   NP = 30;   end
if nargin < 4 || isempty(NIT),  NIT = 250; end
if nargin < 5 || isempty(seed), seed = 1;  end
if nargin < 6 || isempty(alg),  alg = 'pso'; end
alg = lower(alg);
if ~any(strcmp(alg, {'pso','icpso','freq1'}))
    error('shed_pso: alg must be ''pso'', ''icpso'' or ''freq1''');
end
rng(seed, 'twister');
C  = constants();
bs = loadbuses(mpc); n = numel(bs);
ix = arrayfun(@(b) find(mpc.bus(:,C.BUS_I)==b), bs);
Pd = mpc.bus(ix, C.PD); Qd = mpc.bus(ix, C.QD);
capP = mu*sum(mpc.bus(:,C.PD));  capQ = mu*sum(mpc.bus(:,C.QD));
hi = [max(Pd,0); max(Qd,0)];  lo = zeros(2*n,1);   % Eq (7)
D  = 2*n;

rep = @(x) repair(x, n, Pd, Qd, capP, capQ);
fit = @(x) sum(abs(Pd - x(1:n))) + sum(abs(Qd - x(n+1:end)));   % Eq (5)

W = wingbeat();
X = zeros(NP,D);
for i = 1:NP, X(i,:) = rep(lo + rand(D,1).*(hi-lo))'; end
V = zeros(NP,D); vmax = 0.2*(hi-lo)';
F = zeros(NP,1); for i = 1:NP, F(i) = fit(X(i,:)'); end
pb = X; pbF = F; [gbF, gi] = min(pbF); gb = pb(gi,:);
hist = zeros(NIT,1);  w = 0.4;
for k = 1:NIT
    switch alg
        case 'pso'
            w  = 1.0;                                        % [32] gives no schedule
            Fr = ones(NP,1);
        case 'icpso'
            w  = chaotic(w);                                 % [31] Eqs (17)-(18)
            Fr = W.freq(rand(NP,1));                         % [31] Eq (24), per particle
        case 'freq1'
            w  = chaotic(w);
            Fr = ones(NP,1);
    end
    A = 2.05*rand(NP,D).*(pb-X) + 2.05*rand(NP,D).*(repmat(gb,NP,1)-X);
    V = 0.729*(w*V + bsxfun(@times, A, Fr));                 % [31] Eq (22)
    V = max(min(V, repmat(vmax,NP,1)), -repmat(vmax,NP,1));
    Xn = X + bsxfun(@rdivide, V, Fr);                        % [31] Eq (23)
    for i = 1:NP
        X(i,:) = rep(Xn(i,:)')';
        f = fit(X(i,:)');
        if f < pbF(i), pbF(i) = f; pb(i,:) = X(i,:); end
    end
    [m2, gi] = min(pbF);
    if m2 < gbF, gbF = m2; gb = pb(gi,:); end
    hist(k) = gbF;
end
out = struct('buses', bs, 'Pshed', gb(1:n)', 'Qshed', gb(n+1:end)', ...
             'Ptot', sum(gb(1:n)), 'Qtot', sum(gb(n+1:end)), ...
             'capP', capP, 'capQ', capQ, 'f', gbF, 'hist', hist, ...
             'Pd', Pd, 'Qd', Qd, 'alg', alg);
out.onset = find(hist <= hist(end) + 1e-6, 1);
end

function w = chaotic(w)
% Logistic map with mu = 4, [31] Eqs (17)-(18), with the two fixed points of
% the map (0 and 1) escaped so the weight cannot stick.
w = 4.0*w*(1.0-w);
if w <= 1e-12 || w >= 1-1e-12, w = rand(); end
end


function y = repair(x, n, Pd, Qd, capP, capQ)
% Enforce Eq (7) box bounds, then scale onto Eq (6)'s aggregate caps.
p = min(max(x(1:n),     0), max(Pd,0));
q = min(max(x(n+1:end), 0), max(Qd,0));
if sum(p) > capP && sum(p) > 0, p = p*capP/sum(p); end
if sum(q) > capQ && sum(q) > 0, q = q*capQ/sum(q); end
y = [p; q];
end

% =========================================================================
%  [32] PSO-WARM-STARTED LOAD FLOW  (Eqs 10-20)
% =========================================================================

function out = eq10_pso(mpc, NP, NIT, seed, free_slack, scope)
% [32] Eq (10) defines a particle as x = [|V_1| ang d_1, ..., |V_Nb| ang d_Nb]
% -- ALL buses, INCLUDING the slack -- and minimises the summed squared power
% mismatch of Eqs (11)-(19).  Neither reading of that sum is a well-posed load
% flow, which is the finding that explains [32]'s contingency voltages:
%
%   scope = 'all'  sums the P and Q mismatch at every bus.  That needs a
%       SCHEDULED slack injection, which is the one quantity a load flow exists
%       to compute.  Feed it a stale slack P -- as happens the moment load is
%       shed -- and the objective has a non-zero floor no matter what voltages
%       are tried.  Use EQ10_FLOOR to measure that floor.
%
%   scope = 'nr'   sums only the equations Newton-Raphson enforces: P at the PV
%       and PQ buses, Q at the PQ buses.  On IEEE-14 that is 22 equations in 27
%       unknowns (2*Nb voltages less the angle reference), so the system is
%       UNDERDETERMINED by five and has infinitely many zero-residual
%       solutions.  Freeing the generator magnitudes also frees the generator
%       reactive outputs, so the search can settle on a physically pointless
%       point: on the base case it reaches a residual of 7e-5 with line 2 at
%       2.12 pu against the true 0.7561 pu.
%
% Either way, always report out.residual next to the voltages.  A residual that
% is not ~0 means the voltages do not satisfy Kirchhoff's laws, and that is how
% a slack voltage of 0.7214 pu can be printed without the method complaining.
if nargin < 2 || isempty(NP),   NP = 60;  end
if nargin < 3 || isempty(NIT),  NIT = 400; end
if nargin < 4 || isempty(seed), seed = 1;  end
if nargin < 5 || isempty(free_slack), free_slack = true; end
if nargin < 6 || isempty(scope), scope = 'all'; end
rng(seed, 'twister');
C = constants();
nb = size(mpc.bus,1);
[Ybus, ~, ~] = makeYbus(ext2int(mpc));
[Psp, Qsp] = specified(mpc);
[mP, mQ] = scopemask(mpc, scope);
vlo = 0.5; vhi = 1.15; dlo = -60*pi/180; dhi = 20*pi/180;
lo = [repmat(vlo,nb,1); repmat(dlo,nb,1)];
hi = [repmat(vhi,nb,1); repmat(dhi,nb,1)];
ref = find(mpc.bus(:,C.BUS_TYPE) == C.REF, 1);
if ~free_slack, lo(ref) = mpc.bus(ref,C.VM); hi(ref) = mpc.bus(ref,C.VM); end
lo(nb+ref) = 0; hi(nb+ref) = 0;          % angle reference
D = 2*nb;
fit = @(x) mismatch(x, nb, Ybus, Psp, Qsp, mP, mQ);

X = repmat(lo',NP,1) + rand(NP,D).*repmat((hi-lo)',NP,1);
V = zeros(NP,D); vmax = 0.2*(hi-lo)';
F = zeros(NP,1); for i = 1:NP, F(i) = fit(X(i,:)'); end
pb = X; pbF = F; [gbF, gi] = min(pbF); gb = pb(gi,:);
hist = zeros(NIT,1);
for k = 1:NIT
    V = 0.7*V + 1.5*rand(NP,D).*(pb-X) + 1.5*rand(NP,D).*(repmat(gb,NP,1)-X);
    V = max(min(V, repmat(vmax,NP,1)), -repmat(vmax,NP,1));
    X = max(min(X + V, repmat(hi',NP,1)), repmat(lo',NP,1));
    for i = 1:NP
        f = fit(X(i,:)');
        if f < pbF(i), pbF(i) = f; pb(i,:) = X(i,:); end
    end
    [m2, gi] = min(pbF);
    if m2 < gbF, gbF = m2; gb = pb(gi,:); end
    hist(k) = gbF;
end
Vm = gb(1:nb)'; Va = gb(nb+1:end)';
out = struct('Vm', Vm, 'Va', Va*180/pi, 'V', Vm.*exp(1j*Va), ...
             'residual', gbF, 'hist', hist, 'scope', scope, ...
             'neq', sum(mP)+sum(mQ), 'nunk', 2*nb-1);
out.S_pu = branch_S_pu(mpc, out.V);
end

function f = eq10_floor(mpc, V, scope)
% The lowest value [32]'s Eq (10) objective can take at a voltage set that is
% KNOWN to solve the load flow.  This is what separates "the specification is
% inconsistent" from "the optimiser failed", and it is worth reporting because
% the two call for opposite conclusions.
%
% On the base case:            scope 'all' floor = 4.6e-04,  scope 'nr' = 2e-20
% On the outage + 15% shed:    scope 'all' floor = 2.5e-01,  scope 'nr' = 1e-21
%
% So under the 'all' reading the objective cannot reach zero once load has been
% shed, because the scheduled slack injection is then stale by the amount shed.
if nargin < 3 || isempty(scope), scope = 'all'; end
[Ybus, ~, ~] = makeYbus(ext2int(mpc));
[Psp, Qsp] = specified(mpc);
[mP, mQ] = scopemask(mpc, scope);
S = V.*conj(Ybus*V);
f = sum((Psp(mP) - real(S(mP))).^2) + sum((Qsp(mQ) - imag(S(mQ))).^2);
end

function [Psp, Qsp] = specified(mpc)
% Scheduled net injection per bus, in per-unit, from the case data.
C = constants(); nb = size(mpc.bus,1);
Pd = mpc.bus(:,C.PD)/mpc.baseMVA;  Qd = mpc.bus(:,C.QD)/mpc.baseMVA;
Pg = zeros(nb,1); Qg = zeros(nb,1);
for g = 1:size(mpc.gen,1)
    if mpc.gen(g,C.GEN_STATUS) > 0
        b = find(mpc.bus(:,C.BUS_I) == mpc.gen(g,C.GEN_BUS));
        Pg(b) = Pg(b) + mpc.gen(g,C.PG)/mpc.baseMVA;
        Qg(b) = Qg(b) + mpc.gen(g,C.QG)/mpc.baseMVA;
    end
end
Psp = Pg - Pd; Qsp = Qg - Qd;
end

function [mP, mQ] = scopemask(mpc, scope)
C = constants();
switch lower(scope)
    case 'all'
        mP = true(size(mpc.bus,1),1);  mQ = mP;
    case 'nr'
        mP = mpc.bus(:,C.BUS_TYPE) ~= C.REF;
        mQ = mpc.bus(:,C.BUS_TYPE) == C.PQ;
    otherwise
        error('eq10: scope must be ''all'' or ''nr''');
end
end

function e = mismatch(x, nb, Ybus, Psp, Qsp, mP, mQ)
V = x(1:nb).*exp(1j*x(nb+1:end));
S = V.*conj(Ybus*V);
e = sum((Psp(mP) - real(S(mP))).^2) + sum((Qsp(mQ) - imag(S(mQ))).^2);
end


function S = branch_S_pu(mpc, V)
% From-end apparent power per branch, in per-unit, for an arbitrary voltage set.
C = constants(); nl = size(mpc.branch,1); S = zeros(nl,1);
for k = 1:nl
    if mpc.branch(k,C.BR_STATUS) == 0, S(k) = NaN; continue; end
    f = find(mpc.bus(:,C.BUS_I)==mpc.branch(k,C.F_BUS));
    t = find(mpc.bus(:,C.BUS_I)==mpc.branch(k,C.T_BUS));
    ys = 1/(mpc.branch(k,C.BR_R) + 1j*mpc.branch(k,C.BR_X));
    bc = mpc.branch(k,C.BR_B);
    tap = mpc.branch(k,C.TAP); if tap == 0, tap = 1; end
    tap = tap*exp(1j*pi/180*mpc.branch(k,C.SHIFT));
    Vf = V(f)/tap;
    If = (Vf - V(t))*ys + Vf*(1j*bc/2);
    S(k) = abs(Vf*conj(If));
end
end

% =========================================================================
%  [31] ICPSO
% =========================================================================

function W = wingbeat()
% [31] section 3.1: the swarm is Pteropodidae Rousettus (a bat), every particle
% sharing the same wing geometry, so Freq(l) varies only through beta.
W.m_tot = 0.104; W.m_wng = 0.0130; W.b = 0.530; W.s = 0.0465;
W.rho = 1.21;    W.g = 9.81;
fref = @(m) 1.08*( m^(1/3) * W.g^0.5 * W.b^(-1) * W.s^(-0.25) * W.rho^(-1/3) );
W.f_min = fref(W.m_wng);      % 3.0328 Hz
W.f_max = fref(W.m_tot);      % 6.0655 Hz  (exactly 2 x f_min, since m_tot/m_wng = 8)
W.fref  = fref;

% Eq (24) is printed as  Freq = f_min + beta*(f_min - f_max), which DECREASES
% from f_min to exactly zero over beta in [0,1] -- so it spans (0, 3.0328]
% rather than the [3.0328, 6.0655] band [31]'s own text states, and Eq (23)
% divides by it, making the position step unbounded as beta -> 1.  Confirmed
% against the page image, so it is a genuine typesetting error and not a
% text-extraction artefact.  The standard bat-algorithm form is used instead:
W.freq_printed = @(beta) W.f_min + beta.*(W.f_min - W.f_max);   % do not use
W.freq         = @(beta) W.f_min + beta.*(W.f_max - W.f_min);   % [3.0328, 6.0655]

% Eq (21), m_wng = 0.112*m_tot^0.11, returns 0.0873 kg, which is 6.7 times the
% 0.0130 kg stated in the text.  [31] used the measured value from its ref [52],
% so Eq (21) must NOT be used to regenerate m_wng.
W.eq21 = @(mt) 0.112*mt^0.11;
end

function O = orpd_setup(mpc)
% [31]'s IEEE-14 problem: minimise real power loss over 9 control variables.
% Limits are [31] Table 2: V_G 0.95-1.10, tap 0.90-1.10, Q_C 0-20 Mvar.
%
% Constraint handling follows [31] Steps 6, 7 and 10 literally: the penalty is
% on CONTROL variables only, and Step 10 already clamps those to their limits,
% so the penalty never fires and there is NO load-bus voltage constraint.
% That reading is not a liberty -- it is the only one that reproduces [31]'s
% own published PSO figure (12.2777 MW here against 12.275 MW published).
% Imposing 0.94-1.06 instead gives 12.3411 MW and 0.95-1.05 gives 12.4433 MW,
% neither of which matches; and [31]'s own Table 4 ICPSO point drives load-bus
% voltages up to 1.0869 pu, which confirms no such band was enforced.
C = constants();
O.mpc = mpc;
O.lo  = [0.95*ones(1,5), 0.90*ones(1,3),  0];
O.hi  = [1.10*ones(1,5), 1.10*ones(1,3), 20];
O.base = [1.060 1.045 1.010 1.070 1.090, 0.978 0.969 0.932, 19.00];
O.names = {'Vg1','Vg2','Vg3','Vg6','Vg8','T4-7','T4-9','T5-6','Qc9'};
O.apply = @(x) apply_controls(mpc, x);
O.fit   = @(x) orpd_fitness(mpc, x);
end

function m = apply_controls(mpc, x)
C = constants(); m = mpc;
for k = 1:numel(C.VG_BUS)
    b = C.VG_BUS(k);
    g = find(m.gen(:,C.GEN_BUS) == b);
    if ~isempty(g), m.gen(g, C.VG) = x(k); end
    r = find(m.bus(:,C.BUS_I) == b);
    if m.bus(r,C.BUS_TYPE) ~= C.PQ, m.bus(r,C.VM) = x(k); end
end
for k = 1:numel(C.TAP_BR)
    m.branch(C.TAP_BR(k), C.TAP) = x(5+k);
end
m.bus(m.bus(:,C.BUS_I)==C.SHUNT_BUS, C.BS) = x(9);
end

function f = orpd_fitness(mpc, x)
C = constants();
m = apply_controls(mpc, x);
r = solve(m);
if ~r.ok, f = 1e6; return; end
f = r.loss;
% Step 7 penalty on generator reactive limits only (control-variable side).
for g = 1:size(r.res.gen,1)
    if r.res.gen(g,C.GEN_STATUS) > 0
        q = r.res.gen(g,C.QG);
        f = f + (max(0, r.res.gen(g,C.QMIN)-q)^2 + max(0, q-r.res.gen(g,C.QMAX))^2)/100;
    end
end
end

function out = pso(fit, lo, hi, o)
% Constriction PSO with the linearly-decreasing inertia weight of [31] Eq (16).
o = defaults(o, numel(lo));
rng(o.seed, 'twister');
[X,V,pb,pbF,gb,gbF] = swarm_init(fit, lo, hi, o);
umax = 0.5*(hi(:)'-lo(:)');
hist = zeros(o.kmax,1); D = numel(lo);
for k = 1:o.kmax
    w = o.wmax - (o.wmax-o.wmin)*(k-1)/max(o.kmax-1,1);
    A = o.c1*rand(o.np,D).*(pb-X) + o.c2*rand(o.np,D).*(repmat(gb,o.np,1)-X);
    V = 0.729*(w*V + A);
    V = clampr(V, umax);
    X = clampb(X + V, lo, hi);
    [pb,pbF,gb,gbF] = swarm_update(fit, X, pb, pbF, gb, gbF);
    hist(k) = gbF;
end
out = struct('x', gb, 'f', gbF, 'hist', hist);
end

function out = icpso(fit, lo, hi, o)
% [31] ICPSO: chaotic inertia weight (Eqs 17-19) plus multiplicative wingbeat
% frequency (Eqs 22-24).
%
% What the wingbeat term actually does.  Substituting Eq (22) into Eq (23):
%     X(k+1) = X(k) + 0.729*A + 0.729*w*U(k)/Freq
% The attraction term A enters at full strength because Freq cancels, while
% the MOMENTUM term is divided by Freq.  With Freq in [3.03, 6.07] that
% removes most of the inertia, so ICPSO is a momentum-damped PSO.  Setting
% Freq = 1 collapses it to plain chaotic-weight constriction PSO, which is a
% useful control experiment (o.fixed_freq = 1).
%
% Freq is drawn PER PARTICLE per iteration: [31]'s Step 11 loops l over the N
% particles with steps 9 and 10 inside, and the notation Freq(l) carries the
% particle index.  Set o.per_particle = false for one shared draw per iteration.
o = defaults(o, numel(lo));
rng(o.seed, 'twister');
W = wingbeat();
[X,V,pb,pbF,gb,gbF] = swarm_init(fit, lo, hi, o);
umax = 0.5*(hi(:)'-lo(:)');
hist = zeros(o.kmax,1); D = numel(lo); w = o.w0;
for k = 1:o.kmax
    w = 4.0*w*(1.0-w);                                  % Eq (17)-(18), mu = 4
    if w <= 1e-12 || w >= 1-1e-12, w = rand(); end       % escape the fixed points
    if ~isempty(o.fixed_freq)
        Fr = o.fixed_freq;
    elseif o.per_particle
        Fr = W.freq(rand(o.np,1));                       % Eq (24), per particle
    else
        Fr = W.freq(rand());
    end
    A = o.c1*rand(o.np,D).*(pb-X) + o.c2*rand(o.np,D).*(repmat(gb,o.np,1)-X);
    V = 0.729*(w*V + bsxfun(@times, A, Fr));                              % Eq (22)
    V = clampr(V, umax);
    X = clampb(X + bsxfun(@rdivide, V, Fr), lo, hi);                      % Eq (23)
    [pb,pbF,gb,gbF] = swarm_update(fit, X, pb, pbF, gb, gbF);
    hist(k) = gbF;
end
out = struct('x', gb, 'f', gbF, 'hist', hist, 'f_min', W.f_min, 'f_max', W.f_max);
end

function o = defaults(o, D)
if ~isstruct(o), o = struct(); end
d = struct('np',30, 'kmax',250, 'seed',1, 'c1',2.05, 'c2',2.05, ...
           'wmin',0.4, 'wmax',0.9, 'w0',0.4, 'per_particle',true, ...
           'fixed_freq',[], 'seed_point',[]);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(o, fn{i}), o.(fn{i}) = d.(fn{i}); end
end
o.D = D;
end

function [X,V,pb,pbF,gb,gbF] = swarm_init(fit, lo, hi, o)
D = numel(lo); lo = lo(:)'; hi = hi(:)';
X = repmat(lo,o.np,1) + rand(o.np,D).*repmat(hi-lo,o.np,1);
if ~isempty(o.seed_point), X(1,:) = o.seed_point(:)'; end
V = zeros(o.np,D);
o.umax = 0.5*(hi-lo);
pbF = zeros(o.np,1);
for i = 1:o.np, pbF(i) = fit(X(i,:)); end
pb = X; [gbF, gi] = min(pbF); gb = pb(gi,:);
end

function [pb,pbF,gb,gbF] = swarm_update(fit, X, pb, pbF, gb, gbF)
for i = 1:size(X,1)
    f = fit(X(i,:));
    if f < pbF(i), pbF(i) = f; pb(i,:) = X(i,:); end
end
[m, gi] = min(pbF);
if m < gbF, gbF = m; gb = pb(gi,:); end
end

function V = clampr(V, umax)
if isempty(umax), return; end
U = repmat(umax(:)', size(V,1), 1);
V = max(min(V, U), -U);
end

function X = clampb(X, lo, hi)
X = max(min(X, repmat(hi(:)', size(X,1), 1)), repmat(lo(:)', size(X,1), 1));
end

function [x, f, ev] = polish(fit, x0, lo, hi, step0, tol)
% Compass (pattern) search with a shrinking step, respecting the box [lo,hi].
% Deterministic, derivative-free and written out in full so it needs NO
% toolbox -- neither Global Optimization nor Optimization is assumed.
%
% Its purpose is to establish the true optimum of [31]'s ORPD model, so that
% published figures can be measured against the best the model can actually
% attain rather than against another heuristic's output.  Started from a good
% PSO point, or from several random points, it settles on 12.2777 MW.
if nargin < 5 || isempty(step0), step0 = 0.25;  end
if nargin < 6 || isempty(tol),   tol   = 1e-7;  end
lo = lo(:)'; hi = hi(:)';
x  = min(max(x0(:)', lo), hi);
f  = fit(x); ev = 1;
D  = numel(x); s = step0*(hi-lo);
while max(s./(hi-lo)) > tol && ev < 20000
    improved = false;
    for j = 1:D
        for sg = [1 -1]
            t = x; t(j) = min(max(t(j) + sg*s(j), lo(j)), hi(j));
            if t(j) == x(j), continue; end
            ft = fit(t); ev = ev + 1;
            if ft < f - 1e-12, x = t; f = ft; improved = true; break; end
        end
    end
    if ~improved, s = s*0.5; end
end
end

% =========================================================================
%  REPORTING
% =========================================================================

function banner(txt)
fprintf('\n%s\n%s\n', txt, repmat('=', 1, numel(txt)));
end

function cmp(label, computed, published, unit, tol)
% One line of a computed-versus-published comparison table.
if nargin < 4, unit = ''; end
if nargin < 5 || isempty(tol), tol = 5e-3; end
if isnan(published)
    fprintf('  %-34s %12.4f %12s  %s\n', label, computed, '(none)', unit);
    return;
end
d = computed - published;
if abs(d) <= tol*max(1,abs(published)), mark = 'MATCH'; else mark = sprintf('%+.4f', d); end
fprintf('  %-34s %12.4f %12.4f  %-8s %s\n', label, computed, published, unit, mark);
end
