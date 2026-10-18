// Exact optimisation of the stability scores over all N! gene orderings by the
// dynamic program over subsets of genes described in the Methods section of
// Paper II. Vendored from the Paper II reproducibility bundle
// (supplementary/code/core/mu_dp_exact.js, commit 681aeee). Levels n = 0..N-1.
//
// An ordering is the same object as a maximal chain of gene sets
//   {} = S_0 < S_1 < ... < S_N = [N],
// and the scores are additive along the chain:
//   mu_E  is a sum of node weights   E(S) p^(N-|S|),       1 <= |S| <= N-1,
//   mu_A  is a sum of edge weights   A(S,g) p^(N-|S|),     S -> S + {g}, 0 <= |S| <= N-1,
//   mu_I  = N p^N - mu_E - mu_A      (control identity).
// E(S) counts the assignments u of S at which S is not autonomous (the next
// state of the genes in S is not constant over the configurations extending u).
// A(S,g) counts the assignments at which S is autonomous and, in addition, the
// next state of g is constant. The program returns the minimum of mu_E, the
// maximum of mu_A, the number of orderings attaining each, and the orderings.
//
// Boolean networks (p = 2). Configuration c encodes gene i in bit i.
// Input: a CSV with header config_from,config_to and 2^N rows.
//
// Usage:
//   node mu_dp_exact.js <transitions.csv> [genes.csv] [--list-muA]
//   node mu_dp_exact.js <transitions.csv> [genes.csv] --eval 12,3,4,8,5,1,7,2,10,11,9,6,0
//   node mu_dp_exact.js <transitions.csv> --stats   (exact mean and standard deviation of
//                                                  mu_E over all N! orderings)
//
// Examples (from the padic-grn-core/ folder):
//   node tools/mu_dp_exact.js examples/transitions_athaliana_N13_p2.csv
//   node tools/mu_dp_exact.js examples/transitions_celegans_N8_p2.csv
//   node tools/mu_dp_exact.js examples/transitions_drosophila_N6_p2.csv
'use strict';
const fs = require('fs');

const args = process.argv.slice(2);
if (args.length < 1) {
  console.error('usage: node mu_dp_exact.js <transitions.csv> [genes.csv] [--list-muA] [--eval i0,i1,...]');
  process.exit(1);
}
const csvPath = args[0];
const namesPath = args[1] && !args[1].startsWith('--') ? args[1] : null;
const listMuA = args.includes('--list-muA');
const wantStats = args.includes('--stats');
const evalIdx = args.indexOf('--eval');
const evalOrder = evalIdx >= 0 ? args[evalIdx + 1].split(',').map(Number) : null;

// ---------------------------------------------------------------- input
const rows = fs.readFileSync(csvPath, 'utf8').trim().split(/\r?\n/).slice(1);
const SZ = rows.length;
const N = Math.round(Math.log2(SZ));
if (1 << N !== SZ) throw new Error(`expected 2^N rows, got ${SZ}`);
const T = new Int32Array(SZ);
for (const r of rows) { const [a, b] = r.split(',').map(Number); T[a] = b; }

let NAMES = [...Array(N).keys()].map(String);
if (namesPath) {
  const nr = fs.readFileSync(namesPath, 'utf8').trim().split(/\r?\n/).slice(1);
  for (const r of nr) { const [i, nm] = r.split(','); NAMES[+i] = nm; }
}
const label = o => o.map(g => NAMES[g]).join(' > ');
const pc = x => { let n = 0; while (x) { n += x & 1; x >>>= 1; } return n; };
const P = k => 2 ** k;

// ------------------------------------------- node and edge weights
// For every set S, group the configurations by their restriction to S and
// record the OR and AND of the images over each group. A bit that is equal in
// OR and AND is constant on the group.
const wE = new Float64Array(SZ);          // node weight of S for mu_E
const wA = new Float64Array(SZ * N);      // edge weight of S -> S+{g} for mu_A
const orA = new Int32Array(SZ), andA = new Int32Array(SZ), stamp = new Int32Array(SZ).fill(-1);
for (let S = 0; S < SZ - 1; S++) {
  const keys = [];
  for (let c = 0; c < SZ; c++) {
    const k = c & S, im = T[c];
    if (stamp[k] !== S) { stamp[k] = S; orA[k] = im; andA[k] = im; keys.push(k); }
    else { orA[k] |= im; andA[k] &= im; }
  }
  const n = pc(S), wt = P(N - n);
  let e = 0;
  const contracting = new Int32Array(N);
  for (const k of keys) {
    const varying = orA[k] ^ andA[k];
    if (varying & S) { e++; continue; }                  // expanding
    for (let g = 0; g < N; g++) if (!((S >> g) & 1) && !((varying >> g) & 1)) contracting[g]++;
  }
  if (n >= 1) wE[S] = e * wt;
  for (let g = 0; g < N; g++) if (!((S >> g) & 1)) wA[S * N + g] = contracting[g] * wt;
}

// ------------------------------------------------- dynamic program
function solve(better, nodeW, edgeW) {
  const init = better(0, 1) ? Infinity : -Infinity;  // min -> +Inf, max -> -Inf
  const D = new Float64Array(SZ).fill(init), cnt = new Float64Array(SZ);
  D[0] = 0; cnt[0] = 1;
  const order = [...Array(SZ).keys()].sort((a, b) => pc(a) - pc(b));
  for (const S of order) {
    if (!S) continue;
    let bestV = init;
    for (let g = 0; g < N; g++) if ((S >> g) & 1) {
      const R = S & ~(1 << g), v = D[R] + edgeW(R, g);
      if (better(v, bestV)) bestV = v;
    }
    D[S] = bestV + nodeW(S);
    for (let g = 0; g < N; g++) if ((S >> g) & 1) {
      const R = S & ~(1 << g);
      if (D[R] + edgeW(R, g) === bestV) cnt[S] += cnt[R];
    }
  }
  const full = SZ - 1, out = [];
  (function back(S, suffix) {
    if (!S) { out.push(suffix); return; }
    for (let g = 0; g < N; g++) if ((S >> g) & 1) {
      const R = S & ~(1 << g);
      if (D[R] + edgeW(R, g) + nodeW(S) === D[S]) back(R, [g, ...suffix]);
    }
  })(full, []);
  return { value: D[full], count: cnt[full], orders: out };
}
const nodeE = S => (S === SZ - 1 ? 0 : wE[S]);
const zero = () => 0;
const edgeA = (S, g) => wA[S * N + g];

function scores(order) {
  let S = 0, muE = 0, muA = 0;
  for (const g of order) { muA += edgeA(S, g); S |= 1 << g; muE += nodeE(S); }
  return { muE, muA, muI: N * P(N) - muE - muA };
}


// Exact mean and standard deviation of mu_E over the uniform distribution on all N!
// orderings. A set S of size k is the prefix of a uniform ordering with probability
// k!(N-k)!/N!, and two sets are both prefixes only if one contains the other, in which
// case S < T (|S| = k, |T| = m) has probability k!(m-k)!(N-m)!/N!. The second moment is
// therefore a sum over nested pairs, 3^N terms.
function exactStats() {
  const lf = [0]; for (let i = 1; i <= N; i++) lf[i] = lf[i - 1] + Math.log(i);
  const prob = (...blocks) => Math.exp(blocks.reduce((s, b) => s + lf[b], 0) - lf[N]);
  let m1 = 0, m2 = 0;
  for (let T = 1; T < SZ - 1; T++) {
    const wT = wE[T]; if (!wT) continue;
    const m = pc(T);
    m1 += wT * prob(m, N - m);
    m2 += wT * wT * prob(m, N - m);
    for (let S = (T - 1) & T; S > 0; S = (S - 1) & T) {
      const wS = wE[S]; if (!wS) continue;
      const k = pc(S);
      m2 += 2 * wS * wT * prob(k, m - k, N - m);
    }
  }
  return { mean: m1, sd: Math.sqrt(Math.max(0, m2 - m1 * m1)) };
}

// ---------------------------------------------------------- output
console.log(`${csvPath}: N = ${N}, ${SZ} configurations, ${N}! orderings`);
if (evalOrder) {
  const s = scores(evalOrder);
  console.log(`ordering ${label(evalOrder)}`);
  console.log(`  (mu_E, mu_A, mu_I) = (${s.muE}, ${s.muA}, ${s.muI}),  sum = ${s.muE + s.muA + s.muI} = N p^N`);
  process.exit(0);
}
if (wantStats) {
  const st = exactStats(), E0 = solve((a, b) => a < b, nodeE, zero);
  console.log(`mu_E over all ${N}! orderings: mean = ${st.mean.toFixed(2)}, standard deviation = ${st.sd.toFixed(2)}`);
  console.log(`minimum ${E0.value}, z = ${((E0.value - st.mean) / st.sd).toFixed(2)}`);
  process.exit(0);
}
const E = solve((a, b) => a < b, nodeE, zero);
console.log(`min mu_E = ${E.value}, attained by ${E.count} ordering(s)`);
for (const o of E.orders) {
  const s = scores(o);
  console.log(`  [${o.join(',')}]  ${label(o)}   (mu_E, mu_A, mu_I) = (${s.muE}, ${s.muA}, ${s.muI})`);
}
const A = solve((a, b) => a > b, zero, edgeA);
console.log(`max mu_A = ${A.value}, attained by ${A.count} ordering(s)`);
if (listMuA) for (const o of A.orders) {
  const s = scores(o);
  console.log(`  [${o.join(',')}]  ${label(o)}   (mu_E, mu_A, mu_I) = (${s.muE}, ${s.muA}, ${s.muI})`);
}
const shared = E.orders.filter(o => A.orders.some(q => q.join() === o.join())).length;
console.log(`mu_E-minimisers that are also mu_A-maximisers: ${shared}`);
