// Extrae las funciones de cálculo de una demo HTML y las evalúa fuera del navegador.
//
// Uso: node demo_harness.js <demo.html> <casos.json>
// casos.json: [{N, p, tt, orderings: [[...], ...], bnb: bool}, ...]
// Imprime JSON: {networks: {...}, cases: [...]}.
'use strict';
const fs = require('fs');
const vm = require('vm');

const [htmlPath, casesPath] = process.argv.slice(2);
const html = fs.readFileSync(htmlPath, 'utf8');
const cases = JSON.parse(fs.readFileSync(casesPath, 'utf8'));

// Cuerpo completo de `function name(...) {...}` contando llaves fuera de cadenas.
function extractFunction(name) {
  const m = html.match(new RegExp(`(?:async\\s+)?function\\s+${name}\\s*\\(`));
  if (!m) return null;
  let i = html.indexOf('{', m.index + m[0].length);
  let depth = 0, quote = null;
  for (let j = i; j < html.length; j++) {
    const c = html[j];
    if (quote) {
      if (c === '\\') { j++; continue; }
      if (c === quote) quote = null;
      continue;
    }
    if (c === '"' || c === "'") { quote = c; continue; }
    if (c === '{') depth++;
    else if (c === '}' && --depth === 0) return html.slice(m.index, j + 1);
  }
  throw new Error(`llaves sin cerrar en ${name}`);
}

const names = ['configToVector', 'vectorToConfig', 'computeBallDiameter', 'computeScores',
  'classifyBall', 'computeAllLevels', 'computeMuFastBit', 'computeMuFast', 'bnbLevelContrib', 'bnbSearch'];
const found = names.filter(n => extractFunction(n));
const consts = [...html.matchAll(/const (TT_[A-Z]+)\s*=\s*(\[[^\]]*\])/g)].map(m => `const ${m[1]} = ${m[2]};`);
const ctx = vm.createContext({ performance: { now: () => Date.now() }, setTimeout, Math, Array, Int32Array, Uint8Array, Promise });
vm.runInContext(['var _searchStop = false;', ...consts, ...found.map(extractFunction)].join('\n'), ctx);

// Valores esperados declarados en NETWORKS (si la demo los tiene).
const networks = {};
for (const m of html.matchAll(/(\w+):\s*\{\s*name:\s*["'][^"']*["'],\s*N:\s*(\d+),\s*p:\s*(\d+),[\s\S]*?optimal:\s*\[([\d,]+)\],[\s\S]*?tt:\s*(TT_[A-Z]+)/g)) {
  const block = m[0];
  const exp = block.match(/expectedE:(\d+), expectedA:(\d+), expectedI:(\d+)/);
  const tt = vm.runInContext(m[5], ctx);
  const N = +m[2], p = +m[3], A = m[4].split(',').map(Number);
  networks[m[1]] = {
    N, p, ttLength: tt.length, optimal: A,
    expected: exp ? exp.slice(1).map(Number) : null,
    fast: vm.runInContext('computeMuFast', ctx)(tt, N, p, A),
  };
}

(async () => {
  const out = [];
  for (const c of cases) {
    ctx.__tt = c.tt;
    const res = { scores: [], fast: [], levels: [] };
    for (const A of c.orderings) {
      ctx.__A = A;
      if (found.includes('computeScores')) {
        const s = vm.runInContext(`computeScores(${c.N}, ${c.p}, __A, __tt)`, ctx);
        res.scores.push({ mu: [s.muE, s.muA, s.muI], E: s.E, A: s.A, I: s.I });
      }
      if (found.includes('computeAllLevels')) {
        res.levels.push(vm.runInContext(`computeAllLevels(${c.N}, ${c.p}, __A, __tt)`, ctx));
      }
      res.fast.push(vm.runInContext(`computeMuFast(__tt, ${c.N}, ${c.p}, __A)`, ctx));
    }
    if (c.bnb) {
      res.bnb = await new Promise(resolve => {
        vm.runInContext('bnbSearch', ctx)(c.tt, c.N, c.p, () => {}, (optima, best, nodes, prunes) =>
          resolve({ best, prunes, optima: optima.map(o => ({ A: o.A, mu: [o.muE, o.muA, o.muI] })) }));
      });
    }
    out.push(res);
  }
  process.stdout.write(JSON.stringify({ functions: found, networks, cases: out }));
})();
