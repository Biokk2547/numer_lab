const test=require('node:test');
const assert=require('node:assert/strict');
const {ConjugateGradient}=require('../dist/conjugate-gradient.js');
const defaults={a11:4,a12:1,a21:1,a22:3,b1:1,b2:2,initial1:-2,initial2:-1,tolerance:1e-10,maxIterations:20};
const solve=overrides=>new ConjugateGradient({...defaults,...overrides}).solve();
const near=(a,b,tol=1e-9)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
test('CG solves the default SPD system and validates actual residual',()=>{
  const r=solve();assert.ok(r.converged);assert.equal(r.iterations,2);
  near(r.value[0],1/11);near(r.value[1],7/11);near(r.error,Math.hypot(...r.residual(r.value)),1e-15);
  for(const s of r.steps)near(s.residual,Math.hypot(...r.residual([s.x,s.y])),1e-15);
});
test('CG directions are A-conjugate and quadratic energy decreases',()=>{
  const r=solve(),[a,b]=r.steps;near(r.dot(a.search,r.multiply(b.search)),0,1e-9);
  assert.ok(a.objective<r.initial.objective);assert.ok(b.objective<a.objective);
});
test('CG handles different SPD matrices, zero RHS and exact initial solution',()=>{
  const tilted=solve({a11:10,a12:3,a21:3,a22:2});near(tilted.value[0],-4/11);near(tilted.value[1],17/11);
  const zero=solve({b1:0,b2:0});assert.ok(zero.converged);near(zero.value[0],0);near(zero.value[1],0);
  const initial=solve({a11:2,a12:0,a21:0,a22:2,b1:2,b2:4,initial1:1,initial2:2});assert.equal(initial.iterations,0);assert.equal(initial.steps[0].iteration,0);assert.equal(initial.error,0);
  assert.equal(solve({a11:2,a12:0,a21:0,a22:2}).iterations,1);
});
test('CG rejects non-symmetric, indefinite, singular and malformed systems',()=>{
  for(const invalid of [{a21:2},{a11:-4},{a11:1,a12:2,a21:2,a22:1},{a11:1,a12:1,a21:1,a22:1},{b1:NaN},{initial1:Infinity},{maxIterations:0},{maxIterations:1.5},{tolerance:0}])assert.throws(()=>solve(invalid));
});
test('CG reports the iteration limit without claiming convergence',()=>{
  const limited=solve({maxIterations:1});assert.equal(limited.converged,false);assert.equal(limited.iterations,1);assert.ok(limited.error>defaults.tolerance);
});
