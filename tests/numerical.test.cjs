const test=require('node:test');
const assert=require('node:assert/strict');
const {Expression,Bisection,NewtonRaphson,Secant,Quadrature,EulerMethod}=require('../dist/app.js');
const near=(actual,expected,tolerance)=>assert.ok(Math.abs(actual-expected)<tolerance,`${actual} != ${expected}`);
test('Expression precedence, constants, functions and malformed input',()=>{
  near(new Expression('-x^2 + 2^3^2').evaluate(3),503,1e-10);
  near(new Expression('sin(pi/2)+exp(0)+sqrt(4)+log(e)').evaluate(0),5,1e-10);
  near(new Expression('y-x^2+1').evaluate(2,5),2,1e-10);
  for(const expr of ['2x','sin(x','x;alert(1)','foo(x)',''])assert.throws(()=>new Expression(expr));
});
test('Root solvers converge to the cubic root',()=>{
  for(const [Model,params] of [[Bisection,{a:1,b:2}],[NewtonRaphson,{x0:1.5}],[Secant,{x0:1,x1:2}]]){
    const r=new Model(new Expression('x^3-x-2'),{...params,tolerance:1e-8,maxIterations:100}).solve();
    assert.ok(r.converged);near(r.value,1.5213797068,1e-7);assert.ok(r.steps.length>0);
  }
});
test('Invalid bracket, zero slope, equal secant values and iteration limit',()=>{
  assert.throws(()=>new Bisection(new Expression('x^2+1'),{a:-1,b:1,tolerance:1e-6,maxIterations:50}).solve());
  assert.throws(()=>new NewtonRaphson(new Expression('x^2+1'),{x0:0,tolerance:1e-6,maxIterations:50}).solve());
  assert.throws(()=>new Secant(new Expression('x^2+1'),{x0:-1,x1:1,tolerance:1e-6,maxIterations:50}).solve());
  assert.equal(new Bisection(new Expression('x^3-x-2'),{a:1,b:2,tolerance:1e-12,maxIterations:1}).solve().converged,false);
  assert.equal(new Bisection(new Expression('x-1'),{a:1,b:2,tolerance:1e-6,maxIterations:50}).solve().value,1);
});
test('Quadrature area, Simpson exactness and odd interval rejection',()=>{
  near(new Quadrature(new Expression('sin(x)'),{a:0,b:Math.PI,n:100}).solve().value,2,.0002);
  near(new Quadrature(new Expression('x^3'),{a:0,b:2,n:10},true).solve().value,4,1e-10);
  assert.throws(()=>new Quadrature(new Expression('sin(x)'),{a:0,b:Math.PI,n:11},true).solve());
});
test('Euler approaches known ODE solution as interval count increases',()=>{
  const solve=n=>new EulerMethod(new Expression('y'),{x0:0,y0:1,end:1,n}).solve();
  const coarse=solve(20),fine=solve(200);assert.ok(Math.abs(fine.value-Math.E)<Math.abs(coarse.value-Math.E));near(fine.value,Math.E,.007);assert.equal(fine.steps.at(-1).x,1);
});
