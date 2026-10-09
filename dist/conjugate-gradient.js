/* Two-dimensional conjugate gradient model, separate from the visual interface. */
const CGBase = typeof module !== 'undefined' ? require('./app.js').NumericalMethod : NumericalMethod;
class ConjugateGradient extends CGBase {
  constructor(params) { super(null, params); }
  multiply(v) { const p=this.params; return [p.a11*v[0]+p.a12*v[1],p.a21*v[0]+p.a22*v[1]]; }
  dot(a,b) { return a[0]*b[0]+a[1]*b[1]; }
  residual(v) { const av=this.multiply(v); return [this.params.b1-av[0],this.params.b2-av[1]]; }
  objective(v) { return .5*this.dot(v,this.multiply(v))-this.dot([this.params.b1,this.params.b2],v); }
  solve() {
    const p=this.params, names=['a11','a12','a21','a22','b1','b2','initial1','initial2','tolerance','maxIterations'];
    if(names.some(name=>!Number.isFinite(p[name])||Math.abs(p[name])>1e12)) throw new Error('กรุณากรอกสมาชิกเมทริกซ์ เวกเตอร์ และค่าตั้งต้นเป็นจำนวนที่ถูกต้อง');
    if(p.tolerance<1e-12||p.tolerance>1) throw new Error('Tolerance ต้องอยู่ระหว่าง 10⁻¹² และ 1');
    if(!Number.isInteger(p.maxIterations)||p.maxIterations<1||p.maxIterations>500) throw new Error('จำนวนรอบสูงสุดต้องเป็นจำนวนเต็ม 1–500');
    if(p.a12!==p.a21) throw new Error('A ต้องสมมาตร: a₁₂ ต้องเท่ากับ a₂₁');
    const schur=p.a22-p.a12*(p.a12/p.a11);
    if(p.a11<=0||!Number.isFinite(schur)||schur<=0) throw new Error('A ต้องเป็น positive definite: a₁₁ > 0 และ a₂₂ − a₁₂²/a₁₁ > 0');
    const target2=(p.b2-p.a12*(p.b1/p.a11))/schur;
    this.target=[(p.b1-p.a12*target2)/p.a11,target2];
    if(!this.target.every(Number.isFinite)) throw new Error('เมทริกซ์ใกล้ singular เกินไปสำหรับการแสดงกราฟ');
    let x=[p.initial1,p.initial2],r=this.residual(x),direction=[...r],rr=this.dot(r,r);
    this.initial={vector:[...x],residual:Math.hypot(...r),objective:this.objective(x)};
    if(this.initial.residual<=p.tolerance) {
      this.steps.push({iteration:0,x:x[0],y:x[1],alpha:0,beta:0,residual:this.initial.residual,objective:this.initial.objective});
      this.iterations=0;return this.finish([...x],this.initial.residual,true,'ค่าตั้งต้นผ่านเกณฑ์ residual อยู่แล้ว');
    }
    for(let i=0;i<p.maxIterations;i++) {
      const ap=this.multiply(direction),denominator=this.dot(direction,ap);
      if(!Number.isFinite(denominator)||denominator<=0||rr===0) break;
      const alpha=rr/denominator,from=[...x],search=[...direction];
      x=x.map((v,j)=>v+alpha*direction[j]);
      const nextR=r.map((v,j)=>v-alpha*ap[j]),nextRR=this.dot(nextR,nextR),beta=nextRR/rr;
      const actualResidual=Math.hypot(...this.residual(x));
      if(!x.every(Number.isFinite)||!Number.isFinite(actualResidual)||!Number.isFinite(beta)) throw new Error('การคำนวณเกินขอบเขตตัวเลข ลองเปลี่ยนเมทริกซ์หรือค่าตั้งต้น');
      this.steps.push({iteration:i+1,from,search,x:x[0],y:x[1],alpha,beta,residual:actualResidual,objective:this.objective(x)});this.iterations=i+1;
      if(actualResidual<=p.tolerance) return this.finish([...x],actualResidual,true,'‖b − Ax‖₂ ≤ tolerance ที่กำหนด');
      direction=nextR.map((v,j)=>v+beta*direction[j]);r=nextR;rr=nextRR;
    }
    if(!this.steps.length) throw new Error('ไม่สามารถสร้างทิศทางค้นหาที่เสถียรได้ ลองเปลี่ยนเมทริกซ์หรือค่าตั้งต้น');
    return this.finish([...x],Math.hypot(...this.residual(x)),false,this.iterations===p.maxIterations?'ถึงจำนวนรอบสูงสุด':'ถึงขีดจำกัดความแม่นยำของตัวเลข');
  }
}
if(typeof module !== 'undefined') module.exports={ConjugateGradient};

if(typeof document !== 'undefined') {
  const CG_PRESETS=[
    {name:'ระบบพื้นฐาน · A = [4, 1; 1, 3]',values:{a11:4,a12:1,a21:1,a22:3,b1:1,b2:2}},
    {name:'วงรีเอียง · A = [10, 3; 3, 2]',values:{a11:10,a12:3,a21:3,a22:2,b1:1,b2:2}},
    {name:'ทิศทางสมมาตร · A = [6, 4; 4, 6]',values:{a11:6,a12:4,a21:4,a22:6,b1:2,b2:-2}}
  ];
  METHODS.cg={name:'Conjugate Gradient',mode:'linear',model:ConjugateGradient,
    title:'ค้นหาคำตอบ ด้วยทิศทาง conjugate',
    text:'แก้ Ax = b โดยเลือกทิศทางค้นหาที่เป็น A-conjugate กับทิศทางก่อนหน้า เส้นทางบนกราฟเดินสู่จุดต่ำสุดของ φ(x) ระบบนี้ใช้เมทริกซ์สมมาตร positive definite ขนาด 2 × 2',
    formula:'r₀ = b − Ax₀ · p₀ = r₀\nαₖ = (rₖᵀrₖ)/(pₖᵀApₖ)\nxₖ₊₁ = xₖ + αₖpₖ\nrₖ₊₁ = rₖ − αₖApₖ\nβₖ = (rₖ₊₁ᵀrₖ₊₁)/(rₖᵀrₖ)\npₖ₊₁ = rₖ₊₁ + βₖpₖ',
    fields:[['a11','a₁₁',4],['a12','a₁₂',1],['a21','a₂₁',1],['a22','a₂₂',3],['b1','b₁',1],['b2','b₂',2],['initial1','x₁ เริ่มต้น',-2],['initial2','x₂ เริ่มต้น',-1],['tolerance','Tolerance ε',.000001],['maxIterations','จำนวนรอบสูงสุด',20]],
    columns:[['x','x₁'],['y','x₂'],['alpha','α'],['beta','β'],['residual','‖r‖₂'],['objective','φ(x)']]};

  class ConjugateGradientPlot {
    constructor(renderer) { this.renderer=renderer; }
    fit(result) {
      const points=[result.initial.vector,result.target,...result.steps.map(s=>[s.x,s.y])];
      const minX=Math.min(...points.map(v=>v[0])),maxX=Math.max(...points.map(v=>v[0])),minY=Math.min(...points.map(v=>v[1])),maxY=Math.max(...points.map(v=>v[1]));
      const span=Math.max(maxX-minX,maxY-minY,.5),padding=.3*span,cx=(minX+maxX)/2,cy=(minY+maxY)/2;
      this.renderer.bounds={minX:cx-span/2-padding,maxX:cx+span/2+padding,minY:cy-span/2-padding,maxY:cy+span/2+padding};
    }
    draw(index) {
      const renderer=this.renderer,result=renderer.result,canvas=renderer.canvas,ctx=renderer.ctx;
      const width=canvas.clientWidth,height=canvas.clientHeight,dpr=window.devicePixelRatio||1;if(!width||!height)return;
      canvas.width=width*dpr;canvas.height=height*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,width,height);
      const bounds=renderer.bounds,left=48,right=25,top=24,bottom=35;
      const X=x=>left+(x-bounds.minX)/(bounds.maxX-bounds.minX)*(width-left-right),Y=y=>height-bottom-(y-bounds.minY)/(bounds.maxY-bounds.minY)*(height-top-bottom);
      renderer.map=null;ctx.font='9px Manrope, sans-serif';ctx.lineWidth=1;
      for(const axis of ['x','y']) {
        const low=axis==='x'?bounds.minX:bounds.minY,high=axis==='x'?bounds.maxX:bounds.maxY,step=renderer.niceStep(high-low);
        for(let v=Math.ceil(low/step)*step;v<=high;v+=step) {
          ctx.strokeStyle=Math.abs(v)<step*1e-8?'#d4cadf':'#f0edf5';ctx.beginPath();if(axis==='x'){ctx.moveTo(X(v),top);ctx.lineTo(X(v),height-bottom);}else{ctx.moveTo(left,Y(v));ctx.lineTo(width-right,Y(v));}ctx.stroke();
          ctx.fillStyle='#a99db7';ctx.textAlign=axis==='x'?'center':'right';ctx.fillText(Number(v.toPrecision(3)).toString(),axis==='x'?X(v):left-10,axis==='x'?height-bottom+18:Y(v)+3);
        }
      }
      ctx.fillStyle='#9b89b2';ctx.textAlign='right';ctx.fillText('x₁',width-right+13,height-bottom+19);ctx.fillText('x₂',left-11,top-8);
      ctx.save();ctx.beginPath();ctx.rect(left,top,width-left-right,height-top-bottom);ctx.clip();
      const p=result.params,theta=.5*Math.atan2(2*p.a12,p.a11-p.a22),cos=Math.cos(theta),sin=Math.sin(theta);
      const lambda1=(p.a11+p.a22+Math.hypot(p.a11-p.a22,2*p.a12))/2,lambda2=p.a11*(p.a22-p.a12*(p.a12/p.a11))/lambda1;
      const d=result.initial.vector.map((v,j)=>v-result.target[j]),initialEnergy=.5*result.dot(d,result.multiply(d)),span=bounds.maxX-bounds.minX,base=Math.max(initialEnergy,lambda2*span*span*.08,1e-12);
      for(const fraction of [.025,.07,.15,.3,.55,.9,1.35,2]) {
        const energy=base*fraction,r1=Math.sqrt(2*energy/lambda1),r2=Math.sqrt(2*energy/lambda2);ctx.beginPath();
        for(let j=0;j<=180;j++){const t=j/180*2*Math.PI,u=r1*Math.cos(t),v=r2*Math.sin(t),x=result.target[0]+cos*u-sin*v,y=result.target[1]+sin*u+cos*v;if(j===0)ctx.moveTo(X(x),Y(y));else ctx.lineTo(X(x),Y(y));}
        ctx.closePath();ctx.strokeStyle='#cbb9e6';ctx.lineWidth=1;ctx.stroke();
      }
      const path=[result.initial.vector,...result.steps.map(s=>[s.x,s.y])];
      const line=(points,color,dash=[])=>{ctx.beginPath();points.forEach(([x,y],i)=>i?ctx.lineTo(X(x),Y(y)):ctx.moveTo(X(x),Y(y)));ctx.strokeStyle=color;ctx.lineWidth=2;ctx.setLineDash(dash);ctx.stroke();ctx.setLineDash([]);};
      line(path,'#d9d0e9',[4,4]);line(path.slice(0,index+2),'#7858d6');
      const dot=(v,color,radius)=>{ctx.beginPath();ctx.arc(X(v[0]),Y(v[1]),radius+2,0,2*Math.PI);ctx.fillStyle='white';ctx.fill();ctx.beginPath();ctx.arc(X(v[0]),Y(v[1]),radius,0,2*Math.PI);ctx.fillStyle=color;ctx.fill();};
      dot(result.initial.vector,'#db9b54',4);result.steps.slice(0,index+1).forEach(s=>dot([s.x,s.y],'#a98ad5',3));
      const current=result.steps[index];dot([current.x,current.y],'#399d7a',5);
      const tx=X(result.target[0]),ty=Y(result.target[1]);ctx.strokeStyle='#399d7a';ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(tx-7,ty);ctx.lineTo(tx+7,ty);ctx.moveTo(tx,ty-7);ctx.lineTo(tx,ty+7);ctx.stroke();
      ctx.fillStyle='#9783b1';ctx.font='9px Manrope';ctx.textAlign='left';ctx.fillText('เริ่มต้น',X(result.initial.vector[0])+8,Y(result.initial.vector[1])-9);ctx.fillStyle='#399d7a';ctx.fillText('x*',tx+10,ty-10);ctx.restore();
      canvas.setAttribute('aria-label',`Conjugate Gradient รอบ ${current.iteration}: x₁=${current.x}, x₂=${current.y}, residual=${current.residual}`);this.drawResidual(result,index);
    }
    drawResidual(result,index) {
      const c=document.getElementById('cgResidual'),ctx=c.getContext('2d'),w=c.clientWidth,h=c.clientHeight,dpr=window.devicePixelRatio||1;if(!w||!h)return;
      c.width=w*dpr;c.height=h*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
      const series=[{iteration:0,residual:result.initial.residual},...result.steps.filter(s=>s.iteration>0)],floor=Math.min(result.params.tolerance*.01,1e-16),logs=series.map(s=>Math.log10(Math.max(s.residual,floor))),tol=Math.log10(result.params.tolerance);
      const low=Math.min(...logs,tol)-.6,high=Math.max(...logs,tol)+.6,maxIteration=Math.max(1,result.iterations),X=i=>46+i/maxIteration*(w-67),Y=v=>h-23-(v-low)/(high-low)*(h-38);
      ctx.font='8px Manrope';ctx.fillStyle='#a496b6';ctx.textAlign='right';ctx.strokeStyle='#eee8f6';ctx.lineWidth=1;
      const tickStep=Math.max(1,Math.ceil((high-low)/4));for(let v=Math.ceil(low/tickStep)*tickStep;v<=high;v+=tickStep){ctx.beginPath();ctx.moveTo(46,Y(v));ctx.lineTo(w-21,Y(v));ctx.stroke();ctx.fillText('1e'+v,39,Y(v)+3);}
      ctx.strokeStyle='#77b19b';ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(46,Y(tol));ctx.lineTo(w-21,Y(tol));ctx.stroke();ctx.setLineDash([]);ctx.fillStyle='#65a38a';ctx.fillText('ε',w-7,Y(tol)+3);
      ctx.beginPath();series.forEach((s,i)=>i?ctx.lineTo(X(s.iteration),Y(logs[i])):ctx.moveTo(X(s.iteration),Y(logs[i])));ctx.strokeStyle='#b5a1d4';ctx.lineWidth=1.7;ctx.stroke();
      series.forEach((s,i)=>{ctx.beginPath();ctx.arc(X(s.iteration),Y(logs[i]),s.iteration===result.steps[index].iteration?4:2.5,0,2*Math.PI);ctx.fillStyle=s.iteration===result.steps[index].iteration?'#399d7a':'#7858d6';ctx.fill();if(series.length<12||i===0||i===series.length-1){ctx.fillStyle='#a496b6';ctx.textAlign='center';ctx.fillText(s.iteration,X(s.iteration),h-7);}});
      c.setAttribute('aria-label',series.map(s=>`รอบ ${s.iteration}: residual ${s.residual}`).join('; '));
    }
  }

  class ExtendedPlotRenderer extends PlotRenderer {
    fit(){if(this.key!=='cg')return super.fit();this.cgPlot??=new ConjugateGradientPlot(this);this.cgPlot.fit(this.result);this.draw();}
    draw(index=this.stepIndex||0){if(this.key!=='cg')return super.draw(index);this.stepIndex=Math.max(0,Math.min(index,this.result.steps.length-1));this.cgPlot?.draw(this.stepIndex);}
  }

  class ExtendedNumericalLabApp extends NumericalLabApp {
    constructor(){super(ExtendedPlotRenderer);}
    selectMethod(key) {
      this.$('functionFields').hidden=key==='cg';this.$('cgResidualPanel').hidden=key!=='cg';this.$('plotTooltip').hidden=true;this.$('conceptFormula').classList.toggle('cg-formulas',key==='cg');
      document.querySelector('.workspace').classList.toggle('cg-workspace',key==='cg');
      document.querySelector('.heading-art').innerHTML=key==='cg'?'φ<span>(x)</span><i>→ min</i>':'ƒ<span>(x)</span><i>≈ 0</i>';
      if(key!=='cg')return super.selectMethod(key);this.pause();this.key=key;const m=METHODS.cg;
      document.querySelectorAll('[data-mode]').forEach(b=>b.classList.toggle('active',b.dataset.mode==='linear'));
      this.$('crumb').textContent='ระบบสมการเชิงเส้น';this.$('pageTitle').innerHTML='เดินสู่<span>คำตอบของระบบ</span>';this.$('pageSubtitle').textContent='สำรวจ Conjugate Gradient สำหรับ Ax = b ขนาด 2 × 2 ผ่านเส้นชั้นความสูงและ residual';
      this.$('methodTabs').innerHTML='<button class="method-tab active" role="tab" aria-selected="true"><span>01</span>Conjugate Gradient</button>';
      this.$('preset').innerHTML='';CG_PRESETS.forEach((p,i)=>this.$('preset').add(new Option(p.name,String(i))));this.$('preset').add(new Option('กำหนดระบบเอง','custom'));this.$('parameterFields').innerHTML='';
      const groups=[['เมทริกซ์ A · 2 × 2',m.fields.slice(0,4)],['เวกเตอร์ b',m.fields.slice(4,6)],['เวกเตอร์เริ่มต้น x₀',m.fields.slice(6,8)],['เกณฑ์หยุด',m.fields.slice(8)]];
      for(const [title,fields] of groups){const label=document.createElement('div');label.className='cg-group-title';label.textContent=title;this.$('parameterFields').append(label);for(const [name,text,value] of fields){const wrap=document.createElement('div');wrap.className='parameter';const l=document.createElement('label');l.htmlFor='param-'+name;l.textContent=text;const input=document.createElement('input');input.type='number';input.id='param-'+name;input.name=name;input.step=name==='maxIterations'?'1':'any';input.value=value;input.required=true;input.addEventListener('input',()=>{this.$('preset').value='custom';});wrap.append(l,input);this.$('parameterFields').append(wrap);}}
      this.$('conceptTitle').textContent=m.title;this.$('conceptText').textContent=m.text;this.$('conceptFormula').textContent=m.formula;this.$('plotBadge').textContent='CONJUGATE GRADIENT';
      this.$('graphLegend').innerHTML='<span><i class="legend-line dashed"></i> เส้นชั้นความสูง φ(x)</span><span><i class="legend-line"></i> เส้นทางค้นหา</span><span><i class="legend-line green"></i> ค่าปัจจุบัน / x*</span>';this.calculate();
    }
    applyPresetBounds(){if(this.key!=='cg')return super.applyPresetBounds();const preset=CG_PRESETS[Number(this.$('preset').value)];if(preset)for(const [name,value] of Object.entries(preset.values))this.$('param-'+name).value=value;}
    calculate() {
      if(this.key!=='cg')return super.calculate();this.pause();
      try{const result=new ConjugateGradient(this.readParams()).solve();this.result=result;this.$('formError').hidden=true;this.$('plotExpression').textContent='φ(x) = ½xᵀAx − bᵀx';this.$('plotCaption').textContent='เส้นชั้นความสูง · จุดต่ำสุดคือคำตอบของ Ax = b';this.plot.setResult(result,this.key);this.$('stepSlider').max=result.steps.length;this.renderResults();this.renderTable();this.setStep(0);return true;}
      catch(error){this.$('formError').textContent=error.message;this.$('formError').hidden=false;return false;}
    }
    renderResults() {
      if(this.key!=='cg')return super.renderResults();const r=this.result;
      this.$('answerLabel').innerHTML='เวกเตอร์คำตอบ x <span>↗</span>';this.$('answer').textContent=`[${r.value[0].toFixed(5)}, ${r.value[1].toFixed(5)}]`;this.$('answer').classList.add('vector-answer');this.$('answerDetail').textContent='[x₁, x₂] · คำตอบสุดท้าย';this.$('iterations').textContent=String(r.iterations).padStart(2,'0');
      this.$('errorLabel').innerHTML='Residual ‖b − Ax‖₂ <span>◎</span>';this.$('errorValue').textContent=this.format(r.error);this.$('errorDetail').textContent='ใช้ residual จริงตรวจเกณฑ์หยุด';this.$('resultStatus').textContent=r.converged?'ลู่เข้าแล้ว':'ยังไม่ลู่เข้า';this.$('statusDetail').textContent=r.reason;this.$('resultStatus').closest('article').classList.toggle('warn',!r.converged);
    }
    renderTable(){super.renderTable();if(this.key==='cg'){this.$('tableCount').textContent=this.result.iterations+' รอบ';this.$('tableBody').querySelectorAll('.iteration-pill').forEach((cell,i)=>{cell.textContent=this.result.steps[i].iteration===0?'เริ่มต้น':String(this.result.steps[i].iteration).padStart(2,'0');});}else this.$('answer').classList.remove('vector-answer');}
    setStep(index){super.setStep(index);if(this.key==='cg'&&this.result){const s=this.result.steps[this.step];this.$('stepExplanation').textContent=`${s.iteration===0?'ค่าตั้งต้น':`รอบ ${s.iteration}`} → x = [${this.format(s.x)}, ${this.format(s.y)}] · α = ${this.format(s.alpha)} · ‖r‖₂ = ${this.format(s.residual)}`;}}
    tooltip(e){if(this.key==='cg'){this.$('plotTooltip').hidden=true;return;}return super.tooltip(e);}
    registerTools() {
      const cg=METHODS.cg;super.registerTools();
      if(!document.modelContext?.registerTool)return;const fields=Object.fromEntries(cg.fields.map(([name])=>[name,{type:'number'}]));
      try{Promise.resolve(document.modelContext.registerTool({name:'run_conjugate_gradient',title:'Run Conjugate Gradient',description:'Solve a symmetric positive definite 2 × 2 system and update its contour plot, residual chart and iteration table.',inputSchema:{type:'object',properties:fields,required:Object.keys(fields),additionalProperties:false},annotations:{readOnlyHint:false,untrustedContentHint:false},execute:input=>{if(!input||typeof input!=='object'||Object.keys(input).some(k=>!Object.hasOwn(fields,k)))throw new Error('Invalid CG parameters');const result=new ConjugateGradient(input).solve();this.selectMethod('cg');for(const [name,value] of Object.entries(input))this.$('param-'+name).value=value;this.$('preset').value='custom';if(!this.calculate())throw new Error(this.$('formError').textContent);return {value:result.value,residual:result.error,iterations:result.iterations,converged:result.converged};}},{signal:this.toolLifecycle?.signal})).catch(()=>{});}catch{}
    }
  }
  window.numericalLab=new ExtendedNumericalLabApp();
}
