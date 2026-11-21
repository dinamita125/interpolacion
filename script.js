const $ = id => document.getElementById(id);
let lastState = null;
const fmt = n => Number(n).toLocaleString('es-BO',{maximumFractionDigits:8});
const clean = n => Math.abs(n)<1e-12 ? 0 : n;
const names={linear:'Lineal por tramos',lagrange:'Lagrange',newton:'Newton',least:'Mínimos cuadrados'};
function ensureFinite(value){if(!Number.isFinite(value))throw Error('El cálculo produjo un valor no finito. Reduce la magnitud de los datos o evita valores X demasiado cercanos.');return value;}
function parsePoints(){
  const lines=$('points').value.trim().split(/\n+/).filter(Boolean);
  if(lines.length<2) throw Error('Ingresa por lo menos dos pares de datos X,Y.');
  const pts=lines.map((line,i)=>{
    const vals=line.trim().split(/[;,\s]+/).filter(Boolean);
    if(vals.length!==2) throw Error(`La línea ${i+1} debe contener exactamente dos valores: X, Y.`);
    const x=Number(vals[0]), y=Number(vals[1]);
    if(!Number.isFinite(x)||!Number.isFinite(y)) throw Error(`Los valores de la línea ${i+1} deben ser números válidos.`);
    return {x,y};
  });
  const seen=new Set();
  for(const p of pts){if(seen.has(p.x))throw Error(`El valor X = ${p.x} está repetido. Cada X debe ser distinto.`);seen.add(p.x);}
  return pts.sort((a,b)=>a.x-b.x);
}
function linear(pts,x){
  let a,b;
  if(x<=pts[0].x){a=pts[0];b=pts[1];}
  else if(x>=pts.at(-1).x){a=pts.at(-2);b=pts.at(-1);}
  else {let i=0;while(i<pts.length-1 && !(x>=pts[i].x&&x<=pts[i+1].x))i++;a=pts[i];b=pts[i+1];}
  const m=(b.y-a.y)/(b.x-a.x), y=a.y+m*(x-a.x);
  const fn=t=>{let p,q;if(t<=pts[0].x){p=pts[0];q=pts[1]}else if(t>=pts.at(-1).x){p=pts.at(-2);q=pts.at(-1)}else{let k=0;while(k<pts.length-1&&!(t>=pts[k].x&&t<=pts[k+1].x))k++;p=pts[k];q=pts[k+1]}return p.y+(q.y-p.y)/(q.x-p.x)*(t-p.x)};
  return {y,fn,detail:`Se selecciona el tramo formado por (${fmt(a.x)}, ${fmt(a.y)}) y (${fmt(b.x)}, ${fmt(b.y)}).\n\nPendiente:\nm = (y₁ − y₀) / (x₁ − x₀)\nm = (${fmt(b.y)} − ${fmt(a.y)}) / (${fmt(b.x)} − ${fmt(a.x)}) = ${fmt(m)}\n\nFórmula:\nf(x) = y₀ + m(x − x₀)\nf(${fmt(x)}) = ${fmt(a.y)} + ${fmt(m)}(${fmt(x)} − ${fmt(a.x)})\nf(${fmt(x)}) = ${fmt(y)}`};
}
function lagrange(pts,x){
  let y=0, terms=[];
  for(let i=0;i<pts.length;i++){
    let L=1, denom=1, factors=[];
    for(let j=0;j<pts.length;j++)if(i!==j){L*=(x-pts[j].x)/(pts[i].x-pts[j].x);denom*=pts[i].x-pts[j].x;factors.push(`(x − ${fmt(pts[j].x)})`);}
    const part=pts[i].y*L;y+=part;
    terms.push(`L${i}(x) = ${factors.join('') || '1'} / ${fmt(denom)};  y${i}·L${i}(${fmt(x)}) = ${fmt(part)}`);
  }
  const fn=t=>pts.reduce((sum,p,i)=>{let l=1;pts.forEach((q,j)=>{if(i!==j)l*=(t-q.x)/(p.x-q.x)});return sum+p.y*l},0);
  return {y,fn,detail:`Polinomio de Lagrange:\nP(x) = Σ yᵢ Lᵢ(x)\n\nContribuciones evaluadas en X = ${fmt(x)}:\n${terms.join('\n')}\n\nP(${fmt(x)}) = ${fmt(y)}`};
}
function newton(pts,x){
  const n=pts.length, d=pts.map(p=>p.y), coeff=[d[0]], cols=[d.slice()];
  for(let order=1;order<n;order++){
    const next=[];
    for(let i=0;i<n-order;i++)next.push((d[i+1]-d[i])/(pts[i+order].x-pts[i].x));
    coeff.push(next[0]);d.splice(0,d.length,...next);cols.push(next.slice());
  }
  const fn=t=>{let v=coeff[n-1];for(let i=n-2;i>=0;i--)v=coeff[i]+(t-pts[i].x)*v;return v};
  const y=fn(x);
  let polynomial=`P(x) = ${fmt(coeff[0])}`;
  for(let i=1;i<n;i++)polynomial+=` + (${fmt(coeff[i])})`+pts.slice(0,i).map(p=>`(x − ${fmt(p.x)})`).join('');
  const table=cols.map((c,i)=>`Orden ${i}: ${c.map(fmt).join(' | ')}`).join('\n');
  return {y,fn,detail:`Diferencias divididas de Newton:\n${table}\n\nCoeficientes: ${coeff.map(fmt).join(', ')}\n\nPolinomio de Newton:\n${polynomial}\n\nP(${fmt(x)}) = ${fmt(y)}`};
}
function leastSquares(pts,x){
  const n=pts.length,sx=pts.reduce((s,p)=>s+p.x,0),sy=pts.reduce((s,p)=>s+p.y,0),sxy=pts.reduce((s,p)=>s+p.x*p.y,0),sx2=pts.reduce((s,p)=>s+p.x*p.x,0);
  const den=n*sx2-sx*sx;if(Math.abs(den)<1e-14)throw Error('No se puede calcular el ajuste: los valores X no permiten obtener una pendiente.');
  const m=(n*sxy-sx*sy)/den,b=(sy-m*sx)/n,y=m*x+b;
  return {y,fn:t=>m*t+b,detail:`Ajuste lineal por mínimos cuadrados:\ny = mx + b\n\nm = [nΣxy − (Σx)(Σy)] / [nΣx² − (Σx)²]\nm = [${n}(${fmt(sxy)}) − (${fmt(sx)})(${fmt(sy)})] / [${n}(${fmt(sx2)}) − (${fmt(sx)}²)]\nm = ${fmt(m)}\n\nb = (Σy − mΣx) / n\nb = (${fmt(sy)} − ${fmt(m)}·${fmt(sx)}) / ${n} = ${fmt(b)}\n\nRecta aproximada: y = ${fmt(m)}x + ${fmt(b)}\n\nPara X = ${fmt(x)}:\ny = ${fmt(m)}(${fmt(x)}) + ${fmt(b)} = ${fmt(y)}`};
}
function draw(pts,fn,x,y){
  const c=$('plot'),ctx=c.getContext('2d'),W=c.width,H=c.height;
  ctx.clearRect(0,0,W,H);
  const xs=pts.map(p=>p.x).concat([x]),ys=pts.map(p=>p.y).concat([y]);
  let xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
  const xr=xmax-xmin||1,yr=ymax-ymin||1;xmin-=xr*.14;xmax+=xr*.14;ymin-=yr*.2;ymax+=yr*.2;
  const L=65,R=25,T=25,B=52,px=v=>L+(v-xmin)/(xmax-xmin)*(W-L-R),py=v=>H-B-(v-ymin)/(ymax-ymin)*(H-T-B);
  ctx.font='13px Arial';ctx.strokeStyle='#e8edf3';ctx.fillStyle='#667085';ctx.lineWidth=1;
  for(let i=0;i<=5;i++){let xx=xmin+(xmax-xmin)*i/5, yy=ymin+(ymax-ymin)*i/5;
    ctx.beginPath();ctx.moveTo(px(xx),T);ctx.lineTo(px(xx),H-B);ctx.stroke();ctx.fillText(fmt(xx),px(xx)-12,H-B+22);
    ctx.beginPath();ctx.moveTo(L,py(yy));ctx.lineTo(W-R,py(yy));ctx.stroke();ctx.fillText(fmt(yy),8,py(yy)+4);
  }
  ctx.strokeStyle='#344054';ctx.beginPath();ctx.moveTo(L,T);ctx.lineTo(L,H-B);ctx.lineTo(W-R,H-B);ctx.stroke();
  ctx.fillStyle='#344054';ctx.font='bold 14px Arial';ctx.fillText('X',W/2,H-10);ctx.save();ctx.translate(18,H/2);ctx.rotate(-Math.PI/2);ctx.fillText('Y',0,0);ctx.restore();
  ctx.strokeStyle='#1769aa';ctx.lineWidth=3;ctx.beginPath();
  for(let i=0;i<=220;i++){let t=xmin+(xmax-xmin)*i/220,v=fn(t);if(!Number.isFinite(v))continue;i?ctx.lineTo(px(t),py(v)):ctx.moveTo(px(t),py(v));}ctx.stroke();
  pts.forEach(p=>{ctx.beginPath();ctx.arc(px(p.x),py(p.y),6,0,Math.PI*2);ctx.fillStyle='#1769aa';ctx.fill();ctx.strokeStyle='white';ctx.lineWidth=2;ctx.stroke();});
  ctx.beginPath();ctx.arc(px(x),py(y),8,0,Math.PI*2);ctx.fillStyle='#e58b2a';ctx.fill();ctx.strokeStyle='white';ctx.lineWidth=2;ctx.stroke();
  ctx.fillStyle='#344054';ctx.font='12px Arial';ctx.fillText(`Estimación (${fmt(x)}, ${fmt(y)})`,Math.min(px(x)+12,W-190),Math.max(py(y)-12,18));
}
function compute(method,pts,x){if(method==='linear')return linear(pts,x);if(method==='lagrange')return lagrange(pts,x);if(method==='newton')return newton(pts,x);return leastSquares(pts,x)}
function calculate(){
  $('error').style.display='none';$('success').style.display='none';lastState=null;
  try{
    const pts=parsePoints(),x=Number($('target').value),method=$('method').value;
    if(!Number.isFinite(x))throw Error('Indica un valor X válido para realizar la estimación.');
    if(['lagrange','newton'].includes(method)&&pts.length>8)throw Error('Para Lagrange y Newton utiliza un máximo de ocho puntos.');
    const r=compute(method,pts,x),inside=x>=pts[0].x&&x<=pts.at(-1).x;ensureFinite(r.y);
    const types=method==='least'?'Ajuste':inside?'Interpolación':'Extrapolación';
    $('resultValue').textContent=`Y = ${fmt(clean(r.y))}`;$('resultSub').textContent=`${names[method]} · X = ${fmt(x)}`;
    $('countStat').textContent=pts.length;$('methodStat').textContent=names[method];$('typeStat').textContent=types;
    $('statusPill').textContent=types;$('statusPill').classList.toggle('warn',types==='Extrapolación');
    $('procedure').textContent=r.detail;
    $('dataTable').innerHTML=pts.map((p,i)=>`<tr><td>${i+1}</td><td>${fmt(p.x)}</td><td>${fmt(p.y)}</td></tr>`).join('');
    draw(pts,r.fn,x,r.y);
    const methods=['linear','lagrange','newton','least'];
    $('comparison').innerHTML=methods.map(m=>{try{if(['lagrange','newton'].includes(m)&&pts.length>8)throw 0;let z=compute(m,pts,x);return `<div class="compare-item"><small>${names[m]}</small><strong>${fmt(z.y)}</strong></div>`}catch(e){return `<div class="compare-item"><small>${names[m]}</small><strong>—</strong></div>`}}).join('');
    lastState={pts,x,method,r,names,types};$('success').style.display='block';
    const value=$('resultValue');value.classList.remove('pulse');void value.offsetWidth;value.classList.add('pulse');
  }catch(e){$('error').textContent=e.message;$('error').style.display='block';}
}
function load(points,x,method){$('points').value=points;$('target').value=x;$('method').value=method;calculate();}
$('calculate').addEventListener('click',calculate);
$('exampleLinear').addEventListener('click',()=>load('1, 5\n5, 13',3,'linear'));
$('examplePoly').addEventListener('click',()=>load('1, 2\n2, 5\n3, 10',2.5,'lagrange'));
$('exampleLeast').addEventListener('click',()=>load('1, 2\n2, 3\n3, 5',2.5,'least'));
$('clear').addEventListener('click',()=>{$('points').value='';$('target').value='';$('error').style.display='none';$('success').style.display='none';$('resultValue').textContent='Sin calcular';$('resultSub').textContent='Ingresa datos y presiona Calcular';$('procedure').textContent='El procedimiento aparecerá aquí.';$('dataTable').innerHTML='';$('comparison').innerHTML='';$('countStat').textContent='0';$('methodStat').textContent='—';$('typeStat').textContent='—';$('statusPill').textContent='Sin datos';$('plot').getContext('2d').clearRect(0,0,1000,500);lastState=null;});
$('exportCsv').addEventListener('click',()=>{
 if(!lastState){calculate();if(!lastState)return;}
 const s=lastState,rows=[['NUMÉRICALAB - RESULTADOS'],[],['Concepto','Valor'],['Método',s.names[s.method]],['X consultado',s.x],['Y estimado',s.r.y],['Tipo',s.types],[],['Punto','X','Y'],...s.pts.map((p,i)=>[i+1,p.x,p.y]),[],['Comparación de métodos'],...['linear','lagrange','newton','least'].map(m=>{try{return [s.names[m],compute(m,s.pts,s.x).y]}catch{return [s.names[m],'No disponible']}})];
 const csv=rows.map(r=>r.map(v=>`"${String(v??'').replaceAll('"','""')}"`).join(';')).join('\r\n');
 const blob=new Blob(['\ufeff'+csv],{type:'text/csv;charset=utf-8;'});download(blob,'resultados_interpolacion.csv');
});
$('exportPng').addEventListener('click',()=>{if(!lastState){calculate();if(!lastState)return;}const a=document.createElement('a');a.download='grafica_interpolacion.png';a.href=$('plot').toDataURL('image/png');a.click();});
function download(blob,name){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
['points','target','method'].forEach(id=>$(id).addEventListener('input',()=>{lastState=null;$('statusPill').textContent='Datos modificados';$('statusPill').classList.add('warn')}));
document.querySelectorAll('.tab').forEach(tab=>tab.addEventListener('click',()=>{document.querySelectorAll('.tab,.tab-panel').forEach(el=>el.classList.remove('active'));tab.classList.add('active');$(`panel-${tab.dataset.tab}`).classList.add('active')}));
const modal=$('helpModal');function setModal(open){modal.classList.toggle('open',open);modal.setAttribute('aria-hidden',String(!open));document.body.style.overflow=open?'hidden':''}
$('openHelp').addEventListener('click',()=>setModal(true));document.querySelectorAll('[data-close-modal]').forEach(el=>el.addEventListener('click',()=>setModal(false)));document.addEventListener('keydown',e=>{if(e.key==='Escape')setModal(false)});
$('printReport').addEventListener('click',()=>{if(!lastState){calculate();if(!lastState)return}window.print()});$('printManual').addEventListener('click',()=>window.print());
$('themeToggle').addEventListener('click',()=>{document.body.classList.toggle('dark');const dark=document.body.classList.contains('dark');$('themeToggle').textContent=dark?'☀':'☾';localStorage.setItem('numericalab-theme',dark?'dark':'light');if(lastState)draw(lastState.pts,lastState.r.fn,lastState.x,lastState.r.y)});if(localStorage.getItem('numericalab-theme')==='dark'){$('themeToggle').click()}
window.addEventListener('scroll',()=>{$('progress').style.width=`${scrollY/(document.documentElement.scrollHeight-innerHeight)*100||0}%`});
calculate();
