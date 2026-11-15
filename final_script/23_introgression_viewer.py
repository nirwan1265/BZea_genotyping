#!/usr/bin/env python3
"""
23 — Interactive introgression viewer: one HTML file, drop-down menus by taxon, family and pool.

Reads the step-21 calls and writes a single self-contained HTML page (no internet needed):
  - pick a taxon, then a family, then "All pools" (family view) or one pool
  - family view: each chromosome shows one thin track per pool of the family, so shared
    segments line up
  - pool view: the 10 chromosomes of one pool (like the blocks_<GROUP>.pdf pages)
  - hover a segment for its position and state; a table lists the pools of the family
  - "Download SVG" / "Download PNG" save the figure that is on screen

Colors (colorblind-safe, checked for protan/deutan/tritan): B73 blue #3A6FB5, SEG yellow #DDAA33,
DONOR red #BB5566.
States are POOL states: B73 = all plants AA; SEG = the family still segregates, so plants in the
pool are AA, AB or BB (1:2:1); DONOR = all plants BB.

Inputs:  --tracts  tracts_all.tsv.gz (step 21)   --pool_summary pool_summary.tsv (step 21)
         --meta    SAMPLE_META (col 1 sample, col 2 pedigree; family = first two '_' fields)
Output:  --out     introgression_viewer.html
"""
import argparse
import collections
import csv
import gzip
import json

UNIT = 100000      # positions are stored in 100-kb units (the k-mer window)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracts", required=True)
    ap.add_argument("--pool_summary", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="BZea BC2S3 introgressions (k-mer HMM)")
    a = ap.parse_args()

    fam = {}
    with open(a.meta, encoding="utf-8", errors="ignore") as f:
        next(f)
        for line in f:
            p = line.rstrip("\r\n").split("\t")
            if len(p) > 1: fam[p[0]] = "_".join(p[1].split("_")[:2])

    pools = {}
    with open(a.pool_summary) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pools[r["SAMPLE"]] = dict(
                taxon=r["GROUP"], family=fam.get(r["SAMPLE"], "unknown"), flag=r["FLAG"], sigma=float(r["SIGMA"]),
                seg=float(r["FRAC_SEG"]), donor=float(r["FRAC_DONOR"]), share=float(r["DONOR_SHARE"]),
                ntr=int(r["N_TRACTS"]), t=[])

    clen, chroms = {}, []
    code = {"SEG": 1, "DONOR": 2}
    with gzip.open(a.tracts, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            c, s, e = r["CHROM"], int(r["START"]), int(r["END"])
            if c not in clen: chroms.append(c); clen[c] = 0
            clen[c] = max(clen[c], e)
            if r["STATE"] in code and r["SAMPLE"] in pools:
                pools[r["SAMPLE"]]["t"].append([chroms.index(c), s // UNIT, -(-e // UNIT), code[r["STATE"]]])

    tree = collections.defaultdict(lambda: collections.defaultdict(list))
    for s in sorted(pools): tree[pools[s]["taxon"]][pools[s]["family"]].append(s)
    data = dict(title=a.title, unit=UNIT, chroms=chroms, clen=[-(-clen[c] // UNIT) for c in chroms],
                pools=pools, tree={t: dict(sorted(v.items())) for t, v in sorted(tree.items())})
    html = TEMPLATE.replace("__TITLE__", a.title).replace("__DATA__", json.dumps(data, separators=(",", ":")))
    with open(a.out, "w") as f: f.write(html)
    nf = sum(len(v) for v in tree.values())
    print(f"pools: {len(pools)}  families: {nf}  taxa: {len(tree)}  WROTE: {a.out}")


TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Introgression viewer</title>
<style>
:root{--bg:#fcfcfb;--ink:#1d1d1f;--ink2:#555;--line:#ddd;--b73:#3A6FB5;--seg:#DDAA33;--don:#BB5566}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
main{max-width:1200px;margin:0 auto;padding:16px}
h1{font-size:20px;margin:0 0 4px} .sub{color:var(--ink2);margin:0 0 14px}
.bar{display:flex;flex-wrap:wrap;gap:10px;align-items:end;margin-bottom:12px}
label{display:flex;flex-direction:column;font-size:12px;color:var(--ink2);gap:3px}
select,input,button{font:inherit;padding:6px 8px;border:1px solid #bbb;border-radius:6px;background:#fff;color:var(--ink)}
select{min-width:150px} button{cursor:pointer} button:hover{background:#f0f0f0}
.legend{display:flex;flex-wrap:wrap;gap:16px;margin:4px 0 10px;font-size:13px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.sw{width:14px;height:14px;border-radius:3px;display:inline-block}
#fig{overflow-x:auto;border:1px solid var(--line);border-radius:8px;background:#fff}
#tip{position:fixed;pointer-events:none;background:#1d1d1f;color:#fff;padding:6px 8px;border-radius:6px;font-size:12px;display:none;z-index:9;white-space:nowrap}
table{border-collapse:collapse;margin-top:14px;font-size:13px;width:100%}
th,td{padding:5px 8px;border-bottom:1px solid var(--line);text-align:right} th:first-child,td:first-child{text-align:left}
tr.sel td{background:#f3f6fb} tr[data-s]{cursor:pointer}
.noisy{color:#9a3b00;font-weight:600}
.note{color:var(--ink2);font-size:12.5px;margin-top:10px}
</style></head><body><main>
<h1>__TITLE__</h1>
<p class="sub">Pick a taxon and a family. "All pools" shows every pool of the family on each chromosome, so shared segments line up.</p>
<div class="bar">
 <label>Taxon<select id="tax"></select></label>
 <label>Family<select id="fam"></select></label>
 <label>Pool<select id="pool"></select></label>
 <label>Find pool<input id="find" list="allpools" placeholder="e.g. PN12_SID1106" size="16"></label>
 <datalist id="allpools"></datalist>
 <button id="svg">Download SVG</button><button id="png">Download PNG</button>
</div>
<div class="legend">
 <span><i class="sw" style="background:var(--b73)"></i>B73: all plants AA</span>
 <span><i class="sw" style="background:var(--seg)"></i>SEG: family segregates, plants AA, AB or BB (1:2:1)</span>
 <span><i class="sw" style="background:var(--don)"></i>DONOR: all plants BB</span>
</div>
<div id="fig"></div>
<table id="tab"></table>
<p class="note">Each sample is a pool of 6 BC2S3 plants from one family. States are pool states from the k-mer HMM (steps 20–21), in 100-kb windows. SEG is not "het": in a SEG region the plants of the pool are ¼ AA, ½ AB and ¼ BB. Donor share = DONOR + SEG/2. NOISY pools (σ &gt; 0.15) are low-depth and not trusted.</p>
</main><div id="tip"></div>
<script>
const D=__DATA__;
const COL=["#3A6FB5","#DDAA33","#BB5566"], NAME=["B73","SEG","DONOR"];
const $=id=>document.getElementById(id), tax=$("tax"), fam=$("fam"), pool=$("pool"), tip=$("tip");
const NS="http://www.w3.org/2000/svg", maxL=Math.max(...D.clen);
const pct=x=>(100*x).toFixed(1)+"%", mb=u=>(u*D.unit/1e6).toFixed(1);
function opt(sel,vals,lab){sel.innerHTML="";vals.forEach((v,i)=>{const o=document.createElement("option");o.value=v;o.textContent=lab?lab[i]:v;sel.appendChild(o)})}
function el(n,a,p){const e=document.createElementNS(NS,n);for(const k in a)e.setAttribute(k,a[k]);if(p)p.appendChild(e);return e}
opt(tax,Object.keys(D.tree));
Object.keys(D.pools).sort().forEach(s=>{const o=document.createElement("option");o.value=s;$("allpools").appendChild(o)});
function fillFam(){const f=Object.keys(D.tree[tax.value]);opt(fam,f,f.map(x=>`${x} (${D.tree[tax.value][x].length})`))}
function fillPool(){const p=D.tree[tax.value][fam.value];opt(pool,["ALL",...p],["All pools",...p.map(s=>s+(D.pools[s].flag!="ok"?" (NOISY)":""))])}
function seg(g,x,y,w,h,st,s,c,a,b){const r=el("rect",{x,y,width:Math.max(w,0.8),height:h,fill:COL[st]},g);
 if(s)r.addEventListener("mousemove",ev=>{tip.style.display="block";tip.style.left=(ev.clientX+12)+"px";tip.style.top=(ev.clientY+12)+"px";
  tip.textContent=`${s} · ${D.chroms[c]}:${mb(a)}–${mb(b)} Mb · ${NAME[st]} (${(+mb(b-a)).toFixed(1)} Mb)`});
 r.addEventListener("mouseleave",()=>tip.style.display="none");}
function track(g,s,c,x0,y,h,sc){seg(g,x0,y,D.clen[c]*sc,h,0,s,c,0,D.clen[c]);
 D.pools[s].t.forEach(t=>{if(t[0]==c)seg(g,x0+t[1]*sc,y,(t[2]-t[1])*sc,h,t[3],s,c,t[1],t[2])})}
function axis(svg,x0,y,sc){el("line",{x1:x0,x2:x0+maxL*sc,y1:y,y2:y,stroke:"#999"},svg);
 for(let m=0;m*1e6/D.unit<=maxL;m+=50){const x=x0+m*1e6/D.unit*sc;el("line",{x1:x,x2:x,y1:y,y2:y+5,stroke:"#999"},svg);
  const t=el("text",{x,y:y+18,"text-anchor":"middle","font-size":11,fill:"#555"},svg);t.textContent=m}
 const t=el("text",{x:x0+maxL*sc/2,y:y+36,"text-anchor":"middle","font-size":12,fill:"#333"},svg);t.textContent="Position (Mb)"}
function draw(){
 const ps=D.tree[tax.value][fam.value], one=pool.value!="ALL", show=one?[pool.value]:ps;
 const W=1100,x0=one?60:150,sc=(W-x0-30)/maxL, svg=el("svg",{xmlns:NS,width:W,"font-family":"system-ui,sans-serif"});
 el("rect",{x:0,y:0,width:"100%",height:"100%",fill:"#fff"},svg);
 const p0=D.pools[show[0]];
 const ttl=el("text",{x:W/2,y:26,"text-anchor":"middle","font-size":16,"font-weight":700,fill:"#1d1d1f"},svg);
 ttl.textContent=one?`${show[0]} (${p0.taxon}, ${p0.family}) — donor share ${pct(p0.share)}${p0.flag!="ok"?"  [NOISY]":""}`
   :`${fam.value} (${tax.value}) — ${ps.length} pools`;
 let y=50;
 D.chroms.forEach((c,ci)=>{
  if(one){const h=24;const t=el("text",{x:x0-8,y:y+h/2+4,"text-anchor":"end","font-size":12,"font-weight":600,fill:"#333"},svg);t.textContent=c;
   track(svg,show[0],ci,x0,y,h,sc);y+=h+12;}
  else{const h=9;const t=el("text",{x:8,y:y+11,"font-size":12,"font-weight":700,fill:"#333"},svg);t.textContent=c;y+=16;
   show.forEach(s=>{const l=el("text",{x:x0-8,y:y+h-1,"text-anchor":"end","font-size":10,fill:D.pools[s].flag!="ok"?"#9a3b00":"#555"},svg);
    l.textContent=s+(D.pools[s].flag!="ok"?" *":"");track(svg,s,ci,x0,y,h,sc);y+=h+2});y+=10;}
 });
 axis(svg,x0,y+4,sc);
 // legend inside the figure, so downloads carry it
 let lx=x0,ly=y+64;[["B73: all plants AA",0],["SEG: plants AA, AB or BB (1:2:1)",1],["DONOR: all plants BB",2]].forEach(([t,k])=>{
  el("rect",{x:lx,y:ly-10,width:12,height:12,rx:2,fill:COL[k]},svg);const e=el("text",{x:lx+17,y:ly,"font-size":12,fill:"#333"},svg);e.textContent=t;lx+=t.length*6.6+40});
 if(!one&&ps.some(s=>D.pools[s].flag!="ok")){const e=el("text",{x:lx,y:ly,"font-size":12,fill:"#9a3b00"},svg);e.textContent="* NOISY pool"}
 svg.setAttribute("height",ly+16);svg.setAttribute("viewBox",`0 0 ${W} ${ly+16}`);
 $("fig").innerHTML="";$("fig").appendChild(svg);table(ps,one?show[0]:null);
 try{history.replaceState(null,"","#"+encodeURIComponent(tax.value+"|"+fam.value+"|"+pool.value))}catch(e){}
}
function table(ps,sel){const T=$("tab");T.innerHTML="<tr><th>Pool</th><th>Flag</th><th>σ</th><th>Donor share</th><th>SEG</th><th>DONOR</th><th>Tracts</th></tr>";
 ps.forEach(s=>{const p=D.pools[s],r=document.createElement("tr");r.dataset.s=s;if(s==sel)r.className="sel";
  r.innerHTML=`<td>${s}</td><td class="${p.flag!="ok"?"noisy":""}">${p.flag}</td><td>${p.sigma.toFixed(3)}</td><td>${pct(p.share)}</td><td>${pct(p.seg)}</td><td>${pct(p.donor)}</td><td>${p.ntr}</td>`;
  r.onclick=()=>{pool.value=s;draw()};T.appendChild(r)});}
function go(s){const p=D.pools[s];if(!p)return;tax.value=p.taxon;fillFam();fam.value=p.family;fillPool();pool.value=s;draw()}
tax.onchange=()=>{fillFam();fillPool();draw()}; fam.onchange=()=>{fillPool();draw()}; pool.onchange=draw;
$("find").onchange=e=>go(e.target.value.trim());
function name(ext){return (pool.value=="ALL"?fam.value:pool.value)+"_introgression."+ext}
function svgText(){return new XMLSerializer().serializeToString($("fig").querySelector("svg"))}
function dl(href,n){const a=document.createElement("a");a.href=href;a.download=n;document.body.appendChild(a);a.click();a.remove()}
$("svg").onclick=()=>dl(URL.createObjectURL(new Blob([svgText()],{type:"image/svg+xml"})),name("svg"));
$("png").onclick=()=>{const s=$("fig").querySelector("svg"),w=+s.getAttribute("width"),h=+s.getAttribute("height"),img=new Image();
 img.onload=()=>{const c=document.createElement("canvas");c.width=w*2;c.height=h*2;const x=c.getContext("2d");x.scale(2,2);x.drawImage(img,0,0);dl(c.toDataURL("image/png"),name("png"))};
 img.src="data:image/svg+xml;charset=utf-8,"+encodeURIComponent(svgText())};
let h="";try{h=decodeURIComponent(location.hash.slice(1))}catch(e){}
const [ht,hf,hp]=h.split("|");
if(D.tree[ht]){tax.value=ht;fillFam();if(D.tree[ht][hf]){fam.value=hf;fillPool();if(hp)pool.value=hp}else fillPool()}else{fillFam();fillPool()}
draw();
</script></body></html>
"""

if __name__ == "__main__":
    main()
