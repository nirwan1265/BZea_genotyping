import random
random.seed(4)
W,Hh=1100,1110
B73="#2850C0"; SEG="#6A2C9E"; DON="#C0202A"; OK="#2E8B57"; GRAY="#9AA3AD"; INK="#1F2933"; SOFT="#5B6670"
o=[]
def t(x,y,s,size=15,w="normal",fill=INK,anchor="start",style=""):
    o.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{w}" fill="{fill}" text-anchor="{anchor}" {style}>{s}</text>')
def rect(x,y,w,h,fill,rx=0,stroke="none",sw=0,extra=""):
    o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" {extra}/>')
def line(x1,y1,x2,y2,c,w=2,extra=""):
    o.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{c}" stroke-width="{w}" {extra}/>')
def circ(x,y,r,fill,stroke="none",sw=0): o.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
def panel(y,num,title,sub):
    rect(20,y,W-40,Hh and 0 or 0,"none")
    circ(44,y+14,14,INK); t(44,y+19,str(num),14,"bold","#fff","middle")
    t(68,y+19,title,18,"bold"); t(68,y+40,sub,13.5,fill=SOFT)
def stick(x1,x2,y,c,ball_l=True,ball_r=True,w=4):
    line(x1,y,x2,y,c,w,'stroke-linecap="round"')
    if ball_l: circ(x1,y,4.5,c)
    if ball_r: circ(x2,y,4.5,c)
def check(x,y): o.append(f'<path d="M{x-6},{y} l4,5 l9,-10" fill="none" stroke="{OK}" stroke-width="3" stroke-linecap="round"/>')
def cross(x,y,c=DON): o.append(f'<path d="M{x-6},{y-6} l12,12 M{x+6},{y-6} l-12,12" stroke="{c}" stroke-width="3" stroke-linecap="round"/>')

o.append(f'<rect x="0" y="0" width="{W}" height="{Hh}" fill="#FFFFFF"/>')
t(W/2,38,"Calling introgressions with k-mers (no read mapping)",22,"bold",anchor="middle")
t(W/2,62,"one 100-kb window of the B73 genome, one pool of 6 BC2S3 plants",14,fill=SOFT,anchor="middle")

# ---------- panel 1 ----------
y=90; panel(y,1,"Make the B73 k-mer list (step 16)","Cut B73 into every 31-letter word (k-mer). Keep only words that occur ONCE in the whole genome.")
gy=y+95
rect(80,gy-6,940,12,"#E6EAF0",6); t(70,gy+5,"B73",14,"bold",anchor="end")
rect(520,gy-6,190,12,"#C9CED6",6)
t(615,gy+30,"transposon / repeat",12.5,fill=SOFT,anchor="middle")
def row(xs,xe,c,klen=56,gap=14):
    for r_,off in enumerate([0,23,46]):
        x=xs+off
        while x+klen<=xe:
            stick(x,x+klen,gy-20-r_*11,c); x+=klen+gap
row(88,515,B73); row(525,705,GRAY,klen=44,gap=12); row(715,1020,B73)
stick(110,170,gy+58,B73); t(180,gy+63,"single-copy k-mer: kept (680 M of 817 M)",13.5)
stick(560,620,gy+58,GRAY); t(630,gy+63,"repeat k-mer: dropped (present in every genome)",13.5)

# ---------- panel 2 ----------
y=290; panel(y,2,"Look the k-mers up in each pool's reads (step 17)","Count all k-mers in the pool's raw reads (KMC). Reads are never mapped, so donor reads are not lost.")
# B73 region
by=y+128; t(80,by-62,"B73 region",14,"bold",fill=B73)
rect(80,by-6,420,12,"#E6EAF0",6)
for i,(a,b) in enumerate([(90,230),(150,300),(260,410),(330,480),(110,250),(380,495)]):
    ry=by+22+(i%3)*12; line(a,ry,b,ry,"#7C8794",5,'stroke-linecap="round"')
for i in range(6):
    x=100+i*66; stick(x,x+52,by-22,B73); check(x+26,by-42)
t(290,by+80,"every B73 k-mer is found",13.5,"bold",fill=OK,anchor="middle")
# donor region
dx=560; t(dx,by-62,"donor (teosinte) region",14,"bold",fill=DON)
rect(dx,by-6,460,12,"#E6EAF0",6)
snps=[dx+164,dx+290,dx+420]
for i,(a,b) in enumerate([(dx+10,dx+150),(dx+80,dx+230),(dx+190,dx+340),(dx+260,dx+410),(dx+320,dx+455),(dx+40,dx+170)]):
    ry=by+22+(i%3)*12; line(a,ry,b,ry,"#C9A0A3",5,'stroke-linecap="round"')
    for s in snps:
        if a<s<b: circ(s,ry,4,DON)
for s in snps: circ(s,by,5.5,DON)
for i in range(7):
    x=dx+10+i*64; a,b=x,x+52
    hit=any(a-2<=s<=b+2 for s in snps)
    stick(a,b,by-22,DON if hit else B73)
    (cross if hit else check)((a+b)/2,by-42)
t(dx+230,by+80,"k-mers that overlap a donor SNP or indel are missing",13.5,"bold",fill=DON,anchor="middle")
circ(dx+40,by+104,5,DON); t(dx+52,by+109,"teosinte SNP",12.5,fill=SOFT)
line(dx+170,by+104,dx+210,by+104,"#C9A0A3",5,'stroke-linecap="round"'); t(dx+218,by+109,"read",12.5,fill=SOFT)

# ---------- panel 3 ----------
y=568; panel(y,3,"Score the window (step 20)","r = B73 single-copy k-mers seen in the pool ÷ what a B73 pool of the same depth sees (fit of the B73 controls)")
base=y+200; x0=170; h=130
line(x0-20,base,x0+560,base,SOFT,1.5)
line(x0-20,base-h,x0+560,base-h,SOFT,1,'stroke-dasharray="5,5"'); t(x0-28,base-h+5,"1.0",13,fill=SOFT,anchor="end"); t(x0-28,base+5,"0",13,fill=SOFT,anchor="end")
for i,(lab,val,c,desc) in enumerate([("B73",1.0,B73,"all plants AA"),("SEG",0.62,SEG,"family segregates (AA, AB, BB)"),("DONOR",0.22,DON,"all plants BB")]):
    x=x0+30+i*180; rect(x,base-h*val,90,h*val,c,4)
    t(x+45,base-h*val-10,f"r ≈ {val:.1f}" if val!=1.0 else "r ≈ 1",14,"bold",anchor="middle")
    t(x+45,base+22,lab,15,"bold",fill=c,anchor="middle"); t(x+45,base+40,desc,12,fill=SOFT,anchor="middle")
t(x0+620,base-100,"Expected DONOR level per window comes",13.5,fill=INK)
t(x0+620,base-81,"from teosinte genome assemblies",13.5,fill=INK)
t(x0+620,base-62,"(share of B73 k-mers each taxon keeps,",13.5,fill=SOFT)
t(x0+620,base-43,"e.g. Zd ≈ 0.27, Zv ≈ 0.43; step 19).",13.5,fill=SOFT)
t(x0+620,base-18,"SEG lies in between: k-mers the donor lacks",13.5,fill=SOFT)
t(x0+620,base+1,"come only from the B73 plants in the pool.",13.5,fill=SOFT)

# ---------- panel 4 ----------
y=845; panel(y,4,"HMM along each chromosome (steps 20-21)","Windows → B73 / SEG / DONOR. Switches follow the NAM genetic map; switch rate ρ fitted by maximum likelihood.")
px0,px1,pyb,ph=110,1040,y+175,95
line(px0,pyb,px1,pyb,SOFT,1); line(px0,pyb-ph,px1,pyb-ph,SOFT,1,'stroke-dasharray="4,4"')
t(px0-8,pyb-ph+5,"1.0",12,fill=SOFT,anchor="end"); t(px0-8,pyb+5,"0",12,fill=SOFT,anchor="end"); t(px0-40,pyb-ph/2,"r",15,"bold",fill=SOFT)
n=150; states=[]
for i in range(n):
    s=0
    if 40<=i<70: s=1
    if 95<=i<125: s=2
    states.append(s)
lv={0:1.0,1:0.62,2:0.22}
for i,s in enumerate(states):
    v=max(0.02,min(1.25,lv[s]+random.gauss(0,0.09)))
    x=px0+(px1-px0)*i/(n-1); circ(x,pyb-ph*v,2.6,{0:B73,1:SEG,2:DON}[s])
sy=pyb+18; i=0
while i<n:
    j=i
    while j<n and states[j]==states[i]: j+=1
    xa=px0+(px1-px0)*i/(n-1); xb=px0+(px1-px0)*(j-1)/(n-1)
    rect(xa,sy,max(4,xb-xa+6),20,{0:B73,1:SEG,2:DON}[states[i]])
    i=j
t(px0,sy+42,"called states",12.5,fill=SOFT)
for k,(lab,c) in enumerate([("B73 (AA)",B73),("SEG (AB, segregating)",SEG),("DONOR (BB)",DON)]):
    rect(px0+200+k*240,sy+31,14,14,c,3); t(px0+220+k*240,sy+43,lab,13)
svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{Hh}" viewBox="0 0 {W} {Hh}" font-family="Helvetica, Arial, sans-serif">'+"".join(o)+"</svg>"
open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "kmer_method.svg"),"w").write(svg)
print("ok")
