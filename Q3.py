from vpython import *
import numpy as np

win = 500
Natoms = 500
L = 1
gray = color.gray(0.7)
mass = 4E-3/6E23
Ratom = 0.03
k = 1.4E-23
T = 400
dt = 1E-5

# Tapa superior mobil
Ly = L
M = (2/3)*10**(-21)
ptotal = 0
posy = [L,L,0]

def verlet(y, p, M, dt):
    y[2] = 2*y[1] - y[0] + ((2*p)/(M*dt)-9.81)*dt**2
    y[0] = y[1]
    y[1] = y[2]
    return y[2]

def nova_pos_cara():
    global Ly, boxtop, vert1, vert2, vert3, vert4
    d1 = Ly + Ratom
    d = L/2 + Ratom
    boxtop.clear()
    boxtop.append([vector(-d, d1, -d), vector(-d, d1, d),
                   vector(d, d1, d), vector(d, d1, -d), vector(-d, d1, -d)])
    vert1.clear()
    vert1.append([vector(-d, 0, -d), vector(-d, d1, -d)])
    vert2.clear()
    vert2.append([vector(-d, 0, d), vector(-d, d1, d)])
    vert3.clear()
    vert3.append([vector(d, 0, d), vector(d, d1, d)])
    vert4.clear()
    vert4.append([vector(d, 0, -d), vector(d, d1, -d)])

file = open('dades_Ly_isoterm.txt', 'a')

animation = canvas(width=win, height=win, align='left')
animation.range = 3*L
animation.title = 'Proceso isotérmico: Gas de esferas duras'
animation.camera.pos = vector(0, L/2, 0)

d = L/2 + Ratom
d2 = L + Ratom
d0 = -Ratom
r = 0.005

boxbottom = curve(color=gray, radius=r)
boxbottom.append([vector(-d,d0,-d), vector(-d,d0,d), vector(d,d0,d), vector(d,d0,-d), vector(-d,d0,-d)])
boxtop = curve(color=gray, radius=r)
boxtop.append([vector(-d,d2,-d), vector(-d,d2,d), vector(d,d2,d), vector(d,d2,-d), vector(-d,d2,-d)])
vert1 = curve(color=gray, radius=r)
vert2 = curve(color=gray, radius=r)
vert3 = curve(color=gray, radius=r)
vert4 = curve(color=gray, radius=r)
vert1.append([vector(-d,d0,-d), vector(-d,d2,-d)])
vert2.append([vector(-d,d0,d), vector(-d,d2,d)])
vert3.append([vector(d,d0,d), vector(d,d2,d)])
vert4.append([vector(d,d0,-d), vector(d,d2,-d)])

Atoms = []
p = []
apos = []
pavg = sqrt(2*mass*1.5*k*T)

for i in range(Natoms):
    x = L*random()-L/2
    y = L*random()
    z = L*random()-L/2
    if i == 0:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z))
    theta = pi*random()
    phi = 2*pi*random()
    px = pavg*sin(theta)*cos(phi)
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    p.append(vector(px,py,pz))

deltav = 100
def barx(v): return int(v/deltav)
nhisto = int(4500/deltav)
histo = [0.0]*nhisto
histo[barx(pavg/mass)] = Natoms

gg = graph(width=win, height=0.4*win, xmax=3000, align='left',
           xtitle='velocidad (m/s)', ytitle='nº de átomos', ymax=Natoms*deltav/1000)

theory = gcurve(color=color.blue, width=2)
dv = 10
for v in range(0,3001+dv,dv):
    theory.plot(v, (deltav/dv)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) *exp(-0.5*mass*(v**2)/(k*T))*(v**2)*dv )

accum = [[deltav*(i+0.5),0] for i in range(int(3000/deltav))]
vdist = gvbars(color=color.red, delta=deltav)

def interchange(v1, v2):
    barx1 = barx(v1)
    barx2 = barx(v2)
    if barx1 == barx2 or barx1 >= len(histo) or barx2 >= len(histo): return
    histo[barx1] -= 1
    histo[barx2] += 1

def checkCollisions():
    hitlist = []
    r2 = (2*Ratom)**2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i):
            aj = apos[j]
            dr = ai - aj
            if mag2(dr) < r2: hitlist.append([i,j])
    return hitlist

# Aplicació de la font termica
def resample_velocity(i):
    theta = pi*random()
    phi = 2*pi*random()
    speed = np.random.normal(loc=0, scale=np.sqrt(k*T/mass))
    px = speed*sin(theta)*cos(phi)
    py = speed*sin(theta)*sin(phi)
    pz = speed*cos(theta)
    p[i] = vector(px, py, pz)

nhisto = 0

pv_graph = graph(width=win, height=0.4*win, align='left',
                 xtitle='tiempo (s)', ytitle='PV', title='PV vs tiempo (isotérmico)', fast=False)
pv_curve = gcurve(color=color.green)

# Variables grafic PV
t = 0.0
pv_interval = 0.01  
pv_timer = 0.0


while True:
    rate(300)
    for i in range(len(accum)):
        accum[i][1] = (nhisto*accum[i][1] + histo[i])/(nhisto+1)
    if nhisto % 10 == 0:
        vdist.data = accum
    nhisto += 1


    if nhisto % 100 == 0:
        for i in range(Natoms):
            if random() < 0.05:
                resample_velocity(i)

    for i in range(Natoms):
        Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt

    nova_pos_cara()
    hitlist = checkCollisions()

    for ij in hitlist:
        i, j = ij
        ptot = p[i] + p[j]
        posi, posj = apos[i], apos[j]
        vi, vj = p[i]/mass, p[j]/mass
        vrel = vj - vi
        if vrel.mag2 == 0: continue
        rrel = posi - posj
        if rrel.mag > Ratom: continue
        dx = dot(rrel, vrel.hat)
        dy = cross(rrel, vrel.hat).mag
        alpha = asin(dy/(2*Ratom))
        d = (2*Ratom)*cos(alpha) - dx
        deltat = d/vrel.mag
        posi -= vi*deltat
        posj -= vj*deltat
        mtot = 2*mass
        pcmi = p[i] - ptot*mass/mtot
        pcmj = p[j] - ptot*mass/mtot
        rrel = norm(rrel)
        pcmi -= 2*pcmi.dot(rrel)*rrel
        pcmj -= 2*pcmj.dot(rrel)*rrel
        p[i] = pcmi + ptot*mass/mtot
        p[j] = pcmj + ptot*mass/mtot
        apos[i] = posi + (p[i]/mass)*deltat
        apos[j] = posj + (p[j]/mass)*deltat
        interchange(vi.mag, p[i].mag/mass)
        interchange(vj.mag, p[j].mag/mass)
    
    for i in range(Natoms):
        loc = apos[i]
        ptotal = 0
        if abs(loc.x) > L/2:
            p[i].x *= -1
        if abs(loc.z) > L/2:
            p[i].z *= -1
        if loc.y < 0:
            p[i].y = abs(p[i].y)
        if loc.y > Ly:
            apos[i].y = Ly
            p[i].y = -abs(p[i].y)
            ptotal += abs(p[i].y)

        Ly = verlet(posy, ptotal, M, dt)
        

# Per aturar la simuació
        Ly_max =  1.1 
        if Ly > Ly_max:
            Ly = Ly_max
            posy[0] = Ly
            posy[1] = Ly
            posy[2] = Ly


        if i % 600 == 0:
            file.write(str(Ly)+'\n')
        
        # Grafica PV
        t += dt
        pv_timer += dt
        
        V = L * L * Ly
        A = L * L
        P = ptotal / (dt * A)
        
        if pv_timer >= pv_interval:
            pv_curve.plot(t, P * V)
            pv_timer = 0.0