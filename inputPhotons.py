import numpy as n
import pylab as p

nuObs = 4.5e9  # for Radio
# nuObs=4.5e14  #for R filter
# nuObs=4.5e17  #for Xrays
totPh = 1000

Eiso = 1e53
Gamma0 = 300.
ee = 0.3
eB = 0.1
pel = 2.3
thetaJet = 5*n.pi/180
nISM = 1
tEng = 20.

c = 3e10
mp = 1.67e-24
me = mp/1836
kB = 1.38e-16
sT = 6.65e-25
e = 4.8e-10

rES = (3*Eiso/4/n.pi/nISM/mp/c**2/Gamma0**2)**(1/3)
rES2 = (3*Eiso/4/n.pi/nISM/mp/c*tEng)**.25
lrES = n.log10(rES)
r = 10**n.arange(lrES-2, lrES+2.2, .01)
M0 = Eiso/Gamma0/c**2
if rES < rES2:
    M0 = r*0+M0
    jj = n.where(r <= rES2)
    M0[jj] = M0[jj]*(r[jj]/rES2)
effe = 4*n.pi/3*r**3*nISM*mp/M0
Gamma = (Gamma0+effe)/n.sqrt(1+2*Gamma0*effe+effe**2)
beta = n.sqrt(1-1/Gamma**2)
tLab = r.copy()*0
for i in range(r.size):
    tLab[i] = r.min()/beta[0]/c+n.trapz(1/beta[:i+1]/c, x=r[:i+1])
tObs = tLab-r/c

gammaInj = 1+(pel-2)/(pel-1)*ee*(Gamma-1)*1836
B = n.sqrt(32*n.pi*eB*mp*c**2*nISM*(Gamma**2-1))
gammaCool = 15*n.pi*me*c**2*n.sqrt(Gamma**2-1)/sT/B**2/r+1.
tauAbs = 5*e/3.*nISM*r/B/(gammaInj)**5
kk = n.where(gammaCool < gammaInj)
if kk[0].size > 0:
    tauAbs[kk] = 5*e/3.*nISM*r[kk]/B[kk]/(gammaCool[kk])**5
nuPrimeInj = 3*gammaInj**2*e*B/16/me/c
nuPrimeCool = 3*gammaCool**2*e*B/16/me/c
nuPrimeAbs = nuPrimeInj*tauAbs**(3./5.)
if kk[0].size > 0:
    nuPrimeAbs[kk] = nuPrimeCool[kk]*tauAbs[kk]**(3./5.)


# nuPrimeCool=n.maximum(nuPrimeCool,nuPrimeAbs)

gammaAbs = n.sqrt(16*me*c*nuPrimeAbs/3/e/B)

nuPrime = nuObs/2/Gamma
delta = 2*Gamma
IprimeNu = r.copy()*0
dLnu = r.copy()*0

# case (1) abs<inj<cool
# times at which the nus are in standard order
jj = n.where((nuPrimeCool >= nuPrimeInj) & (nuPrimeInj >= nuPrimeAbs))
case1 = jj[0].size
if case1 > 0:
    Pprime = 4./3.*sT*c*B**2/8/n.pi*(gammaInj**2-1)
    IprimePeak = Pprime*nISM*r/nuPrimeInj
    kk = n.where(nuPrime <= nuPrimeAbs)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk] * \
            (nuPrime[kk]/nuPrimeAbs[kk])**2 * \
            (nuPrimeAbs[kk]/nuPrimeInj[kk])**(1./3)
    kk = n.where((nuPrime < nuPrimeInj) & (nuPrime > nuPrimeAbs))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeInj[kk])**(1./3.)
    kk = n.where((nuPrime >= nuPrimeInj) & (nuPrime < nuPrimeCool))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk] * \
            (nuPrime[kk]/nuPrimeInj[kk])**(-(pel-1)/2.)
    kk = n.where(nuPrime >= nuPrimeCool)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrimeInj[kk]/nuPrimeCool[kk])**(
            (pel-1)/2.)*(nuPrime[kk]/nuPrimeCool[kk])**(-pel/2.)
    dLnu[jj] = IprimeNu[jj]*delta[jj]**3*r[jj]**2

# case (2) abs<cool<inj
jj = n.where((nuPrimeAbs <= nuPrimeCool) & (nuPrimeCool <= nuPrimeInj))
case2 = jj[0].size
if case2 > 0:
    Pprime = 4./3.*sT*c*B**2/8/n.pi*(gammaCool**2-1)
    IprimePeak = Pprime*nISM*r/nuPrimeCool
    kk = n.where(nuPrime <= nuPrimeAbs)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk] * \
            (nuPrime[kk]/nuPrimeAbs[kk])**2 * \
            (nuPrimeAbs[kk]/nuPrimeCool[kk])**(1./3)
    kk = n.where((nuPrime <= nuPrimeCool) & (nuPrime > nuPrimeAbs))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeCool[kk])**(1./3.)
    kk = n.where((nuPrime <= nuPrimeInj) & (nuPrime > nuPrimeCool))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeCool[kk])**(-1/2.)
    kk = n.where(nuPrime > nuPrimeInj)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrimeCool[kk]/nuPrimeInj[kk])**(
            1./2.)*(nuPrime[kk]/nuPrimeInj[kk])**(-pel/2.)
    dLnu[jj] = IprimeNu[jj]*delta[jj]**3*r[jj]**2

# case (3) inj<abs<cool
jj = n.where((nuPrimeInj <= nuPrimeAbs) & (nuPrimeAbs <= nuPrimeCool))
case3 = jj[0].size
if case3 > 0:
    Pprime = 4./3.*sT*c*B**2/8/n.pi*(gammaAbs**2-1)
    IprimePeak = Pprime*nISM*r/nuPrimeAbs*(gammaInj/gammaAbs)**(pel-1)
    kk = n.where(nuPrime <= nuPrimeInj)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk] * \
            (nuPrime[kk]/nuPrimeInj[kk])**2 * \
            (nuPrimeInj[kk]/nuPrimeAbs[kk])**(5./2)
    kk = n.where((nuPrime <= nuPrimeAbs) & (nuPrime > nuPrimeInj))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeAbs[kk])**(5./2.)
    kk = n.where((nuPrime <= nuPrimeCool) & (nuPrime > nuPrimeAbs))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk] * \
            (nuPrime[kk]/nuPrimeAbs[kk])**(-(pel-1)/2.)
    kk = n.where(nuPrime > nuPrimeCool)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrimeCool[kk]/nuPrimeAbs[kk]
                                       )**(-(pel-1)/2.)*(nuPrime[kk]/nuPrimeCool[kk])**(-pel/2.)
    dLnu[jj] = IprimeNu[jj]*delta[jj]**3*r[jj]**2

# case (4) inj<cool<abs
jj = n.where((nuPrimeInj <= nuPrimeCool) & (nuPrimeCool <= nuPrimeAbs))
case4 = jj[0].size
if case4 > 0:
    Pprime = 4./3.*sT*c*B**2/8/n.pi*(gammaAbs**2-1)
    IprimePeak = Pprime*nISM*r/nuPrimeAbs*(gammaInj/gammaAbs)**(pel-1)
    kk = n.where(nuPrime <= nuPrimeInj)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk] * \
            (nuPrime[kk]/nuPrimeInj[kk])**2 * \
            (nuPrimeInj[kk]/nuPrimeAbs[kk])**(5./2)
    kk = n.where((nuPrime <= nuPrimeAbs) & (nuPrime > nuPrimeInj))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeAbs[kk])**(5./2.)
    kk = n.where(nuPrime > nuPrimeAbs)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeAbs[kk])**(-pel/2.)
    dLnu[jj] = IprimeNu[jj]*delta[jj]**3*r[jj]**2

# case (5) cool<abs<inj
jj = n.where((nuPrimeCool <= nuPrimeAbs) & (nuPrimeAbs <= nuPrimeInj))
case5 = jj[0].size
if case5 > 0:
    Pprime = 4./3.*sT*c*B**2/8/n.pi*(gammaAbs**2-1)
    IprimePeak = Pprime*nISM*r/nuPrimeAbs
    kk = n.where(nuPrime <= nuPrimeAbs)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeAbs[kk])**2
    kk = n.where((nuPrime <= nuPrimeInj) & (nuPrime > nuPrimeAbs))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeAbs[kk])**(-1/2.)
    kk = n.where(nuPrime > nuPrimeInj)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk] * \
            (nuPrimeAbs[kk]/nuPrimeInj[kk])**(1/2) * \
            (nuPrime[kk]/nuPrimeInj[kk])**(-pel/2.)
    dLnu[jj] = IprimeNu[jj]*delta[jj]**3*r[jj]**2

# case (6) cool<inj<abs
jj = n.where((nuPrimeCool <= nuPrimeInj) & (nuPrimeInj <= nuPrimeAbs))
case6 = jj[0].size
if case6 > 0:
    Pprime = 4./3.*sT*c*B**2/8/n.pi*(gammaAbs**2-1)
    IprimePeak = Pprime*nISM*r/nuPrimeAbs*(gammaInj/gammaAbs)**(pel-1)
    kk = n.where(nuPrime <= nuPrimeInj)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk] * \
            (nuPrime[kk]/nuPrimeInj[kk])**2 * \
            (nuPrimeInj[kk]/nuPrimeAbs[kk])**(5/2)
    kk = n.where((nuPrime <= nuPrimeAbs) & (nuPrime > nuPrimeInj))
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeAbs[kk])**(5/2.)
    kk = n.where(nuPrime > nuPrimeAbs)
    if kk[0].size > 0:
        IprimeNu[kk] = IprimePeak[kk]*(nuPrime[kk]/nuPrimeAbs[kk])**(-pel/2.)
    dLnu[jj] = IprimeNu[jj]*delta[jj]**3*r[jj]**2


totArea = n.trapz(dLnu, x=tLab)

outf = open('testphotons.txt', 'w')
for i in range(tLab.size-1):
    numPh = n.random.poisson(totPh/(tLab.size-1))
    thisArea = (dLnu[i+1]+dLnu[i])*(tObs[i+1]-tObs[i])/2
    weight = thisArea/totArea
    phi = n.random.rand(numPh)*2*n.pi
    xmin = n.cos(thetaJet)
    theta = n.arccos(n.random.rand(numPh)*(1-xmin)+xmin)
    phiv = n.random.rand(numPh)*2*n.pi
    thetav = n.arccos(n.random.rand(numPh)*2-1)
    for j in range(numPh):
        time = (tLab[i+1]+tLab[i])/2
        xPos = (r[i+1]+r[i])/2*n.sin(theta[j])*n.cos(phi[j])
        yPos = (r[i+1]+r[i])/2*n.sin(theta[j])*n.sin(phi[j])
        zPos = (r[i+1]+r[i])/2*n.cos(theta[j])
        pos = (r[i+1]+r[i])/2
        vXp = n.sin(thetav[j])*n.cos(phiv[j])
        vYp = n.sin(thetav[j])*n.sin(phiv[j])
        vZp = n.cos(thetav[j])
        gg = Gamma[i]
        bb = n.sqrt(1-1/gg**2)
        boost = n.array([n.sin(theta[j])*n.cos(phi[j]),
                        n.sin(theta[j])*n.sin(phi[j]), n.cos(theta[j])])*bb
        varr = n.array([vXp, vYp, vZp])
        varrp = (varr/gg+boost+gg/(gg+1)*n.dot(varr, boost)*boost) / \
            (1+n.dot(varr, boost))
        vX = varrp[0]
        vY = varrp[1]
        vZ = varrp[2]
        rnd = n.random.rand()
        if rnd > 0.5:
            zPos = -zPos
            vZ = -vZ
        strout = str(time)+'\t'+str(xPos)+'\t'+str(yPos) + \
            '\t'+str(zPos)+'\t'+str(vX)+'\t'+str(vY)+'\t'+str(vZ) +\
            '\t'+str(nuObs)+'\t'+str(weight)+' \n'
        outf.write(strout)

print('There are ', r.size-case1-case2-case3-case4-case5-case6, 'faulty zones')


outf.close()
p.figure(1)
p.clf()
p.loglog(tObs, nuPrimeAbs, label='Absorption')
p.loglog(tObs, nuPrimeInj, label='Injection')
p.loglog(tObs, nuPrimeCool, label='Cooling')
p.legend()
p.show()


t, x, y, z, vx, vy, vz, nu, w = n.loadtxt('testphotons.txt', unpack=True)


r = n.sqrt(x**2+y**2+z**2)
tObserved = t-n.abs(r)/3e10

tBin = 10**(0.5*(n.log10(tObs[1:])+n.log10(tObs[:-1])))
lCurve = tBin.copy()*0
for i in range(tBin.size):
    jj = n.where((tObserved >= tObs[i]) & (tObserved < tObs[i+1]))
    lCurve[i] = n.sum(w[jj])/(tObs[i+1]-tObs[i])

factor1 = n.sum((tObs[1:]-tObs[:-1])*(dLnu[1:]+dLnu[:-1])/2)
factor2 = n.sum((tBin[1:]-tBin[:-1])*(lCurve[1:]+lCurve[:-1])/2)
factor = factor1/factor2*(1-n.cos(thetaJet))
print('The factor is:', factor)

w = w*factor

outarray = n.zeros([8, t.size])
outarray[0, :] = t
outarray[1, :] = x
outarray[2, :] = y
outarray[3, :] = z
outarray[4, :] = vx
outarray[5, :] = vy
outarray[6, :] = vz
outarray[7, :] = w

n.savetxt('testphotonsWeighs.dat', n.transpose(outarray))

for i in range(tBin.size):
    jj = n.where((tObserved >= tObs[i]) & (tObserved < tObs[i+1]))
    lCurve[i] = n.sum(w[jj])/(tObs[i+1]-tObs[i])

lCurve = lCurve/(1-n.cos(thetaJet))

thAccept = 5
th = n.arctan2(n.sqrt(vx**2+vy**2), n.abs(vz))*180/n.pi
jj = n.where(th < thAccept)
print(jj[0].size, 'out of', t.size)
t = t[jj]
x = x[jj]
y = y[jj]
z = z[jj]
w = w[jj]
nu = nu[jj]

tObserved = t-n.abs(z)/3e10
lCurveObs = tBin.copy()*0
for i in range(tBin.size):
    jj = n.where((tObserved >= tObs[i]) & (tObserved < tObs[i+1]))
    lCurveObs[i] = n.sum(w[jj])/(tObs[i+1]-tObs[i])

lCurveObs = lCurveObs/(1-n.cos(thAccept*n.pi/180))


print(
    "output format: t [s], x [cm], y [cm], z [cm], vx [c], vy [c], vz [c], w")
# p.figure(2)
# p.clf()
# p.loglog(tObs*3,dLnu/3)
# p.loglog(tBin,lCurve)
# p.loglog(tBin, lCurveObs)
# p.show()
