# -*- coding: utf-8 -*-
from numpy import *
import math

g = 9.81

class HFrm(object):
    '''
    Class for apparatus frame hydrodynamic coefficients to be
    calculated; parameters are:
    R,m - frame radius
    H,m - underwater depth of frame center
    sgm,1/sec - circular friequancy'''
        
    def __init__(self,R,H,sgm):
        self.__R = R
        self.__H = H
        self.__sgm = sgm

    @staticmethod
    def rmb(f,a,b,tol=1.0e-6):
        '''
        I = romberg(f,a,b,tol=1.0e-6).
        Romberg intergration of f(x) from x = a to b.
        Returns the integral.
        '''

        def trapezoid(f,a,b,Iold,k):
            '''
            Inew = trapezoid(f,a,b,Iold,k).
            Recursive trapezoidal rule:
            Iold = Integral of f(x) from x = a to b computed by
            trapezoidal rule with 2?(k-1) panels.
            Inew = Same integral computed with 2?k panels.
            '''
            if k == 1:Inew = (f(a) + f(b))*(b - a)/2.0
            else:
                n = 2**(k -2 )   # Number of new points
                h = (b - a)/n    # Spacing of new points
                x = a + h/2.0
                sum = 0.0
                for i in range(n):
                    sum = sum + f(x)
                    x = x + h
                    Inew = (Iold + h*sum)/2.0
            return Inew

        def richardson(r,k):  # function definiting in body of another function
            for j in range(k-1,0,-1):
                const = 4.0**(k-j)
                r[j] = (const*r[j+1] - r[j])/(const - 1.0)
            return r

        r = zeros(21)
        r[1] = trapezoid(f,a,b,0.0,1)
        r_old = r[1]
        for k in range(2,21):
            r[k] = trapezoid(f,a,b,r[k-1],k)
            r = richardson(r,k)
            if abs(r[1]-r_old) < tol*max(abs(r[1]),1.0):
                return r[1]
            r_old = r[1]
        print "Romberg quadrature did not converge"

        # End of function rmb(f,a,b,tol=1.0e-6)

    @property
    def fB(self):
        eps = 1.e-6
        r = self.__R/self.__H
        l = self.__H*sqrt(1.-r*r)
        n = self.__sgm**2*l/g
        if n == 0.: return 0.
        t = r/(1.+sqrt(1.-r*r))
        t0 =-log(t)
        en = eps*exp(n)
        #FIXME:
        N = 10#int(log(1./sqrt(en) + sqrt(1.+1./en))/t0)+1
        #
        b=0.
        em = sg = 1.
        for i in xrange(1,N):
            em *= t
            ep=1./em
            s=.5*(ep-em)
            c=.5*(ep+em)
            s2=s*s
            u = exp(-n*c/s)/s2*sg
            b += u
            sg =-sg
        return pi*n*b

    @property
    def flmb(self):
        r=self.__R/self.__H
        l=self.__H*sqrt(1.-r*r)
        n=self.__sgm**2*l/g
        return sqrt(n)*self.fB**2    # fB()
        
    @property
    def f_minf(self):    # Function to calculate added mass coefs at infinite freq
        eps = 1.e-6
        r = self.__R/self.__H
        if r < eps: r = eps
        t = r/(1.+sqrt(1.-r*r))
        if t < eps : t = eps
        t0 =-log(t)
        N = 100#int(round(log(1./sqrt(eps)+sqrt(1.+1./eps))/t0))+1  
        b = 0.
        em = sg = 1.
        for i in range(1,N):
            em *= t
            ep=1./em
            s=.5*(ep-em)
            s2=s*s
            u=1./s2*sg
            if i != 1: u *= 2.
            b += u
            sg=-sg
        M = pi*b
        return M/pi*sinh(t0)**2    # new value fminf:   /pi*sinh(t0)**2

    def E1(self, x):
            '''Integral exponential function on interval [x,inf)
            with integrand exp(-t)/t'''
            if 0<x<=1.:
                return -.57721566-log(x)+((((.00107857*x-.00976004)* \
                x+.05519968)*x-.24991055)*x+.99999193)*x
            else:
                return exp(-x)*((((x+8.5733287401)*x+18.059016973)*x+8.6347608925 \
                )*x+.2677737343)/x/((((x+9.5733223454)*x+25.6329561486)*x+ \
                21.0996530827)*x+3.9584969228)

    def f1(self, x):
        if x == 0.: return 1.
        return sinh(x)/x

    def F(self, x, c):
        '''Function to calculate singular integral which integrand is
        t^2*exp(-t*c)/(t-x) and 0<=t,n<inf'''
        if x==0. : return 1./c**2
        cx=c*x
        return (1./c+x)/c+x*x*exp(-cx)*(-2.*self.rmb(self.f1,0.,cx) + self.E1(cx))
    
    def vpI(self, r, n):
            '''Function to calculate sigular integral when representing
            added mass coefficient by using Krammers - Kronig relation'''
            eps=1.e-6
            t=r/(1.+sqrt(1.-r*r))   # t=exp(-t0)
            t0=-log(t)
            en=eps*exp(n)
            N=int(log(1./sqrt(en)+sqrt(1.+1./en))/t0)+1
            bi=0.
            emi=1.
            sgi=-1.
            for i in range(1,N):
                emi *= t
                epi=1./emi
                si=.5*(epi-emi)
                ci=.5*(epi+emi)
                cti=ci/si
                s2i=si*si
                bj=0.
                sgj=-1.
                emj=1.
                for j in range(1,N):
                    emj *= t
                    epj=1./emj
                    sj=.5*(epj-emj)
                    cj=.5*(epj+emj)
                    ctj=cj/sj
                    s2j=sj*sj
                    ctij=cti+ctj
                    uj = self.F(n,ctij)/s2j*sgj
                    bj += uj
                    sgj=-sgj
                bi += bj/s2i*sgi
                sgi=-sgi
            return pi*bi

    @property
    def fmu(self):    # Function to calculate added mass coefs at finite freq
        eps=1.e-6
        r=self.__R/self.__H
        l=self.__H*sqrt(1.-r*r)
        n=self.__sgm**2*l/g
        f2 = self.f_minf
        return f2 + self.vpI(r,n)    # fminf()

    # End of class HFrm

class HBd(HFrm):
    '''Class for apparatus body hydrodynamic coefficients to be
        calculated; parameters are:
        L,m - body length
        R,m - max body radius
        H,m - underwater depth of body frame centers
        f - function to define body radii in longitudinal direction
        sgm0,1/sec - incident waves circular friequancy
        v,m/sec - body velocity
        q,rad - angle between body velocity and incident waves velocity'''

    def __init__(self,L,R,f,H,sgm0,v=0.,q=0.):
        self.__L = L
        self.__R = R
        self.__f = f
        self.__H = H
        self.__sgm0 = sgm0
        self.__v = v
        self.__q = q
        #HFrm.__init__(R,H,sgm0)

    def __fsgm(self):    # friequancy of encounter (due to body velocity)
            
            return self.__sgm0+self.__sgm0**2/g*self.__v*cos(self.__q)

    def __fdlv(self):    # dimensionless velocity
        s=self.__fsgm()
        if s==0. : return .0
        return self.__v/s/self.__L

    def __fLmb(self,x):
        #if -.5<x or x>.5: Error
        eps=1.e-6   
        s=self.__fsgm()
        r=self.__R/self.__H*self.__f(x)
        if r<eps: r=eps
        C=HFrm(r,1.,s)
        return C.flmb

    def __fMu(self,x):
        #if -.5<x or x>.5: Error
        eps=1.e-6
        s=self.__fsgm()
        r=self.__R/self.__H*self.__f(x)
        if r<eps: r=eps
        C=HFrm(r,1.,s)
        return C.fmu

    # Calculation of added mass and damping coefficiebts

    @property
    def fLmb0(self):
            return HFrm.rmb(self.__fLmb,-.5,.5,1.0e-3)

    @property
    def fMu0(self):
            return HFrm.rmb(self.__fMu,-.5,.5,1.0e-3)

    def __fLmb1(self,x):
        #if -.5<x or x>.5: Error
        return x*self.__fLmb(x)

    def __fLmb2(self,x):
        #if -.5<x or x>.5: Error
        return x*x*self.__fLmb(x)

    def __fMu1(self,x):
        #if -.5<x or x>.5: Error
        return x*self.__fMu(x)

    def __fMu2(self,x):
        #if -.5<x or x>.5: Error
        return x*x*self.__fMu(x)

    @property
    def fLmb1(self):
            return HFrm.rmb(self.__fLmb1,-.5,.5,1.0e-3)

    @property
    def fMu1(self):
            return HFrm.rmb(self.__fMu1,-.5,.5,1.0e-3)

    @property
    def fLmb2(self):
            return HFrm.rmb(self.__fLmb2,-.5,.5,1.0e-3)

    @property
    def fMu2(self):
            return HFrm.rmb(self.__fMu2 ,-.5,.5,1.0e-3)

    # Ahead velocity of apparatus is absent (zero )

    @property
    def fl22(self):    # Apparatus damping coefficient (l22)
            return self.fLmb0

    @property
    def fm22(self):    # Apparatus added mass coefficient (m22)
            return self.fMu0

    @property
    def fl33(self):    # Apparatus damping coefficient (l33)
            return self.fLmb0

    @property
    def fm33(self):    # Apparatus added mass coefficient (m33)
            return self.fMu0

    @property
    def fl26(self):    # Apparatus damping coefficient (l26)
            return self.fLmb1

    @property
    def fm26(self):    # Apparatus added mass coefficient (m26)
            return self.fMu1

    @property
    def fl35(self):    # Apparatus damping coefficient (l35)
            return self.fLmb1

    @property
    def fm35(self):    # Apparatus added mass coefficient (m35)
            return self.fMu1

    @property
    def fl55(self):    # Apparatus damping coefficient (l55)
            return self.fLmb2

    @property
    def fm55(self):    # Apparatus added mass coefficient (m55)
            return self.fMu2+v*v*self.fMu0

    @property
    def fl66(self):    # Apparatus damping coefficient (l66)
            return self.fLmb2

    @property
    def fm66(self):    # Apparatus added mass coefficient (m66)
            return self.fMu2


    # Ahead velocity of apparatus is not zero (v != 0.)

    @property
    def fL22(self):    # Apparatus damping coefficient (L22)
            return self.fLmb0

    @property
    def fM22(self):    # Apparatus added mass coefficient (M22)
            return self.fMu0

    @property
    def fL33(self):    # Apparatus damping coefficient (L33)
            return self.fLmb0

    @property
    def fM33(self):    # Apparatus added mass coefficient (M33)
            return self.fMu0

    @property
    def fL26(self):    # Apparatus damping coefficient (L26)
            v=self.__fdlv()
            return self.fLmb1-v*self.fMu0

    @property
    def fM26(self):    # Apparatus added mass coefficient (M26)
            v=self.__fdlv()
            return self.fMu1+v*self.fLmb0

    @property
    def fL62(self):    # Apparatus damping coefficient (L62)
            v=self.__fdlv()
            return self.fLmb1+v*self.fMu0

    @property
    def fM62(self):    # Apparatus added mass coefficient (M62)
            v=self.__fdlv()
            return self.fMu1-v*self.fLmb0

    @property
    def fL35(self):    # Apparatus damping coefficient (L35)
            v=self.__fdlv()
            return self.fLmb1-v*self.fMu0

    @property
    def fM35(self):    # Apparatus added mass coefficient (M35)
            v=self.__fdlv()
            return self.fMu1+v*self.fLmb0

    @property
    def fL53(self):    # Apparatus damping coefficient (L53)
            v=self.__fdlv()
            return self.fLmb1+v*self.fMu0

    @property
    def fM53(self):    # Apparatus added mass coefficient (M53)
            v=self.__fdlv()
            return self.fMu1-v*self.fLmb0

    @property
    def fL55(self):    # Apparatus damping coefficient (L55)
            v=self.__fdlv()
            return self.fLmb2+v*v*self.fLmb0

    @property
    def fM55(self):    # Apparatus added mass coefficient (M55)
            v=self.__fdlv()
            return self.fMu2+v*v*self.fMu0

    @property
    def fL66(self):    # Apparatus damping coefficient (L66)
            v=self.__fdlv()
            return self.fLmb2+v*v*self.fLmb0

    @property
    def fM66(self):    # Apparatus added mass coefficient (M66)
            v=self.__fdlv()
            return self.fMu2+v*v*self.fMu0


    # Calculation of excited forces

    def __fBc(self,x):
        #if -.5<x or x>.5: Error
        
        eps=1.e-6
        s0=self.__sgm0 
        s1=s0*sqrt(abs(sin(self.__q)))
        r=self.__R/self.__H*self.__f(x)
        if r<eps: r=eps
        k=s0*s0/g*self.__L*cos(self.__q)
        C=HFrm(r,1.,s1)
        return C.fB*cos(k*x)

    def __fBs(self,x):
        #if -.5<x or x>.5: Error
        
        s0=self.__sgm0 
        s1=s0*sqrt(abs(sin(self.__q)))
        r=self.__R/self.__H*self.__f(x)
        k=s0*s0/g*self.__L*cos(self.__q)
        C=HFrm(r,1.,s1)
        return C.fB*sin(k*x)

    def __fBxc(self,x):
        #if -.5<x or x>.5: Error
        
        eps=1.e-6
        s0=self.__sgm0 
        s1=s0*sqrt(abs(sin(self.__q)))
        r=self.__R/self.__H*self.__f(x)
        if r<eps: r=eps
        k=s0*s0/g*self.__L*cos(self.__q)
        C=HFrm(r,1.,s1)
        return C.fB*x*cos(k*x)

    def __fBxs(self,x):
        #if -.5<x or x>.5: Error
        
        eps=1.e-6
        s0=self.__sgm0 
        s1=s0*sqrt(abs(sin(self.__q)))
        r=self.__R/self.__H*self.__f(x)
        if r<eps: r=eps
        k=s0*s0/g*self.__L*cos(self.__q)
        C=HFrm(r,1.,s1)
        return C.fB*x*sin(k*x)

    # fFc/fFs related (Ar/Ai)2,3 and fFc1/fFs1 - (Ar/Ai)5,6

    @property
    def fAr2(self):    # Real part of excited force amplitude (Ar2)
            return HFrm.rmb(self.__fBc,-.5,.5,1.0e-3)

    @property
    def fAi2(self):    # Imag. part of excited force amplitude (Ai2)
            return HFrm.rmb(self.__fBs,-.5,.5,1.0e-3)

    @property
    def fAr6(self):    # Real part of excited force amplitude (Ar6)
            return HFrm.rmb(self.__fBxc,-.5,.5,1.0e-3)

    @property
    def fAi6(self):    # Imag. part of excited force amplitude (Ai6)
            return HFrm.rmb(self.__fBxs,-.5,.5,1.0e-3)

    @property
    def fAr3(self):    # Real part of excited force amplitude (Ar3)
            return -HFrm.rmb(self.__fBs,-.5,.5,1.0e-3)

    @property
    def fAi3(self):    # Imag. part of excited force amplitude (Ai3)
            return HFrm.rmb(self.__fBc,-.5,.5,1.0e-3)

    @property
    def fAr5(self):    # Real part of excited force amplitude (Ar5)
            return -HFrm.rmb(self.__fBxs,-.5,.5,1.0e-3)

    @property
    def fAi5(self):    # Imag. part of excited force amplitude (Ai5)
            return HFrm.rmb(self.__fBxc,-.5,.5,1.0e-3)


    # End of class HBd

# 2D GRAPHS

def Tb_B(r = 8.5,h = 1.0):
    '''Table of Bj values as function of R/H (col)
    R - circle radius
    H - depth circle center under free fluid surface
    x - frequency parameter x = sigma^2*l/g, l=sqrt(H*H-R*R)


    Bj - define nondimensional amplinude of exciting fofces
    (sum of Krylov's and diffraction parts) applied to circle
    contour due to progressive waves
    '''
    lr = [.1,.2,.3,.4,.5,.6,.7,.8,.9,.95] # for testing
    dx = .05
    x = -dx
    x_l = []
    data_fb = {}
    data_lmb = {}
    data_mu = {}
    navigate_flot = []
    data_fb[r] = []
    data_mu[r] = []
    data_lmb[r] = []

    for i in range(1,128):
        x += dx

    #     for r in lr:
        l = sqrt(1.-r*r)
        s = sqrt(g*x/l)
        C = HFrm(r,h,s)
        data_fb[r].append([float('%.2f'%x),float('%.6f'%C.fB)])
        data_lmb[r].append([float('%.2f'%x),float('%.6f'%C.flmb)])
        data_mu[r].append([float('%.2f'%x),float('%.6f'%C.fmu)])

    #Convert to json for javasript graph
    if data_fb:
        for k,v in data_fb.iteritems():
            navigate_flot.append({"name":"Bj", "data":v})       

    if data_lmb:
        for k,v in data_lmb.iteritems():
            navigate_flot.append({"name":"lmd", "data":v})

    if data_mu:
        for k,v in data_mu.iteritems():
            navigate_flot.append({"name":"mu", "data":v})

    return navigate_flot

# SURFACES

def SrfPrp(frm,prp,sdtDct):
    '''Function to calculate table for surface Prp as a function of
        L/lmb - (rows) and q - (cols) to be constructed, L - apparatus
        length, lmb - incident wave length, q - encounter wave angle.
        Parameters:
        frm - 'E' or 'C' when apparatus is ellipsoid or cone+cylinder
        prp - specified apparatus property
        sdtDct - dictionary of given apparatus data with next keys
        'L' - length, m
        'R' - max hull radius, m
        'H' - depth of underwater centre line, m
        'rb' - max value of L/lmb, dementionless
        'v' - velocity, m/sec
        'dr' - increment of L/lmb, dementionless
        'qmn' - min value of angle q, deg
        'qmx' - max value of angle q, deg
        'dq' - increment of angle q, deg

        Table of mentioned elements' values are written into file
        'Srf_PrpFrm.grd', for example, Srf_M35E.grd'''

    Prps=set(('l22','l26','l66','l33','l35','l55',
    'm22','m26','m66','m33','m35','m55',
        'L22','L26','L62','L66','L33','L35','L53','L55',
    'M22','M26','M62','M66','M33','M35','M53','M55',
    'Ar2','Ai2','Ar6','Ai6','Ar3','Ai3','Ar5','Ai5'))

    eps=1.e-4


    def f1(x):
        ''' Ellipsoid '''
        
        if abs(x)<=.5: y=sqrt(1.-4.*x*x)
        if y<eps: return eps
        return y

    def f2(x):
        ''' Cone & Cylinder '''
        if -.5<=x<=1./6.: return 1.
        if 1./6.<x<=.5: y=3.*(.5-x)
        if y<eps: return eps
        return y

    if prp not in Prps:
        print "Property %s does not exist\n" % prp
        return

    L=sdtDct['L']
    R=sdtDct['R']
    H=sdtDct['H']
    v=sdtDct['v']
    xb=sdtDct['rb']
    dx=sdtDct['dr']
    qmn=sdtDct['qmn']
    qmx=sdtDct['qmx']
    dq=sdtDct['dq']
    

    rad2deg=180./pi
    deg2rad=pi/180.

    xmn=dx 
    xmx=xb
    ymn=qmn*deg2rad 
    ymx=qmx*deg2rad
    dy=dq*deg2rad
    zmn = zmx = 0

    nx = int(round((xmx-xmn)/dx))
    ny = int(round((ymx-ymn)/dy))

    outfile='Srf_'+prp+frm+'.grd'
    p='C.f'+prp

    z=ndarray(shape=(nx,ny), dtype=float, order='C')

    x=0.
    for i in range(0,nx):
        x += dx
        s=sqrt(2.*pi*g/L*x)
        y=0.
        for j in range(0,ny):
            y += dy
        C=HBd(L,R,f2,H,s,v,y)
        z[i,j]=eval(p)
        if (i+j) == 0: zmn=zmx=z[i,j]
        if zmn > z[i,j]: zmn = z[i,j]
        if zmx < z[i,j]: zmx = z[i,j]

    fo=open(outfile,'w')

    print >> fo, "DSAA"                                # info for .grd - file
    print >> fo, '%d   %d' %(ny,nx)                    # quantity of points ny, nx
    print >> fo, '%.6f   %.6f' %(qmn,qmx)              # ymin, ymax
    print >> fo, '%.6f   %.6f' %(dx,xb)                # xmin, xmax
    print >> fo, '%.6f   %.6f' %(zmn,zmx)              # zmin, zmax
    print >> fo
    for i in range(0,nx):
        for j in range(0,ny):
            print >> fo, '%12.6f ' % z[i,j],
    print >> fo

    fo.close()

    # End of function SrfPrp(frm,prp,sdtDct)

def ClcSrfPrp(frmLst,prpLst,sdtDct):
    '''Calculation and write down into files values of apparatus
    property for prp-surfaces to create 
    Parameters:
    frm - 'E' or 'C' when apparatus is ellipsoid or cone+cylinder
    prp - specified apparatus property
    sdtDct - dictionary of given apparatus source data''' 

    for frm in frmLst:
        for prp in prpLst:
            SrfPrp(frm,prp,sdtDct)
        
    # End of function ClcSrfPrp(frmLst,prpLst,sdtDct)

def CrvPrp(frm,prp,sdtDct):
    '''Function to calculate table for curve Prp as a function of
        L/lmb - (rows) and q - (cols) to be constructed, L - apparatus
        length, lmb - incident wave length, q - encounter wave angle.
        Parameters:
        frm - 'E' or 'C' when apparatus is ellipsoid or cone+cylinder
        prp - specified apparatus property
        sdtDct - dictionary of given apparatus data with next keys
        'L' - length, m
        'R' - max hull radius, m
        'H' - depth of underwater centre line, m
        'rb' - max value of L/lmb, dementionless
        'v' - velocity, m/sec
        'dr' - increment of L/lmb, dementionless
        'qmn' - min value of angle q, deg
        'qmx' - max value of angle q, deg
        'dq' - increment of angle q, deg

        Table of mentioned elements' values are written into file
        'Crv_PrpFrm.grd', for example, Crv_M35E.grd'''

    Prps=set(('l22','l26','l66','l33','l35','l55',
    'm22','m26','m66','m33','m35','m55',
        'L22','L26','L62','L66','L33','L35','L53','L55',
    'M22','M26','M62','M66','M33','M35','M53','M55',
    'Ar2','Ai2','Ar6','Ai6','Ar3','Ai3','Ar5','Ai5'))

    def f1(x):
        ''' Ellipsoid '''
        eps=1.e-4
        if abs(x)<=.5: y=sqrt(1.-4.*x*x)
        if y<eps: return eps
        return y

    def f2(x):
        ''' Cone & Cylinder '''
        if -.5<=x<=1./6.: return 1.
        eps=1.e-4
        if 1./6.<x<=.5: y=3.*(.5-x)
        if y<eps: return eps
        return y

    if prp not in Prps:
        print "Property %s does not exist\n" % prp
        return

    L=sdtDct['L']
    R=sdtDct['R']
    H=sdtDct['H']
    v=sdtDct['v']
    xb=sdtDct['rb']
    dx=sdtDct['dr']
    qmn=sdtDct['qmn']
    qmx=sdtDct['qmx']
    dq=sdtDct['dq']
    

    rad2deg=180./pi
    deg2rad=pi/180.

    xmn=dx 
    xmx=xb
    ymn=qmn*deg2rad 
    ymx=qmx*deg2rad
    dy=dq*deg2rad

    nx = int(round((xmx-xmn)/dx))
    ny = int(round((ymx-ymn)/dy))

    outfile='Crv_'+prp+frm+'.grd'
    p='C.f'+prp

    fo=open(outfile,'w')

    x=0.
    for i in range(0,nx+1):
        x += dx
        print >> fo, '%6.3f ' % x,
        s=sqrt(2.*pi*g/L*x)
        y=0.
        for j in range(0,ny):
            y += dy
        C=HBd(L,R,f2,H,s,v,y)
        print >> fo, '%12.6f ' % eval(p),
        print >> fo

    fo.close()

    # End of function CrvPrp(frm,prp,sdtDct)

def ClcCrvPrp(frmLst,prpLst,sdtDct):
    '''Calculation and write down into files values of apparatus
    property for prp-curves to create
    Parameters:
    frm - 'E' or 'C' when apparatus is ellipsoid or cone+cylinder
    prp - specified apparatus property
    sdtDct - dictionary of given apparatus source data''' 

    for frm in frmLst:
        for prp in prpLst:
            CrvPrp(frm,prp,sdtDct)
        
    # End of function ClcCrvPrp(frmLst,prpLst,sdtDct)

def SrfSgm1(sdtDct):
    '''Function to calculate dementionless friequancy of encounter
        (due to body velocity) table for surface Sgm1 as a function of 
        L/lmb - (rows) and q - (cols) to be constructed, L - apparatus
        length, lmb - incident wave length, q - encounter wave angle.
        Parameters:
        sdtDct - dictionary of given apparatus data with next keys
        'L' - length, m
        'R' - max hull radius, m
        'H' - depth of underwater centre line, m
        'rb' - max value of L/lmb, dementionless
        'v' - velocity, m/sec
        'dr' - increment of L/lmb, dementionless
        'qmn' - min value of angle q, deg
        'qmx' - max value of angle q, deg
        'dq' - increment of angle q, deg

        Table of mentioned dementionless friequancy of encounter values 
        (sgm*sqrt(L/g)) are written into file Srf_Sgm1.grd'''

    L=sdtDct['L']
    v=sdtDct['v']
    xb=sdtDct['rb']
    dx=sdtDct['dr']
    qmn=sdtDct['qmn']
    qmx=sdtDct['qmx']
    dq=sdtDct['dq']
    
    Fr=v/sqrt(g*L)  # Froud's number

    rad2deg=180./pi
    deg2rad=pi/180.

    xmn=dx 
    xmx=xb
    ymn=qmn*deg2rad 
    ymx=qmx*deg2rad
    dy=dq*deg2rad

    nx = int(round((xmx-xmn)/dx))
    ny = int(round((ymx-ymn)/dy))

    outfile='Srf_Sgm1.grd'

    z=ndarray(shape=(nx,ny), dtype=float, order='C')

    x=0.
    for i in range(0,nx):
        x += dx
        s=sqrt(2.*pi*x)    # Circular wave friequancy 
        y=0.
        for j in range(0,ny):
            y += dy
        s1=s*(1.+s*cos(y)*Fr)    # Friequancy of encounter 
        z[i,j]=s1
        if (i+j)==0: zmn=zmx=z[i,j]
        if zmn > z[i,j]: zmn = z[i,j]
        if zmx < z[i,j]: zmx = z[i,j]

    fo=open(outfile,'w') 

    print >> fo, "DSAA"                                # info for .grd - file
    print >> fo, '%d   %d' %(ny,nx)                    # quantity of points ny, nx
    print >> fo, '%.6f   %.6f' %(qmn,qmx)              # ymin, ymax
    print >> fo, '%.6f   %.6f' %(dx,xb)                # xmin, xmax
    print >> fo, '%.6f   %.6f' %(zmn,zmx)              # zmin, zmax
    print >> fo
    for i in range(0,nx):
        for j in range(0,ny):
            print >> fo, '%12.6f ' % z[i,j],
    print >> fo

    fo.close()

    # End of function SrfSgm1(sdtDct)


if __name__ == '__main__':
    sdtDct={'L':15.,'R':.75,'H':2.5,'v':5.,'rb':3.5,'dr':.05,'qmn':20.,'qmx':160.,'dq':5.}
    frmLst=['C','E']
    prpLst=['l22','l26','l66','m22','m26','m66','Ar2','Ai2','Ar6','Ai6']
    
    # graph = Tb_B()

    ClcSrfPrp(frmLst,prpLst,sdtDct)
    ClcCrvPrp(frmLst,prpLst,sdtDct)
