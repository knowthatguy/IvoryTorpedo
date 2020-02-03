# -*- coding: utf-8 -*-
import math

from numpy import *

g = 9.81
eps = 1.0e-6


class HFrm(object):
    """
    Class for apparatus frame hydrodynamic coefficients to be
    calculated; parameters are:
    R,m - frame radius
    H,m - underwater depth of frame center
    sgm,1/sec - circular friequancy
    """

    def __init__(self, R, H, sgm):
        self.__R = R
        self.__H = H
        self.__sgm = sgm

        self.r = self.__R / self.__H
        if self.r < eps:
            self.r = eps
        self.l = self.__H * sqrt(1.0 - self.r * self.r)
        self.n = self.__sgm ** 2 * self.l / g
        self.t = self.r / (1.0 + sqrt(1.0 - self.r * self.r))
        self.t0 = -log(self.t)
        self.en = eps * exp(self.n)

    @staticmethod
    def rmb(f, a, b, tol=1.0e-6):
        """
        I = romberg(f,a,b,tol=1.0e-6).
        Romberg intergration of f(x) from x = a to b.
        Returns the integral.
        """

        def trapezoid(f, a, b, Iold, k):
            """
            Inew = trapezoid(f,a,b,Iold,k).
            Recursive trapezoidal rule:
            Iold = Integral of f(x) from x = a to b computed by
            trapezoidal rule with 2?(k-1) panels.
            Inew = Same integral computed with 2?k panels.
            """
            if k == 1:
                Inew = (f(a) + f(b)) * (b - a) / 2.0
            else:
                n = 2 ** (k - 2)  # Number of new points
                h = (b - a) / n  # Spacing of new points
                x = a + h / 2.0
                sum = 0.0
                for i in range(n):
                    sum = sum + f(x)
                    x = x + h
                    Inew = (Iold + h * sum) / 2.0
            return Inew

        def richardson(r, k):
            for j in range(k - 1, 0, -1):
                const = 4.0 ** (k - j)
                r[j] = (const * r[j + 1] - r[j]) / (const - 1.0)
            return r

        r = zeros(21)
        r[1] = trapezoid(f, a, b, 0.0, 1)
        r_old = r[1]
        for k in range(2, 21):
            r[k] = trapezoid(f, a, b, r[k - 1], k)
            r = richardson(r, k)
            if abs(r[1] - r_old) < tol * max(abs(r[1]), 1.0):
                return r[1]
            r_old = r[1]
        print("Romberg quadrature did not converge")

    @property
    def fB(self):
        """
        Function to calculate Bj - const
        """
        if self.n == 0.0:
            return 0.0

        N = int(log(1.0 / sqrt(self.en) + sqrt(1.0 + 1.0 / self.en)) / self.t0) + 1
        b = 0.0
        em = sg = 1.0
        for i in xrange(1, N):
            em *= self.t
            ep = 1.0 / em
            s = 0.5 * (ep - em)
            c = 0.5 * (ep + em)
            s2 = s * s
            u = exp(-self.n * c / s) / s2 * sg
            b += u
            sg = -sg
        return pi * self.n * b

    @property
    def flmb(self):
        """
        Function to calculate damping coefficients
        """
        return sqrt(self.n) * self.fB ** 2

    @property
    def fminf(self):
        """
        Function to calculate added mass coefs at infinite freq
        """
        N = int(round(log(1.0 / sqrt(eps) + sqrt(1.0 + 1.0 / eps)) / self.t0)) + 1
        b = 0.0
        em = sg = 1.0
        for i in range(1, N):
            em *= self.t
            ep = 1.0 / em
            s = 0.5 * (ep - em)
            s2 = s * s
            u = 1.0 / s2 * sg
            if i != 1:
                u *= 2.0
            b += u
            sg = -sg
        M = pi * b
        return M / pi * sinh(self.t0) ** 2

    def E1(self, x):
        """
        Integral exponential function on interval [x,inf)
        with integrand exp(-t)/t
        """
        if 0 < x <= 1.0:
            return (
                -0.57721566
                - log(x)
                + (
                    (((0.00107857 * x - 0.00976004) * x + 0.05519968) * x - 0.24991055)
                    * x
                    + 0.99999193
                )
                * x
            )
        else:
            return (
                exp(-x)
                * (
                    (((x + 8.5733287401) * x + 18.059016973) * x + 8.6347608925) * x
                    + 0.2677737343
                )
                / x
                / (
                    (((x + 9.5733223454) * x + 25.6329561486) * x + 21.0996530827) * x
                    + 3.9584969228
                )
            )

    def f1(self, a):
        return math.sinh(a) / a if a > 0 else 1

    def F(self, x, c):
        """
        Function to calculate singular integral which integrand is
        t^2*exp(-t*c)/(t-x) and 0<=t,n<inf
        """
        if x == 0.0:
            return 1.0 / c ** 2
        cx = c * x
        return (1.0 / c + x) / c + x * x * exp(-cx) * (
            -2.0 * self.rmb(self.f1, 0.0, cx) + self.E1(cx)
        )

    def vpI(self, r, n):
        """
        Function to calculate sigular integral when representing
        added mass coefficient by using Krammers - Kronig relation
        """
        N = int(log(1.0 / sqrt(self.en) + sqrt(1.0 + 1.0 / self.en)) / self.t0) + 1
        bi = 0.0
        emi = 1.0
        sgi = -1.0
        for i in range(1, N):
            emi *= self.t
            epi = 1.0 / emi
            si = 0.5 * (epi - emi)
            ci = 0.5 * (epi + emi)
            cti = ci / si
            s2i = si * si
            bj = 0.0
            sgj = -1.0
            emj = 1.0
            for j in range(1, N):
                emj *= self.t
                epj = 1.0 / emj
                sj = 0.5 * (epj - emj)
                cj = 0.5 * (epj + emj)
                ctj = cj / sj
                s2j = sj * sj
                ctij = cti + ctj
                uj = self.F(self.n, ctij) / s2j * sgj
                bj += uj
                sgj = -sgj
            bi += bj / s2i * sgi
            sgi = -sgi
        return pi * bi

    @property
    def fmu(self):
        """
        Function to calculate added mass coefs at finite freq
        """
        return self.fminf + self.vpI(self.r, self.n)


class HBd(HFrm):
    """
    Class for apparatus body hydrodynamic coefficients to be
    calculated; parameters are:
    L,m - body length
    R,m - max body radius
    H,m - underwater depth of body frame centers
    f - function to define body radii in longitudinal direction
    sgm0,1/sec - incident waves circular friequancy
    v,m/sec - body velocity
    q,rad - angle between body velocity and incident waves velocity
    """

    def __init__(self, L, R, f, H, sgm0, v=0.0, q=0.0):
        self.__L = L
        self.__R = R
        self.__f = f
        self.__H = H
        self.__sgm0 = sgm0
        self.__v = v
        self.__q = q

    def __fsgm(self):
        """
        Friequency of encounter (due to body velocity)
        """
        return self.__sgm0 + self.__sgm0 ** 2 / g * self.__v * cos(self.__q)

    def __fdlv(self):
        """
        Dimensionless velocity
        """
        s = self.__fsgm()
        if s == 0.0:
            return 0.0
        return self.__v / s / self.__L

    def __fLmb(self, x):
        s = self.__fsgm()
        r = self.__R / self.__H * self.__f(x)
        if r < eps:
            r = eps
        C = HFrm(r, 1.0, s)
        return C.flmb

    def __fMu(self, x):
        s = self.__fsgm()
        r = self.__R / self.__H * self.__f(x)
        if r < eps:
            r = eps
        C = HFrm(r, 1.0, s)
        return C.fmu

    # Calculation of added mass and dumping coefficiebts

    @property
    def fLmb0(self):
        return HFrm.rmb(self.__fLmb, -0.5, 0.5, 1.0e-3)

    @property
    def fMu0(self):
        return HFrm.rmb(self.__fMu, -0.5, 0.5, 1.0e-3)

    def __fLmb1(self, x):
        # if -.5<x or x>.5: Error
        return x * self.__fLmb(x)

    def __fLmb2(self, x):
        # if -.5<x or x>.5: Error
        return x * x * self.__fLmb(x)

    def __fMu1(self, x):
        # if -.5<x or x>.5: Error
        return x * self.__fMu(x)

    def __fMu2(self, x):
        # if -.5<x or x>.5: Error
        return x * x * self.__fMu(x)

    @property
    def fLmb1(self):
        return HFrm.rmb(self.__fLmb1, -0.5, 0.5, 1.0e-3)

    @property
    def fMu1(self):
        return HFrm.rmb(self.__fMu1, -0.5, 0.5, 1.0e-3)

    @property
    def fLmb2(self):
        return HFrm.rmb(self.__fLmb2, -0.5, 0.5, 1.0e-3)

    @property
    def fMu2(self):
        return HFrm.rmb(self.__fMu2, -0.5, 0.5, 1.0e-3)

    # Ahead velocity of apparatus is absent (zero )

    @property
    def fl22(self):  # Apparatus damping coefficient (l22)
        return self.fLmb0

    @property
    def fm22(self):  # Apparatus added mass coefficient (m22)
        return self.fMu0

    @property
    def fl33(self):  # Apparatus damping coefficient (l33)
        return self.fLmb0

    @property
    def fm33(self):  # Apparatus added mass coefficient (m33)
        return self.fMu0

    @property
    def fl26(self):  # Apparatus damping coefficient (l26)
        return self.fLmb1

    @property
    def fm26(self):  # Apparatus added mass coefficient (m26)
        return self.fMu1

    @property
    def fl35(self):  # Apparatus damping coefficient (l35)
        return self.fLmb1

    @property
    def fm35(self):  # Apparatus added mass coefficient (m35)
        return self.fMu1

    @property
    def fl55(self):  # Apparatus damping coefficient (l55)
        return self.fLmb2

    @property
    def fm55(self):  # Apparatus added mass coefficient (m55)
        return self.fMu2 + v * v * self.fMu0

    @property
    def fl66(self):  # Apparatus damping coefficient (l66)
        return self.fLmb2

    @property
    def fm66(self):  # Apparatus added mass coefficient (m66)
        return self.fMu2

    # Ahead velocity of apparatus is not zero (v != 0.)

    @property
    def fL22(self):  # Apparatus damping coefficient (L22)
        return self.fLmb0

    @property
    def fM22(self):  # Apparatus added mass coefficient (M22)
        return self.fMu0

    @property
    def fL33(self):  # Apparatus damping coefficient (L33)
        return self.fLmb0

    @property
    def fM33(self):  # Apparatus added mass coefficient (M33)
        return self.fMu0

    @property
    def fL26(self):  # Apparatus damping coefficient (L26)
        v = self.__fdlv()
        return self.fLmb1 - v * self.fMu0

    @property
    def fM26(self):  # Apparatus added mass coefficient (M26)
        v = self.__fdlv()
        return self.fMu1 + v * self.fLmb0

    @property
    def fL62(self):  # Apparatus damping coefficient (L62)
        v = self.__fdlv()
        return self.fLmb1 + v * self.fMu0

    @property
    def fM62(self):  # Apparatus added mass coefficient (M62)
        v = self.__fdlv()
        return self.fMu1 - v * self.fLmb0

    @property
    def fL35(self):  # Apparatus damping coefficient (L35)
        v = self.__fdlv()
        return self.fLmb1 - v * self.fMu0

    @property
    def fM35(self):  # Apparatus added mass coefficient (M35)
        v = self.__fdlv()
        return self.fMu1 + v * self.fLmb0

    @property
    def fL53(self):  # Apparatus damping coefficient (L53)
        v = self.__fdlv()
        return self.fLmb1 + v * self.fMu0

    @property
    def fM53(self):  # Apparatus added mass coefficient (M53)
        v = self.__fdlv()
        return self.fMu1 - v * self.fLmb0

    @property
    def fL55(self):  # Apparatus damping coefficient (L55)
        v = self.__fdlv()
        return self.fLmb2 + v * v * self.fLmb0

    @property
    def fM55(self):  # Apparatus added mass coefficient (M55)
        v = self.__fdlv()
        return self.fMu2 + v * v * self.fMu0

    @property
    def fL66(self):  # Apparatus damping coefficient (L66)
        v = self.__fdlv()
        return self.fLmb2 + v * v * self.fLmb0

    @property
    def fM66(self):  # Apparatus added mass coefficient (M66)
        v = self.__fdlv()
        return self.fMu2 + v * v * self.fMu0

    # Calculation of excited forces

    def __fBc(self, x):
        s0 = self.__sgm0
        s1 = s0 * sqrt(abs(sin(self.__q)))
        r = self.__R / self.__H * self.__f(x)
        if r < eps:
            r = eps
        k = s0 * s0 / g * self.__L * cos(self.__q)
        C = HFrm(r, 1.0, s1)
        return C.fB * cos(k * x)

    def __fBs(self, x):
        # if -.5<x or x>.5: Error

        s0 = self.__sgm0
        s1 = s0 * sqrt(abs(sin(self.__q)))
        r = self.__R / self.__H * self.__f(x)
        k = s0 * s0 / g * self.__L * cos(self.__q)
        C = HFrm(r, 1.0, s1)
        return C.fB * sin(k * x)

    def __fBxc(self, x):
        s0 = self.__sgm0
        s1 = s0 * sqrt(abs(sin(self.__q)))
        r = self.__R / self.__H * self.__f(x)
        if r < eps:
            r = eps
        k = s0 * s0 / g * self.__L * cos(self.__q)
        C = HFrm(r, 1.0, s1)
        return C.fB * x * cos(k * x)

    def __fBxs(self, x):
        s0 = self.__sgm0
        s1 = s0 * sqrt(abs(sin(self.__q)))
        r = self.__R / self.__H * self.__f(x)
        if r < eps:
            r = eps
        k = s0 * s0 / g * self.__L * cos(self.__q)
        C = HFrm(r, 1.0, s1)
        return C.fB * x * sin(k * x)

    # fFc/fFs related (Ar/Ai)2,3 and fFc1/fFs1 - (Ar/Ai)5,6

    @property
    def fAr2(self):  # Real part of excited force amplitude (Ar2)
        return HFrm.rmb(self.__fBc, -0.5, 0.5, 1.0e-3)

    @property
    def fAi2(self):  # Imag. part of excited force amplitude (Ai2)
        return HFrm.rmb(self.__fBs, -0.5, 0.5, 1.0e-3)

    @property
    def fAr6(self):  # Real part of excited force amplitude (Ar6)
        return HFrm.rmb(self.__fBxc, -0.5, 0.5, 1.0e-3)

    @property
    def fAi6(self):  # Imag. part of excited force amplitude (Ai6)
        return HFrm.rmb(self.__fBxs, -0.5, 0.5, 1.0e-3)

    @property
    def fAr3(self):  # Real part of excited force amplitude (Ar3)
        return -HFrm.rmb(self.__fBs, -0.5, 0.5, 1.0e-3)

    @property
    def fAi3(self):  # Imag. part of excited force amplitude (Ai3)
        return HFrm.rmb(self.__fBc, -0.5, 0.5, 1.0e-3)

    @property
    def fAr5(self):  # Real part of excited force amplitude (Ar5)
        return -HFrm.rmb(self.__fBxs, -0.5, 0.5, 1.0e-3)

    @property
    def fAi5(self):  # Imag. part of excited force amplitude (Ai5)
        return HFrm.rmb(self.__fBxc, -0.5, 0.5, 1.0e-3)

    # End of class HBd
