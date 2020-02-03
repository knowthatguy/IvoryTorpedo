# -*- coding: utf-8 -*-
from hydrodynamic_coefficients import *

g = 9.81


class HGraph(object):
    # 2D GRAPHS
    def spline(self, r=0.8, h=1.0):
        """Table of Bj values as function of R/H (col)
        R - circle radius
        H - depth circle center under free fluid surface
        x - frequency parameter x = sigma^2*l/g, l=sqrt(H*H-R*R)


        Bj - define nondimensional amplinude of exciting forces
        (sum of Krylov's and diffraction parts) applied to circle
        contour due to progressive waves
        """
        # lr = [.1,.2,.3,.4,.5,.6,.7,.8,.9,.95] # for testing
        dx = 0.05
        x = -dx
        x_l = []
        data_fb = {}
        data_lmb = {}
        data_mu = {}
        plot_data = []
        data_fb[r] = []
        data_mu[r] = []
        data_lmb[r] = []

        for i in range(1, 128):
            x += dx

            #     for r in lr:
            l = sqrt(1.0 - r * r)
            s = sqrt(g * x / l)
            C = HFrm(r, h, s)
            data_fb[r].append([float("%.2f" % x), float("%.6f" % C.fB)])
            data_lmb[r].append([float("%.2f" % x), float("%.6f" % C.flmb)])
            data_mu[r].append([float("%.2f" % x), float("%.6f" % C.fmu)])

        # Convert to json for javasript graph
        if data_fb:
            for k, v in data_fb.iteritems():
                plot_data.append({"name": "Bj", "data": v})

        if data_lmb:
            for k, v in data_lmb.iteritems():
                plot_data.append({"name": "lmd", "data": v})

        if data_mu:
            for k, v in data_mu.iteritems():
                plot_data.append({"name": "mu", "data": v})

        return plot_data

    # SURFACES

    def SrfPrp(self, frm, prp, sdtDct):
        """Function to calculate table for surface Prp as a function of
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
            'Srf_PrpFrm.grd', for example, Srf_M35E.grd"""

        Prps = set(
            (
                "l22",
                "l26",
                "l66",
                "l33",
                "l35",
                "l55",
                "m22",
                "m26",
                "m66",
                "m33",
                "m35",
                "m55",
                "L22",
                "L26",
                "L62",
                "L66",
                "L33",
                "L35",
                "L53",
                "L55",
                "M22",
                "M26",
                "M62",
                "M66",
                "M33",
                "M35",
                "M53",
                "M55",
                "Ar2",
                "Ai2",
                "Ar6",
                "Ai6",
                "Ar3",
                "Ai3",
                "Ar5",
                "Ai5",
            )
        )

        eps = 1.0e-4

        def f1(x):
            """ Ellipsoid """

            if abs(x) <= 0.5:
                y = sqrt(1.0 - 4.0 * x * x)
            if y < eps:
                return eps
            return y

        def f2(x):
            """ Cone & Cylinder """
            if -0.5 <= x <= 1.0 / 6.0:
                return 1.0
            if 1.0 / 6.0 < x <= 0.5:
                y = 3.0 * (0.5 - x)
            if y < eps:
                return eps
            return y

        if prp not in Prps:
            print ("Property {prp} does not exist\n")
            return

        L = sdtDct["L"]
        R = sdtDct["R"]
        H = sdtDct["H"]
        v = sdtDct["v"]
        xb = sdtDct["rb"]
        dx = sdtDct["dr"]
        qmn = sdtDct["qmn"]
        qmx = sdtDct["qmx"]
        dq = sdtDct["dq"]

        rad2deg = 180.0 / pi
        deg2rad = pi / 180.0

        xmn = dx
        xmx = xb
        ymn = qmn * deg2rad
        ymx = qmx * deg2rad
        dy = dq * deg2rad
        zmn = zmx = 0

        nx = int(round((xmx - xmn) / dx))
        ny = int(round((ymx - ymn) / dy))

        outfile = "Srf_" + prp + frm + ".grd"
        p = "C.f" + prp

        z = ndarray(shape=(nx, ny), dtype=float, order="C")

        x = 0.0
        for i in range(0, nx):
            x += dx
            s = sqrt(2.0 * pi * g / L * x)
            y = 0.0
            for j in range(0, ny):
                y += dy
            C = HBd(L, R, f2, H, s, v, y)
            z[i, j] = eval(p)
            if (i + j) == 0:
                zmn = zmx = z[i, j]
            if zmn > z[i, j]:
                zmn = z[i, j]
            if zmx < z[i, j]:
                zmx = z[i, j]

        fo = open(outfile, "w")

        print >> fo, "DSAA"  # info for .grd - file
        print >> fo, "%d   %d" % (ny, nx)  # quantity of points ny, nx
        print >> fo, "%.6f   %.6f" % (qmn, qmx)  # ymin, ymax
        print >> fo, "%.6f   %.6f" % (dx, xb)  # xmin, xmax
        print >> fo, "%.6f   %.6f" % (zmn, zmx)  # zmin, zmax
        print >> fo
        for i in range(0, nx):
            for j in range(0, ny):
                print >> fo, "%12.6f " % z[i, j],
        print >> fo

        fo.close()

        # End of function SrfPrp(frm,prp,sdtDct)

    def ClcSrfPrp(self, frmLst, prpLst, sdtDct):
        """Calculation and write down into files values of apparatus
        property for prp-surfaces to create
        Parameters:
        frm - 'E' or 'C' when apparatus is ellipsoid or cone+cylinder
        prp - specified apparatus property
        sdtDct - dictionary of given apparatus source data"""

        for frm in frmLst:
            for prp in prpLst:
                self.SrfPrp(frm, prp, sdtDct)

        # End of function ClcSrfPrp(frmLst,prpLst,sdtDct)

    def CrvPrp(self, frm, prp, sdtDct):
        """Function to calculate table for curve Prp as a function of
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
            'Crv_PrpFrm.grd', for example, Crv_M35E.grd"""

        Prps = set(
            (
                "l22",
                "l26",
                "l66",
                "l33",
                "l35",
                "l55",
                "m22",
                "m26",
                "m66",
                "m33",
                "m35",
                "m55",
                "L22",
                "L26",
                "L62",
                "L66",
                "L33",
                "L35",
                "L53",
                "L55",
                "M22",
                "M26",
                "M62",
                "M66",
                "M33",
                "M35",
                "M53",
                "M55",
                "Ar2",
                "Ai2",
                "Ar6",
                "Ai6",
                "Ar3",
                "Ai3",
                "Ar5",
                "Ai5",
            )
        )

        def f1(x):
            """ Ellipsoid """
            eps = 1.0e-4
            if abs(x) <= 0.5:
                y = sqrt(1.0 - 4.0 * x * x)
            if y < eps:
                return eps
            return y

        def f2(x):
            """ Cone & Cylinder """
            if -0.5 <= x <= 1.0 / 6.0:
                return 1.0
            eps = 1.0e-4
            if 1.0 / 6.0 < x <= 0.5:
                y = 3.0 * (0.5 - x)
            if y < eps:
                return eps
            return y

        if prp not in Prps:
            print "Property %s does not exist\n" % prp
            return

        L = sdtDct["L"]
        R = sdtDct["R"]
        H = sdtDct["H"]
        v = sdtDct["v"]
        xb = sdtDct["rb"]
        dx = sdtDct["dr"]
        qmn = sdtDct["qmn"]
        qmx = sdtDct["qmx"]
        dq = sdtDct["dq"]

        rad2deg = 180.0 / pi
        deg2rad = pi / 180.0

        xmn = dx
        xmx = xb
        ymn = qmn * deg2rad
        ymx = qmx * deg2rad
        dy = dq * deg2rad

        nx = int(round((xmx - xmn) / dx))
        ny = int(round((ymx - ymn) / dy))

        outfile = "Crv_" + prp + frm + ".grd"
        p = "C.f" + prp

        fo = open(outfile, "w")

        x = 0.0
        for i in range(0, nx + 1):
            x += dx
            print >> fo, "%6.3f " % x,
            s = sqrt(2.0 * pi * g / L * x)
            y = 0.0
            for j in range(0, ny):
                y += dy
            C = HBd(L, R, f2, H, s, v, y)
            print >> fo, "%12.6f " % eval(p),
            print >> fo

        fo.close()

        # End of function CrvPrp(frm,prp,sdtDct)

    def ClcCrvPrp(self, frmLst, prpLst, sdtDct):
        """Calculation and write down into files values of apparatus
        property for prp-curves to create
        Parameters:
        frm - 'E' or 'C' when apparatus is ellipsoid or cone+cylinder
        prp - specified apparatus property
        sdtDct - dictionary of given apparatus source data"""

        for frm in frmLst:
            for prp in prpLst:
                self.CrvPrp(frm, prp, sdtDct)

        # End of function ClcCrvPrp(frmLst,prpLst,sdtDct)

    def SrfSgm1(self, sdtDct):
        """Function to calculate dementionless friequancy of encounter
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
            (sgm*sqrt(L/g)) are written into file Srf_Sgm1.grd"""

        L = sdtDct["L"]
        v = sdtDct["v"]
        xb = sdtDct["rb"]
        dx = sdtDct["dr"]
        qmn = sdtDct["qmn"]
        qmx = sdtDct["qmx"]
        dq = sdtDct["dq"]

        Fr = v / sqrt(g * L)  # Froud's number

        rad2deg = 180.0 / pi
        deg2rad = pi / 180.0

        xmn = dx
        xmx = xb
        ymn = qmn * deg2rad
        ymx = qmx * deg2rad
        dy = dq * deg2rad

        nx = int(round((xmx - xmn) / dx))
        ny = int(round((ymx - ymn) / dy))

        outfile = "Srf_Sgm1.grd"

        z = ndarray(shape=(nx, ny), dtype=float, order="C")

        x = 0.0
        for i in range(0, nx):
            x += dx
            s = sqrt(2.0 * pi * x)  # Circular wave friequancy
            y = 0.0
            for j in range(0, ny):
                y += dy
            s1 = s * (1.0 + s * cos(y) * Fr)  # Friequancy of encounter
            z[i, j] = s1
            if (i + j) == 0:
                zmn = zmx = z[i, j]
            if zmn > z[i, j]:
                zmn = z[i, j]
            if zmx < z[i, j]:
                zmx = z[i, j]

        fo = open(outfile, "w")

        print >> fo, "DSAA"  # info for .grd - file
        print >> fo, "%d   %d" % (ny, nx)  # quantity of points ny, nx
        print >> fo, "%.6f   %.6f" % (qmn, qmx)  # ymin, ymax
        print >> fo, "%.6f   %.6f" % (dx, xb)  # xmin, xmax
        print >> fo, "%.6f   %.6f" % (zmn, zmx)  # zmin, zmax
        print >> fo
        for i in range(0, nx):
            for j in range(0, ny):
                print >> fo, "%12.6f " % z[i, j],
        print >> fo

        fo.close()

        # End of function SrfSgm1(sdtDct)


if __name__ == "__main__":
    sdtDct = {
        "L": 15.0,
        "R": 0.75,
        "H": 2.5,
        "v": 5.0,
        "rb": 3.5,
        "dr": 0.05,
        "qmn": 20.0,
        "qmx": 160.0,
        "dq": 5.0,
    }
    frmLst = ["C", "E"]
    prpLst = ["l22", "l26", "l66", "m22", "m26", "m66", "Ar2", "Ai2", "Ar6", "Ai6"]

    graph = HGraph()

    graph.ClcSrfPrp(frmLst, prpLst, sdtDct)
    graph.ClcCrvPrp(frmLst, prpLst, sdtDct)
