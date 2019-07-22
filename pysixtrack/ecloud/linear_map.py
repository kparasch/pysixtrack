import numpy as np
from pysixtrack.base_classes import Element

class LinearMap4D(Element):
    """4D Linear Map from tunes and Optics functions"""
    _description = [
        ('mux', 'rad', "Phase advance in the horizontal plane", 0),
        ('muy', 'rad', "Phase advance in the vertical plane", 0),
        ('beta_x_1', 'm', "Beta function of horizontal plane at starting position", 0),
        ('beta_y_1', 'm', "Beta function of vertical plane at starting position", 0),
        ('alfa_x_1', 'm', "Alpha function of horizontal plane at starting position", 0),
        ('alfa_y_1', 'm', "Alpha function of vertical plane at starting position", 0),
        ('beta_x_2', 'm', "Beta function of horizontal plane at final position", 0),
        ('beta_y_2', 'm', "Beta function of vertical plane at final position", 0),
        ('alfa_x_2', 'm', "Alpha function of horizontal plane at final position", 0),
        ('alfa_y_2', 'm', "Alpha function of vertical plane at final position", 0)
#        ('Disp_x_1', 'm', "Horizontal dispersion at starting position", 0),
#        ('Disp_y_1', 'm', "Vertical dispersion at starting position", 0),
#        ('Disp_x_2', 'm', "Horizontal dispersion at final position", 0),
#        ('Disp_y_2', 'm', "Vertical dispersion at final position", 0)
    ]


    def track(self, p):
        cx = p._m.cos(self.mux)
        sx = p._m.sin(self.mux)
        cy = p._m.cos(self.muy)
        sy = p._m.sin(self.muy)

        C00 = p._m.sqrt(self.beta_x_2 / self.beta_x_1)
        C01 = 0.
        C10 = (p._m.sqrt(1. / (self.beta_x_1 * self.beta_x_2)) * (self.alfa_x_1 - self.alfa_x_2))
        C11 = p._m.sqrt(self.beta_x_1 / self.beta_x_2)

        C22 = p._m.sqrt(self.beta_y_2 / self.beta_y_1)
        C23 = 0.
        C32 = (p._m.sqrt(1. / (self.beta_x_1 * self.beta_x_2)) * (self.alfa_x_1 - self.alfa_x_2))
        C33 = p._m.sqrt(self.beta_x_1 / self.beta_x_2)

        S00 = p._m.sqrt(self.beta_x_2 / self.beta_x_1) * self.alpha_x_1
        S01 = p._m.sqrt(self.beta_x_1 * self.beta_x_2)
        S10 = -(p._m.sqrt(1. / (self.beta_x_1 * self.beta_x_2)) * (1. + self.alfa_x_1 * self.alfa_x_2))
        S11 = -(p._m.sqrt(self.beta_x_1 / self.beta_x_2) * self.alfa_x_2)

        S22 = p._m.sqrt(self.beta_y_2 / self.beta_y_1) * self.alfa_y_1
        S23 = p._m.sqrt(self.beta_y_1 * self.beta_y_2)
        S32 = -(p._m.sqrt(1. / (self.beta_y_1 * self.beta_y_2)) * (1. + self.alfa_y_1 * self.alfa_y_2))
        S33 = -(p._m.sqrt(self.beta_y_1 / self.beta_y_2) * self.alfa_y_2)

        M00 = C00 * cx + S00 * sx
        M10 = C10 * cx + S10 * sx
        M01 = C01 * cx + S01 * sx
        M11 = C11 * cx + S11 * sx

        M22 = C22 * cy + S22 * sy
        M32 = C32 * cy + S32 * sy
        M23 = C23 * cy + S23 * sy
        M33 = C33 * cy + S33 * sy

        xp = p.px * p.rpp
        yp = p.py * p.rpp

        x_new = M00 * p.x + M01 * xp
        yp_new = M10 * p.x + M11 * xp
        y_new = M22 * p.y + M23 * yp
        yp_new = M32 * p.y + M33 * yp

        p.x = x_new
        p.px = xp_new / p.rpp
        p.y = y_new
        p.py = yp_new / p.rpp


