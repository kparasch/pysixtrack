import TricubicInterpolation.cTricubic as ti
import numpy as np
import h5py
import myfilemanager as mfm_e

class ElectronCloud(object):

    def __init__(self,fname, scale=1., method='Exact'):
        ob = h5py.File(fname)
        xg = ob['xg'][()]
        yg = ob['yg'][()]
        zg = ob['zg'][()]
        phi = scale*ob['phi'][()]
        
        #ob = mfm_e.myloadmat(fname)
        #xg = ob['xg']
        #yg = ob['yg']
        #zg = ob['zg']
        #phi = scale*ob['phi'].transpose(1,2,0)
        
        print('ecloud:')
        print('\t xg: %f'%xg[-1])
        print('\t yg: %f'%yg[-1])
        print('\t zg: %f'%zg[-1])
        print('\t dx: %f'%(xg[1]-xg[0]))
        print('\t dy: %f'%(yg[1]-yg[0]))
        print('\t dz: %f'%(zg[1]-zg[0]))
        self.TI = ti.Tricubic_Interpolation(A=phi, x0=xg[0], y0=yg[0], z0=zg[0],
                                       dx=xg[1]-xg[0], dy=yg[1]-yg[0], dz=zg[1]-zg[0], method = method)

        #dipolar kicks
        xkick, ykick, zkick = self.TI.kick(0.,0.,0.)
        self.px0 = xkick
        self.py0 = ykick
        self.delta0 = zkick

    def vkick(self, p):
        if len(p.x) == 1:
            xkick, ykick, zkick = self.TI.kick(p.x, p.y, p.zeta)
        else:
            xkick = np.empty([len(p.x)])
            ykick = np.empty_like(xkick)
            zkick = np.empty_like(xkick)
            for i in range(len(p.x)):
                if self.TI.is_inside_box(p.x[i], p.y[i], p.zeta[i]):
                    xkick[i], ykick[i], zkick[i] = self.TI.kick(p.x[i], p.y[i], p.zeta[i])
                else:
                    xkick[i] = float('nan')
                    ykick[i] = float('nan')
                    zkick[i] = float('nan')
        return xkick, ykick, zkick

    def track(self,p):
        xkick, ykick, zkick = self.vkick(p)
        p.px += xkick - self.px0
        p.py += ykick - self.py0
#        p.delta += zkick - self.delta0
