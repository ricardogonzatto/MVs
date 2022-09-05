import numpy as np
import healpy as hp

for i in range (200):
    cls = np.transpose(np.genfromtxt('/home/ricardo/Documentos/EstudosMestrado/Programas/EstudoMest7/datas/CLS.dat'))
    alm = hp.synalm(cls[1])
    np.savetxt(f'/home/ricardo/Documentos/EstudosMestrado/Programas/EstudoMest7/datas/ALS/AL{i}.dat', alm, fmt='%.15f %.15f')
