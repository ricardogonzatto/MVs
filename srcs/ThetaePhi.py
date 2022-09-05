# O retorno da função hp.pix2ang() é theta, phi
# Ambas as listas possuem valores repetidos pq hora mantemos theta fixo para variar phi
# hora variamos phi mantendo theta fixo

import numpy as np
import healpy as hp

nside = 8

coord = []
for i in range(hp.nside2npix(nside)):
    m = hp.pix2ang(8,i)
    coord.append(m)

n = np.transpose(coord)
np.savetxt("../datas/ThetaPhi.dat", np.column_stack([n[0], n[1]]))
