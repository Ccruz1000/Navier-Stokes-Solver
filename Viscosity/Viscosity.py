import numpy as py
from function1d import fir_o_back_y2d
from function1d import fir_o_back_x2d

dudx = fir_o_back_x2d(dx, u, N)
dvdy = fir_o_back_y2d(dy, v, N)
dudy = fir_o_back_y2d(dy, u, N)
dvdx = fir_o_back_x2d(dx, v, N)

miu = 1
lamb = -2 / 3 * miu

def tau(dudx, dvdy, miu, lamb):
    tau = (lamb + 2 * miu) * dudx + lamb * dvdy
    return tau

def tauxy(dudy, dvdx, miu):
    tau = miu * (dudy + dvdx)
    return tau

