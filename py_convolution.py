from scipy.io import loadmat
import numpy as np, matplotlib.pyplot as plt
from IPython import embed as IP
from numpy.random import uniform
from numpy import convolve, cumsum, histogram, linspace, ceil, sqrt

def tmp():
    data_mat = loadmat('./caiso.mat')
    dict_load = dict()
    y0 = 2016
    for i in range(0, data_mat['load'].shape[0]):
        dict_load[str(y0+i)] = data_mat['load'][i, 0]

    pdf_x, bin_edges = np.histogram(dict_load['2016'], density=False, bins= 50)
    # plt.plot(pdf_x/float(dict_load['2016'].size))
    plt.hist(dict_load['2016'])
    plt.show()

def convolution2():
    # From https://stackoverflow.com/a/6526507

    s, e, n= -0.5, 0.5, 1000
    x, y, bins= uniform(s, e, n), uniform(s, e, n), linspace(s, e, n** .75)
    pdf_x= histogram(x, normed= True, bins= bins)[0]
    pdf_y= histogram(y, normed= True, bins= bins)[0]
    c= convolve(pdf_x, pdf_y); c= c/ c.sum()
    bins= linspace(2* s, 2* e, len(c))
    # a simulation
    xpy= uniform(s, e, 10* n)+ uniform(s, e, 10* n)
    c2= histogram(xpy, normed= True, bins= bins)[0]; c2= c2/ c2.sum()

    IP()

    from pylab import grid, plot, show, subplot
    subplot(211), plot(bins, c)
    plot(linspace(xpy.min(), xpy.max(), len(c2)), c2, 'r'), grid(True)
    subplot(212), plot(bins, cumsum(c)), grid(True), show()

def convolution():
    n = 10**2

    x,y = uniform(-0.5,0.5,n),uniform(-0.5,0.5,n)

    bins = int(ceil(sqrt(n)))

    pdf_x = histogram(x,bins=bins,normed=True)
    pdf_y = histogram(y,bins=bins,normed=True)

    s = convolve(pdf_x[0],pdf_y[0])

    plt.plot(s)
    plt.show()


if __name__ == "__main__":
    # tmp()
    # convolution()
    convolution2()