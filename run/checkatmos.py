import h5py
import matplotlib.pyplot as plt
import numpy as np

def relative_gradient(X):

    dx = (X - np.roll(X,1,axis=1)) / X
    dy = (X - np.roll(X,1,axis=0)) / X
    dy[:,0]=0.0
    dr = np.zeros(dx.shape)
    dr = np.copy(np.abs(dx))
    mask = np.where(np.abs(dy) > np.abs(dx)) 
    dr[mask] = np.abs(dy[mask])
    return dr

def refine_mask_gradient(X, thresh = 0.1, psize = 6):

    nx, ny = X.shape
    npx= int(nx / psize)
    npy= int(ny / psize)
    rmask = np.zeros((nx, ny))
    
    for j in range(npy):
        jj = j * psize
        for i in range(npx):
            ii = i * psize
        
            patch = X[jj:jj+psize, ii:ii+psize]

            if patch.max() > thresh:
                rmask[jj:jj+psize, ii:ii+psize] = 1.0

    return rmask

def coarsen_mask_gradient(X, thresh = 0.01, psize = 6):

    nx, ny = X.shape
    npx= int(nx / psize)
    npy= int(ny / psize)
    rmask = np.zeros((nx, ny))
    
    for j in range(npy):
        jj = j * psize
        for i in range(npx):
            ii = i * psize
        
            patch = X[jj:jj+psize, ii:ii+psize]

            if patch.max() < thresh:
                rmask[jj:jj+psize, ii:ii+psize] = 1.0

    return rmask


def refine_mask_minmax(X, vmin,vmax, psize = 6):

    nx, ny = X.shape
    npx= int(nx / psize)
    npy= int(ny / psize)
    rmask = np.zeros((nx, ny))
    
    for j in range(npy):
        jj = j * psize
        for i in range(npx):
            ii = i * psize
            
            patch = X[jj:jj+psize, ii:ii+psize]

            if patch.max() > vmin and patch.min() < vmax:
                rmask[jj:jj+psize, ii:ii+psize] = 1.0

    return rmask

def coarsen_mask_minmax(X, vmin,vmax, psize = 6):

    nx, ny = X.shape
    npx= int(nx / psize)
    npy= int(ny / psize)
    rmask = np.zeros((nx, ny))
    
    for j in range(npy):
        jj = j * psize
        for i in range(npx):
            ii = i * psize
            
            patch = X[jj:jj+psize, ii:ii+psize]

            if patch.max() < vmin or patch.min() > vmax:
                rmask[jj:jj+psize, ii:ii+psize] = 1.0

    return rmask

plt.close("all")

f = h5py.File('inputatmos.hdf5', 'r')

r = np.array(f["mass_density"])
tg = np.array(f["temperature"])
ne = np.array(f["electron_density"])
x = np.array(f["x"])
y = np.array(f["y"])
f.close()

dtg = relative_gradient(tg)
tmask1 = refine_mask_gradient(dtg, thresh=0.1)
tmask2 = refine_mask_minmax(tg, vmin = 0, vmax = 1e5)
tmask = np.logical_and(tmask1>0, tmask2>0)

tmask1 = coarsen_mask_gradient(dtg, thresh=0.01)
tmask2 = coarsen_mask_minmax(tg, vmin = 0, vmax = 1e5)
tmask_coarsen = np.logical_or(tmask1>0, tmask2>0)


dr = relative_gradient(r)
rmask1 = refine_mask_gradient(dr, thresh=0.2)
rmask2 = refine_mask_minmax(r, vmin = 1e-13, vmax = 1e-6)
rmask = np.logical_and(rmask1>0, rmask2>0)



dne = relative_gradient(ne)
nemask1 = refine_mask_gradient(dne, thresh=0.2)
nemask2 = refine_mask_minmax(ne, vmin = 1e10, vmax = 1e14)
nemask = np.logical_and(nemask1>0, nemask2>0)

mask = np.logical_or(nemask>0, np.logical_or(rmask > 0, tmask > 0) > 0)



f1, axarr = plt.subplots(3, 3, sharex='col', sharey='row')
axarr[0,0].imshow(tmask1,cmap="Greys_r")
axarr[0,1].imshow(tmask2,cmap="Greys_r")
axarr[0,2].imshow(np.log(tg), cmap="Greys_r")
axarr[0,1].set_title("temperature conditions")

axarr[1,0].imshow(rmask1,cmap="Greys_r")
axarr[1,1].imshow(rmask2,cmap="Greys_r")
axarr[1,2].imshow(np.log(r), cmap="Greys_r")
axarr[1,1].set_title("density conditions")

axarr[2,0].imshow(nemask1,cmap="Greys_r")
axarr[2,1].imshow(nemask2,cmap="Greys_r")
axarr[2,2].imshow(np.log(ne), cmap="Greys_r")
axarr[2,1].set_title("electron density conditions")

f1.tight_layout()
f1.show()

k=np.zeros((768,768))

a=np.where(mask>0)
k[a] = 1
a=np.where(mask==0)
k[a] -= tmask_coarsen[a]

f2,ax = plt.subplots(1,1)
ax.imshow(k,cmap="bwr")
ax.set_title("combined refinement map")
f2.show()

        
