import mdtraj as md
from msmbuilder.dataset import dataset
from msmbuilder.decomposition import tICA
from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.cluster import KCenters
from msmbuilder.cluster import MiniBatchKMeans
from msmbuilder.msm import MarkovStateModel
from msmbuilder.utils import dump
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import numpy as np
from scipy.spatial import distance
import timeit

# nsteps = 100000
z = np.zeros(100000)
t, x, y = np.loadtxt('trajectories/xy_traj_100k.trj', usecols=(0,1,2), unpack=True)
t = t[1:]
x = x[1:]
y = y[1:]

trajs = np.array([])

for i in xrange(10):
    trajs = [np.column_stack((x, y, z))]
    # len(trajs) = 1, len(trajs[0]) = nsteps + 1

# Featurize trajectory
# featurizer = RawPositionsFeaturizer(atom_indices=None, ref_traj=None)
# sequences = featurizer.transform(trajs)

sequences = trajs
lag_time = 100000
tICA = tICA(n_components=3, lag_time=lag_time)
transformed = tICA.fit_transform(sequences)    
transformed = list([np.concatenate(transformed)])

# Visualize TICA
def draw_arrow(a, v, color):
    plt.arrow(0, 0, a*v[0], a*v[1], color=color, width=0.02, linewidth=3)

transformed_c = np.concatenate(transformed)
v = tICA.eigenvectors_
eigval = tICA.eigenvalues_
print v
print eigval
print v[0]
print v[:,0]
# plt.hist2d(transformed_c[:,0], transformed_c[:,1], bins=100, cmap='viridis')
plt.figure()
plt.hexbin(transformed_c[:,0], transformed_c[:,1], bins='log', mincnt=1)
draw_arrow(2*eigval[0], v[0], color='red')
plt.text(-1.5, 1.7, 'TICA', color='red', fontsize = 20, fontweight='bold') # rotation='vertical')
plt.xlabel('TICA #1')
plt.ylabel('TICA #2')
plt.colorbar()
# plt.figure()
# plt.hexbin(transformed_c[:,0], transformed_c[:,1], bins='log', mincnt=1)
# plt.xlabel('TICA #1')
# plt.ylabel('TICA #2')
# plt.figure()
# plt.hexbin(transformed_c[:,1], transformed_c[:,2], bins='log', mincnt=1)
# plt.xlabel('TICA #2')
# plt.ylabel('TICA #3')
# plt.savefig('figures/msmb_tica_t%d' % (lag_time), dpi=1080)
# plt.show()
plt.close()


# Clustering in reduced dimension
n_clusters = 80
# nframe = len(transformed[0])
cluster = MiniBatchKMeans(n_clusters=n_clusters, init_size=n_clusters*3)
xtransformed = cluster.fit_transform(transformed)  
# xtransformed = transformed.fit_transform(clusterer)              
xtransformed_centers = cluster.cluster_centers_
print len(xtransformed[0]), len(xtransformed_centers)
print


# MSM
msm = MarkovStateModel(lag_time=5)
msm.fit(xtransformed)
eig = msm.eigenvalues_
timescales = msm.timescales_
print
# print 'transmat = ', '\n', msm.transmat_, '\n', len(msm.transmat_), len(msm.transmat_[0]), '\n'
print 'eigenvalues = ', '\n', msm.eigenvalues_, '\n', len(msm.eigenvalues_), '\n'
print 'timescales = ', '\n', msm.timescales_, '\n', len(msm.timescales_), '\n'

zipped = zip(eig, eig[1:], timescales)
# np.savetxt('2dwell_msm.txt', zipped)

