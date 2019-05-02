import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.grid(False)
ax.set_axis_off()
ax.set_facecolor((.95,.95,.95))

s_pts = np.array([
    [1,0,1],
    [1,0,-1.5],
    [-1,0,1],
    [-1,0,-1.5],
    ])

u_pts = np.array([
    [0,0,1],
    [2,0,1],
    [-2,0,1],
    ])

v_pts = np.array([
    [1,1.5,1],
    [1,-1.5,1],
    [-1,1.5,1],
    [-1,-1.5,1],
    ])

w_pts = np.array([
    [1,0,-.25],
    [1,0,2.25],
    [-1,0,-.25],
    [-1,0,2.25],
    ])

lines=[
    # x lines
    [[-10,0,1],[10,0,1]],
    [[-10,0,-1.5],[10,0,-1.5]],
    # y lines
    [[1,-10,1],[1,10,1]],
    [[1,-10,-1.5],[1,10,-1.5]],
    [[-1,-10,1],[-1,10,1]],
    #[[-1,-5,-1.5],[-1,5,-1.5]],
    # z lines
    [[-1,0,-10],[-1,0,10]],
    [[1,0,-10],[1,0,10]],
    ]

ax.add_collection(Line3DCollection(lines, colors='grey'))
ax.scatter(s_pts[:,0],s_pts[:,1],s_pts[:,2],s=50,color='k',depthshade=False)
ax.scatter3D(u_pts[:,0],u_pts[:,1],u_pts[:,2],s=100,color='b',marker='|',depthshade=False)
ax.scatter3D(v_pts[:,0],v_pts[:,1],v_pts[:,2],s=100,color='r',marker='_',depthshade=False)
ax.scatter3D(w_pts[:,0],w_pts[:,1],w_pts[:,2],s=100,color='g',marker='_',depthshade=False)

#ax.text(3,0,0,s='x',color='grey')
#ax.text(0,3,0,s='y',color='grey')
#ax.text(0,0,3,s='z',color='grey')
# s labels
ax.text(.33,0,.7,s=r'$s_{i,j,k}$',color='k')
ax.text(1.25,0,-1.4,s=r'$s_{i,j,k-1}$',color='k')
ax.text(-.75,0,1.05,s=r'$s_{i-1,j,k}$',color='k')

# u labels
ax.text(2.1 ,0,1.05,s=r'$u_{i+1,j,k}$',color='b')
ax.text(.1,0,1.05,s=r'$u_{i,j,k}$',color='b')
#ax.text(-1.9,0,1.1,s=r'$u_{i-1,j,k}$',color='b')

# v labels
ax.text(1.2,1.5,1,s=r'$v_{i,j,k}$',color='r')
ax.text(1.2,-1.5,1,s=r'$v_{i,j,k-1}$',color='r')
#ax.text(-.8,1.5,1,s=r'$v_{i-1,j,k}$',color='r')

# w labels
ax.text(1.15,0,2.25,s=r'$w_{i,j,k}$',color='g')
ax.text(1.15,0,-.25,s=r'$w_{i,j,k-1}$',color='g')
#ax.text(-.8,0,2,s=r'$w_{i-1,j,k}$',color='g')

lims=(-2,2)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_zlim(lims)

# Save to pdf
pdf = PdfPages('/home/sdomke/workspace/thesis/masters/data/arakawaC.pdf')
pdf.savefig(fig)
pdf.close()