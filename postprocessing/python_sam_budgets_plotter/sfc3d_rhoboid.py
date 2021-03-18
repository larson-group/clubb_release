import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Setup
ax.grid(False)
ax.set_axis_off()
#ax.set_axis_on()
#ax.set_frame_on(True)

# Define list of point coordinates
points = np.array([
                # bottom, top
                [0, 0, -2],
                [0, 0, 2],
                # layer 1
                [0, 0, 1],
                [1, 0, 1],
                [0, 1, 1],
                [-1, 0, 1],
                [0, -1, 1],
                # layer -1
                [0, 0, -1],
                [1, 0, -1],
                [0, 1, -1],
                [-1, 0, -1],
                [0, -1, -1],
                # layer 0
                # outer cross points
                [2, 0, 0],
                [0, 2, 0],
                [-2, 0, 0],
                [0, -2, 0],
                # inner cross points
                #[0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [-1, 0, 0],
                [0, -1, 0],
                # dagonal points
                [1, 1, 0],
                [1, -1, 0],
                [-1, 1, 0],
                [-1, -1, 0]])

# list of sides' polygons of figure
verts = [[[0,0,2],[2,0,0],[0,2,0]],
         [[0,0,2],[2,0,0],[0,-2,0]],
         [[0,0,2],[-2,0,0],[0,2,0]],
         [[0,0,2],[-2,0,0],[0,-2,0]],
         [[0,0,-2],[2,0,0],[0,2,0]],
         [[0,0,-2],[2,0,0],[0,-2,0]],
         [[0,0,-2],[-2,0,0],[0,2,0]],
         [[0,0,-2],[-2,0,0],[0,-2,0]],
         ]

# Axis limits
lims=(-1.5,1.6)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot vertices
ax.scatter3D(points[:, 0], points[:, 1], points[:, 2], color='k')
ax.scatter3D(0,0,0,s=100,color='r')
#lines = [[[0,0,2],[2,0,0]],
         #[[0,0,2],[-2,0,0]],
         #[[0,0,2],[0,2,0]],
         #[[0,0,2],[0,-2,0]],
         #[[0,0,-2],[2,0,0]],
         #[[0,0,-2],[-2,0,0]],
         #[[0,0,-2],[0,2,0]],
         #[[0,0,-2],[0,-2,0]],
         #[[2,0,0],[0,2,0]],
         #[[2,0,0],[0,-2,0]],
         #[[-2,0,0],[0,2,0]],
         #[[-2,0,0],[0,-2,0]],
        #]


# plot sides
ax.add_collection3d(Poly3DCollection(verts, facecolors=(.8,.8,.8,.2), linewidths=1, edgecolors='k'))

#ax.plot([-10,10],[0,0],zs=[0,0],color='k')
#ax.plot([0,0],[-10,10],zs=[0,0],color='k')
#ax.plot([0,0],[0,0],zs=[-10,10],color='k')
ax.quiver([0,0,0],[0,0,0],[0,0,0],[6,0,0],[0,6,0],[0,0,6],color='k', pivot='middle', arrow_length_ratio=.05)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_zlim(lims)

# Axis labels
ax.text(2.9,-0.6,0,'x')
ax.text(.2,2.8,0,'y')
ax.text(0.2,0,2.8,'z')
# Point descriptions
ax.text(1.9,0.1,0.1,'(i+2,j,k)')
ax.text(.2,1.9,0,'(i,j+2,k)')
ax.text(0.1,0.1,1.9,'(i,j,k+2)')
# Center
ax.text(0.1,.2,0,'(i,j,k)',color='r',fontsize=11)
