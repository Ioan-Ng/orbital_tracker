import planet_data as pd
import numpy as np
import pyvista as pv
from pyvista import examples
import tools as tools
# sphere = pv.Sphere()
# pl = pv.Plotter(shape = (1,1))
# pl.subplot(0,0)
# pl.add_mesh(sphere)
# pl.show()
# globe = examples.load_globe()
# plane = pv.examples.load_airplane()
# pl = pv.Plotter()


# pl.add_mesh(plane)
# pl.add_points(plane.points,color = "red", render_points_as_spheres=True)
# pl.show_bounds()
# pl.show_axes()
# pl.show()
cb = pd.earth
earth = examples.planets.load_earth(radius=cb["radius"])
earth_texture = examples.load_globe_texture()


pl = pv.Plotter()

pl.show_axes()
pl.show_bounds()

point = np.array([0,0,3])



earth_point = pv.MultiBlock([earth,point])
earth_point.rotate_y(23.44,inplace=True)
pl.add_mesh(earth,texture=earth_texture)

rs = tools.orbitsPropagate(cb = pd.earth, file = "SpaceStationData.csv")
r = rs[0]
rGround  = np.array(r)
rGround *=.97
print(rGround)
zeros = np.zeros((86400,3))
#line_mesh = pv.MultipleLines(rs)
pl.add_points(point)
pl.show_axes()
pl.show_bounds()


zero_line_mesh = pv.MultipleLines(zeros)
zero_groundline_mesh = pv.MultipleLines(zeros)
earth_ground = pv.MultiBlock([earth, zero_groundline_mesh])
pl.add_mesh(zero_line_mesh, color = "red")
pl.add_mesh(zero_groundline_mesh, color = "pink")
pl.open_gif("earth_rotate.gif",fps = 60, iterations = 1)
frames = 24*60*60

for i in range(frames):
    earth_ground.rotate_z(0.00417*60, inplace=True)
    zero_line_mesh.points[i] = r[i*60]
    zero_groundline_mesh.points[i] = rGround[i*60]
    pl.write_frame()
    
