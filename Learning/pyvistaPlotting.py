import numpy as np
import pyvista as pv
from pyvista import examples
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

earth = pv.Sphere(radius =3, center = (0,0,0), direction = (0,0,1))
sphere = pv.Sphere(radius =3, center = (0,0,0), direction = (0,0,1))

pl = pv.Plotter(shape = (2,1))
pl.subplot(0,0)
pl.show_axes()
pl.show_bounds()

point = np.array([0,0,3])
theta = np.linspace(0, 2 * np.pi, 100)
phi = np.linspace(0, 2 * np.pi, 100)
x = 3 * np.sin(phi) 
y = 3 * np.cos(phi) 
z = np.full_like(phi, 1)
points = np.column_stack((x, y,z))
spline = pv.Spline(points, 100)
pl.add_mesh(spline)
earth_point = pv.MultiBlock([earth,point])
earth_point.rotate_y(23.44,inplace=True)
pl.add_mesh(earth_point)

pl.subplot(1,0)
pl.add_mesh(sphere)
pl.add_points(point)
pl.show_axes()
pl.show_bounds()
pl.show()