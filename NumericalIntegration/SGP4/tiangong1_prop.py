import numpy as np
import math

from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv

line1 = ('1 37820U 11053A   18012.22140694  .00063269  99835-5  13071-3 0  9990')
line2 = ('2 37820  42.7537 344.4268 0017667 147.3056 342.3989 15.98674657361189')

satellite = twoline2rv(line1, line2, wgs84)
print(satellite.error)
print(satellite.error_message)
position, velocity = satellite.propagate(2018, 8, 6, 03, 8, 54)
print(position)
print(velocity)
print(np.linalg.norm(position) - 6378.137)
