# The Algorithm

[Multilateration](https://en.wikipedia.org/wiki/Multilateration) is the underlying principle and in general is very simple.  Imagine a simplified experiment with receivers placed in an infinite, homogeneous, acoustic medium.  There is one source and the observed arrival at all receivers is a single impulse.  It is then straightforward to identify the arrival time at any station and compute the differential travel-time for any station pair.  Realize that an infinite number of sources can produce the observed differential travel-time.  The consequence of the homogeneous model is that the candidate sources will be distributed along a hyperbola.  Hence, the location algorithm proceeds by overlaying and summing many hyperbolas.  The intersection of the hyperbolas generated from the differential travel-times in turn produces the event location.  

Of course, differential travel-times and migration in homogeneous media is a bit too simplistic.  Instead, the travel-time migration of processed correlograms is performed in arbitrarily heterogeneous media.  Thus, what is really done under the hood is

  1.  Phase-correlograms <a href="https://www.codecogs.com/eqnedit.php?latex=C_{j,i}(t)&space;=&space;\mathcal{F}^{-1}&space;\left&space;\{&space;\frac{U_i(\omega)&space;U_j^*(\omega)}{|U_i(\omega)&space;U_j(\omega)|}&space;\right&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?C_{j,i}(t)&space;=&space;\mathcal{F}^{-1}&space;\left&space;\{&space;\frac{U_i(\omega)&space;U_j^*(\omega)}{|U_i(\omega)&space;U_j(\omega)|}&space;\right&space;\}" title="C_{j,i}(t) = \mathcal{F}^{-1} \left \{ \frac{U_i(\omega) U_j^*(\omega)}{|U_i(\omega) U_j(\omega)|} \right \}" /></a> are computed for the all the differential arrival times.  
  2.  Then, the RMS averaging window or envelope of the phase-correlograms are computed.  
  3.  Finally, the processed phase-correlograms are migrated in the provided travel-time models.  

While this algorithm is quite simple and only requires a few lines of Python code the actual implementation is much more complicated.  Particularly because of parallelism and performance are essential.   

