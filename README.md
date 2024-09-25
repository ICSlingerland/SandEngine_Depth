# SandEngine_Depth

Estimating depth in 2D imagery


In this project we will test some deep learning methods. First we will choose and apply some
methods on aerial imagery obtained by a drone, together with a LiDAR dataset acquired at the same
time, in order to validate the estimation. As a drone is still expensive, and aerial imagery might be less
applicable for depth estimation, we will also use so-called ARGUS imagery to test the methods. 

We will design a fieldwork setup to calibrate the imagery from the ARGUS tower using our own ground
control points (GCPs). Then use the obtained transformation matrix to transform a 3D point
cloud into the image plane of the ARGUS imagery. Afterward we will estimate depth in these images
and validate this using the reprojected point clouds.


