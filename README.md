# SandEngine_Depth

Estimating depth in 2D imagery


In this project you will test some of these deep learning methods. First you will choose and apply some
methods on aerial imagery obtained by a drone, together with a LiDAR dataset acquired at the same
time, in order to validate the estimation. As a drone is still expensive, and aerial imagery might be less
applicable for depth estimation, you will also use so-called ARGUS imagery to test the methods. First
you will design a fieldwork setup to calibrate the imagery from the ARGUS tower using your own ground
control points (GCPs). Then you will use your obtained transformation matrix to transform a 3D point
cloud into the image plane of the ARGUS imagery. Afterward you will estimate depth in these images
and validate this using the reprojected point clouds.
