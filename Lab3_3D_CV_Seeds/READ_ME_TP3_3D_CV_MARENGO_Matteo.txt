Student: MARENGO Matteo
Date: 24/10/2023

The objective of Lab Exercise #3 was to compute the disparity map associated to a pair of images.

A disparity map is a representation of the difference in horizontal positions of corresponding points in two stereo images, which are images of the same scene taken from slightly different viewpoints. This difference, known as disparity, is inversely proportional to the depth of objects in the scene. In other words, objects that are closer to the camera will have a larger disparity, while objects that are farther away will have a smaller disparity. By calculating the disparity between corresponding points in the stereo images, it is possible to create a depth map, which can be used for 3D reconstruction of the scene. In addition,  seeds typically refer to initial disparity values that have been estimated with high confidence. These seed disparities are then used as starting points to propagate disparity information to neighboring pixels in the image.

In this lab:

    Part 1 involves computing disparity map from image 1 to 2 of all points by highest NCC score. This means that we compute the NCC score between two images and keep the one with the highest disparity. The NCC value ranges from -1 to 1, with 1 indicating a perfect match between the two patches, 0 indicating no correlation, and -1 indicating a perfect inverse correlation.

    Part 2 focuses on keeping disparity where NCC is sufficiently high so we can consider them as seeds. 

    Part 3, the propagation of seeds, we keep the seeds for the neighbor that has the highest NCC score among dP-1, dP and dP+1. This propagation process helps to fill in the disparity map, creating a more complete and accurate representation of the scene's depth.

