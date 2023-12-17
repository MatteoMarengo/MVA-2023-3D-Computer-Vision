Student: MARENGO Matteo
Date: 10/10/2023

The objective of Lab Exercise #1 was to create a panorama using the technique of homography.

Homography is a mathematical transformation used to relate corresponding points between two images taken from slightly different viewpoints. It is represented as a 3x3 matrix (H) and is employed to map points from one image to another in a way that accounts for perspective transformations, rotations, and translations.

In this lab:

    Part 1 involves manually selecting corresponding points in the two images and obtaining their coordinates (x, y). These correspondences are used to establish a set of point pairs, which are crucial for calculating the homography matrix.

    Part 2 focuses on the mathematical foundation. The goal is to construct matrices A and B, which are essential components of a linear system that can be solved to find the elements of the homography matrix. The linear system can be represented as Ax = B, where A is constructed based on the coordinates of the corresponding points and B is determined from the coordinates of the points in the second image.

    Part 3, the assembly phase, is where the magic happens. Using the homography matrix derived from Part 2, the two images are combined. The homography matrix is used to warp one image so that it aligns correctly with the other. Pixels in the overlapping regions are blended together to ensure a seamless transition between the two images. The resulting composite image forms the panorama.

In conclusion, this lab provides experience with homographies, illustrating their mathematical strength in creating panoramas by aligning and blending multiple images, taking into account the complex geometric transformations that occur when capturing scenes from different perspectives.