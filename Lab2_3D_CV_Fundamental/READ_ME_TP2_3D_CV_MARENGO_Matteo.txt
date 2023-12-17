Student: MARENGO Matteo
Date: 17/10/2023

The objective of Lab Exercise #2 was to compute the fundamental matrix with RANSAC Algorithm.

The RANdom SAmple Consensus (RANSAC) algorithm stands out as a method tailored to estimate model parameters when data countains outliers. In the context of deriving the Fundamental matrix, RANSAC operates by iteratively selecting 8 random point correspondences, estimating the matrix via the 8-point method, and then evaluating inliers based on their alignment with the projected matrix. It identifies the best matrix estimation with the maximum inliers. 

In this lab:

	Part 1: In the RANSAC algorithm, the first step is about randomly selecting a minimal subset of 8 data points required to estimate the model. This is achieved by randomly selecting 8 matches from the pool of matches detected by the SIFT algorithm. The randomness ensures a broad exploration of potential inliers.

	Part 2: Once the 8 matches are selected, the next step is to compute the Fundamental matrix using these matches. The 8-point algorithm is the method that's used here. 

	Part 3: With a candidate F, the script determines how well this matrix fits the entire dataset. A match is considered an "inlier" if its corresponding epipolar line (determined by F) is close enough to the actual point in the other image. The threshold for this closeness is defined by distMax. By counting the number of inliers, the script can assess the quality of the estimated F.

	Part 4: RANSAC is an iterative algorithm. However, instead of using a fixed number of iterations, the script dynamically adjusts the number of iterations based on the proportion of inliers found. This is an optimization to make RANSAC run faster while still maintaining its robustness. The more inliers it finds in an iteration, the fewer iterations it believes it needs to continue.

	Part 5: Once F is computed, it can be used to derive epipolar lines for any point in one image on the other image. By visualizing these lines, one can visually assess how well the estimated F aligns with the actual geometry of the scene.