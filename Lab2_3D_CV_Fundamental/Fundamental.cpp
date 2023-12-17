// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

// Student: MARENGO Matteo
// Mail: matteo.marengo@ens-paris-saclay.fr
// Date: 17/10/2023

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>

#include <cstdlib>
#include <ctime>

using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// #1: SIFT ----------------------------------------------- //
//--------------------------------------------------------- //
// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

// #2: RANSAC ----------------------------------------------- //
// ##################################################
// #2.1: Function to compute the fundamental matrix F
// ##################################################

FMatrix<float, 3, 3> computeFundamentalMatrix(vector<Match>& matches) {
    // Matrix for linear equation
    FMatrix<float, 9, 9> linearSystemMatrix;
    for (int idx = 0; idx < 9; idx++) {
        if (idx < 8) {
            Match matchData = matches[idx];

            // Extract coordinates
            DoublePoint3 firstPoint;
            DoublePoint3 secondPoint;
            // Normalized points
            float x1n,y1n,x2n,y2n;
            x1n = 0.001 * matchData.x1;
            y1n = 0.001 * matchData.y1;
            x2n = 0.001 * matchData.x2;
            y2n = 0.001 * matchData.y2;
            firstPoint = { x1n, y1n, 1 };
            secondPoint = { x2n, y2n, 1 };

            // Fill in linear system with 9 unknowns
            linearSystemMatrix(idx, 0) = firstPoint[0] * secondPoint[0];
            linearSystemMatrix(idx, 1) = firstPoint[0] * secondPoint[1];
            linearSystemMatrix(idx, 2) = firstPoint[0];
            linearSystemMatrix(idx, 3) = firstPoint[1] * secondPoint[0];
            linearSystemMatrix(idx, 4) = firstPoint[1] * secondPoint[1];
            linearSystemMatrix(idx, 5) = firstPoint[1];
            linearSystemMatrix(idx, 6) = secondPoint[0];
            linearSystemMatrix(idx, 7) = secondPoint[1];
            linearSystemMatrix(idx, 8) = 1;
        } else {
            for (int j = 0; j < 9; j++) { // We use the Hint
                linearSystemMatrix(8, j) = 0;
            }
        }
    }

    // Determine fundamental matrix using SVD
    // A can be diagonalised as A = USV
    // S is a diagonal matrix, which diagonal is the eigenvalues of A
    FVector<float, 9> singularValues; // S
    FMatrix<float, 9, 9> leftMatrix, rightMatrixTranspose; // U & V
    svd(linearSystemMatrix, leftMatrix, singularValues, rightMatrixTranspose);

    FMatrix<float, 3, 3> fundamentalMatrix;
    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            fundamentalMatrix(m, n) = rightMatrixTranspose.getRow(8)[3 * m + n];
        }
    }

    // Enforce rank-2 constraint
    FVector<float, 3> constrainedSingularValues;
    FMatrix<float, 3, 3> leftSubMatrix, rightSubMatrixTranspose;
    // 1: We compute its svd
    svd(fundamentalMatrix, leftSubMatrix, constrainedSingularValues, rightSubMatrixTranspose);
    // 2: We put sigma3 = 0
    constrainedSingularValues[2] = 0; // sigma3 = 0
    // 3: We recompute F using the modified constraint
    fundamentalMatrix = leftSubMatrix * Diagonal(constrainedSingularValues) * rightSubMatrixTranspose;

    // Setup normalization matrix
    FMatrix<float, 3, 3> normalizationMatrix(0.f);
    normalizationMatrix(0, 0) = 0.001;
    normalizationMatrix(1, 1) = 0.001;
    normalizationMatrix(2, 2) = 1;

    // Denormalize resulting matrix
    fundamentalMatrix = normalizationMatrix * fundamentalMatrix * normalizationMatrix;
    fundamentalMatrix = fundamentalMatrix / norm(fundamentalMatrix);
    return fundamentalMatrix;
}

// ############################################################################################
// #2.1 bis: Function to refine resulting F with least square minimization based on all inliers
// ############################################################################################
FMatrix<float,3,3> computeFundamentalMatrixLS(vector<Match>& matches) {

    // Matrix for linear equation
    int sizemat = matches.size();
    // Dynamic allocation of the matrix
    Matrix<float> linearSystemMatrix(sizemat,9);
    // Fill the linear matrix A
    for (int idx = 0; idx < sizemat; idx++) {

        Match matchData = matches[idx];
        // Extract coordinates
        DoublePoint3 firstPoint;
        DoublePoint3 secondPoint;
        // Normalized points
        float x1n,y1n,x2n,y2n;
        x1n = 0.001 * matchData.x1;
        y1n = 0.001 * matchData.y1;
        x2n = 0.001 * matchData.x2;
        y2n = 0.001 * matchData.y2;
        firstPoint = { x1n, y1n, 1 };
        secondPoint = { x2n, y2n, 1 };

        // Fill in linear system with 9 unknowns
        linearSystemMatrix(idx, 0) = firstPoint[0] * secondPoint[0];
        linearSystemMatrix(idx, 1) = firstPoint[0] * secondPoint[1];
        linearSystemMatrix(idx, 2) = firstPoint[0];
        linearSystemMatrix(idx, 3) = firstPoint[1] * secondPoint[0];
        linearSystemMatrix(idx, 4) = firstPoint[1] * secondPoint[1];
        linearSystemMatrix(idx, 5) = firstPoint[1];
        linearSystemMatrix(idx, 6) = secondPoint[0];
        linearSystemMatrix(idx, 7) = secondPoint[1];
        linearSystemMatrix(idx, 8) = 1;
    }

    // Determine fundamental matrix using SVD
    // F can be diagonalised as A = USV
    // S is a diagonal matrix, which diagonal is the eigenvalues of A
    //linearSystemMatrix = transpose(linearSystemMatrix)*linearSystemMatrix;
    Vector<float> singularValues(sizemat); // S
    Matrix<float> leftMatrix(sizemat,sizemat), rightMatrixTranspose(9,9); // U & V
    svd(linearSystemMatrix, leftMatrix, singularValues, rightMatrixTranspose);

    FMatrix<float,3,3> fundamentalMatrix;
    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            fundamentalMatrix(m, n) = rightMatrixTranspose.getRow(8)[3 * m + n];
        }
    }

    // Enforce rank-2 constraint
    FVector<float,3> constrainedSingularValues;
    FMatrix<float,3,3> leftSubMatrix, rightSubMatrixTranspose;
    // 1: We compute its svd
    svd(fundamentalMatrix, leftSubMatrix, constrainedSingularValues, rightSubMatrixTranspose);
    // 2: We put sigma3 = 0
    constrainedSingularValues[2] = 0; // sigma3 = 0
    // 3: We recompute F using the modified constraint
    fundamentalMatrix = leftSubMatrix * Diagonal(constrainedSingularValues) * rightSubMatrixTranspose;

    // Setup normalization matrix
    FMatrix<float,3,3> normalizationMatrix(0.f);
    normalizationMatrix(0, 0) = 0.001;
    normalizationMatrix(1, 1) = 0.001;
    normalizationMatrix(2, 2) = 1;

    // Denormalize resulting matrix
    fundamentalMatrix = normalizationMatrix * fundamentalMatrix * normalizationMatrix;
    fundamentalMatrix = fundamentalMatrix / norm(fundamentalMatrix);
    return fundamentalMatrix;
}

// #########################
// #2.2: ComputeF algorithm
//##########################
// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    // Define some variables
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF; // The matrix F with the most inliers
    vector<int> bestInliers;

    // --------------- TODO ------------ //
    // DO NOT FORGET NORMALIZATION OF POINTS
    FMatrix<float,3,3> F; // The matrix F we want to compute
    vector<int> inliers;
    int NbMatches = matches.size(); // The number of matches computed by the SIFT algorithm
    int const SAMPLE_SIZE = 8; // Size of the sample (8 for the 8-point method)
    int currentSampleNum = 0; // Number of the sample being studied
    // Check that there are enough matches
    if(SAMPLE_SIZE > NbMatches){
        cerr << "Insufficient matches for sample creation" << endl;
    }
    while(currentSampleNum < Niter){
        // Create the matches
        vector<Match> SVDMatches;
        vector<int> selectedIndices;
        selectedIndices.reserve(SAMPLE_SIZE);// Allocate memory of size 8
        while (selectedIndices.size() < SAMPLE_SIZE) {
            int randomIndex = rand() % NbMatches;// It returns a value between 0 and NbMatches (exclusive)
            // Check whether the randomindex is not already selected
            if (find(selectedIndices.begin(), selectedIndices.end(), randomIndex) == selectedIndices.end()) {
                selectedIndices.push_back(randomIndex);
                SVDMatches.push_back(matches[randomIndex]);
            }
        }

        // Compute Matrix F with the matches that were extracted
        F = computeFundamentalMatrix(SVDMatches);
        // Compute the number of inliers
        vector<int> inliers;
        for (int idx = 0; idx < NbMatches; idx++) {
            Match matchItem = matches[idx];
            DoublePoint3 pt1 = { matchItem.x1, matchItem.y1, 1 };
            DoublePoint3 pt2 = { matchItem.x2, matchItem.y2, 1 };
            // distance of xi' to epipolar line associated to xi (F^Txi)
            FVector<float, 3> computedLine = F * pt2;
            float magnitude = sqrt(computedLine[0] * computedLine[0] + computedLine[1] * computedLine[1]);
            computedLine = computedLine / magnitude;
            float dist = fabs(pt1 * computedLine);
            if (dist < distMax) {
                inliers.push_back(idx);
            }
        }

        if(inliers.size() > bestInliers.size() && inliers.size()>50){
            bestInliers = inliers;
            bestF = F;
            cout <<"Best number of Inliers: " << bestInliers.size()<<endl;
            int nIterNb;
            float m = bestInliers.size();
            float n = NbMatches;
            float k = SAMPLE_SIZE;
            nIterNb = log(BETA)/log(1- pow( m/n, k));
            if (nIterNb<Niter){
                Niter = nIterNb;
            }
        }
        currentSampleNum = currentSampleNum + 1;
    }
    // ----------- END TO DO ---------- //////
    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);

    // Refine resulting F with least square minimization based on all inliers
    FMatrix<float,3,3> bestFLS;
    bestFLS = computeFundamentalMatrixLS(matches);

    // You can return either bestF or bestFLS to compare the difference between before and after lS minimization on all inliers
    return bestFLS;
}

// #3: EPIPOLAR----------------------------------------------- //
// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    int nbepipoleft = 0,nbepiporight = 0;
    while(true) {
        int x,y;
        // Loop until the right mouse button is clicked
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------
        // Determine which image was clicked, compute the associated epipolar line
        // in the alternative image, and showcase it
        DoublePoint3 clickedPoint;
        int imageWidth = I1.width();
        // If the x-coordinate of the clicked point is greater than the width of the first image it belongs to the
        // second image
        if (x > imageWidth) { // Click in RIGHT - Display in LEFT
            cout <<"Right Image clicked"<<endl;
            cout <<"Number of epipolar lines in Left Image: "<<nbepiporight<<endl;
            clickedPoint = { (double)(x - imageWidth), (double)y, 1 };
            FVector<float, 3> computedLine = F * clickedPoint;
            // Equation of line ax + by + c = 0
            // a = computedLine[0]
            // b = computedLine[1]
            // c = computedLine[2]
            // y1 is the coordinate when x = 0
            // y2 is the coordinate when x = imageWidth
            float y1 = (-1) * computedLine[2] / computedLine[1];
            float y2 = (-1) * (computedLine[2] + computedLine[0] * imageWidth) / computedLine[1];
            drawLine(0,(int)y1 , imageWidth, (int)y2, RED);
            nbepiporight = nbepiporight +1;
        } else { // Click in LEFT - Display in RIGHT
            cout <<"Left Image clicked"<<endl;
            cout <<"Number of epipolar lines in Right Image: "<<nbepipoleft<<endl;
            clickedPoint = { (double)x, (double)y, 1 };
            FVector<float, 3> computedLine = transpose(F) * clickedPoint;
            float y1 = (-1) * computedLine[2] / computedLine[1];
            float y2 = (-1) * (computedLine[2] + computedLine[0] * imageWidth) / computedLine[1];
            drawLine(imageWidth, (int)y1, 2 * imageWidth, (int)y2, BLUE);
            nbepipoleft = nbepipoleft + 1;
        }
        // ------------- END TO DO --------------
    }
}

// #4: MAIN FUNCTION ----------------------------------------------- //
int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    // SIFT matches
    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100,20,std::to_string(n)+ " matches",RED);
    click();
    
    // RANSAC
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    drawString(100, 20, to_string(matches.size())+"/"+to_string(n)+" inliers", BLUE);
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
