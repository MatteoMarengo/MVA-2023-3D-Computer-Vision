// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Date:     2013/10/08

// Student : MARENGO Matteo - MVA (matteo.marengo@ens-paris-saclay.fr)
// Date: 2023/10/10

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

//// ----------------------------------------------------------------////
// Record clicks in two images, until right button click
void getClicks(Window w1, Window w2,
               vector<IntPoint2>& pts1, vector<IntPoint2>& pts2) {

    // ------------- TODO/A completer ----------

    int x, y,nbpoints = 0;
    while (true){
        setActiveWindow(w1); //Image 1 is on
        showWindow(w1); //Image 1 is shown
        getMouse(x,y);
        pts1.push_back(IntPoint2(x,y));
        drawCircle(x,y,10,RED);

        setActiveWindow(w2); //Image 2 is on
        showWindow(w2); //Image 2 is shown
        int click = getMouse(x,y);
        pts2.push_back(IntPoint2(x,y));
        drawCircle(x,y,10,BLUE);

        nbpoints = nbpoints + 1;
        cout << "There are "<<nbpoints<<" point selected !"<< endl;
        if (click == 3 && nbpoints >=4){
            break;//Leave this function
        }
    }

    // ------ END TO DO ------
}
//// ----------------------------------------------------------------////
// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2>& pts1,
                            const vector<IntPoint2>& pts2) {
    size_t n = min(pts1.size(), pts2.size());
    if(n<4) {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }
    Matrix<double> A(2*n,8);
    Vector<double> B(2*n);
    // ------------- TODO/A completer ----------
    IntPoint2 pt1,pt2;
    for (int i = 0; i < n; i ++){
        int x1,y1,x2,y2;
        pt1 = pts1[i];
        pt2 = pts2[i];

        x1 = pt1[0];
        y1 = pt1[1];
        x2 = pt2[0];
        y2 = pt2[1];

        // B is determined from the coordinates of the points in the second image
        B[2*i] = x2;
        B[2*i+1] = y2;

        // A is constructed based on the coordinates of the corresponding points
        A(2*i, 0) = x1;
        A(2*i, 1) = y1;
        A(2*i, 2) = 1;
        A(2*i, 3) = 0;
        A(2*i, 4) = 0;
        A(2*i, 5) = 0;
        A(2*i, 6) = -x1*x2;
        A(2*i, 7) = -x2*y1;

        A(2*i+1, 0) = 0;
        A(2*i+1, 1) = 0;
        A(2*i+1, 2) = 0;
        A(2*i+1, 3) = x1;
        A(2*i+1, 4) = y1;
        A(2*i+1, 5) = 1;
        A(2*i+1, 6) = -x1*y2;
        A(2*i+1, 7) = -y2*y1;
    }

    // ------ FIN TO DO ------

    B = linSolve(A, B); // Solve the linear system
    // Define the homography matrix
    Matrix<float> H(3, 3);
    H(0,0)=B[0]; H(0,1)=B[1]; H(0,2)=B[2];
    H(1,0)=B[3]; H(1,1)=B[4]; H(1,2)=B[5];
    H(2,0)=B[6]; H(2,1)=B[7]; H(2,2)=1;

    // Sanity check
    for(size_t i=0; i<n; i++) {
        float v1[]={(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[]={(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1,3);
        Vector<float> x2(v2,3);
        x1 = H*x1;
        cout << x1[1]*x2[2]-x1[2]*x2[1] << ' '
             << x1[2]*x2[0]-x1[0]*x2[2] << ' '
             << x1[0]*x2[1]-x1[1]*x2[0] << endl;
    }
    return H;
}

//// ----------------------------------------------------------------////
// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void growTo(float& x0, float& y0, float& x1, float& y1, float x, float y) {
    if(x<x0) x0=x;
    if(x>x1) x1=x;
    if(y<y0) y0=y;
    if(y>y1) y1=y;
}

//// ----------------------------------------------------------------////
// Panorama construction
void panorama(const Image<Color,2>& I1, const Image<Color,2>& I2,
              Matrix<float> H) {
    Vector<float> v(3);
    float x0=0, y0=0, x1=I2.width(), y1=I2.height();

    v[0]=0; v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=0; v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1<<endl;
    Image<Color> I(int(x1-x0), int(y1-y0));
    setActiveWindow( openWindow(I.width(), I.height()) );
    I.fill(WHITE);

    // ------------- TODO/A completer ----------
    Color col1,col2;
    int I2H = I2.height();
    int I1H = I1.height();
    int I1W = I1.width();
    int I2W = I2.width();
    int xim,yim;
    // Invert the matrix H
    Matrix<float> Hin = inverse(H);
    for (xim = int(x0); xim < int(x1); xim ++){
        for (yim = int(y0); yim < int(y1); yim++){

            if (( xim>=0 && xim < I2W) && (yim>=0 && yim < I2H)){
                col2 = I2(xim,yim);
                v[0] = xim;
                v[1] = yim;
                v[2] = 1;
                v = Hin*v;
                v/=v[2];
                int xfus = int(v[0]);
                int yfus = int(v[1]);
                if ((xfus>=0 && xfus < I1W) && (yfus>=0 && yfus < I1H)) {
                    col1 = I1(xfus, yfus);
                    I(xim-x0,yim-y0) = col1/byte(2)+col2/byte(2);
                } else{
                    I(xim-x0,yim-y0) = col2;
                }
            } else{
                v[0] = xim;
                v[1] = yim;
                v[2] = 1;
                v = Hin*v;
                v/=v[2];
                int xfus = int(v[0]);
                int yfus = int(v[1]);
                if ((xfus>=0 && xfus < I1W) && (yfus>=0 && yfus < I1H)) {
                    col1 = I1(xfus, yfus);
                    I(xim-x0,yim-y0) = col1;
                }
            }
        }
    }

    // ------ END TO DO ------
    display(I,0,0);
}

//// ----------------------------------------------------------------////
// Main function
int main(int argc, char* argv[]) {
    const char* s1 = argc>1? argv[1]: srcPath("image0006.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("image0007.jpg");

    // Load and display images
    Image<Color> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1,0,0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2,0,0);

    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1="<<endl;
    for(it=pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2="<<endl;
    for(it=pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H/H(2,2);

    // Apply homography
    panorama(I1, I2, H);

    endGraphics();
    return 0;
}
