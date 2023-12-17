// Imagine++ project
// Project:  Seeds
// Author:   Pascal Monasse
// Student: MARENGO Matteo
// Mail: matteo.marengo@ens-paris-saclay.fr
// Date: 24/10/2023

#include <Imagine/Images.h>
#include <queue>
#include <string>
#include <iostream>
using namespace Imagine;
using namespace std;

// Default data
// Images path
// For Release Mode
const char *DEF_im1=srcPath("im1.jpg"), *DEF_im2=srcPath("im2.jpg");
// Min and max disparities
static int dmin=-30, dmax=-7;

/// Min NCC for a seed
static const float nccSeed=0.95f;

/// Radius of patch for correlation
static const int win=(9-1)/2;
/// To avoid division by 0 for constant patch
static const float EPS=0.1f;

/// A seed
struct Seed {
    Seed(int x0, int y0, int d0, float ncc0)
    : x(x0), y(y0), d(d0), ncc(ncc0) {}
    int x,y, d;
    float ncc;
};

/// Order by NCC
bool operator<(const Seed& s1, const Seed& s2) {
    return (s1.ncc<s2.ncc);
}

/// 4-neighbors
static const int dx[]={+1,  0, -1,  0};
static const int dy[]={ 0, -1,  0, +1};

// SCRIPT 1 ##########################################################################
/// Display disparity map
static Image<Color> displayDisp(const Image<int>& disp, Window W, int subW) {
    Image<Color> im(disp.width(), disp.height());
    for(int j=0; j<disp.height(); j++)
        for(int i=0; i<disp.width(); i++) {
            // If the disparity value is invalid
            if(disp(i,j)<dmin || disp(i,j)>dmax)
                im(i,j) = CYAN;
            // If the disparity value is valid
            else {
                int g = 255*(disp(i,j)-dmin)/(dmax-dmin);
                im(i,j)= Color(g,g,g);
            }
        }
    setActiveWindow(W,subW);
    display(im);
    showWindow(W,subW);
    return im;
}

// SCRIPT 2##########################################################################
/// Show 3D window
static void show3D(const Image<Color>& im, const Image<int>& disp) {
#ifdef IMAGINE_OPENGL // Imagine++ must have been built with OpenGL support...
    // Intrinsic parameters given by Middlebury website
    const float f=3740;
    const float d0=-200; // Doll images cropped by this amount
    const float zoom=2; // Half-size images, should double measured disparity
    const float B=0.160; // Baseline in m
    FMatrix<float,3,3> K(0.0f);
    K(0,0)=-f/zoom; K(0,2)=disp.width()/2;
    K(1,1)= f/zoom; K(1,2)=disp.height()/2;
    K(2,2)=1.0f;
    K = inverse(K);
    K /= K(2,2);
    std::vector<FloatPoint3> pts;
    std::vector<Color> col;
    for(int j=0; j<disp.height(); j++)
        for(int i=0; i<disp.width(); i++)
            if(dmin<=disp(i,j) && disp(i,j)<=dmax) {
                float z = B*f/(zoom*disp(i,j)+d0);
                FloatPoint3 pt((float)i,(float)j,1.0f);
                pts.push_back(K*pt*z);
                col.push_back(im(i,j));
            }
    Mesh mesh(&pts[0], pts.size(), 0,0,0,0,VERTEX_COLOR);
    mesh.setColors(VERTEX, &col[0]);
    Window W = openWindow3D(512,512,"3D");
    setActiveWindow(W);
    showMesh(mesh);
#else
    std::cout << "No 3D: Imagine++ not built with OpenGL support" << std::endl;
#endif
}

// SCRIPT 3 ##########################################################################
/// Correlation between patches centered on (i1,j1) and (i2,j2). The values
/// m1 or m2 are subtracted from each pixel value.
/// m1 and m2 are the mean average pixel values of the patches im1 and im2 respectively
/// Computation of the NCC distance for a defined patch
static float correl(const Image<byte>& im1, int i1,int j1,float m1,
                    const Image<byte>& im2, int i2,int j2,float m2) {
    float dist=0.0f;
    // ------------- TODO -------------
    // #####################
    // SCRIPT IS HERE
    float numerator=0.0f, denumeratora=0.0f, denumeratorb=0.0f, denumerator=0.0f;
    // Loop through all the pixels of Left and Right Images
    for (int pixx=-win; pixx<=win;pixx++){
        for (int pixy=-win; pixy<= win; pixy++){
            denumeratora += pow(im1(i1+pixx, j1+pixy)-m1, 2);
            denumeratorb += pow(im2(i2+pixx, j2+pixy) -m2, 2);
            numerator += (im1(i1+pixx, j1+pixy)-m1) * (im2(i2+pixx, j2+pixy)-m2);
        }
    }
    denumerator = sqrt(denumeratora*denumeratorb);
    dist = numerator/denumerator;
    // #####################
    return dist;
}

// SCRIPT 4 ##########################################################################
/// Sum of pixel values in patch centered on (i,j).
/// Will be used to compute m1 and m2 (mean)
static float sum(const Image<byte>& im, int i, int j) {
    float s=0.0f;
    // ------------- TODO -------------
    // ################################
    // SCRIPT IS HERE
    // Loop through all the image
    for (int a = -win; a < win +1; a++){
        for (int b = -win; b < win + 1; b++){
            s += im(i + a, j+b);
        }
    }
    // ################################
    return s;
}

// SCRIPT 5 ##########################################################################
/// Centered correlation of patches of size 2*win+1.
/// Return NCC distance value
static float ccorrel(const Image<byte>& im1,int i1,int j1,
                     const Image<byte>& im2,int i2,int j2) {
    int w = 2*win+1; // width of the patch (including the central pixel)
    float m1 = sum(im1,i1,j1)/(w*w); // mean value of first patch
    float m2 = sum(im2,i2,j2)/(w*w); // mean value of second patch
    return correl(im1,i1,j1,m1, im2,i2,j2,m2);
}

// SCRIPT 6 ##########################################################################
/// Compute disparity map from im1 to im2, but only at points where NCC is
/// above nccSeed (0.95). Set to true the seeds and put them in Q.
static void find_seeds(Image<byte> im1, Image<byte> im2,
                       float nccSeed,
                       Image<int>& disp, Image<bool>& seeds,
                       std::priority_queue<Seed>& Q) {
    disp.fill(dmin-1); // fill the disparity map
    seeds.fill(false); // fill the seed map
    // Empty the priority queue and make it ready to store new elements
    while(! Q.empty()) // clear the priority queue
        Q.pop();

    const int maxy = std::min(im1.height(),im2.height()); // minimum height of the two images
    const int refreshStep = (maxy-2*win)*5/100;
    for(int y=win; y+win<maxy; y++) {
        if((y-win-1)/refreshStep != (y-win)/refreshStep)
            std::cout << "Seeds: " << 5*(y-win)/refreshStep <<"%\r"<<std::flush;
        for(int x=win; x+win<im1.width(); x++) {
            // ------------- TODO -------------
            // ################################
            // Hint: just ignore windows that are not fully in image
            // SCRIPT IS HERE
            int topD; // Best disparity value
            float topNCC=-5.0f; // Low value so it will be updated as lowest value is -1
            float tempNcc;
            int newx;
            int _d = dmin;
            while (_d <= dmax) {
                int newX = x + _d;
                if (newX >= win) {
                    float _ncc = ccorrel(im1, x, y, im2, newX, y);
                    if (_ncc > topNCC) {
                        topNCC = _ncc;
                        topD = _d;
                    }
                }
                _d++;
            }
            // Compute the NCC with all points on the horizontal line
            // Keep the ones where Ncc > nccSeed
            if (topNCC >= nccSeed){
                Seed seed(x, y, topD, topNCC);
                Q.push(seed);
                disp(x,y) = topD; // update the disparity map
                seeds(x,y) = true; // update the seed map
            }
            // #################################
        }
    }
    std::cout << std::endl;
}

// SCRIPT 7 ##########################################################################
/// Propagate seeds
static void propagate(Image<byte> im1, Image<byte> im2,
                      Image<int>& disp, Image<bool>& seeds,
                      std::priority_queue<Seed>& Q) {
    const int maxy = std::min(im1.height(),im2.height());
    while(! Q.empty()) { // While Queue is not empty
        Seed s=Q.top(); // gets the seed with the highest NCC Value
        Q.pop(); // Remove it
        for(int i=0; i<4; i++) { // Loop through its 4 neighbors
            int x=s.x+dx[i], y=s.y+dy[i];
            if(0<=x-win && x+win<im1.width() && 0<=y-win && y+win<maxy && ! seeds(x,y)) {
                // ################################
                // ------------- TODO -------------
                // SCRIPT IS HERE
                int topD;
                float topNCC = -5.0f; // low value so it can be updated
                float ncc;
                int newX;
                // set dQ by highest NCC score among dp-1 and dp+1
                int diff = -1;
                while (diff <= 2) {
                    newX = x + s.d + diff;
                    if (newX >= win) {
                        ncc = ccorrel(im1, x, y, im2, newX, y);
                        if (ncc > topNCC) {
                            topNCC = ncc;
                            topD = s.d + diff;
                        }
                    }
                    diff++;
                }
                if (topNCC >= -4){
                    if(topD < dmin){topD = dmin;}
                    else if (topD > dmax){topD = dmax;}
                    disp(x,y) = topD;
                    seeds(x,y) = true;
                    Seed propagateSeed(x,y,topD,topNCC);
                    Q.push(propagateSeed);
                }

                // ################################
            }
        }
    }
}

// ############################################
// MAIN #######################################
// ############################################

int main(int argc, char* argv[]) {
    if(argc!=1 && argc!=5) {
        cerr << "Usage: " << argv[0] << " im1 im2 dmin dmax" << endl;
        return 1;
    }
    const char *im1=DEF_im1, *im2=DEF_im2;
    if(argc>1) {
        im1 = argv[1]; im2=argv[2]; dmin=stoi(argv[3]); dmax=stoi(argv[4]);
    }
    // Load and display images
    Image<Color> I1, I2;
    if(!load(I1,im1) || !load(I2,im2)) {
        cerr<< "Error loading image files" << endl;
        return 1;
    }
    std::string names[5]={"image 1","image 2","dense","seeds","propagation"};
    Window W = openComplexWindow(I1.width(), I1.height(), "Seeds propagation",
                                 5, names);
    setActiveWindow(W,0);
    display(I1,0,0);
    setActiveWindow(W,1);
    display(I2,0,0);

    Image<int> disp(I1.width(), I1.height());
    Image<bool> seeds(I1.width(), I1.height());
    std::priority_queue<Seed> Q;

    // Dense disparity
    find_seeds(I1, I2, -1.0f, disp, seeds, Q);
    save(displayDisp(disp,W,2), srcPath("0dense.png"));

    // Only seeds
    find_seeds(I1, I2, nccSeed, disp, seeds, Q);
    save(displayDisp(disp,W,3), srcPath("1seeds.png"));

    // Propagation of seeds
    propagate(I1, I2, disp, seeds, Q);
    save(displayDisp(disp,W,4), srcPath("2final.png"));

    // Show 3D (use shift click to animate)
    show3D(I1,disp);

    endGraphics();
    return 0;
}
