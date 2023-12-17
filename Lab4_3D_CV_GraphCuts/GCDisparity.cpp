// Imagine++ project
// Project:  GraphCutsDisparity
// Author:   Renaud Marlet/Pascal Monasse
// Student: MARENGO Matteo
// Email: matteo.marengo@ens-paris-saclay.fr
// Date : 07/11/2023

#include <Imagine/Images.h>
#include <iostream>
#include <algorithm>
#include <string>
#include "maxflow/graph.h"
using namespace Imagine;
using namespace std;

typedef Image<byte> byteImage;
typedef Image<double> doubleImage;

// Default data
const char *DEF_im1=srcPath("face0.png"), *DEF_im2=srcPath("face1.png");
static int dmin=10, dmax=55; // Min and max disparities


// Parameters of the algorithm
// OPTIMIZATION: to make the program faster, a zoom factor is used to
// down-sample the input images on the fly. You will
// only look at pixels (win+zoom*i,win+zoom*j) with win the radius of patch.
const int win = (7-1)/2;    // Correlation patches of size (2n+1)*(2n+1)
const float lambdaf = 0.1;  // Weight of regularization (smoothing) term
const int zoom = 2;         // Zoom factor (to speedup computations)
const float sigma = 6/zoom; // Gaussian blur parameter for disparity
// Energy discretization precision (as we build a graph with 'int' weights)
const int wcc = std::max(1+int(1/lambdaf),20);
const int lambda = lambdaf*wcc; // Regularization term (must be >= 1)

// ################
// ### SCRIPT 1 ###
// ################
#ifdef IMAGINE_OPENGL
void doc() {
    cout << "***** 3D mesh renderings *****" << endl;
    cout << "- Button 1: toggle textured or gray rendering" << endl;
    cout << "- SHIFT+Button 1: rotate" << endl;
    cout << "- SHIFT+Button 3: translate" << endl;
    cout << "- Mouse wheel: zoom" << endl;
    cout << "- SHIFT+a: zoom out" << endl;
    cout << "- SHIFT+z: zoom in" << endl;
    cout << "- SHIFT+r: recenter camera" << endl;
    cout << "- SHIFT+m: toggle solid/wire/points mode" << endl;
    cout << "- Button 3: exit" << endl;
}

// ################
// ### SCRIPT 2 ###
// Display 3D mesh renderings
void control3D(const Mesh& Mt, const Mesh& Mg) {
    setActiveWindow( openWindow3D(512, 512, "3D") );
    showMesh(Mt);
    bool textured = true;
    while (true) {
        Event evt;
        getEvent(5, evt);
        // On mouse button 1
        if (evt.type == EVT_BUT_ON && evt.button == 1) {
            // Toggle textured rendering and gray rendering
            if (textured) {
                hideMesh(Mt,false);
                showMesh(Mg,false);
            } else {
                hideMesh(Mg,false);
                showMesh(Mt,false);
            }
            textured = !textured;
        }
        // On mouse button 3
        if (evt.type == EVT_BUT_ON && evt.button == 3)
            break;
    }
}
#endif

// ################
// ### SCRIPT 3 ###
// ################
void show3D(const byteImage& I, const doubleImage& D, int zoom) {
#ifdef IMAGINE_OPENGL
    cout << "Click to compute depth map and 3D mesh renderings... " << flush;
    click();

    // Compute 3D point cloud: magic constants depending on camera pose
    const float f  = 750;  // Focal
    const float d0 = 100;  // Disparity for infinite depth (due to crop)
    const float B  = -0.20;// Baseline

    const int nx=D.width(), ny=D.height();
    FMatrix<float,3,3> K(0.0f);
    K(0,0)= -f/zoom; K(0,2)=nx/2;
    K(1,1)=  f/zoom; K(1,2)=ny/2;
    K(2,2)=1.0f;
    K = inverse(K);
    K /= K(2,2);
    Array<FloatPoint3> p(nx*ny);
    Array<Color> pcol(nx*ny);
    for(int j=0; j<ny; j++)
        for(int i=0; i<nx; i++) {
            float z = f*B/(d0+D(i,j));
            FloatPoint3 pt((float)i,(float)j,1.0f);
            p[i+nx*j] = K*pt*z;
            pcol[i+nx*j] = Color(I(i*zoom,j*zoom));
        }
    // Create mesh from 3D point cloud
    Array<Triangle> t(2*(nx-1)*(ny-1));
    Array<Color> tcol(2*(nx-1)*(ny-1));
    for(int j=0; j+1<ny; j++)
        for(int i=0; i+1<nx; i++) {
            // Create triangles with next pixels in line/column
            t[2*(i+j*(nx-1))]   = Triangle(i+nx*j,   i+1+nx*j,     i+nx*(j+1));
            t[2*(i+j*(nx-1))+1] = Triangle(i+1+nx*j, i+1+nx*(j+1), i+nx*(j+1));
            tcol[2*(i+j*(nx-1))]   = pcol[i+nx*j];
            tcol[2*(i+j*(nx-1))+1] = pcol[i+nx*j];
        }
    // Mesh with texture from original image
    Mesh Mt(p.data(), nx*ny, t.data(), 2*(nx-1)*(ny-1), 0, 0, FACE_COLOR);
    Mt.setColors(TRIANGLE, tcol.data());
    // Mesh with artificial light
    Mesh Mg(p.data(), nx*ny, t.data(), 2*(nx-1)*(ny-1), 0, 0,
            CONSTANT_COLOR, SMOOTH_SHADING);
    cout << "done" << endl;

    doc();
    control3D(Mt, Mg);
#else
    cout << "No 3D: Imagine++ not built with OpenGL support" << endl;
#endif
}

// ################
// ### SCRIPT 4 ###
// Return image of mean intensity value over (2win+1)x(2win+1) patch
doubleImage meanImage(const byteImage& I) {
    int w = I.width(), h = I.height();
    doubleImage IM(w,h);
    double area = (2*win+1)*(2*win+1);
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++) {
            if (j-win<0 || j+win>=h || i-win<0 ||i+win>=w) {
                IM(i,j)=0;
                continue;
            }
            double sum = 0;
            for(int y=-win; y<=win; y++)
                for(int x=-win; x<=win; x++)
                    sum += I(x+i,y+j);
            IM(i,j) = sum / area;
        }
    return IM;
}

// ################
// ### SCRIPT 5 ###
// Compute correlation between two pixels in images 1 and 2
double correl(const byteImage& I1,    // Image 1
              const doubleImage& I1M, // Image of mean value over patch
              const byteImage& I2,    // Image2
              const doubleImage& I2M, // Image of mean value over patch
              int u1, int v1,         // Pixel of interest in image 1
              int u2, int v2) {       // Pixel of interest in image 2
    double c=0;
    for(int y=-win; y<=win; y++)
        for(int x=-win; x<=win; x++) 
            c += (I1(u1+x,v1+y)-I1M(u1,v1)) * (I2(u2+x,v2+y)-I2M(u2,v2));
    return c / ((2*win+1)*(2*win+1));
}

// ################
// ### SCRIPT 6 ###
// Compute ZNCC between two patches in images 1 and 2
double zncc(const byteImage& I1,    // Image 1
            const doubleImage& I1M, // Image of mean intensity value over patch
            const byteImage& I2,    // Image2
            const doubleImage& I2M, // Image of mean intensity value over patch
            int u1, int v1,         // Pixel of interest in image 1
            int u2, int v2) {       // Pixel of interest in image 2
    double var1 = correl(I1, I1M, I1, I1M, u1, v1, u1, v1);
    if(var1 == 0)
        return 0;
    double var2 = correl(I2, I2M, I2, I2M, u2, v2, u2, v2);
    if(var2 == 0)
        return 0;
    return correl(I1, I1M, I2, I2M, u1, v1, u2, v2) / sqrt(var1 * var2);
}

// Define a new structure
struct FourFloats {
    float a, b, c, d;
};

// ### ZOOM COORDINATES MODIFICATION  ###
FourFloats zoom_coordinates(int x,int y,int zoom,int win, int dmin, int d){
    float zx1 = (x*zoom)+win; float zy1 = (y*zoom)+win;
    float zx2 = zx1 + dmin + d; float zy2 = zy1;
    FourFloats result {zx1,zy1,zx2,zy2};
    return result;
}

// ### CALCULATE RHO ###
// As seen in the slides
double rho(double c){
    if (c<0){
        return 1;
    }
    if (c>= 0 && c <= 1){
        return sqrt(1-c);
    }
    else{
        std::cerr<< "Error with the value of the cost ! "<< std::endl;
        return 0;
    }
}

// ################
// ### SCRIPT 7 ###
/// 1) Create graph
/// The graph library works with node numbers.
/// 2) To clarify the setting, create a formula to associate a unique node number to a triplet (x,y,d) of pixel
/// coordinates and disparity.
/// 3) The library assumes an edge consists of a pair of oriented edges, one in each direction.
/// Put correct weights to the edges, such as 0, INF, or an intermediate weight.
void build_graph(Graph<int,int,int>& G,
                 const byteImage& I1, const byteImage& I2,
                 int nx, int ny, int nd) {
    const int INF=1000000; // "Infinite" value for edge impossible to cut
    // Precompute images of mean intensity value over patch
    doubleImage I1M = meanImage(I1), I2M = meanImage(I2);
    // ------------- TODO -------------
    // nx and ny have already been modified (zoom)
    G.add_node(nx*ny*nd); // graph node creation
    // Define a large penalty for discrepancies in pairwise potential
    double K = 1 + nd*4*lambda; // As seen in Slide 89 with 4 neighbours per node
    // Calculate horizontal, vertical and terminal edges in a single loop
    for (int node = 0; node < nx * ny * nd; ++node) {
        int x = node % nx;
        int y = (node / nx) % ny;
        int d = node / (nx * ny);
        int NodeA = d*nx*ny+y*nx+x;

        // 1 - Assign weights for each of neighboring sites
        // The node's right neighbor is the next node in the same row and disparity level
        if (x + 1 < nx) {
            int rNode = NodeA + 1;
            G.add_edge(NodeA, rNode, lambda, lambda);
        }

        // The node's bottom neighbor is nx positions away in the same column and disparity level
        if (y + 1 < ny) {
            int bNode = NodeA + nx;
            G.add_edge(NodeA, bNode, lambda, lambda);
        }

        // 2 - Assign weights for terminal edges
        // Adjusting coordinates for zoom
        FourFloats ZoomC = zoom_coordinates(x,y,zoom,win,dmin,d);
        float zx1 = ZoomC.a; float zy1 = ZoomC.b; float zx2 = ZoomC.c; float zy2 = ZoomC.d;

        if (d==0){
            // SOURCE weights | d = 0
            double src_weight = K + wcc*rho(zncc(I1,I1M,I2,I2M,zx1,zy1,zx2,zy2));
            G.add_tweights(x+y*nx,src_weight,0);
        }

        if (d != nd-1){
            // 3 - Assign weights for each location
            int node1 = x + y*nx + d*nx*ny;
            int node2 = node1 + nx*ny;
            double edge_weight = K + wcc*rho(zncc(I1,I1M,I2,I2M,zx1,zy1,zx2,zy2));
            G.add_edge(node1,node2,edge_weight,INF);
        }
        else{
            // SINK weigths | d = (nd-1)
            double sink_weight = K + wcc*rho(zncc(I1,I1M,I2,I2M,zx1,zy1,zx1 + dmin + nd-1 ,zy2));
            G.add_tweights((nd-1)*nx*ny + x+y*nx ,0,sink_weight);
        }
    }
    // END TO DO

}

/// Extract disparity from minimum cut
doubleImage decode_graph(Graph<int,int,int>& G, int nx, int ny, int nd) {
    doubleImage D(nx,ny);

    // ------------- TODO -------------
    for (int node = 0; node < nx * ny; ++node){
        int x = node % nx;
        int y = (node / nx) % ny;
        float actuDisparity = -1;
        for(int d=0; d<nd; d++) {
            if(G.what_segment(d*nx*ny+y*nx+x) == Graph<int,int,int>::SINK) {
                actuDisparity = d + dmin;
                break;
            }
        }
        D(x,y) = actuDisparity;
    }
    // END TO DO
    return D;
}
    
// Load two rectified images.
// Compute the disparity of image 2 w.r.t. image 1.
// Display disparity map.
// Display 3D mesh of corresponding depth map.
int main(int argc, char* argv[]) {
    if(argc!=1 && argc!=5) {
        cerr << "Usage: " << argv[0] << " im1 im2 dmin dmax" << endl;
        return 1;
    }
    const char *im1=DEF_im1, *im2=DEF_im2;
    if(argc>1) {
        im1 = argv[1]; im2=argv[2]; dmin=stoi(argv[3]); dmax=stoi(argv[4]);
    }
    cout << "Loading images... " << flush;
    byteImage I1,I2;
    if(!load(I1, im1) || !load(I2,im2)) {
        cerr << "Error loading image files" << endl;
        return 1;
    }
    cout << "done" << endl;

    cout << "Parameters: " << "d=" << dmin << "..." << dmax
         << ", win="<<win << ", lambda="<<lambdaf << ", sigma="<<sigma
         << ", zoom="<<zoom << endl;

    cout << "Displaying images... " << flush;
    int w1=I1.width(), w2=I2.width(), h=I1.height();
    openWindow(w1+w2, h);
    display(I1); display(I2,w1,0);
    cout << "done" << endl;

    // Zoomed image dim, disregarding borders (strips of width the patch radius)
    const int nx=(w1-2*win)/zoom, ny=(h-2*win)/zoom;
    const int nd=dmax-dmin; // Disparity range

    cout << "Constructing graph (be patient)... " << flush;
    Graph<int,int,int> G(nx*ny*nd,2*nx*ny*nd);
    build_graph(G, I1, I2, nx, ny, nd);
    cout << "done" << endl;

    cout << "Computing minimum cut... " << flush;
    int f = G.maxflow();
    cout << "done" << endl << "  max flow = " << f << endl;

    cout << "Extracting disparity map from minimum cut... " << flush;
    doubleImage D=decode_graph(G, nx, ny, nd);
    cout << "done" << endl;

    cout << "Displaying disparity map... " << flush;
    fillRect(0,0,w1,h,CYAN);
    display(enlarge(grey(D),zoom),win,win);
    cout << "done" << endl;
    cout << "Click to compute and display blured disparity map... " << flush;
    click();
    D=blur(D,sigma);
    display(enlarge(grey(D),zoom),win,win);
    cout << "done" << endl;

    show3D(I1.getSubImage(win,win,w1-2*win,h-2*win), D, zoom);
    endGraphics();
    return 0;
}
