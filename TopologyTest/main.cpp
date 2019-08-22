
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"

#include <TPZRefPattern.h>

#include "TPZMaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"

#include "pzgengrid.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"
#include "pzshapetriang.h"
using namespace std;
//typedef REAL RND_MAX;
REAL RND_MAX = 4294967295.0;


int main(){
    pztopology::TPZTriangle trian;
    double val = (((double) arc4random())/RND_MAX);
    TPZTransform<> tr = trian.SideToSideTransform(6, 3);
    TPZFMatrix<REAL> MULT = tr.Mult();
    TPZFMatrix<REAL> SUM =tr.Sum();
    TPZVec<double> vectorin(2,1.0/10.0);
 //   vectorin[1]=1.0/5.0;
    TPZVec<double> vectorout(1);
    tr.Apply(vectorin, vectorout);
    std::cout<<vectorout[0]<<std::endl;
   // std::cout<<vectorout[1]<<std::endl;
   
    pzshape::TPZShapeTriang triang;
    REAL out;
    triang.ProjectPoint2dTriangToRib(0, vectorin, out);
    
    
//    TPZManVector<REAL> pt(3,0.0);
//    quad.RandomPoint(pt);
//    std::cout<<pt[0]<<std::endl;
//    std::cout<<pt[1]<<std::endl;
//    std::cout<<pt[2]<<std::endl;
    return 0;
}
