
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
#include "pzshapequad.h"

void TestTriangleTransforms(TPZVec<double> vectorin);
void TestQuadrilateralTransforms(TPZVec<double> vectorin);
void TestPyramidTransforms(TPZVec<double> vectorin);
void TestTetrahedronTransforms(TPZVec<double> vectorin);

using namespace std;
//typedef REAL RND_MAX;
REAL RND_MAX = 4294967295.0;


int main(){
    TPZVec<double> vectorin(3,0.1);
//  TestTriangleTransforms();
//  TestQuadrilateralTransforms();
    TestTetrahedronTransforms(vectorin);
    return 0;
}
void TestTriangleTransforms(TPZVec<double> vectorin){
    std::cout<<"************************ "<<std::endl;
    std::cout<<"TransformElementToSide TEST->Triangle"<<std::endl;
    std::cout<<"Parametric Domain 2D {xi,0,1}, {eta,0,1-xi} "<<std::endl;
    std::cout<<"Parametric Domain 1D {xi,-1,1}"<<std::endl;
    std::cout<<"************************ "<<std::endl;
  
    TPZVec<double> vectorout(1);
    
    std::cout<<"Point to transform: "<<vectorin<<std::endl;
    
    pztopology::TPZTriangle trian;
    int npoints = trian.NCornerNodes;
    int nfaces = trian.NFaces;
    for (int iside = npoints; iside<npoints + nfaces; iside++) {
        TPZTransform<> transform = trian.TransformElementToSide(iside);
        transform.Apply(vectorin, vectorout);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
    }
    
 
    return;
    
    
}
void TestQuadrilateralTransforms(TPZVec<double> vectorin){
    
    std::cout<<"************************ "<<std::endl;
    std::cout<<"TransformElementToSide TEST->Quadrilateral"<<std::endl;
    std::cout<<"Parametric Domain 2D {xi,-1,1}, {eta,-1,1} "<<std::endl;
    std::cout<<"Parametric Domain 1D {xi,-1,1}"<<std::endl;
    std::cout<<"************************ "<<std::endl;
  
    TPZVec<double> vectorout(1);
    
    std::cout<<"Point to transform: "<<vectorin<<std::endl;
    
    pztopology::TPZQuadrilateral trian;
    int npoints = trian.NCornerNodes;
    int nfaces = trian.NFaces;
    for (int iside = npoints; iside<npoints + nfaces; iside++) {
        TPZTransform<> transform = trian.TransformElementToSide(iside);
        transform.Apply(vectorin, vectorout);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
    }
    
    
    return;
    
}
void TestPyramidTransforms(TPZVec<double> vectorin){
    
    std::cout<<"************************ "<<std::endl;
    std::cout<<"TransformElementToSide TEST->Pyramid"<<std::endl;
    std::cout<<"Parametric Domain 3D {xi,-1,1}, {eta,-1,1}, {zeta,0,1} "<<std::endl;
    std::cout<<"Parametric Domain 2D {xi,-1,1}, {eta,-1,1} "<<std::endl;
    std::cout<<"Parametric Domain 1D {xi,-1,1}"<<std::endl;
    std::cout<<"************************ "<<std::endl;
  
    TPZVec<double> vectorout(1);
    
    std::cout<<"Point to transform: "<<vectorin<<std::endl;
    
    pztopology::TPZPyramid pyram;
    int npoints = pyram.NCornerNodes;
    int nfaces = pyram.NFaces;
    int nribs =  pyram.NSides - npoints - nfaces -1;
    for (int iside = npoints; iside<npoints + nfaces + nribs; iside++) {
        if (iside > 12){
            vectorout.resize(2);
        }
        TPZTransform<> transform = pyram.TransformElementToSide(iside);
        transform.Apply(vectorin, vectorout);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
    }
    
    
    return;
    
}
void TestTetrahedronTransforms(TPZVec<double> vectorin){
    std::cout<<"************************ "<<std::endl;
    std::cout<<"TransformElementToSide TEST->Hexahedra"<<std::endl;
    std::cout<<"Parametric Domain 3D {xi,0,1}, {eta,0,1}, {zeta,0,1} "<<std::endl;
    std::cout<<"Parametric Domain 2D {xi,-1,1}, {eta,-1,1} "<<std::endl;
    std::cout<<"Parametric Domain 1D {xi,-1,1}"<<std::endl;
    std::cout<<"************************ "<<std::endl;
    
    vectorin[1]=0.1;
    TPZVec<double> vectorout(1);
    
    std::cout<<"Point to transform: "<<vectorin<<std::endl;
    
    pztopology::TPZTetrahedron tetra;
    int npoints = tetra.NCornerNodes;
    int nfaces = tetra.NFaces;
    int nribs =  tetra.NSides - npoints - nfaces -1;
    for (int iside = npoints; iside<npoints + nfaces + nribs; iside++) {
        if (iside > 9){
            vectorout.resize(2);
        }
        TPZTransform<> transform = tetra.TransformElementToSide(iside);
        transform.Apply(vectorin, vectorout);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
    }
    
    
    return;
}


