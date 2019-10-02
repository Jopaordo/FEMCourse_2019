
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
#include "TPZRefPatternDataBase.h"
void ShapeTriangle();
void HexahedronTest();
int main(){
   // gRefDBase.InitializeUniformRefPattern(EPrisma);
    gRefDBase.InitializeAllUniformRefPatterns();
 
  //  HexahedronTest();
    
    
    return 0;
}
void ShapeTriangle(){
    pzshape::TPZShapeTriang triangle;
    TPZVec<REAL> pt(2,0.3);
    TPZVec<int64_t> id(3,0);
    id[1]=1;
    id[2]=2;
    TPZVec<int> order(4,2);
    TPZFMatrix<REAL> phi(3,1);
    TPZFMatrix<REAL> dphi(2,3);
    triangle.Shape(pt, id, order, phi, dphi);
    phi.Print(std::cout);
}
void HexahedronTest(){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(1);
    TPZVec<REAL> nodeCoord(3);
    nodeCoord[0] =  0.;
    nodeCoord[1] =  0.;
    nodeCoord[2] = -1.;
    TPZGeoNode node(0, nodeCoord, *gmesh);
    gmesh->NodeVec()[0]=node;
    
    nodeCoord[0] =  1.;
    nodeCoord[1] =  0.;
    nodeCoord[2] = -1.;
    gmesh->NodeVec().Resize(2);
    TPZGeoNode node1(1, nodeCoord, *gmesh);
    gmesh->NodeVec()[1]=node1;
    
    nodeCoord[0] =  0.;
    nodeCoord[1] =  1.;
    nodeCoord[2] = -1.;
gmesh->NodeVec().Resize(3);
    TPZGeoNode node2(2, nodeCoord, *gmesh);
    gmesh->NodeVec()[2]=node2;
    nodeCoord[0] = 0.;
    nodeCoord[1] = 0.;
    nodeCoord[2] = 1.;
    gmesh->NodeVec().Resize(4);
    TPZGeoNode node3(3, nodeCoord, *gmesh);
    gmesh->NodeVec()[3]=node3;
    
    nodeCoord[0] = 0.;
    nodeCoord[1] = 0.;
    nodeCoord[2] = 1.;
    gmesh->NodeVec().Resize(5);
    TPZGeoNode node4(4, nodeCoord, *gmesh);
    gmesh->NodeVec()[4]=node4;
    
    nodeCoord[0] = 0.;
    nodeCoord[1] = 1.;
    nodeCoord[2] = 1.;
    gmesh->NodeVec().Resize(6);
    TPZGeoNode node5(5, nodeCoord, *gmesh);
    gmesh->NodeVec()[5]=node5;
    
//    nodeCoord[0] = 1.;
//    nodeCoord[1] = 1.;
//    nodeCoord[2] = 1.;
//
//    gmesh->NodeVec().Resize(7);
//    TPZGeoNode node6(6, nodeCoord, *gmesh);
//    gmesh->NodeVec()[6]=node6;
//    nodeCoord[0] = -1.;
//    nodeCoord[1] =  1.;
//    nodeCoord[2] =  1.;
//
//    TPZGeoNode node7(7, nodeCoord, *gmesh);
//      gmesh->NodeVec().Resize(8);
//    gmesh->NodeVec()[7]=node7;
    
    TPZVec<int64_t> cornerindexes(6);
    cornerindexes[0]=0;
    cornerindexes[1]=1;
    cornerindexes[2]=2;
    cornerindexes[3]=3;
    cornerindexes[4]=4;
    cornerindexes[5]=5;
//    cornerindexes[6]=6;
//    cornerindexes[7]=7;
    
   // cornerindexes[3]=3;
    int64_t index=0;
    gmesh->CreateGeoElement(EPrisma, cornerindexes, 1, index);
    
    TPZVec<TPZGeoEl *> sons;
    gmesh->BuildConnectivity();
    gmesh->Element(0)->Divide(sons);
    
    std::ofstream file("tetra.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
}
