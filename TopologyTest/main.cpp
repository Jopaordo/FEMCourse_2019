
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


#include "pzrefquad.h"
TPZGeoMesh *GenerateMeshRef(double l, double h, int side, int nivel);
void TestTriangleTransforms();

using namespace std;
//typedef REAL RND_MAX;
REAL RND_MAX = 4294967295.0;


int main(){

    TPZVec<double> vectorin(3,0.1);
//  TestTriangleTransforms();
//  TestQuadrilateralTransforms();
    TestTetrahedronTransforms(vectorin);

    GenerateMeshRef(2,2, 4, 4);
    
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




TPZGeoMesh *GenerateMeshRef(double l, double h, int side, int nivel){
    TPZVec<REAL> p0(3,0.0);
    TPZVec<REAL> p1(3,0.0);
    p1[0]= l;
    p1[1]= h;
    TPZVec<int> nx(2,1);
    
    TPZGenGrid gen(nx, p0, p1);
    TPZGeoMesh *gmsh = new TPZGeoMesh;
    gen.Read(gmsh);
    gmsh->BuildConnectivity();
    double val = h;
    int count = 0;
    if (side == 4) {
        for (count=0; count<nivel; count++) {
            val= val/2.0;
            int nels = gmsh->NElements();
            for (int j=0; j<nels; j++) {
                TPZGeoEl *gel = gmsh->Element(j);
                if (gel->HasSubElement()) {
                    continue;
                }
                TPZFMatrix<REAL> cooridnates;
                gel->NodesCoordinates(cooridnates);
                cooridnates.Print(std::cout);
                if (cooridnates(1,0)< val) {
                    TPZVec<TPZGeoEl *> filios;
                    gel->Divide(filios);
                }
                
//                TPZVec<TPZGeoEl *> filhos;
//                gel->Divide(filhos);
            }
        }
   
    }
 
//    for (int i=0; i<nivel; i++) {
//        int nels = gmsh->NElements();
//        for (int j=0; j<nels; j++) {
//            TPZGeoEl *gel = gmsh->Element(j);
//            if (gel->HasSubElement()) {
//                continue;
//            }
//            TPZVec<TPZGeoEl *> filhos;
//            gel->Divide(filhos);
//        }
//
//    }
    
    std::cout<<gmsh->NElements()<<endl;
    std::ofstream file("2dGeomesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmsh, file);
    return gmsh;
    
}

void Divide(TPZGeoEl *gel, int side, double val){
    TPZFMatrix<REAL> cooridnates;
    TPZGeoMesh *gmesh = gel->Mesh();
    gel->NodesCoordinates(cooridnates);
    TPZVec<int64_t> nodeindices;
    gel->GetNodeIndices(nodeindices);
    int64_t index1;
    pzrefine::TPZRefQuad::NewMidSideNode(gel, 5, index1);
    int64_t index2;
    pzrefine::TPZRefQuad::NewMidSideNode(gel, 7, index2);
    
    
    
//    gel->MidSideNodeIndex(5, )
//    TPZGeoEl *subel = gel->Mesh()->CreateGeoElement(EQuadrilateral,cornerindexes,matid,index,0);
//    tpzrefquad
//    gel->SetSubElement(i , subel);
}

TPZGeoNode *Geonode(TPZGeoMesh *gmsh, int side){
    
//    MidSideNodeIndex(gel,side,index);
//    if(index < 0) {
//        TPZGeoElSide gelside = gel->Neighbour(side);
//        if(gelside.Element()) {
//            while(gelside.Element() != gel) {
//                gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
//                if(index!=-1) return;
//                gelside = gelside.Neighbour();
//            }
//        }
//        TPZVec<REAL> par(3,0.);
//        TPZVec<REAL> coord(3,0.);
//        if(side < TPZShapeQuad::NCornerNodes) {
//            index = gel->NodeIndex(side);
//            return;
//        }
//        //aqui side = 8 a 26
//        side-=TPZShapeQuad::NCornerNodes;//0,1,..,18
//        par[0] = MidCoord[side][0];
//        par[1] = MidCoord[side][1];
//        gel->X(par,coord);
//        index = gel->Mesh()->NodeVec().AllocateNewElement();
//        gel->Mesh()->NodeVec()[index].Initialize(coord,*gel->Mesh());
    
    return 0;
};

