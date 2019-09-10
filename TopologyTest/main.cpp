
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
#include "pzrefquad.h"
TPZGeoMesh *GenerateMeshRef(double l, double h, int side, int nivel);
void TestTriangleTransforms();
using namespace std;
//typedef REAL RND_MAX;
REAL RND_MAX = 4294967295.0;


int main(){
    GenerateMeshRef(2,2, 4, 4);
    
//  TestTriangleTransforms();
    
    
//    pztopology::TPZQuadrilateral trian;
//    double val = (((double) arc4random())/RND_MAX);
//    TPZTransform<> tr = trian.SideToSideTransform(8, 6);
//    TPZFMatrix<REAL> MULT = tr.Mult();
//    TPZFMatrix<REAL> SUM =tr.Sum();
//  //  TPZVec<double> vectorin(2,1/3);
//    //vectorin[1]=0.3;
//    TPZVec<double> vectorout(1);
//    tr.Apply(vectorin, vectorout);
//    std::cout<<vectorout[0]<<std::endl;
//   // std::cout<<vectorout[1]<<std::endl;
//    TPZVec<int64_t> id(4,0);
//    id[0]=4;
//    id[1]=3;
//    id[2]=2;
//    id[3]=1;
//    int valret = trian.GetTransformId(8,id);
//    TPZVec<int> permute(9,0);
//    trian.GetGatherPermute(valret, permute);
//
//    for (int i=0; i<9; i++) {
//        std::cout<<"id "<<i<<":"<<permute[i]<<std::endl;
//    }
//
//
//
//    pzshape::TPZShapeQuad quad;
//    REAL out;
//    quad.ProjectPoint2dQuadToRib(2, vectorin, out);
//
//
//
//    TPZManVector<REAL> pt(3,0.0);
//    quad.RandomPoint(pt);
//    std::cout<<pt[0]<<std::endl;
//    std::cout<<pt[1]<<std::endl;
//    std::cout<<pt[2]<<std::endl;
    return 0;
}
void TestTriangleTransforms(){
    
    
    int sidefrom=6;
    int sideto=3;
    
    TPZVec<double> vectorin(2,1.0/3.0);
   // vectorin[1]=0.3;
    
    pztopology::TPZTriangle trian;
    TPZTransform<> tr = trian.SideToSideTransform(sidefrom, sideto);
    TPZTransform<> trinv = trian.SideToSideTransform(sideto, sidefrom);

    TPZFMatrix<REAL> mult = trinv.Mult();
    TPZFMatrix<REAL> sum  = trinv.Sum();
    
    mult.Print(std::cout);
    sum.Print(std::cout);
    
    TPZVec<double> vectorout(1);
    tr.Apply(vectorin, vectorout);
    std::cout<<vectorout[0]<<std::endl;
    
    REAL val = (vectorout[0]+1.0)/2.0;
    TPZVec<double> vectoroutinv(2);
    TPZVec<double> vectorinInv(1,vectorout[0]);
    trinv.Apply(vectorinInv, vectoroutinv);
    vectoroutinv[1]=vectorin[1];
    std::cout<<vectoroutinv[0]<<std::endl;
    std::cout<<vectoroutinv[1]<<std::endl;
    
    
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
