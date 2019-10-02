
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
#include "pzshapelinear.h"
#include "pzshapecube.h"
#include "pzshapepiram.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"

void TestTriangleTransforms(TPZVec<double> vectorin);
void TestQuadrilateralTransforms(TPZVec<double> vectorin);
void TestPyramidTransforms(TPZVec<double> vectorin);
void TestTetrahedronTransforms(TPZVec<double> vectorin);
void TestCubeTransforms(TPZVec<double> vectorin);
void TestPrismaTransforms(TPZVec<double> vectorin);

#include "pzrefquad.h"
TPZGeoMesh *GenerateMeshRef(double l, double h, int side, int nivel);
void TestTriangleTransforms();

using namespace std;
//typedef REAL RND_MAX;
// REAL RND_MAX = 4294967295.0;

template <class TSHAPE>
TPZTransform<REAL> GetSideTransform(int side, int trans_id) {
    
    MElementType type_side = TSHAPE::Type(side);
    TPZTransform<REAL> TransElToSide = TSHAPE::TransformElementToSide(side);
    
    TPZTransform<REAL> TransParametric(1,1);
    switch (type_side) {
        case EOned:
        {
            TransParametric = pzshape::TPZShapeLinear::ParametricTransform(trans_id);
        }
            break;
        case EQuadrilateral:
        {
            TransParametric = pzshape::TPZShapeQuad::ParametricTransform(trans_id);
        }
            break;
        case ETriangle:
        {
            TransParametric = pzshape::TPZShapeTriang::ParametricTransform(trans_id);
        }
            break;
        default:
            TPZFMatrix<REAL> Ident(TSHAPE::SideDimension(side),TSHAPE::SideDimension(side));
            break;
    }
    
    TPZFMatrix<double> resul_mult;
    
    if (side == TSHAPE::NSides - 1) {
        return TransElToSide;
    }
    else{
        TransParametric.Mult().Multiply(TransElToSide.Mult(), resul_mult);
        TransElToSide.Mult() = resul_mult;
    }
    
    // TransElToSide.Mult().Print(std::cout);
    
    return TransElToSide;
    
}
template <class TSHAPE>
void ComputeTransforms(TPZVec<int64_t> &id,  TPZVec<TPZTransform<REAL> > &transvec) {
    // int dim = TSHAPE::Dimension;
    int NSides = TSHAPE::NSides;
    int NCorners = TSHAPE::NCornerNodes;
    transvec.resize(NSides - NCorners);
    for (int iside = NCorners; iside< NSides ; iside++) {
        int pos = iside - NCorners;
        int trans_id = TSHAPE::GetTransformId(iside, id); // Foi criado
        std::cout<<"SIDE= "<<iside<<" TransID= "<<trans_id<<std::endl;
        transvec[pos] = GetSideTransform<TSHAPE>(iside, trans_id); // Foi criado
    }
    
}

template <class TSHAPE>
void Shape(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> phi,TPZFMatrix<REAL> dphi) {
    
    int dim = TSHAPE::Dimension;
    int NSides = TSHAPE::NSides;
    int NCorners = TSHAPE::NCornerNodes;
    
    TPZFNMatrix<100> phiblend(NSides,1),dphiblend(dim,NSides);
    
    TSHAPE::ShapeCorner(pt,phi,dphi);
    
    phiblend = phi;
    phiblend.Resize(NSides, 1);
    dphiblend = dphi;
    dphiblend.Resize(dim, NSides);
    
    TSHAPE::ShapeGenerating(pt,phiblend, dphiblend);
    int shape = NCorners;
    for (int side = NCorners; side<NSides ; side++)
    {
        TPZTransform<REAL> transform = transvec[side - NCorners];
        //  int numshape = orders[side - NCorners]-1;
        int numshape =TSHAPE::NConnectShapeF( side, orders[side - NCorners]);
        TPZFNMatrix<20,REAL> phin(numshape,1),dphin(TSHAPE::Dimension,numshape);
        
        int sidedim = TSHAPE::SideDimension(side);
        TPZManVector<REAL,1> outvec(sidedim);
        //  transform.Mult().Print(std::cout);
        transform.Apply(pt, outvec);
        
        
        dphin.Zero();
        TSHAPE::ShapeInternal(side, outvec,orders[side - NCorners], phin, dphin);
        
        int ndphin = phin.Rows();
        TPZFMatrix<REAL> dphiaux(dphin);
        dphiaux.Zero();
        transform.Mult().Print(std::cout);
        for(int ish=0; ish<ndphin; ish++)
        {
            for (int d=0; d<dim; d++) {
                REAL val = 0.;
                for (int sd=0; sd<sidedim; sd++) {
                    val+= transform.Mult()(sd,d)*dphin(sd,ish);
                }
                dphiaux(d,ish) = val;
            }
        }
        
        dphin.Print("dphiN= ",std::cout);
        dphiaux.Print("dphiAux= ",std::cout);
        //        std::cout<<"Side: "<<side<<std::endl;
        //        std::cout<<"Rows: "<<transform.Mult().Rows()<<" Cols: "<<transform.Mult().Cols()<<std::endl;
        for (int i = 0; i < numshape; i++) {
            phi(shape,0) = phiblend(side,0)*phin(i,0);
            for(int xj=0;xj<TSHAPE::Dimension;xj++) {
                dphi(xj,shape) = dphiblend(xj,side)*phin(i,0)+phiblend(side,0)*dphiaux(xj,i);
            }
            shape++;
        }
    }
    
   
    
    
    
    dphi.Print("phi=",std::cout,EMathematicaInput);
    //    TSHAPE::ShapeGeneratingInternal();
    
    
}

int main(){
    
    TPZVec<REAL> pt(2, 0.2);
    pt[1]=0.1;
    TPZVec<int64_t> id(4,0);
    for (int i=0; i<4; i++) {
        id[i]=i;
    }
    //    id[0]=4;
    //    id[1]=40;
    //    id[2]=20;
    //    id[3]=10;
    //    id[4]=100;
    //    id[5]=31;
    
    TPZVec<int> order(5,4);
    TPZFMatrix<REAL> phi(25,1,0.0);
    TPZFMatrix<REAL> dphi(2,25,0.0);
    
    TPZVec<TPZTransform<>> transvec;
    ComputeTransforms<pzshape::TPZShapeQuad>(id, transvec);
    Shape<pzshape::TPZShapeQuad>(pt, order, transvec, phi, dphi);
    //    pzshape::TPZShapeTriang trian;
    //    TPZVec<REAL> pt(2, 0.2);
    //    pt[1]=0.1;
    //    TPZVec<int64_t> id(3,0);
    //    id[1]=1;
    //    id[2]=2;
    //    TPZVec<int> order(4,3);
    //    TPZFMatrix<REAL> phi(10,1);
    //    TPZFMatrix<REAL> dphi(2,10);
    //    trian.Shape(pt, id, order, phi, dphi);
    //    phi.Print(std::cout);
    //    int aka=0;
    //   TestQuadrilateralTransforms(in);
    
    //  pzshape::TPZShapeCube quad;
    //    TPZVec<REAL> in(3,1.0/3.0);
    ////    in[1] = -0.45;
    //    TPZVec<int64_t> id(5,0);
    //    id[1]=1;
    //    id[2]=2;
    //    id[3]=3;
    //    id[4]=4;
    ////    id[5]=5;
    ////    id[6]=6;
    ////    id[7]=7;
    //
    //    TPZVec<TPZTransform<REAL> > transvec;
    //    ComputeTransforms<pzshape::TPZShapePiram>(id, transvec);
    //    TPZVec<double> vectorout(1);
    //    std::cout<<"N transfor: "<<transvec.size()<<std::endl;
    //    transvec[5].Apply(in, vectorout);
    //    std::cout<<"TRANSFORM: "<<vectorout[0]<<std::endl;
    //
    //
    //    return 0;
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
    static REAL psiTriangle[3]{1,1,0};
    static REAL etaTriangle[3]{0,1,1};
    int counterror=0;
    for (int iside = npoints; iside<npoints + nfaces; iside++) {
        TPZTransform<> transform = trian.TransformElementToSide(iside);
        TPZTransform<> transformInv = trian.TransformSideToElement(iside);
        transform.Apply(vectorin, vectorout);
        REAL valximapReal =0.0;
        REAL valetamapReal =0.0;
        switch (iside) {
            case 3:
                valximapReal =vectorin[0]+vectorin[1]/2;
                valetamapReal = 0.0;
                break;
            case 4:
                valximapReal=(vectorin[0] - vectorin[1] + 1)/2;
                valetamapReal = -0.5*(-vectorin[0] + vectorin[1] + 1) + 1;
                break;
            case 5:
                valximapReal = 0.0;
                valetamapReal = 1.0 - (1 - vectorin[1] - (vectorin[0])/2);
                break;
                
        }
        
        TPZVec<double> vectorout_inv(2);
        transformInv.Apply(vectorout, vectorout_inv);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
        std::cout<<"TransformSideToElement("<<iside<<") = "<<vectorout_inv<<std::endl;
        REAL xierror =vectorout_inv[0] - valximapReal;
        REAL etaerror =vectorout_inv[1] - valetamapReal;
        REAL tol = 1.0e-9;
        if (xierror<tol && etaerror<tol) {
            std::cout<<"Transformation is ok!..."<<std::endl;
            std::cout<<"************************ "<<std::endl;
        }
        else{
            std::cout<<"Error!..."<<std::endl;
            counterror++;
        }
    }
    
    std::cout<<"Errors number: "<<counterror<<std::endl;
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
    
    pztopology::TPZQuadrilateral quad;
    int npoints = quad.NCornerNodes;
    int nfaces = quad.NFaces;
    int counterror=0;
    for (int iside = npoints; iside<npoints + nfaces; iside++) {
        TPZTransform<> transform = quad.TransformElementToSide(iside);
        TPZTransform<> transformInv = quad.TransformSideToElement(iside);
        transform.Apply(vectorin, vectorout);
        REAL valximapReal =0.0;
        REAL valetamapReal =0.0;
        switch (iside) {
            case 4:
                valximapReal =vectorin[0];
                valetamapReal = -1.0;
                break;
            case 5:
                valximapReal=1.0;
                valetamapReal = vectorin[1] ;
                break;
            case 6:
                valximapReal = vectorin[0];
                valetamapReal = 1.0 ;
                break;
            case 7:
                valximapReal = -1.0;
                valetamapReal = vectorin[1];
                break;
                
        }
        
        TPZVec<double> vectorout_inv(2);
        transformInv.Apply(vectorout, vectorout_inv);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
        std::cout<<"TransformSideToElement("<<iside<<") = "<<vectorout_inv<<std::endl;
        REAL xierror =vectorout_inv[0] - valximapReal;
        REAL etaerror =vectorout_inv[1] - valetamapReal;
        REAL tol = 1.0e-9;
        if (xierror<tol && etaerror<tol) {
            std::cout<<"Transformation is ok!..."<<std::endl;
            std::cout<<"************************ "<<std::endl;
        }
        else{
            std::cout<<"Error!..."<<std::endl;
            counterror++;
        }
    }
    
    std::cout<<"Errors number: "<<counterror<<std::endl;
    return;
    
}

void TestCubeTransforms(TPZVec<double> vectorin){
    
    std::cout<<"************************ "<<std::endl;
    std::cout<<"TransformElementToSide TEST->Cube"<<std::endl;
    std::cout<<"Parametric Domain 2D {xi,-1,1}, {eta,-1,1} "<<std::endl;
    std::cout<<"Parametric Domain 1D {xi,-1,1}"<<std::endl;
    std::cout<<"************************ "<<std::endl;
    
    TPZVec<double> vectorout(1);
    
    std::cout<<"Point to transform: "<<vectorin<<std::endl;
    int counterror=0;
    pztopology::TPZCube cube;
    int npoints = cube.NCornerNodes;
    int nfaces = cube.NFaces;
    int nribs =  cube.NSides - npoints - nfaces -1;
    TPZStack<int> smallsides;
    for (int iside = npoints; iside<npoints + nfaces + nribs; iside++) {
        if (iside>=20) {
            vectorout.resize(2);
            
            
        }
        TPZTransform<> transform = cube.TransformElementToSide(iside);
        TPZTransform<> transformInv = cube.TransformSideToElement(iside);
        
        
        
        transform.Apply(vectorin, vectorout);
        REAL valximapReal =0.0;
        REAL valetamapReal =0.0;
        REAL valzetamapReal =0.0;
        switch (iside) {
            case 8:
                valximapReal =vectorin[0];
                valetamapReal = -1.0;
                valzetamapReal = -1.0;
                break;
            case 9:
                valximapReal =  1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = -1.0;
                break;
            case 10:
                valximapReal =vectorin[0];
                valetamapReal =  1.0;
                valzetamapReal = -1.0;
                break;
            case 11:
                valximapReal =  -1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = -1.0;
                break;
            case 12:
                valximapReal =  -1.0;
                valetamapReal = -1.0;
                valzetamapReal = vectorin[2];
                break;
            case 13:
                valximapReal =  1.0;
                valetamapReal = -1.0;
                valzetamapReal = vectorin[2];
                break;
            case 14:
                valximapReal =  1.0;
                valetamapReal = 1.0;
                valzetamapReal = vectorin[2];
                break;
            case 15:
                valximapReal =  -1.0;
                valetamapReal =  1.0;
                valzetamapReal = vectorin[2];
                break;
            case 16:
                valximapReal =vectorin[0];
                valetamapReal = -1.0;
                valzetamapReal = 1.0;
                break;
            case 17:
                valximapReal =  1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = 1.0;
                break;
            case 18:
                valximapReal =vectorin[0];
                valetamapReal =  1.0;
                valzetamapReal = 1.0;
                break;
            case 19:
                valximapReal =  -1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = 1.0;
                break;
            case 20:
                valximapReal =  vectorin[0];
                valetamapReal = vectorin[1];
                valzetamapReal = -1.0;
                break;
            case 21:
                valximapReal =  vectorin[0];
                valetamapReal = -1.0;
                valzetamapReal = vectorin[2];
                break;
            case 22:
                valximapReal =  1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = vectorin[2];
                break;
            case 23:
                valximapReal =  vectorin[0];
                valetamapReal = 1.0;
                valzetamapReal = vectorin[2];
                break;
            case 24:
                valximapReal =  -1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = vectorin[2];
                break;
            case 25:
                valximapReal =  vectorin[0];
                valetamapReal = vectorin[1];
                valzetamapReal = 1.0;
                break;
                
        }
        
        TPZVec<double> vectorout_inv(3);
        transformInv.Apply(vectorout, vectorout_inv);
        
        
        //
        
        int nsmallsides = smallsides.size();
        TPZTransform<> transformRibToSide;
        TPZTransform<> transformElementToRib;
        
        
        
        
        
        
        //
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
        std::cout<<"TransformSideToElement("<<iside<<") = "<<vectorout_inv<<std::endl;
        REAL xierror =abs(vectorout_inv[0] - valximapReal);
        REAL etaerror =abs(vectorout_inv[1] - valetamapReal);
        REAL zetaerror =abs(vectorout_inv[2] - valzetamapReal);
        REAL tol = 1.0e-6;
        if (xierror<tol && etaerror<tol && zetaerror<tol) {
            std::cout<<"Transformation is ok!..."<<std::endl;
            std::cout<<"************************ "<<std::endl;
        }
        else{
            std::cout<<"Error!..."<<std::endl;
            counterror++;
        }
    }
    std::cout<<"Errors number: "<<counterror<<std::endl;
    
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
    int counterror =0;
    pztopology::TPZPyramid pyramid;
    int npoints = pyramid.NCornerNodes;
    int nfaces = pyramid.NFaces;
    int nribs =  pyramid.NSides - npoints - nfaces -1;
    for (int iside = npoints; iside<npoints + nfaces + nribs; iside++) {
        if (iside>=13) {
            vectorout.resize(2);
        }
        TPZTransform<> transform = pyramid.TransformElementToSide(iside);
        TPZTransform<> transformInv = pyramid.TransformSideToElement(iside);
        transform.Apply(vectorin, vectorout);
        REAL valximapReal =0.0;
        REAL valetamapReal =0.0;
        REAL valzetamapReal =0.0;
        switch (iside) {
            case 5:
                valximapReal = vectorin[0];
                valetamapReal = -1.0;
                valzetamapReal = 0.0;
                break;
            case 6:
                valximapReal =  1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = 0.0;
                break;
            case 7:
                valximapReal = vectorin[0];
                valetamapReal = 1.0;
                valzetamapReal = 0.0;
                break;
            case 8:
                valximapReal =  -1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = 0.0;
                break;
            case 9:
                valximapReal =  vectorin[2]+0.5;
                valetamapReal = 0.0 ;
                valzetamapReal =-0.5*vectorin[0] +0.5*vectorin[2] + 0.5 ;
                break;
            case 100:
                valximapReal = 0.0;
                valetamapReal =  0.5*vectorin[1]-0.5*vectorin[2]+0.5;
                valzetamapReal = -0.5*vectorin[1]+0.5*vectorin[2]+0.5;
                break;
            case 10:
                valximapReal =  vectorin[0];
                valetamapReal = vectorin[1];
                valzetamapReal = 0.0;
                break;
            case 11:
                valximapReal =  vectorin[0];
                valetamapReal =0.0;
                valzetamapReal = vectorin[2];
                break;
            case 12:
                valximapReal =  (1.0/3.0)*(-vectorin[0]+2*vectorin[1]-vectorin[2]+1);
                valetamapReal =  (1.0/3.0)*(-vectorin[0]-vectorin[1]+2*vectorin[2]+1);
                valzetamapReal = 1.0 - valximapReal -valetamapReal;
                break;
            case 13:
                valximapReal =  0.0;
                valetamapReal =vectorin[1];
                valzetamapReal = vectorin[2];
                break;
                
        }
        
        TPZVec<double> vectorout_inv(3);
        transformInv.Apply(vectorout, vectorout_inv);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
        std::cout<<"TransformSideToElement("<<iside<<") = "<<vectorout_inv<<std::endl;
        REAL xierror =abs(vectorout_inv[0] - valximapReal);
        REAL etaerror =abs(vectorout_inv[1] - valetamapReal);
        REAL zetaerror =abs(vectorout_inv[2] - valzetamapReal);
        REAL tol = 1.0e-6;
        if (xierror<tol && etaerror<tol && zetaerror<tol) {
            std::cout<<"Transformation is ok!..."<<std::endl;
            std::cout<<"************************ "<<std::endl;
        }
        else{
            std::cout<<"Error!..."<<std::endl;
            counterror++;
        }
    }
    std::cout<<"Errors number: "<<counterror<<std::endl;
    
    return;
    
}

void TestTetrahedronTransforms(TPZVec<double> vectorin){
    std::cout<<"************************ "<<std::endl;
    std::cout<<"TransformElementToSide TEST->Hexahedra"<<std::endl;
    std::cout<<"Parametric Domain 3D {xi,0,1}, {eta,0,1}, {zeta,0,1} "<<std::endl;
    std::cout<<"Parametric Domain 2D {xi,-1,1}, {eta,-1,1} "<<std::endl;
    std::cout<<"Parametric Domain 1D {xi,-1,1}"<<std::endl;
    std::cout<<"************************ "<<std::endl;
    
    TPZVec<double> vectorout(1);
    
    std::cout<<"Point to transform: "<<vectorin<<std::endl;
    int counterror=0;
    pztopology::TPZTetrahedron tetra;
    int npoints = tetra.NCornerNodes;
    int nfaces = tetra.NFaces;
    int nribs =  tetra.NSides - npoints - nfaces -1;
    for (int iside = npoints; iside<npoints + nfaces + nribs; iside++) {
        if (iside>=10) {
            vectorout.resize(2);
        }
        TPZTransform<> transform = tetra.TransformElementToSide(iside);
        TPZTransform<> transformInv = tetra.TransformSideToElement(iside);
        transform.Apply(vectorin, vectorout);
        REAL valximapReal =0.0;
        REAL valetamapReal =0.0;
        REAL valzetamapReal =0.0;
        switch (iside) {
            case 4:
                valximapReal =vectorin[0]+0.5*vectorin[1]+0.5*vectorin[2];
                valetamapReal = 0.0;
                valzetamapReal = 0.0;
                break;
            case 5:
                valximapReal =  (vectorin[0]-vectorin[1]+1)/2.0;
                valetamapReal = 1 + (-0.5*vectorin[0]+0.5*vectorin[1]-0.5);
                valzetamapReal = 0.0;
                break;
            case 6:
                valximapReal = 0.0;
                valetamapReal = 1.0 - 0.5*vectorin[0]-vectorin[1]-0.5*vectorin[2];
                valzetamapReal = 0.0;
                break;
            case 7:
                valximapReal =  0.0;
                valetamapReal = 0.0;
                valzetamapReal = 0.5*vectorin[0]+0.5*vectorin[1]+vectorin[2];
                break;
            case 8:
                valximapReal =  0.5*vectorin[0]-0.5*vectorin[2]+0.5;
                valetamapReal = 0.0 ;
                valzetamapReal =-0.5*vectorin[0] +0.5*vectorin[2] + 0.5 ;
                break;
            case 9:
                valximapReal = 0.0;
                valetamapReal =  0.5*vectorin[1]-0.5*vectorin[2]+0.5;
                valzetamapReal = -0.5*vectorin[1]+0.5*vectorin[2]+0.5;
                break;
            case 10:
                valximapReal =  vectorin[0];
                valetamapReal = vectorin[1];
                valzetamapReal = 0.0;
                break;
            case 11:
                valximapReal =  vectorin[0];
                valetamapReal =0.0;
                valzetamapReal = vectorin[2];
                break;
            case 12:
                valximapReal =  (1.0/3.0)*(-vectorin[0]+2*vectorin[1]-vectorin[2]+1);
                valetamapReal =  (1.0/3.0)*(-vectorin[0]-vectorin[1]+2*vectorin[2]+1);
                valzetamapReal = 1.0 - valximapReal -valetamapReal;
                break;
            case 13:
                valximapReal =  0.0;
                valetamapReal =vectorin[1];
                valzetamapReal = vectorin[2];
                break;
                
        }
        
        TPZVec<double> vectorout_inv(3);
        transformInv.Apply(vectorout, vectorout_inv);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
        std::cout<<"TransformSideToElement("<<iside<<") = "<<vectorout_inv<<std::endl;
        REAL xierror =abs(vectorout_inv[0] - valximapReal);
        REAL etaerror =abs(vectorout_inv[1] - valetamapReal);
        REAL zetaerror =abs(vectorout_inv[2] - valzetamapReal);
        REAL tol = 1.0e-6;
        if (xierror<tol && etaerror<tol && zetaerror<tol) {
            std::cout<<"Transformation is ok!..."<<std::endl;
            std::cout<<"************************ "<<std::endl;
        }
        else{
            std::cout<<"Error!..."<<std::endl;
            counterror++;
        }
    }
    std::cout<<"Errors number: "<<counterror<<std::endl;
    
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

void TestPrismaTransforms(TPZVec<double> vectorin){
    
    std::cout<<"************************ "<<std::endl;
    std::cout<<"TransformElementToSide TEST->Cube"<<std::endl;
    std::cout<<"Parametric Domain 2D {xi,-1,1}, {eta,-1,1} "<<std::endl;
    std::cout<<"Parametric Domain 1D {xi,-1,1}"<<std::endl;
    std::cout<<"************************ "<<std::endl;
    
    TPZVec<double> vectorout(1);
    
    std::cout<<"Point to transform: "<<vectorin<<std::endl;
    int counterror=0;
    pztopology::TPZPrism prism;
    int npoints = prism.NCornerNodes;
    int nfaces = prism.NFaces;
    int nribs =  prism.NSides - npoints - nfaces -1;
    for (int iside = npoints; iside<npoints + nfaces + nribs; iside++) {
        if (iside>=15) {
            vectorout.resize(2);
        }
        TPZTransform<> transform = prism.TransformElementToSide(iside);
        TPZTransform<> transformInv = prism.TransformSideToElement(iside);
        transform.Apply(vectorin, vectorout);
        REAL valximapReal =0.0;
        REAL valetamapReal =0.0;
        REAL valzetamapReal =0.0;
        switch (iside) {
            case 8:
                valximapReal =vectorin[0];
                valetamapReal = -1.0;
                valzetamapReal = -1.0;
                break;
            case 9:
                valximapReal =  1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = -1.0;
                break;
            case 10:
                valximapReal =vectorin[0];
                valetamapReal =  1.0;
                valzetamapReal = -1.0;
                break;
            case 11:
                valximapReal =  -1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = -1.0;
                break;
            case 12:
                valximapReal =  -1.0;
                valetamapReal = -1.0;
                valzetamapReal = vectorin[2];
                break;
            case 13:
                valximapReal =  1.0;
                valetamapReal = -1.0;
                valzetamapReal = vectorin[2];
                break;
            case 14:
                valximapReal =  1.0;
                valetamapReal = 1.0;
                valzetamapReal = vectorin[2];
                break;
            case 15:
                valximapReal =  -1.0;
                valetamapReal =  1.0;
                valzetamapReal = vectorin[2];
                break;
            case 16:
                valximapReal =vectorin[0];
                valetamapReal = -1.0;
                valzetamapReal = 1.0;
                break;
            case 17:
                valximapReal =  1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = 1.0;
                break;
            case 18:
                valximapReal =vectorin[0];
                valetamapReal =  1.0;
                valzetamapReal = 1.0;
                break;
            case 19:
                valximapReal =  -1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = 1.0;
                break;
            case 20:
                valximapReal =  vectorin[0];
                valetamapReal = vectorin[1];
                valzetamapReal = -1.0;
                break;
            case 21:
                valximapReal =  vectorin[0];
                valetamapReal = -1.0;
                valzetamapReal = vectorin[2];
                break;
            case 22:
                valximapReal =  1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = vectorin[2];
                break;
            case 23:
                valximapReal =  vectorin[0];
                valetamapReal = 1.0;
                valzetamapReal = vectorin[2];
                break;
            case 24:
                valximapReal =  -1.0;
                valetamapReal = vectorin[1];
                valzetamapReal = vectorin[2];
                break;
            case 25:
                valximapReal =  vectorin[0];
                valetamapReal = vectorin[1];
                valzetamapReal = 1.0;
                break;
                
        }
        
        TPZVec<double> vectorout_inv(3);
        transformInv.Apply(vectorout, vectorout_inv);
        std::cout<<"TransformElementToSide("<<iside<<") = "<<vectorout<<std::endl;
        std::cout<<"TransformSideToElement("<<iside<<") = "<<vectorout_inv<<std::endl;
        REAL xierror =abs(vectorout_inv[0] - valximapReal);
        REAL etaerror =abs(vectorout_inv[1] - valetamapReal);
        REAL zetaerror =abs(vectorout_inv[2] - valzetamapReal);
        REAL tol = 1.0e-6;
        if (xierror<tol && etaerror<tol && zetaerror<tol) {
            std::cout<<"Transformation is ok!..."<<std::endl;
            std::cout<<"************************ "<<std::endl;
        }
        else{
            std::cout<<"Error!..."<<std::endl;
            counterror++;
        }
    }
    std::cout<<"Errors number: "<<counterror<<std::endl;
    
    return;
    
}

template void Shape<pzshape::TPZShapeLinear>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> phi,TPZFMatrix<REAL> dphi);
template void Shape<pzshape::TPZShapeTriang>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> phi,TPZFMatrix<REAL> dphi);
template void Shape<pzshape::TPZShapeQuad>(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> phi,TPZFMatrix<REAL> dphi);

//void Shape(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> phi,TPZFMatrix<REAL> dphi);
