//const boundBox& bounds = mesh.bounds();

vector vec1(bounds.max().x(), bounds.min().y(), bounds.min().z());
vector vec5(bounds.max().x(), bounds.min().y(), bounds.max().z());

vector vec2(bounds.max().x(), bounds.max().y(), bounds.min().z());
vector vec6(bounds.max().x(), bounds.max().y(), bounds.max().z());

vector midBottom(
                    (bounds.max().x() - bounds.min().x())/2,
                    (bounds.max().y() - bounds.min().y())/2,
                    bounds.min().z()
                );

scalar height = bounds.max().z() - bounds.min().z();
scalar length = bounds.max().x() - bounds.min().x();
scalar width = bounds.max().y() - bounds.min().y();
//Info<< height << endl;

label n1 = mesh.findCell(vec1);
label n5 = mesh.findCell(vec5);

label n2 = mesh.findCell(vec2);
label n6 = mesh.findCell(vec6);

label nBottom = mesh.findCell(midBottom);

vector D1 = D[n1];
vector D5 = D[n5];
vector delta15 = D1 - D5;

vector D2 = D[n2];
vector D6 = D[n6];
vector delta26 = D2 - D6;

vector DBottom = D[nBottom];

//Info << "delta12 = " << delta12 << endl;
Info<< "shrinkage15 = "<< delta15.component(2)/height*100 << "%" << endl;

//Info << "delta34 = " << delta34 << endl;
Info<< "shrinkage26 = "<< delta26.component(2)/height*100 << "%" << endl;

Info<< "deltaZ(bottom) = " << DBottom.component(2) << endl;

Info<< "deltaZ1 = " << D1.component(2) << endl;

Info<< "deltaZ2 = " << D2.component(2) << endl;

Info<< "curve1 = " << (DBottom.component(2) - D1.component(2))/length*100 << "%" << endl;

Info<< "curve2 = " << (DBottom.component(2) - D2.component(2))/length*100 << "%" << endl;

Info<< "deltaX2 = " << D2.component(0)/length*100 << "%" << endl;

Info<< "deltaX1 = " << D1.component(0)/length*100 << "%" << endl;

Info<< "deltaY1 = " << D1.component(1)/width*100 << "%" << endl;

Info<< "deltaY2 = " << D2.component(1)/width*100 << "%" << endl;



