const scalarField& diag = DEqn.diag();
const vectorField& source = DEqn.source();
const labelUList& neighbor = DEqn.lduAddr().lowerAddr();
const labelUList& owner = DEqn.lduAddr().upperAddr();
const scalarField& lower = DEqn.lower();
const scalarField& upper = DEqn.upper();

Info << "Data = [" << endl;
for(label i=0; i<diag.size(); i++)
{
Info << i+1 << "\t" << diag[i] << endl;
}
for(label f=0; f<upper.size(); f++)
{
Info << owner[f] << "\t" << neighbor[f] << "\t" << lower[f] << endl;
}
Info << "];\n " << "Matrix = sparse(Data(:,1), Data(:,2), Data(:,3));\n";
Info << "Source = [" << endl;
for(label i=0; i<source.size(); i++)
{
Info << i+1 << "\t" << source[i][0] << "\t" << source[i][1] << "\t" << source[i][2] << endl;
}
Info << "];\n";

Info<< "Source size = " << source.size() << endl;
