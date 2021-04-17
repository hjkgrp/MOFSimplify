
/* Functions to handle symmetry operatons and equivalent atom detection

   Spacegroups taken from http://it.iucr.org/A/

   added by M. Haranczyk, 11/11/11

*/

vector <XYZ> GetEquivalentPositions(int spacegroup,XYZ *pt);
bool IsUniqueVertex(XYZ *p, ATOM_NETWORK &cell);
int get_sym_ID(string s);

