// OMS, a library for determining the exposure grade of a molecule
//
// OMS.hh
//
// Routines for single molecule analysis
//
// Author   : Ismael Gómez García
// Email    : 
// Date     : November 27th 2015

#include "OMS.h"


void TransformToSpherical(double *TargetMolecule)
{

	double x_coord, y_coord, z_coord;
	double Aux;
	
	// Save the coordinates in order to avoid overwritting
	x_coord = TargetMolecule[X_COORD];
	y_coord = TargetMolecule[Y_COORD];
	z_coord = TargetMolecule[Z_COORD];

	

	// Calculate the phi coordinate (interval [0, 2*PI] )
	TargetMolecule[PHI_COORD] = atan2(y_coord, x_coord);
	if (TargetMolecule[PHI_COORD] < 0)
		TargetMolecule[PHI_COORD] = 2*PI + TargetMolecule[PHI_COORD]; 

	// Calcuate the theta coordinate (interval [0, PI] )
	Aux = sqrt( x_coord*x_coord + y_coord*y_coord );
	TargetMolecule[THETA_COORD] = atan2(Aux, z_coord);
	if (TargetMolecule[THETA_COORD] < 0)
		TargetMolecule[THETA_COORD] = PI + TargetMolecule[THETA_COORD]; 
	
	

	// Calculate the radius (not really necessary)
	//SphericalCoordinates[R_COORD] = sqrt

}

void Sort(	double **MoleculeCoordinates, 
		int Criterion,
		int Size)
{

	int i,j;
	double *Swap;

	for (i = 0 ; i < Size-1; i++)
		for (j = 0 ; j < Size - i - 1; j++)
			if (MoleculeCoordinates[j][Criterion] > MoleculeCoordinates[j+1][Criterion]) 
			{
				Swap = MoleculeCoordinates[j];
				MoleculeCoordinates[j] = MoleculeCoordinates[j+1];
				MoleculeCoordinates[j+1] = Swap;
			}
}


/*
 *	Checks if two vectors are linearly dependent between them
 */
bool LinearlyDependent (double *Vector1, double *Vector2)
{

	Vector3d v1(Vector1[X_COORD], Vector1[Y_COORD], Vector1[Z_COORD]);
	Vector3d v2(Vector2[X_COORD], Vector2[Y_COORD], Vector2[Z_COORD]);

	v1 = v1.cross(v2);

	//cout << "Checking LD: \n" << Vector1[0] << Vector1[1] << Vector1[2] << "\n" << Vector2[0] << Vector2[1] << Vector2[2] << "\n" << v1 << "\n";	

	//if (v1(0) == 0 && v1(1) == 0 && v2(0) == 0)
	if ( v1(0)*v1(0)+v1(1)*v1(1)+v1(2)*v1(2) == 0)
	{
		//cout << "Vectors are lin dep " << v1(0)*v1(0)+v1(1)*v1(1)+v1(2)*v1(2) << "\n";
		return true;
	}

	
	return false;



}


/* For a plane given through two vectors (departing from the origin), and a point, returns:
 * 0 if the point is in the plane
 * -1 if the point is in the "bottom" of the plane (normal sign negative)
 * +1 if the point is in the "top" of the plane (normal sign positive)
 */
double PlaneSide (double *Vector1, double *Vector2, double *Point)
{

	// Elements for system resolution
	Matrix3d A;
	Vector3d b(Point[X_COORD], Point[Y_COORD], Point[Z_COORD]);
	Vector3d x;

	// Vectors defining the plane
	Vector3d v1(Vector1[X_COORD], Vector1[Y_COORD], Vector1[Z_COORD]);
	Vector3d v2(Vector2[X_COORD], Vector2[Y_COORD], Vector2[Z_COORD]);
	Vector3d n(Vector1[X_COORD], Vector1[Y_COORD], Vector1[Z_COORD]);
	
	n = n.cross(v2);


	/*cout << "Vector 1: \n" << v1 << "\n";
	cout << "Vector 2: \n" << v2 << "\n";
	cout << "Normal: \n" << n << "\n";

	cout << "B: \n" << b << "\n";*/

	/**** Fill the matrix ****/
	// Each column contains one vector, so the system solved returns the coordinates of the point in the system v1,v2,n
	for (int i = 0; i < DIM; i++)
		for (int j=0; j<DIM; j++)
			if (i == 0)			
				A(j,i) = v1(j);
			else if (i == 1) 
				A(j,i) = v2(j);
			else 	
				A(j,i) = n(j);

	//cout << "A\n" << A << "\n";

	/**** Compute the system ****/
	x = A.colPivHouseholderQr().solve(b);


	//cout << "x: \n" << x << "\n";

	/**** Return the value of the third element (lambda parameter of the normal) ****/
	return x(Z_COORD);

}

int Sign(double Lambda)
{
	if (Lambda > 0)
		return 1;

	if (Lambda == 0)
		return 0;

	return -1;
}


bool IsExposedMolecule (vector < vector <double> > MoleculeCoordinates  )
{

	int NPoints = (int) MoleculeCoordinates.size();
	//vector < vector <double> > MoleculeSphericalCoordinates(Size, vector<double>(DIM) );
	//vector < vector <double> > MoleculeSphericalCoordinates(Size);

	double **TranslatedCoordinates = (double **) malloc((NPoints-1)*sizeof(double *));
	double *OriginCoordinates = (double *) malloc (DIM*sizeof(double));
	
	int i,j,k;

	double Aux;

	double Lambda = 0.;
	int S_lambda_old = 0, S_lambda;

	bool PlaneFound;


	// Assign initial coordinates
	for (i = 0; i < DIM; i++)
		OriginCoordinates[i] = MoleculeCoordinates.at(0).at(i);


	// Displace all the coordinates in order to set the metal molecule at the center
	for(i = 1; i < NPoints; i++)
	{

		TranslatedCoordinates[i-1] = (double *) malloc(sizeof(double)*DIM);
		// Assign values from input to the spherical coordinates vector and translate it to origin
		for (j = 0; j < DIM; j++)
		{
			Aux = MoleculeCoordinates.at(i).at(j);
			TranslatedCoordinates[i-1][j] = Aux - OriginCoordinates[j];
		}


	}

	/* For every two non-metal atoms in the molecule, form a vector (M-O), then a plane, then decide whether the rest of the poins are
	 * - Either in the same plane
 	 * - Either at the same side of the plane
	 */
	for (i = 0; i < NPoints-1; i++)
		for (j = i+1; j < NPoints-1; j++)
		{

			// Reinitialize the variables to avoid error
			S_lambda_old = 0;
			S_lambda = 0;
			Lambda = 0;
			PlaneFound = true;
			// For every point except the two composing the plane, compute
			if (!LinearlyDependent(TranslatedCoordinates[i], TranslatedCoordinates[j]))
			{		
				for (k = 0; k < NPoints-1; k++)
					if (k != i && k != j ) 
					{
						Lambda = PlaneSide (TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);

						//cout << "k: " << k << "\nLambda: " << Lambda << "\n";
			
						S_lambda = Sign(Lambda);

						// Update S_lambda_old only the first time that the point is not in the plane				
						if (Lambda != 0 && S_lambda_old == 0)
							S_lambda_old = Sign(Lambda);

						// If the sign changes anytime, stop this loop (the plane has points on both sides)
						if (S_lambda != S_lambda_old)
						{
							PlaneFound = false;
							//cout << "This plane separates the points \n";
							break;
						}
					}
			}
			else 	// If the vectors are linearly dependent the iteration is ignored
			{
				PlaneFound = false;
				//cout << "Setting PlaneFound false \n";
			}

			// Check if a plane with all the points on one side was found
			if (PlaneFound)
				return true;
		}
						
					

	return false;

}


double PlaneAngle (double *Vector1, double *Vector2, double *Point)
{

	
	double Angle;

	// Elements for system resolution
	Matrix3d A;
	Vector3d b(Point[X_COORD], Point[Y_COORD], Point[Z_COORD]);
	Vector3d x;

	// Vectors defining the plane
	Vector3d v1(Vector1[X_COORD], Vector1[Y_COORD], Vector1[Z_COORD]);
	Vector3d v2(Vector2[X_COORD], Vector2[Y_COORD], Vector2[Z_COORD]);
	Vector3d n(Vector1[X_COORD], Vector1[Y_COORD], Vector1[Z_COORD]);
	
	n = n.cross(v2);

	/**** Fill the matrix ****/
	// Each column contains one vector, so the system solved returns the coordinates of the point in the system v1,v2,n
	for (int i = 0; i < DIM; i++)
		for (int j=0; j<DIM; j++)
			if (i == 0)			
				A(j,i) = v1(j);
			else if (i == 1) 
				A(j,i) = v2(j);
			else 	
				A(j,i) = n(j);

	//cout << "A\n" << A << "\n";

	/**** Compute the system ****/
	x = A.colPivHouseholderQr().solve(b);


	// Multiply normal by its coordinate
	n = x(Z_COORD)*n;


//cout << "Compute angle from " << n.norm() << " and " << b.norm() << "\n";
	Angle = asin(n.norm()/b.norm());
	
	/**** Return the value of the third element (lambda parameter of the normal) ****/
	return Angle;	



	// Elements for system resolution
	/*Vector3d VectorToPoint(Point[X_COORD], Point[Y_COORD], Point[Z_COORD]);

	// Variable for angle computation
	double Angle;
	double PVMod, PNMod;

	// Vectors defining the plane
	Vector3d v2(Vector2[X_COORD], Vector2[Y_COORD], Vector2[Z_COORD]);
	Vector3d PlaneNormal(Vector1[X_COORD], Vector1[Y_COORD], Vector1[Z_COORD]);
	
	PlaneNormal = PlaneNormal.cross(v2);

	PVMod = sqrt(	VectorToPoint(0)*VectorToPoint(0) + VectorToPoint(1)*VectorToPoint(1) + VectorToPoint(2)*VectorToPoint(2) );
	PNMod = sqrt(	PlaneNormal(0)*PlaneNormal(0) + PlaneNormal(1)*PlaneNormal(1) + PlaneNormal(2)*PlaneNormal(2) );


cout << "Compute angle from " << PNMod/PVMod << "\n";
	Angle = asin(PNMod/PVMod);


	return Angle;*/

}

// Degree of exposure of the molecule (in radians)
/*bool IsExposedMoleculeThreshold (vector < vector <double> > MoleculeCoordinates, double Threshold  )
{

	int NPoints = (int) MoleculeCoordinates.size();

	double **TranslatedCoordinates = (double **) malloc((NPoints-1)*sizeof(double *));
	double *OriginCoordinates = (double *) malloc (DIM*sizeof(double));
	
	int i,j,k;

	// Angles
	double Theta = 0.;
	double Exposure = 0.;

	// Plane side
	double Lambda = 0.;
	int S_lambda_old = 0, S_lambda = 0;

	bool PlaneFound;


	// Assign initial coordinates
	for (i = 0; i < DIM; i++)
		OriginCoordinates[i] = MoleculeCoordinates.at(0).at(i);


	// Displace all the coordinates in order to set the metal molecule at the center
	for(i = 1; i < NPoints; i++)
	{
		TranslatedCoordinates[i-1] = (double *) malloc(sizeof(double)*DIM);
		// Assign values from input to the spherical coordinates vector and translate it to origin
		for (j = 0; j < DIM; j++)
			TranslatedCoordinates[i-1][j] = MoleculeCoordinates.at(i).at(j) - OriginCoordinates[j];
	}

	
	for (i = 0; i < NPoints-1; i++)
		for (j = i+1; j < NPoints-1; j++)
		{

			//cout << "Iteration " << i << j << "\n";

			// Reinitialize the variables to avoid error
			S_lambda_old = 0;
			S_lambda = 0;
			Lambda = 0;
			PlaneFound = true;
			Theta = 0.;
			Exposure = 0.;
			// For every point except the two composing the plane, compute
			if (!LinearlyDependent(TranslatedCoordinates[i], TranslatedCoordinates[j]))
			{		
				for (k = 0; k < NPoints-1; k++)
					if (k != i && k != j ) 
					{
						// Compute the side of the plane
						Lambda = PlaneSide (TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);
									
						S_lambda = Sign(Lambda);

						//cout << "k: " << k << " Lambda " << Lambda << "\n";

						// Update S_lambda_old only the first time that the point is not in the plane				
						if (Lambda != 0 && S_lambda_old == 0)
							S_lambda_old = Sign(Lambda);

						// If the sign changes anytime, check the angle 
						if (S_lambda != S_lambda_old)
						{
							Theta = PlaneAngle(TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);

							//cout << "\nTheta: " << Theta << "\tThreshold: "  << Threshold << "\n";
							
							// If the absolute value exceeds the threshold, stop this iteration
							if (abs(Theta) > Threshold)
							{
								PlaneFound = false;
								//cout << "This plane separates the points \n";
								break;
							}
							else if (abs(Theta) > Exposure)
							{
								Exposure = abs(Theta);
							}
						}
					}
			}
			else 	// If the vectors are linearly dependent the iteration is ignored
			{
				PlaneFound = false;
				//cout << "Setting PlaneFound false \n";
			}

			// Check if a plane with all the points on one side was found
			if (PlaneFound)
				return true;
				//return Exposure;
		}
						
					

	//return -Threshold;
	return false;
	//return Exposure;

}*/

bool IsExposedMoleculeThreshold (vector < vector <double> > MoleculeCoordinates, double Threshold  )
{

	int NPoints = (int) MoleculeCoordinates.size();

	double **TranslatedCoordinates = (double **) malloc((NPoints-1)*sizeof(double *));
	double *OriginCoordinates = (double *) malloc (DIM*sizeof(double));
	
	int i,j,k;

	// Angles
	double Theta = 0.;
	double Exposure_pos = 0., Exposure_neg = 0., Exposure_plane = 0.;

	// Plane side
	double Lambda = 0.;
	int S_lambda_old = 0, S_lambda = 0;


	// Assign initial coordinates
	for (i = 0; i < DIM; i++)
		OriginCoordinates[i] = MoleculeCoordinates.at(0).at(i);


	// Displace all the coordinates in order to set the metal molecule at the center
	for(i = 1; i < NPoints; i++)
	{
		TranslatedCoordinates[i-1] = (double *) malloc(sizeof(double)*DIM);
		// Assign values from input to the spherical coordinates vector and translate it to origin
		for (j = 0; j < DIM; j++)
			TranslatedCoordinates[i-1][j] = MoleculeCoordinates.at(i).at(j) - OriginCoordinates[j];
	}

	/* For every two non-metal atoms in the molecule, form a vector (M-O), then a plane, then decide whether the rest of the poins are
	 * - Either in the same plane
 	 * - Either at the same side of the plane
	 */
	for (i = 0; i < NPoints-1; i++)
		for (j = i+1; j < NPoints-1; j++)
		{

			//cout << "Plane " << i << j << "\n";

			// Reinitialize the variables to avoid error
			S_lambda_old = 0;
			S_lambda = 0;
			Lambda = 0;
			Theta = 0.;

			Exposure_pos = 0.;
			Exposure_neg = 0.;

			// For every point except the two composing the plane, compute
			if (!LinearlyDependent(TranslatedCoordinates[i], TranslatedCoordinates[j]))
			{	
				// POINTS ITERATION: Compare the points with the plane selected in previous iterations
				for (k = 0; k < NPoints-1; k++)
					if (k != i && k != j ) 
					{
						// Compute the side of the plane
						Lambda = PlaneSide (TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);
						//cout << "k: " << k << " Lambda " << Lambda << "\n";

						S_lambda = Sign(Lambda);

						// Store the angles at both sides of the plane 
						// The exposure angle is the minimum of both
						if (S_lambda > 0)
						{
							Theta = PlaneAngle(TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);
							if (Theta > Exposure_pos)							
								Exposure_pos = Theta;
						}
						if (S_lambda < 0)
						{
							Theta = PlaneAngle(TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);
							if (Theta > Exposure_neg)							
								Exposure_neg = Theta;
						}

					}

				// Compute exposures (only if the vector forms a plane)
				Exposure_plane = min(Exposure_pos, Exposure_neg);
				if (Exposure_plane < Threshold)
					return true;

			}

			

		}
						
					
	// Maximum angle formed by a point not contained within the same side of the rest of the points
	return false;




}

double DegreeOfExposure(vector < vector <double> > MoleculeCoordinates)
{

	
	int NPoints = (int) MoleculeCoordinates.size();

	double **TranslatedCoordinates = (double **) malloc((NPoints-1)*sizeof(double *));
	double *OriginCoordinates = (double *) malloc (DIM*sizeof(double));
	
	int i,j,k;

	// Angles
	double Theta = 0.;
	double Exposure_pos = 0., Exposure_neg = 0., Exposure = PI/2., Exposure_plane = 0.;

	// Plane side
	double Lambda = 0.;
	int S_lambda_old = 0, S_lambda = 0;

	// Assign initial coordinates
	for (i = 0; i < DIM; i++)
		OriginCoordinates[i] = MoleculeCoordinates.at(0).at(i);


	// Displace all the coordinates in order to set the metal molecule at the center
	for(i = 1; i < NPoints; i++)
	{
		TranslatedCoordinates[i-1] = (double *) malloc(sizeof(double)*DIM);
		// Assign values from input to the spherical coordinates vector and translate it to origin
		for (j = 0; j < DIM; j++)
			TranslatedCoordinates[i-1][j] = MoleculeCoordinates.at(i).at(j) - OriginCoordinates[j];
	}

	/* For every two non-metal atoms in the molecule, form a vector (M-O), then a plane, then decide whether the rest of the poins are
	 * - Either in the same plane
 	 * - Either at the same side of the plane
	 */
	for (i = 0; i < NPoints-1; i++)
		for (j = i+1; j < NPoints-1; j++)
		{

			//cout << "Plane " << i << j << "\n";

			// Reinitialize the variables to avoid error
			S_lambda_old = 0;
			S_lambda = 0;
			Lambda = 0;
			Theta = 0.;

			Exposure_pos = 0.;
			Exposure_neg = 0.;
				
			// For every point except the two composing the plane, compute
			if (!LinearlyDependent(TranslatedCoordinates[i], TranslatedCoordinates[j]))
			{	
				// POINTS ITERATION: Compare the points with the plane selected in previous iterations
				for (k = 0; k < NPoints-1; k++)
					if (k != i && k != j ) 
					{
						// Compute the side of the plane
						Lambda = PlaneSide (TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);
						//cout << "k: " << k << " Lambda " << Lambda << "\n";

						S_lambda = Sign(Lambda);

						// Store the angles at both sides of the plane 
						// The exposure angle is the minimum of both
						if (S_lambda > 0)
						{
							Theta = PlaneAngle(TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);
							if (Theta > Exposure_pos)							
								Exposure_pos = Theta;
						}
						if (S_lambda < 0)
						{
							Theta = PlaneAngle(TranslatedCoordinates[i], TranslatedCoordinates[j], TranslatedCoordinates[k]);
							if (Theta > Exposure_neg)							
								Exposure_neg = Theta;
						}

					}

				// Compute exposures (only if the vector forms a plane)
				Exposure_plane = min(Exposure_pos, Exposure_neg);
				/*cout << "Exposure plane: \t" << Exposure_plane << "\n";
				cout << "Exposure: \t" << Exposure << "\n";*/

				if (Exposure_plane < Exposure)
				{
					//cout << "Exposure updated: " << Exposure << "\n";			
					Exposure = Exposure_plane;
				}
			}

			

		}
						
					
	// Maximum angle formed by a point not contained within the same side of the rest of the points
	return Exposure;


}






