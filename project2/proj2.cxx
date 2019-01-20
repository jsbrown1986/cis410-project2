/*=========================================================================

			Program 2 for CIS 410 (W18)
			BoundingBoxForCell, CountNumberOfStraddlingCells, 
			EvaluateFieldAtLocation, and Interpolation functions
			implemented by Jacob Brown 1/18/2018

===========================================================================*/

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// Interpolation Function
float Interpolation(float value, float a, float b, float fa, float fb) {
	return (b - a) != 0 ? fa + ((value - a) / (b - a)) * (fb - fa) : 0;
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float *pt, const int *dims,
	const float *X, const float *Y, const float *F) {

	int bottomLeftPointGet[2];
	int bottomRightPointGet[2];
	int topLeftPointGet[2];
	int topRightPointGet[2];

	// Checks to see if values are outside arrays before searching
	if (pt[0] < X[0] || pt[0] >= X[dims[0] - 1] || pt[1] < Y[0] || pt[1] >= Y[dims[1] - 1])
		return 0;

	// Run through X array, seeing if its between two values
	for (int i = 0; i < dims[0] - 1; i++) {
		if (X[i] <= pt[0] && X[i + 1] > pt[0]) {
			bottomLeftPointGet[0] = i;
			break;
		}
	}

	// Run through Y-array, seeing if its between two values
	for (int i = 0; i < dims[1] - 1; i++) {
		if (Y[i] <= pt[1] && Y[i + 1] > pt[1]) {
			bottomLeftPointGet[1] = i;
			break;
		}
	}

	// Uses logical point index to find the other indices of the cell
	bottomRightPointGet[0] = bottomLeftPointGet[0] + 1;
	bottomRightPointGet[1] = bottomLeftPointGet[1];

	topLeftPointGet[0] = bottomLeftPointGet[0];
	topLeftPointGet[1] = bottomLeftPointGet[1] + 1;

	topRightPointGet[0] = bottomLeftPointGet[0] + 1;
	topRightPointGet[1] = bottomLeftPointGet[1] + 1;

	// Uses indices to find all four F-values of cell
	float bottom_left_f = F[GetPointIndex(bottomLeftPointGet, dims)];
	float bottom_right_f = F[GetPointIndex(bottomRightPointGet, dims)];
	float top_left_f = F[GetPointIndex(topLeftPointGet, dims)];
	float top_right_f = F[GetPointIndex(topRightPointGet, dims)];

	// Interpolates between bottom left & right and top left & right
	float bottom_interp = Interpolation(pt[0], X[bottomLeftPointGet[0]], X[bottomRightPointGet[0]], bottom_left_f, bottom_right_f);
	float top_interp = Interpolation(pt[0], X[topLeftPointGet[0]], X[topRightPointGet[0]], top_left_f, top_right_f);

	// Returns final interpolation
    return Interpolation(pt[1], Y[bottomLeftPointGet[1]], Y[topLeftPointGet[1]], bottom_interp, top_interp);
}

// ****************************************************************************
//  Function: BoundingBoxForCell
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     cellId: a cellIndex (I.e., between 0 and GetNumberOfCells(dims))
//     bbox (output): the bounding box of cellId.  Format should be
//                     bbox[0]: the minimum X value in cellId.
//                     bbox[1]: the maximum X value in cellId.
//                     bbox[2]: the minimum Y value in cellId.
//                     bbox[3]: the maximum Y value in cellId.
//
//  Returns:  None (argument bbox is output)
//
// ****************************************************************************

void
BoundingBoxForCell(const float *X, const float *Y, const int *dims,
	int cellId, float *bbox){

	int indexGet[2];

	// Gets logical cell index for cellId
	GetLogicalCellIndex(indexGet, cellId, dims);

	// Verifies cellId validity, assigns accordingly
	if (cellId < GetNumberOfCells(dims) && cellId >= 0) {
		bbox[0] = X[indexGet[0]];
		bbox[1] = X[indexGet[0] + 1];
		bbox[2] = Y[indexGet[1]];
		bbox[3] = Y[indexGet[1] + 1];
	} else {
		bbox[0] = -100;
		bbox[1] = +100;
		bbox[2] = -100;
		bbox[3] = +100;
	}
}

// ****************************************************************************
//  Function: CountNumberOfStraddingCells
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//  Returns:  the number of cells that straddle 0, i.e., the number of cells
//            that contains points that have F>0 and also have points with F<0.
//
// ****************************************************************************

int
CountNumberOfStraddlingCells(const float *X, const float *Y, const int *dims,
	const float *F){

	int counter = 0, indexGet[2], pointGet[2];
	bool positive = false, negative = false;

	// Runs through each cell
	for (int i = 0; i < GetNumberOfCells(dims); i++) {

		// Gets logical cell index at each cell
		GetLogicalCellIndex(indexGet, i, dims);

		// Gets point index at each corner of cell
		// (Bottom Left (xmin, ymin))
		pointGet[0] = indexGet[0];
		pointGet[1] = indexGet[1];

		if (F[GetPointIndex(pointGet, dims)] > 0) {
			positive = true;
		} else if (F[GetPointIndex(pointGet, dims)] < 0) {
			negative = true;
		}

		// (Bottom Right (xmax, ymin))
		pointGet[0] = indexGet[0] + 1;
		pointGet[1] = indexGet[1];

		if (F[GetPointIndex(pointGet, dims)] > 0) {
			positive = true;
		} else if (F[GetPointIndex(pointGet, dims)] < 0) {
			negative = true;
		}

		// (Top Left (xmin, ymax))
		pointGet[0] = indexGet[0];
		pointGet[1] = indexGet[1] + 1;

		if (F[GetPointIndex(pointGet, dims)] > 0) {
			positive = true;
		} else if (F[GetPointIndex(pointGet, dims)] < 0) {
			negative = true;
		}

		// (Top Right (xmax, ymax))
		pointGet[0] = indexGet[0] + 1;
		pointGet[1] = indexGet[1] + 1;

		if (F[GetPointIndex(pointGet, dims)] > 0) {
			positive = true;
		} else if (F[GetPointIndex(pointGet, dims)] < 0) {
			negative = true;
		}

		// If there's a straddle, count the cell
		if (positive && negative)
			counter++;

		positive = false;
		negative = false;
	}

    return counter;
}

int main()
{
    int  i;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj2_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int numCells = CountNumberOfStraddlingCells(X, Y, dims, F);
    cerr << "The number of cells straddling zero is " << numCells << endl;

    float bbox[4];
    const int ncells = 5;
    int cellIds[ncells] = { 0, 50, 678, 1000, 1200 };
    for (i = 0 ; i < ncells ; i++)
    {
        BoundingBoxForCell(X, Y, dims, cellIds[i], bbox);
        cerr << "The bounding box for cell " << cellIds[i] << " is " 
             << bbox[0] << "->" << bbox[1] << ", " << bbox[2] << "->" << bbox[3]
             << endl;
    }

    const int npts = 10;
    float pt[npts][3] = 
         {
            {1.01119, 0.122062, 0},
            {0.862376, 1.33839, 0},
            {0.155026, 0.126123, 0},
            {0.69736, 0.0653565, 0},
            {0.2, 0.274117, 0},
            {0.893699, 1.04111, 0},
            {0.608791, -0.0533753, 0},
            {1.00543, 0.138024, 0},
            {0.384128, -0.0768977, 0},
            {0.666757, 0.60259, 0},
         };

    

    for (i = 0 ; i < npts ; i++)
    {
        float f = EvaluateFieldAtLocation(pt[i], dims, X, Y, F);
        cerr << "Evaluated field at (" << pt[i][0] <<"," << pt[i][1] << ") as "
             << f << endl;
    }
    
   
    cerr << "Infinite loop here, else Windows people may have the terminal "
         << "disappear before they see the output."
         << " Remove these lines if they annoy you." << endl;
    cerr << "(press Ctrl-C to exit program)" << endl;
    while (1) ; 
}