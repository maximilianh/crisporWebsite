/*
 * A namespace for functions that perform additional analysis on structures and RNA objects.
 *
 * (c) 2017 Mathews Lab, University of Rochester Medical Center.
 * Written by Richard M. Watson
 */


#include "StructureTools.h"
#include <limits>
#include <vector>
#include "loop_utils.h"
#include "../RNA_class/RNA.h"

namespace StructureTools {

// ------------- nuclayout class -------------------
nuclayout::nuclayout(const int seqlen) :
	seqlen(seqlen), 
	labelCount(seqlen/10+1),
	motifs(new int[seqlen]),
	coords(new geometry2D::vec2D[seqlen]),
	labels(new geometry2D::vec2D[labelCount]) {
		for(int i = 0; i < seqlen; i++) {
			coords[i].set(0,0);
			motifs[i] = 0;
		}
}
nuclayout::~nuclayout() {
	delete[] coords;
	delete[] motifs;
	delete[] labels;
}

int nuclayout::getCoords(float * const arr, const int arrSize) const {
	if (arrSize < 2*seqlen) return 2*seqlen;
	for(int i = 0; i < seqlen; i++) {
		arr[2*i]   = (float)coords[i].x;
		arr[2*i+1] = (float)coords[i].y;
	}
	return 0;
}
int nuclayout::getNumberLabels(float * const arr, const int arrSize) const {
	if (arrSize < 2*labelCount) return 2*labelCount;
	for(int i = 0; i < labelCount; i++) {
		arr[2*i]   = (float)labels[i].x;
		arr[2*i+1] = (float)labels[i].y;
	}
	return 0;
}

// ------------- ^ nuclayout class ^ -------------------

//StructureTools::StructureTools() : rna(NULL), ownsRNA(false) {
//	
//}
//
//StructureTools::~StructureTools() {
//	if (ownsRNA & rna != NULL)
//		delete rna;
//}

//	//Determine the coordinates for drawing a secondary structure.
//int RNA::DetermineDrawingCoordinates(const int height, const int width, const int structurenumber) {
//	int diameter = (int)ceil(sqrt((double)((width)*width + (height)*height)));
//	return DetermineDrawingCoordinates(diameter, (int)ceil(diameter * 0.4), (int)ceil(diameter * 1/8.0f), (int)ceil(diameter * 1/8.0f), structurenumber);
//}
//int StructureTools::
nuclayout * layoutAngular(RNA &rna, const int structurenumber, const float nucRadius, const float bondLength, const float helixSpacing, const float loopSpacing) {
	structure &ct = *rna.GetStructure();
	const int seqlen = ct.GetSequenceLength();
	
	nuclayout &layout = *new nuclayout(seqlen);

	//check to make sure that a sequence has been read
	if (seqlen==0) return &(layout.setError(20)); //No sequence has been read.

	//verify the structure number
	if (structurenumber<0||structurenumber>seqlen) return &(layout.setError(3)); //Structure number out of range

	return &layout;
}

}