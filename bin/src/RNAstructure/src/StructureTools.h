/*
 * A namespace for functions that perform additional analysis on structures and RNA objects.
 *
 * (c) 2017 Mathews Lab, University of Rochester Medical Center.
 * Written by Richard M. Watson
 */

#ifndef STRUCTURE_TOOLS_H
#define STRUCTURE_TOOLS_H

#include "../RNA_class/RNA.h"
#include "geometry2D.h"

using namespace std;

namespace StructureTools {
	class nuclayout {
	public:
		nuclayout(const int seqlen);
		~nuclayout();
		int getCoords(float* const arr, const int arrSize) const;
		int getNumberLabels(float* const arr, const int arrSize) const;
		//int getCoords(vector<float>&arr) const;
		inline int size() const { return seqlen; }
		inline geometry2D::vec2D& operator[](const int index) const {
			return coords[index];
		}
		inline nuclayout& setError(const int code) { errorCode=code; return *this; }
		inline int error() const { return errorCode; }
#ifdef SWIG
	private:
#endif
		geometry2D::vec2D* coords;
		geometry2D::vec2D* labels;
		int* motifs;
		const int seqlen;
		const int labelCount;
	private:
		int errorCode;
	};

	//public:
	//	StructureTools();
	//	~StructureTools();
		nuclayout * layoutAngular(RNA &rna, const int structurenumber, const float nucRadius, const float bondLength, const float helixSpacing, const float loopSpacing);
	//private:
}

#endif /* STRUCTURE_TOOLS_H */
