#include "nn_model.hpp"
#include "nn_model_efn2.hpp"

#include "../RNA_class/thermodynamics.h"

#include <vector>
#include <string>
#include <iostream>
#include <cstring>

using namespace std;


const string path = "../data_tables";


int main() {

	Thermodynamics th;

	th.ReadThermodynamic("../data_tables/");

	datatable *dt = th.GetDatatable();
	NNModel nnm(dt, "GGGCCC");
	nnm.SetRNA("AUGCUAGCUGUGACGGAAACC");

	cout << nnm.TwoLoop(0, 2, 7, 8) << endl;
	cout << nnm.TwoLoop(2, 3, 11, 13) << endl;
	cout << nnm.TwoLoop(14, 15, 19, 20) << endl;
	cout << nnm.OneLoop(2, 7) << endl;
	cout << nnm.MismatchCoax(0, 8, 2, 7) << endl;
	cout << nnm.ClosingFiveDangle(15, 19) << endl;

	NNModelLinear nnml(dt, "");
	cout << nnml.MLClosure(4, 5) << endl;

	NNModelEFN2 nnme(dt, "");
	cout << nnme.MLClosure(1.333333, 4, 6) << endl;
	cout << nnme.MLClosure(1.5673, 4, 11) << endl;
	cout << nnme.MLClosure(1.333333, 3, 1) << endl;
}
