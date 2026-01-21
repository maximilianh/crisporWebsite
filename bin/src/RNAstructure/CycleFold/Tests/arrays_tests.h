#include "../arrays.h"
#include "../NCM_parameters.h"
arrays<real_t> arr(10);

BOOST_AUTO_TEST_CASE(test_v){
    for(int i=0;i<num_pairtypes;i++){
        arr.set_v(1.0,i,0,0);
    }
    for(int i=0;i<(int)NCM::all_NCMs().size();i++){
        arr.set_v(1.0,i,0,0);
    }
}
