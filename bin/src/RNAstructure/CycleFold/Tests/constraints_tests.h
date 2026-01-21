#include "../constraints.h"
#include "../NCM_parameters.h"
#include <vector>
using std::vector;

constraints setup_constraints(){
    std::vector<int> v;
    for(int i=0;i<12;i++){
        v.push_back(i);
    }
    v[2] = 8;
    v[8] = 2;
    constraints cst = constraints(v);
    return cst;
}

constraints setup_unpair_constraints(){
    std::vector<bool> v;
    for(int i=0;i<12;i++){
        v.push_back(false);
    }
    v[2] = true;
    v[8] = true;
    constraints cst = constraints(v);
    return cst;
}

BOOST_AUTO_TEST_CASE(allowed_pair_works_correctly){
    const constraints cst = setup_constraints();
    BOOST_CHECK(allowed_pair(2,8,cst));
    BOOST_CHECK(allowed_pair(1,7,cst));
    BOOST_CHECK(!allowed_pair(2,7,cst));

	const constraints cst2 = setup_unpair_constraints();
	BOOST_CHECK(allowed_pair(3,7,cst2));
	BOOST_CHECK(!allowed_pair(2,7,cst2));
	BOOST_CHECK(!allowed_pair(1,8,cst2));
}

BOOST_AUTO_TEST_CASE(unpaired_works_correctly){
    const constraints cst = setup_constraints();
    BOOST_CHECK(cst.unpaired(1));
    BOOST_CHECK(cst.unpaired(3));
    BOOST_CHECK(!cst.unpaired(2));
    BOOST_CHECK(!cst.unpaired(8));
	const constraints cst2 = setup_unpair_constraints();
	BOOST_CHECK(cst2.unpaired(3));
}

BOOST_AUTO_TEST_CASE(partner_works_correctly){
    const constraints cst = setup_constraints();
    BOOST_CHECK_EQUAL(cst.partner(2),8);
    BOOST_CHECK_EQUAL(cst.partner(8),2);
}

BOOST_AUTO_TEST_CASE(conflicts_works_correctly){
    const constraints cst = setup_constraints();
    NCM_type a = std::string("232");
    NCM_type b = std::string("223");
    NCM_type h = std::string("14");
    //at the pair
    BOOST_CHECK(!conflicts(2,8,a,cst));
    BOOST_CHECK(!conflicts(2,8,h,cst));
    //bulge on 5' side
    BOOST_CHECK(conflicts(1,9,a,cst));
    BOOST_CHECK(!conflicts(0,9,a,cst));
    //bulge on 3' side
    BOOST_CHECK(conflicts(1,9,b,cst));
    BOOST_CHECK(!conflicts(1,10,b,cst));
    //hairpin
    BOOST_CHECK(conflicts(1,4,h,cst));
    BOOST_CHECK(!conflicts(3,6,h,cst));
	//check exterior case
    //at the pair
    BOOST_CHECK(!conflicts(8,2,a,cst));
    BOOST_CHECK(!conflicts(8,2,h,cst));
    //bulge on 5' side
    BOOST_CHECK(conflicts(7,3,a,cst));
    BOOST_CHECK(!conflicts(7,4,a,cst));
    //bulge on 3' side
    BOOST_CHECK(conflicts(7,1,b,cst));
    BOOST_CHECK(conflicts(7,2,b,cst));
    BOOST_CHECK(!conflicts(6,3,b,cst));
    BOOST_CHECK(!conflicts(8,2,b,cst));

	const constraints cst2 = setup_unpair_constraints();
	BOOST_CHECK(conflicts(2,8,a,cst2));
	BOOST_CHECK(conflicts(2,9,a,cst2));
    BOOST_CHECK(conflicts(5,8,h,cst2));
    BOOST_CHECK(!conflicts(1,4,h,cst2));
	//exterior case
	BOOST_CHECK(conflicts(8,2,a,cst2));
	BOOST_CHECK(conflicts(9,4,a,cst2));
	//five' bulge
	BOOST_CHECK(!conflicts(9,3,a,cst2));
	//three' bulge
	BOOST_CHECK(!conflicts(7,1,b,cst2));



	const constraints empty = constraints();
    BOOST_CHECK(!conflicts(1,10,b,empty));
    BOOST_CHECK(!conflicts(5,8,h,empty));
    BOOST_CHECK(!conflicts(1,4,h,empty));
}
