#include "../NCM_parameters.h"
#include<iostream>
namespace parameters_tests{
//put setup code here
}

BOOST_AUTO_TEST_CASE(static_parameters){
	BOOST_CHECK(NCM::fivep_length("222")==2);
	BOOST_CHECK(NCM::fivep_length("232")==3);
	BOOST_CHECK(NCM::threep_length("222")==2);
	BOOST_CHECK(NCM::fivep_length("16")==6);
}

BOOST_AUTO_TEST_CASE(test_NCM_ID){
	NCM_id foo = NCM_id("222","AGCU");
	NCM_id bar = NCM_id(3,9,"222","GGGAGGGGCUGGG");
	BOOST_CHECK(foo==bar);
}

BOOST_AUTO_TEST_CASE(test_junction_id){
	junction_id foo ("16","222");
	junction_id bar ("16","222");
	//test id generation works as expected
	BOOST_CHECK(foo.getID() == std::string("16222"));
	//test operator== works as expected
	BOOST_CHECK(foo==bar);
	//test std::less works as expected
	//BOOST_CHECK(std::less<junction_id>(junction_id("15","222"));
}
BOOST_AUTO_TEST_CASE(test_hinge_id){
	hinge_id foo("222","232","cWW");
	hinge_id bar("222","232","cWW");
	hinge_id baz("16","222","cWW");
	BOOST_CHECK(foo.getID()==std::string("222232cWW"));
	BOOST_CHECK(foo==bar);
	BOOST_CHECK(!(foo==baz));

}
BOOST_AUTO_TEST_CASE(test_pair_id){
	pair_id foo("cWW","GC");
	pair_id bar = foo;
	pair_id baz("cHS","GC");
	BOOST_CHECK(foo.getID()==std::string("cWWGC"));
	BOOST_CHECK(foo==bar);
	BOOST_CHECK(!(bar==baz));
}

BOOST_AUTO_TEST_CASE(junction_energy_is_symmetrical){
	//add inner onto outer
	NCM_id outer = NCM_id("222","AGCU");
	NCM_id inner = NCM_id("222","GCGC");
	pairtype p = "cWW";
	double e1 = double(par.junction_energy(outer,inner,false));
	double e2 = double(par.junction_energy(inner,outer,true));
	BOOST_CHECK_CLOSE(e1, e2, 0.000001);

	NCM_id outer2 = NCM_id("232","AAGCU");
	NCM_id inner2 = NCM_id("222","GCGC");
	p="cWW";
	e1 = double(par.junction_energy(outer2,inner2,false));
	e2 = double(par.junction_energy(inner2,outer2,true));
	BOOST_CHECK_CLOSE(e1, e2, 0.000001);
}

BOOST_AUTO_TEST_CASE(NCMs_have_expected_sequence){
    const std::string seq = std::string("GACAAAGUC");
    //test here that "232" and "223" interior and exterior NCMs have correct sequence

}

BOOST_AUTO_TEST_CASE(test_NCM_id_exterior_equals_interior){
    const std::string seq = std::string("GACAAAGUC");
    NCM_type stk = std::string("222");
    NCM_type blg = std::string("232");
	pairtype p = "cWW";
    NCM_id inner = NCM_id(1,7,stk,seq);
    NCM_id outer = NCM_id(6,2,stk,seq);
    BOOST_CHECK_EQUAL(inner.getID(),outer.getID());
    NCM_id inner2 = NCM_id(0,8,blg,seq);
    NCM_id outer2 = NCM_id(7,2,blg,seq);
    BOOST_CHECK_EQUAL(inner2.getID(),outer2.getID());
    //two problems here: sequence is off by one and runs off and not equal
}
/*
BOOST_AUTO_TEST_CASE(test_inner_equals_outer_for_simple_stack){
    //simple hairpin structure
    const std::string seq = std::string("GACAAAGUC");
    //going inner
    NCM_type stk = std::string("222");
	pairtype p = "cWW";
    NCM_id outer = NCM_id(0,8,stk,seq);
    NCM_id inner = NCM_id(1,7,stk,seq);
    std::cout<<"outer "<<outer.getID()<<std::endl;
    std::cout<<"innner "<<inner.getID()<<std::endl;
    double inner_result = double(par.energy(outer,inner,p,false));
    NCM_id outer2 = NCM_id(6,2,stk,seq);
    NCM_id inner2 = NCM_id(7,1,stk,seq);
    std::cout<<"outer "<<outer2.getID()<<std::endl;
    std::cout<<"innner "<<inner2.getID()<<std::endl;
    double outer_result = double(par.energy(outer2,inner2,p,true));
    BOOST_CHECK_CLOSE(inner_result,outer_result,1e-12);

}
*/

