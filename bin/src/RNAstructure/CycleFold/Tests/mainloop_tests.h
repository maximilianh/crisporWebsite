#include "../mainloop.h"

BOOST_AUTO_TEST_CASE(wrap_works_correctly){
	sequence s = sequence(std::string("foo"),std::string("AAAGGGUUUCCC"));
	const int length = s.getLength();
	BOOST_CHECK(wrap(2,length)==2);
	BOOST_CHECK(wrap(length,length)==0);
	BOOST_CHECK(wrap(length+2,length)==2);
}

BOOST_AUTO_TEST_CASE(ncm_length_works_correctly){
    BOOST_CHECK_EQUAL(ncm_length(std::string("222")),4);
    BOOST_CHECK_EQUAL(ncm_length(std::string("13")),3);
}

BOOST_AUTO_TEST_CASE(too_close_works_correctly){
    const int length = 9;
    NCM_type outer = std::string("222");
    NCM_type inner = std::string("222");
    NCM_type hairpin = std::string("13");
    //test too_close for stacks
    //adding 222 onto 222 on the interior
    BOOST_CHECK(too_close(length,3,5,outer,inner));
    BOOST_CHECK(!too_close(length,1,7,outer,inner));
    //adding 222 onto 222 on the exterior
    BOOST_CHECK(too_close(length,8,0,outer,inner));
    BOOST_CHECK(!too_close(length,6,2,outer,inner));
    //adding 222 onto 13 on the interior
    BOOST_CHECK(too_close(length,3,5,hairpin,inner));
    BOOST_CHECK(!too_close(length,2,6,hairpin,inner));

    //test too_close for multibranch loops
    //double NCM in inretior calculation case
    BOOST_CHECK(too_close(length,3,7,inner,std::string("mb")));
    BOOST_CHECK(!too_close(length,2,7,inner,std::string("mb")));
    //double NCM in exterior calculation case
    BOOST_CHECK(too_close(length,7,1,inner,std::string("mb")));
    BOOST_CHECK(!too_close(length,6,2,inner,std::string("mb")));
    //single_NCM case
    //too far apart
    BOOST_CHECK(too_close(length,3,6,hairpin,std::string("mb")));
    //doesn't leave a fragment on the end
    BOOST_CHECK(too_close(length,0,2,hairpin,std::string("mb")));
    //fragment on the end and the correct distance apart
    BOOST_CHECK(!too_close(length,4,6,hairpin,std::string("mb")));

    //test too_close for exterior loops
    BOOST_CHECK(too_close(length,8,0,std::string("221"),std::string("ext")));
    BOOST_CHECK(too_close(10,8,0,std::string("221"),std::string("ext")));
    BOOST_CHECK(!too_close(length,7,1,std::string("222"),std::string("ext")));
}

BOOST_AUTO_TEST_CASE(set_di_dj_works_correctly){
    const int i = 1;
    const int j = 7;
    int di = 0;
    int dj = 0;
    set_di_dj(i,j,di,dj,std::string("232"));
    BOOST_CHECK_EQUAL(di,2);
    BOOST_CHECK_EQUAL(dj,1);
    set_di_dj(i,j,di,dj,std::string("223"));
    BOOST_CHECK_EQUAL(di,1);
    BOOST_CHECK_EQUAL(dj,2);
    set_di_dj(j,i,di,dj,std::string("223"));
    BOOST_CHECK_EQUAL(di,2);
    BOOST_CHECK_EQUAL(dj,1);
}

BOOST_AUTO_TEST_CASE(folding_works){
	sequence testSeq(std::string("test"), std::string("GCGAAACGC"));
	parameters<int> par = parameters<int>(options());
	arrays<int> arr = fill_arrays(testSeq,par,constraints());
    std::vector<std::pair<int,int>> pairs = traceback<int>(testSeq, arr, par, options());
}

BOOST_AUTO_TEST_CASE(probs_works){
	sequence testSeq(std::string("test"), std::string("GCGAAACGC"));
	parameters<real_t> par = parameters<real_t>(options());
	table_t probs = calculate_pairing_probabilities(testSeq,par,constraints());
	for(int i=0;i<testSeq.getLength();i++)
	for(int j=0;j<testSeq.getLength();j++){
		BOOST_CHECK_CLOSE(probs[i][j], probs[j][i], 0.0000001);
		BOOST_ASSERT(probs[i][j] <= 1.0);
	}
}
