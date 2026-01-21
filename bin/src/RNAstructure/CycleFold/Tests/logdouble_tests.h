#include "../logdouble.h"

BOOST_AUTO_TEST_CASE(log_math){
	log_double foo = log_double(5.0);
	log_double bar = log_double(2.0);
	log_double zero = log_double(0.0);
	log_double bazz = log_double(0.25);

    BOOST_CHECK_CLOSE((zero*foo).asDouble(),0.0,1e-20);
    BOOST_CHECK_CLOSE((foo*zero).asDouble(),0.0,1e-20);

//5+2=7
	BOOST_CHECK_CLOSE((foo+bar).asDouble(),7.0,1e-6);
//5-2=3
	BOOST_CHECK_CLOSE((foo-bar).asDouble(),3.0,1e-6);
//5*2=10
	BOOST_CHECK_CLOSE((foo*bar).asDouble(),10.0,1e-6);
//5/2=2.5
	BOOST_CHECK_CLOSE((foo/bar).asDouble(),2.5,1e-6);
//5+0=5
	BOOST_CHECK_CLOSE((foo+zero).asDouble(),foo.asDouble(),1e-6);
//5*0=0
	BOOST_CHECK_CLOSE((foo*zero).asDouble(),zero.asDouble(),1e-6);

//2^3
	BOOST_CHECK_CLOSE((bar^3).asDouble(),8.0,1e-6);
//0.25 ^ 2
	BOOST_CHECK_CLOSE((bazz^2).asDouble(),0.25*0.25,1e-6);

//test in-place mutations
	foo += bar;
	BOOST_CHECK_CLOSE(foo.asDouble(),7.0,1e-6);
	foo *= log_double(5.0);
	BOOST_CHECK_CLOSE(foo.asDouble(), 35.0,1e-6);
	foo -= log_double(5.0);
	BOOST_CHECK_CLOSE(foo.asDouble(), 30.0,1e-6);
	foo /= log_double(6.0);
	BOOST_CHECK_CLOSE(foo.asDouble(), 5.0,1e-6);

//test initialization
    log_double baz = 0.0;
    log_double thud = log_double(0.0);
    BOOST_CHECK_CLOSE(baz.asDouble(),thud.asDouble(),1e-6);
    thud = 2.0;
    BOOST_ASSERT(thud > baz);
    BOOST_CHECK_CLOSE((zero+thud).asDouble(),2.0,1e-6);

}

BOOST_AUTO_TEST_CASE(exponents){
	log_double five = log_double(5.0);
	log_double two = log_double(2.0);
	log_double one = log_double(1.0);
	log_double zero = log_double(0.0);
	BOOST_CHECK_CLOSE(pow(zero, one).asDouble(), 0.0, 1e-20);
	BOOST_CHECK_CLOSE(pow(one, one).asDouble(), 1.0, 1e-20);
	BOOST_CHECK_CLOSE(pow(five, zero).asDouble(), 1.0, 1e-20);
	BOOST_CHECK_CLOSE(pow(five, two).asDouble(), 25.0, 1e-7);
}
