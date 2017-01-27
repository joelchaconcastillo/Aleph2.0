#include "math/StepFunction.hh"

#include "tests/Base.hh" 

using namespace aleph::math;

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Step function: Basic properties" );

  StepFunction<T> f;
  f.add( 0, 1, 1   );
  f.add( 2, 3, 1   );
  f.add( 3, 4, 2 );

  ALEPH_ASSERT_THROW( f(0)   == 1 );
  ALEPH_ASSERT_THROW( f(1)   == 1 );
  ALEPH_ASSERT_THROW( f(1.5) == 0 );
  ALEPH_ASSERT_THROW( f(2)   == 1 );
  ALEPH_ASSERT_THROW( f(3)   == 1 ); // TODO: Incorrect
  ALEPH_ASSERT_THROW( f(3.5) == 2 );

  ALEPH_TEST_END();
}

int main()
{
  test<double>();
  test<float> ();
}