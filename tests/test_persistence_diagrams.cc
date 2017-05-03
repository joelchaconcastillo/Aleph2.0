#include "tests/Base.hh"

#include "distances/Wasserstein.hh"

#include "persistenceDiagrams/Mean.hh"
#include "persistenceDiagrams/MultiScaleKernel.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <limits>
#include <random>
#include <vector>

#include <cmath>

template <class T> aleph::PersistenceDiagram<T> createRandomPersistenceDiagram( unsigned n )
{
  std::random_device rd;
  std::default_random_engine rng( rd() );
  std::uniform_real_distribution<T> distribution( T(0), T( std::nextafter( T(1), std::numeric_limits<T>::max() ) ) );

  using PersistenceDiagram = aleph::PersistenceDiagram<T>;
  PersistenceDiagram D;

  for( unsigned i = 0; i < n; i++ )
  {
    auto x = distribution( rng );
    auto y = distribution( rng );

    if( x > y )
      std::swap( x,y );

    D.add( x,y );
  }

  return D;
}

template <class T> void testMultiScaleKernel()
{
  ALEPH_TEST_BEGIN( "Multi-scale kernel" );

  auto D1 = createRandomPersistenceDiagram<T>( 50 );
  auto D2 = createRandomPersistenceDiagram<T>( 50 );

  auto d1 = aleph::multiScalePseudoMetric( D1, D2, 1.0 );
  auto d2 = aleph::multiScalePseudoMetric( D1, D2, 2.0 );
  auto d3 = aleph::distances::wassersteinDistance( D1, D2, T(1) );

  // Non-negativity
  ALEPH_ASSERT_THROW( d1 > 0.0 );
  ALEPH_ASSERT_THROW( d2 > 0.0 );
  ALEPH_ASSERT_THROW( d3 > 0.0 );

  // Symmetry
  ALEPH_ASSERT_EQUAL( aleph::multiScalePseudoMetric(D1, D1, 1.0), aleph::multiScalePseudoMetric(D2, D2, 1.0) );
  ALEPH_ASSERT_EQUAL( aleph::multiScalePseudoMetric(D1, D1, 1.0), 0.0 );

  // Stability
  ALEPH_ASSERT_THROW( d1 < 1.0 / ( 1.0 / ( 1.0 * std::sqrt( 8.0 * M_PI ) ) * d3 ) );
  ALEPH_ASSERT_THROW( d2 < 1.0 / ( 1.0 / ( 2.0 * std::sqrt( 8.0 * M_PI ) ) * d3 ) );

  ALEPH_TEST_END();
}

template <class T> void testFrechetMean()
{
  using PersistenceDiagram = aleph::PersistenceDiagram<T>;

  ALEPH_TEST_BEGIN( "Persistence diagram mean");

  unsigned n = 20;

  std::vector<PersistenceDiagram> diagrams;
  diagrams.reserve( n );

  for( decltype(n) i = 0; i < n; i++ )
    diagrams.emplace_back( createRandomPersistenceDiagram<T>( 50 ) );

  auto D = aleph::mean( diagrams.begin(), diagrams.end() );

  ALEPH_ASSERT_THROW( D.size() > 0 );

  std::cerr << D << "\n";

  ALEPH_TEST_END();
}

int main(int, char**)
{
  testMultiScaleKernel<float> ();
  testMultiScaleKernel<double>();

  testFrechetMean<float> ();
  testFrechetMean<double>();
}