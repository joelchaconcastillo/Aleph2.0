#ifndef ALEPH_GEOMETRY_SPHERE_SAMPLING_HH__
#define ALEPH_GEOMETRY_SPHERE_SAMPLING_HH__

#include <cmath>

#include <random>
#include <vector>

namespace aleph
{

namespace geometry
{

/**
  Samples $n$ points from a sphere such that the expected number of
  points per area of the sphere is uniform. Only the angular values
  of the sampled points ($\theta$, $\phi$) will be returned.

  @param n Number of samples to draw

  @returns Vector of angle values, i.e. $\theta$ and $\phi$, which are
  sufficient to describe the sphere. Use aleph::geometry::makeSphere()
  for creating a point cloud from the resulting angles.
*/

template <class T>
std::vector< std::pair<T, T> > sphereSampling( unsigned n )
{
  std::random_device rd;

  std::mt19937 rngU( rd() );
  std::mt19937 rngV( rd() );

  std::vector< std::pair<T, T> > angles;
  angles.reserve( n );

  std::uniform_real_distribution<T> uDistribution( std::nextafter( T(0), std::numeric_limits<T>::max() ), T(1) );
  std::uniform_real_distribution<T> vDistribution( std::nextafter( T(0), std::numeric_limits<T>::max() ), T(1) );

  for( unsigned i = 0; i < n; i++ )
  {
    auto u = uDistribution( rngU );
    auto v = uDistribution( rngV );

    T theta = T( 2*M_PI*u );
    T phi   = std::acos( 2*v - 1 );

    angles.push_back( std::make_pair( theta, phi ) );
  }

  return angles;
}

/**
  Converts a vector of angles into a point cloud that contains samples
  from a sphere of a given radius.
*/

template <class T> aleph::containers::PointCloud<T> makeSphere( const std::vector< std::pair<T, T> >& angles, T r )
{
  using PointCloud = aleph::containers::PointCloud<T>;

  PointCloud pc( angles.size(), 3 );

  unsigned index = 0;
  for( auto&& pair : angles )
  {
    auto theta = pair.first;
    auto phi   = pair.second;

    auto x = r * std::sin( theta ) * std::cos( phi );
    auto y = r * std::sin( theta ) * std::sin( phi );
    auto z = r * std::cos( theta );

    pc.set( index++, {x,y,z} );
  }

  return pc;
}

} // namespace geometry

} // namespace aleph

#endif