/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to calculate a *witness complex* from an
  unstructured point cloud (using Euclidean distances) and calculate its
  persistent homology.

  Demonstrated classes:

    - aleph::PersistenceDiagram
    - aleph::containers::PointCloud
    - aleph::geometry::FLANN
    - aleph::geometry::BruteForce

  Demonstrated functions:

    - aleph::geometry::buildWitnessComplex
    - aleph::calculatePersistenceDiagrams
    - aleph::utilities::convert

  Original author: Bastian Rieck
*/

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/DowkerComplexWLGraph.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/utilities/String.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <iostream>
#include <string>
#include <vector>

#include <getopt.h>

using namespace aleph;
using namespace containers;
using namespace geometry;
using namespace distances;

void usage()
{
  std::cerr << "Usage: witness_complex FILE [DIMENSION]\n"
            << "\n"
            << "Calculates the witness complex of an unstructured point cloud,\n"
            << "stored in FILE. Euclidean distances are used for the expansion\n"
            << "process. Other optional parameters can be adjusted in order to\n"
            << "change the complex that is built. An optional second argument,\n"
            << "indicating the DIMENSION, can be used to truncate the complex,\n"
            << "making it easier to handle.\n"
            << "\n";
}
typedef double T;
int main( int argc, char** argv )
{
  srand(1);

using Matrix = std::vector< std::vector<T> >;
  //1) set up graph
  Matrix X( 3, std::vector<T>( 3 ) );

  X[0][1] = T(6);
  X[0][2] = T(4);
  X[1][0] = T(1);
  X[1][2] = T(5);
  X[2][0] = T(3);
  X[2][1] = T(3);

  //2) find landmarks if are not given

  //3) compute admissible pairs
  auto X_pairs = aleph::geometry::admissiblePairs( X, T(6) );

  //4) compute dowker complexes
  auto dowkerComplexes = aleph::geometry::buildDowkerComplexGraph<unsigned, T>( X_pairs );

  //5) persistence Diagrams..
  auto diagrams = aleph::calculatePersistenceDiagrams(dowkerComplexes);
  for( auto&& D : diagrams )
    D.removeDiagonal();


  using DataType   = double;
  using Distance   = Euclidean<DataType>;

  double  landmarksFraction = 0.9;
  DataType radius           = DataType();
  bool randomLandmarks      = false;

  std::string input = argv[optind++];

  // This loads the point cloud from an unstructured file. The point
  // cloud loader is smart enough to handle things such as different
  // separators in a file.
  auto pointCloud   = aleph::containers::load<DataType>( input );
  auto dimension    = static_cast<unsigned>( pointCloud.dimension() + 1 );
  auto numLandmarks = static_cast<std::size_t>( static_cast<double>( pointCloud.size() ) * landmarksFraction );

  if( ( argc - optind ) >= 1 )
    dimension = static_cast<unsigned>( std::stoul( argv[optind++] ) );

  std::vector<std::size_t> landmarks;

  if( randomLandmarks ){
    std::cerr << "* Generating landmarks using random strategy...";

    generateRandomLandmarks( pointCloud.size(),
                             numLandmarks,
                             std::back_inserter( landmarks ) );

    std::cerr << "finished\n";
  }
  else
  {
    std::cerr << "* Generating landmarks using max--min strategy...";

    generateMaxMinLandmarks( pointCloud,
                             numLandmarks,
                             std::back_inserter( landmarks ), Distance() );

    std::cerr << "finished\n";
  }

  std::cerr << "* Calculating witness complex with " << ", R=" << radius << ", and d=" << dimension << "...";

  
//  auto K = buildDowkerComplexEuclidean<Distance>( pointCloud, landmarks.begin(), landmarks.end(), dimension, radius );
//
//
//  std::cerr << "finished\n"
//            << "* Obtained simplicial complex with " << K.size() << " simplices\n";
//  std::cerr << "* Calculating persistence diagrams...";
//
//  // Finally, this function will calculate all persistence diagrams of
//  // the simplicial complex. Again, this is a convenience function which
//  // assumes that the complex is already in filtration order.
//  auto diagrams = aleph::calculatePersistenceDiagrams( K );
//
//  std::cerr << "finished\n"
//            << "* Obtained " << diagrams.size() << " persistence diagrams\n";
//
//  for( auto&& D : diagrams )
//  {
//    // Removes all features of zero persistence. They only clutter up
//    // the diagonal.
//    D.removeDiagonal();
//    // This output contains a sort of header (in gnuplot style) so that
//    // it is possible to store multiple persistence diagrams in the same
//    // file.
//    //
//    // Note that for the output, it would also be possible just to loop
//    // over the individual points of the persistence diagram.
//    std::cout << "# Persistence diagram <" << input << ">\n"
//              << "#\n"
//              << "# Dimension: " << D.dimension() << "\n"
//              << "# Entries  : " << D.size() << "\n"
//              << D << "\n\n";
//  }
}
