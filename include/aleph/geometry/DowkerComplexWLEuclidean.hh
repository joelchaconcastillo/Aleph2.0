#ifndef ALEPH_GEOMETRY_WITNESS_COMPLEXLW_HH__
#define ALEPH_GEOMETRY_WITNESS_COMPLEXLW_HH__
#include <algorithm>
#include <limits>
#include <iterator>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>
#include <iostream> //TODO: temporal its for developing..
#include <aleph/math/Combinations.hh>
#include <aleph/geometry/RipsExpander.hh>

#include <aleph/geometry/distances/Traits.hh>

#include <aleph/math/SymmetricMatrix.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>

namespace aleph
{

  namespace geometry
  {
     namespace detail{
        // FIXME: the weight type of an edge should be configurable as an
        // additional template parameter
        using EdgeWeightProperty = boost::property<boost::edge_weight_t, double>;
        
        // FIXME: the data type of the graph should be configurable as an
        // additional template parameter.
        using Graph = boost::adjacency_list<
          boost::vecS,
          boost::vecS,
          boost::directedS,
          boost::no_property,
          EdgeWeightProperty
        >;

        using VertexDescriptor = boost::graph_traits<Graph>::vertex_descriptor;
        
        template <class T, class I = std::size_t> struct Pair
        {
          I p; // first index
          I q; // second index
          T w; // weight
        };
        
        template <class D, class V> struct Vertex
        {
          V p; // vertex index
          D w; // weight
        };
        
        template <class D, class V> void swap( Vertex<D, V>& first, Vertex<D, V>& second )
        {
          std::swap( first.p, second.p );
          std::swap( first.w, second.w );
        }
     }

     template < class Distance, class Container, class InputIterator > 
     topology::SimplicialComplex< topology::Simplex<typename Distance::ResultType, typename std::iterator_traits<InputIterator>::value_type > >
     buildDowkerComplexEuclidean( const Container& container, InputIterator begin, InputIterator end, unsigned dimension = 0, typename Distance::ResultType R = typename Distance::ResultType(), Distance = Distance() )
     {
	     using namespace detail;
             using IndexType         = typename std::iterator_traits<InputIterator>::value_type;
             using VertexType        = IndexType;
             using DataType          = typename Distance::ResultType;
             using Simplex           = topology::Simplex<DataType, VertexType>;
             using SimplicialComplex = topology::SimplicialComplex<Simplex>;
	     using namespace std;
             using Vertex = Vertex<DataType, VertexType>;

	     //  1) First create a set of admissible pairs for a Euclidean space, which can be O(N x n) where N is for the number of points, and n is for the number of landmarks..
	     vector<IndexType> landmarksIndices(begin, end);
	     auto n = landmarksIndices.size();
	     auto N = container.size(); //witness...
	     auto d = container.dimension();
	     Distance dist;
             std::unordered_map< VertexType, std::vector<Vertex> > BasePointMap; //Adjacency  list

	     for( size_t j = 0; j < N ; j++){ // for each witness point...
		vector<DataType> distances;
		distances.reserve(n);
		auto && point = container[j];
		for( size_t i = 0 ; i < n; i++){ //for each landmark..
		   auto && landmark =  container[landmarksIndices.at(i)];
		   auto distancePoints = dist(landmark.begin(), point.begin(), d );
		   if(distancePoints  <= R){
		      BasePointMap[VertexType(j)].push_back({VertexType(i), distancePoints});
		   }
		}
	     }
	     //
	     //
	     //  2) create a combination for each admisible pair, where each admisible pair consists of (w, l), stands for w (witness), l (landmarks)
	     //     --- for each combination make sure of removing witness points
	     //     --
	   using Iterator = typename std::vector<Vertex>::const_iterator;
          auto getWeight = [] ( Iterator begin, Iterator end ){
              DataType weight = std::numeric_limits<DataType>::lowest();
          
              for( auto it = begin; it != end; ++it )
                weight = std::max( weight, it->w);
          
              return weight;
            };
           auto makeSimplices = [&dimension, &getWeight] ( const std::unordered_map< VertexType, std::vector<Vertex> >& map ){
               std::vector<Simplex> simplices;
               std::unordered_map<Simplex, DataType> simplex_to_weight;
               for( auto&& pair : map ){ //for each pair-shortest path
                   auto vertices = pair.second;
                   std::size_t maxDimension = 0;
                   if( dimension == 0 )
                     maxDimension = vertices.size();
                   else
                     maxDimension = dimension + 1;
                   using DifferenceType = typename decltype(vertices)::difference_type;
       
                   for( std::size_t d = std::min( vertices.size(), maxDimension ); d >= 1; d-- ) //for each dimension...
                   {
                       math::for_each_combination( vertices.begin(), vertices.begin() + DifferenceType(d), vertices.end(),
       	                 [&simplex_to_weight, &getWeight] ( Iterator first, Iterator last )
       	                 {
       	                   std::vector<VertexType> vertices_;
       	                   vertices_.reserve( typename std::vector<VertexType>::size_type( std::distance( first, last ) ) );

       	                   for( auto it = first; it != last; ++it )
       	                     vertices_.push_back(it->p);

       	                   Simplex s( vertices_.begin(), vertices_.end() );

       	                   if( simplex_to_weight.find(s) == simplex_to_weight.end() )
       	                     simplex_to_weight[s] = getWeight(first, last);
       	                   else
       	                     simplex_to_weight[s] = std::min( simplex_to_weight[s], getWeight(first, last) );

       	                   return false;
       	                 }
       	               );

             } //end dimenstion.
           } //end shortest path
           for( auto&& pair : simplex_to_weight )
           {
	    // std::cout<<pair.second<<"--->";
	    // for(auto i:pair.first)std::cout<<i<<",";
	    // std::cout<<"\n";
             auto s = pair.first;
             s.setData( pair.second );
             simplices.push_back( s );
           }
           return simplices;
         }; 
             std::vector<Simplex> simplices;
             auto Edges = makeSimplices( BasePointMap );
             SimplicialComplex dowkerComplex( Edges.begin(), Edges.end() );
             dowkerComplex.sort( topology::filtrations::Data<Simplex>() );
	 return dowkerComplex;
     }

    template <class T, class OutputIterator> void generateRandomLandmarks( T n, T k, OutputIterator result )
    {
      std::random_device rd;
      std::mt19937 rng( rd() );
    
      std::vector<T> indices( n );
    
      using DifferenceType = typename std::vector<T>::difference_type;
    
      std::iota( indices.begin(), indices.end(), T() );
      std::shuffle( indices.begin(), indices.end(), rng );
      std::copy( indices.begin(), indices.begin() + static_cast<DifferenceType>(k), result );
    }

    template <
      class Distance,
      class Container,
      class OutputIterator
    > void generateMaxMinLandmarks( const Container& container, std::size_t n, OutputIterator result, Distance distance = Distance() )
    {
      if( n > container.size() )
        throw std::out_of_range( "Number of landmarks is out of range" );
    
      using SizeType = decltype( container.size() );
    
      std::random_device rd;
      std::mt19937 rng( rd() );
    
      std::uniform_int_distribution<SizeType> distribution( SizeType(0), container.size() - 1 );
    
      std::vector<SizeType> indices;
      indices.reserve( n );
    
//      indices.emplace_back( distribution( rng ) );
      indices.emplace_back(0);
    
      using DataType = typename Distance::ResultType;
      auto N         = container.size();
      auto d         = container.dimension();
    
      while( indices.size() < n )
      {
        auto index = SizeType(0);
        auto max   = std::numeric_limits<DataType>::lowest();
    
        for( SizeType i = 0; i < N; i++ )
        {
          auto min = std::numeric_limits<DataType>::max();
    
          for( auto&& landmarkIndex : indices )
          {
            auto dist = distance( container[i].begin(), container[landmarkIndex].begin(), d );
            min       = std::min( min, dist );
          }
    
          if( min > max )
          {
            max   = min;
            index = i;
          }
        }
    
        indices.push_back( index );
      }
    
      std::copy( indices.begin(), indices.end(), result );
    }
  }
}
#endif