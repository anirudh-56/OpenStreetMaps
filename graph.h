// graph.h <Starter Code>
// application.cpp <Starter Code>
// Anirudh Yallapragada (ayall2@uic.edu)
//
// Adam T Koehler, PhD
// University of Illinois Chicago
// CS 251, Fall 2023
//
// Project Original Variartion By:
// Joe Hummel, PhD
// University of Illinois at Chicago
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>

using namespace std;


// each edge is a pair (vertex to, weight)
template<typename VertexT, typename WeightT>
class graph {
 private:

  static constexpr int ListSize = 100;
  // adjacency_list for a weighted directed graph:
  // maps each vertex to the collection of edges radiating out of that vertex
  map<VertexT, map<VertexT, WeightT>> adjList;


 public:
  //
  // constructor:
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  graph(){
      // Initialize an empty graph
  }

  //
  // NumVertices
  // Returns the # of vertices currently in the graph.
  int NumVertices() const {
    //Count number of vertices within the map

    return adjList.size();
  }

  //
  // NumEdges
  // Returns the # of edges currently in the graph.
  int NumEdges() const {
    int count = 0;

    //
    // loop through the adjacency list and count how many
    // edges currently exist:

    for(const auto& vertex: adjList){
      count += vertex.second.size();
    }

    return count;

  }


  // addVertex
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  // This needs to be an O(lgN)
  bool addVertex(VertexT v) {
    // is the vertex already in the graph?  If so, we do not
    // insert again otherwise Vertices may fill with duplicates:
    if(adjList.find(v) != adjList.end()){
      return false;
    }

    // if we get here, vertex does not exist so insert.  Where
    // we insert becomes the rows and col position for this
    // vertex in the adjacency matrix.
    adjList[v] = map<VertexT, WeightT>(); //Adding an empty map as the value to the key Vertex
    return true;
  }


  // addEdge
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    //
    // Search the Vertices and finding the position
    // of each vertex; this will denote the row and col to
    // access in the adjacency list:

    bool fromExists = adjList.find(from) != adjList.end();
    bool toExists = adjList.find(to) != adjList.end();

    if(fromExists && toExists){ //If the vertex is already in the list
      adjList[from][to] = weight; //Add associataed weight to it
      return true;
    }
    else{
      return false; //Not adding a new element 
    }
  }


  // getWeight
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    //
    // we need to search the Vertices and find the position
    // of each vertex; this will denote the row and col to
    // access in the adjacency matrix:
    //
    bool fromExists = adjList.find(from) != adjList.end();
    bool toExists = adjList.find(to) != adjList.end();
    if(!fromExists || !toExists){ //If the vertex is already in the list
      return false;
    }

    auto itFrom = adjList.find(from); //itFrom will be the outer map 
    const auto& vertex = itFrom->second; //vertex will be the inner map
    auto neighbor = vertex.find(to); //Searching for the 'to' vertex in the inner map

    if (neighbor == vertex.end()) {  // no 'to' vertex is found, so there's no weight is found:
      return false;
    }

    // Okay, the edge exists, return the weight via the
    // reference parameter:
    weight = neighbor->second; //The weight will be the value from the inner map
    return true;
  }


  // neighbors
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT>  S;


    // we need to search the Vertices and find the position of v
    auto vertices = adjList.find(v);
    if (vertices == adjList.end()) {  // not found:
      return S;
    }

    // we found the row, so loop along the row and for every
    // edge that exists, add the column vertex to V:
    // NOTE: how many columns are there?  The # of vertices.
    for(const auto& neighbor: vertices->second){
      S.insert(neighbor.first);
    }

    return S;
  }


  // getVertices
  vector<VertexT> getVertices() const {

    vector<VertexT> tempV;


    for(const auto& vertex: adjList){
      tempV.push_back(vertex.first); 
    }

    return tempV; //Returns a vector containing all the vertices currently in the graph.
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;

    int i = 0;
    for(const auto& vertex: adjList){
      output << " " << i << "." << vertex.first << endl;
      i++;
    }

    output << endl;
    output << "**Edges:" << endl;
    for(const auto& vertex: adjList){
      output << " row " << vertex.first << ": ";

      for(const auto& neighbor: vertex.second){
        if(vertex.second.find(neighbor.first) != vertex.second.end()){
          output << "(T," << neighbor.second << ") ";
        }
        else{
          output << "F ";
        }
      }
      output << endl;
    }
    output << "**************************************************" << endl;
  }
};
