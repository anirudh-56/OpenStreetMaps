// application.cpp <Starter Code>
// Anirudh Yallapragada (ayall2@uic.edu)
//
//
// Adam T Koehler, PhD
// University of Illinois Chicago
// CS 251, Fall 2023
//
// Project Original Variartion By:
// Joe Hummel, PhD
// University of Illinois at Chicago
//
// 
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <map>
#include <queue>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "tinyxml2.h"
#include "dist.h"
#include "graph.h"
#include "osm.h"


using namespace std;
using namespace tinyxml2;

//
// Implement your standard application here
//
const double INF = numeric_limits<double>::max();

class prioritize  // you could also use a struct
{
public:
  bool operator()(const pair<long long, double>& p1, const pair<long long, double>& p2) const
  {
    return p1.second > p2.second; 
  }
};

vector<long long> findPath(map<long long, long long>& predecessor, long long endVertex){
    long long currentVertex = endVertex;
    stack<long long> path;
    vector<long long> finalPath;

    while(predecessor[currentVertex] != -1){
      //Vertices are going into a stack so the order can be reveresed
      path.push(currentVertex);
      currentVertex = predecessor[currentVertex];
    }

    path.push(currentVertex);

    //Reverse the order of vertices
    while(!path.empty()){
      currentVertex = path.top();
      path.pop();
      finalPath.push_back(currentVertex);
    }

    return finalPath; //Returning a vector of the path
}

//Startvertex (long long), graph, two maps(passed by referencs, distances, and predessors
void DijkstraShortestPath(graph<long long, double>&G, long long& nearestNode, map<long long, double>&nodeDistances, map<long long, long long>& predecessors){
  priority_queue<pair<long long, double>, vector<pair<long long, double>>, prioritize> unvisitedQueue;
  set<long long> visitedNodes;

  //each vertex currentV in graph
   for (auto& vertex: G.getVertices()) {
      nodeDistances[vertex] = INF;
      predecessors[vertex] = -1;
      //Enqueue currentV in unvisitedQueue
      if (vertex != nearestNode) {
        unvisitedQueue.push(pair<long long, double>(vertex, INF));
      } 
   }

   //startV has a distance of 0 from itself
   //startV⇢distance = 0
   nodeDistances[nearestNode] = 0.00;
   unvisitedQueue.push(pair<long long, double>(nearestNode, 0));

   while (!unvisitedQueue.empty()) {
      // Visit vertex with minimum distance from startV
      pair<long long, double> currentVertex = unvisitedQueue.top();
      unvisitedQueue.pop();
      //currentV = DequeueMin unvisitedQueue
      if(nearestNode == INF){
        break;
      }
      
      if(visitedNodes.count(currentVertex.first)){
        continue;
      }
      
      visitedNodes.insert(currentVertex.first);

      //Check if we have the shortest path to V
      //for each vertex adjV adjacent to currentV
      for (auto& neighbor: G.neighbors(currentVertex.first)) {
         double edgeWeight;
         G.getWeight(currentVertex.first, neighbor, edgeWeight);
         double alternativePathDistance = nodeDistances[currentVertex.first] + edgeWeight;
            
         // If shorter path from startV to adjV is found,
         // update adjV's distance and predecessor
         //cout << alternativePathDistance << endl;
         if (alternativePathDistance < nodeDistances[neighbor]) {
            nodeDistances[neighbor] = alternativePathDistance;
            predecessors[neighbor] = currentVertex.first;

            unvisitedQueue.push(pair<long long, double>(neighbor, alternativePathDistance));
         }
      }
   }
}

Coordinates nearestBuilding(BuildingInfo building, map<long long, Coordinates> Nodes, vector<FootwayInfo> Footways){
  double min = INF;

  Coordinates closestCoord;

  for(const auto& footways: Footways){ 
    //Looping through the nodes of the current footway
    for(long unsigned i = 0; i < footways.Nodes.size(); i++){
      long long node1 = footways.Nodes[i];

      double distance = distBetween2Points(building.Coords.Lat, building.Coords.Lon, Nodes[node1].Lat, Nodes[node1].Lon);
      //Distance between building's coordinates and the current Cordinates

      if(distance < min){ 
        min = distance;
        closestCoord.Lat = Nodes[node1].Lat; 
        closestCoord.Lon = Nodes[node1].Lon;
        closestCoord.ID = node1;
      }
      //Getting the cloeset coordinates
    }
  } 

  return closestCoord; 
}

bool searchBuilding(vector<BuildingInfo>& Buildings, string query, BuildingInfo& info){
  for(int i = 0; i < Buildings.size(); i++){ //Searching the vector of Buildings
    if(query == Buildings[i].Abbrev){ //If the full string input is the title of the Building 
      info = Buildings[i]; //Info is now the Building
      return true;
    }
    else if (Buildings[i].Fullname.find(query) != string::npos){ //If they entered the accrynonym
      info = Buildings[i]; //Info is now the Building
      return true;
    }
  }

  return false; //If the building was not found
}

void application(
    map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings, graph<long long, double>& G) {
  string person1Building, person2Building;

  map<long long, double> nodeDistancesOne; 
  map<long long, double> nodeDistancesTwo;
  map<long long, long long> predecessorsOne; 
  map<long long, long long> predecessorsTwo; 
  stack<long long> personOne;
  stack<long long> personTwo;
  set<string> unreachableBuilding;
  //For Dijkstra's

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);

    bool foundBuilding1 = false;
    bool foundBuilding2 = false;
   
    //Milestone 7: Search Buildings 1 and 2
    //Searching for the buildings, using a boolean to determine if they are found
    BuildingInfo info1;
    BuildingInfo info2;
    foundBuilding1 = searchBuilding(Buildings, person1Building, info1);
    foundBuilding2 = searchBuilding(Buildings, person2Building, info2);
    bool exitSearch = false;
    //A boolean function to determine if the Building exists, if it is then the Building info is transferred to the info

    while(!foundBuilding1 || !foundBuilding2){ //Nested while loop to constatntly take in inputs if the inputed Buildings were incorrect
      if(!foundBuilding1){
        cout << "Person 1's building not found" << endl;
      }
      if(!foundBuilding2){
        cout << "Person 2's building not found" << endl;
      }

      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      foundBuilding1 = searchBuilding(Buildings, person1Building, info1);

      if(person1Building == "#"){ //Breaks out of current while loop
        exitSearch = true;
        break;
      }

      cout << "Enter person 2's building (partial name or abbreviation)> ";
      getline(cin, person2Building);
      foundBuilding2 = searchBuilding(Buildings, person2Building, info2);  
    }

    if(exitSearch){ //Breaks out of main while loop 
      break;
    }

    cout << "Person 1's point:" << endl;
    cout << " " << info1.Fullname << endl;
    cout << " (" << info1.Coords.Lat << ", " << info1.Coords.Lon << ")" << endl;

    cout << "Person 2's point:" << endl;
    cout << " " << info2.Fullname << endl;
    cout << " (" << info2.Coords.Lat << ", " << info2.Coords.Lon << ")" << endl;
    //Prining the Building information

    bool pathFound = false;
    

    if(foundBuilding1 && foundBuilding2){
        while(!pathFound){
        
        //Milestone 8 (Locate the Center Building)
        Coordinates midpoint;
        BuildingInfo midBuilding;
      
        midpoint = centerBetween2Points(info1.Coords.Lat, info1.Coords.Lon, info2.Coords.Lat, info2.Coords.Lon); //The center coordinates

        double min = INF;
        for(int i = 0; i < Buildings.size(); i++){
          double distance = distBetween2Points(midpoint.Lat, midpoint.Lon, Buildings[i].Coords.Lat, Buildings[i].Coords.Lon); 

          //Checking the building is unreachable
          if(unreachableBuilding.find(Buildings[i].Fullname) != unreachableBuilding.end()){
            continue;
          }

          //Getting the minimum distance for the closest building
          if(distance < min){
            min = distance;
            midBuilding = Buildings[i];
          }
        }

        cout << "Destination Building:" << endl;
        cout << " " << midBuilding.Fullname << endl;
        cout << " (" << midBuilding.Coords.Lat << ", " << midBuilding.Coords.Lon << ")" << endl;        
        
        //Milestone 9 (Find Nearest Nodes from buildings 1, 2, & Center)
        Coordinates info1Coords = nearestBuilding(info1, Nodes, Footways);
        Coordinates info2Coords = nearestBuilding(info2, Nodes, Footways);
        Coordinates centerCoords = nearestBuilding(midBuilding, Nodes, Footways);
        //Finding nearest nodes from the buildings
        
        cout << endl;
        cout << "Nearest P1 node:" << endl;
        cout << " " << info1Coords.ID << endl; 
        cout << "(" << info1Coords.Lat << ", " << info1Coords.Lon << ")" << endl; //Cordinates 

        cout << "Nearest P2 node:" << endl;
        cout << " " << info2Coords.ID << endl; 
        cout << "(" << info2Coords.Lat << ", " << info2Coords.Lon << ")" << endl; //Cordinates 

        cout << "Nearest destination node:" << endl;
        cout << " " << centerCoords.ID << endl; 
        cout << "(" << centerCoords.Lat << ", " << centerCoords.Lon << ")" << endl; //Cordinates       

        //Milestone 10 (run Dijkstra's alg from each start)
        DijkstraShortestPath(G, info1Coords.ID, nodeDistancesOne, predecessorsOne); 
        DijkstraShortestPath(G, info2Coords.ID, nodeDistancesTwo, predecessorsTwo);
        //Using Dijkstra's algorithm to find then shortest path between maps

        
        //Milestone 11 (Print path, (path found! break))
        if(nodeDistancesOne[info2Coords.ID] >= INF){ 
          cout << endl;
          cout << "Sorry, destination unreachable" << endl;
          break;
        } 
        //Determining if the destination is unreachable

        else if(nodeDistancesOne[centerCoords.ID] >= INF || nodeDistancesTwo[centerCoords.ID] >= INF){ //If either one is greater than Infinity
          unreachableBuilding.insert(midBuilding.Fullname);
          cout << "At least one person was unable to reach the destination building" << endl;
          continue;
        }
        //If one of the people's distance was greater than Infinity, then one of them are unreachable

        cout << endl;
        cout << "Person 1's distance to dest: " << nodeDistancesOne[centerCoords.ID] << " miles" << endl;
        cout << "Path: ";
        vector<long long> pathOne = findPath(predecessorsOne, centerCoords.ID); //Function returns a vector of IDs for me to loop through
        for(long long vertex: pathOne){
          cout << vertex;
          if(vertex != pathOne.back()){
            cout << "->";
          }
        }
        cout << endl;
        //Outputting the person 1's path to the distnation 

        cout << endl;

        cout << "Person 2's distance to dest: " << nodeDistancesTwo[centerCoords.ID] << " miles" << endl;
        cout << "Path: ";
        vector<long long> pathTwo = findPath(predecessorsTwo, centerCoords.ID);
        for(long long vertex: pathTwo){
          cout << vertex;
          if(vertex != pathTwo.back()){
            cout << "->";
          }
        }
        cout << endl;
        //Outputting the person 2's path to the distnation 

        pathFound = true; //Path was found, we can exit the loop
        break;
      }
    }

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }    
}

int main() {
  graph<long long, double> G;

  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }


  // Load XML-based map file
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }


  // Read the nodes, which are the various known positions on the map:
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  // Read the footways, which are the walking paths:
  int footwayCount = ReadFootways(xmldoc, Footways);

  // Read the university buildings:
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  // Stats
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;



  //Milestone 5 (Add the Vertices)
  for(const auto& vertex: Nodes){ 
    G.addVertex(vertex.first); //Making use of funnction from graph.h
  }

  //Milestone 6 (Add the Edges)
  for(const auto& Footway: Footways){
    for(int i = 0; i < Footway.Nodes.size() - 1; i++){
      long long node1 = Footway.Nodes[i];
      long long node2 = Footway.Nodes[i + 1];
      double distance = distBetween2Points(Nodes[node1].Lat, Nodes[node1].Lon, Nodes[node2].Lat, Nodes[node2].Lon);

      G.addEdge(Footway.Nodes.at(i), Footway.Nodes.at(i + 1), distance); //Making use of funnction from graph.h
      G.addEdge(Footway.Nodes.at(i + 1), Footway.Nodes.at(i), distance);
      //bidirectional edge
    }
  }

  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  // Execute Application
  application(Nodes, Footways, Buildings, G);


  // done:
  cout << "** Done **" << endl;
  return 0;
}
