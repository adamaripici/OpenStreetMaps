// application.cpp <Starter Code>
// <Ada Pici>
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
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

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stack>
#include <string>
#include <vector>
#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;
const double INF = numeric_limits<double>::max();

class prioritize {
 public:
  bool operator()(const pair<long long, double>& p1,
                  const pair<long long, double>& p2) const {
    if (p1.second == p2.second) return p1.first > p2.first;
    return p1.second > p2.second;
  }
};

// dijstra algorithm
vector<long long> Dijkstra(graph<long long, double>& G, long long startV,
                           map<long long, double>& distances,
                           map<long long, long long>& prev) {
  long long currentV;

  double alternativePath;
  vector<long long> visited;

  priority_queue<pair<long long, double>, vector<pair<long long, double>>,
                 prioritize>
      unvisitedQueue;

  vector<long long> vertex = G.getVertices();

  for (long long x : vertex) {
    unvisitedQueue.push(make_pair(x, INF));
    distances[x] = INF;
    prev[x] = -1;
  }

  unvisitedQueue.push(make_pair(startV, 0));
  distances[startV] = 0;
  set<long long> visitedSet;

  while (!unvisitedQueue.empty()) {
    currentV = unvisitedQueue.top().first;

    double weight = unvisitedQueue.top().second;
    unvisitedQueue.pop();

    if (distances[currentV] == INF) {
      break;
    } else if (weight == INF) {
      break;
    } else if (visitedSet.count(currentV) > 0) {
      continue;
    } else {
      visited.push_back(currentV);
      visitedSet.insert(currentV);
    }

    set<long long> neighbors = G.neighbors(currentV);

    for (auto& y : neighbors) {
      double edgeWeight;
      G.getWeight(currentV, y, edgeWeight);
      alternativePath = distances[currentV] + edgeWeight;

      if (alternativePath < distances[y]) {
        unvisitedQueue.push(make_pair(y, alternativePath));
        prev[y] = currentV;
        distances[y] = alternativePath;
      }
    }
  }
  return visited;
}

//
// Implement your creative component application here
// TO DO: add arguments
//
void creative() {}

void getPath(map<long long, long long> predecessors, long long start,
             long long destination) {
  stack<long long> mystack;

  long long curr = destination;
  while (curr != start) {
    mystack.push(curr);
    curr = predecessors[curr];
  }

  cout << start;
  while (!mystack.empty()) {
    cout << "->" << mystack.top();
    mystack.pop();
  }

  cout << endl;
}
bool search(string building, vector<BuildingInfo>& Buildings,
            BuildingInfo& info) {
  // input can be abbreviation or partial name
  // first search through abbreviation
  // if none found, search for partial names
  for (auto const& b : Buildings) {
    if (b.Abbrev == building) {
      info = b;
      return true;
    }
    if (b.Fullname.find(building) != string::npos) {
      info = b;
      return true;
    }
  }

  return false;
}

//
// Implement your standard application here
// TO DO: add a parameter for the graph you make.
//
void application(map<long long, Coordinates>& Nodes,
                 vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
                 graph<long long, double> G) {
  string person1Building, person2Building;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);
  bool searching1 = true;
  bool searching2 = true;
  BuildingInfo info1, info2;
  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    // TODO MILESTONE 7: Search Buildings 1 and 2
    // look for building
    searching1 = search(person1Building, Buildings, info1);
    searching2 = search(person2Building, Buildings, info2);
    if (searching1 == false) {
      cout << "Person 1's building not found" << endl;
      cout
          << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    } else if (searching2 == false) {
      cout << "Person 2's building not found" << endl;
      cout
          << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    } else {
      cout << "Person 1's point:" << endl;
      cout << " " << info1.Fullname << endl;
      cout << " (" << info1.Coords.Lat << ", " << info1.Coords.Lon << ")"
           << endl;

      cout << "Person 2's point:" << endl;
      cout << " " << info2.Fullname << endl;
      cout << " (" << info2.Coords.Lat << ", " << info2.Coords.Lon << ")"
           << endl;

      // find the destination building
      Coordinates midpoint;

      // loop through all the buildings
      midpoint = centerBetween2Points(info1.Coords.Lat, info1.Coords.Lon,
                                      info2.Coords.Lat, info2.Coords.Lon);

      // search the building vector for the building
      // w coordinates closest to the center of that line
      double min = INF;
      double distance;
      BuildingInfo destinationBuilding;
      for (auto const& building : Buildings) {
        // calculate distance between midpoint
        // and each building
        distance = distBetween2Points(midpoint.Lat, midpoint.Lon,
                                      building.Coords.Lat, building.Coords.Lon);

        if (distance < min) {
          destinationBuilding = building;
          min = distance;
        }
      }
      cout << "Destination Building:" << endl;
      cout << " " << destinationBuilding.Fullname << endl;
      cout << " (" << destinationBuilding.Coords.Lat << ", "
           << destinationBuilding.Coords.Lon << ")" << endl;
      // calculate distance between midpoint and
      // each building

      // calculate distance between midpoint and
      // each building

      // find nearest nodes
      long long firstStart, firstDest;
      double startMinDist, endMinDist;
      firstStart = Footways[0].Nodes[0];
      firstDest = Footways[0].Nodes[0];
      startMinDist = distBetween2Points(Nodes[Footways[0].Nodes[0]].Lat,
                                        Nodes[Footways[0].Nodes[0]].Lon,
                                        info1.Coords.Lat, info1.Coords.Lon);
      endMinDist = distBetween2Points(Nodes[Footways[0].Nodes[0]].Lat,
                                      Nodes[Footways[0].Nodes[0]].Lon,
                                      info2.Coords.Lat, info2.Coords.Lon);
      for (auto const& x : Footways) {
        for (auto const& y : x.Nodes) {
          if (distBetween2Points(Nodes[y].Lat, Nodes[y].Lon, info1.Coords.Lat,
                                 info1.Coords.Lon) < startMinDist) {
            firstStart = y;
            startMinDist = distBetween2Points(
                Nodes[y].Lat, Nodes[y].Lon, info1.Coords.Lat, info1.Coords.Lon);
          }

          if (distBetween2Points(Nodes[y].Lat, Nodes[y].Lon, info2.Coords.Lat,
                                 info2.Coords.Lon) < endMinDist) {
            firstDest = y;
            endMinDist = distBetween2Points(Nodes[y].Lat, Nodes[y].Lon,
                                            info2.Coords.Lat, info2.Coords.Lon);
          }
        }
      }
      cout << endl;
      cout << "Nearest P1 node:" << endl;
      cout << " " << firstStart << endl;
      cout << " (" << Nodes[firstStart].Lat << ","
           << " " << Nodes[firstStart].Lon << ")" << endl;

      cout << "Nearest P2 node:" << endl;
      cout << " " << firstDest << endl;
      cout << " (" << Nodes[firstDest].Lat << ","
           << " " << Nodes[firstDest].Lon << ")" << endl;

      // find the nearest node to the destination building
      // this becomes the "Dest" node
      // long long destinationNode = nearestNode(BuildingInfo
      // destinationBuilding);
      double minimum = INF;
      long long dest = Footways[0].Nodes[0];
      // BuildingInfo destinationNode;
      for (auto const& x : Footways) {
        for (auto const& y : x.Nodes) {
          distance = distBetween2Points(Nodes[y].Lat, Nodes[y].Lon,
                                        destinationBuilding.Coords.Lat,
                                        destinationBuilding.Coords.Lon);
          if (distance < minimum) {
            dest = y;
            minimum = distance;
          }
        }
      }
      map<long long, double> distances1;
      map<long long, long long> predecessors;
      cout << "Nearest destination node:" << endl;
      cout << " " << dest << endl;
      cout << " (" << Nodes[dest].Lat << ","
           << " " << Nodes[dest].Lon << ")" << endl;
      cout << endl;
      // TODO MILESTONE 10: Run Dijkstra’s Algorithm
      vector<long long> answer =
          Dijkstra(G, firstStart, distances1, predecessors);
      if (distances1[firstDest] >= INF) {
        cout << "Sorry, destination unreachable." << endl;

        cout << "Enter person 1's building (partial name or abbreviation), or "
                "#> ";
        getline(cin, person1Building);
        continue;
      }

      map<long long, double> distances2;

      if (distances1[dest] >= INF || distances2[dest] >= INF) {
        cout << "At least one person was unable to reach the destination "
                "building.";
        cout << " Finding next closest building..." << endl;
        // find the 2nd closest building to the center
        set<string>
            unreachableBuildings;  // storing the unreachableBuildings names
        unreachableBuildings.insert(destinationBuilding.Fullname);
        BuildingInfo newDestinationBuilding;
        double newMin = INF;
        double newDistance;
        for (auto const& building : Buildings) {
          // cout << building.Fullname << endl;

          if (unreachableBuildings.count(building.Fullname) != 0) {
            // FOund in set
            continue;
          }
          newDistance =
              distBetween2Points(midpoint.Lat, midpoint.Lon,
                                 building.Coords.Lat, building.Coords.Lon);

          if (newDistance < newMin) {
            newDestinationBuilding = building;
            newMin = newDistance;
          }
        }
        cout << "New destination building:" << endl;
        cout << " " << newDestinationBuilding.Fullname << endl;
        cout << " (" << newDestinationBuilding.Coords.Lat << ", "
             << newDestinationBuilding.Coords.Lon << ")" << endl;
        // find the nearest destination node:
        newMin = INF;
        long long newDestNode = Footways[0].Nodes[0];
        for (auto const& a : Footways) {
          for (auto const& b : a.Nodes) {
            newDistance = distBetween2Points(Nodes[b].Lat, Nodes[b].Lon,
                                             newDestinationBuilding.Coords.Lat,
                                             newDestinationBuilding.Coords.Lon);
            if (newDistance < newMin) {
              newDestNode = b;
              newMin = newDistance;
            }
          }
        }
        map<long long, double> newDistances1;
        map<long long, double> newDistances2;
        map<long long, long long> newPredecessors;
        cout << "Nearest destination node:" << endl;
        cout << " " << newDestNode << endl;
        cout << " (" << Nodes[newDestNode].Lat << ","
             << " " << Nodes[newDestNode].Lon << ")" << endl;
        cout << endl;
        vector<long long> newAnswer =
            Dijkstra(G, firstStart, newDistances1, newPredecessors);

        cout << "Person 1's distance to dest: " << newDistances1[newDestNode]
             << " miles" << endl;

        cout << "Path: ";
        getPath(newPredecessors, firstStart, newDestNode);

        cout << endl;
        newAnswer = Dijkstra(G, firstDest, newDistances2, newPredecessors);
        cout << "Person 2's distance to dest: " << newDistances2[newDestNode]
             << " miles" << endl;
        cout << "Path: ";
        getPath(newPredecessors, firstDest, newDestNode);
        cout << endl;

        cout << "Enter person 1's building (partial name or abbreviation), or "
                "#> ";
        getline(cin, person1Building);
        continue;
      }

      cout << "Person 1's distance to dest: " << distances1[dest] << " miles"
           << endl;
      // get path
      cout << "Path: ";
      getPath(predecessors, firstStart, dest);
      cout << endl;
      answer = Dijkstra(G, firstDest, distances2, predecessors);

      cout << "Person 2's distance to dest: " << distances2[dest] << " miles"
           << endl;
      // TODO MILESTONE 11: Print path using predecessors (if MS10 succeeded)
      cout << "Path: ";
      getPath(predecessors, firstDest, dest);

      // TODO MILESTONE 11: Find Second Nearest Destination (if MS10 failed)
      //
      // another navigation?
      //

      //
      // TO DO: lookup buildings, find nearest start and dest nodes, find center
      // run Dijkstra's alg from each start, output distances and paths to
      // destination:
      //

      // cout << "Person 1's building not found" << endl;
      // cout << "Person 2's building not found" << endl;

      //
      // another navigation?
      //
      cout << endl;
      cout
          << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
    }
  }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates> Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo> Footways;
  // info about each building, in no particular order
  vector<BuildingInfo> Buildings;
  XMLDocument xmldoc;

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

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  graph<long long, double> G;
  // TODO MILESTONE 5: add vertices
  // add each node to the graph
  for (auto x : Nodes) {
    G.addVertex(x.first);
  }
  // TODO MILESTONE 6: add edges
  // add edges based on footways
  for (auto const c : Footways) {
    for (int i = 0; i < (int)c.Nodes.size() - 1; i++) {
      G.addEdge(
          c.Nodes[i], c.Nodes[i + 1],
          distBetween2Points(Nodes.at(c.Nodes[i]).Lat, Nodes.at(c.Nodes[i]).Lon,
                             Nodes.at(c.Nodes[i + 1]).Lat,
                             Nodes.at(c.Nodes[i + 1]).Lon));
      G.addEdge(
          c.Nodes[i + 1], c.Nodes[i],
          distBetween2Points(Nodes.at(c.Nodes[i]).Lat, Nodes.at(c.Nodes[i]).Lon,
                             Nodes.at(c.Nodes[i + 1]).Lat,
                             Nodes.at(c.Nodes[i + 1]).Lon));
    }
  }

  //
  // TO DO: build the graph, output stats:
  //

  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
       << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    // TO DO: add argument for the graph you make.
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    // TO DO: add arguments
    creative();
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
