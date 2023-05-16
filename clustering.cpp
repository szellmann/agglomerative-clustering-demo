// ======================================================================== //
// Copyright 2023-2023 Stefan Zellmann                                      //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include <algorithm>
#include <iostream>
#include <cassert>
#include <cfloat>
#include "clustering.h"
using namespace math;

Dendrogram::Dendrogram(const std::vector<vec2f> &points)
{
  reset(points);
}

void Dendrogram::reset(const std::vector<vec2f> &points)
{
  float step = 1.f/points.size();

  input = points;
  clusters.clear();
  levels.clear();
  currentLevel = 0;
  // at the bottom level, each point becomes its own cluster
  for (size_t i=0; i<points.size(); ++i) {
    AggCluster cluster;
    cluster.pointIDs.push_back(i);
    cluster.bounds = {points[i],points[i]};
    cluster.level = currentLevel;
    cluster.centroid = step/2+step*i;
    cluster.id1 = -1;
    cluster.id2 = -1;
    clusters.push_back(cluster);
  }
}

bool Dendrogram::step()
{
  if (clusters.size() == 1)
    return false;

  else if (clusters.size() > 1) { // only if we're not done yet

    // take a "snapshot" of the current clusters
    // we do this primarily so we can later draw the clusters
    levels.push_back(clusters);

    // find the two most similar clusters
    // this is the costly part of this algorithm!
    // in the next lecture we will learn how to optimize
    // this for better/real-time performance
    size_t id1, id2;
    float bestSimilarity = FLT_MAX;
    for (size_t i=0; i<clusters.size(); ++i) {
      for (size_t j=i+1; j<clusters.size(); ++j) {
        float dist = length(clusters[i].bounds.center()-clusters[j].bounds.center());
        if (dist < bestSimilarity) {
          id1 = i;
          id2 = j;
          bestSimilarity = dist;
        }
      }
    }

    assert(bestSimilarity < FLT_MAX);

    // std::cout << "combining clusters " << id1 << ',' << id2
    //           << ", bounds: " << clusters[id1].bounds << ','
    //           << clusters[id2].bounds << '\n';

    // combine the two clusters we've found
    AggCluster newCluster;
    box2f bounds{FLT_MAX,-FLT_MAX};
    for (size_t i=0; i<clusters[id1].pointIDs.size(); ++i) {
      newCluster.pointIDs.push_back(clusters[id1].pointIDs[i]);
      bounds.extend(input[clusters[id1].pointIDs[i]]);
    }

    for (size_t i=0; i<clusters[id2].pointIDs.size(); ++i) {
      newCluster.pointIDs.push_back(clusters[id2].pointIDs[i]);
      bounds.extend(input[clusters[id2].pointIDs[i]]);
    }

    newCluster.bounds = bounds;
    newCluster.level = currentLevel+1;
    newCluster.centroid = (clusters[id1].centroid+clusters[id2].centroid)/2;
    newCluster.id1 = id1;
    newCluster.id2 = id2;

    // remove the two clusters
    clusters[id1].level = -1; // mark for deletion
    clusters[id2].level = -1; //  "
    clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                  [](const AggCluster &c) { return c.level == -1; }),
                   clusters.end());

    // add the new cluster
    clusters.push_back(newCluster);

    currentLevel++;
  }

  return true;
}

Points Dendrogram::getPoints() const
{
  return input;
}

Clusters Dendrogram::getClusters(int level) const
{
  if (level == -1)
    return clusters;
  else
    return levels[level];
}
// vim: sw=2:expandtab:softtabstop=2:ts=2:cino=\:0g0t0

