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

#pragma once

#include <vector>
#include "math.h"

enum class SimilarityMetric
{
  ManhattanDistance,
  EuclideanDistance,
  ChebyshevDistance,
  SurfaceArea,
};

enum class Linkage
{
  Single,
  Complete,
  Centroid,
  Average,
};

struct AggCluster
{
  std::vector<size_t> pointIDs;
  math::box2f bounds;
  int level;

  // [0,1] centroid in X to draw dendograms from
  // (don't confuse with cluster centroid!)
  float centroid;

  // for convenience, we also store the cluster IDs
  // from the previous level that we combined
  int id1, id2;
};

typedef std::vector<math::vec2f> Points;
typedef std::vector<AggCluster> Clusters;
typedef std::vector<math::box1f> Ranges;

class Dendrogram
{
public:
  Dendrogram() = default;
  Dendrogram(const std::vector<math::vec2f> &points);
  void reset(const std::vector<math::vec2f> &points);
  bool step();

  Points getPoints() const;
  Clusters getClusters(int level = -1) const;
  Ranges getRanges(int level = -1) const;

  int currentLevel = 0;

  SimilarityMetric metric = SimilarityMetric::EuclideanDistance;
  Linkage linkage = Linkage::Centroid;
private:
  std::vector<math::vec2f> input;
  Clusters clusters;
  Ranges hranges;
  std::vector<Clusters> levels;
};
// vim: sw=2:expandtab:softtabstop=2:ts=2:cino=\:0g0t0

