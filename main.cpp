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

#include <cfloat>
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <QApplication>
#include <QMainWindow>
#include <QLayout>
#include <QScrollArea>
#include <QPainter>
#include <QPainterPath>
#include <QLabel>
#include <QPushButton>
#include <QCheckBox>
#include <QComboBox>
#include <QMouseEvent>
#include "clustering.h"
#include "math.h"

using namespace math;

int g_initialDist = 2;
bool g_showClusters = true;
bool g_colorizeByCluster = false;
float g_cutLevel = 0.f;
float g_maxCutLevel = 0.f;

// [0,1] to cool/warm color map
static QColor mapColor(float val01)
{
  vec3f cols[] = {{0.231f, 0.298f, 0.752f},
                  {0.552f, 0.690f, 0.996f},
                  {0.866f, 0.866f, 0.866f},
                  {0.956f, 0.603f, 0.486f},
                  {0.705f, 0.015f, 0.149f}};
  int numCols = sizeof(cols)/sizeof(cols[0]);
  int col1 = std::max(0,std::min((int)ceilf(val01*numCols),numCols-1));
  int col2 = std::max(0,std::min((int)floorf(val01*numCols),numCols-1));
  float frac = val01*numCols-floorf(val01*numCols);
  vec3f col = lerp(cols[col1],cols[col2],frac);
  return QColor::fromRgbF(col.x,col.y,col.z);
}

inline QColor randomColor(size_t idx)
{
  unsigned int r = (unsigned int)(idx*13*17 + 0x234235);
  unsigned int g = (unsigned int)(idx*7*3*5 + 0x773477);
  unsigned int b = (unsigned int)(idx*11*19 + 0x223766);
  return QColor::fromRgbF((r&255)/255.f,
                          (g&255)/255.f,
                          (b&255)/255.f);
}

std::vector<vec2f> randomPoints(size_t n)
{
  std::uniform_real_distribution<float> dist(0,1);
  std::default_random_engine rng;

  std::vector<vec2f> data(n);

  for (size_t i=0; i<n; ++i) {
    data[i] = {dist(rng),dist(rng)};
  }

  return data;
}

// lifted/adapted from: https://www.naftaliharris.com/blog/visualizing-k-means-clustering/
std::vector<vec2f> smiley(size_t n)
{
  std::uniform_real_distribution<float> dist(0,1);
  std::default_random_engine rng;

  std::vector<vec2f> data(n);

  // scale from [-10,+10] to [0,1] (and flip!)
  auto scale = [](vec2f p) { return vec2f{(p.x+10.f)/20.f,1.f-(p.y+10.f)/20.f}; };

  for (size_t i=0; i<n; ) {
    float x = dist(rng) * 20 - 10;
    float y = dist(rng) * 20 - 10;

    // Smiley params
    float a = 6;
    float b = 8;
    float d = 4;
    float c = 0.15;
    float C = (c - a / (b*b));

    // The border
    if(x*x + y*y < 100 && x*x + y*y > 81) {
      data[i] = scale({x,y});
      i += 1;
    }
    // Left eye
    else if((x+4)*(x+4) + (y-4)*(y-4) < 1) {
      data[i] = scale({x,y});
      i += 1;
    }
    // Right eye
    else if((x-4)*(x-4) + (y-4)*(y-4) < 1) {
      data[i] = scale({x,y});
      i += 1;
    }
    // Smile
    else if((x > -b && x < b) && (y > c*x*x - a) && (y < C*x*x - d)) {
      data[i] = scale({x,y});
      i += 1;
    }
  }

  return data;
}

std::vector<vec2f> preset1()
{
  std::vector<vec2f> data(16);

  data[ 0] = {0.13f, 0.12f};
  data[ 1] = {0.21f, 0.29f};
  data[ 2] = {0.18f, 0.18f};
  data[ 3] = {0.08f, 0.05f};
  data[ 4] = {0.03f, 0.07f};
  data[ 5] = {0.13f, 0.22f};
  data[ 6] = {0.03f, 0.13f};
  data[ 7] = {0.15f, 0.14f};
  data[ 8] = {0.11f, 0.03f};
  data[ 9] = {0.02f, 0.03f};
  data[10] = {0.25f, 0.13f};

  data[11] = {0.78f, 0.75f};
  data[12] = {0.74f, 0.71f};
  data[13] = {0.81f, 0.69f};
  data[14] = {0.61f, 0.73f};
  data[15] = {0.80f, 0.91f};

  return data;
}

class ClusterWidget : public QWidget
{
public:
  Q_OBJECT

public:
  ClusterWidget(QWidget *parent) : QWidget(parent) {}

  void setDendrogram(const Dendrogram &dendro)
  {
    this->dendro = &dendro;
  }

  void paintEvent(QPaintEvent *)
  {
    QPainter painter(this);

    QBrush brush;
    brush.setColor(Qt::gray);
    brush.setStyle(Qt::SolidPattern);

    painter.setBrush(brush);
    painter.drawRect(0, 0, width(), height());

    // scale from [0,1] to [W,H), and center at [W/2,H/2]
    auto scale = [this](vec2f p) {
      int dim = std::min(width(),height());
      int offX = dim==width()?0:(width()-dim)/2;
      int offY = dim==height()?0:(height()-dim)/2;
      return vec2f{
        offX+p.x*dim,
        offY+p.y*dim
      };
    };

    constexpr float R = 5.f;

    const Points &points = dendro->getPoints();
    for (size_t i=0; i<points.size(); ++i) {
      QPen pen;
      pen.setColor(Qt::black);
      painter.setPen(pen);

      QBrush cbrush;
      if (g_colorizeByCluster) {
        const Clusters &clusters = dendro->getClusters(g_cutLevel);
        for (size_t clusterID=0;clusterID<clusters.size();++clusterID) {
          bool found=false;
          for (size_t j=0; j<clusters[clusterID].pointIDs.size(); ++j) {
            if (clusters[clusterID].pointIDs[j] == i) {
              cbrush.setColor(randomColor(clusterID));
              found = true;
              break;
            }
          }
          if (found) break;
        }
      }
      else
        cbrush.setColor(mapColor(i/float(points.size()-1)));
      cbrush.setStyle(Qt::SolidPattern);
      painter.setBrush(cbrush);

      vec2f center = scale(points[i]);
      painter.drawEllipse(QPointF(center.x,center.y),R,R);
    }

    if (g_showClusters) {
      const Clusters &clusters = dendro->getClusters(g_cutLevel);
      for (size_t i=0; i<clusters.size(); ++i) {
        QPen pen;
        pen.setColor(Qt::black);
        if (clusters.size() < points.size() && i==clusters.size()-1) {
          // that's the cluster we just added
          pen.setWidth(5);
        }
        painter.setPen(pen);

        QBrush cbrush;
        cbrush.setStyle(Qt::NoBrush);
        painter.setBrush(cbrush);

        box2f bounds = clusters[i].bounds;
        vec2f pt1 = scale(bounds.lower)-R-1;
        vec2f pt2 = scale(bounds.upper)+R+2;
        painter.drawRect(QRect(pt1.x,pt1.y,pt2.x-pt1.x,pt2.y-pt1.y));
      }
    }

    QFont f = font();
    f.setPointSize(28);
    painter.setFont(f);

    QPen pen;
    pen.setColor(Qt::white);
    painter.setPen(pen);
#ifdef __APPLE__
    painter.drawText(QPointF(0.f,22.f),"2D points");
#else
    painter.drawText(QPointF(0.f,28.f),"2D points");
#endif
  }

private:
  const Dendrogram *dendro = nullptr;
};

class DendroWidget : public QWidget
{
public:
  Q_OBJECT

signals:
  void cutUpdated();

public:
  DendroWidget(QWidget *parent) : QWidget(parent)
  {
    setMouseTracking(true);
  }

  void setDendrogram(const Dendrogram &dendro)
  {
    this->dendro = &dendro;
  }

  float computeXCoord(size_t pointID)
  {
    const Points &inputPoints = dendro->getPoints();
    float step = 1.f/inputPoints.size();
    return step/2+step*pointID;
  }

  float computeYCoord(float level)
  {
    constexpr float paddingBottom = 0.05f;
    constexpr float paddingTop = 0.1f;
    float result = 1.f-paddingBottom;
    result -= level*((1.f-paddingTop-paddingBottom)/dendro->getPoints().size());
    return result;
  }

  float computeLevel(float y)
  {
    constexpr float paddingBottom = 0.05f;
    constexpr float paddingTop = 0.1f;
    y = 1.f - y;
    y -= paddingBottom;
    y /= ((1.f-paddingTop-paddingBottom)/dendro->getPoints().size());
    return y;
  }

  // leave some room for dendrogram-cut slider
  int canvasWidth()
  {
    return width()-18;
  }

  bool sliderDragging=false;
  void mousePressEvent(QMouseEvent *event)
  {
    QPoint pos = event->pos();
    if (pos.x() <= canvasWidth())
      return;

    if (!sliderDragging) {
      float cutY = computeYCoord(g_cutLevel)*height();
      sliderDragging = true;
      float y = pos.y();
      g_cutLevel = std::max(0.f,std::min(computeLevel(y/height()),g_maxCutLevel));
      emit cutUpdated();
    }
  }

  void mouseReleaseEvent(QMouseEvent *event)
  {
    QPoint pos = event->pos();
    if (pos.x() <= canvasWidth()) {
      sliderDragging = false;
      emit cutUpdated();
      return;
    }

    float y = pos.y();
    g_cutLevel = std::max(0.f,std::min(computeLevel(y/height()),g_maxCutLevel));
    sliderDragging = false;
    emit cutUpdated();
  }

  void mouseMoveEvent(QMouseEvent *event)
  {
    if (sliderDragging) {
      float y = event->pos().y();
      g_cutLevel = std::max(0.f,std::min(computeLevel(y/height()),g_maxCutLevel));
      emit cutUpdated();
    }
  }

  void paintEvent(QPaintEvent *)
  {
    QPainter painter(this);

    QBrush brush;
    brush.setColor(Qt::white);
    brush.setStyle(Qt::SolidPattern);

    painter.setBrush(brush);
    painter.drawRect(0, 0, canvasWidth(), height());

    brush.setColor(Qt::gray);

    painter.setBrush(brush);
    painter.drawRect(canvasWidth(), 0, width()-canvasWidth(), height());

    // scale from [0,1] to [W,H]
    auto scale = [this](vec2f p) {
      return vec2f{
        p.x*canvasWidth(),
        p.y*height()
      };
    };

    const Points &inputPoints = dendro->getPoints();
    Points points(dendro->getPoints().size());
    for (size_t i=0; i<points.size(); ++i) {
      points[i] = {computeXCoord(i),computeYCoord(0)};
    }

    for (int l=1; l<=dendro->currentLevel; ++l) {
      const Clusters &prev = dendro->getClusters(l-1);
      const Clusters &curr = dendro->getClusters(l<dendro->currentLevel? l: -1);

      int id1 = curr.back().id1;
      int id2 = curr.back().id2;

      float y1 = computeYCoord(prev[id1].level)*height();
      float y2 = computeYCoord(prev[id2].level)*height();
      float y = computeYCoord(l)*height();

      float x1 = prev[id1].centroid*canvasWidth();
      float x2 = prev[id2].centroid*canvasWidth();

      //constexpr float minDist = 7.f;
      //if (fabsf(x1-x2) < minDist) {
      //  if (x1>x2) {
      //    x2-=minDist/2;
      //    x1+=minDist/2;
      //  } else {
      //    x1-=minDist/2;
      //    x2+=minDist/2;
      //  }
      //}

      QPen pen;
      pen.setColor(Qt::black);
      if (g_showClusters && l==int(g_cutLevel)) {
        pen.setWidth(3);
      }
      painter.setPen(pen);
      painter.drawLine(QPoint(x1,y1), QPoint(x1,y));
      painter.drawLine(QPoint(x2,y2), QPoint(x2,y));
      painter.drawLine(QPoint(x1,y), QPoint(x2,y));
    }

    for (size_t i=0; i<points.size(); ++i) {
      QPen pen;
      pen.setColor(Qt::black);
      painter.setPen(pen);

      QBrush cbrush;
      if (g_colorizeByCluster) {
        const Clusters &clusters = dendro->getClusters(g_cutLevel);
        for (size_t clusterID=0;clusterID<clusters.size();++clusterID) {
          bool found=false;
          for (size_t j=0; j<clusters[clusterID].pointIDs.size(); ++j) {
            if (clusters[clusterID].pointIDs[j] == i) {
              cbrush.setColor(randomColor(clusterID));
              found = true;
              break;
            }
          }
          if (found) break;
        }
      }
      else
        cbrush.setColor(mapColor(i/float(points.size()-1)));
      cbrush.setStyle(Qt::SolidPattern);
      painter.setBrush(cbrush);

      vec2f center = scale(points[i]);
      constexpr float R = 5.f;
      painter.drawEllipse(QPointF(center.x,center.y),R,R);
    }

    // cut-slider
    float cutY = computeYCoord(g_cutLevel)*height();
    QPainterPath path;
    path.moveTo(canvasWidth()+1,cutY);
    path.lineTo(canvasWidth()+1+6,cutY+6);
    path.lineTo(canvasWidth()+1+14,cutY+6);
    path.lineTo(canvasWidth()+1+14,cutY-6);
    path.lineTo(canvasWidth()+1+6,cutY-6);
    path.lineTo(canvasWidth()+1,cutY);

    QPen pen;
    pen.setColor(Qt::black);
    painter.setPen(pen);

    brush.setColor(Qt::white);
    painter.setBrush(brush);

    painter.drawPath(path);

    if (1) { // dashed line
      pen.setStyle(Qt::DashLine);
      pen.setWidth(sliderDragging? 3: 1);
      pen.setColor(Qt::gray);

      painter.setPen(pen);
      painter.drawLine(QPoint(0.f,cutY), QPoint(canvasWidth(),cutY));

      pen.setStyle(Qt::SolidLine);
      pen.setWidth(1);
      pen.setColor(Qt::black);
    }

    QFont f = font();
    f.setPointSize(28);
    painter.setFont(f);

    pen.setColor(Qt::black);
    painter.setPen(pen);
#ifdef __APPLE__
    painter.drawText(QPointF(0.f,22.f),"Dendrogram");
#else
    painter.drawText(QPointF(0.f,28.f),"Dendrogram");
#endif
  }

private:
  const Dendrogram *dendro = nullptr;
};

int main(int argc, char **argv)
{
  QApplication app(argc,argv);

  QMainWindow win;
  win.setWindowTitle("Agglomerative Clustering Interactive Demo");
  win.resize(1280,800);

  // sort the points so the dendrogram has fewer crossings
  // (merely for aesthetic reasons)
  auto doSort = [](const std::vector<vec2f> &pts) {
    std::vector<vec2f> out = pts;
    std::stable_sort(out.begin(),out.end(),[](vec2f a, vec2f b) { return a.x<b.x; });
    std::stable_sort(out.begin(),out.end(),[](vec2f a, vec2f b) { return a.y<b.y; });
    return out;
  };

  Dendrogram dendro;

  std::vector<vec2f> pts;

  auto reset = [&]() {
    if (g_initialDist == 0)
      pts = doSort(randomPoints(100));
    else if (g_initialDist == 1)
      pts = doSort(smiley(150));
    else if (g_initialDist == 2)
      pts = doSort(preset1());
    dendro.reset(pts);
    g_cutLevel = 0.f;
    g_maxCutLevel = 0.f;
  };

  reset();

  // cluster and dendroid widgets
  ClusterWidget *clusterWidget = new ClusterWidget(&win);
  clusterWidget->setDendrogram(dendro);

  DendroWidget *dendroWidget = new DendroWidget(&win);
  dendroWidget->setDendrogram(dendro);

  QObject::connect(dendroWidget, &DendroWidget::cutUpdated,
    [&]() {
      clusterWidget->repaint();
      dendroWidget->repaint();
    });

  QHBoxLayout *hlayout = new QHBoxLayout;
  hlayout->addWidget(clusterWidget);
  hlayout->addWidget(dendroWidget);

  QVBoxLayout *vlayout = new QVBoxLayout;
  vlayout->addLayout(hlayout);

  // buttons
  QHBoxLayout *buttonLayout = new QHBoxLayout;
  vlayout->addLayout(buttonLayout);

  QPushButton *resetRandomButton = new QPushButton(&win);
  resetRandomButton->setText("Reset (random/uniform)");
  resetRandomButton->setMaximumWidth(160);
  QObject::connect(resetRandomButton, &QPushButton::pressed,
    [&]() {
      g_initialDist = 0;
      reset();
      clusterWidget->repaint();
      dendroWidget->repaint();
    });
  buttonLayout->addWidget(resetRandomButton);

  QPushButton *resetSmileyButton = new QPushButton(&win);
  resetSmileyButton->setText("Reset (smiley)");
  resetSmileyButton->setMaximumWidth(160);
  QObject::connect(resetSmileyButton, &QPushButton::pressed,
    [&]() {
      g_initialDist = 1;
      reset();
      clusterWidget->repaint();
      dendroWidget->repaint();
    });
  buttonLayout->addWidget(resetSmileyButton);

  QPushButton *resetPreset1Button = new QPushButton(&win);
  resetPreset1Button->setText("Reset (preset 1)");
  resetPreset1Button->setMaximumWidth(160);
  QObject::connect(resetPreset1Button, &QPushButton::pressed,
    [&]() {
      g_initialDist = 2;
      reset();
      clusterWidget->repaint();
      dendroWidget->repaint();
    });
  buttonLayout->addWidget(resetPreset1Button);

  QHBoxLayout *metricLayout = new QHBoxLayout;
  QLabel *metricLabel = new QLabel(&win);
  metricLabel->setText("Similarity Metric:");
  metricLabel->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  QComboBox *metricBox = new QComboBox(&win);
  metricBox->addItem("Manhattan Distance");
  metricBox->addItem("Euclidean Distance");
  metricBox->addItem("Chebyshev Distance");
  metricBox->addItem("Surface Area (2D)");
  if (dendro.metric == SimilarityMetric::ManhattanDistance)
    metricBox->setCurrentIndex(0);
  else if (dendro.metric == SimilarityMetric::EuclideanDistance)
    metricBox->setCurrentIndex(1);
  else if (dendro.metric == SimilarityMetric::ChebyshevDistance)
    metricBox->setCurrentIndex(2);
  else if (dendro.metric == SimilarityMetric::SurfaceArea)
    metricBox->setCurrentIndex(3);
  QObject::connect(metricBox, qOverload<int>(&QComboBox::currentIndexChanged),
    [&](int index) {
      reset();
      if (index == 0)
        dendro.metric = SimilarityMetric::ManhattanDistance;
      else if (index == 1)
        dendro.metric = SimilarityMetric::EuclideanDistance;
      else if (index == 2)
        dendro.metric = SimilarityMetric::ChebyshevDistance;
      else if (index == 3)
        dendro.metric = SimilarityMetric::SurfaceArea;
      clusterWidget->repaint();
      dendroWidget->repaint();
    });
  metricLayout->addWidget(metricLabel);
  metricLayout->addWidget(metricBox);
  buttonLayout->addLayout(metricLayout);

  QCheckBox *showClustersBox = new QCheckBox(&win);
  showClustersBox->setText("Show Clusters");
  showClustersBox->setMaximumWidth(160);
  showClustersBox->setCheckState(g_showClusters? Qt::Checked: Qt::Unchecked);
  QObject::connect(showClustersBox, &QCheckBox::clicked,
    [&](bool value) {
      g_showClusters = value;
      clusterWidget->repaint();
      dendroWidget->repaint();
    });
  buttonLayout->addWidget(showClustersBox);

  QCheckBox *colorizeBox = new QCheckBox(&win);
  colorizeBox->setText("Color by Cluster ID");
  colorizeBox->setMaximumWidth(160);
  colorizeBox->setCheckState(g_colorizeByCluster? Qt::Checked: Qt::Unchecked);
  QObject::connect(colorizeBox, &QCheckBox::clicked,
    [&](bool value) {
      g_colorizeByCluster = value;
      clusterWidget->repaint();
      dendroWidget->repaint();
    });
  buttonLayout->addWidget(colorizeBox);

  QPushButton *computeButton = new QPushButton(&win);
  computeButton->setText("Compute Dendrogram");
  computeButton->setMaximumWidth(160);
  QObject::connect(computeButton, &QPushButton::pressed,
    [&]() {
      for (;;) {
        bool done = !dendro.step();
        if (done) break;
        g_cutLevel = dendro.currentLevel;
        g_maxCutLevel = dendro.currentLevel;
      }
      clusterWidget->repaint();
      dendroWidget->repaint();
    });
  buttonLayout->addWidget(computeButton);

  QPushButton *stepButton = new QPushButton(&win);
  stepButton->setText("Step");
  stepButton->setMaximumWidth(160);
  QObject::connect(stepButton, &QPushButton::pressed,
    [&]() {
      dendro.step();
      g_cutLevel = dendro.currentLevel;
      g_maxCutLevel = dendro.currentLevel;
      clusterWidget->repaint();
      dendroWidget->repaint();
    });
  buttonLayout->addWidget(stepButton);

  // central widget/scroll area
  QWidget *centralWidget = new QWidget(&win);
  centralWidget->setLayout(vlayout);

  QScrollArea *scrollArea = new QScrollArea(&win);
  scrollArea->setWidget(centralWidget);
  scrollArea->setWidgetResizable(true);
  win.setCentralWidget(scrollArea);

  win.show();
  app.exec();
}

#include "main.moc"
// vim: sw=2:expandtab:softtabstop=2:ts=2:cino=\:0g0t0

