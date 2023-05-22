Agglomerative Clustering Demo App
=================================

Demo App for "Real-time computer graphics, Lecture 8"

Dependencies
------------

- cmake
- C++-Compiler
- QT5 or QT6

Build instructions
------------------

Tested on macOS and Ubuntu Linux, for example on macOS/terminal (QT 6.5.0 via homebrew):
```
brew install qt6
mkdir build && cd build
CMAKE_PREFIX_PATH=/opt/homebrew/Cellar/qt/6.5.0/lib/cmake/Qt6Widgets cmake <source-dir>
make
open agglomerativeClustering.app/
```

Code organization
-----------------

The clustering implementation can be found in [clustering.cpp](/clustering.cpp).
Definitions for the dendrogram are found in [clustering.h](/clustering.h).
The `Dendrogram` class stores the sets of clusters generated during (naive)
construction, allowing us to easily compute cuts, but prohibits scaling to
serious sizes beyond the toy examples. Besides, the data gets initialized
in [main.cpp](/main.cpp) in the `main` function.
