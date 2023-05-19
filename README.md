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
