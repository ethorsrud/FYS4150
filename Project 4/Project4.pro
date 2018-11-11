TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/Cellar/armadillo/9.100.5_1/include
LIBS += -L/usr/local/Cellar/armadillo/9.100.5_1/lib -larmadillo -llapack -lblas
INCLUDEPATH += /usr/local/Cellar/open-mpi/3.1.3/include
LIBS += -L/usr/local/Cellar/open-mpi/3.1.3/lib

SOURCES += \
        main.cpp \
    ising.cpp

HEADERS += \
    ising.h


# MPI Settings
QMAKE_CXX = /usr/local/bin/mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = /usr/local/bin/mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
