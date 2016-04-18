#-------------------------------------------------
#
# Project created by QtCreator 2013-09-12T16:39:14
#
#-------------------------------------------------

QT       += core

QT       -= gui  
TARGET = Random
CONFIG   += console
CONFIG   -= app_bundle

DESTDIR = ../bin

TEMPLATE = app
QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS+=-march=corei7 -msse4.2

SOURCES += main.cpp \
    ../Lib/PhysMtr.cpp \
    ../Lib/matrix.cpp \
    ../Lib/Geometry.cpp \
    ../Lib/Crystal.cpp \
    ../Lib/Mueller.cpp \
    ../Lib/particle.cpp \
    ../Lib/trajectory.cpp \
    ../Lib/compl.cpp \
    ../Lib/Scattering.cpp \
    ../Lib/beam.cpp

HEADERS += \
    ../Lib/Beam.hpp \
    ../Lib/compl.hpp \
    ../Lib/Crystal.hpp \
    ../Lib/Geometry.hpp \
    ../Lib/matrix.hpp \
    ../Lib/Mueller.hpp \
    ../Lib/particle.hpp \
    ../Lib/PhysMtr.hpp \
    ../Lib/service.hpp \
    ../Lib/trajectory.hpp \
    ../Lib/Intersection.h
