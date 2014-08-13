#-------------------------------------------------
#
# Project created by QtCreator 2013-09-12T16:39:14
#
#-------------------------------------------------

QT       += core

QT       -= gui
  bnvbnvbnvbnvbnv
TARGET = Random
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app
QMAKE_CXXFLAGS += -std=gnu++11

QMAKE_CFLAGS += -march=corei7-avx -mtune=corei7-avx
QMAKE_CXXFLAGS += -march=corei7-avx -mtune=corei7-avx -mavx

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

HEADERS +=
