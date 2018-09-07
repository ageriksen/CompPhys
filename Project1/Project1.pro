TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp


LIBS += -LC:\armadillo\newblas -llibblas
LIBS += -LC:\armadillo\newblas -lliblapack

INCLUDEPATH += C:\armadillo\include
DEPENDPATH += C:\armadillo\include

HEADERS += \
    project1h.h
