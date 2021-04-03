#include "mainwindow.h"
#include <QApplication>
#include <QLibrary>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <dlfcn.h>
#include <gmp.h>
#include <mpfr.h>
#include "mpreal.h"


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();



    return a.exec();
}
