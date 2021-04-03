#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QDebug>
#include <QLibrary>
#include <gmp.h>
#include <mpfr.h>
#include "mpreal.h"

QString libname;

typedef mpfr::mpreal(*MyPrototype)(mpfr::mpreal);
typedef std::string(*MyPrototype2)();



MyPrototype func;
MyPrototype dfunc;
MyPrototype d2func;


QString NewtonRaphson(mpfr::mpreal x, mpfr::mpreal eps, int it);

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QWidget::setWindowTitle("Wyznaczanie pierwiastków metodą Newtona-Raphsona");

    ui->FileDisplay->setText(QString(". . ."));
}

MainWindow::~MainWindow()
{
    delete ui;

}

void MainWindow::on_pushButton_clicked()
{
    QString file_name = QFileDialog::getOpenFileName(this,"Choose a file", QDir::homePath(),tr("Shared Library Files(*.so)"));
    QMessageBox::information(this,"Wybrano plik",file_name);
    std::string base = file_name.toUtf8().constData();
    base = base.substr(base.find_last_of("/") + 1);
    QString Qbase = QString::fromStdString(base);
    ui->FileDisplay->setText(Qbase);
    libname = file_name;
}

void MainWindow::on_pushButton_2_clicked()
{

    double x = QString(ui->initialx->text()).toDouble();
    int mit = QString(ui->mitlabel->text()).toInt();
    mpfr::mpreal x2 = x;

    qDebug() << "x: " << QString::fromStdString(x2.toString());

    QLibrary myLib(libname);


    MyPrototype function = (MyPrototype) myLib.resolve("func");
    MyPrototype dfunction = (MyPrototype) myLib.resolve("dfunc");
    MyPrototype d2function = (MyPrototype) myLib.resolve("d2func");
    MyPrototype2 funcname = (MyPrototype2) myLib.resolve("funcname");



    if(function && dfunction && d2function && funcname){
        std::string funcName = funcname();
        qDebug() <<  QString::fromStdString(funcName);
        qDebug() << NewtonRaphson(x2,0.0000000000000001,mit);
        ui->results->setText(NewtonRaphson(x2,0.0000000000000001,mit));

    }
    else{
        QMessageBox::information(this,"Niekompatybilna biblioteka","Biblioteka nie jest zgodna ze wzorem z dokumentacji!");
    }
}
QString NewtonRaphson(mpfr::mpreal x, mpfr::mpreal eps, int mit = 1){
    int st = 0;
    mpfr::mpreal x_1;

    //Formatowanie przedstawienia danych
    QString result;
    result.push_back("Dane:\n x = ");
    result.push_back(QString::fromStdString(x.toString(3)));
    result.push_back(", mit = ");
    result.push_back(QString::number(mit));
    result.push_back(", eps = 1e-16\nRównanie:\n");

    QString initialx = QString::fromStdString(x.toString());
    QLibrary myLib(libname);
    MyPrototype func = (MyPrototype) myLib.resolve("func");
    MyPrototype dfunc = (MyPrototype) myLib.resolve("dfunc");
    MyPrototype d2func = (MyPrototype) myLib.resolve("d2func");
    MyPrototype2 funcname = (MyPrototype2) myLib.resolve("funcname");

    std::string funcName = funcname();
    qDebug() <<  QString::fromStdString(funcName);
    result.push_back(QString::fromStdString(funcName));
    result.push_back("\nWyniki:");


    //mpfr::mpreal h = func(x) / dfunc(x);
    mpfr::mpreal hplus = (dfunc(x) + mpfr::sqrt(mpfr::pow(dfunc(x),2)-2*func(x)*d2func(x))) / d2func(x);
    mpfr::mpreal hminus = (dfunc(x) - mpfr::sqrt(mpfr::pow(dfunc(x),2)-2*func(x)*d2func(x))) / d2func(x);
    //mpfr::mpreal x1 = x;

    //Wartość wykonanych iteracji it
    int count = 0;

   do
    {
        if(mit < 1) {
            st = 1;
            break;
        }
        if (d2func(x) == 0){
            st = 2; break;
        }
        hplus = (dfunc(x) + mpfr::sqrt(mpfr::pow(dfunc(x),2)-2*func(x)*d2func(x))) / d2func(x);
        hminus = (dfunc(x) - mpfr::sqrt(mpfr::pow(dfunc(x),2)-2*func(x)*d2func(x))) / d2func(x);

        //hplus = (dfunc(x1) + mpfr::sqrt(mpfr::pow(dfunc(x1),2)-2*func(x1)*d2func(x1))) / d2func(x1);
        //hminus = (dfunc(x1) - mpfr::sqrt(mpfr::pow(dfunc(x1),2)-2*func(x1)*d2func(x1))) / d2func(x1);
        if (mpfr::pow(dfunc(x),2)-2*func(x)*d2func(x) < 0){
            st = 4; break;
        }
        x_1 = x;
        x = x - hplus;
        //x1 = x1 - hminus;

        count ++;
        qDebug() << "NewtonRaphson(x,eps,it) = "  << QString::fromStdString(x.toString());
        qDebug() << "fatx = " << QString::fromStdString(func(x).toString(30));
        qDebug() << "It : " <<QString::number(count);

        //qDebug() << "NewtonRaphson(x1,eps,it) = "  << QString::fromStdString(x1.toString());
        //qDebug() << "fatx1 = " << QString::fromStdString(func(x1).toString(30));
        //qDebug() << "It : " <<QString::number(count);
        if(func(x) == 0) {
            st = 0;
            break;
        }
        else if(mpfr::abs(x-x_1)/mpfr::max(mpfr::abs(x),mpfr::abs(x_1)) < eps){
            break;
        }
        else{
            st = 3;
            continue;
        }
    }
    while(count < mit);

    //Formatowanie wyniku
    if(st == 0 || st == 3){
        result.push_back("\nNewtonRaphson(x,eps,it) = ");
        result.push_back(QString::fromStdString(x.toString()));
        result.push_back("\nfatx = ");
        result.push_back(QString::fromStdString(func(x).toString()));
        result.push_back("\n");
        result.push_back("it = ");
        result.push_back(QString::number(count));
        result.push_back(",st = ");
        result.push_back(QString::number(st)); // st
    return result;
    }
    else{
        result.clear();
        result.push_back("Dane:\n x = ");
        result.push_back(QString::fromStdString(x.toString(3)));
        result.push_back(", mit = ");
        result.push_back(QString::number(mit));
        result.push_back(", eps = 1e-16\n Równanie:\n");
        std::string funcName = funcname();
        qDebug() <<  QString::fromStdString(funcName);
        result.push_back(QString::fromStdString(funcName));
        result.push_back("\nWyniki:");
        result.push_back("st = ");
        result.push_back(QString::number(st));
        return result;
    }
}


