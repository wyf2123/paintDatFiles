#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    int index = 0;
    index = ui->lineEdit->text().toInt();
    int k = 0;
    for(int i = index;(i < index+24)&&i < 105;i++)
    {
        paintdat[i] = new qaintdat(k,i);
        k++;
        paintdat[i]->show();
    }
}

void MainWindow::on_pushButton_2_clicked()
{
    int index = 0;
    index = ui->lineEdit->text().toInt();
    for (int i = index;(i < index+24)&&i < 105;i++) {
        paintdat[i]->close();
        delete paintdat[i];
    }
}
