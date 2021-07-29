#ifndef QAINTDAT_H
#define QAINTDAT_H

#include <QWidget>

namespace Ui {
class qaintdat;
}

class qaintdat : public QWidget
{
    Q_OBJECT

public:
    explicit qaintdat(int index,int Nub,QWidget *parent = nullptr);
    ~qaintdat();
    bool ssCalibrate(int Nub);
    void transfer(char *buffer, double *bufTr, long size, int sampleScale);
private:
    Ui::qaintdat *ui;
};

#endif // QAINTDAT_H
