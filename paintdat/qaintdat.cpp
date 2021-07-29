#include "qaintdat.h"
#include "ui_qaintdat.h"
#include <QDebug>

extern const int sampleScale = 2; // 降采样倍率
extern const int posDev1 = 0;
extern const int posDev2 = 1174528/sampleScale;
extern const int size1 = 1174528/sampleScale;
extern const int size2 = 922624/sampleScale;

// 角膜后表面
extern const int posDev3 = 17600/sampleScale;
extern const int size3 = 47000/sampleScale;
extern const int dataPointSrc = 8*1024*1024; // (原数据)数据点
extern const int dataPoint = 8*1024*1024/4/sampleScale*4; // (降采样后)数据点
extern const int dataSize = 16*1024*1024; // (原数据)数据量

char *dataAcq = new char [dataSize];//dataSize
double *dataTrans = new double [dataPoint/4];
unsigned short *dataShort = new unsigned short [dataPoint/8];
double *dataPro = new double [dataPoint/4];//dataPoint/4

qaintdat::qaintdat(int index,int Nub,QWidget *parent) :
    QWidget(parent),
    ui(new Ui::qaintdat)
{
    ui->setupUi(this);

    setGeometry((index%6)*320,(index/6)*270,320,250);

    ui->widget->setGeometry(0,0,320,250);
    ssCalibrate(Nub);

    QVector<double> x(1024*1024), y(1024*1024); // initialize with entries 0..100
    for (int i=0; i<1024*1024; ++i)
    {
        x[i] = i; // x goes from -1 to 1
        y[i] = dataPro[i];  // let's plot a quadratic function
    }
    // create graph and assign data to it:
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x, y);
    // give the axes some labels:
    ui->widget->xAxis->setLabel("xx");
    ui->widget->yAxis->setLabel("yy");
    // set axes ranges, so we see all data:
    ui->widget->xAxis->setRange(0, 1024*1024);
    ui->widget->yAxis->setRange(-0.9, 0.9);
}

bool qaintdat::ssCalibrate(int Nub)
{

    qDebug() << "dataSize" << dataSize << "dataPoint" << dataPoint << "size1" << size1 << "dataPointSrc" << dataPointSrc;

    double front[104] = {0.0};
    double back[104] = {0.0};
    double y[104] = {0.0};

    int frontStartIndex = -1;
    int frontSize = 0;
    int backStartIndex = -1;
    int backSize = 0;

    QTime prev,cur;
    for(int i=Nub; i<Nub+1; i++)
    {
        cur = QTime::currentTime();
        // 数据采集
        memset(dataAcq,0,dataSize);
#if 1
        // 保存
        QString fileName = QString("C:\\Users\\Administrator\\Desktop\\9hao\\%1.dat").arg(i);
        QFile file(fileName);
        file.open(QFile::ReadOnly);
        file.read(dataAcq,dataSize);
        file.close();

        prev = cur;
        cur = QTime::currentTime();
//        qDebug() << "trans data cost" << prev.msecsTo(cur) << "ms";
#endif
        // 转换
        transfer(dataAcq,dataPro,dataPointSrc/4,sampleScale);
        // 找最大值
        double frontMax = 0.0;
        int frontPos = 0;
        double backMax = 0.0;
        int backPos = 0;

        // generate some data:
        QVector<double> x(1024*1024), y(1024*1024); // initialize with entries 0..100
        for (int i=0; i<1024*1024; ++i)
        {
            x[i] = i; // x goes from -1 to 1
            y[i] = dataPro[i];  // let's plot a quadratic function
        }
//        // create graph and assign data to it:
//        ui->centralWidget->addGraph();
//        ui->centralWidget->graph(0)->setData(x, y);
//        // give the axes some labels:
//        ui->centralWidget->xAxis->setLabel("xx");
//        ui->centralWidget->yAxis->setLabel("yy");
//        // set axes ranges, so we see all data:
//        ui->centralWidget->xAxis->setRange(0, 1024*1024);
//        ui->centralWidget->yAxis->setRange(-0.4, 0.4);

#if 0
        transfer1(dataAcq,dataShort,dataPointSrc/4,sampleScale);


        unsigned short frontMaxShort = 0;
        for(;;)
        {
            for(int j=0; j<dataPoint/8; j++)//size1
            {
                if(dataShort[j] > frontMaxShort)
                {
                    frontMaxShort = dataShort[j];
                    frontPos = j;
                }
            }
            if(frontMaxShort > 0.1)
            {
                front[104-i] = frontPos;
                frontStartIndex = 104-i;
                frontSize++;
            }
            qDebug() << "i>45&&i<104 i frontMax fontPos frontStartIndex frontSize" << i << frontMaxShort << frontPos << frontStartIndex << frontSize;
        }
#endif
        // 带通滤波
//        OpticsAlgorithm::filtfilt(dataTrans,dataPro,dataPoint/4,a1,b1,11);
        // 取绝对值
        for(int i=0; i<dataPoint/4; i++)
            dataTrans[i] = fabs(dataPro[i]);

        for(int j=0; j<1024*1024; j++)//size1
        {
            if(dataTrans[j] > frontMax)
            {
                frontMax = dataTrans[j];
                frontPos = j;
            }
        }
        if(frontMax > 0.1)
        {
            front[104-i] = frontPos;
            frontStartIndex = 104-i;
            frontSize++;
        }
        qDebug() << "i frontMax frontPos" << i << frontMax << "\t" << frontPos;
//        qDebug() << "i>45&&i<104 i frontMax fontPos frontStartIndex frontSize" << i << frontMax << frontPos << frontStartIndex << frontSize;

        //        for(int i=1 ;i < 1024*512;i+2)
        //        {
        //            int y1 = (int)(*(dataTrans+i)*1000);
        //            int y2 = (int)(*(dataTrans+i+1)*1000);
        //            painter.drawLine(i,y1,i-1,y2);
        //        }

//        update();

        // 低通滤波
//        OpticsAlgorithm::filtfilt(dataTrans,dataPro,dataPoint/4,a2,b2,5);

/*
        if(i>=45 && i<=104)
        {
            for(int j=0; j<1024*1024; j++)//size1
            {
                if(dataPro[j] > frontMax)
                {
                    frontMax = dataPro[j];
                    frontPos = j;
                }
            }
            if(frontMax > 0.1)
            {
                front[104-i] = frontPos;
                frontStartIndex = 104-i;
                frontSize++;
            }
            qDebug() << "i>45&&i<104 i frontMax fontPos frontStartIndex frontSize" << i << frontMax << frontPos << frontStartIndex << frontSize;
        }
        if(i>=1 && i<=60)
        {
            for(int j=0; j<1024*1024; j++)//size1    dataPoint/4
            {
                if(dataPro[j] > backMax)
                {
                    backMax = dataPro[j];
                    backPos = j;
                }
            }
            if(backMax > 0.1)
            {
                back[104-i] = backPos;
                backStartIndex = 104-i;
                backSize++;
            }
            qDebug() << "i>1&&i<60 i backMax backPos backStartIndex backSize" << i << backMax << backPos << backStartIndex << backSize;
        }

//        emit sigProgress(i);
        qApp->processEvents();

        if(i == 104)
            break;

        QThread::sleep(5);
    }

    //  剔除掉不符合的点
    int totalPoint = 104;
    qDebug() << "1剔除掉不符合的点 frontStartIndex" << frontStartIndex << endl;
    for (int k = 0 ; k < 104;k++)
    {
        qDebug() << front[k];
    }
    for(int i=frontStartIndex; i<totalPoint; i++)
        y[i] = 500*(i-frontStartIndex);
    for(int i=0; i<104; i++)
        qDebug() << y[i];

    for(int i=1; i<totalPoint; i++)
    {
        int removeIndex = -1;
        if((front[i]>0 && front[i]<=front[i-1]) || (i+1<totalPoint && front[i]==0 && front[i+1]!=0))
        {
            removeIndex = i;
        }
        else if((back[i]>0 && back[i]<=back[i-1]) || (i+1<totalPoint && back[i]==0 && back[i+1]!=0))
        {
            removeIndex = i;
        }
        if(back[i]==0 && front[i]==0)
            removeIndex = i;
        if(removeIndex >= 0)
        {
            // 后面的数据前移
            for(int j=removeIndex; j<totalPoint-1; j++)
            {
                front[j] = front[j+1];
                back[j] = back[j+1];
                y[j] = y[j+1];
            }
            if(removeIndex>=frontStartIndex && removeIndex<frontStartIndex+frontSize)
                frontSize--;
            if(removeIndex>=backStartIndex && removeIndex<backStartIndex+backSize)
                backSize--;
            if(removeIndex < frontStartIndex)
                frontStartIndex--;
            if(removeIndex < backStartIndex)
                backStartIndex--;
            totalPoint--;
            i--;
        }
    }

    for(int i=0; i<totalPoint; i++)
        qDebug() << "front back y"<< front[i] << back[i] << y[i];
//        for(int i=0; i<104; i++)
//            qDebug() << y[i];

    double ps1[4] = {0.0};
    double ps2[4] = {0.0};
    qDebug() << "front" << frontStartIndex << frontStartIndex+frontSize;
    qDebug() << "back" << backStartIndex << backStartIndex+backSize;
    if(frontSize<10 || backSize<10)
        return false;
    OpticsAlgorithm::polyfit(front+frontStartIndex,y+frontStartIndex,frontSize,3,ps1);
    OpticsAlgorithm::polyfit(back+backStartIndex,y+backStartIndex,backSize,3,ps2);

    for(int i=0; i<4; i++)
    {
        qDebug() << "p1" << ps1[i] << "p2" << ps2[i];
        p1[i] = ps1[3-i];
        p2[i] = ps2[3-i];

    }
    // save to file
    QFile file("./positionCalibrate1.txt");
    file.open(QFile::ReadWrite);
    QByteArray arr = file.readAll();
    QList<QByteArray> list = arr.split('\n');
    if(arr.contains("p1") && arr.contains("p2"))
    {
        for(int i=list.size()-1; i>=0; i--)
        {
            const QByteArray &a = list.at(i);
            if(a.contains("p1") || a.contains("p2"))
                list.removeAt(i);
        }
    }
    list.append(QString("p1=%1 %2 %3 %4").arg(p1[0]).arg(p1[1]).arg(p1[2]).arg(p1[3]).toLocal8Bit());
    list.append(QString("p2=%1 %2 %3 %4").arg(p2[0]).arg(p2[1]).arg(p2[2]).arg(p2[3]).toLocal8Bit());
    file.resize(0);
    foreach(QByteArray a,list)
    {
        file.write(a);
        file.write("\n");
    }
    file.flush();
    file.close();*/
}
//    delete [] dataAcq;
//    delete [] dataTrans;
//    delete [] dataPro;

    return true;
}

void qaintdat::transfer(char *buffer, double *bufTr, long size, int sampleScale) // size 数据点
{
//    for(int i=0; i<size; i++)
//    {
//        unsigned char chLow = buffer[i*2];
//        unsigned char chHigh = buffer[i*2+1];
//        unsigned short vol = (chHigh<<8)|chLow;
//        bufTr[i] = ((double)((vol+0x8000)&0xFFFF)/0xFFFF)*10-5;
//    }

    unsigned short *buf = (unsigned short *)buffer;
    for(int i=0; i<size/sampleScale/2; i++)
    {
        bufTr[i] = ((double)((buf[i*sampleScale]+0x8000)&0xFFFF)/0xFFFF)*10-5;
    }
}

qaintdat::~qaintdat()
{
    delete ui;
}
