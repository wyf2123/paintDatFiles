#include "algorithm.h"
#include <QDebug>
#include <QImage>
#include <vector>
#include <cmath>
#include <QFileInfo>
#include <QDir>
#include <QTime>

#define EPS 0.000001

const double diopt_n  = 1.3375;
const double Pi = 3.1415926;
const double &PI = Pi;
//const double scaleRate = 0.0102825;
const double scaleRate = 0.0069;
extern const int b0 = 82;
extern const int h_in = 21;
extern const int h_out = 26;

extern double position[4];
extern double cct_c[2]; // c_in c_out
extern double cct_k[2]; // k1 k2

extern int posDevX, posDevY;

//void saveImage(const QString &fileName, const Mat &img)
//{
//    // qDebug() << fileName;

//    QFileInfo info(fileName);
//    QDir dir(info.absolutePath());
//    // qDebug() << info.absolutePath();
//    if(!dir.exists())
//        system(QString("mkdir -p %1").arg(info.absolutePath()).toStdString().c_str());
//    imwrite(fileName.toStdString(),img);
//}

int OpticsAlgorithm::getIterativeBestThreshold(int histGram[])
{
    int X, Iter = 0;
    int MeanValueOne, MeanValueTwo, SumOne, SumTwo, SumIntegralOne, SumIntegralTwo;
    int MinValue, MaxValue;
    int Threshold, NewThreshold;

    for (MinValue = 0; MinValue < 256 && histGram[MinValue] == 0; MinValue++) ;
    for (MaxValue = 255; MaxValue > MinValue && histGram[MaxValue] == 0; MaxValue--) ;

    if (MaxValue == MinValue) return MaxValue;          // 图像中只有一个颜色
    if (MinValue + 1 == MaxValue) return MinValue;      // 图像中只有二个颜色

    Threshold = MinValue;
    NewThreshold = (MaxValue + MinValue) >> 1;
    while (Threshold != NewThreshold)    // 当前后两次迭代的获得阈值相同时，结束迭代
    {
        SumOne = 0; SumIntegralOne = 0;
        SumTwo = 0; SumIntegralTwo = 0;
        Threshold = NewThreshold;
        for (X = MinValue; X <= Threshold; X++)         //根据阈值将图像分割成目标和背景两部分，求出两部分的平均灰度值
        {
            SumIntegralOne += histGram[X] * X;
            SumOne += histGram[X];
        }

        MeanValueOne = SumIntegralOne / SumOne;
        for (X = Threshold + 1; X <= MaxValue; X++)
        {
            SumIntegralTwo += histGram[X] * X;
            SumTwo += histGram[X];
        }

        MeanValueTwo = SumIntegralTwo / SumTwo;
        NewThreshold = (MeanValueOne + MeanValueTwo) >> 1;     //求出新的阈值
        Iter++;
        if (Iter >= 1000) return -1;
    }
    //qDebug() << "threshold" << Threshold;
    return Threshold;
}

//void OpticsAlgorithm::contoursFilter(vector< vector<Point> > &contours,double minArea)
//{
//    vector< vector<Point> >::iterator it;
//    for(it=contours.begin();it!=contours.end();)
//    {
//        if(contourArea(*it) < minArea)
//            it = contours.erase(it);
//        else
//            ++it;
//    }
//}

//bool OpticsAlgorithm::calculateCornealCurvature(const Mat &img, CornealCurvature &curve, QString saveFileName, int index)
//{
//    qDebug() << __FUNCTION__ << "in";
//    int x = (img.cols-350)/2+posDevX;
//    int y = (img.rows-350)/2+posDevY;
//    Mat image;
//    if(x>=0 && x+350<=img.cols && y>=0 && y+350<=img.rows)
//        image = img(Rect((img.cols-350)/2+posDevX,(img.rows-350)/2+posDevY,350,350));
//    else
//        image = img(Rect((img.cols-350)/2,(img.rows-350)/2,350,350));

//    // 旋转180度
//    flip(image, image, -1);
//    // 保存图片
//    saveImage(saveFileName,image);
//    qDebug() << "save image";
//    // 眼轴长数据找不到前表面
//    double pos = position[index];
//    if(pos == 0)
//        return false;

//    // 像素灰度统计
//    //    int pixelCount[256] = {0};
//    //    for(int i=0; i<image.rows; i++)
//    //    {
//    //        for(int j=0; j<image.cols; j++)
//    //        {
//    //            pixelCount[image.at<uchar>(i,j)]++;
//    //        }
//    //    }
//    //    for(int i=0; i<256; i++)
//    //        qDebug() << i << pixelCount[i];

//    // 二值化--迭代最佳阈值
//    // int thresh = getIterativeBestThreshold(pixelCount);
//    Mat binary;
//    threshold(image,binary,255,255.0,CV_THRESH_OTSU);
//    qDebug() << "binary";
//    // 查找轮廓
//    vector< vector<Point> > contours;
//    findContours(binary, contours, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);

//    // 筛选轮廓，剔除噪点
//    contoursFilter(contours,10);

//    qDebug() << "find " << contours.size() << "contours";

//    // 轮廓不足30个--眼睑下垂或睫毛遮挡严重
//    if(contours.size() < 30)
//    {
//        qDebug() << "eyelid droop";
//        return false;
//    }

//    // 获得轮廓中心
//    vector<PointDistance> contourMidpoints;
//    int fitCnt = 0;
//    for(size_t i=0; i<contours.size(); i++)
//    {
//        Rect rec = boundingRect(contours[i]);
//        PointDistance pd;
//        pd.point = Point(rec.x+rec.width/2,rec.y+rec.height/2);
//        contourMidpoints.push_back(pd);

//        double whr = (double)rec.width/rec.height;
//        if(whr>0.77 && whr<1.3)
//            fitCnt++;
//    }
//    // 判断轮廓的宽高比是否符合预期
//    if(fitCnt < 16)
//    {
//        qDebug() << "符合宽高比[0.9,1.11]的轮廓不足16个，当前" << fitCnt;
//        return false;
//    }
//    // 计算所有点的算术中心
//    int totalX=0, totalY=0;
//    for(size_t i=0; i<contourMidpoints.size(); i++)
//    {
//        totalX += contourMidpoints[i].point.x;
//        totalY += contourMidpoints[i].point.y;
//    }
//    Point MidPoint(totalX/contourMidpoints.size(),totalY/contourMidpoints.size());
//    // 若中心点离图片中心点距离超过15，则取它们的中点作为中心点
//    if(abs(MidPoint.x-image.cols/2)>15 || abs(MidPoint.y-image.rows/2)>15)
//    {
//        MidPoint.x = (MidPoint.x+image.cols/2)/2;
//        MidPoint.y = (MidPoint.y+image.rows/2)/2;
//    }
//    // 计算离轮廓中心的距离
//    for(size_t i=0; i<contourMidpoints.size(); i++)
//    {
//        contourMidpoints[i].distance = twoPointsDistance(contourMidpoints[i].point,MidPoint);
//    }
//    // 内外圈分组--按距离排序
//    std::sort(contourMidpoints.begin(),contourMidpoints.end(),sortPoint);

//    // 剔除中心固视灯
//    contourMidpoints.erase(contourMidpoints.begin());

//    vector<Point> innerCircle; // 内圈
//    vector<Point> outerCircle; // 外圈
//    for(size_t i=0; i<16; i++)
//    {
//        innerCircle.push_back(contourMidpoints[i].point);
//        if(i+16 < contourMidpoints.size())
//            outerCircle.push_back(contourMidpoints[i+16].point);
//    }
//    // 椭圆拟合
//    RotatedRect innerEllipse = fitEllipse(innerCircle);
//    RotatedRect outerEllipse = fitEllipse(outerCircle);

//    //    qDebug() << "angle " << innerEllipse.angle << outerEllipse.angle;
//    //    qDebug() << "size " << innerEllipse.size.width << innerEllipse.size.height << outerEllipse.size.width << outerEllipse.size.height;

//    double innerAngleL = innerEllipse.angle;
//    if(innerAngleL>135) innerAngleL -= 180;
//    double innerAngleS = innerAngleL>=90?innerAngleL-90:innerAngleL+90;

//    int angleL = round(innerAngleL);
//    int angleS = round(innerAngleS);

//    // 长轴轴位为水平，短轴轴位为垂直
//    curve.horizontalMeridianDegree = angleL;
//    curve.verticalMeridianDegree = angleS;
//    curve.astigmatismDegree = angleS;

//    float innerAxisL = MAX(innerEllipse.size.width,innerEllipse.size.height)/2; // 内圈长轴半径
//    float innerAxisS = MIN(innerEllipse.size.width,innerEllipse.size.height)/2; // 内圈短轴半径
//    float outerAxisL = MAX(outerEllipse.size.width,outerEllipse.size.height)/2; // 外圈长轴半径
//    float outerAxisS = MIN(outerEllipse.size.width,outerEllipse.size.height)/2; // 外圈短轴半径

//    qDebug() << "轴长" << innerAxisL << innerAxisS << outerAxisL << outerAxisS;

//    double h0Lin = 2.2*innerAxisL/1000;
//    double h0Sin = 2.2*innerAxisS/1000;
//    double h0Lout = 2.2*outerAxisL/1000;
//    double h0Sout = 2.2*outerAxisS/1000;

//    double c_in = cct_c[0];
//    double c_out = cct_c[1];
//    double beta_in = 4e-6*pos*pos-8e-4*pos+c_in;
//    double beta_out = 4e-6*pos*pos-8e-4*pos+c_out;

//    double innerRL = 2*(b0+pos)*h0Lin/(beta_in*h_in-2*h0Lin);  // 内圈长轴角膜曲率半径
//    double innerRS = 2*(b0+pos)*h0Sin/(beta_in*h_in-2*h0Sin);  // 内圈短轴角膜曲率半径
//    double outerRL = 2*(b0+pos)*h0Lout/(beta_out*h_out-2*h0Lout);  // 外圈长轴角膜曲率半径
//    double outerRS = 2*(b0+pos)*h0Sout/(beta_out*h_out-2*h0Sout);  // 外圈短轴角膜曲率半径

//    qDebug() << "position" << pos;
//    qDebug() << "半径" << innerRL << innerRS << outerRL << outerRS;
//    qDebug() << "c_in" << c_in << "c_out" << c_out;
//    qDebug() << "k1 k2" << cct_k[0] << cct_k[1];
//    qDebug() << "beta_in" << beta_in << "beta_out" << beta_out;
//    qDebug() << "h0" << h0Lin << h0Sin << h0Lout << h0Sout;

//    innerRL = cct_k[0]*innerRL+cct_k[1];
//    innerRS = cct_k[0]*innerRS+cct_k[1];
//    outerRL = cct_k[0]*outerRL+cct_k[1];
//    outerRS = cct_k[0]*outerRS+cct_k[1];

//    double innerKL = 1000*(diopt_n-1)/innerRL; // 内圈长轴屈光度
//    double innerKS = 1000*(diopt_n-1)/innerRS; // 内圈短轴屈光度
//    double outerKL = 1000*(diopt_n-1)/outerRL; // 外圈长轴屈光度
//    double outerKS = 1000*(diopt_n-1)/outerRS; // 外圈短轴屈光度

//    qDebug() << "屈光度" << innerKL << innerKS << outerKL << outerKS;

//    double KL = (innerKL*2+outerKL)/3;
//    double KS = (innerKS*2+outerKS)/3;

//    curve.horizontalMeridianDiopter = KL;
//    curve.verticalMeridianDiopter = KS;
//    curve.astigmatismDiopter = KS-KL;

//    qDebug() << "角膜曲率" << curve.horizontalMeridianDiopter << curve.verticalMeridianDiopter;

//    // 35-56D筛除
//    double hor = curve.horizontalMeridianDiopter;
//    double ver = curve.verticalMeridianDiopter;
//    if(hor<35 || hor>56 || ver<35 || ver>56)
//        return false;

//    // 释放矩阵资源
//    binary.release();

//    return true;
//}

//short getGray(Mat img, Point center, double degree, double radius)
//{
//    int row = round(center.y-sin(degree/180*Pi)*radius);
//    int col = round(center.x+cos(degree/180*Pi)*radius);
//    return img.at<uchar>(row,col);
//}

//bool OpticsAlgorithm::roundDetection(const Mat &img, const QString &saveFileName, CalculateType type, double &diameter, double &barycenterDX, double &barycenterDY, double &kappaX, double &kappaY)
//{
//    Mat image;
//    cv::resize(img,image,Size(324,243));

//    QTime t1,t2;
//    t1 = QTime::currentTime();

//    // 旋转180度
//    // transpose(image, image);
//    flip(image, image, -1);

//    // 亮度大于220的点坐标累加
//    int posXSum = 0;
//    int posYSum = 0;
//    int posCnt = 0;

//    uchar brightest = 0;
//    for(int row=0; row<image.rows; row++)
//    {
//        for(int col=0; col<image.cols; col++)
//        {
//            // 亮度大于220的点坐标累加
//            if(image.at<uchar>(Point(col,row)) > 220)
//            {
//                posXSum += col;
//                posYSum += row;
//                posCnt++;
//            }
//            // 找最大灰度值
//            uchar t = image.at<uchar>(Point(col,row));
//            if(t > brightest)
//            {
//                brightest = t;
//            }
//        }
//    }

//    // 图片增亮
//    // 灰度对数平均值
//    if(type == TypePupilDiameter)
//    {
//        double Iw = 0;
//        for(int row=0; row<image.rows; row++)
//        {
//            for(int col=0; col<image.cols; col++)
//            {
//                Iw += log(image.at<uchar>(Point(col,row))+0.0001);
//            }
//        }
//        Iw = exp(Iw/(image.rows*image.rows));
//        qDebug() << "brightest" << brightest << "Iw" << Iw;
//        for(int row=0; row<image.rows; row++)
//        {
//            for(int col=0; col<image.cols; col++)
//            {
//                uchar &pos = image.at<uchar>(Point(col,row));
//                pos = log(pos/Iw+1)/log(brightest/Iw+1)*255;
//            }
//        }
//    }

//    t2 = QTime::currentTime();
//    qDebug() << "图片调亮" << t1.msecsTo(t2) << "ms";
//    t1 = t2;

//    saveImage(saveFileName,image);

//    t2 = QTime::currentTime();
//    qDebug() << "保存图片" << t1.msecsTo(t2) << "ms";
//    t1 = t2;

//    //获取自定义核
//    Mat element = getStructuringElement(MORPH_RECT, Size(9, 9));
//    //高级形态学处理--闭运算
//    Mat imageClose;
//    morphologyEx(image, imageClose, MORPH_CLOSE, element);

//    Mat imageGradient(imageClose.rows,imageClose.cols,CV_8UC1); // 梯度图

//    for(int i=0; i<imageClose.rows; i++)
//    {
//        for(int j=0; j<imageClose.cols; j++)
//        {
//            short diff;
//            if(j==0 || j==imageClose.cols-1 || i==0 || i==imageClose.rows-1)
//            {
//                continue;
//            }
//            else
//            {
//                short l = imageClose.at<uchar>(Point(j-1,i));
//                short r = imageClose.at<uchar>(Point(j+1,i));
//                short u = imageClose.at<uchar>(Point(j,i-1));
//                short d = imageClose.at<uchar>(Point(j,i+1));
//                diff = sqrt((l-r)*(l-r)+(u-d)*(u-d));
//                // qDebug() << "diff" << diff;
//            }
//            imageGradient.at<uchar>(i,j) = diff;
//        }
//    }

//    t2 = QTime::currentTime();
//    qDebug() << "生成梯度图" << t1.msecsTo(t2) << "ms";
//    t1 = t2;

//    Point bestCenter;
//    double bestRadius = 65;
//    double maxSum = 0;
//    double d1 = -20;
//    double d2 = 40;
//    double d3 = 140;
//    double d4 = 200;
//    double gap = 2;
//    double r1 = 42;
//    double r2 = 95;

//    Point center;
//    if(type == TypeWhiteToWhite)
//    {
//        d1 = -40;
//        d2 = 0;
//        d3 = 180;
//        d4 = 220;
//        r1 = 52;
//        r2 = 125;
//        center = Point(image.cols/2+posDevX/8,image.rows/2+posDevY/8);
//    }
//    else
//    {
//        // 灰度值大于220的点位置平均
//        if(posCnt > 0)
//            center = Point(posXSum/posCnt,posYSum/posCnt);
//        else
//            center = Point(image.cols/2+posDevX/4,image.rows/2+posDevY/4);
//    }

//    // 生成环形模板
//    Mat annular = Mat::zeros(image.rows, image.cols, CV_8UC1);
//    circle(annular, center, 130, Scalar(1), -1);
//    if(type == TypeWhiteToWhite)
//        circle(annular, center, 50, Scalar(0), -1);
//    else
//        circle(annular, center, 40, Scalar(0), -1);

//    // 梯度图与模板相乘
//    for(int i=0; i<imageGradient.rows; i++)
//    {
//        for(int j=0; j<imageGradient.cols; j++)
//        {
//            Point pt(j,i);
//            imageGradient.at<uchar>(pt) *= annular.at<uchar>(pt);
//        }
//    }

//    for(int x=center.x-15; x<=center.x+15; x++)
//    {
//        for(int y=center.y-10; y<=center.y+10; y++)
//        {
//            Point centerTmp(x,y);
//            for(double r=r1; r<r2; r+=1)
//            {
//                if(centerTmp.x-r<0 || centerTmp.x+r>image.cols-1 || centerTmp.y-r<0 || centerTmp.y+r>image.rows-1)
//                    break;

//                int sumL = 0;
//                for(double d=d1; d<=d2; d+=gap)
//                {
//                    sumL += getGray(imageGradient,centerTmp,d,r);
//                }

//                int sumR = 0;
//                for(double d=d3; d<=d4; d+=gap)
//                {
//                    sumR += getGray(imageGradient,centerTmp,d,r);
//                }

//                if(sumL>150 && sumR>150 && abs(sumL-sumR)<200 && sumL+sumR>maxSum)
//                {
//                    maxSum = sumL+sumR;
//                    bestCenter = centerTmp;
//                    bestRadius = r;
//                }
//            }
//        }
//    }

//    t2 = QTime::currentTime();
//    qDebug() << "找圆" << t1.msecsTo(t2) << "ms";

//    qDebug() << bestRadius << bestCenter.x << bestCenter.y;
//    diameter = bestRadius*2*8;
//    barycenterDX = (bestCenter.x-image.cols/2)*8;
//    barycenterDY = (bestCenter.y-image.rows/2)*8;

//    if(type == TypePupilDiameter)
//    {
//        // 计算眼轴光源中心坐标
//        Mat imageBin;
//        threshold(image,imageBin,200,255,THRESH_BINARY);
//        vector< vector<Point> > contours;
//        findContours(imageBin,contours, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);
//        int ledPosX = 162;
//        int ledPosY = 121;
//        if(contours.size() > 0)
//        {
//            // 轮廓坐标平均值
//            Point avg(0,0);
//            for(size_t i=0; i<contours.size(); i++)
//            {
//                Point &p = contours.at(i).at(0);
//                avg.x += p.x;
//                avg.y += p.y;
//            }
//            avg.x /= contours.size();
//            avg.y /= contours.size();
//            int nearestIdx = 0;
//            int nearestDst = 255;

//            for(size_t i=0; i<contours.size(); i++)
//            {
//                Point &p = contours.at(i).at(0);
//                int dist = (avg.x-p.x)*(avg.x-p.x)+(avg.y-p.y)*(avg.y-p.y);
//                if(dist < nearestDst)
//                {
//                    nearestIdx = i;
//                    nearestDst = dist;
//                }
//            }
//            Rect rect = boundingRect(contours.at(nearestIdx));
//            ledPosX = rect.x+rect.width/2;
//            ledPosY = rect.y+rect.height/2;
//        }

//        kappaX = (bestCenter.x-ledPosX)*8;
//        kappaY = (bestCenter.y-ledPosY)*8;

//        imageBin.release();
//    }

//    image.release();
//    element.release();
//    imageClose.release();
//    imageGradient.release();
//    annular.release();

//    return true;
//}

//bool OpticsAlgorithm::calculateWhiteToWhite(const Mat &img, WhiteToWhite &wtw, QString saveFileName)
//{
//    double diameter=0,barycenterDX=0, barycenterDY=0, ignore=0;
//    roundDetection(img,saveFileName,TypeWhiteToWhite,diameter,barycenterDX,barycenterDY,ignore,ignore);

//    wtw.whiteToWhite = diameter*scaleRate;
//    wtw.irisBarycenterDX = barycenterDX*scaleRate;
//    wtw.irisBarycenterDY = barycenterDY*scaleRate;
//    qDebug() << "wtw" << wtw.whiteToWhite;

//    return true;
//}

//bool OpticsAlgorithm::calculatePupilDiameter(const Mat &img, PupilData &pupil, QString saveFileName)
//{
//    double diameter=0,barycenterDX=0, barycenterDY=0, kappaX=0, kappaY=0;
//    roundDetection(img,saveFileName,TypePupilDiameter,diameter,barycenterDX,barycenterDY,kappaX,kappaY);

//    pupil.pupilDiameter = diameter*scaleRate/2;
//    pupil.pupilBarycenterDX = barycenterDX*scaleRate/2;
//    pupil.pupilBarycenterDY = barycenterDY*scaleRate/2;
//    pupil.kappaX = kappaX*scaleRate/2;
//    pupil.kappaY = kappaY*scaleRate/2;
//    qDebug() << "pupil" << pupil.pupilDiameter;

//    return true;
//}

//矩阵乘法  m为a的行数，n为a的列数数，k为b的行数，第一个矩阵列数必须要和第二个矩阵的行数相同
void OpticsAlgorithm::trmul(double *a,double *b,double *c,int m,int n,int k)
{
    int i,j,l,u;
    for (i=0; i<=m-1; i++)
        for (j=0; j<=k-1; j++)
        { u=i*k+j; c[u]=0.0;
            for (l=0; l<=n-1; l++)
                c[u]=c[u]+a[i*n+l]*b[l*k+j];
        }
    return;
}

//求逆矩阵，当返回值为0时成功，a变为逆矩阵
int OpticsAlgorithm::rinv(double *a,int n)
{ int *is,*js,i,j,k,l,u,v;
    double d,p;
    is=(int *)malloc(n*sizeof(int));
    js=(int *)malloc(n*sizeof(int));
    for (k=0; k<=n-1; k++)
    { d=0.0;
        for (i=k; i<=n-1; i++)
            for (j=k; j<=n-1; j++)
            { l=i*n+j; p=fabs(a[l]);
                if (p>d) { d=p; is[k]=i; js[k]=j;}
            }
        if (d+1.0==1.0)
        { free(is); free(js); printf("err**not inv\n");
            return(0);
        }
        if (is[k]!=k)
            for (j=0; j<=n-1; j++)
            { u=k*n+j; v=is[k]*n+j;
                p=a[u]; a[u]=a[v]; a[v]=p;
            }
        if (js[k]!=k)
            for (i=0; i<=n-1; i++)
            { u=i*n+k; v=i*n+js[k];
                p=a[u]; a[u]=a[v]; a[v]=p;
            }
        l=k*n+k;
        a[l]=1.0/a[l];
        for (j=0; j<=n-1; j++)
            if (j!=k)
            { u=k*n+j; a[u]=a[u]*a[l];}
        for (i=0; i<=n-1; i++)
            if (i!=k)
                for (j=0; j<=n-1; j++)
                    if (j!=k)
                    { u=i*n+j;
                        a[u]=a[u]-a[i*n+k]*a[k*n+j];
                    }
        for (i=0; i<=n-1; i++)
            if (i!=k)
            { u=i*n+k; a[u]=-a[u]*a[l];}
    }
    for (k=n-1; k>=0; k--)
    { if (js[k]!=k)
            for (j=0; j<=n-1; j++)
            { u=k*n+j; v=js[k]*n+j;
                p=a[u]; a[u]=a[v]; a[v]=p;
            }
        if (is[k]!=k)
            for (i=0; i<=n-1; i++)
            { u=i*n+k; v=i*n+is[k];
                p=a[u]; a[u]=a[v]; a[v]=p;
            }
    }
    free(is); free(js);
    return(1);
}

void OpticsAlgorithm::filter(const double* x, double* y, int xlen, double* a, double* b, int nfilt, double* zi)
{
    double tmp;
    int i,j;

    //normalization
    if( (*a-1.0>EPS) || (*a-1.0<-EPS) )
    {
        tmp=*a;
        for(i=0;i<nfilt;i++)
        {
            b[i]/=tmp;
            a[i]/=tmp;
        }
    }

    memset(y,0,xlen*sizeof(double));//将y清零，以双浮点为单位

    a[0]=0.0;
    for(i=0;i<xlen;i++)
    {
        for(j=0;i>=j&&j<nfilt;j++)
        {
            y[i] += (b[j]*x[i-j]-a[j]*y[i-j]);
        }
        if(zi&&i<nfilt-1) y[i] += zi[i];
    }
    a[0]=1.0;
}

int OpticsAlgorithm::filtfilt(double* x, double* y, int xlen, double* a, double* b, int nfilt){
    int nfact;
    int tlen;    //length of tx
    int i;
    double *tx,*tx1,*p,*t,*end;
    double *sp,*tvec,*zi;
    double tmp,tmp1;

    nfact=nfilt-1;    //3*nfact: length of edge transients

    if(xlen<=3*nfact || nfilt<2) return -1;
    //too short input x or a,b
    //Extrapolate beginning and end of data sequence using a "reflection
    //method". Slopes of original and extrapolated sequences match at
    //the end points.
    //This reduces end effects.
    tlen=6*nfact+xlen;
    tx=(double *)malloc(tlen*sizeof(double));
    tx1=(double *)malloc(tlen*sizeof(double));

    sp=(double *)malloc( sizeof(double) * nfact * nfact );
    tvec=(double *)malloc( sizeof(double) * nfact );
    zi=(double *)malloc( sizeof(double) * nfact );

    if( !tx || !tx1 || !sp || !tvec || !zi ){
        free(tx);
        free(tx1);
        free(sp);
        free(tvec);
        free(zi);
        return 1;
    }

    tmp=x[0];
    for(p=x+3*nfact,t=tx;p>x;--p,++t) *t=2.0*tmp-*p;
    for(end=x+xlen;p<end;++p,++t) *t=*p;
    tmp=x[xlen-1];
    for(end=tx+tlen,p-=2;t<end;--p,++t) *t=2.0*tmp-*p;
    //now tx is ok.

    end = sp + nfact*nfact;
    p=sp;
    while(p<end) *p++ = 0.0L; //clear sp
    sp[0]=1.0+a[1];
    for(i=1;i<nfact;i++){
        sp[i*nfact]=a[i+1];
        sp[i*nfact+i]=1.0L;
        sp[(i-1)*nfact+i]=-1.0L;
    }

    for(i=0;i<nfact;i++){
        tvec[i]=b[i+1]-a[i+1]*b[0];
    }

    if(rinv(sp,nfact)){
        free(zi);
        zi=NULL;
    }
    else{
        trmul(sp,tvec,zi,nfact,nfact,1);
    }//zi is ok

    free(sp);free(tvec);

    //filtering tx, save it in tx1
    tmp1=tx[0];
    if(zi)
        for( p=zi,end=zi+nfact; p<end;) *(p++) *= tmp1;
    filter(tx,tx1,tlen,a,b,nfilt,zi);

    //reverse tx1
    for( p=tx1,end=tx1+tlen-1; p<end; p++,end--){
        tmp = *p;
        *p = *end;
        *end = tmp;
    }

    //filter again
    tmp1 = (*tx1)/tmp1;
    if(zi)
        for( p=zi,end=zi+nfact; p<end;) *(p++) *= tmp1;
    filter(tx1,tx,tlen,a,b,nfilt,zi);

    //reverse to y
    end = y+xlen;
    p = tx+3*nfact+xlen-1;
    while(y<end)
    {
        *y++ = *p--;
    }

    free(zi);
    free(tx);
    free(tx1);

    return 0;
}

int OpticsAlgorithm::polyfit(const double *const independentVariables, const double *const dependentVariables, unsigned int countOfElements, unsigned int order, double *coefficients)
{
    // Declarations...
    // ----------------------------------
    enum {maxOrder = 5};

    double B[maxOrder+1] = {0.0f};
    double P[((maxOrder+1) * 2)+1] = {0.0f};
    double A[(maxOrder + 1)*2*(maxOrder + 1)] = {0.0f};

    double x, y, powx;

    unsigned int ii, jj, kk;

    // Verify initial conditions....
    // ----------------------------------

    // This method requires that the countOfElements >
    // (order+1)
    if (countOfElements <= order)
        return -1;

    // This method has imposed an arbitrary bound of
    // order <= maxOrder.  Increase maxOrder if necessary.
    if (order > maxOrder)
        return -1;

    // Begin Code...
    // ----------------------------------

    // Identify the column vector
    for (ii = 0; ii < countOfElements; ii++)
    {
        x    = independentVariables[ii];
        y    = dependentVariables[ii];
        powx = 1;

        for (jj = 0; jj < (order + 1); jj++)
        {
            B[jj] = B[jj] + (y * powx);
            powx  = powx * x;
        }
    }

    // Initialize the PowX array
    P[0] = countOfElements;

    // Compute the sum of the Powers of X
    for (ii = 0; ii < countOfElements; ii++)
    {
        x    = independentVariables[ii];
        powx = independentVariables[ii];

        for (jj = 1; jj < ((2 * (order + 1)) + 1); jj++)
        {
            P[jj] = P[jj] + powx;
            powx  = powx * x;
        }
    }

    // Initialize the reduction matrix
    //
    for (ii = 0; ii < (order + 1); ii++)
    {
        for (jj = 0; jj < (order + 1); jj++)
        {
            A[(ii * (2 * (order + 1))) + jj] = P[ii+jj];
        }

        A[(ii*(2 * (order + 1))) + (ii + (order + 1))] = 1;
    }

    // Move the Identity matrix portion of the redux matrix
    // to the left side (find the inverse of the left side
    // of the redux matrix
    for (ii = 0; ii < (order + 1); ii++)
    {
        x = A[(ii * (2 * (order + 1))) + ii];
        if (x != 0)
        {
            for (kk = 0; kk < (2 * (order + 1)); kk++)
            {
                A[(ii * (2 * (order + 1))) + kk] =
                        A[(ii * (2 * (order + 1))) + kk] / x;
            }

            for (jj = 0; jj < (order + 1); jj++)
            {
                if ((jj - ii) != 0)
                {
                    y = A[(jj * (2 * (order + 1))) + ii];
                    for (kk = 0; kk < (2 * (order + 1)); kk++)
                    {
                        A[(jj * (2 * (order + 1))) + kk] =
                                A[(jj * (2 * (order + 1))) + kk] -
                                y * A[(ii * (2 * (order + 1))) + kk];
                    }
                }
            }
        }
        else
        {
            // Cannot work with singular matrices
            return -1;
        }
    }

    // Calculate and Identify the coefficients
    for (ii = 0; ii < (order + 1); ii++)
    {
        for (jj = 0; jj < (order + 1); jj++)
        {
            x = 0;
            for (kk = 0; kk < (order + 1); kk++)
            {
                x = x + (A[(ii * (2 * (order + 1))) + (kk + (order + 1))] *
                        B[kk]);
            }
            coefficients[ii] = x;
        }
    }

    return 0;
}
