#ifndef ALGORITHM_H
#define ALGORITHM_H

//#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/core/core.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
//#include "mysqlitepatients.h"
#include <stdio.h>

extern const double scaleRate;

//using namespace cv;
//using std::sort;
//using std::vector;

//struct PointDistance
//{
//    Point point;
//    double distance;
//};

enum CalculateType {TypeCornealCurvature,TypeWhiteToWhite,TypePupilDiameter};

class OpticsAlgorithm
{
public:
    OpticsAlgorithm() {}
    ~OpticsAlgorithm() {}

    // 交换
    template<class T>
    inline static void swap(T &v1, T &v2)
    {
         T vt = v1;
         v1 = v2;
         v2 = vt;
    }

    // 计算平均值
//    inline static double average(const double *data, int size)
//    {
//        if(data==NULL || size<=0)
//            return 0;
//        double sum = 0;
//        for(int i=0; i<size; i++)
//            sum += data[i];
//        return sum/size;
//    }

//    inline static double average(double d1,double d2, double d3=0, int size=3)
//    {
//        double data[3] = {0};
//        data[0] = d1;
//        data[1] = d2;
//        data[2] = d3;
//        return average(data,size);
//    }

    // 计算标准偏差
//    inline static double standardDeviation(int size, double n1, double n2, double n3=0)
//    {
//        double data[3];
//        data[0] = n1;
//        data[1] = n2;
//        data[2] = n3;

//        return standardDeviation(data,size);
//    }

//    inline static double standardDeviation(const double *data, int size)
//    {
//        if(size <= 1)
//            return 0;
//        double sum = 0;
//        for(int i=0; i<size; i++)
//            sum += data[i];
//        double avg = sum/size;
//        double sumOfSquare = 0;
//        for(int i=0; i<size; i++)
//            sumOfSquare += (data[i]-avg)*(data[i]-avg);
//        return sqrt(sumOfSquare/(size-1));
//    }

    // 角度有效性检测
//    inline static void degreeIdentity(double *data, int size)
//    {
//        // 是否有相差大于90度
//        bool diff = false;
//        for(int i=0; i<size-1; i++)
//        {
//            if(abs(data[i]-data[i+1]) > 90)
//            {
//                diff = true;
//                break;
//            }
//        }
//        // 如果是，大于90度的值减去180
//        if(diff)
//        {
//            for(int i=0; i<size; i++)
//                if(data[i] > 90)
//                    data[i] -= 180;
//        }
//    }

    // 计算两点距离
//    inline static double twoPointsDistance(Point p1, Point p2)
//    {
//        return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2));
//    }

    // 将点按照距离排序的比较函数
//    inline static bool sortPoint(PointDistance pd1, PointDistance pd2)
//    {
//        return pd1.distance<pd2.distance;
//    }

    // 迭代最佳阈值法
    static int getIterativeBestThreshold(int histGram[]);

    // 轮廓大小过滤
//    static void contoursFilter(vector< vector<Point> > &contours,double minArea=50);

    // 角膜曲率算法
//    static bool calculateCornealCurvature(const Mat &image, CornealCurvature &curve, QString saveFileName, int index);

    // 圆检测算法，计算直径
//    static bool roundDetection(const Mat &img, const QString &saveFileName, CalculateType type, double &diameter, double &barycenterDX, double &barycenterDY, double &kappaX, double &kappaY);

    // 白到白
//    static bool calculateWhiteToWhite(const Mat &img, WhiteToWhite &wtw, QString saveFileName);

    // 计算瞳孔直径
//    static bool calculatePupilDiameter(const Mat &img, PupilData &pupil, QString saveFileName);

    inline static double DS(double AL, double K1, double K2)
    {
        return 81.915 - 1.953*AL - 0.569*K1 - 0.28*K2;
    }

    // 滤波函数
    static void trmul(double *a,double *b,double *c,int m,int n,int k);
    static int rinv(double *a,int n);
    static void filter(const double *x, double *y, int xlen, double *a, double *b, int nfilt, double *zi);
    static int filtfilt(double* x, double* y, int xlen, double* a, double* b, int nfilt);

    // 曲线拟合
    static int polyfit(const double *const independentVariables, const double *const dependentVariables, unsigned int countOfElements, unsigned int order, double *coefficients);
};

#endif // ALGORITHM_H
