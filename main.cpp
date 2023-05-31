#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <cstring>

using namespace std;
typedef double MyFloat;
const int   M=51,              // Кол-во фиктивных источников по трубе
            N=5,                // Max Кол-во фиктивных источников по аноду
            NN = 4*M-1;     // Число неизвестных/уравнений системы

//NN = 3;

const MyFloat pi = 3.141592653589793238462643383279,
            Rt = 0.557/2.0,     // Внешний радиус трубы
            Ht = 0.008,         // Толщина стенки трубы
            Zt = 1.5,           // Глубина залегания трубы
            Sms = pi*(Rt*Rt - (Rt-Ht)*(Rt-Ht)),
                                // Площадь металлического сечения трубы
            Ya = 200.0,         // Расстояние от анода до трубы
            Za = 2.5,           // Глубина точечного анода под поверхностью земли
            ro_t = 2.45e-7,     // 2.45e-7 Удельное сопр стали Ом*м
            ro_g = 35.0,        // Удельное сопр грунта Ом*м
            Ct = 5000.0,        // Удельное сопр изоляции трубы Ом*м2
            Sigma_t = 1/ro_t,   // Удельная электропроводность металла трубы
            Sigma_g = 1/ro_g,   // Удельная электропроводность грунта

            point_povr_Ct = 0.75, /// Относительная коор-та центра участка повреждения изоляции (0,1)
            koef_povr_Ct = 0.2,   /// Коэффициент поврежденности изоляции (0,1)
                                  /// 0 - полное ровреждение, 1 - без повреждений
            L_povr_Ct = 1000;       /// Длина участка с поврежденной изоляцией, м
/// на участке трубы на расстоянии L*point_povr_Ct от начала защищаемого участка задается участок
/// трубы длиной L_povr_Ct м, состояние изоляции трубы которого таково, что среднее значение
/// удельного электрического сопротивления изолции на данном участке равно Ct*koef_povr_Ct

/// !!! Пока не проработан механизм учитывающий, что задать повреждение мы можем только
/// для участка, равного по длине целому числу ФИ!!!


/// Длина ФИ м/б достаточно велика, и задать повреждение можем лишь для целого числа ФИ,
/// поэтому для обеспечения заданного условия вычислим кол-во ФИ с повреждением и
/// вычислим скорректированный коэфф-т неповрежденности изоляции
/// count_povr_FI = int(L_povr_Ct/Li);  - кол-во ФИ, для которых изменяется Ct.
/// составим пропорцию count_povr_FI*Li*x_Ct*Ct = L_povr_Ct*koef_povr_Ct*Ct,
/// откуда  x_Ct = (L_povr_Ct*koef_povr_Ct) / count_povr_FI*Li
/// Все это делаем в Init_FIt()
/// int((M-1)*point_povr_Ct - номер ФИ - центра повреждения
/// найдем номер первого ФИ на участке повреждения
/// first_povr_FI = int((M-1)*point_povr_Ct - int(count_povr_FI/2)
/// НЕОБХОДИМО ПРОАНАЛИЗИРОВАТЬ И ОБРАБОТАТЬ СЛУЧАИ, КОГДА
/// first_povr_FI ПОЛУЧИТСЯ ОТРИЦАТЕЛЬНЫМ, ИЛИ (first_povr_FI + count_povr_FI) > M
///    for (i=0;i < count_povr_FI;i++)
///    {
///        FI_t[first_povr_FI+i][3]= Ct*koef_povr_Ct);
///    }

MyFloat     //I0,               // Ток катодной станции
            Li,                 // Длина участка трубы соотв-го одному ФИ
            St,                 // Площадь боковой поверхности, соответствующая одному ФИ
            L,                  // Длина защищаемого участка трубы
            varR1,             // Выражение для потенциала цилиндрического электрода в ПЭА
//            C_t[M],           // удельные электропроводности изоляции трубы Ом*м2
            FI_t[M][4],         // Координаты центров ФИ по трубе
            FI_a[N][5],         // Координаты центров ФИ анодов
            Pi[3],              // коор-ты точек-центров ФИ
            Pj[3],              // коор-ты точек-центров ФИ
            A[NN][NN+1],        // Расширенная матрица системы
//            B[NN],              /// Столбец свободных членов - ИСПОЛЬЗУЕТСЯ??
            Nv;                 // Норма невязки системы
//            Pks_x,              /// координата Х точки подключения к трубе катодной станции - ИСПОЛЬЗУЕТСЯ??
// R_self; /// - ИСПОЛЬЗУЕТСЯ??

int         i_ks;               /// номер ФИ трубы, к центру которого подкючен контакт КС - ИСПОЛЬЗУЕТСЯ??

//////////// ГАУСС ////////////////////

int         num[NN];            // Массив для учета перестановки столбцов
                                // при выборе главного элемента в методе Гаусса
MyFloat     Ais[NN][NN+1],      // Массив для хранения исходной копии матрицы системы
                                // необходим для вычисления невязки
            X[NN],              // Вектор неизвестных системы
            Nev[NN];            // Вектор невязки системы


//////////// ГАУСС ////////////////////

/// //////// Импорт данных  /////////////////
const int D1_0 = 100;   // max кол-во строк данных о точках замера
const int D2_0 = 6;     // max кол-во столбцов данных для одной точки замера
int Dt,                 // реальное кол-во строк даных о точках замера
    Dskz;               // реальное кол-во строк даных о СКЗ (кол-во анодов)
MyFloat DataSKZ[D1_0][D2_0]; /// массив для импорта данных из файла

/// //////// Импорт данных  /////////////////

void MatrPrn_Nx2N(MyFloat A[NN][2*NN],int M,int N,char* fn)
/// Вывод в файл точечной матрицы Nx2N
{
    FILE *fp;
    int i,j;
    if((fp=fopen(fn, "w")) == NULL)
    {
        printf("%s %s \n","Not opened file",fn);
    }
    else
    {
        printf("%s %s \n","Opened file ",fn);
        fprintf(fp,"Расширенная матрица системы;\n");

        for (i=0;i<M;i++)
        {
            for (j=0;j<N;j++)
            {
                fprintf(fp,"%24.18f%s",A[i][j],"  ");
            }
            fprintf(fp,"\n");
        }
    }
}

void MatrPrn_NxNp1(MyFloat A[NN][NN+1],int M,int N,char* fn)
/// Вывод в файл точечной матрицы NxN
{
    FILE *fp;
    int i,j;
    if((fp=fopen(fn, "w")) == NULL)
    {
        printf("%s %s \n","Not opened file ",fn);
    }
    else
    {
        printf("%s %s \n","Opened file ",fn);
        fprintf(fp,"Расширенная матрица системы;\n");

        for (i=0;i<M;i++)
        {
            for (j=0;j<N;j++)
            {
                fprintf(fp,"%24.18f%s",A[i][j],"  ");
            }
            fprintf(fp,"\n");
        }
    }
}
/// //////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////
/// Интервалы begin
/// //////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////
#include <fenv.h>
typedef MyFloat Ifloat[2];
Ifloat  i_AB[NN][NN+1],   /// "Рабочая" расширенная матрица, приводимая к треугольной для м Гаусса
        i_ABis[NN][NN+1], /// Исходная расширенная матрица для вычисления невязки
        i_X[NN],          /// "Рабочий" вектор неизвестниых, компоненты переставляются при выборе главного эл-та
        i_Xis[NN];        /// Вектор неизвестниых для проверки невязки по исходной м-це

//        i_nev_is[NN],     /// Вектор невязки по исходной м-це
//        i_nev_tr[NN];     /// Вектор невязки по треугольной м-це

MyFloat i_AB_iw[NN][NN+1], /// Массив, содержащий значения ширины интервалов коэф-в треуг м-цы - для контроля
        i_nev_is[NN],     /// Вектор невязки по исходной м-це
        i_nev_tr[NN],     /// Вектор невязки по треугольной м-це
        i_dist_sol[NN];   /// Вектор для оценивани алгебраического решения



int i_num[NN];            /// массив для учета перенумерации переменных при выборе главного эл-та
char i_format[24] = "%s%24.18f%s%24.18f%s"; /// Формат вывода интервальных коэф-в в общий файл
char i_format1[28] = "%16.8f%s%24.18f%s%24.18f\n"; /// Формат вывода интервальных результатов в отдельные файлы для построения графиков


MyFloat disp_ro_g = 50, /// дисперсия, разброс знач уд сопр. грунта в % +/- к принятому выше значению
        disp_Ct = 20;  /// дисперсия, разброс знач уд сопр. изоляции трубы в % +/- к принятому выше значению

Ifloat i_ro_g = {(ro_g*(100-disp_ro_g)/100),(ro_g*(100+disp_ro_g)/100)},
       i_sig = {1/i_ro_g[1],1/i_ro_g[0]},
       i_Ct = {Ct*(100-disp_Ct)/100,Ct*(100+disp_Ct)/100}, /// одно значение Ст для всех ФИ, исп в первых версиях

       i_Ct_FI[M];   /// Для каждого ФИ свое Ст, интервальный вариант FI_t[i][3]

MyFloat AE[NN][2*NN], /// Матрица вида АЕ для нахождения обратной к А
        E[NN][NN+1];  /// изначально единичная м-ца


void i_movp(Ifloat r,MyFloat a)
{
    r[0] = a;
    r[1] = a;
}

void i_movi(Ifloat r,Ifloat a)
{
    r[0] = a[0];
    r[1] = a[1];
}

void i_add(Ifloat r,Ifloat a, Ifloat b)
{
    r[0] = a[0] + b[0];
    r[1] = -a[1] - b[1];
    r[1] = - r[1];
}
void i_sub(Ifloat r,Ifloat a, Ifloat b)
{
    r[0] = a[0] - b[1];
    r[1] = - a[1] + b[0];
    r[1] = - r[1];
}
void i_mul(Ifloat r,Ifloat a, Ifloat b)
{
    int i;
    MyFloat i_tl[4],i_tr[4],i_t_min,i_t_max; /// временная переменная,
    i_tl[0] = a[0] * b[0];
    i_tl[1] = a[0] * b[1];
    i_tl[2] = a[1] * b[0];
    i_tl[3] = a[1] * b[1];

    i_tr[0] = - a[0] * b[0];
    i_tr[1] = - a[0] * b[1];
    i_tr[2] = - a[1] * b[0];
    i_tr[3] = - a[1] * b[1];

    i_tr[0] = - i_tr[0];
    i_tr[1] = - i_tr[1];
    i_tr[2] = - i_tr[2];
    i_tr[3] = - i_tr[3];

    i_t_min = i_tl[0];
    i_t_max = i_tr[0];
    for (i=0;i<3;i++)
    {
        if (i_tl[i] < i_t_min) i_t_min = i_tl[i];
        if (i_tr[i] > i_t_max) i_t_max = i_tr[i];
    }
    r[0] = i_t_min;
    r[1] = i_t_max;
}

bool i_incl0(Ifloat a)
{
    if (a[0]*a[1] <= 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void i_div(Ifloat r,Ifloat a,Ifloat b)
{
    int i;
    if (i_incl0(b))
    {
        cout<<"division by zero included interval";
        cin>>i;
    }
    else
    {
        MyFloat i_tl[4],i_tr[4],i_t_min,i_t_max; /// временная переменная,
        i_tl[0] = a[0] / b[0];
        i_tl[1] = a[0] / b[1];
        i_tl[2] = a[1] / b[0];
        i_tl[3] = a[1] / b[1];

        i_tr[0] = - a[0] / b[0];
        i_tr[1] = - a[0] / b[1];
        i_tr[2] = - a[1] / b[0];
        i_tr[3] = - a[1] / b[1];

        i_tr[0] = - i_tr[0];
        i_tr[1] = - i_tr[1];
        i_tr[2] = - i_tr[2];
        i_tr[3] = - i_tr[3];

        i_t_min = i_tl[0];
        i_t_max = i_tr[0];
        for (i=0;i<3;i++)
        {
            if (i_tl[i] < i_t_min) i_t_min = i_tl[i];
            if (i_tr[i] > i_t_max) i_t_max = i_tr[i];
        }
        r[0] = i_t_min;
        r[1] = i_t_max;
    }


}

MyFloat i_abs(Ifloat a)
{
    MyFloat abt;
    abt = abs(a[0]);
    if (abt < abs(a[1]))
    {
        abt = abs(a[1]);
    }
    return abt;
}

MyFloat i_dist(Ifloat a, Ifloat b)
{
    Ifloat t;
    t[0] = abs(a[0] - b[0]);
    t[1] = - a[1] + b[1];
    t[1] = - t[1];
    return i_abs(t);
}


MyFloat i_wid(Ifloat a)
{
    return a[1]-a[0];
}


void i_Gauss()
{
    int l,i,j,t,k,max_num;
    Ifloat Li1,tmp,sum;
    for (l=0;l<NN;l++)
{
/// поиск макс элемента
    max_num = l;
    for (k=l+1;k<NN;k++)
    {
///  МЫСЛИ О МОДЕРНИЗАЦИИ ВЫБОРА ГЛАВНОГО ЭЛЕМЕНТА С УЧЕТОМ ШИРИНЫ ИНТЕРВАЛА ВЕДУЩЕГО ЭЛЕМЕНТА
///        cout<<width_int(i_AB[l][max_num])<<endl;

        if (i_abs(i_AB[l][k]) > i_abs(i_AB[l][max_num])) max_num = k;
//        if ((i_abs(i_AB[l][k]) > i_abs(i_AB[l][max_num])) and (!i_incl0(i_AB[l][k]))) max_num = k;
    }
        if (max_num != l)
        {
//            меняем столбцы
            t = i_num[l];
            i_num[l] = i_num[max_num];
            i_num[max_num] = t;

            for (t = 0;t<NN;t++)
            {
                tmp[0] = i_AB[t][l][0];
                tmp[1] = i_AB[t][l][1];

                i_AB[t][l][0] = i_AB[t][max_num][0];
                i_AB[t][l][1] = i_AB[t][max_num][1];

                i_AB[t][max_num][0] = tmp[0];
                i_AB[t][max_num][1] = tmp[1];
            }
        }

/// поиск макс элемента
        for (i=l+1;i<NN;i++)
        {
            i_div(Li1,i_AB[i][l],i_AB[l][l]);

            for (j=0;j<NN+1;j++)
            {
                i_mul(tmp,Li1,i_AB[l][j]);
                i_sub(i_AB[i][j],i_AB[i][j],tmp);
            }
        }
}

///        AB[NN][NN+1] = AB[NN][NN+1]/AB[NN][NN];
        i_div(i_X[NN-1],i_AB[NN-1][NN],i_AB[NN-1][NN-1]);
        for (i=NN-2;i>=0;i--)
        {
            sum[0] = 0;
            sum[1] = 0;
///            for (j=NN;j=i+1;j--) sum = sum + AB[i][j]*AB[j][NN+1];
            for (j=NN-1;j>=i;j--)
            {
                i_mul(tmp,i_AB[i][j],i_X[j]);
                i_add(sum,sum,tmp);
            }
///            AB[i][NN+1] = (AB[i][NN+1] - sum)/(AB[i][i]);
            i_sub(tmp,i_AB[i][NN],sum);
            i_div(i_X[i],tmp,i_AB[i][i]);
        }
    for (i=0;i<NN;i++)
    {
       i_Xis[i_num[i]][0] = i_X[i][0];
       i_Xis[i_num[i]][1] = i_X[i][1];
    }

        /// ВЫЧИСЛЯЕМ ШИРИНУ ИНТЕРВАЛОВ-КОЭФФИЦИЕНТОВ ТРЕУГОЛЬНОЙ МАТРИЦЫ
    fesetround (FE_UPWARD);
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            i_AB_iw[i][j] = i_wid(i_AB[i][j]);
        }
    }
    fesetround (FE_DOWNWARD);
 }


void i_printResult()
{
    FILE *fp;
//        dirName[] = "\\1111\\",
char    dirName[9] = "",
    fileType[9] = ".txt",
    fileName_Res[9] = "i_Resall",
    fileName_Itx[9] = "i_Itx",
    fileName_Itg[9] = "i_Itg",
    fileName_Utg[9] = "i_Utg",
    fileName_Utm[9] = "i_Utm",
    fileName_Upr[9] = "i_Upr",
    fileName_Jtg[9] = "i_Jtg",
    curFileName[9] = "";
    int i,j;
    Ifloat t;
//    char filename[100] = __TIME__;
//    cout<<filename<<endl;
//    ofstream fout(filename);
    strcat(curFileName,dirName);
    strcat(curFileName,fileName_Res);
    strcat(curFileName,fileType);
    if((fp=fopen(curFileName, "w")) == NULL)
    {
        printf("%s %s \n","File not opened ",curFileName);
    }
    else
    {
        printf("%s%s\n","File opened ",curFileName);

        fprintf(fp,"%s%f%s%f%s\n","Удельная электропроводность грунта (интервальная) i_ro = (",i_ro_g[0],";",i_ro_g[1],")");
        fprintf(fp,"%s%f%s%f%s\n","Удельное сопротивление грунта (интервальное) i_sig = (",i_sig[0],";",i_sig[1],")");
        fprintf(fp,"%s%f%s%f%s\n","Удельное сопротивление изоляции трубы (интервальное) i_Ct = (",i_Ct[0],";",i_Ct[1],")");
        fprintf(fp,"%s\n","Исходная матрица системы;");
        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,i_format,"(",i_ABis[i][j][0],";",i_ABis[i][j][1],") ");
            }
        fprintf(fp,"\n");
        }

        fprintf(fp,"%s\n","Решение для исходной системы;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_Xis[i][0],";",i_Xis[i][1],") ");
        }
        fprintf(fp,"\n");

        fprintf(fp,"%s\n","Треугольная матрица системы;");

        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,i_format,"(",i_AB[i][j][0],";",i_AB[i][j][1],") ");
            }
        fprintf(fp,"\n");
        }

        fprintf(fp,"%s\n","Решение для треугольной системы;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_X[i][0],";",i_X[i][1],") ");
        }
        fprintf(fp,"\n");

        fprintf(fp,"%s\n","Массив перенумерации неизв системы;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,"%i%s",i_num[i],"  ");
        }
        fprintf(fp,"\n");

        fprintf(fp,"%s\n","Ширина интервалов треугольной матрицы;");

        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,"%s%24.18f%s"," (",i_AB_iw[i][j],")");
            }
        fprintf(fp,"\n");
        }

        fprintf(fp,"%s\n","Типа Невязка для исходной системы;");

        for (i=0;i<NN;i++)
        {
//            fprintf(fp,i_format,"(",i_nev_is[i][0],";",i_nev_is[i][1],") ");
            fprintf(fp,i_format,"(",i_nev_is[i],";",i_nev_is[i],") ");

        }
        fprintf(fp,"\n");
        fprintf(fp,"%s\n","Правая часть для исходной системы;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_ABis[i][NN][0],";",i_ABis[i][NN][1],") ");
        }
        fprintf(fp,"\n");

        fprintf(fp,"%s\n","Типа Невязка для треугольной системы;");

        for (i=0;i<NN;i++)
        {
//            fprintf(fp,i_format,"(",i_nev_tr[i][0],";",i_nev_tr[i][1],") ");
            fprintf(fp,i_format,"(",i_nev_tr[i],";",i_nev_tr[i],") ");
        }
        fprintf(fp,"%\n");
        fprintf(fp,"%s\n","Правая часть для треугольной системы;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_AB[i][NN][0],";",i_AB[i][NN][1],") ");
        }

        fprintf(fp,"%\n");
        fprintf(fp,"%s\n","Соответствие алгебраическому решению;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_dist_sol[i],";",i_dist_sol[i],") ");
        }


        fprintf(fp,"%\n");


        fprintf(fp,"Ток вдоль трубы Itx_i;\n");

        for (i=0;i<M-1;i++) fprintf(fp,i_format,"(",i_Xis[i][0],";",i_Xis[i][1],") ");
        fprintf(fp,"%\n");

        fprintf(fp,"Ток на границе грунт-труба Itg_i, А;\n");
//    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[M-1+i],";");
        for (i=0;i<M;i++) fprintf(fp,i_format,"(",i_Xis[M-1+i][0],";",i_Xis[M-1+i][1],") ");

        fprintf(fp,"%\n");
/*
    fprintf(fp,"Проверка: Суммарный ток на границе грунт-труба Itg_i;\n");
    Nv = 0;
    for (i=0;i<M;i++) Nv = Nv + X[M-1+i];
    fprintf(fp,outputFormat,Nv);
    fprintf(fp,"%\n");
*/
        fprintf(fp,"Эл. потенциал на границе грунт-туба Utg_i;\n");
//    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[2*M-1+i],";");
        for (i=0;i<M;i++) fprintf(fp,i_format,"(",i_Xis[2*M-1+i][0],";",i_Xis[2*M-1+i][1],") ");
        fprintf(fp,"%\n");


        fprintf(fp,"Эл. потенциал металла трубы Utm_i;\n");
//        for (i=0;i<M;i++) fprintf(fp,outputFormat,X[3*M-1+i],";");
        for (i=0;i<M;i++) fprintf(fp,i_format,"(",i_Xis[3*M-1+i][0],";",i_Xis[3*M-1+i][1],") ");

        fprintf(fp,"%\n");
        fprintf(fp,"Защитный потенциал Utg_i-Utm_i;\n");
//    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[2*M-1+i]-X[3*M-1+i],";");
///НАДО ВЫЧИСЛЯТЬ ИНТЕРВ РАЗНОСТЬ по правилам интерв арифм
/// но мы вычислим по простому
///    Ifloat t;

        for (i=0;i<M;i++)
        {
            i_sub(t,i_Xis[2*M-1+i],i_Xis[3*M-1+i]);
            fprintf(fp,i_format,"(",t[0],";",t[1],") ");
        }
/*
        for (i=0;i<M;i++)
        {
            t[0] = i_Xis[2*M-1+i][0] - i_Xis[3*M-1+i][0];
            t[1] = i_Xis[2*M-1+i][1] - i_Xis[3*M-1+i][1];
            fprintf(fp,i_format,"(",t[0],";",t[1],") ");
        }
*/


        fprintf(fp,"%\n");
        fclose(fp);
        }

  //  char outputFormat1[28] = "%24.16f%s%24.16f%s";
        strcpy(curFileName,"");
        strcat(curFileName,dirName);
        strcat(curFileName,"i_Upr");
        strcat(curFileName,fileType);
        if((fp=fopen(curFileName, "w")) == NULL)
        {
            printf("%s %s \n","Невозможно открыть файл ",curFileName);
        }
        else
        {
            printf("%s %s \n","Возможно открыть файл ",curFileName);
//        fprintf(fp,"Величина защитного потенциала;\n");
            fprintf(fp,"Защитный потенциал Utg_i-Utm_i;\n");
            for (i=0;i<M;i++)
            {
                i_sub(t,i_Xis[2*M-1+i],i_Xis[3*M-1+i]);

                t[0] = i_Xis[2*M-1+i][0] - i_Xis[3*M-1+i][0];
                t[1] = i_Xis[2*M-1+i][1] - i_Xis[3*M-1+i][1];
                fprintf(fp,i_format1,FI_t[i][0],";",t[0],";",t[1]);
            }
            fprintf(fp,"%\n");
            fclose(fp);
        }

        strcpy(curFileName,"");
        strcat(curFileName,dirName);
        strcat(curFileName,"i_Jtg");
        strcat(curFileName,fileType);
        if((fp=fopen(curFileName, "w")) == NULL)
        {
            printf("%s %s \n","Невозможно открыть файл ",curFileName);
        }
        else
        {
            printf("%s%s \n","Возможно открыть файл ",curFileName);
            fprintf(fp,"Плотность тока на границе грунт-труба, мА/м2;\n");

///        i_format1[] = "%16.8f%s24.18f%s%24.18f\n";
            for (i=0;i<M;i++) fprintf(fp,i_format1,FI_t[i][0],";",(1000/St)*i_Xis[M-1+i][0],";",(1000/St)*i_Xis[M-1+i][1]);
//        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",1000*X[M-1+i]/St,";\n");
            fprintf(fp,"%\n");
            fclose(fp);
        }

}


void i_ProvAB()
{
    int i,j;
    Ifloat sum,tmp;
    for (i=0;i<NN;i++)
        {
        sum[0] = 0;
        sum[1] = 0;
        for (j=0;j<NN;j++)
        {
            i_mul(tmp,i_AB[i][j],i_X[j]);
            i_add(sum,sum,tmp);
        }

//        i_movi(i_nev_tr[i],sum);
/// Не будем вычиатать праую часть ,покажем отдельно!!!

        i_nev_tr[i] = i_dist(i_AB[i][NN],sum);
//        i_sub(i_nev_tr[i],i_AB[i][NN],sum);
        }
}

void det_alg_sol()
/// Определение соответствия алгебраическому решению
/// на основе интервальных i_ABis, и i_Xis вычисляем точечное i_dist_sol
{
    int i,j;
    Ifloat t,t1,sum;
    /// З-ны Кирхгофа
    i_dist_sol[0] = i_dist(i_Xis[M-1+0],i_Xis[0]); /// Itg_0 = Itx_0
//i_dist_sol[0] = i_dist_sol[0] +100;
    for (i=1;i<M-1;i++)
    {
        i_add(t,i_Xis[M-1+i],i_Xis[i-1]);    /// Itg_i + Itx_(i-1)
        i_add(t1,i_Xis[i],i_ABis[i][NN]); /// Itx_i + Ia Берем интенсивность анодов Ia из правой части СЛАУ

        i_dist_sol[i] = i_dist(t,t1); /// Itg_i + Itx_(i-1) = Itx_i + Ia
//i_dist_sol[i] = i_dist_sol[i] +100;
    }

    i_movp(t,-1);
    i_mul(t,t,i_Xis[M-2]);
    i_dist_sol[M-1] = i_dist(i_Xis[M-1+M-1],t); /// Itg_M = - Itx_(M-1)
//i_dist_sol[M-1] = i_dist_sol[2*M-2] +100;

//ERROOR!!! HERE!
    /// Граничн услов 3 рода
    for (i=0;i<M;i++)
    {

///        A[i+M][M+i-1]=-FI_t[i][3]/St; Из Init_A

        i_mul(t,i_ABis[i+M][M+i-1],i_Xis[M-1+i]); /// (Ct/St)*Itg_i - коэф в матрице с минусом!
        i_add(t,i_Xis[2*M-1+i],t);                /// Utg_i - (Ct/St)*Itg_i поэтому add!
        i_dist_sol[M+i] = i_dist(i_Xis[3*M-1+i],t); /// Utm_i = Utg_i - (Ct/St)*Itg_i
//i_dist_sol[M+i] = i_dist_sol[M+i] +300;
    }
    /// Законы Ома
    for (i=0;i<M-1;i++)
    {
//        i_movp(t,-ro_t*Li/Sms); /// вОЗЬМЕМ ИЗ МАТРИЦЫ i_Ais
        i_movp(t,-1);
        i_mul(t,t,i_ABis[i+2*M][i]); /// ro_t*Li/Sms вОЗЬМЕМ ИЗ МАТРИЦЫ i_Ais
        i_mul(t,t,i_Xis[i]);
        i_sub(t1,i_Xis[3*M-1+i+1],i_Xis[3*M-1+i]);

/// /////////////// Попробовали вычислить разность по неправильному - ОЧЕНЬ ТОЧНО!!
///        t1[0] = i_Xis[3*M-1+i+1][0] -i_Xis[3*M-1+i][0];
///        t1[1] = i_Xis[3*M-1+i+1][1] -i_Xis[3*M-1+i][1];
/// ///////////////
        i_dist_sol[2*M+i] = i_dist(t,t1); /// -Rt*Itx_i = Utm_(i+1) - Utm_i
//i_dist_sol[2*M+i] = i_dist_sol[2*M+i] +200;
    }

    /// Выражения для ПЭА
    /// сначала сделаем по простому, не перенося в знаменатели сигма
    /// 4Pi*Utg_i = -Sum(Itg/(R*i_sig))  + Sum(Ia/(R*i_sig)) Попробовать и сравнить с переносом sig влево
    /// A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;         /// коэффициент при Utg_i
    /// A[i+3*M-1][4*M-1]=A[i+3*M-1][4*M-1]+FI_a[j][3]*(funR1(Pi,Pj)+funR2(Pi,Pj)); /// коэфф правой части
    /// A[i+3*M-1][M-1+j] = varR1 + funR2(Pi,Pj); /// коэф при Itg_i


    for (i=0;i<M;i++)
    {
        i_movi(sum,i_ABis[i+3*M-1][NN]);  /// Sum(Ia/(R*i_sig))
    cout<<"i_ABis[i+3*M-1][NN][0] = "<<i_ABis[i+3*M-1][NN][0]<<endl;
        for (j=0;j<M;j++)
        {
            i_mul(t,i_Xis[M-1+j],i_ABis[i+3*M-1][M-1+j]);
            i_sub(sum,sum,t); /// -Sum(Itg/(R*i_sig))
        }
    i_mul(t,i_ABis[i+3*M-1][2*M+i-1],i_Xis[2*M-1+i]);     /// 4Pi*Utg_i = -Sum(Itg/(R*i_sig))  + Sum(Ia/(R*i_sig))
    i_dist_sol[3*M-1+i] = i_dist(t,sum);
//i_dist_sol[3*M-1+i] = i_dist_sol[3*M-1+i] +400;
    }
    /// Ток вдоль трубы                     Itx_i for M-1, X[i]
    /// Ток на границе грунт-труба          Itg_i for M, X[M-1+i]
    /// Эл. потенциал на границе грунт-туба Utg_i for M, X[2*M-1+i]
    /// Эл. потенциал металла трубы         Utm_i for M, X[3*M-1+i]
}

void i_ProvABis()
{
    int i,j;
    Ifloat sum,tmp;
    for (i=0;i<NN;i++)
        {
        sum[0] = 0;
        sum[1] = 0;
        for (j=0;j<NN;j++)
        {
            i_mul(tmp,i_ABis[i][j],i_Xis[j]);
            i_add(sum,sum,tmp);
        }
//        i_movi(i_nev_is[i],sum);
/// Не будем вычиатать праую часть ,покажем отдельно!!!

        i_nev_is[i] = i_dist(i_ABis[i][NN],sum);
//        i_sub(i_nev_is[i],i_ABis[i][NN],sum);
        }
}

void no_simply_int()
{
    fesetround (FE_DOWNWARD);

    int i,j;
    Ifloat Itg[M],
           Itx[M],
           Utg[M],
           Utm[M],
           t;

    Ifloat t1[M],t2[M];
    FILE* fp;


printf("%24.18f%\n",Sigma_g);
printf(i_format,"(",i_sig[0],";",i_sig[1],")\n ");

printf("%24.18f%\n",Ct);
printf(i_format,"(",i_Ct[0],";",i_Ct[1],")\n ");

    for (i=0;i<M;i++)
    {
        i_movp(Itg[i],X[M-1+i]);
        i_movp(Itx[i],X[i]);
        i_movp(Utg[i],0);
        i_movp(Utm[i],0);

/// СООТВЕТСТВИЕ ПАРАМЕТРОВ СИСТЕМЫ
///    Ток вдоль трубы Itx_i  X[i]
///    Ток на границе грунт-труба Itg_i X[M-1+i]
///    Эл. потенциал на границе грунт-туба Utg_i X[2*M-1+i]
///    Эл. потенциал металла трубы Utm_i X[3*M-1+i]
    }

int d;
for (d=0;d<1;d++)
{

/// ПРОБУЕМ ПОСЧИТАТЬ ИНТЕРВАЛЬНЫЕ ПАРАМЕТРЫ "ПО НЕ ПРОСТОМУ 1"

/// Вычисляем Utg из ПЭА для точ Itg и инт Sigma_grunta
    for (i=0;i<M;i++)
    {
        i_movp(Utg[i],Ais[i+3*M-1][4*M-1]);
        for (j=0;j<M;j++)
        {
            i_movp(t,Ais[i+3*M-1][M-1+j]);
            i_mul(t,Itg[j],t);
            i_sub(Utg[i],Utg[i],t); /// КОЭФ-ТЫ (1/Rij) ПРИ Itg ИЗ ИСХОДНОЙ М-ЦЫ ИЗ УР-НИЙ ДЛЯ ПЭА
        }
        i_movp(t,4*pi);
        i_mul(t,t,i_sig);
        i_div(Utg[i],Utg[i],t);
    }

/// Вычисляем  инт Utm из гран.услов. 3 рода

    i_movp(t,-1/St);
    i_mul(t,t,i_Ct);

    for (i=0;i<M;i++)
    {
//        Utm[i] = Utg[i] - i_Ct*Itg[i]/St;

        i_movi(Utm[i],Itg[i]);
        i_mul(Utm[i],Utm[i],t);
        i_add(Utm[i],Utm[i],Utg[i]);
    }


/// Вычисляем Itg снова из граничн. условий 3 рода принимая Ct, Utg Utm интервальными

//        Utm[i] = Utg[i] - i_Ct*Itg[i]/St;
///       Itg[i] = St*(Utg[i] - Utm[i])/Ct

    i_movp(t,St);
    for (i=0;i<M;i++)
    {
        i_movi(Itg[i],Utg[i]);
        i_sub(Itg[i],Itg[i],Utm[i]);
        i_mul(Itg[i],Itg[i],t);
        i_div(Itg[i],Itg[i],i_Ct);
    }

/// ПЕРЕНОСИМ ИНТЕРВАЛЬНЫЕ РЕЗУЛЬТАТЫ В ОДИН ВЕКТОР i_Xis

    for (i=0;i<M-1;i++)
    {
        i_movi(i_Xis[i],Itx[i]);
    }
    for (i=0;i<M;i++)
    {
        i_movi(i_Xis[M+i-1],Itg[i]);
        i_movi(i_Xis[2*M+i-1],Utg[i]);
        i_movi(i_Xis[3*M+i-1],Utm[i]);
    }
/// СООТВЕТСТВИЕ ПАРАМЕТРОВ СИСТЕМЫ
///    Ток вдоль трубы Itx_i  X[i]
///    Ток на границе грунт-труба Itg_i X[M-1+i]
///    Эл. потенциал на границе грунт-туба Utg_i X[2*M-1+i]
///    Эл. потенциал металла трубы Utm_i X[3*M-1+i]

/*
    Ifloat sum,tmp;
    for (i=0;i<NN;i++)
        {
        sum[0] = 0;
        sum[1] = 0;
        for (j=0;j<NN;j++)
        {
            i_mul(tmp,i_ABis[i][j],i_Xis[j]);
            i_add(sum,sum,tmp);
        }
        i_sub(i_nev_is[i],i_ABis[i][NN],sum);
        }
   if((fp=fopen("nosimply_Nev.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","ii_Utg.txt");
//        fprintf(fp,"Величина защитного потенциала;\n");
        fprintf(fp,"Невязка no_simply;\n");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_nev_is[i][0],";",i_nev_is[i][1],") \n");
        }

        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fclose(fp);
    }
}

///

   if((fp=fopen("nosimply_Itg.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","ii_Utg.txt");
//        fprintf(fp,"Величина защитного потенциала;\n");
        fprintf(fp,"Ток Itg, интервальное значение;\n");
        for (i=0;i<M-1;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",Itg[i][0],";",Itg[i][1]);
        }
        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fprintf(fp,"Ток Itg, точечное значение;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[i+M-1],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }




   if((fp=fopen("nosimply_Utg.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","ii_Utg.txt");
        fprintf(fp,"Потенциал Utg, интервальное значение;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",Utg[i][0],";",Utg[i][1]);
        }
        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fprintf(fp,"Потенциал Utg, точечное значение;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[2*M-1+i],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

   if((fp=fopen("nosimply_Utm.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","ii_Utg.txt");
//        fprintf(fp,"Величина защитного потенциала;\n");
        fprintf(fp,"Потенциал Utm, интервальное значение;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",Utm[i][0],";",Utm[i][1]);
        }
        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fprintf(fp,"Потенциал Utg, точечное значение;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[3*M-1+i],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

/// Вычислим итервальное Upr

    for (i=0;i<M;i++)
    {
        i_sub(t1[i],Utg[i],Utm[i]);
    }
   if((fp=fopen("nosimply_Upr.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","ii_Utg.txt");
//        fprintf(fp,"Величина защитного потенциала;\n");
        fprintf(fp,"Защитный потенциал Upr, интервальное значение;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",t1[i][0],";",t1[i][1]);
        }
        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fprintf(fp,"Потенциал Upr, точечное значение;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[2*M-1+i]-X[3*M-1+i],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }
*/
}
fesetround (FE_TONEAREST);
}
void simply_int()
{
    fesetround (FE_DOWNWARD);
    int i,j;
    Ifloat Itx[M],
           Itg[M],
           Utg[M],
           Utm[M],
           Upr[M];

    Ifloat t1[M],t2[M],t;
    FILE* fp;

    for (i=0;i<M;i++)
    {
        i_movp(Itg[i],X[M-1+i]);
        i_movp(Itx[i],X[i]);
        i_movp(Utg[i],X[2*M-1+i]);
        i_movp(Utm[i],X[3*M-1+i]);
        i_sub(Upr[i],Utg[i],Utm[i]);
    }

/// ПОСЧИТАТЛИ ИНТЕРВАЛЬНЫЕ ПАРАМЕТРЫ "ПО ПРОСТОМУ"
    i_movp(t,St);
    for (i=0;i<M;i++)
    {
        i_movp(t1[i],St);
        i_mul(t1[i],t1[i],Upr[i]);
        i_div(t1[i],t1[i],i_Ct); /// находим Itg = St*Upr/Ct

        i_movp(t2[i],1/St);
        i_mul(t2[i],t2[i],i_Ct); /// находим Upr = Ct*Itg/St
        i_mul(t2[i],t2[i],Itg[i]);

        i_sub(Utg[i],t2[i],Utm[i]); /// Utg = Upr + Utm
    }
    for (i=0;i<M;i++)
    {
        i_movi(Itg[i],t1[i]);
        i_movi(Itg[i],t1[i]);
    }


/// ПЕРЕНОСИМ ИНТЕРВАЛЬНЫЕ РЕЗУЛЬТАТЫ В ОДИН ВЕКТОР i_Xis

    for (i=0;i<M-1;i++)
    {
        i_movi(i_Xis[i],Itx[i]);
    }
    for (i=0;i<M;i++)
    {
        i_movi(i_Xis[M+i-1],Itg[i]);
        i_movi(i_Xis[2*M+i-1],Utg[i]);
        i_movi(i_Xis[3*M+i-1],Utm[i]);
    }
/// СООТВЕТСТВИЕ ПАРАМЕТРОВ СИСТЕМЫ
///    Ток вдоль трубы Itx_i  X[i]
///    Ток на границе грунт-труба Itg_i X[M-1+i]
///    Эл. потенциал на границе грунт-туба Utg_i X[2*M-1+i]
///    Эл. потенциал металла трубы Utm_i X[3*M-1+i]




/// ПОПЫТКА "ПО ХИТРОМУ" ПОСЧИТАТЬ ИНТЕРВАЛЬНЫЕ ПАРАМЕТРЫ - ИТЕР СХ-СЯ К НУЛЮ
/*    for (i=0;i<M;i++)
    {
//        for (j=0;j<1;j++)
        {
            i_div(dUdn[i],Jtg[i],i_sig);
            i_mul(Upr[i],i_Ct,dUdn[i]);
            i_mul(Upr[i],i_sig,Upr[i]);
            i_div(Jtg[i],Upr[i],i_Ct);
        }
    }
    if((fp=fopen("file1.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","file.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","file.txt");
//        fprintf(fp,"Величина защитного потенциала;\n");
        fprintf(fp,"Защитный потенциал Utg_i-Utm_i;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",Upr[i][0],";",Upr[i][1]);
        }
        fprintf(fp,"%\n");
        fclose(fp);
    }

    if((fp=fopen("file2.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","file1.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","file1.txt");
//        fprintf(fp,"Величина защитного потенциала;\n");
        fprintf(fp,"Плотность тока Jtg;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",1000*Jtg[i][0],";",1000*Jtg[i][1]);
        }
        fprintf(fp,"%\n");
        fclose(fp);
    }
*/
fesetround (FE_TONEAREST);
}

void fox_simply_int()
{
    fesetround (FE_DOWNWARD);
    int i,j;
    Ifloat Itx[M],
           Itg[M],
           Utg[M],
           Utm[M],
           Upr[M],
           dUdn[M];

    Ifloat t1[M],t2[M],t;
    FILE* fp;

    for (i=0;i<M;i++)
    {
        i_movp(Itg[i],X[M-1+i]);
        i_movp(Itx[i],X[i]);
        i_movp(Utg[i],X[2*M-1+i]);
        i_movp(Utm[i],X[3*M-1+i]);
        i_sub(Upr[i],Utg[i],Utm[i]);
    }

/// CЧИТАЕМ dUdn = Itg/(Ct*i_sig_g*St)

    for (i=0;i<M;i++)
    {
        i_movp(dUdn[i],1/St);
        i_div(dUdn[i],dUdn[i],i_sig);
        i_div(dUdn[i],dUdn[i],i_Ct);
        i_mul(dUdn[i],dUdn[i],Itg[i]);
    }
/// НЕ ПОНЯЛ ПОКА КАК ТУТ БЫТЬ
    for (i=0;i<M;i++)
    {
//        i_div(dUdn[i],Jtg[i],i_sig);
        i_mul(Upr[i],i_Ct,dUdn[i]);
        i_mul(Upr[i],i_sig,Upr[i]);
//        i_div(Jtg[i],Upr[i],i_Ct);
    }
/// ПЕРЕНОСИМ ИНТЕРВАЛЬНЫЕ РЕЗУЛЬТАТЫ В ОДИН ВЕКТОР i_Xis

    for (i=0;i<M-1;i++)
    {
        i_movi(i_Xis[i],Itx[i]);
    }
    for (i=0;i<M;i++)
    {
        i_movi(i_Xis[M+i-1],Itg[i]);
        i_movi(i_Xis[2*M+i-1],Utg[i]);
        i_movi(i_Xis[3*M+i-1],Utm[i]);
    }
/// СООТВЕТСТВИЕ ПАРАМЕТРОВ СИСТЕМЫ
///    Ток вдоль трубы Itx_i  X[i]
///    Ток на границе грунт-труба Itg_i X[M-1+i]
///    Эл. потенциал на границе грунт-туба Utg_i X[2*M-1+i]
///    Эл. потенциал металла трубы Utm_i X[3*M-1+i]

fesetround (FE_TONEAREST);
}







/// Интервалы end


MyFloat funR1(MyFloat P_i[3],MyFloat P_j[3])
{
    return      1/sqrt((P_i[0]-P_j[0])*(P_i[0]-P_j[0])+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]-P_j[2])*(P_i[2]-P_j[2]))
                +
                1/sqrt((P_i[0]-P_j[0])*(P_i[0]-P_j[0])+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]+P_j[2])*(P_i[2]+P_j[2]));
}

MyFloat funR2(MyFloat P_i[3],MyFloat P_j[3])
{
    return      1/sqrt((P_i[0]+P_j[0])*(P_i[0]+P_j[0])+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]-P_j[2])*(P_i[2]-P_j[2]))        /// ФИ зерк от-но x=0
                +
                1/sqrt((P_i[0]-(2*L-P_j[0]))*(P_i[0]-(2*L-P_j[0]))+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]-P_j[2])*(P_i[2]-P_j[2]))        /// ФИ зерк от-но x=L
                +
                 1/sqrt((P_i[0]+P_j[0])*(P_i[0]+P_j[0])+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]+P_j[2])*(P_i[2]+P_j[2]))        /// ФИ зерк от-но z=0 x=0
                +
                1/sqrt((P_i[0]-(2*L-P_j[0]))*(P_i[0]-(2*L-P_j[0]))+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]+P_j[2])*(P_i[2]+P_j[2]));       /// ФИ зерк от-но z=0 x=L
}



void Init_FI_t()
{
    int i;
    L = DataSKZ[Dt-1][0]-DataSKZ[0][0] + DataSKZ[1][0]-DataSKZ[0][0]; // Общая длина защищаемого участка трубы
    cout<<"L= "<<L<<endl;
    Li = L / M;                         // Длина участка трубы соотв-го одному ФИ
    St = 2*pi*Rt*Li;                    // Площадь боковой поверхности, соответствующая ФИ
    cout<<"Li= "<<Li<<endl;
    varR1 = 2.0*log((sqrt(Rt*Rt+Li*Li)+Li)/Rt)/Li+        // Выражение для потенциала
            2.0*log((sqrt(4*Zt*Zt+Li*Li)+Li)/(2*Zt))/Li;  // цилиндрического электрода
    cout<<"varR1 = "<<varR1<<endl;
    for (i=0;i<M;i++)
    {
        FI_t[i][0]=i*Li+(Li/2);
        FI_t[i][1]=0;
        FI_t[i][2]=Zt;
        FI_t[i][3]= Ct;
        cout<<"Parametri truby  "<<FI_t[i][0]<<"  "<<FI_t[i][1]<<"  "<<FI_t[i][2]<<"  "<<FI_t[i][3]<<endl;
    }

/// ЗАДАЕМ ПОВРЕЖДЕНИЕ
/// point_povr_Ct = 0.75, /// Относительная коор-та центра участка повреждения изоляции (0,1)
/// koef_povr_Ct = 0.2,   /// Коэффициент поврежденности изоляции (0,1)
///                             0 - полное ровреждение, 1 - без повреждений
/// L_povr_Ct = 100;       /// Длина участка с поврежденной изоляцией, м
/// на участке трубы на расстоянии L*point_povr_Ct от начала защищаемого участка задается участок
/// трубы длиной L_povr_Ct м, состояние изоляции трубы которого таково, что среднее значение
/// удельного электрического сопротивления изолции на данном участке равно Ct*koef_povr_Ct
/// Длина ФИ м/б достаточно велика, и задать повреждение можем лишь для целого числа ФИ,
/// поэтому для обеспечения заданного условия вычислим кол-во ФИ с повреждением и
/// вычислим скорректированный коэфф-т неповрежденности изоляции
/// count_povr_FI = int(L_povr_Ct/Li)+1;  - кол-во ФИ, для которых изменяется Ct.
/// составим пропорцию count_povr_FI*Li*x_Ct*Ct = L_povr_Ct*koef_povr_Ct*Ct,
/// откуда  x_Ct = (L_povr_Ct*koef_povr_Ct) / count_povr_FI*Li
/// Все это делаем в Init_FI_t()

/// int((M-1)*point_povr_Ct - номер ФИ - центра повреждения
/// найдем номер первого ФИ на участке повреждения
/// first_povr_FI = int((M-1)*point_povr_Ct - int(count_povr_FI/2)


/// НЕОБХОДИМО ПРОАНАЛИЗИРОВАТЬ И ОБРАБОТАТЬ СЛУЧАИ, КОГДА
/// first_povr_FI ПОЛУЧИТСЯ ОТРИЦАТЕЛЬНЫМ, ИЛИ (first_povr_FI + count_povr_FI) > M

/// Корректировка коэффициента koef_povr_Ct в соответствии с длиной участка
/// (L_povr_Ct*koef_povr_Ct) / count_povr_FI*Li
/// по целому кол-ву ФИ не продумана, при koef_povr_Ct = 1 мы получается будем
/// корректировать Ст!!!
/// поэтому пока без этих ухищрений

int count_povr_FI = int(L_povr_Ct/Li)+1;
int first_povr_FI = int((M-1)*point_povr_Ct) - int(count_povr_FI/2);
    for (i=0;i < count_povr_FI;i++)
    {
        FI_t[first_povr_FI+i][3]= Ct*koef_povr_Ct;
    }
}

void Init_FI_a()
//void InitKS()
{
    int i,j;
    /// Кол-во точек Dt и СКЗ Dskz посчитано в процедуре LoadDataSKz()
    j = 0;
    for (i=0;i<Dt;i++)
    {
        if (DataSKZ[i][3] != 0) /// Столбец силы тока содержит не ноль (это СКЗ)
        {
//            FI_a[j][0] = DataSKZ[i][0] - DataSKZ[0][0]+0.5*Li;
            FI_a[j][0] = DataSKZ[i][0] - DataSKZ[0][0]+(DataSKZ[1][0] - DataSKZ[0][0])/2;
cout<<"FI_a[j][0] = "<<FI_a[j][0]<<endl;
            FI_a[j][1] = Ya;
            FI_a[j][2] = Za;
            FI_a[j][3] = DataSKZ[i][3];         // сила тока СКЗ

cout<<"i = "<<i<<endl;
cout<<"j = "<<j<<endl;
cout<<"FI_a[j][3] = "<<FI_a[j][3]<<endl;
cout<<"DataSKZ[i][3] = "<<DataSKZ[i][3]<<endl;

            FI_a[j][4] = int(FI_a[j][0]/Li);  // номер ФИ, подключенного к КС
            if (FI_a[j][4] > M-1) i_ks=M-1;
cout<<"FI_a[j][4] = "<<FI_a[j][4]<<endl;
            j++;
        }
    }
}





void Init_A()
{
    int i,j;

    for (i=0;i<NN;i++)
        for (j=0;j<NN+1;j++)
        {
            A[i][j] = 0;
        }
    //Законы Киргофа для ФИ по трубе
    A[0][0]=-1;         /// коэффициент при Itx_0
    A[0][M-1]=1;        /// коэффициент при Itg_0

    for (i=1;i<=M-2;i++)
    {
        A[i][i-1]=1;    /// коэффициент при Itx_(i-1)
        A[i][i]=-1;     /// коэффициент при Itx_i
        A[i][M+i-1]=1;  /// коэффициент при Itg_i
    }
    A[M-1][M-2]=1;      /// коэффициент при Itx_(M-1)
    A[M-1][2*M-2]=1;    /// коэффициент при Itg_M
    for (i=0;i<Dskz;i++)
    {
        cout<<"FI_a[i][3] = "<<FI_a[i][3]<<endl;
        cout<<"int(FI_a[i][4]) = "<<int(FI_a[i][4])<<endl;

        A[int(FI_a[i][4])][NN] = FI_a[i][3];  /// Интенсивности анодов
                                              /// в правую часть системы

        cout<<"A[int(FI_a[i][4])][NN] = "<<A[int(FI_a[i][4])][NN]<<endl;
    }


// Граничные условия 3 рода
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1]=-FI_t[i][3]/St; /// коэффициент при Itg_i Ct/St
        A[i+M][i+2*M-1]=1;            /// коэффициент при Utg_i
        A[i+M][i+3*M-1]=-1;           /// коэффициент при Utm_i
    }

     // Закон Ома между соседними фиктивными источниками
    for (i=0;i<M-1;i++)
    {
        A[i+2*M][i]=ro_t*Li/Sms;    /// коэффициент при Itx_(i) +++!!!
        A[i+2*M][i+3*M-1]=-1;       /// коэффициент при Utx_(i)
        A[i+2*M][i+3*M]=1;          /// коэффициент при Utx_(i+1)
    }


    // Выражения для электростатической аналогии
    for (i=0;i<M;i++)
        {
        A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;   /// коэффициент при Utg_i

        Pi[0] = FI_t[i][0];
        Pi[1] = FI_t[i][1];
        Pi[2] = FI_t[i][2];

        for (j=0;j<Dskz;j++)
        {
                Pj[0] = FI_a[j][0];
                Pj[1] = FI_a[j][1];
                Pj[2] = FI_a[j][2];

                // Вычисляем коэффициент правой части
                A[i+3*M-1][4*M-1]=A[i+3*M-1][4*M-1]+FI_a[j][3]*(funR1(Pi,Pj)+funR2(Pi,Pj));
        }

        for (j=0;j<M;j++)
            {
                Pj[0] = FI_t[j][0];
                Pj[1] = FI_t[j][1];
                Pj[2] = FI_t[j][2];

                /// Вычисляем коэффициент при Itg_i
                if (i==j)
                {
                    A[i+3*M-1][M-1+j] = varR1 + funR2(Pi,Pj);
                }
                else
                {
                    A[i+3*M-1][M-1+j] = funR1(Pi,Pj) + funR2(Pi,Pj);
                }
             }
         }

    for (i=0;i<NN;i++)
    {
        num[i] = i;              // Иходный порядок переменных
        for (j=0;j<NN+1;j++)
        {
            Ais[i][j] = A[i][j]; // Исходная матрица для вычисления невязки
        }
    }
}

//////////// ГАУСС ////////////////////
//////////// ГАУСС ////////////////////


void gl_el(int k)
{
    int i,j,gl;
    MyFloat gll;
//    double gll;
    gl = 0;

    // ищем главный элемент в строке к
    for (i=0;i<NN;i++)
    {
        if (abs(A[k][i])>abs(A[k][gl]))
            {
                gl=i;
            }
    }

    if (gl!=k) // если главный элемент не на диагонали
    {
        for (i=0;i<NN;i++) // меняем столбцы местами
        {
            gll=A[i][k];
            A[i][k] = A[i][gl];
            A[i][gl]=gll;
        }
        //фиксируем изменение порядка записи переменных
        i=num[k];
        num[k]=num[gl];
        num[gl]=i;
    }

    // получаем 1 для ведущего элемента строки к
    gll = A[k][k];
    for (i=0;i<NN+1;i++)
    {
        A[k][i] = A[k][i]/gll;
    }
    // получаем нули в столбце вне ведущей строки
    for (i=0;i<NN;i++)
        if ((i!=k) & (A[i][k]!=0)) // не трогаем сам столбец и те, где уже 0
        {
            gll = A[i][k];
            for (j=0;j<NN+1;j++)
            {
                A[i][j] = A[i][j]/gll;
                A[i][j] = A[i][j] - A[k][j];
            }
        }
}

void mygauss()
{
    int i,j;
    MyFloat gll;

    for (i=0;i<NN;i++)
    {
        gl_el(i);
    }

     for (i=0;i<NN;i++)
     {
         A[i][NN] = A[i][NN]/A[i][i];
         A[i][i] = 1;
     }
     for (i=0;i<NN;i++)
    {
        X[num[i]] = A[i][NN];
    }
     for (i=0;i<NN;i++)
        {
        Nev[i] = - Ais[i][NN];
        for (j=0;j<NN;j++)
            {
                Nev[i] = Nev[i] + Ais[i][j]*X[j];
            }
        }
    Nv = 0;
    for (i=0;i<NN;i++)
    {
       Nv = Nv + abs(Nev[i]);
    }
}




void step_gl_el(int k)
{
    int i,j,gl;
    MyFloat gll;
//    double gll;
    gl = 0;

    // ищем главный элемент в строке к
    for (i=0;i<NN;i++)
    {
        if (abs(A[k][i])>abs(A[k][gl]))
            {
                gl=i;
            }
    }

    if (gl!=k) // если главный элемент не на диагонали
    {
        for (i=0;i<NN;i++) // меняем столбцы местами
        {
            gll=A[i][k];
            A[i][k] = A[i][gl];
            A[i][gl]=gll;
        }
        //фиксируем изменение порядка записи переменных
        i=num[k];
        num[k]=num[gl];
        num[gl]=i;
    }

    // получаем 1 для ведущего элемента строки к
    gll = A[k][k];
    for (i=0;i<NN+1;i++)
    {
        A[k][i] = A[k][i]/gll;
    }
    // получаем нули в столбце вне ведущей строки
    for (i=0;i<NN;i++)
        if ((i!=k) & (A[i][k]!=0)) // не трогаем сам столбец и те, где уже 0
        {
            gll = A[i][k];
            for (j=0;j<NN+1;j++)
            {
                A[i][j] = A[i][j]/gll;
                A[i][j] = A[i][j] - A[k][j];
            }
        }
}

void step_Gauss()
{
    int i,j;
    MyFloat gll;

    for (i=0;i<NN;i++)
    {
        step_gl_el(i);
    }

     for (i=0;i<NN;i++)
     {
         A[i][NN] = A[i][NN]/A[i][i];
         A[i][i] = 1;
     }
     for (i=0;i<NN;i++)
    {
        X[num[i]] = A[i][NN];
    }
}






//void printResult(char* fname,char* comment)
void printResult(char fname[50],char comment[100])
{
    FILE *fp;
    char outputFormat[12] = "%24.18f%s";
    int i,j;
//    strcat(fileName_Res,fname);
//    char filename[100] = __TIME__;
    if((fp=fopen(fname, "w")) == NULL)
    {
        printf("%s %s \n","Not opened file ",fname);
    }
    else
    {
        printf("%s %s \n","Opened file ",fname);

        fprintf(fp,"%s\n",comment);
        fprintf(fp,"%s%f\n","Удельная электропроводность грунта (точечная) ro_g = ",ro_g);
        fprintf(fp,"%s%f\n","Удельное сопротивление грунта (точечное) Sigma_g = ",Sigma_g);
        fprintf(fp,"%s%f%s%f%s\n","Удельное сопротивление изоляции трубы (интервальное) Ct = ",i_Ct[0],";",i_Ct[1],")");
        fprintf(fp,"%s%f\n","Удельное сопротивление изоляции трубы (точечное) Ct = ",Ct);
        fprintf(fp,"%s\n","Исходная матрица системы;");

        fprintf(fp,"Ais - Расширенная матрица системы;\n");
        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,outputFormat,Ais[i][j],";");
            }
            fprintf(fp,"%\n");
        }

        fprintf(fp,"A - Преобразованная матрица системы;\n");
        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,outputFormat,A[i][j],";");
            }
            fprintf(fp,"%\n");
        }
    fprintf(fp,"X - Решение системы;\n");
    for (i=0;i<NN;i++) fprintf(fp,outputFormat,X[i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"Вектор невязки системы;\n");
    for (i=0;i<NN;i++) fprintf(fp,outputFormat,Nev[i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"Норма вектора невязки системы;\n");
    fprintf(fp,outputFormat,Nv);
    fprintf(fp,"%\n");
    fprintf(fp,"Ток вдоль трубы Itx_i;\n");
    for (i=0;i<M-1;i++) fprintf(fp,outputFormat,X[i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"Ток на границе грунт-труба Itg_i;\n");
    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[M-1+i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"Проверка: Суммарный ток на границе грунт-труба Itg_i;\n");
    Nv = 0;
    for (i=0;i<M;i++) Nv = Nv + X[M-1+i];
    fprintf(fp,outputFormat,Nv);
    fprintf(fp,"%\n");
    fprintf(fp,"Эл. потенциал на границе грунт-туба Utg_i;\n");
    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[2*M-1+i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"Эл. потенциал металла трубы Utm_i;\n");
    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[3*M-1+i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"Защитный потенциал Utg_i-Utm_i;\n");
    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[2*M-1+i]-X[3*M-1+i],";");
    fprintf(fp,"%\n");
    fclose(fp);
    }
    char outputFormat1[28] = "%24.16f%s%24.16f%s";
    if((fp=fopen("Upr.txt", "w")) == NULL)
    {
        printf("%s %s \n","Not opened file ","Upr.txt");
    }
    else
    {
        printf("%s %s \n","Opened file ","Upr.txt");
        fprintf(fp,"Защитный потенциал Utg_i-Utm_i;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[2*M-1+i]-X[3*M-1+i],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

    char outputFormat2[28] = "%24.16f%s%24.16f%s";
    if((fp=fopen("Jtg.txt", "w")) == NULL)
    {
        printf("%s %s \n","Not opened file ","Jtg.txt");
    }
    else
    {
        printf("%s %s \n","Opened file ","Jtg.txt");
        fprintf(fp,"Плотность тока на границе грунт-труба;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",1000*X[M-1+i]/St,";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }


        if((fp=fopen("Ct.txt", "w")) == NULL)
    {
        printf("%s %s \n","Not opened file ","Ct.txt");
    }
    else
    {
        printf("%s %s \n","Opened file ","Ct.txt");
        fprintf(fp,"Удельное электрическое сопротивление изоляции трубы;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",FI_t[i][3],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

}

//////////// ГАУСС ////////////////////
//////////// ГАУСС ////////////////////

/// //////// Импорт данных  /////////////////

void LoadDataSKZ()
{
    setlocale(LC_ALL, "rus");
    ifstream fin("Data_SKZ.csv");
    if (!fin.is_open())
        cout << "File Data_SKZ.csv not found!\n";
    else
{
    char buff[500],bf[100],*s[1];
    int i,row,col,lenbuff;
    row = 0;
    col = 0;
    Dskz = 0;
    Dt = 0;
    strcpy(bf,"");
    *s = ";";
    while (fin.getline(buff, 500))
    {
        strncat(buff,s[0],1);                       /// добавляем ";" в конец строки
        cout<<"buff["<<col<<"] = ";
        lenbuff = strlen(buff);                     /// длина строки
        for (i=0;i<lenbuff;i++)
        {
            if (buff[i] != *s[0])                   /// если символ - не разделитель столбцов
            {
                cout<<buff[i];
                strncat(bf,buff+i,1);               /// продолжаем собирать символи очередного числа
            }
            else                                    /// иначе
            {
                DataSKZ[row][col] = atof(bf);       /// преобразуем собранное число к типу double
                                                    /// и заносим в массив
                cout<<"="<<DataSKZ[row][col]<<"; ";
                col++;                              /// переходим к следующему стлбцу
                strcpy(bf,"");                      /// очищаем строку для сбора следующего числа
            }
        }
        row++;
        col = 0;
        cout<<endl;
    }
   for (i=0;i<row;i++) if (DataSKZ[i][3] != 0) Dskz++;     /// считаем кол-во СКЗ
   Dt = row;
//    cout<<"Strok count   "<<Dt<<endl;
//    cout<<"SKZ count   "<<Dskz<<endl;
    fin.close();
//* посмотрели содержимое массива
    for (row=0;row<Dt;row++)
    {
        for (col=0;col<6;col++)
    {
        cout<<DataSKZ[row][col]<<"  ";
    }
    cout<<endl;
    }
//*/
    }
}

/// //////// Импорт данных  /////////////////

void i_InitAB()
/// для случая Ст одинакового для всех ФИ
/// ВВОДИМ ИНТЕРВАЛЬНУЮ МАТРИЦУ i_AB и i_ABis
/// на основе точечной матрицы Ais и интервальных i_sig i_Ct
/// ПОКА ДЕЛАЕМ НА ОСНОВЕ ГОТОВОЙ ТОЧЕЧНОЙ
/// ПРАВИЛЬНЕЕ БУДЕТ ВЫЧИСЛЯТЬ КОЭФ-ТЫ ЗАНОВО
/// С НАПРАВЛЕННЫМИ ОКРУГЛЕНИЯМИ
{
    int i,j;
    Ifloat t;
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            i_movp(i_AB[i][j],Ais[i][j]);
            i_movp(i_ABis[i][j],Ais[i][j]);
        }
        i_num[i] = i;
    }
/// Вводим интервальность по уд сопр грунта в ПЭА
    for (i=0;i<M;i++)
    {
//        A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;
        i_movp(t,Sigma_g);
        i_div(t,i_sig,t);
        i_mul(i_AB[i+3*M-1][2*M+i-1],i_AB[i+3*M-1][2*M+i-1],t);

        i_movi(i_ABis[i+3*M-1][2*M+i-1],i_AB[i+3*M-1][2*M+i-1]); /// копируем результат в матр i_ABis
    }

/// Вводим интервальность по уд сопр изоляции трубы в граничные условия 3 рода

/*
    Фрагмент из проц Init_A
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1]=-FI_t[i][3]/St; /// коэффициент при Itg_i Ct/St
        A[i+M][i+2*M-1]=1;            /// коэффициент при Utg_i
        A[i+M][i+3*M-1]=-1;           /// коэффициент при Utm_i
    }
*/
    for (i=0;i<M;i++)
    {
  /*      i_movp(t,Ct); /// здесь ошибки нет, знак уже "сидит" ЗДЕСЬ БЫЛА ОШИБКА!!!
        i_div(t,i_Ct,t);
        i_mul(i_AB[i+M][M+i-1],i_AB[i+M][M+i-1],t); /// коэффициент при Itg_i Ct/St


*/
        /// заново!
        i_movp(i_AB[i+M][M+i-1],-1/St);
        i_mul(i_AB[i+M][M+i-1],i_AB[i+M][M+i-1],i_Ct);

        i_movi(i_ABis[i+M][M+i-1],i_AB[i+M][M+i-1]);/// копируем результат в матр i_ABis
    }
}


void i_main_Gauss()
{
    int i,j;
    fesetround (FE_DOWNWARD);

    i_InitAB();  /// строим интервальную матрицу
    i_Gauss();
    i_ProvABis();
    i_ProvAB();
    i_printResult();

    fesetround (FE_TONEAREST);


}

/// ////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////// ПРОЕКТ ЛЕВЫХ ПРАВЫХ ИНТЕРВАЛОВ - НАЧАЛО ////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////

void i_InitABl()
/// ВВОДИМ ТОЧЕЧНУЮ МАТРИЦУ ДЛЯ ЛЕВОЙ ГРАНИЦЫ ИНТЕРВАЛА
{
    int i,j;
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            A[i][j] = i_ABis[i][j][0];
            Ais[i][j] = i_ABis[i][j][0];
        }
        num[i] = i;              // Иходный порядок переменных
    }
}

void i_InitABr()
/// ВВОДИМ ТОЧЕЧНУЮ МАТРИЦУ ДЛЯ ПРАВОЙ ГРАНИЦЫ ИНТЕРВАЛА
{
    int i,j;
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            A[i][j] = i_ABis[i][j][1];
            Ais[i][j] = i_ABis[i][j][1];
        }
        num[i] = i;              // Иходный порядок переменных
    }
}


void i_main_LR()
{
//    cur_h_Ct = i_wid(i_Ct_h[i]) // i_Ct_N; - в процедуру перебора интервалов

    int i,j;
//    fesetround (FE_DOWNWARD);

    i_InitAB();  /// строим интервальную матрицу

    i_InitABl(); /// заносим в А левую точечную м-цу
    mygauss();   /// решаем задачу для левой границы интервала

    printResult("resall_l.txt","Результаты точечной левой матрицы"); /// Результаты решения левой системы

/// запомнить результаты для левой матрицы
    for (i=0;i<NN;i++)
    {
        i_Xis[i][0] = X[i];
//        Xtmp[i][0] = X[i];
    }


    i_InitABr(); /// заносим в АВ правую точечную м-цу
    mygauss();   /// решаем задачу для левой границы интервала

   printResult("resall_r.txt","Результаты точечной правой системы"); /// Результаты решения правой системы


/// запомнить результаты для левой матрицы
/// собрать интервальное решение
    for (i=0;i<NN;i++)
    {
        if (X[i] < i_Xis[i][0])
        {
            i_Xis[i][1] = i_Xis[i][0];
            i_Xis[i][0] = X[i];
        }
        else
        {
            i_Xis[i][1] = X[i];
        }
    }

/// посчитать невязку
    i_ProvABis();
    det_alg_sol();
    i_printResult();

 //   fesetround (FE_TONEAREST);
}
/// ////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////// ПРОЕКТ ЛЕВЫХ ПРАВЫХ ИНТЕРВАЛОВ - КОНЕЦ ////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////



/// ////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////// ПРОЕКТ ПЕРЕБОРА ВСЕХ ИНТЕРВАЛОВ - НАЧАЛО ////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////



const unsigned int i_Ct_N = 1,     /// число шагов разбиения интервала i_Ct_FI[i]
                   i_sig_N = 1;   /// число шагов разбиения интервала i_sig

const unsigned long i_Ct_cmb_cnt = 1024*64; /// Число комбинаций по i_Ct (i_Ct_N+1+1)^M 9765625  390625 3125
const int  i_sig_cmb_cnt = i_sig_N+1; /// Число комбинаций по i_sig (i_sig_N+1)


unsigned long step;    ///    Текущий шаг перебора комбинаций - будет содержать число комбинаций для которых решеныточечные задачи

MyFloat cur_Ct[M],  /// для хранения текущих точечных Ст при пробегании интервалов
        cur_sig,    /// для хранения текущего точечного sigma_g
        h_Ct[M],       /// шаг перебора
        h_sig,      /// шаг перебора
        all_X[4][NN+2];  /// Все решения для всех комбинаций - только для проекта lim
//        all_Ct[i_Ct_cmb_cnt][M];  /// Все комбинации Ct
//        all_X[i_Ct_N+2][NN+2];

Ifloat i_Ext_Sol[NN+M]; /// Вектор для получения внешней оценки решения all_X[] + защитного потенциала.
float eps_Ct = 1,       /// Величина, малая отностиельно Ct
      eps_sig = 0.0001; /// Величина, малая отностиельно sigma_g
                        /// необходимы для точного "попадания" в границу интервала в циклах
    FILE* fn_step;


void i_Init_Ct()
/// Получаем интервальные Ct для каждого ФИ, изменения в отдельных ФИ можно вносить здесь или в Init_FI_t()
{
    int i;
//    i_Ct = {Ct*(100-disp_Ct)/100,Ct*(100+disp_Ct)/100}, /// одно значение Ст для всех ФИ, исп в первых вериях
//    i_Ct_FI[i][0] = FI_t[3]+disp
    for (i=0;i<M;i++)
    {
        i_Ct_FI[i][0] = FI_t[i][3]*(100-disp_Ct)/100;
        i_Ct_FI[i][1] = FI_t[i][3]*(100+disp_Ct)/100;
        h_Ct[i] = i_wid(i_Ct_FI[i])/i_Ct_N;
// normal        cout<<i_Ct_FI[i][0]<<"; "<<i_Ct_FI[i][1]<<endl;
//        i_movi(i_Ct_FI[i],{FI_t[i][3]*(100-disp_Ct)/100,FI_t[i][3]*(100+disp_Ct)/100});
    }
}
void i_Init_ABis()
/// ВВОДИМ ИНТЕРВАЛЬНУЮ МАТРИЦУ i_AB И i_ABis для проверок
/// на основе точечной матрицы Ais и интервальных i_sig i_Ct_FI[M]
{
    int i,j;
    Ifloat t;
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            i_movp(i_AB[i][j],Ais[i][j]);
//            i_movp(i_AB[i][j],i_ABis[i][j]);
        }
        i_num[i] = i;
    }
/// Вводим интервальность по уд сопр грунта в ПЭА
    for (i=0;i<M;i++)
    {
//        A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;
        i_movp(t,4*pi);
        i_mul(i_AB[i+3*M-1][2*M+i-1],t,i_sig);

    }

/// Вводим интервальность по уд сопр изоляции трубы в граничные условия 3 рода

/*
    Фрагмент из проц Init_A
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1]=-FI_t[i][3]/St; /// коэффициент при Itg_i Ct/St
        A[i+M][i+2*M-1]=1;            /// коэффициент при Utg_i
        A[i+M][i+3*M-1]=-1;           /// коэффициент при Utm_i
    }
*/
    for (i=0;i<M;i++)
    {
        i_movp(i_AB[i+M][M+i-1],-1/St);
        i_mul(i_AB[i+M][M+i-1],i_AB[i+M][M+i-1],i_Ct_FI[i]); /// коэффициент при Itg_i -Ct/St
    }

/// Копируем i_AB[i][j] в i_ABis[i][j]
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            i_movi(i_ABis[i][j],i_AB[i][j]);
        }
    }
}

void i_InitStepAB(MyFloat cur_Ct[M],MyFloat cur_sig)
/// cur_Ct[i] содержит точечное значение Ст для i-го ФИ
/// ВВОДИМ ТОЧЕЧНУЮ МАТРИЦУ A
/// ДЛЯ ЗАДАННЫХ ТОЧЕЧНЫХ Ct[i] и sigma_g
/// надо бы и i_Ais изменить для выч-я точности!!! в отдельной процедуре
{
    unsigned int i,j;
    Ifloat t;
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            A[i][j] = Ais[i][j]; /// берем коэф-ты из начальной, "средней" точечной м-цы
        }
        num[i] = i;
    }

    for (i=0;i<M;i++)
    {
        A[i+3*M-1][2*M+i-1]=4*pi*cur_sig; /// уд сопр грунта в ПЭА A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;
        A[i+M][M+i-1] = -cur_Ct[i]/St; /// коэффициент при Itg_i Ct/St  гран усл 3 рода
    }

/// Вносим текущее точечное уд сопр изоляции трубы в граничные условия 3 рода
/*
    Фрагмент из проц Init_A
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1]=-FI_t[i][3]/St; /// коэффициент при Itg_i Ct/St
        A[i+M][i+2*M-1]=1;            /// коэффициент при Utg_i
        A[i+M][i+3*M-1]=-1;           /// коэффициент при Utm_i
    }
*/
/*
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1] = -cur_Ct[i]/St; /// коэффициент при Itg_i Ct/St  гран усл 3 рода
    }

    for (i=0;i<NN;i++)
    {
        num[i] = i;
        for (j=0;j<NN+1;j++)
        {
            Ais[i][j] = A[i][j]; ///
        }

    }
*/


/// normal    MatrPrn_NxNp1(A,NN,NN+1,"stA.txt");
}

void print_cur_step()
{
    int i,j;
    fprintf(fn_step,"%8.5f%s",cur_sig,"; ");
    for (i=0;i<M;i++)
    {
        fprintf(fn_step,"%12.9f%s",cur_Ct[i],"; ");
    }
//    fprintf(fn_step,"\n");
//    fprintf(fn_step,"%8.5f%s",cur_sig,"; ");
    for (i=0;i<M;i++)
    {
        fprintf(fn_step,"%12.9f%s",X[M-1+i],"; ");
    }
    fprintf(fn_step,"\n");
}

void Iter(int p)
/// эта процедура только формирует Ct, решает точечные СЛАУ и ищет макс и мин решений покомпонентно
/// НУЖНО ВИДЕТЬ ВНЕШНИЕ ПЕРЕМЕННЫЕ
/// all_X[NN]
/// cur_Ct[M];
/// cur_sig
/// h_Ct
/// step
{
unsigned int w;
    if (p < M-1)
    {
        Iter(p+1);
        while (cur_Ct[p+1] < i_Ct_FI[p+1][1] - 2*h_Ct[p+1] + eps_Ct)
        {
            cur_Ct[p+1] = cur_Ct[p+1] + h_Ct[p+1]; /// При разных Ст будет ли правильно, подумать!
            Iter(p+1);

        }

        cur_Ct[p+1] = i_Ct_FI[p+1][1]; /// Верхнюю границу не вычисляем а берем из исходного интервала
                                       /// для обеспечения  большей точночти и совпадения с lim
        Iter(p+1);

//        printf("%s%f\n","cur_Ct[p+1] = ",cur_Ct[p+1]-h_Ct);
        cur_Ct[p+1] = i_Ct_FI[p+1][0];
    }
    else
    {
/*
        printf("%s%f","cur_sig = ",cur_sig);
        for (w=0;w<M;w++)
        {

            printf("%s%f%s%f",",   cur_Ct[w] = ",cur_Ct[w]);
        }
        printf("\n");
*/
        /// Формируем точечную матрицу для фикс точечных cur_Ct и cur_sig
        i_InitStepAB(cur_Ct,cur_sig);

        /// Решаем систему с текущими точечными cur_Ct и cur_sig
        step_Gauss(); /// Немного оптимизированный, также не считает невязку для ускорения
//        mygauss();


        /// Анализируем полученное точечное решение на макс и мин
        /// i_Ext_Sol будет содержать покомпонентные макс и мин получаемых решений X
        for (w=0;w<NN;w++)
        {
            if (X[w] < i_Ext_Sol[w][0])
            {
                i_Ext_Sol[w][0] = X[w];
/// Была мысль запомнить сигма и Ст при которых будет макс и мин - они разные для разных компонент!!!
///                i_Ext_sig = cur_sig;/// Запомним Сигма
/// В общем пока не понял возможно ли это! 17.40 17.08
                /// Запомним комбинацию Сигма

///                i_Ext_Ct[] = cur_Ct[];

            }
            if (X[w] > i_Ext_Sol[w][1])
            {
                i_Ext_Sol[w][1] = X[w];
            }
        }
        for (w=0;w<M;w++)
        {
            if (X[2*M-1+w] - X[3*M-1+w] < i_Ext_Sol[NN+w][0])
            {
                i_Ext_Sol[NN+w][0] = X[2*M-1+w] - X[3*M-1+w];
            }

            if (X[2*M-1+w] - X[3*M-1+w] > i_Ext_Sol[NN+w][1])
            {
                i_Ext_Sol[NN+w][1] = X[2*M-1+w] - X[3*M-1+w];
            }
        }

print_cur_step(); ///  Выводим результаты текущего расчета в файл.
        step++;
        printf("%u%s%u\n",step," из ",i_sig_cmb_cnt*i_Ct_cmb_cnt);


    }
}

void i_print_Comb()
{
    unsigned int i,j,k;
    /// Вывести построчн в отдельный файл результат
    FILE *fp;
  /*  if((fp=fopen("all_X.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file");
    }
    else
    {
        printf("%s \n","Opened file ");
        fprintf(fp,"Решения all_X \n");

        for (i=0;i<step*(i_sig_cmb_cnt);i++)
        {
            for (j=0;j<NN;j++)
            {
                fprintf(fp,"%24.18f%s",all_X[i][j],";  ");
            }
            fprintf(fp,"\n");
        }
    }
    fclose(fp);

    if((fp=fopen("all_Ct_Upr.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file");
    }
    else
    {
        printf("%s \n","Opened file ");
        fprintf(fp,"комбинации sigma, Ct и соответствующий защитный потенциал;\n");

        cur_sig = i_sig[1];    /// ВОЗМОЖНО НУЖНО ПЕРЕБИРАТЬ В ОБРАТНОМ НАПРАВЛЕНИИ
                      ///  ЧТОБЫ МИН И МАКС БЫЛИ В НАЧАЛЕ И В КОНЦЕ СПИСКА КОМПБИНАЦИЙ
        for (k=0;k<i_sig_N+1;k++) /// цикл перебора точечных sigma_g
        {
            for (j=0;j<step;j++) /// цикл перебора комбинаций Ct по ФИ
            {
                fprintf(fp,"%24.18f%s",cur_sig,";  ");
                for (i=0;i<M;i++) /// Цикл перебора Ct ФИ
                {
                    fprintf(fp,"%24.18f%s",all_Ct[j][i],";  ");
                }
                fprintf(fp,"\n");
                fprintf(fp,"%24.18f%s",0,";  ");
                for (i=0;i<M;i++) /// Цикл перебора Ct ФИ
                {
                    fprintf(fp,"%24.18f%s",all_X[j+k*step][2*M-1+i]-all_X[j+k*step][3*M-1+i],";  ");
                }
                fprintf(fp,"\n");
            }
            cur_sig = cur_sig + h_sig;
        }
    }
    fclose(fp);
*/
    if((fp=fopen("Ext_Sol.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file");
    }
    else
    {
        printf("%s \n","Opened file ");
        fprintf(fp,"Внешняя оценка решения Ext_Sol[NN+M] \n");

        for (i=0;i<NN+M;i++) ///
        {
            fprintf(fp,"%24.18f%s",i_Ext_Sol[i][0],";  ");
        }
        fprintf(fp,"\n");
        for (i=0;i<NN+M;i++) ///
        {
            fprintf(fp,"%24.18f%s",i_Ext_Sol[i][1],";  ");
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}
/*
void i_Create_Sol()
/// Использует all_X и строит брус - интервальную оболочку приближенного решения
/// методом нахожденя макс и мин значений для каждого столбца (ФИ)
/// результат заносит в Ifloat i_Ext_Sol[M] - внешняя оценка решения
/// Подумать про защитный потенциал, его в all_X нет!!!
{
    unsigned int i,j;



    for (i=0;i<NN;i++) /// цикл перебора компонент вектора all_X - неизвестных системы
    {
        i_Ext_Sol[i][0] = all_X[0][i];
        i_Ext_Sol[i][1] = all_X[0][i];
        for (j=1;j<(i_sig_cmb_cnt)*step;j++) /// цикл перебора решений (строк массива all_X)
        {
            if (all_X[j][i] < i_Ext_Sol[i][0])
            {
                i_Ext_Sol[i][0] = all_X[j][i];
            }
            if (all_X[j][i] > i_Ext_Sol[i][1])
            {
                i_Ext_Sol[i][1] = all_X[j][i];
            }
        }
    }


    for (i=0;i<M;i++) /// цикл перебора компонент вектора all_X - неизвестных системы
    {
        i_Ext_Sol[NN+i][0] = all_X[0][2*M-1+NN+i] - all_X[0][3*M-1+NN+i];
        i_Ext_Sol[NN+i][1] = i_Ext_Sol[NN+i][0];

        for (j=1;j<(i_sig_cmb_cnt)*step;j++) /// цикл перебора решений (строк массива all_X)
        {

            if (all_X[j][2*M-1+i] - all_X[j][3*M-1+i] < i_Ext_Sol[NN+i][0])
            {
                i_Ext_Sol[NN+i][0] = all_X[j][2*M-1+i] - all_X[j][3*M-1+i];
            }

            if (all_X[j][2*M-1+i] - all_X[j][3*M-1+i] > i_Ext_Sol[NN+i][1])
            {
                i_Ext_Sol[NN+i][1] = all_X[j][2*M-1+i] - all_X[j][3*M-1+i];
            }
        }
    }

    for (i=0;i<NN+M;i++) ///
    {
        cout<<i_Ext_Sol[i][0]<<";  ";
    }
    cout<<endl;
    for (i=0;i<NN+M;i++) ///
    {
        cout<<i_Ext_Sol[i][1]<<";  ";
    }
    cout<<endl;

    for (i=0;i<M;i++) ///
    {
        cout<<i_Ext_Sol[NN+i][0]<<";  ";
    }
    cout<<endl;
    for (i=0;i<M;i++) ///
    {
        cout<<i_Ext_Sol[NN+i][1]<<";  ";
    }
    cout<<endl;



}
*/

void i_main_step()
/// СУЩЕСТВЕННОЕ ЗАМЕЧАНИЕ
/// МЫ ИМЕЕМ ФАКТИЧЕСКИ М ИНТЕРВАЛЬНЫХ ПАРАМЕТРОВ Ст,
/// И ПО КАЖДОМУ ИЗ НИХ НУЖНО ПРОБЕГАТЬ ИНТЕРВАЛ ДЛЯ ВРЕРБОРА ВСЕХ
/// ВОЗМОЖНЫХ ВАРИАНТОВ - СОЧЕТАНИЙ. ЭТО ВОЗМОЖЕО СДЕЛАТЬ С ПОМОЩЬЮ РЕКУРСИИ?

/// ПРИБЛИЖЕННОЕ ВНЕШНЕЕ ОЦЕНИВАНИЕ ПУТЕМ РЕШЕНИЯ
/// ТОЧЕЧНЫХ ЗАДАЧ С ШАГОМ ПО ИНТЕРВАЛАМ

/// Использует
/// i_Ct_FI[M] - интервальные уд сопр изол для каждого ФИ
/// i_sig - интервальное уд сопр грунта
/// i_Ct_N кол-во шагов по интервалу Ct
/// i_sig_N кол-во шагов по интервалу sigma
{
    unsigned int i,j,k;

    i_Init_Ct();  /// в массив i_Ct_FI заносим интеррвальные Ct для каждого ФИ
                  /// на отснове i_Ct_FI и i_sig

    i_Init_ABis(); /// строим интервальную матрицу i_ABis для проверок решения
/// normal i_printResult();
    /// Будем из i_Ct_FI ,брать с шагом h_Ct точечные вектора cur_Ct,
    /// h_Ct свой для каждого i h_Ct = (i_wid(i_Ct_FI[i]))/i_Ct_N

    step = 0;
    for (i=0;i<M;i++)
    {
        cur_Ct[i] = i_Ct_FI[i][0]; /// начинаем с левой границы Ct
    }

    for (i=0;i<NN+M;i++) ///
    {
        i_Ext_Sol[i][0] = 100;  /// задаем начальные заведомо большие и малые значения
        i_Ext_Sol[i][1] = -100; /// для последующего поиска
    }


    cur_sig = i_sig[0];    /// ВОЗМОЖНО НУЖНО ПЕРЕБИРАТЬ В ОБРАТНОМ НАПРАВЛЕНИИ
                           ///  ЧТОБЫ МИН И МАКС БЫЛИ В НАЧАЛЕ И В КОНЦЕ СПИСКА КОМПБИНАЦИЙ НЕТ!!! ВСЕ НОРМАЛЬНО!!!
    h_sig = i_wid(i_sig) / i_sig_N;

/// Открываем файл для записи fn_step

    if((fn_step=fopen("i_Itg_step.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file i_Itg_step.txt");
    }
    else
    {
        printf("%s \n","Opened file i_Itg_step.txt");
        fprintf(fn_step,"Комбинации для тока Itg грунт-труба \n");
    }

    for (k=0;k<i_sig_N;k++) /// цикл перебора точечных sigma_g
    {
        Iter(-1);
//        printf("%f\n",cur_sig);
        cur_sig = cur_sig + h_sig;
    }

    cur_sig = i_sig[1]; /// Аналогично Ct для большей точности и совпадения с lim правую границу не вычисляем а берем из i_sig
    Iter(-1);



    printf("%s%f%s%f \n","i_sig = ",i_sig[0],"  ",i_sig[1]);

    fclose(fn_step);
    i_print_Comb(); /// нужно передавать имя файла?


}




void i_Create_Sol_lim()
/// Использует all_X и строит брус - интервальную оболочку приближенного решения
/// методом нахожденя макс и мин значений для каждого столбца (ФИ)
/// результат заносит в Ifloat i_Ext_Sol[M] - внешняя оценка решения

{
    unsigned int i,j;



    for (i=0;i<NN;i++) /// цикл перебора компонент вектора all_X - неизвестных системы
    {
        i_Ext_Sol[i][0] = all_X[0][i];
        i_Ext_Sol[i][1] = all_X[0][i];

        for (j=1;j<4;j++) /// цикл перебора решений (строк массива all_X)
        {
            if (all_X[j][i] < i_Ext_Sol[i][0])
            {
                i_Ext_Sol[i][0] = all_X[j][i];
            }
            if (all_X[j][i] > i_Ext_Sol[i][1])
            {
                i_Ext_Sol[i][1] = all_X[j][i];
            }
        }
    }

    for (i=0;i<M;i++) /// цикл перебора компонент вектора all_X - неизвестных системы
    {
//        i_Ext_Sol[NN+i][0] = all_X[0][2*M-1+NN+i] - all_X[0][3*M-1+NN+i];
        i_Ext_Sol[NN+i][0] = all_X[0][2*M-1+i] - all_X[0][3*M-1+i];
        i_Ext_Sol[NN+i][1] = i_Ext_Sol[NN+i][0];

        for (j=1;j<4;j++) /// цикл перебора решений (строк массива all_X)
        {

            if (all_X[j][2*M-1+i] - all_X[j][3*M-1+i] < i_Ext_Sol[NN+i][0])
            {
                i_Ext_Sol[NN+i][0] = all_X[j][2*M-1+i] - all_X[j][3*M-1+i];
            }

            if (all_X[j][2*M-1+i] - all_X[j][3*M-1+i] > i_Ext_Sol[NN+i][1])
            {
                i_Ext_Sol[NN+i][1] = all_X[j][2*M-1+i] - all_X[j][3*M-1+i];
            }
        }
    }
}

void i_print_Comb_lim()
{
    unsigned int i,j,k,s;
    /// Вывести построчн в отдельный файл результат
    FILE *fp;
    if((fp=fopen("all_X_lim.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file all_X_lim.txt");
    }
    else
    {
        printf("%s \n","Opened file all_X_lim.txt");
        fprintf(fp,"Решения all_X \n");

        for (i=0;i<step;i++)
        {
            for (j=0;j<NN;j++)
            {
                fprintf(fp,"%24.18f%s",all_X[i][j],";  ");
            }
            fprintf(fp,"\n");
        }
    }
    fclose(fp);

    if((fp=fopen("all_Ct_Upr_lim.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file all_Ct_Upr_lim.txt");
    }
    else
    {
        printf("%s \n","Opened file all_Ct_Upr_lim.txt");
        fprintf(fp,"комбинации sigma, Ct и соответствующий защитный потенциал;\n");
        s = 0;
        for (j=0;j<2;j++)
        {
            for (i=0;i<M;i++)
            {
                cur_Ct[i] = i_Ct_FI[i][j]; /// начинаем с левой границы Ct
            }
            for (k=0;k<2;k++)
            {
                cur_sig = i_sig[k];

            /// Делаем вывод в файл

/// ////////////////////////////////
                fprintf(fp,"%24.18f%s",cur_sig,";  ");
                for (i=0;i<M;i++) /// Цикл перебора Ct ФИ
                {
                    fprintf(fp,"%24.18f%s",cur_Ct[i],";  ");
                }
                fprintf(fp,"\n");
                fprintf(fp,"%24.18f%s",0,";  ");
                for (i=0;i<M;i++) /// Цикл перебора Ct ФИ
                {
///                    !!!fprintf(fp,"%24.18f%s",all_X[(k+1)*(j+1)-1][2*M-1+i]-all_X[(k+1)*(j+1)-1][3*M-1+i],";  ");
                    fprintf(fp,"%24.18f%s",all_X[s][2*M-1+i]-all_X[s][3*M-1+i],";  ");
                }
                s++;
                fprintf(fp,"\n");

/// /////////////////////////////////

            }
        }
        fclose(fp);
    }

    if((fp=fopen("Ext_Sol_lim.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file Ext_Sol_lim.txt");
    }
    else
    {
        printf("%s \n","Opened file Ext_Sol_lim.txt");
        fprintf(fp,"Внешняя оценка решения Ext_Sol[NN+M] \n");

        for (i=0;i<NN+M;i++) ///
        {
            fprintf(fp,"%24.18f%s",i_Ext_Sol[i][0],";  ");
        }
        fprintf(fp,"\n");
        for (i=0;i<NN+M;i++) ///
        {
            fprintf(fp,"%24.18f%s",i_Ext_Sol[i][1],";  ");
        }
        fprintf(fp,"\n");
        fclose(fp);
    }

///        i_format1[] = "%16.8f%s24.18f%s%24.18f\n";

/// Вывод Upr в формате для графера
        if((fp=fopen("i_Upr_lim.txt", "w")) == NULL)
        {
            printf("%s %s \n","Not opened file ","i_Upr_lim.txt");
        }
        else
        {
            printf("%s %s \n","Opened file ","i_Upr_lim.txt");
//        fprintf(fp,"Величина защитного потенциала;\n");
            fprintf(fp,"Защитный потенциал Utg_i-Utm_i;\n");
            for (i=0;i<M;i++)
            {
                fprintf(fp,i_format1,FI_t[i][0],";",i_Ext_Sol[NN+i][0],";",i_Ext_Sol[NN+i][1]);
            }
            fprintf(fp,"%\n");
            fclose(fp);
        }

/// Вывод Jtg в формате для графера
        if((fp=fopen("i_Jtg_lim.txt", "w")) == NULL)
        {
            printf("%s %s \n","Not opened file ","i_Jtg_lim.txt");
        }
        else
        {
            printf("%s%s \n","Opened file ","i_Jtg_lim.txt");
            fprintf(fp,"Плотность тока на границе грунт-труба, мА/м2;\n");


            for (i=0;i<M;i++) fprintf(fp,i_format1,FI_t[i][0],";",(1000/St)*i_Ext_Sol[M-1+i][0],";",(1000/St)*i_Ext_Sol[M-1+i][1]);
//        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",1000*X[M-1+i]/St,";\n");
            fprintf(fp,"%\n");
            fclose(fp);
        }

/// Вывод Ct в формате для графера
        if((fp=fopen("i_Ct_lim.txt", "w")) == NULL)
        {
            printf("%s %s \n","Not opened file ","i_Ct_lim.txt");
        }
        else
        {
            printf("%s%s \n","Opened file ","i_Ct_lim.txt");
            fprintf(fp,"Удельное сопротивление изоляции трубы;\n");


            for (i=0;i<M;i++) fprintf(fp,i_format1,FI_t[i][0],";",i_Ct_FI[i][0],";",i_Ct_FI[i][1]);
//        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",1000*X[M-1+i]/St,";\n");
            fprintf(fp,"%\n");
            fclose(fp);
        }





/*
    char outputFormat1[28] = "%24.16f%s%24.16f%s";
    if((fp=fopen("Upr_lim.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","Upr.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","Upr.txt");
        fprintf(fp,"Защитный потенциал Utg_i-Utm_i;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[2*M-1+i]-X[3*M-1+i],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

    char outputFormat2[28] = "%24.16f%s%24.16f%s";
    if((fp=fopen("Jtg.txt", "w")) == NULL)
    {
        printf("%s %s \n","Невозможно открыть файл ","Jtg.txt");
    }
    else
    {
        printf("%s %s \n","Возможно открыть файл ","Jtg.txt");
        fprintf(fp,"Плотность тока на границе грунт-труба;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",1000*X[M-1+i]/St,";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }
*/


}

void i_Itg_comb()
/// ПРОЦЕДУРА ЗАНОСИТ РЕЗУЛЬТАТЫ В i_Ext_Sol, ПОЭТОМУ ЗАПУСКАТЬ ПОСЛЕ НАХОЖДЕНИЯ Upr
/// Будем решать системы, перебирая комбинации
/// вида (Ct_0, Ct_1,...,Ct_i-1,Ct^i,Ct_i+1,...Ct_M-1)
/// и (Ct^0, Ct^1,...,Ct^i-1,Ct_i,Ct^i+1,...Ct^M-1)
/// Всего будет 4*М комбинаций
{
    int i,j,k,w;
    FILE* fp;
    step = 0;
    if((fp=fopen("i_Ig_comb.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file i_Ig_comb.txt");
    }
    else
    {
        printf("%s \n","Opened file i_Ig_comb.txt");
    }
    fprintf(fp,"Все комбинации Itg для Ct вида (Ct_0, Ct_1,...,Ct_i-1,Ct^i,Ct_i+1,...Ct_M-1) и (Ct^0, Ct^1,...,Ct^i-1,Ct_i,Ct^i+1,...Ct^M-1) \n");

    for (i=0;i<M;i++) /// Перебор ФИ для которого будем искать мин и макс Itg
    {

        i_Ext_Sol[M-1+i][0] = 100; /// Можно оставить мин и макс найденные вместе с Upr?
        i_Ext_Sol[M-1+i][1] = -100;


        for (w=0;w<2;w++) /// Цикл перебора гранмчных значений Ct,i
        {
            /// Формируем вектор Ст требуемого вида
            for (j=0;j<M;j++)
            {
                cur_Ct[j] = i_Ct_FI[j][w];
            }
//            cur_Ct[i] = i_Ct_FI[i][((w+1) mod 2)]; /// если w =0 то берем 1 и наоборот
            cur_Ct[i] = i_Ct_FI[i][1-w]; /// если w =0 то берем 1 и наоборот

            /// Цикл перебора граничных значений сигма
            for (k=0;k<2;k++)
            {
                cur_sig = i_sig[k];

                /// Формируем точечную систему и решаем

                i_InitStepAB(cur_Ct,cur_sig);
                step_Gauss();

                /// Необходимо найти макс и мин по каждому ФИ  и объединить с общими рез-ми
                /// Будем заносить в i_Ext_Sol

                if (X[M-1+i] < i_Ext_Sol[M-1+i][0])
                {
                    i_Ext_Sol[M-1+i][0] = X[M-1+i];
                }
                if (X[M-1+i] > i_Ext_Sol[M-1+i][1])
                {
                    i_Ext_Sol[M-1+i][1] = X[M-1+i];
                }

                /// сохраним результат
                for (j=0;j<M;j++)
                {
                    fprintf(fp,"%8.6f%s",X[M-1+j],"; ");
                }

                /// Сохраним комбинацию Ст
                for (j=0;j<M;j++)
                {
                    fprintf(fp,"%7.1f%s",cur_Ct[j],"; ");
                }

                ///  и значение сигма
                fprintf(fp,"%7.5f%s",cur_sig,"; ");
                fprintf(fp,"\n");
                step++;
                printf("%i%s%i\n",step," из ",4*M);
            }

        }

    }
    fclose(fp);

}



void i_main_lim()
/// ПРИБЛИЖЕННОЕ ВНЕШНЕЕ ОЦЕНИВАНИЕ ПУТЕМ РЕШЕНИЯ
/// КОМБИНАЦИЙ ТОЧЕЧНЫХ ЗАДАЧ ДЛЯ ГРАНИЧНЫХ ЗНАЧЕНИЙ ИНТЕРВАЛОВ

/// Использует
/// i_Ct_FI[M] - интервальные уд сопр изол для каждого ФИ
/// i_sig - интервальное уд сопр грунта
/// i_Ct_N кол-во шагов по интервалу Ct
/// i_sig_N кол-во шагов по интервалу sigma
{
    unsigned int i,j,k;

    i_Init_Ct();  /// в массив i_Ct_FI заносим интервальные Ct для каждого ФИ
                  /// на основе i_Ct_FI и i_sig

    i_Init_ABis(); /// строим интервальную матрицу i_ABis для проверок решения
/// normal i_printResult();
    step = 0;
    for (j=0;j<2;j++)
    {
        for (i=0;i<M;i++)
        {
            cur_Ct[i] = i_Ct_FI[i][j]; /// начинаем с левой границы Ct
        }
        for (k=0;k<2;k++)
        {
            cur_sig = i_sig[k];

            /// РЕШАЕМ ТОЧЕЧНУЮ ЗАДАЧУ!!!
            i_InitStepAB(cur_Ct,cur_sig);
            mygauss();

            /// Записываем текущее точечное решение в массив
            for (i=0;i<NN;i++)
            {
                all_X[step][i] = X[i];
            }
            step++;
        }
    }

    i_Create_Sol_lim();

    i_Itg_comb();  /// Находит внешнюю оценку для Itg,i и заносит в i_Ex_Sol[M-1+i]

//    i_print_Comb_lim(); /// нужно передавать имя файла? Да!
}



/// СООТВЕТСТВИЕ ПАРАМЕТРОВ СИСТЕМЫ
///    Ток вдоль трубы Itx_i  X[i]
///    Ток на границе грунт-труба Itg_i X[M-1+i]
///    Эл. потенциал на границе грунт-туба Utg_i X[2*M-1+i]
///    Эл. потенциал металла трубы Utm_i X[3*M-1+i]


/// ////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////// ПРОЕКТ ПЕРЕБОРА ВСЕХ ИНТЕРВАЛОВ - КОНЕЦ /////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////



int main()
{
int i;

    LoadDataSKZ(); /// Загружает данные из файла DataKSZ.csv

    /// использует константы
    /// const int D1_0 = 100; max кол-во строк данных о точках замера
    /// const int D2_0 = 6;   max кол-во столбцов данных одной точки замера

    /// вносит данные в глобальные переменные
    /// int Dt реальное кол-во строк даных о точках замера
    /// int Dskz  реальное кол-во строк даных о СКЗ (кол-во анодов)
    /// float DataSKZ[D1_0][D2_0] массив, содержащий исходную таблицу в числовом формате
    /// столбцы:
    ///  0 - абсолютная коор-та точк замера в км
    ///  1 - измеренный Potential V
    ///  2 - измеренная разность потенциалов - защитный потенциал
    ///  3 - сила тока катодной станции
    ///  4 - напряжение создаваемое СКЗ
    ///  5 - сопотивление грунта растеканию тока

    Init_FI_t(); /// Инициализирует данные по трубе
    /// использует константы
    /// M - принятое кол-во ФИ по трубе
    /// Dt - кол-во точек замера в исходных данных

    /// Ct Ом*м2 удельное сопротивление изоляции ПОКА КОНСТАНТА!!!

    /// использует глобальные переменные
    /// DataSKZ[][] массив с исходными данными (см. LoadDataSKZ())

    /// вносит данные в глобальные переменные
    /// L - общая длина рассматриваемого участка трубы
    /// Li = L / M Длина участка трубы соотв-го одному ФИ
    /// St = 2*pi*Rt*Li;  Площадь боковой поверхности, соответствующая одному ФИ
    /// FI_t[i][1..3] координаты центра ФИ с номером i
    /// FI_tg[i][1..3] координаты точки на грнанице ФИ с номером i
    /// C_t[i]=Ct/St;  сопротивление изоляции учстка трубы, соот-го ФИ с номером i


    Init_FI_a();
    /// Инициализирует данные об анодах (катодных станциях)
    /// вносит изменения в массивы и переменные
    /// FI_a[i][0..2] - коор-ты анода номером k (i=1 кол-во ФИ по аноду!!!)
    /// FI_a[i][3] - сила тока СКЗ с номером i
    /// FI_a[i][4] содержит номер ФИ, к центру которого подключена i-я СКЗ
/// ВВЕСТИ РАССТОЯНИЯ ДО ТРУБЫ И ГЛУБИНУ ДЛЯ КАЖДОГО АНОДА В ИСХОДНЫЕ ДАННЫЕ???

    Init_A();   /// строит расширенную матрицу А для СЛАУ и ее копию Ais для проверок

    step_Gauss();
    //mygauss(); /// решает СЛАУ, заданную расширенной матрицей А, результат содержится в X

    printResult("resall.txt","Результаты решения точечной системы");


/// ////////////////////////////////////////////////////////////////


    i_main_lim();
//    i_main_step();
//  i_main_LR(); /// Интервальные вычисления

///i_main_Gauss();
//    no_simply_int();
//    simply_int();

//    i_invA();
//    i_prov_inv();


    return 0;
}
