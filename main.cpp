#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <cstring>

using namespace std;
typedef double MyFloat;
const int   M=51,              // ���-�� ��������� ���������� �� �����
            N=5,                // Max ���-�� ��������� ���������� �� �����
            NN = 4*M-1;     // ����� �����������/��������� �������

//NN = 3;

const MyFloat pi = 3.141592653589793238462643383279,
            Rt = 0.557/2.0,     // ������� ������ �����
            Ht = 0.008,         // ������� ������ �����
            Zt = 1.5,           // ������� ��������� �����
            Sms = pi*(Rt*Rt - (Rt-Ht)*(Rt-Ht)),
                                // ������� �������������� ������� �����
            Ya = 200.0,         // ���������� �� ����� �� �����
            Za = 2.5,           // ������� ��������� ����� ��� ������������ �����
            ro_t = 2.45e-7,     // 2.45e-7 �������� ���� ����� ��*�
            ro_g = 35.0,        // �������� ���� ������ ��*�
            Ct = 5000.0,        // �������� ���� �������� ����� ��*�2
            Sigma_t = 1/ro_t,   // �������� ������������������ ������� �����
            Sigma_g = 1/ro_g,   // �������� ������������������ ������

            point_povr_Ct = 0.75, /// ������������� ����-�� ������ ������� ����������� �������� (0,1)
            koef_povr_Ct = 0.2,   /// ����������� �������������� �������� (0,1)
                                  /// 0 - ������ �����������, 1 - ��� �����������
            L_povr_Ct = 1000;       /// ����� ������� � ������������ ���������, �
/// �� ������� ����� �� ���������� L*point_povr_Ct �� ������ ����������� ������� �������� �������
/// ����� ������ L_povr_Ct �, ��������� �������� ����� �������� ������, ��� ������� ��������
/// ��������� �������������� ������������� ������� �� ������ ������� ����� Ct*koef_povr_Ct

/// !!! ���� �� ���������� �������� �����������, ��� ������ ����������� �� ����� ������
/// ��� �������, ������� �� ����� ������ ����� ��!!!


/// ����� �� �/� ���������� ������, � ������ ����������� ����� ���� ��� ������ ����� ��,
/// ������� ��� ����������� ��������� ������� �������� ���-�� �� � ������������ �
/// �������� ����������������� �����-� ���������������� ��������
/// count_povr_FI = int(L_povr_Ct/Li);  - ���-�� ��, ��� ������� ���������� Ct.
/// �������� ��������� count_povr_FI*Li*x_Ct*Ct = L_povr_Ct*koef_povr_Ct*Ct,
/// ������  x_Ct = (L_povr_Ct*koef_povr_Ct) / count_povr_FI*Li
/// ��� ��� ������ � Init_FIt()
/// int((M-1)*point_povr_Ct - ����� �� - ������ �����������
/// ������ ����� ������� �� �� ������� �����������
/// first_povr_FI = int((M-1)*point_povr_Ct - int(count_povr_FI/2)
/// ���������� ���������������� � ���������� ������, �����
/// first_povr_FI ��������� �������������, ��� (first_povr_FI + count_povr_FI) > M
///    for (i=0;i < count_povr_FI;i++)
///    {
///        FI_t[first_povr_FI+i][3]= Ct*koef_povr_Ct);
///    }

MyFloat     //I0,               // ��� �������� �������
            Li,                 // ����� ������� ����� �����-�� ������ ��
            St,                 // ������� ������� �����������, ��������������� ������ ��
            L,                  // ����� ����������� ������� �����
            varR1,             // ��������� ��� ���������� ��������������� ��������� � ���
//            C_t[M],           // �������� ������������������ �������� ����� ��*�2
            FI_t[M][4],         // ���������� ������� �� �� �����
            FI_a[N][5],         // ���������� ������� �� ������
            Pi[3],              // ����-�� �����-������� ��
            Pj[3],              // ����-�� �����-������� ��
            A[NN][NN+1],        // ����������� ������� �������
//            B[NN],              /// ������� ��������� ������ - ������������??
            Nv;                 // ����� ������� �������
//            Pks_x,              /// ���������� � ����� ����������� � ����� �������� ������� - ������������??
// R_self; /// - ������������??

int         i_ks;               /// ����� �� �����, � ������ �������� �������� ������� �� - ������������??

//////////// ����� ////////////////////

int         num[NN];            // ������ ��� ����� ������������ ��������
                                // ��� ������ �������� �������� � ������ ������
MyFloat     Ais[NN][NN+1],      // ������ ��� �������� �������� ����� ������� �������
                                // ��������� ��� ���������� �������
            X[NN],              // ������ ����������� �������
            Nev[NN];            // ������ ������� �������


//////////// ����� ////////////////////

/// //////// ������ ������  /////////////////
const int D1_0 = 100;   // max ���-�� ����� ������ � ������ ������
const int D2_0 = 6;     // max ���-�� �������� ������ ��� ����� ����� ������
int Dt,                 // �������� ���-�� ����� ����� � ������ ������
    Dskz;               // �������� ���-�� ����� ����� � ��� (���-�� ������)
MyFloat DataSKZ[D1_0][D2_0]; /// ������ ��� ������� ������ �� �����

/// //////// ������ ������  /////////////////

void MatrPrn_Nx2N(MyFloat A[NN][2*NN],int M,int N,char* fn)
/// ����� � ���� �������� ������� Nx2N
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
        fprintf(fp,"����������� ������� �������;\n");

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
/// ����� � ���� �������� ������� NxN
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
        fprintf(fp,"����������� ������� �������;\n");

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
/// ��������� begin
/// //////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////
#include <fenv.h>
typedef MyFloat Ifloat[2];
Ifloat  i_AB[NN][NN+1],   /// "�������" ����������� �������, ���������� � ����������� ��� � ������
        i_ABis[NN][NN+1], /// �������� ����������� ������� ��� ���������� �������
        i_X[NN],          /// "�������" ������ ������������, ���������� �������������� ��� ������ �������� ��-��
        i_Xis[NN];        /// ������ ������������ ��� �������� ������� �� �������� �-��

//        i_nev_is[NN],     /// ������ ������� �� �������� �-��
//        i_nev_tr[NN];     /// ������ ������� �� ����������� �-��

MyFloat i_AB_iw[NN][NN+1], /// ������, ���������� �������� ������ ���������� ����-� ����� �-�� - ��� ��������
        i_nev_is[NN],     /// ������ ������� �� �������� �-��
        i_nev_tr[NN],     /// ������ ������� �� ����������� �-��
        i_dist_sol[NN];   /// ������ ��� ��������� ��������������� �������



int i_num[NN];            /// ������ ��� ����� ������������� ���������� ��� ������ �������� ��-��
char i_format[24] = "%s%24.18f%s%24.18f%s"; /// ������ ������ ������������ ����-� � ����� ����
char i_format1[28] = "%16.8f%s%24.18f%s%24.18f\n"; /// ������ ������ ������������ ����������� � ��������� ����� ��� ���������� ��������


MyFloat disp_ro_g = 50, /// ���������, ������� ���� �� ����. ������ � % +/- � ��������� ���� ��������
        disp_Ct = 20;  /// ���������, ������� ���� �� ����. �������� ����� � % +/- � ��������� ���� ��������

Ifloat i_ro_g = {(ro_g*(100-disp_ro_g)/100),(ro_g*(100+disp_ro_g)/100)},
       i_sig = {1/i_ro_g[1],1/i_ro_g[0]},
       i_Ct = {Ct*(100-disp_Ct)/100,Ct*(100+disp_Ct)/100}, /// ���� �������� �� ��� ���� ��, ��� � ������ �������

       i_Ct_FI[M];   /// ��� ������� �� ���� ��, ������������ ������� FI_t[i][3]

MyFloat AE[NN][2*NN], /// ������� ���� �� ��� ���������� �������� � �
        E[NN][NN+1];  /// ���������� ��������� �-��


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
    MyFloat i_tl[4],i_tr[4],i_t_min,i_t_max; /// ��������� ����������,
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
        MyFloat i_tl[4],i_tr[4],i_t_min,i_t_max; /// ��������� ����������,
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
/// ����� ���� ��������
    max_num = l;
    for (k=l+1;k<NN;k++)
    {
///  ����� � ������������ ������ �������� �������� � ������ ������ ��������� �������� ��������
///        cout<<width_int(i_AB[l][max_num])<<endl;

        if (i_abs(i_AB[l][k]) > i_abs(i_AB[l][max_num])) max_num = k;
//        if ((i_abs(i_AB[l][k]) > i_abs(i_AB[l][max_num])) and (!i_incl0(i_AB[l][k]))) max_num = k;
    }
        if (max_num != l)
        {
//            ������ �������
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

/// ����� ���� ��������
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

        /// ��������� ������ ����������-������������� ����������� �������
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

        fprintf(fp,"%s%f%s%f%s\n","�������� ������������������ ������ (������������) i_ro = (",i_ro_g[0],";",i_ro_g[1],")");
        fprintf(fp,"%s%f%s%f%s\n","�������� ������������� ������ (������������) i_sig = (",i_sig[0],";",i_sig[1],")");
        fprintf(fp,"%s%f%s%f%s\n","�������� ������������� �������� ����� (������������) i_Ct = (",i_Ct[0],";",i_Ct[1],")");
        fprintf(fp,"%s\n","�������� ������� �������;");
        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,i_format,"(",i_ABis[i][j][0],";",i_ABis[i][j][1],") ");
            }
        fprintf(fp,"\n");
        }

        fprintf(fp,"%s\n","������� ��� �������� �������;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_Xis[i][0],";",i_Xis[i][1],") ");
        }
        fprintf(fp,"\n");

        fprintf(fp,"%s\n","����������� ������� �������;");

        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,i_format,"(",i_AB[i][j][0],";",i_AB[i][j][1],") ");
            }
        fprintf(fp,"\n");
        }

        fprintf(fp,"%s\n","������� ��� ����������� �������;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_X[i][0],";",i_X[i][1],") ");
        }
        fprintf(fp,"\n");

        fprintf(fp,"%s\n","������ ������������� ����� �������;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,"%i%s",i_num[i],"  ");
        }
        fprintf(fp,"\n");

        fprintf(fp,"%s\n","������ ���������� ����������� �������;");

        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,"%s%24.18f%s"," (",i_AB_iw[i][j],")");
            }
        fprintf(fp,"\n");
        }

        fprintf(fp,"%s\n","���� ������� ��� �������� �������;");

        for (i=0;i<NN;i++)
        {
//            fprintf(fp,i_format,"(",i_nev_is[i][0],";",i_nev_is[i][1],") ");
            fprintf(fp,i_format,"(",i_nev_is[i],";",i_nev_is[i],") ");

        }
        fprintf(fp,"\n");
        fprintf(fp,"%s\n","������ ����� ��� �������� �������;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_ABis[i][NN][0],";",i_ABis[i][NN][1],") ");
        }
        fprintf(fp,"\n");

        fprintf(fp,"%s\n","���� ������� ��� ����������� �������;");

        for (i=0;i<NN;i++)
        {
//            fprintf(fp,i_format,"(",i_nev_tr[i][0],";",i_nev_tr[i][1],") ");
            fprintf(fp,i_format,"(",i_nev_tr[i],";",i_nev_tr[i],") ");
        }
        fprintf(fp,"%\n");
        fprintf(fp,"%s\n","������ ����� ��� ����������� �������;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_AB[i][NN][0],";",i_AB[i][NN][1],") ");
        }

        fprintf(fp,"%\n");
        fprintf(fp,"%s\n","������������ ��������������� �������;");

        for (i=0;i<NN;i++)
        {
            fprintf(fp,i_format,"(",i_dist_sol[i],";",i_dist_sol[i],") ");
        }


        fprintf(fp,"%\n");


        fprintf(fp,"��� ����� ����� Itx_i;\n");

        for (i=0;i<M-1;i++) fprintf(fp,i_format,"(",i_Xis[i][0],";",i_Xis[i][1],") ");
        fprintf(fp,"%\n");

        fprintf(fp,"��� �� ������� �����-����� Itg_i, �;\n");
//    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[M-1+i],";");
        for (i=0;i<M;i++) fprintf(fp,i_format,"(",i_Xis[M-1+i][0],";",i_Xis[M-1+i][1],") ");

        fprintf(fp,"%\n");
/*
    fprintf(fp,"��������: ��������� ��� �� ������� �����-����� Itg_i;\n");
    Nv = 0;
    for (i=0;i<M;i++) Nv = Nv + X[M-1+i];
    fprintf(fp,outputFormat,Nv);
    fprintf(fp,"%\n");
*/
        fprintf(fp,"��. ��������� �� ������� �����-���� Utg_i;\n");
//    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[2*M-1+i],";");
        for (i=0;i<M;i++) fprintf(fp,i_format,"(",i_Xis[2*M-1+i][0],";",i_Xis[2*M-1+i][1],") ");
        fprintf(fp,"%\n");


        fprintf(fp,"��. ��������� ������� ����� Utm_i;\n");
//        for (i=0;i<M;i++) fprintf(fp,outputFormat,X[3*M-1+i],";");
        for (i=0;i<M;i++) fprintf(fp,i_format,"(",i_Xis[3*M-1+i][0],";",i_Xis[3*M-1+i][1],") ");

        fprintf(fp,"%\n");
        fprintf(fp,"�������� ��������� Utg_i-Utm_i;\n");
//    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[2*M-1+i]-X[3*M-1+i],";");
///���� ��������� ������ �������� �� �������� ������ �����
/// �� �� �������� �� ��������
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
            printf("%s %s \n","���������� ������� ���� ",curFileName);
        }
        else
        {
            printf("%s %s \n","�������� ������� ���� ",curFileName);
//        fprintf(fp,"�������� ��������� ����������;\n");
            fprintf(fp,"�������� ��������� Utg_i-Utm_i;\n");
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
            printf("%s %s \n","���������� ������� ���� ",curFileName);
        }
        else
        {
            printf("%s%s \n","�������� ������� ���� ",curFileName);
            fprintf(fp,"��������� ���� �� ������� �����-�����, ��/�2;\n");

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
/// �� ����� ��������� ����� ����� ,������� ��������!!!

        i_nev_tr[i] = i_dist(i_AB[i][NN],sum);
//        i_sub(i_nev_tr[i],i_AB[i][NN],sum);
        }
}

void det_alg_sol()
/// ����������� ������������ ��������������� �������
/// �� ������ ������������ i_ABis, � i_Xis ��������� �������� i_dist_sol
{
    int i,j;
    Ifloat t,t1,sum;
    /// �-�� ��������
    i_dist_sol[0] = i_dist(i_Xis[M-1+0],i_Xis[0]); /// Itg_0 = Itx_0
//i_dist_sol[0] = i_dist_sol[0] +100;
    for (i=1;i<M-1;i++)
    {
        i_add(t,i_Xis[M-1+i],i_Xis[i-1]);    /// Itg_i + Itx_(i-1)
        i_add(t1,i_Xis[i],i_ABis[i][NN]); /// Itx_i + Ia ����� ������������� ������ Ia �� ������ ����� ����

        i_dist_sol[i] = i_dist(t,t1); /// Itg_i + Itx_(i-1) = Itx_i + Ia
//i_dist_sol[i] = i_dist_sol[i] +100;
    }

    i_movp(t,-1);
    i_mul(t,t,i_Xis[M-2]);
    i_dist_sol[M-1] = i_dist(i_Xis[M-1+M-1],t); /// Itg_M = - Itx_(M-1)
//i_dist_sol[M-1] = i_dist_sol[2*M-2] +100;

//ERROOR!!! HERE!
    /// ������� ����� 3 ����
    for (i=0;i<M;i++)
    {

///        A[i+M][M+i-1]=-FI_t[i][3]/St; �� Init_A

        i_mul(t,i_ABis[i+M][M+i-1],i_Xis[M-1+i]); /// (Ct/St)*Itg_i - ���� � ������� � �������!
        i_add(t,i_Xis[2*M-1+i],t);                /// Utg_i - (Ct/St)*Itg_i ������� add!
        i_dist_sol[M+i] = i_dist(i_Xis[3*M-1+i],t); /// Utm_i = Utg_i - (Ct/St)*Itg_i
//i_dist_sol[M+i] = i_dist_sol[M+i] +300;
    }
    /// ������ ���
    for (i=0;i<M-1;i++)
    {
//        i_movp(t,-ro_t*Li/Sms); /// ������� �� ������� i_Ais
        i_movp(t,-1);
        i_mul(t,t,i_ABis[i+2*M][i]); /// ro_t*Li/Sms ������� �� ������� i_Ais
        i_mul(t,t,i_Xis[i]);
        i_sub(t1,i_Xis[3*M-1+i+1],i_Xis[3*M-1+i]);

/// /////////////// ����������� ��������� �������� �� ������������� - ����� �����!!
///        t1[0] = i_Xis[3*M-1+i+1][0] -i_Xis[3*M-1+i][0];
///        t1[1] = i_Xis[3*M-1+i+1][1] -i_Xis[3*M-1+i][1];
/// ///////////////
        i_dist_sol[2*M+i] = i_dist(t,t1); /// -Rt*Itx_i = Utm_(i+1) - Utm_i
//i_dist_sol[2*M+i] = i_dist_sol[2*M+i] +200;
    }

    /// ��������� ��� ���
    /// ������� ������� �� ��������, �� �������� � ����������� �����
    /// 4Pi*Utg_i = -Sum(Itg/(R*i_sig))  + Sum(Ia/(R*i_sig)) ����������� � �������� � ��������� sig �����
    /// A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;         /// ����������� ��� Utg_i
    /// A[i+3*M-1][4*M-1]=A[i+3*M-1][4*M-1]+FI_a[j][3]*(funR1(Pi,Pj)+funR2(Pi,Pj)); /// ����� ������ �����
    /// A[i+3*M-1][M-1+j] = varR1 + funR2(Pi,Pj); /// ���� ��� Itg_i


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
    /// ��� ����� �����                     Itx_i for M-1, X[i]
    /// ��� �� ������� �����-�����          Itg_i for M, X[M-1+i]
    /// ��. ��������� �� ������� �����-���� Utg_i for M, X[2*M-1+i]
    /// ��. ��������� ������� �����         Utm_i for M, X[3*M-1+i]
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
/// �� ����� ��������� ����� ����� ,������� ��������!!!

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

/// ������������ ���������� �������
///    ��� ����� ����� Itx_i  X[i]
///    ��� �� ������� �����-����� Itg_i X[M-1+i]
///    ��. ��������� �� ������� �����-���� Utg_i X[2*M-1+i]
///    ��. ��������� ������� ����� Utm_i X[3*M-1+i]
    }

int d;
for (d=0;d<1;d++)
{

/// ������� ��������� ������������ ��������� "�� �� �������� 1"

/// ��������� Utg �� ��� ��� ��� Itg � ��� Sigma_grunta
    for (i=0;i<M;i++)
    {
        i_movp(Utg[i],Ais[i+3*M-1][4*M-1]);
        for (j=0;j<M;j++)
        {
            i_movp(t,Ais[i+3*M-1][M-1+j]);
            i_mul(t,Itg[j],t);
            i_sub(Utg[i],Utg[i],t); /// ����-�� (1/Rij) ��� Itg �� �������� �-�� �� ��-��� ��� ���
        }
        i_movp(t,4*pi);
        i_mul(t,t,i_sig);
        i_div(Utg[i],Utg[i],t);
    }

/// ���������  ��� Utm �� ����.�����. 3 ����

    i_movp(t,-1/St);
    i_mul(t,t,i_Ct);

    for (i=0;i<M;i++)
    {
//        Utm[i] = Utg[i] - i_Ct*Itg[i]/St;

        i_movi(Utm[i],Itg[i]);
        i_mul(Utm[i],Utm[i],t);
        i_add(Utm[i],Utm[i],Utg[i]);
    }


/// ��������� Itg ����� �� �������. ������� 3 ���� �������� Ct, Utg Utm �������������

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

/// ��������� ������������ ���������� � ���� ������ i_Xis

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
/// ������������ ���������� �������
///    ��� ����� ����� Itx_i  X[i]
///    ��� �� ������� �����-����� Itg_i X[M-1+i]
///    ��. ��������� �� ������� �����-���� Utg_i X[2*M-1+i]
///    ��. ��������� ������� ����� Utm_i X[3*M-1+i]

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
        printf("%s %s \n","���������� ������� ���� ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","ii_Utg.txt");
//        fprintf(fp,"�������� ��������� ����������;\n");
        fprintf(fp,"������� no_simply;\n");

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
        printf("%s %s \n","���������� ������� ���� ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","ii_Utg.txt");
//        fprintf(fp,"�������� ��������� ����������;\n");
        fprintf(fp,"��� Itg, ������������ ��������;\n");
        for (i=0;i<M-1;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",Itg[i][0],";",Itg[i][1]);
        }
        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fprintf(fp,"��� Itg, �������� ��������;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[i+M-1],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }




   if((fp=fopen("nosimply_Utg.txt", "w")) == NULL)
    {
        printf("%s %s \n","���������� ������� ���� ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","ii_Utg.txt");
        fprintf(fp,"��������� Utg, ������������ ��������;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",Utg[i][0],";",Utg[i][1]);
        }
        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fprintf(fp,"��������� Utg, �������� ��������;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[2*M-1+i],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

   if((fp=fopen("nosimply_Utm.txt", "w")) == NULL)
    {
        printf("%s %s \n","���������� ������� ���� ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","ii_Utg.txt");
//        fprintf(fp,"�������� ��������� ����������;\n");
        fprintf(fp,"��������� Utm, ������������ ��������;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",Utm[i][0],";",Utm[i][1]);
        }
        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fprintf(fp,"��������� Utg, �������� ��������;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[3*M-1+i],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

/// �������� ����������� Upr

    for (i=0;i<M;i++)
    {
        i_sub(t1[i],Utg[i],Utm[i]);
    }
   if((fp=fopen("nosimply_Upr.txt", "w")) == NULL)
    {
        printf("%s %s \n","���������� ������� ���� ","ii_Utg.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","ii_Utg.txt");
//        fprintf(fp,"�������� ��������� ����������;\n");
        fprintf(fp,"�������� ��������� Upr, ������������ ��������;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",t1[i][0],";",t1[i][1]);
        }
        fprintf(fp,"%\n");
        char outputFormat1[28] = "%24.18f%s%24.18f%s";
        fprintf(fp,"��������� Upr, �������� ��������;\n");
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

/// ���������� ������������ ��������� "�� ��������"
    i_movp(t,St);
    for (i=0;i<M;i++)
    {
        i_movp(t1[i],St);
        i_mul(t1[i],t1[i],Upr[i]);
        i_div(t1[i],t1[i],i_Ct); /// ������� Itg = St*Upr/Ct

        i_movp(t2[i],1/St);
        i_mul(t2[i],t2[i],i_Ct); /// ������� Upr = Ct*Itg/St
        i_mul(t2[i],t2[i],Itg[i]);

        i_sub(Utg[i],t2[i],Utm[i]); /// Utg = Upr + Utm
    }
    for (i=0;i<M;i++)
    {
        i_movi(Itg[i],t1[i]);
        i_movi(Itg[i],t1[i]);
    }


/// ��������� ������������ ���������� � ���� ������ i_Xis

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
/// ������������ ���������� �������
///    ��� ����� ����� Itx_i  X[i]
///    ��� �� ������� �����-����� Itg_i X[M-1+i]
///    ��. ��������� �� ������� �����-���� Utg_i X[2*M-1+i]
///    ��. ��������� ������� ����� Utm_i X[3*M-1+i]




/// ������� "�� �������" ��������� ������������ ��������� - ���� ��-�� � ����
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
        printf("%s %s \n","���������� ������� ���� ","file.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","file.txt");
//        fprintf(fp,"�������� ��������� ����������;\n");
        fprintf(fp,"�������� ��������� Utg_i-Utm_i;\n");
        for (i=0;i<M;i++)
        {
            fprintf(fp,i_format1,FI_t[i][0],";",Upr[i][0],";",Upr[i][1]);
        }
        fprintf(fp,"%\n");
        fclose(fp);
    }

    if((fp=fopen("file2.txt", "w")) == NULL)
    {
        printf("%s %s \n","���������� ������� ���� ","file1.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","file1.txt");
//        fprintf(fp,"�������� ��������� ����������;\n");
        fprintf(fp,"��������� ���� Jtg;\n");
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

/// C������ dUdn = Itg/(Ct*i_sig_g*St)

    for (i=0;i<M;i++)
    {
        i_movp(dUdn[i],1/St);
        i_div(dUdn[i],dUdn[i],i_sig);
        i_div(dUdn[i],dUdn[i],i_Ct);
        i_mul(dUdn[i],dUdn[i],Itg[i]);
    }
/// �� ����� ���� ��� ��� ����
    for (i=0;i<M;i++)
    {
//        i_div(dUdn[i],Jtg[i],i_sig);
        i_mul(Upr[i],i_Ct,dUdn[i]);
        i_mul(Upr[i],i_sig,Upr[i]);
//        i_div(Jtg[i],Upr[i],i_Ct);
    }
/// ��������� ������������ ���������� � ���� ������ i_Xis

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
/// ������������ ���������� �������
///    ��� ����� ����� Itx_i  X[i]
///    ��� �� ������� �����-����� Itg_i X[M-1+i]
///    ��. ��������� �� ������� �����-���� Utg_i X[2*M-1+i]
///    ��. ��������� ������� ����� Utm_i X[3*M-1+i]

fesetround (FE_TONEAREST);
}







/// ��������� end


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
                (P_i[1]-P_j[1])+(P_i[2]-P_j[2])*(P_i[2]-P_j[2]))        /// �� ���� ��-�� x=0
                +
                1/sqrt((P_i[0]-(2*L-P_j[0]))*(P_i[0]-(2*L-P_j[0]))+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]-P_j[2])*(P_i[2]-P_j[2]))        /// �� ���� ��-�� x=L
                +
                 1/sqrt((P_i[0]+P_j[0])*(P_i[0]+P_j[0])+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]+P_j[2])*(P_i[2]+P_j[2]))        /// �� ���� ��-�� z=0 x=0
                +
                1/sqrt((P_i[0]-(2*L-P_j[0]))*(P_i[0]-(2*L-P_j[0]))+(P_i[1]-P_j[1])*
                (P_i[1]-P_j[1])+(P_i[2]+P_j[2])*(P_i[2]+P_j[2]));       /// �� ���� ��-�� z=0 x=L
}



void Init_FI_t()
{
    int i;
    L = DataSKZ[Dt-1][0]-DataSKZ[0][0] + DataSKZ[1][0]-DataSKZ[0][0]; // ����� ����� ����������� ������� �����
    cout<<"L= "<<L<<endl;
    Li = L / M;                         // ����� ������� ����� �����-�� ������ ��
    St = 2*pi*Rt*Li;                    // ������� ������� �����������, ��������������� ��
    cout<<"Li= "<<Li<<endl;
    varR1 = 2.0*log((sqrt(Rt*Rt+Li*Li)+Li)/Rt)/Li+        // ��������� ��� ����������
            2.0*log((sqrt(4*Zt*Zt+Li*Li)+Li)/(2*Zt))/Li;  // ��������������� ���������
    cout<<"varR1 = "<<varR1<<endl;
    for (i=0;i<M;i++)
    {
        FI_t[i][0]=i*Li+(Li/2);
        FI_t[i][1]=0;
        FI_t[i][2]=Zt;
        FI_t[i][3]= Ct;
        cout<<"Parametri truby  "<<FI_t[i][0]<<"  "<<FI_t[i][1]<<"  "<<FI_t[i][2]<<"  "<<FI_t[i][3]<<endl;
    }

/// ������ �����������
/// point_povr_Ct = 0.75, /// ������������� ����-�� ������ ������� ����������� �������� (0,1)
/// koef_povr_Ct = 0.2,   /// ����������� �������������� �������� (0,1)
///                             0 - ������ �����������, 1 - ��� �����������
/// L_povr_Ct = 100;       /// ����� ������� � ������������ ���������, �
/// �� ������� ����� �� ���������� L*point_povr_Ct �� ������ ����������� ������� �������� �������
/// ����� ������ L_povr_Ct �, ��������� �������� ����� �������� ������, ��� ������� ��������
/// ��������� �������������� ������������� ������� �� ������ ������� ����� Ct*koef_povr_Ct
/// ����� �� �/� ���������� ������, � ������ ����������� ����� ���� ��� ������ ����� ��,
/// ������� ��� ����������� ��������� ������� �������� ���-�� �� � ������������ �
/// �������� ����������������� �����-� ���������������� ��������
/// count_povr_FI = int(L_povr_Ct/Li)+1;  - ���-�� ��, ��� ������� ���������� Ct.
/// �������� ��������� count_povr_FI*Li*x_Ct*Ct = L_povr_Ct*koef_povr_Ct*Ct,
/// ������  x_Ct = (L_povr_Ct*koef_povr_Ct) / count_povr_FI*Li
/// ��� ��� ������ � Init_FI_t()

/// int((M-1)*point_povr_Ct - ����� �� - ������ �����������
/// ������ ����� ������� �� �� ������� �����������
/// first_povr_FI = int((M-1)*point_povr_Ct - int(count_povr_FI/2)


/// ���������� ���������������� � ���������� ������, �����
/// first_povr_FI ��������� �������������, ��� (first_povr_FI + count_povr_FI) > M

/// ������������� ������������ koef_povr_Ct � ������������ � ������ �������
/// (L_povr_Ct*koef_povr_Ct) / count_povr_FI*Li
/// �� ������ ���-�� �� �� ���������, ��� koef_povr_Ct = 1 �� ���������� �����
/// �������������� ��!!!
/// ������� ���� ��� ���� ���������

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
    /// ���-�� ����� Dt � ��� Dskz ��������� � ��������� LoadDataSKz()
    j = 0;
    for (i=0;i<Dt;i++)
    {
        if (DataSKZ[i][3] != 0) /// ������� ���� ���� �������� �� ���� (��� ���)
        {
//            FI_a[j][0] = DataSKZ[i][0] - DataSKZ[0][0]+0.5*Li;
            FI_a[j][0] = DataSKZ[i][0] - DataSKZ[0][0]+(DataSKZ[1][0] - DataSKZ[0][0])/2;
cout<<"FI_a[j][0] = "<<FI_a[j][0]<<endl;
            FI_a[j][1] = Ya;
            FI_a[j][2] = Za;
            FI_a[j][3] = DataSKZ[i][3];         // ���� ���� ���

cout<<"i = "<<i<<endl;
cout<<"j = "<<j<<endl;
cout<<"FI_a[j][3] = "<<FI_a[j][3]<<endl;
cout<<"DataSKZ[i][3] = "<<DataSKZ[i][3]<<endl;

            FI_a[j][4] = int(FI_a[j][0]/Li);  // ����� ��, ������������� � ��
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
    //������ ������� ��� �� �� �����
    A[0][0]=-1;         /// ����������� ��� Itx_0
    A[0][M-1]=1;        /// ����������� ��� Itg_0

    for (i=1;i<=M-2;i++)
    {
        A[i][i-1]=1;    /// ����������� ��� Itx_(i-1)
        A[i][i]=-1;     /// ����������� ��� Itx_i
        A[i][M+i-1]=1;  /// ����������� ��� Itg_i
    }
    A[M-1][M-2]=1;      /// ����������� ��� Itx_(M-1)
    A[M-1][2*M-2]=1;    /// ����������� ��� Itg_M
    for (i=0;i<Dskz;i++)
    {
        cout<<"FI_a[i][3] = "<<FI_a[i][3]<<endl;
        cout<<"int(FI_a[i][4]) = "<<int(FI_a[i][4])<<endl;

        A[int(FI_a[i][4])][NN] = FI_a[i][3];  /// ������������� ������
                                              /// � ������ ����� �������

        cout<<"A[int(FI_a[i][4])][NN] = "<<A[int(FI_a[i][4])][NN]<<endl;
    }


// ��������� ������� 3 ����
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1]=-FI_t[i][3]/St; /// ����������� ��� Itg_i Ct/St
        A[i+M][i+2*M-1]=1;            /// ����������� ��� Utg_i
        A[i+M][i+3*M-1]=-1;           /// ����������� ��� Utm_i
    }

     // ����� ��� ����� ��������� ���������� �����������
    for (i=0;i<M-1;i++)
    {
        A[i+2*M][i]=ro_t*Li/Sms;    /// ����������� ��� Itx_(i) +++!!!
        A[i+2*M][i+3*M-1]=-1;       /// ����������� ��� Utx_(i)
        A[i+2*M][i+3*M]=1;          /// ����������� ��� Utx_(i+1)
    }


    // ��������� ��� ������������������ ��������
    for (i=0;i<M;i++)
        {
        A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;   /// ����������� ��� Utg_i

        Pi[0] = FI_t[i][0];
        Pi[1] = FI_t[i][1];
        Pi[2] = FI_t[i][2];

        for (j=0;j<Dskz;j++)
        {
                Pj[0] = FI_a[j][0];
                Pj[1] = FI_a[j][1];
                Pj[2] = FI_a[j][2];

                // ��������� ����������� ������ �����
                A[i+3*M-1][4*M-1]=A[i+3*M-1][4*M-1]+FI_a[j][3]*(funR1(Pi,Pj)+funR2(Pi,Pj));
        }

        for (j=0;j<M;j++)
            {
                Pj[0] = FI_t[j][0];
                Pj[1] = FI_t[j][1];
                Pj[2] = FI_t[j][2];

                /// ��������� ����������� ��� Itg_i
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
        num[i] = i;              // ������� ������� ����������
        for (j=0;j<NN+1;j++)
        {
            Ais[i][j] = A[i][j]; // �������� ������� ��� ���������� �������
        }
    }
}

//////////// ����� ////////////////////
//////////// ����� ////////////////////


void gl_el(int k)
{
    int i,j,gl;
    MyFloat gll;
//    double gll;
    gl = 0;

    // ���� ������� ������� � ������ �
    for (i=0;i<NN;i++)
    {
        if (abs(A[k][i])>abs(A[k][gl]))
            {
                gl=i;
            }
    }

    if (gl!=k) // ���� ������� ������� �� �� ���������
    {
        for (i=0;i<NN;i++) // ������ ������� �������
        {
            gll=A[i][k];
            A[i][k] = A[i][gl];
            A[i][gl]=gll;
        }
        //��������� ��������� ������� ������ ����������
        i=num[k];
        num[k]=num[gl];
        num[gl]=i;
    }

    // �������� 1 ��� �������� �������� ������ �
    gll = A[k][k];
    for (i=0;i<NN+1;i++)
    {
        A[k][i] = A[k][i]/gll;
    }
    // �������� ���� � ������� ��� ������� ������
    for (i=0;i<NN;i++)
        if ((i!=k) & (A[i][k]!=0)) // �� ������� ��� ������� � ��, ��� ��� 0
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

    // ���� ������� ������� � ������ �
    for (i=0;i<NN;i++)
    {
        if (abs(A[k][i])>abs(A[k][gl]))
            {
                gl=i;
            }
    }

    if (gl!=k) // ���� ������� ������� �� �� ���������
    {
        for (i=0;i<NN;i++) // ������ ������� �������
        {
            gll=A[i][k];
            A[i][k] = A[i][gl];
            A[i][gl]=gll;
        }
        //��������� ��������� ������� ������ ����������
        i=num[k];
        num[k]=num[gl];
        num[gl]=i;
    }

    // �������� 1 ��� �������� �������� ������ �
    gll = A[k][k];
    for (i=0;i<NN+1;i++)
    {
        A[k][i] = A[k][i]/gll;
    }
    // �������� ���� � ������� ��� ������� ������
    for (i=0;i<NN;i++)
        if ((i!=k) & (A[i][k]!=0)) // �� ������� ��� ������� � ��, ��� ��� 0
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
        fprintf(fp,"%s%f\n","�������� ������������������ ������ (��������) ro_g = ",ro_g);
        fprintf(fp,"%s%f\n","�������� ������������� ������ (��������) Sigma_g = ",Sigma_g);
        fprintf(fp,"%s%f%s%f%s\n","�������� ������������� �������� ����� (������������) Ct = ",i_Ct[0],";",i_Ct[1],")");
        fprintf(fp,"%s%f\n","�������� ������������� �������� ����� (��������) Ct = ",Ct);
        fprintf(fp,"%s\n","�������� ������� �������;");

        fprintf(fp,"Ais - ����������� ������� �������;\n");
        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,outputFormat,Ais[i][j],";");
            }
            fprintf(fp,"%\n");
        }

        fprintf(fp,"A - ��������������� ������� �������;\n");
        for (i=0;i<NN;i++)
        {
            for (j=0;j<NN+1;j++)
            {
                fprintf(fp,outputFormat,A[i][j],";");
            }
            fprintf(fp,"%\n");
        }
    fprintf(fp,"X - ������� �������;\n");
    for (i=0;i<NN;i++) fprintf(fp,outputFormat,X[i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"������ ������� �������;\n");
    for (i=0;i<NN;i++) fprintf(fp,outputFormat,Nev[i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"����� ������� ������� �������;\n");
    fprintf(fp,outputFormat,Nv);
    fprintf(fp,"%\n");
    fprintf(fp,"��� ����� ����� Itx_i;\n");
    for (i=0;i<M-1;i++) fprintf(fp,outputFormat,X[i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"��� �� ������� �����-����� Itg_i;\n");
    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[M-1+i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"��������: ��������� ��� �� ������� �����-����� Itg_i;\n");
    Nv = 0;
    for (i=0;i<M;i++) Nv = Nv + X[M-1+i];
    fprintf(fp,outputFormat,Nv);
    fprintf(fp,"%\n");
    fprintf(fp,"��. ��������� �� ������� �����-���� Utg_i;\n");
    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[2*M-1+i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"��. ��������� ������� ����� Utm_i;\n");
    for (i=0;i<M;i++) fprintf(fp,outputFormat,X[3*M-1+i],";");
    fprintf(fp,"%\n");
    fprintf(fp,"�������� ��������� Utg_i-Utm_i;\n");
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
        fprintf(fp,"�������� ��������� Utg_i-Utm_i;\n");
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
        fprintf(fp,"��������� ���� �� ������� �����-�����;\n");
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
        fprintf(fp,"�������� ������������� ������������� �������� �����;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",FI_t[i][3],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

}

//////////// ����� ////////////////////
//////////// ����� ////////////////////

/// //////// ������ ������  /////////////////

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
        strncat(buff,s[0],1);                       /// ��������� ";" � ����� ������
        cout<<"buff["<<col<<"] = ";
        lenbuff = strlen(buff);                     /// ����� ������
        for (i=0;i<lenbuff;i++)
        {
            if (buff[i] != *s[0])                   /// ���� ������ - �� ����������� ��������
            {
                cout<<buff[i];
                strncat(bf,buff+i,1);               /// ���������� �������� ������� ���������� �����
            }
            else                                    /// �����
            {
                DataSKZ[row][col] = atof(bf);       /// ����������� ��������� ����� � ���� double
                                                    /// � ������� � ������
                cout<<"="<<DataSKZ[row][col]<<"; ";
                col++;                              /// ��������� � ���������� ������
                strcpy(bf,"");                      /// ������� ������ ��� ����� ���������� �����
            }
        }
        row++;
        col = 0;
        cout<<endl;
    }
   for (i=0;i<row;i++) if (DataSKZ[i][3] != 0) Dskz++;     /// ������� ���-�� ���
   Dt = row;
//    cout<<"Strok count   "<<Dt<<endl;
//    cout<<"SKZ count   "<<Dskz<<endl;
    fin.close();
//* ���������� ���������� �������
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

/// //////// ������ ������  /////////////////

void i_InitAB()
/// ��� ������ �� ����������� ��� ���� ��
/// ������ ������������ ������� i_AB � i_ABis
/// �� ������ �������� ������� Ais � ������������ i_sig i_Ct
/// ���� ������ �� ������ ������� ��������
/// ���������� ����� ��������� ����-�� ������
/// � ������������� ������������
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
/// ������ �������������� �� �� ���� ������ � ���
    for (i=0;i<M;i++)
    {
//        A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;
        i_movp(t,Sigma_g);
        i_div(t,i_sig,t);
        i_mul(i_AB[i+3*M-1][2*M+i-1],i_AB[i+3*M-1][2*M+i-1],t);

        i_movi(i_ABis[i+3*M-1][2*M+i-1],i_AB[i+3*M-1][2*M+i-1]); /// �������� ��������� � ���� i_ABis
    }

/// ������ �������������� �� �� ���� �������� ����� � ��������� ������� 3 ����

/*
    �������� �� ���� Init_A
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1]=-FI_t[i][3]/St; /// ����������� ��� Itg_i Ct/St
        A[i+M][i+2*M-1]=1;            /// ����������� ��� Utg_i
        A[i+M][i+3*M-1]=-1;           /// ����������� ��� Utm_i
    }
*/
    for (i=0;i<M;i++)
    {
  /*      i_movp(t,Ct); /// ����� ������ ���, ���� ��� "�����" ����� ���� ������!!!
        i_div(t,i_Ct,t);
        i_mul(i_AB[i+M][M+i-1],i_AB[i+M][M+i-1],t); /// ����������� ��� Itg_i Ct/St


*/
        /// ������!
        i_movp(i_AB[i+M][M+i-1],-1/St);
        i_mul(i_AB[i+M][M+i-1],i_AB[i+M][M+i-1],i_Ct);

        i_movi(i_ABis[i+M][M+i-1],i_AB[i+M][M+i-1]);/// �������� ��������� � ���� i_ABis
    }
}


void i_main_Gauss()
{
    int i,j;
    fesetround (FE_DOWNWARD);

    i_InitAB();  /// ������ ������������ �������
    i_Gauss();
    i_ProvABis();
    i_ProvAB();
    i_printResult();

    fesetround (FE_TONEAREST);


}

/// ////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////// ������ ����� ������ ���������� - ������ ////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////

void i_InitABl()
/// ������ �������� ������� ��� ����� ������� ���������
{
    int i,j;
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            A[i][j] = i_ABis[i][j][0];
            Ais[i][j] = i_ABis[i][j][0];
        }
        num[i] = i;              // ������� ������� ����������
    }
}

void i_InitABr()
/// ������ �������� ������� ��� ������ ������� ���������
{
    int i,j;
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            A[i][j] = i_ABis[i][j][1];
            Ais[i][j] = i_ABis[i][j][1];
        }
        num[i] = i;              // ������� ������� ����������
    }
}


void i_main_LR()
{
//    cur_h_Ct = i_wid(i_Ct_h[i]) // i_Ct_N; - � ��������� �������� ����������

    int i,j;
//    fesetround (FE_DOWNWARD);

    i_InitAB();  /// ������ ������������ �������

    i_InitABl(); /// ������� � � ����� �������� �-��
    mygauss();   /// ������ ������ ��� ����� ������� ���������

    printResult("resall_l.txt","���������� �������� ����� �������"); /// ���������� ������� ����� �������

/// ��������� ���������� ��� ����� �������
    for (i=0;i<NN;i++)
    {
        i_Xis[i][0] = X[i];
//        Xtmp[i][0] = X[i];
    }


    i_InitABr(); /// ������� � �� ������ �������� �-��
    mygauss();   /// ������ ������ ��� ����� ������� ���������

   printResult("resall_r.txt","���������� �������� ������ �������"); /// ���������� ������� ������ �������


/// ��������� ���������� ��� ����� �������
/// ������� ������������ �������
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

/// ��������� �������
    i_ProvABis();
    det_alg_sol();
    i_printResult();

 //   fesetround (FE_TONEAREST);
}
/// ////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////// ������ ����� ������ ���������� - ����� ////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////



/// ////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////// ������ �������� ���� ���������� - ������ ////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////



const unsigned int i_Ct_N = 1,     /// ����� ����� ��������� ��������� i_Ct_FI[i]
                   i_sig_N = 1;   /// ����� ����� ��������� ��������� i_sig

const unsigned long i_Ct_cmb_cnt = 1024*64; /// ����� ���������� �� i_Ct (i_Ct_N+1+1)^M 9765625  390625 3125
const int  i_sig_cmb_cnt = i_sig_N+1; /// ����� ���������� �� i_sig (i_sig_N+1)


unsigned long step;    ///    ������� ��� �������� ���������� - ����� ��������� ����� ���������� ��� ������� �������������� ������

MyFloat cur_Ct[M],  /// ��� �������� ������� �������� �� ��� ���������� ����������
        cur_sig,    /// ��� �������� �������� ��������� sigma_g
        h_Ct[M],       /// ��� ��������
        h_sig,      /// ��� ��������
        all_X[4][NN+2];  /// ��� ������� ��� ���� ���������� - ������ ��� ������� lim
//        all_Ct[i_Ct_cmb_cnt][M];  /// ��� ���������� Ct
//        all_X[i_Ct_N+2][NN+2];

Ifloat i_Ext_Sol[NN+M]; /// ������ ��� ��������� ������� ������ ������� all_X[] + ��������� ����������.
float eps_Ct = 1,       /// ��������, ����� ������������ Ct
      eps_sig = 0.0001; /// ��������, ����� ������������ sigma_g
                        /// ���������� ��� ������� "���������" � ������� ��������� � ������
    FILE* fn_step;


void i_Init_Ct()
/// �������� ������������ Ct ��� ������� ��, ��������� � ��������� �� ����� ������� ����� ��� � Init_FI_t()
{
    int i;
//    i_Ct = {Ct*(100-disp_Ct)/100,Ct*(100+disp_Ct)/100}, /// ���� �������� �� ��� ���� ��, ��� � ������ ������
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
/// ������ ������������ ������� i_AB � i_ABis ��� ��������
/// �� ������ �������� ������� Ais � ������������ i_sig i_Ct_FI[M]
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
/// ������ �������������� �� �� ���� ������ � ���
    for (i=0;i<M;i++)
    {
//        A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;
        i_movp(t,4*pi);
        i_mul(i_AB[i+3*M-1][2*M+i-1],t,i_sig);

    }

/// ������ �������������� �� �� ���� �������� ����� � ��������� ������� 3 ����

/*
    �������� �� ���� Init_A
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1]=-FI_t[i][3]/St; /// ����������� ��� Itg_i Ct/St
        A[i+M][i+2*M-1]=1;            /// ����������� ��� Utg_i
        A[i+M][i+3*M-1]=-1;           /// ����������� ��� Utm_i
    }
*/
    for (i=0;i<M;i++)
    {
        i_movp(i_AB[i+M][M+i-1],-1/St);
        i_mul(i_AB[i+M][M+i-1],i_AB[i+M][M+i-1],i_Ct_FI[i]); /// ����������� ��� Itg_i -Ct/St
    }

/// �������� i_AB[i][j] � i_ABis[i][j]
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            i_movi(i_ABis[i][j],i_AB[i][j]);
        }
    }
}

void i_InitStepAB(MyFloat cur_Ct[M],MyFloat cur_sig)
/// cur_Ct[i] �������� �������� �������� �� ��� i-�� ��
/// ������ �������� ������� A
/// ��� �������� �������� Ct[i] � sigma_g
/// ���� �� � i_Ais �������� ��� ���-� ��������!!! � ��������� ���������
{
    unsigned int i,j;
    Ifloat t;
    for (i=0;i<NN;i++)
    {
        for (j=0;j<NN+1;j++)
        {
            A[i][j] = Ais[i][j]; /// ����� ����-�� �� ���������, "�������" �������� �-��
        }
        num[i] = i;
    }

    for (i=0;i<M;i++)
    {
        A[i+3*M-1][2*M+i-1]=4*pi*cur_sig; /// �� ���� ������ � ��� A[i+3*M-1][2*M+i-1]=4*pi*Sigma_g;
        A[i+M][M+i-1] = -cur_Ct[i]/St; /// ����������� ��� Itg_i Ct/St  ���� ��� 3 ����
    }

/// ������ ������� �������� �� ���� �������� ����� � ��������� ������� 3 ����
/*
    �������� �� ���� Init_A
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1]=-FI_t[i][3]/St; /// ����������� ��� Itg_i Ct/St
        A[i+M][i+2*M-1]=1;            /// ����������� ��� Utg_i
        A[i+M][i+3*M-1]=-1;           /// ����������� ��� Utm_i
    }
*/
/*
    for (i=0;i<M;i++)
    {
        A[i+M][M+i-1] = -cur_Ct[i]/St; /// ����������� ��� Itg_i Ct/St  ���� ��� 3 ����
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
/// ��� ��������� ������ ��������� Ct, ������ �������� ���� � ���� ���� � ��� ������� �������������
/// ����� ������ ������� ����������
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
            cur_Ct[p+1] = cur_Ct[p+1] + h_Ct[p+1]; /// ��� ������ �� ����� �� ���������, ��������!
            Iter(p+1);

        }

        cur_Ct[p+1] = i_Ct_FI[p+1][1]; /// ������� ������� �� ��������� � ����� �� ��������� ���������
                                       /// ��� �����������  ������� �������� � ���������� � lim
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
        /// ��������� �������� ������� ��� ���� �������� cur_Ct � cur_sig
        i_InitStepAB(cur_Ct,cur_sig);

        /// ������ ������� � �������� ��������� cur_Ct � cur_sig
        step_Gauss(); /// ������� ����������������, ����� �� ������� ������� ��� ���������
//        mygauss();


        /// ����������� ���������� �������� ������� �� ���� � ���
        /// i_Ext_Sol ����� ��������� �������������� ���� � ��� ���������� ������� X
        for (w=0;w<NN;w++)
        {
            if (X[w] < i_Ext_Sol[w][0])
            {
                i_Ext_Sol[w][0] = X[w];
/// ���� ����� ��������� ����� � �� ��� ������� ����� ���� � ��� - ��� ������ ��� ������ ���������!!!
///                i_Ext_sig = cur_sig;/// �������� �����
/// � ����� ���� �� ����� �������� �� ���! 17.40 17.08
                /// �������� ���������� �����

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

print_cur_step(); ///  ������� ���������� �������� ������� � ����.
        step++;
        printf("%u%s%u\n",step," �� ",i_sig_cmb_cnt*i_Ct_cmb_cnt);


    }
}

void i_print_Comb()
{
    unsigned int i,j,k;
    /// ������� �������� � ��������� ���� ���������
    FILE *fp;
  /*  if((fp=fopen("all_X.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file");
    }
    else
    {
        printf("%s \n","Opened file ");
        fprintf(fp,"������� all_X \n");

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
        fprintf(fp,"���������� sigma, Ct � ��������������� �������� ���������;\n");

        cur_sig = i_sig[1];    /// �������� ����� ���������� � �������� �����������
                      ///  ����� ��� � ���� ���� � ������ � � ����� ������ �����������
        for (k=0;k<i_sig_N+1;k++) /// ���� �������� �������� sigma_g
        {
            for (j=0;j<step;j++) /// ���� �������� ���������� Ct �� ��
            {
                fprintf(fp,"%24.18f%s",cur_sig,";  ");
                for (i=0;i<M;i++) /// ���� �������� Ct ��
                {
                    fprintf(fp,"%24.18f%s",all_Ct[j][i],";  ");
                }
                fprintf(fp,"\n");
                fprintf(fp,"%24.18f%s",0,";  ");
                for (i=0;i<M;i++) /// ���� �������� Ct ��
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
        fprintf(fp,"������� ������ ������� Ext_Sol[NN+M] \n");

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
/// ���������� all_X � ������ ���� - ������������ �������� ������������� �������
/// ������� ��������� ���� � ��� �������� ��� ������� ������� (��)
/// ��������� ������� � Ifloat i_Ext_Sol[M] - ������� ������ �������
/// �������� ��� �������� ���������, ��� � all_X ���!!!
{
    unsigned int i,j;



    for (i=0;i<NN;i++) /// ���� �������� ��������� ������� all_X - ����������� �������
    {
        i_Ext_Sol[i][0] = all_X[0][i];
        i_Ext_Sol[i][1] = all_X[0][i];
        for (j=1;j<(i_sig_cmb_cnt)*step;j++) /// ���� �������� ������� (����� ������� all_X)
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


    for (i=0;i<M;i++) /// ���� �������� ��������� ������� all_X - ����������� �������
    {
        i_Ext_Sol[NN+i][0] = all_X[0][2*M-1+NN+i] - all_X[0][3*M-1+NN+i];
        i_Ext_Sol[NN+i][1] = i_Ext_Sol[NN+i][0];

        for (j=1;j<(i_sig_cmb_cnt)*step;j++) /// ���� �������� ������� (����� ������� all_X)
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
/// ������������ ���������
/// �� ����� ���������� � ������������ ���������� ��,
/// � �� ������� �� ��� ����� ��������� �������� ��� �������� ����
/// ��������� ��������� - ���������. ��� �������� ������� � ������� ��������?

/// ������������ ������� ���������� ����� �������
/// �������� ����� � ����� �� ����������

/// ����������
/// i_Ct_FI[M] - ������������ �� ���� ���� ��� ������� ��
/// i_sig - ������������ �� ���� ������
/// i_Ct_N ���-�� ����� �� ��������� Ct
/// i_sig_N ���-�� ����� �� ��������� sigma
{
    unsigned int i,j,k;

    i_Init_Ct();  /// � ������ i_Ct_FI ������� ������������� Ct ��� ������� ��
                  /// �� ������� i_Ct_FI � i_sig

    i_Init_ABis(); /// ������ ������������ ������� i_ABis ��� �������� �������
/// normal i_printResult();
    /// ����� �� i_Ct_FI ,����� � ����� h_Ct �������� ������� cur_Ct,
    /// h_Ct ���� ��� ������� i h_Ct = (i_wid(i_Ct_FI[i]))/i_Ct_N

    step = 0;
    for (i=0;i<M;i++)
    {
        cur_Ct[i] = i_Ct_FI[i][0]; /// �������� � ����� ������� Ct
    }

    for (i=0;i<NN+M;i++) ///
    {
        i_Ext_Sol[i][0] = 100;  /// ������ ��������� �������� ������� � ����� ��������
        i_Ext_Sol[i][1] = -100; /// ��� ������������ ������
    }


    cur_sig = i_sig[0];    /// �������� ����� ���������� � �������� �����������
                           ///  ����� ��� � ���� ���� � ������ � � ����� ������ ����������� ���!!! ��� ���������!!!
    h_sig = i_wid(i_sig) / i_sig_N;

/// ��������� ���� ��� ������ fn_step

    if((fn_step=fopen("i_Itg_step.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file i_Itg_step.txt");
    }
    else
    {
        printf("%s \n","Opened file i_Itg_step.txt");
        fprintf(fn_step,"���������� ��� ���� Itg �����-����� \n");
    }

    for (k=0;k<i_sig_N;k++) /// ���� �������� �������� sigma_g
    {
        Iter(-1);
//        printf("%f\n",cur_sig);
        cur_sig = cur_sig + h_sig;
    }

    cur_sig = i_sig[1]; /// ���������� Ct ��� ������� �������� � ���������� � lim ������ ������� �� ��������� � ����� �� i_sig
    Iter(-1);



    printf("%s%f%s%f \n","i_sig = ",i_sig[0],"  ",i_sig[1]);

    fclose(fn_step);
    i_print_Comb(); /// ����� ���������� ��� �����?


}




void i_Create_Sol_lim()
/// ���������� all_X � ������ ���� - ������������ �������� ������������� �������
/// ������� ��������� ���� � ��� �������� ��� ������� ������� (��)
/// ��������� ������� � Ifloat i_Ext_Sol[M] - ������� ������ �������

{
    unsigned int i,j;



    for (i=0;i<NN;i++) /// ���� �������� ��������� ������� all_X - ����������� �������
    {
        i_Ext_Sol[i][0] = all_X[0][i];
        i_Ext_Sol[i][1] = all_X[0][i];

        for (j=1;j<4;j++) /// ���� �������� ������� (����� ������� all_X)
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

    for (i=0;i<M;i++) /// ���� �������� ��������� ������� all_X - ����������� �������
    {
//        i_Ext_Sol[NN+i][0] = all_X[0][2*M-1+NN+i] - all_X[0][3*M-1+NN+i];
        i_Ext_Sol[NN+i][0] = all_X[0][2*M-1+i] - all_X[0][3*M-1+i];
        i_Ext_Sol[NN+i][1] = i_Ext_Sol[NN+i][0];

        for (j=1;j<4;j++) /// ���� �������� ������� (����� ������� all_X)
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
    /// ������� �������� � ��������� ���� ���������
    FILE *fp;
    if((fp=fopen("all_X_lim.txt", "w")) == NULL)
    {
        printf("%s \n","Not opened file all_X_lim.txt");
    }
    else
    {
        printf("%s \n","Opened file all_X_lim.txt");
        fprintf(fp,"������� all_X \n");

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
        fprintf(fp,"���������� sigma, Ct � ��������������� �������� ���������;\n");
        s = 0;
        for (j=0;j<2;j++)
        {
            for (i=0;i<M;i++)
            {
                cur_Ct[i] = i_Ct_FI[i][j]; /// �������� � ����� ������� Ct
            }
            for (k=0;k<2;k++)
            {
                cur_sig = i_sig[k];

            /// ������ ����� � ����

/// ////////////////////////////////
                fprintf(fp,"%24.18f%s",cur_sig,";  ");
                for (i=0;i<M;i++) /// ���� �������� Ct ��
                {
                    fprintf(fp,"%24.18f%s",cur_Ct[i],";  ");
                }
                fprintf(fp,"\n");
                fprintf(fp,"%24.18f%s",0,";  ");
                for (i=0;i<M;i++) /// ���� �������� Ct ��
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
        fprintf(fp,"������� ������ ������� Ext_Sol[NN+M] \n");

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

/// ����� Upr � ������� ��� �������
        if((fp=fopen("i_Upr_lim.txt", "w")) == NULL)
        {
            printf("%s %s \n","Not opened file ","i_Upr_lim.txt");
        }
        else
        {
            printf("%s %s \n","Opened file ","i_Upr_lim.txt");
//        fprintf(fp,"�������� ��������� ����������;\n");
            fprintf(fp,"�������� ��������� Utg_i-Utm_i;\n");
            for (i=0;i<M;i++)
            {
                fprintf(fp,i_format1,FI_t[i][0],";",i_Ext_Sol[NN+i][0],";",i_Ext_Sol[NN+i][1]);
            }
            fprintf(fp,"%\n");
            fclose(fp);
        }

/// ����� Jtg � ������� ��� �������
        if((fp=fopen("i_Jtg_lim.txt", "w")) == NULL)
        {
            printf("%s %s \n","Not opened file ","i_Jtg_lim.txt");
        }
        else
        {
            printf("%s%s \n","Opened file ","i_Jtg_lim.txt");
            fprintf(fp,"��������� ���� �� ������� �����-�����, ��/�2;\n");


            for (i=0;i<M;i++) fprintf(fp,i_format1,FI_t[i][0],";",(1000/St)*i_Ext_Sol[M-1+i][0],";",(1000/St)*i_Ext_Sol[M-1+i][1]);
//        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",1000*X[M-1+i]/St,";\n");
            fprintf(fp,"%\n");
            fclose(fp);
        }

/// ����� Ct � ������� ��� �������
        if((fp=fopen("i_Ct_lim.txt", "w")) == NULL)
        {
            printf("%s %s \n","Not opened file ","i_Ct_lim.txt");
        }
        else
        {
            printf("%s%s \n","Opened file ","i_Ct_lim.txt");
            fprintf(fp,"�������� ������������� �������� �����;\n");


            for (i=0;i<M;i++) fprintf(fp,i_format1,FI_t[i][0],";",i_Ct_FI[i][0],";",i_Ct_FI[i][1]);
//        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",1000*X[M-1+i]/St,";\n");
            fprintf(fp,"%\n");
            fclose(fp);
        }





/*
    char outputFormat1[28] = "%24.16f%s%24.16f%s";
    if((fp=fopen("Upr_lim.txt", "w")) == NULL)
    {
        printf("%s %s \n","���������� ������� ���� ","Upr.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","Upr.txt");
        fprintf(fp,"�������� ��������� Utg_i-Utm_i;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat1,FI_t[i][0],";",X[2*M-1+i]-X[3*M-1+i],";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }

    char outputFormat2[28] = "%24.16f%s%24.16f%s";
    if((fp=fopen("Jtg.txt", "w")) == NULL)
    {
        printf("%s %s \n","���������� ������� ���� ","Jtg.txt");
    }
    else
    {
        printf("%s %s \n","�������� ������� ���� ","Jtg.txt");
        fprintf(fp,"��������� ���� �� ������� �����-�����;\n");
        for (i=0;i<M;i++) fprintf(fp,outputFormat2,FI_t[i][0],";",1000*X[M-1+i]/St,";\n");
        fprintf(fp,"%\n");
        fclose(fp);
    }
*/


}

void i_Itg_comb()
/// ��������� ������� ���������� � i_Ext_Sol, ������� ��������� ����� ���������� Upr
/// ����� ������ �������, ��������� ����������
/// ���� (Ct_0, Ct_1,...,Ct_i-1,Ct^i,Ct_i+1,...Ct_M-1)
/// � (Ct^0, Ct^1,...,Ct^i-1,Ct_i,Ct^i+1,...Ct^M-1)
/// ����� ����� 4*� ����������
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
    fprintf(fp,"��� ���������� Itg ��� Ct ���� (Ct_0, Ct_1,...,Ct_i-1,Ct^i,Ct_i+1,...Ct_M-1) � (Ct^0, Ct^1,...,Ct^i-1,Ct_i,Ct^i+1,...Ct^M-1) \n");

    for (i=0;i<M;i++) /// ������� �� ��� �������� ����� ������ ��� � ���� Itg
    {

        i_Ext_Sol[M-1+i][0] = 100; /// ����� �������� ��� � ���� ��������� ������ � Upr?
        i_Ext_Sol[M-1+i][1] = -100;


        for (w=0;w<2;w++) /// ���� �������� ��������� �������� Ct,i
        {
            /// ��������� ������ �� ���������� ����
            for (j=0;j<M;j++)
            {
                cur_Ct[j] = i_Ct_FI[j][w];
            }
//            cur_Ct[i] = i_Ct_FI[i][((w+1) mod 2)]; /// ���� w =0 �� ����� 1 � ��������
            cur_Ct[i] = i_Ct_FI[i][1-w]; /// ���� w =0 �� ����� 1 � ��������

            /// ���� �������� ��������� �������� �����
            for (k=0;k<2;k++)
            {
                cur_sig = i_sig[k];

                /// ��������� �������� ������� � ������

                i_InitStepAB(cur_Ct,cur_sig);
                step_Gauss();

                /// ���������� ����� ���� � ��� �� ������� ��  � ���������� � ������ ���-��
                /// ����� �������� � i_Ext_Sol

                if (X[M-1+i] < i_Ext_Sol[M-1+i][0])
                {
                    i_Ext_Sol[M-1+i][0] = X[M-1+i];
                }
                if (X[M-1+i] > i_Ext_Sol[M-1+i][1])
                {
                    i_Ext_Sol[M-1+i][1] = X[M-1+i];
                }

                /// �������� ���������
                for (j=0;j<M;j++)
                {
                    fprintf(fp,"%8.6f%s",X[M-1+j],"; ");
                }

                /// �������� ���������� ��
                for (j=0;j<M;j++)
                {
                    fprintf(fp,"%7.1f%s",cur_Ct[j],"; ");
                }

                ///  � �������� �����
                fprintf(fp,"%7.5f%s",cur_sig,"; ");
                fprintf(fp,"\n");
                step++;
                printf("%i%s%i\n",step," �� ",4*M);
            }

        }

    }
    fclose(fp);

}



void i_main_lim()
/// ������������ ������� ���������� ����� �������
/// ���������� �������� ����� ��� ��������� �������� ����������

/// ����������
/// i_Ct_FI[M] - ������������ �� ���� ���� ��� ������� ��
/// i_sig - ������������ �� ���� ������
/// i_Ct_N ���-�� ����� �� ��������� Ct
/// i_sig_N ���-�� ����� �� ��������� sigma
{
    unsigned int i,j,k;

    i_Init_Ct();  /// � ������ i_Ct_FI ������� ������������ Ct ��� ������� ��
                  /// �� ������ i_Ct_FI � i_sig

    i_Init_ABis(); /// ������ ������������ ������� i_ABis ��� �������� �������
/// normal i_printResult();
    step = 0;
    for (j=0;j<2;j++)
    {
        for (i=0;i<M;i++)
        {
            cur_Ct[i] = i_Ct_FI[i][j]; /// �������� � ����� ������� Ct
        }
        for (k=0;k<2;k++)
        {
            cur_sig = i_sig[k];

            /// ������ �������� ������!!!
            i_InitStepAB(cur_Ct,cur_sig);
            mygauss();

            /// ���������� ������� �������� ������� � ������
            for (i=0;i<NN;i++)
            {
                all_X[step][i] = X[i];
            }
            step++;
        }
    }

    i_Create_Sol_lim();

    i_Itg_comb();  /// ������� ������� ������ ��� Itg,i � ������� � i_Ex_Sol[M-1+i]

//    i_print_Comb_lim(); /// ����� ���������� ��� �����? ��!
}



/// ������������ ���������� �������
///    ��� ����� ����� Itx_i  X[i]
///    ��� �� ������� �����-����� Itg_i X[M-1+i]
///    ��. ��������� �� ������� �����-���� Utg_i X[2*M-1+i]
///    ��. ��������� ������� ����� Utm_i X[3*M-1+i]


/// ////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////// ������ �������� ���� ���������� - ����� /////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////



int main()
{
int i;

    LoadDataSKZ(); /// ��������� ������ �� ����� DataKSZ.csv

    /// ���������� ���������
    /// const int D1_0 = 100; max ���-�� ����� ������ � ������ ������
    /// const int D2_0 = 6;   max ���-�� �������� ������ ����� ����� ������

    /// ������ ������ � ���������� ����������
    /// int Dt �������� ���-�� ����� ����� � ������ ������
    /// int Dskz  �������� ���-�� ����� ����� � ��� (���-�� ������)
    /// float DataSKZ[D1_0][D2_0] ������, ���������� �������� ������� � �������� �������
    /// �������:
    ///  0 - ���������� ����-�� ���� ������ � ��
    ///  1 - ���������� Potential V
    ///  2 - ���������� �������� ����������� - �������� ���������
    ///  3 - ���� ���� �������� �������
    ///  4 - ���������� ����������� ���
    ///  5 - ������������ ������ ���������� ����

    Init_FI_t(); /// �������������� ������ �� �����
    /// ���������� ���������
    /// M - �������� ���-�� �� �� �����
    /// Dt - ���-�� ����� ������ � �������� ������

    /// Ct ��*�2 �������� ������������� �������� ���� ���������!!!

    /// ���������� ���������� ����������
    /// DataSKZ[][] ������ � ��������� ������� (��. LoadDataSKZ())

    /// ������ ������ � ���������� ����������
    /// L - ����� ����� ���������������� ������� �����
    /// Li = L / M ����� ������� ����� �����-�� ������ ��
    /// St = 2*pi*Rt*Li;  ������� ������� �����������, ��������������� ������ ��
    /// FI_t[i][1..3] ���������� ������ �� � ������� i
    /// FI_tg[i][1..3] ���������� ����� �� �������� �� � ������� i
    /// C_t[i]=Ct/St;  ������������� �������� ������ �����, ����-�� �� � ������� i


    Init_FI_a();
    /// �������������� ������ �� ������ (�������� ��������)
    /// ������ ��������� � ������� � ����������
    /// FI_a[i][0..2] - ����-�� ����� ������� k (i=1 ���-�� �� �� �����!!!)
    /// FI_a[i][3] - ���� ���� ��� � ������� i
    /// FI_a[i][4] �������� ����� ��, � ������ �������� ���������� i-� ���
/// ������ ���������� �� ����� � ������� ��� ������� ����� � �������� ������???

    Init_A();   /// ������ ����������� ������� � ��� ���� � �� ����� Ais ��� ��������

    step_Gauss();
    //mygauss(); /// ������ ����, �������� ����������� �������� �, ��������� ���������� � X

    printResult("resall.txt","���������� ������� �������� �������");


/// ////////////////////////////////////////////////////////////////


    i_main_lim();
//    i_main_step();
//  i_main_LR(); /// ������������ ����������

///i_main_Gauss();
//    no_simply_int();
//    simply_int();

//    i_invA();
//    i_prov_inv();


    return 0;
}
